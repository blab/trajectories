"""
Extract per-tip trajectory FASTA files from phylogenetic tree data.

Each trajectory contains sequences from root to tip with branch Hamming distances.
"""

import argparse
import subprocess
import sys
sys.setrecursionlimit(100000)
import csv
import json
import os
import re
import statistics
from Bio import SeqIO
from tqdm import tqdm
import zstandard as zstd

from shard_writer import ShardWriter

def get_git_commit():
    """Return short git commit hash, with '-dirty' suffix if working tree has changes."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        try:
            subprocess.check_call(
                ["git", "diff", "--quiet", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            commit += "-dirty"
        return commit
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


# Pre-compile regex for filename sanitization
UNSAFE_CHARS_RE = re.compile(r'[/\\:*?"<>| ]')


def sanitize_filename(name):
    """
    Sanitize tip name for use as filename.

    Strips 'hCoV-19' prefix and removes unsafe filesystem characters.
    """
    # Strip hCoV-19 prefix
    if name.startswith("hCoV-19"):
        name = name[7:]  # len("hCoV-19") = 7

    return UNSAFE_CHARS_RE.sub('', name)


def parse_branches(branches_path):
    """
    Parse branches.tsv to build parent-child relationships, hamming distances, and train/test labels.

    Returns:
        parent_of: dict mapping child -> parent
        hamming_of: dict mapping (parent, child) -> hamming distance
        train_test_of: dict mapping node -> train/test label
    """
    parent_of = {}
    hamming_of = {}
    train_test_of = {}

    with open(branches_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            parent = row['parent']
            child = row['child']
            hamming = row['hamming']
            train_test = row.get('train_test', '')

            parent_of[child] = parent
            # Handle missing hamming values (marked as '?')
            if hamming != '?':
                hamming_of[(parent, child)] = int(hamming)
            else:
                hamming_of[(parent, child)] = 0

            if train_test:
                train_test_of[child] = train_test

    return parent_of, hamming_of, train_test_of


def load_sequences(alignment_path):
    """Load sequences from FASTA file into a dictionary."""
    sequences = {}
    for record in SeqIO.parse(alignment_path, 'fasta'):
        sequences[record.id] = str(record.seq)
    return sequences


def find_tips(parent_of):
    """
    Find tip nodes (nodes that are children but never parents).
    """
    children = set(parent_of.keys())
    parents = set(parent_of.values())
    tips = children - parents
    return tips


def get_path_to_root(tip, parent_of):
    """
    Trace path from tip back to root.

    Returns list of nodes from tip to root.
    """
    path = [tip]
    current = tip
    while current in parent_of:
        current = parent_of[current]
        path.append(current)
    return path


def find_test_boundary(path, train_test_of):
    """
    Find the index of the first test node in a root-to-tip path.

    Returns the index where test nodes begin, or None if all nodes are train.
    """
    for i, node in enumerate(path):
        if train_test_of.get(node) == 'test':
            return i
    return None


def format_sequence(seq, line_width=60):
    """Format sequence with line breaks, optimized for speed."""
    return '\n'.join(seq[i:i+line_width] for i in range(0, len(seq), line_width))


def hamming_distance(seq1, seq2):
    """Hamming distance ignoring positions where either sequence has a gap or N."""
    d = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 in '-N' or c2 in '-N':
            continue
        if c1 != c2:
            d += 1
    return d


def trim_path_gaps(path, sequences):
    """
    Per-trajectory gap trimming.

    Drop columns where every node on this root-to-tip path has '-' or 'N'.
    Rationale: columns that are structural gaps throughout this lineage carry
    no signal for this trajectory — they only exist because MAFFT aligned
    against OTHER lineages' insertions. Removing them preserves every real
    insertion/deletion event (any column where at least one node has a base
    is kept) while shortening trajectories 2-5x in practice.
    Hamming distances are unchanged (all-'-' columns contribute 0 to Hamming).

    Returns a new dict with trimmed sequences for just the nodes on this path.
    """
    path_seqs = [sequences[n] for n in path if n in sequences and sequences[n]]
    if not path_seqs:
        return sequences
    L = len(path_seqs[0])
    # For each column, keep it if ANY node on this path has a real base
    keep = bytearray(L)
    for s in path_seqs:
        s_str = str(s)
        for i, c in enumerate(s_str):
            if c not in '-N':
                keep[i] = 1
    keep_idx = [i for i in range(L) if keep[i]]
    if len(keep_idx) == L:
        return sequences  # nothing to trim
    trimmed = {}
    for node, s in sequences.items():
        if not s:
            trimmed[node] = s
            continue
        s_str = str(s)
        trimmed[node] = type(s)(''.join(s_str[i] for i in keep_idx))
    return trimmed


def build_trajectory_content(path, sequences, hamming_of):
    """
    Build trajectory FASTA content for a single tip.

    Path should be in root-to-tip order.
    Each header includes the branch Hamming distance from the previous emitted
    node and the direct Hamming distance from the start node:
    ``>{node}|{branch_distance}|{direct_distance}``
    Skips nodes with zero branch distance. If the tip is skipped, the last
    emitted node is relabeled with the tip's name so the trajectory always
    ends with the tip.

    Alignment columns that are gap ('-'/'N') in every node on this path are
    dropped (per-trajectory trim), so trajectories only contain positions
    biologically relevant to this lineage.

    Returns:
        tuple: (content, tip_distance, path_depth) where content is the FASTA string,
               tip_distance is cumulative Hamming distance from root, and path_depth is
               number of frames written.
    """
    # Per-trajectory gap trim: drop columns that are always gap on this path
    sequences = trim_path_gaps(path, sequences)

    cumulative_distance = 0
    tip_node = path[-1] if path else None
    frames_written = 0
    last_emitted_seq = None
    # First emitted sequence, used to compute direct Hamming distance
    start_seq = None

    # Build content
    parts = []
    for i, node in enumerate(path):
        # Calculate cumulative distance
        if i > 0:
            parent = path[i - 1]
            branch_dist = hamming_of.get((parent, node), 0)
            cumulative_distance += branch_dist

            # Skip nodes with no distance increase
            if branch_dist == 0:
                continue

        # Get sequence
        seq = sequences.get(node)
        if not seq:
            continue

        # Compute header distance from the last emitted sequence directly,
        # rather than using the tree-edge Hamming. Gap-ignoring Hamming is
        # not additive along tree paths when gap patterns change at skipped
        # intermediate nodes.
        if last_emitted_seq is not None:
            emitted_dist = hamming_distance(last_emitted_seq, seq)
        else:
            emitted_dist = 0

        # Set start_seq on the first emitted sequence.
        if start_seq is None:
            start_seq = seq

        # Direct Hamming distance from the start node
        direct_dist = hamming_distance(start_seq, seq)

        # Add FASTA entry
        parts.append(f">{node}|{emitted_dist}|{direct_dist}\n{format_sequence(seq)}\n")
        last_emitted_seq = seq
        frames_written += 1

    # Ensure the trajectory ends with the tip name. If the tip was
    # skipped (zero branch distance), relabel the last emitted frame.
    if parts and tip_node is not None:
        last_header = parts[-1].split('\n', 1)[0]
        last_name = last_header[1:].split('|')[0]
        if last_name != tip_node:
            parts[-1] = parts[-1].replace(f">{last_name}|", f">{tip_node}|", 1)

    content = ''.join(parts)
    return content, cumulative_distance, frames_written


def write_trajectory(path, sequences, hamming_of, output_path, compress=False, compressor=None):
    """
    Write trajectory FASTA file for a single tip.

    Path should be in root-to-tip order.
    Each header includes the branch Hamming distance from the previous emitted node.
    Skips nodes with zero branch distance. If the tip is skipped, the last
    emitted node is relabeled with the tip's name.

    Returns:
        tuple: (tip_distance, path_depth) where tip_distance is cumulative
               Hamming distance from root and path_depth is number of frames written.
    """
    content, cumulative_distance, frames_written = build_trajectory_content(
        path, sequences, hamming_of
    )

    # Write to file (compressed or plain)
    if compress:
        cctx = compressor or zstd.ZstdCompressor()
        with open(output_path, 'wb') as f:
            f.write(cctx.compress(content.encode('utf-8')))
    else:
        with open(output_path, 'w') as f:
            f.write(content)

    return cumulative_distance, frames_written


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-tip trajectory FASTA files from phylogenetic tree data."
    )
    parser.add_argument(
        "--branches", required=True,
        help="Path to branches.tsv file"
    )
    parser.add_argument(
        "--alignment", required=True,
        help="Path to alignment.fasta file"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for trajectory files"
    )
    parser.add_argument(
        "--shard-size", type=int, default=10000,
        help="Number of trajectories per shard (default: 10000)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for shuffling (default: 42)"
    )
    parser.add_argument(
        "--summary",
        help="Path to write summary statistics JSON file"
    )
    parser.add_argument(
        "--dataset",
        help="Dataset name (used as key in summary JSON)"
    )
    parser.add_argument(
        "--url",
        help="Dataset source URL (included in summary JSON)"
    )
    args = parser.parse_args()

    # Parse input files
    print("Loading branches...")
    parent_of, hamming_of, train_test_of = parse_branches(args.branches)

    print("Loading sequences...")
    sequences = load_sequences(args.alignment)

    # Find tips
    tips = find_tips(parent_of)
    print(f"Found {len(tips)} tips")

    # Get sequence length from first sequence
    seq_length = len(next(iter(sequences.values()))) if sequences else 0

    # Compute branch statistics
    branch_distances = list(hamming_of.values())
    zero_distance_branches = sum(1 for d in branch_distances if d == 0)

    # Check if we have train/test labels
    has_train_test = bool(train_test_of)

    # Process each tip and collect statistics
    tip_distances = []
    path_depths = []
    train_tips = 0
    test_tips = 0

    print(f"Writing trajectory shards (shard_size={args.shard_size})...")

    with ShardWriter(args.output_dir, "forwards-train", args.shard_size, shuffle=True, seed=args.seed) as train_writer, \
         ShardWriter(args.output_dir, "forwards-test", args.shard_size, shuffle=True, seed=args.seed) as test_writer:

        for tip in tqdm(tips, desc="Processing tips"):
            # Get path from tip to root, then reverse to root-to-tip
            path = get_path_to_root(tip, parent_of)
            path.reverse()

            # Sanitize filename
            safe_name = sanitize_filename(tip)
            filename = f"{safe_name}.fasta"

            # Determine which writer to use and potentially truncate path for test tips
            is_test = False
            if has_train_test:
                tip_label = train_test_of.get(tip, 'train')
                if tip_label == 'test':
                    # Find where test clade begins and truncate path
                    boundary_idx = find_test_boundary(path, train_test_of)
                    if boundary_idx is not None:
                        path = path[boundary_idx:]  # Start from first test node
                    writer = test_writer
                    is_test = True
                else:
                    writer = train_writer
            else:
                writer = train_writer

            # Build trajectory content and collect stats
            content, tip_dist, path_depth = build_trajectory_content(
                path, sequences, hamming_of
            )

            # Skip trajectories with fewer than 2 sequences
            if path_depth < 2:
                continue

            # Count tips after filtering
            if is_test:
                test_tips += 1
            else:
                train_tips += 1

            writer.add(filename, content)
            tip_distances.append(tip_dist)
            path_depths.append(path_depth)

    # Report results
    train_shards, train_files, train_bytes = train_writer.stats
    test_shards, test_files, test_bytes = test_writer.stats

    print(f"Done! Wrote {train_tips} train trajectories ({train_shards} shards, {train_bytes / 1024 / 1024:.1f} MB)")
    print(f"      Wrote {test_tips} test trajectories ({test_shards} shards, {test_bytes / 1024 / 1024:.1f} MB)")

    # Write summary statistics if requested
    if args.summary and args.dataset:
        # Build stats for this dataset
        dataset_summary = {
            "git_commit": get_git_commit(),
            "url": args.url,
            "num_tips": len(tips),
            "num_nodes": len(sequences),
            "sequence_length": seq_length,
            "hamming_from_root": {
                "min": min(tip_distances),
                "max": max(tip_distances),
                "mean": round(statistics.mean(tip_distances), 2)
            },
            "path_depth": {
                "min": min(path_depths),
                "max": max(path_depths),
                "mean": round(statistics.mean(path_depths), 2)
            },
            "total_branches": len(branch_distances),
            "zero_distance_branches": zero_distance_branches,
            "per_branch_hamming": {
                "min": min(branch_distances) if branch_distances else 0,
                "max": max(branch_distances) if branch_distances else 0,
                "mean": round(statistics.mean(branch_distances), 2) if branch_distances else 0
            }
        }

        # Add train/test counts if available
        if has_train_test:
            dataset_summary["train_tips"] = train_tips
            dataset_summary["test_tips"] = test_tips

        # Load existing summary or start fresh
        if os.path.exists(args.summary):
            with open(args.summary, 'r') as f:
                all_summaries = json.load(f)
        else:
            all_summaries = {}

        # Update this dataset's entry
        all_summaries[args.dataset] = dataset_summary

        # Write back
        with open(args.summary, 'w') as f:
            json.dump(all_summaries, f, indent=2)
        print(f"Wrote summary for '{args.dataset}' to {args.summary}")


if __name__ == "__main__":
    main()
