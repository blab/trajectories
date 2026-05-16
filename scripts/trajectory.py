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
import numpy as np
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


def load_div_map(auspice_path):
    """
    Load the auspice JSON and build a {node_name: div} map from node_attrs.div.

    ``div`` is the auspice tree divergence (cumulative subs/site augur
    computed; root ~0). Walks the whole ``tree`` object. The parsed JSON is
    freed before returning so only the compact float map is retained.
    """
    with open(auspice_path) as f:
        auspice = json.load(f)

    div_map = {}
    stack = [auspice.get("tree")]
    while stack:
        node = stack.pop()
        if not node:
            continue
        name = node.get("name")
        node_attrs = node.get("node_attrs") or {}
        div = node_attrs.get("div")
        if name is not None and div is not None:
            div_map[name] = float(div)
        children = node.get("children")
        if children:
            stack.extend(children)

    return div_map


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
    a = np.frombuffer(str(seq1).encode("ascii"), dtype=np.uint8)
    b = np.frombuffer(str(seq2).encode("ascii"), dtype=np.uint8)
    valid = (a != ord("-")) & (a != ord("N")) & (b != ord("-")) & (b != ord("N"))
    return int(((a != b) & valid).sum())


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
    Uses numpy for fast column reduction.
    """
    path_seqs = [str(sequences[n]) for n in path if n in sequences and sequences[n]]
    if not path_seqs:
        return sequences
    L = len(path_seqs[0])

    # OR together the "has-a-base" bitmaps across all nodes on the path.
    # Keep columns where at least one node has a real base (not '-' or 'N').
    path_mask = np.zeros(L, dtype=bool)
    for s in path_seqs:
        arr = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        path_mask |= (arr != ord("-")) & (arr != ord("N"))

    if path_mask.all():
        return sequences  # nothing to trim

    keep_idx = np.nonzero(path_mask)[0]
    trimmed = {}
    for node in path:
        s = sequences.get(node)
        if not s:
            continue
        arr = np.frombuffer(str(s).encode("ascii"), dtype=np.uint8)
        new_str = arr[keep_idx].tobytes().decode("ascii")
        trimmed[node] = type(s)(new_str)
    return trimmed


def build_trajectory_content(path, sequences, hamming_of,
                             max_divergence=None, min_nodes=3,
                             div_map=None):
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

    Divergence-based basal truncation (when ``max_divergence`` is set):
        The basal end of the trajectory is truncated so the most-basal kept
        node ``A_k`` satisfies ``divergence(A_k -> tip) <= max_divergence``,
        where ``divergence(ancestor -> tip) = div[tip] - div[ancestor]`` and
        ``div`` is the auspice tree divergence (``node_attrs.div``, cumulative
        subs/site augur computed, root ~0) passed in via ``div_map``.
        Walking from the tip toward the root, ancestors are dropped once
        ``div[tip] - div[ancestor]`` would exceed the threshold. ``min_nodes``
        is a DROP FILTER, not an extend guard: if the truncated trajectory has
        fewer than ``min_nodes`` emitted frames the whole trajectory is
        dropped (the path is NEVER extended basally past the cutoff). When
        ``max_divergence`` is None, the full root-to-tip path is emitted
        (fully backward-compatible) and no drop filter applies.

    Returns:
        tuple: (content, tip_distance, path_depth, trimmed_seq_length) where
               content is the FASTA string, tip_distance is the cumulative
               tree-edge Hamming distance from the (possibly truncated) origin,
               path_depth is the number of frames written, and
               trimmed_seq_length is the per-trajectory trimmed column count.
               When ``max_divergence`` is set and the divergence-truncated
               trajectory has fewer than ``min_nodes`` frames, the trajectory
               is dropped (content ""). path_depth is -1 when the drop is
               attributable to truncation -- the untruncated trajectory would
               have been emitted (>= 2 frames and >= min_nodes) -- and 0 when
               the trajectory was too short to emit regardless of truncation.
    """
    # Per-trajectory gap trim: drop columns that are always gap on this path
    sequences = trim_path_gaps(path, sequences)

    tip_node = path[-1] if path else None

    # First pass: collect every emitted frame as
    # (node, emitted_dist, branch_dist, seq).
    #   emitted_dist - gap-ignoring Hamming from the previous emitted node
    #                  (the per-branch substitution count in the header)
    #   branch_dist  - the tree-edge Hamming from hamming_of; cumulative
    #                  tip_distance is the sum of these (unchanged from the
    #                  original behaviour) so the summary JSON is stable.
    frames = []
    last_emitted_seq = None
    for i, node in enumerate(path):
        branch_dist = 0
        if i > 0:
            parent = path[i - 1]
            branch_dist = hamming_of.get((parent, node), 0)
            # Skip nodes with no tree-edge distance increase
            if branch_dist == 0:
                continue

        seq = sequences.get(node)
        if not seq:
            continue

        # Header distance computed directly from the last emitted sequence;
        # gap-ignoring Hamming is not additive across skipped intermediates.
        if last_emitted_seq is not None:
            emitted_dist = hamming_distance(last_emitted_seq, seq)
        else:
            emitted_dist = 0

        frames.append((node, emitted_dist, branch_dist, seq))
        last_emitted_seq = seq

    # Divergence-based basal truncation. divergence(ancestor -> tip) is the
    # auspice tree divergence difference div[tip] - div[ancestor]; div is
    # cumulative-from-root subs/site so the difference is path divergence.
    if max_divergence is not None and frames:
        div_map = div_map or {}
        div_tip = div_map.get(tip_node)
        # Number of emittable frames BEFORE truncation. A tip is only counted
        # as a divergence drop if it would otherwise have been emitted -- the
        # untruncated trajectory had >= 2 frames (the emit threshold). Tips
        # already below 2 frames are skipped regardless of --max-divergence,
        # so truncation is not their cause.
        n_untruncated = len(frames)
        # Walk from the tip toward the root; keep the most-basal frame whose
        # node is still within max_divergence of the tip.
        # start_idx = index of the most-basal kept frame.
        start_idx = len(frames) - 1
        if div_tip is not None:
            for j in range(len(frames) - 1, -1, -1):
                node_j = frames[j][0]
                div_j = div_map.get(node_j)
                if div_j is None:
                    # No divergence for this node: stop, do not include it.
                    break
                if (div_tip - div_j) > max_divergence:
                    # Including frame j exceeds the threshold; stop above it.
                    break
                start_idx = j
        frames = frames[start_idx:]
        # min_nodes is a DROP FILTER: if fewer than min_nodes frames remain
        # after truncation, drop the whole trajectory (never extend basally).
        # path_depth is returned as -1 to mark a divergence drop -- but only
        # when the untruncated trajectory would have been emitted (>= 2
        # frames); otherwise path_depth 0 lets the normal "< 2 frames" skip
        # handle it (truncation is not the cause of those drops).
        if len(frames) < min_nodes:
            if n_untruncated >= 2:
                return "", 0, -1, 0
            return "", 0, 0, 0
        # The most-basal kept frame is the new origin: zero its branch dists
        # so the |0|0 header and cumulative tip_distance start fresh there.
        node0, _, _, seq0 = frames[0]
        frames[0] = (node0, 0, 0, seq0)

    # Second pass: emit FASTA. cumulative_distance is the tree-edge Hamming
    # sum (matches the original, unchanged when max_divergence is None).
    cumulative_distance = 0
    frames_written = 0
    start_seq = frames[0][3] if frames else None
    parts = []
    for node, emitted_dist, branch_dist, seq in frames:
        cumulative_distance += branch_dist
        direct_dist = hamming_distance(start_seq, seq)
        parts.append(f">{node}|{emitted_dist}|{direct_dist}\n{format_sequence(seq)}\n")
        frames_written += 1

    # Ensure the trajectory ends with the tip name. If the tip was
    # skipped (zero branch distance), relabel the last emitted frame.
    if parts and tip_node is not None:
        last_header = parts[-1].split('\n', 1)[0]
        last_name = last_header[1:].split('|')[0]
        if last_name != tip_node:
            parts[-1] = parts[-1].replace(f">{last_name}|", f">{tip_node}|", 1)

    content = ''.join(parts)
    trimmed_seq_length = len(start_seq) if start_seq else 0
    return content, cumulative_distance, frames_written, trimmed_seq_length


def write_trajectory(path, sequences, hamming_of, output_path, compress=False,
                     compressor=None, max_divergence=None, min_nodes=3,
                     div_map=None):
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
    content, cumulative_distance, frames_written, _ = build_trajectory_content(
        path, sequences, hamming_of, max_divergence=max_divergence,
        min_nodes=min_nodes, div_map=div_map
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
    parser.add_argument(
        "--auspice",
        help="Path to the auspice JSON (required when --max-divergence is "
             "set). Used to read node_attrs.div, the auspice tree divergence."
    )
    parser.add_argument(
        "--max-divergence", type=float, default=None,
        help="If set, truncate each trajectory's basal end so the most-basal "
             "kept node is within this auspice tree divergence "
             "(node_attrs.div, cumulative subs/site) of the tip. "
             "Unset = full root-to-tip path. Requires --auspice."
    )
    parser.add_argument(
        "--min-nodes", type=int, default=3,
        help="Minimum number of nodes per trajectory. When --max-divergence "
             "is set this is a DROP FILTER: trajectories with fewer than "
             "this many nodes after truncation are dropped entirely (default: 3)"
    )
    args = parser.parse_args()

    if args.max_divergence is not None and not args.auspice:
        parser.error("--auspice is required when --max-divergence is set")

    # Parse input files
    print("Loading branches...")
    parent_of, hamming_of, train_test_of = parse_branches(args.branches)

    print("Loading sequences...")
    sequences = load_sequences(args.alignment)

    # Load auspice divergence map for divergence-based truncation.
    div_map = None
    if args.max_divergence is not None:
        print(f"Loading auspice divergence map from {args.auspice}...")
        div_map = load_div_map(args.auspice)
        print(f"Loaded div for {len(div_map)} nodes")

    # Find tips
    tips = find_tips(parent_of)
    print(f"Found {len(tips)} tips")

    # Get alignment length from first sequence (before per-trajectory trimming)
    alignment_length = len(next(iter(sequences.values()))) if sequences else 0

    # Compute branch statistics
    branch_distances = list(hamming_of.values())
    zero_distance_branches = sum(1 for d in branch_distances if d == 0)

    # Check if we have train/test labels
    has_train_test = bool(train_test_of)

    # Process each tip and collect statistics
    tip_distances = []
    path_depths = []
    trimmed_lengths = []
    train_tips = 0
    test_tips = 0
    dropped_tips = 0

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
            content, tip_dist, path_depth, trimmed_len = build_trajectory_content(
                path, sequences, hamming_of,
                max_divergence=args.max_divergence, min_nodes=args.min_nodes,
                div_map=div_map
            )

            # Divergence drop filter: when --max-divergence is set, a tip
            # whose truncated trajectory has < min_nodes frames is dropped
            # (build_trajectory_content signals this with path_depth == -1).
            if path_depth == -1:
                dropped_tips += 1
                continue

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
            trimmed_lengths.append(trimmed_len)

    # Report results
    train_shards, train_files, train_bytes = train_writer.stats
    test_shards, test_files, test_bytes = test_writer.stats

    print(f"Done! Wrote {train_tips} train trajectories ({train_shards} shards, {train_bytes / 1024 / 1024:.1f} MB)")
    print(f"      Wrote {test_tips} test trajectories ({test_shards} shards, {test_bytes / 1024 / 1024:.1f} MB)")
    if args.max_divergence is not None:
        print(f"      Dropped {dropped_tips} tips "
              f"(< {args.min_nodes} nodes within divergence {args.max_divergence})")

    # Write summary statistics if requested
    if args.summary and args.dataset:
        # min/max/mean over a list, tolerant of an empty list (e.g. when the
        # divergence drop filter dropped every tip for this marker).
        def _stats(values):
            if not values:
                return {"min": 0, "max": 0, "mean": 0}
            return {
                "min": min(values),
                "max": max(values),
                "mean": round(statistics.mean(values), 2),
            }

        # Build stats for this dataset
        dataset_summary = {
            "git_commit": get_git_commit(),
            "url": args.url,
            "num_tips": len(tips),
            "num_nodes": len(sequences),
            "alignment_length": alignment_length,
            "trimmed_length": _stats(trimmed_lengths),
            "hamming_from_root": _stats(tip_distances),
            "path_depth": _stats(path_depths),
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

        # Record divergence-truncation drop count when active
        if args.max_divergence is not None:
            dataset_summary["max_divergence"] = args.max_divergence
            dataset_summary["min_nodes"] = args.min_nodes
            dataset_summary["dropped_tips"] = dropped_tips

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
