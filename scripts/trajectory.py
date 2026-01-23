"""
Extract per-tip trajectory FASTA files from phylogenetic tree data.

Each trajectory contains sequences from root to tip with cumulative Hamming distances.
"""

import argparse
import csv
import json
import os
import re
import statistics
from Bio import SeqIO
from tqdm import tqdm
import zstandard as zstd


def sanitize_filename(name):
    """
    Sanitize tip name for use as filename.

    Strips 'hCoV-19' prefix and removes unsafe filesystem characters.
    """
    # Strip hCoV-19 prefix
    if name.startswith("hCoV-19"):
        name = name[7:]  # len("hCoV-19") = 7

    # Remove unsafe characters: / \ : * ? " < > | and space
    unsafe_chars = r'[/\\:*?"<>| ]'
    return re.sub(unsafe_chars, '', name)


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


def write_trajectory(path, sequences, hamming_of, output_path, compress=False):
    """
    Write trajectory FASTA file for a single tip.

    Path should be in root-to-tip order.
    Each header includes cumulative Hamming distance from root.
    Skips intermediate nodes with zero branch distance.

    Returns:
        tuple: (tip_distance, path_depth) where tip_distance is cumulative
               Hamming distance and path_depth is number of frames written.
    """
    cumulative_distance = 0
    last_node = path[-1] if path else None  # tip node
    frames_written = 0

    # Build content
    lines = []
    for i, node in enumerate(path):
        # Calculate cumulative distance
        if i > 0:
            parent = path[i - 1]
            branch_dist = hamming_of.get((parent, node), 0)
            cumulative_distance += branch_dist

            # Skip if no distance increase (except for tip)
            if branch_dist == 0 and node != last_node:
                continue

        # Get sequence
        seq = sequences.get(node, '')
        if not seq:
            continue

        # Add FASTA entry
        lines.append(f">{node}|{cumulative_distance}\n")
        # Add sequence in 60-char lines
        for j in range(0, len(seq), 60):
            lines.append(seq[j:j+60] + '\n')
        frames_written += 1

    content = ''.join(lines)

    # Write to file (compressed or plain)
    if compress:
        cctx = zstd.ZstdCompressor()
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
        "--compress", action="store_true",
        help="Compress output files with zstd (.fasta.zst)"
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

    # Create output subdirectories if we have train/test labels
    if has_train_test:
        train_dir = os.path.join(args.output_dir, 'forwards-train')
        test_dir = os.path.join(args.output_dir, 'forwards-test')
        os.makedirs(train_dir, exist_ok=True)
        os.makedirs(test_dir, exist_ok=True)
    else:
        os.makedirs(args.output_dir, exist_ok=True)

    # Process each tip and collect statistics
    tip_distances = []
    path_depths = []
    train_tips = 0
    test_tips = 0
    ext = ".fasta.zst" if args.compress else ".fasta"
    print(f"Writing trajectory files{' (compressed)' if args.compress else ''}...")
    for tip in tqdm(tips, desc="Processing tips"):
        # Get path from tip to root, then reverse to root-to-tip
        path = get_path_to_root(tip, parent_of)
        path.reverse()

        # Sanitize filename
        safe_name = sanitize_filename(tip)

        # Determine output path and potentially truncate path for test tips
        if has_train_test:
            tip_label = train_test_of.get(tip, 'train')
            if tip_label == 'test':
                # Find where test clade begins and truncate path
                boundary_idx = find_test_boundary(path, train_test_of)
                if boundary_idx is not None:
                    path = path[boundary_idx:]  # Start from first test node
                output_path = os.path.join(test_dir, f"{safe_name}{ext}")
                test_tips += 1
            else:
                output_path = os.path.join(train_dir, f"{safe_name}{ext}")
                train_tips += 1
        else:
            output_path = os.path.join(args.output_dir, f"{safe_name}{ext}")

        # Write trajectory and collect stats
        tip_dist, path_depth = write_trajectory(
            path, sequences, hamming_of, output_path, compress=args.compress
        )
        tip_distances.append(tip_dist)
        path_depths.append(path_depth)

    if has_train_test:
        print(f"Done! Wrote {train_tips} train and {test_tips} test trajectory files to {args.output_dir}")
    else:
        print(f"Done! Wrote {len(tips)} trajectory files to {args.output_dir}")

    # Write summary statistics if requested
    if args.summary and args.dataset:
        # Build stats for this dataset
        dataset_summary = {
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
