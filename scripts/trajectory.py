"""
Extract per-tip trajectory FASTA files from phylogenetic tree data.

Each trajectory contains sequences from root to tip with cumulative Hamming distances.
"""

import argparse
import csv
import os
import re
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
    Parse branches.tsv to build parent-child relationships and hamming distances.

    Returns:
        parent_of: dict mapping child -> parent
        hamming_of: dict mapping (parent, child) -> hamming distance
    """
    parent_of = {}
    hamming_of = {}

    with open(branches_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            parent = row['parent']
            child = row['child']
            hamming = row['hamming']

            parent_of[child] = parent
            # Handle missing hamming values (marked as '?')
            if hamming != '?':
                hamming_of[(parent, child)] = int(hamming)
            else:
                hamming_of[(parent, child)] = 0

    return parent_of, hamming_of


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


def write_trajectory(path, sequences, hamming_of, output_path, compress=False):
    """
    Write trajectory FASTA file for a single tip.

    Path should be in root-to-tip order.
    Each header includes cumulative Hamming distance from root.
    Skips intermediate nodes with zero branch distance.
    """
    cumulative_distance = 0
    last_node = path[-1] if path else None  # tip node

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

    content = ''.join(lines)

    # Write to file (compressed or plain)
    if compress:
        cctx = zstd.ZstdCompressor()
        with open(output_path, 'wb') as f:
            f.write(cctx.compress(content.encode('utf-8')))
    else:
        with open(output_path, 'w') as f:
            f.write(content)


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
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse input files
    print("Loading branches...")
    parent_of, hamming_of = parse_branches(args.branches)

    print("Loading sequences...")
    sequences = load_sequences(args.alignment)

    # Find tips
    tips = find_tips(parent_of)
    print(f"Found {len(tips)} tips")

    # Process each tip
    ext = ".fasta.zst" if args.compress else ".fasta"
    print(f"Writing trajectory files{' (compressed)' if args.compress else ''}...")
    for tip in tqdm(tips, desc="Processing tips"):
        # Get path from tip to root, then reverse to root-to-tip
        path = get_path_to_root(tip, parent_of)
        path.reverse()

        # Sanitize filename
        safe_name = sanitize_filename(tip)
        output_path = os.path.join(args.output_dir, f"{safe_name}{ext}")

        # Write trajectory
        write_trajectory(path, sequences, hamming_of, output_path, compress=args.compress)

    print(f"Done! Wrote {len(tips)} trajectory files to {args.output_dir}")


if __name__ == "__main__":
    main()
