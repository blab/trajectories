"""
Generate pairwise trajectory FASTA files from phylogenetic tree data.

Each file contains two tip sequences with their Hamming distance in headers.
"""

import argparse
import itertools
import json
import os
import random
import statistics
from tqdm import tqdm

# Reuse utilities from trajectory.py
from trajectory import sanitize_filename, parse_branches, load_sequences, find_tips


def calculate_hamming_distance(seq1, seq2):
    """
    Calculate Hamming distance between two sequences.

    Ignores positions with gaps (-) or ambiguous bases (N).
    """
    distance = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 in '-N' or c2 in '-N':
            continue
        if c1 != c2:
            distance += 1
    return distance


def find_test_clade_roots(parent_of, train_test_of):
    """
    Find test clade roots (test nodes with train parent).

    Returns list of node names that are the roots of test clades.
    """
    clade_roots = []
    for child, parent in parent_of.items():
        child_label = train_test_of.get(child, 'train')
        parent_label = train_test_of.get(parent, 'train')
        if child_label == 'test' and parent_label == 'train':
            clade_roots.append(child)
    return clade_roots


def get_clade_membership(tips, parent_of, train_test_of):
    """
    For each test tip, determine which clade root it belongs to.

    Returns dict mapping tip_name -> clade_root_name.
    """
    clade_roots = set(find_test_clade_roots(parent_of, train_test_of))
    membership = {}

    for tip in tips:
        if train_test_of.get(tip) != 'test':
            continue
        # Walk up from tip to find the clade root
        current = tip
        while current in parent_of:
            if current in clade_roots:
                membership[tip] = current
                break
            current = parent_of[current]

    return membership


def generate_pairs(items, limit=None, seed=42):
    """
    Generate all pairs or random sample if limit specified.

    For N items, there are N*(N-1)/2 possible pairs.
    """
    items = list(items)
    all_pairs = list(itertools.combinations(items, 2))

    if limit is None or limit >= len(all_pairs):
        return all_pairs

    random.seed(seed)
    return random.sample(all_pairs, limit)


def generate_test_pairs_by_clade(clade_membership, limit=None, seed=42):
    """
    Generate test pairs only within same clade.

    Returns list of (tip1, tip2) tuples.
    """
    # Group tips by clade
    clades = {}
    for tip, clade_root in clade_membership.items():
        if clade_root not in clades:
            clades[clade_root] = []
        clades[clade_root].append(tip)

    # Generate pairs within each clade
    all_pairs = []
    for clade_tips in clades.values():
        if len(clade_tips) >= 2:
            all_pairs.extend(itertools.combinations(clade_tips, 2))

    if limit is None or limit >= len(all_pairs):
        return all_pairs

    random.seed(seed)
    return random.sample(all_pairs, limit)


def write_pairwise_fasta(tip1, tip2, sequences, output_path):
    """
    Write pairwise FASTA file.

    First sequence gets |0, second gets |{hamming_distance}.
    Returns the Hamming distance, or None if sequences missing.
    """
    seq1 = sequences.get(tip1, '')
    seq2 = sequences.get(tip2, '')

    if not seq1 or not seq2:
        return None

    hamming = calculate_hamming_distance(seq1, seq2)

    lines = []
    # First sequence with |0
    lines.append(f">{tip1}|0\n")
    for j in range(0, len(seq1), 60):
        lines.append(seq1[j:j+60] + '\n')

    # Second sequence with |{hamming}
    lines.append(f">{tip2}|{hamming}\n")
    for j in range(0, len(seq2), 60):
        lines.append(seq2[j:j+60] + '\n')

    with open(output_path, 'w') as f:
        f.write(''.join(lines))

    return hamming


def main():
    parser = argparse.ArgumentParser(
        description="Generate pairwise trajectory FASTA files."
    )
    parser.add_argument("--branches", required=True,
                        help="Path to branches.tsv file")
    parser.add_argument("--alignment", required=True,
                        help="Path to alignment FASTA file")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for pairwise files")
    parser.add_argument("--train-limit", type=int, default=None,
                        help="Max training pairs (default: all)")
    parser.add_argument("--test-limit", type=int, default=None,
                        help="Max test pairs (default: all)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for sampling")
    parser.add_argument("--summary", help="Path to summary JSON file")
    parser.add_argument("--dataset", help="Dataset name for summary")
    parser.add_argument("--url", help="Dataset URL for summary")

    args = parser.parse_args()

    # Parse inputs
    print("Loading branches and alignment...")
    parent_of, hamming_of, train_test_of = parse_branches(args.branches)
    sequences = load_sequences(args.alignment)
    tips = find_tips(parent_of)

    # Separate train/test tips
    train_tips = [t for t in tips if train_test_of.get(t, 'train') == 'train']
    test_tips = [t for t in tips if train_test_of.get(t) == 'test']

    print(f"Found {len(train_tips)} train tips and {len(test_tips)} test tips")

    # Create output directories
    train_dir = os.path.join(args.output_dir, 'pairwise-train')
    test_dir = os.path.join(args.output_dir, 'pairwise-test')
    os.makedirs(train_dir, exist_ok=True)
    os.makedirs(test_dir, exist_ok=True)

    # Generate train pairs
    train_pairs = generate_pairs(train_tips, args.train_limit, args.seed)
    train_distances = []

    print(f"Writing {len(train_pairs)} train pairs...")
    for tip1, tip2 in tqdm(train_pairs, desc="Train pairs"):
        safe1 = sanitize_filename(tip1)
        safe2 = sanitize_filename(tip2)
        output_path = os.path.join(train_dir, f"{safe1}__{safe2}.fasta")
        dist = write_pairwise_fasta(tip1, tip2, sequences, output_path)
        if dist is not None:
            train_distances.append(dist)

    # Generate test pairs (within clades only)
    clade_membership = get_clade_membership(test_tips, parent_of, train_test_of)
    test_pairs = generate_test_pairs_by_clade(clade_membership, args.test_limit, args.seed)
    test_distances = []

    # Count clades
    unique_clades = set(clade_membership.values())
    print(f"Found {len(unique_clades)} test clades")

    print(f"Writing {len(test_pairs)} test pairs...")
    for tip1, tip2 in tqdm(test_pairs, desc="Test pairs"):
        safe1 = sanitize_filename(tip1)
        safe2 = sanitize_filename(tip2)
        output_path = os.path.join(test_dir, f"{safe1}__{safe2}.fasta")
        dist = write_pairwise_fasta(tip1, tip2, sequences, output_path)
        if dist is not None:
            test_distances.append(dist)

    print(f"Done! Wrote {len(train_distances)} train and {len(test_distances)} test pairs")

    # Update summary if requested
    if args.summary and args.dataset:
        # Load existing summary or create new
        if os.path.exists(args.summary):
            with open(args.summary, 'r') as f:
                summary = json.load(f)
        else:
            summary = {}

        # Ensure dataset entry exists
        if args.dataset not in summary:
            summary[args.dataset] = {}

        # Add pairwise statistics
        summary[args.dataset]['pairwise_train_pairs'] = len(train_distances)
        summary[args.dataset]['pairwise_test_pairs'] = len(test_distances)
        summary[args.dataset]['pairwise_test_clades'] = len(unique_clades)

        if train_distances:
            summary[args.dataset]['pairwise_train_hamming'] = {
                'min': min(train_distances),
                'max': max(train_distances),
                'mean': round(statistics.mean(train_distances), 2)
            }

        if test_distances:
            summary[args.dataset]['pairwise_test_hamming'] = {
                'min': min(test_distances),
                'max': max(test_distances),
                'mean': round(statistics.mean(test_distances), 2)
            }

        with open(args.summary, 'w') as f:
            json.dump(summary, f, indent=2)

        print(f"Updated summary at {args.summary}")


if __name__ == "__main__":
    main()
