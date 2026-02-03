#!/usr/bin/env python3
"""
Train/test split for phylogenetic branches.

This script reads branches_raw.tsv (which contains parent-child relationships)
and marks entire clades as "test" data. Test clades are selected by randomly
choosing tip nodes, walking back a specified number of mutations, and marking
all descendants of that ancestor as test data.

The tree structure is built directly from the branches file, avoiding the need
to parse a separate tree file.
"""

import argparse
import random
from collections import defaultdict


def build_tree_from_branches(branches_dict):
    """
    Build tree structure from branches dictionary.

    Args:
        branches_dict: dict mapping child -> (parent, hamming, train_test)

    Returns:
        children: dict mapping node_name -> list of child names
        parents: dict mapping child_name -> parent_name
        root: name of root node
    """
    children = defaultdict(list)
    parents = {}

    all_nodes = set()
    child_nodes = set()

    for child, (parent, _, _) in branches_dict.items():
        parents[child] = parent
        children[parent].append(child)
        all_nodes.add(parent)
        all_nodes.add(child)
        child_nodes.add(child)

    # Root is the node that has no parent (not a child of any other node)
    root_candidates = all_nodes - child_nodes
    if len(root_candidates) == 1:
        root = root_candidates.pop()
    elif len(root_candidates) > 1:
        # Multiple roots - pick the one with most descendants
        root = max(root_candidates, key=lambda n: len(children.get(n, [])))
    else:
        # All nodes are children - shouldn't happen with valid data
        root = None

    return dict(children), parents, root


def get_all_tips_iterative(children, root):
    """Get all tip (leaf) node names using iterative traversal."""
    if root is None:
        return []

    tips = []
    stack = [root]

    while stack:
        node = stack.pop()
        if node not in children or not children[node]:
            tips.append(node)
        else:
            stack.extend(children[node])

    return tips


def get_clade_tips_iterative(children, node):
    """Get all tips in a clade using iterative traversal."""
    tips = []
    stack = [node]

    while stack:
        current = stack.pop()
        if current not in children or not children[current]:
            tips.append(current)
        else:
            stack.extend(children[current])

    return tips


def get_all_descendants_iterative(children, node):
    """Get all descendants (tips + internal) using iterative traversal."""
    descendants = []
    stack = [node]

    while stack:
        current = stack.pop()
        descendants.append(current)
        if current in children:
            stack.extend(children[current])

    return descendants


def count_branch_mutations(branches_dict, node):
    """Count mutations on branch leading to this node (using Hamming distance)."""
    if node in branches_dict:
        _, hamming, _ = branches_dict[node]
        return hamming if hamming != "?" else 0
    return 0


def walk_back_mutations(node, parents, branches_dict, target_mutations):
    """
    Walk back from node counting mutations, return ancestor name.

    Returns the name of the ancestor node reached after accumulating
    at least target_mutations mutations.
    """
    current = node
    accumulated = 0

    while current in parents:
        branch_muts = count_branch_mutations(branches_dict, current)
        accumulated += branch_muts

        if accumulated >= target_mutations:
            return current

        current = parents[current]

    # Reached root
    return current


def iterative_test_selection(
    children, parents, root, branches_dict,
    test_proportion, mutations_back, max_clade_proportion, rng
):
    """
    Select test nodes by walking back from random tips and marking clades.

    Uses iterative traversal for all operations to handle large trees.
    """
    all_tips = get_all_tips_iterative(children, root)
    total_tips = len(all_tips)

    if total_tips == 0:
        print("Warning: No tips found in tree")
        return set(), set()

    target_test_count = max(1, int(total_tips * test_proportion))
    max_clade_tips = max(1, int(total_tips * max_clade_proportion))

    print(f"Total tips: {total_tips}")
    print(f"Target test tips: {target_test_count} ({test_proportion:.1%})")
    print(f"Max clade size: {max_clade_tips} tips")

    test_nodes = set()
    test_tips = set()
    tried_seeds = set()

    # Shuffle tips for random selection
    available_tips = list(all_tips)
    rng.shuffle(available_tips)

    while len(test_tips) < target_test_count and available_tips:
        # Get next seed tip that hasn't been tried and isn't already test
        seed_tip = None
        while available_tips:
            candidate = available_tips.pop()
            if candidate not in tried_seeds and candidate not in test_tips:
                seed_tip = candidate
                tried_seeds.add(candidate)
                break

        if seed_tip is None:
            break

        # Walk back from seed tip to find ancestor
        ancestor = walk_back_mutations(
            seed_tip, parents, branches_dict, mutations_back
        )

        # Check clade size
        clade_tips = get_clade_tips_iterative(children, ancestor)
        if len(clade_tips) > max_clade_tips:
            # Clade too large, skip this seed
            continue

        # Mark all descendants as test
        descendants = get_all_descendants_iterative(children, ancestor)
        for desc in descendants:
            test_nodes.add(desc)
            if desc in all_tips:
                test_tips.add(desc)

    return test_nodes, test_tips


def load_branches_raw(branches_path):
    """
    Load branches_raw.tsv and return dict mapping child -> (parent, hamming, train_test).
    """
    branches = {}

    with open(branches_path, "r") as f:
        header = f.readline()  # Skip header

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            parent = parts[0]
            child = parts[1]
            hamming = parts[2]

            try:
                hamming = int(hamming)
            except ValueError:
                hamming = "?"

            train_test = parts[3] if len(parts) > 3 else ""

            branches[child] = (parent, hamming, train_test)

    return branches


def write_branches_with_split(branches_raw, test_nodes, output_path):
    """Write branches.tsv with train_test labels populated."""
    print(f"Writing branches with train/test labels to {output_path}...")

    with open(output_path, "w") as f:
        f.write("parent\tchild\thamming\ttrain_test\n")

        for child, (parent, hamming, _) in branches_raw.items():
            label = "test" if child in test_nodes else "train"
            hamming_str = "?" if hamming == "?" else str(hamming)
            f.write(f"{parent}\t{child}\t{hamming_str}\t{label}\n")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--branches-raw", required=True,
        help="Input branches_raw.tsv (without train_test labels)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output branches.tsv with train_test labels"
    )
    parser.add_argument(
        "--test-proportion", type=float, default=0.1,
        help="Target proportion of tips as test (0.0-1.0, default: 0.1)"
    )
    parser.add_argument(
        "--mutations-back", type=int, default=5,
        help="Mutations to walk back from seed tip (default: 5)"
    )
    parser.add_argument(
        "--max-clade-proportion", type=float, default=0.01,
        help="Max size of any single test clade as proportion of total tips (default: 0.01)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility"
    )

    args = parser.parse_args()

    # Validate arguments
    if not 0.0 < args.test_proportion <= 1.0:
        raise ValueError("Test proportion must be between 0.0 and 1.0")
    if not 0.0 < args.max_clade_proportion <= 1.0:
        raise ValueError("Max clade proportion must be between 0.0 and 1.0")

    # Set up RNG
    rng = random.Random(args.seed)

    # Load branches
    print(f"Loading branches from {args.branches_raw}...")
    branches_raw = load_branches_raw(args.branches_raw)
    print(f"Loaded {len(branches_raw)} branches")

    # Build tree structure from branches
    print("Building tree structure from branches...")
    children, parents, root = build_tree_from_branches(branches_raw)
    print(f"Tree has {len(parents) + 1} nodes, root: {root}")

    # Select test nodes
    print("Selecting test clades...")
    test_nodes, test_tips = iterative_test_selection(
        children, parents, root, branches_raw,
        args.test_proportion,
        args.mutations_back,
        args.max_clade_proportion,
        rng
    )

    # Write output
    write_branches_with_split(branches_raw, test_nodes, args.output)

    # Report results
    all_tips = get_all_tips_iterative(children, root)
    total_tips = len(all_tips)
    test_tip_count = len(test_tips)
    train_tip_count = total_tips - test_tip_count
    achieved_proportion = test_tip_count / total_tips if total_tips > 0 else 0

    print(f"\nResults:")
    print(f"  Total tips: {total_tips}")
    print(f"  Test tips: {test_tip_count} ({achieved_proportion:.1%})")
    print(f"  Train tips: {train_tip_count} ({1-achieved_proportion:.1%})")
    print(f"  Target proportion: {args.test_proportion:.1%}")
    print(f"  Output written to: {args.output}")


if __name__ == "__main__":
    main()
