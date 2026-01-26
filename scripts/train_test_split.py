#!/usr/bin/env python3
"""
Mark entire clades as "test" data in an Auspice JSON file for train/test splits.

This enables proper train/test splits where test data is truly out-of-sample
(not sharing evolutionary paths with training data).
"""

import argparse
import json
import random


def count_mutations_on_branch(node_dict, gene="nuc", trim_begin=None, trim_end=None):
    """Count mutations in branch_attrs.mutations[gene] for this node.

    If trim_begin and trim_end are specified, only count mutations within that range.
    Positions are 1-indexed and inclusive on both ends.
    """
    branch_attrs = node_dict.get("branch_attrs", {})
    mutations = branch_attrs.get("mutations", {})
    gene_mutations = mutations.get(gene, [])

    if trim_begin is None or trim_end is None:
        return len(gene_mutations)

    count = 0
    for mut_str in gene_mutations:
        try:
            position = int(mut_str[1:-1])
            if trim_begin <= position <= trim_end:
                count += 1
        except (ValueError, IndexError):
            continue

    return count


def build_node_map(tree_dict):
    """Build dict mapping node name -> node dict."""
    node_map = {}

    def traverse(node):
        name = node.get("name")
        if name:
            node_map[name] = node
        for child in node.get("children", []):
            traverse(child)

    traverse(tree_dict)
    return node_map


def build_parent_map(tree_dict):
    """Build dict mapping child name -> parent name."""
    parent_map = {}

    def traverse(node, parent_name=None):
        name = node.get("name")
        if name and parent_name:
            parent_map[name] = parent_name
        for child in node.get("children", []):
            traverse(child, name)

    traverse(tree_dict)
    return parent_map


def get_all_tips(tree_dict):
    """Get list of all tip (leaf) node names."""
    tips = []

    def traverse(node):
        if "children" not in node or len(node.get("children", [])) == 0:
            name = node.get("name")
            if name:
                tips.append(name)
        else:
            for child in node.get("children", []):
                traverse(child)

    traverse(tree_dict)
    return tips


def walk_back_mutations(tip_name, parent_map, node_map, target_mutations, gene, trim_begin=None, trim_end=None):
    """
    Walk back from tip counting mutations, return ancestor name.

    Returns the name of the ancestor node reached after accumulating
    at least target_mutations mutations. If we reach the root before
    accumulating enough mutations, returns the root.

    If trim_begin and trim_end are specified, only mutations within that
    position range are counted.
    """
    current_name = tip_name
    accumulated_mutations = 0

    while current_name in parent_map:
        parent_name = parent_map[current_name]
        parent_node = node_map[parent_name]

        # Count mutations on the branch leading to current node
        current_node = node_map[current_name]
        branch_mutations = count_mutations_on_branch(current_node, gene, trim_begin, trim_end)
        accumulated_mutations += branch_mutations

        if accumulated_mutations >= target_mutations:
            return current_name

        current_name = parent_name

    # Reached root without accumulating enough mutations
    return current_name


def get_all_descendants(node_dict):
    """Recursively get all descendant node names (tips + internal)."""
    descendants = []

    def traverse(node):
        name = node.get("name")
        if name:
            descendants.append(name)
        for child in node.get("children", []):
            traverse(child)

    traverse(node_dict)
    return descendants


def get_clade_tip_count(node_dict):
    """Count tips in a clade (for size checking)."""
    count = 0

    def traverse(node):
        nonlocal count
        if "children" not in node or len(node.get("children", [])) == 0:
            count += 1
        else:
            for child in node.get("children", []):
                traverse(child)

    traverse(node_dict)
    return count


def iterative_test_selection(tree_dict, test_proportion, mutations_back, max_clade_proportion, gene, rng, trim_begin=None, trim_end=None):
    """
    Main loop: select seed tips, find ancestors, skip oversized clades, mark until target reached.

    Returns a set of node names marked as test.

    If trim_begin and trim_end are specified, only mutations within that
    position range are counted when walking back.
    """
    all_tips = get_all_tips(tree_dict)
    total_tips = len(all_tips)
    target_test_count = max(1, int(total_tips * test_proportion))
    max_clade_tips = int(total_tips * max_clade_proportion)

    node_map = build_node_map(tree_dict)
    parent_map = build_parent_map(tree_dict)

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
        ancestor_name = walk_back_mutations(seed_tip, parent_map, node_map, mutations_back, gene, trim_begin, trim_end)
        ancestor_node = node_map[ancestor_name]

        # Check clade size
        clade_tip_count = get_clade_tip_count(ancestor_node)
        if clade_tip_count > max_clade_tips:
            # Clade too large, skip this seed
            continue

        # Mark all descendants as test
        descendants = get_all_descendants(ancestor_node)
        for desc in descendants:
            test_nodes.add(desc)
            # Track if this is a tip
            if desc in all_tips:
                test_tips.add(desc)

    return test_nodes, test_tips


def find_monophyletic_clade(tree_dict, target_size, tolerance=0.5):
    """Find a single clade closest to target_size for test set.

    JSON traversal version - no BioPython dependency.

    Args:
        tree_dict: Root node of tree as dict with 'children' key
        target_size: Target number of tips in the clade
        tolerance: Acceptable deviation from target (0.5 = ±50%)

    Returns:
        Tuple of (node_dict, tip_count, depth) for the best matching clade
    """
    candidates = []
    min_size = target_size * (1 - tolerance)
    max_size = target_size * (1 + tolerance)

    def search(node, depth=0):
        tip_count = get_clade_tip_count(node)
        if min_size <= tip_count <= max_size:
            candidates.append((node, tip_count, depth))
        for child in node.get("children", []):
            search(child, depth + 1)

    search(tree_dict)

    if not candidates:
        raise ValueError(f"No clade found with size near {target_size} (tolerance: ±{tolerance*100:.0f}%)")

    # Return clade closest to target size
    return min(candidates, key=lambda x: abs(x[1] - target_size))


def monophyletic_test_selection(tree_dict, test_proportion, tolerance=0.5):
    """Select test nodes using single monophyletic clade strategy.

    Finds a single contiguous clade in the tree that is closest to the
    target test size. This ensures the test set is phylogenetically
    coherent (one contiguous subtree).

    Args:
        tree_dict: Root node of tree as dict
        test_proportion: Target proportion of tips for test set
        tolerance: Acceptable deviation from target size (default: 0.5)

    Returns:
        Tuple of (test_nodes set, test_tips set) like iterative_test_selection
    """
    all_tips = get_all_tips(tree_dict)
    total_tips = len(all_tips)
    target_test_count = max(1, int(total_tips * test_proportion))

    # Find the best monophyletic clade
    clade_node, clade_tip_count, depth = find_monophyletic_clade(
        tree_dict, target_test_count, tolerance
    )

    # Get all descendants of the selected clade
    test_nodes = set(get_all_descendants(clade_node))
    test_tips = set(tip for tip in all_tips if tip in test_nodes)

    print(f"Selected monophyletic clade: {clade_tip_count} tips at depth {depth}")

    return test_nodes, test_tips


def add_train_test_coloring(auspice_json):
    """Add train_test to meta.colorings with categorical scale."""
    if "meta" not in auspice_json:
        auspice_json["meta"] = {}
    if "colorings" not in auspice_json["meta"]:
        auspice_json["meta"]["colorings"] = []

    # Check if train_test coloring already exists
    for coloring in auspice_json["meta"]["colorings"]:
        if coloring.get("key") == "train_test":
            # Update existing
            coloring["title"] = "Train/Test Split"
            coloring["type"] = "categorical"
            coloring["scale"] = [["train", "#4C78A8"], ["test", "#E45756"]]
            return

    # Add new coloring
    train_test_coloring = {
        "key": "train_test",
        "title": "Train/Test Split",
        "type": "categorical",
        "scale": [["train", "#4C78A8"], ["test", "#E45756"]]
    }
    auspice_json["meta"]["colorings"].append(train_test_coloring)


def annotate_nodes_train_test(tree_dict, test_node_names):
    """Set node_attrs.train_test.value for all nodes."""

    def traverse(node):
        name = node.get("name")
        if name:
            if "node_attrs" not in node:
                node["node_attrs"] = {}
            value = "test" if name in test_node_names else "train"
            node["node_attrs"]["train_test"] = {"value": value}

        for child in node.get("children", []):
            traverse(child)

    traverse(tree_dict)


def main():
    parser = argparse.ArgumentParser(
        description="Mark entire clades as test data in an Auspice JSON file"
    )
    parser.add_argument(
        "--json", required=True,
        help="Input Auspice JSON file"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output Auspice JSON file"
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
        "--max-clade-proportion", type=float, default=0.1,
        help="Max size of any single test clade as proportion of total tips (default: 0.1)"
    )
    parser.add_argument(
        "--gene", type=str, default="nuc",
        help="Gene key for counting mutations (default: nuc)"
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--trim-begin", type=str, default="",
        help="Start position for filtering mutations (1-indexed, inclusive). Empty string means no filtering."
    )
    parser.add_argument(
        "--trim-end", type=str, default="",
        help="End position for filtering mutations (1-indexed, inclusive). Empty string means no filtering."
    )
    parser.add_argument(
        "--strategy", type=str, choices=["random-clades", "monophyletic"], default="random-clades",
        help="Test selection strategy: random-clades (multiple random clades) or monophyletic (single contiguous clade)"
    )
    parser.add_argument(
        "--tolerance", type=float, default=0.5,
        help="Size tolerance for monophyletic strategy (default: 0.5, meaning ±50%% of target size)"
    )

    args = parser.parse_args()

    # Validate arguments
    if not 0.0 < args.test_proportion <= 1.0:
        raise ValueError("Test proportion must be between 0.0 and 1.0")
    if not 0.0 < args.max_clade_proportion <= 1.0:
        raise ValueError("Max clade proportion must be between 0.0 and 1.0")
    if args.mutations_back < 0:
        raise ValueError("Mutations back must be non-negative")
    if not 0.0 < args.tolerance <= 1.0:
        raise ValueError("Tolerance must be between 0.0 and 1.0")

    # Set up random number generator
    rng = random.Random(args.seed)

    # Convert trim args to int or None
    trim_begin = int(args.trim_begin) if args.trim_begin else None
    trim_end = int(args.trim_end) if args.trim_end else None

    # Load Auspice JSON
    with open(args.json, 'r') as f:
        auspice_data = json.load(f)

    tree_dict = auspice_data.get("tree")
    if not tree_dict:
        raise ValueError("No 'tree' key found in Auspice JSON")

    # Get total tips for reporting
    all_tips = get_all_tips(tree_dict)
    total_tips = len(all_tips)

    # Select test nodes based on strategy
    if args.strategy == "monophyletic":
        test_nodes, test_tips = monophyletic_test_selection(
            tree_dict,
            args.test_proportion,
            args.tolerance
        )
    else:  # random-clades (default)
        test_nodes, test_tips = iterative_test_selection(
            tree_dict,
            args.test_proportion,
            args.mutations_back,
            args.max_clade_proportion,
            args.gene,
            rng,
            trim_begin,
            trim_end
        )

    # Annotate nodes with train/test labels
    annotate_nodes_train_test(tree_dict, test_nodes)

    # Add coloring metadata
    add_train_test_coloring(auspice_data)

    # Write output
    with open(args.output, 'w') as f:
        json.dump(auspice_data, f, separators=(',', ':'))

    # Report results
    test_tip_count = len(test_tips)
    train_tip_count = total_tips - test_tip_count
    achieved_proportion = test_tip_count / total_tips if total_tips > 0 else 0

    print(f"Strategy: {args.strategy}")
    print(f"Total tips: {total_tips}")
    print(f"Test tips: {test_tip_count} ({achieved_proportion:.1%})")
    print(f"Train tips: {train_tip_count} ({1-achieved_proportion:.1%})")
    print(f"Target proportion: {args.test_proportion:.1%}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
