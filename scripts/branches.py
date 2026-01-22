#!/usr/bin/env python3
"""
Extract parent-child relationships from Auspice JSON tree with Hamming distances.
Creates TSV file with columns: parent, child, hamming, train_test
"""

import argparse
import json
import Bio.Phylo
from Bio import SeqIO
from tqdm import tqdm

def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary."""
    # Check for v2 JSON which has combined metadata and tree data
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute
    if hasattr(node, "attr"):
        node.numdate = node.attr.get("num_date")
        node.branch_length = node.attr.get("div")

        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):
        node.branch_length = node.node_attrs.get("div")

    if root:
        node = annotate_parents_for_tree(node)

    return node

def annotate_parents_for_tree(tree):
    """Annotate each node in the given tree with its parent."""
    tree.root.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node
    return tree

def calculate_hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences."""
    if len(seq1) != len(seq2):
        print(f"Warning: Sequences have different lengths ({len(seq1)} vs {len(seq2)})")
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]

    # Count differences, ignoring gaps and ambiguous bases
    distance = 0
    for c1, c2 in zip(seq1, seq2):
        # Skip if either position has gap or ambiguous base
        if c1 in '-N' or c2 in '-N':
            continue
        if c1 != c2:
            distance += 1

    return distance

def extract_branches_with_hamming(tree, sequences):
    """
    Extract parent-child relationships with Hamming distances and train/test labels.

    Args:
        tree: Bio.Phylo tree
        sequences: Dictionary mapping sequence IDs to sequences

    Returns:
        List of (parent, child, hamming, train_test) tuples
    """
    branches = []

    # Get all nodes
    nodes = list(tree.find_clades(order="postorder"))

    for node in tqdm(nodes, desc="Processing nodes"):
        if node.parent is not None:
            # Clean node names (remove prefixes)
            child_name = node.name.removeprefix('hCoV-19/')
            parent_name = node.parent.name.removeprefix('hCoV-19/')

            # Calculate Hamming distance if sequences available
            if child_name in sequences and parent_name in sequences:
                child_seq = sequences[child_name]
                parent_seq = sequences[parent_name]
                hamming = calculate_hamming_distance(parent_seq, child_seq)
            else:
                # One or both sequences missing (common for internal nodes)
                hamming = -1  # Flag for missing sequence(s)

            # Extract train_test label from node_attrs
            train_test = ""
            if hasattr(node, 'node_attrs') and node.node_attrs:
                train_test_attr = node.node_attrs.get('train_test', {})
                train_test = train_test_attr.get('value', '')

            branches.append((parent_name, child_name, hamming, train_test))

    return branches

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", required=True,
                       help="Path to the auspice.json file")
    parser.add_argument("--alignment", required=True,
                       help="Path to the FASTA alignment file")
    parser.add_argument("--output", default="branches.tsv",
                       help="Output TSV file")

    args = parser.parse_args()

    # Load auspice JSON file
    print("Loading tree from JSON...")
    with open(args.json, 'r') as f:
        tree_json = json.load(f)
    tree = json_to_tree(tree_json)

    # Load sequences
    print("Loading sequences...")
    sequences = {}
    for record in SeqIO.parse(args.alignment, "fasta"):
        # Store both with and without common prefixes
        seq_str = str(record.seq).upper()
        sequences[record.id] = seq_str
        # Also store without prefix
        clean_id = record.id.removeprefix('hCoV-19/')
        sequences[clean_id] = seq_str

    print(f"Loaded {len(sequences)} sequences")

    # Extract branches with Hamming distances
    print("Extracting branches and computing Hamming distances...")
    branches = extract_branches_with_hamming(tree, sequences)

    # Keep all branches, including those with missing sequences
    valid_branches = branches
    missing_branches = sum(1 for _, _, h, _ in branches if h == -1)

    print(f"Found {len(branches)} total branches")
    if missing_branches > 0:
        print(f"  {missing_branches} branches have missing sequences (marked with '?' for hamming)")

    # Write output
    print(f"Writing to {args.output}...")
    with open(args.output, 'w') as f:
        # Write header
        f.write("parent\tchild\thamming\ttrain_test\n")

        # Write branches
        for parent, child, hamming, train_test in valid_branches:
            # Use '?' for missing sequences (hamming = -1)
            hamming_str = '?' if hamming == -1 else str(hamming)
            f.write(f"{parent}\t{child}\t{hamming_str}\t{train_test}\n")

    # Print statistics
    if valid_branches:
        # Only compute statistics for branches with valid hamming values
        hamming_values = [h for _, _, h, _ in valid_branches if h >= 0]

        if hamming_values:
            print(f"\nHamming distance statistics (for branches with sequences):")
            print(f"  Mean: {sum(hamming_values)/len(hamming_values):.1f}")
            print(f"  Min: {min(hamming_values)}")
            print(f"  Max: {max(hamming_values)}")
            print(f"  Branches with sequences: {len(hamming_values)}/{len(valid_branches)}")
        else:
            print(f"\nNo branches with valid sequences found")
            print(f"  Total branches: {len(valid_branches)}")

        # Distribution (only for valid hamming values)
        if hamming_values:
            print(f"\nDistance distribution:")
            bins = [0, 1, 5, 10, 20, 50, 100, 200, 500, 1000]
            for i in range(len(bins)-1):
                count = sum(1 for h in hamming_values if bins[i] <= h < bins[i+1])
                if count > 0:
                    pct = 100.0 * count / len(hamming_values)
                    print(f"  [{bins[i]:4d}-{bins[i+1]:4d}): {count:5d} ({pct:5.1f}%)")

if __name__ == "__main__":
    main()