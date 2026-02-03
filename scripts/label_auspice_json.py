#!/usr/bin/env python3
"""
Add train/test labels from branches.tsv to Auspice JSON for visualization.

This script reads an Auspice JSON and a branches.tsv file (with train_test column),
and adds the train_test labels to each node in the JSON for visualization in Auspice.
"""

import argparse
import json


def load_train_test_labels(branches_path):
    """Load train/test labels from branches.tsv."""
    labels = {}
    with open(branches_path, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                child = parts[1]
                train_test = parts[3]
                labels[child] = train_test
    return labels


def annotate_nodes_train_test(tree_dict, labels):
    """Set node_attrs.train_test.value for all nodes based on labels dict."""
    def traverse(node):
        name = node.get("name")
        if name:
            if "node_attrs" not in node:
                node["node_attrs"] = {}
            # Strip hCoV-19/ prefix to match branches.tsv format
            clean_name = name.removeprefix("hCoV-19/")
            # Use label from branches.tsv, default to "train" if not found
            value = labels.get(clean_name, "train")
            node["node_attrs"]["train_test"] = {"value": value}

        for child in node.get("children", []):
            traverse(child)

    traverse(tree_dict)


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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--json", required=True,
        help="Input Auspice JSON file (without train/test labels)"
    )
    parser.add_argument(
        "--branches", required=True,
        help="Input branches.tsv file (with train_test column)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output Auspice JSON file (with train/test labels)"
    )

    args = parser.parse_args()

    # Load train/test labels from branches.tsv
    print(f"Loading train/test labels from {args.branches}...")
    labels = load_train_test_labels(args.branches)
    print(f"  Found labels for {len(labels)} nodes")

    # Count train vs test
    train_count = sum(1 for v in labels.values() if v == "train")
    test_count = sum(1 for v in labels.values() if v == "test")
    print(f"  Train: {train_count}, Test: {test_count}")

    # Load Auspice JSON
    print(f"Loading Auspice JSON from {args.json}...")
    with open(args.json, 'r') as f:
        auspice_data = json.load(f)

    tree_dict = auspice_data.get("tree")
    if not tree_dict:
        raise ValueError("No 'tree' key found in Auspice JSON")

    # Annotate nodes with train/test labels
    print("Annotating nodes with train/test labels...")
    annotate_nodes_train_test(tree_dict, labels)

    # Add coloring metadata
    add_train_test_coloring(auspice_data)

    # Write output
    print(f"Writing labeled JSON to {args.output}...")
    with open(args.output, 'w') as f:
        json.dump(auspice_data, f, separators=(',', ':'))

    print("Done!")


if __name__ == "__main__":
    main()
