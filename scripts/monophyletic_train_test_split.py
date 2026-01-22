#!/usr/bin/env python3
"""
Create train/test split from trajectories based on monophyletic clades.

This ensures test set is phylogenetically coherent (one contiguous subtree),
which is useful for evaluating generalization to held-out clades.

Usage:
    python scripts/monophyletic_train_test_split.py \
        --auspice data/spike-xs/auspice.json \
        --trajectories results/spike-xs \
        --output results/spike-xs/train_test_split \
        --test-fraction 0.2
"""

import argparse
import json
import shutil
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade


def auspice_to_newick(auspice_path, output_path):
    """Convert auspice JSON tree to Newick format."""
    with open(auspice_path) as f:
        data = json.load(f)
    
    def to_clade(node):
        div = node.get('node_attrs', {}).get('div', 0)
        name = node.get('name')
        is_tip = 'children' not in node
        
        clade = Clade(branch_length=div, name=name if is_tip else None)
        
        if 'children' in node:
            for child in node['children']:
                child_clade = to_clade(child)
                child_div = child.get('node_attrs', {}).get('div', 0)
                child_clade.branch_length = child_div - div
                clade.clades.append(child_clade)
        
        return clade
    
    root = to_clade(data['tree'])
    root.branch_length = 0
    tree = Tree(root=root)
    Phylo.write(tree, str(output_path), 'newick')
    return tree


def find_monophyletic_clade(tree, target_size, tolerance=0.5):
    """Find a clade closest to target_size for test set."""
    candidates = []
    
    def search(clade, depth=0):
        n_tips = len(clade.get_terminals())
        if target_size * (1-tolerance) <= n_tips <= target_size * (1+tolerance):
            candidates.append((clade, n_tips, depth))
        if hasattr(clade, 'clades') and clade.clades:
            for child in clade.clades:
                search(child, depth + 1)
    
    search(tree.root)
    
    if not candidates:
        raise ValueError(f"No clade found with size near {target_size}")
    
    # Return clade closest to target size
    return min(candidates, key=lambda x: abs(x[1] - target_size))


def tip_to_filename(tip_name):
    """Convert tree tip name to trajectory filename."""
    name = tip_name.replace("hCoV-19/", "").replace("/", "")
    return f"{name}.fasta"


def collect_nodes_from_trajectories(traj_dir, tips):
    """Collect all node names that appear in trajectories for given tips."""
    nodes = set()
    for tip in tips:
        src = traj_dir / tip_to_filename(tip)
        if src.exists():
            with open(src) as f:
                for line in f:
                    if line.startswith(">"):
                        # Header format: >NODE_XXXXXXX|position or >tip_name|position
                        node_name = line[1:].split("|")[0].strip()
                        nodes.add(node_name)
    return nodes


def filter_trajectory(src_path, dest_path, exclude_nodes):
    """Copy trajectory file, excluding nodes in exclude_nodes set."""
    with open(src_path) as f:
        content = f.read()

    # Parse FASTA entries
    entries = []
    current_header = None
    current_seq = []

    for line in content.strip().split("\n"):
        if line.startswith(">"):
            if current_header is not None:
                entries.append((current_header, "".join(current_seq)))
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

    if current_header is not None:
        entries.append((current_header, "".join(current_seq)))

    # Filter out excluded nodes
    filtered = []
    for header, seq in entries:
        node_name = header[1:].split("|")[0]
        if node_name not in exclude_nodes:
            filtered.append((header, seq))

    # Write filtered trajectory
    with open(dest_path, "w") as f:
        for header, seq in filtered:
            f.write(f"{header}\n{seq}\n")


def create_auspice_with_split(auspice_path, train_tips, test_tips, output_path):
    """Add train/test split coloring to auspice JSON."""
    with open(auspice_path) as f:
        auspice = json.load(f)
    
    train_set = set(train_tips)
    test_set = set(test_tips)
    
    def add_split_attr(node):
        name = node.get("name", "")
        if "node_attrs" not in node:
            node["node_attrs"] = {}
        
        if name in train_set:
            node["node_attrs"]["split"] = {"value": "train"}
        elif name in test_set:
            node["node_attrs"]["split"] = {"value": "test"}
        
        if "children" in node:
            for child in node["children"]:
                add_split_attr(child)
    
    add_split_attr(auspice["tree"])
    
    # Add coloring definition
    auspice["meta"]["colorings"] = [
        c for c in auspice["meta"].get("colorings", []) 
        if c.get("key") != "split"
    ]
    auspice["meta"]["colorings"].insert(0, {
        "key": "split",
        "title": "Train/Test Split",
        "type": "categorical",
        "scale": [["train", "#4C78A8"], ["test", "#E45756"]]
    })
    
    with open(output_path, "w") as f:
        json.dump(auspice, f)


def main():
    parser = argparse.ArgumentParser(
        description="Create train/test split based on monophyletic clades"
    )
    parser.add_argument("--auspice", required=True, help="Input auspice JSON")
    parser.add_argument("--trajectories", required=True, help="Directory with trajectory FASTA files")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--test-fraction", type=float, default=0.2, help="Fraction for test set (default: 0.2)")
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    traj_dir = Path(args.trajectories)
    
    # Convert auspice to newick and load tree
    tree_path = output_dir / "tree.nwk"
    print(f"Extracting tree from {args.auspice}...")
    tree = auspice_to_newick(args.auspice, tree_path)
    
    total_tips = len(tree.get_terminals())
    print(f"Total tips: {total_tips}")
    
    # Find monophyletic clade for test set
    target_test_size = int(total_tips * args.test_fraction)
    print(f"Finding monophyletic clade with ~{target_test_size} tips for test set...")
    
    test_clade, test_size, depth = find_monophyletic_clade(tree, target_test_size)
    test_tips = [t.name for t in test_clade.get_terminals()]
    train_tips = [t.name for t in tree.get_terminals() if t.name not in set(test_tips)]
    
    print(f"Selected test clade: {test_size} tips at depth {depth}")
    print(f"Train: {len(train_tips)}, Test: {len(test_tips)}")

    # Save tip lists
    with open(output_dir / "train_tips.txt", "w") as f:
        f.write("\n".join(sorted(train_tips)) + "\n")
    with open(output_dir / "test_tips.txt", "w") as f:
        f.write("\n".join(sorted(test_tips)) + "\n")
    
    # Copy trajectory files
    train_dir = output_dir / "train"
    test_dir = output_dir / "test"
    train_dir.mkdir(exist_ok=True)
    test_dir.mkdir(exist_ok=True)
    
    print("Copying train trajectory files...")
    for tip in train_tips:
        src = traj_dir / tip_to_filename(tip)
        if src.exists():
            shutil.copy(src, train_dir / src.name)

    # Collect all nodes from train trajectories
    print("Collecting nodes from train trajectories...")
    train_nodes = collect_nodes_from_trajectories(traj_dir, train_tips)
    print(f"Found {len(train_nodes)} unique nodes in train trajectories")

    # Filter test trajectories to exclude train nodes
    print("Filtering test trajectory files (removing shared ancestors)...")
    for tip in test_tips:
        src = traj_dir / tip_to_filename(tip)
        if src.exists():
            filter_trajectory(src, test_dir / src.name, train_nodes)

    print(f"Copied {len(list(train_dir.glob('*.fasta')))} train, {len(list(test_dir.glob('*.fasta')))} test files")
    
    # Create auspice JSON with split coloring
    auspice_out = output_dir / "auspice_with_split.json"
    print(f"Creating {auspice_out}...")
    create_auspice_with_split(args.auspice, train_tips, test_tips, auspice_out)
    
    print(f"\nDone! Output in {output_dir}")
    print(f"  - train/: {len(train_tips)} trajectory files")
    print(f"  - test/: {len(test_tips)} trajectory files (monophyletic)")
    print(f"  - auspice_with_split.json: view in auspice.us")
    print(f"\nTo package for S3:")
    print(f"  python scripts/package.py --input-dir {output_dir}/train --output-dir {output_dir}/export_train --shuffle")
    print(f"  python scripts/package.py --input-dir {output_dir}/test --output-dir {output_dir}/export_test --shuffle")


if __name__ == "__main__":
    main()
