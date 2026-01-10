# code from Trevor Bedford and Katie Kistler
# MIT license
import argparse
import json
import os
import Bio.Phylo
from tqdm import tqdm

# Metadata fields to extract from node_attrs
# Note: We'll dynamically extract all categorical fields found in the tree
METADATA_FIELDS = ["clade_membership", "clade", "S1_mutations", "family_name"]

def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.
    Assigns links back to parent nodes for the root of the tree.
    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute.
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract metadata from Auspice JSON as TSV")
    parser.add_argument("--json", required=True, help="Path to the auspice.json file")
    parser.add_argument("--output", type=str, default="metadata.tsv", help="Output TSV file")

    args = parser.parse_args()

    # Load auspice JSON file
    with open(args.json, 'r') as f:
        tree_json = json.load(f)
    tree = json_to_tree(tree_json)

    data = []
    # Track which metadata fields are present in the data
    fields_present = set()
    all_fields_discovered = set()

    nodes = list(tree.find_clades(order="postorder"))

    # First pass: discover all categorical fields
    print("Discovering metadata fields...")
    for n in nodes[:100]:  # Sample first 100 nodes to discover fields
        if hasattr(n, 'node_attrs'):
            for field, value in n.node_attrs.items():
                if isinstance(value, dict) and 'value' in value:
                    # Check if it's likely a categorical field (has string value)
                    if isinstance(value['value'], str):
                        all_fields_discovered.add(field)

    # Combine discovered fields with known fields
    all_metadata_fields = list(set(METADATA_FIELDS) | all_fields_discovered)

    # Process all nodes
    for n in tqdm(nodes, desc="Processing nodes", unit="node"):
        node_elements = {"name": n.name.removeprefix('hCoV-19/')}
        node_elements["parent"] = n.parent.name.removeprefix('hCoV-19/') if n.parent else None

        # Extract all metadata fields dynamically
        if hasattr(n, 'node_attrs'):
            for field in all_metadata_fields:
                if field in n.node_attrs and 'value' in n.node_attrs[field]:
                    node_elements[field] = n.node_attrs[field]["value"]
                    fields_present.add(field)

        data.append(node_elements)

    # Determine column headers dynamically
    headers = ["name", "parent"]
    # Add metadata fields that were actually found in the data
    # Sort fields for consistent ordering
    for field in sorted(fields_present):
        headers.append(field)

    with open(args.output, 'w', encoding='utf-8') as handle:
        print(*headers, sep='\t', file=handle)
        for elem in data:
            row = [elem.get(header, "") for header in headers]
            print(*row, sep='\t', file=handle)
