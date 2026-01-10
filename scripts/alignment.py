# Code from Katie Kistler
# MIT license
"""
Given an auspice.json file, finds the sequences of each node in the tree
and outputs a FASTA file with these node sequences.

The root sequence can be either:
1. Embedded in the auspice.json file under 'root_sequence'
2. In a sidecar file named {auspice}_root-sequence.json

If a gene is specified, the sequences will be the AA sequence of that gene
at that node. If 'nuc' is specified, the whole genome nucleotide sequence
at the node will be output. (this is default if no gene is specified).
The FASTA header is the node's name in the tree.
"""
import argparse
import json
import os
from augur.utils import json_to_tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

def apply_muts_to_root(root_seq, list_of_muts):
    """
    Apply a list of mutations to the root sequence
    to find the sequence at a given node. The list of mutations
    is ordered from root to node, so multiple mutations at the
    same site will correctly overwrite each other
    """
    # make the root sequence mutable
    root_plus_muts = MutableSeq(root_seq)

    # apply all mutations to root sequence
    for mut in list_of_muts:
        # subtract 1 to deal with biological numbering vs python
        mut_site = int(mut[1:-1])-1
        # get the nucleotide that the site was mutated to
        mutation = mut[-1]
        # apply mutation
        root_plus_muts[mut_site] = mutation

    return root_plus_muts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--gene", default="nuc",
                        help="Name of gene to return AA sequences for. 'nuc' will return full genome nucleotide seq")
    parser.add_argument("--json", required=True,
                        help="Path to the auspice.json file")
    parser.add_argument("--output", type=str, default="alignment.fasta",
                        help="Output FASTA file for sequences")
    parser.add_argument("--tips-only", action="store_true",
                        help="If set, only include tip sequences (leaf nodes)")

    args = parser.parse_args()

    # Load auspice JSON file
    with open(args.json, 'r') as f:
        auspice_json = json.load(f)

    # Auto-detect root sequence source
    root_seq = None

    # First, check for sidecar _root-sequence.json file
    json_dir = os.path.dirname(args.json)
    json_base = os.path.basename(args.json).replace('.json', '')
    sidecar_path = os.path.join(json_dir, f"{json_base}_root-sequence.json")

    if os.path.exists(sidecar_path):
        # Use sidecar file
        with open(sidecar_path, 'r') as f:
            root_json = json.load(f)
        root_seq = root_json.get(args.gene)
        if not root_seq:
            raise ValueError(f"Root sequence for gene '{args.gene}' not found in sidecar file {sidecar_path}")
    else:
        # No sidecar file, check for embedded root sequence
        if 'root_sequence' in auspice_json:
            root_seq = auspice_json['root_sequence'].get(args.gene)
            if not root_seq:
                raise ValueError(f"Root sequence for gene '{args.gene}' not found in auspice.json's root_sequence field")
        else:
            raise ValueError(f"No root sequence found. Expected either:\n"
                           f"1. Sidecar file at {sidecar_path}, or\n"
                           f"2. 'root_sequence' field in auspice.json")

    # Convert tree JSON to Bio.phylo format
    tree = json_to_tree(auspice_json)

    # Initialize list to store sequence records for each node
    sequence_records = []

    # Iterate over tip nodes only if --tips-only is set, otherwise iterate over all nodes
    nodes = list(tree.get_terminals() if args.tips_only else tree.find_clades())

    for node in tqdm(nodes, desc="Processing nodes", unit="node"):

        # Get path back to the root
        path = tree.get_path(node)

        # Get all mutations relative to root
        muts = [branch.branch_attrs['mutations'].get(args.gene, []) for branch in path]
        # Flatten the list of mutations
        muts = [item for sublist in muts for item in sublist]
        # Get sequence at node
        node_seq = apply_muts_to_root(root_seq, muts)
        # Strip trailing stop codons
        stripped_seq = Seq(str(node_seq).rstrip('*'))
        # Strip hCoV-19/ from beginning of strain name
        strain = node.name.removeprefix('hCoV-19/')
        # Only keep records without stop codons (*)
        if '*' not in stripped_seq:
            sequence_records.append(SeqRecord(stripped_seq, strain, '', ''))

    # Write sequences to output FASTA file
    SeqIO.write(sequence_records, args.output, "fasta")
