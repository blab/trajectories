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

Sequence reconstruction is linear in the total number of mutations. A single
depth-first traversal carries the parent sequence down to each child and
applies only that child's own branch mutations -- the parent already encodes
every ancestral mutation. A mutable buffer is mutated on descent and undone
on ascent, so the same node sequence is produced as re-walking the full
root->node path, but in O(total mutations) instead of O(nodes * depth).
"""
import argparse
import json
import os
# augur sets the Python recursion limit on import, from the AUGUR_RECURSION_LIMIT
# env var (default 10000). json_to_tree() below recurses to tree depth, so it
# relies on that -- set AUGUR_RECURSION_LIMIT for very deep trees. No manual
# sys.setrecursionlimit() here: it would just be overwritten by this import.
from augur.utils import json_to_tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def reconstruct_all_sequences(tree, root_seq, gene):
    """
    Reconstruct the sequence of every node in the tree in a single pass.

    Walks the tree once with an explicit stack (no recursion -- trees can be
    very deep) carrying one mutable sequence buffer. Descending into a node
    applies that node's own branch mutations; ascending undoes them, exactly
    restoring the parent state. For every node this yields the root sequence
    with every mutation on the root->node path applied in root->tip order
    (a later mutation at the same site overwrites an earlier one) -- identical
    to flattening tree.get_path(node) and applying each mutation to the root.

    Returns a dict mapping each Bio.Phylo clade object to its sequence string.
    """
    # make the root sequence mutable; this single buffer is reused for every node
    buf = MutableSeq(root_seq)
    sequences = {}

    # The root keeps the given root sequence unchanged: its own branch_attrs
    # are never applied (tree.get_path() excludes the root in the original).
    sequences[tree.root] = str(buf)

    # Stack entries are ("descend", node) or ("ascend", undo_list).
    stack = [("descend", child) for child in reversed(tree.root.clades)]

    pbar = tqdm(desc="Reconstructing sequences", unit="node")
    while stack:
        kind, payload = stack.pop()

        if kind == "ascend":
            # Revert this node's mutations, last-applied first, so that
            # repeated mutations at the same site unwind to the exact
            # parent base.
            for site, old_base in reversed(payload):
                buf[site] = old_base
            continue

        node = payload
        # Apply only this node's own branch mutations to the parent sequence.
        muts = node.branch_attrs['mutations'].get(gene, [])
        undo = []
        for mut in muts:
            # subtract 1 to deal with biological numbering vs python
            mut_site = int(mut[1:-1]) - 1
            # get the nucleotide that the site was mutated to
            mutation = mut[-1]
            # record the pre-change base so the descent can be undone exactly
            undo.append((mut_site, buf[mut_site]))
            # apply mutation
            buf[mut_site] = mutation
        sequences[node] = str(buf)

        # The undo runs only after the whole subtree below this node is done.
        stack.append(("ascend", undo))
        for child in reversed(node.clades):
            stack.append(("descend", child))
        pbar.update(1)
    pbar.close()

    return sequences


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
    parser.add_argument("--trim-begin", type=int, default=None,
                        help="Start position for trimming (1-indexed, inclusive)")
    parser.add_argument("--trim-end", type=int, default=None,
                        help="End position for trimming (1-indexed, inclusive)")
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

    # Reconstruct every node's sequence in a single linear traversal
    sequences = reconstruct_all_sequences(tree, root_seq, args.gene)

    # Initialize list to store sequence records for each node
    sequence_records = []

    # Iterate over tip nodes only if --tips-only is set, otherwise iterate over all nodes
    nodes = list(tree.get_terminals() if args.tips_only else tree.find_clades())

    for node in tqdm(nodes, desc="Processing nodes", unit="node"):

        # Get sequence at node (precomputed in the linear traversal above).
        # pop() frees each reconstructed sequence as soon as it is consumed,
        # so the dict and the output records are not both fully held at once.
        node_seq = sequences.pop(node)
        # Strip trailing stop codons
        stripped_seq = Seq(str(node_seq).rstrip('*'))
        # Apply trimming if specified
        if args.trim_begin is not None or args.trim_end is not None:
            start = (args.trim_begin - 1) if args.trim_begin else 0
            end = args.trim_end if args.trim_end else len(stripped_seq)
            stripped_seq = stripped_seq[start:end]
        # Strip hCoV-19/ from beginning of strain name
        strain = node.name.removeprefix('hCoV-19/')
        # Only keep records without stop codons (*)
        if '*' not in stripped_seq:
            sequence_records.append(SeqRecord(stripped_seq, strain, '', ''))

    # Write sequences to output FASTA file
    SeqIO.write(sequence_records, args.output, "fasta")
