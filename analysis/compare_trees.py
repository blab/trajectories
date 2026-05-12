#!/usr/bin/env python3
"""
Compare UShER and Viridian mutation-annotated trees for data quality.

Subsamples both trees to a common number of tips, extracts mutations via
matUtils, and compares reversion-to-reference and indel rates. A cleaner
tree should have fewer reversions and fewer indels per branch.

Reuses helper functions from the main pipeline's provisioning script.
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

# Add scripts directory to path so we can import pipeline utilities
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from usher_provision_alignment_branches import (
    find_matutils,
    run_matutils_extract,
    parse_newick_iterative,
    parse_mutations_file,
    fetch_reference_sequence,
)


def decompress_pb(pb_path, work_dir):
    """Decompress .pb.gz to .pb, returning path to decompressed file."""
    if not pb_path.endswith(".gz"):
        return pb_path
    decompressed = os.path.join(work_dir, "tree.pb")
    print(f"  Decompressing {os.path.basename(pb_path)}...")
    with gzip.open(pb_path, "rb") as f_in, open(decompressed, "wb") as f_out:
        while True:
            chunk = f_in.read(1024 * 1024)
            if not chunk:
                break
            f_out.write(chunk)
    return decompressed


def extract_tree_data(pb_path, work_dir, matutils_path, subsample=None, seed=42):
    """Extract tree structure and mutations from a protobuf file.

    Optionally subsamples to a target number of tips using matUtils -z flag
    in the same extraction step.

    Returns (children, parents, root, mutations) dicts.
    """
    pb_file = decompress_pb(pb_path, work_dir)

    abs_output_dir = os.path.abspath(work_dir)
    abs_pb_file = os.path.abspath(pb_file)
    tree_path = os.path.join(abs_output_dir, "tree.nwk")
    mutations_path = os.path.join(abs_output_dir, "all_mutations.txt")

    cmd = [
        matutils_path, "extract",
        "-i", abs_pb_file,
        "-t", "tree.nwk",
        "-A", "all_mutations.txt",
    ]
    if subsample:
        cmd.extend(["-z", str(subsample)])
        print(f"  Extracting with subsample to {subsample} tips...")
    else:
        print(f"  Extracting tree and mutations...")

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=abs_output_dir)
    if result.returncode != 0:
        print(f"  matUtils stderr: {result.stderr[:500]}")
        result.check_returncode()

    with open(tree_path) as f:
        newick_str = f.read().strip()
    children, parents, root = parse_newick_iterative(newick_str)
    mutations = parse_mutations_file(mutations_path)

    n_tips = sum(1 for n in parents if n not in children or not children.get(n))
    n_internal = len(parents) + 1 - n_tips
    print(f"  Tree: {n_tips} tips, {n_internal} internal nodes, {len(mutations)} annotated")

    return children, parents, root, mutations


def analyze_mutations(mutations, reference_seq, trim_begin=None, trim_end=None):
    """Analyze mutations for reversions and indels.

    A reversion is a mutation where the derived base matches the reference.
    An indel involves N as the ancestral or derived base (UShER represents
    deletions/missing data as mutations to N).

    Returns dict with counts and per-branch distributions.
    """
    trim_start = (trim_begin - 1) if trim_begin else 0
    trim_stop = trim_end if trim_end else len(reference_seq)

    total_mutations = 0
    total_reversions = 0
    total_indels = 0
    total_branches = 0
    branches_with_reversions = 0
    branches_with_indels = 0

    per_branch_mutations = []
    per_branch_reversions = []
    per_branch_indels = []

    for node, muts in mutations.items():
        if not muts:
            continue

        branch_muts = 0
        branch_revs = 0
        branch_indels = 0

        for mut in muts:
            if len(mut) < 3:
                continue
            try:
                ancestral = mut[0]
                derived = mut[-1]
                pos = int(mut[1:-1])  # 1-indexed
            except (ValueError, IndexError):
                continue

            # Filter to trimmed region
            if pos < trim_start + 1 or pos > trim_stop:
                continue

            branch_muts += 1

            # Check for indel (N as ancestral or derived)
            if ancestral == "N" or derived == "N":
                branch_indels += 1
                continue

            # Check for reversion to reference
            ref_pos = pos - 1  # 0-indexed
            if ref_pos < len(reference_seq) and derived == reference_seq[ref_pos]:
                branch_revs += 1

        if branch_muts > 0:
            total_branches += 1
            total_mutations += branch_muts
            total_reversions += branch_revs
            total_indels += branch_indels
            per_branch_mutations.append(branch_muts)
            per_branch_reversions.append(branch_revs)
            per_branch_indels.append(branch_indels)
            if branch_revs > 0:
                branches_with_reversions += 1
            if branch_indels > 0:
                branches_with_indels += 1

    return {
        "total_mutations": total_mutations,
        "total_reversions": total_reversions,
        "total_indels": total_indels,
        "total_branches": total_branches,
        "branches_with_reversions": branches_with_reversions,
        "branches_with_indels": branches_with_indels,
        "reversion_rate": total_reversions / total_mutations if total_mutations else 0,
        "indel_rate": total_indels / total_mutations if total_mutations else 0,
        "per_branch_mutations": per_branch_mutations,
        "per_branch_reversions": per_branch_reversions,
        "per_branch_indels": per_branch_indels,
    }


def print_comparison(usher_stats, viridian_stats):
    """Print side-by-side comparison table."""
    def pct(n, d):
        return f"{100 * n / d:.2f}%" if d else "N/A"

    def mean(vals):
        return f"{sum(vals) / len(vals):.2f}" if vals else "N/A"

    print("\n" + "=" * 70)
    print(f"{'Metric':<35} {'UShER':>15} {'Viridian':>15}")
    print("=" * 70)

    rows = [
        ("Branches with mutations", usher_stats["total_branches"], viridian_stats["total_branches"]),
        ("Total mutations", usher_stats["total_mutations"], viridian_stats["total_mutations"]),
        ("Mean mutations/branch", mean(usher_stats["per_branch_mutations"]), mean(viridian_stats["per_branch_mutations"])),
    ]
    for label, u, v in rows:
        print(f"{label:<35} {str(u):>15} {str(v):>15}")

    print("-" * 70)

    print(f"{'Total reversions':<35} {usher_stats['total_reversions']:>15} {viridian_stats['total_reversions']:>15}")
    print(f"{'Reversion rate (% of mutations)':<35} {pct(usher_stats['total_reversions'], usher_stats['total_mutations']):>15} {pct(viridian_stats['total_reversions'], viridian_stats['total_mutations']):>15}")
    print(f"{'Branches with reversions':<35} {usher_stats['branches_with_reversions']:>15} {viridian_stats['branches_with_reversions']:>15}")
    print(f"{'% branches with reversions':<35} {pct(usher_stats['branches_with_reversions'], usher_stats['total_branches']):>15} {pct(viridian_stats['branches_with_reversions'], viridian_stats['total_branches']):>15}")
    print(f"{'Mean reversions/branch':<35} {mean(usher_stats['per_branch_reversions']):>15} {mean(viridian_stats['per_branch_reversions']):>15}")

    print("-" * 70)

    print(f"{'Total indels (N mutations)':<35} {usher_stats['total_indels']:>15} {viridian_stats['total_indels']:>15}")
    print(f"{'Indel rate (% of mutations)':<35} {pct(usher_stats['total_indels'], usher_stats['total_mutations']):>15} {pct(viridian_stats['total_indels'], viridian_stats['total_mutations']):>15}")
    print(f"{'Branches with indels':<35} {usher_stats['branches_with_indels']:>15} {viridian_stats['branches_with_indels']:>15}")
    print(f"{'% branches with indels':<35} {pct(usher_stats['branches_with_indels'], usher_stats['total_branches']):>15} {pct(viridian_stats['branches_with_indels'], viridian_stats['total_branches']):>15}")
    print(f"{'Mean indels/branch':<35} {mean(usher_stats['per_branch_indels']):>15} {mean(viridian_stats['per_branch_indels']):>15}")

    print("=" * 70)

    # Summary
    rev_change = viridian_stats["reversion_rate"] / usher_stats["reversion_rate"] if usher_stats["reversion_rate"] else float("inf")

    print(f"\nViridian reversion rate is {rev_change:.2f}x UShER ({('lower' if rev_change < 1 else 'higher')})")

    if usher_stats["indel_rate"] and viridian_stats["indel_rate"]:
        indel_change = viridian_stats["indel_rate"] / usher_stats["indel_rate"]
        print(f"Viridian indel rate is {indel_change:.2f}x UShER ({('lower' if indel_change < 1 else 'higher')})")
    elif usher_stats["total_indels"] == 0 and viridian_stats["total_indels"] == 0:
        print("No indels detected in either tree (UShER MAT format may not encode indels as nucleotide mutations)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--usher", required=True, help="Path to UShER .pb.gz file")
    parser.add_argument("--viridian", required=True, help="Path to Viridian .pb.gz file")
    parser.add_argument("--subsample", type=int, default=10000, help="Subsample to N tips (default: 10000)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--trim-begin", type=int, default=None, help="Start position for spike S1 trimming (1-indexed)")
    parser.add_argument("--trim-end", type=int, default=None, help="End position for spike S1 trimming (1-indexed)")
    parser.add_argument("--reference-cache", default=None, help="Path to cache reference FASTA")
    args = parser.parse_args()

    matutils_path = find_matutils()
    if not matutils_path:
        print("ERROR: matUtils not found. Install via: conda create -n usher -c conda-forge -c bioconda usher")
        sys.exit(1)
    print(f"Using matUtils: {matutils_path}")

    reference_seq = fetch_reference_sequence(args.reference_cache)
    print(f"Reference length: {len(reference_seq)}")
    if args.trim_begin and args.trim_end:
        print(f"Analyzing positions {args.trim_begin}-{args.trim_end} (spike S1)")

    # Process UShER tree
    print("\n--- UShER tree ---")
    with tempfile.TemporaryDirectory(prefix="usher_cmp_") as usher_dir:
        _, _, _, usher_mutations = extract_tree_data(
            args.usher, usher_dir, matutils_path,
            subsample=args.subsample, seed=args.seed,
        )
        usher_stats = analyze_mutations(
            usher_mutations, reference_seq,
            trim_begin=args.trim_begin, trim_end=args.trim_end,
        )

    # Process Viridian tree
    print("\n--- Viridian tree ---")
    with tempfile.TemporaryDirectory(prefix="viridian_cmp_") as viridian_dir:
        _, _, _, viridian_mutations = extract_tree_data(
            args.viridian, viridian_dir, matutils_path,
            subsample=args.subsample, seed=args.seed,
        )
        viridian_stats = analyze_mutations(
            viridian_mutations, reference_seq,
            trim_begin=args.trim_begin, trim_end=args.trim_end,
        )

    print_comparison(usher_stats, viridian_stats)


if __name__ == "__main__":
    main()
