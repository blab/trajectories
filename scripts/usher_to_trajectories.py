#!/usr/bin/env python3
"""
Convert UShER mutation-annotated tree (protobuf) to trajectory pipeline format.

This script:
1. Runs matUtils extract to get tree structure and mutations
2. Fetches SARS-CoV-2 reference sequence from NCBI
3. Reconstructs sequences for all nodes by applying mutations from root
4. Outputs alignment.fasta and branches_raw.tsv

The output files are compatible with the downstream trajectory pipeline.

Generated with Claude Code assistance.
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

import requests
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


# NCBI reference for SARS-CoV-2 (Wuhan-Hu-1)
REFERENCE_ACCESSION = "NC_045512.2"
REFERENCE_URL = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={REFERENCE_ACCESSION}&rettype=fasta&retmode=text"


def fetch_reference_sequence(cache_path=None):
    """Fetch SARS-CoV-2 reference sequence from NCBI."""
    # Check cache first
    if cache_path and os.path.exists(cache_path):
        print(f"Loading cached reference from {cache_path}")
        record = SeqIO.read(cache_path, "fasta")
        return str(record.seq).upper()

    print(f"Fetching reference sequence {REFERENCE_ACCESSION} from NCBI...")
    response = requests.get(REFERENCE_URL, timeout=30)
    response.raise_for_status()

    # Parse FASTA from response
    lines = response.text.strip().split("\n")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))

    # Cache if path provided
    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "w") as f:
            f.write(f">{REFERENCE_ACCESSION}\n{sequence}\n")
        print(f"Cached reference to {cache_path}")

    return sequence.upper()


def find_matutils():
    """Find matUtils executable, checking common locations."""
    # Check if matUtils is in PATH
    try:
        result = subprocess.run(
            ["which", "matUtils"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except FileNotFoundError:
        pass

    # Check common conda locations (without requiring conda activation)
    conda_paths = [
        "/opt/homebrew/Caskroom/miniconda/base/envs/usher/bin/matUtils",
        os.path.expanduser("~/miniconda3/envs/usher/bin/matUtils"),
        os.path.expanduser("~/anaconda3/envs/usher/bin/matUtils"),
        os.path.expanduser("~/miniforge3/envs/usher/bin/matUtils"),
    ]

    for path in conda_paths:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path

    return None


def check_docker_available():
    """Check if Docker is available and running."""
    try:
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def run_matutils_extract(pb_file, output_dir):
    """Run matUtils extract to get tree and mutations.

    Tries local matUtils first (including conda envs), falls back to Docker if not found.
    """
    # Use absolute paths for tracking, but run matUtils from output_dir with relative paths
    # (matUtils has a bug where it prepends cwd to absolute paths)
    abs_output_dir = os.path.abspath(output_dir)
    abs_pb_file = os.path.abspath(pb_file)
    tree_path = os.path.join(abs_output_dir, "tree.nwk")
    mutations_path = os.path.join(abs_output_dir, "all_mutations.txt")

    # Find matUtils executable
    matutils_path = find_matutils()

    if matutils_path:
        print(f"Running matUtils extract using: {matutils_path}")
        # Run from output directory with relative output paths to work around matUtils path bug
        cmd = [
            matutils_path, "extract",
            "-i", abs_pb_file,  # Input can be absolute
            "-t", "tree.nwk",   # Output must be relative
            "-A", "all_mutations.txt"
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=abs_output_dir)
            return tree_path, mutations_path
        except subprocess.CalledProcessError as e:
            print(f"ERROR: matUtils failed: {e.stderr}")
            sys.exit(1)

    # Fall back to Docker
    print("Local matUtils not found, trying Docker...")
    if not check_docker_available():
        print("ERROR: Neither matUtils nor Docker is available.")
        print("Please install one of:")
        print("  - UShER via conda: conda create -n usher -c conda-forge -c bioconda usher")
        print("  - Docker: docker pull pathogengenomics/usher:latest")
        sys.exit(1)

    # Use absolute paths for Docker volume mounting
    abs_pb_file = os.path.abspath(pb_file)
    abs_output_dir = os.path.abspath(output_dir)

    # Docker command with volume mounts
    # Use --workdir / to ensure absolute paths are interpreted correctly
    docker_cmd = [
        "docker", "run", "--rm",
        "--workdir", "/",
        "-v", f"{os.path.dirname(abs_pb_file)}:/input:ro",
        "-v", f"{abs_output_dir}:/output",
        "pathogengenomics/usher:latest",
        "matUtils", "extract",
        "-i", f"/input/{os.path.basename(abs_pb_file)}",
        "-t", "/output/tree.nwk",
        "-A", "/output/all_mutations.txt"
    ]

    print(f"Running via Docker: pathogengenomics/usher:latest")
    try:
        result = subprocess.run(docker_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Docker matUtils failed: {e.stderr}")
        sys.exit(1)

    return tree_path, mutations_path


def parse_newick_iterative(newick_str):
    """
    Parse Newick string iteratively to avoid stack overflow on large trees.

    Returns:
        dict mapping node_name -> list of child names
        dict mapping node_name -> parent name
        str: root node name
    """
    children = defaultdict(list)
    parents = {}

    # Tokenize: split on special characters but keep them
    tokens = re.findall(r"[^(),;:]+|[(),;:]", newick_str)

    # State machine parsing
    stack = []  # Stack of (parent_name, current_children)
    current_children = []
    node_counter = 0

    i = 0
    while i < len(tokens):
        token = tokens[i]

        if token == "(":
            # Start new clade - push current state
            stack.append((None, current_children))
            current_children = []
            i += 1

        elif token == ",":
            i += 1

        elif token == ")":
            # End clade - get name (might be empty for internal nodes)
            clade_children = current_children

            # Look ahead for name and branch length
            i += 1
            name = ""
            while i < len(tokens) and tokens[i] not in "(),;":
                if tokens[i] == ":":
                    # Skip branch length
                    i += 1
                    if i < len(tokens) and tokens[i] not in "(),;":
                        i += 1  # Skip the actual length value
                else:
                    name = tokens[i]
                    i += 1

            # Generate name for unnamed internal nodes
            if not name:
                name = f"node_{node_counter}"
                node_counter += 1

            # Set up parent-child relationships
            for child in clade_children:
                children[name].append(child)
                parents[child] = name

            # Pop state and add this node as child of parent
            _, current_children = stack.pop()
            current_children.append(name)

        elif token == ";":
            i += 1

        elif token == ":":
            # Branch length follows - skip it
            i += 1
            if i < len(tokens) and tokens[i] not in "(),;":
                i += 1

        else:
            # Leaf node name
            name = token.strip()
            if name:
                current_children.append(name)
                # Check for branch length
                i += 1
                if i < len(tokens) and tokens[i] == ":":
                    i += 1
                    if i < len(tokens) and tokens[i] not in "(),;":
                        i += 1
            else:
                i += 1

    # Root is the last remaining node in current_children
    root = current_children[0] if current_children else None

    return dict(children), parents, root


def parse_mutations_file(mutations_path):
    """
    Parse matUtils mutations file.

    Format: node_name: mut1,mut2,mut3...
    (colon and space separated, not tab)

    Mutations are like: A123T (nucleotide) or S:D614G (amino acid)

    Returns dict mapping node_name -> list of nucleotide mutations
    """
    mutations = {}

    print(f"Parsing mutations from {mutations_path}...")
    with open(mutations_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Format is "node_name: mut1,mut2,..." with colon-space separator
            if ": " in line:
                node_name, mut_str = line.split(": ", 1)
            elif line.endswith(":"):
                # Handle "node_name:" with no mutations
                node_name = line[:-1]
                mut_str = ""
            else:
                # Fallback: treat whole line as node name with no mutations
                node_name = line
                mut_str = ""

            # Parse mutations, keeping only nucleotide mutations
            # Nucleotide mutations look like: A123T, G456C
            # Skip amino acid mutations which look like: S:D614G
            node_muts = []
            if mut_str:
                for mut in mut_str.split(","):
                    mut = mut.strip()
                    # Skip amino acid mutations (contain colon like S:D614G)
                    if ":" in mut:
                        continue
                    # Valid nucleotide mutation: starts with base, ends with base
                    if mut and mut[0] in "ACGTN" and mut[-1] in "ACGTN":
                        node_muts.append(mut)

            mutations[node_name] = node_muts

    return mutations


def apply_mutations(sequence, mutations):
    """Apply list of mutations to a sequence, returning new sequence."""
    seq = MutableSeq(sequence)

    for mut in mutations:
        # Parse mutation like A123T
        if len(mut) < 3:
            continue
        try:
            pos = int(mut[1:-1]) - 1  # Convert to 0-indexed
            new_base = mut[-1]
            if 0 <= pos < len(seq):
                seq[pos] = new_base
        except (ValueError, IndexError):
            continue

    return str(seq)


def calculate_hamming(seq1, seq2):
    """Calculate Hamming distance between two sequences, ignoring N and gaps."""
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]

    distance = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 in "-N" or c2 in "-N":
            continue
        if c1 != c2:
            distance += 1
    return distance


def extract_year_from_name(name):
    """Extract year from a sequence name, returning None if not found.

    Looks for SARS-CoV-2 date patterns. Valid years are 2019-2026.
    """
    import re
    # Primary: Look for full date pattern YYYY-MM-DD (most reliable)
    match = re.search(r'(20(?:19|2[0-6]))-\d{2}-\d{2}', name)
    if match:
        return int(match.group(1))
    # Secondary: Look for /YYYY| pattern (e.g., "England/sample/2020|accession")
    match = re.search(r'/(20(?:19|2[0-6]))\|', name)
    if match:
        return int(match.group(1))
    # Tertiary: Look for /YYYY at end before any suffix
    match = re.search(r'/(20(?:19|2[0-6]))(?:[|\s]|$)', name)
    if match:
        return int(match.group(1))
    return None


def subsample_tips_stratified(tips, n_samples, rng):
    """
    Subsample tips with stratification by year for temporal diversity.

    Args:
        tips: list of tip names
        n_samples: target number of samples
        rng: random.Random instance

    Returns:
        set of sampled tip names
    """
    from collections import defaultdict

    # Group tips by year
    tips_by_year = defaultdict(list)
    tips_no_year = []

    for tip in tips:
        year = extract_year_from_name(tip)
        if year:
            tips_by_year[year].append(tip)
        else:
            tips_no_year.append(tip)

    years = sorted(tips_by_year.keys())

    if not years:
        # No years found, fall back to random sampling
        print("  No date information found, using random sampling")
        return set(rng.sample(tips, min(n_samples, len(tips))))

    # Report year distribution before sampling
    print(f"  Tips by year before sampling:")
    for year in years:
        print(f"    {year}: {len(tips_by_year[year])}")
    if tips_no_year:
        print(f"    unknown: {len(tips_no_year)}")

    # Calculate samples per year (equal distribution)
    # Reserve 10% for tips without dates if any exist
    if tips_no_year:
        n_for_dated = int(n_samples * 0.9)
        n_for_undated = n_samples - n_for_dated
    else:
        n_for_dated = n_samples
        n_for_undated = 0

    samples_per_year = n_for_dated // len(years)
    extra_samples = n_for_dated % len(years)

    sampled_tips = set()

    # Sample from each year
    for i, year in enumerate(years):
        year_tips = tips_by_year[year]
        # Distribute extra samples to earlier years (more diverse)
        n_year = samples_per_year + (1 if i < extra_samples else 0)
        n_year = min(n_year, len(year_tips))

        if n_year > 0:
            sampled_tips.update(rng.sample(year_tips, n_year))

    # Sample from undated tips
    if tips_no_year and n_for_undated > 0:
        n_undated = min(n_for_undated, len(tips_no_year))
        sampled_tips.update(rng.sample(tips_no_year, n_undated))

    # If we didn't get enough samples (some years had fewer tips), sample more from larger years
    if len(sampled_tips) < n_samples:
        remaining = n_samples - len(sampled_tips)
        available = [t for t in tips if t not in sampled_tips]
        if available:
            extra = rng.sample(available, min(remaining, len(available)))
            sampled_tips.update(extra)

    # Report year distribution after sampling
    sampled_by_year = defaultdict(int)
    sampled_no_year = 0
    for tip in sampled_tips:
        year = extract_year_from_name(tip)
        if year:
            sampled_by_year[year] += 1
        else:
            sampled_no_year += 1

    print(f"  Tips by year after sampling ({len(sampled_tips)} total):")
    for year in sorted(sampled_by_year.keys()):
        pct = 100 * sampled_by_year[year] / len(sampled_tips)
        print(f"    {year}: {sampled_by_year[year]} ({pct:.1f}%)")
    if sampled_no_year:
        pct = 100 * sampled_no_year / len(sampled_tips)
        print(f"    unknown: {sampled_no_year} ({pct:.1f}%)")

    return sampled_tips


def reconstruct_sequences_streaming(
    reference_seq, children, parents, root, mutations,
    output_fasta, subsample=None, seed=42, trim_begin=None, trim_end=None
):
    """
    Reconstruct sequences by depth-first traversal, streaming to disk.

    Returns dict mapping node_name -> sequence (trimmed if specified).
    For memory efficiency with large trees, we store sequences in a dict
    but write to disk as we go.
    """
    import random
    rng = random.Random(seed)

    # Get all node names
    all_nodes = set(mutations.keys())
    if root and root not in all_nodes:
        all_nodes.add(root)
    for parent in parents.values():
        all_nodes.add(parent)
    for child_list in children.values():
        for child in child_list:
            all_nodes.add(child)

    # Subsample if requested
    if subsample and subsample < len(all_nodes):
        print(f"Subsampling {subsample} tips from {len(all_nodes)} total nodes...")
        # We need to keep tree structure intact, so subsample tips and include their ancestors
        tips = [n for n in all_nodes if n not in children or not children[n]]
        print(f"  Found {len(tips)} tips")

        if subsample < len(tips):
            # Use stratified sampling by year for better temporal diversity
            sampled_tips = subsample_tips_stratified(tips, subsample, rng)
        else:
            sampled_tips = set(tips)

        # Include all ancestors of sampled tips
        nodes_to_keep = set(sampled_tips)
        for tip in sampled_tips:
            current = tip
            while current in parents:
                nodes_to_keep.add(current)
                current = parents[current]
            nodes_to_keep.add(current)  # Add root

        print(f"After including ancestors: {len(nodes_to_keep)} nodes")
    else:
        nodes_to_keep = all_nodes

    # Compute trimming indices
    trim_start = (trim_begin - 1) if trim_begin else 0
    trim_stop = trim_end if trim_end else len(reference_seq)

    sequences = {}

    # Iterative depth-first traversal using explicit stack
    # Stack contains (node_name, parent_sequence)
    root_seq = reference_seq
    if root in mutations:
        root_seq = apply_mutations(reference_seq, mutations[root])

    stack = [(root, root_seq)]
    processed = set()

    print("Reconstructing sequences...")
    with open(output_fasta, "w") as fasta_out:
        with tqdm(total=len(nodes_to_keep), desc="Nodes") as pbar:
            while stack:
                node, parent_seq = stack.pop()

                if node in processed:
                    continue
                processed.add(node)

                if node not in nodes_to_keep:
                    continue

                # Apply this node's mutations to parent sequence
                node_muts = mutations.get(node, [])
                if node == root:
                    node_seq = parent_seq
                else:
                    node_seq = apply_mutations(parent_seq, node_muts)

                # Trim sequence
                trimmed_seq = node_seq[trim_start:trim_stop]

                # Store and write
                sequences[node] = trimmed_seq
                fasta_out.write(f">{node}\n{trimmed_seq}\n")
                pbar.update(1)

                # Add children to stack
                for child in children.get(node, []):
                    if child not in processed and child in nodes_to_keep:
                        stack.append((child, node_seq))

    return sequences


def write_branches_tsv(sequences, parents, output_path):
    """Write branches.tsv with parent-child relationships and Hamming distances."""
    print(f"Writing branches to {output_path}...")

    with open(output_path, "w") as f:
        f.write("parent\tchild\thamming\ttrain_test\n")

        for child, parent in tqdm(parents.items(), desc="Branches"):
            if child not in sequences or parent not in sequences:
                continue

            hamming = calculate_hamming(sequences[parent], sequences[child])
            # train_test will be filled in by usher_train_test_split.py
            f.write(f"{parent}\t{child}\t{hamming}\t\n")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pb", required=True,
        help="Path to UShER protobuf file (.pb or .pb.gz)"
    )
    parser.add_argument(
        "--output-fasta", required=True,
        help="Output FASTA file for reconstructed sequences"
    )
    parser.add_argument(
        "--output-branches", required=True,
        help="Output TSV file for branch information"
    )
    parser.add_argument(
        "--output-tree", required=True,
        help="Output Newick tree file"
    )
    parser.add_argument(
        "--reference-cache",
        help="Path to cache reference sequence (optional)"
    )
    parser.add_argument(
        "--subsample", type=int, default=None,
        help="Subsample to this many tips (for testing)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for subsampling"
    )
    parser.add_argument(
        "--trim-begin", type=int, default=None,
        help="Start position for trimming (1-indexed, inclusive)"
    )
    parser.add_argument(
        "--trim-end", type=int, default=None,
        help="End position for trimming (1-indexed, inclusive)"
    )
    parser.add_argument(
        "--work-dir",
        help="Working directory for intermediate files (default: temp dir)"
    )

    args = parser.parse_args()

    # Set up working directory
    if args.work_dir:
        work_dir = args.work_dir
        os.makedirs(work_dir, exist_ok=True)
        cleanup = False
    else:
        work_dir = tempfile.mkdtemp(prefix="usher_")
        cleanup = True

    try:
        # Decompress if needed
        pb_file = args.pb
        if args.pb.endswith(".gz"):
            print(f"Decompressing {args.pb}...")
            decompressed = os.path.join(work_dir, "tree.pb")
            with gzip.open(args.pb, "rb") as f_in:
                with open(decompressed, "wb") as f_out:
                    # Read in chunks to handle large files
                    while True:
                        chunk = f_in.read(1024 * 1024)  # 1MB chunks
                        if not chunk:
                            break
                        f_out.write(chunk)
            pb_file = decompressed

        # Run matUtils extract
        tree_path, mutations_path = run_matutils_extract(pb_file, work_dir)

        # Copy tree to output location
        import shutil
        shutil.copy(tree_path, args.output_tree)
        print(f"Tree written to {args.output_tree}")

        # Fetch reference sequence
        reference_seq = fetch_reference_sequence(args.reference_cache)
        print(f"Reference sequence length: {len(reference_seq)}")

        # Parse tree structure
        print(f"Parsing Newick tree from {tree_path}...")
        with open(tree_path, "r") as f:
            newick_str = f.read().strip()

        children, parents, root = parse_newick_iterative(newick_str)
        print(f"Tree has {len(parents) + 1} nodes, root: {root}")

        # Parse mutations
        mutations = parse_mutations_file(mutations_path)
        print(f"Parsed mutations for {len(mutations)} nodes")

        # Reconstruct sequences
        sequences = reconstruct_sequences_streaming(
            reference_seq, children, parents, root, mutations,
            args.output_fasta,
            subsample=args.subsample,
            seed=args.seed,
            trim_begin=args.trim_begin,
            trim_end=args.trim_end
        )
        print(f"Wrote {len(sequences)} sequences to {args.output_fasta}")

        # Filter parents to only include nodes we kept
        filtered_parents = {c: p for c, p in parents.items()
                          if c in sequences and p in sequences}

        # Write branches
        write_branches_tsv(sequences, filtered_parents, args.output_branches)

        print("Done!")

    finally:
        # Cleanup temp directory if we created it
        if cleanup and os.path.exists(work_dir):
            import shutil
            shutil.rmtree(work_dir)


if __name__ == "__main__":
    main()
