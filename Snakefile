configfile: "defaults/config.yaml"

# Additional dataset configs (bac120-*.yaml, trellis.yaml, odb-fungi.yaml)
# are opt-in via --configfile defaults/{name}.yaml on the command line.

# Auto-discover rdrp subtrees (only if include_subtrees=true)
import glob
import os
import re
import shutil


def _file_gb(path):
    """Return file size in GB, or 0 if the file doesn't exist yet."""
    try:
        return os.path.getsize(path) / 1024**3
    except OSError:
        return 0


def _mem_gb(analysis, *filenames, multiplier=3, base=0):
    """Estimate mem_gb from input file sizes within an analysis data dir.

    Default multiplier of 3 covers most rules (Python dicts use ~2-3x raw
    file size). trajectory.py and pairwise need higher multipliers — they
    hold full alignment + branch data in memory simultaneously.
    """
    total = base + sum(_file_gb(f"data/{analysis}/{f}") for f in filenames) * multiplier
    return max(1, int(total) + 1)


def _alignment_mem_gb(analysis):
    """Estimate peak memory for alignment.py from JSON size and seq_length.

    alignment.py reconstructs full sequences for every tree node (~2x tips).
    Peak memory ≈ 3 * n_nodes * seq_length bytes (dict + mutable seqs).
    n_nodes estimated from JSON size (each node ≈ 300 bytes in compact JSON).
    """
    json_bytes = _file_gb(f"data/{analysis}/auspice_raw.json") * 1024**3
    if json_bytes == 0:
        return 1
    n_nodes_est = int(json_bytes / 300)
    seq_length = config["analysis"].get(analysis, {}).get("seq_length", 5000)
    peak_bytes = 3 * n_nodes_est * seq_length
    return max(1, int(peak_bytes / 1024**3) + 1)

# Helper function to detect UShER datasets
def is_usher_dataset(analysis_name):
    """Check if a dataset uses UShER protobuf format instead of Auspice JSON."""
    return "usher_pb" in config["analysis"].get(analysis_name, {})

if config.get("include_subtrees", False):
    RDRP_BASE = "../rdrp/phylogenetic/auspice"
    RDRP_FAMILIES = {
        "paramyxoviridae": 1653,
        "flaviviridae": 1884,
        "picornaviridae": 1386,
    }

    for family, seq_length in RDRP_FAMILIES.items():
        subtree_pattern = f"{RDRP_BASE}/{family}/subtrees/{family}_*.json"
        for json_path in sorted(glob.glob(subtree_pattern)):
            # Extract subtree ID (e.g., "001" from "paramyxoviridae_001.json")
            filename = os.path.basename(json_path)
            subtree_id = filename.replace(f"{family}_", "").replace(".json", "")
            analysis_name = f"rdrp-{family}-xs_{subtree_id}"

            # Skip if already in config (allows manual overrides)
            if analysis_name not in config["analysis"]:
                config["analysis"][analysis_name] = {
                    "dataset": json_path,
                    "gene": "nuc",
                    "seq_length": seq_length,
                }

# Get all analyses from config, or use target_analyses if specified on command line
ANALYSES = config.get("target_analyses", list(config["analysis"].keys()))

# Constrain {analysis} wildcard to valid analysis names only
wildcard_constraints:
    analysis = "|".join(config["analysis"].keys())

# Separate wildcard constraints for UShER vs Auspice datasets
# This ensures rules only match the appropriate dataset types
def _usher_analyses():
    return [a for a in config["analysis"].keys() if "usher_pb" in config["analysis"][a]]

def _auspice_analyses():
    return [a for a in config["analysis"].keys() if "usher_pb" not in config["analysis"][a]]

# Split analyses into UShER and non-UShER for different rule requirements
USHER_ANALYSES = [a for a in ANALYSES if is_usher_dataset(a)]
AUSPICE_ANALYSES = [a for a in ANALYSES if not is_usher_dataset(a)]
RDRP_ANALYSES = [a for a in ANALYSES if re.match(r'^rdrp-[a-z]+-xs$', a)]

rule all:
    input:
        # Auspice datasets have metadata and labeled JSON, UShER datasets don't
        expand("data/{analysis}/metadata.tsv", analysis=AUSPICE_ANALYSES),
        expand("data/{analysis}/auspice.json", analysis=AUSPICE_ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        expand("results/{analysis}/.trajectories.done", analysis=ANALYSES),
        expand("results/{analysis}/.pairwise.done", analysis=ANALYSES),

rule results:
    input:
        # Auspice datasets have metadata and labeled JSON, UShER datasets don't
        expand("data/{analysis}/metadata.tsv", analysis=AUSPICE_ANALYSES),
        expand("data/{analysis}/auspice.json", analysis=AUSPICE_ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        expand("results/{analysis}/.trajectories.done", analysis=ANALYSES),
        expand("results/{analysis}/.pairwise.done", analysis=ANALYSES),

rule download_auspice_json:
    wildcard_constraints:
        analysis = "|".join(_auspice_analyses()) if _auspice_analyses() else "NOMATCH"
    output:
        tree = "data/{analysis}/auspice_raw.json"
    params:
        dataset = lambda wildcards: config["analysis"][wildcards.analysis]["dataset"]
    shell:
        """
        mkdir -p data/{wildcards.analysis}
        dataset={params.dataset:q}
        if [[ "$dataset" == file://* ]]; then
            # Local file with file:// prefix - copy it
            cp "${{dataset#file://}}" {output.tree:q}
        elif [[ "$dataset" != http* && "$dataset" != s3://* && "$dataset" != *nextstrain.org* ]]; then
            # Relative or absolute local path - copy it
            cp "$dataset" {output.tree:q}
        else
            nextstrain remote download "$dataset" {output.tree:q}
        fi
        """

rule provision_alignment:
    wildcard_constraints:
        analysis = "|".join(_auspice_analyses()) if _auspice_analyses() else "NOMATCH"
    input:
        auspice = "data/{analysis}/auspice_raw.json"
    output:
        alignment = "data/{analysis}/alignment.fasta"
    resources:
        mem_gb = lambda wc: _alignment_mem_gb(wc.analysis)
    params:
        gene = lambda wildcards: config["analysis"][wildcards.analysis]["gene"],
        trim_begin_arg = lambda wildcards: f"--trim-begin {config['analysis'][wildcards.analysis]['trim_begin']}" if "trim_begin" in config['analysis'][wildcards.analysis] else "",
        trim_end_arg = lambda wildcards: f"--trim-end {config['analysis'][wildcards.analysis]['trim_end']}" if "trim_end" in config['analysis'][wildcards.analysis] else ""
    shell:
        """
        # Copy root-sequence sidecar file if it exists (some datasets have inline root sequence instead)
        if [ -f "data/{wildcards.analysis}/auspice_raw_root-sequence.json" ]; then
            cp "data/{wildcards.analysis}/auspice_raw_root-sequence.json" "data/{wildcards.analysis}/auspice_root-sequence.json"
        fi

        python scripts/alignment.py \
            --json {input.auspice:q} \
            --output {output.alignment:q} \
            --gene {params.gene:q} \
            {params.trim_begin_arg} \
            {params.trim_end_arg}
        """

rule provision_metadata:
    wildcard_constraints:
        analysis = "|".join(_auspice_analyses()) if _auspice_analyses() else "NOMATCH"
    input:
        auspice = "data/{analysis}/auspice_raw.json"
    output:
        metadata = "data/{analysis}/metadata.tsv"
    resources:
        mem_gb = lambda wc: _mem_gb(wc.analysis, "auspice_raw.json")
    shell:
        """
        python scripts/metadata.py \
            --json {input.auspice:q} \
            --output {output.metadata:q}
        """

rule provision_branches:
    wildcard_constraints:
        analysis = "|".join(_auspice_analyses()) if _auspice_analyses() else "NOMATCH"
    input:
        auspice = "data/{analysis}/auspice_raw.json",
        alignment = "data/{analysis}/alignment.fasta"
    output:
        branches = "data/{analysis}/branches_raw.tsv"
    resources:
        mem_gb = lambda wc: _mem_gb(wc.analysis, "auspice_raw.json", "alignment.fasta")
    shell:
        """
        python scripts/branches.py \
            --json {input.auspice:q} \
            --alignment {input.alignment:q} \
            --output {output.branches:q}
        """

rule provision_colors:
    wildcard_constraints:
        analysis = "|".join(_auspice_analyses()) if _auspice_analyses() else "NOMATCH"
    input:
        auspice = "data/{analysis}/auspice_raw.json"
    output:
        colors = "data/{analysis}/colors.json"
    shell:
        """
        python scripts/colors.py \
            --json {input.auspice:q} \
            --output {output.colors:q}
        """

rule train_test_split:
    """
    Unified train/test split for both Auspice and UShER datasets.
    Works on branches_raw.tsv and outputs branches.tsv with train_test labels.
    """
    input:
        branches_raw = "data/{analysis}/branches_raw.tsv"
    output:
        branches = "data/{analysis}/branches.tsv"
    resources:
        mem_gb = lambda wc: _mem_gb(wc.analysis, "branches_raw.tsv")
    params:
        test_proportion = lambda wildcards: config["analysis"][wildcards.analysis].get("test_proportion", 0.1),
        mutations_back = lambda wildcards: config["analysis"][wildcards.analysis].get("mutations_back", 5),
        max_clade_proportion = lambda wildcards: config["analysis"][wildcards.analysis].get("max_clade_proportion", 0.01),
        seed = lambda wildcards: config["analysis"][wildcards.analysis].get("seed", 42)
    shell:
        """
        python scripts/train_test_split.py \
            --branches-raw {input.branches_raw:q} \
            --output {output.branches:q} \
            --test-proportion {params.test_proportion} \
            --mutations-back {params.mutations_back} \
            --max-clade-proportion {params.max_clade_proportion} \
            --seed {params.seed}
        """

rule label_auspice_json:
    """
    Add train/test labels to Auspice JSON for visualization.
    """
    wildcard_constraints:
        analysis = "|".join(_auspice_analyses()) if _auspice_analyses() else "NOMATCH"
    input:
        auspice_raw = "data/{analysis}/auspice_raw.json",
        branches = "data/{analysis}/branches.tsv"
    output:
        auspice = "data/{analysis}/auspice.json"
    resources:
        mem_gb = lambda wc: _mem_gb(wc.analysis, "auspice_raw.json")
    shell:
        """
        python scripts/label_auspice_json.py \
            --json {input.auspice_raw:q} \
            --branches {input.branches:q} \
            --output {output.auspice:q}
        """

rule sample:
    input:
        alignment = "data/{analysis}/alignment.fasta"
    output:
        sampled = "data/{analysis}/sample.fasta"
    shell:
        """
        python scripts/sample.py \
            --input {input.alignment:q} \
            --output {output.sampled:q} \
            --fraction 0.2
        """

rule trajectories:
    input:
        branches = "data/{analysis}/branches.tsv",
        alignment = "data/{analysis}/alignment.fasta"
    output:
        done = "results/{analysis}/.trajectories.done"
    resources:
        mem_gb = lambda wc: _mem_gb(wc.analysis, "alignment.fasta", multiplier=10, base=5)
    params:
        output_dir = "results/{analysis}",
        summary = "results/summary.json",
        url = lambda wildcards: config["analysis"][wildcards.analysis].get("dataset", config["analysis"][wildcards.analysis].get("usher_pb", "")),
        shard_size = lambda wildcards: config["analysis"][wildcards.analysis].get("shard_size", 10000),
        seed = lambda wildcards: config["analysis"][wildcards.analysis].get("seed", 42)
    shell:
        """
        python scripts/trajectory.py \
            --branches {input.branches:q} \
            --alignment {input.alignment:q} \
            --output-dir {params.output_dir:q} \
            --shard-size {params.shard_size} \
            --seed {params.seed} \
            --summary {params.summary:q} \
            --dataset {wildcards.analysis} \
            --url {params.url:q}
        touch {output.done}
        """

rule pairwise:
    input:
        branches = "data/{analysis}/branches.tsv",
        alignment = "data/{analysis}/alignment.fasta"
    output:
        done = "results/{analysis}/.pairwise.done"
    resources:
        mem_gb = lambda wc: _mem_gb(wc.analysis, "alignment.fasta", multiplier=6)
    params:
        output_dir = "results/{analysis}",
        train_limit = lambda wildcards: config["analysis"][wildcards.analysis].get("pairwise_train_limit", 100000),
        test_limit = lambda wildcards: config["analysis"][wildcards.analysis].get("pairwise_test_limit", 50000),
        shard_size = lambda wildcards: config["analysis"][wildcards.analysis].get("shard_size", 10000),
        seed = lambda wildcards: config["analysis"][wildcards.analysis].get("seed", 42),
        summary = "results/summary.json",
        url = lambda wildcards: config["analysis"][wildcards.analysis].get("dataset", config["analysis"][wildcards.analysis].get("usher_pb", ""))
    shell:
        """
        python scripts/pairwise_trajectory.py \
            --branches {input.branches:q} \
            --alignment {input.alignment:q} \
            --output-dir {params.output_dir:q} \
            --train-limit {params.train_limit} \
            --test-limit {params.test_limit} \
            --shard-size {params.shard_size} \
            --seed {params.seed} \
            --summary {params.summary:q} \
            --dataset {wildcards.analysis} \
            --url {params.url:q}
        touch {output.done}
        """

# ============================================================================
# UShER-specific rules for processing mutation-annotated trees (protobuf)
# ============================================================================

rule download_usher_pb:
    """Download UShER protobuf file from UCSC."""
    wildcard_constraints:
        analysis = "|".join(_usher_analyses()) if _usher_analyses() else "NOMATCH"
    output:
        pb = "data/{analysis}/tree.pb.gz"
    params:
        url = lambda wildcards: config["analysis"][wildcards.analysis]["usher_pb"]
    shell:
        """
        mkdir -p data/{wildcards.analysis}
        curl -L -o {output.pb:q} {params.url:q}
        """

rule usher_provision_alignment_branches:
    """
    Extract tree and mutations from UShER protobuf, reconstruct sequences.

    This rule:
    1. Runs matUtils extract to get tree.nwk and mutations
    2. Fetches SARS-CoV-2 reference from NCBI
    3. Reconstructs sequences by applying mutations from root
    4. Outputs alignment.fasta (trimmed) and branches_raw.tsv
    """
    wildcard_constraints:
        analysis = "|".join(_usher_analyses()) if _usher_analyses() else "NOMATCH"
    input:
        pb = "data/{analysis}/tree.pb.gz"
    output:
        alignment = "data/{analysis}/alignment.fasta",
        branches_raw = "data/{analysis}/branches_raw.tsv",
        tree = "data/{analysis}/tree.nwk"
    params:
        subsample = lambda wildcards: config["analysis"][wildcards.analysis].get("subsample"),
        seed = lambda wildcards: config["analysis"][wildcards.analysis].get("seed", 42),
        trim_begin = lambda wildcards: config["analysis"][wildcards.analysis].get("trim_begin"),
        trim_end = lambda wildcards: config["analysis"][wildcards.analysis].get("trim_end"),
        subsample_arg = lambda wildcards: f"--subsample {config['analysis'][wildcards.analysis]['subsample']}" if config['analysis'][wildcards.analysis].get('subsample') else "",
        trim_begin_arg = lambda wildcards: f"--trim-begin {config['analysis'][wildcards.analysis]['trim_begin']}" if config['analysis'][wildcards.analysis].get('trim_begin') else "",
        trim_end_arg = lambda wildcards: f"--trim-end {config['analysis'][wildcards.analysis]['trim_end']}" if config['analysis'][wildcards.analysis].get('trim_end') else "",
        reference_cache = "data/{analysis}/reference.fasta"
    shell:
        """
        python scripts/usher_provision_alignment_branches.py \
            --pb {input.pb:q} \
            --output-fasta {output.alignment:q} \
            --output-branches {output.branches_raw:q} \
            --output-tree {output.tree:q} \
            --reference-cache {params.reference_cache:q} \
            --seed {params.seed} \
            {params.subsample_arg} \
            {params.trim_begin_arg} \
            {params.trim_end_arg}
        """

# ============================================================================
# End of UShER-specific rules
# ============================================================================

rule upload:
    input:
        expand("results/{analysis}/.trajectories.done", analysis=ANALYSES),
        expand("results/{analysis}/.pairwise.done", analysis=ANALYSES)
    params:
        analyses = ANALYSES,
        s3_prefix = config.get("s3_prefix", "trajectories")
    shell:
        """
        python scripts/upload-to-s3.py \
            --upload-dir results \
            --prefix {params.s3_prefix} \
            --analyses {params.analyses}
        """

rule upload_rdrp:
    input:
        expand("results/{analysis}/.trajectories.done", analysis=RDRP_ANALYSES),
        expand("results/{analysis}/.pairwise.done", analysis=RDRP_ANALYSES)
    params:
        analyses = RDRP_ANALYSES,
        s3_prefix = config.get("s3_prefix", "trajectories")
    shell:
        """
        python scripts/upload-to-s3.py \
            --upload-dir results \
            --prefix {params.s3_prefix} \
            --analyses {params.analyses}
        """
