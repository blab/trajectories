configfile: "defaults/config.yaml"

import os


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


# Get all analyses from config, or use target_analyses if specified on command line
ANALYSES = config.get("target_analyses", list(config["analysis"].keys()))

# Constrain {analysis} wildcard to valid analysis names only
wildcard_constraints:
    analysis = "|".join(config["analysis"].keys()) or "NOMATCH"

# Gate which trajectory pipelines (forwards, pairwise, or both) are required
# by the default consumer rules (all, results).
TRAJECTORY_MODE = config.get("trajectory_mode", "both")
if TRAJECTORY_MODE not in ("forwards", "pairwise", "both"):
    raise ValueError(
        f"trajectory_mode must be 'forwards', 'pairwise', or 'both'; got {TRAJECTORY_MODE!r}"
    )

def trajectory_targets(analyses):
    targets = []
    if TRAJECTORY_MODE in ("forwards", "both"):
        targets += expand("results/{analysis}/.trajectories.done", analysis=analyses)
    if TRAJECTORY_MODE in ("pairwise", "both"):
        targets += expand("results/{analysis}/.pairwise.done", analysis=analyses)
    return targets

rule all:
    input:
        expand("data/{analysis}/metadata.tsv", analysis=ANALYSES),
        expand("data/{analysis}/auspice.json", analysis=ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        trajectory_targets(ANALYSES),

rule results:
    input:
        expand("data/{analysis}/metadata.tsv", analysis=ANALYSES),
        expand("data/{analysis}/auspice.json", analysis=ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        trajectory_targets(ANALYSES),

rule download_auspice_json:
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

rule train_test_split:
    """
    Split branches into train/test by clade selection, output branches.tsv.
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
        url = lambda wildcards: config["analysis"][wildcards.analysis]["dataset"],
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
        url = lambda wildcards: config["analysis"][wildcards.analysis]["dataset"]
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
