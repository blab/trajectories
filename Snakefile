configfile: "defaults/config.yaml"

# Auto-discover rdrp subtrees
import glob
import os

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
        analysis_name = f"rdrp-{family}_{subtree_id}"

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

rule all:
    input:
        expand("data/{analysis}/metadata.tsv", analysis=ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        expand("results/{analysis}/forwards-train", analysis=ANALYSES),
        expand("results/{analysis}/forwards-test", analysis=ANALYSES),
        expand("results/{analysis}/pairwise-train", analysis=ANALYSES),
        expand("results/{analysis}/pairwise-test", analysis=ANALYSES),

rule results:
    input:
        expand("data/{analysis}/metadata.tsv", analysis=ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        expand("results/{analysis}/forwards-train", analysis=ANALYSES),
        expand("results/{analysis}/forwards-test", analysis=ANALYSES),
        expand("results/{analysis}/pairwise-train", analysis=ANALYSES),
        expand("results/{analysis}/pairwise-test", analysis=ANALYSES),

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

rule train_test_split:
    input:
        auspice = "data/{analysis}/auspice_raw.json"
    output:
        auspice = "data/{analysis}/auspice.json"
    params:
        gene = lambda wildcards: config["analysis"][wildcards.analysis]["gene"],
        test_proportion = lambda wildcards: config["analysis"][wildcards.analysis].get("test_proportion", 0.1),
        mutations_back = lambda wildcards: config["analysis"][wildcards.analysis].get("mutations_back", 5),
        max_clade_proportion = lambda wildcards: config["analysis"][wildcards.analysis].get("max_clade_proportion", 0.01),
        seed = lambda wildcards: config["analysis"][wildcards.analysis].get("seed", 42),
        trim_begin_arg = lambda wildcards: f"--trim-begin {config['analysis'][wildcards.analysis]['trim_begin']}" if "trim_begin" in config['analysis'][wildcards.analysis] else "",
        trim_end_arg = lambda wildcards: f"--trim-end {config['analysis'][wildcards.analysis]['trim_end']}" if "trim_end" in config['analysis'][wildcards.analysis] else ""
    shell:
        """
        python scripts/train_test_split.py \
            --json {input.auspice:q} \
            --output {output.auspice:q} \
            --gene {params.gene:q} \
            --test-proportion {params.test_proportion} \
            --mutations-back {params.mutations_back} \
            --max-clade-proportion {params.max_clade_proportion} \
            --seed {params.seed} \
            {params.trim_begin_arg} \
            {params.trim_end_arg}
        """

rule provision_alignment:
    input:
        auspice = "data/{analysis}/auspice.json"
    output:
        alignment = "data/{analysis}/full_alignment.fasta"
    params:
        gene = lambda wildcards: config["analysis"][wildcards.analysis]["gene"]
    shell:
        """
        # Copy root-sequence sidecar file if it exists (some datasets have inline root sequence instead)
        if [ -f "data/{wildcards.analysis}/auspice_raw_root-sequence.json" ]; then
            cp "data/{wildcards.analysis}/auspice_raw_root-sequence.json" "data/{wildcards.analysis}/auspice_root-sequence.json"
        fi

        python scripts/alignment.py \
            --json {input.auspice:q} \
            --output {output.alignment:q} \
            --gene {params.gene:q}
        """

rule provision_metadata:
    input:
        auspice = "data/{analysis}/auspice.json"
    output:
        metadata = "data/{analysis}/metadata.tsv"
    shell:
        """
        python scripts/metadata.py \
            --json {input.auspice:q} \
            --output {output.metadata:q}
        """

rule provision_branches:
    input:
        auspice = "data/{analysis}/auspice.json",
        alignment = "data/{analysis}/alignment.fasta"
    output:
        branches = "data/{analysis}/branches.tsv"
    shell:
        """
        python scripts/branches.py \
            --json {input.auspice:q} \
            --alignment {input.alignment:q} \
            --output {output.branches:q}
        """

rule provision_colors:
    input:
        auspice = "data/{analysis}/auspice.json"
    output:
        colors = "data/{analysis}/colors.json"
    shell:
        """
        python scripts/colors.py \
            --json {input.auspice:q} \
            --output {output.colors:q}
        """

rule trim:
    input:
        alignment = "data/{analysis}/full_alignment.fasta"
    output:
        alignment = "data/{analysis}/alignment.fasta"
    params:
        begin = lambda wildcards: config["analysis"][wildcards.analysis].get("trim_begin"),
        end = lambda wildcards: config["analysis"][wildcards.analysis].get("trim_end"),
        begin_arg = lambda wildcards: f"--begin {config['analysis'][wildcards.analysis]['trim_begin']}" if "trim_begin" in config['analysis'][wildcards.analysis] else "",
        end_arg = lambda wildcards: f"--end {config['analysis'][wildcards.analysis]['trim_end']}" if "trim_end" in config['analysis'][wildcards.analysis] else ""
    shell:
        """
        python scripts/trim.py \
            --input-alignment {input.alignment:q} \
            --output-alignment {output.alignment:q} \
            {params.begin_arg} \
            {params.end_arg}
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
        train_dir = directory("results/{analysis}/forwards-train"),
        test_dir = directory("results/{analysis}/forwards-test")
    params:
        output_dir = "results/{analysis}",
        summary = "results/summary.json",
        url = lambda wildcards: config["analysis"][wildcards.analysis]["dataset"]
    shell:
        """
        python scripts/trajectory.py \
            --branches {input.branches:q} \
            --alignment {input.alignment:q} \
            --output-dir {params.output_dir:q} \
            --summary {params.summary:q} \
            --dataset {wildcards.analysis} \
            --url {params.url:q}
        """

rule pairwise:
    input:
        branches = "data/{analysis}/branches.tsv",
        alignment = "data/{analysis}/alignment.fasta"
    output:
        train_dir = directory("results/{analysis}/pairwise-train"),
        test_dir = directory("results/{analysis}/pairwise-test")
    params:
        output_dir = "results/{analysis}",
        train_limit = lambda wildcards: config["analysis"][wildcards.analysis].get("pairwise_train_limit", 100000),
        test_limit = lambda wildcards: config["analysis"][wildcards.analysis].get("pairwise_test_limit", 50000),
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
            --seed {params.seed} \
            --summary {params.summary:q} \
            --dataset {wildcards.analysis} \
            --url {params.url:q}
        """

rule package:
    input:
        forwards_train = "results/{analysis}/forwards-train",
        forwards_test = "results/{analysis}/forwards-test",
        pairwise_train = "results/{analysis}/pairwise-train",
        pairwise_test = "results/{analysis}/pairwise-test"
    output:
        sharddir = directory("export/{analysis}")
    params:
        input_dir = "results/{analysis}"
    shell:
        """
        python scripts/package.py \
            --input-dir {params.input_dir:q} \
            --output-dir {output.sharddir:q} \
            --shuffle
        """

rule copy_summary:
    input:
        forwards = expand("results/{analysis}/forwards-train", analysis=ANALYSES),
        pairwise = expand("results/{analysis}/pairwise-train", analysis=ANALYSES)
    output:
        "export/summary.json"
    shell:
        "cp results/summary.json {output}"

rule export:
    input:
        expand("export/{analysis}", analysis=ANALYSES),
        "export/summary.json"

rule upload:
    input:
        expand("export/{analysis}", analysis=ANALYSES),
        "export/summary.json"
    shell:
        """
        python scripts/upload-to-s3.py \
            --export-dir export \
            --prefix trajectories
        """
