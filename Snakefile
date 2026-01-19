configfile: "defaults/config.yaml"

# Get all analyses from config, or use target_analyses if specified on command line
ANALYSES = config.get("target_analyses", list(config["analysis"].keys()))

# Constrain {analysis} wildcard to valid analysis names only
wildcard_constraints:
    analysis = "|".join(config["analysis"].keys())

rule all:
    input:
        expand("data/{analysis}/metadata.tsv", analysis=ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        expand("results/{analysis}", analysis=ANALYSES),

rule results:
    input:
        expand("data/{analysis}/metadata.tsv", analysis=ANALYSES),
        expand("data/{analysis}/branches.tsv", analysis=ANALYSES),
        expand("results/{analysis}", analysis=ANALYSES),

rule download_auspice_json:
    output:
        tree = "data/{analysis}/auspice.json"
    params:
        dataset = lambda wildcards: config["analysis"][wildcards.analysis]["dataset"]
    shell:
        """
        mkdir -p data/{wildcards.analysis}
        nextstrain remote download {params.dataset:q} {output.tree:q}
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
        outdir = directory("results/{analysis}")
    params:
        summary = "results/summary.json",
        url = lambda wildcards: config["analysis"][wildcards.analysis]["dataset"]
    shell:
        """
        python scripts/trajectory.py \
            --branches {input.branches:q} \
            --alignment {input.alignment:q} \
            --output-dir {output.outdir:q} \
            --summary {params.summary:q} \
            --dataset {wildcards.analysis} \
            --url {params.url:q}
        """

rule package:
    input:
        trajdir = "results/{analysis}"
    output:
        sharddir = directory("export/{analysis}")
    shell:
        """
        python scripts/package.py \
            --input-dir {input.trajdir:q} \
            --output-dir {output.sharddir:q} \
            --shuffle
        """

rule copy_summary:
    input:
        summary = "results/summary.json",
        trajdirs = expand("results/{analysis}", analysis=ANALYSES)
    output:
        "export/summary.json"
    shell:
        "cp {input.summary} {output}"

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
