# Trajectories

This is a simple repo to provision evolutionary sequence trajectories from Nextstrain trees

# Installation

Install various dependencies with pip:

```
pip install -r requirements.txt
```

Nextstrain CLI should be installed [following its docs](https://docs.nextstrain.org/en/latest/install.html#install-nextstrain-cli) and looks like

```
# Mac
curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/mac | bash
# Linux
curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/linux | bash
```

# Workflow

The workflow is managed through Snakemake with three main targets:

```
snakemake --cores 1 -p results   # Generate trajectory FASTA files
snakemake --cores 1 -p export    # Package into tar.zst shards
snakemake --cores 1 -p upload    # Upload to S3
```

Running `snakemake` with no target defaults to `results`.

## Dataset-specific outputs

To provision specific datasets, use the `target_analyses` config:
```
# Provision spike-xs dataset only
snakemake --cores 1 -p results --config target_analyses='["spike-xs"]'

# Provision multiple datasets
snakemake --cores 1 -p results --config target_analyses='["spike-xs","cytb-xs"]'
```

## Available datasets

The following datasets are pre-configured and can be used with the above commands:

- `spike-xs`: [SARS-CoV-2 spike S1 sequences](https://nextstrain.org/groups/blab/ncov/trajectories/10k/2025-09-15) (10,195 sequences x 2055 nucleotides)
- `spike-sm`: [SARS-CoV-2 spike S1 sequences](https://nextstrain.org/groups/blab/ncov/trajectories/30k/2026-01-01) (34,707 sequences x 2055 nucleotides)
- `cytb-xs`: [Mammalian cytochrome B sequences](https://nextstrain.org/groups/blab/cytb) (5059 sequences x 1140 nucleotides)
- `n450-xs`: [Measles N450 sequences](https://nextstrain.org/measles/N450@2025-10-01) (2429 sequences x 450 nucleotides)

Dataset names include a size suffix indicating the number of tips:
- `xs`: 1k - 10k tips
- `sm`: 10k - 100k tips
- `md`: 100k - 1M tips
- `lg`: 1M - 10M tips

## S3 upload

The `upload` target packages trajectories and uploads to S3:

```
snakemake --cores 1 -p upload
```

This uploads `export/` to `s3://{bucket}/trajectories/`. Requires `S3_BUCKET`, `AWS_ACCESS_KEY_ID`, and `AWS_SECRET_ACCESS_KEY` environment variables to be set.

# Outputs

## Intermediate data files

For each dataset, the workflow generates intermediate files in `data/{dataset}/`:

- `auspice.json` - Original Nextstrain tree data
- `alignment.fasta` - Sequences for all nodes (tips and internal)
- `metadata.tsv` - Phylogenetic metadata with parent relationships
- `branches.tsv` - Parent-child relationships with Hamming distances

## Trajectory files

The main output is one trajectory file per tip in `results/{dataset}/`:

```
results/
├── spike-xs/
│   ├── USACA-CDC-STM-A1234562021.fasta
│   ├── NigeriaISTH-E02312020.fasta
│   └── ...
├── cytb-xs/
│   └── ...
└── n450-xs/
    └── ...
```

Each trajectory is a FASTA file containing the evolutionary path from root to tip:

```
>NODE_0000000|0
ATGTTCGTTTTT...
>NODE_0001234|15
ATGTTCGTTTTT...
>TipName|42
ATGTTCGTTTTT...
```

Where each header contains `>{node_name}|{cumulative_hamming_distance}`. Intermediate nodes with zero mutations are skipped (root and tip are always included).

## Export shards

For distribution and efficient data loading, trajectories are packaged into sharded tar.zst archives in `export/`:

```
export/
├── spike-xs/
│   └── trajectories-000.tar.zst
├── cytb-xs/
│   └── trajectories-000.tar.zst
├── n450-xs/
│   └── trajectories-000.tar.zst
└── summary.json
```

Each shard contains up to 10,000 trajectories. Files are shuffled before sharding. Larger datasets will have multiple shards (e.g., `trajectories-000.tar.zst`, `trajectories-001.tar.zst`, etc.).

## Summary statistics

A consolidated `results/summary.json` file contains statistics for all processed datasets:

```json
{
  "cytb-xs": {
    "url": "nextstrain.org/groups/blab/cytb",
    "num_tips": 5059,
    "num_nodes": 10117,
    "sequence_length": 1140,
    "hamming_from_root": { "min": 154, "max": 602, "mean": 377.81 },
    "path_depth": { "min": 6, "max": 43, "mean": 22.28 },
    "total_branches": 10116,
    "zero_distance_branches": 373,
    "per_branch_hamming": { "min": 0, "max": 253, "mean": 29.24 }
  },
  "spike-xs": { ... },
  "spike-sm": { ... },  
  "n450-xs": { ... }
}
```

Each dataset entry is added or updated when its trajectories are generated.

# License

This repository is licensed under the MIT License. See the LICENSE file for
details.

## Important Disclaimer About Copyright and AI-Generated Code

Some portions of the code in this repository were generated with the assistance
of large language models (LLMs), primarily Claude Code. Individual scripts are
commented to state their provenance. While I have reviewed, modified, and
integrated these contributions, the copyright status of LLM-generated code is
uncertain and may vary depending on jurisdiction.

As a result:

- Human-Authored Contributions: Code written by me (the repository owner) is
  explicitly licensed under the MIT License and is subject to the terms outlined
  in the LICENSE file.
- LLM-Generated Contributions: For any portions of the code generated by LLMs, I
  do not assert copyright ownership and disclaim any responsibility for the
  originality or copyright status of such code.
- User Responsibility: Users of this repository are encouraged to independently
  verify the legal status of any LLM-generated portions of the code before reuse
  or redistribution.
