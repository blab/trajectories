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
snakemake --cores 1 -p upload    # After export, upload to S3
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

- `cytb-xs`: [Mammalian cytochrome B sequences](https://next.nextstrain.org/groups/trajectories/cytb-xs) (5059 sequences x 1140 nucleotides)
- `n450-xs`: [Measles N450 sequences](https://next.nextstrain.org/groups/trajectories/n450-xs) (2429 sequences x 450 nucleotides)
- `spike-xs`: [SARS-CoV-2 spike S1 sequences](https://next.nextstrain.org/groups/trajectories/spike-xs) (10,195 sequences x 2055 nucleotides)
- `spike-sm`: [SARS-CoV-2 spike S1 sequences](https://next.nextstrain.org/groups/trajectories/spike-sm) (34,707 sequences x 2055 nucleotides)

**RdRp datasets (local):**
- `rdrp-paramyxoviridae`: Paramyxoviridae L Domain V (3,985 sequences x 1,653 nucleotides)
- `rdrp-flaviviridae`: Flaviviridae NS5 RdRp (4,785 sequences x 1,884 nucleotides)
- `rdrp-picornaviridae`: Picornaviridae 3D polymerase (2,627 sequences x 1,386 nucleotides)
- Plus subtree datasets (e.g., `rdrp-paramyxoviridae_001`, `rdrp-flaviviridae_001`)

Dataset names include a size suffix indicating the number of tips:
- `xs`: 1k - 10k tips
- `sm`: 10k - 100k tips
- `md`: 100k - 1M tips
- `lg`: 1M - 10M tips

## Train/test split

The workflow automatically splits tips into train and test sets by marking entire clades as test data. This ensures test trajectories represent independent evolutionary lineages. See [notes/train_test.md](notes/train_test.md) for details on the algorithm and how it affects trajectory construction.

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
- `branches.tsv` - Parent-child relationships with Hamming distances and train/test labels

## Trajectory files

The main output is trajectory files in `results/{dataset}/`, organized into forwards (root-to-tip) and pairwise (tip-to-tip) subdirectories:

```
results/
├── spike-xs/
│   ├── forwards-train/
│   │   ├── USACA-CDC-STM-A1234562021.fasta
│   │   └── ...
│   ├── forwards-test/
│   │   ├── Algeria10372023.fasta
│   │   └── ...
│   ├── pairwise-train/
│   │   ├── TipA__TipB.fasta
│   │   └── ...
│   └── pairwise-test/
│       ├── TipX__TipY.fasta
│       └── ...
├── cytb-xs/
│   └── ...
└── n450-xs/
    └── ...
```

### Forwards trajectories

Each forwards trajectory is a FASTA file containing the evolutionary path from root to tip:

```
>NODE_0000000|0
ATGTTCGTTTTT...
>NODE_0001234|15
ATGTTCGTTTTT...
>TipName|42
ATGTTCGTTTTT...
```

Where each header contains `>{node_name}|{cumulative_hamming_distance}`. Intermediate nodes with zero mutations are skipped (root and tip are always included).

**Training trajectories** contain the full root-to-tip path. **Test trajectories** are truncated to start at the test clade boundary, ensuring they contain only evolutionary history unseen during training.

### Pairwise trajectories

Each pairwise trajectory is a FASTA file containing two tip sequences with their Hamming distance:

```
>TipA|0
ATGTTCGTTTTT...
>TipB|23
ATGTTCGTTTAT...
```

The first sequence gets `|0` and the second gets `|{hamming_distance}` representing the distance between the two sequences. File naming uses double underscore separator: `{tip1}__{tip2}.fasta`.

**Training pairs** are random samples from all training tips (default limit: 100K pairs). **Test pairs** are only generated within the same test clade to avoid overlap with training branches (default limit: 50K pairs). Limits can be configured via `pairwise_train_limit` and `pairwise_test_limit` in config.

## Export shards

For distribution and efficient data loading, trajectories are packaged into sharded tar.zst archives in `export/`:

```
export/
├── spike-xs/
│   ├── forwards-train-000.tar.zst
│   ├── forwards-test-000.tar.zst
│   ├── pairwise-train-000.tar.zst
│   ├── pairwise-train-001.tar.zst
│   ├── ...
│   └── pairwise-test-000.tar.zst
├── cytb-xs/
│   └── ...
├── n450-xs/
│   └── ...
└── summary.json
```

Each shard contains up to 10,000 trajectories. Files are shuffled before sharding. Larger datasets will have multiple shards (e.g., `pairwise-train-000.tar.zst`, `pairwise-train-001.tar.zst`, etc.).

## Summary statistics

A consolidated `results/summary.json` file contains statistics for all processed datasets:

```json
{
  "spike-xs": {
    "url": "next.nextstrain.org/groups/trajectories/spike-xs",
    "num_tips": 10195,
    "num_nodes": 19960,
    "sequence_length": 2055,
    "hamming_from_root": { "min": 0, "max": 80, "mean": 27.18 },
    "path_depth": { "min": 1, "max": 24, "mean": 10.57 },
    "total_branches": 19959,
    "zero_distance_branches": 15095,
    "per_branch_hamming": { "min": 0, "max": 37, "mean": 0.35 },
    "train_tips": 9172,
    "test_tips": 1023,
    "pairwise_train_pairs": 100000,
    "pairwise_test_pairs": 8500,
    "pairwise_test_clades": 25,
    "pairwise_train_hamming": { "min": 0, "max": 80, "mean": 35.2 },
    "pairwise_test_hamming": { "min": 0, "max": 45, "mean": 12.3 }
  },
  "cytb-xs": { ... },
  "spike-sm": { ... },
  "n450-xs": { ... }
}
```

Each dataset entry is added or updated when its trajectories are generated. The `train_tips` and `test_tips` fields indicate the number of forwards trajectories in each split. The `pairwise_*` fields show pairwise pair counts, number of test clades, and Hamming distance statistics.

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
