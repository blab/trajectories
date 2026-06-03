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

`defaults/config.yaml` defines a handful of example datasets out of the box. The default invocation runs them all:

```
snakemake --cores 1 -p results
```

Running `snakemake` with no target defaults to `results`.

## Dataset-specific outputs

To narrow to specific datasets, use `target_analyses`:

```
# Provision flu-h3-xs only
snakemake --cores 1 -p results --config target_analyses='["flu-h3-xs"]'

# Provision multiple datasets
snakemake --cores 1 -p results --config target_analyses='["flu-h3-xs","cytb-xs"]'
```

## Available datasets

The four example datasets in `defaults/config.yaml`:

- `flu-h3-xs`: [H3N2 HA1 sequences](https://nextstrain.org/groups/trajectories/flu-h3-xs) (10,263 sequences x 987 nucleotides)
- `cytb-xs`: [Vertebrate cytochrome b](https://nextstrain.org/groups/trajectories/cytb-xs) (1,140 nucleotides)
- `spike-sm`: [SARS-CoV-2 spike S1](https://nextstrain.org/groups/trajectories/spike-sm) (2,055 nucleotides)
- `trellis-18aa-KEVT`: [Trellis KEVT 18-amino-acid ligand phylogeny](https://nextstrain.org/groups/trajectories/trellis/18aa/KEVT) (54 nucleotides)

Dataset names include a size suffix indicating the number of tips:
- `xs`: 1k - 10k tips
- `sm`: 10k - 100k tips
- `md`: 100k - 1M tips
- `lg`: 1M - 10M tips

To add your own dataset, append an entry under `analysis:` in `defaults/config.yaml` with at minimum a `dataset` URL (or local path), `gene`, and `seq_length`. Optional `trim_begin` / `trim_end` restrict the alignment to a window.

## Train/test split

The workflow automatically splits tips into train and test sets by marking entire clades as test data. This ensures test trajectories represent independent evolutionary lineages. See [notes/train_test.md](notes/train_test.md) for details on the algorithm and how it affects trajectory construction.

# Outputs

## Intermediate data files

For each dataset, the workflow generates intermediate files in `data/{dataset}/`:

- `auspice.json` - Original Nextstrain tree data
- `alignment.fasta` - Sequences for all nodes (tips and internal)
- `metadata.tsv` - Phylogenetic metadata with parent relationships
- `branches.tsv` - Parent-child relationships with Hamming distances and train/test labels

## Trajectory shards

The main output is sharded tar.zst archives in `results/{dataset}/`:

```
results/
├── flu-h3-xs/
│   ├── forwards-train-000.tar.zst
│   ├── forwards-test-000.tar.zst
│   ├── pairwise-train-000.tar.zst
│   ├── pairwise-train-001.tar.zst
│   ├── ...
│   └── pairwise-test-000.tar.zst
├── cytb-xs/
│   └── ...
└── spike-sm/
    └── ...
```

Each shard contains up to 10,000 trajectories (configurable via `shard_size` in config). Files are shuffled within each shard before writing. Larger datasets will have multiple shards (e.g., `pairwise-train-000.tar.zst`, `pairwise-train-001.tar.zst`, etc.).

To inspect shard contents:

```bash
# List files in a shard
zstd -d -c results/flu-h3-xs/forwards-train-000.tar.zst | tar -tf -

# View first trajectory
zstd -d -c results/flu-h3-xs/forwards-train-000.tar.zst | tar -xOf - | head -50

# Extract a specific file
zstd -d -c results/flu-h3-xs/forwards-train-000.tar.zst | tar -xOf - SomeFile.fasta

# Extract all files to current directory
zstd -d -c results/flu-h3-xs/forwards-train-000.tar.zst | tar -xf -
```

See [notes/data_format.md](notes/data_format.md) for a detailed worked example with a small tree illustrating both trajectory formats.

### Forwards trajectories

Each forwards trajectory is a FASTA file containing the evolutionary path from root to tip:

```
>NODE_0000000|0|0
ATGTTCGTTTTT...
>NODE_0001234|15|15
ATGTTCGTTTTT...
>TipName|14|27
ATGTTCGTTTTT...
```

Where each header contains `>{node_name}|{branch_hamming_distance}|{direct_hamming_distance}` — the branch distance from the previous emitted node (0 for root) and the direct Hamming distance from the start node. All Hamming distances ignore positions where either sequence has a gap (`-`) or ambiguous base (`N`). Intermediate nodes with zero mutations are skipped. If the tip has zero branch distance from the last emitted node, the last emitted node is relabeled with the tip's name rather than adding a zero-distance frame.

**Training trajectories** contain the full root-to-tip path. **Test trajectories** are truncated to start at the test clade boundary, ensuring they contain only evolutionary history unseen during training.

### Pairwise trajectories

Each pairwise trajectory is a FASTA file containing two tip sequences with their Hamming distance:

```
>TipA|0|0
ATGTTCGTTTTT...
>TipB|23|23
ATGTTCGTTTAT...
```

Headers use the same three-field format as forwards trajectories: `>{name}|{branch_distance}|{direct_distance}`. The first sequence gets `|0|0` and the second gets `|{hamming}|{hamming}` (branch and direct are always identical for pairwise). File naming uses double underscore separator: `{tip1}__{tip2}.fasta`.

**Training pairs** are random samples from all training tips (default limit: 100K pairs). **Test pairs** are only generated within the same test clade to avoid overlap with training branches (default limit: 50K pairs). Limits can be configured via `pairwise_train_limit` and `pairwise_test_limit` in config.

## Summary statistics

A consolidated `results/summary.json` file contains statistics for all processed datasets:

```json
{
  "flu-h3-xs": {
    "git_commit": "d7c62d4",
    "url": "nextstrain.org/groups/trajectories/flu-h3-xs",
    "num_tips": 10195,
    "num_nodes": 19960,
    "alignment_length": 987,
    "trimmed_length": { "min": 987, "max": 987, "mean": 987.0 },
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
    "pairwise_trimmed_length": { "min": 987, "max": 987, "mean": 987.0 },
    "pairwise_train_hamming": { "min": 0, "max": 80, "mean": 35.2 },
    "pairwise_test_hamming": { "min": 0, "max": 45, "mean": 12.3 }
  },
  "cytb-xs": { ... },
  "spike-sm": { ... }
}
```

Each dataset entry is added or updated when its trajectories are generated. `alignment_length` is the full alignment width; `trimmed_length` shows the per-trajectory length after dropping columns that are all-gap on each path (identical to `alignment_length` for viral datasets with no insertions, shorter for diverse-phylum alignments). The `train_tips` and `test_tips` fields indicate the number of forwards trajectories in each split. The `pairwise_*` fields show pairwise pair counts, number of test clades, trimmed lengths, and Hamming distance statistics.

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
