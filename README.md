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

For UShER-based datasets (like `spike-lg`), matUtils is required. Install via miniconda:

```bash
# Install miniconda via Homebrew (Mac)
brew install --cask miniconda

# Create isolated usher environment (no need to run conda init)
/opt/homebrew/Caskroom/miniconda/base/bin/conda create -n usher -c conda-forge -c bioconda usher -y
```

The workflow automatically detects matUtils in common conda locations. To verify installation:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/usher/bin/matUtils --version
```

On Linux, miniconda installs to `~/miniconda3/` instead of `/opt/homebrew/Caskroom/miniconda/base/`.

# Workflow

Before executing the workflow, please run `nextstrain login https://nextstrain.org` to access `trajectories-private` datasets, or alternatively, remove those datasets from your chosen config.

`defaults/config.yaml` holds only system-wide settings (`s3_prefix`, `trajectory_mode`) and no analyses. Every invocation must explicitly select a dataset config via `--configfile`:

```
snakemake --configfile defaults/viral.yaml --cores 1 -p results   # Generate trajectory shards
snakemake --configfile defaults/viral.yaml --cores 1 -p upload    # Upload results to S3
```

Running `snakemake` with no target defaults to `results`. Running with no `--configfile` is a no-op (no analyses defined).

## Dataset config files

Pre-defined dataset bundles live in `defaults/`. Pick one with `--configfile`:

```
snakemake --configfile defaults/viral.yaml --cores 1 -p results
snakemake --configfile defaults/bac120-cyano.yaml --cores 1 -p results
snakemake --configfile defaults/trellis.yaml --cores 1 -p results
snakemake --configfile defaults/odb-fungi.yaml --cores 1 -p results
```

Multiple `--configfile` flags can be stacked; later files override earlier keys. Each dataset config can also set its own `s3_prefix` (e.g. `defaults/trellis.yaml` uploads to `s3://{bucket}/trellis-trajectories/` to keep phylogenetic outputs separate from the base `trajectories/` prefix).

## Dataset-specific outputs

To narrow to specific datasets within a config, use `target_analyses`:
```
# Provision n450-xs only from viral.yaml
snakemake --configfile defaults/viral.yaml --cores 1 -p results --config target_analyses='["n450-xs"]'

# Provision multiple datasets
snakemake --configfile defaults/viral.yaml --cores 1 -p results --config target_analyses='["n450-xs","flu-h3-xs"]'
```

## Available datasets

**Viral datasets** (`defaults/viral.yaml`):

- `n450-xs`: [Measles N450 sequences](https://nextstrain.org/groups/trajectories/n450-xs) (2429 sequences x 450 nucleotides)
- `flu-h3-xs`: [H3N2 HA1 sequences](https://nextstrain.org/groups/trajectories/flu-h3-xs) (10,263 sequences x 987 nucleotides)
- `spike-lg`: SARS-CoV-2 full spike from [Viridian](https://www.nature.com/articles/s41592-025-02947-1) global tree (~4.5M sequences x 2055 nucleotides)

**RdRp datasets** (`defaults/viral.yaml`):
- `rdrp-paramyxoviridae-xs`: Paramyxoviridae L Domain V (3,985 sequences x 1,653 nucleotides)
- `rdrp-flaviviridae-xs`: Flaviviridae NS5 RdRp (4,785 sequences x 1,884 nucleotides)
- `rdrp-picornaviridae-xs`: Picornaviridae 3D polymerase (2,627 sequences x 1,386 nucleotides)

**RdRp subtrees (opt-in):**

Subtree datasets are auto-discovered from `../rdrp/phylogenetic/auspice/*/subtrees/` but must be explicitly enabled (and depend on the rdrp parents in `defaults/viral.yaml`):
```bash
snakemake --configfile defaults/viral.yaml --cores 8 -p results --config include_subtrees=true
```
This adds datasets like `rdrp-paramyxoviridae-xs_001`, `rdrp-flaviviridae-xs_001`, etc.

**bac120 marker gene datasets:**

Bacterial marker gene datasets from the [GTDB bac120](https://gtdb.ecogenomic.org/) set are configured in separate per-phylum YAML files in `defaults/`. Select one explicitly with `--configfile defaults/bac120-{phylum}.yaml`. Currently available:

- `defaults/bac120-cyano.yaml` — 123 Cyanobacteria marker genes (~2,500 genomes per marker, 389–11,348 nucleotides)
- `defaults/bac120-bacteroidota.yaml` — 124 Bacteroidota marker genes (~22,000 genomes per marker, 658–41,943 nucleotides)

This gives 247 datasets totaling ~3M sequences. Auspice JSONs are hosted at `nextstrain.org/groups/trajectories-private/bac120/{cyano,bacteroidota}/<marker>`. New phyla can be added by dropping a `defaults/bac120-*.yaml` config file.

**Suffixing:**

Dataset names include a size suffix indicating the number of tips:
- `xs`: 1k - 10k tips
- `sm`: 10k - 100k tips
- `md`: 100k - 1M tips
- `lg`: 1M - 10M tips

## UShER datasets

Some large datasets use UShER mutation-annotated trees (protobuf format) instead of Auspice JSON. These datasets have a `usher_pb` field in the config instead of `dataset`. The workflow automatically detects and processes them differently:

1. Downloads the protobuf file
2. Uses `matUtils extract` to get tree structure and mutations
3. Fetches reference sequence from NCBI
4. Reconstructs all node sequences by applying mutations from root
5. Applies train/test split using tree structure

The `spike-lg` dataset uses the [Viridian](https://www.nature.com/articles/s41592-025-02947-1) global SARS-CoV-2 tree (~4.5M genomes), which produces cleaner sequences than the standard UShER tree due to its amplicon-aware consensus calling pipeline. By default it processes all sequences. To subsample for testing, add a `subsample` line to the config.

## Train/test split

The workflow automatically splits tips into train and test sets by marking entire clades as test data. This ensures test trajectories represent independent evolutionary lineages. See [notes/train_test.md](notes/train_test.md) for details on the algorithm and how it affects trajectory construction.

## S3 upload

The `upload` target uploads trajectory shards to S3:

```
snakemake --configfile defaults/viral.yaml --cores 1 -p upload
```

This uploads `results/` to `s3://{bucket}/{s3_prefix}/`, where `s3_prefix` defaults to `trajectories` (set in `defaults/config.yaml`) and can be overridden per dataset config (e.g. `defaults/trellis.yaml` sets `s3_prefix: trellis-trajectories`). Requires `S3_BUCKET`, `AWS_ACCESS_KEY_ID`, and `AWS_SECRET_ACCESS_KEY` environment variables to be set.

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
├── n450-xs/
│   ├── forwards-train-000.tar.zst
│   ├── forwards-test-000.tar.zst
│   ├── pairwise-train-000.tar.zst
│   ├── pairwise-train-001.tar.zst
│   ├── ...
│   └── pairwise-test-000.tar.zst
├── flu-h3-xs/
│   └── ...
└── spike-lg/
    └── ...
```

Each shard contains up to 10,000 trajectories (configurable via `shard_size` in config). Files are shuffled within each shard before writing. Larger datasets will have multiple shards (e.g., `pairwise-train-000.tar.zst`, `pairwise-train-001.tar.zst`, etc.).

To inspect shard contents:

```bash
# List files in a shard
zstd -d -c results/n450-xs/forwards-train-000.tar.zst | tar -tf -

# View first trajectory
zstd -d -c results/n450-xs/forwards-train-000.tar.zst | tar -xOf - | head -50

# Extract a specific file
zstd -d -c results/n450-xs/forwards-train-000.tar.zst | tar -xOf - SomeFile.fasta

# Extract all files to current directory
zstd -d -c results/n450-xs/forwards-train-000.tar.zst | tar -xf -
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
  "n450-xs": {
    "git_commit": "d7c62d4",
    "url": "nextstrain.org/groups/trajectories/n450-xs",
    "num_tips": 10195,
    "num_nodes": 19960,
    "alignment_length": 2055,
    "trimmed_length": { "min": 2055, "max": 2055, "mean": 2055.0 },
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
    "pairwise_trimmed_length": { "min": 2055, "max": 2055, "mean": 2055.0 },
    "pairwise_train_hamming": { "min": 0, "max": 80, "mean": 35.2 },
    "pairwise_test_hamming": { "min": 0, "max": 45, "mean": 12.3 }
  },
  "flu-h3-xs": { ... },
  "spike-lg": { ... }
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
