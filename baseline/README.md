# Baseline evaluation

Random-mutation baseline predictor for sequence trajectory prediction. For each source-target pair, flips N random positions in the source sequence (Jukes-Cantor model) and scores predictions using the mutation accuracy metric defined in `notes/mutation_accuracy_metric.md`.

## Usage

```bash
cd baseline
snakemake --cores 1 -p
```

By default this runs on `spike-xs` and `cytb-xs`. To target other datasets:

```bash
snakemake --cores 1 -p --config target_analyses='["spike-sm"]'
```
