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

## Results

![Mutation accuracy histogram for spike-xs](figures/accuracy_histogram.png)

Distribution of per-prediction mutation accuracy for the spike-xs dataset. The random baseline scores near zero on both forwards (mean = 0.0002) and pairwise (mean = 0.0005) trajectories, as expected.
