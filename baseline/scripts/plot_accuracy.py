"""
Plot mutation accuracy histograms from baseline evaluation detail TSVs.

Usage:
    python scripts/plot_accuracy.py \
        --results-dir results \
        --output figures/accuracy_histogram.png
"""

import argparse
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def discover_detail_tsvs(results_dir):
    """Find all *-detail.tsv files and return sorted (dataset, traj_type, path) tuples."""
    pattern = os.path.join(results_dir, "*", "*-detail.tsv")
    entries = []
    for path in sorted(glob.glob(pattern)):
        basename = os.path.basename(path)
        traj_type = basename.replace("-detail.tsv", "")
        dataset = os.path.basename(os.path.dirname(path))
        entries.append((dataset, traj_type, path))
    return entries


def read_accuracy_values(path):
    """Read mutation_accuracy column from a detail TSV, skipping NA rows."""
    values = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc_str = row["mutation_accuracy"]
            if acc_str != "NA":
                values.append(float(acc_str))
    return values


def main():
    parser = argparse.ArgumentParser(description="Plot mutation accuracy histograms")
    parser.add_argument("--results-dir", required=True,
                        help="Path to results directory containing dataset subdirs")
    parser.add_argument("--output", required=True,
                        help="Output PNG path")
    args = parser.parse_args()

    entries = discover_detail_tsvs(args.results_dir)
    if not entries:
        print(f"No *-detail.tsv files found in {args.results_dir}")
        return

    n_panels = len(entries)
    fig, axes = plt.subplots(n_panels, 1, figsize=(8, 2.5 * n_panels),
                             squeeze=False)

    for i, (dataset, traj_type, path) in enumerate(entries):
        ax = axes[i, 0]
        values = read_accuracy_values(path)

        if values:
            mean_val = sum(values) / len(values)
            ax.hist(values, bins=50, edgecolor="white", linewidth=0.3)
            ax.axvline(mean_val, color="red", linestyle="--", linewidth=1)
            ax.text(0.98, 0.92, f"mean = {mean_val:.4f}",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=9, color="red")
        else:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")

        ax.set_title(f"{dataset} {traj_type}", fontsize=11, loc="left")
        ax.set_xlabel("Mutation accuracy")
        ax.set_ylabel("Count")

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
