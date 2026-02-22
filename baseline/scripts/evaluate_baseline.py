"""
Baseline random-mutation predictor evaluation.

For each source-target pair, generates a single random prediction by flipping
N random positions in the source to a uniformly random different nucleotide
(Jukes-Cantor model), then scores the prediction using mutation accuracy.
Accuracy is aggregated across all pairs.

Two modes:
  Evaluate: read test shards, run baseline predictor, write detail TSV
  Summarize: aggregate detail TSVs into summary JSON

Usage:
  python evaluate_baseline.py \
      --shards shard1.tar.zst shard2.tar.zst \
      --trajectory-type forwards \
      --dataset spike-xs \
      --output baseline/results/spike-xs/forwards-detail.tsv

  python evaluate_baseline.py \
      --summarize \
      --results-dir baseline/results \
      --output baseline/results/summary.json
"""

import argparse
import csv
import io
import json
import os
import random
import tarfile

import zstandard as zstd
from tqdm import tqdm

VALID_BASES = set("ACGT")


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta_frames(content):
    """Parse FASTA content into list of (name, generation_token, sequence).

    Handles both header formats:
      - Forwards: >name|branch_distance|direct_distance  (3 fields)
      - Pairwise: >name|distance                         (2 fields)

    In both cases, parts[1] is the branch/pairwise distance used as the
    generation token (N for the prediction task).
    Sequences may be wrapped across multiple lines.
    """
    frames = []
    current_name = None
    current_token = None
    seq_parts = []

    for line in content.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                frames.append((current_name, current_token, "".join(seq_parts)))
            header = line[1:]
            parts = header.split("|")
            current_name = parts[0]
            current_token = int(parts[1])
            seq_parts = []
        else:
            seq_parts.append(line)

    if current_name is not None:
        frames.append((current_name, current_token, "".join(seq_parts)))

    return frames


# ---------------------------------------------------------------------------
# Prediction and scoring
# ---------------------------------------------------------------------------

def random_flip(seq, n, rng):
    """Generate a random prediction by flipping n positions to a random different base.

    Selects n random ACGT positions in the source and flips each to one of
    the 3 alternative nucleotides with equal probability (Jukes-Cantor model).
    Returns the predicted sequence as a string.
    """
    bases = list("ACGT")
    mutated = list(seq)
    valid_positions = [i for i in range(len(seq)) if seq[i] in VALID_BASES]
    if n == 0 or not valid_positions:
        return "".join(mutated)
    chosen = rng.sample(valid_positions, min(n, len(valid_positions)))
    for pos in chosen:
        original = mutated[pos]
        alternatives = [b for b in bases if b != original]
        mutated[pos] = rng.choice(alternatives)
    return "".join(mutated)


def score_prediction(source, target, predicted):
    """Score a predicted sequence against the true target.

    Returns (N, L, M, D, mutation_accuracy) where:
      N = true mutation count (source vs target, ACGT positions only)
      L = number of positions where both source and target are ACGT
      M = predicted mutation count (source vs predicted, ACGT positions only)
      D = Hamming distance between predicted and target (ACGT positions only)
      mutation_accuracy = (correct - |M - N|) / N, or None if N = 0
    """
    N = 0
    L = 0
    M = 0
    D = 0
    correct = 0

    for s, t, p in zip(source, target, predicted):
        if s not in VALID_BASES or t not in VALID_BASES:
            continue
        L += 1
        if s != t:
            N += 1
            if p == t:
                correct += 1
        if p in VALID_BASES and p != s:
            M += 1
        if p in VALID_BASES and p != t:
            D += 1

    if N == 0:
        return N, L, M, D, None
    accuracy = (correct - abs(M - N)) / N
    return N, L, M, D, accuracy


# ---------------------------------------------------------------------------
# Shard reading
# ---------------------------------------------------------------------------

def iter_fasta_from_shards(shard_paths):
    """Yield (filename, fasta_content_str) from tar.zst shards."""
    dctx = zstd.ZstdDecompressor()
    for shard_path in shard_paths:
        with open(shard_path, "rb") as fh:
            with dctx.stream_reader(fh) as reader:
                with tarfile.open(fileobj=reader, mode="r|") as tar:
                    for member in tar:
                        if not member.isfile():
                            continue
                        f = tar.extractfile(member)
                        if f is None:
                            continue
                        content = f.read().decode("utf-8")
                        yield member.name, content


# ---------------------------------------------------------------------------
# Pair extraction
# ---------------------------------------------------------------------------

def extract_forwards_pairs(frames):
    """From a forwards trajectory, yield consecutive (source, target, N) tuples.

    For consecutive frames (frame_i, frame_i+1):
      source = frame_i sequence
      N = frame_i+1 generation_token (per-edge branch distance from frame_i)
      target = frame_i+1 sequence
    """
    for i in range(len(frames) - 1):
        source_name, source_token, source_seq = frames[i]
        target_name, target_token, target_seq = frames[i + 1]
        n = target_token
        yield source_name, target_name, source_seq, target_seq, n


def extract_pairwise_pair(frames):
    """From a pairwise trajectory, yield single (source, target, N) tuple.

    First tip has |0, second tip has |N.
    """
    if len(frames) != 2:
        return None
    source_name, source_token, source_seq = frames[0]
    target_name, target_token, target_seq = frames[1]
    return source_name, target_name, source_seq, target_seq, target_token


# ---------------------------------------------------------------------------
# Evaluate mode
# ---------------------------------------------------------------------------

def run_evaluate(args):
    """Process shards: for each pair, generate one random prediction and score it."""
    rng = random.Random(args.seed)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Count files for progress bar
    total_files = 0
    dctx = zstd.ZstdDecompressor()
    for shard_path in args.shards:
        with open(shard_path, "rb") as fh:
            with dctx.stream_reader(fh) as reader:
                with tarfile.open(fileobj=reader, mode="r|") as tar:
                    for member in tar:
                        if member.isfile():
                            total_files += 1

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "dataset", "trajectory", "source_node", "target_node",
            "N", "L", "M", "D", "mutation_accuracy",
        ])

        for filename, content in tqdm(
            iter_fasta_from_shards(args.shards),
            total=total_files,
            desc=f"{args.dataset} {args.trajectory_type}",
        ):
            frames = parse_fasta_frames(content)
            trajectory_name = os.path.splitext(os.path.basename(filename))[0]

            if args.trajectory_type == "forwards":
                pairs = list(extract_forwards_pairs(frames))
            else:
                pair = extract_pairwise_pair(frames)
                pairs = [pair] if pair else []

            for source_name, target_name, source_seq, target_seq, n in pairs:
                predicted_seq = random_flip(source_seq, n, rng)
                N, L, M, D, accuracy = score_prediction(
                    source_seq, target_seq, predicted_seq,
                )

                writer.writerow([
                    args.dataset, trajectory_name, source_name, target_name,
                    N, L, M, D,
                    f"{accuracy:.6f}" if accuracy is not None else "NA",
                ])


# ---------------------------------------------------------------------------
# Summarize mode
# ---------------------------------------------------------------------------

def run_summarize(args):
    """Read all detail TSVs and write summary JSON."""
    summary = {}

    results_dir = args.results_dir
    for dataset_dir in sorted(os.listdir(results_dir)):
        dataset_path = os.path.join(results_dir, dataset_dir)
        if not os.path.isdir(dataset_path):
            continue

        dataset_summary = {}

        for traj_type in ["forwards", "pairwise"]:
            detail_path = os.path.join(dataset_path, f"{traj_type}-detail.tsv")
            if not os.path.exists(detail_path):
                continue

            rows = []
            with open(detail_path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    rows.append(row)

            if not rows:
                continue

            n_predictions = len(rows)
            n_values = []
            d_values = []
            acc_values = []
            by_n = {}

            for row in rows:
                n = int(row["N"])
                n_values.append(n)
                d = int(row["D"])
                d_values.append(d)

                acc_str = row["mutation_accuracy"]
                if acc_str != "NA":
                    acc = float(acc_str)
                    acc_values.append(acc)
                else:
                    acc = None

                n_str = str(n)
                if n_str not in by_n:
                    by_n[n_str] = {"count": 0, "d_values": [], "acc_values": []}
                by_n[n_str]["count"] += 1
                by_n[n_str]["d_values"].append(d)
                if acc is not None:
                    by_n[n_str]["acc_values"].append(acc)

            # Compute aggregates
            mean_n = sum(n_values) / len(n_values) if n_values else 0
            mean_d = sum(d_values) / len(d_values) if d_values else 0
            mean_acc = sum(acc_values) / len(acc_values) if acc_values else 0

            # Build by_N summary
            by_n_summary = {}
            for n_str in sorted(by_n.keys(), key=int):
                entry = by_n[n_str]
                mean_d_n = sum(entry["d_values"]) / len(entry["d_values"])
                mean_acc_n = (
                    sum(entry["acc_values"]) / len(entry["acc_values"])
                    if entry["acc_values"]
                    else None
                )
                by_n_entry = {
                    "count": entry["count"],
                    "mean_D": round(mean_d_n, 4),
                }
                if mean_acc_n is not None:
                    by_n_entry["mean_accuracy"] = round(mean_acc_n, 6)
                else:
                    by_n_entry["mean_accuracy"] = None
                by_n_summary[n_str] = by_n_entry

            dataset_summary[traj_type] = {
                "n_predictions": n_predictions,
                "mean_N": round(mean_n, 2),
                "mean_mutation_accuracy": round(mean_acc, 6),
                "mean_D": round(mean_d, 4),
                "by_N": by_n_summary,
            }

        if dataset_summary:
            summary[dataset_dir] = dataset_summary

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote summary to {args.output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Baseline random-mutation predictor evaluation")

    # Mode selection
    parser.add_argument("--summarize", action="store_true",
                        help="Run in summarize mode (aggregate detail TSVs)")

    # Evaluate mode args
    parser.add_argument("--shards", nargs="+",
                        help="Path(s) to test shard tar.zst files")
    parser.add_argument("--trajectory-type", choices=["forwards", "pairwise"],
                        help="Type of trajectory to evaluate")
    parser.add_argument("--dataset", help="Dataset name for TSV output")
    parser.add_argument("--output", required=True,
                        help="Output path (detail TSV or summary JSON)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")

    # Summarize mode args
    parser.add_argument("--results-dir",
                        help="Path to results directory containing dataset subdirs")

    args = parser.parse_args()

    if args.summarize:
        if not args.results_dir:
            parser.error("--results-dir required in summarize mode")
        run_summarize(args)
    else:
        if not args.shards:
            parser.error("--shards required in evaluate mode")
        if not args.trajectory_type:
            parser.error("--trajectory-type required in evaluate mode")
        if not args.dataset:
            parser.error("--dataset required in evaluate mode")
        run_evaluate(args)


if __name__ == "__main__":
    main()
