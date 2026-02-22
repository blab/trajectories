"""
Baseline random-mutation predictor evaluation.

Two modes:
  Evaluate: read test shards, run Monte Carlo random-flip predictions, write detail TSV
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
N_TRIALS = 100


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta_frames(content):
    """Parse FASTA content into list of (name, generation_token, sequence).

    Header format: >name|generation_token
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
            parts = header.rsplit("|", 1)
            current_name = parts[0]
            current_token = int(parts[1])
            seq_parts = []
        else:
            seq_parts.append(line)

    if current_name is not None:
        frames.append((current_name, current_token, "".join(seq_parts)))

    return frames


# ---------------------------------------------------------------------------
# Mutation helpers
# ---------------------------------------------------------------------------

def find_mutation_positions(source, target):
    """Return list of indices where source differs from target (both ACGT)."""
    positions = []
    for i in range(len(source)):
        if source[i] in VALID_BASES and target[i] in VALID_BASES and source[i] != target[i]:
            positions.append(i)
    return positions


def count_valid_positions(source, target):
    """Count positions where both source and target are ACGT."""
    count = 0
    for i in range(len(source)):
        if source[i] in VALID_BASES and target[i] in VALID_BASES:
            count += 1
    return count


def random_flip(seq, n, valid_positions, rng):
    """Flip n random positions to a uniformly random different nucleotide.

    Returns mutated sequence as a list of characters.
    """
    bases = list("ACGT")
    mutated = list(seq)
    if n == 0 or not valid_positions:
        return mutated
    chosen = rng.sample(valid_positions, min(n, len(valid_positions)))
    for pos in chosen:
        original = mutated[pos]
        alternatives = [b for b in bases if b != original]
        mutated[pos] = rng.choice(alternatives)
    return mutated


def hamming_distance(seq_a, seq_b):
    """Hamming distance counting only positions where both are ACGT."""
    d = 0
    for i in range(len(seq_a)):
        a, b = seq_a[i], seq_b[i]
        if a in VALID_BASES and b in VALID_BASES and a != b:
            d += 1
    return d


# ---------------------------------------------------------------------------
# Monte Carlo evaluation
# ---------------------------------------------------------------------------

def evaluate_prediction(source, target, n_mutations, n_trials, rng):
    """Run n_trials random-flip predictions and return (mean_D, mean_M, mean_accuracy).

    mean_D: mean Hamming distance between prediction and target
    mean_M: mean number of predicted mutations (positions where prediction != source)
    mean_accuracy: mean mutation_accuracy = (correct - |M - N|) / N
    """
    source_list = list(source)
    target_list = list(target)

    # Valid positions in source for flipping
    valid_positions = [i for i in range(len(source)) if source[i] in VALID_BASES]

    # Real mutation positions (source != target, both ACGT)
    mutation_positions = set(find_mutation_positions(source, target))
    n_real_mutations = len(mutation_positions)

    total_d = 0
    total_m = 0
    total_accuracy = 0

    for _ in range(n_trials):
        predicted = random_flip(source_list, n_mutations, valid_positions, rng)

        # Hamming distance: prediction vs target (only valid positions)
        d = 0
        correct_mutations = 0
        m_predicted = 0
        for i in range(len(predicted)):
            p, t, s = predicted[i], target_list[i], source_list[i]
            if p in VALID_BASES and t in VALID_BASES:
                if p != t:
                    d += 1
            # Count predicted mutations (positions where prediction != source)
            if p in VALID_BASES and s in VALID_BASES and p != s:
                m_predicted += 1
            # Check if this mutation was correctly predicted
            if i in mutation_positions and p == t:
                correct_mutations += 1

        total_d += d
        total_m += m_predicted
        if n_real_mutations > 0:
            total_accuracy += (correct_mutations - abs(m_predicted - n_real_mutations)) / n_real_mutations

    mean_d = total_d / n_trials
    mean_m = total_m / n_trials
    mean_accuracy = total_accuracy / n_trials if n_real_mutations > 0 else float("nan")
    return mean_d, mean_m, mean_accuracy


# ---------------------------------------------------------------------------
# Analytical expectations
# ---------------------------------------------------------------------------

def expected_random_hamming(N, L):
    """E[D] = 2N - (4/3)(N^2/L) for random Jukes-Cantor flipping."""
    if L == 0:
        return 0.0
    return 2 * N - (4 / 3) * (N * N / L)


def expected_random_accuracy(N, L):
    """E[accuracy] = N/(3L) for random Jukes-Cantor flipping."""
    if L == 0 or N == 0:
        return 0.0
    return N / (3 * L)


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
    """Process shards and write detail TSV."""
    rng = random.Random(42)

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
            "N", "L", "M_predicted", "D_predicted", "mutation_accuracy",
            "D_analytical", "accuracy_analytical",
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
                L = count_valid_positions(source_seq, target_seq)

                if n <= 0:
                    # No mutations to predict — D=0, M=0, accuracy is N/A
                    writer.writerow([
                        args.dataset, trajectory_name, source_name, target_name,
                        0, L, 0.0, 0.0, "NA",
                        0.0, "NA",
                    ])
                    continue

                mean_d, mean_m, mean_acc = evaluate_prediction(
                    source_seq, target_seq, n, N_TRIALS, rng,
                )

                d_analytical = expected_random_hamming(n, L)
                acc_analytical = expected_random_accuracy(n, L)

                writer.writerow([
                    args.dataset, trajectory_name, source_name, target_name,
                    n, L,
                    f"{mean_m:.4f}", f"{mean_d:.4f}", f"{mean_acc:.6f}",
                    f"{d_analytical:.4f}", f"{acc_analytical:.6f}",
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
                d = float(row["D_predicted"])
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
