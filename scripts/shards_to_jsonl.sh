#!/usr/bin/env bash
# Convert all trajectory shards for one (analysis, mode, split) into a single
# concatenated JSONL in json_sdiff format.
#
# Inputs:
#   $ANALYSIS_DIR   directory containing {mode}-{split}-NNN.tar.zst shards
#                   (e.g. results/n450-xs)
#   $MODE           "forwards" or "pairwise"
#   $SPLIT          "train" or "test"
#   $OUT_JSONL      path to write the concatenated per-(mode,split) JSONL
#   $SHARDS_DIR     directory for per-shard intermediates (resumable cache)
#   $PY             python interpreter to invoke fasta_to_jsonl.py with
#   $FASTA_TO_JSONL path to fasta_to_jsonl.py
#   $PARALLEL       worker count for parallel shard processing
#   $MAX_RAW_LEN    max assistant-content char length (passed through)
#   $MAX_LEN        optional max trajectory length (pass empty string to skip)
#   $CONTEXT_SIZE   left-context bp per hunk
#   $NOPROMPT       "1" to pass --noprompt
#
# Resumability:
#   - Per-shard outputs live under $SHARDS_DIR and are NOT deleted between runs.
#   - Each shard writes to *.jsonl.partial then atomically renames, so an
#     interrupted run never leaves a half-written file marked complete.
#   - Re-running the script skips any shard whose final .jsonl already exists.
#   - To force a redo, delete $SHARDS_DIR (or the specific shard files).

set -uo pipefail   # no -e: still concat surviving shards if some fail

: "${ANALYSIS_DIR:?ANALYSIS_DIR required}"
: "${MODE:?MODE required}"
: "${SPLIT:?SPLIT required}"
: "${OUT_JSONL:?OUT_JSONL required}"
: "${SHARDS_DIR:?SHARDS_DIR required}"
: "${PY:?PY required}"
: "${FASTA_TO_JSONL:?FASTA_TO_JSONL required}"

PARALLEL="${PARALLEL:-16}"
MAX_RAW_LEN="${MAX_RAW_LEN:-14000}"
MAX_LEN="${MAX_LEN:-}"
CONTEXT_SIZE="${CONTEXT_SIZE:-16}"
NOPROMPT="${NOPROMPT:-1}"

mkdir -p "$(dirname "$OUT_JSONL")" "$SHARDS_DIR"
export SHARDS_DIR PY FASTA_TO_JSONL MAX_RAW_LEN MAX_LEN CONTEXT_SIZE NOPROMPT

shard_pattern="${MODE}-${SPLIT}-*.tar.zst"
n_total=$(find "$ANALYSIS_DIR" -maxdepth 1 -name "$shard_pattern" | wc -l)
n_done=$(find "$SHARDS_DIR" -maxdepth 1 -name "${MODE}-${SPLIT}-*.jsonl" -not -name '*.partial' 2>/dev/null | wc -l)

echo "Analysis:  $ANALYSIS_DIR"
echo "Mode:      $MODE   Split: $SPLIT"
echo "Output:    $OUT_JSONL"
echo "Shards:    $SHARDS_DIR  (persists across runs)"
echo "Parallel:  $PARALLEL"
echo "Progress:  $n_done / $n_total shards already complete"

if [ "$n_total" -eq 0 ]; then
    echo "No ${MODE}-${SPLIT} shards in $ANALYSIS_DIR; writing empty $OUT_JSONL"
    : > "$OUT_JSONL"
    exit 0
fi

# Clean up any stale .partial files from a previous interrupted run.
find "$SHARDS_DIR" -maxdepth 1 -name "${MODE}-${SPLIT}-*.partial" -delete 2>/dev/null || true

find "$ANALYSIS_DIR" -maxdepth 1 -name "$shard_pattern" -print0 \
  | xargs -0 -I{} -P "$PARALLEL" bash -c '
      shard="$1"
      base=$(basename "$shard" .tar.zst)
      outfile="$SHARDS_DIR/${base}.jsonl"

      if [ -s "$outfile" ]; then
          exit 0
      fi

      workdir="$SHARDS_DIR/${base}.work"
      rm -rf "$workdir"
      mkdir -p "$workdir"

      if ! tar --use-compress-program=unzstd -xf "$shard" -C "$workdir"; then
          echo "FAIL extract: $base" >&2
          rm -rf "$workdir"
          exit 1
      fi

      extra=()
      [ "$NOPROMPT" = "1" ] && extra+=(--noprompt)
      [ -n "$MAX_LEN" ] && extra+=(--max-len "$MAX_LEN")

      if ! "$PY" "$FASTA_TO_JSONL" \
              --input "$workdir" \
              --output "${outfile}.partial" \
              --max-raw-len "$MAX_RAW_LEN" \
              --context-size "$CONTEXT_SIZE" \
              "${extra[@]}"; then
          echo "FAIL convert: $base" >&2
          rm -rf "$workdir" "${outfile}.partial"
          exit 1
      fi

      mv "${outfile}.partial" "$outfile"
      rm -rf "$workdir"
    ' _ {}
xargs_status=$?

if [ "$xargs_status" -ne 0 ]; then
    echo "WARNING: at least one shard failed (xargs status $xargs_status). Re-run to retry." >&2
fi

# Concatenate surviving per-shard outputs in stable (sorted) order.
shopt -s nullglob
shards=( "$SHARDS_DIR"/${MODE}-${SPLIT}-*.jsonl )
if [ "${#shards[@]}" -eq 0 ]; then
    echo "No per-shard JSONL files produced; writing empty $OUT_JSONL"
    : > "$OUT_JSONL"
    exit "$xargs_status"
fi

# Sort for deterministic ordering, then concat.
printf '%s\n' "${shards[@]}" | sort | xargs cat > "$OUT_JSONL"

n_lines=$(wc -l < "$OUT_JSONL")
echo "Wrote $n_lines lines to $OUT_JSONL"

exit "$xargs_status"
