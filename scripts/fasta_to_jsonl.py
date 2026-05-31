#!/usr/bin/env python3
"""Convert FASTA trajectory files to JSON-structured SDIFF JSONL format.

Vendored from pegasus-evals (scripts/fasta_to_jsonl.py). Self-contained —
only external dep is numpy. Invoked by the `jsonl_mode_split` rule in
Snakefile via scripts/shards_to_jsonl.sh.


Standalone script (no dFactory imports) for converting FASTA files into
the json_sdiff format used by pegasus-evals for eval and inference.

Output format (with --noprompt):
    {"messages":[{"role":"user","content":""},{"role":"assistant","content":"REFSEQ\n@@N\n{\"h\":[...]}"}]}

Usage:
    python scripts/fasta_to_jsonl.py \
        --input /path/to/file.fasta \
        --output /path/to/output.jsonl \
        --max-raw-len 7000 \
        --noprompt
"""

import argparse
import json
import glob as glob_mod
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Core functions (from dFactory/scripts/preprocess_data/fasta_to_json_sdiff.py)
# ---------------------------------------------------------------------------


def parse_fasta_file(filepath: str) -> List[Tuple[str, str]]:
    """Parse a FASTA file into a list of (header, sequence) tuples, in order."""
    sequences: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_sequence: List[str] = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_sequence)))
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line.upper())

        if current_header is not None:
            sequences.append((current_header, "".join(current_sequence)))

    return sequences


def find_unique_context_size(
    ref_seq: str,
    pos: int,
    ref_allele: str,
    min_context: int = 8,
    max_context: int = 16,
) -> int:
    """Find minimum left-context size that uniquely identifies position."""
    for ctx_size in range(min_context, max_context + 1):
        left_start = max(0, pos - ctx_size)
        left = ref_seq[left_start:pos]

        search_pattern = left + ref_allele
        occurrences = ref_seq.count(search_pattern)

        if occurrences == 1:
            return ctx_size

    return max_context


def compute_variants_aligned(
    ref_seq: str,
    alt_seq: str,
) -> List[Tuple[int, str, str, str]]:
    """Compute variants between two aligned sequences (same length)."""
    if len(ref_seq) != len(alt_seq):
        raise ValueError("Sequences must be same length for aligned comparison")

    ref_arr = np.array(list(ref_seq), dtype="U1")
    alt_arr = np.array(list(alt_seq), dtype="U1")

    diff_mask = ref_arr != alt_arr
    diff_positions = np.where(diff_mask)[0]

    variants = []
    i = 0
    while i < len(diff_positions):
        start_pos = int(diff_positions[i])

        end_pos = start_pos
        while i + 1 < len(diff_positions) and diff_positions[i + 1] == end_pos + 1:
            i += 1
            end_pos = int(diff_positions[i])

        ref_allele = ref_seq[start_pos : end_pos + 1]
        alt_allele = alt_seq[start_pos : end_pos + 1]

        var_type = "SNP" if len(ref_allele) == 1 else "MNP"
        variants.append((start_pos, ref_allele, alt_allele, var_type))
        i += 1

    return variants


def compute_variants(
    ref_seq: str,
    alt_seq: str,
) -> List[Tuple[int, str, str, str]]:
    """Compute variants between reference and alternate sequence."""
    return compute_variants_aligned(ref_seq, alt_seq)


def format_example(content: str, prompt: str) -> Dict:
    """Format content into SFT-style conversation JSON."""
    return {
        "messages": [
            {"role": "user", "content": prompt},
            {"role": "assistant", "content": content},
        ],
    }


def format_json_sdiff_hunk(
    pos: int,
    ref_allele: str,
    alt_allele: str,
    var_type: str,
    ref_seq: str,
    context_size: int = 8,
    adaptive_context: bool = True,
) -> Dict[str, str]:
    """Format a single mutation as a JSON hunk dict."""
    if adaptive_context and ref_allele:
        ctx = find_unique_context_size(
            ref_seq, pos, ref_allele, min_context=context_size, max_context=12
        )
    else:
        ctx = context_size

    left_start = max(0, pos - ctx)
    left = ref_seq[left_start:pos]

    return {"c": left, "r": ref_allele, "a": alt_allele}


def format_json_sdiff_generation(
    variants: List[Tuple[int, str, str, str]],
    ref_seq: str,
    context_size: int = 8,
    adaptive_context: bool = True,
) -> str:
    """Format all variants for a generation as a JSON string."""
    hunks = []
    for pos, ref_allele, alt_allele, var_type in variants:
        hunk = format_json_sdiff_hunk(
            pos, ref_allele, alt_allele, var_type,
            ref_seq, context_size, adaptive_context,
        )
        hunks.append(hunk)

    return json.dumps({"h": hunks}, ensure_ascii=False, separators=(",", ":"))


def process_fasta_to_json_sdiff(
    sequences: List[Tuple[str, str]],
    max_len: Optional[int] = None,
    max_raw_len: Optional[int] = None,
    context_size: int = 8,
    adaptive_context: bool = True,
) -> Optional[str]:
    """Process sequences into JSON-structured SDIFF format.

    Invariant: every @@N block's hunks are computed against the root sequence
    (sequences[start_idx]), not against the previous generation. Downstream
    consumers must apply hunks to the root to reconstruct any generation.
    """
    if len(sequences) < 2:
        return None

    k = len(sequences)

    if max_len is not None and k > max_len:
        start_idx = k - max_len
    else:
        start_idx = 0

    if start_idx >= k - 1:
        return None

    ref_header, ref_seq = sequences[start_idx]

    if max_raw_len is not None and len(ref_seq) > max_raw_len:
        return None

    current_len = len(ref_seq)
    parts = [ref_seq]
    gen_count = 0

    for idx in range(start_idx + 1, k):
        alt_header, alt_seq = sequences[idx]

        variants = compute_variants(ref_seq, alt_seq)

        if not variants:
            continue

        gen_body = format_json_sdiff_generation(
            variants, ref_seq, context_size, adaptive_context
        )

        gen_id = len(variants)
        header_str = f"@@{gen_id}"
        added_len = 1 + len(header_str) + 1 + len(gen_body)

        if max_raw_len is not None and current_len + added_len > max_raw_len:
            break

        parts.append(header_str)
        parts.append(gen_body)
        current_len += added_len
        gen_count += 1

    if gen_count == 0:
        return None

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert FASTA trajectory files to JSON-structured SDIFF JSONL format."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Single FASTA file or directory of FASTA files",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output JSONL path (one JSON line per FASTA file)",
    )
    parser.add_argument(
        "--max-raw-len",
        type=int,
        default=14000,
        help="Max output string length in characters (default: 14000)",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=None,
        help="Max trajectory length (number of sequences); triggers truncation",
    )
    parser.add_argument(
        "--noprompt",
        action="store_true",
        help='Set user content to "" (default for inference)',
    )
    parser.add_argument(
        "--context-size",
        type=int,
        default=16,
        help="Left-context bases for each mutation (default: 16)",
    )
    parser.add_argument(
        "--file-pattern",
        type=str,
        default=None,
        help="Glob pattern when input is a directory (default: *.fasta + *.fa)",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    prompt = "" if args.noprompt else "evodiff-json-sdiff"

    # Collect FASTA files
    if input_path.is_file():
        fasta_files = [input_path]
    elif input_path.is_dir():
        if args.file_pattern:
            fasta_files = sorted(input_path.glob(args.file_pattern))
        else:
            fasta_files = sorted(
                list(input_path.glob("*.fasta")) + list(input_path.glob("*.fa"))
            )
    else:
        raise FileNotFoundError(f"Input not found: {input_path}")

    if not fasta_files:
        print(f"No FASTA files found in {input_path}")
        return

    os.makedirs(Path(args.output).parent, exist_ok=True)

    total = 0
    skipped = 0

    with open(args.output, "w", encoding="utf-8") as out_f:
        for fasta_file in fasta_files:
            try:
                sequences = parse_fasta_file(str(fasta_file))
                if not sequences:
                    skipped += 1
                    continue

                sdiff_content = process_fasta_to_json_sdiff(
                    sequences,
                    max_len=args.max_len,
                    max_raw_len=args.max_raw_len,
                    context_size=args.context_size,
                )
                if not sdiff_content:
                    skipped += 1
                    continue

                example = format_example(sdiff_content, prompt)
                out_f.write(
                    json.dumps(example, ensure_ascii=False, sort_keys=True) + "\n"
                )
                total += 1

            except Exception as e:
                print(f"Error processing {fasta_file.name}: {e}")
                skipped += 1

    print(f"Converted {total} files ({skipped} skipped) -> {args.output}")


if __name__ == "__main__":
    main()
