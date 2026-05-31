"""Streaming JSONL filter for json_sdiff pretrain records.

Vendored from pegasus-datasets (src/clean.py). Self-contained — stdlib only.
Invoked via scripts/run_clean.py and the `clean_jsonl` rule in Snakefile.


Applies four record-level quality filters (independent — any failure drops):

  1. max_hunk_len:    drop if any hunk has ``len(r) > T`` or ``len(a) > T``
  2. gap_allele_frac: drop if ``>T`` fraction of hunks have ``-`` in r or a
  3. ref_gap_frac:    drop if ``ref.count('-') / len(ref) > T``
  4. mut_density:     drop if ``n_hunks / ungapped_ref_len > T``

Input/output are JSONL files (one record per line). Workers parse and
filter chunks in parallel; the parent process writes survivors and
accumulates a manifest CSV of drop counts by reason.
"""
from __future__ import annotations

import csv
import json
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

HUNK_RE = re.compile(r'\{"c":"([^"]*)","r":"([^"]*)","a":"([^"]*)"\}')
GEN_HEADER_RE = re.compile(r'^@@(\d+)$', re.MULTILINE)


@dataclass
class CleanConfig:
    input_path: Path
    output_path: Path
    manifest_path: Path | None = None

    max_hunk_len: int = 100
    """Threshold on the longest single-hunk ref/alt allele (chars)."""

    gap_allele_frac: float = 0.7
    """Threshold on fraction of hunks that contain '-' in r or a."""

    ref_gap_frac: float = 0.3
    """Threshold on ref's gap fraction (chars of '-' / len(ref))."""

    mut_density: float = 0.5
    """Threshold on max-per-generation hunk count / ungapped_ref_len.
    For pairwise (one generation) this equals total_hunks / ungapped_ref_len.
    For forwards, using max-per-gen avoids double-counting positions
    re-mutated across generations."""

    parallel: int = 8
    chunk_size: int = 2000

    def __post_init__(self) -> None:
        self.input_path = Path(self.input_path)
        self.output_path = Path(self.output_path)
        if self.manifest_path is not None:
            self.manifest_path = Path(self.manifest_path)


DROP_REASONS = (
    "max_hunk_len",
    "gap_allele_frac",
    "ref_gap_frac",
    "mut_density",
    "parse_error",
)


def _evaluate(
    line: str,
    *,
    max_hunk_len: int,
    gap_allele_frac: float,
    ref_gap_frac: float,
    mut_density: float,
) -> tuple[bool, str | None]:
    """Return (keep, drop_reason). drop_reason is None when keep=True.

    Filter order is fixed (ref_gap -> max_hunk -> gap_allele -> mut_density)
    so manifest counts reflect the *first* failing check.
    """
    try:
        rec = json.loads(line)
        content = rec["messages"][1]["content"]
    except Exception:
        return False, "parse_error"

    nl = content.find("\n")
    if nl < 0:
        ref = content
        rest = ""
    else:
        ref = content[:nl]
        rest = content[nl + 1 :]

    ref_len = len(ref)
    gap_count = ref.count("-")

    # 3. ref_gap_frac (cheap — check first)
    if ref_len > 0 and gap_count / ref_len > ref_gap_frac:
        return False, "ref_gap_frac"

    hunks = HUNK_RE.findall(rest)
    n_hunks = len(hunks)
    if n_hunks == 0:
        # No mutations -> no quality signal to filter on; keep.
        return True, None

    # 4. max_hunk_len
    for _c, r, a in hunks:
        if len(r) > max_hunk_len or len(a) > max_hunk_len:
            return False, "max_hunk_len"

    # 2. gap_allele_frac
    n_gap_allele = 0
    for _c, r, a in hunks:
        if "-" in r or "-" in a:
            n_gap_allele += 1
    if n_gap_allele / n_hunks > gap_allele_frac:
        return False, "gap_allele_frac"

    # 5. mut_density: use max-per-generation hunk count (the @@N header value
    # IS the variant count for that generation). For pairwise this equals
    # n_hunks; for forwards it picks the worst single step.
    ungapped = ref_len - gap_count
    if ungapped > 0:
        gen_counts = [int(x) for x in GEN_HEADER_RE.findall(rest)]
        max_gen_muts = max(gen_counts) if gen_counts else n_hunks
        if max_gen_muts / ungapped > mut_density:
            return False, "mut_density"

    return True, None


def _process_chunk(args: tuple) -> tuple[list[str], dict[str, int]]:
    lines, thresholds = args
    kept: list[str] = []
    drops: dict[str, int] = {r: 0 for r in DROP_REASONS}
    for line in lines:
        keep, reason = _evaluate(line, **thresholds)
        if keep:
            kept.append(line)
        else:
            drops[reason] = drops.get(reason, 0) + 1
    return kept, drops


def _chunked_lines(path: Path, chunk_size: int):
    buf: list[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            buf.append(line)
            if len(buf) >= chunk_size:
                yield buf
                buf = []
    if buf:
        yield buf


def run_clean(cfg: CleanConfig) -> dict:
    cfg.output_path.parent.mkdir(parents=True, exist_ok=True)
    if cfg.manifest_path is not None:
        cfg.manifest_path.parent.mkdir(parents=True, exist_ok=True)

    thresholds = {
        "max_hunk_len": cfg.max_hunk_len,
        "gap_allele_frac": cfg.gap_allele_frac,
        "ref_gap_frac": cfg.ref_gap_frac,
        "mut_density": cfg.mut_density,
    }

    n_in = 0
    n_kept = 0
    drops: dict[str, int] = {r: 0 for r in DROP_REASONS}

    print(f"[clean] input:   {cfg.input_path}")
    print(f"[clean] output:  {cfg.output_path}")
    print(f"[clean] thresholds: {thresholds}")
    print(f"[clean] parallel={cfg.parallel} chunk_size={cfg.chunk_size}")

    chunks = ((lines, thresholds) for lines in _chunked_lines(cfg.input_path, cfg.chunk_size))

    with cfg.output_path.open("w", encoding="utf-8") as out:
        if cfg.parallel <= 1:
            for chunk in chunks:
                kept, chunk_drops = _process_chunk(chunk)
                for line in kept:
                    out.write(line)
                n_in += len(chunk[0])
                n_kept += len(kept)
                for k, v in chunk_drops.items():
                    drops[k] += v
                if n_in % 200000 == 0:
                    print(f"[clean] {n_in:>10}  kept={n_kept}  drops={sum(drops.values())}")
        else:
            with ProcessPoolExecutor(max_workers=cfg.parallel) as pool:
                # imap_unordered for memory; chunksize=1 since each chunk is already big
                for kept, chunk_drops in pool.map(_process_chunk, chunks, chunksize=1):
                    for line in kept:
                        out.write(line)
                    n_kept += len(kept)
                    chunk_total = sum(chunk_drops.values()) + len(kept)
                    n_in += chunk_total
                    for k, v in chunk_drops.items():
                        drops[k] += v
                    if n_in // cfg.chunk_size % 100 == 0:
                        print(f"[clean] {n_in:>10}  kept={n_kept}  drops={sum(drops.values())}")

    if cfg.manifest_path is not None:
        with cfg.manifest_path.open("w", newline="", encoding="utf-8") as mf:
            w = csv.writer(mf)
            w.writerow(["reason", "count", "fraction_of_input"])
            for r in DROP_REASONS:
                frac = drops[r] / n_in if n_in else 0.0
                w.writerow([r, drops[r], f"{frac:.6f}"])
            w.writerow(["kept", n_kept, f"{(n_kept / n_in if n_in else 0.0):.6f}"])
            w.writerow(["total", n_in, "1.000000"])

    summary = {
        "input": str(cfg.input_path),
        "output": str(cfg.output_path),
        "n_in": n_in,
        "n_kept": n_kept,
        "kept_fraction": (n_kept / n_in) if n_in else 0.0,
        "drops": drops,
    }
    print(f"[clean] done: kept {n_kept}/{n_in} ({100 * (n_kept / n_in if n_in else 0):.2f}%)")
    for r in DROP_REASONS:
        if drops[r] > 0:
            print(f"  drop {r}: {drops[r]} ({100 * drops[r] / n_in:.2f}%)")
    return summary
