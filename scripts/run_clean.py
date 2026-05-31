#!/usr/bin/env python3
"""Thin invoker for pegasus-datasets `clean` filter.

Loads ``clean.py`` from a configurable location and runs ``run_clean`` against
a single input JSONL with the thresholds passed on the CLI. Used by the
`clean_jsonl` rule in the Snakefile so the cleaning step does not require
pegasus-datasets to be installed as a package -- only the source path must be
resolvable.

The default datasets root matches the layout on this machine
(/home/ubuntu/datasets); override via --datasets-root or the
DATASETS_SRC environment variable.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _load_clean(datasets_root: Path):
    src_dir = datasets_root / "src"
    if not src_dir.is_dir():
        sys.exit(f"clean.py source dir not found: {src_dir}")
    sys.path.insert(0, str(src_dir))
    try:
        from clean import CleanConfig, run_clean  # type: ignore
    except ImportError as e:
        sys.exit(f"Failed to import clean from {src_dir}: {e}")
    return CleanConfig, run_clean


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument(
        "--datasets-root",
        type=Path,
        default=Path(os.environ.get("DATASETS_SRC", "/home/ubuntu/datasets")),
        help="Root of the pegasus-datasets repo (contains src/clean.py)",
    )
    p.add_argument("--max-hunk-len", type=int, default=100)
    p.add_argument("--gap-allele-frac", type=float, default=0.7)
    p.add_argument("--ref-gap-frac", type=float, default=0.3)
    p.add_argument("--mut-density", type=float, default=0.5)
    p.add_argument("--parallel", type=int, default=8)
    p.add_argument("--chunk-size", type=int, default=2000)
    args = p.parse_args(argv)

    CleanConfig, run_clean = _load_clean(args.datasets_root)

    cfg = CleanConfig(
        input_path=args.input,
        output_path=args.output,
        manifest_path=args.manifest,
        max_hunk_len=args.max_hunk_len,
        gap_allele_frac=args.gap_allele_frac,
        ref_gap_frac=args.ref_gap_frac,
        mut_density=args.mut_density,
        parallel=args.parallel,
        chunk_size=args.chunk_size,
    )
    stats = run_clean(cfg)
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
