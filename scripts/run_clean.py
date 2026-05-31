#!/usr/bin/env python3
"""CLI wrapper around scripts/clean_lib.py for the `clean_jsonl` Snakefile rule.

Imports the vendored ``clean_lib`` module (sibling file in scripts/) and runs
``run_clean`` against a single input JSONL with the thresholds passed on the
CLI.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# clean_lib.py is a sibling of this file under scripts/.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from clean_lib import CleanConfig, run_clean  # type: ignore  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument("--max-hunk-len", type=int, default=100)
    p.add_argument("--gap-allele-frac", type=float, default=0.7)
    p.add_argument("--ref-gap-frac", type=float, default=0.3)
    p.add_argument("--mut-density", type=float, default=0.5)
    p.add_argument("--parallel", type=int, default=8)
    p.add_argument("--chunk-size", type=int, default=2000)
    args = p.parse_args(argv)

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
