#!/usr/bin/env python3
"""Step 06: mtDNA breakpoint 分布统计。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from numts_frequency.io_utils import ensure_dir, read_delim, write_tsv
from numts_frequency.stats_utils import compute_mt_breakpoint_tables

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def run(catalog_tsv: str, output_dir: str) -> int:
    ensure_dir(output_dir)
    catalog = read_delim(catalog_tsv, sep="\t")
    tables = compute_mt_breakpoint_tables(catalog)
    for name, df in tables.items():
        write_tsv(df, str(Path(output_dir) / f"{name}.tsv"))
    LOG.info("MT breakpoint results completed")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute mtDNA breakpoint statistics")
    parser.add_argument("--catalog-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
