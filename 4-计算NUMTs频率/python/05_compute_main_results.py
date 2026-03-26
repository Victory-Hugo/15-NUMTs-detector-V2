#!/usr/bin/env python3
"""Step 05: 主结果统计。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from numts_frequency.io_utils import ensure_dir, read_delim, write_tsv
from numts_frequency.stats_utils import compute_main_results

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def run(
    catalog_tsv: str,
    sample_map_tsv: str,
    sample_events_tsv: str,
    sample_inclusion_tsv: str,
    output_dir: str,
) -> int:
    ensure_dir(output_dir)
    catalog = read_delim(catalog_tsv, sep="\t")
    sample_map = read_delim(sample_map_tsv, sep="\t")
    sample_events = read_delim(sample_events_tsv, sep="\t")
    sample_inclusion = read_delim(sample_inclusion_tsv, sep="\t")

    results = compute_main_results(catalog, sample_map, sample_events, sample_inclusion)
    for name, df in results.items():
        write_tsv(df, str(Path(output_dir) / f"{name}.tsv"))
    LOG.info("Main results completed")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute main NUMT frequency results")
    parser.add_argument("--catalog-tsv", required=True)
    parser.add_argument("--sample-map-tsv", required=True)
    parser.add_argument("--sample-events-tsv", required=True)
    parser.add_argument("--sample-inclusion-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
