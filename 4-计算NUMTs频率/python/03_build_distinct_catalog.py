#!/usr/bin/env python3
"""Step 03: 构建 distinct NUMT catalog。"""

from __future__ import annotations

import argparse
import logging

import pandas as pd

from numts_frequency.catalog_utils import build_distinct_catalog
from numts_frequency.io_utils import read_delim, write_tsv

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def run(
    sample_events_tsv: str,
    sample_inclusion_tsv: str,
    nuclear_merge_gap_bp: int,
    mt_merge_gap_bp: int,
    distinct_id_prefix: str,
    prevalent_min: float,
    common_min: float,
    rare_min: float,
    catalog_tsv: str,
    sample_map_tsv: str,
) -> int:
    sample_events = read_delim(sample_events_tsv, sep="\t")
    retained_sample_count = int(read_delim(sample_inclusion_tsv, sep="\t", usecols=["sample_id"])["sample_id"].nunique())
    catalog, sample_map = build_distinct_catalog(
        sample_events=sample_events,
        retained_sample_count=retained_sample_count,
        nuclear_merge_gap_bp=nuclear_merge_gap_bp,
        mt_merge_gap_bp=mt_merge_gap_bp,
        distinct_id_prefix=distinct_id_prefix,
        prevalent_min=prevalent_min,
        common_min=common_min,
        rare_min=rare_min,
    )
    write_tsv(catalog, catalog_tsv)
    write_tsv(sample_map, sample_map_tsv)
    LOG.info("Built %d distinct NUMTs", len(catalog))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build distinct NUMT catalog")
    parser.add_argument("--sample-events-tsv", required=True)
    parser.add_argument("--sample-inclusion-tsv", required=True)
    parser.add_argument("--nuclear-merge-gap-bp", required=True, type=int)
    parser.add_argument("--mt-merge-gap-bp", required=True, type=int)
    parser.add_argument("--distinct-id-prefix", required=True)
    parser.add_argument("--prevalent-min", required=True, type=float)
    parser.add_argument("--common-min", required=True, type=float)
    parser.add_argument("--rare-min", required=True, type=float)
    parser.add_argument("--catalog-tsv", required=True)
    parser.add_argument("--sample-map-tsv", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
