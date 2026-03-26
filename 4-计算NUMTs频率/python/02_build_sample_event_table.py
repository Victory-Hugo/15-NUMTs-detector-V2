#!/usr/bin/env python3
"""Step 02: 构建样本级 NUMT 记录表。"""

from __future__ import annotations

import argparse
import logging

import pandas as pd

from numts_frequency.catalog_utils import build_sample_event_table
from numts_frequency.io_utils import read_breakpoints, read_clusters, read_delim, read_length_summary, write_tsv

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def run(
    breakpoint_tsv: str,
    cluster_tsv: str,
    length_summary_tsv: str,
    sample_inclusion_tsv: str,
    output_tsv: str,
) -> int:
    breakpoints = read_breakpoints(breakpoint_tsv)
    clusters = read_clusters(cluster_tsv)
    length_summary = read_length_summary(length_summary_tsv)
    retained_samples = set(read_delim(sample_inclusion_tsv, sep="\t", usecols=["sample_id"])["sample_id"].astype(str))
    sample_events = build_sample_event_table(clusters, breakpoints, retained_samples, length_summary)
    write_tsv(sample_events, output_tsv)
    LOG.info("Built %d sample events", len(sample_events))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build sample-level NUMT event table")
    parser.add_argument("--breakpoint-tsv", required=True)
    parser.add_argument("--cluster-tsv", required=True)
    parser.add_argument("--length-summary-tsv", required=True)
    parser.add_argument("--sample-inclusion-tsv", required=True)
    parser.add_argument("--output-tsv", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
