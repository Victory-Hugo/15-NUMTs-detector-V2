#!/usr/bin/env python3
"""Step 04: 注释样本事件与 distinct catalog。"""

from __future__ import annotations

import argparse
import logging

import pandas as pd

from numts_frequency.annotation_utils import annotate_catalog, annotate_sample_events
from numts_frequency.io_utils import read_delim, read_length_summary, read_mt_gene_table, write_tsv

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def run(
    sample_events_tsv: str,
    catalog_tsv: str,
    sample_map_tsv: str,
    length_summary_tsv: str,
    mt_gene_table: str,
    sample_events_output_tsv: str,
    catalog_output_tsv: str,
) -> int:
    sample_events = read_delim(sample_events_tsv, sep="\t")
    catalog = read_delim(catalog_tsv, sep="\t")
    sample_map = read_delim(sample_map_tsv, sep="\t")
    length_summary = read_length_summary(length_summary_tsv)
    mt_genes = read_mt_gene_table(mt_gene_table)

    sample_events_annotated = annotate_sample_events(sample_events, length_summary)
    catalog_annotated = annotate_catalog(catalog, sample_map, sample_events_annotated, mt_genes)

    write_tsv(sample_events_annotated, sample_events_output_tsv)
    write_tsv(catalog_annotated, catalog_output_tsv)
    LOG.info("Annotation completed")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Annotate sample events and distinct catalog")
    parser.add_argument("--sample-events-tsv", required=True)
    parser.add_argument("--catalog-tsv", required=True)
    parser.add_argument("--sample-map-tsv", required=True)
    parser.add_argument("--length-summary-tsv", required=True)
    parser.add_argument("--mt-gene-table", required=True)
    parser.add_argument("--sample-events-output-tsv", required=True)
    parser.add_argument("--catalog-output-tsv", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
