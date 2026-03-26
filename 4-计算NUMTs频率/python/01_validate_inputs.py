#!/usr/bin/env python3
"""Step 01: 输入一致性检查。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from numts_frequency.io_utils import (
    ensure_dir,
    read_breakpoints,
    read_clusters,
    read_length_summary,
    read_meta,
    subset_length_summary,
    summarize_length_bed,
    write_tsv,
)
from numts_frequency.qc_utils import (
    build_event_traceability,
    build_input_summary,
    build_sample_filters,
    build_validation_report,
    write_qc_outputs,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def run(
    breakpoint_tsv: str,
    cluster_tsv: str,
    length_bed_gz: str,
    meta_tsv: str,
    meta_id_col: str,
    meta_qc_col: str,
    meta_qc_pass_value: str,
    output_dir: str,
    length_summary_cache_tsv: str | None = None,
) -> int:
    ensure_dir(output_dir)
    breakpoints = read_breakpoints(breakpoint_tsv)
    clusters = read_clusters(cluster_tsv)
    meta = read_meta(meta_tsv, meta_id_col)
    inclusion, exclusion = build_sample_filters(breakpoints, clusters, meta, meta_qc_col, meta_qc_pass_value)
    retained_ids = set(inclusion["sample_id"].astype(str))
    retained_clusters = clusters[clusters["sample_id"].isin(retained_ids)].copy()
    region_lookup = {
        row.region_key: {
            "sample_id": row.sample_id,
            "chr": row.chr,
            "start": int(row.start),
            "end": int(row.end),
        }
        for row in retained_clusters.itertuples(index=False)
    }
    if length_summary_cache_tsv and Path(length_summary_cache_tsv).exists():
        LOG.info("Reusing cached full length summary: %s", length_summary_cache_tsv)
        full_length_summary = read_length_summary(length_summary_cache_tsv)
    else:
        LOG.info("Building full length summary cache from: %s", length_bed_gz)
        full_length_summary = summarize_length_bed(length_bed_gz)
        if length_summary_cache_tsv:
            write_tsv(full_length_summary, length_summary_cache_tsv)
            LOG.info("Wrote full length summary cache: %s", length_summary_cache_tsv)
    length_summary = subset_length_summary(full_length_summary, region_lookup)
    input_summary = build_input_summary(breakpoints, clusters, meta, length_summary, meta_qc_col, meta_qc_pass_value)
    traceability = build_event_traceability(breakpoints, clusters, retained_ids)
    write_qc_outputs(output_dir, input_summary, inclusion, exclusion, traceability, length_summary)
    build_validation_report(
        str(Path(output_dir) / "05-input_validation_report.md"),
        input_summary,
        inclusion,
        exclusion,
        traceability,
        breakpoints,
        clusters,
        meta,
    )
    LOG.info("Validation completed: retained %d samples", len(inclusion))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate NUMTs frequency pipeline inputs")
    parser.add_argument("--breakpoint-tsv", required=True)
    parser.add_argument("--cluster-tsv", required=True)
    parser.add_argument("--length-bed-gz", required=True)
    parser.add_argument("--meta-tsv", required=True)
    parser.add_argument("--meta-id-col", required=True)
    parser.add_argument("--meta-qc-col", required=True)
    parser.add_argument("--meta-qc-pass-value", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--length-summary-cache-tsv", required=False, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
