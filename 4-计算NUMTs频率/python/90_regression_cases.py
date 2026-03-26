#!/usr/bin/env python3
"""NUMT 频率流程关键回归用例。"""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from numts_frequency.catalog_utils import build_distinct_catalog, finalize_sample_event_table
from numts_frequency.io_utils import summarize_length_bed
from numts_frequency.stats_utils import compute_main_results, compute_mt_breakpoint_tables


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def case_discontinuous_mt_intervals() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        bed_path = Path(tmpdir) / "length.tsv.gz"
        with gzip.open(bed_path, "wt") as handle:
            handle.write("S1_chr1.100.200|a\t10\t20\n")
            handle.write("S1_chr1.100.200|b\t40\t50\n")
            handle.write("S1_chr1.100.200|c\t51\t55\n")
        summary = summarize_length_bed(
            str(bed_path),
            region_lookup={
                "S1_chr1.100.200": {
                    "sample_id": "S1",
                    "chr": "chr1",
                    "start": 100,
                    "end": 200,
                }
            },
        )
    row = summary.iloc[0]
    _assert(int(row["mt_source_envelope_span_bp"]) == 46, "envelope span mismatch")
    _assert(int(row["mt_source_span_bp"]) == 27, "union span mismatch")
    _assert(bool(row["complex_mt_interval"]) is True, "complex interval flag mismatch")
    _assert(int(row["mt_primary_fragment_span_bp"]) == 16, "primary fragment span mismatch")


def case_finalize_sample_events() -> None:
    sample_events = pd.DataFrame(
        [
            {
                "sample_event_id": "S1|E1",
                "sample_id": "S1",
                "cluster_id_raw": "E1",
                "chr": "chr1",
                "cluster_start": 100,
                "cluster_end": 200,
                "cluster_width_bp": 101,
                "cluster_support_n": 2,
                "cluster_reads_raw": "[]",
                "cluster_mt_positions_raw": "[]",
                "nuclear_core_start": 100,
                "nuclear_core_end": 200,
                "nuclear_core_midpoint": 150.0,
                "mt_core_start": 300,
                "mt_core_end": 350,
                "mt_core_midpoint": 325.0,
                "mt_core_span_bp": 51,
                "nuclear_left_breakpoint_pos": 110.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "nuclear_anchor_midpoint": np.nan,
                "mt_breakpoint_start": 305.0,
                "mt_breakpoint_end": np.nan,
                "supporting_breakpoint_rows": 1,
                "length_bed_region_key": "S1_chr1.100.200",
                "length_bed_match_status": "matched",
            },
            {
                "sample_event_id": "S2|E2",
                "sample_id": "S2",
                "cluster_id_raw": "E2",
                "chr": "chr1",
                "cluster_start": 300,
                "cluster_end": 400,
                "cluster_width_bp": 101,
                "cluster_support_n": 2,
                "cluster_reads_raw": "[]",
                "cluster_mt_positions_raw": "[]",
                "nuclear_core_start": 300,
                "nuclear_core_end": 400,
                "nuclear_core_midpoint": 350.0,
                "mt_core_start": 500,
                "mt_core_end": 580,
                "mt_core_midpoint": 540.0,
                "mt_core_span_bp": 81,
                "nuclear_left_breakpoint_pos": 310.0,
                "nuclear_right_breakpoint_pos": 320.0,
                "nuclear_anchor_midpoint": 315.0,
                "mt_breakpoint_start": 505.0,
                "mt_breakpoint_end": 515.0,
                "supporting_breakpoint_rows": 2,
                "length_bed_region_key": "S2_chr1.300.400",
                "length_bed_match_status": "matched",
            },
        ]
    )
    length_summary = pd.DataFrame(
        [
            {
                "region_key": "S1_chr1.100.200",
                "mt_source_min": 300,
                "mt_source_max": 350,
                "mt_source_envelope_span_bp": 51,
                "mt_source_span_bp": 51,
                "length_bed_read_count": 1,
                "merged_interval_count": 1,
                "complex_mt_interval": False,
                "mt_primary_fragment_start": 300,
                "mt_primary_fragment_end": 350,
                "mt_primary_fragment_span_bp": 51,
                "mt_merged_intervals": "300-350",
            },
            {
                "region_key": "S2_chr1.300.400",
                "mt_source_min": 500,
                "mt_source_max": 580,
                "mt_source_envelope_span_bp": 81,
                "mt_source_span_bp": 41,
                "length_bed_read_count": 2,
                "merged_interval_count": 2,
                "complex_mt_interval": True,
                "mt_primary_fragment_start": 500,
                "mt_primary_fragment_end": 520,
                "mt_primary_fragment_span_bp": 21,
                "mt_merged_intervals": "500-520;560-580",
            },
        ]
    )
    finalized = finalize_sample_event_table(sample_events, length_summary)
    row_single = finalized.loc[finalized["sample_event_id"] == "S1|E1"].iloc[0]
    row_complex = finalized.loc[finalized["sample_event_id"] == "S2|E2"].iloc[0]
    _assert(bool(row_single["eligible_for_main_analysis"]) is True, "single-sided breakpoint event should be eligible")
    _assert(int(row_single["mt_length_for_main_bp"]) == 51, "single-sided event should fall back to union span")
    _assert(bool(row_complex["eligible_for_main_analysis"]) is False, "complex event should be excluded")
    _assert(row_complex["main_analysis_exclusion_reason"] == "complex_mt_interval", "complex reason mismatch")


def case_distinct_interval_clustering_and_stats() -> None:
    sample_events = pd.DataFrame(
        [
            {
                "sample_event_id": "S1|A",
                "sample_id": "S1",
                "chr": "chr1",
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 100.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 500.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 31.0,
                "mt_length_for_main_bp": 31.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S2|B",
                "sample_id": "S2",
                "chr": "chr1",
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 130.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 530.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 35.0,
                "mt_length_for_main_bp": 35.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S1|C",
                "sample_id": "S1",
                "chr": "chr1",
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 5000.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 9000.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 20.0,
                "mt_length_for_main_bp": 20.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S3|D",
                "sample_id": "S3",
                "chr": "chr1",
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 6000.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 9050.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 22.0,
                "mt_length_for_main_bp": 22.0,
                "eligible_for_main_analysis": True,
            },
        ]
    )
    catalog, sample_map = build_distinct_catalog(
        sample_events=sample_events,
        retained_sample_count=3,
        nuclear_merge_gap_bp=50,
        mt_merge_gap_bp=50,
        distinct_id_prefix="NUMT",
        prevalent_min=0.10,
        common_min=0.01,
        rare_min=0.001,
    )
    _assert(len(catalog) == 3, "interval clustering should yield three distinct NUMTs")
    _assert(int(catalog["carrier_count"].max()) == 2, "carrier_count max mismatch")
    _assert(int(sample_map["distinct_numt_id"].nunique()) == 3, "sample_map distinct mismatch")
    merged_pair = catalog.sort_values("carrier_count", ascending=False).iloc[0]
    _assert(float(merged_pair["merged_nuclear_breakpoint_start"]) == 100.0, "merged nuclear start mismatch")
    _assert(float(merged_pair["merged_nuclear_breakpoint_end"]) == 100.0, "merged nuclear end mismatch")
    _assert(float(merged_pair["nuclear_cluster_range_end"]) == 130.0, "nuclear cluster range end mismatch")
    _assert(float(merged_pair["distinct_nuclear_start"]) == float(merged_pair["merged_nuclear_breakpoint_start"]), "distinct start should follow merged breakpoint start")
    _assert(float(merged_pair["distinct_nuclear_end"]) == float(merged_pair["merged_nuclear_breakpoint_end"]), "distinct end should follow merged breakpoint end")

    catalog = catalog.assign(
        mt_primary_gene=["ND1", "ND2", "ND3"][: len(catalog)],
        mt_primary_gene_class=["protein_coding"] * len(catalog),
        mt_overlapping_genes=["ND1", "ND2", "ND3"][: len(catalog)],
    )
    sample_events = sample_events.assign(
        main_analysis_exclusion_reason=["", "", "", ""],
    )
    retained_samples = pd.DataFrame({"sample_id": ["S1", "S2", "S3"]})
    results = compute_main_results(catalog, sample_map, sample_events, retained_samples)
    freq = results["frequency_class_summary"]
    carrier_ratio = results["carrier_ratio_by_frequency_class"]
    _assert((freq["carrier_fraction"] <= 1).all(), "carrier_fraction must be <= 1")
    _assert((carrier_ratio["sample_ratio"] <= 1).all(), "sample_ratio must be <= 1")

    mt_tables = compute_mt_breakpoint_tables(catalog)
    cumulative = mt_tables["mt_breakpoint_cumulative"]
    overall = cumulative[["mt_position", "overall_cumulative_fraction"]].drop_duplicates()
    _assert(abs(float(overall["overall_cumulative_fraction"].iloc[-1]) - 1.0) < 1e-9, "overall cumulative fraction must end at 1")


def main() -> int:
    case_discontinuous_mt_intervals()
    case_finalize_sample_events()
    case_distinct_interval_clustering_and_stats()
    print("All regression cases passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
