#!/usr/bin/env python3
"""NUMT 频率流程关键回归用例。"""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from numts_frequency.catalog_utils import build_distinct_catalog, build_sample_event_table, finalize_sample_event_table
from numts_frequency.io_utils import summarize_length_bed
from numts_frequency.stats_utils import compute_main_results, compute_mt_breakpoint_tables
from numts_frequency.annotation_utils import annotate_catalog


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
            {
                "sample_event_id": "S3|E3",
                "sample_id": "S3",
                "cluster_id_raw": "E3",
                "chr": "chr2",
                "cluster_start": 500,
                "cluster_end": 620,
                "cluster_width_bp": 121,
                "cluster_support_n": 2,
                "cluster_reads_raw": "[]",
                "cluster_mt_positions_raw": "[]",
                "nuclear_core_start": 500,
                "nuclear_core_end": 620,
                "nuclear_core_midpoint": 560.0,
                "mt_core_start": 1000,
                "mt_core_end": 9000,
                "mt_core_midpoint": 5000.0,
                "mt_core_span_bp": 8001,
                "nuclear_left_breakpoint_pos": 510.0,
                "nuclear_right_breakpoint_pos": 610.0,
                "nuclear_anchor_midpoint": 560.0,
                "mt_breakpoint_start": 1500.0,
                "mt_breakpoint_end": 1401.0,
                "supporting_breakpoint_rows": 2,
                "length_bed_region_key": "S3_chr2.500.620",
                "length_bed_match_status": "not_matched",
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
    row_reverse_mt = finalized.loc[finalized["sample_event_id"] == "S3|E3"].iloc[0]
    _assert(bool(row_single["eligible_for_main_analysis"]) is True, "single-sided breakpoint event should be eligible")
    _assert(int(row_single["mt_length_for_main_bp"]) == 51, "single-sided event should fall back to union span")
    _assert(bool(row_complex["eligible_for_main_analysis"]) is False, "complex event should be excluded")
    _assert(row_complex["main_analysis_exclusion_reason"] == "complex_mt_interval", "complex reason mismatch")
    _assert(int(row_reverse_mt["mt_breakpoint_span_bp"]) == 100, "reverse-strand mt breakpoint span should be positive")
    _assert(int(row_reverse_mt["mt_source_span_bp"]) == 100, "not_matched event should fall back to mt breakpoint span")
    _assert(int(row_reverse_mt["mt_length_for_main_bp"]) == 100, "main length should use positive mt breakpoint span")


def case_distinct_interval_clustering_and_stats() -> None:
    sample_events = pd.DataFrame(
        [
            {
                "sample_event_id": "S1|A",
                "sample_id": "S1",
                "chr": "chr1",
                "cluster_start": 100,
                "cluster_end": 180,
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 100.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 500.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 31.0,
                "mt_primary_fragment_span_bp": 31.0,
                "mt_length_for_main_bp": 31.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S2|B",
                "sample_id": "S2",
                "chr": "chr1",
                "cluster_start": 150,
                "cluster_end": 230,
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 220.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 530.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 35.0,
                "mt_primary_fragment_span_bp": 35.0,
                "mt_length_for_main_bp": 35.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S1|C",
                "sample_id": "S1",
                "chr": "chr1",
                "cluster_start": 5000,
                "cluster_end": 5080,
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 5000.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 9000.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 20.0,
                "mt_primary_fragment_span_bp": 20.0,
                "mt_length_for_main_bp": 20.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S3|D",
                "sample_id": "S3",
                "chr": "chr1",
                "cluster_start": 6000,
                "cluster_end": 6080,
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 6000.0,
                "nuclear_right_breakpoint_pos": np.nan,
                "mt_breakpoint_start": 9050.0,
                "mt_breakpoint_end": np.nan,
                "mt_source_span_bp": 22.0,
                "mt_primary_fragment_span_bp": 22.0,
                "mt_length_for_main_bp": 22.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S4|E",
                "sample_id": "S4",
                "chr": "chr1",
                "cluster_start": 7000,
                "cluster_end": 7080,
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 7000.0,
                "nuclear_right_breakpoint_pos": 7010.0,
                "mt_breakpoint_start": 1000.0,
                "mt_breakpoint_end": 1099.0,
                "mt_source_span_bp": 100.0,
                "mt_primary_fragment_span_bp": 100.0,
                "mt_length_for_main_bp": 100.0,
                "eligible_for_main_analysis": True,
            },
            {
                "sample_event_id": "S5|F",
                "sample_id": "S5",
                "chr": "chr1",
                "cluster_start": 7020,
                "cluster_end": 7100,
                "cluster_width_bp": 100,
                "nuclear_left_breakpoint_pos": 7025.0,
                "nuclear_right_breakpoint_pos": 7035.0,
                "mt_breakpoint_start": 1005.0,
                "mt_breakpoint_end": 1104.0,
                "mt_source_span_bp": 100.0,
                "mt_primary_fragment_span_bp": 101.0,
                "mt_length_for_main_bp": 101.0,
                "eligible_for_main_analysis": True,
            },
        ]
    )
    catalog, sample_map = build_distinct_catalog(
        sample_events=sample_events,
        retained_sample_count=5,
        nuclear_merge_gap_bp=50,
        mt_merge_gap_bp=50,
        distinct_id_prefix="NUMT",
        prevalent_min=0.10,
        common_min=0.01,
        rare_min=0.001,
    )
    _assert(len(catalog) == 4, "interval clustering should yield four distinct NUMTs")
    _assert(int(catalog["carrier_count"].max()) == 2, "carrier_count max mismatch")
    _assert(int(sample_map["distinct_numt_id"].nunique()) == 4, "sample_map distinct mismatch")
    merged_pair = catalog.sort_values("carrier_count", ascending=False).iloc[0]
    _assert(float(merged_pair["nuclear_cluster_range_start"]) == 100.0, "nuclear cluster range start mismatch")
    _assert(float(merged_pair["nuclear_cluster_range_end"]) == 230.0, "nuclear cluster range end mismatch")
    _assert(float(merged_pair["merged_nuclear_breakpoint_start"]) == 100.0, "merged nuclear start mismatch")
    _assert(float(merged_pair["merged_nuclear_breakpoint_end"]) == 100.0, "merged nuclear end mismatch")
    _assert(
        float(catalog.loc[catalog["nuclear_cluster_range_start"] == 7000.0, "median_mt_length_for_main_bp"].iloc[0]) == 100.0,
        "representative length should remain integer after aggregation",
    )

    catalog = catalog.assign(
        mt_primary_gene=["ND1", "ND2", "ND3", "COX1"][: len(catalog)],
        mt_primary_gene_class=["protein_coding"] * len(catalog),
        mt_overlapping_genes=["ND1", "ND2", "ND3", "COX1"][: len(catalog)],
    )
    sample_events = sample_events.assign(
        main_analysis_exclusion_reason=["", "", "", "", "", ""],
        mt_breakpoint_span_bp=[np.nan, np.nan, np.nan, np.nan, 100.0, 100.0],
        nuclear_breakpoint_span_bp=[np.nan, np.nan, np.nan, np.nan, 11.0, 11.0],
    )
    retained_samples = pd.DataFrame({"sample_id": ["S1", "S2", "S3", "S4", "S5"]})
    results = compute_main_results(catalog, sample_map, sample_events, retained_samples)
    freq = results["frequency_class_summary"]
    carrier_ratio = results["carrier_ratio_by_frequency_class"]
    _assert((freq["carrier_fraction"] <= 1).all(), "carrier_fraction must be <= 1")
    _assert((carrier_ratio["sample_ratio"] <= 1).all(), "sample_ratio must be <= 1")

    mt_tables = compute_mt_breakpoint_tables(catalog)
    cumulative = mt_tables["mt_breakpoint_cumulative"]
    overall = cumulative[["mt_position", "overall_cumulative_fraction"]].drop_duplicates()
    _assert(abs(float(overall["overall_cumulative_fraction"].iloc[-1]) - 1.0) < 1e-9, "overall cumulative fraction must end at 1")

    mt_genes = pd.DataFrame(
        [
            {"GENENAME": "ND1", "GENE_CLASS": "protein_coding", "POSSTART": 1, "POSEND": 2000},
            {"GENENAME": "COX1", "GENE_CLASS": "protein_coding", "POSSTART": 2001, "POSEND": 5000},
            {"GENENAME": "ND2", "GENE_CLASS": "protein_coding", "POSSTART": 5001, "POSEND": 9000},
        ]
    )
    annotated_catalog = annotate_catalog(catalog, sample_map, sample_events, mt_genes)
    _assert(
        float(annotated_catalog.loc[annotated_catalog["nuclear_cluster_range_start"] == 7000.0, "median_mt_length_for_main_bp"].iloc[0]) == 100.0,
        "annotation stage should preserve integer representative length",
    )


def case_breakpoint_source_region_matching() -> None:
    clusters = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "chr": "chr1",
                "cluster_id_raw": "C1",
                "start": 100,
                "end": 200,
                "Cluster_No": 3,
                "reads": "[120, 140, 160]",
                "mt_positions": "[500, 500, 500]",
            },
            {
                "sample_id": "S1",
                "chr": "chr1",
                "cluster_id_raw": "C2",
                "start": 110,
                "end": 210,
                "Cluster_No": 3,
                "reads": "[130, 150, 170]",
                "mt_positions": "[500, 500, 500]",
            },
        ]
    )
    breakpoints = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "chr": "chr1",
                "pointGroup": "nu_Tend_Bleft",
                "pos": 125.0,
                "readsCount": 5,
                "source_region_key": "S1_chr1.100.200",
            },
            {
                "sample_id": "S1",
                "chr": "chrM",
                "pointGroup": "mt_Tstart",
                "pos": 500.0,
                "readsCount": 5,
                "source_region_key": "S1_chr1.100.200",
            },
        ]
    )
    length_summary = pd.DataFrame(
        columns=[
            "region_key",
            "sample_id",
            "chr",
            "start",
            "end",
            "mt_source_min",
            "mt_source_max",
            "mt_source_envelope_span_bp",
            "mt_source_span_bp",
            "length_bed_read_count",
            "merged_interval_count",
            "complex_mt_interval",
            "mt_primary_fragment_start",
            "mt_primary_fragment_end",
            "mt_primary_fragment_span_bp",
            "mt_merged_intervals",
        ]
    )
    sample_events = build_sample_event_table(clusters, breakpoints, retained_samples={"S1"}, length_summary=length_summary)
    row_c1 = sample_events.loc[sample_events["cluster_id_raw"] == "C1"].iloc[0]
    row_c2 = sample_events.loc[sample_events["cluster_id_raw"] == "C2"].iloc[0]
    _assert(row_c1["breakpoint_match_mode"] == "source_region_key", "exact region key should be preferred when available")
    _assert(pd.notna(row_c1["nuclear_left_breakpoint_pos"]), "matching cluster should recover nuclear breakpoint")
    _assert(pd.isna(row_c2["nuclear_left_breakpoint_pos"]), "neighboring cluster must not borrow breakpoint from another region")
    _assert(pd.isna(row_c2["mt_breakpoint_start"]), "neighboring cluster must not borrow mt breakpoint from another region")


def main() -> int:
    case_discontinuous_mt_intervals()
    case_finalize_sample_events()
    case_distinct_interval_clustering_and_stats()
    case_breakpoint_source_region_matching()
    print("All regression cases passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
