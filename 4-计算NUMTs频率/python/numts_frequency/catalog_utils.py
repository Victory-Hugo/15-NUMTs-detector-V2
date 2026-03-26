#!/usr/bin/env python3
"""样本事件表与 distinct catalog 构建。"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .parse_utils import (
    LEFT_GROUP,
    MT_END_GROUP,
    MT_START_GROUP,
    RIGHT_GROUP,
    build_region_key,
    classify_frequency,
    parse_int_list,
    safe_midpoint,
    weighted_median,
)

LOG = logging.getLogger(__name__)


def build_breakpoint_index(breakpoints: pd.DataFrame) -> dict:
    """按样本与染色体建立二分检索索引。"""
    index: dict = {}
    for key, sub in breakpoints.groupby(["sample_id", "chr"], sort=False):
        ordered = sub.sort_values("pos").reset_index(drop=True)
        index[key] = {
            "positions": ordered["pos"].to_numpy(dtype=float),
            "weights": ordered["readsCount"].to_numpy(dtype=float),
            "point_groups": ordered["pointGroup"].astype(str).to_numpy(),
        }
    return index


def query_interval(index: dict, sample_id: str, chrom: str, start: float, end: float) -> dict:
    data = index.get((sample_id, chrom))
    if data is None:
        return {"positions": np.array([]), "weights": np.array([]), "point_groups": np.array([])}
    left = int(np.searchsorted(data["positions"], float(start), side="left"))
    right = int(np.searchsorted(data["positions"], float(end), side="right"))
    return {
        "positions": data["positions"][left:right],
        "weights": data["weights"][left:right],
        "point_groups": data["point_groups"][left:right],
    }


def _weighted_point(interval: dict, point_group: str) -> float | None:
    mask = interval["point_groups"] == point_group
    if mask.sum() == 0:
        return None
    return weighted_median(interval["positions"][mask], interval["weights"][mask])


def _normalize_bool_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("false").str.lower().isin({"true", "1", "yes"})


def _join_reason(reasons: list[str]) -> str:
    return "|".join(reasons)


def _single_point_or_interval(start: float | None, end: float | None) -> tuple[float | None, float | None]:
    start_is_na = pd.isna(start)
    end_is_na = pd.isna(end)
    if start_is_na and end_is_na:
        return np.nan, np.nan
    if start_is_na:
        return float(end), float(end)
    if end_is_na:
        return float(start), float(start)
    return float(min(start, end)), float(max(start, end))


def _interval_single_linkage_labels(df: pd.DataFrame, start_col: str, end_col: str, gap: int) -> list[int]:
    labels: list[int] = []
    current_label = -1
    current_max_end: float | None = None
    for row in df.itertuples(index=False):
        start = float(getattr(row, start_col))
        end = float(getattr(row, end_col))
        if current_max_end is None or start > current_max_end + gap:
            current_label += 1
            current_max_end = end
        else:
            current_max_end = max(current_max_end, end)
        labels.append(current_label)
    return labels


def _mode_smallest(values: pd.Series) -> float | None:
    clean = values.dropna()
    if clean.empty:
        return None
    rounded = clean.round().astype(int)
    counts = rounded.value_counts()
    max_count = counts.max()
    return float(min(counts[counts == max_count].index.tolist()))


def finalize_sample_event_table(sample_events: pd.DataFrame, length_summary: pd.DataFrame) -> pd.DataFrame:
    summary_cols = [
        "region_key",
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
    drop_cols = [col for col in summary_cols if col in sample_events.columns]
    prepared = sample_events.drop(columns=drop_cols, errors="ignore").copy()
    merged = prepared.merge(
        length_summary[summary_cols],
        left_on="length_bed_region_key",
        right_on="region_key",
        how="left",
    ).copy()
    merged["nuclear_breakpoint_span_bp"] = merged["nuclear_right_breakpoint_pos"] - merged["nuclear_left_breakpoint_pos"] + 1
    merged.loc[
        merged["nuclear_left_breakpoint_pos"].isna() | merged["nuclear_right_breakpoint_pos"].isna(),
        "nuclear_breakpoint_span_bp",
    ] = np.nan
    merged["mt_breakpoint_span_bp"] = merged["mt_breakpoint_end"] - merged["mt_breakpoint_start"] + 1
    merged.loc[
        merged["mt_breakpoint_start"].isna() | merged["mt_breakpoint_end"].isna(),
        "mt_breakpoint_span_bp",
    ] = np.nan
    merged["mt_source_span_bp"] = merged["mt_source_span_bp"].fillna(merged["mt_core_span_bp"])
    merged["mt_source_min"] = merged["mt_source_min"].fillna(merged["mt_core_start"])
    merged["mt_source_max"] = merged["mt_source_max"].fillna(merged["mt_core_end"])
    merged["mt_source_envelope_span_bp"] = merged["mt_source_envelope_span_bp"].fillna(merged["mt_core_span_bp"])
    merged["length_bed_read_count"] = merged["length_bed_read_count"].fillna(0).astype(int)
    merged["merged_interval_count"] = merged["merged_interval_count"].fillna(1).astype(int)
    merged["complex_mt_interval"] = _normalize_bool_series(merged["complex_mt_interval"])
    merged["mt_primary_fragment_start"] = merged["mt_primary_fragment_start"].fillna(merged["mt_source_min"])
    merged["mt_primary_fragment_end"] = merged["mt_primary_fragment_end"].fillna(merged["mt_source_max"])
    merged["mt_primary_fragment_span_bp"] = merged["mt_primary_fragment_span_bp"].fillna(merged["mt_source_span_bp"])
    merged["mt_merged_intervals"] = merged["mt_merged_intervals"].fillna(
        merged.apply(
            lambda row: (
                f"{int(row.mt_source_min)}-{int(row.mt_source_max)}"
                if pd.notna(row.mt_source_min) and pd.notna(row.mt_source_max)
                else ""
            ),
            axis=1,
        )
    )
    merged["length_source"] = np.where(merged["region_key"].notna(), "length_bed", "cluster_mt_positions")
    merged["mt_length_for_main_bp"] = merged["mt_breakpoint_span_bp"]
    fallback_mask = merged["mt_length_for_main_bp"].isna() & (~merged["complex_mt_interval"])
    merged.loc[fallback_mask, "mt_length_for_main_bp"] = merged.loc[fallback_mask, "mt_source_span_bp"]

    reasons: list[str] = []
    reason_rows = []
    for row in merged.itertuples(index=False):
        reasons = []
        has_nuclear = pd.notna(row.nuclear_left_breakpoint_pos) or pd.notna(row.nuclear_right_breakpoint_pos)
        has_mt = pd.notna(row.mt_breakpoint_start) or pd.notna(row.mt_breakpoint_end)
        if bool(row.complex_mt_interval):
            reasons.append("complex_mt_interval")
        if not has_nuclear:
            reasons.append("missing_nuclear_breakpoint")
        if not has_mt:
            reasons.append("missing_mt_breakpoint")
        reason_rows.append(_join_reason(reasons))
    merged["main_analysis_exclusion_reason"] = reason_rows
    merged["eligible_for_main_analysis"] = merged["main_analysis_exclusion_reason"].eq("")
    merged = merged.drop(columns=["region_key"], errors="ignore")
    return merged


def build_sample_event_table(
    clusters: pd.DataFrame,
    breakpoints: pd.DataFrame,
    retained_samples: set[str],
    length_summary: pd.DataFrame,
) -> pd.DataFrame:
    """从 cluster 主表和 breakpoint 证据构建样本事件表。"""
    index = build_breakpoint_index(breakpoints[breakpoints["sample_id"].isin(retained_samples)].copy())
    length_keys = set(length_summary["region_key"].astype(str))
    length_lookup = (
        length_summary.set_index("region_key")[["mt_source_min", "mt_source_max"]].to_dict("index")
        if not length_summary.empty
        else {}
    )
    records = []

    for row in clusters[clusters["sample_id"].isin(retained_samples)].itertuples(index=False):
        reads = parse_int_list(row.reads)
        mt_positions = parse_int_list(row.mt_positions)
        region_key = build_region_key(row.sample_id, row.chr, row.start, row.end)
        mt_lookup = length_lookup.get(region_key, {})
        mt_query_start = mt_lookup.get("mt_source_min")
        mt_query_end = mt_lookup.get("mt_source_max")
        if pd.isna(mt_query_start) or pd.isna(mt_query_end):
            if mt_positions:
                mt_query_start = min(mt_positions)
                mt_query_end = max(mt_positions)
            else:
                mt_query_start = np.nan
                mt_query_end = np.nan
        nuclear_interval = query_interval(index, row.sample_id, row.chr, row.start, row.end)
        if pd.notna(mt_query_start) and pd.notna(mt_query_end):
            mt_interval = query_interval(index, row.sample_id, "chrM", mt_query_start, mt_query_end)
        else:
            mt_interval = {"positions": np.array([]), "weights": np.array([]), "point_groups": np.array([])}
        left_pos = _weighted_point(nuclear_interval, LEFT_GROUP)
        right_pos = _weighted_point(nuclear_interval, RIGHT_GROUP)
        mt_start = _weighted_point(mt_interval, MT_START_GROUP)
        mt_end = _weighted_point(mt_interval, MT_END_GROUP)
        cluster_width_bp = int(row.end - row.start + 1)
        records.append(
            {
                "sample_event_id": f"{row.sample_id}|{row.cluster_id_raw}",
                "sample_id": row.sample_id,
                "cluster_id_raw": row.cluster_id_raw,
                "chr": row.chr,
                "cluster_start": int(row.start),
                "cluster_end": int(row.end),
                "cluster_width_bp": cluster_width_bp,
                "cluster_support_n": int(row.Cluster_No) if not pd.isna(row.Cluster_No) else len(reads),
                "cluster_reads_raw": row.reads,
                "cluster_mt_positions_raw": row.mt_positions,
                "nuclear_core_start": min(reads) if reads else np.nan,
                "nuclear_core_end": max(reads) if reads else np.nan,
                "nuclear_core_midpoint": safe_midpoint(min(reads), max(reads)) if reads else np.nan,
                "mt_core_start": min(mt_positions) if mt_positions else np.nan,
                "mt_core_end": max(mt_positions) if mt_positions else np.nan,
                "mt_core_midpoint": safe_midpoint(min(mt_positions), max(mt_positions)) if mt_positions else np.nan,
                "mt_core_span_bp": (max(mt_positions) - min(mt_positions) + 1) if mt_positions else np.nan,
                "nuclear_left_breakpoint_pos": left_pos,
                "nuclear_right_breakpoint_pos": right_pos,
                "nuclear_anchor_midpoint": safe_midpoint(left_pos, right_pos),
                "mt_breakpoint_start": mt_start,
                "mt_breakpoint_end": mt_end,
                "supporting_breakpoint_rows": int(nuclear_interval["positions"].shape[0] + mt_interval["positions"].shape[0]),
                "length_bed_region_key": region_key,
                "length_bed_match_status": "matched" if region_key in length_keys else "not_matched",
            }
        )
    result = pd.DataFrame.from_records(records)
    return finalize_sample_event_table(result, length_summary)


def build_distinct_catalog(
    sample_events: pd.DataFrame,
    retained_sample_count: int,
    nuclear_merge_gap_bp: int,
    mt_merge_gap_bp: int,
    distinct_id_prefix: str,
    prevalent_min: float,
    common_min: float,
    rare_min: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """按核断点和 mt 断点双层聚类形成 distinct NUMT。"""
    catalog_columns = [
        "distinct_numt_id",
        "chr",
        "nuclear_cluster_range_start",
        "nuclear_cluster_range_end",
        "distinct_nuclear_start",
        "distinct_nuclear_end",
        "distinct_nuclear_midpoint",
        "distinct_mt_start",
        "distinct_mt_end",
        "distinct_mt_midpoint",
        "merged_nuclear_breakpoint_start",
        "merged_nuclear_breakpoint_end",
        "representative_nuclear_left_breakpoint",
        "representative_nuclear_right_breakpoint",
        "representative_mt_breakpoint",
        "carrier_count",
        "carrier_frequency",
        "sample_event_count",
        "median_cluster_width_bp",
        "median_nuclear_breakpoint_span_bp",
        "median_mt_breakpoint_span_bp",
        "median_mt_source_union_span_bp",
        "median_mt_source_span_bp",
        "median_mt_length_for_main_bp",
        "frequency_class",
    ]
    events = sample_events.copy()
    if "eligible_for_main_analysis" not in events.columns:
        raise ValueError("sample_events missing eligible_for_main_analysis")
    events["eligible_for_main_analysis"] = _normalize_bool_series(events["eligible_for_main_analysis"])
    events = events[events["eligible_for_main_analysis"]].copy()
    if events.empty:
        return pd.DataFrame(columns=catalog_columns), pd.DataFrame(
            columns=["distinct_numt_id", "sample_id", "sample_event_id", "chr"]
        )

    nuclear_bounds = [
        _single_point_or_interval(start, end)
        for start, end in zip(events["nuclear_left_breakpoint_pos"], events["nuclear_right_breakpoint_pos"])
    ]
    mt_bounds = [
        _single_point_or_interval(start, end)
        for start, end in zip(events["mt_breakpoint_start"], events["mt_breakpoint_end"])
    ]
    events["nuclear_merge_start"] = [item[0] for item in nuclear_bounds]
    events["nuclear_merge_end"] = [item[1] for item in nuclear_bounds]
    events["mt_merge_start"] = [item[0] for item in mt_bounds]
    events["mt_merge_end"] = [item[1] for item in mt_bounds]

    combined_parts = []
    for chrom, chrom_df in events.groupby("chr", sort=True):
        chrom_df = chrom_df.sort_values(["nuclear_merge_start", "nuclear_merge_end", "sample_event_id"]).reset_index(drop=True)
        chrom_df["nuclear_cluster_label"] = _interval_single_linkage_labels(
            chrom_df,
            "nuclear_merge_start",
            "nuclear_merge_end",
            nuclear_merge_gap_bp,
        )
        for nuc_label, nuc_df in chrom_df.groupby("nuclear_cluster_label", sort=True):
            nuc_df = nuc_df.sort_values(["mt_merge_start", "mt_merge_end", "sample_event_id"]).reset_index(drop=True)
            nuc_df["mt_cluster_label"] = _interval_single_linkage_labels(
                nuc_df,
                "mt_merge_start",
                "mt_merge_end",
                mt_merge_gap_bp,
            )
            nuc_df["combined_cluster_key"] = [f"{chrom}|{nuc_label}|{sub}" for sub in nuc_df["mt_cluster_label"]]
            combined_parts.append(nuc_df)
    events = pd.concat(combined_parts, ignore_index=True)

    grouped = []
    for key, sub in events.groupby("combined_cluster_key", sort=False):
        grouped.append(
            (
                sub["chr"].iloc[0],
                float(sub["nuclear_merge_start"].min()),
                float(sub["mt_merge_start"].min()),
                key,
                sub.copy(),
            )
        )
    grouped.sort(key=lambda item: (item[0], item[1], item[2]))

    distinct_rows = []
    sample_map_rows = []
    for idx, (_, _, _, _, sub) in enumerate(grouped, start=1):
        distinct_id = f"{distinct_id_prefix}{idx:07d}"
        carrier_count = int(sub["sample_id"].nunique())
        carrier_frequency = float(carrier_count) / float(max(retained_sample_count, 1))
        nuclear_span = sub["nuclear_right_breakpoint_pos"] - sub["nuclear_left_breakpoint_pos"] + 1
        nuclear_span = nuclear_span.where(
            sub["nuclear_right_breakpoint_pos"].notna() & sub["nuclear_left_breakpoint_pos"].notna(),
            np.nan,
        )
        mt_span = sub["mt_breakpoint_end"] - sub["mt_breakpoint_start"] + 1
        mt_span = mt_span.where(sub["mt_breakpoint_end"].notna() & sub["mt_breakpoint_start"].notna(), np.nan)
        frequency_class = classify_frequency(
            carrier_count=carrier_count,
            carrier_frequency=carrier_frequency,
            prevalent_min=prevalent_min,
            common_min=common_min,
            rare_min=rare_min,
        )
        representative_mt_values = pd.concat(
            [sub["mt_breakpoint_start"], sub["mt_breakpoint_end"]],
            ignore_index=True,
        )
        nuclear_cluster_range_start = float(sub["nuclear_merge_start"].min())
        nuclear_cluster_range_end = float(sub["nuclear_merge_end"].max())
        representative_nuclear_left = _mode_smallest(sub["nuclear_left_breakpoint_pos"])
        representative_nuclear_right = _mode_smallest(sub["nuclear_right_breakpoint_pos"])
        if representative_nuclear_left is None and representative_nuclear_right is None:
            merged_nuclear_breakpoint_start = nuclear_cluster_range_start
            merged_nuclear_breakpoint_end = nuclear_cluster_range_end
        elif representative_nuclear_left is None:
            merged_nuclear_breakpoint_start = float(representative_nuclear_right)
            merged_nuclear_breakpoint_end = float(representative_nuclear_right)
        elif representative_nuclear_right is None:
            merged_nuclear_breakpoint_start = float(representative_nuclear_left)
            merged_nuclear_breakpoint_end = float(representative_nuclear_left)
        else:
            merged_nuclear_breakpoint_start = float(min(representative_nuclear_left, representative_nuclear_right))
            merged_nuclear_breakpoint_end = float(max(representative_nuclear_left, representative_nuclear_right))
        distinct_nuclear_start = merged_nuclear_breakpoint_start
        distinct_nuclear_end = merged_nuclear_breakpoint_end
        distinct_mt_start = float(sub["mt_merge_start"].min())
        distinct_mt_end = float(sub["mt_merge_end"].max())
        distinct_rows.append(
            {
                "distinct_numt_id": distinct_id,
                "chr": sub["chr"].iloc[0],
                "nuclear_cluster_range_start": nuclear_cluster_range_start,
                "nuclear_cluster_range_end": nuclear_cluster_range_end,
                "distinct_nuclear_start": distinct_nuclear_start,
                "distinct_nuclear_end": distinct_nuclear_end,
                "distinct_nuclear_midpoint": safe_midpoint(distinct_nuclear_start, distinct_nuclear_end),
                "distinct_mt_start": distinct_mt_start,
                "distinct_mt_end": distinct_mt_end,
                "distinct_mt_midpoint": safe_midpoint(distinct_mt_start, distinct_mt_end),
                "merged_nuclear_breakpoint_start": merged_nuclear_breakpoint_start,
                "merged_nuclear_breakpoint_end": merged_nuclear_breakpoint_end,
                "representative_nuclear_left_breakpoint": representative_nuclear_left,
                "representative_nuclear_right_breakpoint": representative_nuclear_right,
                "representative_mt_breakpoint": _mode_smallest(representative_mt_values),
                "carrier_count": carrier_count,
                "carrier_frequency": carrier_frequency,
                "sample_event_count": int(len(sub)),
                "median_cluster_width_bp": float(sub["cluster_width_bp"].median()),
                "median_nuclear_breakpoint_span_bp": float(nuclear_span.median()) if nuclear_span.notna().any() else np.nan,
                "median_mt_breakpoint_span_bp": float(mt_span.median()) if mt_span.notna().any() else np.nan,
                "median_mt_source_union_span_bp": (
                    float(sub["mt_source_span_bp"].median()) if sub["mt_source_span_bp"].notna().any() else np.nan
                ),
                "median_mt_source_span_bp": (
                    float(sub["mt_source_span_bp"].median()) if sub["mt_source_span_bp"].notna().any() else np.nan
                ),
                "median_mt_length_for_main_bp": (
                    float(sub["mt_length_for_main_bp"].median()) if sub["mt_length_for_main_bp"].notna().any() else np.nan
                ),
                "frequency_class": frequency_class,
            }
        )
        for event_row in sub.itertuples(index=False):
            sample_map_rows.append(
                {
                    "distinct_numt_id": distinct_id,
                    "sample_id": event_row.sample_id,
                    "sample_event_id": event_row.sample_event_id,
                    "chr": event_row.chr,
                }
            )

    catalog = pd.DataFrame.from_records(distinct_rows)
    sample_map = pd.DataFrame.from_records(sample_map_rows).drop_duplicates(
        subset=["distinct_numt_id", "sample_id", "sample_event_id"]
    )
    return catalog, sample_map
