#!/usr/bin/env python3
"""统计汇总函数。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .parse_utils import classify_frequency


def _normalize_bool_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("false").str.lower().isin({"true", "1", "yes"})


def _build_exclusion_summary(sample_events: pd.DataFrame) -> pd.DataFrame:
    excluded = sample_events[~_normalize_bool_series(sample_events["eligible_for_main_analysis"])].copy()
    excluded = excluded[excluded["main_analysis_exclusion_reason"].fillna("").astype(str).ne("")]
    if excluded.empty:
        return pd.DataFrame(columns=["main_analysis_exclusion_reason", "event_n", "sample_n"])
    return (
        excluded.groupby("main_analysis_exclusion_reason", as_index=False)
        .agg(event_n=("sample_event_id", "nunique"), sample_n=("sample_id", "nunique"))
        .sort_values(["event_n", "sample_n", "main_analysis_exclusion_reason"], ascending=[False, False, True])
        .reset_index(drop=True)
    )


def compute_main_results(
    catalog: pd.DataFrame,
    sample_map: pd.DataFrame,
    sample_events: pd.DataFrame,
    retained_samples: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    catalog = catalog.copy()
    sample_map = sample_map.copy()
    sample_events = sample_events.copy()
    retained_samples = retained_samples.copy()
    catalog["distinct_numt_id"] = catalog["distinct_numt_id"].astype(str)
    sample_map["distinct_numt_id"] = sample_map["distinct_numt_id"].astype(str)
    sample_map["sample_id"] = sample_map["sample_id"].astype(str)
    retained_samples["sample_id"] = retained_samples["sample_id"].astype(str)
    sample_events["sample_id"] = sample_events["sample_id"].astype(str)
    sample_events["sample_event_id"] = sample_events["sample_event_id"].astype(str)
    retained_ids = retained_samples["sample_id"].astype(str).tolist()

    sample_counts = (
        sample_map[["sample_id", "distinct_numt_id"]]
        .drop_duplicates()
        .groupby("sample_id", as_index=False)
        .agg(distinct_numt_count=("distinct_numt_id", "nunique"))
    )
    sample_counts = pd.DataFrame({"sample_id": retained_ids}).merge(sample_counts, on="sample_id", how="left").fillna(
        {"distinct_numt_count": 0}
    )
    sample_counts["distinct_numt_count"] = sample_counts["distinct_numt_count"].astype(int)

    distinct_lengths = catalog[["distinct_numt_id", "median_mt_length_for_main_bp"]].copy()
    distinct_lengths["length_level"] = "distinct"
    distinct_lengths = distinct_lengths.rename(columns={"median_mt_length_for_main_bp": "length_bp", "distinct_numt_id": "entity_id"})
    event_lengths = sample_events[["sample_event_id", "mt_length_for_main_bp"]].copy()
    event_lengths["length_level"] = "sample_event"
    event_lengths = event_lengths.rename(columns={"sample_event_id": "entity_id", "mt_length_for_main_bp": "length_bp"})
    length_frames = [frame for frame in [distinct_lengths, event_lengths] if not frame.empty]
    if length_frames:
        length_distribution = pd.concat(length_frames, ignore_index=True).dropna(subset=["length_bp"])
    else:
        length_distribution = pd.DataFrame(columns=["entity_id", "length_bp", "length_level"])

    carrier_ratio = (
        sample_map.merge(catalog[["distinct_numt_id", "frequency_class"]], on="distinct_numt_id", how="left")
        .drop_duplicates(["sample_id", "frequency_class"])
        .groupby("frequency_class", as_index=False)
        .agg(carrier_sample_n=("sample_id", "nunique"))
    )
    carrier_ratio["carrier_fraction"] = carrier_ratio["carrier_sample_n"] / max(len(retained_ids), 1)
    carrier_ratio["sample_count"] = carrier_ratio["carrier_sample_n"]
    carrier_ratio["sample_ratio"] = carrier_ratio["carrier_fraction"]

    carrier_sum = (
        catalog.groupby("frequency_class", as_index=False)
        .agg(
            distinct_numt_n=("distinct_numt_id", "nunique"),
            carrier_count_sum=("carrier_count", "sum"),
        )
        .copy()
    )
    frequency_summary = carrier_sum.merge(
        carrier_ratio[["frequency_class", "carrier_sample_n", "carrier_fraction"]],
        on="frequency_class",
        how="left",
    ).fillna({"carrier_sample_n": 0, "carrier_fraction": 0.0})
    frequency_summary["carrier_sample_n"] = frequency_summary["carrier_sample_n"].astype(int)
    frequency_summary["distinct_fraction"] = frequency_summary["distinct_numt_n"] / max(len(catalog), 1)
    frequency_summary["mean_numt_count_per_sample"] = frequency_summary["carrier_count_sum"] / max(len(retained_ids), 1)

    spectrum = (
        catalog.groupby("carrier_count", as_index=False)
        .agg(distinct_count=("distinct_numt_id", "nunique"))
        .sort_values("carrier_count")
    )
    spectrum["carrier_frequency"] = spectrum["carrier_count"] / max(len(retained_ids), 1)

    exclusion_summary = _build_exclusion_summary(sample_events)

    return {
        "sample_numt_counts": sample_counts,
        "distinct_length_distribution": length_distribution,
        "frequency_class_summary": frequency_summary[
            [
                "frequency_class",
                "distinct_numt_n",
                "distinct_fraction",
                "carrier_sample_n",
                "carrier_fraction",
                "carrier_count_sum",
                "mean_numt_count_per_sample",
            ]
        ],
        "carrier_ratio_by_frequency_class": carrier_ratio[
            ["frequency_class", "carrier_sample_n", "carrier_fraction", "sample_count", "sample_ratio"]
        ],
        "distinct_frequency_spectrum": spectrum,
        "main_analysis_exclusion_summary": exclusion_summary,
    }


def compute_mt_breakpoint_tables(catalog: pd.DataFrame) -> dict[str, pd.DataFrame]:
    catalog = catalog.copy()
    if catalog.empty:
        empty_cumulative = pd.DataFrame(
            columns=[
                "mt_position",
                "overall_count",
                "overall_cumulative_count",
                "overall_cumulative_fraction",
                "frequency_class",
                "count",
                "cumulative_count",
                "cumulative_fraction",
            ]
        )
        return {
            "mt_gene_breakpoint_counts": pd.DataFrame(columns=["mt_primary_gene", "distinct_numt_n"]),
            "mt_gene_breakpoint_counts.by_freq": pd.DataFrame(columns=["mt_primary_gene", "frequency_class", "distinct_numt_n"]),
            "mt_breakpoint_cumulative": empty_cumulative,
        }

    gene_counts = (
        catalog.groupby("mt_primary_gene", as_index=False)
        .agg(distinct_numt_n=("distinct_numt_id", "nunique"))
        .sort_values("distinct_numt_n", ascending=False)
    )
    by_freq = (
        catalog.groupby(["mt_primary_gene", "frequency_class"], as_index=False)
        .agg(distinct_numt_n=("distinct_numt_id", "nunique"))
        .sort_values(["mt_primary_gene", "frequency_class"])
    )

    positions = pd.DataFrame({"mt_position": np.arange(1, 16569 + 1, dtype=int)})
    breakpoint_df = catalog[["distinct_numt_id", "representative_mt_breakpoint", "frequency_class"]].copy()
    breakpoint_df = breakpoint_df.dropna(subset=["representative_mt_breakpoint"]).copy()
    if breakpoint_df.empty:
        positions["overall_count"] = 0
        positions["overall_cumulative_count"] = 0
        positions["overall_cumulative_fraction"] = 0.0
        empty_by_class = pd.DataFrame(columns=["mt_position", "frequency_class", "count", "cumulative_count", "cumulative_fraction"])
        cumulative = positions.merge(empty_by_class, on="mt_position", how="left")
        return {
            "mt_gene_breakpoint_counts": gene_counts,
            "mt_gene_breakpoint_counts.by_freq": by_freq,
            "mt_breakpoint_cumulative": cumulative,
        }

    breakpoint_df["mt_position_rounded"] = breakpoint_df["representative_mt_breakpoint"].round().astype(int).clip(1, 16569)
    overall_counts = breakpoint_df.groupby("mt_position_rounded").size().rename("overall_count").reset_index()
    positions = positions.merge(overall_counts, left_on="mt_position", right_on="mt_position_rounded", how="left").drop(
        columns=["mt_position_rounded"]
    )
    positions["overall_count"] = positions["overall_count"].fillna(0).astype(int)
    positions["overall_cumulative_count"] = positions["overall_count"].cumsum()
    overall_total = int(breakpoint_df["distinct_numt_id"].nunique())
    positions["overall_cumulative_fraction"] = positions["overall_cumulative_count"] / max(overall_total, 1)

    class_tables = []
    for freq_class, sub in breakpoint_df.groupby("frequency_class", sort=False):
        counts = sub.groupby("mt_position_rounded").size().rename("count").reset_index()
        merged = pd.DataFrame({"mt_position": np.arange(1, 16569 + 1, dtype=int)}).merge(
            counts,
            left_on="mt_position",
            right_on="mt_position_rounded",
            how="left",
        )
        merged["count"] = merged["count"].fillna(0).astype(int)
        merged["cumulative_count"] = merged["count"].cumsum()
        class_total = int(sub["distinct_numt_id"].nunique())
        merged["cumulative_fraction"] = merged["cumulative_count"] / max(class_total, 1)
        merged["frequency_class"] = freq_class
        class_tables.append(merged[["mt_position", "frequency_class", "count", "cumulative_count", "cumulative_fraction"]])
    by_class = pd.concat(class_tables, ignore_index=True) if class_tables else pd.DataFrame(
        columns=["mt_position", "frequency_class", "count", "cumulative_count", "cumulative_fraction"]
    )
    cumulative = positions.merge(by_class, on="mt_position", how="left")
    return {
        "mt_gene_breakpoint_counts": gene_counts,
        "mt_gene_breakpoint_counts.by_freq": by_freq,
        "mt_breakpoint_cumulative": cumulative,
    }


def recompute_catalog_for_subset(
    catalog: pd.DataFrame,
    sample_map: pd.DataFrame,
    subgroup_samples: list[str],
    subgroup_sample_n: int,
    prevalent_min: float,
    common_min: float,
    rare_min: float,
) -> pd.DataFrame:
    catalog = catalog.copy()
    sample_map = sample_map.copy()
    if catalog.empty:
        return pd.DataFrame(columns=list(catalog.columns))
    catalog["distinct_numt_id"] = catalog["distinct_numt_id"].astype(str)
    sample_map["distinct_numt_id"] = sample_map["distinct_numt_id"].astype(str)
    sample_map["sample_id"] = sample_map["sample_id"].astype(str)
    subgroup_samples = [str(sample_id) for sample_id in subgroup_samples]
    subset_map = sample_map[sample_map["sample_id"].isin(subgroup_samples)].copy()
    if subset_map.empty:
        return pd.DataFrame(columns=list(catalog.columns))
    counts = subset_map.groupby("distinct_numt_id", as_index=False).agg(carrier_count=("sample_id", "nunique"))
    merged = catalog.merge(counts, on="distinct_numt_id", how="inner", suffixes=("", "_sub"))
    merged["carrier_count"] = merged["carrier_count_sub"]
    merged = merged.drop(columns=["carrier_count_sub"])
    merged["carrier_frequency"] = merged["carrier_count"] / max(subgroup_sample_n, 1)
    merged["frequency_class"] = [
        classify_frequency(count, freq, prevalent_min, common_min, rare_min)
        for count, freq in zip(merged["carrier_count"], merged["carrier_frequency"])
    ]
    return merged
