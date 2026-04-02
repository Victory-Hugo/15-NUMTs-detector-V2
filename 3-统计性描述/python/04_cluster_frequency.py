#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""基于多样本 GroupID 聚类结果输出 distinct NUMTs 统计。"""

from __future__ import annotations

import argparse
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List

import pandas as pd

from groupNumtCluster_fromMultipleSamples import (
    build_cluster_tables,
    filter_cluster_tables,
    read_input_table,
    write_cluster_tables,
)

LOG = logging.getLogger(__name__)


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="按 GroupID 聚类结果统计 distinct NUMTs。")
    parser.add_argument("--cluster-input-file", required=True)
    parser.add_argument("--cluster-prefix", required=True)
    parser.add_argument("--cluster-out", required=True)
    parser.add_argument("--class-summary-out", required=True)
    parser.add_argument("--support-summary-out", required=True)
    parser.add_argument("--top-out", required=True)
    parser.add_argument("--cluster-gap-bp", type=int, required=True)
    parser.add_argument("--denominator", type=int, required=True)
    parser.add_argument("--min-supports", required=True)
    parser.add_argument("--primary-min-support", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    return parser


def classify_frequency(sample_count: int, denominator: int) -> str:
    frequency = sample_count / denominator
    if frequency >= 0.05:
        return "common"
    if frequency >= 0.01:
        return "low-frequency"
    if frequency >= 0.001:
        return "rare"
    return "ultra-rare"


def parse_min_supports(raw_value: str) -> List[int]:
    values: List[int] = []
    for token in raw_value.split(","):
        token = token.strip()
        if not token:
            continue
        value = int(token)
        if value < 1:
            raise ValueError("--min-supports 里的值必须 >= 1")
        values.append(value)
    if not values:
        raise ValueError("--min-supports 不能为空")
    return sorted(set(values))


def with_suffix_before_extension(path_str: str, suffix: str) -> str:
    path = Path(path_str)
    return str(path.with_name(f"{path.stem}{suffix}{path.suffix}"))


def build_frequency_cluster_table(summary_df: pd.DataFrame, denominator: int) -> pd.DataFrame:
    if summary_df.empty:
        return pd.DataFrame(
            columns=[
                "merged_cluster_id",
                "GroupID",
                "chr",
                "cluster_min_midpoint",
                "cluster_max_midpoint",
                "cluster_midpoint_mean",
                "POS",
                "sample_count",
                "chosen_region_length",
                "frequency",
                "frequency_class",
                "is_singleton",
                "is_doubleton",
                "is_tripleton",
            ]
        )

    working_df = summary_df.copy()
    if "mean" not in working_df.columns:
        working_df["mean"] = (
            pd.to_numeric(working_df["min"], errors="coerce") + pd.to_numeric(working_df["max"], errors="coerce")
        ) / 2

    freq_df = working_df.rename(
        columns={
            "mergedClusterID": "merged_cluster_id",
            "CHR": "chr",
            "PositionCount": "POS",
            "SampleCount": "sample_count",
            "mean": "cluster_midpoint_mean",
            "min": "cluster_min_midpoint",
            "max": "cluster_max_midpoint",
            "chosen_region_length": "chosen_region_length",
        }
    ).copy()
    freq_df["sample_count"] = freq_df["sample_count"].astype(int)
    freq_df["frequency"] = freq_df["sample_count"] / denominator
    freq_df["frequency_class"] = freq_df["sample_count"].map(lambda value: classify_frequency(int(value), denominator))
    freq_df["is_singleton"] = freq_df["sample_count"] == 1
    freq_df["is_doubleton"] = freq_df["sample_count"] == 2
    freq_df["is_tripleton"] = freq_df["sample_count"] == 3
    return freq_df[
        [
            "merged_cluster_id",
            "GroupID",
            "chr",
            "cluster_min_midpoint",
            "cluster_max_midpoint",
            "cluster_midpoint_mean",
            "POS",
            "sample_count",
            "chosen_region_length",
            "frequency",
            "frequency_class",
            "is_singleton",
            "is_doubleton",
            "is_tripleton",
        ]
    ]


def build_class_summary(freq_df: pd.DataFrame) -> pd.DataFrame:
    class_order = ["common", "low-frequency", "rare", "ultra-rare"]
    counts = freq_df["frequency_class"].value_counts().to_dict() if not freq_df.empty else {}
    return pd.DataFrame(
        [{"category": category, "cluster_count": int(counts.get(category, 0))} for category in class_order]
    )


def build_support_summary_row(min_sample_support: int, freq_df: pd.DataFrame) -> Dict[str, int]:
    sample_counts = freq_df["sample_count"].astype(int) if not freq_df.empty else pd.Series(dtype=int)
    return {
        "min_sample_support": min_sample_support,
        "distinct_numt_count": int(len(freq_df)),
        "singleton_clusters": int((sample_counts == 1).sum()),
        "doubleton_clusters": int((sample_counts == 2).sum()),
        "tripleton_clusters": int((sample_counts == 3).sum()),
        "clusters_ge_2_samples": int((sample_counts >= 2).sum()),
        "clusters_ge_3_samples": int((sample_counts >= 3).sum()),
        "clusters_ge_4_samples": int((sample_counts >= 4).sum()),
        "clusters_ge_5_samples": int((sample_counts >= 5).sum()),
        "clusters_ge_10_samples": int((sample_counts >= 10).sum()),
        "clusters_ge_20_samples": int((sample_counts >= 20).sum()),
        "clusters_ge_50_samples": int((sample_counts >= 50).sum()),
        "clusters_ge_100_samples": int((sample_counts >= 100).sum()),
    }


def write_threshold_outputs(
    min_sample_support: int,
    detail_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    cluster_prefix: str,
    cluster_out: str,
    class_summary_out: str,
    top_out: str,
    denominator: int,
    is_primary: bool,
) -> Dict[str, object]:
    filtered_detail_df, filtered_summary_df = filter_cluster_tables(
        detail_df=detail_df,
        summary_df=summary_df,
        min_sample_support=min_sample_support,
    )

    detail_out = f"{cluster_prefix}.min-support-{min_sample_support}.allCluster.tsv"
    summary_out = f"{cluster_prefix}.min-support-{min_sample_support}.allCluster.sum.tsv"
    write_cluster_tables(filtered_detail_df, filtered_summary_df, detail_out=detail_out, summary_out=summary_out)

    freq_df = build_frequency_cluster_table(filtered_summary_df, denominator=denominator)
    class_df = build_class_summary(freq_df)
    top_df = freq_df.sort_values(
        ["sample_count", "frequency", "chr", "cluster_min_midpoint"],
        ascending=[False, False, True, True],
    ).head(20)

    per_support_cluster_out = with_suffix_before_extension(cluster_out, f".min-support-{min_sample_support}")
    per_support_class_out = with_suffix_before_extension(class_summary_out, f".min-support-{min_sample_support}")
    per_support_top_out = with_suffix_before_extension(top_out, f".min-support-{min_sample_support}")
    Path(per_support_cluster_out).parent.mkdir(parents=True, exist_ok=True)
    freq_df.to_csv(per_support_cluster_out, sep="\t", index=False)
    class_df.to_csv(per_support_class_out, sep="\t", index=False)
    top_df.to_csv(per_support_top_out, sep="\t", index=False)

    if is_primary:
        freq_df.to_csv(cluster_out, sep="\t", index=False)
        class_df.to_csv(class_summary_out, sep="\t", index=False)
        top_df.to_csv(top_out, sep="\t", index=False)

    return {
        "support_row": build_support_summary_row(min_sample_support=min_sample_support, freq_df=freq_df),
        "cluster_count": int(len(freq_df)),
    }


def run(
    cluster_input_file: str,
    cluster_prefix: str,
    cluster_out: str,
    class_summary_out: str,
    support_summary_out: str,
    top_out: str,
    cluster_gap_bp: int,
    denominator: int,
    min_supports: str,
    primary_min_support: int,
    threads: int,
) -> int:
    support_values = parse_min_supports(min_supports)
    if primary_min_support not in support_values:
        raise ValueError("--primary-min-support 必须包含在 --min-supports 内")

    input_df = read_input_table(cluster_input_file)
    observed_sample_count = int(input_df["sampleID"].astype(str).nunique()) if not input_df.empty else 0
    effective_denominator = max(int(denominator), observed_sample_count)
    if effective_denominator != int(denominator):
        LOG.warning(
            "配置 denominator=%s 小于聚类输入观测样本数=%s，自动使用 %s 以避免 frequency > 1",
            denominator,
            observed_sample_count,
            effective_denominator,
        )
    detail_df, summary_df = build_cluster_tables(input_df=input_df, max_gap_bp=cluster_gap_bp, threads=threads)

    results: List[Dict[str, object]] = []
    max_workers = max(1, min(int(threads), len(support_values)))
    if max_workers > 1 and len(support_values) > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(
                    write_threshold_outputs,
                    min_sample_support=value,
                    detail_df=detail_df,
                    summary_df=summary_df,
                    cluster_prefix=cluster_prefix,
                    cluster_out=cluster_out,
                    class_summary_out=class_summary_out,
                    top_out=top_out,
                    denominator=effective_denominator,
                    is_primary=value == primary_min_support,
                )
                for value in support_values
            ]
            for future in futures:
                results.append(future.result())
    else:
        for value in support_values:
            results.append(
                write_threshold_outputs(
                    min_sample_support=value,
                    detail_df=detail_df,
                    summary_df=summary_df,
                    cluster_prefix=cluster_prefix,
                    cluster_out=cluster_out,
                    class_summary_out=class_summary_out,
                    top_out=top_out,
                    denominator=effective_denominator,
                    is_primary=value == primary_min_support,
                )
            )

    support_rows = sorted(
        [item["support_row"] for item in results],
        key=lambda row: int(row["min_sample_support"]),
    )
    pd.DataFrame(support_rows).to_csv(support_summary_out, sep="\t", index=False)

    primary_cluster_count = next(
        int(item["cluster_count"]) for value, item in zip(support_values, results) if value == primary_min_support
    )
    LOG.info(
        "distinct NUMTs 聚类完成: total_group_ids=%s primary_min_support=%s primary_distinct_numts=%s",
        len(summary_df),
        primary_min_support,
        primary_cluster_count,
    )
    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        cluster_input_file=args.cluster_input_file,
        cluster_prefix=args.cluster_prefix,
        cluster_out=args.cluster_out,
        class_summary_out=args.class_summary_out,
        support_summary_out=args.support_summary_out,
        top_out=args.top_out,
        cluster_gap_bp=args.cluster_gap_bp,
        denominator=args.denominator,
        min_supports=args.min_supports,
        primary_min_support=args.primary_min_support,
        threads=args.threads,
    )


if __name__ == "__main__":
    raise SystemExit(main())
