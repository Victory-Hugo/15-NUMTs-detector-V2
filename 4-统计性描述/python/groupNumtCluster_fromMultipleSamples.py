#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""对多样本 NUMT breakpoint 输入做跨样本聚类。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Sequence, Tuple

import pandas as pd

LOG = logging.getLogger(__name__)

INPUT_COLUMNS = ["sampleID", "Cluster_No", "disFile", "splitFile", "wgsBAM", "chr", "start", "end"]
DETAIL_COLUMNS = ["GroupID", "Index", "POS", "CHR", "mergedClusterID"]
SUMMARY_COLUMNS = ["mergedClusterID", "GroupID", "CHR", "PositionCount", "SampleCount", "min", "max"]


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="对多样本 breakpoint 输入做跨样本聚类。")
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--detail-out")
    parser.add_argument("--summary-out")
    parser.add_argument("--max-gap-bp", type=int, default=1000)
    parser.add_argument("--min-sample-support", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    return parser


def resolve_output_paths(input_file: str, detail_out: str | None, summary_out: str | None) -> Tuple[Path, Path]:
    input_path = Path(input_file)
    default_prefix = Path(str(input_path))
    detail_path = Path(detail_out) if detail_out else Path(str(default_prefix) + ".allCluster.tsv")
    summary_path = Path(summary_out) if summary_out else Path(str(default_prefix) + ".allCluster.sum.tsv")
    return detail_path, summary_path


def read_input_table(input_file: str) -> pd.DataFrame:
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file}")
    if input_path.stat().st_size == 0:
        return pd.DataFrame(columns=INPUT_COLUMNS)

    df = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        names=INPUT_COLUMNS,
        dtype={"sampleID": str, "chr": str},
        compression="infer",
        engine="c",
    )
    if df.empty:
        return df

    df = df.drop_duplicates(["sampleID", "chr", "start"]).copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["sampleID", "chr", "start", "end"]).copy()
    if df.empty:
        return df

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["position"] = ((df["start"] + df["end"]) / 2).astype(int)
    df = df.sort_values(["chr", "position"]).copy()
    df = df.drop_duplicates(["sampleID", "chr", "position"]).copy()
    return df


def cluster_positions(positions: Sequence[int], max_gap_bp: int) -> List[List[int]]:
    if not positions:
        return []

    sorted_positions = sorted(int(item) for item in positions)
    clusters: List[List[int]] = [[sorted_positions[0]]]
    for position in sorted_positions[1:]:
        if position - clusters[-1][-1] <= max_gap_bp:
            clusters[-1].append(position)
        else:
            clusters.append([position])
    return clusters


def build_cluster_tables(input_df: pd.DataFrame, max_gap_bp: int, threads: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if input_df.empty:
        return pd.DataFrame(columns=DETAIL_COLUMNS), pd.DataFrame(columns=SUMMARY_COLUMNS)

    detail_rows: List[dict[str, object]] = []
    summary_rows: List[dict[str, object]] = []
    grouped_df = input_df.groupby("chr", sort=False)
    for chrom, chrom_df in grouped_df:
        positions = sorted(chrom_df["position"].astype(int).unique().tolist())
        sub_clusters = cluster_positions(positions, max_gap_bp=max_gap_bp)
        for group_id, sub_cluster in enumerate(sub_clusters):
            member_rows = chrom_df.loc[chrom_df["position"].astype(int).isin(sub_cluster)].copy()
            sample_count = int(member_rows["sampleID"].astype(str).nunique())
            merged_cluster_id = f"{group_id}_{chrom}"
            for index, position in enumerate(sub_cluster):
                detail_rows.append(
                    {
                        "GroupID": group_id,
                        "Index": index,
                        "POS": int(position),
                        "CHR": str(chrom),
                        "mergedClusterID": merged_cluster_id,
                    }
                )
            summary_rows.append(
                {
                    "mergedClusterID": merged_cluster_id,
                    "GroupID": group_id,
                    "CHR": str(chrom),
                    "PositionCount": len(sub_cluster),
                    "SampleCount": sample_count,
                    "min": int(min(sub_cluster)),
                    "max": int(max(sub_cluster)),
                }
            )

    detail_df = pd.DataFrame(detail_rows, columns=DETAIL_COLUMNS)
    summary_df = pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
    return detail_df, summary_df


def filter_cluster_tables(
    detail_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    min_sample_support: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if min_sample_support < 1:
        raise ValueError("--min-sample-support 必须 >= 1")
    if summary_df.empty:
        return pd.DataFrame(columns=DETAIL_COLUMNS), pd.DataFrame(columns=SUMMARY_COLUMNS)

    filtered_summary = summary_df.loc[summary_df["SampleCount"].astype(int) >= min_sample_support].copy()
    if filtered_summary.empty:
        return pd.DataFrame(columns=DETAIL_COLUMNS), pd.DataFrame(columns=SUMMARY_COLUMNS)

    valid_group_ids = set(filtered_summary["GroupID"].astype(int).tolist())
    filtered_detail = detail_df.loc[detail_df["GroupID"].astype(int).isin(valid_group_ids)].copy()
    filtered_detail = filtered_detail[DETAIL_COLUMNS]
    filtered_summary = filtered_summary[SUMMARY_COLUMNS]
    return filtered_detail, filtered_summary


def write_cluster_tables(detail_df: pd.DataFrame, summary_df: pd.DataFrame, detail_out: str, summary_out: str) -> None:
    detail_path = Path(detail_out)
    summary_path = Path(summary_out)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    detail_df.to_csv(detail_path, sep="\t", index=False)
    summary_df.to_csv(summary_path, sep="\t", index=False)


def run(
    input_file: str,
    detail_out: str | None = None,
    summary_out: str | None = None,
    max_gap_bp: int = 1000,
    min_sample_support: int = 1,
    threads: int = 1,
) -> int:
    resolved_detail_out, resolved_summary_out = resolve_output_paths(input_file, detail_out, summary_out)
    input_df = read_input_table(input_file)
    detail_df, summary_df = build_cluster_tables(input_df=input_df, max_gap_bp=max_gap_bp, threads=threads)
    filtered_detail_df, filtered_summary_df = filter_cluster_tables(
        detail_df=detail_df,
        summary_df=summary_df,
        min_sample_support=min_sample_support,
    )
    write_cluster_tables(
        detail_df=filtered_detail_df,
        summary_df=filtered_summary_df,
        detail_out=str(resolved_detail_out),
        summary_out=str(resolved_summary_out),
    )
    LOG.info(
        "NUMT 聚类完成: total_clusters=%s retained_clusters=%s min_sample_support=%s",
        len(summary_df),
        len(filtered_summary_df),
        min_sample_support,
    )
    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        input_file=args.input_file,
        detail_out=args.detail_out,
        summary_out=args.summary_out,
        max_gap_bp=args.max_gap_bp,
        min_sample_support=args.min_sample_support,
        threads=args.threads,
    )


if __name__ == "__main__":
    raise SystemExit(main())
