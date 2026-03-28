#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""按染色体对样本级事件进行 1 kb 聚类并计算人群频率。"""

from __future__ import annotations

import argparse
import csv
import logging
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

LOG = logging.getLogger(__name__)


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="对事件表做染色体内 1 kb 聚类并输出频率表。")
    parser.add_argument("--events-tsv", required=True)
    parser.add_argument("--cluster-out", required=True)
    parser.add_argument("--class-summary-out", required=True)
    parser.add_argument("--support-summary-out", required=True)
    parser.add_argument("--top-out", required=True)
    parser.add_argument("--cluster-gap-bp", type=int, required=True)
    parser.add_argument("--denominator", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    return parser


def chrom_sort_key(chrom: str) -> Tuple[int, str]:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return (int(suffix), chrom)
        if suffix == "X":
            return (23, chrom)
        if suffix == "Y":
            return (24, chrom)
    return (10_000, chrom)


def classify_frequency(sample_count: int, denominator: int) -> str:
    frequency = sample_count / denominator
    if frequency >= 0.05:
        return "common"
    if frequency >= 0.01:
        return "low-frequency"
    if frequency >= 0.001:
        return "rare"
    return "ultra-rare"


def cluster_one_chromosome(task: Tuple[str, List[Tuple[int, str]], int, int]) -> List[Dict[str, object]]:
    chrom, rows, cluster_gap_bp, denominator = task
    rows = sorted(rows, key=lambda item: item[0])
    if not rows:
        return []

    result = []
    cluster_index = 1

    current_points = [rows[0][0]]
    current_unique_points = {rows[0][0]}
    current_samples = {rows[0][1]}
    previous_midpoint = rows[0][0]

    for midpoint, sample_id in rows[1:]:
        if midpoint - previous_midpoint <= cluster_gap_bp:
            current_points.append(midpoint)
            current_unique_points.add(midpoint)
            current_samples.add(sample_id)
        else:
            sample_count = len(current_samples)
            result.append(
                {
                    "merged_cluster_id": f"{cluster_index}_{chrom}",
                    "chr": chrom,
                    "cluster_min_midpoint": min(current_points),
                    "cluster_max_midpoint": max(current_points),
                    "cluster_midpoint_mean": sum(current_points) / len(current_points),
                    "POS": len(current_unique_points),
                    "sample_count": sample_count,
                    "frequency": sample_count / denominator,
                    "frequency_class": classify_frequency(sample_count, denominator),
                    "is_singleton": sample_count == 1,
                    "is_doubleton": sample_count == 2,
                    "is_tripleton": sample_count == 3,
                }
            )
            cluster_index += 1
            current_points = [midpoint]
            current_unique_points = {midpoint}
            current_samples = {sample_id}
        previous_midpoint = midpoint

    sample_count = len(current_samples)
    result.append(
        {
            "merged_cluster_id": f"{cluster_index}_{chrom}",
            "chr": chrom,
            "cluster_min_midpoint": min(current_points),
            "cluster_max_midpoint": max(current_points),
            "cluster_midpoint_mean": sum(current_points) / len(current_points),
            "POS": len(current_unique_points),
            "sample_count": sample_count,
            "frequency": sample_count / denominator,
            "frequency_class": classify_frequency(sample_count, denominator),
            "is_singleton": sample_count == 1,
            "is_doubleton": sample_count == 2,
            "is_tripleton": sample_count == 3,
        }
    )
    return result


def run(
    events_tsv: str,
    cluster_out: str,
    class_summary_out: str,
    support_summary_out: str,
    top_out: str,
    cluster_gap_bp: int,
    denominator: int,
    threads: int,
) -> int:
    event_df = pd.read_csv(events_tsv, sep="\t", dtype={"sampleID": str, "region_chr": str})
    event_df["midpoint"] = event_df["midpoint"].astype(int)

    events_by_chr: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    for row in event_df.itertuples(index=False):
        events_by_chr[row.region_chr].append((int(row.midpoint), str(row.sampleID)))

    tasks = [(chrom, rows, cluster_gap_bp, denominator) for chrom, rows in events_by_chr.items()]
    cluster_rows: List[Dict[str, object]] = []
    if max(1, threads) > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(threads, len(tasks))) as pool:
            for rows in pool.map(cluster_one_chromosome, tasks):
                cluster_rows.extend(rows)
    else:
        for task in tasks:
            cluster_rows.extend(cluster_one_chromosome(task))

    cluster_rows.sort(key=lambda row: (chrom_sort_key(str(row["chr"])), int(row["cluster_min_midpoint"])))

    cluster_df = pd.DataFrame(cluster_rows)
    cluster_df.to_csv(cluster_out, sep="\t", index=False)

    class_order = ["common", "low-frequency", "rare", "ultra-rare"]
    class_counts = Counter(cluster_df["frequency_class"].tolist())
    with Path(class_summary_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "cluster_count"])
        for category in class_order:
            writer.writerow([category, class_counts.get(category, 0)])

    sample_counts = cluster_df["sample_count"].astype(int)
    support_rows = [
        ("all_clusters", int(len(cluster_df))),
        ("singleton_clusters", int((sample_counts == 1).sum())),
        ("doubleton_clusters", int((sample_counts == 2).sum())),
        ("tripleton_clusters", int((sample_counts == 3).sum())),
        ("clusters_ge_2_samples", int((sample_counts >= 2).sum())),
        ("clusters_ge_3_samples", int((sample_counts >= 3).sum())),
        ("clusters_ge_4_samples", int((sample_counts >= 4).sum())),
        ("clusters_ge_5_samples", int((sample_counts >= 5).sum())),
        ("clusters_ge_10_samples", int((sample_counts >= 10).sum())),
        ("clusters_ge_20_samples", int((sample_counts >= 20).sum())),
        ("clusters_ge_50_samples", int((sample_counts >= 50).sum())),
        ("clusters_ge_100_samples", int((sample_counts >= 100).sum())),
    ]
    with Path(support_summary_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(support_rows)

    top_df = cluster_df.sort_values(["sample_count", "frequency", "chr", "cluster_min_midpoint"], ascending=[False, False, True, True]).head(20)
    top_df.to_csv(top_out, sep="\t", index=False)

    LOG.info(
        "频率聚类完成，总 cluster 数 %s；singleton %s；>=2 samples %s；>=3 samples %s；>=4 samples %s",
        len(cluster_df),
        int((sample_counts == 1).sum()),
        int((sample_counts >= 2).sum()),
        int((sample_counts >= 3).sum()),
        int((sample_counts >= 4).sum()),
    )
    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        events_tsv=args.events_tsv,
        cluster_out=args.cluster_out,
        class_summary_out=args.class_summary_out,
        support_summary_out=args.support_summary_out,
        top_out=args.top_out,
        cluster_gap_bp=args.cluster_gap_bp,
        denominator=args.denominator,
        threads=args.threads,
    )


if __name__ == "__main__":
    raise SystemExit(main())
