#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""过滤质控通过样本并输出过滤后的大表（v2：新增 cluster input QC 过滤）。"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, Set

import pandas as pd

LOG = logging.getLogger(__name__)


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="过滤 PASS 样本并输出过滤后的 TSV/GZ 文件（v2）。")
    parser.add_argument("--meta-tsv", required=True)
    parser.add_argument("--breakpoint-tsv-gz", required=True)
    parser.add_argument("--merge-bed-tsv-gz", required=True)
    # v2 新增：cluster input 过滤
    parser.add_argument("--cluster-input-gz", required=True,
                        help="原始 cluster input TSV.GZ（无表头，第一列为 sampleID）")
    parser.add_argument("--cluster-out-gz", required=True,
                        help="输出：QC 过滤后的 cluster input TSV.GZ")
    parser.add_argument("--meta-id-col", required=True)
    parser.add_argument("--meta-qc-col", required=True)
    parser.add_argument("--meta-qc-pass-value", required=True)
    parser.add_argument("--pass-sample-tsv", required=True)
    parser.add_argument("--breakpoint-out-gz", required=True)
    parser.add_argument("--merge-bed-out-gz", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--threads", type=int, default=2)
    return parser


def write_pass_samples(pass_samples: Iterable[str], output_path: str) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sampleID"])
        for sample_id in sorted(pass_samples):
            writer.writerow([sample_id])


def parse_merge_sample_id(record_key: str) -> str:
    region_key = record_key.split("|", 1)[0]
    return region_key.split("_", 1)[0]


def filter_breakpoint_file(input_path: str, output_path: str, pass_samples: Set[str]) -> Dict[str, int]:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    total_rows = 0
    kept_rows = 0
    seen_samples: Set[str] = set()
    kept_samples: Set[str] = set()

    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as src, gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"空文件: {input_path}")
        dst.write(header)
        for line in src:
            total_rows += 1
            sample_id = line.split("\t", 1)[0]
            seen_samples.add(sample_id)
            if sample_id in pass_samples:
                dst.write(line)
                kept_rows += 1
                kept_samples.add(sample_id)
            if total_rows % 1_000_000 == 0:
                LOG.info("breakpoint 已处理 %s 行", total_rows)

    return {
        "total_rows": total_rows,
        "kept_rows": kept_rows,
        "total_samples": len(seen_samples),
        "kept_samples": len(kept_samples),
    }


def filter_merge_bed_file(input_path: str, output_path: str, pass_samples: Set[str]) -> Dict[str, int]:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    total_rows = 0
    kept_rows = 0
    seen_samples: Set[str] = set()
    kept_samples: Set[str] = set()

    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as src, gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as dst:
        for line in src:
            total_rows += 1
            record_key = line.split("\t", 1)[0]
            sample_id = parse_merge_sample_id(record_key)
            seen_samples.add(sample_id)
            if sample_id in pass_samples:
                dst.write(line)
                kept_rows += 1
                kept_samples.add(sample_id)
            if total_rows % 2_000_000 == 0:
                LOG.info("merge_bed 已处理 %s 行", total_rows)

    return {
        "total_rows": total_rows,
        "kept_rows": kept_rows,
        "total_samples": len(seen_samples),
        "kept_samples": len(kept_samples),
    }


def filter_cluster_input_file(input_path: str, output_path: str, pass_samples: Set[str]) -> Dict[str, int]:
    """过滤 cluster input 文件（无表头，第一列为 sampleID）。

    cluster input 格式（无表头，制表符分隔）：
      sampleID  Cluster_No  disFile  splitFile  wgsBAM  chr  start  end
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    total_rows = 0
    kept_rows = 0
    seen_samples: Set[str] = set()
    kept_samples: Set[str] = set()

    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as src, gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as dst:
        for line in src:
            total_rows += 1
            sample_id = line.split("\t", 1)[0]
            seen_samples.add(sample_id)
            if sample_id in pass_samples:
                dst.write(line)
                kept_rows += 1
                kept_samples.add(sample_id)
            if total_rows % 1_000_000 == 0:
                LOG.info("cluster_input 已处理 %s 行", total_rows)

    return {
        "total_rows": total_rows,
        "kept_rows": kept_rows,
        "total_samples": len(seen_samples),
        "kept_samples": len(kept_samples),
    }


def write_summary(
    output_path: str,
    meta_total: int,
    meta_pass: int,
    breakpoint_summary: Dict[str, int],
    merge_summary: Dict[str, int],
    cluster_summary: Dict[str, int],
) -> None:
    rows = [
        ("meta_total_samples", meta_total),
        ("meta_pass_samples", meta_pass),
        ("meta_fail_samples", meta_total - meta_pass),
        ("breakpoint_total_rows", breakpoint_summary["total_rows"]),
        ("breakpoint_kept_rows", breakpoint_summary["kept_rows"]),
        ("breakpoint_total_samples", breakpoint_summary["total_samples"]),
        ("breakpoint_pass_observed_samples", breakpoint_summary["kept_samples"]),
        ("breakpoint_dropped_samples", breakpoint_summary["total_samples"] - breakpoint_summary["kept_samples"]),
        ("merge_bed_total_rows", merge_summary["total_rows"]),
        ("merge_bed_kept_rows", merge_summary["kept_rows"]),
        ("merge_bed_total_samples", merge_summary["total_samples"]),
        ("merge_bed_pass_observed_samples", merge_summary["kept_samples"]),
        ("merge_bed_dropped_samples", merge_summary["total_samples"] - merge_summary["kept_samples"]),
        # v2 新增：cluster input 统计
        ("cluster_input_total_rows", cluster_summary["total_rows"]),
        ("cluster_input_kept_rows", cluster_summary["kept_rows"]),
        ("cluster_input_total_samples", cluster_summary["total_samples"]),
        ("cluster_input_pass_observed_samples", cluster_summary["kept_samples"]),
        ("cluster_input_dropped_samples", cluster_summary["total_samples"] - cluster_summary["kept_samples"]),
    ]
    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def run(
    meta_tsv: str,
    breakpoint_tsv_gz: str,
    merge_bed_tsv_gz: str,
    cluster_input_gz: str,
    cluster_out_gz: str,
    meta_id_col: str,
    meta_qc_col: str,
    meta_qc_pass_value: str,
    pass_sample_tsv: str,
    breakpoint_out_gz: str,
    merge_bed_out_gz: str,
    summary_out: str,
    threads: int = 2,
) -> int:
    meta_df = pd.read_csv(meta_tsv, sep="\t", dtype=str, usecols=[meta_id_col, meta_qc_col])
    meta_df[meta_id_col] = meta_df[meta_id_col].fillna("")
    meta_df[meta_qc_col] = meta_df[meta_qc_col].fillna("")

    pass_samples = set(meta_df.loc[meta_df[meta_qc_col] == meta_qc_pass_value, meta_id_col])
    meta_total = int(meta_df[meta_id_col].nunique())
    meta_pass = len(pass_samples)

    LOG.info("meta 总样本数: %s", meta_total)
    LOG.info("meta PASS 样本数: %s", meta_pass)

    write_pass_samples(pass_samples, pass_sample_tsv)

    # 三个过滤任务并行：breakpoint、merge_bed、cluster_input
    max_workers = max(1, min(threads, 3))
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        future_bp = pool.submit(filter_breakpoint_file, breakpoint_tsv_gz, breakpoint_out_gz, pass_samples)
        future_merge = pool.submit(filter_merge_bed_file, merge_bed_tsv_gz, merge_bed_out_gz, pass_samples)
        future_cluster = pool.submit(filter_cluster_input_file, cluster_input_gz, cluster_out_gz, pass_samples)
        breakpoint_summary = future_bp.result()
        merge_summary = future_merge.result()
        cluster_summary = future_cluster.result()

    write_summary(
        summary_out,
        meta_total,
        meta_pass,
        breakpoint_summary,
        merge_summary,
        cluster_summary,
    )
    LOG.info("过滤完成，summary 已写出: %s", summary_out)
    return 0


def main(argv: list[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        meta_tsv=args.meta_tsv,
        breakpoint_tsv_gz=args.breakpoint_tsv_gz,
        merge_bed_tsv_gz=args.merge_bed_tsv_gz,
        cluster_input_gz=args.cluster_input_gz,
        cluster_out_gz=args.cluster_out_gz,
        meta_id_col=args.meta_id_col,
        meta_qc_col=args.meta_qc_col,
        meta_qc_pass_value=args.meta_qc_pass_value,
        pass_sample_tsv=args.pass_sample_tsv,
        breakpoint_out_gz=args.breakpoint_out_gz,
        merge_bed_out_gz=args.merge_bed_out_gz,
        summary_out=args.summary_out,
        threads=args.threads,
    )


if __name__ == "__main__":
    raise SystemExit(main())
