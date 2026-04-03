#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""基于 merge_bed 的 mtDNA 全覆盖跨度（max(mt_end)-min(mt_start) per locus）
修正 NUMTs 长度，输出含真实 mtDNA 插入长度的频率表及散点图。"""

from __future__ import annotations

import argparse
import logging
import math
import re
from pathlib import Path
from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)

# 用于解析 read_name 末尾 locus 坐标的正则：sampleID_chrN.start.end
LOCUS_RE = re.compile(r"_(chr[^.]+)\.(\d+)\.(\d+)$")

FREQ_CLASS_ORDER: List[str] = ["common", "low-frequency", "rare", "ultra-rare"]
FREQ_CLASS_COLORS = {
    "common": "#d62728",
    "low-frequency": "#ff7f0e",
    "rare": "#2ca02c",
    "ultra-rare": "#1f77b4",
}
FREQ_BOUNDARIES = {
    "5%": 0.05,
    "1%": 0.01,
    "0.1%": 0.001,
}


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="基于 mtDNA 全覆盖跨度修正 NUMTs 长度，输出修正频率表和散点图。"
    )
    parser.add_argument(
        "--merge-bed-gz",
        required=True,
        help="过滤后的 merge_bed 文件（3-merge_bed.pass.tsv.gz）",
    )
    parser.add_argument(
        "--cluster-detail-tsv",
        required=True,
        help="min-support-1 的 allCluster.tsv（含 GroupID/POS/CHR/mergedClusterID）",
    )
    parser.add_argument(
        "--freq-cluster-tsv",
        required=True,
        help="4-numt-frequency-by-cluster.tsv",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="修正后的频率+长度 TSV 输出路径",
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help="散点图 PDF 输出路径",
    )
    parser.add_argument(
        "--output-pdf-revised",
        default=None,
        help="过滤 D-loop circular artifact 后的散点图 PDF（可选）",
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=500_000,
        help="merge_bed 分块行数（默认 500000）",
    )
    return parser


def build_locus_coverage_span(merge_bed_gz: str, chunk_rows: int) -> pd.DataFrame:
    """按 locus（sampleID_chr.start.end）计算 mtDNA 全覆盖跨度。

    对每个 locus 的所有 reads 取 max(mt_end) - min(mt_start)，
    再对同一 (chr, nDNA midpoint) 的多个样本取中位数，作为该基因组位置的代表长度。

    这与 chosen_length（单条 read 最大 span）的关键区别在于：
    当多条 read 分布在 mtDNA 的不同区段时，覆盖跨度能正确还原完整插入长度。
    """
    LOG.info("解析 merge_bed 构建 mtDNA 覆盖跨度查找表: %s", merge_bed_gz)

    chunk_agg_list: List[pd.DataFrame] = []
    total_rows = 0

    for chunk_idx, chunk in enumerate(
        pd.read_csv(
            merge_bed_gz,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            dtype={0: str, 1: np.int32, 2: np.int32},
            compression="infer",
            chunksize=chunk_rows,
            engine="c",
        )
    ):
        chunk.columns = ["read_name", "mt_start", "mt_end"]
        # locus_id = read_name 去掉末尾的 |readID 部分
        chunk["locus_id"] = chunk["read_name"].str.rsplit("|", n=1).str[0]

        # 每块内按 locus_id 做 min/max 聚合
        agg = chunk.groupby("locus_id", sort=False).agg(
            chunk_min=("mt_start", "min"),
            chunk_max=("mt_end", "max"),
        )
        chunk_agg_list.append(agg)
        total_rows += len(chunk)

        if (chunk_idx + 1) % 20 == 0:
            LOG.info("已处理 %d 个分块（%d 行）", chunk_idx + 1, total_rows)

    LOG.info("merge_bed 读取完成，共 %d 行，开始合并聚合结果...", total_rows)

    # 跨块合并：对同一 locus_id 取全局 min/max
    combined = pd.concat(chunk_agg_list)
    final = (
        combined.groupby("locus_id", sort=False)
        .agg(min_mt_start=("chunk_min", "min"), max_mt_end=("chunk_max", "max"))
        .reset_index()
    )
    final["coverage_span"] = final["max_mt_end"] - final["min_mt_start"]

    # 从 locus_id 解析 (chr, nDNA midpoint)
    locus_parsed = final["locus_id"].str.extract(LOCUS_RE.pattern)
    locus_parsed.columns = ["chr", "nDNA_start", "nDNA_end"]
    valid = locus_parsed[["nDNA_start", "nDNA_end"]].notna().all(axis=1)

    locus_df = pd.DataFrame(
        {
            "chr": locus_parsed.loc[valid, "chr"],
            "midpoint": (
                (
                    locus_parsed.loc[valid, "nDNA_start"].astype(np.int64)
                    + locus_parsed.loc[valid, "nDNA_end"].astype(np.int64)
                )
                // 2
            ).astype(np.int64),
            "coverage_span": final.loc[valid, "coverage_span"].values,
            "min_mt_start": final.loc[valid, "min_mt_start"].values,
            "max_mt_end": final.loc[valid, "max_mt_end"].values,
        }
    )

    # 同一 nDNA 位置跨样本取中位数，得到每个 locus 的代表长度及 D-loop 检测所需的端点信息
    locus_agg = (
        locus_df.groupby(["chr", "midpoint"], sort=False)
        .agg(
            mt_span_bp=("coverage_span", "median"),
            median_min_mt_start=("min_mt_start", "median"),
            median_max_mt_end=("max_mt_end", "median"),
        )
        .reset_index()
    )
    locus_agg["mt_span_bp"] = locus_agg["mt_span_bp"].round().astype(int)
    locus_span = locus_agg
    LOG.info(
        "查找表构建完成，共 %d 个唯一 (chr, midpoint) locus，最大跨度 %d bp",
        len(locus_span),
        int(locus_span["mt_span_bp"].max()),
    )
    return locus_span


def compute_cluster_mt_span(
    cluster_detail_tsv: str, locus_span: pd.DataFrame
) -> pd.DataFrame:
    """对每个 cluster 取各 position 的 mtDNA 覆盖跨度中位数作为代表长度。"""
    LOG.info("加载 cluster detail: %s", cluster_detail_tsv)
    detail = pd.read_csv(
        cluster_detail_tsv,
        sep="\t",
        dtype={"CHR": str, "POS": np.int64, "mergedClusterID": str},
    )
    merged = detail.merge(
        locus_span,
        left_on=["CHR", "POS"],
        right_on=["chr", "midpoint"],
        how="left",
    )
    cluster_span = (
        merged.groupby("mergedClusterID", sort=False)
        .agg(
            numt_length_bp=("mt_span_bp", "median"),
            cluster_min_mt_start=("median_min_mt_start", "median"),
            cluster_max_mt_end=("median_max_mt_end", "median"),
        )
        .reset_index()
        .rename(columns={"mergedClusterID": "merged_cluster_id"})
    )
    cluster_span["numt_length_bp"] = cluster_span["numt_length_bp"].round()
    # D-loop 环状 artifact：reads 同时聚集在 mtDNA 线性表示的两端
    cluster_span["is_dloop_artifact"] = (
        (cluster_span["cluster_min_mt_start"] < 200)
        & (cluster_span["cluster_max_mt_end"] > 16350)
    )
    na_count = int(cluster_span["numt_length_bp"].isna().sum())
    if na_count:
        LOG.warning(
            "有 %d 个 cluster 未匹配到 mtDNA 长度（输出中标为 NA，散点图中跳过）",
            na_count,
        )
    dloop_count = int(cluster_span["is_dloop_artifact"].sum())
    LOG.info(
        "cluster 长度计算完成，共 %d 个 cluster，最大长度 %s bp，其中 D-loop 环状 artifact %d 个",
        len(cluster_span),
        cluster_span["numt_length_bp"].max(),
        dloop_count,
    )
    return cluster_span


def build_result_table(
    freq_tsv: str, cluster_span: pd.DataFrame
) -> pd.DataFrame:
    """合并频率表与 mtDNA 长度，添加 neg_log2_frequency 列。"""
    freq_df = pd.read_csv(freq_tsv, sep="\t")
    result = freq_df.merge(cluster_span, on="merged_cluster_id", how="left")
    result["neg_log2_frequency"] = -np.log2(result["frequency"])
    return result


def plot_scatter(result_df: pd.DataFrame, output_pdf: str) -> None:
    """绘制 -log2(frequency) vs NUMTs 长度散点图，颜色按 frequency_class。"""
    LOG.info("绘制散点图 → %s", output_pdf)
    plot_df = result_df.dropna(subset=["numt_length_bp", "neg_log2_frequency"]).copy()
    plot_df["numt_length_bp"] = plot_df["numt_length_bp"].astype(float)

    fig, ax = plt.subplots(figsize=(8, 6))

    # 由稀到常绘制，确保 common 点在最上层
    for cls in reversed(FREQ_CLASS_ORDER):
        sub = plot_df[plot_df["frequency_class"] == cls]
        if sub.empty:
            continue
        ax.scatter(
            sub["neg_log2_frequency"],
            sub["numt_length_bp"],
            c=FREQ_CLASS_COLORS[cls],
            label=cls,
            s=6,
            alpha=0.45,
            linewidths=0,
            rasterized=True,
        )

    # 频率类别分界参考竖线
    y_max = plot_df["numt_length_bp"].max()
    for label, freq_val in FREQ_BOUNDARIES.items():
        xval = -math.log2(freq_val)
        ax.axvline(x=xval, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.text(
            xval + 0.1,
            y_max,
            label,
            fontsize=7,
            color="gray",
            va="top",
            rotation=90,
        )

    ax.set_yscale("log")
    ax.set_xlabel(r"$-\log_2$(Frequency of NUMTs)", fontsize=12)
    ax.set_ylabel("Size of NUMTs (bp)", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)

    handles = [
        mpatches.Patch(color=FREQ_CLASS_COLORS[cls], label=cls)
        for cls in FREQ_CLASS_ORDER
        if cls in plot_df["frequency_class"].values
    ]
    ax.legend(
        handles=handles,
        title="Frequency class",
        fontsize=9,
        title_fontsize=9,
        loc="upper left",
        framealpha=0.85,
    )

    plt.tight_layout()
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    LOG.info("散点图已保存: %s", output_pdf)


def run(
    merge_bed_gz: str,
    cluster_detail_tsv: str,
    freq_cluster_tsv: str,
    output_tsv: str,
    output_pdf: str,
    chunk_rows: int,
    output_pdf_revised: Optional[str] = None,
) -> int:
    locus_span = build_locus_coverage_span(merge_bed_gz, chunk_rows)
    cluster_span = compute_cluster_mt_span(cluster_detail_tsv, locus_span)
    result_df = build_result_table(freq_cluster_tsv, cluster_span)

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_tsv, sep="\t", index=False)
    LOG.info(
        "修正后的频率+长度表已保存: %s，共 %d 行", output_tsv, len(result_df)
    )

    plot_scatter(result_df, output_pdf)

    if output_pdf_revised:
        filtered_df = result_df[~result_df["is_dloop_artifact"].fillna(False)]
        n_removed = int(len(result_df) - len(filtered_df))
        LOG.info("过滤 D-loop circular artifact：移除 %d 个 cluster，剩余 %d 个", n_removed, len(filtered_df))
        plot_scatter(filtered_df, output_pdf_revised)

    return 0


def main(argv: Optional[List[str]] = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        merge_bed_gz=args.merge_bed_gz,
        cluster_detail_tsv=args.cluster_detail_tsv,
        freq_cluster_tsv=args.freq_cluster_tsv,
        output_tsv=args.output_tsv,
        output_pdf=args.output_pdf,
        chunk_rows=args.chunk_rows,
        output_pdf_revised=args.output_pdf_revised,
    )


if __name__ == "__main__":
    raise SystemExit(main())
