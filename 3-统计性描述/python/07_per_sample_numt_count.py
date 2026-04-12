#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""统计每个个体携带的非参考 NUMT 数量，并生成分布图。

数据来源: breakpoint TSV.gz
算法:
  1. 读取 breakpoint TSV.gz，提取所有唯一 sampleID 作为样本宇宙；
  2. 按 (sampleID, source_region_key) 去重，每个唯一 insertion region 算 1 个 NUMT；
  3. 统计每人 NUMT 数量，对无检出的样本补 0；
  4. 输出每样本计数表、分布汇总表，以及双面板分布图。
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd

LOG = logging.getLogger(__name__)

FILL_COLOR = "#0072b2"
EDGE_COLOR = "#30203e"
GRID_COLOR = "#d9d9d9"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="统计每个个体携带的非参考 NUMT 数量，并生成分布图。"
    )
    parser.add_argument("--breakpoint-tsv-gz", required=True,
                        help="breakpoint TSV.gz，含 sampleID 和 source_region_key 列")
    parser.add_argument("--output-per-sample-tsv", required=True,
                        help="每样本 NUMT 数量输出表（sampleID, numt_count）")
    parser.add_argument("--output-distribution-tsv", required=True,
                        help="分布汇总表（numt_count, sample_count）")
    parser.add_argument("--out-prefix", required=True,
                        help="图形输出前缀（不含扩展名），自动生成 .pdf / .svg / .png")
    return parser


# ---------------------------------------------------------------------------
# 核心计算
# ---------------------------------------------------------------------------

def compute_per_sample_counts(
    breakpoint_tsv_gz: str,
) -> pd.DataFrame:
    """返回 DataFrame，列为 [sampleID, numt_count]，涵盖 breakpoint 文件中所有样本（含 0）。"""

    LOG.info("读取 breakpoint 文件: %s", breakpoint_tsv_gz)
    bp_df = pd.read_csv(
        breakpoint_tsv_gz, sep="\t",
        dtype={"sampleID": str, "source_region_key": str},
        usecols=["sampleID", "source_region_key"],
    )
    LOG.info("总 breakpoint 行数: %d", len(bp_df))

    # 从 breakpoint 文件中提取全部唯一样本 ID 作为样本宇宙
    all_ids: set[str] = set(bp_df["sampleID"].tolist())
    LOG.info("breakpoint 文件中唯一样本数: %d", len(all_ids))

    # 每个 (sampleID, source_region_key) 唯一组合 = 1 个 NUMT
    bp_df = bp_df.drop_duplicates(["sampleID", "source_region_key"])

    per_sample: pd.DataFrame = (
        bp_df.groupby("sampleID", as_index=False)
        .size()
        .rename(columns={"size": "numt_count"})
    )

    # 对出现在 breakpoint 文件但去重后无检出的样本补 0
    detected_ids = set(per_sample["sampleID"].tolist())
    zero_ids = sorted(all_ids - detected_ids)
    if zero_ids:
        zero_df = pd.DataFrame({"sampleID": zero_ids, "numt_count": 0})
        per_sample = pd.concat([per_sample, zero_df], ignore_index=True)

    per_sample["numt_count"] = per_sample["numt_count"].astype(int)
    per_sample = per_sample.sort_values(["numt_count", "sampleID"]).reset_index(drop=True)

    n_with = int((per_sample["numt_count"] > 0).sum())
    max_val = int(per_sample["numt_count"].max()) if not per_sample.empty else 0
    median_val = float(per_sample["numt_count"].median()) if not per_sample.empty else 0.0
    LOG.info(
        "每样本 NUMT 计数完成: 样本总数=%d  有NUMT=%d  max=%d  median=%.1f",
        len(per_sample), n_with, max_val, median_val,
    )
    return per_sample


def build_distribution_df(per_sample_df: pd.DataFrame) -> pd.DataFrame:
    counts = per_sample_df["numt_count"].value_counts().sort_index()
    return pd.DataFrame({
        "numt_count": counts.index.astype(int),
        "sample_count": counts.values.astype(int),
    })


# ---------------------------------------------------------------------------
# 绘图
# ---------------------------------------------------------------------------

def configure_matplotlib() -> None:
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False


def _add_stats_box(ax: plt.Axes, lines: list[str]) -> None:
    ax.text(
        0.98, 0.98, "\n".join(lines),
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9.5, color=EDGE_COLOR,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f8f8f8", "edgecolor": "#d9d9d9"},
    )


def _draw_bar_panel(
    ax: plt.Axes,
    dist_df: pd.DataFrame,
    title: str,
    stats_lines: list[str],
) -> None:
    if not dist_df.empty:
        ax.bar(
            dist_df["numt_count"], dist_df["sample_count"],
            color=FILL_COLOR, edgecolor=EDGE_COLOR, linewidth=0.7, width=0.75,
        )
        max_x = int(dist_df["numt_count"].max())
        ax.set_xticks(range(0, max_x + 1))
    ax.set_xlabel("NUMTs per individual")
    ax.set_ylabel("Number of individuals")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)
    _add_stats_box(ax, stats_lines)


def plot_per_sample_distribution(per_sample_df: pd.DataFrame, out_prefix: Path) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor("white")

    # Panel A: 全部 QC-pass 样本（含 0-NUMT）
    dist_all = build_distribution_df(per_sample_df)
    n_total = len(per_sample_df)
    n_zero = int((per_sample_df["numt_count"] == 0).sum())
    _draw_bar_panel(
        axes[0], dist_all,
        title="A. All QC-pass individuals",
        stats_lines=[
            f"n = {n_total:,}",
            f"n(0 NUMT) = {n_zero:,}",
            f"mean = {per_sample_df['numt_count'].mean():.1f}",
            f"median = {float(per_sample_df['numt_count'].median()):.1f}",
            f"max = {int(per_sample_df['numt_count'].max())}",
        ],
    )

    # Panel B: 仅含 ≥1 NUMT 的样本
    pos_df = per_sample_df.loc[per_sample_df["numt_count"] >= 1].copy()
    dist_pos = build_distribution_df(pos_df)
    _draw_bar_panel(
        axes[1], dist_pos,
        title="B. Individuals with \u2265\u20091 NUMT",
        stats_lines=(
            [
                f"n = {len(pos_df):,}",
                f"mean = {pos_df['numt_count'].mean():.1f}",
                f"median = {float(pos_df['numt_count'].median()):.1f}",
                f"max = {int(pos_df['numt_count'].max())}",
            ]
            if not pos_df.empty
            else ["n = 0"]
        ),
    )

    fig.suptitle("Per-individual NUMT Count", fontsize=14, fontweight="bold")
    fig.tight_layout()

    for suffix, extra in [(".pdf", {}), (".svg", {}), (".png", {"dpi": 300})]:
        fig.savefig(str(out_prefix) + suffix, bbox_inches="tight", **extra)
    plt.close(fig)
    LOG.info("分布图已保存: %s.[pdf/svg/png]", out_prefix)


# ---------------------------------------------------------------------------
# run / main
# ---------------------------------------------------------------------------

def run(
    breakpoint_tsv_gz: str,
    output_per_sample_tsv: str,
    output_distribution_tsv: str,
    out_prefix: str,
) -> int:
    per_sample_df = compute_per_sample_counts(
        breakpoint_tsv_gz=breakpoint_tsv_gz,
    )

    out_per_sample = Path(output_per_sample_tsv)
    out_per_sample.parent.mkdir(parents=True, exist_ok=True)
    per_sample_df.to_csv(out_per_sample, sep="\t", index=False)
    LOG.info("每样本计数表已保存: %s", out_per_sample)

    dist_df = build_distribution_df(per_sample_df)
    out_dist = Path(output_distribution_tsv)
    out_dist.parent.mkdir(parents=True, exist_ok=True)
    dist_df.to_csv(out_dist, sep="\t", index=False)
    LOG.info("分布汇总表已保存: %s", out_dist)

    out_prefix_path = Path(out_prefix)
    out_prefix_path.parent.mkdir(parents=True, exist_ok=True)
    plot_per_sample_distribution(per_sample_df, out_prefix_path)

    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        breakpoint_tsv_gz=args.breakpoint_tsv_gz,
        output_per_sample_tsv=args.output_per_sample_tsv,
        output_distribution_tsv=args.output_distribution_tsv,
        out_prefix=args.out_prefix,
    )


if __name__ == "__main__":
    raise SystemExit(main())
