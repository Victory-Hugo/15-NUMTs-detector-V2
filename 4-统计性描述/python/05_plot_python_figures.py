#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成单一多面板统计图。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)

DISCRETE_PALETTE = ["#0072b2", "#56b4e9", "#009e73", "#f0e442", "#e69f00", "#d55e00"]
CONTINUOUS_PALETTE = [
    "#0b0405",
    "#30203e",
    "#3e4d93",
    "#366b9f",
    "#3488a6",
    "#36a4ab",
    "#49c1ad",
    "#60ceac",
    "#84d8b0",
    "#c4e9cf",
    "#def5e5",
]
EDGE_COLOR = "#30203e"
GRID_COLOR = "#d9d9d9"
CLASS_ORDER = ["common", "low-frequency", "rare", "ultra-rare"]
CLASS_COLORS = dict(zip(CLASS_ORDER, DISCRETE_PALETTE[: len(CLASS_ORDER)]))


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def configure_matplotlib() -> None:
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="生成 NUMT 单一多面板统计图。")
    parser.add_argument("--length-tsv", required=True)
    parser.add_argument("--chr-count-tsv", required=True)
    parser.add_argument("--freq-cluster-tsv", required=True)
    parser.add_argument("--freq-class-tsv", required=True)
    parser.add_argument("--support-summary-tsv", required=True)
    parser.add_argument("--out-prefix", required=True)
    return parser


def save_figure(fig: plt.Figure, out_prefix: Path) -> None:
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def choose_histogram_bins(lengths: np.ndarray) -> np.ndarray:
    if lengths.size == 0:
        raise ValueError("长度数组为空")
    if np.all(lengths == lengths[0]):
        return np.array([lengths[0] - 0.5, lengths[0] + 0.5])
    bins = np.histogram_bin_edges(lengths, bins="fd")
    bin_count = max(1, len(bins) - 1)
    if bin_count < 30:
        bins = np.histogram_bin_edges(lengths, bins="rice")
    elif bin_count > 120:
        bins = np.histogram_bin_edges(lengths, bins=120)
    return np.unique(np.round(bins).astype(int))


def chrom_sort_value(chrom: str) -> int:
    chrom = str(chrom)
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return int(suffix)
        if suffix == "X":
            return 23
        if suffix == "Y":
            return 24
    return 10_000


def apply_axis_style(ax: plt.Axes) -> None:
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)


def add_panel_stats_box(ax: plt.Axes, lines: list[str]) -> None:
    ax.text(
        0.02,
        0.98,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        color=EDGE_COLOR,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f8f8f8", "edgecolor": "#d9d9d9"},
    )


def choose_finer_histogram_bins(lengths: np.ndarray, max_bins: int = 180) -> np.ndarray:
    if lengths.size == 0:
        raise ValueError("长度数组为空")
    if np.all(lengths == lengths[0]):
        return np.array([lengths[0] - 0.5, lengths[0] + 0.5])
    bins = np.histogram_bin_edges(lengths, bins="sqrt")
    bin_count = max(1, len(bins) - 1)
    if bin_count < 45:
        bins = np.histogram_bin_edges(lengths, bins="rice")
    elif bin_count > max_bins:
        bins = np.histogram_bin_edges(lengths, bins=max_bins)
    return np.unique(np.round(bins).astype(int))


def plot_length_histogram(ax: plt.Axes, freq_cluster_df: pd.DataFrame) -> None:
    lengths = freq_cluster_df["chosen_region_length"].astype(float).to_numpy()
    bins = choose_histogram_bins(lengths)
    if bins.size < 2:
        bins = np.array([lengths.min() - 0.5, lengths.max() + 0.5])
    ax.hist(lengths, bins=bins, color=CONTINUOUS_PALETTE[6], edgecolor=EDGE_COLOR, linewidth=0.7)
    ax.set_xlabel("Distinct NUMT locus length (bp)")
    ax.set_ylabel("Distinct locus count")
    ax.set_title("A. Distinct locus length distribution", loc="left", fontweight="bold")
    apply_axis_style(ax)
    add_panel_stats_box(
        ax,
        [
            f"n = {len(lengths):,}",
            f"median = {np.median(lengths):.1f} bp",
            f"max = {np.max(lengths):.0f} bp",
        ],
    )


def plot_region_length_histogram(ax: plt.Axes, length_df: pd.DataFrame) -> None:
    lengths = length_df["chosen_length"].astype(float).to_numpy()
    bins = choose_finer_histogram_bins(lengths)
    if bins.size < 2:
        bins = np.array([lengths.min() - 0.5, lengths.max() + 0.5])
    ax.hist(lengths, bins=bins, color=CONTINUOUS_PALETTE[8], edgecolor=EDGE_COLOR, linewidth=0.7)
    ax.set_xlabel("Per-NUMT length (bp)")
    ax.set_ylabel("NUMT count")
    ax.set_title("B. Per-NUMT length distribution", loc="left", fontweight="bold")
    apply_axis_style(ax)
    add_panel_stats_box(
        ax,
        [
            f"n = {len(lengths):,}",
            f"median = {np.median(lengths):.1f} bp",
            f"max = {np.max(lengths):.0f} bp",
        ],
    )


def plot_region_length_histogram_zoom(ax: plt.Axes, length_df: pd.DataFrame, x_max: int = 500) -> None:
    lengths = length_df["chosen_length"].astype(float).to_numpy()
    zoom_lengths = lengths[lengths <= x_max]
    bins = choose_finer_histogram_bins(zoom_lengths, max_bins=120)
    if bins.size < 2:
        bins = np.array([zoom_lengths.min() - 0.5, zoom_lengths.max() + 0.5])
    ax.hist(zoom_lengths, bins=bins, color=CONTINUOUS_PALETTE[9], edgecolor=EDGE_COLOR, linewidth=0.7)
    ax.set_xlim(0, x_max)
    ax.set_xlabel("Per-NUMT length (bp)")
    ax.set_ylabel("NUMT count")
    ax.set_title("C. Per-NUMT length distribution (0-500 bp)", loc="left", fontweight="bold")
    apply_axis_style(ax)
    add_panel_stats_box(
        ax,
        [
            f"n <=500 bp = {len(zoom_lengths):,}",
            f"share = {len(zoom_lengths) / len(lengths):.1%}",
            f"median = {np.median(zoom_lengths):.1f} bp",
        ],
    )


def plot_chr_counts(ax: plt.Axes, chr_df: pd.DataFrame) -> None:
    plot_df = chr_df.copy()
    plot_df["sort_key"] = plot_df["chr"].map(chrom_sort_value)
    plot_df = plot_df.sort_values(["sort_key", "chr"])
    colors = [DISCRETE_PALETTE[idx % len(DISCRETE_PALETTE)] for idx in range(len(plot_df))]
    ax.bar(plot_df["chr"], plot_df["event_count"], color=colors, edgecolor=EDGE_COLOR, linewidth=0.7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Event count")
    ax.set_title("D. Events by chromosome", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    apply_axis_style(ax)


def plot_frequency_class(ax: plt.Axes, freq_class_df: pd.DataFrame) -> None:
    plot_df = (
        freq_class_df.set_index("category")
        .reindex(CLASS_ORDER)
        .fillna({"cluster_count": 0})
        .reset_index()
    )
    ax.bar(
        plot_df["category"],
        plot_df["cluster_count"].astype(int),
        color=[CLASS_COLORS[item] for item in plot_df["category"]],
        edgecolor=EDGE_COLOR,
        linewidth=0.7,
    )
    ax.set_xlabel("Frequency class")
    ax.set_ylabel("Cluster count")
    ax.set_title("E. Frequency class distribution", loc="left", fontweight="bold")
    apply_axis_style(ax)


def plot_support_sensitivity(ax: plt.Axes, support_df: pd.DataFrame) -> None:
    plot_df = support_df.sort_values("min_sample_support").copy()
    ax.plot(
        plot_df["min_sample_support"],
        plot_df["distinct_numt_count"],
        color=DISCRETE_PALETTE[4],
        marker="o",
        linewidth=1.8,
        markersize=4.5,
    )
    ax.set_xlabel("Minimum supporting samples")
    ax.set_ylabel("Distinct NUMT count")
    ax.set_title("F. Support threshold sensitivity", loc="left", fontweight="bold")
    ax.set_xticks(plot_df["min_sample_support"].tolist())
    ax.grid(color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)


def plot_genomewide_recurrence(ax: plt.Axes, freq_cluster_df: pd.DataFrame) -> None:
    plot_df = freq_cluster_df.copy()
    plot_df["chr_order"] = plot_df["chr"].map(chrom_sort_value)
    plot_df = plot_df.loc[plot_df["chr_order"] < 10_000].copy()
    plot_df["CHR"] = plot_df["chr"]
    plot_df = plot_df.sort_values(["chr_order", "cluster_midpoint_mean"])

    chr_sizes = (
        plot_df.groupby(["CHR", "chr_order"], as_index=False)["cluster_max_midpoint"]
        .max()
        .rename(columns={"cluster_max_midpoint": "chr_max"})
        .sort_values("chr_order")
    )
    chr_sizes["offset"] = chr_sizes["chr_max"].cumsum().shift(fill_value=0)
    plot_df = plot_df.merge(chr_sizes, on=["CHR", "chr_order"], how="left")
    plot_df["x"] = plot_df["offset"] + plot_df["cluster_midpoint_mean"].astype(float)
    chr_centers = (
        plot_df.groupby("CHR", as_index=False)["x"].mean().sort_values(by="x")
    )

    class_colors = plot_df["frequency_class"].map(CLASS_COLORS).fillna(DISCRETE_PALETTE[5])
    point_sizes = np.clip(np.sqrt(plot_df["POS"].astype(float).to_numpy()) * 0.8, 8, 60)
    ax.scatter(
        plot_df["x"],
        plot_df["sample_count"].astype(float),
        c=class_colors,
        s=point_sizes,
        alpha=0.6,
        linewidths=0,
    )
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Supporting sample count")
    ax.set_title("G. Genome-wide recurrent loci", loc="left", fontweight="bold")
    ax.set_xticks(chr_centers["x"])
    ax.set_xticklabels(chr_centers["CHR"], rotation=90)
    ax.grid(color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)


def build_summary_text(length_df: pd.DataFrame, freq_cluster_df: pd.DataFrame, support_df: pd.DataFrame) -> str:
    locus_lengths = freq_cluster_df["chosen_region_length"].astype(float)
    sample_counts = freq_cluster_df["sample_count"].astype(int)
    primary_support = support_df.sort_values("min_sample_support").iloc[0]
    return "\n".join(
        [
            f"n regions = {len(length_df):,}",
            f"distinct NUMTs = {len(freq_cluster_df):,}",
            f"distinct length median = {locus_lengths.median():.1f} bp",
            f"median support = {sample_counts.median():.1f} samples",
            f"support>= {int(primary_support['min_sample_support'])}: {int(primary_support['distinct_numt_count']):,}",
        ]
    )


def run(
    length_tsv: str,
    chr_count_tsv: str,
    freq_cluster_tsv: str,
    freq_class_tsv: str,
    support_summary_tsv: str,
    out_prefix: str,
) -> int:
    configure_matplotlib()
    out_path = Path(out_prefix)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    length_df = pd.read_csv(length_tsv, sep="\t")
    chr_df = pd.read_csv(chr_count_tsv, sep="\t")
    freq_cluster_df = pd.read_csv(freq_cluster_tsv, sep="\t")
    freq_class_df = pd.read_csv(freq_class_tsv, sep="\t")
    support_df = pd.read_csv(support_summary_tsv, sep="\t")

    fig, axes = plt.subplots(4, 2, figsize=(18, 21))
    fig.patch.set_facecolor("white")

    plot_length_histogram(axes[0, 0], freq_cluster_df)
    plot_region_length_histogram(axes[0, 1], length_df)
    plot_region_length_histogram_zoom(axes[1, 0], length_df)
    plot_chr_counts(axes[1, 1], chr_df)
    plot_frequency_class(axes[2, 0], freq_class_df)
    plot_support_sensitivity(axes[2, 1], support_df)
    plot_genomewide_recurrence(axes[3, 0], freq_cluster_df)
    axes[3, 1].axis("off")

    summary_text = build_summary_text(length_df, freq_cluster_df, support_df)
    fig.suptitle("NUMT Statistical Description", fontsize=18, fontweight="bold", y=0.995)
    fig.text(
        0.012,
        0.985,
        summary_text,
        ha="left",
        va="top",
        fontsize=10.5,
        color=EDGE_COLOR,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#f8f8f8", "edgecolor": "#d9d9d9"},
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    save_figure(fig, out_path)
    LOG.info("单一多面板图输出完成: %s", out_path)
    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        length_tsv=args.length_tsv,
        chr_count_tsv=args.chr_count_tsv,
        freq_cluster_tsv=args.freq_cluster_tsv,
        freq_class_tsv=args.freq_class_tsv,
        support_summary_tsv=args.support_summary_tsv,
        out_prefix=args.out_prefix,
    )


if __name__ == "__main__":
    raise SystemExit(main())
