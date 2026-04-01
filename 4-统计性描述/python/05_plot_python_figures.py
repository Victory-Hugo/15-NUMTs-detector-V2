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
CLASS_LABELS = {
    "common": "Common",
    "low-frequency": "Low-frequency",
    "rare": "Rare",
    "ultra-rare": "Ultra-rare",
    "private": "Private",
}
PRIVATE_CATEGORY = "private"
PRIVATE_COLOR = DISCRETE_PALETTE[4]
CATEGORY_COLORS = {**CLASS_COLORS, PRIVATE_CATEGORY: PRIVATE_COLOR}


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
    parser.add_argument("--scatter-tsv")
    parser.add_argument("--relative-frequency-tsv")
    parser.add_argument("--out-prefix", required=True)
    return parser


def save_figure(fig: plt.Figure, out_prefix: Path) -> None:
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_scatter_table(plot_df: pd.DataFrame, out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    export_df = plot_df.rename(
        columns={
            "midpoint": "region_midpoint",
            "chosen_length": "numt_size_bp",
            "frequency": "numt_frequency",
            "neg_log2_frequency": "neg_log2_numt_frequency",
        }
    )[
        [
            "sampleID",
            "region_id",
            "region_chr",
            "region_start",
            "region_end",
            "region_midpoint",
            "numt_size_bp",
            "numt_frequency",
            "neg_log2_numt_frequency",
            "frequency_class",
            "cluster_min_midpoint",
            "cluster_max_midpoint",
        ]
    ].copy()
    export_df.to_csv(out_tsv, sep="\t", index=False)


def save_relative_frequency_table(plot_df: pd.DataFrame, out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    export_df = plot_df.rename(
        columns={
            "category": "frequency_category",
            "label": "frequency_category_label",
            "mean_frequency": "mean_numt_frequency",
            "relative_percentage": "relative_percentage_vs_common",
        }
    )[
        [
            "sort_order",
            "frequency_category",
            "frequency_category_label",
            "locus_count",
            "mean_numt_frequency",
            "relative_percentage_vs_common",
        ]
    ].copy()
    export_df["common_reference"] = "Common mean frequency = 100%"
    export_df["private_definition"] = export_df["frequency_category"].map(
        lambda value: "singleton distinct NUMT loci (sample_count == 1)" if value == PRIVATE_CATEGORY else ""
    )
    export_df["relationship_to_ultra_rare"] = export_df["frequency_category"].map(
        lambda value: "subset_of_ultra_rare" if value == PRIVATE_CATEGORY else ""
    )
    export_df.to_csv(out_tsv, sep="\t", index=False)


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


def build_relative_frequency_percentage_df(freq_cluster_df: pd.DataFrame) -> pd.DataFrame:
    working_df = freq_cluster_df.copy()
    working_df["sample_count"] = pd.to_numeric(working_df["sample_count"], errors="coerce")
    working_df["frequency"] = pd.to_numeric(working_df["frequency"], errors="coerce")
    working_df = working_df.dropna(subset=["sample_count", "frequency", "frequency_class"]).copy()
    if working_df.empty:
        raise ValueError("频率聚类表为空，无法绘制频率百分比子图")

    category_order = ["common", "low-frequency", "rare", "ultra-rare", PRIVATE_CATEGORY]
    rows: list[dict[str, float | int | str]] = []
    for category in CLASS_ORDER:
        class_df = working_df.loc[working_df["frequency_class"] == category].copy()
        rows.append(
            {
                "category": category,
                "mean_frequency": float(class_df["frequency"].mean()) if not class_df.empty else 0.0,
                "locus_count": int(len(class_df)),
            }
        )

    private_df = working_df.loc[working_df["sample_count"].astype(int) == 1].copy()
    rows.append(
        {
            "category": PRIVATE_CATEGORY,
            "mean_frequency": float(private_df["frequency"].mean()) if not private_df.empty else 0.0,
            "locus_count": int(len(private_df)),
        }
    )

    plot_df = pd.DataFrame(rows)
    common_reference = float(
        plot_df.loc[plot_df["category"] == "common", "mean_frequency"].iloc[0]
    )
    if common_reference <= 0:
        raise ValueError("Common 频率均值 <= 0，无法将其归一化为 100%")

    plot_df["relative_percentage"] = plot_df["mean_frequency"] / common_reference * 100.0
    plot_df["label"] = plot_df["category"].map(CLASS_LABELS)
    plot_df["color"] = plot_df["category"].map(CATEGORY_COLORS)
    plot_df["sort_order"] = plot_df["category"].map({category: idx for idx, category in enumerate(category_order)})
    return plot_df.sort_values("sort_order").reset_index(drop=True)


def plot_relative_frequency_percentage(ax: plt.Axes, freq_cluster_df: pd.DataFrame) -> None:
    plot_df = build_relative_frequency_percentage_df(freq_cluster_df)
    y_pos = np.arange(len(plot_df))
    bar_values = plot_df["relative_percentage"].astype(float).to_numpy()

    ax.barh(
        y_pos,
        bar_values,
        color=plot_df["color"].tolist(),
        edgecolor=EDGE_COLOR,
        linewidth=0.7,
        height=0.68,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df["label"])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Percentage (Common mean frequency = 100%)")
    ax.set_ylabel("Frequency class")
    ax.set_title("I. Relative frequency scale across classes", loc="left", fontweight="bold")
    ax.grid(axis="x", color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)

    min_positive = max(float(bar_values[bar_values > 0].min()), 1e-4)
    ax.set_xlim(min_positive / 1.6, 140)

    for idx, (_, row) in enumerate(plot_df.iterrows()):
        value = float(row["relative_percentage"])
        label_text = f"{value:.2f}%"
        x_pos = min(value * 1.15, 132.0)
        ha = "left"
        if value >= 85:
            x_pos = value / 1.15
            ha = "right"
        ax.text(x_pos, idx, label_text, va="center", ha=ha, fontsize=9.5, color=EDGE_COLOR)

    private_count = int(plot_df.loc[plot_df["category"] == PRIVATE_CATEGORY, "locus_count"].iloc[0])
    ultra_rare_count = int(plot_df.loc[plot_df["category"] == "ultra-rare", "locus_count"].iloc[0])
    add_panel_stats_box(
        ax,
        [
            "Private = singleton distinct loci",
            "Private is a subset of Ultra-rare",
            f"n(private) = {private_count:,}; n(ultra-rare) = {ultra_rare_count:,}",
        ],
    )


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


def build_region_frequency_length_df(length_df: pd.DataFrame, freq_cluster_df: pd.DataFrame) -> pd.DataFrame:
    region_df = length_df.copy()
    region_df["region_start"] = pd.to_numeric(region_df["region_start"], errors="coerce")
    region_df["region_end"] = pd.to_numeric(region_df["region_end"], errors="coerce")
    region_df["chosen_length"] = pd.to_numeric(region_df["chosen_length"], errors="coerce")
    region_df = region_df.dropna(subset=["region_chr", "region_start", "region_end", "chosen_length"]).copy()
    region_df["midpoint"] = ((region_df["region_start"] + region_df["region_end"]) / 2).astype(int)

    cluster_df = freq_cluster_df.copy()
    cluster_df["cluster_min_midpoint"] = pd.to_numeric(cluster_df["cluster_min_midpoint"], errors="coerce")
    cluster_df["cluster_max_midpoint"] = pd.to_numeric(cluster_df["cluster_max_midpoint"], errors="coerce")
    cluster_df["frequency"] = pd.to_numeric(cluster_df["frequency"], errors="coerce")
    cluster_df = cluster_df.dropna(
        subset=["chr", "cluster_min_midpoint", "cluster_max_midpoint", "frequency", "frequency_class"]
    ).copy()

    merged_df = region_df.merge(
        cluster_df[["chr", "cluster_min_midpoint", "cluster_max_midpoint", "frequency", "frequency_class"]],
        left_on="region_chr",
        right_on="chr",
        how="left",
    )
    plot_df = merged_df.loc[
        (merged_df["midpoint"] >= merged_df["cluster_min_midpoint"])
        & (merged_df["midpoint"] <= merged_df["cluster_max_midpoint"])
        & (merged_df["frequency"] > 0)
    ].copy()
    if plot_df.empty:
        raise ValueError("未找到可用于频率-长度散点图的 region 与 cluster 映射")

    matched_counts = plot_df.groupby("region_id").size()
    duplicated_regions = matched_counts[matched_counts > 1]
    if not duplicated_regions.empty:
        raise ValueError("存在 region 映射到多个 frequency cluster，无法唯一绘图")

    plot_df["neg_log2_frequency"] = -np.log2(plot_df["frequency"])
    return plot_df


def plot_frequency_length_scatter(ax: plt.Axes, length_df: pd.DataFrame, freq_cluster_df: pd.DataFrame) -> None:
    plot_df = build_region_frequency_length_df(length_df=length_df, freq_cluster_df=freq_cluster_df)
    class_colors = plot_df["frequency_class"].map(CLASS_COLORS).fillna(DISCRETE_PALETTE[5])
    ax.scatter(
        plot_df["neg_log2_frequency"],
        plot_df["chosen_length"],
        c=class_colors,
        s=10,
        alpha=0.55,
        edgecolors="none",
    )
    ax.set_xlabel("-log2(frequency of NUMTs)")
    ax.set_ylabel("Size of NUMTs (bp)")
    ax.set_title("H. NUMT size versus frequency", loc="left", fontweight="bold")
    ax.grid(color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    ax.set_axisbelow(True)

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=6,
            markerfacecolor=CLASS_COLORS[class_name],
            markeredgecolor="none",
            label=class_name,
        )
        for class_name in CLASS_ORDER
    ]
    ax.legend(
        handles=legend_handles,
        title="Frequency class",
        frameon=False,
        loc="upper right",
    )
    add_panel_stats_box(
        ax,
        [
            f"n = {len(plot_df):,}",
            f"median X = {plot_df['neg_log2_frequency'].median():.2f}",
            f"median Y = {plot_df['chosen_length'].median():.1f} bp",
        ],
    )


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
    scatter_tsv: str | None,
    relative_frequency_tsv: str | None,
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
    scatter_plot_df = build_region_frequency_length_df(length_df=length_df, freq_cluster_df=freq_cluster_df)

    if scatter_tsv:
        scatter_out_path = Path(scatter_tsv)
    else:
        scatter_out_path = out_path.with_name(f"{out_path.name}-frequency-length-scatter.tsv")
    save_scatter_table(scatter_plot_df, scatter_out_path)
    relative_frequency_plot_df = build_relative_frequency_percentage_df(freq_cluster_df)
    if relative_frequency_tsv:
        relative_frequency_out_path = Path(relative_frequency_tsv)
    else:
        relative_frequency_out_path = out_path.with_name(f"{out_path.name}-relative-frequency-percentage.tsv")
    save_relative_frequency_table(relative_frequency_plot_df, relative_frequency_out_path)

    fig = plt.figure(figsize=(18, 25))
    gridspec = fig.add_gridspec(5, 2, height_ratios=[1, 1, 1, 1, 0.95])
    axes = np.empty((5, 2), dtype=object)
    axes[0, 0] = fig.add_subplot(gridspec[0, 0])
    axes[0, 1] = fig.add_subplot(gridspec[0, 1])
    axes[1, 0] = fig.add_subplot(gridspec[1, 0])
    axes[1, 1] = fig.add_subplot(gridspec[1, 1])
    axes[2, 0] = fig.add_subplot(gridspec[2, 0])
    axes[2, 1] = fig.add_subplot(gridspec[2, 1])
    axes[3, 0] = fig.add_subplot(gridspec[3, 0])
    axes[3, 1] = fig.add_subplot(gridspec[3, 1])
    axes[4, 0] = fig.add_subplot(gridspec[4, :])
    axes[4, 1] = None
    fig.patch.set_facecolor("white")

    plot_length_histogram(axes[0, 0], freq_cluster_df)
    plot_region_length_histogram(axes[0, 1], length_df)
    plot_region_length_histogram_zoom(axes[1, 0], length_df)
    plot_chr_counts(axes[1, 1], chr_df)
    plot_frequency_class(axes[2, 0], freq_class_df)
    plot_support_sensitivity(axes[2, 1], support_df)
    plot_genomewide_recurrence(axes[3, 0], freq_cluster_df)
    class_colors = scatter_plot_df["frequency_class"].map(CLASS_COLORS).fillna(DISCRETE_PALETTE[5])
    axes[3, 1].scatter(
        scatter_plot_df["neg_log2_frequency"],
        scatter_plot_df["chosen_length"],
        c=class_colors,
        s=10,
        alpha=0.55,
        edgecolors="none",
    )
    axes[3, 1].set_xlabel("-log2(frequency of NUMTs)")
    axes[3, 1].set_ylabel("Size of NUMTs (bp)")
    axes[3, 1].set_title("H. NUMT size versus frequency", loc="left", fontweight="bold")
    axes[3, 1].grid(color=GRID_COLOR, linewidth=0.8, alpha=0.5)
    axes[3, 1].set_axisbelow(True)
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=6,
            markerfacecolor=CLASS_COLORS[class_name],
            markeredgecolor="none",
            label=class_name,
        )
        for class_name in CLASS_ORDER
    ]
    axes[3, 1].legend(handles=legend_handles, title="Frequency class", frameon=False, loc="upper right")
    add_panel_stats_box(
        axes[3, 1],
        [
            f"n = {len(scatter_plot_df):,}",
            f"median X = {scatter_plot_df['neg_log2_frequency'].median():.2f}",
            f"median Y = {scatter_plot_df['chosen_length'].median():.1f} bp",
        ],
    )
    plot_relative_frequency_percentage(axes[4, 0], freq_cluster_df)

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
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    save_figure(fig, out_path)
    LOG.info("单一多面板图输出完成: %s", out_path)
    LOG.info("频率-长度散点图数据表输出完成: %s", scatter_out_path)
    LOG.info("频率百分比子图数据表输出完成: %s", relative_frequency_out_path)
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
        scatter_tsv=args.scatter_tsv,
        relative_frequency_tsv=args.relative_frequency_tsv,
        out_prefix=args.out_prefix,
    )


if __name__ == "__main__":
    raise SystemExit(main())
