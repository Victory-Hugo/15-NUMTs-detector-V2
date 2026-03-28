#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""绘制长度分布、染色体计数和频率分级图。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)

BAR_FILL = "#49c1ad"
BAR_EDGE = "#30203e"
FREQ_COLORS = {
    "common": "#0072b2",
    "low-frequency": "#56b4e9",
    "rare": "#009e73",
    "ultra-rare": "#f0e442",
}
TAIL_COLORS = ["#0072b2", "#56b4e9", "#009e73", "#f0e442"]


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def configure_matplotlib() -> None:
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="生成长度图、染色体计数图和频率图。")
    parser.add_argument("--length-tsv", required=True)
    parser.add_argument("--chr-count-tsv", required=True)
    parser.add_argument("--freq-cluster-tsv", required=True)
    parser.add_argument("--freq-class-tsv", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser


def save_figure(fig: plt.Figure, base_path: Path) -> None:
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def choose_histogram_bins(lengths: np.ndarray) -> np.ndarray:
    if lengths.size == 0:
        raise ValueError("长度数组为空")
    if np.all(lengths == lengths[0]):
        return np.array([lengths[0] - 0.5, lengths[0] + 0.5])
    bins = np.histogram_bin_edges(lengths, bins="fd")
    bin_count = max(1, len(bins) - 1)
    if bin_count < 15 or bin_count > 80:
        bins = np.histogram_bin_edges(lengths, bins="sturges")
    return np.unique(np.round(bins).astype(int))


def chrom_sort_value(chrom: str) -> int:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return int(suffix)
        if suffix == "X":
            return 23
        if suffix == "Y":
            return 24
    return 10_000


def plot_length_histogram(length_df: pd.DataFrame, out_dir: Path) -> None:
    lengths = length_df["chosen_length"].astype(float).to_numpy()
    bins = choose_histogram_bins(lengths)
    if bins.size < 2:
        bins = np.array([lengths.min() - 0.5, lengths.max() + 0.5])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(lengths, bins=bins, color=BAR_FILL, edgecolor=BAR_EDGE, linewidth=0.8)
    ax.set_xlabel("NUMT length (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title("NUMT length distribution")
    ax.grid(axis="y", alpha=0.2)
    save_figure(fig, out_dir / "2-numt-length-histogram")


def plot_length_ecdf(length_df: pd.DataFrame, out_dir: Path) -> None:
    lengths = np.sort(length_df["chosen_length"].astype(float).to_numpy())
    y_values = np.arange(1, len(lengths) + 1) / len(lengths)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(lengths, y_values, where="post", color=BAR_EDGE, linewidth=1.5)
    ax.fill_between(lengths, y_values, step="post", color=BAR_FILL, alpha=0.25)
    ax.set_xlabel("NUMT length (bp)")
    ax.set_ylabel("ECDF")
    ax.set_title("NUMT length ECDF")
    ax.grid(alpha=0.2)
    save_figure(fig, out_dir / "2-numt-length-ecdf")


def plot_chr_counts(chr_df: pd.DataFrame, out_dir: Path) -> None:
    plot_df = chr_df.copy()
    plot_df["sort_key"] = plot_df["chr"].map(chrom_sort_value)
    plot_df = plot_df.sort_values(["sort_key", "chr"])

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(plot_df["chr"], plot_df["event_count"], color=BAR_FILL, edgecolor=BAR_EDGE, linewidth=0.8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Event count")
    ax.set_title("NUMT event counts by chromosome")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.2)
    save_figure(fig, out_dir / "3-numt-chromosome-counts")


def plot_frequency_class(freq_class_df: pd.DataFrame, out_dir: Path) -> None:
    class_order = ["common", "low-frequency", "rare", "ultra-rare"]
    plot_df = (
        freq_class_df.set_index("category")
        .reindex(class_order)
        .fillna({"cluster_count": 0})
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(
        plot_df["category"],
        plot_df["cluster_count"].astype(int),
        color=[FREQ_COLORS[item] for item in plot_df["category"]],
        edgecolor=BAR_EDGE,
        linewidth=0.8,
    )
    ax.set_xlabel("Frequency class")
    ax.set_ylabel("Cluster count")
    ax.set_title("NUMT frequency class distribution")
    ax.grid(axis="y", alpha=0.2)
    save_figure(fig, out_dir / "4-numt-frequency-class-barplot")


def plot_frequency_tail(freq_cluster_df: pd.DataFrame, out_dir: Path) -> None:
    sample_counts = freq_cluster_df["sample_count"].astype(int)
    tail_series = pd.Series(
        {
            "singleton": int((sample_counts == 1).sum()),
            "doubleton": int((sample_counts == 2).sum()),
            "tripleton": int((sample_counts == 3).sum()),
            "4-8 samples": int(((sample_counts >= 4) & (sample_counts <= 8)).sum()),
        }
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(tail_series.index, tail_series.values, color=TAIL_COLORS, edgecolor=BAR_EDGE, linewidth=0.8)
    ax.set_xlabel("Rare tail class")
    ax.set_ylabel("Cluster count")
    ax.set_title("NUMT rare-tail distribution")
    ax.grid(axis="y", alpha=0.2)
    save_figure(fig, out_dir / "4-numt-frequency-tail-barplot")


def run(length_tsv: str, chr_count_tsv: str, freq_cluster_tsv: str, freq_class_tsv: str, out_dir: str) -> int:
    configure_matplotlib()
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    length_df = pd.read_csv(length_tsv, sep="\t")
    chr_df = pd.read_csv(chr_count_tsv, sep="\t")
    freq_cluster_df = pd.read_csv(freq_cluster_tsv, sep="\t")
    freq_class_df = pd.read_csv(freq_class_tsv, sep="\t")

    plot_length_histogram(length_df, out_path)
    plot_length_ecdf(length_df, out_path)
    plot_chr_counts(chr_df, out_path)
    plot_frequency_class(freq_class_df, out_path)
    plot_frequency_tail(freq_cluster_df, out_path)
    LOG.info("Python 图形输出完成")
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
        out_dir=args.out_dir,
    )


if __name__ == "__main__":
    raise SystemExit(main())
