#!/usr/bin/env python3
"""统一绘图工具。"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, to_hex


DISCRETE_PRIMARY = ["#0072b2", "#56b4e9", "#009e73", "#f0e442", "#e69f00", "#d55e00"]
DISCRETE_SECONDARY = ["#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500"]
FREQUENCY_COLORS = {
    "prevalent": "#0072b2",
    "common": "#56b4e9",
    "rare": "#009e73",
    "ultra_rare": "#e69f00",
    "private": "#d55e00",
}


def enable_editable_vector_fonts(font_family: str = "DejaVu Sans") -> None:
    plt.rcParams["font.sans-serif"] = [font_family]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"


def get_discrete_palette(n: int) -> list[str]:
    if n <= 6:
        return DISCRETE_PRIMARY[:n]
    if n <= 11:
        return DISCRETE_PRIMARY + DISCRETE_SECONDARY[: n - 6]
    anchors = DISCRETE_PRIMARY + DISCRETE_SECONDARY
    cmap = LinearSegmentedColormap.from_list("custom_discrete_interp", anchors)
    return [to_hex(cmap(x)) for x in np.linspace(0, 1, n)]


def save_figure(fig: plt.Figure, output_base: Path, format_main: str, format_aux: str, dpi: int) -> None:
    fig.tight_layout()
    fig.savefig(output_base.with_suffix(f".{format_main}"), dpi=dpi)
    if format_aux and format_aux != format_main:
        fig.savefig(output_base.with_suffix(f".{format_aux}"), dpi=dpi)
    plt.close(fig)


def plot_histogram(series: pd.Series, title: str, xlabel: str, output_base: Path, format_main: str, format_aux: str, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(series.dropna(), bins=50, color=DISCRETE_PRIMARY[0], alpha=0.85)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    save_figure(fig, output_base, format_main, format_aux, dpi)


def plot_bar(df: pd.DataFrame, x_col: str, y_col: str, title: str, xlabel: str, ylabel: str, output_base: Path, format_main: str, format_aux: str, dpi: int, use_frequency_colors: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    x_values = df[x_col].astype(str).tolist()
    if use_frequency_colors:
      colors = [FREQUENCY_COLORS.get(value, DISCRETE_PRIMARY[0]) for value in x_values]
    else:
      colors = get_discrete_palette(len(x_values))
    ax.bar(x_values, df[y_col], color=colors)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=45)
    save_figure(fig, output_base, format_main, format_aux, dpi)


def plot_line(df: pd.DataFrame, x_col: str, y_col: str, group_col: str | None, title: str, xlabel: str, ylabel: str, output_base: Path, format_main: str, format_aux: str, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    if group_col is None:
        ax.plot(df[x_col], df[y_col], color=DISCRETE_PRIMARY[0], linewidth=2)
    else:
        groups = sorted(df[group_col].dropna().astype(str).unique())
        palette = {group: FREQUENCY_COLORS.get(group, color) for group, color in zip(groups, get_discrete_palette(len(groups)))}
        for group in groups:
            sub = df[df[group_col].astype(str) == group]
            ax.plot(sub[x_col], sub[y_col], label=group, color=palette[group], linewidth=2)
        ax.legend(frameon=False)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    save_figure(fig, output_base, format_main, format_aux, dpi)
