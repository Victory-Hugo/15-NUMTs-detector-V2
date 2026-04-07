#!/usr/bin/env python3
"""Plot enrichment summary heatmap and barplot."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

log = logging.getLogger(__name__)

CLASS_ORDER = ["all", "common", "low-frequency", "rare", "ultra-rare"]
CLASS_COLORS = {
    "all": "#0072b2",
    "common": "#56b4e9",
    "low-frequency": "#009e73",
    "rare": "#f0e442",
    "ultra-rare": "#e69f00",
}
REGION_LABELS = {
    "1-启动子区域.bed": "Promoter",
    "2-蛋白编码区域.bed": "Protein-coding",
    "3-大规模片段重复区域.bed": "Segmental duplication",
    "4-去除片段重复区域的人类STR区域.bed": "STR (segdup removed)",
    "05-CpG_islands.bed": "CpG islands",
    "06-Microsat.bed": "Microsat",
    "07-Start_codon.bed": "Start codon",
    "08-snRNA.bed": "snRNA",
    "09-FuncElems.bed": "FuncElems",
    "10-Retroposon.bed": "Retroposon",
    "11-CDS.bed": "CDS",
    "12-Stop_codon.bed": "Stop codon",
    "13-SINE.bed": "SINE",
    "14-rmsk-DNA.bed": "rmsk-DNA",
    "15-Centromeres.bed": "Centromeres",
    "16-Satellite.bed": "Satellite",
    "17-LTR.bed": "LTR",
    "18-Exon.bed": "Exon",
    "19-UTR.bed": "UTR",
    "20-LINE.bed": "LINE",
    "21-Simple_Repeats.bed": "Simple repeats",
    "22-Genomic_superdups.bed": "Genomic superdups",
    "23-snoRNA_miRNA.bed": "snoRNA/miRNA",
    "24-Gene.bed": "Gene",
    "25-Intron.bed": "Intron",
    "26-Regulatory_elements.bed": "Regulatory elements",
}
SINGLE_HUE_COLORS = [
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


def display_region_name(region_name: str) -> str:
    """Return an ASCII label suitable for manuscript figures."""
    return REGION_LABELS.get(str(region_name), str(region_name))


def save_figure(fig: plt.Figure, output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(output_prefix.with_suffix(suffix), dpi=300, bbox_inches="tight")


def run(summary_tsv: str | Path, output_dir: str | Path) -> dict[str, str]:
    summary_path = Path(summary_tsv)
    output_path = Path(output_dir)
    if not summary_path.is_file():
        raise FileNotFoundError(f"Summary TSV not found: {summary_path}")

    df = pd.read_csv(summary_path, sep="\t")
    if df.empty:
        raise ValueError(f"Summary TSV is empty: {summary_path}")
    df["frequency_class"] = pd.Categorical(df["frequency_class"], categories=CLASS_ORDER, ordered=True)
    df["region_label"] = df["region_name"].map(display_region_name)
    df = df.sort_values(["region_name", "frequency_class"]).copy()

    heatmap_df = df.pivot(index="region_label", columns="frequency_class", values="p_value_greater")
    heatmap_df = heatmap_df.reindex(columns=[item for item in CLASS_ORDER if item in heatmap_df.columns])
    heatmap_values = -np.log10(heatmap_df.clip(lower=1e-300).astype(float))

    cmap = LinearSegmentedColormap.from_list("continuous_single_2", SINGLE_HUE_COLORS)
    fig, ax = plt.subplots(figsize=(max(7, 1.2 * heatmap_values.shape[1]), max(4, 0.6 * heatmap_values.shape[0] + 2)))
    image = ax.imshow(heatmap_values.values, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(heatmap_values.columns)))
    ax.set_xticklabels(heatmap_values.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(heatmap_values.index)))
    ax.set_yticklabels(heatmap_values.index)
    ax.set_title("NUMT nuclear breakpoint flank enrichment")
    ax.set_xlabel("NUMT frequency class")
    ax.set_ylabel("Target region")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("-log10(P value for enrichment)")
    fig.tight_layout()
    heatmap_prefix = output_path / "enrichment_pvalue_heatmap"
    save_figure(fig, heatmap_prefix)
    plt.close(fig)

    bar_df = df.copy()
    regions = list(dict.fromkeys(bar_df["region_label"].tolist()))
    classes = [item for item in CLASS_ORDER if item in set(bar_df["frequency_class"].astype(str))]
    x = np.arange(len(regions))
    width = 0.8 / max(1, len(classes))
    fig, ax = plt.subplots(figsize=(max(8, 1.6 * len(regions)), 5))
    for idx, frequency_class in enumerate(classes):
        class_values = []
        for region in regions:
            match = bar_df.loc[
                (bar_df["region_label"] == region) & (bar_df["frequency_class"].astype(str) == frequency_class),
                "observed_percentage",
            ]
            class_values.append(float(match.iloc[0]) if not match.empty else 0.0)
        ax.bar(
            x + (idx - (len(classes) - 1) / 2) * width,
            class_values,
            width=width,
            label=frequency_class,
            color=CLASS_COLORS.get(frequency_class, "#0072b2"),
        )
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=35, ha="right")
    ax.set_ylabel("Observed overlap percentage")
    ax.set_xlabel("Target region")
    ax.set_title("Observed NUMT breakpoint flank overlap")
    ax.legend(frameon=False)
    fig.tight_layout()
    barplot_prefix = output_path / "enrichment_observed_percentage_barplot"
    save_figure(fig, barplot_prefix)
    plt.close(fig)

    log.info("Wrote figures to %s", output_path)
    return {"heatmap": str(heatmap_prefix), "barplot": str(barplot_prefix)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot enrichment summary figures.")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
