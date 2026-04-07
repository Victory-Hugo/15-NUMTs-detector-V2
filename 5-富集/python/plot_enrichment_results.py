#!/usr/bin/env python3
"""Plot enrichment summary barplot."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

    log.info("Wrote barplot to %s", output_path)
    return {"barplot": str(barplot_prefix)}


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
