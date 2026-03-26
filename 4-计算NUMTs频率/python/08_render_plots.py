#!/usr/bin/env python3
"""Step 08: 统一绘图。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

from numts_frequency.io_utils import ensure_dir, read_delim
from numts_frequency.plot_utils import enable_editable_vector_fonts, plot_bar, plot_histogram, plot_line

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)
logging.getLogger("fontTools").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)


def _read_tsv_or_empty(path: Path, columns: list[str]) -> pd.DataFrame:
    try:
        return read_delim(path, sep="\t")
    except EmptyDataError:
        return pd.DataFrame(columns=columns)


def run(
    main_results_dir: str,
    mt_results_dir: str,
    stratified_dir: str,
    format_main: str,
    format_aux: str,
    dpi: int,
    font_family: str,
) -> int:
    enable_editable_vector_fonts(font_family=font_family)
    main_tables = Path(main_results_dir) / "tables"
    main_plots = ensure_dir(Path(main_results_dir) / "plots")
    mt_tables = Path(mt_results_dir) / "tables"
    mt_plots = ensure_dir(Path(mt_results_dir) / "plots")

    sample_counts = _read_tsv_or_empty(main_tables / "sample_numt_counts.tsv", ["sample_id", "distinct_numt_count"])
    length_distribution = _read_tsv_or_empty(main_tables / "distinct_length_distribution.tsv", ["entity_id", "length_level", "length_bp"])
    frequency_summary = _read_tsv_or_empty(
        main_tables / "frequency_class_summary.tsv",
        [
            "frequency_class",
            "distinct_numt_n",
            "distinct_fraction",
            "carrier_sample_n",
            "carrier_fraction",
            "carrier_count_sum",
            "mean_numt_count_per_sample",
        ],
    )
    carrier_ratio = _read_tsv_or_empty(
        main_tables / "carrier_ratio_by_frequency_class.tsv",
        ["frequency_class", "carrier_sample_n", "carrier_fraction", "sample_count", "sample_ratio"],
    )
    spectrum = _read_tsv_or_empty(main_tables / "distinct_frequency_spectrum.tsv", ["carrier_count", "distinct_count", "carrier_frequency"])
    mt_gene_counts = _read_tsv_or_empty(mt_tables / "mt_gene_breakpoint_counts.tsv", ["mt_primary_gene", "distinct_numt_n"])
    mt_cumulative = _read_tsv_or_empty(
        mt_tables / "mt_breakpoint_cumulative.tsv",
        [
            "mt_position",
            "overall_count",
            "overall_cumulative_count",
            "overall_cumulative_fraction",
            "frequency_class",
            "count",
            "cumulative_count",
            "cumulative_fraction",
        ],
    )

    if not sample_counts.empty:
        plot_histogram(
            sample_counts["distinct_numt_count"],
            "Per-sample distinct NUMT count distribution",
            "Distinct NUMT count",
            Path(main_plots) / "sample_numt_count_distribution",
            format_main,
            format_aux,
            dpi,
        )
    if not length_distribution.empty:
        plot_histogram(
            length_distribution[length_distribution["length_level"] == "distinct"]["length_bp"],
            "Distinct NUMT length distribution",
            "Length (bp)",
            Path(main_plots) / "distinct_numt_length_distribution",
            format_main,
            format_aux,
            dpi,
        )
    if not frequency_summary.empty:
        plot_bar(
            frequency_summary.sort_values("frequency_class"),
            "frequency_class",
            "distinct_numt_n",
            "Frequency class summary",
            "Frequency class",
            "Distinct NUMT count",
            Path(main_plots) / "frequency_class_summary",
            format_main,
            format_aux,
            dpi,
            use_frequency_colors=True,
        )
    if not carrier_ratio.empty:
        plot_bar(
            carrier_ratio.sort_values("frequency_class"),
            "frequency_class",
            "carrier_fraction",
            "Carrier ratio by frequency class",
            "Frequency class",
            "Carrier fraction",
            Path(main_plots) / "carrier_ratio_by_frequency_class",
            format_main,
            format_aux,
            dpi,
            use_frequency_colors=True,
        )
    if not spectrum.empty:
        plot_bar(
            spectrum,
            "carrier_count",
            "distinct_count",
            "Distinct NUMT site frequency spectrum",
            "Carrier count",
            "Distinct NUMT count",
            Path(main_plots) / "distinct_frequency_spectrum",
            format_main,
            format_aux,
            dpi,
        )
    if not mt_gene_counts.empty:
        plot_bar(
            mt_gene_counts.head(20),
            "mt_primary_gene",
            "distinct_numt_n",
            "MT gene breakpoint count",
            "MT gene",
            "Distinct NUMT count",
            Path(mt_plots) / "mt_gene_breakpoint_counts",
            format_main,
            format_aux,
            dpi,
        )
    if not mt_cumulative.empty:
        plot_line(
            mt_cumulative[["mt_position", "overall_cumulative_fraction"]].drop_duplicates(),
            "mt_position",
            "overall_cumulative_fraction",
            None,
            "MT breakpoint cumulative distribution",
            "mtDNA position",
            "Cumulative fraction",
            Path(mt_plots) / "mt_breakpoint_cumulative_overall",
            format_main,
            format_aux,
            dpi,
        )
    if "frequency_class" in mt_cumulative.columns:
        mt_class = mt_cumulative.dropna(subset=["frequency_class", "cumulative_fraction"]).copy()
        if not mt_class.empty:
            plot_line(
                mt_class,
                "mt_position",
                "cumulative_fraction",
                "frequency_class",
                "MT breakpoint cumulative distribution by frequency class",
                "mtDNA position",
                "Cumulative fraction",
                Path(mt_plots) / "mt_breakpoint_cumulative_by_frequency_class",
                format_main,
                format_aux,
                dpi,
            )

    for dimension in ["region", "population", "sex"]:
        dim_dir = Path(stratified_dir) / dimension
        plots_dir = ensure_dir(dim_dir / "plots")
        length_df = _read_tsv_or_empty(dim_dir / "subgroup_length_distribution.tsv", ["subgroup", "entity_id", "length_level", "length_bp"])
        freq_df = _read_tsv_or_empty(dim_dir / "subgroup_frequency_summary.tsv", ["subgroup", "frequency_class", "distinct_numt_n"])
        mt_df = _read_tsv_or_empty(
            dim_dir / "subgroup_mt_breakpoint_cumulative.tsv",
            [
                "mt_position",
                "overall_count",
                "overall_cumulative_count",
                "overall_cumulative_fraction",
                "frequency_class",
                "count",
                "cumulative_count",
                "cumulative_fraction",
                "subgroup",
            ],
        )
        chr_df = _read_tsv_or_empty(dim_dir / "subgroup_nuclear_chr_distribution.tsv", ["subgroup", "chr", "distinct_numt_n"])

        if not length_df.empty:
            for subgroup, sub in length_df.groupby("subgroup", sort=True):
                plot_histogram(
                    sub["length_bp"],
                    f"{dimension}:{subgroup} length distribution",
                    "Length (bp)",
                    plots_dir / f"{subgroup}_length_distribution",
                    format_main,
                    format_aux,
                    dpi,
                )
        if not freq_df.empty:
            for subgroup, sub in freq_df.groupby("subgroup", sort=True):
                plot_bar(
                    sub.sort_values("frequency_class"),
                    "frequency_class",
                    "distinct_numt_n",
                    f"{dimension}:{subgroup} frequency summary",
                    "Frequency class",
                    "Distinct NUMT count",
                    plots_dir / f"{subgroup}_frequency_summary",
                    format_main,
                    format_aux,
                    dpi,
                    use_frequency_colors=True,
                )
        if not mt_df.empty:
            for subgroup, sub in mt_df.groupby("subgroup", sort=True):
                overall_col = "overall_cumulative_fraction"
                plot_line(
                    sub.drop_duplicates(subset=["mt_position", overall_col]),
                    "mt_position",
                    overall_col,
                    None,
                    f"{dimension}:{subgroup} mt breakpoint cumulative",
                    "mtDNA position",
                    "Cumulative fraction",
                    plots_dir / f"{subgroup}_mt_breakpoint_cumulative",
                    format_main,
                    format_aux,
                    dpi,
                )
        if not chr_df.empty:
            for subgroup, sub in chr_df.groupby("subgroup", sort=True):
                plot_bar(
                    sub.sort_values("chr"),
                    "chr",
                    "distinct_numt_n",
                    f"{dimension}:{subgroup} nuclear chromosome distribution",
                    "Chromosome",
                    "Distinct NUMT count",
                    plots_dir / f"{subgroup}_nuclear_chr_distribution",
                    format_main,
                    format_aux,
                    dpi,
                )
    LOG.info("Plots rendered")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render plots from NUMT frequency result tables")
    parser.add_argument("--main-results-dir", required=True)
    parser.add_argument("--mt-results-dir", required=True)
    parser.add_argument("--stratified-dir", required=True)
    parser.add_argument("--format-main", required=True)
    parser.add_argument("--format-aux", required=True)
    parser.add_argument("--dpi", required=True, type=int)
    parser.add_argument("--font-family", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
