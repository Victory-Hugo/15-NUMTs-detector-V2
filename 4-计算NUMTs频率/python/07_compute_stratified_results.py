#!/usr/bin/env python3
"""Step 07: 分层统计。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from numts_frequency.io_utils import ensure_dir, read_delim, read_meta, write_tsv
from numts_frequency.stats_utils import compute_mt_breakpoint_tables, recompute_catalog_for_subset

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def _write_dimension_outputs(
    dimension_dir: Path,
    summary_rows: list[dict],
    subgroup_metrics_rows: list[dict],
    subgroup_length_rows: list[dict],
    subgroup_frequency_rows: list[dict],
    subgroup_mt_rows: list[pd.DataFrame],
    subgroup_chr_rows: list[dict],
) -> None:
    write_tsv(
        pd.DataFrame(summary_rows, columns=["subgroup", "sample_n", "warning"]),
        str(dimension_dir / "summary.tsv"),
    )
    write_tsv(
        pd.DataFrame(subgroup_metrics_rows, columns=["subgroup", "sample_n", "mean_numt_per_sample", "distinct_numt_n"]),
        str(dimension_dir / "subgroup_metrics.tsv"),
    )
    write_tsv(
        pd.DataFrame(subgroup_length_rows, columns=["subgroup", "entity_id", "length_level", "length_bp"]),
        str(dimension_dir / "subgroup_length_distribution.tsv"),
    )
    write_tsv(
        pd.DataFrame(subgroup_frequency_rows, columns=["subgroup", "frequency_class", "distinct_numt_n"]),
        str(dimension_dir / "subgroup_frequency_summary.tsv"),
    )
    mt_df = (
        pd.concat(subgroup_mt_rows, ignore_index=True)
        if subgroup_mt_rows
        else pd.DataFrame(
            columns=[
                "mt_position",
                "overall_count",
                "overall_cumulative_count",
                "overall_cumulative_fraction",
                "frequency_class",
                "count",
                "cumulative_count",
                "cumulative_fraction",
                "subgroup",
            ]
        )
    )
    write_tsv(mt_df, str(dimension_dir / "subgroup_mt_breakpoint_cumulative.tsv"))
    write_tsv(
        pd.DataFrame(subgroup_chr_rows, columns=["subgroup", "chr", "distinct_numt_n"]),
        str(dimension_dir / "subgroup_nuclear_chr_distribution.tsv"),
    )


def run(
    catalog_tsv: str,
    sample_map_tsv: str,
    sample_events_tsv: str,
    sample_inclusion_tsv: str,
    meta_tsv: str,
    meta_id_col: str,
    geography_col: str,
    ethnicity_col: str,
    sex_col: str,
    prevalent_min: float,
    common_min: float,
    rare_min: float,
    min_group_n_warn: int,
    output_dir: str,
) -> int:
    ensure_dir(output_dir)
    catalog = read_delim(catalog_tsv, sep="\t")
    sample_map = read_delim(sample_map_tsv, sep="\t")
    sample_events = read_delim(sample_events_tsv, sep="\t")
    sample_inclusion = read_delim(sample_inclusion_tsv, sep="\t")
    meta = read_meta(meta_tsv, meta_id_col)
    catalog["distinct_numt_id"] = catalog["distinct_numt_id"].astype(str)
    sample_map["distinct_numt_id"] = sample_map["distinct_numt_id"].astype(str)
    sample_map["sample_id"] = sample_map["sample_id"].astype(str)
    sample_map["sample_event_id"] = sample_map["sample_event_id"].astype(str)
    sample_events["sample_event_id"] = sample_events["sample_event_id"].astype(str)
    sample_inclusion["sample_id"] = sample_inclusion["sample_id"].astype(str)
    meta["sample_id"] = meta["sample_id"].astype(str)
    retained_meta = sample_inclusion[["sample_id"]].merge(meta, on="sample_id", how="left")

    dimensions = {
        "region": geography_col,
        "population": ethnicity_col,
        "sex": sex_col,
    }

    for dim_name, meta_col in dimensions.items():
        dim_dir = ensure_dir(Path(output_dir) / dim_name)
        summary_rows = []
        subgroup_metrics_rows = []
        subgroup_length_rows = []
        subgroup_frequency_rows = []
        subgroup_mt_rows = []
        subgroup_chr_rows = []
        retained_meta[meta_col] = retained_meta[meta_col].fillna("Unknown").astype(str)

        for subgroup, subgroup_meta in retained_meta.groupby(meta_col, sort=True):
            subgroup_samples = subgroup_meta["sample_id"].astype(str).tolist()
            subgroup_catalog = recompute_catalog_for_subset(
                catalog,
                sample_map,
                subgroup_samples,
                len(subgroup_samples),
                prevalent_min,
                common_min,
                rare_min,
            )
            subgroup_map = sample_map[sample_map["sample_id"].isin(subgroup_samples)].copy()
            subgroup_event_ids = subgroup_map["sample_event_id"].astype(str).unique().tolist()
            subgroup_events = sample_events[sample_events["sample_event_id"].astype(str).isin(subgroup_event_ids)].copy()
            sample_level_counts = (
                subgroup_map[["sample_id", "distinct_numt_id"]]
                .drop_duplicates()
                .groupby("sample_id")
                .size()
                .reindex(subgroup_samples, fill_value=0)
            )
            warning = "sample_n_lt_threshold" if len(subgroup_samples) < min_group_n_warn else ""
            summary_rows.append(
                {
                    "subgroup": subgroup,
                    "sample_n": len(subgroup_samples),
                    "warning": warning,
                }
            )
            subgroup_metrics_rows.append(
                {
                    "subgroup": subgroup,
                    "sample_n": len(subgroup_samples),
                    "mean_numt_per_sample": float(sample_level_counts.mean()) if len(sample_level_counts) > 0 else 0.0,
                    "distinct_numt_n": subgroup_catalog["distinct_numt_id"].nunique(),
                }
            )
            for row in subgroup_catalog[["distinct_numt_id", "median_mt_length_for_main_bp"]].itertuples(index=False):
                subgroup_length_rows.append(
                    {
                        "subgroup": subgroup,
                        "entity_id": row.distinct_numt_id,
                        "length_level": "distinct",
                        "length_bp": row.median_mt_length_for_main_bp,
                    }
                )
            for row in subgroup_events[["sample_event_id", "mt_length_for_main_bp"]].itertuples(index=False):
                subgroup_length_rows.append(
                    {
                        "subgroup": subgroup,
                        "entity_id": row.sample_event_id,
                        "length_level": "sample_event",
                        "length_bp": row.mt_length_for_main_bp,
                    }
                )
            frequency_rows = (
                subgroup_catalog.groupby("frequency_class", as_index=False)
                .agg(distinct_numt_n=("distinct_numt_id", "nunique"))
                .copy()
            )
            for row in frequency_rows.itertuples(index=False):
                subgroup_frequency_rows.append(
                    {
                        "subgroup": subgroup,
                        "frequency_class": row.frequency_class,
                        "distinct_numt_n": row.distinct_numt_n,
                    }
                )
            mt_tables = compute_mt_breakpoint_tables(subgroup_catalog) if not subgroup_catalog.empty else {}
            if mt_tables:
                mt_df = mt_tables["mt_breakpoint_cumulative"].copy()
                mt_df["subgroup"] = subgroup
                subgroup_mt_rows.append(mt_df)
            chr_distribution = (
                subgroup_catalog.groupby("chr", as_index=False)
                .agg(distinct_numt_n=("distinct_numt_id", "nunique"))
                .copy()
            )
            for row in chr_distribution.itertuples(index=False):
                subgroup_chr_rows.append(
                    {
                        "subgroup": subgroup,
                        "chr": row.chr,
                        "distinct_numt_n": row.distinct_numt_n,
                    }
                )
        _write_dimension_outputs(
            dim_dir,
            summary_rows,
            subgroup_metrics_rows,
            subgroup_length_rows,
            subgroup_frequency_rows,
            subgroup_mt_rows,
            subgroup_chr_rows,
        )

    LOG.info("Stratified results completed")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute stratified NUMT frequency statistics")
    parser.add_argument("--catalog-tsv", required=True)
    parser.add_argument("--sample-map-tsv", required=True)
    parser.add_argument("--sample-events-tsv", required=True)
    parser.add_argument("--sample-inclusion-tsv", required=True)
    parser.add_argument("--meta-tsv", required=True)
    parser.add_argument("--meta-id-col", required=True)
    parser.add_argument("--geography-col", required=True)
    parser.add_argument("--ethnicity-col", required=True)
    parser.add_argument("--sex-col", required=True)
    parser.add_argument("--prevalent-min", required=True, type=float)
    parser.add_argument("--common-min", required=True, type=float)
    parser.add_argument("--rare-min", required=True, type=float)
    parser.add_argument("--min-group-n-warn", required=True, type=int)
    parser.add_argument("--output-dir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
