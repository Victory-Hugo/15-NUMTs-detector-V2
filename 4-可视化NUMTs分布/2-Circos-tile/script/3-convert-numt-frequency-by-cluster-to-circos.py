#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Convert distinct NUMT cluster summary to circos text with real mtDNA ranges."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ["chr", "cluster_min_midpoint", "cluster_max_midpoint"]
RAW_CLUSTER_REQUIRED_COLUMNS = ["Cluster_ID"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert NUMT cluster frequency TSV to six-column circos text format."
    )
    parser.add_argument("--input", required=True, help="Input 4-numt-frequency-by-cluster.tsv")
    parser.add_argument("--output", required=True, help="Output circos_SLE.txt")
    parser.add_argument("--raw-cluster-file", required=True, help="Input 7-mt_disc_cluster.tsv.gz")
    parser.add_argument("--mt-chr", default="hsM", help="Mitochondrial chromosome label. Default: hsM")
    parser.add_argument(
        "--dedup-window-bp",
        type=int,
        default=1000,
        help="Remove a row when both nuclear Start and End are within this bp distance of the previous kept row on the same chromosome. Default: 1000",
    )
    return parser


def normalize_chr(value: str) -> str:
    value = str(value).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return f"hs{value}"


def parse_cluster_id(cluster_id: str) -> tuple[str, int, int, int, int]:
    parts = str(cluster_id).split("_")
    if len(parts) < 6:
        raise ValueError(f"Unexpected Cluster_ID format: {cluster_id}")
    nu_chr = parts[0]
    nu_start = int(parts[1])
    nu_end = int(parts[2])
    mt_start = int(parts[4])
    mt_end = int(parts[5])
    return nu_chr, nu_start, nu_end, mt_start, mt_end


def load_raw_cluster_table(raw_cluster_file: str) -> pd.DataFrame:
    raw_df = pd.read_csv(raw_cluster_file, sep="\t", usecols=RAW_CLUSTER_REQUIRED_COLUMNS)
    parsed = raw_df["Cluster_ID"].map(parse_cluster_id)
    parsed_df = pd.DataFrame(
        parsed.tolist(),
        columns=["chr", "nu_start", "nu_end", "mt_start", "mt_end"],
    )
    parsed_df["chr"] = parsed_df["chr"].astype(str)
    parsed_df["nu_midpoint"] = ((parsed_df["nu_start"] + parsed_df["nu_end"]) / 2).astype(int)
    return parsed_df


def build_raw_cluster_index(raw_cluster_df: pd.DataFrame) -> dict[str, dict[str, np.ndarray]]:
    index: dict[str, dict[str, np.ndarray]] = {}
    for chrom, chrom_df in raw_cluster_df.groupby("chr", sort=False):
        chrom_df = chrom_df.sort_values("nu_midpoint", kind="stable").reset_index(drop=True)
        index[str(chrom)] = {
            "midpoints": chrom_df["nu_midpoint"].to_numpy(dtype=int),
            "mt_start": chrom_df["mt_start"].to_numpy(dtype=int),
            "mt_end": chrom_df["mt_end"].to_numpy(dtype=int),
        }
    return index


def deduplicate_nearby_rows(out_df: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    if out_df.empty or window_bp < 0:
        return out_df

    kept_rows: list[dict[str, object]] = []
    for chrom, chrom_df in out_df.groupby("Chr", sort=False):
        last_start: int | None = None
        last_end: int | None = None
        for row in chrom_df.itertuples(index=False):
            start = int(row.Start)
            end = int(row.End)
            if last_start is not None and abs(start - last_start) <= window_bp and abs(end - last_end) <= window_bp:
                continue
            kept_rows.append(
                {
                    "Chr": chrom,
                    "Start": start,
                    "End": end,
                    "MtChr": row.MtChr,
                    "MtStart": int(row.MtStart),
                    "MtEnd": int(row.MtEnd),
                }
            )
            last_start = start
            last_end = end

    return pd.DataFrame(kept_rows, columns=["Chr", "Start", "End", "MtChr", "MtStart", "MtEnd"])


def main() -> None:
    args = build_parser().parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", dtype={"chr": str})
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    raw_cluster_df = load_raw_cluster_table(args.raw_cluster_file)
    raw_cluster_index = build_raw_cluster_index(raw_cluster_df)
    out_rows: list[dict[str, object]] = []

    for row in df.itertuples(index=False):
        chrom = str(row.chr)
        start = int(row.cluster_min_midpoint)
        end = int(row.cluster_max_midpoint)
        chrom_index = raw_cluster_index.get(chrom)
        if chrom_index is None:
            raise ValueError(f"No raw chromosome records matched distinct locus: {chrom}:{start}-{end}")

        midpoints = chrom_index["midpoints"]
        left = int(np.searchsorted(midpoints, start, side="left"))
        right = int(np.searchsorted(midpoints, end, side="right"))
        if left >= right:
            raise ValueError(f"No raw mtDNA coordinates matched distinct locus: {chrom}:{start}-{end}")

        mt_start_values = chrom_index["mt_start"][left:right]
        mt_end_values = chrom_index["mt_end"][left:right]

        out_rows.append(
            {
                "Chr": normalize_chr(chrom),
                "Start": start,
                "End": end,
                "MtChr": args.mt_chr,
                "MtStart": int(mt_start_values.min()),
                "MtEnd": int(mt_end_values.max()),
            }
        )

    out_df = pd.DataFrame(out_rows, columns=["Chr", "Start", "End", "MtChr", "MtStart", "MtEnd"])

    out_df = out_df.sort_values(["Chr", "Start", "End"], kind="stable")
    out_df = deduplicate_nearby_rows(out_df, window_bp=args.dedup_window_bp)
    out_df.to_csv(output_path, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
