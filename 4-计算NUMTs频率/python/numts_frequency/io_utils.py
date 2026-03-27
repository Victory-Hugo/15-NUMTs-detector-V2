#!/usr/bin/env python3
"""输入输出与标准化读取。"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Iterable

import pandas as pd

from .parse_utils import build_region_key, normalize_chr_name, parse_length_region_key

LOG = logging.getLogger(__name__)
FAST_READ_ENGINE = "pyarrow"


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted((int(start), int(end)) for start, end in intervals)
    merged: list[tuple[int, int]] = []
    for start, end in ordered:
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
            continue
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def _primary_interval(intervals: list[tuple[int, int]]) -> tuple[int, int]:
    return min(intervals, key=lambda item: (-(item[1] - item[0] + 1), item[0], item[1]))


def _empty_length_summary() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "region_key",
            "sample_id",
            "chr",
            "start",
            "end",
            "mt_source_min",
            "mt_source_max",
            "mt_source_envelope_span_bp",
            "mt_source_span_bp",
            "length_bed_read_count",
            "merged_interval_count",
            "complex_mt_interval",
            "mt_primary_fragment_start",
            "mt_primary_fragment_end",
            "mt_primary_fragment_span_bp",
            "mt_merged_intervals",
        ]
    )


def _finalize_length_summary(summary_map: dict[str, dict]) -> pd.DataFrame:
    if not summary_map:
        return _empty_length_summary()
    records = []
    for value in summary_map.values():
        merged_intervals = _merge_intervals(value.pop("raw_intervals"))
        primary_start, primary_end = _primary_interval(merged_intervals)
        value["mt_source_envelope_span_bp"] = value["mt_source_max"] - value["mt_source_min"] + 1
        value["mt_source_span_bp"] = sum(end - start + 1 for start, end in merged_intervals)
        value["merged_interval_count"] = len(merged_intervals)
        value["complex_mt_interval"] = len(merged_intervals) > 1
        value["mt_primary_fragment_start"] = primary_start
        value["mt_primary_fragment_end"] = primary_end
        value["mt_primary_fragment_span_bp"] = primary_end - primary_start + 1
        value["mt_merged_intervals"] = ";".join(f"{start}-{end}" for start, end in merged_intervals)
        records.append(value)
    return pd.DataFrame.from_records(records)


def ensure_dir(path: str | Path) -> Path:
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return target


def read_delim(path: str | Path, sep: str = "\t", **kwargs) -> pd.DataFrame:
    """Prefer pyarrow for faster local TSV/CSV reads, fall back to pandas default engine."""
    try:
        return pd.read_csv(path, sep=sep, engine=FAST_READ_ENGINE, **kwargs)
    except Exception:
        return pd.read_csv(path, sep=sep, **kwargs)


def read_breakpoints(path: str) -> pd.DataFrame:
    header = read_delim(path, sep="\t", nrows=0)
    available_cols = set(header.columns.astype(str).tolist())
    usecols = ["sampleID", "chr", "pointGroup", "pos", "readsCount"]
    if "source_region_key" in available_cols:
        usecols.append("source_region_key")
    elif {"source_sample_id", "source_region_chr", "source_region_start", "source_region_end"}.issubset(available_cols):
        usecols.extend(["source_sample_id", "source_region_chr", "source_region_start", "source_region_end"])

    dtype_map = {"sampleID": str, "chr": str, "pointGroup": str}
    if "source_region_key" in usecols:
        dtype_map["source_region_key"] = str
    if "source_sample_id" in usecols:
        dtype_map["source_sample_id"] = str
    if "source_region_chr" in usecols:
        dtype_map["source_region_chr"] = str

    df = read_delim(
        path,
        sep="\t",
        usecols=usecols,
        dtype=dtype_map,
    )
    df = df.rename(columns={"sampleID": "sample_id"}).copy()
    df["chr"] = df["chr"].map(normalize_chr_name)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    if "readsCount" in df.columns:
        df["readsCount"] = pd.to_numeric(df["readsCount"], errors="coerce").fillna(1).astype(int)
    else:
        df["readsCount"] = 1
    if "source_region_chr" in df.columns:
        df["source_region_chr"] = df["source_region_chr"].map(normalize_chr_name)
        df["source_region_start"] = pd.to_numeric(df["source_region_start"], errors="coerce")
        df["source_region_end"] = pd.to_numeric(df["source_region_end"], errors="coerce")
    if "source_region_key" not in df.columns and {"source_sample_id", "source_region_chr", "source_region_start", "source_region_end"}.issubset(df.columns):
        df["source_region_key"] = [
            build_region_key(sample_id, chrom, start, end)
            if pd.notna(start) and pd.notna(end)
            else pd.NA
            for sample_id, chrom, start, end in zip(
                df["source_sample_id"],
                df["source_region_chr"],
                df["source_region_start"],
                df["source_region_end"],
            )
        ]
    if "source_region_key" in df.columns:
        df["source_region_key"] = df["source_region_key"].astype("string").str.strip()
        df.loc[df["source_region_key"].eq(""), "source_region_key"] = pd.NA
    return df


def read_clusters(path: str) -> pd.DataFrame:
    df = read_delim(
        path,
        sep="\t",
        usecols=["IndividualID", "chr", "Cluster_ID", "start", "end", "Cluster_No", "reads", "mt_positions"],
        dtype={"IndividualID": str, "chr": str, "Cluster_ID": str},
    )
    df = df.rename(columns={"IndividualID": "sample_id", "Cluster_ID": "cluster_id_raw"}).copy()
    df["chr"] = df["chr"].map(normalize_chr_name)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["Cluster_No"] = pd.to_numeric(df["Cluster_No"], errors="coerce")
    df["region_key"] = [
        build_region_key(sample_id, chrom, start, end)
        for sample_id, chrom, start, end in zip(df["sample_id"], df["chr"], df["start"], df["end"])
    ]
    return df


def read_meta(path: str, id_col: str) -> pd.DataFrame:
    df = read_delim(path, sep="\t", dtype={id_col: str})
    df = df.rename(columns={id_col: "sample_id"}).copy()
    return df


def summarize_length_bed(path: str, region_lookup: dict[str, dict] | None = None) -> pd.DataFrame:
    """把长度 bed 汇总成区域级表。"""
    summary_map: dict[str, dict] = {}
    allowed_region_keys = set(region_lookup.keys()) if region_lookup is not None else None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            region_key = parts[0].split("|", 1)[0]
            if allowed_region_keys is not None and region_key not in allowed_region_keys:
                continue
            mt_start = pd.to_numeric(parts[1], errors="coerce")
            mt_end = pd.to_numeric(parts[2], errors="coerce")
            if pd.isna(mt_start) or pd.isna(mt_end):
                continue
            if region_key not in summary_map:
                if region_lookup is None:
                    info = parse_length_region_key(parts[0])
                else:
                    info = region_lookup[region_key]
                summary_map[region_key] = {
                    "region_key": region_key,
                    "sample_id": info["sample_id"],
                    "chr": info["chr"],
                    "start": info["start"],
                    "end": info["end"],
                    "mt_source_min": mt_start,
                    "mt_source_max": mt_end,
                    "length_bed_read_count": 1,
                    "raw_intervals": [(int(mt_start), int(mt_end))],
                }
                continue
            summary_map[region_key]["mt_source_min"] = min(summary_map[region_key]["mt_source_min"], mt_start)
            summary_map[region_key]["mt_source_max"] = max(summary_map[region_key]["mt_source_max"], mt_end)
            summary_map[region_key]["length_bed_read_count"] += 1
            summary_map[region_key]["raw_intervals"].append((int(mt_start), int(mt_end)))
    summary = _finalize_length_summary(summary_map)
    LOG.info("Loaded %d summarized region rows from length bed", len(summary))
    return summary


def subset_length_summary(length_summary: pd.DataFrame, region_lookup: dict[str, dict]) -> pd.DataFrame:
    if length_summary.empty:
        return _empty_length_summary()
    region_keys = set(region_lookup.keys())
    subset = length_summary[length_summary["region_key"].astype(str).isin(region_keys)].copy()
    if subset.empty:
        return _empty_length_summary()
    return subset


def read_length_summary(path: str) -> pd.DataFrame:
    return read_delim(
        path,
        sep="\t",
        dtype={
            "sample_id": str,
            "chr": str,
            "region_key": str,
            "mt_merged_intervals": str,
        },
    )


def read_mt_gene_table(path: str) -> pd.DataFrame:
    df = read_delim(path, sep=",")
    df["POSSTART"] = pd.to_numeric(df["POSSTART"], errors="coerce").astype(int)
    df["POSEND"] = pd.to_numeric(df["POSEND"], errors="coerce").astype(int)
    return df


def write_tsv(df: pd.DataFrame, path: str) -> None:
    ensure_dir(Path(path).parent)
    df.to_csv(path, sep="\t", index=False)
