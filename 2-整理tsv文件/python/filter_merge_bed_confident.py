#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Set


@dataclass
class FileStats:
    confident_keys: int = 0
    input_rows: int = 0
    kept_rows: int = 0
    dropped_rows: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter 8-merge_bed.tsv.gz by source_region_key values from 3-confident_breakpoints.tsv.gz."
    )
    parser.add_argument("--confident-gz", required=True, help="Path to 3-confident_breakpoints.tsv.gz")
    parser.add_argument("--merge-bed-gz", required=True, help="Path to 8-merge_bed.tsv.gz")
    parser.add_argument("--output-gz", required=True, help="Path to output 9-merge_bed_confident.tsv.gz")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    return parser.parse_args()


def split_tsv(line: str):
    return line.rstrip("\n").split("\t")


def open_gzip_text(path: Path, mode: str):
    kwargs = {"encoding": "utf-8", "newline": ""}
    if "w" in mode:
        kwargs["compresslevel"] = 1
    return gzip.open(path, mode, **kwargs)


def build_header_index(header_line: str) -> Dict[str, int]:
    fields = split_tsv(header_line)
    return {name: idx for idx, name in enumerate(fields)}


def ensure_output_path(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists, use --overwrite to replace: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def atomic_tmp_path(path: Path, overwrite: bool) -> Path:
    ensure_output_path(path, overwrite)
    return path.with_suffix(path.suffix + ".tmp")


def finalize_atomic_write(tmp_path: Path, final_path: Path) -> None:
    os.replace(tmp_path, final_path)


def load_confident_region_keys(confident_gz: Path) -> Set[str]:
    region_keys: Set[str] = set()
    with open_gzip_text(confident_gz, "rt") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty input: {confident_gz}")
        header_idx = build_header_index(header)
        if "source_region_key" not in header_idx:
            raise ValueError(f"Missing column source_region_key in {confident_gz}")

        key_idx = header_idx["source_region_key"]
        for line_no, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = split_tsv(line)
            try:
                region_key = parts[key_idx]
            except IndexError as exc:
                raise ValueError(f"Missing source_region_key at {confident_gz}:{line_no}") from exc
            if not region_key:
                raise ValueError(f"Empty source_region_key at {confident_gz}:{line_no}")
            region_keys.add(region_key)
    return region_keys


def filter_merge_bed(
    merge_bed_gz: Path,
    output_gz: Path,
    confident_region_keys: Set[str],
    overwrite: bool,
) -> FileStats:
    stats = FileStats(confident_keys=len(confident_region_keys))
    tmp_path = atomic_tmp_path(output_gz, overwrite)

    with open_gzip_text(merge_bed_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.input_rows += 1
            parts = split_tsv(line)
            if len(parts) < 1:
                raise ValueError(f"Expected at least 1 column at {merge_bed_gz}:{line_no}")

            query_name = parts[0]
            if "|" not in query_name:
                raise ValueError(f"Missing region delimiter in query_name at {merge_bed_gz}:{line_no}")

            region_key = query_name.split("|", 1)[0]
            if region_key in confident_region_keys:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, output_gz)
    return stats


def main() -> int:
    args = parse_args()
    confident_gz = Path(args.confident_gz).resolve()
    merge_bed_gz = Path(args.merge_bed_gz).resolve()
    output_gz = Path(args.output_gz).resolve()

    if not confident_gz.exists():
        raise FileNotFoundError(f"Confident breakpoint file not found: {confident_gz}")
    if not merge_bed_gz.exists():
        raise FileNotFoundError(f"Merge bed file not found: {merge_bed_gz}")

    confident_region_keys = load_confident_region_keys(confident_gz)
    stats = filter_merge_bed(merge_bed_gz, output_gz, confident_region_keys, args.overwrite)

    print(
        f"confident_keys={stats.confident_keys}\t"
        f"input_rows={stats.input_rows}\t"
        f"kept_rows={stats.kept_rows}\t"
        f"dropped_rows={stats.dropped_rows}\t"
        f"output={output_gz}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
