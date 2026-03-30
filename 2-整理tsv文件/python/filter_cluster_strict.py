#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple

# 上游默认使用 discordant cluster 表中的宽松阈值作为过滤标准。
# 此处可以针对已经生成的 tsv.gz 文件进行过滤，保留 Cluster_No >= 5 的记录，生成2套文件。
# 生成/覆盖的输出文件（位于 output_dir）
# - 7-mt_disc_cluster.tsv.gz
# - 4-mt_disc_breakpoint_input.tsv.gz
# - 6-mt_disc_cluster_summary.tsv.gz
# - 5-mt_disc_cluster_old.tsv.gz
# - 1-all_breakpoints.tsv.gz
# - 3-confident_breakpoints.tsv.gz
# - 2-breakpoints_old.tsv.gz
# - 8-merge_bed.tsv.gz

# 返回值
# - int: 程序退出码，成功时返回 0。

# 依赖（调用的外部函数/模块）
# - parse_args, require_inputs, build_kept_region_keys, filter_cluster_table, filter_breakpoint_input, filter_cluster_summary, filter_cluster_old, filter_breakpoints_by_region_key, rebuild_breakpoints_old_from_confident, filter_merge_bed, print_stats

REQUIRED_INPUTS = (
    "1-all_breakpoints.tsv.gz",
    "2-breakpoints_old.tsv.gz",
    "3-confident_breakpoints.tsv.gz",
    "4-mt_disc_breakpoint_input.tsv.gz",
    "5-mt_disc_cluster_old.tsv.gz",
    "6-mt_disc_cluster_summary.tsv.gz",
    "7-mt_disc_cluster.tsv.gz",
    "8-merge_bed.tsv.gz",
)


@dataclass
class FileStats:
    path: Path
    total_rows: int = 0
    kept_rows: int = 0
    dropped_rows: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter NUMTs summary tables by Cluster_No threshold."
    )
    parser.add_argument("--input-dir", required=True, help="Directory with the 8 input TSV.GZ files.")
    parser.add_argument("--output-dir", required=True, help="Directory for filtered TSV.GZ outputs.")
    parser.add_argument(
        "--min-cluster-no",
        type=int,
        default=5,
        help="Minimum Cluster_No to keep. Default: 5",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs in the output directory.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    print(message, file=sys.stderr)


def open_gzip_text(path: Path, mode: str):
    kwargs = {"encoding": "utf-8", "newline": ""}
    if "w" in mode:
        kwargs["compresslevel"] = 1
    return gzip.open(path, mode, **kwargs)


def ensure_output_path(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists, use --overwrite to replace: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def atomic_writer(path: Path, overwrite: bool):
    ensure_output_path(path, overwrite)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    return tmp_path


def finalize_atomic_write(tmp_path: Path, final_path: Path) -> None:
    os.replace(tmp_path, final_path)


def require_inputs(input_dir: Path) -> Dict[str, Path]:
    paths: Dict[str, Path] = {}
    missing: List[str] = []
    for name in REQUIRED_INPUTS:
        path = input_dir / name
        if path.exists():
            paths[name] = path
        else:
            missing.append(str(path))
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    return paths


def split_tsv(line: str) -> List[str]:
    return line.rstrip("\n").split("\t")


def region_key(sample_id: str, chrom: str, start: str, end: str) -> str:
    return f"{sample_id}_{chrom}.{start}.{end}"


def build_header_index(header_line: str) -> Dict[str, int]:
    fields = split_tsv(header_line)
    return {name: idx for idx, name in enumerate(fields)}


def build_kept_region_keys(cluster_tsv_gz: Path, min_cluster_no: int) -> Tuple[Set[str], Set[str], FileStats]:
    kept: Set[str] = set()
    dropped: Set[str] = set()
    stats = FileStats(path=cluster_tsv_gz)

    with open_gzip_text(cluster_tsv_gz, "rt") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty cluster table: {cluster_tsv_gz}")
        header_idx = build_header_index(header)
        required = ("chr", "start", "end", "Cluster_No", "IndividualID")
        for field in required:
            if field not in header_idx:
                raise ValueError(f"Missing column {field} in {cluster_tsv_gz}")

        for line_no, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            try:
                cluster_no = int(parts[header_idx["Cluster_No"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid Cluster_No at {cluster_tsv_gz}:{line_no}") from exc

            key = region_key(
                parts[header_idx["IndividualID"]],
                parts[header_idx["chr"]],
                parts[header_idx["start"]],
                parts[header_idx["end"]],
            )
            if cluster_no >= min_cluster_no:
                kept.add(key)
                stats.kept_rows += 1
            else:
                dropped.add(key)
                stats.dropped_rows += 1

    return kept, dropped, stats


def filter_cluster_table(
    in_gz: Path,
    out_gz: Path,
    kept_region_keys: Set[str],
    min_cluster_no: int,
    overwrite: bool,
) -> FileStats:
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"Empty input: {in_gz}")
        dst.write(header)
        idx = build_header_index(header)
        for field in ("chr", "start", "end", "Cluster_No", "IndividualID"):
            if field not in idx:
                raise ValueError(f"Missing column {field} in {in_gz}")

        for line_no, line in enumerate(src, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            try:
                cluster_no = int(parts[idx["Cluster_No"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid Cluster_No at {in_gz}:{line_no}") from exc
            key = region_key(
                parts[idx["IndividualID"]],
                parts[idx["chr"]],
                parts[idx["start"]],
                parts[idx["end"]],
            )
            if key in kept_region_keys and cluster_no >= min_cluster_no:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_breakpoint_input(
    in_gz: Path,
    out_gz: Path,
    kept_region_keys: Set[str],
    min_cluster_no: int,
    overwrite: bool,
) -> Tuple[FileStats, int]:
    stats = FileStats(path=out_gz)
    mismatch_count = 0
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            if len(parts) < 8:
                raise ValueError(f"Expected 8 columns at {in_gz}:{line_no}")
            sample_id, cluster_no_raw, _, _, _, chrom, start, end = parts[:8]
            try:
                cluster_no = int(cluster_no_raw)
            except ValueError as exc:
                raise ValueError(f"Invalid Cluster_No at {in_gz}:{line_no}") from exc
            key = region_key(sample_id, chrom, start, end)
            keep_by_key = key in kept_region_keys
            keep_by_value = cluster_no >= min_cluster_no
            if keep_by_key != keep_by_value:
                mismatch_count += 1
            if keep_by_key:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats, mismatch_count


def filter_cluster_summary(
    in_gz: Path,
    out_gz: Path,
    min_cluster_no: int,
    overwrite: bool,
) -> FileStats:
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            if len(parts) < 6:
                raise ValueError(f"Expected 6 columns at {in_gz}:{line_no}")
            try:
                cluster_no = int(parts[3])
            except ValueError as exc:
                raise ValueError(f"Invalid Cluster_No at {in_gz}:{line_no}") from exc
            if cluster_no >= min_cluster_no:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_cluster_old(
    in_gz: Path,
    out_gz: Path,
    min_cluster_no: int,
    overwrite: bool,
) -> FileStats:
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"Empty input: {in_gz}")
        dst.write(header)
        idx = build_header_index(header)
        if "Cluster_No" not in idx:
            raise ValueError(f"Missing Cluster_No column in {in_gz}")

        for line_no, line in enumerate(src, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            try:
                cluster_no = int(parts[idx["Cluster_No"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid Cluster_No at {in_gz}:{line_no}") from exc
            if cluster_no >= min_cluster_no:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_breakpoints_by_region_key(
    in_gz: Path,
    out_gz: Path,
    kept_region_keys: Set[str],
    overwrite: bool,
) -> FileStats:
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"Empty input: {in_gz}")
        dst.write(header)
        idx = build_header_index(header)
        if "source_region_key" not in idx:
            raise ValueError(f"Missing source_region_key in {in_gz}")

        for line_no, line in enumerate(src, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            try:
                source_key = parts[idx["source_region_key"]]
            except IndexError as exc:
                raise ValueError(f"Missing source_region_key at {in_gz}:{line_no}") from exc
            if source_key in kept_region_keys:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_merge_bed(
    in_gz: Path,
    out_gz: Path,
    kept_region_keys: Set[str],
    overwrite: bool,
) -> FileStats:
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            if len(parts) < 4:
                raise ValueError(f"Expected 4 columns at {in_gz}:{line_no}")
            query_name = parts[0]
            if "|" not in query_name:
                raise ValueError(f"Missing region delimiter in query_name at {in_gz}:{line_no}")
            key = query_name.split("|", 1)[0]
            if key in kept_region_keys:
                dst.write(line)
                stats.kept_rows += 1
            else:
                stats.dropped_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def rebuild_breakpoints_old_from_confident(
    filtered_confident_gz: Path,
    out_gz: Path,
    overwrite: bool,
) -> FileStats:
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(filtered_confident_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"Empty input: {filtered_confident_gz}")
        idx = build_header_index(header)
        required = ("sampleID", "chr", "strand", "pointGroup", "Group", "readsCount", "pos")
        for field in required:
            if field not in idx:
                raise ValueError(f"Missing column {field} in {filtered_confident_gz}")

        for line_no, line in enumerate(src, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)
            point_group = parts[idx["pointGroup"]]
            sample_id = parts[idx["sampleID"]]
            chrom = parts[idx["chr"]]
            strand = parts[idx["strand"]]
            reads_count = parts[idx["readsCount"]]
            try:
                pos_int = str(int(float(parts[idx["pos"]])))
            except ValueError as exc:
                raise ValueError(f"Invalid pos at {filtered_confident_gz}:{line_no}") from exc

            if point_group == "nu_Tend_Bleft":
                out_parts = [point_group, "nuleft", chrom, pos_int, strand, reads_count, "-1", sample_id]
            elif point_group == "nu_Tstart_Bright":
                out_parts = [point_group, "nuright", chrom, "-1", strand, reads_count, pos_int, sample_id]
            elif point_group == "mt_Tend":
                out_parts = [point_group, "mtRight", "MT", pos_int, strand, reads_count, "-1", sample_id]
            elif point_group == "mt_Tstart":
                out_parts = [point_group, "mtRight", "MT", "-1", strand, reads_count, pos_int, sample_id]
            else:
                raise ValueError(f"Unsupported pointGroup {point_group} at {filtered_confident_gz}:{line_no}")

            dst.write("\t".join(out_parts) + "\n")
            stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def print_stats(label: str, stats: FileStats) -> None:
    print(
        f"{label}\tinput_rows={stats.total_rows}\tkept_rows={stats.kept_rows}"
        f"\tdropped_rows={stats.dropped_rows}\toutput={stats.path}"
    )


def main() -> int:

    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    min_cluster_no = args.min_cluster_no

    if min_cluster_no < 1:
        raise ValueError("--min-cluster-no must be >= 1")

    paths = require_inputs(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    kept_keys, dropped_keys, cluster_key_stats = build_kept_region_keys(
        paths["7-mt_disc_cluster.tsv.gz"], min_cluster_no
    )

    print(
        f"kept_region_keys={len(kept_keys)}\tdropped_region_keys={len(dropped_keys)}"
        f"\tmin_cluster_no={min_cluster_no}"
    )

    stats_cluster = filter_cluster_table(
        paths["7-mt_disc_cluster.tsv.gz"],
        output_dir / "7-mt_disc_cluster.tsv.gz",
        kept_keys,
        min_cluster_no,
        args.overwrite,
    )
    stats_break_input, mismatch_count = filter_breakpoint_input(
        paths["4-mt_disc_breakpoint_input.tsv.gz"],
        output_dir / "4-mt_disc_breakpoint_input.tsv.gz",
        kept_keys,
        min_cluster_no,
        args.overwrite,
    )
    stats_cluster_summary = filter_cluster_summary(
        paths["6-mt_disc_cluster_summary.tsv.gz"],
        output_dir / "6-mt_disc_cluster_summary.tsv.gz",
        min_cluster_no,
        args.overwrite,
    )
    stats_cluster_old = filter_cluster_old(
        paths["5-mt_disc_cluster_old.tsv.gz"],
        output_dir / "5-mt_disc_cluster_old.tsv.gz",
        min_cluster_no,
        args.overwrite,
    )
    stats_all = filter_breakpoints_by_region_key(
        paths["1-all_breakpoints.tsv.gz"],
        output_dir / "1-all_breakpoints.tsv.gz",
        kept_keys,
        args.overwrite,
    )
    stats_confident = filter_breakpoints_by_region_key(
        paths["3-confident_breakpoints.tsv.gz"],
        output_dir / "3-confident_breakpoints.tsv.gz",
        kept_keys,
        args.overwrite,
    )
    stats_break_old = rebuild_breakpoints_old_from_confident(
        output_dir / "3-confident_breakpoints.tsv.gz",
        output_dir / "2-breakpoints_old.tsv.gz",
        args.overwrite,
    )
    stats_merge_bed = filter_merge_bed(
        paths["8-merge_bed.tsv.gz"],
        output_dir / "8-merge_bed.tsv.gz",
        kept_keys,
        args.overwrite,
    )

    print_stats("7-mt_disc_cluster.tsv.gz", stats_cluster)
    print_stats("4-mt_disc_breakpoint_input.tsv.gz", stats_break_input)
    print(f"4-mt_disc_breakpoint_input.tsv.gz\tmismatch_count={mismatch_count}")
    print_stats("6-mt_disc_cluster_summary.tsv.gz", stats_cluster_summary)
    print_stats("5-mt_disc_cluster_old.tsv.gz", stats_cluster_old)
    print_stats("1-all_breakpoints.tsv.gz", stats_all)
    print_stats("3-confident_breakpoints.tsv.gz", stats_confident)
    print_stats("2-breakpoints_old.tsv.gz", stats_break_old)
    print_stats("8-merge_bed.tsv.gz", stats_merge_bed)
    print_stats("region_key_source(7)", cluster_key_stats)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        log(f"ERROR: {exc}")
        raise
