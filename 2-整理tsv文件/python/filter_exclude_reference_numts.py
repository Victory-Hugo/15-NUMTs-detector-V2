#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 过滤排除参考NUMTs：对已经过严格Cluster_No阈值过滤的8个文件，
# 进一步排除核侧坐标与参考NUMTs重叠（上下游各扩展extension bp）的记录。

import argparse
import bisect
import gzip
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

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
        description="Filter NUMTs data by excluding reference NUMTs overlaps."
    )
    parser.add_argument("--input-dir", required=True, help="Directory with the 8 input TSV.GZ files.")
    parser.add_argument("--output-dir", required=True, help="Directory for filtered TSV.GZ outputs.")
    parser.add_argument("--ref-file", required=True, help="Reference NUMTs TSV file.")
    parser.add_argument(
        "--extension", type=int, default=500,
        help="Extension in bp around nuclear coordinates for overlap check. Default: 500"
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
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


def atomic_writer(path: Path, overwrite: bool) -> Path:
    ensure_output_path(path, overwrite)
    return path.with_suffix(path.suffix + ".tmp")


def finalize_atomic_write(tmp_path: Path, final_path: Path) -> None:
    os.replace(tmp_path, final_path)


def split_tsv(line: str) -> List[str]:
    return line.rstrip("\n").split("\t")


def build_header_index(header_line: str) -> Dict[str, int]:
    fields = split_tsv(header_line)
    return {name: idx for idx, name in enumerate(fields)}


def normalize_chr(chr_val: str) -> str:
    """统一去除chr前缀：'chr1' -> '1', '1' -> '1'"""
    val = chr_val.strip()
    if val.startswith("chr") or val.startswith("Chr"):
        val = val[3:]
    return val


def load_reference_numts(ref_file: Path) -> Dict[str, Tuple[List[int], List[int], List[int]]]:
    """加载参考NUMTs，按染色体分组，并预处理为二分查找结构。

    返回：{chr_norm: (starts, ends, cummax_ends)}
    - starts: 按start排序的起始坐标列表
    - ends: 对应的终止坐标列表
    - cummax_ends[i] = max(ends[0..i])，用于快速判断最大覆盖范围

    重叠判断可降至 O(log n)：
    若 cummax_ends[bisect_right(starts, qe)-1] >= qs，则存在重叠。
    """
    raw: Dict[str, List[Tuple[int, int]]] = {}
    with open(ref_file, "r", encoding="utf-8") as f:
        header = f.readline()
        idx = build_header_index(header)
        for required in ("Chr", "Nuc Start", "Nuc End"):
            if required not in idx:
                raise ValueError(f"Missing column '{required}' in {ref_file}")

        for line in f:
            parts = split_tsv(line)
            chr_val = parts[idx["Chr"]].strip()
            start_val = parts[idx["Nuc Start"]].strip()
            end_val = parts[idx["Nuc End"]].strip()

            if not chr_val or not start_val or not end_val:
                continue

            try:
                start = int(start_val)
                end = int(end_val)
            except ValueError:
                continue

            chr_norm = normalize_chr(chr_val)
            raw.setdefault(chr_norm, []).append((start, end))

    # 构建排序+累积最大值索引
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]] = {}
    for chr_norm, intervals in raw.items():
        intervals.sort(key=lambda x: x[0])
        starts = [iv[0] for iv in intervals]
        ends = [iv[1] for iv in intervals]
        cummax: List[int] = []
        cur_max = -1
        for e in ends:
            cur_max = max(cur_max, e)
            cummax.append(cur_max)
        ref_numts[chr_norm] = (starts, ends, cummax)

    total = sum(len(v[0]) for v in ref_numts.values())
    log(f"Loaded {total} reference NUMTs across {len(ref_numts)} chromosomes from {ref_file}")
    return ref_numts


def overlaps_reference(
    chr_val: str,
    start: int,
    end: int,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    extension: int,
) -> bool:
    """O(log n) 判断 [start-extension, end+extension] 是否与参考NUMTs中任一区间重叠。

    算法：
    1. 扩展查询区间 [qs, qe]
    2. 用 bisect_right(starts, qe) 找到所有 start <= qe 的区间（前 i 个）
    3. 若这 i 个区间中的最大 end（cummax_ends[i-1]）>= qs，则存在重叠
    """
    chr_norm = normalize_chr(chr_val)
    if chr_norm not in ref_numts:
        return False

    starts, ends, cummax = ref_numts[chr_norm]
    qs = max(0, start - extension)
    qe = end + extension

    i = bisect.bisect_right(starts, qe)
    if i == 0:
        return False
    return cummax[i - 1] >= qs


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


def print_stats(label: str, stats: FileStats) -> None:
    print(
        f"{label}\ttotal={stats.total_rows}\tkept={stats.kept_rows}"
        f"\tdropped={stats.dropped_rows}\toutput={stats.path}"
    )


# --- File filter functions ---

def filter_with_header_cols(
    in_gz: Path,
    out_gz: Path,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    chr_col: str,
    start_col: str,
    end_col: str,
    extension: int,
    overwrite: bool,
) -> FileStats:
    """通用：有表头，使用指定列名过滤。
    用于文件1、3（source_region_chr/start/end）和文件7（chr/start/end）。
    """
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"Empty input: {in_gz}")
        dst.write(header)
        idx = build_header_index(header)

        for col in (chr_col, start_col, end_col):
            if col not in idx:
                raise ValueError(f"Missing column '{col}' in {in_gz}")

        for line_no, line in enumerate(src, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)

            try:
                chr_val = parts[idx[chr_col]]
                start = int(parts[idx[start_col]])
                end = int(parts[idx[end_col]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Parse error at {in_gz}:{line_no}") from exc

            if overlaps_reference(chr_val, start, end, ref_numts, extension):
                stats.dropped_rows += 1
            else:
                dst.write(line)
                stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_breakpoints_old(
    in_gz: Path,
    out_gz: Path,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    extension: int,
    overwrite: bool,
) -> FileStats:
    """文件2: 2-breakpoints_old.tsv.gz
    无表头，列：point_group | direction | chr | pos_left | strand | reads | pos_right | sample_id
    - nu_* 行：核侧chr=col[2]，坐标取col[3]和col[6]中不为"-1"的值（单点）
    - mt_* 行：线粒体坐标，直接保留（不做参考NUMTs比较）
    """
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)

            if len(parts) < 7:
                dst.write(line)
                stats.kept_rows += 1
                continue

            point_group = parts[0]
            if point_group.startswith("mt_"):
                dst.write(line)
                stats.kept_rows += 1
                continue

            chr_val = parts[2]
            pos_left = parts[3]
            pos_right = parts[6]

            try:
                pos = int(pos_left) if pos_left != "-1" else int(pos_right)
            except ValueError:
                dst.write(line)
                stats.kept_rows += 1
                continue

            if overlaps_reference(chr_val, pos, pos, ref_numts, extension):
                stats.dropped_rows += 1
            else:
                dst.write(line)
                stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_breakpoint_input(
    in_gz: Path,
    out_gz: Path,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    extension: int,
    overwrite: bool,
) -> FileStats:
    """文件4: 4-mt_disc_breakpoint_input.tsv.gz
    无表头，列：sample_id | cluster_no | mt_disc_sam | mt_split_sam | cram | chr | start | end | ...
    核侧坐标：col[5]=chr, col[6]=start, col[7]=end
    """
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)

            if len(parts) < 8:
                raise ValueError(f"Expected >= 8 columns at {in_gz}:{line_no}, got {len(parts)}")

            try:
                chr_val = parts[5]
                start = int(parts[6])
                end = int(parts[7])
            except ValueError as exc:
                raise ValueError(f"Parse error at {in_gz}:{line_no}") from exc

            if overlaps_reference(chr_val, start, end, ref_numts, extension):
                stats.dropped_rows += 1
            else:
                dst.write(line)
                stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_cluster_old(
    in_gz: Path,
    out_gz: Path,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    extension: int,
    overwrite: bool,
) -> FileStats:
    """文件5: 5-mt_disc_cluster_old.tsv.gz
    有表头（SAM扩展格式）。RNAME=核侧染色体，POS=核侧位置（单点）。
    """
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        header = src.readline()
        if not header:
            raise ValueError(f"Empty input: {in_gz}")
        dst.write(header)
        idx = build_header_index(header)

        for col in ("RNAME", "POS"):
            if col not in idx:
                raise ValueError(f"Missing column '{col}' in {in_gz}")

        for line_no, line in enumerate(src, start=2):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)

            try:
                chr_val = parts[idx["RNAME"]]
                pos = int(parts[idx["POS"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Parse error at {in_gz}:{line_no}") from exc

            if overlaps_reference(chr_val, pos, pos, ref_numts, extension):
                stats.dropped_rows += 1
            else:
                dst.write(line)
                stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_cluster_summary(
    in_gz: Path,
    out_gz: Path,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    extension: int,
    overwrite: bool,
) -> FileStats:
    """文件6: 6-mt_disc_cluster_summary.tsv.gz
    无表头，第3列（index=2）复合字段：chr11_49861843_49862261_chrMboth_16089_16477
    核侧坐标：tokens[0]=chr, tokens[1]=start, tokens[2]=end
    """
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)

            if len(parts) < 3:
                raise ValueError(f"Expected >= 3 columns at {in_gz}:{line_no}")

            try:
                composite = parts[2]
                tokens = composite.split("_")
                if len(tokens) < 3:
                    raise ValueError(f"Composite field too short: {composite!r}")
                chr_val = tokens[0]
                start = int(tokens[1])
                end = int(tokens[2])
            except ValueError as exc:
                raise ValueError(f"Parse error at {in_gz}:{line_no}: {exc}") from exc

            if overlaps_reference(chr_val, start, end, ref_numts, extension):
                stats.dropped_rows += 1
            else:
                dst.write(line)
                stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


def filter_merge_bed(
    in_gz: Path,
    out_gz: Path,
    ref_numts: Dict[str, Tuple[List[int], List[int], List[int]]],
    extension: int,
    overwrite: bool,
) -> FileStats:
    """文件8: 8-merge_bed.tsv.gz
    无表头，第1列格式：sample_chr.start.end|read_name
    例：1000434_chr1.118141557.118142709|E100028702L1...
    核侧坐标：取最后一个'_'后的部分 chr.start.end，按'.'分割解析。
    """
    stats = FileStats(path=out_gz)
    tmp_path = atomic_writer(out_gz, overwrite)

    with open_gzip_text(in_gz, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
        for line_no, line in enumerate(src, start=1):
            if not line.strip():
                continue
            stats.total_rows += 1
            parts = split_tsv(line)

            try:
                query_name = parts[0]
                region_key = query_name.split("|", 1)[0]
                # region_key: "1000434_chr1.118141557.118142709"
                region_part = region_key.rsplit("_", 1)[-1]
                # region_part: "chr1.118141557.118142709"
                region_tokens = region_part.split(".")
                if len(region_tokens) < 3:
                    raise ValueError(f"Cannot parse region from: {region_part!r}")
                chr_val = region_tokens[0]
                start = int(region_tokens[1])
                end = int(region_tokens[2])
            except ValueError as exc:
                raise ValueError(f"Parse error at {in_gz}:{line_no}: {exc}") from exc

            if overlaps_reference(chr_val, start, end, ref_numts, extension):
                stats.dropped_rows += 1
            else:
                dst.write(line)
                stats.kept_rows += 1

    finalize_atomic_write(tmp_path, out_gz)
    return stats


# --- Main ---

def main() -> int:
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    ref_file = Path(args.ref_file).resolve()
    extension = args.extension

    if not ref_file.exists():
        raise FileNotFoundError(f"Reference NUMTs file not found: {ref_file}")

    ref_numts = load_reference_numts(ref_file)
    paths = require_inputs(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"Input dir:  {input_dir}")
    log(f"Output dir: {output_dir}")
    log(f"Extension:  {extension} bp")

    log("Processing 1-all_breakpoints.tsv.gz...")
    stats_1 = filter_with_header_cols(
        paths["1-all_breakpoints.tsv.gz"],
        output_dir / "1-all_breakpoints.tsv.gz",
        ref_numts, "source_region_chr", "source_region_start", "source_region_end",
        extension, args.overwrite,
    )

    log("Processing 2-breakpoints_old.tsv.gz...")
    stats_2 = filter_breakpoints_old(
        paths["2-breakpoints_old.tsv.gz"],
        output_dir / "2-breakpoints_old.tsv.gz",
        ref_numts, extension, args.overwrite,
    )

    log("Processing 3-confident_breakpoints.tsv.gz...")
    stats_3 = filter_with_header_cols(
        paths["3-confident_breakpoints.tsv.gz"],
        output_dir / "3-confident_breakpoints.tsv.gz",
        ref_numts, "source_region_chr", "source_region_start", "source_region_end",
        extension, args.overwrite,
    )

    log("Processing 4-mt_disc_breakpoint_input.tsv.gz...")
    stats_4 = filter_breakpoint_input(
        paths["4-mt_disc_breakpoint_input.tsv.gz"],
        output_dir / "4-mt_disc_breakpoint_input.tsv.gz",
        ref_numts, extension, args.overwrite,
    )

    log("Processing 5-mt_disc_cluster_old.tsv.gz...")
    stats_5 = filter_cluster_old(
        paths["5-mt_disc_cluster_old.tsv.gz"],
        output_dir / "5-mt_disc_cluster_old.tsv.gz",
        ref_numts, extension, args.overwrite,
    )

    log("Processing 6-mt_disc_cluster_summary.tsv.gz...")
    stats_6 = filter_cluster_summary(
        paths["6-mt_disc_cluster_summary.tsv.gz"],
        output_dir / "6-mt_disc_cluster_summary.tsv.gz",
        ref_numts, extension, args.overwrite,
    )

    log("Processing 7-mt_disc_cluster.tsv.gz...")
    stats_7 = filter_with_header_cols(
        paths["7-mt_disc_cluster.tsv.gz"],
        output_dir / "7-mt_disc_cluster.tsv.gz",
        ref_numts, "chr", "start", "end",
        extension, args.overwrite,
    )

    log("Processing 8-merge_bed.tsv.gz...")
    stats_8 = filter_merge_bed(
        paths["8-merge_bed.tsv.gz"],
        output_dir / "8-merge_bed.tsv.gz",
        ref_numts, extension, args.overwrite,
    )

    print()
    print("=== Statistics ===")
    print_stats("1-all_breakpoints.tsv.gz", stats_1)
    print_stats("2-breakpoints_old.tsv.gz", stats_2)
    print_stats("3-confident_breakpoints.tsv.gz", stats_3)
    print_stats("4-mt_disc_breakpoint_input.tsv.gz", stats_4)
    print_stats("5-mt_disc_cluster_old.tsv.gz", stats_5)
    print_stats("6-mt_disc_cluster_summary.tsv.gz", stats_6)
    print_stats("7-mt_disc_cluster.tsv.gz", stats_7)
    print_stats("8-merge_bed.tsv.gz", stats_8)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
