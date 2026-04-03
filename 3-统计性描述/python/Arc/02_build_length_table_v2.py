#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""基于 merge_bed 构建样本区域级 NUMT 长度表（v2：chosen_length 改为 mt_start 位置跨度）。

v1 问题：chosen_length = max(mt_end - mt_start)，即单条 read 最长比对长度，
对于 ~800 bp 的 NUMT，每条 read 只跨越断点附近 ~150 bp，导致严重低估。

v2 修复：
  - 收集每个 source_region_key 下所有 read 的 mt_start 位置
  - 计算 mt_start 的跨度（max_mt_start - min_mt_start），即 reads 铺开的 mtDNA 范围
  - 对环形 mtDNA（reads 同时出现在序列头部和尾部）做环形坐标矫正
  - chosen_length = mt_start_span（矫正后）
  - 同时保留原有的单 read 比对统计列（min/median/max read span）供参考
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
import statistics
import zlib
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

LOG = logging.getLogger(__name__)
REGION_PATTERN = re.compile(r"^(?P<sample>.+)_(?P<chrom>chr[^.]+)\.(?P<start>\d+)\.(?P<end>\d+)$")

# 环形检测阈值：若同一 region 的 mt_start 既有 < MT_CIRCULAR_LOW 又有 > MT_CIRCULAR_HIGH，
# 则判定为跨环形边界，需做坐标变换。
# 人类 mtDNA = 16569 bp，以中点 8284 为界分高低两半。
MT_CIRCULAR_LOW = 1000
MT_CIRCULAR_HIGH = 15569


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="从过滤后的 merge_bed 生成区域级长度统计表（v2）。")
    parser.add_argument("--input-gz", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--chunk-rows", type=int, default=1_000_000)
    parser.add_argument("--temp-dir", required=True)
    parser.add_argument(
        "--mt-genome-length",
        type=int,
        default=16569,
        help="mtDNA 全长（bp），用于环形坐标矫正，默认 16569（人类 mtDNA）",
    )
    return parser


def parse_region_id(region_id: str) -> Dict[str, str]:
    match = REGION_PATTERN.match(region_id)
    if match is None:
        raise ValueError(f"无法解析 region_id: {region_id}")
    info = match.groupdict()
    return {
        "sampleID": info["sample"],
        "region_chr": info["chrom"],
        "region_start": info["start"],
        "region_end": info["end"],
    }


def bucket_index(region_id: str, bucket_count: int) -> int:
    return zlib.crc32(region_id.encode("utf-8")) % bucket_count


def compute_mt_start_span(mt_starts: List[int], mt_genome_length: int) -> Tuple[int, bool]:
    """计算 mt_start 位置跨度，自动处理环形 mtDNA。

    返回：(span_bp, is_circular)
      - span_bp：矫正后的跨度（bp），代表该 region 对应的 NUMT 在 mtDNA 上的覆盖范围
      - is_circular：是否检测到环形跨越

    环形检测逻辑：
      若 mt_start 中既有 < MT_CIRCULAR_LOW 又有 > MT_CIRCULAR_HIGH，
      则判定为跨 mtDNA 环形边界（reads 同时覆盖序列起点和终点附近），
      对 > 半长 的值做 -= mt_genome_length 变换后再计算跨度。
    """
    if not mt_starts:
        return 0, False

    if len(mt_starts) == 1:
        return 0, False

    has_low = any(s < MT_CIRCULAR_LOW for s in mt_starts)
    has_high = any(s > MT_CIRCULAR_HIGH for s in mt_starts)
    is_circular = has_low and has_high

    if is_circular:
        # 将高端位置变换为负值，使所有点落在同一线性区间
        half = mt_genome_length // 2
        transformed = [s - mt_genome_length if s > half else s for s in mt_starts]
        span = max(transformed) - min(transformed)
    else:
        span = max(mt_starts) - min(mt_starts)

    return max(0, span), is_circular


def first_pass(
    input_gz: str,
    temp_dir: str,
    chunk_rows: int,
    bucket_count: int = 64,
) -> Dict[str, Dict[str, object]]:
    """第一遍：流式读取 merge_bed，分桶写出 (region_id, mt_start, read_span)。

    同时在内存维护每个 region 的基础统计（n_hits、read span 的 min/max）。
    """
    temp_path = Path(temp_dir)
    temp_path.mkdir(parents=True, exist_ok=True)
    for old_bucket in temp_path.glob("bucket_*.tsv"):
        old_bucket.unlink()

    bucket_handles = []
    for bucket_id in range(bucket_count):
        bucket_file = temp_path / f"bucket_{bucket_id:03d}.tsv"
        bucket_handles.append(bucket_file.open("w", encoding="utf-8", newline=""))

    region_stats: Dict[str, Dict[str, object]] = {}
    total_rows = 0

    try:
        with gzip.open(input_gz, "rt", encoding="utf-8", newline="") as handle:
            for line in handle:
                total_rows += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(f"merge_bed 行列数不足: {line[:200]}")
                record_key = fields[0]
                region_id = record_key.split("|", 1)[0]
                mt_start = int(fields[1])
                mt_end = int(fields[2])
                read_span = mt_end - mt_start
                if read_span < 0:
                    raise ValueError(f"检测到负 read_span: {line[:200]}")

                bucket_id = bucket_index(region_id, bucket_count)
                # 桶文件格式：region_id \t mt_start \t read_span
                bucket_handles[bucket_id].write(f"{region_id}\t{mt_start}\t{read_span}\n")

                stats = region_stats.get(region_id)
                if stats is None:
                    region_info = parse_region_id(region_id)
                    region_stats[region_id] = {
                        **region_info,
                        "region_id": region_id,
                        "n_hits": 1,
                        "min_read_span": read_span,
                        "max_read_span": read_span,
                    }
                else:
                    stats["n_hits"] = int(stats["n_hits"]) + 1
                    stats["min_read_span"] = min(int(stats["min_read_span"]), read_span)
                    stats["max_read_span"] = max(int(stats["max_read_span"]), read_span)

                if total_rows % chunk_rows == 0:
                    LOG.info("first pass 已处理 %s 行，唯一 region 数 %s", total_rows, len(region_stats))
    finally:
        for bucket_handle in bucket_handles:
            bucket_handle.close()

    LOG.info("first pass 完成，总行数 %s，唯一 region 数 %s", total_rows, len(region_stats))
    return region_stats


def second_pass_and_write(
    region_stats: Dict[str, Dict[str, object]],
    temp_dir: str,
    output_tsv: str,
    mt_genome_length: int,
) -> Tuple[List[int], int]:
    """第二遍：从桶文件计算 mt_start 跨度，写出完整长度表。

    返回：(chosen_lengths, circular_region_count)
    """
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    chosen_lengths: List[int] = []
    circular_count = 0

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sampleID",
                "region_chr",
                "region_start",
                "region_end",
                "region_id",
                "n_hits",
                # 单 read 比对跨度统计（保留供参考，不再作为 chosen_length）
                "min_read_span",
                "median_read_span",
                "max_read_span",
                # v2 新增：mt_start 位置跨度统计（核心修复）
                "mt_start_span",       # mt_start 位置跨度（环形矫正后），即 NUMT 在 mtDNA 上的覆盖范围
                "is_circular",         # 是否检测到环形跨越（True/False）
                # chosen_length = mt_start_span，与 v1 接口保持一致
                "chosen_length",
            ]
        )

        for bucket_file in sorted(Path(temp_dir).glob("bucket_*.tsv")):
            # 从桶文件收集：mt_start 列表 和 read_span 列表（按 region 分组）
            bucket_mt_starts: Dict[str, List[int]] = defaultdict(list)
            bucket_read_spans: Dict[str, List[int]] = defaultdict(list)

            with bucket_file.open("r", encoding="utf-8", newline="") as bucket_handle:
                for line in bucket_handle:
                    parts = line.rstrip("\n").split("\t")
                    region_id = parts[0]
                    mt_start = int(parts[1])
                    read_span = int(parts[2])
                    bucket_mt_starts[region_id].append(mt_start)
                    bucket_read_spans[region_id].append(read_span)

            for region_id in sorted(bucket_mt_starts):
                mt_starts = bucket_mt_starts[region_id]
                read_spans = bucket_read_spans[region_id]

                median_read_span = int(statistics.median(sorted(read_spans)))
                mt_start_span, is_circular = compute_mt_start_span(mt_starts, mt_genome_length)

                if is_circular:
                    circular_count += 1

                chosen_length = mt_start_span
                chosen_lengths.append(chosen_length)

                stats = region_stats[region_id]
                writer.writerow(
                    [
                        stats["sampleID"],
                        stats["region_chr"],
                        stats["region_start"],
                        stats["region_end"],
                        region_id,
                        stats["n_hits"],
                        stats["min_read_span"],
                        median_read_span,
                        stats["max_read_span"],
                        mt_start_span,
                        is_circular,
                        chosen_length,
                    ]
                )
            bucket_file.unlink()

    return chosen_lengths, circular_count


def write_summary(
    summary_tsv: str,
    chosen_lengths: List[int],
    region_count: int,
    circular_count: int,
) -> None:
    if not chosen_lengths:
        raise ValueError("未获得任何区域长度结果")
    rows = [
        ("region_count", region_count),
        ("circular_region_count", circular_count),
        # chosen_length = mt_start_span
        ("chosen_length_min", min(chosen_lengths)),
        ("chosen_length_max", max(chosen_lengths)),
        ("chosen_length_mean", f"{statistics.mean(chosen_lengths):.6f}"),
        ("chosen_length_median", statistics.median(chosen_lengths)),
    ]
    with Path(summary_tsv).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def run(
    input_gz: str,
    output_tsv: str,
    summary_tsv: str,
    chunk_rows: int,
    temp_dir: str,
    mt_genome_length: int = 16569,
) -> int:
    region_stats = first_pass(input_gz=input_gz, temp_dir=temp_dir, chunk_rows=chunk_rows)
    chosen_lengths, circular_count = second_pass_and_write(
        region_stats=region_stats,
        temp_dir=temp_dir,
        output_tsv=output_tsv,
        mt_genome_length=mt_genome_length,
    )
    write_summary(
        summary_tsv=summary_tsv,
        chosen_lengths=chosen_lengths,
        region_count=len(region_stats),
        circular_count=circular_count,
    )
    LOG.info(
        "长度统计完成（v2），输出文件: %s，环形 region 数: %s/%s",
        output_tsv,
        circular_count,
        len(region_stats),
    )
    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        input_gz=args.input_gz,
        output_tsv=args.output_tsv,
        summary_tsv=args.summary_tsv,
        chunk_rows=args.chunk_rows,
        temp_dir=args.temp_dir,
        mt_genome_length=args.mt_genome_length,
    )


if __name__ == "__main__":
    raise SystemExit(main())
