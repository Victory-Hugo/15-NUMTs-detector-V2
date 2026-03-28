#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""基于 merge_bed 构建样本区域级 NUMT 长度表。"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
import re
import statistics
import zlib
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

LOG = logging.getLogger(__name__)
REGION_PATTERN = re.compile(r"^(?P<sample>.+)_(?P<chrom>chr[^.]+)\.(?P<start>\d+)\.(?P<end>\d+)$")


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="从过滤后的 merge_bed 生成区域级长度统计表。")
    parser.add_argument("--input-gz", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--chunk-rows", type=int, default=1_000_000)
    parser.add_argument("--temp-dir", required=True)
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


def first_pass(
    input_gz: str,
    temp_dir: str,
    chunk_rows: int,
    bucket_count: int = 64,
) -> Dict[str, Dict[str, object]]:
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
                start = int(fields[1])
                end = int(fields[2])
                length = end - start
                if length < 0:
                    raise ValueError(f"检测到负长度: {line[:200]}")

                bucket_id = bucket_index(region_id, bucket_count)
                bucket_handles[bucket_id].write(f"{region_id}\t{length}\n")

                stats = region_stats.get(region_id)
                if stats is None:
                    region_info = parse_region_id(region_id)
                    region_stats[region_id] = {
                        **region_info,
                        "region_id": region_id,
                        "n_hits": 1,
                        "min_length": length,
                        "max_length": length,
                    }
                else:
                    stats["n_hits"] = int(stats["n_hits"]) + 1
                    stats["min_length"] = min(int(stats["min_length"]), length)
                    stats["max_length"] = max(int(stats["max_length"]), length)

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
) -> List[int]:
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    chosen_lengths: List[int] = []

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
                "min_length",
                "median_length",
                "max_length",
                "chosen_length",
            ]
        )

        for bucket_file in sorted(Path(temp_dir).glob("bucket_*.tsv")):
            bucket_lengths: Dict[str, List[int]] = defaultdict(list)
            with bucket_file.open("r", encoding="utf-8", newline="") as bucket_handle:
                for line in bucket_handle:
                    region_id, length_text = line.rstrip("\n").split("\t")
                    bucket_lengths[region_id].append(int(length_text))

            for region_id in sorted(bucket_lengths):
                lengths = bucket_lengths[region_id]
                median_length = statistics.median(sorted(lengths))
                stats = region_stats[region_id]
                chosen_length = int(stats["max_length"])
                chosen_lengths.append(chosen_length)
                writer.writerow(
                    [
                        stats["sampleID"],
                        stats["region_chr"],
                        stats["region_start"],
                        stats["region_end"],
                        region_id,
                        stats["n_hits"],
                        stats["min_length"],
                        median_length,
                        stats["max_length"],
                        chosen_length,
                    ]
                )
            bucket_file.unlink()

    return chosen_lengths


def write_summary(summary_tsv: str, chosen_lengths: List[int], region_count: int) -> None:
    if not chosen_lengths:
        raise ValueError("未获得任何区域长度结果")
    rows = [
        ("region_count", region_count),
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
) -> int:
    region_stats = first_pass(input_gz=input_gz, temp_dir=temp_dir, chunk_rows=chunk_rows)
    chosen_lengths = second_pass_and_write(region_stats=region_stats, temp_dir=temp_dir, output_tsv=output_tsv)
    write_summary(summary_tsv=summary_tsv, chosen_lengths=chosen_lengths, region_count=len(region_stats))
    LOG.info("长度统计完成，输出文件: %s", output_tsv)
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
    )


if __name__ == "__main__":
    raise SystemExit(main())
