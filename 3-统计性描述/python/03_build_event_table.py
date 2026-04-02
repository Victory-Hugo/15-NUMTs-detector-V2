#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从过滤后的 breakpoint 表构建样本级事件表。"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from collections import Counter
from pathlib import Path
from typing import Iterable, List, Tuple

LOG = logging.getLogger(__name__)


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="从 breakpoint 表去重生成样本级 NUMT 事件表。")
    parser.add_argument("--input-gz", required=True)
    parser.add_argument("--events-out", required=True)
    parser.add_argument("--ideogram-out", required=True)
    parser.add_argument("--chr-count-out", required=True)
    parser.add_argument("--chunk-rows", type=int, default=1_000_000)
    return parser


def chrom_sort_key(chrom: str) -> Tuple[int, str]:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return (int(suffix), chrom)
        if suffix == "X":
            return (23, chrom)
        if suffix == "Y":
            return (24, chrom)
    return (10_000, chrom)


def run(
    input_gz: str,
    events_out: str,
    ideogram_out: str,
    chr_count_out: str,
    chunk_rows: int,
) -> int:
    unique_events = set()
    total_rows = 0

    with gzip.open(input_gz, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_cols = {"sampleID", "source_region_chr", "source_region_start", "source_region_end"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"breakpoint 缺少必要列: {sorted(missing)}")

        for row in reader:
            total_rows += 1
            key = (
                row["sampleID"],
                row["source_region_chr"],
                int(row["source_region_start"]),
                int(row["source_region_end"]),
            )
            unique_events.add(key)
            if total_rows % chunk_rows == 0:
                LOG.info("event 去重已处理 %s 行，当前唯一事件数 %s", total_rows, len(unique_events))

    sorted_events = sorted(unique_events, key=lambda item: (chrom_sort_key(item[1]), item[2], item[3], item[0]))
    chr_counter = Counter(event[1] for event in sorted_events)

    Path(events_out).parent.mkdir(parents=True, exist_ok=True)
    with Path(events_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sampleID", "region_chr", "region_start", "region_end", "midpoint", "region_id"])
        for sample_id, chrom, start, end in sorted_events:
            midpoint = (start + end) // 2
            region_id = f"{sample_id}_{chrom}.{start}.{end}"
            writer.writerow([sample_id, chrom, start, end, midpoint, region_id])

    with Path(ideogram_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Chr", "Start", "End"])
        for _, chrom, start, end in sorted_events:
            writer.writerow([chrom, start, end])

    with Path(chr_count_out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chr", "event_count"])
        for chrom in sorted(chr_counter, key=chrom_sort_key):
            writer.writerow([chrom, chr_counter[chrom]])

    LOG.info("事件表构建完成，总断点行数 %s，唯一样本级事件数 %s", total_rows, len(sorted_events))
    return 0


def main(argv: List[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        input_gz=args.input_gz,
        events_out=args.events_out,
        ideogram_out=args.ideogram_out,
        chr_count_out=args.chr_count_out,
        chunk_rows=args.chunk_rows,
    )


if __name__ == "__main__":
    raise SystemExit(main())
