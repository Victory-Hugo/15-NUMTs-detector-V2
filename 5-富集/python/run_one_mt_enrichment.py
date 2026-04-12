#!/usr/bin/env python3
"""Run one mtDNA target-region by breakpoint-stratum enrichment task."""

from __future__ import annotations

import argparse
import hashlib
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def stable_seed(base_seed: int, task_id: str) -> int:
    digest = hashlib.sha256(f"{base_seed}:{task_id}".encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def read_bed3(path: str | Path) -> pd.DataFrame:
    bed_path = Path(path)
    if not bed_path.is_file():
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    df = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    df["chr"] = df["chr"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df.loc[df["end"] > df["start"]].copy()


def target_to_mt_linear(target_bed: str | Path) -> pd.DataFrame:
    """Convert mtDNA target BED to linear coordinate DataFrame (offset=0, chrM only)."""
    df = read_bed3(target_bed)
    df = df.loc[df["chr"] == "chrM"].copy()
    if df.empty:
        raise ValueError(f"No chrM intervals found in target BED: {target_bed}")
    return pd.DataFrame({"Start": df["start"], "End": df["end"]})


def merge_intervals(interval_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    sorted_df = interval_df.sort_values(["Start", "End"]).reset_index(drop=True)
    merged: list[tuple[int, int]] = []
    current_start = int(sorted_df.loc[0, "Start"])
    current_end = int(sorted_df.loc[0, "End"])
    for row in sorted_df.itertuples(index=False):
        start = int(row.Start)
        end = int(row.End)
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start = start
            current_end = end
    merged.append((current_start, current_end))
    starts = np.array([item[0] for item in merged], dtype=np.int64)
    ends = np.array([item[1] for item in merged], dtype=np.int64)
    return starts, ends


def count_vectorized_overlaps(
    random_start: np.ndarray, random_end: np.ndarray, starts: np.ndarray, ends: np.ndarray
) -> int:
    idx = np.searchsorted(starts, random_end, side="left") - 1
    valid = idx >= 0
    if not np.any(valid):
        return 0
    overlap = np.zeros(len(random_start), dtype=bool)
    overlap[valid] = ends[idx[valid]] > random_start[valid]
    return int(overlap.sum())


def count_observed_overlap(breakpoint_bed: str | Path, target_bed: str | Path, bedtools_bin: str) -> int:
    cmd = [bedtools_bin, "intersect", "-a", str(breakpoint_bed), "-b", str(target_bed), "-u"]
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    if not result.stdout:
        return 0
    return len(result.stdout.rstrip("\n").splitlines())


def run(
    task_id: str,
    region_id: str,
    region_name: str,
    target_bed: str,
    frequency_class: str,
    breakpoint_bed: str,
    output_dir: str,
    bedtools_bin: str,
    simulation_runs: int,
    random_seed: int,
    genome_length: int,
    flank_bp: int,
) -> dict[str, str | int | float]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    log_path = output_path / "task.log"
    handler = logging.FileHandler(log_path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    breakpoint_df = read_bed3(breakpoint_bed)
    breakpoints_total = int(len(breakpoint_df))
    if breakpoints_total <= 0:
        raise ValueError(f"No breakpoint flanks found for task {task_id}: {breakpoint_bed}")

    flank_length = int(flank_bp) * 2
    observed_overlap = count_observed_overlap(
        breakpoint_bed=breakpoint_bed,
        target_bed=target_bed,
        bedtools_bin=bedtools_bin,
    )
    observed_percentage = observed_overlap / breakpoints_total

    target_starts, target_ends = merge_intervals(target_to_mt_linear(target_bed))
    rng = np.random.default_rng(stable_seed(int(random_seed), task_id))
    detail_rows: list[dict[str, int | float]] = []

    max_start = int(genome_length) - flank_length - 1
    if max_start <= 0:
        raise ValueError(f"Genome length {genome_length} is too small for flank length {flank_length}")

    replace = breakpoints_total >= max_start
    for iteration in range(1, int(simulation_runs) + 1):
        random_start = rng.choice(max_start, size=breakpoints_total, replace=replace) + 1
        random_end = random_start + flank_length
        overlap_count = count_vectorized_overlaps(random_start, random_end, target_starts, target_ends)
        detail_rows.append(
            {
                "iteration": iteration,
                "overlap_count": int(overlap_count),
                "percentage": overlap_count / breakpoints_total,
            }
        )

    detail_df = pd.DataFrame(detail_rows)
    p_value_less = float((detail_df["percentage"] <= observed_percentage).mean())
    p_value_greater = float((detail_df["percentage"] >= observed_percentage).mean())

    detail_file = output_path / "simulation_detail.tsv"
    summary_file = output_path / "summary.tsv"
    detail_df.to_csv(detail_file, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "task_id": task_id,
                "region_id": region_id,
                "region_name": region_name,
                "target_bed": target_bed,
                "frequency_class": frequency_class,
                "breakpoint_bed": breakpoint_bed,
                "breakpoints_total": breakpoints_total,
                "breakpoints_target": int(observed_overlap),
                "flank_bp": int(flank_bp),
                "flank_length": int(flank_length),
                "simulation_runs": int(simulation_runs),
                "observed_percentage": float(observed_percentage),
                "p_value_less": p_value_less,
                "p_value_greater": p_value_greater,
                "significance_less_than_random": "Yes" if p_value_less < 0.05 else "No",
                "significance_greater_than_random": "Yes" if p_value_greater < 0.05 else "No",
                "detail_file": str(detail_file),
            }
        ]
    ).to_csv(summary_file, sep="\t", index=False)

    log.info(
        "Finished %s: class=%s region=%s observed=%s/%s p_less=%.6g p_greater=%.6g",
        task_id,
        frequency_class,
        region_name,
        observed_overlap,
        breakpoints_total,
        p_value_less,
        p_value_greater,
    )
    log.removeHandler(handler)
    handler.close()
    return {
        "summary_file": str(summary_file),
        "detail_file": str(detail_file),
        "p_value_less": p_value_less,
        "p_value_greater": p_value_greater,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run one mtDNA breakpoint flank enrichment task.")
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--region-id", required=True)
    parser.add_argument("--region-name", required=True)
    parser.add_argument("--target-bed", required=True)
    parser.add_argument("--frequency-class", required=True)
    parser.add_argument("--breakpoint-bed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bedtools-bin", required=True)
    parser.add_argument("--simulation-runs", required=True, type=int)
    parser.add_argument("--random-seed", required=True, type=int)
    parser.add_argument("--genome-length", required=True, type=int)
    parser.add_argument("--flank-bp", required=True, type=int)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
