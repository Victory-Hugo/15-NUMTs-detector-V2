#!/usr/bin/env python3
"""Prepare distinct NUMT BED files and shell task manifests.

Example:
    python prepare_frequency_numts.py --frequency-cluster-tsv input.tsv \
      --regions-dir data/regions --prepared-dir output/1-prepared-numts \
      --enrichment-dir output/2-enrichment --task-manifest tasks.tsv
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

REQUIRED_COLUMNS = {
    "merged_cluster_id",
    "chr",
    "cluster_min_midpoint",
    "cluster_max_midpoint",
    "sample_count",
    "chosen_region_length",
    "frequency",
    "frequency_class",
}


def sanitize_label(value: str) -> str:
    label = re.sub(r"\.bed$", "", str(value))
    label = re.sub(r"[^0-9A-Za-z._\-\u4e00-\u9fff]+", "_", label)
    return label.strip("._-") or "item"


def parse_classes(raw_value: str) -> list[str]:
    classes = [item.strip() for item in raw_value.split(",") if item.strip()]
    if not classes:
        raise ValueError("--frequency-classes must not be empty.")
    return classes


def build_bed(df: pd.DataFrame) -> pd.DataFrame:
    working = df.copy()
    for column in ["cluster_min_midpoint", "cluster_max_midpoint", "sample_count", "chosen_region_length", "frequency"]:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    working = working.dropna(
        subset=["chr", "cluster_min_midpoint", "cluster_max_midpoint", "chosen_region_length", "sample_count"]
    ).copy()
    if working.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "merged_cluster_id", "sample_count", "frequency_class"])

    start = (working["cluster_min_midpoint"] - working["chosen_region_length"] / 2).astype(int)
    end = (working["cluster_max_midpoint"] + working["chosen_region_length"] / 2).astype(int)
    start = start.mask(start < 0, 0)
    end = end.mask(end <= start, start + 1)
    bed = pd.DataFrame(
        {
            "chr": working["chr"].astype(str),
            "start": start.astype(int),
            "end": end.astype(int),
            "merged_cluster_id": working["merged_cluster_id"].astype(str),
            "sample_count": working["sample_count"].astype(int),
            "frequency_class": working["frequency_class"].astype(str),
            "frequency": working["frequency"],
            "chosen_region_length": working["chosen_region_length"].astype(int),
        }
    )
    return bed


def discover_target_beds(regions_dir: Path, target_bed_glob: str) -> list[Path]:
    beds = sorted(path for path in regions_dir.glob(target_bed_glob) if path.is_file())
    if not beds:
        raise FileNotFoundError(f"No target BED files matched {regions_dir / target_bed_glob}")
    return beds


def run(
    frequency_cluster_tsv: str | Path,
    regions_dir: str | Path,
    target_bed_glob: str,
    frequency_classes: str,
    prepared_dir: str | Path,
    enrichment_dir: str | Path,
    task_manifest: str | Path,
) -> dict[str, str | int]:
    input_path = Path(frequency_cluster_tsv)
    regions_path = Path(regions_dir)
    prepared_path = Path(prepared_dir)
    enrichment_path = Path(enrichment_dir)
    task_manifest_path = Path(task_manifest)

    if not input_path.is_file():
        raise FileNotFoundError(f"Frequency cluster table not found: {input_path}")
    if not regions_path.is_dir():
        raise NotADirectoryError(f"Target BED directory not found: {regions_path}")

    prepared_path.mkdir(parents=True, exist_ok=True)
    enrichment_path.mkdir(parents=True, exist_ok=True)
    task_manifest_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t")
    missing = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {input_path}: {', '.join(missing)}")

    requested_classes = parse_classes(frequency_classes)
    target_beds = discover_target_beds(regions_path, target_bed_glob)

    manifest_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for frequency_class in requested_classes:
        if frequency_class == "all":
            class_df = df.copy()
        else:
            class_df = df.loc[df["frequency_class"].astype(str) == frequency_class].copy()
        bed_df = build_bed(class_df)
        class_label = sanitize_label(frequency_class)
        bed_path = prepared_path / f"distinct_numts.{class_label}.bed"
        bed_df.to_csv(bed_path, sep="\t", header=False, index=False)
        summary_rows.append(
            {
                "frequency_class": frequency_class,
                "numt_count": int(len(bed_df)),
                "mean_numt_length": float((bed_df["end"] - bed_df["start"]).mean()) if not bed_df.empty else 0.0,
                "bed_path": str(bed_path),
            }
        )

        for target_bed in target_beds:
            region_name = target_bed.name
            region_id = sanitize_label(region_name)
            task_id = f"{region_id}__{class_label}"
            output_dir = enrichment_path / region_id / class_label
            manifest_rows.append(
                {
                    "task_id": task_id,
                    "region_id": region_id,
                    "region_name": region_name,
                    "target_bed": str(target_bed),
                    "frequency_class": frequency_class,
                    "numt_bed": str(bed_path),
                    "output_dir": str(output_dir),
                }
            )

    prepared_summary = prepared_path / "prepared_numts_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(prepared_summary, sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(task_manifest_path, sep="\t", index=False)
    log.info("Prepared %d NUMT strata and %d enrichment tasks.", len(summary_rows), len(manifest_rows))
    log.info("Prepared summary: %s", prepared_summary)
    log.info("Task manifest: %s", task_manifest_path)
    return {
        "prepared_summary": str(prepared_summary),
        "task_manifest": str(task_manifest_path),
        "task_count": int(len(manifest_rows)),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare distinct NUMT BED files and enrichment task manifest.")
    parser.add_argument("--frequency-cluster-tsv", required=True)
    parser.add_argument("--regions-dir", required=True)
    parser.add_argument("--target-bed-glob", default="*.bed")
    parser.add_argument("--frequency-classes", default="all,common,low-frequency,rare,ultra-rare")
    parser.add_argument("--prepared-dir", required=True)
    parser.add_argument("--enrichment-dir", required=True)
    parser.add_argument("--task-manifest", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
