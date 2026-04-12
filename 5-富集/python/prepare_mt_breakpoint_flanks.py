#!/usr/bin/env python3
"""Prepare mtDNA breakpoint flank BED files and enrichment task manifest."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_breakpoint_flanks import (
    build_event_to_cluster,
    load_frequency_table,
    parse_classes,
    discover_target_beds,
    sanitize_label,
)

log = logging.getLogger(__name__)

MT_GENOME_LENGTH = 16569


def load_mt_breakpoint_candidates(
    breakpoint_tsv_gz: Path,
    event_map: pd.DataFrame,
    frequency_df: pd.DataFrame,
    chunksize: int,
) -> pd.DataFrame:
    # mt_confident rows: chr=="chrM", pos=mtDNA position;
    # source_region_chr/start/end = nuclear insertion site (used to link to clusters)
    usecols = [
        "sampleID",
        "Group",
        "readsCount",
        "pos",
        "source_region_chr",
        "source_region_start",
        "source_region_end",
    ]
    candidate_frames: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        breakpoint_tsv_gz,
        sep="\t",
        compression="infer",
        usecols=usecols,
        dtype={"sampleID": str, "Group": str, "source_region_chr": str},
        chunksize=chunksize,
    ):
        chunk = chunk.loc[chunk["Group"] == "mt_confident"].copy()
        if chunk.empty:
            continue
        for col in ["readsCount", "pos", "source_region_start", "source_region_end"]:
            chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
        chunk = chunk.dropna(subset=["pos", "source_region_start", "source_region_end"]).copy()
        chunk["source_region_start"] = chunk["source_region_start"].astype(int)
        chunk["source_region_end"] = chunk["source_region_end"].astype(int)
        chunk["mt_pos"] = chunk["pos"].round().astype(int)
        # event_key uses nuclear insertion coords to match cluster_input (same as nuclear version)
        chunk["event_key"] = (
            chunk["sampleID"].astype(str)
            + "|"
            + chunk["source_region_chr"].astype(str)
            + "|"
            + chunk["source_region_start"].astype(str)
            + "|"
            + chunk["source_region_end"].astype(str)
        )
        mapped = chunk.merge(event_map, on="event_key", how="inner")
        if mapped.empty:
            continue
        mapped = mapped.merge(frequency_df, left_on="mergedClusterID", right_on="merged_cluster_id", how="inner")
        if mapped.empty:
            continue
        mapped["readsCount"] = pd.to_numeric(mapped["readsCount"], errors="coerce").fillna(0).astype(int)
        candidate_frames.append(
            mapped[["merged_cluster_id", "frequency_class", "mt_pos", "readsCount"]]
        )

    if candidate_frames:
        return pd.concat(candidate_frames, ignore_index=True)
    return pd.DataFrame(columns=["merged_cluster_id", "frequency_class", "mt_pos", "readsCount"])


def collapse_mt_cluster_breakpoints(candidates: pd.DataFrame, flank_bp: int) -> pd.DataFrame:
    if candidates.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "breakpoint_id", "merged_cluster_id", "frequency_class"])
    candidates = candidates.sort_values(
        ["merged_cluster_id", "readsCount", "mt_pos"],
        ascending=[True, False, True],
    )
    selected = candidates.drop_duplicates("merged_cluster_id", keep="first").copy()
    selected["start"] = (selected["mt_pos"] - int(flank_bp)).clip(lower=0)
    selected["end"] = (selected["mt_pos"] + int(flank_bp)).clip(upper=MT_GENOME_LENGTH)
    selected["chr"] = "chrM"
    selected["breakpoint_id"] = selected["merged_cluster_id"].astype(str) + "__mt"
    bed = selected[["chr", "start", "end", "breakpoint_id", "merged_cluster_id", "frequency_class"]].copy()
    bed = bed.loc[bed["end"] > bed["start"]].copy()
    return bed


def run(
    breakpoint_tsv_gz: str | Path,
    cluster_input_tsv_gz: str | Path,
    cluster_detail_tsv: str | Path,
    frequency_cluster_tsv: str | Path,
    mt_regions_dir: str | Path,
    target_bed_glob: str,
    frequency_classes: str,
    flank_bp: int,
    prepared_dir: str | Path,
    enrichment_dir: str | Path,
    task_manifest: str | Path,
    chunksize: int = 500000,
) -> None:
    breakpoint_path = Path(breakpoint_tsv_gz)
    cluster_input_path = Path(cluster_input_tsv_gz)
    cluster_detail_path = Path(cluster_detail_tsv)
    frequency_path = Path(frequency_cluster_tsv)
    mt_regions_path = Path(mt_regions_dir)
    prepared_path = Path(prepared_dir)
    enrichment_path = Path(enrichment_dir)
    task_manifest_path = Path(task_manifest)

    for p in [breakpoint_path, cluster_input_path, cluster_detail_path, frequency_path]:
        if not p.is_file():
            raise FileNotFoundError(f"Input file not found: {p}")
    if not mt_regions_path.is_dir():
        raise NotADirectoryError(f"mt_regions_dir not found: {mt_regions_path}")

    prepared_path.mkdir(parents=True, exist_ok=True)
    enrichment_path.mkdir(parents=True, exist_ok=True)
    task_manifest_path.parent.mkdir(parents=True, exist_ok=True)

    requested_classes = parse_classes(frequency_classes)
    target_beds = discover_target_beds(mt_regions_path, target_bed_glob)
    event_map = build_event_to_cluster(cluster_input_path, cluster_detail_path)
    frequency_df = load_frequency_table(frequency_path)
    candidates = load_mt_breakpoint_candidates(breakpoint_path, event_map, frequency_df, int(chunksize))
    bed_df = collapse_mt_cluster_breakpoints(candidates, flank_bp)
    log.info("mtDNA: %d distinct cluster breakpoints after collapsing.", len(bed_df))

    manifest_rows: list[dict] = []
    summary_rows: list[dict] = []

    for frequency_class in requested_classes:
        if frequency_class == "all":
            class_df = bed_df.copy()
        else:
            class_df = bed_df.loc[bed_df["frequency_class"].astype(str) == frequency_class].copy()
        class_label = sanitize_label(frequency_class)
        bed_path = prepared_path / f"mt_breakpoint_flanks.{class_label}.bed"
        class_df.to_csv(bed_path, sep="\t", header=False, index=False)
        summary_rows.append(
            {
                "frequency_class": frequency_class,
                "breakpoint_count": int(len(class_df)),
                "cluster_count": int(class_df["merged_cluster_id"].nunique()) if not class_df.empty else 0,
                "bed_path": str(bed_path),
            }
        )
        for target_bed in target_beds:
            region_name = target_bed.name
            region_id = sanitize_label(region_name)
            task_id = f"mt_{region_id}__{class_label}"
            output_dir = enrichment_path / region_id / class_label
            manifest_rows.append(
                {
                    "task_id": task_id,
                    "region_id": region_id,
                    "region_name": region_name,
                    "target_bed": str(target_bed),
                    "frequency_class": frequency_class,
                    "breakpoint_bed": str(bed_path),
                    "output_dir": str(output_dir),
                }
            )

    prepared_summary = prepared_path / "prepared_breakpoint_flanks_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(prepared_summary, sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(task_manifest_path, sep="\t", index=False)
    log.info("Wrote %d tasks to %s", len(manifest_rows), task_manifest_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare mtDNA breakpoint flank BEDs and task manifest.")
    parser.add_argument("--breakpoint-tsv-gz", required=True)
    parser.add_argument("--cluster-input-tsv-gz", required=True)
    parser.add_argument("--cluster-detail-tsv", required=True)
    parser.add_argument("--frequency-cluster-tsv", required=True)
    parser.add_argument("--mt-regions-dir", required=True)
    parser.add_argument("--target-bed-glob", default="*.bed")
    parser.add_argument("--frequency-classes", default="all,common,low-frequency,rare,ultra-rare")
    parser.add_argument("--flank-bp", type=int, default=100)
    parser.add_argument("--prepared-dir", required=True)
    parser.add_argument("--enrichment-dir", required=True)
    parser.add_argument("--task-manifest", required=True)
    parser.add_argument("--chunksize", type=int, default=500000)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
