#!/usr/bin/env python3
"""Prepare paper-aligned NUMT nuclear breakpoint flank BED files."""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

CLUSTER_INPUT_COLUMNS = ["sampleID", "Cluster_No", "disFile", "splitFile", "wgsBAM", "chr", "start", "end"]
PRIMARY_CHROMS = {f"chr{idx}" for idx in range(1, 23)} | {"chrX", "chrY"}


def sanitize_label(value: str) -> str:
    label = re.sub(r"\.bed$", "", str(value))
    label = re.sub(r"[^0-9A-Za-z._\-\u4e00-\u9fff]+", "_", label)
    return label.strip("._-") or "item"


def parse_classes(raw_value: str) -> list[str]:
    classes = [item.strip() for item in raw_value.split(",") if item.strip()]
    if not classes:
        raise ValueError("--frequency-classes must not be empty.")
    return classes


def discover_target_beds(regions_dir: Path, target_bed_glob: str) -> list[Path]:
    beds = sorted(path for path in regions_dir.glob(target_bed_glob) if path.is_file())
    if not beds:
        raise FileNotFoundError(f"No target BED files matched {regions_dir / target_bed_glob}")
    return beds


def build_event_to_cluster(cluster_input_tsv_gz: Path, cluster_detail_tsv: Path) -> pd.DataFrame:
    cluster_input = pd.read_csv(
        cluster_input_tsv_gz,
        sep="\t",
        header=None,
        names=CLUSTER_INPUT_COLUMNS,
        compression="infer",
        dtype={"sampleID": str, "chr": str},
    )
    cluster_input["start"] = pd.to_numeric(cluster_input["start"], errors="coerce")
    cluster_input["end"] = pd.to_numeric(cluster_input["end"], errors="coerce")
    cluster_input = cluster_input.dropna(subset=["sampleID", "chr", "start", "end"]).copy()
    cluster_input["start"] = cluster_input["start"].astype(int)
    cluster_input["end"] = cluster_input["end"].astype(int)
    cluster_input["POS"] = ((cluster_input["start"] + cluster_input["end"]) / 2).astype(int)
    cluster_input["event_key"] = (
        cluster_input["sampleID"].astype(str)
        + "|"
        + cluster_input["chr"].astype(str)
        + "|"
        + cluster_input["start"].astype(str)
        + "|"
        + cluster_input["end"].astype(str)
    )

    cluster_detail = pd.read_csv(cluster_detail_tsv, sep="\t", dtype={"CHR": str})
    required = {"POS", "CHR", "mergedClusterID"}
    missing = required - set(cluster_detail.columns)
    if missing:
        raise ValueError(f"Cluster detail table missing columns: {sorted(missing)}")
    cluster_detail["POS"] = pd.to_numeric(cluster_detail["POS"], errors="coerce")
    cluster_detail = cluster_detail.dropna(subset=["POS", "CHR", "mergedClusterID"]).copy()
    cluster_detail["POS"] = cluster_detail["POS"].astype(int)
    cluster_detail = cluster_detail.drop_duplicates(["CHR", "POS"])

    event_map = cluster_input.merge(
        cluster_detail[["CHR", "POS", "mergedClusterID"]],
        left_on=["chr", "POS"],
        right_on=["CHR", "POS"],
        how="left",
    )
    event_map = event_map.dropna(subset=["mergedClusterID"]).copy()
    event_map = event_map[["event_key", "mergedClusterID"]].drop_duplicates()
    log.info("Mapped %d sample-level NUMT events to merged clusters.", len(event_map))
    return event_map


def load_frequency_table(frequency_cluster_tsv: Path) -> pd.DataFrame:
    freq = pd.read_csv(frequency_cluster_tsv, sep="\t", dtype={"merged_cluster_id": str})
    required = {"merged_cluster_id", "frequency_class"}
    missing = required - set(freq.columns)
    if missing:
        raise ValueError(f"Frequency cluster table missing columns: {sorted(missing)}")
    return freq[["merged_cluster_id", "frequency_class"]].drop_duplicates()


def load_nuclear_breakpoint_candidates(
    breakpoint_tsv_gz: Path,
    event_map: pd.DataFrame,
    frequency_df: pd.DataFrame,
    chunksize: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    candidate_frames: list[pd.DataFrame] = []
    total_rows = 0
    nuclear_rows = 0
    mapped_rows = 0
    usecols = [
        "sampleID",
        "chr",
        "Group",
        "readsCount",
        "pos",
        "source_region_chr",
        "source_region_start",
        "source_region_end",
    ]
    for chunk in pd.read_csv(
        breakpoint_tsv_gz,
        sep="\t",
        compression="infer",
        usecols=usecols,
        dtype={"sampleID": str, "chr": str, "Group": str, "source_region_chr": str},
        chunksize=chunksize,
    ):
        total_rows += len(chunk)
        chunk = chunk.loc[chunk["Group"].isin(["nuleft", "nuright"])].copy()
        nuclear_rows += len(chunk)
        if chunk.empty:
            continue
        for column in ["readsCount", "pos", "source_region_start", "source_region_end"]:
            chunk[column] = pd.to_numeric(chunk[column], errors="coerce")
        chunk = chunk.dropna(subset=["pos", "source_region_start", "source_region_end"]).copy()
        chunk["source_region_start"] = chunk["source_region_start"].astype(int)
        chunk["source_region_end"] = chunk["source_region_end"].astype(int)
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
        mapped_rows += len(mapped)
        if mapped.empty:
            continue
        mapped = mapped.merge(frequency_df, left_on="mergedClusterID", right_on="merged_cluster_id", how="inner")
        if mapped.empty:
            continue
        mapped["readsCount"] = pd.to_numeric(mapped["readsCount"], errors="coerce").fillna(0).astype(int)
        mapped["pos"] = pd.to_numeric(mapped["pos"], errors="coerce").round().astype("Int64")
        mapped = mapped.dropna(subset=["pos"]).copy()
        mapped["pos"] = mapped["pos"].astype(int)
        candidate_frames.append(
            mapped[
                [
                    "merged_cluster_id",
                    "frequency_class",
                    "chr",
                    "Group",
                    "pos",
                    "readsCount",
                ]
            ]
        )

    if candidate_frames:
        candidates = pd.concat(candidate_frames, ignore_index=True)
    else:
        candidates = pd.DataFrame(columns=["merged_cluster_id", "frequency_class", "chr", "Group", "pos", "readsCount"])
    qc = {
        "breakpoint_rows_total": int(total_rows),
        "nuclear_breakpoint_rows": int(nuclear_rows),
        "mapped_nuclear_breakpoint_rows": int(mapped_rows),
        "candidate_rows_after_frequency_join": int(len(candidates)),
    }
    return candidates, qc


def collapse_to_distinct_cluster_breakpoints(candidates: pd.DataFrame, flank_bp: int) -> pd.DataFrame:
    if candidates.empty:
        return pd.DataFrame(
            columns=["chr", "start", "end", "breakpoint_id", "merged_cluster_id", "breakpoint_side", "frequency_class"]
        )
    candidates = candidates.loc[candidates["chr"].isin(PRIMARY_CHROMS)].copy()
    candidates = candidates.sort_values(
        ["merged_cluster_id", "Group", "readsCount", "pos"],
        ascending=[True, True, False, True],
    )
    selected = candidates.drop_duplicates(["merged_cluster_id", "Group"], keep="first").copy()
    selected["start"] = selected["pos"].astype(int) - int(flank_bp)
    selected["start"] = selected["start"].mask(selected["start"] < 0, 0)
    selected["end"] = selected["pos"].astype(int) + int(flank_bp)
    selected["breakpoint_id"] = selected["merged_cluster_id"].astype(str) + "__" + selected["Group"].astype(str)
    bed = selected[
        ["chr", "start", "end", "breakpoint_id", "merged_cluster_id", "Group", "frequency_class"]
    ].rename(columns={"Group": "breakpoint_side"})
    bed = bed.loc[bed["end"] > bed["start"]].copy()
    return bed


def write_outputs(
    bed_df: pd.DataFrame,
    requested_classes: list[str],
    target_beds: list[Path],
    prepared_dir: Path,
    enrichment_dir: Path,
    task_manifest_path: Path,
) -> dict[str, int | str]:
    manifest_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for frequency_class in requested_classes:
        if frequency_class == "all":
            class_df = bed_df.copy()
        else:
            class_df = bed_df.loc[bed_df["frequency_class"].astype(str) == frequency_class].copy()
        class_label = sanitize_label(frequency_class)
        bed_path = prepared_dir / f"nuclear_breakpoint_flanks.{class_label}.bed"
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
            task_id = f"{region_id}__{class_label}"
            output_dir = enrichment_dir / region_id / class_label
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

    prepared_summary = prepared_dir / "prepared_breakpoint_flanks_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(prepared_summary, sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(task_manifest_path, sep="\t", index=False)
    return {
        "prepared_summary": str(prepared_summary),
        "task_manifest": str(task_manifest_path),
        "task_count": int(len(manifest_rows)),
    }


def run(
    breakpoint_tsv_gz: str | Path,
    cluster_input_tsv_gz: str | Path,
    cluster_detail_tsv: str | Path,
    frequency_cluster_tsv: str | Path,
    regions_dir: str | Path,
    target_bed_glob: str,
    frequency_classes: str,
    flank_bp: int,
    prepared_dir: str | Path,
    enrichment_dir: str | Path,
    task_manifest: str | Path,
    chunksize: int = 500000,
) -> dict[str, str | int]:
    breakpoint_path = Path(breakpoint_tsv_gz)
    cluster_input_path = Path(cluster_input_tsv_gz)
    cluster_detail_path = Path(cluster_detail_tsv)
    frequency_path = Path(frequency_cluster_tsv)
    regions_path = Path(regions_dir)
    prepared_path = Path(prepared_dir)
    enrichment_path = Path(enrichment_dir)
    task_manifest_path = Path(task_manifest)

    for input_path in [breakpoint_path, cluster_input_path, cluster_detail_path, frequency_path]:
        if not input_path.is_file():
            raise FileNotFoundError(f"Input file not found: {input_path}")
    if not regions_path.is_dir():
        raise NotADirectoryError(f"Target BED directory not found: {regions_path}")
    if flank_bp <= 0:
        raise ValueError("--flank-bp must be positive.")

    prepared_path.mkdir(parents=True, exist_ok=True)
    enrichment_path.mkdir(parents=True, exist_ok=True)
    task_manifest_path.parent.mkdir(parents=True, exist_ok=True)

    requested_classes = parse_classes(frequency_classes)
    target_beds = discover_target_beds(regions_path, target_bed_glob)
    event_map = build_event_to_cluster(cluster_input_path, cluster_detail_path)
    frequency_df = load_frequency_table(frequency_path)
    candidates, qc = load_nuclear_breakpoint_candidates(breakpoint_path, event_map, frequency_df, int(chunksize))
    bed_df = collapse_to_distinct_cluster_breakpoints(candidates, flank_bp)
    all_bed_path = prepared_path / "nuclear_breakpoint_flanks.all_records.bed"
    bed_df.to_csv(all_bed_path, sep="\t", header=False, index=False)

    qc.update(
        {
            "distinct_cluster_breakpoint_flanks": int(len(bed_df)),
            "distinct_clusters_with_breakpoints": int(bed_df["merged_cluster_id"].nunique()) if not bed_df.empty else 0,
            "flank_bp": int(flank_bp),
        }
    )
    qc_path = prepared_path / "prepared_breakpoint_flanks_qc.tsv"
    pd.DataFrame([qc]).to_csv(qc_path, sep="\t", index=False)

    result = write_outputs(bed_df, requested_classes, target_beds, prepared_path, enrichment_path, task_manifest_path)
    result["qc_file"] = str(qc_path)
    result["all_breakpoint_bed"] = str(all_bed_path)
    log.info("Prepared %d breakpoint flank records and %d tasks.", len(bed_df), result["task_count"])
    log.info("QC: %s", qc_path)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare paper-aligned nuclear breakpoint flank BEDs.")
    parser.add_argument("--breakpoint-tsv-gz", required=True)
    parser.add_argument("--cluster-input-tsv-gz", required=True)
    parser.add_argument("--cluster-detail-tsv", required=True)
    parser.add_argument("--frequency-cluster-tsv", required=True)
    parser.add_argument("--regions-dir", required=True)
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
