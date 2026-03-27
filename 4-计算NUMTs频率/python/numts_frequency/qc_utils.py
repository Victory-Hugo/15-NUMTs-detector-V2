#!/usr/bin/env python3
"""输入检查与样本过滤逻辑。"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .catalog_utils import build_breakpoint_index, query_interval, query_interval_by_key
from .io_utils import write_tsv
from .parse_utils import build_region_key
from .parse_utils import LEFT_GROUP, RIGHT_GROUP

LOG = logging.getLogger(__name__)


def build_input_summary(
    breakpoints: pd.DataFrame,
    clusters: pd.DataFrame,
    meta: pd.DataFrame,
    length_summary: pd.DataFrame,
    meta_qc_col: str,
    meta_qc_pass_value: str,
) -> pd.DataFrame:
    summary = [
        {"dataset": "breakpoints", "rows": len(breakpoints), "unique_samples": breakpoints["sample_id"].nunique()},
        {"dataset": "clusters", "rows": len(clusters), "unique_samples": clusters["sample_id"].nunique()},
        {"dataset": "meta", "rows": len(meta), "unique_samples": meta["sample_id"].nunique()},
        {"dataset": "length_bed_summary", "rows": len(length_summary), "unique_samples": length_summary["sample_id"].nunique()},
        {
            "dataset": "meta_pass",
            "rows": int(meta[meta_qc_col].astype(str).eq(meta_qc_pass_value).sum()),
            "unique_samples": meta.loc[meta[meta_qc_col].astype(str).eq(meta_qc_pass_value), "sample_id"].nunique(),
        },
    ]
    return pd.DataFrame(summary)


def build_sample_filters(
    breakpoints: pd.DataFrame,
    clusters: pd.DataFrame,
    meta: pd.DataFrame,
    meta_qc_col: str,
    meta_qc_pass_value: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    bp_ids = set(breakpoints["sample_id"].astype(str))
    cluster_ids = set(clusters["sample_id"].astype(str))
    pass_ids = set(meta.loc[meta[meta_qc_col].astype(str).eq(meta_qc_pass_value), "sample_id"].astype(str))
    retained_ids = sorted(pass_ids & bp_ids & cluster_ids)
    all_ids = sorted(set(meta["sample_id"].astype(str)) | bp_ids | cluster_ids)

    inclusion = pd.DataFrame(
        {
            "sample_id": retained_ids,
            "in_meta": True,
            "meta_qc_pass": True,
            "in_breakpoint": True,
            "in_cluster": True,
            "retained": True,
        }
    )

    exclusion_records = []
    for sample_id in all_ids:
        in_meta = sample_id in set(meta["sample_id"].astype(str))
        meta_qc_pass = sample_id in pass_ids
        in_breakpoint = sample_id in bp_ids
        in_cluster = sample_id in cluster_ids
        retained = sample_id in retained_ids
        reasons = []
        if not in_meta:
            reasons.append("missing_in_meta")
        if in_meta and not meta_qc_pass:
            reasons.append("meta_qc_fail")
        if not in_breakpoint:
            reasons.append("missing_in_breakpoint")
        if not in_cluster:
            reasons.append("missing_in_cluster")
        if retained:
            continue
        exclusion_records.append(
            {
                "sample_id": sample_id,
                "in_meta": in_meta,
                "meta_qc_pass": meta_qc_pass,
                "in_breakpoint": in_breakpoint,
                "in_cluster": in_cluster,
                "retained": retained,
                "exclude_reason": ",".join(reasons),
            }
        )
    exclusion = pd.DataFrame.from_records(exclusion_records)
    return inclusion, exclusion


def build_event_traceability(breakpoints: pd.DataFrame, clusters: pd.DataFrame, retained_ids: set[str]) -> pd.DataFrame:
    bp = breakpoints[breakpoints["sample_id"].isin(retained_ids)].copy()
    nuclear_bp = bp[bp["pointGroup"].isin([LEFT_GROUP, RIGHT_GROUP])].copy()
    fallback_index = build_breakpoint_index(nuclear_bp)
    region_index = {}
    use_region_key_match = False
    if "source_region_key" in nuclear_bp.columns:
        region_bp = nuclear_bp[nuclear_bp["source_region_key"].notna()].copy()
        if not region_bp.empty:
            region_index = build_breakpoint_index(region_bp, group_cols=("source_region_key", "chr"))
            use_region_key_match = True
    results = []
    for row in clusters[clusters["sample_id"].isin(retained_ids)].itertuples(index=False):
        region_key = build_region_key(row.sample_id, row.chr, row.start, row.end)
        if use_region_key_match:
            sub = query_interval_by_key(region_index, (region_key, row.chr), row.start, row.end)
            traceability_mode = "source_region_key"
        else:
            sub = query_interval(fallback_index, row.sample_id, row.chr, row.start, row.end)
            traceability_mode = "sample_chr_interval"
        results.append(
            {
                "sample_id": row.sample_id,
                "cluster_id_raw": row.cluster_id_raw,
                "chr": row.chr,
                "cluster_start": int(row.start),
                "cluster_end": int(row.end),
                "breakpoint_region_key": region_key if use_region_key_match else pd.NA,
                "traceability_mode": traceability_mode,
                "supporting_nuclear_breakpoint_rows": int(sub["positions"].shape[0]),
                "traceability_status": "ok" if sub["positions"].shape[0] > 0 else "event_without_confident_nuclear_breakpoint",
            }
        )
    return pd.DataFrame.from_records(results)


def build_validation_report(
    output_path: str,
    input_summary: pd.DataFrame,
    inclusion: pd.DataFrame,
    exclusion: pd.DataFrame,
    traceability: pd.DataFrame,
    breakpoints: pd.DataFrame,
    clusters: pd.DataFrame,
    meta: pd.DataFrame,
) -> None:
    duplicated_bp_cols = ["sample_id", "chr", "pointGroup", "pos"]
    if "source_region_key" in breakpoints.columns and breakpoints["source_region_key"].notna().any():
        duplicated_bp_cols.append("source_region_key")
    duplicated_bp = int(breakpoints.duplicated(duplicated_bp_cols).sum())
    duplicated_cluster = int(clusters.duplicated(["sample_id", "cluster_id_raw"]).sum())
    duplicated_meta = int(meta.duplicated(["sample_id"]).sum())
    inconsistent_bp_chr = int((~breakpoints["chr"].astype(str).str.startswith("chr")).sum())
    inconsistent_cluster_chr = int((~clusters["chr"].astype(str).str.startswith("chr")).sum())
    traceability_fail = int((traceability["traceability_status"] != "ok").sum())

    lines = [
        "# 输入一致性检查报告",
        "",
        "## 数据规模",
    ]
    for row in input_summary.itertuples(index=False):
        lines.append(f"- {row.dataset}: rows={row.rows}, unique_samples={row.unique_samples}")
    lines.extend([
        "",
        "## 重复与命名",
        f"- breakpoint 重复行数: {duplicated_bp}",
        f"- cluster 重复行数: {duplicated_cluster}",
        f"- meta 重复样本数: {duplicated_meta}",
        f"- breakpoint 非 chr 命名行数: {inconsistent_bp_chr}",
        f"- cluster 非 chr 命名行数: {inconsistent_cluster_chr}",
        "",
        "## 样本过滤",
        f"- 保留样本数: {len(inclusion)}",
        f"- 排除样本数: {len(exclusion)}",
        "",
        "## 事件可追踪性",
        f"- 无 confident nuclear breakpoint 的事件数: {traceability_fail}",
    ])
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")


def write_qc_outputs(
    output_dir: str,
    input_summary: pd.DataFrame,
    inclusion: pd.DataFrame,
    exclusion: pd.DataFrame,
    traceability: pd.DataFrame,
    length_summary: pd.DataFrame,
) -> None:
    write_tsv(input_summary, str(Path(output_dir) / "01-input_summary.tsv"))
    write_tsv(inclusion, str(Path(output_dir) / "02-sample_inclusion.tsv"))
    write_tsv(exclusion, str(Path(output_dir) / "03-sample_exclusion.tsv"))
    write_tsv(traceability, str(Path(output_dir) / "04-event_traceability.tsv"))
    write_tsv(length_summary, str(Path(output_dir) / "_length_bed_region_summary.tsv"))
