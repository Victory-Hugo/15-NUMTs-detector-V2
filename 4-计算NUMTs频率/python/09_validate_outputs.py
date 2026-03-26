#!/usr/bin/env python3
"""Step 09: 输出一致性自检。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

from numts_frequency.io_utils import read_delim

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger(__name__)


def _series_in_unit_interval(series: pd.Series) -> bool:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    if clean.empty:
        return True
    return bool(((clean >= 0) & (clean <= 1)).all())


def _is_monotonic_non_decreasing(series: pd.Series) -> bool:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    if clean.empty:
        return True
    return bool(clean.is_monotonic_increasing)


def _final_value_is_one(series: pd.Series) -> bool:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    if clean.empty:
        return True
    return abs(float(clean.iloc[-1]) - 1.0) < 1e-9


def _read_tsv_or_empty(path: Path, columns: list[str]) -> pd.DataFrame:
    try:
        return read_delim(path, sep="\t")
    except EmptyDataError:
        return pd.DataFrame(columns=columns)


def _validate_mt_table(mt_df: pd.DataFrame, label: str, failures: list[str]) -> None:
    required_overall_cols = {"mt_position", "overall_cumulative_fraction"}
    if not required_overall_cols.issubset(mt_df.columns):
        failures.append(f"{label}: missing required overall fraction columns")
        return
    if not _series_in_unit_interval(mt_df["overall_cumulative_fraction"]):
        failures.append(f"{label}: overall_cumulative_fraction outside [0, 1]")
    overall = mt_df[["mt_position", "overall_cumulative_fraction"]].drop_duplicates().sort_values("mt_position")
    if not _is_monotonic_non_decreasing(overall["overall_cumulative_fraction"]):
        failures.append(f"{label}: overall_cumulative_fraction not monotonic")
    if not _final_value_is_one(overall["overall_cumulative_fraction"]):
        failures.append(f"{label}: final overall_cumulative_fraction is not 1.0")
    if "frequency_class" in mt_df.columns and "cumulative_fraction" in mt_df.columns:
        class_df = mt_df.dropna(subset=["frequency_class"]).copy()
        for freq_class, sub in class_df.groupby("frequency_class", sort=True):
            sub = sub.sort_values("mt_position")
            if not _series_in_unit_interval(sub["cumulative_fraction"]):
                failures.append(f"{label}:{freq_class} cumulative_fraction outside [0, 1]")
            if not _is_monotonic_non_decreasing(sub["cumulative_fraction"]):
                failures.append(f"{label}:{freq_class} cumulative_fraction not monotonic")
            if not _final_value_is_one(sub["cumulative_fraction"]):
                failures.append(f"{label}:{freq_class} final cumulative_fraction is not 1.0")
    elif "frequency_class" in mt_df.columns:
        failures.append(f"{label}: missing cumulative_fraction column")


def run(
    main_results_dir: str,
    mt_results_dir: str,
    stratified_dir: str,
    sample_events_tsv: str,
    report_path: str,
) -> int:
    failures: list[str] = []
    checks: list[str] = []

    main_tables = Path(main_results_dir)
    mt_tables = Path(mt_results_dir)
    stratified_root = Path(stratified_dir)

    frequency_summary = _read_tsv_or_empty(main_tables / "frequency_class_summary.tsv", ["carrier_fraction"])
    carrier_ratio = _read_tsv_or_empty(main_tables / "carrier_ratio_by_frequency_class.tsv", ["sample_ratio"])
    sample_counts = _read_tsv_or_empty(main_tables / "sample_numt_counts.tsv", ["distinct_numt_count"])
    exclusion_summary_path = main_tables / "main_analysis_exclusion_summary.tsv"
    sample_events = _read_tsv_or_empty(Path(sample_events_tsv), ["eligible_for_main_analysis", "main_analysis_exclusion_reason"])
    mt_cumulative = _read_tsv_or_empty(
        mt_tables / "mt_breakpoint_cumulative.tsv",
        ["mt_position", "overall_cumulative_fraction", "frequency_class", "cumulative_fraction"],
    )

    if not _series_in_unit_interval(frequency_summary["carrier_fraction"]):
        failures.append("frequency_class_summary.tsv: carrier_fraction outside [0, 1]")
    else:
        checks.append("frequency_class_summary.tsv: carrier_fraction within [0, 1]")

    if not _series_in_unit_interval(carrier_ratio["sample_ratio"]):
        failures.append("carrier_ratio_by_frequency_class.tsv: sample_ratio outside [0, 1]")
    else:
        checks.append("carrier_ratio_by_frequency_class.tsv: sample_ratio within [0, 1]")

    if (pd.to_numeric(sample_counts["distinct_numt_count"], errors="coerce").fillna(-1) < 0).any():
        failures.append("sample_numt_counts.tsv: distinct_numt_count contains negative values")
    else:
        checks.append("sample_numt_counts.tsv: distinct_numt_count non-negative")

    if not exclusion_summary_path.exists():
        failures.append("main_analysis_exclusion_summary.tsv is missing")
    else:
        checks.append("main_analysis_exclusion_summary.tsv exists")

    if "eligible_for_main_analysis" not in sample_events.columns:
        failures.append("sample_numt_events.annotated.tsv missing eligible_for_main_analysis")
    elif "main_analysis_exclusion_reason" not in sample_events.columns:
        failures.append("sample_numt_events.annotated.tsv missing main_analysis_exclusion_reason")
    else:
        eligible = sample_events["eligible_for_main_analysis"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
        excluded = sample_events.loc[~eligible]
        if not excluded.empty and excluded["main_analysis_exclusion_reason"].fillna("").astype(str).eq("").any():
            failures.append("sample_numt_events.annotated.tsv has excluded events without exclusion reason")
        else:
            checks.append("sample_numt_events.annotated.tsv exclusion reasons populated")

    _validate_mt_table(mt_cumulative, "mt_breakpoint_cumulative.tsv", failures)
    if not any(item.startswith("mt_breakpoint_cumulative.tsv") for item in failures):
        checks.append("mt_breakpoint_cumulative.tsv fractions valid")

    for dimension in ["region", "population", "sex"]:
        dim_file = stratified_root / dimension / "subgroup_mt_breakpoint_cumulative.tsv"
        if not dim_file.exists():
            failures.append(f"{dimension}/subgroup_mt_breakpoint_cumulative.tsv is missing")
            continue
        dim_df = _read_tsv_or_empty(
            dim_file,
            ["mt_position", "overall_cumulative_fraction", "frequency_class", "cumulative_fraction", "subgroup"],
        )
        if dim_df.empty:
            checks.append(f"{dimension}/subgroup_mt_breakpoint_cumulative.tsv empty")
            continue
        for subgroup, sub in dim_df.groupby("subgroup", sort=True):
            _validate_mt_table(sub, f"{dimension}:{subgroup}", failures)
        if not any(item.startswith(f"{dimension}:") for item in failures):
            checks.append(f"{dimension}/subgroup_mt_breakpoint_cumulative.tsv fractions valid")

    lines = ["# 输出一致性自检报告", "", "## 结果"]
    if failures:
        lines.append("- 状态: FAILED")
        lines.append(f"- 失败项数量: {len(failures)}")
    else:
        lines.append("- 状态: PASSED")
        lines.append(f"- 通过项数量: {len(checks)}")
    lines.extend(["", "## 检查项"])
    for item in checks:
        lines.append(f"- PASS: {item}")
    if failures:
        lines.extend(["", "## 失败项"])
        for item in failures:
            lines.append(f"- FAIL: {item}")

    Path(report_path).write_text("\n".join(lines), encoding="utf-8")
    if failures:
        LOG.error("Output validation failed with %d issues", len(failures))
        return 1
    LOG.info("Output validation passed")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate NUMT frequency outputs")
    parser.add_argument("--main-results-dir", required=True)
    parser.add_argument("--mt-results-dir", required=True)
    parser.add_argument("--stratified-dir", required=True)
    parser.add_argument("--sample-events-tsv", required=True)
    parser.add_argument("--report-path", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
