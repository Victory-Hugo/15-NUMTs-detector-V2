#!/usr/bin/env python3
"""QC visualization pipeline for Ancestry/CopyNumber/VerifyBam tables.

Outputs:
- img/summary/*.csv + summary.md
- img/data/*.csv (plot-ready tables)
- img/plots/*.png + *.pdf

Extra rule:
- Samples used for analysis and visualization are restricted to ID intersections
  with IDs from 总表.tsv.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


@dataclass
class Paths:
    base_img: Path
    summary: Path
    data: Path
    plots: Path


ANCESTRY_REQUIRED = ["ID", "优先纳入", "PC", "ContaminatingSample", "IntendedSample"]
COPY_REQUIRED = [
    "ID",
    "Sex",
    "优先纳入",
    "NUMT断点状态",
    "质控状态",
    "mean_mtDNA_depth",
    "mean_autosomal_depth",
    "mtDNA_copy_number",
]
VERIFY_REQUIRED = [
    "ID",
    "优先纳入",
    "#SNPS",
    "#READS",
    "AVG_DP",
    "FREEMIX",
    "FREELK1",
    "FREELK0",
    "status",
]
TOTAL_REQUIRED = ["ID"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC三表归纳与可视化")
    parser.add_argument(
        "--input-dir",
        required=True,
        help="包含 Ancestry.tsv / CopyNumber.tsv / VerifyBam.tsv / 总表.tsv 的目录",
    )
    parser.add_argument("--img-dir", required=True, help="输出目录（会创建 summary/data/plots 子目录）")
    parser.add_argument("--total-name", default="总表.tsv", help="总表文件名，默认: 总表.tsv")
    return parser.parse_args()


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="talk", palette="colorblind")
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"


def ensure_dirs(img_dir: Path) -> Paths:
    summary = img_dir / "summary"
    data = img_dir / "data"
    plots = img_dir / "plots"
    summary.mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    return Paths(base_img=img_dir, summary=summary, data=data, plots=plots)


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8")


def require_columns(df: pd.DataFrame, required: Iterable[str], file_name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{file_name} missing required columns: {missing}")


def clean_id_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip()


def filter_by_total_ids(df: pd.DataFrame, total_ids: set[str], id_col: str = "ID") -> pd.DataFrame:
    out = df.copy()
    out[id_col] = clean_id_series(out[id_col])
    return out[out[id_col].isin(total_ids)].copy()


def save_plot(fig: plt.Figure, plot_path_base: Path) -> None:
    fig.savefig(plot_path_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(plot_path_base.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_counts(series: pd.Series, prefix: str) -> List[Dict[str, str]]:
    vc = series.value_counts(dropna=False)
    return [{"metric": f"{prefix}_{idx}", "value": str(val)} for idx, val in vc.items()]


def summarize_quantiles(values: pd.Series, prefix: str) -> List[Dict[str, str]]:
    q_points = [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]
    quants = values.quantile(q_points)
    records: List[Dict[str, str]] = []
    for q, v in quants.items():
        records.append({"metric": f"{prefix}_q{q:.2f}", "value": f"{v:.10g}"})
    return records


def save_summary_csv(records: List[Dict[str, str]], out_csv: Path) -> pd.DataFrame:
    df = pd.DataFrame(records)
    df.to_csv(out_csv, index=False, encoding="utf-8")
    return df


def process_ancestry(df_raw: pd.DataFrame, total_ids: set[str], paths: Paths) -> Dict[str, str]:
    df = filter_by_total_ids(df_raw, total_ids)
    df["PC_num"] = pd.to_numeric(df["PC"], errors="coerce").astype("Int64")
    df["ContaminatingSample_num"] = pd.to_numeric(df["ContaminatingSample"], errors="coerce")
    df["IntendedSample_num"] = pd.to_numeric(df["IntendedSample"], errors="coerce")

    records: List[Dict[str, str]] = [
        {"metric": "rows_before_intersection", "value": str(len(df_raw))},
        {"metric": "rows_after_intersection", "value": str(len(df))},
        {"metric": "unique_id_after_intersection", "value": str(df["ID"].nunique())},
    ]
    records.extend(summarize_counts(df["优先纳入"], "priority_count"))
    records.extend(summarize_counts(df["PC_num"].astype(str), "pc_count"))
    records.extend(summarize_quantiles(df["ContaminatingSample_num"].dropna(), "contaminating"))
    records.extend(summarize_quantiles(df["IntendedSample_num"].dropna(), "intended"))
    save_summary_csv(records, paths.summary / "ancestry_summary.csv")

    # A1 data and plot
    a1 = df[["ID", "优先纳入", "PC_num", "ContaminatingSample_num", "IntendedSample_num"]].copy()
    a1 = a1.rename(columns={"PC_num": "PC"})
    a1_long = a1.melt(
        id_vars=["ID", "优先纳入", "PC"],
        value_vars=["ContaminatingSample_num", "IntendedSample_num"],
        var_name="role",
        value_name="value",
    )
    a1_long["role"] = a1_long["role"].replace(
        {
            "ContaminatingSample_num": "Contaminating",
            "IntendedSample_num": "Intended",
        }
    )
    a1_long = a1_long.dropna(subset=["PC", "value"]).copy()
    a1_long.to_csv(paths.data / "ancestry_role_distribution_by_pc.csv", index=False, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=a1_long, x="PC", y="value", hue="role", ax=ax, fliersize=0)
    plot_points = a1_long.sample(n=min(6000, len(a1_long)), random_state=42)
    sns.stripplot(
        data=plot_points,
        x="PC",
        y="value",
        hue="role",
        dodge=True,
        size=1.8,
        alpha=0.35,
        ax=ax,
        linewidth=0,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, title="Role")
    ax.set_title("Ancestry Role Distribution by PC")
    ax.set_xlabel("PC")
    ax.set_ylabel("Coordinate Value")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "ancestry_role_distribution_by_pc")

    # A2 data and plot
    a2 = df[["ID", "优先纳入", "PC_num", "ContaminatingSample_num", "IntendedSample_num"]].copy()
    a2 = a2.rename(
        columns={
            "PC_num": "PC",
            "ContaminatingSample_num": "ContaminatingSample",
            "IntendedSample_num": "IntendedSample",
        }
    )
    a2 = a2.dropna(subset=["PC", "ContaminatingSample", "IntendedSample"]).copy()
    a2.to_csv(paths.data / "ancestry_contam_vs_intended_scatter.csv", index=False, encoding="utf-8")

    g = sns.FacetGrid(a2, col="PC", col_wrap=2, height=4.2, sharex=False, sharey=False)
    g.map_dataframe(
        sns.scatterplot,
        x="ContaminatingSample",
        y="IntendedSample",
        alpha=0.45,
        s=12,
        linewidth=0,
    )
    g.set_axis_labels("ContaminatingSample", "IntendedSample")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Contaminating vs Intended by PC")
    for ax in g.axes.flatten():
        sns.despine(ax=ax)
    save_plot(g.fig, paths.plots / "ancestry_contam_vs_intended_scatter")

    # A3 data and plot
    a3 = df[["ID", "优先纳入", "PC_num", "ContaminatingSample_num"]].copy()
    a3["abs_contam"] = a3["ContaminatingSample_num"].abs()
    a3 = a3.dropna(subset=["abs_contam"]).sort_values("abs_contam", ascending=False).head(50).copy()
    a3["rank"] = np.arange(1, len(a3) + 1)
    a3 = a3.rename(columns={"PC_num": "PC", "ContaminatingSample_num": "ContaminatingSample"})
    a3.to_csv(paths.data / "ancestry_abs_contam_top50.csv", index=False, encoding="utf-8")

    a3_plot = a3.copy().sort_values("abs_contam", ascending=True)
    a3_plot["ID_PC"] = a3_plot["ID"] + " | PC" + a3_plot["PC"].astype(str)
    fig, ax = plt.subplots(figsize=(11, 12))
    sns.barplot(data=a3_plot, x="abs_contam", y="ID_PC", hue="PC", dodge=False, ax=ax)
    ax.set_title("Top 50 |ContaminatingSample| Values")
    ax.set_xlabel("Absolute ContaminatingSample")
    ax.set_ylabel("ID | PC")
    ax.legend(frameon=False, title="PC")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "ancestry_abs_contam_top50")

    return {
        "rows_before": str(len(df_raw)),
        "rows_after": str(len(df)),
        "unique_id_after": str(df["ID"].nunique()),
    }


def process_copynumber(df_raw: pd.DataFrame, total_ids: set[str], paths: Paths) -> Dict[str, str]:
    df = filter_by_total_ids(df_raw, total_ids)

    for col in ["mean_mtDNA_depth", "mean_autosomal_depth", "mtDNA_copy_number"]:
        df[f"{col}_num"] = pd.to_numeric(df[col], errors="coerce")

    numeric_cols = ["mean_mtDNA_depth_num", "mean_autosomal_depth_num", "mtDNA_copy_number_num"]
    anomalies = df[df[numeric_cols].isna().any(axis=1)].copy()
    clean = df[df[numeric_cols].notna().all(axis=1)].copy()

    records: List[Dict[str, str]] = [
        {"metric": "rows_before_intersection", "value": str(len(df_raw))},
        {"metric": "rows_after_intersection", "value": str(len(df))},
        {"metric": "unique_id_after_intersection", "value": str(df["ID"].nunique())},
        {"metric": "numeric_valid_rows", "value": str(len(clean))},
        {"metric": "numeric_invalid_rows", "value": str(len(anomalies))},
    ]
    records.extend(summarize_counts(df["质控状态"], "qc_state_count"))
    records.extend(summarize_counts(df["Sex"], "sex_count"))
    records.extend(summarize_quantiles(clean["mean_mtDNA_depth_num"], "mean_mtDNA_depth"))
    records.extend(summarize_quantiles(clean["mean_autosomal_depth_num"], "mean_autosomal_depth"))
    records.extend(summarize_quantiles(clean["mtDNA_copy_number_num"], "mtDNA_copy_number"))
    save_summary_csv(records, paths.summary / "copynumber_summary.csv")

    # anomalies csv
    anomaly_cols = [
        "ID",
        "Sex",
        "优先纳入",
        "NUMT断点状态",
        "质控状态",
        "mean_mtDNA_depth",
        "mean_autosomal_depth",
        "mtDNA_copy_number",
    ]
    anomalies[anomaly_cols].to_csv(paths.data / "copynumber_numeric_anomalies.csv", index=False, encoding="utf-8")

    # C1 data and plot
    c1 = clean[["ID", "Sex", "优先纳入", "质控状态", "mean_autosomal_depth_num", "mean_mtDNA_depth_num"]].copy()
    c1 = c1.rename(
        columns={
            "mean_autosomal_depth_num": "mean_autosomal_depth",
            "mean_mtDNA_depth_num": "mean_mtDNA_depth",
        }
    )
    c1.to_csv(paths.data / "copynumber_depth_scatter.csv", index=False, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 6.5))
    sns.scatterplot(
        data=c1,
        x="mean_autosomal_depth",
        y="mean_mtDNA_depth",
        hue="质控状态",
        style="质控状态",
        alpha=0.65,
        s=28,
        ax=ax,
    )
    ax.set_title("CopyNumber: Autosomal Depth vs mtDNA Depth")
    ax.set_xlabel("Mean Autosomal Depth")
    ax.set_ylabel("Mean mtDNA Depth")
    ax.legend(frameon=False, title="QC")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "copynumber_depth_scatter")

    # C2 data and plot
    c2 = clean[["ID", "Sex", "优先纳入", "质控状态", "mtDNA_copy_number_num"]].copy()
    c2 = c2.rename(columns={"mtDNA_copy_number_num": "mtDNA_copy_number"})
    c2.to_csv(paths.data / "copynumber_distribution_by_qc.csv", index=False, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sns.violinplot(data=c2, x="质控状态", y="mtDNA_copy_number", inner=None, cut=0, ax=ax)
    sns.boxplot(
        data=c2,
        x="质控状态",
        y="mtDNA_copy_number",
        width=0.25,
        showcaps=True,
        boxprops={"facecolor": "white", "alpha": 0.9},
        showfliers=False,
        ax=ax,
    )
    ax.set_title("CopyNumber Distribution by QC Status")
    ax.set_xlabel("QC Status")
    ax.set_ylabel("mtDNA Copy Number")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "copynumber_distribution_by_qc")

    # C3 data and plot
    c3 = pd.DataFrame(
        {
            "group": ["numeric_valid", "numeric_invalid"],
            "n": [len(clean), len(anomalies)],
        }
    )
    c3.to_csv(paths.data / "copynumber_anomaly_count.csv", index=False, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(6.5, 5))
    sns.barplot(data=c3, x="group", y="n", ax=ax)
    ax.set_title("CopyNumber Numeric Data Validity")
    ax.set_xlabel("Group")
    ax.set_ylabel("Sample Count")
    for p in ax.patches:
        ax.annotate(f"{int(p.get_height())}", (p.get_x() + p.get_width() / 2, p.get_height()), ha="center", va="bottom")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "copynumber_anomaly_count")

    return {
        "rows_before": str(len(df_raw)),
        "rows_after": str(len(df)),
        "unique_id_after": str(df["ID"].nunique()),
        "numeric_invalid": str(len(anomalies)),
    }


def process_verifybam(df_raw: pd.DataFrame, total_ids: set[str], paths: Paths) -> Dict[str, str]:
    df = filter_by_total_ids(df_raw, total_ids)

    df["FREEMIX_num"] = pd.to_numeric(df["FREEMIX"], errors="coerce")
    df["AVG_DP_num"] = pd.to_numeric(df["AVG_DP"], errors="coerce")
    df["READS_num"] = pd.to_numeric(df["#READS"], errors="coerce")

    id_consistency_mismatch = 0
    freemix_consistency_mismatch = 0
    if "ID.1" in df.columns:
        id_consistency_mismatch = int((clean_id_series(df["ID.1"]) != clean_id_series(df["ID"])).sum())
    if "freemix" in df.columns:
        freemix_alt = pd.to_numeric(df["freemix"], errors="coerce")
        diff = (df["FREEMIX_num"] - freemix_alt).abs()
        freemix_consistency_mismatch = int((diff > 1e-6).fillna(True).sum())

    bins = [0, 0.001, 0.005, 0.01, 0.02, 0.05, np.inf]
    labels = ["<=0.001", "0.001-0.005", "0.005-0.01", "0.01-0.02", "0.02-0.05", ">0.05"]
    df["freemix_bin"] = pd.cut(df["FREEMIX_num"], bins=bins, labels=labels, include_lowest=True)

    records: List[Dict[str, str]] = [
        {"metric": "rows_before_intersection", "value": str(len(df_raw))},
        {"metric": "rows_after_intersection", "value": str(len(df))},
        {"metric": "unique_id_after_intersection", "value": str(df["ID"].nunique())},
        {"metric": "id_consistency_mismatch_count", "value": str(id_consistency_mismatch)},
        {"metric": "freemix_consistency_mismatch_count", "value": str(freemix_consistency_mismatch)},
    ]
    records.extend(summarize_counts(df["status"], "status_count"))
    records.extend(summarize_counts(df["优先纳入"], "priority_count"))
    records.extend(summarize_quantiles(df["FREEMIX_num"].dropna(), "freemix"))
    records.extend(summarize_quantiles(df["AVG_DP_num"].dropna(), "avg_dp"))

    bin_counts = df["freemix_bin"].value_counts().reindex(labels).fillna(0).astype(int)
    for k, v in bin_counts.items():
        records.append({"metric": f"freemix_bin_count_{k}", "value": str(v)})

    # acceptance check
    fail_below_or_equal_threshold = int(((df["status"] == "FAIL") & (df["FREEMIX_num"] <= 0.02 + 1e-12)).sum())
    records.append({"metric": "fail_with_freemix_le_0.02", "value": str(fail_below_or_equal_threshold)})

    save_summary_csv(records, paths.summary / "verifybam_summary.csv")

    # V1 data and plot
    v1 = df[["ID", "优先纳入", "FREEMIX_num", "status"]].copy()
    v1 = v1.rename(columns={"FREEMIX_num": "FREEMIX"})
    v1 = v1.dropna(subset=["FREEMIX"]).copy()
    v1 = v1[v1["FREEMIX"] > 0].copy()
    v1["log10_FREEMIX"] = np.log10(v1["FREEMIX"])
    v1.to_csv(paths.data / "verifybam_freemix_hist.csv", index=False, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(8.5, 6))
    sns.histplot(data=v1, x="FREEMIX", bins=60, ax=ax, color="#4C72B0")
    ax.set_xscale("log")
    ax.axvline(0.02, color="#D62728", linestyle="--", linewidth=1.5, label="Threshold=0.02")
    ax.set_title("VerifyBam FREEMIX Distribution")
    ax.set_xlabel("FREEMIX (log scale)")
    ax.set_ylabel("Count")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "verifybam_freemix_hist")

    # V2 data and plot
    v2 = df[["ID", "优先纳入", "AVG_DP_num", "FREEMIX_num", "status"]].copy()
    v2 = v2.rename(columns={"AVG_DP_num": "AVG_DP", "FREEMIX_num": "FREEMIX"})
    v2 = v2.dropna(subset=["AVG_DP", "FREEMIX"]).copy()
    v2 = v2[v2["FREEMIX"] > 0].copy()
    v2.to_csv(paths.data / "verifybam_depth_vs_freemix.csv", index=False, encoding="utf-8")

    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    sns.scatterplot(data=v2, x="AVG_DP", y="FREEMIX", hue="status", alpha=0.68, s=24, ax=ax)
    ax.set_yscale("log")
    ax.axhline(0.02, color="#D62728", linestyle="--", linewidth=1.3)
    ax.set_title("AVG_DP vs FREEMIX")
    ax.set_xlabel("AVG_DP")
    ax.set_ylabel("FREEMIX (log scale)")
    ax.legend(frameon=False, title="status")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "verifybam_depth_vs_freemix")

    # V3 data and plot
    v3 = (
        df.groupby(["freemix_bin", "status"], dropna=False, observed=False)
        .size()
        .reset_index(name="n")
        .rename(columns={"freemix_bin": "freemix_bin"})
    )
    v3 = v3[v3["freemix_bin"].notna()].copy()
    v3["freemix_bin"] = pd.Categorical(v3["freemix_bin"], categories=labels, ordered=True)
    v3 = v3.sort_values(["freemix_bin", "status"]).copy()
    total_per_bin = v3.groupby("freemix_bin", observed=False)["n"].transform("sum")
    v3["proportion"] = (v3["n"] / total_per_bin).fillna(0)
    v3.to_csv(paths.data / "verifybam_freemix_bin_status.csv", index=False, encoding="utf-8")

    pivot = v3.pivot(index="freemix_bin", columns="status", values="n").fillna(0)
    fig, ax = plt.subplots(figsize=(9, 6))
    bottom = np.zeros(len(pivot), dtype=float)
    for status in pivot.columns:
        vals = pivot[status].values
        ax.bar(pivot.index.astype(str), vals, bottom=bottom, label=status)
        bottom += vals
    ax.set_title("FREEMIX Bin by Status")
    ax.set_xlabel("FREEMIX Bin")
    ax.set_ylabel("Count")
    ax.legend(frameon=False, title="status")
    sns.despine(ax=ax)
    save_plot(fig, paths.plots / "verifybam_freemix_bin_status")

    return {
        "rows_before": str(len(df_raw)),
        "rows_after": str(len(df)),
        "unique_id_after": str(df["ID"].nunique()),
        "fail_le_0.02": str(fail_below_or_equal_threshold),
    }


def write_summary_md(paths: Paths, total_id_count: int, ancestry_meta: Dict[str, str], copy_meta: Dict[str, str], verify_meta: Dict[str, str]) -> None:
    lines = [
        "# QC三表归纳总结",
        "",
        f"- 总表ID集合大小: {total_id_count}",
        "- 规则: 三个输入TSV均先与总表ID取交集，再进入归纳与可视化。",
        "",
        "## Ancestry.tsv",
        f"- rows_before_intersection: {ancestry_meta['rows_before']}",
        f"- rows_after_intersection: {ancestry_meta['rows_after']}",
        f"- unique_id_after_intersection: {ancestry_meta['unique_id_after']}",
        "",
        "## CopyNumber.tsv",
        f"- rows_before_intersection: {copy_meta['rows_before']}",
        f"- rows_after_intersection: {copy_meta['rows_after']}",
        f"- unique_id_after_intersection: {copy_meta['unique_id_after']}",
        f"- numeric_invalid_rows: {copy_meta['numeric_invalid']}",
        "",
        "## VerifyBam.tsv",
        f"- rows_before_intersection: {verify_meta['rows_before']}",
        f"- rows_after_intersection: {verify_meta['rows_after']}",
        f"- unique_id_after_intersection: {verify_meta['unique_id_after']}",
        f"- fail_with_freemix_le_0.02: {verify_meta['fail_le_0.02']}",
        "",
        "详细统计见 summary/*.csv；每张图的绘图配套数据见 data/*.csv。",
    ]
    (paths.summary / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup_style()

    input_dir = Path(args.input_dir)
    img_dir = Path(args.img_dir)
    paths = ensure_dirs(img_dir)

    ancestry_path = input_dir / "Ancestry.tsv"
    copy_path = input_dir / "CopyNumber.tsv"
    verify_path = input_dir / "VerifyBam.tsv"
    total_path = input_dir / args.total_name

    ancestry_raw = read_tsv(ancestry_path)
    copy_raw = read_tsv(copy_path)
    verify_raw = read_tsv(verify_path)
    total_raw = read_tsv(total_path)

    require_columns(ancestry_raw, ANCESTRY_REQUIRED, ancestry_path.name)
    require_columns(copy_raw, COPY_REQUIRED, copy_path.name)
    require_columns(verify_raw, VERIFY_REQUIRED, verify_path.name)
    require_columns(total_raw, TOTAL_REQUIRED, total_path.name)

    total_ids = set(clean_id_series(total_raw["ID"]).dropna().tolist())

    ancestry_meta = process_ancestry(ancestry_raw, total_ids, paths)
    copy_meta = process_copynumber(copy_raw, total_ids, paths)
    verify_meta = process_verifybam(verify_raw, total_ids, paths)

    write_summary_md(paths, len(total_ids), ancestry_meta, copy_meta, verify_meta)

    print("Done.")
    print(f"Summary dir: {paths.summary}")
    print(f"Data dir: {paths.data}")
    print(f"Plots dir: {paths.plots}")


if __name__ == "__main__":
    main()
