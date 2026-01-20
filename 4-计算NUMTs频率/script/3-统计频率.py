#!/usr/bin/env python3
import os
import argparse

import matplotlib
import pandas as pd
import re

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def load_cluster_raw(raw_path):
    raw = pd.read_csv(raw_path, sep="\t", header=0)
    if not {"chr", "start", "end"}.issubset(raw.columns):
        raw = pd.read_csv(
            raw_path,
            sep="\t",
            header=None,
            names=["sampleID", "Cluster_No", "disFile", "splitFile", "wgsBAM", "chr", "start", "end"],
        )
    else:
        raw = raw.rename(columns={"IndividualID": "sampleID"})
    total_samples = raw["sampleID"].nunique()
    return raw, total_samples

def parse_positions_span(value):
    if pd.isna(value):
        return None
    nums = [int(x) for x in re.findall(r"\d+", str(value))]
    if not nums:
        return None
    return max(nums) - min(nums) + 1


def build_distinct_numts(raw, total_samples, window_size=1000):
    raw = raw.copy()
    raw["mid_pos"] = (raw["start"].astype(int) + raw["end"].astype(int)) / 2.0
    nuc_len = raw["end"].astype(int) - raw["start"].astype(int) + 1
    if "mt_positions" in raw.columns:
        mt_span = raw["mt_positions"].map(parse_positions_span)
        raw["event_length"] = mt_span.fillna(nuc_len).astype(int)
    else:
        raw["event_length"] = nuc_len
    raw["window"] = (raw["mid_pos"] // window_size).astype(int)
    raw["numt_id"] = raw["chr"].astype(str) + ":" + raw["window"].astype(str)

    grouped = (
        raw.groupby(["chr", "window", "numt_id"], as_index=False)
        .agg(
            min_start=("start", "min"),
            max_end=("end", "max"),
            sample_count=("sampleID", "nunique"),
            median_length=("event_length", "median"),
        )
    )
    grouped["length_bp"] = grouped["median_length"].astype(float).round().astype(int)
    grouped["freq"] = grouped["sample_count"].astype(float) / float(total_samples)
    return grouped


def load_breakpoints(confident_path):
    raw = pd.read_csv(confident_path, sep="\t", header=0)
    raw = raw.rename(columns={"sampleID": "sampleID", "chr": "chr", "pos": "pos"})
    total_samples = raw["sampleID"].nunique()
    return raw, total_samples


def build_distinct_breakpoints(raw, total_samples):
    distinct = raw[["chr", "pos"]].drop_duplicates().copy()
    distinct["length_bp"] = 1
    sample_counts = (
        raw.groupby(["chr", "pos"])["sampleID"]
        .nunique()
        .reset_index()
        .rename(columns={"sampleID": "sample_count"})
    )
    distinct = distinct.merge(sample_counts, on=["chr", "pos"], how="left")
    distinct["freq"] = distinct["sample_count"].astype(float) / float(total_samples)
    return distinct


def plot_length_hist(df, output_path, tsv_path, highlight_short=False, short_cutoff=400):
    out_df = df.copy()
    out_df["is_short"] = out_df["length_bp"] < short_cutoff
    out_df.to_csv(tsv_path, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    if highlight_short:
        short = df[df["length_bp"] < short_cutoff]["length_bp"]
        long = df[df["length_bp"] >= short_cutoff]["length_bp"]
        ax.hist(long, bins=60, color="#6aaed6", alpha=0.8, label=f">= {short_cutoff} bp")
        ax.hist(short, bins=60, color="#d43e3e", alpha=0.8, label=f"< {short_cutoff} bp")
        ax.legend(frameon=False)
        ax.set_title(f"NUMT length distribution (highlight < {short_cutoff} bp)")
    else:
        ax.hist(df["length_bp"], bins=60, color="#4ca1e8", alpha=0.85)
        ax.set_title("NUMT length distribution")
    ax.set_xlabel("Length (bp)")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def classify_frequency(row):
    if row["sample_count"] <= 1:
        return "private"
    if row["freq"] >= 0.01:
        return "common"
    if row["freq"] >= 0.001:
        return "rare"
    return "ultra-rare"


def plot_freq_pie(df, output_path, tsv_path):
    merged = df["freq_class"].replace({"private": "ultra-rare"})
    counts = merged.value_counts().reindex(
        ["common", "rare", "ultra-rare"], fill_value=0
    )
    counts.to_frame(name="count").to_csv(tsv_path, sep="\t")
    labels = counts.index.tolist()
    values = counts.values.tolist()
    colors = ["#3dcd8b", "#4ca1e8", "#dd6a3b"]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.pie(values, labels=labels, autopct="%1.1f%%", colors=colors, startangle=90)
    ax.set_title("NUMT population frequency categories")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_sample_proportion_from_raw(raw, distinct, total_samples, output_path, tsv_path, mode, window_size=1000):
    if mode == "cluster":
        assigned = raw[["sampleID", "chr", "start", "end"]].drop_duplicates()
        assigned["mid_pos"] = (assigned["start"].astype(int) + assigned["end"].astype(int)) / 2.0
        assigned["window"] = (assigned["mid_pos"] // window_size).astype(int)
        assigned["numt_id"] = assigned["chr"].astype(str) + ":" + assigned["window"].astype(str)
        assigned = assigned.merge(
            distinct[["numt_id", "freq_class"]],
            on="numt_id",
            how="left",
        )
    else:
        assigned = raw[["sampleID", "chr", "pos"]].drop_duplicates()
        assigned = assigned.merge(
            distinct[["chr", "pos", "freq_class"]],
            on=["chr", "pos"],
            how="left",
        )
    counts = (
        assigned.groupby("freq_class")["sampleID"]
        .nunique()
        .reindex(["common", "rare", "ultra-rare", "private"], fill_value=0)
    )
    proportions = counts / float(total_samples)
    pd.DataFrame(
        {
            "freq_class": counts.index.tolist(),
            "sample_count": counts.values.tolist(),
            "proportion": proportions.values.tolist(),
        }
    ).to_csv(tsv_path, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(proportions.index.tolist(), proportions.values.tolist(), color="#5a86ee")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion of samples")
    ax.set_title("Samples carrying ≥1 NUMT per frequency class")
    for i, v in enumerate(proportions.values.tolist()):
        ax.text(i, v + 0.01, f"{v:.2%}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="NUMT frequency summary plots (cluster + breakpoint)")
    parser.add_argument(
        "--cluster",
        default="/mnt/l/13_SLE_NUMT/1-所有的NUMTs/data/all_individuals_mt.disc.sam.cluster.tsv",
        help="Path to per-sample cluster TSV",
    )
    parser.add_argument(
        "--confident",
        default="/mnt/l/13_SLE_NUMT/1-所有的NUMTs/data/all_individuals_ConfidentBreakpoints.tsv",
        help="Path to confident breakpoints TSV",
    )
    parser.add_argument(
        "--outdir",
        default="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-计算NUMTs频率/output",
        help="Base output directory",
    )
    parser.add_argument(
        "--outdir-cluster",
        default=None,
        help="Output directory for cluster-based results",
    )
    parser.add_argument(
        "--outdir-breakpoint",
        default=None,
        help="Output directory for breakpoint-based results",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=1000,
        help="Window size in bp for cross-sample clustering (cluster mode)",
    )
    args = parser.parse_args()

    outdir_base = os.path.abspath(args.outdir)
    if not args.outdir_cluster or not args.outdir_breakpoint:
        raise SystemExit("Error: --outdir-cluster and --outdir-breakpoint are required")

    # A: per-sample cluster events
    cluster_outdir = os.path.abspath(args.outdir_cluster)
    ensure_dir(cluster_outdir)
    cluster_raw_path = os.path.abspath(args.cluster)
    cluster_raw, cluster_total_samples = load_cluster_raw(cluster_raw_path)
    cluster_df = build_distinct_numts(cluster_raw, cluster_total_samples, window_size=args.window)
    cluster_df["freq_class"] = cluster_df.apply(classify_frequency, axis=1)

    plot_length_hist(
        cluster_df,
        os.path.join(cluster_outdir, "1-numt_length_hist.png"),
        os.path.join(cluster_outdir, "1-numt_length_hist.tsv"),
        highlight_short=False,
    )
    plot_length_hist(
        cluster_df,
        os.path.join(cluster_outdir, "2-numt_length_hist_short_lt400.png"),
        os.path.join(cluster_outdir, "2-numt_length_hist_short_lt400.tsv"),
        highlight_short=True,
        short_cutoff=400,
    )
    plot_freq_pie(
        cluster_df,
        os.path.join(cluster_outdir, "3-numt_frequency_pie.png"),
        os.path.join(cluster_outdir, "3-numt_frequency_pie.tsv"),
    )
    plot_sample_proportion_from_raw(
        cluster_raw,
        cluster_df,
        cluster_total_samples,
        os.path.join(cluster_outdir, "4-sample_proportion_by_freq.png"),
        os.path.join(cluster_outdir, "4-sample_proportion_by_freq.tsv"),
        mode="cluster",
        window_size=args.window,
    )

    # B: confident breakpoints
    break_outdir = os.path.abspath(args.outdir_breakpoint)
    ensure_dir(break_outdir)
    break_raw_path = os.path.abspath(args.confident)
    break_raw, break_total_samples = load_breakpoints(break_raw_path)
    break_df = build_distinct_breakpoints(break_raw, break_total_samples)
    break_df["freq_class"] = break_df.apply(classify_frequency, axis=1)

    plot_length_hist(
        break_df,
        os.path.join(break_outdir, "1-numt_length_hist.png"),
        os.path.join(break_outdir, "1-numt_length_hist.tsv"),
        highlight_short=False,
    )
    plot_length_hist(
        break_df,
        os.path.join(break_outdir, "2-numt_length_hist_short_lt400.png"),
        os.path.join(break_outdir, "2-numt_length_hist_short_lt400.tsv"),
        highlight_short=True,
        short_cutoff=400,
    )
    plot_freq_pie(
        break_df,
        os.path.join(break_outdir, "3-numt_frequency_pie.png"),
        os.path.join(break_outdir, "3-numt_frequency_pie.tsv"),
    )
    plot_sample_proportion_from_raw(
        break_raw,
        break_df,
        break_total_samples,
        os.path.join(break_outdir, "4-sample_proportion_by_freq.png"),
        os.path.join(break_outdir, "4-sample_proportion_by_freq.tsv"),
        mode="breakpoint",
    )


if __name__ == "__main__":
    main()
