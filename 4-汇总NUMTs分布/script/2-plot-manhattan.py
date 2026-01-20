#!/usr/bin/env python3
import argparse
import os
import re

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def chr_to_order(chr_name):
    name = chr_name.strip()
    name = re.sub(r"^chr", "", name, flags=re.IGNORECASE)
    if name.isdigit():
        return int(name)
    if name.upper() == "X":
        return 23
    if name.upper() == "Y":
        return 24
    if name.upper() in ("M", "MT"):
        return 25
    return None


def infer_paths(sum_path, all_path, raw_path):
    if all_path is None and sum_path.endswith(".allCluster.sum.tsv"):
        all_path = sum_path.replace(".allCluster.sum.tsv", ".allCluster.tsv")
    if raw_path is None and sum_path.endswith(".allCluster.sum.tsv"):
        raw_path = sum_path.replace(".allCluster.sum.tsv", "")
    return all_path, raw_path


def output_with_suffix(output_path, suffix):
    base, ext = os.path.splitext(output_path)
    if not ext:
        ext = ".png"
    return f"{base}{suffix}{ext}"

def output_with_suffix_and_ext(output_path, suffix, new_ext):
    base, _ = os.path.splitext(output_path)
    if not new_ext.startswith("."):
        new_ext = f".{new_ext}"
    return f"{base}{suffix}{new_ext}"

def default_report_path(output_path):
    out_dir = os.path.dirname(output_path)
    report_dir = os.path.normpath(os.path.join(out_dir, "..", "markdown"))
    return os.path.join(report_dir, "2-manhattan-report.md")

def write_report(report_path, sum_df, sample_df):
    lines = []
    lines.append("# NUMTs Manhattan summary report")
    lines.append("")
    lines.append("## Summary")
    lines.append(f"- total_clusters: {len(sum_df)}")
    lines.append(f"- chromosomes: {', '.join(sorted(sum_df['CHR'].dropna().unique(), key=lambda x: chr_to_order(x) or 99))}")
    lines.append(f"- POS_mean: {sum_df['POS'].mean():.2f}")
    lines.append(f"- POS_median: {sum_df['POS'].median():.2f}")
    lines.append(f"- POS_max: {sum_df['POS'].max()}")
    top_pos = sum_df.loc[sum_df['POS'].idxmax()]
    lines.append(f"- top_cluster_by_POS: {top_pos['mergedClusterID']} {top_pos['CHR']}:{int(top_pos['min'])}-{int(top_pos['max'])} POS={int(top_pos['POS'])}")
    if sample_df is not None and not sample_df.empty:
        lines.append(f"- sample_count_mean: {sample_df['sample_count'].mean():.2f}")
        lines.append(f"- sample_count_median: {sample_df['sample_count'].median():.2f}")
        lines.append(f"- sample_count_max: {int(sample_df['sample_count'].max())}")
        top_sample = sample_df.loc[sample_df['sample_count'].idxmax()]
        lines.append(f"- top_cluster_by_samples: {top_sample['mergedClusterID']} {top_sample['CHR']}:{int(top_sample['min'])}-{int(top_sample['max'])} samples={int(top_sample['sample_count'])}")
        if sample_df['POS'].notna().all():
            corr = sample_df['POS'].corr(sample_df['sample_count'], method='spearman')
            lines.append(f"- spearman_POS_vs_samples: {corr:.3f}")
    lines.append("")
    lines.append("## Top 10 clusters by POS")
    lines.append("")
    lines.append("|mergedClusterID|CHR|min|max|POS|sample_count|")
    lines.append("|---|---|---|---|---|---|")
    if sample_df is not None and not sample_df.empty:
        top10 = sample_df.sort_values("POS", ascending=False).head(10)
        for _, row in top10.iterrows():
            lines.append(f"|{row['mergedClusterID']}|{row['CHR']}|{int(row['min'])}|{int(row['max'])}|{int(row['POS'])}|{int(row['sample_count'])}|")
    else:
        top10 = sum_df.sort_values("POS", ascending=False).head(10)
        for _, row in top10.iterrows():
            lines.append(f"|{row['mergedClusterID']}|{row['CHR']}|{int(row['min'])}|{int(row['max'])}|{int(row['POS'])}|NA|")
    lines.append("")
    lines.append("## Top 10 clusters by sample count")
    lines.append("")
    lines.append("|mergedClusterID|CHR|min|max|sample_count|POS|")
    lines.append("|---|---|---|---|---|---|")
    if sample_df is not None and not sample_df.empty:
        top10 = sample_df.sort_values("sample_count", ascending=False).head(10)
        for _, row in top10.iterrows():
            lines.append(f"|{row['mergedClusterID']}|{row['CHR']}|{int(row['min'])}|{int(row['max'])}|{int(row['sample_count'])}|{int(row['POS'])}|")
    else:
        lines.append("|NA|NA|NA|NA|NA|NA|")

    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

def build_genome_offsets(df, chr_col="CHR", max_col="max"):
    chr_sizes = (
        df.groupby("chr_order")[max_col]
        .max()
        .sort_index()
        .to_dict()
    )
    offsets = {}
    running = 0
    for chr_order, size in chr_sizes.items():
        offsets[chr_order] = running
        running += int(size)
    return chr_sizes, offsets, running


def add_chr_order(df, chr_col="CHR"):
    df["chr_order"] = df[chr_col].map(chr_to_order)
    df = df.dropna(subset=["chr_order"]).copy()
    df["chr_order"] = df["chr_order"].astype(int)
    return df


def add_genome_x(df, offsets, pos_col):
    df["x"] = df.apply(lambda r: offsets[r["chr_order"]] + r[pos_col], axis=1)
    return df


def manhattan_plot(df, chr_sizes, offsets, running, output_path, y_col, title, y_label):
    df = df.sort_values(["chr_order", "mid_pos"])
    fig, ax = plt.subplots(figsize=(14, 5))
    colors = ["#1f77b4", "#ff7f0e"]
    for i, (chr_order, sub) in enumerate(df.groupby("chr_order")):
        ax.scatter(
            sub["x"],
            sub[y_col],
            s=10,
            color=colors[i % 2],
            alpha=0.7,
            edgecolors="none",
        )

    xticks = []
    xlabels = []
    for chr_order, size in chr_sizes.items():
        start = offsets[chr_order]
        end = start + int(size)
        xticks.append((start + end) / 2.0)
        if chr_order <= 22:
            label = str(chr_order)
        elif chr_order == 23:
            label = "X"
        elif chr_order == 24:
            label = "Y"
        else:
            label = "M"
        xlabels.append(label)
        ax.axvline(end, color="#dddddd", linewidth=0.8, zorder=0)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_xlim(0, running)
    ax.set_ylabel(y_label)
    ax.set_xlabel("Chromosome")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def manhattan_plot_with_sample(df, chr_sizes, offsets, running, output_path):
    df = df.sort_values(["chr_order", "mid_pos"])
    fig, ax = plt.subplots(figsize=(14, 5))
    colors = ["#1f77b4", "#ff7f0e"]

    for i, (chr_order, sub) in enumerate(df.groupby("chr_order")):
        sizes = 15 + (sub["sample_count"].astype(float) * 4.0)
        ax.scatter(
            sub["x"],
            sub["POS"],
            s=sizes,
            c=sub["sample_count"],
            cmap="viridis",
            alpha=0.7,
            edgecolors="none",
        )

    xticks = []
    xlabels = []
    for chr_order, size in chr_sizes.items():
        start = offsets[chr_order]
        end = start + int(size)
        xticks.append((start + end) / 2.0)
        if chr_order <= 22:
            label = str(chr_order)
        elif chr_order == 23:
            label = "X"
        elif chr_order == 24:
            label = "Y"
        else:
            label = "M"
        xlabels.append(label)
        ax.axvline(end, color="#dddddd", linewidth=0.8, zorder=0)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_xlim(0, running)
    ax.set_ylabel("Cluster size (unique positions)")
    ax.set_xlabel("Chromosome")
    ax.set_title("NUMTs cluster Manhattan plot (size=sample count)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.4)

    cbar = fig.colorbar(ax.collections[0], ax=ax)
    cbar.set_label("Samples per cluster")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def density_heatmap(all_df, chr_sizes, offsets, running, output_path, bin_size=1000000):
    all_df = add_chr_order(all_df, chr_col="CHR")
    all_df["mid_pos"] = all_df["POS"].astype(int)
    all_df = add_genome_x(all_df, offsets, pos_col="mid_pos")
    x_bins = max(int(running / bin_size), 1)
    y_bins = list(range(1, 27))

    fig, ax = plt.subplots(figsize=(14, 5))
    hist = ax.hist2d(
        all_df["x"],
        all_df["chr_order"],
        bins=[x_bins, y_bins],
        cmap="viridis",
    )
    fig.colorbar(hist[3], ax=ax, label="Point density")

    xticks = []
    xlabels = []
    for chr_order, size in chr_sizes.items():
        start = offsets[chr_order]
        end = start + int(size)
        xticks.append((start + end) / 2.0)
        if chr_order <= 22:
            label = str(chr_order)
        elif chr_order == 23:
            label = "X"
        elif chr_order == 24:
            label = "Y"
        else:
            label = "M"
        xlabels.append(label)
        ax.axvline(end, color="#dddddd", linewidth=0.8, zorder=0)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_xlim(0, running)
    ax.set_ylabel("Chromosome")
    ax.set_xlabel("Chromosome position (binned)")
    ax.set_yticks([1, 5, 9, 13, 17, 21, 23, 24, 25])
    ax.set_yticklabels(["1", "5", "9", "13", "17", "21", "X", "Y", "M"])
    ax.set_title("NUMTs point density heatmap")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def sample_count_by_cluster(raw_path, sum_df):
    col_names = ["sampleID", "Cluster_No", "disFile", "splitFile", "wgsBAM", "chr", "start", "end"]
    raw = pd.read_csv(raw_path, sep="\t", header=None, names=col_names)
    raw["mid_pos"] = (raw["start"].astype(int) + raw["end"].astype(int)) / 2.0
    raw = raw[["sampleID", "chr", "mid_pos"]]

    sum_df = sum_df.copy()
    sample_rows = []
    for chr_name, sub_sum in sum_df.groupby("CHR"):
        sub_raw = raw[raw["chr"] == chr_name]
        if sub_raw.empty:
            continue
        intervals = pd.IntervalIndex.from_arrays(
            sub_sum["min"].astype(int),
            sub_sum["max"].astype(int),
            closed="both",
        )
        idx = intervals.get_indexer(sub_raw["mid_pos"].astype(int))
        matched = sub_raw[idx >= 0].copy()
        if matched.empty:
            continue
        matched["mergedClusterID"] = sub_sum.iloc[idx[idx >= 0]]["mergedClusterID"].values
        sample_rows.append(matched[["mergedClusterID", "sampleID"]])

    if not sample_rows:
        return pd.DataFrame(columns=["mergedClusterID", "sample_count"])

    sample_df = pd.concat(sample_rows, ignore_index=True)
    sample_counts = sample_df.groupby("mergedClusterID")["sampleID"].nunique().reset_index()
    sample_counts = sample_counts.rename(columns={"sampleID": "sample_count"})
    return sample_counts


def main():
    parser = argparse.ArgumentParser(description="Plot genome-wide Manhattan/scatter for NUMTs clusters")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to merge.tsv.allCluster.sum.tsv",
    )
    parser.add_argument(
        "--all",
        default=None,
        help="Path to merge.tsv.allCluster.tsv (optional)",
    )
    parser.add_argument(
        "--raw",
        default=None,
        help="Path to merge.tsv (optional)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output image path (png)",
    )
    parser.add_argument(
        "--output-density",
        default=None,
        help="Output density heatmap path (png)",
    )
    parser.add_argument(
        "--output-sample",
        default=None,
        help="Output sample-count Manhattan path (png)",
    )
    parser.add_argument(
        "--output-combined",
        default=None,
        help="Output combined Manhattan path (png)",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Output markdown report path",
    )
    parser.add_argument(
        "--output-tsv",
        default=None,
        help="Output TSV for cluster-size Manhattan points",
    )
    parser.add_argument(
        "--output-sample-tsv",
        default=None,
        help="Output TSV for sample-count Manhattan points",
    )
    args = parser.parse_args()

    all_path, raw_path = infer_paths(args.input, args.all, args.raw)
    density_output = args.output_density or output_with_suffix(args.output, ".density")
    sample_output = args.output_sample or output_with_suffix(args.output, ".samplecount")
    combined_output = args.output_combined or output_with_suffix(args.output, ".combined")
    tsv_output = args.output_tsv or output_with_suffix_and_ext(args.output, ".points", ".tsv")
    sample_tsv_output = args.output_sample_tsv or output_with_suffix_and_ext(args.output, ".samplecount.points", ".tsv")
    report_output = args.report or default_report_path(args.output)

    sum_df = pd.read_csv(args.input, sep="\t")
    sum_df = add_chr_order(sum_df, chr_col="CHR")
    sum_df["mid_pos"] = (sum_df["min"].astype(int) + sum_df["max"].astype(int)) / 2.0
    chr_sizes, offsets, running = build_genome_offsets(sum_df, chr_col="CHR", max_col="max")
    sum_df = add_genome_x(sum_df, offsets, pos_col="mid_pos")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    manhattan_plot(
        sum_df,
        chr_sizes,
        offsets,
        running,
        args.output,
        y_col="POS",
        title="NUMTs cluster Manhattan plot",
        y_label="Cluster size (unique positions)",
    )
    sum_df[
        ["mergedClusterID", "CHR", "chr_order", "min", "max", "mid_pos", "x", "POS"]
    ].to_csv(tsv_output, sep="\t", index=False)

    if all_path and os.path.exists(all_path):
        all_df = pd.read_csv(all_path, sep="\t")
        density_heatmap(all_df, chr_sizes, offsets, running, density_output)

    if raw_path and os.path.exists(raw_path):
        sample_counts = sample_count_by_cluster(raw_path, sum_df)
        if not sample_counts.empty:
            sample_df = sum_df.merge(sample_counts, on="mergedClusterID", how="left").fillna(0)
            manhattan_plot(
                sample_df,
                chr_sizes,
                offsets,
                running,
                sample_output,
                y_col="sample_count",
                title="NUMTs cluster Manhattan plot (sample count)",
                y_label="Samples per cluster",
            )
            sample_df[
                ["mergedClusterID", "CHR", "chr_order", "min", "max", "mid_pos", "x", "sample_count", "POS"]
            ].to_csv(sample_tsv_output, sep="\t", index=False)
            manhattan_plot_with_sample(
                sample_df,
                chr_sizes,
                offsets,
                running,
                combined_output,
            )
            write_report(report_output, sum_df, sample_df)
        else:
            write_report(report_output, sum_df, None)
    else:
        write_report(report_output, sum_df, None)


if __name__ == "__main__":
    main()
