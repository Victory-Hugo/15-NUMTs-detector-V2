#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Conservative breakpoint pairing to derive NUMT lengths.

Input: all_individuals_ConfidentBreakpoints.tsv
Output: paired NUMT intervals (start/end) for downstream frequency plots.
"""

import argparse
import os
import sys
import pandas as pd

LEFT_GROUP = "nu_Tend_Bleft"
RIGHT_GROUP = "nu_Tstart_Bright"


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_breakpoints(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=0)
    # If no header, fall back to legacy columns
    if "pointGroup" not in df.columns or "pos" not in df.columns:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["pointGroup", "Group", "chr", "Tend", "strand", "readsCount", "Tstart", "sampleID"],
        )
        df["pos"] = df["Tend"].where(df["Tend"].notna(), df["Tstart"])
        df = df.rename(columns={"sampleID": "sampleID"})
    return df


def cluster_positions(df: pd.DataFrame, max_gap: int) -> list:
    positions = df[["pos", "readsCount"]].sort_values("pos").values.tolist()
    clusters = []
    current = []
    for pos, count in positions:
        if not current:
            current = [(pos, count)]
            continue
        if pos - current[-1][0] <= max_gap:
            current.append((pos, count))
        else:
            clusters.append(current)
            current = [(pos, count)]
    if current:
        clusters.append(current)
    return clusters


def summarize_cluster(cluster: list) -> dict:
    # Weighted mean of positions by readsCount
    total = sum(c for _, c in cluster)
    if total <= 0:
        total = len(cluster)
    weighted = sum(p * c for p, c in cluster) / float(total)
    return {
        "pos": int(round(weighted)),
        "reads": int(total),
    }


def pair_clusters(lefts: list, rights: list, max_len: int) -> list:
    pairs = []
    used_right = set()
    rights_sorted = sorted(enumerate(rights), key=lambda x: x[1]["pos"])
    for left in sorted(lefts, key=lambda x: x["pos"]):
        best_idx = None
        best_dist = None
        for idx, right in rights_sorted:
            if idx in used_right:
                continue
            if right["pos"] <= left["pos"]:
                continue
            dist = right["pos"] - left["pos"]
            if dist > max_len:
                # rights are sorted; further ones will be even larger
                break
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_idx = idx
        if best_idx is not None:
            used_right.add(best_idx)
            right = rights[best_idx]
            pairs.append((left, right))
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(description="Pair conservative NUMT breakpoints and derive lengths")
    parser.add_argument("--input", required=True, help="Path to all_individuals_ConfidentBreakpoints.tsv")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--cluster-gap", type=int, default=1000, help="Max gap for clustering breakpoints")
    parser.add_argument("--max-length", type=int, default=1000, help="Max length allowed for pairing (bp)")
    args = parser.parse_args()

    df = load_breakpoints(args.input)
    if "readsCount" not in df.columns:
        df["readsCount"] = 1

    # Keep only nuclear left/right breakpoints
    df = df[df["pointGroup"].isin([LEFT_GROUP, RIGHT_GROUP])].copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["readsCount"] = pd.to_numeric(df["readsCount"], errors="coerce").fillna(1).astype(int)
    df = df.dropna(subset=["sampleID", "chr", "pos"])

    results = []
    grouped = df.groupby(["sampleID", "chr"], sort=False)
    for (sample_id, chrom), sub in grouped:
        left_df = sub[sub["pointGroup"] == LEFT_GROUP]
        right_df = sub[sub["pointGroup"] == RIGHT_GROUP]
        if left_df.empty or right_df.empty:
            continue

        left_clusters = [summarize_cluster(c) for c in cluster_positions(left_df, args.cluster_gap)]
        right_clusters = [summarize_cluster(c) for c in cluster_positions(right_df, args.cluster_gap)]
        for left, right in pair_clusters(left_clusters, right_clusters, args.max_length):
            start = int(left["pos"])
            end = int(right["pos"])
            length_bp = end - start + 1
            results.append({
                "sampleID": sample_id,
                "chr": chrom,
                "start": start,
                "end": end,
                "length_bp": length_bp,
                "left_pos": left["pos"],
                "right_pos": right["pos"],
                "left_reads": left["reads"],
                "right_reads": right["reads"],
            })

    columns = [
        "sampleID",
        "chr",
        "start",
        "end",
        "length_bp",
        "left_pos",
        "right_pos",
        "left_reads",
        "right_reads",
    ]
    out_df = pd.DataFrame(results, columns=columns)
    ensure_dir(os.path.dirname(args.output))
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out_df)} paired NUMTs -> {args.output}")


if __name__ == "__main__":
    main()
