#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NUMT clustering pipeline

This script reads discordant mtDNA alignments, extracts clusters,
produces summary statistics, and defines breakpoint regions.
"""

import sys
import csv
import json
import os
import numpy as np
import pandas as pd

# Increase CSV field size limit
csv.field_size_limit(sys.maxsize)

# Constants (modify as needed)
INPUT_DISC = sys.argv[1]
SAMPLE_ID = sys.argv[2]
WGS_BAM = sys.argv[3]

# Column names for SAM file
COLUMNS = [
    'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT',
    'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'SM', 'RG', 'NM', 'BC',
    'OC', 'ZX', 'ZY', 'SA'
]

def legacy_mt_names(mt_chr: str, aliases: dict) -> set:
    """Collect MT chromosome names for legacy filtering."""
    mt_names = {mt_chr} if mt_chr else set()
    if aliases:
        mt_names.update({k for k, v in aliases.items() if v == mt_chr})
    return mt_names

def load_config():
    cfg_path = os.environ.get("NUMTS_CONFIG")
    if not cfg_path:
        conf_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "confs")
        )
        preferred = os.path.join(conf_dir, "GRCH38.json")
        fallback = os.path.join(conf_dir, "numts_config.json")
        cfg_path = preferred if os.path.exists(preferred) else fallback
    try:
        with open(cfg_path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        sys.exit(f"❌ 找不到配置文件: {cfg_path}")
    except json.JSONDecodeError as exc:
        sys.exit(f"❌ 配置文件 JSON 解析失败: {cfg_path} ({exc})")


cfg = load_config().get("cluster", {})
VALID_CHROMS = set(cfg.get("valid_chroms", []))
MT_CHR = cfg.get("mt_chrom")
CHROM_ALIASES = cfg.get("chrom_aliases", {})

if not VALID_CHROMS or not MT_CHR:
    sys.exit("❌ 配置缺少 cluster.valid_chroms 或 cluster.mt_chrom")
if MT_CHR not in VALID_CHROMS:
    sys.exit("❌ 配置错误：cluster.mt_chrom 不在 cluster.valid_chroms 中")


def read_and_filter(path: str) -> pd.DataFrame:
    """
    Read SAM-like file and apply initial filtering of RNAME.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        names=COLUMNS,
        comment="@",
        dtype={'RNAME': str, 'RNEXT': str},
        engine="python",
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )
    if CHROM_ALIASES:
        df['RNAME'] = df['RNAME'].map(CHROM_ALIASES).fillna(df['RNAME'])
        df['RNEXT'] = df['RNEXT'].map(CHROM_ALIASES).fillna(df['RNEXT'])
    # Keep only valid chromosomes
    df = df[df['RNAME'].isin(VALID_CHROMS)]
    # 与旧版保持一致：旧版在聚类前过滤 MAPQ > 0，排除多重比对reads
    df = df[df['MAPQ'].astype(int) > 0]
    return df


def fix_rnext(df: pd.DataFrame) -> pd.DataFrame:
    mask_eq  = df['RNEXT'] == '='
    mask_chr = df['RNAME'].isin(VALID_CHROMS - {MT_CHR})
    # 两种子掩码
    mask1 = mask_eq & mask_chr
    mask2 = mask_eq & ~mask_chr

    # 对子掩码分别赋值
    df.loc[mask1, 'RNEXT'] = df.loc[mask1, 'RNAME']
    df.loc[mask2, 'RNEXT'] = MT_CHR
    return df



def select_discordant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only reads where one end maps to MT and the other to a nuclear chromosome.
    """
    mask = (
        (df['RNAME'] != MT_CHR) & (df['RNEXT'] == MT_CHR)
    ) | (
        (df['RNAME'] == MT_CHR) & (df['RNEXT'] != MT_CHR)
    )
    return df[mask].copy()

def read_legacy_df(path: str) -> pd.DataFrame:
    """
    Read SAM-like file using the legacy filtering rules.
    """
    df0 = pd.read_csv(
        path,
        sep="\t",
        names=COLUMNS,
        comment="@",
        dtype={'RNAME': str, 'RNEXT': str},
        engine="python",
        quoting=csv.QUOTE_NONE,
        escapechar="\\"
    )
    df0['RNAME'] = df0['RNAME'].fillna('')
    # Legacy: drop Un_/random/./ contigs
    df = df0[~df0.RNAME.str.contains(r"Un_|random|\.", regex=True)]
    return df

def legacy_cluster_outputs(path: str, sample_id: str) -> pd.DataFrame:
    """
    Reproduce the legacy read-level cluster output.
    """
    df = read_legacy_df(path)
    mt_names = legacy_mt_names(MT_CHR, CHROM_ALIASES)

    df1 = df[
        (df["MAPQ"].astype(int) > 0) &
        (df['RNAME'].isin(mt_names) | df['RNEXT'].isin(mt_names))
    ]
    df3 = df1[~df1.RNAME.isin(mt_names)]
    df3 = df3.sort_values(['RNAME', 'POS'])

    output1 = pd.DataFrame([])
    df_chr = df3.groupby(['RNAME'])
    for clusterID, myclusters in df_chr:
        myclusters['POS'] = myclusters['POS'].astype(int)
        sub_cluster = cluster_positions(myclusters['POS'].tolist(), maxgap=500)
        for x in sub_cluster:
            if len(x) >= 2:
                mycluster = x
                df_cluster = df3[df3.POS.isin(mycluster)]
                df_cluster_pairMT = df[df.QNAME.isin(df_cluster['QNAME'])]
                mt_cluster = cluster_positions(df_cluster['PNEXT'].tolist(), maxgap=500)
                output = pd.DataFrame()
                for y in mt_cluster:
                    df_cluster_pairMT_out = df_cluster_pairMT[df_cluster_pairMT.PNEXT.isin(y)]
                    df_cluster_pairMT_out = df_cluster_pairMT_out.assign(
                        subCluster_No=len(y),
                        Cluster_No=len(x),
                        Cluster_ID=f"{clusterID}_{min(df_cluster['POS'])}_{max(df_cluster['POS'])}_MTboth_{min(df_cluster_pairMT_out['PNEXT'])}_{max(df_cluster_pairMT_out['PNEXT'])}"
                    )
                    output = pd.concat([pd.DataFrame(), df_cluster_pairMT_out], ignore_index=True)

                output1 = pd.concat([output1, output])
                output1["IndividualID"] = sample_id
                if len(output1) != 0:
                    output1['clusterID'] = output1['RNAME'].astype(str) + "_" + output1['Cluster_No'].astype(str)
                else:
                    continue

    return output1

#! 这里指定相邻读段最远的距离，默认上下各自500bp。
def cluster_positions(positions: list, maxgap: int = 500) -> list:
    """
    Cluster sorted positions by max gap.
    """
    if not positions:
        return []
    arr = np.sort(np.array(positions, dtype=int))
    breaks = np.where(np.diff(arr) > maxgap)[0] + 1
    groups = np.split(arr, breaks)
    return [group.tolist() for group in groups if group.size]

#! 这里指定相邻读段最远的距离，默认上下各自500bp。
def extract_regions(df: pd.DataFrame) -> list:
    """
    Generate candidate clusters for nuclear->MT discordant reads.
    """
    regions = []
    for chrom in sorted(df['RNAME'].unique()):
        if chrom == MT_CHR:
            continue
        sub = df[df['RNAME'] == chrom]
        clusters = cluster_positions(sub['POS'].astype(int).tolist())
        for pts in clusters:
            # 文献严格阈值为 >=5 对 discordant reads；当前 >=2 作为初筛，
            # 若需严格标准请将下行改为 if len(pts) < 5
            if len(pts) < 2:
                continue
            start = min(pts) - 500
            # end = max(pts) + 500  # 新版
            end = max(pts) + 500 + 150  # 与旧版保持一致：旧版区间右端额外扩展150bp（read长度）
            regions.append({
                'chr': chrom,
                'start': start,
                'end': end,
                'Cluster_No': len(pts),
                'reads': pts
            })
    return regions


def annotate_mt_positions(df: pd.DataFrame, regions: list) -> None:
    """
    For each region, attach the list of MT positions from discordant pairs.
    """
    for r in regions:
        subset = df[
            (df['RNAME'] == r['chr']) &
            (df['POS'].between(r['start'], r['end'])) &
            (df['RNEXT'] == MT_CHR)
        ]
        r['mt_positions'] = subset['PNEXT'].astype(int).tolist()
        r['subCluster_No'] = len(r['mt_positions'])
        if r['mt_positions']:
            r['Cluster_ID'] = f"{r['chr']}_{r['start']+500}_{r['end']-500}" \
                              f"_{MT_CHR}both_{min(r['mt_positions'])}_{max(r['mt_positions'])}"
        else:
            r['Cluster_ID'] = f"{r['chr']}_{r['start']+500}_{r['end']-500}_{MT_CHR}both_None_None"


def make_breakpoint_df(regions: list) -> pd.DataFrame:
    """
    Prepare breakpoint input DataFrame.
    """
    df = pd.DataFrame(regions)
    df['IndividualID'] = SAMPLE_ID
    df['disFile'] = INPUT_DISC
    df['splitFile'] = INPUT_DISC.replace('disc', 'split')
    df['wgsBAM'] = WGS_BAM
    return df[[
        'IndividualID', 'Cluster_No', 'disFile',
        'splitFile', 'wgsBAM', 'chr', 'start', 'end'
    ]]


def make_summary_df(regions: list) -> pd.DataFrame:
    """
    Generate summary counts for each cluster.
    """
    df = pd.DataFrame(regions)
    keys = ['IndividualID', 'Cluster_ID', 'Cluster_No', 'subCluster_No']
    return df.groupby(keys).size().to_frame('size').reset_index()

#! 输出各个结果文件。
def main():
    # Step 1: Read and filter
    df = read_and_filter(INPUT_DISC)
    df = fix_rnext(df)
    df = select_discordant(df)

    # Step 2: Extract nuclear regions
    regions = extract_regions(df)
    annotate_mt_positions(df, regions)
    for r in regions:
        r['IndividualID'] = SAMPLE_ID
    if not regions:
        legacy_out = legacy_cluster_outputs(INPUT_DISC, SAMPLE_ID)
        if legacy_out.empty:
            legacy_cols = COLUMNS + ['subCluster_No', 'Cluster_No', 'Cluster_ID', 'IndividualID', 'clusterID']
            pd.DataFrame(columns=legacy_cols).to_csv(f"{INPUT_DISC}.cluster.old.tsv",
                                                    sep='\t', header=True, index=False)
        else:
            legacy_out.to_csv(f"{INPUT_DISC}.cluster.old.tsv",
                              sep='\t', header=True, index=False)

        pd.DataFrame(columns=[
            'IndividualID', 'Cluster_No', 'disFile',
            'splitFile', 'wgsBAM', 'chr', 'start', 'end'
        ]).to_csv(f"{INPUT_DISC}.breakpointINPUT.tsv",
                  sep='\t', header=False, index=False)
        pd.DataFrame(columns=[
            'IndividualID', 'Cluster_ID', 'Cluster_No',
            'subCluster_No', 'size'
        ]).to_csv(f"{INPUT_DISC}.cluster.summary.tsv",
                  sep='\t', header=False, index=True)
        pd.DataFrame(columns=[
            'chr', 'start', 'end', 'Cluster_No', 'reads',
            'mt_positions', 'subCluster_No', 'Cluster_ID',
            'IndividualID'
        ]).to_csv(f"{INPUT_DISC}.cluster.tsv",
                  sep='\t', header=True, index=False)

        print("输出完成：")
        print(f"  - {INPUT_DISC}.cluster.tsv")
        print(f"  - {INPUT_DISC}.cluster.summary.tsv")
        print(f"  - {INPUT_DISC}.breakpointINPUT.tsv")
        print(f"  - {INPUT_DISC}.cluster.old.tsv")
        return
    # Step 3: Build DataFrames
    df_break = make_breakpoint_df(regions)
    df_all = pd.DataFrame(regions)
    df_summary = make_summary_df(regions)

    # Step 4: Write outputs
    df_break.to_csv(f"{INPUT_DISC}.breakpointINPUT.tsv",
                    sep='\t', header=False, index=False)
    df_all.to_csv(f"{INPUT_DISC}.cluster.tsv",
                  sep='\t', header=True, index=False)
    df_summary.to_csv(f"{INPUT_DISC}.cluster.summary.tsv",
                      sep='\t', header=False, index=True)
    legacy_out = legacy_cluster_outputs(INPUT_DISC, SAMPLE_ID)
    if legacy_out.empty:
        legacy_cols = COLUMNS + ['subCluster_No', 'Cluster_No', 'Cluster_ID', 'IndividualID', 'clusterID']
        pd.DataFrame(columns=legacy_cols).to_csv(f"{INPUT_DISC}.cluster.old.tsv",
                                                sep='\t', header=True, index=False)
    else:
        legacy_out.to_csv(f"{INPUT_DISC}.cluster.old.tsv",
                          sep='\t', header=True, index=False)

    print("输出完成：")
    print(f"  - {INPUT_DISC}.cluster.tsv")
    print(f"  - {INPUT_DISC}.cluster.summary.tsv")
    print(f"  - {INPUT_DISC}.breakpointINPUT.tsv")
    print(f"  - {INPUT_DISC}.cluster.old.tsv")


if __name__ == '__main__':
    main()
