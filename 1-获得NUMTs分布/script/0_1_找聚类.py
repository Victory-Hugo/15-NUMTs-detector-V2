#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NUMT clustering pipeline

This script reads discordant mtDNA alignments, extracts clusters,
produces summary statistics, and defines breakpoint regions.
"""

import sys
import csv
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

# Valid chromosomes
VALID_CHROMS = {str(i) for i in range(1, 23)} | {'X', 'Y', 'MT'}


def read_and_filter(path: str) -> pd.DataFrame:
    """
    Read SAM-like file and apply initial filtering of RNAME.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        names=COLUMNS,
        comment="@",
        low_memory=False,
        dtype={'RNAME': str, 'RNEXT': str}
    )
    # Keep only valid chromosomes
    return df[df['RNAME'].isin(VALID_CHROMS)]


def fix_rnext(df: pd.DataFrame) -> pd.DataFrame:
    mask_eq  = df['RNEXT'] == '='
    mask_chr = df['RNAME'].isin(VALID_CHROMS - {'MT'})
    # 两种子掩码
    mask1 = mask_eq & mask_chr
    mask2 = mask_eq & ~mask_chr

    # 对子掩码分别赋值
    df.loc[mask1, 'RNEXT'] = df.loc[mask1, 'RNAME']
    df.loc[mask2, 'RNEXT'] = 'MT'
    return df



def select_discordant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only reads where one end maps to MT and the other to a nuclear chromosome.
    """
    mask = (
        (df['RNAME'] != 'MT') & (df['RNEXT'] == 'MT')
    ) | (
        (df['RNAME'] == 'MT') & (df['RNEXT'] != 'MT')
    )
    return df[mask].copy()

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
        if chrom == 'MT':
            continue
        sub = df[df['RNAME'] == chrom]
        clusters = cluster_positions(sub['POS'].astype(int).tolist())
        for pts in clusters:
            if len(pts) < 2:
                continue
            start = min(pts) - 500
            end = max(pts) + 500
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
            (df['RNEXT'] == 'MT')
        ]
        r['mt_positions'] = subset['PNEXT'].astype(int).tolist()
        r['subCluster_No'] = len(r['mt_positions'])
        if r['mt_positions']:
            r['Cluster_ID'] = f"{r['chr']}_{r['start']+500}_{r['end']-500}" \
                              f"_MTboth_{min(r['mt_positions'])}_{max(r['mt_positions'])}"
        else:
            r['Cluster_ID'] = f"{r['chr']}_{r['start']+500}_{r['end']-500}_MTboth_None_None"


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

    print("输出完成：")
    print(f"  - {INPUT_DISC}.cluster.tsv")
    print(f"  - {INPUT_DISC}.cluster.summary.tsv")
    print(f"  - {INPUT_DISC}.breakpointINPUT.tsv")


if __name__ == '__main__':
    main()
