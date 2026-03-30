#!/usr/bin/env python3

from Bio import AlignIO
import os
import pandas as pd
import sys

aln, numts = sys.argv[1:]

mynumts = pd.read_csv(numts, names=['numtID', 'numt_start', 'numt_end', 'comment'], sep="\t")
mynumts['numtID'] = mynumts['numtID'].astype(str)

mysequence = AlignIO.read(aln, 'fasta')
if len(mysequence) < 2:
    raise ValueError(f"alignment序列数不足(需要至少2条): {aln}")

basename = os.path.basename(aln)
sampleID = basename.split('.')[0] if '.' in basename else basename


def build_tables(ref_seq, alt_seq, alt_id):
    """优化版：避免不必要的内存复制，使用向量化操作"""
    # 直接转字符串，避免中间列表
    ref_str = str(ref_seq)
    alt_str = str(alt_seq)

    # 直接构建DataFrame
    mydf1 = pd.DataFrame({
        'ref': list(ref_str),
        'alt': list(alt_str),
        'pos': range(1, len(ref_str) + 1)
    })

    # 向量化查找起止位置（替代iterrows）
    non_gap_mask = mydf1['alt'] != '-'
    non_gap_pos = mydf1.loc[non_gap_mask, 'pos']

    if len(non_gap_pos) == 0:
        start_pos = 1
        end_pos = 1
    else:
        start_pos = non_gap_pos.iloc[0]
        end_pos = non_gap_pos.iloc[-1]

    mydf2 = mydf1[
        (mydf1['pos'] >= start_pos) &
        (mydf1['pos'] <= end_pos) &
        (mydf1['pos'] <= 16569)
    ].copy()

    mydf2['sampleID'] = sampleID
    mydf2['numtID'] = str(alt_id)
    mydf2['contigID'] = str(alt_id)
    mydf2['numtLength'] = len(mydf2)
    mydf2['numtID'] = mydf2['numtID'].astype(str)

    mydf2 = pd.merge(mydf2, mynumts, on='numtID', how='left')

    mydf3 = mydf2[~(mydf2['ref'] == mydf2['alt'])].copy()
    mydf3_dedup = mydf3.drop_duplicates(subset=['ref', 'alt', 'pos'])
    mydf3['varCount'] = len(mydf3_dedup)
    mydf2['varCount'] = len(mydf3_dedup)

    mydf4 = mydf3[
        (mydf3['pos'] > mydf3['numt_start']) &
        (mydf3['pos'] < mydf3['numt_end'])
    ].copy()

    return mydf2, mydf3_dedup, mydf4


ref_seq = mysequence[0].seq
alt_records = mysequence[1:]

all_mydf2 = []
all_mydf3_dedup = []
all_mydf4 = []

for alt in alt_records:
    df2, df3_dedup, df4 = build_tables(ref_seq, alt.seq, alt.id)
    all_mydf2.append(df2)
    all_mydf3_dedup.append(df3_dedup)
    all_mydf4.append(df4)

# 保持原版的concat逻辑，确保输出完全一致
pd.concat(all_mydf2, ignore_index=True).to_csv(aln + '.numt.tsv', index=0, sep='\t')
pd.concat(all_mydf3_dedup, ignore_index=True).to_csv(aln + '.numtVar.tsv', index=0, sep='\t')
pd.concat(all_mydf4, ignore_index=True).to_csv(aln + '.numtVarFilterPos.tsv', index=0, sep='\t')
