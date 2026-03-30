#!/usr/bin/env python3

from Bio import AlignIO
import os
import pandas as pd
import sys

aln, numts = sys.argv[1:]

mynumts = pd.read_csv(numts, names=['numtID', 'numt_start', 'numt_end', 'comment'], sep="\t")
mynumts['numtID'] = mynumts['numtID'].astype(str)

mysequence = AlignIO.read(aln, 'fasta')
if len(mysequence) < 3:
    raise ValueError(f"alignment序列数不足(需要至少3条): {aln}")

basename = os.path.basename(aln)
sampleID = basename.split('.')[0] if '.' in basename else basename


def build_tables(human_seq, chimp_seq, alt_seq, alt_id):
    """优化版：避免不必要的内存复制，使用向量化操作"""
    human_str = str(human_seq)
    chimp_str = str(chimp_seq)
    alt_str = str(alt_seq)

    mydf1 = pd.DataFrame({
        'human': list(human_str),
        'chimp': list(chimp_str),
        'alt': list(alt_str),
        'pos': range(1, len(human_str) + 1)
    })

    # 向量化查找起止位置
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

    mydf3 = mydf2[
        (mydf2['pos'] > mydf2['numt_start']) &
        (mydf2['pos'] < mydf2['numt_end'])
    ].copy()

    mydf4 = mydf3[mydf3['human'] != mydf3['chimp']].copy()
    mydf4['category'] = 'none'
    mydf4.loc[mydf4['human'] == mydf4['alt'], 'category'] = 'humanM'
    mydf4.loc[mydf4['chimp'] == mydf4['alt'], 'category'] = 'chimpM'
    mydf4 = mydf4.drop_duplicates(subset=['human', 'chimp', 'alt', 'pos'])

    # 向量化计算
    humanN = (mydf4['human'] == mydf4['alt']).sum()
    chimpN = (mydf4['chimp'] == mydf4['alt']).sum()
    total_vars = len(mydf4)
    denom = chimpN + humanN

    summary_stats = pd.DataFrame({
        'numtID': [str(alt_id)],
        'sampleID': [sampleID],
        'contigID': [str(alt_id)],
        'total': [total_vars],
        'humanN': [humanN],
        'chimpN': [chimpN],
        'pct_human': [humanN / denom if denom > 0 else 0],
        'age': [6 * (chimpN / denom) if denom > 0 else 0]
    })

    return mydf2, mydf3, mydf4, summary_stats


human_seq = mysequence[0].seq
chimp_seq = mysequence[1].seq
alt_records = mysequence[2:]

all_mydf2 = []
all_mydf3 = []
all_mydf4 = []
all_summary = []

for alt in alt_records:
    df2, df3, df4, summary = build_tables(human_seq, chimp_seq, alt.seq, alt.id)
    all_mydf2.append(df2)
    all_mydf3.append(df3)
    all_mydf4.append(df4)
    all_summary.append(summary)

# 保持原版的concat逻辑，确保输出完全一致
pd.concat(all_mydf2, ignore_index=True).to_csv(aln + '.full.tsv', index=0, sep='\t')
pd.concat(all_mydf3, ignore_index=True).to_csv(aln + '.numt.tsv', index=0, sep='\t')
pd.concat(all_mydf4, ignore_index=True).to_csv(aln + '.numtDhumanChimp.tsv', index=0, sep='\t')
pd.concat(all_summary, ignore_index=True).to_csv(aln + '.numtDhumanChimp.sum.tsv', index=False, sep='\t')
