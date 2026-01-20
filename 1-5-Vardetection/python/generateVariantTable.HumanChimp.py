#!/usr/bin/env python3

from Bio import AlignIO
import os
import pandas as pd
import sys

# input aln - multiple alignment file from NUMT sequences
# input numts - target numts contains numtID, numt_start and numt_end
aln, numts = sys.argv[1:]

mynumts = pd.read_csv(numts, names=['numtID', 'numt_start', 'numt_end', 'comment'], sep="\t")
mynumts['numtID'] = mynumts['numtID'].astype(str)

mysequence = AlignIO.read(aln, 'fasta')
if len(mysequence) < 3:
    raise ValueError(f"alignment序列数不足(需要至少3条): {aln}")

basename = os.path.basename(aln)
sampleID = basename.split('.')[0] if '.' in basename else basename


def build_tables(human_seq, chimp_seq, alt_seq, alt_id):
    myseq = [human_seq, chimp_seq, alt_seq]
    myout_list = []
    myseq_zip = list(zip(*myseq))
    for i in range(0, (len(myseq_zip))):
        myout_list.append(myseq_zip[i])

    mydf1 = pd.DataFrame(myout_list)
    mydf1.columns = ['human', 'chimp', 'alt']
    mydf1['pos'] = mydf1.index + 1

    # find start/end by first/last non-gap in alt
    start_pos = 1
    for _, row in mydf1.iterrows():
        if row['alt'] != '-':
            start_pos = row['pos']
            break

    end_pos = len(mydf1)
    mydf1_rev = mydf1.sort_values(by='pos', ascending=False)
    for _, row in mydf1_rev.iterrows():
        if row['alt'] != '-':
            end_pos = row['pos']
            break

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
    mydf4['category'] = ''
    mydf4.loc[mydf4['human'] == mydf4['alt'], 'category'] = 'humanM'
    mydf4.loc[mydf4['chimp'] == mydf4['alt'], 'category'] = 'chimpM'
    mydf4.loc[mydf4['category'] == '', 'category'] = 'none'
    mydf4 = mydf4.drop_duplicates(subset=['human', 'chimp', 'alt', 'pos'])

    humanN = len(mydf4[mydf4['human'] == mydf4['alt']])
    chimpN = len(mydf4[mydf4['chimp'] == mydf4['alt']])
    total_vars = len(mydf4)

    summary_stats = pd.DataFrame({
        'numtID': [str(alt_id)],
        'sampleID': [sampleID],
        'contigID': [str(alt_id)],
        'total': [total_vars],
        'humanN': [humanN],
        'chimpN': [chimpN],
        'pct_human': [humanN / (chimpN + humanN) if (chimpN + humanN) > 0 else 0],
        'age': [6 * (chimpN / (chimpN + humanN)) if (chimpN + humanN) > 0 else 0]
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

pd.concat(all_mydf2, ignore_index=True).to_csv(aln + '.full.tsv', index=0, sep='\t')
pd.concat(all_mydf3, ignore_index=True).to_csv(aln + '.numt.tsv', index=0, sep='\t')
pd.concat(all_mydf4, ignore_index=True).to_csv(aln + '.numtDhumanChimp.tsv', index=0, sep='\t')
pd.concat(all_summary, ignore_index=True).to_csv(aln + '.numtDhumanChimp.sum.tsv', index=False, sep='\t')
