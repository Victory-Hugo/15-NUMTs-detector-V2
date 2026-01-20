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
if len(mysequence) < 2:
    raise ValueError(f"alignment序列数不足(需要至少2条): {aln}")

basename = os.path.basename(aln)
sampleID = basename.split('.')[0] if '.' in basename else basename


def build_tables(ref_seq, alt_seq, alt_id):
    myseq = [ref_seq, alt_seq]
    myout_list = []
    myseq_zip = list(zip(*myseq))
    for i in range(0, (len(myseq_zip))):
        myout_list.append(myseq_zip[i])

    mydf1 = pd.DataFrame(myout_list)
    mydf1.columns = ['ref', 'alt']
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

pd.concat(all_mydf2, ignore_index=True).to_csv(aln + '.numt.tsv', index=0, sep='\t')
pd.concat(all_mydf3_dedup, ignore_index=True).to_csv(aln + '.numtVar.tsv', index=0, sep='\t')
pd.concat(all_mydf4, ignore_index=True).to_csv(aln + '.numtVarFilterPos.tsv', index=0, sep='\t')
