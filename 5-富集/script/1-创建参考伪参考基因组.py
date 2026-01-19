#!/usr/bin/env python
#!这个脚本的作用是创建一个用于富集分析的参考基因组文件。
#!简单来说，就是把人类基因组的24条染色体（1-22号染色体 + X和Y染色体）重新排列成一条连续的"虚拟染色体"。
#?输入文件格式：一个BED格式的文件（通常是.bed或.tsv文件）
#?文件内容：包含三列数据，分别是：
# 第1列：染色体名称（如chr1, chr2, chr3...chrX, chrY）
# 第2列：起始位置
# 第3列：结束位置
#?示例输入文件内容：
# chr1    1000    2000
# chr2    5000    6000
# chr3    3000    4000
#?期待的输出文件
# 输出文件名：原输入文件名 + .convert.noGap.tsv
# 文件内容：在原有三列基础上，新增两列：
# refStart_noChr：在虚拟连续基因组中的新起始位置
# refEnd_noChr：在虚拟连续基因组中的新结束位置

import fileinput
import sys, os
import pandas as pd
import glob
import scipy.stats as stats
import numpy as np

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Usage: python 1-创建参考伪参考基因组.py <input_bed_file> [output_file]")
    print("Example: python 1-创建参考伪参考基因组.py genes.bed")
    print("Example: python 1-创建参考伪参考基因组.py genes.bed output.tsv")
    sys.exit(1)

input1 = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) == 3 else input1 + '.convert.noGap.tsv'

if not os.path.exists(input1):
    print(f"Error: Input file '{input1}' does not exist!")
    sys.exit(1) 

df_raw = pd.read_csv(input1, names=['chr','start','end'],sep="\t")
df = df_raw[df_raw.chr.isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])].copy()

refStart_noChr = []
refEnd_noChr = []
for index, row in df.iterrows():
        if row['chr'] == "chr1":
                refStart_noChr.append(row['start'])
                refEnd_noChr.append(row['end'])
        elif row['chr'] == "chr2":
                refStart_noChr.append(row['start'] + 248956422)
                refEnd_noChr.append(row['end'] + 248956422)
        elif row['chr'] == "chr3":
                refStart_noChr.append(row['start'] + 491149951)
                refEnd_noChr.append(row['end'] + 491149951)
        elif row['chr'] == "chr4":
                refStart_noChr.append(row['start'] + 689445510)
                refEnd_noChr.append(row['end'] + 689445510)
        elif row['chr'] == "chr5":
                refStart_noChr.append(row['start'] + 879660065)
                refEnd_noChr.append(row['end'] + 879660065)
        elif row['chr'] == "chr6":
                refStart_noChr.append(row['start'] + 1061198324)
                refEnd_noChr.append(row['end'] + 1061198324)
        elif row['chr'] == "chr7":
                refStart_noChr.append(row['start'] + 1232004303)
                refEnd_noChr.append(row['end'] + 1232004303)
        elif row['chr'] == "chr8":
                refStart_noChr.append(row['start'] + 1391350276)
                refEnd_noChr.append(row['end'] + 1391350276)
        elif row['chr'] == "chr9":
                refStart_noChr.append(row['start'] + 1536488912)
                refEnd_noChr.append(row['end'] + 1536488912)
        elif row['chr'] == "chr10":
                refStart_noChr.append(row['start'] + 1674883629)
                refEnd_noChr.append(row['end'] + 1674883629)
        elif row['chr'] == "chr11":
                refStart_noChr.append(row['start'] + 1808681051)
                refEnd_noChr.append(row['end'] + 1808681051)
        elif row['chr'] == "chr12":
                refStart_noChr.append(row['start'] + 1943767673)
                refEnd_noChr.append(row['end'] + 1943767673)
        elif row['chr'] == "chr13":
                refStart_noChr.append(row['start'] + 2077042982)
                refEnd_noChr.append(row['end'] + 2077042982)
        elif row['chr'] == "chr14":
                refStart_noChr.append(row['start'] + 2191407310)
                refEnd_noChr.append(row['end'] + 2191407310)
        elif row['chr'] == "chr15":
                refStart_noChr.append(row['start'] + 2298451028)
                refEnd_noChr.append(row['end'] + 2298451028)
        elif row['chr'] == "chr16":
                refStart_noChr.append(row['start'] + 2400442217)
                refEnd_noChr.append(row['end'] + 2400442217)
        elif row['chr'] == "chr17":
                refStart_noChr.append(row['start'] + 2490780562)
                refEnd_noChr.append(row['end'] + 2490780562)
        elif row['chr'] == "chr18":
                refStart_noChr.append(row['start'] + 2574038003)
                refEnd_noChr.append(row['end'] + 2574038003)
        elif row['chr'] == "chr19":
                refStart_noChr.append(row['start'] + 2654411288)
                refEnd_noChr.append(row['end'] + 2654411288)
        elif row['chr'] == "chr20":
                refStart_noChr.append(row['start'] + 2713028904)
                refEnd_noChr.append(row['end'] + 2713028904)
        elif row['chr'] == "chr21":
                refStart_noChr.append(row['start'] + 2777473071)
                refEnd_noChr.append(row['end'] + 2777473071)
        elif row['chr'] == "chr22":
                refStart_noChr.append(row['start'] + 2824183054)
                refEnd_noChr.append(row['end'] + 2824183054)
        elif row['chr'] == "chrX":
                refStart_noChr.append(row['start'] + 2875001522)
                refEnd_noChr.append(row['end'] + 2875001522)
        elif row['chr'] == "chrY":
                refStart_noChr.append(row['start'] + 3031042417)
                refEnd_noChr.append(row['end'] + 3031042417)
        else:
                print("Something is wrong!")

df['refStart_noChr'] = refStart_noChr
df['refEnd_noChr'] = refEnd_noChr
df.to_csv(output_file, index=False, sep='\t')
print(f"转换完成，输出文件: {output_file}")
print(f"处理了 {len(df)} 个区域")
