#!/usr/bin/env python
################################################################################
## This script does the enrichment permutation test
## Example of run: python enrichment_simulation.py inputFile 100 10 200 outputFile
################################################################################

import fileinput
import sys, os
import pandas as pd
import glob
import scipy.stats as stats
import numpy as np
import random
try:
    import pybedtools 
    from pybedtools import BedTool
    USE_PYBEDTOOLS = True
except ImportError:
    print("Warning: pybedtools is not installed. Using pandas-based implementation.")
    USE_PYBEDTOOLS = False

if len(sys.argv) < 6 or len(sys.argv) > 7:
    print("Usage: python enrichment_simulation.py inputFile numtsNo numtsTargetNo numtsLength outputPrefix [outputDir]")
    print("Example: python enrichment_simulation.py input.tsv 100 10 200 output")
    print("Example: python enrichment_simulation.py input.tsv 100 10 200 output /path/to/output")
    sys.exit(1)

input1,numtsNo,numtsTargetNo,numtsLength,output = sys.argv[1:6]
output_dir = sys.argv[6] if len(sys.argv) == 7 else os.path.dirname(input1)

numtsNo = int(numtsNo) #total number of NUMTs
numtsTargetNo = int(numtsTargetNo) #number of NUMTs in the target/testing region 
numtsLength = int(numtsLength) #length of flaking region of NUMTs  

def overlap_count_pandas(df1, df2):
    """
    使用pandas计算两个BED格式数据框的重叠数量
    df1: 随机NUMTs DataFrame (columns: chr, start, end)
    df2: 目标区域 DataFrame (columns: chr, start, end)
    """
    overlap_count = 0
    for _, row1 in df1.iterrows():
        for _, row2 in df2.iterrows():
            if (row1['chr'] == row2['chr'] and 
                not (row1['randomEnd'] <= row2['refStart_noChr'] or 
                     row1['randomStart'] >= row2['refEnd_noChr'])):
                overlap_count += 1
                break  # 一个NUMT只计算一次重叠
    return overlap_count

simuRuns = 1001 #number of simulation
numtsTargetPerc = numtsTargetNo/numtsNo
df_ref = pd.read_csv(input1, sep="\t") #bedfile of target/testing regions
df_ref['chr'] = 'chr'
df_ref = df_ref[['chr','refStart_noChr','refEnd_noChr']]

if USE_PYBEDTOOLS:
    df_ref = pybedtools.BedTool.from_dataframe(df_ref)

freq_list=list()
for i in range(1, simuRuns):
    # 使用核基因组范围进行随机采样
    randomStart = random.sample(range(1, 2937639397), numtsNo) # nuclear genome    
    # 注释：如果要使用线粒体基因组，取消下一行注释并注释上一行
    # randomStart = random.sample(range(1, 16570), numtsNo) # for mtDNA genome

    randomNUMTs = pd.DataFrame(randomStart, columns=['randomStart'])
    randomNUMTs['randomStart'] = randomNUMTs['randomStart'].astype(int)
    randomNUMTs['randomEnd'] = randomNUMTs['randomStart'] + numtsLength
    randomNUMTs['chr'] = "chr"
    randomNUMTs = randomNUMTs[['chr','randomStart','randomEnd']]
    
    if USE_PYBEDTOOLS:
        randomNUMTs_bed = pybedtools.BedTool.from_dataframe(randomNUMTs)
        overlap_count = (randomNUMTs_bed + df_ref).count()
    else:
        overlap_count = overlap_count_pandas(randomNUMTs, df_ref)
    
    overlap_freq = overlap_count / numtsNo
    freq_list.append(overlap_freq)

df_pro = pd.DataFrame(freq_list, columns=["Percentage"])
Pvalue_less = 1 - (len(df_pro[df_pro['Percentage'] > numtsTargetPerc ]) / simuRuns)
Pvalue_greater = 1 - (len(df_pro[df_pro['Percentage'] < numtsTargetPerc ]) / simuRuns)

df_pro.columns = ["Pvalue_less(" + str(Pvalue_less) + ") & Pvalue_greater(" + str(Pvalue_greater) + ")"]

print("Pvalue(less):", Pvalue_less)
print("Pvalue(greater):", Pvalue_greater)

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 构造输出文件路径
detail_file = os.path.join(output_dir, f"{os.path.basename(input1)}.enrichmentOutput.{output}.tsv")
summary_file = os.path.join(output_dir, f"enrichment_summary_{output}.tsv")

# 保存详细的模拟数据（原格式）
df_pro.to_csv(detail_file, index=True, sep='\t')

# 创建用户友好的摘要文件
summary_data = {
    'Parameter': [
        'Total_NUMTs', 
        'Target_Region_NUMTs', 
        'NUMT_Length_bp',
        'Simulation_Runs',
        'Observed_Percentage',
        'P_value_less',
        'P_value_greater',
        'Significance_less_than_random',
        'Significance_greater_than_random'
    ],
    'Value': [
        numtsNo,
        numtsTargetNo, 
        numtsLength,
        simuRuns-1,
        f"{numtsTargetPerc:.4f}",
        f"{Pvalue_less:.6f}",
        f"{Pvalue_greater:.6f}",
        "Yes" if Pvalue_less < 0.05 else "No",
        "Yes" if Pvalue_greater < 0.05 else "No"
    ],
    'Description': [
        'Total number of NUMTs in the genome',
        'Number of NUMTs in target regions',
        'Average length of NUMTs',
        'Number of permutation simulations',
        'Observed percentage of NUMTs in target regions',
        'Probability that random distribution has fewer overlaps',
        'Probability that random distribution has more overlaps', 
        'Is target region depleted of NUMTs?',
        'Is target region enriched for NUMTs?'
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(summary_file, index=False, sep='\t')

print(f"富集分析完成！")
print(f"  - 详细结果: {detail_file}")
print(f"  - 摘要结果: {summary_file}")

# 显示统计摘要
print(f"\n=== 模拟结果统计 ===")
simulation_percentages = df_pro.iloc[:, 0].values
print(f"随机模拟中重叠百分比统计:")
print(f"  平均值: {np.mean(simulation_percentages):.6f}")
print(f"  标准差: {np.std(simulation_percentages):.6f}")
print(f"  最小值: {np.min(simulation_percentages):.6f}")
print(f"  最大值: {np.max(simulation_percentages):.6f}")
print(f"  观察值: {numtsTargetPerc:.6f}")
print(f"=====================")
