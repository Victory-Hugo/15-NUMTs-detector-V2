#!/usr/bin/env python

################################################################################
## This script takes the cluster sum files generated from searchNumtCluster_fromDiscordantReads.py
## to look for shared NUMT clusters across different individuals
################################################################################

import fileinput
import sys, os
import numpy as np
import pandas as pd
import glob
import scipy.stats as stats

def cluster(data, maxgap):
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


input1 = sys.argv[1]

def check_and_fix_file_format(input_file):
    """检查文件格式并修复不规范的染色体格式"""
    
    # 读取文件的第一行来检查格式
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()
    
    # 检查染色体列（第6列）是否为元组格式
    columns = first_line.split('\t')
    if len(columns) >= 6:
        chr_column = columns[5]
        # 检查是否包含元组格式的特征：包含括号和引号
        if '(' in chr_column and ')' in chr_column and "'" in chr_column:
            print("检测到格式不规范：染色体列为元组格式")
            print("格式不规范，备份原文件，生成规范的文件！")
            
            # 备份原文件
            backup_file = input_file + '.bak'
            os.rename(input_file, backup_file)
            print(f"原文件已备份为: {backup_file}")
            
            # 读取备份文件并修复格式
            column_names = ['sampleID', 'Cluster_No', 'disFile', 'splitFile', 'wgsBAM', 'chr', 'start', 'end']
            df_temp = pd.read_csv(backup_file, sep="\t", engine='python', header=None, names=column_names)
            
            # 清理染色体列格式，处理类似 ('chr1',) 的元组字符串
            df_temp['chr'] = df_temp['chr'].astype(str).str.replace(r"[()']", "", regex=True).str.replace(",", "").str.strip()
            
            # 保存修复后的文件
            df_temp.to_csv(input_file, sep="\t", header=False, index=False)
            print(f"规范化文件已生成: {input_file}")
            print(f"染色体格式已从元组格式修复为标准格式")
            
            return True
        else:
            print("检测到文件格式规范，直接处理")
            return False
    else:
        print("警告: 文件格式可能有问题，列数不足")
        return False

# 检查并修复文件格式
format_fixed = check_and_fix_file_format(input1)

if format_fixed:
    print("重新读取修复后的规范文件...")
else:
    print("直接读取原始规范文件...")

# 读取新版本的breakpointINPUT.tsv文件（无表头，8列）
# 列顺序：sampleID, Cluster_No, disFile, splitFile, wgsBAM, chr, start, end
column_names = ['sampleID', 'Cluster_No', 'disFile', 'splitFile', 'wgsBAM', 'chr', 'start', 'end']
df0 = pd.read_csv(input1, sep="\t", engine='python', header=None, names=column_names)

print(f"成功读取数据，共 {len(df0)} 行")

# 去重并筛选需要的列
df = df0[['sampleID', 'chr', 'start', 'end']].drop_duplicates(["sampleID","chr","start"])
#df = df0[~df0.RNAME.str.contains("chrUn_")]

##### order by chromosome and pos #####
df['position'] = (df['start'].astype(int) + df['end'].astype(int))/2
df.sort_values(['chr','position'], inplace=True)
df = df.drop_duplicates(['chr','position'])
##### look for the clutser by mapgap on pos #####
output1_list = []
df1 = df.groupby(['chr'])
for clusterID, myclusters in df1:
    # 提取真实的染色体名称（从元组中取第一个元素）
    chr_name = clusterID[0] if isinstance(clusterID, tuple) else clusterID
    
    myclusters['position'] = myclusters['position'].astype(int)
    sub_cluster = cluster(myclusters['position'].tolist(), maxgap=1000)
    
    # 处理聚类结果
    cluster_data = []
    for group_idx, positions in enumerate(sub_cluster):
        for pos in positions:
            cluster_data.append({
                'GroupID': group_idx,
                'Index': 0,  # 可以根据需要调整
                'POS': pos,
                'CHR': chr_name  # 使用提取的染色体名称
            })
    
    if cluster_data:
        df2 = pd.DataFrame(cluster_data)
        output1_list.append(df2)

# 合并所有结果
if output1_list:
    output1 = pd.concat(output1_list, ignore_index=True)
else:
    output1 = pd.DataFrame(columns=['GroupID','Index','POS','CHR'])
output1["mergedClusterID"] = output1['GroupID'].astype(str) + "_" +  output1['CHR'].astype(str)

output2_1 = output1.groupby(['mergedClusterID']).count().reset_index()
output2_2 = output1.groupby(['mergedClusterID']).POS.agg(['min','max']).reset_index()
output2_1 = pd.DataFrame(output2_1)
output2_2 = pd.DataFrame(output2_2)
output2 = pd.merge(output2_1, output2_2, how='outer', left_on=["mergedClusterID"], right_on=["mergedClusterID"])

output1.to_csv(input1 + '.allCluster.tsv', sep="\t", header=True, index=False)
output2.to_csv(input1 + '.allCluster.sum.tsv', sep="\t", header=True, index=True)

