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

# 读取新版本的breakpointINPUT.tsv文件（无表头，8列）
# 列顺序：sampleID, Cluster_No, disFile, splitFile, wgsBAM, chr, start, end
column_names = ['sampleID', 'Cluster_No', 'disFile', 'splitFile', 'wgsBAM', 'chr', 'start', 'end']
df0 = pd.read_csv(input1, sep="\t", engine='python', header=None, names=column_names)

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
	myclusters['position'] = myclusters['position'].astype(int)
	sub_cluster = cluster(myclusters['position'].tolist(), maxgap=1000)
	df1_temp = pd.DataFrame.from_records(sub_cluster)
	df2 = df1_temp.stack()
	df2 = df2.reset_index()
	df2.columns = ['GroupID','Index','POS']
	df2['CHR'] = clusterID
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

