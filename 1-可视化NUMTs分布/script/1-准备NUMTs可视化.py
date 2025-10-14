#!/usr/bin/env python3
import sys
import pandas as pd

if len(sys.argv) != 3:
    print("用法: python script.py <输入文件> <输出文件>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# 加载数据
data = pd.read_csv(input_file, sep='\t')

# 检查 'Cluster_ID' 是否缺失
missing_count = data['Cluster_ID'].isna().sum()
if missing_count > 0:
    print(f"警告：'Cluster_ID' 列中有 {missing_count} 个缺失值，将删除这些行。")
    data = data.dropna(subset=['Cluster_ID'])

# 提取列
cluster_id_col = data['Cluster_ID']

# 转换函数
def convert_to_circos_format(cluster_id):
    try:
        cluster_id = str(cluster_id)
        parts = cluster_id.split('_')
        if len(parts) < 6:
            raise ValueError(f"Unexpected format: {cluster_id}")
        hs1 = parts[0].replace("('", "").replace("',)", "")
        start1, end1, startM, endM = parts[1], parts[2], parts[4], parts[5]
        return f"{hs1} {start1} {end1} hsM {startM} {endM}"
    except Exception as e:
        print(f"Error processing Cluster_ID '{cluster_id}': {e}")
        return None

# 应用转换
circos_format_data = cluster_id_col.apply(convert_to_circos_format).dropna()

# 转为 DataFrame
circos_format_df = circos_format_data.str.split(expand=True)
circos_format_df.columns = ['hs1', 'start1', 'end1', 'hsM', 'startM', 'endM']

# 添加 'hs' 前缀
def add_hs_prefix(value):
    return f"hs{value}" if value.isdigit() else value

circos_format_df['hs1'] = circos_format_df['hs1'].apply(add_hs_prefix)

# 替换 'chr' -> 'hs' 并去重
circos_format_df = circos_format_df.replace({'chr': 'hs'}, regex=True).drop_duplicates()

# 保存结果
circos_format_df.to_csv(output_file, sep=' ', index=False, header=False)

print(f"文件已保存至 {output_file}")
