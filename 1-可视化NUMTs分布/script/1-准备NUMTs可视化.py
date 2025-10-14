import pandas as pd

# 加载上传的文件
file_path = '/mnt/f/12_甘肃藏族_NUMT/1-所有的NUMTs/data/all_individuals_mt.disc.sam.cluster.tsv'
data = pd.read_csv(file_path, sep='\t')

# 检查 'Cluster_ID' 列是否存在缺失值
missing_count = data['Cluster_ID'].isna().sum()
if missing_count > 0:
    print(f"警告：'Cluster_ID' 列中有 {missing_count} 个缺失值。将删除这些行。")
    data = data.dropna(subset=['Cluster_ID'])

# 提取相关列
cluster_id_col = data['Cluster_ID']

# 定义函数，将 Cluster_ID 列转换为所需的格式
def convert_to_circos_format(cluster_id):
    try:
        # 确保 cluster_id 是字符串
        cluster_id = str(cluster_id)
        parts = cluster_id.split('_')
        if len(parts) < 6:
            raise ValueError(f"Unexpected format: {cluster_id}")
        
        hs1 = parts[0].replace("('", "").replace("',)", "")  # 提取第一个字段并清理
        start1 = parts[1]  # 提取第二个字段
        end1 = parts[2]  # 提取第三个字段
        hsM = 'hsM'  # 固定值 'hsM'
        startM = parts[4]  # 提取第五个字段
        endM = parts[5]  # 提取第六个字段
        return f"{hs1} {start1} {end1} {hsM} {startM} {endM}"
    except Exception as e:
        print(f"Error processing Cluster_ID '{cluster_id}': {e}")
        return None  # 或者返回一个默认值，例如 "NA NA NA hsM NA NA"

# 将函数应用到 Cluster_ID 列
circos_format_data = cluster_id_col.apply(convert_to_circos_format)

# 删除转换失败的行
circos_format_data = circos_format_data.dropna()

# 将结果转换为 DataFrame
circos_format_df = circos_format_data.str.split(expand=True)
circos_format_df.columns = ['hs1', 'start1', 'end1', 'hsM', 'startM', 'endM']

# 定义函数，如果第一列是数字，则添加 'hs' 前缀
def add_hs_prefix(value):
    return f"hs{value}" if value.isdigit() else value

# 将函数应用到第一列
circos_format_df['hs1'] = circos_format_df['hs1'].apply(add_hs_prefix)

# 去除文件头并将 'chr' 替换为 'hs'
circos_format_df_no_header = circos_format_df.replace({'chr': 'hs'}, regex=True)

# 对 DataFrame 去重，保留完全相同的行中的一个
circos_format_df_no_header = circos_format_df_no_header.drop_duplicates()

# 将结果保存到文件中，不包含表头
output_file_path = '/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布/output/circos.txt'
circos_format_df_no_header.to_csv(output_file_path, sep=' ', index=False, header=False)

print(f"文件已保存至 {output_file_path}")
