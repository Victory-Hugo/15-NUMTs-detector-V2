#!/bin/bash

################################################################################
## 准备基因启动子区域BED文件 - 修正版
## 功能：从GENCODE注释中提取蛋白编码基因的启动子区域（TSS上游2kb，下游500bp）
################################################################################

set -euo pipefail

# 设置工作目录
WORK_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-富集/meta"
cd "$WORK_DIR"

echo "=== 准备基因启动子区域BED文件 ==="

# 检查输入文件
if [[ ! -f "gencode.v44.annotation.gtf.gz" ]]; then
    echo "错误：未找到 gencode.v44.annotation.gtf.gz 文件"
    exit 1
fi

echo "1. 提取蛋白编码基因的转录起始位点..."

# 提取蛋白编码基因信息
zcat gencode.v44.annotation.gtf.gz | \
awk '$3=="gene"' | \
grep 'gene_type "protein_coding"' | \
awk -F'\t' '{
    # 解析属性字段以获取基因名称
    gene_name = "unknown";
    if (match($9, /gene_name "([^"]+)"/, arr)) {
        gene_name = arr[1];
    }
    
    # 计算启动子区域（TSS上游2kb，下游500bp）
    if ($7 == "+") {
        # 正链：TSS = start
        start = $4 - 2000;
        end = $4 + 500;
    } else {
        # 负链：TSS = end  
        start = $5 - 500;
        end = $5 + 2000;
    }
    
    # 确保坐标不为负数
    if (start < 1) start = 1;
    
    # 输出BED格式：chr, start, end, gene_name, score, strand
    print $1 "\t" start "\t" end "\t" gene_name "\t" "1000" "\t" $7;
}' > promoters_with_names.bed

echo "2. 检查提取结果..."
lines=$(wc -l < promoters_with_names.bed)
echo "   提取到 $lines 个基因的启动子区域"

if [[ $lines -eq 0 ]]; then
    echo "错误：未提取到任何基因信息，请检查GTF文件格式"
    exit 1
fi

echo "3. 清理和标准化染色体名称..."

# 只保留主要染色体，排除scaffold和patch
awk '$1 ~ /^chr[0-9XY]+$/ || $1 ~ /^chr[12]?[0-9]$/' promoters_with_names.bed > promoters_main_chr.bed

echo "4. 按位置排序..."
sort -k1,1 -k2,2n promoters_main_chr.bed > promoters_sorted.bed

echo "5. 去除重叠区域..."

# 使用bedtools merge去除重叠（如果可用）
if command -v bedtools &> /dev/null; then
    echo "   使用bedtools合并重叠区域..."
    # 只使用前3列进行合并
    cut -f1-3 promoters_sorted.bed | bedtools merge -i - > promoters_merged.bed
    mv promoters_merged.bed promoters_final.bed
else
    echo "   bedtools未安装，使用awk去重..."
    # 简单的重叠去除，只保留前3列
    awk '{
        if (prev_chr == $1 && prev_end > $2) {
            # 重叠区域，合并
            prev_end = ($3 > prev_end) ? $3 : prev_end;
        } else {
            # 非重叠区域，输出前一个区域
            if (NR > 1) print prev_chr "\t" prev_start "\t" prev_end;
            prev_chr = $1; prev_start = $2; prev_end = $3;
        }
    } END {
        if (NR > 0) print prev_chr "\t" prev_start "\t" prev_end;
    }' promoters_sorted.bed > promoters_final.bed
fi

echo "6. 生成统计报告..."

# 统计信息
if [[ -s promoters_final.bed ]]; then
    total_regions=$(wc -l < promoters_final.bed)
    total_length=$(awk '{sum += ($3-$2)} END {print sum}' promoters_final.bed)
    avg_length=$(awk '{sum += ($3-$2); count++} END {printf "%.0f", sum/count}' promoters_final.bed)
    chr_count=$(cut -f1 promoters_final.bed | sort -u | wc -l)
    
    # 计算基因组覆盖度（基于GRCh38）
    genome_size=3088269832  # GRCh38 总长度
    coverage_percent=$(awk -v total=$total_length -v genome=$genome_size 'BEGIN {printf "%.3f", total/genome*100}')
    
    echo "=== 启动子区域BED文件生成完成 ==="
    echo "输出文件: promoters_final.bed"
    echo ""
    echo "=== 统计信息 ==="
    echo "区域总数: $total_regions"
    echo "染色体数: $chr_count"
    echo "总覆盖长度: $total_length bp"
    echo "平均区域长度: $avg_length bp"
    echo "基因组覆盖度: $coverage_percent%"
    echo ""
    
    # 显示染色体分布
    echo "=== 各染色体区域分布 ==="
    cut -f1 promoters_final.bed | sort | uniq -c | sort -k1,1nr | head -10
    
    echo ""
    echo "=== 前10个区域示例 ==="
    head -10 promoters_final.bed
    
    # 创建最终的标准BED文件（3列格式）
    cp promoters_final.bed promoters.bed
    
    echo ""
    echo "=== 最终输出文件 ==="
    echo "启动子区域BED文件: promoters.bed"
    
else
    echo "错误：最终文件为空"
    exit 1
fi

# 清理临时文件
rm -f promoters_with_names.bed promoters_main_chr.bed promoters_sorted.bed

echo "完成！"
echo ""
echo "=== 使用建议 ==="
echo "此BED文件适合用于NUMTs富集分析："
echo "- 区域类型: 蛋白编码基因启动子（TSS±2kb/500bp）"
echo "- 建议用于研究NUMTs对基因转录调控的影响"
echo "- 如果覆盖度太小，可以考虑扩大flanking区域"