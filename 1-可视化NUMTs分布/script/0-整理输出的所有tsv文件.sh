#!/bin/bash

# 使用示例 - 0-打包.py
# 这个脚本展示了如何使用 0-打包.py 工具

PYTHON3="/home/luolintao/miniconda3/envs/pyg/bin/python3" #* 请根据实际情况修改为你的Python3路径
PYTHON_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/2-汇总NUMTs分布/script/0-打包.py"
PREPARE_CIRCOS_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布/script/1-准备NUMTs可视化.py" #* 用于转换为circos输入格式的脚本路径
INPUT_CSV_DIR="/mnt/c/Users/Administrator/Desktop/download" #* 里面存放了所有的tsv文件: *.mt.disc.sam.breakpointINPUT.tsv / *.mt.disc.sam.cluster.summary.tsv / *.mt.disc.sam.cluster.tsv / *.Breakpoints.tsv
OUTPUT_DIR="/mnt/c/Users/Administrator/Desktop/藏族_甘肃完成_output" #* 输出目录，存放合并后的tsv文件
TAR_GZ_PATH="/mnt/c/Users/Administrator/Desktop/藏族_甘肃完成.tar.gz" #* 输出的tar.gz压缩文件路径

${PYTHON3} \
    ${PYTHON_SCRIPT} \
    ${INPUT_CSV_DIR} \
    ${OUTPUT_DIR} \
    ${TAR_GZ_PATH}

${PYTHON3} \
    ${PREPARE_CIRCOS_SCRIPT} \
    ${OUTPUT_DIR}/all_individuals_mt.disc.sam.cluster.tsv \
    ${OUTPUT_DIR}/circos.txt

echo "==============================================="
echo "   TSV文件合并与压缩工具使用示例"
echo "==============================================="

# 示例1: 基本用法
echo "示例1: 基本用法"
echo "python 1-打包.py <输入目录> <输出目录> <压缩文件路径>"
echo ""

# 示例2: 处理藏族完成1数据
echo "示例2: 处理藏族完成1数据"
echo "python 1-打包.py ./藏族完成1/tsv ./藏族完成1/merged_output ./藏族完成1/original_files.tar.gz"
echo ""

# 示例3: 处理其他民族数据
echo "示例3: 处理其他民族数据"
echo "python 1-打包.py ./甘肃完成1/tsv ./甘肃完成1/merged_output ./甘肃完成1/original_files.tar.gz"
echo ""

# 功能说明
echo "==============================================="
echo "功能说明:"
echo "1. 自动合并4种类型的TSV文件:"
echo "   - *.Breakpoints.tsv"
echo "   - *.mt.disc.sam.breakpointINPUT.tsv"
echo "   - *.mt.disc.sam.cluster.summary.tsv"
echo "   - *.mt.disc.sam.cluster.tsv"
echo ""
echo "2. 智能处理空文件，避免解析错误"
echo ""
echo "3. 压缩原始文件为 tar.gz 格式"
echo ""
echo "4. 验证压缩文件完整性"
echo ""
echo "5. 安全删除原始文件（仅在验证通过后）"
echo ""
echo "6. 彩色输出，便于查看进度"
echo ""
echo "7. 详细的统计报告和磁盘空间释放信息"
echo ""

# 安全特性说明
echo "==============================================="
echo "安全特性:"
echo "1. 压缩文件完整性验证"
echo "   - 检查压缩文件是否可以正常打开"
echo "   - 抽样验证文件内容可读性"
echo "   - 只有验证通过才会删除原始文件"
echo ""
echo "2. 详细的操作日志"
echo "   - 记录每个步骤的执行状态"
echo "   - 显示文件数量和大小统计"
echo "   - 失败时立即停止，保护数据安全"
echo ""
echo "3. 磁盘空间管理"
echo "   - 显示释放的磁盘空间大小"
echo "   - 统计处理的文件数量"
echo "==============================================="