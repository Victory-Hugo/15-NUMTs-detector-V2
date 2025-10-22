#!/bin/bash

# 使用示例 - 0-打包.py / 0-打包_LDY.py
# 这个脚本会自动检测文件版本，使用对应的Python脚本处理

PYTHON3="/home/luolintao/miniconda3/envs/pyg/bin/python3" #* 请根据实际情况修改为你的Python3路径
PYTHON_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布/script/0-打包.py"
PYTHON_SCRIPT_LDY="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布/script/0-打包_LDY.py" #* 新版本脚本 (LDY)
PREPARE_CIRCOS_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布/script/1-准备NUMTs可视化.py" #* 用于转换为circos输入格式的脚本路径
INPUT_CSV_DIR="/mnt/f/13_SLE_NUMT/share/home/grp-wangyf/luolintao/1-luolintao/2-NUMTs/1-整理/output" #* 里面存放了所有的tsv文件
OUTPUT_DIR="/mnt/f/13_SLE_NUMT/SLE_完成_output" #* 输出目录，存放合并后的tsv文件
TAR_GZ_PATH="/mnt/f/13_SLE_NUMT/SLE完成.tar.gz" #* 输出的tar.gz压缩文件路径

# ============================================================================
# 版本检测函数
# ============================================================================

detect_file_version() {
    """
    检测输入目录中的文件版本
    返回值: 
      - 0 表示旧版本 (*.Breakpoints.tsv)
      - 1 表示新版本 (*ConfidentBreakpoints.tsv, *AllBreakpoints.tsv)
      - 2 表示混合版本 (两种格式都有)
    """
    local input_dir="$1"
    
    # 检查旧版本文件
    old_breakpoints=$(find "$input_dir" -maxdepth 1 -type f -name "*.Breakpoints.tsv" 2>/dev/null | wc -l)
    
    # 检查新版本文件
    new_confident=$(find "$input_dir" -maxdepth 1 -type f -name "*ConfidentBreakpoints.tsv" 2>/dev/null | wc -l)
    new_all=$(find "$input_dir" -maxdepth 1 -type f -name "*AllBreakpoints.tsv" 2>/dev/null | wc -l)
    
    # 判断版本
    if [ "$new_confident" -gt 0 ] || [ "$new_all" -gt 0 ]; then
        if [ "$old_breakpoints" -gt 0 ]; then
            echo 2  # 混合版本
        else
            echo 1  # 新版本
        fi
    else
        echo 0  # 旧版本
    fi
}

# ============================================================================
# 主程序
# ============================================================================

echo "==============================================="
echo "   TSV文件合并与压缩工具 (智能版本检测)"
echo "==============================================="
echo ""

# 检查输入目录是否存在
if [ ! -d "$INPUT_CSV_DIR" ]; then
    echo "错误: 输入目录不存在: $INPUT_CSV_DIR"
    exit 1
fi

# 检测文件版本
echo "正在检测文件版本..."
VERSION=$(detect_file_version "$INPUT_CSV_DIR")

echo ""
echo "==============================================="

case $VERSION in
    0)
        echo "检测结果: 旧版本文件"
        echo "文件格式: *.Breakpoints.tsv, *.mt.disc.sam.breakpointINPUT.tsv 等"
        echo "使用脚本: 0-打包.py"
        echo "==============================================="
        echo ""
        SELECTED_SCRIPT="$PYTHON_SCRIPT"
        SCRIPT_NAME="0-打包.py"
        ;;
    1)
        echo "检测结果: 新版本文件 (LDY版本)"
        echo "文件格式: *ConfidentBreakpoints.tsv, *AllBreakpoints.tsv 等"
        echo "使用脚本: 0-打包_LDY.py"
        echo "==============================================="
        echo ""
        SELECTED_SCRIPT="$PYTHON_SCRIPT_LDY"
        SCRIPT_NAME="0-打包_LDY.py"
        ;;
    2)
        echo "检测结果: 混合版本文件"
        echo "警告: 目录中同时包含旧版本和新版本的文件!"
        echo "旧版本: *.Breakpoints.tsv"
        echo "新版本: *ConfidentBreakpoints.tsv, *AllBreakpoints.tsv"
        echo ""
        echo "请选择处理方式:"
        echo "1) 使用旧版本脚本 (0-打包.py) - 仅处理旧格式"
        echo "2) 使用新版本脚本 (0-打包_LDY.py) - 仅处理新格式"
        echo "3) 取消处理"
        echo ""
        read -p "请输入选择 (1/2/3): " choice
        echo "==============================================="
        echo ""
        
        case $choice in
            1)
                SELECTED_SCRIPT="$PYTHON_SCRIPT"
                SCRIPT_NAME="0-打包.py"
                echo "已选择: 使用旧版本脚本 (0-打包.py)"
                echo ""
                ;;
            2)
                SELECTED_SCRIPT="$PYTHON_SCRIPT_LDY"
                SCRIPT_NAME="0-打包_LDY.py"
                echo "已选择: 使用新版本脚本 (0-打包_LDY.py)"
                echo ""
                ;;
            3)
                echo "已取消处理"
                exit 0
                ;;
            *)
                echo "无效选择，已取消"
                exit 1
                ;;
        esac
        ;;
esac

# 检查选择的脚本是否存在
if [ ! -f "$SELECTED_SCRIPT" ]; then
    echo "错误: 选择的脚本不存在: $SELECTED_SCRIPT"
    exit 1
fi

# 执行选择的Python脚本
echo "开始执行 $SCRIPT_NAME ..."
echo ""

${PYTHON3} \
    ${SELECTED_SCRIPT} \
    ${INPUT_CSV_DIR} \
    ${OUTPUT_DIR} \
    ${TAR_GZ_PATH}

# 检查Python脚本的执行结果
if [ $? -ne 0 ]; then
    echo ""
    echo "错误: Python脚本执行失败"
    exit 1
fi

# 执行circos准备脚本
echo ""
echo "准备circos输入文件..."

# 根据版本选择合适的输出文件
if [ "$VERSION" -eq 1 ] || ([ "$VERSION" -eq 2 ] && [ "$SCRIPT_NAME" = "0-打包_LDY.py" ]); then
    # 新版本可能没有cluster.summary，检查cluster文件是否存在
    CLUSTER_FILE="${OUTPUT_DIR}/all_individuals_mt.disc.sam.cluster.tsv"
    if [ ! -f "$CLUSTER_FILE" ]; then
        echo "警告: 未找到聚类文件，跳过circos准备步骤"
    else
        ${PYTHON3} \
            ${PREPARE_CIRCOS_SCRIPT} \
            ${CLUSTER_FILE} \
            ${OUTPUT_DIR}/circos.txt
    fi
else
    # 旧版本
    ${PYTHON3} \
        ${PREPARE_CIRCOS_SCRIPT} \
        ${OUTPUT_DIR}/all_individuals_mt.disc.sam.cluster.tsv \
        ${OUTPUT_DIR}/circos.txt
fi

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