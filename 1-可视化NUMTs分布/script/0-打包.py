#!/usr/bin/env python3
"""
TSV文件合并与压缩工具 - 1-打包.py

该脚本将合并指定目录下的TSV文件并压缩原始文件：
1. *.Breakpoints.tsv - 所有个体的可信断点文件 (输出为: all_individuals_ConfidentBreakpoints.tsv)
2. *.mt.disc.sam.breakpointINPUT.tsv - 所有个体的断点输入文件  
3. *.mt.disc.sam.cluster.summary.tsv - 所有个体的聚类摘要文件
4. *.mt.disc.sam.cluster.tsv - 所有个体的聚类详细文件

用法:
python 1-打包.py <输入目录> <输出目录> <压缩文件路径>

示例:
python 1-打包.py ./完成1/tsv ./完成1/merged_output ./完成1/original_files.tar.gz
"""

import os
import sys
import glob
import pandas as pd
import tarfile
import gzip
import shutil
from pathlib import Path
from datetime import datetime

# 定义颜色代码
class Colors:
    """终端颜色定义"""
    # 基本颜色
    RED = '\033[91m'      # 错误信息
    GREEN = '\033[92m'    # 成功信息
    YELLOW = '\033[93m'   # 警告信息
    BLUE = '\033[94m'     # 信息提示
    MAGENTA = '\033[95m'  # 重要信息
    CYAN = '\033[96m'     # 标题
    WHITE = '\033[97m'    # 普通文本
    
    # 样式
    BOLD = '\033[1m'      # 粗体
    UNDERLINE = '\033[4m' # 下划线
    
    # 重置
    RESET = '\033[0m'     # 重置所有样式

def print_colored(message, color=Colors.WHITE):
    """打印彩色文本"""
    print(f"{color}{message}{Colors.RESET}")

def print_header(title):
    """打印标题"""
    print_colored("=" * 50, Colors.CYAN)
    print_colored(f"{Colors.BOLD}{title}", Colors.CYAN)
    print_colored("=" * 50, Colors.CYAN)

def print_success(message):
    """打印成功信息"""
    print_colored(f"✓ {message}", Colors.GREEN)

def print_warning(message):
    """打印警告信息"""
    print_colored(f"⚠ {message}", Colors.YELLOW)

def print_error(message):
    """打印错误信息"""
    print_colored(f"✗ {message}", Colors.RED)

def print_info(message):
    """打印信息"""
    print_colored(f"ℹ {message}", Colors.BLUE)

def print_progress(message):
    """打印进度信息"""
    print_colored(f"→ {message}", Colors.MAGENTA)

def is_file_empty_or_header_only(file_path, has_header=True):
    """
    检查文件是否为空或只有标题行
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # 去除空行
        non_empty_lines = [line.strip() for line in lines if line.strip()]
        
        if len(non_empty_lines) == 0:
            return True
        
        if has_header and len(non_empty_lines) == 1:
            return True
            
        if not has_header and len(non_empty_lines) == 0:
            return True
            
        return False
    except Exception:
        return True

def verify_compressed_file(compress_path):
    """
    验证压缩文件的完整性
    
    Args:
        compress_path (str): 压缩文件路径
        
    Returns:
        bool: 验证成功返回True，失败返回False
    """
    print_progress("验证压缩文件完整性...")
    
    try:
        # 尝试打开并列出压缩文件内容
        with tarfile.open(compress_path, "r:gz") as tar:
            # 获取压缩文件中的所有文件
            members = tar.getmembers()
            print_info(f"压缩文件包含 {len(members)} 个条目")
            
            # 验证每个文件是否可以正常读取
            for member in members[:10]:  # 验证前10个文件作为抽样检查
                if member.isfile():
                    try:
                        tar.extractfile(member).read(1024)  # 尝试读取前1KB
                    except Exception as e:
                        print_error(f"验证文件 {member.name} 时出错：{e}")
                        return False
            
            print_success("压缩文件完整性验证通过")
            return True
            
    except Exception as e:
        print_error(f"压缩文件验证失败：{e}")
        return False

def delete_original_files(input_dir):
    """
    删除原始文件目录
    
    Args:
        input_dir (str): 要删除的目录路径
        
    Returns:
        bool: 删除成功返回True，失败返回False
    """
    print_progress("删除原始文件...")
    
    try:
        # 计算原始目录大小
        total_size = 0
        file_count = 0
        for root, dirs, files in os.walk(input_dir):
            for file in files:
                file_path = os.path.join(root, file)
                if os.path.exists(file_path):
                    total_size += os.path.getsize(file_path)
                    file_count += 1
        
        print_info(f"准备删除 {file_count} 个文件，总大小 {total_size:,} bytes ({total_size/1024/1024:.2f} MB)")
        
        # 删除整个目录
        shutil.rmtree(input_dir)
        
        print_success(f"成功删除原始文件目录：{input_dir}")
        print_success(f"释放磁盘空间：{total_size/1024/1024:.2f} MB")
        
        return True
        
    except Exception as e:
        print_error(f"删除原始文件失败：{e}")
        return False

def compress_original_files(input_dir, compress_path):
    """
    压缩原始文件目录为tar.gz格式
    
    Args:
        input_dir (str): 输入目录路径
        compress_path (str): 压缩文件输出路径
    """
    
    print_progress("开始压缩原始文件...")
    
    try:
        # 确保压缩文件的目录存在
        compress_dir = os.path.dirname(compress_path)
        if compress_dir:
            os.makedirs(compress_dir, exist_ok=True)
        
        # 创建tar.gz压缩文件
        with tarfile.open(compress_path, "w:gz") as tar:
            # 添加整个输入目录到压缩文件
            tar.add(input_dir, arcname=os.path.basename(input_dir))
        
        # 获取压缩文件大小
        compress_size = os.path.getsize(compress_path)
        
        print_success(f"压缩完成：{compress_path}")
        print_info(f"压缩文件大小：{compress_size:,} bytes ({compress_size/1024/1024:.2f} MB)")
        
        return True
        
    except Exception as e:
        print_error(f"压缩失败：{e}")
        return False

def merge_files_by_type(input_dir, output_dir):
    """
    按文件类型合并TSV文件
    
    Args:
        input_dir (str): 输入目录路径
        output_dir (str): 输出目录路径
    """
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 定义文件类型模式和输出文件名
    file_patterns = {
        'breakpoints': {
            'pattern': '*.Breakpoints.tsv',
            'output': 'all_individuals_ConfidentBreakpoints.tsv',
            'has_header': False
        },
        'breakpoint_input': {
            'pattern': '*.mt.disc.sam.breakpointINPUT.tsv',
            'output': 'all_individuals_mt.disc.sam.breakpointINPUT.tsv',
            'has_header': False
        },
        'cluster_summary': {
            'pattern': '*.mt.disc.sam.cluster.summary.tsv',
            'output': 'all_individuals_mt.disc.sam.cluster.summary.tsv',
            'has_header': False
        },
        'cluster': {
            'pattern': '*.mt.disc.sam.cluster.tsv',
            'output': 'all_individuals_mt.disc.sam.cluster.tsv',
            'has_header': True
        }
    }
    
    # 处理每种文件类型
    for file_type, config in file_patterns.items():
        print_progress(f"处理 {file_type} 文件类型...")
        
        # 搜索匹配的文件
        pattern = os.path.join(input_dir, config['pattern'])
        files = glob.glob(pattern)
        
        if not files:
            print_warning(f"未找到匹配模式 {config['pattern']} 的文件")
            continue
            
        print_info(f"找到 {len(files)} 个文件")
        
        # 过滤掉空文件
        valid_files = []
        for file_path in sorted(files):
            if not is_file_empty_or_header_only(file_path, config['has_header']):
                valid_files.append(file_path)
            else:
                print_warning(f"跳过空文件: {os.path.basename(file_path)}")
        
        if not valid_files:
            print_warning("所有文件都是空文件")
            continue
            
        print_info(f"有效文件数量: {len(valid_files)}")
        
        # 合并文件
        merged_data = []
        
        for i, file_path in enumerate(valid_files):
            print_colored(f"  处理文件 {i+1}/{len(valid_files)}: {os.path.basename(file_path)}", Colors.WHITE)
            
            try:
                # 读取文件
                if config['has_header']:
                    df = pd.read_csv(file_path, sep='\t')
                else:
                    df = pd.read_csv(file_path, sep='\t', header=None)
                
                if len(df) == 0:
                    print_warning(f"文件 {os.path.basename(file_path)} 数据为空")
                    continue
                
                merged_data.append(df)
                
            except Exception as e:
                print_error(f"读取文件 {os.path.basename(file_path)} 时出错: {e}")
                continue
        
        # 合并所有数据
        if merged_data:
            try:
                combined_df = pd.concat(merged_data, ignore_index=True)
                
                # 输出文件路径
                output_file = os.path.join(output_dir, config['output'])
                
                # 保存合并后的文件
                combined_df.to_csv(output_file, sep='\t', index=False, header=config['has_header'])
                
                print_success(f"成功合并 {len(merged_data)} 个文件，共 {len(combined_df):,} 行数据")
                print_info(f"输出文件: {output_file}")
                
            except Exception as e:
                print_error(f"合并文件时出错: {e}")
        else:
            print_warning("没有有效数据可合并")

def show_output_summary(output_dir):
    """显示输出文件摘要"""
    if not os.path.exists(output_dir):
        print_warning("输出目录不存在")
        return
    
    print_header("输出文件摘要")
    
    files = [f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))]
    
    if not files:
        print_warning("输出目录为空")
        return
    
    total_size = 0
    for file in sorted(files):
        file_path = os.path.join(output_dir, file)
        file_size = os.path.getsize(file_path)
        total_size += file_size
        
        # 统计行数
        try:
            with open(file_path, 'r') as f:
                line_count = sum(1 for _ in f)
            print_info(f"{file}: {file_size:,} bytes, {line_count:,} 行")
        except:
            print_info(f"{file}: {file_size:,} bytes")
    
    print_success(f"总大小: {total_size:,} bytes ({total_size/1024/1024:.2f} MB)")

def print_usage():
    """打印使用说明"""
    print_header("使用说明")
    print_colored("用法:", Colors.BOLD)
    print_colored("  python 1-打包.py <输入目录> <输出目录> <压缩文件路径>", Colors.WHITE)
    print_colored("\n参数说明:", Colors.BOLD)
    print_colored("  输入目录     - 包含TSV文件的目录路径", Colors.WHITE)
    print_colored("  输出目录     - 合并文件的输出目录", Colors.WHITE)
    print_colored("  压缩文件路径 - 原始文件压缩包的保存路径", Colors.WHITE)
    print_colored("\n示例:", Colors.BOLD)
    print_colored("  python 1-打包.py ./藏族完成1/tsv ./藏族完成1/merged_output ./藏族完成1/original_files.tar.gz", Colors.WHITE)

def main():
    """主函数"""
    
    # 检查命令行参数
    if len(sys.argv) != 4:
        print_error("参数数量不正确")
        print_usage()
        sys.exit(1)
    
    input_directory = sys.argv[1]
    output_directory = sys.argv[2] 
    compress_file_path = sys.argv[3]
    
    # 转换为绝对路径
    input_directory = os.path.abspath(input_directory)
    output_directory = os.path.abspath(output_directory)
    compress_file_path = os.path.abspath(compress_file_path)
    
    print_header("TSV文件合并与压缩工具")
    print_info(f"输入目录: {input_directory}")
    print_info(f"输出目录: {output_directory}")
    print_info(f"压缩文件: {compress_file_path}")
    print_info(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 检查输入目录是否存在
    if not os.path.exists(input_directory):
        print_error(f"输入目录不存在: {input_directory}")
        sys.exit(1)
    
    # 步骤1: 执行文件合并
    print_header("步骤1: 合并TSV文件")
    try:
        merge_files_by_type(input_directory, output_directory)
        print_success("文件合并完成")
    except Exception as e:
        print_error(f"文件合并失败: {e}")
        sys.exit(1)
    
    # 步骤2: 压缩原始文件
    print_header("步骤2: 压缩原始文件")
    try:
        success = compress_original_files(input_directory, compress_file_path)
        if not success:
            print_error("压缩失败")
            sys.exit(1)
        print_success("文件压缩完成")
    except Exception as e:
        print_error(f"压缩过程出错: {e}")
        sys.exit(1)
    
    # 步骤3: 验证压缩文件
    print_header("步骤3: 验证压缩文件")
    try:
        verify_success = verify_compressed_file(compress_file_path)
        if not verify_success:
            print_error("压缩文件验证失败，不会删除原始文件")
            sys.exit(1)
        print_success("压缩文件验证成功")
    except Exception as e:
        print_error(f"验证过程出错: {e}")
        sys.exit(1)
    
    # 步骤4: 删除原始文件
    print_header("步骤4: 删除原始文件")
    try:
        delete_success = delete_original_files(input_directory)
        if not delete_success:
            print_error("删除原始文件失败")
            sys.exit(1)
        print_success("原始文件删除完成")
    except Exception as e:
        print_error(f"删除过程出错: {e}")
        sys.exit(1)
    
    # 显示结果摘要
    show_output_summary(output_directory)
    
    print_header("任务完成")
    print_success("所有操作已成功完成!")
    print_info(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()