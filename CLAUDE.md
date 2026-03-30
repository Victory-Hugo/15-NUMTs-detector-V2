# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

NUMTs Detector V2 是一个用于识别和分析 **NUMTs（Nuclear Mitochondrial DNA segments，核线粒体 DNA 片段）** 的生物信息学流水线。项目采用模块化设计，支持从 BAM 文件检测 NUMTs、变异分析、统计描述和可视化。

## 核心架构

### 目录结构规范

本项目遵循科研流水线标准目录结构：

- `conf/` - 配置文件（Config.yaml）
- `pipe/` - 主控流程脚本（bash）
- `script/` - 流程调用的辅助脚本
- `python/` - Python 分析模块
- `src/` - 核心源码或可复用组件
- `meta/` - 元数据文件（样本信息表等）
- `tmp/` - 临时文件和日志
- `output/` - 结果输出目录

### 主要分析模块

项目按功能划分为多个独立模块，每个模块都是独立的子目录：

1. **0-GRCH37-liftover** - 坐标系转换（GRCh37/hg19）
2. **1-1-获得NUMTs分布** - 从 BAM 文件检测 NUMTs 候选区域和断点
3. **1-5-Vardetection** - NUMTs 变异检测（与人类/黑猩猩 mtDNA 比对）
4. **2-整理tsv文件** - 结果文件整理和格式化
5. **3-可视化NUMTs分布** - Circos 和 Ideogram 可视化
6. **4-统计性描述** - 频率统计、聚类分析、图表生成
7. **5-富集** - 功能富集分析
8. **9-参考NUMTs** - 参考 NUMTs 数据库

## 关键技术栈

### 生物信息学工具
- **samtools** - BAM/SAM 文件处理
- **bwa-mem** - 序列比对
- **BLAT** - 快速序列比对（用于 mtDNA 比对）
- **CAP3** - 序列拼接（可选）

### Python 依赖
- pandas - 数据处理
- numpy - 数值计算
- matplotlib/seaborn - 可视化
- pysam - BAM/SAM 文件读写（部分模块）

### 运行环境
- Conda 环境管理
- Bash shell 脚本
- 支持多线程并行处理

## 常用命令

### 1. NUMTs 检测（单样本）

```bash
cd 1-1-获得NUMTs分布/pipe
bash 0_NUMTs_detection_2.0_单样本.sh
```

需要在配置文件中设置：
- `INPUT_BAM` - 输入 BAM 文件路径
- `OUTPUT_DIR` - 输出目录
- `REF_GRCh` - 参考基因组路径
- `THREADS` - 线程数

### 2. NUMTs 检测（多样本并行）

```bash
cd 1-1-获得NUMTs分布/pipe
bash 0_NUMTs_detection_2.0_多样本_并行.sh
```

### 3. 统计性描述

```bash
cd 4-统计性描述/pipe
bash 1-run-statistical-description.sh
```

需要配置 `conf/Config.yaml`：
- `input.breakpoint_tsv_gz` - 断点结果文件
- `input.distinct_numt_cluster_input_tsv_gz` - 聚类输入文件
- `input.meta_tsv` - 样本元数据
- `analysis.distinct_numt_min_sample_supports` - 最小样本支持数阈值

### 4. 变异检测

```bash
cd 1-5-Vardetection/script
# 查看具体脚本的使用说明
```

## 配置文件说明

### Config.yaml 结构

配置文件采用 YAML 格式，主要包含以下部分：

```yaml
input:           # 输入文件路径
meta:            # 元数据配置（列名、质控标准）
reference:       # 参考数据路径
analysis:        # 分析参数（阈值、窗口大小等）
runtime:         # 运行时配置（conda 环境、线程数）
output:          # 输出目录结构
```

**重要原则**：
- 所有软件路径必须在配置文件中定义，禁止硬编码
- 包括 conda 可执行文件路径、conda 环境名/前缀、各类工具路径
- 模块脚本通过命名 CLI 参数接收配置，不直接读取配置文件

## 输出文件说明

### NUMTs 检测主要输出

1. **聚类结果**
   - `*.cluster.tsv` - 候选 NUMT 区域总表（带表头）
   - `*.cluster.summary.tsv` - 聚类汇总统计
   - `*.breakpointINPUT.tsv` - 断点分析输入区域清单（无表头）

2. **断点结果**
   - `*.AllBreakpoints.tsv` - 所有检测到的断点
   - `*.ConfidentBreakpoints.tsv` - 高可信度断点（核侧和线粒体侧证据互相对应）

3. **比对文件**
   - `*.mt.disc.sam` - discordant reads
   - `*.mt.split.sam` - split reads
   - `*_all_regions.psl` - BLAT 比对结果

### 变异检测输出

位于 `1-5-Vardetection/output/<SAMPLE_ID>/`：

- `alnHuman/*.numtVarFilterPos.tsv` - 可信 NUMT 变异集合
- `alnHumanChimp/*.numtDhumanChimp.sum.tsv` - NUMT 年龄推断结果

## 开发注意事项

### 代码规范

1. **模块化设计**
   - 主控脚本（pipe/）负责流程编排和参数传递
   - 功能模块（script/、python/）实现具体算法
   - 禁止在单个文件内混写多种语言

2. **路径管理**
   - 禁止在模块中硬编码绝对路径
   - 所有路径通过配置文件或命令行参数传入
   - 使用相对路径时基于 `$SCRIPT_DIR` 或 `$PROJECT_ROOT`

3. **日志和断点续跑**
   - 使用基于日志的断点续跑机制
   - 按样本写入 `success.log` 和 `fail.log`
   - 禁止通过扫描输出目录判断完成状态

4. **Python 脚本**
   - 支持命令行参数解析（argparse）
   - 包含 `if __name__ == "__main__"` 入口
   - 可作为模块导入或独立运行

### 测试和调试

- 单样本测试：使用 `0_NUMTs_detection_2.0_单样本.sh`
- 查看日志：检查 `tmp/logs/` 目录
- 中间结果：保留在输出目录中，可用于调试

### 常见问题

1. **聚类结果为空**
   - 检查输入 BAM 是否包含 chrM 比对
   - 确认 discordant reads 提取是否成功

2. **断点检测失败**
   - 查看 PSL 文件是否生成
   - 检查候选区域是否有足够的 reads 支持

3. **配置文件错误**
   - 使用 `set -euo pipefail` 确保错误及时暴露
   - 检查必需变量是否定义（`: "${VAR:?missing VAR}"`）

## 数据流概览

```
BAM 文件
  ↓ [提取 discordant/split reads]
SAM 文件 (*.mt.disc.sam, *.mt.split.sam)
  ↓ [聚类脚本 0_1_找聚类.py]
候选区域 (*.cluster.tsv, *.breakpointINPUT.tsv)
  ↓ [提取序列 + bwa-mem 重比对]
PSL 文件 (*_all_regions.psl)
  ↓ [断点脚本 0_3_找断点.py]
断点结果 (*.AllBreakpoints.tsv, *.ConfidentBreakpoints.tsv)
  ↓ [统计分析]
频率表、可视化图表
```

## 参考文档

- 各模块的 README.md 包含详细的输出文件解读
- `1-1-获得NUMTs分布/README.md` - 检测流程和结果文件说明
- `1-5-Vardetection/README.md` - 变异检测结果解读
- `4-统计性描述/README.md` - 统计分析说明
