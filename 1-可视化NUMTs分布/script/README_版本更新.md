# TSV文件合并工具 - 版本更新说明

## 概述

由于朋友修改了TSV文件的后缀，现在有两个版本的脚本来处理不同格式的文件：

- **旧版本脚本**: `0-打包.py` - 处理原始格式
- **新版本脚本**: `0-打包_LDY.py` - 处理新格式 (LDY版本)
- **智能Shell脚本**: `0-整理输出的所有tsv文件.sh` - 自动检测并使用合适的脚本

## 文件格式对比

### 旧版本文件格式
```
*.Breakpoints.tsv
*.mt.disc.sam.breakpointINPUT.tsv
*.mt.disc.sam.cluster.summary.tsv
*.mt.disc.sam.cluster.tsv
```

### 新版本文件格式 (LDY)
```
*ConfidentBreakpoints.tsv
*AllBreakpoints.tsv
*.mt.disc.sam.breakpointINPUT.tsv
*.mt.disc.sam.cluster.summary.tsv
*.mt.disc.sam.cluster.tsv
```

## 使用方法

### 方法1: 自动版本检测 (推荐)

直接执行智能Shell脚本，它会自动检测文件版本并选择合适的Python脚本：

```bash
bash 0-整理输出的所有tsv文件.sh
```

脚本会自动:
- ✓ 检测输入目录中的文件格式
- ✓ 选择对应的处理脚本
- ✓ 如果混合版本，提示用户选择
- ✓ 执行文件合并和压缩
- ✓ 验证压缩文件完整性
- ✓ 准备circos输入文件

### 方法2: 手动指定旧版本脚本

如果你确定是旧版本文件，可以直接使用旧脚本：

```bash
python 0-打包.py <输入目录> <输出目录> <压缩文件路径>

# 示例
python 0-打包.py ./tsv_files ./output ./backup.tar.gz
```

### 方法3: 手动指定新版本脚本 (LDY)

如果你确定是新版本文件，可以直接使用新脚本：

```bash
python 0-打包_LDY.py <输入目录> <输出目录> <压缩文件路径>

# 示例
python 0-打包_LDY.py ./tsv_files ./output ./backup.tar.gz
```

## 脚本配置

编辑Shell脚本中的以下变量以适应你的环境：

```bash
PYTHON3="/home/luolintao/miniconda3/envs/pyg/bin/python3"  # Python路径
PYTHON_SCRIPT="...0-打包.py"                                # 旧版本脚本
PYTHON_SCRIPT_LDY="...0-打包_LDY.py"                        # 新版本脚本
PREPARE_CIRCOS_SCRIPT="...1-准备NUMTs可视化.py"            # Circos准备脚本
INPUT_CSV_DIR="/path/to/input"                              # 输入目录
OUTPUT_DIR="/path/to/output"                                # 输出目录
TAR_GZ_PATH="/path/to/backup.tar.gz"                        # 压缩文件路径
```

## 版本检测逻辑

Shell脚本使用以下逻辑来检测文件版本：

| 检测到的文件 | 版本 | 使用脚本 |
|-----------|------|--------|
| `*.Breakpoints.tsv` 存在 | 旧版本 (0) | `0-打包.py` |
| `*ConfidentBreakpoints.tsv` 或 `*AllBreakpoints.tsv` 存在 | 新版本 (1) | `0-打包_LDY.py` |
| 两种格式都存在 | 混合版本 (2) | 用户选择 |

## 处理流程

### 旧版本处理流程
```
输入目录
    ├── sample1.Breakpoints.tsv
    ├── sample1.mt.disc.sam.breakpointINPUT.tsv
    ├── sample1.mt.disc.sam.cluster.summary.tsv
    ├── sample1.mt.disc.sam.cluster.tsv
    └── ...
         ↓
    0-打包.py
         ↓
输出目录
    ├── all_individuals_Breakpoints.tsv
    ├── all_individuals_mt.disc.sam.breakpointINPUT.tsv
    ├── all_individuals_mt.disc.sam.cluster.summary.tsv
    └── all_individuals_mt.disc.sam.cluster.tsv
```

### 新版本处理流程
```
输入目录
    ├── sample1.ConfidentBreakpoints.tsv
    ├── sample1.AllBreakpoints.tsv
    ├── sample1.mt.disc.sam.breakpointINPUT.tsv
    ├── sample1.mt.disc.sam.cluster.tsv
    └── ...
         ↓
    0-打包_LDY.py
         ↓
输出目录
    ├── all_individuals_ConfidentBreakpoints.tsv
    ├── all_individuals_AllBreakpoints.tsv
    ├── all_individuals_mt.disc.sam.breakpointINPUT.tsv
    ├── all_individuals_mt.disc.sam.cluster.summary.tsv (如果存在)
    └── all_individuals_mt.disc.sam.cluster.tsv
```

## 新版本脚本的主要特性

✓ **自动版本检测**: 无需手动指定，脚本自动判断文件格式
✓ **支持新的后缀格式**: `*ConfidentBreakpoints.tsv` 和 `*AllBreakpoints.tsv`
✓ **保持兼容性**: 保留对其他文件格式的支持
✓ **智能文件合并**: 自动识别并合并相同类型的文件
✓ **完整性验证**: 验证压缩文件的完整性
✓ **详细日志**: 彩色输出，清晰显示处理进度
✓ **安全删除**: 只在验证通过后才删除原始文件
✓ **空文件处理**: 智能跳过空文件或仅含标题的文件

## 故障排查

### 问题1: 无法找到Python脚本
**解决方案**: 检查Shell脚本中的 `PYTHON_SCRIPT` 和 `PYTHON_SCRIPT_LDY` 路径是否正确

### 问题2: 混合版本无法自动处理
**解决方案**: 按照提示选择要使用的脚本版本（1或2），或者分别在不同的目录中处理

### 问题3: 压缩文件验证失败
**解决方案**: 检查磁盘空间是否充足，或手动检查输出目录中的文件

### 问题4: 找不到circos准备脚本
**解决方案**: 检查 `PREPARE_CIRCOS_SCRIPT` 路径，确保脚本存在

## 文件清单

脚本目录中应该包含以下文件：

```
/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布/script/
├── 0-打包.py                          # 旧版本合并脚本
├── 0-打包_LDY.py                      # 新版本合并脚本 (保留旧代码)
├── 0-整理输出的所有tsv文件.sh         # 智能检测脚本 (已更新)
├── 0-整理输出的所有tsv文件_智能版.sh  # 智能版本备份
├── 1-准备NUMTs可视化.py               # Circos准备脚本
├── 2-circos.sh                        # Circos绘图脚本
└── README_版本更新.md                 # 本文件
```

## 示例使用场景

### 场景1: 处理旧版本文件
```bash
# 自动检测并处理
bash 0-整理输出的所有tsv文件.sh
# 或手动指定
python 0-打包.py ./old_tsv_dir ./merged_output ./backup.tar.gz
```

### 场景2: 处理新版本文件 (LDY)
```bash
# 自动检测并处理
bash 0-整理输出的所有tsv文件.sh
# 或手动指定
python 0-打包_LDY.py ./new_tsv_dir ./merged_output ./backup.tar.gz
```

### 场景3: 处理混合版本文件
```bash
# 使用智能脚本，会提示选择
bash 0-整理输出的所有tsv文件.sh
# 选择1: 使用旧版本脚本处理旧格式
# 选择2: 使用新版本脚本处理新格式
```

## 技术细节

### 旧版本脚本 (0-打包.py)
- 识别文件模式: `*.Breakpoints.tsv`, `*.mt.disc.sam.*`
- 处理4种文件类型
- 适用于原始数据格式

### 新版本脚本 (0-打包_LDY.py)
- 识别文件模式: `*ConfidentBreakpoints.tsv`, `*AllBreakpoints.tsv`, `*.mt.disc.sam.*`
- 处理5种文件类型 (包括新的断点类型)
- 适用于修改后缀的数据格式
- 代码结构完全保留，便于维护

### 智能Shell脚本
- 使用 `find` 命令扫描输入目录
- 基于文件数量判断版本
- 支持交互式选择
- 自动调用合适的Python脚本

## 更新历史

- **v1.0** (原始版本): 支持旧格式文件的合并
- **v2.0** (LDY新版本): 
  - 新增 `0-打包_LDY.py` 脚本处理新格式
  - 更新Shell脚本支持自动版本检测
  - 保留原始代码以确保向后兼容

## 联系与支持

如有问题或建议，请检查:
1. 脚本配置变量是否正确
2. 输入目录路径是否存在
3. Python环境是否正确安装依赖包
4. 磁盘空间是否充足

---

**最后更新**: 2025-10-22
**脚本版本**: v2.0
