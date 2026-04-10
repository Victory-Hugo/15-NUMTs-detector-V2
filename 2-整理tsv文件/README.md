# 整理 TSV 文件流程

本目录包含 NUMTs 检测结果的整理、过滤和筛选流程。通过多次调用 `pipe/` 下的脚本，可以从原始结果产出不同阈值条件和样本子集的结果。

## 流程概览

```
原始 per-sample TSV (NUMTs Detector 输出)
        │
        ▼
  1-1-merge-tsv.sh          合并所有样本的 TSV 文件
        │
        ▼
  1-2-strict-tsv.sh         宽松阈值 → 严格阈值 (Cluster_No ≥ 5)
        │
        ▼
  3-exclude-reference-NUMT.sh  排除已知参考 NUMTs (±500bp)
        │
        ├──▶ 2-1-filter-merge-bed-confident.sh  提取 confident merge_bed
        │
        └──▶ 2-5-filter-new-sample.sh           筛选新样本子集
```

## 各步骤说明

| 步骤 | 脚本 | 配置文件 | 功能 |
|------|------|----------|------|
| 1-1 | `pipe/1-1-merge-tsv.sh` | `conf/1-1-merge.yaml` | 将各样本的 TSV 按类型合并为 7 个汇总文件 |
| 1-2 | `pipe/1-2-strict-tsv.sh` | `conf/1-2-strict.yaml` | 对宽松阈值结果施加严格过滤 (Cluster_No ≥ 5) |
| 3 | `pipe/3-exclude-reference-NUMT.sh` | `conf/3-exclude-reference.yaml` | 去除与已知参考 NUMTs 重叠的片段 (两侧扩展 500bp) |
| 2-1 | `pipe/2-1-filter-merge-bed-confident.sh` | `conf/2-1-filter-merge-bed-confident.yaml` | 从 merge_bed 中仅保留 confident 区域 |
| 2-5 | `pipe/2-5-filter-new-sample.sh` | `conf/2-5-filter-new-sample.yaml` | 从所有样本中按 ID 列表提取新测序样本 |

## 结果目录说明

SLE 队列的结果按以下 8 种条件组织，每个条件对应一个目录：

### 条件组合

| 编号 | 阈值 | 参考 NUMTs | 样本范围 | 目录 |
|------|------|------------|----------|------|
| 1 | 宽松 | 未排除 | 所有样本 | `1-原始结果/data/3-SLE/1-宽松阈值/1-所有样本` |
| 2 | 宽松 | 未排除 | 新样本 | `1-原始结果/data/3-SLE/1-宽松阈值/2-新样本` |
| 3 | 严格 | 未排除 | 所有样本 | `1-原始结果/data/3-SLE/2-严格阈值/1-所有样本` |
| 4 | 严格 | 未排除 | 新样本 | `1-原始结果/data/3-SLE/2-严格阈值/2-新样本` |
| 5 | 严格 | 已排除 | 所有样本 | `1-原始结果/data/3-SLE/3-严格阈值-排除参考NUMTs/1-所有样本` |
| 6 | 严格 | 已排除 | 新样本 | `1-原始结果/data/3-SLE/3-严格阈值-排除参考NUMTs/2-新样本` |
| 7 | 宽松 | 已排除 | 所有样本 | `1-原始结果/data/3-SLE/4-宽松阈值-排除参考NUMTs/1-所有样本` |
| 8 | 宽松 | 已排除 | 新样本 | `1-原始结果/data/3-SLE/4-宽松阈值-排除参考NUMTs/2-新样本` |

## 使用方式

### 运行单个步骤

```bash
# 以步骤 1-1 为例，合并所有样本 TSV
bash pipe/1-1-merge-tsv.sh --config conf/1-1-merge.yaml

# 步骤 1-2，宽松阈值转为严格阈值
bash pipe/1-2-strict-tsv.sh --config conf/1-2-strict.yaml

# 步骤 3，排除参考 NUMTs
bash pipe/3-exclude-reference-NUMT.sh --config conf/3-exclude-reference.yaml

# 步骤 2-1，提取 confident merge_bed
bash pipe/2-1-filter-merge-bed-confident.sh --config conf/2-1-filter-merge-bed-confident.yaml

# 步骤 2-5，筛选新样本
bash pipe/2-5-filter-new-sample.sh --config conf/2-5-filter-new-sample.yaml
```

### 产出 8 种结果条件

每种结果条件通过修改对应步骤配置文件中的 `input_dir` 和 `output_dir` 路径，对同一流程多次调用获得。以严格阈值-排除参考 NUMTs-所有样本为例：

```bash
# 1. 修改 conf/1-2-strict.yaml
#    input_dir  指向 宽松阈值/1-所有样本
#    output_dir 指向 严格阈值/1-所有样本
bash pipe/1-2-strict-tsv.sh --config conf/1-2-strict.yaml

# 2. 修改 conf/3-exclude-reference.yaml
#    input_dir  指向 严格阈值/1-所有样本
#    output_dir 指向 严格阈值-排除参考NUMTs/1-所有样本
bash pipe/3-exclude-reference-NUMT.sh --config conf/3-exclude-reference.yaml
```

其余条件通过类似方式，将 `input_dir` / `output_dir` 指向目标目录即可。每个脚本都支持 `force: true` 配置项跳过断点续跑检查，或使用 `overwrite: true` 覆盖已有输出。

### 新样本筛选 (步骤 2-5)

新样本筛选依赖 `meta/1-新测序数据ID.tsv` 中的样本 ID 列表。修改配置文件中的 `input_dir` 为任意"所有样本"目录，`output_dir` 为对应的"新样本"目录即可：

```bash
# 修改 conf/2-5-filter-new-sample.yaml 中的 input_dir 和 output_dir
bash pipe/2-5-filter-new-sample.sh --config conf/2-5-filter-new-sample.yaml
```

## 配置文件关键参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `input_dir` | 输入结果目录 | `.../1-宽松阈值/1-所有样本` |
| `output_dir` | 输出结果目录 | `.../2-严格阈值/1-所有样本` |
| `force` | 强制重新运行（忽略断点续跑日志） | `true` / `false` |
| `overwrite` | 覆盖已有输出文件 | `true` / `false` |
| `min_cluster_no` | 严格阈值的最小 Cluster_No | `5` |
| `extension` | 参考 NUMTs 两侧扩展碱基数 | `500` |
| `jobs` | 并行任务数 | `9` |
