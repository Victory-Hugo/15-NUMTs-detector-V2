# NUMTs 核断点富集分析 pipeline

本流程用于复现论文 Figure 7D 思路：以 confident NUMT nuclear breakpoints 上下游各 `100 bp` 的 flank 作为分析区间，在多个核基因组功能区域中做 permutation 富集分析。频率分层保留本项目当前分类：

```text
all
common
low-frequency
rare
ultra-rare
```

## 输入

主流程需要四类输入：

```text
/mnt/l/20-NUMTs/6-NUMTs频率分布描述/2-严格阈值-v2/output/01-qc/2-confident_breakpoints.pass.tsv.gz
/mnt/l/20-NUMTs/6-NUMTs频率分布描述/2-严格阈值-v2/output/01-qc/4-mt_disc_breakpoint_input.pass.tsv.gz
/mnt/l/20-NUMTs/6-NUMTs频率分布描述/2-严格阈值-v2/output/02-tables/4-mt_disc_breakpoint_input.min-support-1.allCluster.tsv
/mnt/l/20-NUMTs/6-NUMTs频率分布描述/2-严格阈值-v2/output/02-tables/4-numt-frequency-by-cluster.tsv
```

其中 `4-numt-frequency-by-cluster.tsv` 只用于给 breakpoint flank 赋予 `common/low-frequency/rare/ultra-rare` 分层；核基因组位置主输入来自 confident nuclear breakpoint 表。

目标区域默认读取：

```text
data/regions/*.bed
```

## BED 资源准备

`data/regions/` 是长期可复用资源，不接入主流程。当前目录同时包含论文 Figure 7D 的 22 个区域 BED 和本项目原有的 4 个自定义区域 BED。需要重新构建文章区域 BED 时，单独运行：

```bash
bash script/prepare_article_regions_bed.sh \
  --output-dir data/regions \
  --source-dir data/source/article_regions \
  --tmp-dir data/tmp/article_regions \
  --raw-dir data/raw \
  --bedtools-bin bedtools
```

当前应包含 26 个 BED：本项目原有的 `01-Promoter.bed` 到 `04-STR_without_segmental_duplications.bed`，以及 `.agent/补齐bed.md` 中列出的 22 个文章区域 BED。

## 运行

```bash
bash pipe/run_enrichment_analysis.sh
```

可指定配置：

```bash
bash pipe/run_enrichment_analysis.sh conf/Config.yaml
```

并行由 shell 后台任务队列调度；Python 模块只处理单一任务，不在 Python 内部做并行。

## 输出

```text
output/1-prepared-breakpoint-flanks/  # 分层 nuclear breakpoint flank BED 与任务清单
output/2-enrichment/                  # 每个 BED x 频率层级的单任务富集结果
output/3-summary/                     # 合并后的 enrichment_summary.tsv
output/4-figures/                     # 热图与观测比例条形图
output/logs/                          # pipeline 和模块日志
report/                               # Markdown 结果报告
```

最终报告：

```text
report/enrichment_report.md
```
