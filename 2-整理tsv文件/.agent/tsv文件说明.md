# TSV 文件说明与下游使用建议

本文档用于解释目录下各个 TSV/GZ 文件的来源、含义、上下游关系，以及后续分析时应该优先使用哪些文件。

## 1. 先给结论

原始文件目录下的文件不是同一层级的结果，主要分成两条链路：

1. `1-1-获得NUMTs分布` 产生单样本 NUMT 候选区域和断点结果。
2. `2-整理tsv文件` 把这些单样本 TSV 按文件类型跨样本合并，得到 `1` 到 `7` 号文件。
3. `1-5-Vardetection` 另外基于每个样本的 `*_all_regions.fasta / *.psl` 生成每个样本的 `*_numts.bed`，`8-merge_bed.tsv.gz` 是这类 BED 风格表的跨样本汇总产物，不属于 `2-整理tsv文件` 那条合并脚本直接输出。

如果后续是做 NUMT 位点/断点层面的主分析，通常应优先使用：

- `7-mt_disc_cluster.tsv`：候选 NUMT 区域主表
- `3-confident_breakpoints.tsv`：高可信断点主表

如果需要保留全量候选证据，再辅助使用：

- `1-all_breakpoints.tsv`

通常不建议把以下文件作为主分析输入：

- `2-breakpoints_old.tsv`
- `5-mt_disc_cluster_old.tsv`
- `4-mt_disc_breakpoint_input.tsv`
- `6-mt_disc_cluster_summary.tsv`

`8-merge_bed.tsv.gz` 也不适合作为 NUMT 事件主结果表，它更适合服务于 `1-5-Vardetection` 阶段的线粒体区间注释或变异分析。

## 2. 文件来源链路

### 2.1 上游脚本与目录

本次确认的核心目录如下：

- `/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-1-获得NUMTs分布`
- `/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-5-Vardetection`
- `/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/2-整理tsv文件`

### 2.2 `1-1-获得NUMTs分布` 的逻辑

该步骤从样本 BAM/CRAM 中提取 mtDNA 相关 discordant reads，再聚类得到候选 NUMT 区域，并对候选区域做断点检测。

其核心输出分两类：

1. 候选区域级：
- `*.mt.disc.sam.cluster.tsv`
- `*.mt.disc.sam.cluster.summary.tsv`
- `*.mt.disc.sam.breakpointINPUT.tsv`
- `*.mt.disc.sam.cluster.old.tsv`

2. 断点级：
- `*.AllBreakpoints.tsv`
- `*.ConfidentBreakpoints.tsv`
- `*.Breakpoints.old.tsv`

### 2.3 `2-整理tsv文件` 的逻辑

`2-整理tsv文件/conf/merge_config.yaml` 明确指定了 7 类待合并文件：

- `*.AllBreakpoints.tsv` -> `1-all_breakpoints.tsv`
- `*.Breakpoints.old.tsv` -> `2-breakpoints_old.tsv`
- `*.ConfidentBreakpoints.tsv` -> `3-confident_breakpoints.tsv`
- `*.mt.disc.sam.breakpointINPUT.tsv` -> `4-mt_disc_breakpoint_input.tsv`
- `*.mt.disc.sam.cluster.old.tsv` -> `5-mt_disc_cluster_old.tsv`
- `*.mt.disc.sam.cluster.summary.tsv` -> `6-mt_disc_cluster_summary.tsv`
- `*.mt.disc.sam.cluster.tsv` -> `7-mt_disc_cluster.tsv`

因此，原始文件目录里的 `1` 到 `7` 号文件，本质上都是 `1-1` 单样本结果的跨样本拼接版。

其中：

- `1-all_breakpoints.tsv`
- `3-confident_breakpoints.tsv`

在合并时额外添加了来源列：

- `source_file_basename`
- `source_sample_id`
- `source_region_chr`
- `source_region_start`
- `source_region_end`
- `source_region_key`

这说明它们在跨样本合并后，仍然保留了“来自哪个样本、哪个候选区域”的追溯信息。

### 2.4 `1-5-Vardetection` 的逻辑

该步骤对每个样本的候选 NUMT 序列做 mtDNA 参考比对，并生成：

- `<sample>_numts.bed`
- Human-only 变异表
- Human-Chimp 比较表

其中 `<sample>_numts.bed` 由 `python/psl_to_numts_bed.py` 从样本级 PSL 文件中提取“每条 query 在 mtDNA 上的最佳命中区间”生成，格式为 4 列：

- query_name
- numt_start
- numt_end
- 标签 `best_hit_to_mtDNA`

`8-merge_bed.tsv.gz` 的内容格式与 `<sample>_numts.bed` 完全一致，因此可以判断它是多个样本 `*_numts.bed` 的合并结果或等价汇总结果，而不是 `2-整理tsv文件` 直接产物。

## 3. 原始文件 目录下每个文件是什么意思

目录内容如下：

- `1-all_breakpoints.tsv`
- `2-breakpoints_old.tsv`
- `3-confident_breakpoints.tsv`
- `4-mt_disc_breakpoint_input.tsv`
- `5-mt_disc_cluster_old.tsv`
- `6-mt_disc_cluster_summary.tsv`
- `7-mt_disc_cluster.tsv`
- `8-merge_bed.tsv.gz`

### 3.1 `1-all_breakpoints.tsv`

含义：
每个候选区域中“所有统计到的断点”总表，为全量断点结果。

来源：
`1-1-获得NUMTs分布/script/0_3_找断点.py` 产生单区域 `*.AllBreakpoints.tsv`，再由 `2-整理tsv文件` 跨样本合并。

特点：

- 有表头
- 每行对应一个断点统计条目
- `type=all`
- 包含核侧和线粒体侧断点
- 合并后带有来源区域信息

用途：

- 适合做探索性分析
- 适合保留“所有候选断点”证据

限制：

- 包含低可信或单边证据
- 不应默认当作最终高可信断点结果

### 3.2 `2-breakpoints_old.tsv`

含义：
旧版兼容断点结果，实质上是旧格式的 confident breakpoint 输出。

来源：
`0_3_找断点.py` 输出的 `*.Breakpoints.old.tsv` 合并版。

特点：

- 无表头
- 列结构老旧
- 同时保留 `Tstart/Tend`
- 不如新版表直观

用途：

- 仅用于兼容旧脚本或回溯旧分析

限制：

- 不建议作为新分析主输入

### 3.3 `3-confident_breakpoints.tsv`

含义：
高可信断点主表。只保留核侧和线粒体侧可由同一批 reads 互相对应上的断点。

来源：
`0_3_找断点.py` 输出的 `*.ConfidentBreakpoints.tsv` 合并版。

特点：

- 有表头
- `type=confident`
- 合并后保留来源样本和候选区域信息
- 是断点层面的主结果文件

用途：

- 断点位置统计
- 高可信 NUMT 事件支持证据筛选
- 跨样本复现分析

建议：

- 若后续需要做“断点级主分析”，优先使用这个文件

### 3.4 `4-mt_disc_breakpoint_input.tsv`

含义：
断点检测阶段的输入区域清单，不是最终结果。

来源：
`1-1` 聚类脚本 `0_1_找聚类.py` 生成的 `*.mt.disc.sam.breakpointINPUT.tsv` 合并版。

内容：

- 样本 ID
- cluster 支持数
- `disc.sam` 路径
- `split.sam` 路径
- 原始 BAM/CRAM 路径
- 候选区域坐标

用途：

- 重新跑断点分析
- 追踪上游输入来源

限制：

- 这是流程输入表，不代表已验证的 NUMT 事件
- 不应用作最终统计结果

### 3.5 `5-mt_disc_cluster_old.tsv`

含义：
旧版兼容的 cluster 明细表，更接近 read 级记录。

来源：
`0_1_找聚类.py` 生成的 `*.mt.disc.sam.cluster.old.tsv` 合并版。

特点：

- 有表头
- 一行更接近一条 read 证据
- 包含大量 SAM 字段
- 数据量通常很大

用途：

- 追踪某个 cluster 由哪些 reads 支持
- 检查原始证据

限制：

- 不适合直接做样本间主统计
- 不建议作为主结果输入

### 3.6 `6-mt_disc_cluster_summary.tsv`

含义：
cluster 汇总简表。

来源：
`0_1_找聚类.py` 生成的 `*.mt.disc.sam.cluster.summary.tsv` 合并版。

特点：

- 无正式表头
- 第一列常是 pandas 自动索引
- 本质是 cluster 计数摘要，不是完整事件表

用途：

- 快速 QC
- 粗略查看每个 cluster 的支持数概况

限制：

- 信息不完整
- 不适合作为正式下游分析输入

### 3.7 `7-mt_disc_cluster.tsv`

含义：
候选 NUMT 区域主表，是 cluster 层面的核心结果。

来源：
`0_1_找聚类.py` 生成单样本 `*.mt.disc.sam.cluster.tsv`，再由 `2-整理tsv文件` 合并。

字段核心含义：

- `chr/start/end`：候选 NUMT 区域
- `Cluster_No`：核侧支持位置数
- `reads`：核侧位置列表
- `mt_positions`：对应的 mtDNA 位置列表
- `subCluster_No`：线粒体侧位置数
- `Cluster_ID`：候选事件唯一标识
- `IndividualID`：样本 ID

用途：

- NUMT 候选事件统计
- 跨样本共享位点归并
- 候选区域可视化
- 作为后续断点结果的事件锚点

建议：

- 若后续需要做“NUMT 位点/事件级主分析”，优先使用这个文件

### 3.8 `8-merge_bed.tsv.gz`

含义：
样本级 `*_numts.bed` 的汇总结果，记录每条 query 在 mtDNA 上的最佳命中区间。

来源：
推断来自 `1-5-Vardetection` 中每个样本的 `<sample>_numts.bed` 合并后再压缩得到。

理由：

- 文件内容与 `psl_to_numts_bed.py` 的输出格式完全一致
- 每行形如 `sample_region|readID  start  end  best_hit_to_mtDNA`
- `2-整理tsv文件` 配置中并未定义这个文件的合并规则

用途：

- 作为变异检测阶段的 mtDNA 命中区间参考
- 追踪每条序列在 mtDNA 上的最佳比对范围

限制：

- 这是 read/query 命中区间表，不是 NUMT 事件主表
- 不适合直接代替 `cluster.tsv` 或 `confident_breakpoints.tsv`

## 4. 哪些文件属于上游，哪些属于最终结果

### 4.1 上游输入/中间文件

以下文件更偏向流程输入、兼容输出或中间摘要：

- `4-mt_disc_breakpoint_input.tsv`
- `5-mt_disc_cluster_old.tsv`
- `6-mt_disc_cluster_summary.tsv`
- `2-breakpoints_old.tsv`
- `8-merge_bed.tsv.gz`

其中：

- `4` 是后续断点检测输入
- `5` 和 `2` 是 legacy 兼容输出
- `6` 是简表/QC 表
- `8` 是变异注释链路使用的区间表

### 4.2 最终主结果文件

真正适合作为主结果进行下游分析的文件是：

- `7-mt_disc_cluster.tsv`
- `3-confident_breakpoints.tsv`

可作为补充结果使用：

- `1-all_breakpoints.tsv`

理解方式：

- `7` 负责回答“有哪些 NUMT 候选事件”
- `3` 负责回答“这些事件的高可信断点在哪里”
- `1` 负责保留“所有候选断点，不论是否高可信”

## 5. 后续分析应该用哪些文件

### 5.1 如果目标是统计 NUMT 位点/事件

应优先使用：

- `7-mt_disc_cluster.tsv`

原因：

- 一行对应一个候选 NUMT 区域
- 粒度清晰
- 可直接按 `Cluster_ID` 或 `chr/start/end` 做事件级分析

典型用途：

- 每个样本 NUMT 事件数
- NUMT 染色体分布
- 跨样本共享/特异 NUMT 位点

### 5.2 如果目标是统计断点位置

应优先使用：

- `3-confident_breakpoints.tsv`

必要时补充：

- `1-all_breakpoints.tsv`

原因：

- `3` 是高可信断点集合，更适合正式统计和作图
- `1` 可用于探索未进入 confident 集的边缘证据

典型用途：

- 核侧左/右断点分布
- mtDNA 断点热点
- 每个候选事件的断点支持强度

### 5.3 如果目标是做变异分析或 mtDNA 命中区间辅助注释

可以使用：

- `8-merge_bed.tsv.gz`

但应明确：

- 它不是 NUMT 事件表
- 更适合给 `1-5-Vardetection` 相关分析做辅助，而不是替代 `7` 或 `3`

## 6. 不建议直接作为主分析输入的文件

以下文件通常不建议直接用于正式下游主分析：

- `2-breakpoints_old.tsv`
- `4-mt_disc_breakpoint_input.tsv`
- `5-mt_disc_cluster_old.tsv`
- `6-mt_disc_cluster_summary.tsv`

原因分别是：

- `2`：旧版格式，不如新版清晰
- `4`：流程输入，不是结果
- `5`：read 级明细，过细且冗余
- `6`：摘要过简，信息不足

## 7. 推荐的下游分析组合

最推荐的组合是：

1. 用 `7-mt_disc_cluster.tsv` 定义 NUMT 候选事件集合。
2. 用 `3-confident_breakpoints.tsv` 为这些事件补充高可信断点信息。
3. 如需保留全量证据，再额外参考 `1-all_breakpoints.tsv`。
4. 若进入变异检测链路，再根据需要使用 `8-merge_bed.tsv.gz`。

## 8. 一句话版判断

如果只能记一句话：

- `7-mt_disc_cluster.tsv` 是 NUMT 事件主表。
- `3-confident_breakpoints.tsv` 是高可信断点主表。
- `1-all_breakpoints.tsv` 是全量断点补充表。
- 其余文件大多属于上游输入、兼容格式、QC 摘要或变异链路辅助文件。
