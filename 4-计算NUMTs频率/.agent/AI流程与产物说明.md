# NUMTs 频率流程简明说明

本文档用于让后续协作者快速理解这个项目的代码结构、数据流、每个脚本的职责、每个主要产物文件的意义，以及**哪些地方最容易误读**。

目标不是介绍背景知识，而是避免修改时犯错。

---

## 1. 先记住的几条原则

### 1.1 主入口

正式流程入口是：

```bash
bash pipe/01-run_numts_frequency.sh --config conf/Config.yaml
```

不要调用 `script/3-统计频率.py` 来生成正式结果。那个是旧脚本，保留作历史参考，不是当前主流程。

### 1.2 当前主分析口径

当前流程的主分析已经按较新的修正规则运行：

1. `length_bed` 的长度不是外包 span，而是**区间并集长度**。
2. `complex_mt_interval=True` 的事件会被标记，且**默认不进入主 distinct 统计**。
3. distinct NUMT 的跨样本合并是基于**真实 breakpoint 区间**，不是旧版中点聚类。
4. `carrier_fraction` 表示“至少携带 1 个该类 NUMT 的样本比例”，不是 carrier 总和除以样本数。
5. mt breakpoint 累积图默认看的是 **cumulative fraction**，不是单纯 cumulative count。

### 1.3 不要混淆的几个字段

- `mt_source_span_bp`
  当前语义是 **union span**，不是旧版 envelope span。
- `mt_source_envelope_span_bp`
  只是辅助 QC 字段，表示最外层包络长度。
- `mt_length_for_main_bp`
  这是主分析使用的样本事件长度。优先用 mt breakpoint span，其次在非复杂事件里回退到 union span。
- `eligible_for_main_analysis`
  这是是否进入主 distinct 统计的硬条件，不是参考字段。
- `distinct_nuclear_start/end`
  当前应理解为**merged nuclear breakpoint 对**，不是简单簇范围。
- `nuclear_cluster_range_start/end`
  这是簇范围信息，和 `distinct_nuclear_start/end` 不同。

---

## 2. 目录职责

### 2.1 `pipe/`

- `01-run_numts_frequency.sh`
  总控脚本。负责：
  - 读 `conf/Config.yaml`
  - 激活 conda 环境
  - 组织 1 到 9 步
  - 管理 `step 1` 的缓存复用
  - 把 `length_bed` 的全量汇总缓存写到 `output/cache/`

### 2.2 `python/`

这里的 `01` 到 `09` 是流程步骤脚本，负责调度对应模块。

### 2.3 `python/numts_frequency/`

这里是业务逻辑模块，真正的算法和表结构定义在这里。

### 2.4 `script/`

旧脚本和辅助脚本目录。**不要把这里的历史脚本当正式流程输出来源。**

### 2.5 `output/`

当前所有正式产物都写在这里。后续 AI 应优先读取这里的表，而不是去拼装中间逻辑。

---

## 3. 每个代码文件的作用

## 3.1 `pipe/01-run_numts_frequency.sh`

作用：

1. 读配置。
2. 激活 conda。
3. 跑 9 个步骤。
4. 如果 `output/00-qc/` 中关键 QC 文件都已存在且非空，则跳过 `step 1`。
5. 复用或生成 `output/cache/length_bed.full_region_summary.tsv`。

注意：

- 这里的 `step 1 cache` 只是流程缓存，不代表科学逻辑变化。
- 如果你删掉 `output/00-qc/`，流程会重跑 `step 1`。
- 如果你删掉 `output/cache/length_bed.full_region_summary.tsv`，流程会重新扫描整个 `length_bed.tsv.gz`。

## 3.2 `python/01_validate_inputs.py`

作用：

1. 读取原始 `breakpoints`、`clusters`、`meta`。
2. 生成样本纳入/排除结果。
3. 构建或复用 `length_bed` 全量缓存。
4. 从全量缓存中筛出当前保留样本对应的 `_length_bed_region_summary.tsv`。
5. 生成 QC 报告和 traceability 表。

这是整个流程里最慢的一步之一。

## 3.3 `python/02_build_sample_event_table.py`

作用：

把 cluster 表和 breakpoint 证据整合为 `sample_event` 层表。

这是样本内事件层，不是 distinct 层。

## 3.4 `python/03_build_distinct_catalog.py`

作用：

把 `sample_event` 跨样本合并成 distinct NUMT catalog，并输出 sample-event 到 distinct 的映射表。

这个阶段会过滤 `eligible_for_main_analysis=False` 的事件。

## 3.5 `python/04_annotate_catalog.py`

作用：

1. 重新补全/统一 sample_event 的长度与主分析字段。
2. 给 distinct catalog 增加长度摘要和 mt gene 注释。

## 3.6 `python/05_compute_main_results.py`

作用：

生成主结果表：

- 每样本 distinct NUMT 数
- 主分析长度分布
- 频率分类摘要
- carrier ratio
- frequency spectrum
- 主分析排除摘要

## 3.7 `python/06_compute_mt_breakpoint_results.py`

作用：

生成 mt breakpoint 相关表，包括：

- gene-level breakpoint counts
- 按频率类别分层的 gene counts
- breakpoint cumulative 表

## 3.8 `python/07_compute_stratified_results.py`

作用：

按 `Region`、`Population`、`Sex` 重算 distinct 频率和长度分布，并输出分层结果表。

## 3.9 `python/08_render_plots.py`

作用：

读取主结果表和分层结果表，统一绘图。

注意：

- mt breakpoint 图默认用的是 `overall_cumulative_fraction` / `cumulative_fraction`
- 不是旧版 `cumulative_count`

## 3.10 `python/09_validate_outputs.py`

作用：

对结果表做一致性检查，输出：

- `output/00-qc/06-output_validation_report.md`

这是结果级校验，不是单元测试。

## 3.11 `python/90_regression_cases.py`

作用：

最小回归测试脚本。主要验证：

- 不连续 mt 区间的 union span
- complex 标记
- distinct 合并逻辑
- carrier_fraction 逻辑
- cumulative_fraction 逻辑

修改核心逻辑后，应该先跑这个脚本。

## 3.12 `python/numts_frequency/io_utils.py`

作用：

- 读写 TSV/CSV
- 优先使用 `pyarrow` 快速读表
- 汇总 `length_bed`
- 构建 `length_bed` 全量缓存
- 从全量缓存切片当前 `region_key`

这是长度口径的上游真源之一。

## 3.13 `python/numts_frequency/catalog_utils.py`

作用：

- 建 breakpoint 检索索引
- 构建 sample_event 表
- 计算主分析纳入标记
- 基于 breakpoint 区间合并 distinct NUMT

这里定义了最核心的 distinct 合并逻辑。

## 3.14 `python/numts_frequency/annotation_utils.py`

作用：

- sample_event 与 length_summary 的对齐
- distinct 长度摘要
- mt gene 注释

## 3.15 `python/numts_frequency/stats_utils.py`

作用：

- 主统计
- exclusion summary
- mt breakpoint cumulative 表
- 分层 subgroup catalog 的重算

`carrier_fraction`、`cumulative_fraction` 的正式定义都在这里。

## 3.16 `python/numts_frequency/qc_utils.py`

作用：

- step 1 的输入摘要
- 样本纳入/排除
- traceability
- QC markdown 报告

## 3.17 `python/numts_frequency/plot_utils.py`

作用：

统一的绘图包装。

---

## 4. 主要产物文件的意义

下面只列正式流程中最重要的文件。

## 4.1 `output/cache/length_bed.full_region_summary.tsv`

意义：

- 对整个 `length_bed.tsv.gz` 做的一次性全量 region 汇总缓存。
- 后续 `step 1` 直接读取它，再按当前保留样本过滤。

不要误解：

- 这不是当前分析样本的最终长度表。
- 这是**全量缓存**。

## 4.2 `output/00-qc/_length_bed_region_summary.tsv`

意义：

- 从全量缓存切出来的、只属于当前 retained samples 的 region 长度汇总表。

重要字段：

- `mt_source_span_bp`：union span
- `mt_source_envelope_span_bp`：外包 span
- `complex_mt_interval`
- `mt_primary_fragment_span_bp`

## 4.3 `output/00-qc/02-sample_inclusion.tsv`

意义：

- 进入主流程的样本列表。

这是 retained sample 的权威来源。

## 4.4 `output/00-qc/03-sample_exclusion.tsv`

意义：

- 未进入 retained set 的样本及其排除原因。

## 4.5 `output/00-qc/04-event_traceability.tsv`

意义：

- 每个 retained cluster 是否能追溯到 confident nuclear breakpoint 的证据。

## 4.6 `output/01-sample-events/sample_numt_events.tsv`

意义：

- 样本内事件表。

这是 sample-event 层的核心表，不是 distinct 表。

重要字段：

- `nuclear_left_breakpoint_pos`
- `nuclear_right_breakpoint_pos`
- `mt_breakpoint_start`
- `mt_breakpoint_end`
- `mt_source_span_bp`
- `mt_breakpoint_span_bp`
- `mt_length_for_main_bp`
- `complex_mt_interval`
- `eligible_for_main_analysis`
- `main_analysis_exclusion_reason`

最重要的判断：

- 如果某个 event `eligible_for_main_analysis=False`，它不会进入主 distinct catalog。

## 4.7 `output/02-catalog/distinct_numt_catalog.tsv`

意义：

- 主分析 distinct NUMT catalog，已经是跨样本合并后的位点级表。

重要字段：

- `merged_nuclear_breakpoint_start`
- `merged_nuclear_breakpoint_end`
- `distinct_nuclear_start`
- `distinct_nuclear_end`
- `nuclear_cluster_range_start`
- `nuclear_cluster_range_end`
- `representative_mt_breakpoint`
- `carrier_count`
- `carrier_frequency`
- `frequency_class`

避免误读：

- `distinct_nuclear_start/end` 当前应理解为 merged nuclear breakpoint 对。
- `nuclear_cluster_range_start/end` 才是簇覆盖范围。

## 4.8 `output/02-catalog/distinct_numt_sample_map.tsv`

意义：

- sample_event 到 distinct_numt 的映射表。

这是样本计数、carrier 统计和 subgroup 重算的桥表。

## 4.9 `output/03-annotation/sample_numt_events.annotated.tsv`

意义：

- sample_event 层的补全版本，供后续统计与输出验证使用。

通常优先读这个，而不是未注释版。

## 4.10 `output/03-annotation/distinct_numt_catalog.annotated.tsv`

意义：

- distinct catalog 的补全版本，供主统计和 mt breakpoint 分析使用。

通常主统计读取的是这个文件。

## 4.11 `output/04-main-results/tables/sample_numt_counts.tsv`

意义：

- 每个 retained 样本拥有多少个主分析 distinct NUMT。

## 4.12 `output/04-main-results/tables/distinct_length_distribution.tsv`

意义：

- 主分析长度分布表。

注意：

- `length_level=distinct` 用的是 `median_mt_length_for_main_bp`
- `length_level=sample_event` 用的是 `mt_length_for_main_bp`

## 4.13 `output/04-main-results/tables/frequency_class_summary.tsv`

意义：

- 各频率类别的摘要。

重要字段：

- `carrier_sample_n`
- `carrier_fraction`
- `carrier_count_sum`
- `mean_numt_count_per_sample`

避免误读：

- `carrier_fraction` 现在是真正的样本携带比例。
- `carrier_count_sum` 才是所有位点 carrier 数的总和。

## 4.14 `output/04-main-results/tables/carrier_ratio_by_frequency_class.tsv`

意义：

- 供绘图和兼容旧接口的 carrier ratio 表。

语义和 `frequency_class_summary.tsv` 中的 `carrier_fraction` 一致。

## 4.15 `output/04-main-results/tables/main_analysis_exclusion_summary.tsv`

意义：

- 被排除出主分析的事件按 reason 聚合后的汇总表。

注意：

- 这是汇总，不是复杂事件逐条清单。

## 4.16 `output/05-mt-breakpoint/tables/mt_breakpoint_cumulative.tsv`

意义：

- mt breakpoint 的累计分布表。

重要字段：

- `overall_cumulative_fraction`
- `cumulative_fraction`

避免误读：

- 画图和解释时优先看 fraction，不要只看 count。

## 4.17 `output/06-stratified/...`

意义：

- 三个维度的 subgroup 重算结果：
  - `region`
  - `population`
  - `sex`

每个子目录下最重要的文件：

- `summary.tsv`
- `subgroup_metrics.tsv`
- `subgroup_length_distribution.tsv`
- `subgroup_frequency_summary.tsv`
- `subgroup_mt_breakpoint_cumulative.tsv`
- `subgroup_nuclear_chr_distribution.tsv`

## 4.18 `output/00-qc/06-output_validation_report.md`

意义：

- 最终输出是否满足基本逻辑约束的报告。

这是快速判断“这次跑出来的结果表是否明显坏掉”的第一入口。

---

## 5. 未来最容易犯的错

### 5.1 把旧脚本当正式流程

不要根据 `script/3-统计频率.py` 推断当前主流程逻辑。

### 5.2 把 `mt_source_span_bp` 误认为外包长度

当前不是。

当前：

- `mt_source_span_bp` = union span
- `mt_source_envelope_span_bp` = 外包 span

### 5.3 把 `distinct_nuclear_start/end` 误认为簇范围

当前不是。

当前：

- `distinct_nuclear_start/end` = merged breakpoint 对
- `nuclear_cluster_range_start/end` = 簇范围

### 5.4 忽略 `eligible_for_main_analysis`

这是主分析入口的硬条件。

如果你在别处直接重算 distinct，而没有先过滤它，就会把复杂事件重新混回主结果。

### 5.5 用 `carrier_count_sum / N` 解释为携带比例

这是错误的。

真正的携带比例看：

- `carrier_fraction`
- `carrier_ratio_by_frequency_class.tsv`

### 5.6 拿 mt cumulative count 直接比组间富集

当前应该优先比较：

- `overall_cumulative_fraction`
- `cumulative_fraction`

### 5.7 忘记 step 1 缓存

如果你改了 `step 1` 逻辑，但 `00-qc` 已经存在，主控会跳过 `step 1`。

修改这部分代码后，记得：

1. 删 `output/00-qc/` 里的缓存文件，或者
2. 至少删掉 `output/cache/length_bed.full_region_summary.tsv`

否则你会误以为代码没生效。

---

## 6. 推荐的检查顺序

当未来需要改代码时，推荐顺序：

1. 先看 `pipe/01-run_numts_frequency.sh`
2. 再看 `python/01` 到 `python/09`
3. 真正改逻辑时再看 `python/numts_frequency/*.py`
4. 修改核心逻辑后先跑：

```bash
PYTHONPATH=python python/90_regression_cases.py
```

5. 再跑正式流程
6. 最后检查：

- `output/00-qc/06-output_validation_report.md`


