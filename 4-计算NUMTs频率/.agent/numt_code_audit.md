# NUMT 频率与分布流程代码核查说明

## 1. 目的

本文档用于把**文献中的原始做法**与**当前代码实现**逐项对照，定位最可能导致结果异常的环节，并给出：

1. 需要先核查什么
2. 核查时应看什么证据
3. 若确认有误，代码应如何修改

当前重点只放在以下 4 类问题：

- 长度定义
- distinct NUMT 合并规则
- 频率与携带者统计
- mtDNA 断点累计曲线

不涉及外部注释、known 和 novel、cohort-specific、富集分析。

---

## 2. 文献中的原做法

### 2.1 Wei 2022 与 NyuWa 2025 的上游检测与合并框架

两篇文章在 NUMT 检测和跨样本合并上基本是一条线。NyuWa 明确说明其 NUMT 分析直接沿用了 Wei 2022 的公共 pipeline。

**原始规则要点：**

1. 先从 WGS 中提取 discordant reads 和 split reads。
2. 将 **500 bp 内的 discordant reads 作为一个 cluster**。
3. **少于 5 对 discordant reads 的 cluster 去除**。
4. 在 cluster 两侧 **500 bp flanks** 内搜索 split reads，并用 **BLAT** 重新比对。
5. 用 split reads 来确定 **nuclear breakpoint 和 mtDNA breakpoint**。
6. **来自不同样本、相距 1000 bp 内的 NUMT** 合并为同一个 NUMT。
7. NyuWa 明确说明：**不考虑 concatenated NUMTs 和其他复杂情况**。

这几点在 NyuWa 方法部分写得很清楚，且其代码来源直接指向 Wei 2022 的公开 pipeline。  

### 2.2 文献对 distinct NUMT 的理解

文献中的 distinct NUMT 不是“按窗口粗分箱”的概念，而是：

- 先由读段证据确定一个样本内的 insertion 和断点；
- 再跨样本基于**断点位置相近**进行聚合；
- 最终得到一个“跨样本去重后的插入位点”。

因此，**断点本身**比“简单中点”更接近文献中的主对象。

### 2.3 文献中的 mtDNA 断点和长度

Wei 2022 正文明确写到：NUMT 的长度是基于 **两端 mtDNA breakpoints** 定义出来的。文中给出的 germline NUMT 长度分布为：中位数 156 bp，77.8% 小于 500 bp。NyuWa 也给出类似结论：中位数 119 bp，89.15% 小于 500 bp。  

因此，文献中的长度本质上应尽量对应于：

- 一个插入事件对应的 mtDNA 片段长度；
- 或一个 distinct NUMT 中该片段长度的代表值；
- 而不是把多个不连续 mt 片段的最小起点和最大终点直接外包成一整段。

### 2.4 Zhu 2026 对跨样本合并的更明确写法

Zhu 2026 把“跨样本合并”写得更具体：

- **核断点**：1000 bp 内聚类；
- **合并后的核断点**：取簇内插入的第一和第二 start positions 作为 merged nuclear breakpoint 的起止；
- **mtDNA 断点**：在簇内统计每个 mtDNA breakpoint 的出现频率，**取最高频的 breakpoint** 作为 merged mt breakpoint。

这个写法比“直接取中点再求中位数”更接近文献语义。

### 2.5 文献中的频率分层与携带比例

Wei 2022 与 NyuWa 都使用：

- common：$F \ge 1\%$
- rare：$0.1\% \le F < 1\%$
- ultra-rare：$F < 0.1\%$
- private：单独统计，NyuWa 为仅 1 个个体，Wei 2022 为仅 1 个家系

并且文中的“某类 NUMT 的携带比例”是指：

- **至少携带 1 个该类 NUMT 的样本数 / 总样本数**

而不是把该类所有 NUMT 的 carrier_count 直接相加后再除以样本总数。

### 2.6 文献中的 mtDNA 累积分布

NyuWa 方法部分明确说：对 mtDNA 断点分布使用 **cumulative distribution** 来观察不同位置的斜率，反映断点富集程度。  

这意味着它更接近：

- 沿 mtDNA 坐标的**累计比例曲线**，或
- 至少是不同组之间可比较的标准化累计曲线。

如果只输出累计计数，不同类别之间就不方便直接比较。

---

## 3. 当前代码的实际做法

下面只写与结果偏差最相关的部分。

### 3.1 `io_utils.py::summarize_length_bed()`

当前实现对同一 `region_key`：

- 收集所有 bed 行；
- 取 `mt_start` 的最小值；
- 取 `mt_end` 的最大值；
- 计算 `mt_source_span_bp = max_end - min_start + 1`。

**含义：**
它实际上把一个 region_key 下面出现的所有 mt 片段，直接用一个“最外层包络区间”表示。

**风险：**
如果一个 region_key 下有多个不连续 mt 片段，这种写法会显著放大长度。

---

### 3.2 `catalog_utils.py::build_distinct_catalog()`

当前实现先为每个 sample event 定义两个用于合并的位置：

- `merge_nuclear_pos = nuclear_anchor_midpoint`，缺失时退回 `nuclear_core_midpoint`
- `merge_mt_pos = (mt_breakpoint_start + mt_breakpoint_end) / 2`，缺失时退回 `mt_core_midpoint`

然后：

1. 先在每条染色体内按 `merge_nuclear_pos` 做 1000 bp 单链聚类；
2. 再在核簇内部按 `merge_mt_pos` 做 1000 bp 单链聚类；
3. 用 `combined_cluster_key` 形成 distinct NUMT；
4. distinct 的起止坐标再由这些中点的最小值、最大值、中位数来汇总。

**含义：**
当前 distinct catalog 的核心不是“断点”，而是“中点的两层聚类结果”。

**风险：**
这与文献中“基于 breakpoints 合并”的语义并不完全一致。

---

### 3.3 `annotation_utils.py::annotate_catalog()`

当前 distinct 长度是：

- 先把 sample event 的 `mt_source_span_bp` 合并到 distinct；
- 然后对每个 distinct 取 `median_mt_source_span_bp`。

**含义：**
如果 sample event 层的长度本来已经被高估，那么 distinct 层只是把这个高估再做一次中位数汇总。

---

### 3.4 `stats_utils.py::compute_main_results()`

当前 `frequency_class_summary` 中：

- `carrier_total = sum(carrier_count)`
- `carrier_fraction = carrier_total / N`

这实际上不是“携带比例”，而是“该类 NUMT 的载体总和对样本数的比值”。

例如一个样本同时带多个 rare NUMT，会被重复计入多次。

但同一个文件里，`carrier_ratio_by_frequency_class` 又是按：

- 每类中至少带 1 个 NUMT 的样本数 / 总样本数

这个才是文献常用语义。

因此，当前同一套输出中混入了两种不同定义。

---

### 3.5 `stats_utils.py::compute_mt_breakpoint_tables()`

当前 `mt_breakpoint_cumulative`：

- 先按位置计数；
- 再直接对 count 做 `cumsum()`。

输出是：

- `overall_cumulative_count`
- `cumulative_count`

**含义：**
这是累计计数，不是累计比例。

**风险：**
不同频率类别之间不能直接比较斜率，和文献中用 cumulative distribution 比较 enrichment 的口径并不一致。

---

## 4. 文献做法与当前代码的关键差异

| 主题 | 文献原做法 | 当前代码做法 | 主要风险 |
| --- | --- | --- | --- |
| 样本内事件基础 | 基于 discordant 和 split reads 定断点 | 部分回退到 reads 和 mt_positions 的范围中点 | 断点信息被中点化 |
| 跨样本合并 | 基于核与 mt 断点相近来合并 | 基于核中点和 mt 中点做两层单链聚类 | 可能把真实事件拆碎，或把不同事件并错 |
| merged mt breakpoint | Zhu 2026 取簇内最高频 breakpoint | 当前代码取 breakpoint 中点，再取中位数 | 代表性不足 |
| 长度定义 | 尽量由两端 mt 断点定义片段长度 | 对 bed 多行取最外层 span | 长度可能被显著拉长 |
| 复杂 NUMT | NyuWa 明确不纳入主分析 | 当前代码未单独标记复杂候选 | 复杂事件混入主统计 |
| 携带比例 | 至少携带 1 个该类事件的样本比例 | 一处按样本去重，一处按 carrier_count 求和 | 容易误读 |
| mt 累积分布 | cumulative distribution | cumulative count | 组间不可直接比较 |

---

## 5. 需要优先核查的地方

下面按优先级排序。

### P0. 长度定义是否被错误放大

#### 为什么先查

如果长度定义错了，会直接把：

- 长度分布
- 短片段比例
- common 和 rare 的长度比较
- 各分层长度图

全部带偏。

#### 应核查什么

对 20 到 50 个随机 region_key 检查：

1. length bed 中该 region_key 是否对应多行；
2. 这些行在 mtDNA 上是否连续；
3. 当前输出的 `mt_source_span_bp` 是否等于“最外层包络长度”；
4. 是否存在明显不连续但被算成一整段的情况；
5. 是否出现大量小数长度，或大量远大于文献常见范围的长度。

#### 若确认有误，建议修改

**修改原则：不要再用最小起点到最大终点的外包 span 作为默认长度。**

建议分两种情况：

**情况 A：同一 region_key 下的 mt 区间彼此重叠或相邻**

- 先对区间做 merge；
- 再用合并后区间的**并集长度**作为 `mt_source_span_bp`。

**情况 B：同一 region_key 下存在明显不连续区间**

- 标记为 `complex_mt_interval = True`；
- 默认不直接进入主长度统计；
- 若主分析必须保留，则使用**主片段长度**，例如支持最多读段的片段，或与 breakpoint 最一致的片段；
- 同时单独输出复杂候选清单供人工核查。

如果你的主分析目标是尽量贴近 Wei 2022 和 NyuWa，那么更保守的做法是：

- **复杂和疑似 concatenated 事件单独标记，不纳入主长度分布。**

---

### P0. distinct 合并是否被“中点化”得过头

#### 为什么先查

如果合并规则错了，会直接造成：

- distinct NUMT 总数虚高或虚低；
- 每样本 distinct NUMT 数异常；
- private、rare、ultra-rare 比例异常；
- 核和 mt 断点图失真。

#### 应核查什么

针对 20 到 30 个 carrier 数较高、以及 20 到 30 个 private 的 distinct NUMT：

1. 回到原始 `sample_events`，检查其 `nuclear_left_breakpoint_pos`、`nuclear_right_breakpoint_pos`、`mt_breakpoint_start`、`mt_breakpoint_end`；
2. 看同一个 distinct 中，是否只是“中点接近”，但真实 breakpoints 差别很大；
3. 看同一个 1000 bp 核窗口内是否存在多个稳定的 mt breakpoint 簇；
4. 看当前 distinct 的 `distinct_nuclear_start` 与 `distinct_nuclear_end` 是否只是中点的范围，而非真实 breakpoint 边界；
5. 看同一样本是否有多个 sample event 被并进同一个 distinct，且本质上像不同插入。

#### 若确认有误，建议修改

**修改原则：用 breakpoints 合并，不要用 breakpoints 的中点替代 breakpoints。**

建议改法：

1. **核侧代表位置**
   - 优先用 `nuclear_left_breakpoint_pos` 和 `nuclear_right_breakpoint_pos`；
   - 若两端都存在，可把两端分别保留；
   - 跨样本聚类时，至少应基于真实 breakpoint，而不是 `nuclear_anchor_midpoint`。

2. **mt 侧代表位置**
   - 优先用 `mt_breakpoint_start` 和 `mt_breakpoint_end`；
   - 若一个核簇内有多个 mt breakpoint 候选，参考 Zhu 2026，取**簇内最常见的 mt breakpoint** 作为 merged mt breakpoint；
   - 不建议直接取 `(start + end) / 2` 的中点作为代表。

3. **distinct 起止坐标**
   - `distinct_nuclear_start` 与 `distinct_nuclear_end` 应尽量来源于真实断点边界；
   - 若采用 Zhu 2026 风格，可在核簇内使用第一和第二 start positions 作为 merged nuclear breakpoint 边界；
   - 至少不要再把“中点的最小值和最大值”当成断点边界。

4. **复杂候选单独标记**
   - 若同一核簇内存在多个稳定 mt 簇，建议增加字段：
     - `multi_mt_cluster_in_nuclear_window`
     - `complex_event_flag`
   - 默认不进入主 distinct 统计，或单独列出。

---

### P0. `carrier_fraction` 是否定义错误

#### 为什么先查

这会直接影响你对“某类 NUMT 在多少样本中出现”的解释。

#### 应核查什么

检查 `frequency_class_summary.tsv` 中：

- `carrier_fraction`
- `carrier_total`

并与 `carrier_ratio_by_frequency_class.tsv` 对照：

- 若二者差很多，说明两个文件的统计语义不同；
- 若 `carrier_fraction > 1`，那它就不应再叫“fraction”。

#### 若确认有误，建议修改

二选一：

**方案 A：保留当前求和结果，但改名**

把：

- `carrier_total`
- `carrier_fraction`

改为：

- `carrier_count_sum`
- `mean_numt_count_per_sample`

这样至少语义不误导。

**方案 B：改成真正的携带比例，推荐**

对每个 `frequency_class`：

1. 从 `sample_map` 与 `catalog` 联表；
2. 对 `sample_id + frequency_class` 去重；
3. 统计至少携带 1 个该类 NUMT 的样本数；
4. 除以总样本数。

然后把这列保留为：

- `carrier_fraction`

同时保留 `carrier_count_sum` 作为辅助描述。

---

### P1. mt 断点累计曲线是否缺乏标准化

#### 为什么要查

如果你后面要比较：

- common 与 rare
- 各地区
- 各性别
- 各民族

那么累计计数曲线没有标准化会很难解释。

#### 应核查什么

查看 `mt_breakpoint_cumulative.tsv`：

- 是否只有累计计数；
- 是否没有 `cumulative_fraction` 或 `ecdf`；
- 不同频率类的总断点数差很多时，曲线是否只是“数量高低”，而不是“分布形状”。

#### 若确认有误，建议修改

在现有累计计数基础上再增加：

- `overall_cumulative_fraction = overall_cumulative_count / overall_total`
- `cumulative_fraction = cumulative_count / class_total`

其中 `overall_total` 与 `class_total` 都应是各自纳入统计的 distinct NUMT 总数或 breakpoint 总数。

后续绘图时优先用 `cumulative_fraction`，而不是 `cumulative_count`。

---

### P1. 旧脚本中的固定窗口法是否仍在被使用

压缩包中的旧脚本 `script/3-统计频率.py` 用的是：

- `mid_pos // 1000`
- `chr:window`

来定义 distinct NUMT。

这属于**固定 1000 bp 分箱**，不是文献中的断点聚类逻辑。

#### 应核查什么

1. 当前正式流程是否仍调用这份旧脚本；
2. 旧脚本输出是否混入最终结果目录；
3. 是否存在“新旧脚本结果混用”的情况。

#### 若确认有误，建议修改

- 明确把旧脚本移出正式流程；
- 在主脚本中只保留当前新版模块化流程；
- 在 README 中标注旧脚本为 deprecated。

---

## 6. 建议的修改顺序

### 第一步：先修长度

优先修改：

- `numts_frequency/io_utils.py::summarize_length_bed()`

建议：

1. 收集同一 `region_key` 的所有 mt 区间；
2. 先判断是否重叠或近邻；
3. 对连续区间求并集长度；
4. 对不连续区间加 `complex_mt_interval` 标志；
5. 主分析优先使用“并集长度”或“主片段长度”，不要再用最外层包络长度。

### 第二步：再修 distinct 合并

优先修改：

- `numts_frequency/catalog_utils.py::build_distinct_catalog()`

建议：

1. 用真实核断点代替 `merge_nuclear_pos`；
2. 用真实 mt breakpoint 代替 `merge_mt_pos`；
3. 核簇内 mt breakpoint 用“最高频 breakpoint”表示；
4. 增加复杂事件标记；
5. 输出 distinct 时保留：
   - merged nuclear breakpoint start 和 end
   - merged mt breakpoint
   - cluster 内事件数
   - cluster 内样本数
   - 是否复杂

### 第三步：统一频率与携带者统计

优先修改：

- `numts_frequency/stats_utils.py::compute_main_results()`

建议：

1. `carrier_fraction` 统一改为“至少携带 1 个该类 NUMT 的样本比例”；
2. 若需要现有求和结果，则单独命名为 `carrier_count_sum`；
3. 确保主文图和补充图都使用同一语义。

### 第四步：补标准化累计曲线

优先修改：

- `numts_frequency/stats_utils.py::compute_mt_breakpoint_tables()`

建议：

1. 保留累计计数；
2. 新增累计比例；
3. 画图时默认用累计比例。

---

## 7. 最小核查清单

在真正改代码前，建议先做一个最小人工核查。

### 7.1 长度核查

随机抽 30 个 distinct NUMT：

- 看其对应 sample event 的 `mt_source_span_bp`
- 对照原始 length bed 行
- 判断是否存在“多个不连续片段被包成一整段”

### 7.2 合并核查

随机抽：

- 10 个 prevalent 或 common
- 10 个 rare
- 20 个 private

检查其：

- 核侧 breakpoints 是否真的相近
- mt 侧 breakpoints 是否真的属于同一簇
- 是否只是中点接近

### 7.3 统计核查

检查：

- `frequency_class_summary.tsv`
- `carrier_ratio_by_frequency_class.tsv`

确认两个文件中的“比例”定义是否一致。

### 7.4 结果合理性核查

把修正前后的结果与文献做数量级比对：

- 平均每样本 distinct NUMT 数
- distinct NUMT 总数
- 长度中位数
- 小于 500 bp 的比例
- ultra-rare 和 private 的携带率

如果修正后这些指标明显向 Wei 2022 和 NyuWa 靠拢，说明方向是对的。

---

## 8. 结论

当前代码最可能影响结果真实性的，不是画图，而是以下三点：

1. **长度汇总把多段 mt 区间外包成一个大 span**；
2. **distinct 合并过度依赖中点，而不是断点本身**；
3. **携带比例与累计曲线的统计口径不够统一**。

因此，建议先按如下顺序处理：

1. 修 `summarize_length_bed()`
2. 修 `build_distinct_catalog()`
3. 修 `compute_main_results()`
4. 修 `compute_mt_breakpoint_tables()`

在这 4 处修完并复核之前，**不建议把当前结果当作最终 NUMT landscape 结果解释。**

---

## 9. 参考文献

1. Wei W, et al. Nature, 2022. Nuclear-embedded mitochondrial DNA sequences in 66,083 human genomes.
2. Wang Y, et al. bioRxiv, 2025. Mitochondrial Genome Variants and Nuclear Mitochondrial DNA Segments in 7331 Individuals from NyuWa and 1KGP.
3. Zhu Q, et al. Science Advances, 2026. Introgressed mitochondrial fragments from archaic hominins alter nuclear genome function in modern humans.
