# AGENTS.md

## 目的

本文档用于告诉后续人工智能或脚本：本流程不同结果文件后缀分别表示什么生物学含义，避免编码时误读。

如文档与代码不一致，以以下脚本为准：

- `pipe/0_NUMTs_detection_2.0_单样本.sh`
- `script/0_1_找聚类.py`
- `script/0_2_sam2psl.py`
- `script/0_3_找断点.py`

## 总原则

1. 候选区域不等于已确认 NUMT。
2. `cluster` 类文件表示候选区域及支持证据，不表示最终断点。
3. `AllBreakpoints` 是全量断点信号，包含较弱证据。
4. `ConfidentBreakpoints` 是更高可信断点，优先用于结果汇总。
5. `Cluster_No` 是支持数，不是编号。
6. `strand` 是比对方向，不能直接当作插入方向。
7. `start/end` 是扩展后的候选区间，不是精确插入位点。
8. `old` 文件是兼容旧格式，非首选结果。

## 各类文件的生物学意义

### `*.mt.disc.sam`

线粒体相关的 不一致 读段 原始证据。

只能表示“可能存在线粒体-核基因组连接信号”，不能直接当作 NUMT 结果。

### `*.mt.split.sam`

线粒体相关的分裂读段原始证据。

只能作为辅助证据，不能单独视为已确认插入。

### `*.mt.disc.sam.cluster.tsv`

主候选区域表。

每行表示一个疑似 NUMT 相关核基因组区域，汇总了：

- 核侧染色体位置
- 核侧支持位点
- 线粒体侧对应位置
- 支持数量

编码注意：

- 用它统计候选区域。
- 不要把每一行写成“一个已确认 NUMT”。
- `Cluster_No` 是支持位点数量，不是流水号。

### `*.mt.disc.sam.cluster.summary.tsv`

候选区域简表，用于快速查看和计数。

第一列只是表格行号，没有生物学意义。需要详细解释时优先看 `cluster.tsv`。

### `*.mt.disc.sam.breakpointINPUT.tsv`

断点分析输入区域列表。

每行表示一个候选区域将被送去做后续断点识别。它是流程输入清单，不是最终断点结果。

### `*_all_regions.fasta`

所有候选区域提取出的序列集合。

中间文件，无直接最终生物学结论。

### `*_all_regions.bwa.sam`

所有候选区域重新比对后的 SAM 结果。

中间文件，用于后续断点判断。

### `*_all_regions.psl`

重新比对结果转成的 PSL 文件。

中间文件，用于断点识别逻辑。

### `PREFIX.psl`

单个候选区域对应的 PSL 文件。

一个文件对应一个候选区域，用于该区域的断点识别。

### `PREFIX.chr*.AllBreakpoints.tsv`

单个候选区域内检测到的全部断点信号。

包括：

- 核侧断点信号
- 线粒体侧断点信号

常见理解：

- `nu_` 开头：核基因组侧断点
- `mt_` 开头：线粒体侧断点
- `读段Count`：支持该断点信号的 读段 数
- `pos`：统一后的断点位置

编码注意：

- 这是全量结果，不等于高可信结果。
- 同一区域多行不等于有多个真实 NUMT。

### `PREFIX.chr*.ConfidentBreakpoints.tsv`

单个候选区域的高可信断点结果。

只有当核侧证据和线粒体侧证据能通过同一批 读段 或匹配证据对应起来时，才会保留下来。

编码注意：

- 结果汇总优先使用这个文件。
- 如果该文件为空，不能把 `AllBreakpoints` 强行当作高可信结论。

### `*.mt.disc.sam.cluster.old.tsv`

聚类阶段旧版兼容输出。

更接近 read 级记录，不适合作为首选候选区域总表。优先使用 `*.mt.disc.sam.cluster.tsv`。

### `PREFIX.Breakpoints.old.tsv`

断点阶段旧版兼容输出。

优先使用 `PREFIX.chr*.ConfidentBreakpoints.tsv`，不要和新版结果直接混用。

## 编码时的优先顺序

1. 候选区域：`*.mt.disc.sam.cluster.tsv`
2. 流程输入区间：`*.mt.disc.sam.breakpointINPUT.tsv`
3. 全量断点：`PREFIX.chr*.AllBreakpoints.tsv`
4. 高可信断点：`PREFIX.chr*.ConfidentBreakpoints.tsv`
5. 旧版兼容结果：仅在必须兼容旧流程时使用

## 推荐表述

建议使用：

- 候选 NUMT 区域
- 疑似断点
- 高可信断点
- 支持 读段 数
- 核侧断点信号
- 线粒体侧断点信号

避免使用：

- 已确认 NUMT
- 精确插入位点
- 插入方向

除非代码或外部验证明确支持这些结论。

## HG001 示例

`output/HG001` 中共有 3 个候选区域：

- `chr1:633396-634413`
- `chr11:49861375-49862857`
- `chr5:32337638-32339095`

其中：

- `chr1` 有 `AllBreakpoints`，但没有 `ConfidentBreakpoints`，应视为只有候选信号，没有高可信断点支持。
- `chr11` 有高可信断点，证据强于 `chr1`。
- `chr5` 有高可信断点，证据强于仅出现在 `AllBreakpoints` 的区域。

## 最终约束

后续人工智能写代码时，禁止把所有后缀文件都当成同一种“NUMT结果”。

至少必须区分：

- 原始证据
- 候选区域
- 全量断点
- 高可信断点
- 旧版兼容结果
