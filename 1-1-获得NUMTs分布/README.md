# README.md

本文件夹负责从`BAM`文件里提取出 NUMTs 。

输入：样本的`BAM`文件。
输出：多个`TSV`文件、`psl`、`sam`、`fasta`文件，记录了NUMTs检测流程中不同阶段的结果。


# NUMTs 检测流程结果解读

本文根据以下代码整理而成，不修改任何代码，只解释结果文件是什么意思：

- `pipe/0_NUMTs_detection_2.0_单样本.sh`
- `script/0_1_找聚类.py`
- `script/0_2_sam2psl.py`
- `script/0_3_找断点.py`

重点是说明这个流程产出的多个 `tsv` 文件分别记录了什么，以及每一列在生物学上或流程上代表什么。

> 说明：当前 `example/output` 目录里的文件可能不齐全。本文以代码里实际定义的输出为准，而不是只以示例目录中“目前看得到的文件”为准。

按代码实际会出现的 `old` 类 TSV，至少包括下面两种：

- `*.mt.disc.sam.cluster.old.tsv`：样本级、聚类阶段的旧版兼容输出
- `PREFIX.Breakpoints.old.tsv`：候选区域级、断点阶段的旧版兼容输出

---

## 1. 先用一句话理解这套流程

这套流程做的事情可以概括成：

1. 先从 WGS 的 BAM 里找出“一个末端像在线粒体、另一个末端像在核基因组”的 reads。
2. 把这些 reads 在核基因组上的落点聚成一簇，得到“疑似 NUMT 插入区域”。
3. 提取这些候选区域里的 reads 序列，重新做一次全基因组比对。
4. 根据重新比对的结果，统计核基因组断点和线粒体断点。
5. 输出“所有断点”和“高可信断点”两套结果。

如果只想快速看结果，通常优先看这 4 类文件：

- `*.cluster.tsv`：候选 NUMT 区域总表
- `*.breakpointINPUT.tsv`：后续断点分析的输入区域
- `*.AllBreakpoints.tsv`：这个区域里所有检测到的断点
- `*.ConfidentBreakpoints.tsv`：可信度更高的断点

---

## 2. 文件是怎么一步步生成的

```text
输入：样本 BAM
  -> 提取线粒体相关 discordant / split reads
  -> 生成：
     *.mt.disc.sam
     *.mt.split.sam

  -> 聚类脚本 0_1_找聚类.py
  -> 生成：
     *.mt.disc.sam.cluster.tsv
     *.mt.disc.sam.cluster.summary.tsv
     *.mt.disc.sam.breakpointINPUT.tsv
     *.mt.disc.sam.cluster.old.tsv

  -> 提取所有候选区域序列
  -> 生成：
     *_all_regions.fasta

  -> bwa mem 重新比对
  -> 生成：
     *_all_regions.bwa.sam

  -> 0_2_sam2psl.py 转 PSL
  -> 生成：
     *_all_regions.psl

  -> 按每个候选区域拆分 PSL，再跑 0_3_找断点.py
  -> 每个候选区域各自产生：
     PREFIX.psl
     PREFIX.Breakpoints.old.tsv
     PREFIX.chr*.AllBreakpoints.tsv
     PREFIX.chr*.ConfidentBreakpoints.tsv
```

其中 `PREFIX` 的格式是：

```text
{sampleID}_{chr}.{start}.{end}
例如：
HG001_chr1.633396.634413
```

注意：

- `cluster.old.tsv` 是按样本输出 1 个
- `Breakpoints.old.tsv` 是按候选区域输出多个
- 如果某个候选区域拆分出的 `PSL` 没有真实命中，shell 脚本会直接跳过该区域，此时对应的 `PREFIX.Breakpoints.old.tsv` 可能不会生成

---

## 3. 这些文件可以分成哪几类

### 3.1 候选区域层面的文件

这类文件回答的问题是：

- 样本里找到了哪些疑似 NUMT 区域？
- 这些区域是由哪些 discordant reads 支持的？

对应文件：

- `*.cluster.tsv`
- `*.cluster.summary.tsv`
- `*.breakpointINPUT.tsv`
- `*.cluster.old.tsv`

### 3.2 断点层面的文件

这类文件回答的问题是：

- 在某一个候选区域里，具体可能断在什么位置？
- 哪些断点只是“都看到了”，哪些断点是“核基因组 hit 和线粒体 hit 能对上同一批 reads”，因此更可信？

对应文件：

- `*.AllBreakpoints.tsv`
- `*.ConfidentBreakpoints.tsv`
- `*.Breakpoints.old.tsv`

---

## 4. 每个 TSV 文件的详细解释

## 4.1 `*.breakpointINPUT.tsv`

### 这个文件是什么

这是断点分析的“区域清单”。  
每一行就是一个候选 NUMT 区域，后面的 shell 脚本会按这一行给出的染色体和坐标去提取序列、拆分 PSL、调用断点脚本。

这个文件没有表头。

### 每列是什么意思

| 列序号 | 列名 | 通俗解释 |
|---|---|---|
| 1 | `IndividualID` | 样本 ID，例如 `HG001`。 |
| 2 | `Cluster_No` | 这个候选区域里，聚在一起的核基因组侧 discordant read 落点有多少个。它不是“第几个 cluster 编号”，而是“这个 cluster 里有多少个位置点”。 |
| 3 | `disFile` | discordant SAM 文件路径，即 `*.mt.disc.sam`。 |
| 4 | `splitFile` | split-read SAM 文件路径，即 `*.mt.split.sam`。当前新版主流程里这个字段主要是作为记录信息保留下来。 |
| 5 | `wgsBAM` | 原始 WGS BAM 文件路径。后面提取候选区域 reads 时会用到。 |
| 6 | `chr` | 候选 NUMT 区域所在的核染色体，例如 `chr1`。 |
| 7 | `start` | 候选区域起点。代码里是在该 cluster 最小 read 位置的基础上再往左扩 500 bp。 |
| 8 | `end` | 候选区域终点。代码里是在该 cluster 最大 read 位置的基础上再往右扩 500 bp。 |

### 例子

```text
HG001    3    .../HG001.mt.disc.sam    .../HG001.mt.split.sam    .../HG001.bam    chr1    633396    634413
```

这行的意思是：

- 样本 `HG001`
- 在 `chr1:633396-634413` 这个候选区域
- 有 3 个核侧 discordant 落点支持它
- 后续分析要用到对应的 `disc`、`split` 和原始 `bam`

---

## 4.2 `*.cluster.tsv`

### 这个文件是什么

这是新版聚类结果的主表。  
每一行对应一个候选 NUMT 区域，比 `breakpointINPUT.tsv` 更详细，因为它把支持这个区域的核侧位置和线粒体侧位置都列出来了。

这个文件有表头。

### 每列是什么意思

| 列名 | 通俗解释 |
|---|---|
| `chr` | 这个候选区域位于哪条核染色体。 |
| `start` | 候选区域左边界，等于 cluster 最小位置减 500。 |
| `end` | 候选区域右边界，等于 cluster 最大位置加 500。 |
| `Cluster_No` | 这个候选区域内，聚在一起的核基因组侧 read 落点数。名字里虽然有 `No`，但它实际是“数量”，不是流水号。 |
| `reads` | 支持这个候选区域的核基因组侧 read 位置列表，也就是这些 discordant reads 在核基因组上的 `POS`。 |
| `mt_positions` | 与这些核侧 reads 配对的线粒体侧位置列表，来自 SAM 中的 `PNEXT`。可以理解为“同一批证据在 chrM 上对应落在哪里”。 |
| `subCluster_No` | `mt_positions` 里有多少个线粒体位置。当前代码里它本质上等于 `mt_positions` 的长度。 |
| `Cluster_ID` | 给这个候选区域拼出来的唯一字符串 ID，格式大致为：`核染色体_核侧最小原始位置_核侧最大原始位置_chrMboth_线粒体最小位置_线粒体最大位置`。 |
| `IndividualID` | 样本 ID。 |

### 例子

```text
chr1    633396    634413    3    [633896, 633900, 633913]    [9735, 9872, 9818]    3    chr1_633896_633913_chrMboth_9735_9872    HG001
```

可以直接读成：

- 在 `chr1:633396-634413` 发现 1 个候选 NUMT 区域
- 核侧支持位置是 `633896, 633900, 633913`
- 这几条 reads 在线粒体上的对应位置是 `9735, 9872, 9818`
- 因为核侧这些位置彼此很近，所以被归为同一个候选插入区域

### 一个容易混淆的点

`start/end` 是“扩展后的候选区域边界”，而 `Cluster_ID` 里放的是“扩展前的原始 cluster 核心范围”。  
所以你会看到：

- `start=633396, end=634413`
- 但 `Cluster_ID` 里是 `633896_633913`

这是正常的，不是写错了。

---

## 4.3 `*.cluster.summary.tsv`

### 这个文件是什么

这是 `cluster.tsv` 的一个简化汇总表。  
它没有表头，而且第一列是 pandas 自动写出来的行号。

这个文件最容易误读，因为第一列看起来像“统计值”，其实不是，它只是行索引。

### 每列是什么意思

| 列序号 | 对应内容 | 通俗解释 |
|---|---|---|
| 1 | pandas index | 纯粹是表格行号，例如 `0, 1, 2`，没有生物学意义。 |
| 2 | `IndividualID` | 样本 ID。 |
| 3 | `Cluster_ID` | 候选区域 ID。 |
| 4 | `Cluster_No` | 核侧聚类中包含的 read 落点数。 |
| 5 | `subCluster_No` | 线粒体侧位置数。 |
| 6 | `size` | 对相同 `IndividualID + Cluster_ID + Cluster_No + subCluster_No` 组合做 `groupby` 后的行数。按当前新版输出逻辑，通常每个 cluster 只占 1 行，所以这里常常是 `1`。它不是 reads 数，也不是 cluster 编号。 |

### 例子

```text
0    HG001    chr11_49861875_49862357_chrMboth_16089_16340    2    2    1
```

这行的意思是：

- 第 0 行
- 样本 `HG001`
- 这个 cluster 的 ID 是 `chr11_49861875_49862357_chrMboth_16089_16340`
- 核侧有 2 个支持位置
- 线粒体侧有 2 个位置
- 这个组合在汇总表里出现了 1 次

---

## 4.4 `*.cluster.old.tsv`

### 这个文件是什么

这是旧版兼容格式的聚类结果。  
和 `cluster.tsv` 不同，它不是“一行一个候选区域”，而是更接近“一行一条 read 记录”。  
前面大部分列直接来自 `disc SAM`，后面再附加几个 cluster 字段。

这个文件有表头。

### 每列是什么意思

| 列名 | 通俗解释 |
|---|---|
| `QNAME` | read 名称。 |
| `FLAG` | SAM flag，描述配对状态、方向、是否 supplementary 等。 |
| `RNAME` | 这条记录当前比对到的参考序列名称。这里通常是核染色体。 |
| `POS` | 这条 read 在 `RNAME` 上的起始比对位置。 |
| `MAPQ` | 比对质量。数值越大通常代表越可信。 |
| `CIGAR` | 比对模式，比如匹配、软剪切、插入缺失等。 |
| `RNEXT` | mate read 比对到的参考序列名称。这里常见是 `chrM`。 |
| `PNEXT` | mate read 在 `RNEXT` 上的比对位置。这里常作为线粒体对应位置。 |
| `TLEN` | 模板长度字段。 |
| `SEQ` | read 序列。 |
| `QUAL` | 碱基质量字符串。 |
| `SM` | 这里是按“第 12 列可选标签”读进来的占位列，常见情况下可能存放某个 SAM optional tag，也可能为空。 |
| `RG` | 同上，按固定列位读取的可选标签占位列，常见可表示 read group，但是否真是 `RG` 取决于输入 SAM。 |
| `NM` | 同上，常见可表示编辑距离 tag，但这里本质上仍是按列位读取的原始可选字段。 |
| `BC` | 同上，常见可表示 barcode tag。 |
| `OC` | 同上，常见可表示 original CIGAR tag。 |
| `ZX` | 同上，自定义或其他 optional tag 占位列。 |
| `ZY` | 同上，自定义或其他 optional tag 占位列。 |
| `SA` | 同上，常见可表示 supplementary alignment tag。 |
| `subCluster_No` | 线粒体侧子聚类大小。旧版逻辑里通常表示该条 read 所属线粒体位置簇包含多少条记录。 |
| `Cluster_No` | 核侧聚类大小，即该条 read 所属核侧 cluster 里有多少个位置点。 |
| `Cluster_ID` | 该条 read 所属 cluster 的 ID。 |
| `IndividualID` | 样本 ID。 |
| `clusterID` | 旧版简化 cluster 名，一般是 `染色体_Cluster_No` 这样的格式。 |

### 怎么理解它

如果你想看“某个候选区域到底是由哪些原始 reads 支持的”，`cluster.old.tsv` 会比 `cluster.tsv` 更细。  
如果你只是想看每个候选区域的概况，优先看 `cluster.tsv` 即可。

---

## 4.5 `*.AllBreakpoints.tsv`

### 这个文件是什么

这是某一个候选区域的“全部断点统计表”。  
只要某个断点类型在该区域的 PSL 结果中被统计到了，就会写进来，因此它偏“全量结果”。

这个文件有表头。

### 每列是什么意思

| 列名 | 通俗解释 |
|---|---|
| `sampleID` | 样本 ID。 |
| `chr` | 断点所在参考序列。可能是核染色体，也可能是 `chrM`。 |
| `strand` | 该命中的比对方向，`+` 或 `-`。这是比对方向，不要直接把它理解成 NUMT 插入方向。 |
| `pointGroup` | 断点类别，是脚本根据 PSL 的 `Qstart/Qend/Tstart/Tend/strand` 规则分出来的标签。 |
| `Group` | 对 `pointGroup` 的简化分组，方便汇总阅读。 |
| `readsCount` | 支持这个断点位置的 reads 数量。 |
| `pos` | 统一后的断点位置坐标。脚本会根据断点类型从 `Tstart` 或 `Tend` 里取一个作为最终位置。 |
| `type` | 结果类型，这个文件里固定为 `all`。 |

### `pointGroup` 应该怎么读

最常见的 4 类是：

| `pointGroup` | 通俗理解 |
|---|---|
| `nu_Tend_Bleft` | 核基因组侧左边界断点。 |
| `nu_Tstart_Bright` | 核基因组侧右边界断点。 |
| `mt_Tend` | 线粒体侧一个方向的断点。 |
| `mt_Tstart` | 线粒体侧另一个方向的断点。 |

可以简单理解成：  
`nu_` 开头的是核基因组端，`mt_` 开头的是线粒体端。

### `Group` 应该怎么读

常见值包括：

| `Group` | 通俗理解 |
|---|---|
| `nuleft` | 核侧左断点 |
| `nuright` | 核侧右断点 |
| `mt_Tend` | 线粒体侧一种端点 |
| `mt_Tstart` | 线粒体侧另一种端点 |

### 例子

```text
HG001    chr11    +    nu_Tend_Bleft      nuleft     1    49862079.0    all
HG001    chrM     -    mt_Tstart          mt_Tstart  1    16088.0       all
```

可以读成：

- 在样本 `HG001` 的这个候选区域里
- 有 1 条 read 支持核侧左断点在 `chr11:49862079`
- 另有 1 条 read 支持线粒体侧断点在 `chrM:16088`

### 这个文件适合做什么

适合初步浏览“这个区域里可能有哪些断点”。  
但它包含的是全量统计，可信度不一定都一样。

---

## 4.6 `*.ConfidentBreakpoints.tsv`

### 这个文件是什么

这是高可信断点表。  
它和 `AllBreakpoints.tsv` 的列完全一样，但不是所有断点都会进来。

更严格地说，代码只保留这类证据：

- 核侧断点：这些 reads 的 `Qname` 同时能在线粒体命中集合里找到
- 线粒体断点：这些 reads 的 `Qname` 同时能在当前核染色体命中集合里找到

也就是说，它要求核侧证据和线粒体侧证据能通过同一批 reads 对上，因此比 `AllBreakpoints.tsv` 更可信。

### 每列是什么意思

和 `AllBreakpoints.tsv` 完全相同：

| 列名 | 通俗解释 |
|---|---|
| `sampleID` | 样本 ID。 |
| `chr` | 断点所在参考序列。 |
| `strand` | 比对方向。 |
| `pointGroup` | 断点类别。 |
| `Group` | 简化分组。 |
| `readsCount` | 支持该断点的 reads 数。 |
| `pos` | 断点位置。 |
| `type` | 在这个文件里固定为 `confident`。 |

### 例子

```text
HG001    chr5    +    nu_Tend_Bleft    nuleft        1    32338477.0    confident
HG001    chrM    +    mt_Tend          mt_confident  1    12821.0       confident
```

可以理解为：

- 这不是“单边看到”的断点
- 而是核基因组侧和线粒体侧能互相对应上的、更可信的断点

### 一个常见现象

有些区域 `AllBreakpoints.tsv` 里有很多行，但 `ConfidentBreakpoints.tsv` 可能是空的。  
这通常表示：

- 这个区域里虽然统计到了不少可能断点
- 但没有足够证据证明核侧 hit 和线粒体侧 hit 来自同一批 read 名称

示例中的 `chr1` 候选区域就是这种情况。

---

## 4.7 `*.Breakpoints.old.tsv`

### 这个文件是什么

这是旧版兼容格式的断点结果。  
它没有表头，列顺序也和新版不一样。

它的内容本质上来自 `ConfidentBreakpoints` 的旧输出方式，因此如果你在做新分析，优先看 `ConfidentBreakpoints.tsv` 更直观。

### 每列是什么意思

| 列序号 | 列名 | 通俗解释 |
|---|---|---|
| 1 | `pointGroup` | 断点类别，例如 `nu_Tend_Bleft`、`mt_Tstart`。 |
| 2 | `Group` | 旧版分组标签。核侧常见 `nuleft`，线粒体侧在旧版里常见 `mtLeft` 或 `mtRight`。 |
| 3 | `chr` | 染色体名称。旧版里线粒体有时会写成 `MT`。 |
| 4 | `Tend` | 目标序列末端坐标。只有部分断点类型这里有意义；没有意义时常被填成 `-1`。 |
| 5 | `strand` | 比对方向。 |
| 6 | `readsCount` | 支持该断点的 reads 数。 |
| 7 | `Tstart` | 目标序列起始坐标。只有部分断点类型这里有意义；没有意义时常被填成 `-1`。 |
| 8 | `sampleID` | 样本 ID。 |

### 为什么会同时有 `Tstart` 和 `Tend`

因为旧版格式没有新版里的统一坐标列 `pos`。  
所以脚本把 `Tstart` 和 `Tend` 都保留下来：

- 某些断点真正该看的位置在 `Tstart`
- 另一些断点真正该看的位置在 `Tend`
- 不适用的那个位置会被填成 `-1`

因此这个旧文件读起来不如新版方便。

---

## 5. 哪些文件最值得优先看

如果你是第一次看这套结果，建议按下面顺序：

### 第一步：看候选区域

先看 `*.cluster.tsv`：

- 有几个候选 NUMT 区域
- 每个区域在什么染色体、什么范围
- 有多少核侧支持位置和线粒体支持位置

### 第二步：看这个区域有没有明确断点

再看每个区域对应的 `*.AllBreakpoints.tsv`：

- 核侧左断点在哪里
- 核侧右断点在哪里
- 线粒体侧断点在哪里

### 第三步：看哪些断点更可信

最后看 `*.ConfidentBreakpoints.tsv`：

- 如果这里也有结果，说明核侧和线粒体侧证据能通过 read 名称互相对上
- 这通常比只看 `AllBreakpoints.tsv` 更可靠

---

## 6. 示例结果怎么读

以示例里的 `HG001_chr11.49861375.49862857` 为例：

### 在 `HG001.mt.disc.sam.cluster.tsv` 里

```text
chr11    49861375    49862857    2    [49861875, 49862357]    [16089, 16340]    2    chr11_49861875_49862357_chrMboth_16089_16340    HG001
```

说明：

- 在 `chr11` 上有一个候选区域
- 扩展后的候选区间是 `49861375-49862857`
- 由 2 个核侧位置支持
- 对应的线粒体位置是 `16089` 和 `16340`

### 在 `HG001_chr11.49861375.49862857.chr11.AllBreakpoints.tsv` 里

```text
HG001    chr11    +    nu_Tend_Bleft      nuleft      1    49862079.0    all
HG001    chrM     -    mt_Tstart          mt_Tstart   1    16088.0       all
```

说明：

- 核基因组端检测到一个左侧断点候选
- 线粒体端也检测到一个断点候选

### 在 `HG001_chr11.49861375.49862857.chr11.ConfidentBreakpoints.tsv` 里

```text
HG001    chr11    +    nu_Tend_Bleft    nuleft        1    49862079.0    confident
HG001    chrM     -    mt_Tstart        mt_confident  1    16088.0       confident
```

说明：

- 这两个断点不是“单独统计到”
- 而是能够由互相对应的 read 证据支持
- 因此这组断点比单纯出现在 `AllBreakpoints.tsv` 中更可信

---

## 7. 最后给一个简化版记忆方法

如果只想记住最关键的区别，可以这么记：

- `cluster.tsv`：告诉你“哪里像是 NUMT 区域”
- `breakpointINPUT.tsv`：告诉流程“后面要分析哪些区域”
- `AllBreakpoints.tsv`：告诉你“这个区域里所有可能的断点”
- `ConfidentBreakpoints.tsv`：告诉你“其中哪些断点证据更硬”
- `*.old.tsv`：旧版兼容输出，能看，但不如新版直观

---

## 8. 一个实用提醒

看到 `Cluster_No`、`subCluster_No` 这种字段时，不要先入为主把它当“编号”。  
在这套代码里，它们更接近“数量”：

- `Cluster_No`：核侧支持位置数量
- `subCluster_No`：线粒体侧支持位置数量

另外，`cluster.summary.tsv` 第一列只是行号，不是统计结果；  
`Breakpoints.old.tsv` 里的 `Tstart/Tend` 也不是同时都有效，通常只会有一个是真正该看的断点坐标。
