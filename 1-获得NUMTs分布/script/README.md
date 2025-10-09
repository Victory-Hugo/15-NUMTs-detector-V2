# 0_1_找聚类.py
该代码输出3个结果文件。
## cluster.tsv
“每一个核→MT 不匹配聚类”的原始记录。每一列的内容如下所示：
1. chr：核染色体名称（1–22, X 或 Y），表示这组读段所在的核染色体。
2. start：核染色体聚类区间的起始位置（以 bp 为单位），等于 min(POS) - 500。
3. end：核染色体聚类区间的结束位置，等于 max(POS) + 500。
4. Cluster_No：该聚簇中原始核染色体端读段的数量（聚簇大小）。
5. reads：列表形式，包含参与该聚簇的所有核端 POS 值。
6. mt_positions：列表形式，包含与这些核端读段配对的 MT 端 PNEXT 值（即落在线粒体上的配对位置）。
7. subCluster_No：对上面 mt_positions 再聚簇（maxgap=500）后，保留下来的 MT 子聚簇大小。
8. Cluster_ID：一个文本 ID，格式为：用来唯一标识这一对“核区间↔MT区间”关系。
```
{chr}_{原POS最小}_{原POS最大}_MTboth_{mtPOS最小}_{mtPOS最大}
```
9. IndividualID：样本标识（脚本的第二个参数），例如 Y21100000492606.deduped。

## cluster.summary.tsv
这是对上面 .cluster.tsv 的汇总：按照
```
[IndividualID, Cluster_ID, Cluster_No, subCluster_No]
```
1. 行索引：脚本自动生成的行号（无实际意义，只是为了区分）。
2. IndividualID：同上。
3. Cluster_ID：同上，用来唯一标记一个“核→MT”对。
4. Cluster_No：同上。
5. subCluster_No：同上。
6. size：在 regions 列表中，对应这一组合的记录数目（理论上通常为 1，除非后续逻辑合并了相同 ID 的多行）。

## breakpointINPUT.tsv
这是供下游断点搜索脚本使用的“区间输入”表，列顺序固定（写入时 header=False, index=False）：
1. IndividualID：样本标识，同上。
2. Cluster_No：核聚簇大小，同上。
3. disFile：用于聚类的 discordant SAM 文件路径（${OUTPUT}.mt.disc.sam）。
4. splitFile：对应的 split SAM 文件路径（${OUTPUT}.mt.split.sam）。
5. wgsBAM：原始全基因组 BAM 文件路径。
6. chr：核染色体名称。
7. start：起始坐标（bp），同上。
8. end：结束坐标（bp），同上。
这个表会被断点脚本读取，用来对每个核区间在 MT 和 genome 上精准定位潜在的 NUMT 断点。

## 总计
先把那些乱七八糟、染色体名字对不上号的读段过滤掉。
挑出那些一端落在核上、另一端落在线粒体上的读段。
把这些嫌疑按照核基因组的位置——用 500 bp 的“邻居距离”——分成一簇簇。
聚好簇后，再看看每簇里连到线粒体的读段都跑到线粒体哪个位置，同样按“500 bp 规则”再分子簇。
给每一簇贴上一个独一无二的 ID，像 “chr3_10000_10500_MTboth_500_800” 这样，标明核区间、线粒体区间。
最后把这些详细的「侦查大表」、按簇统计的「嫌疑概要表」和给下一步“精确破案”用的「断点输入表」统统打包输出。

# 0_2_sam2psl.py
该代码功能单一，仅仅将sam文件转为psl文件。

# 0_3_找断点.py
改代码输出精确定位的断点位置。
## 代码筛选标准总结
1. 初始过滤条件
- 比对长度过滤: matchLEN < 150（比对到参考序列的长度小于150 bp）
- 错配数过滤: misMatch <= 5（允许最多5个错配）
- 断点位置过滤: Tend >= 147 或 Tend <= 3（比对末端在参考序列的末端附近）
2. 核序列（nu）分类规则
- 正向链 (+):
- Qend >= 145（readLEN - mismatchLEN = 150 - 5） → nu_Tstart_Bright（右端断点）
- Qstart <= 5 → nu_Tend_Bleft（左端断点）
- 负向链 (-): 直接标记为 nu_NegStrand（负链断点）
3. 线粒体序列（mt）分类规则
- 正向链 (+):
- Qend >= 145 → mt_Tstart（线粒体起始断点）
- Qstart <= 5 → mt_Tend（线粒体末端断点）
- 负向链 (-):
- Qend >= 145 → mt_Tend
- Qstart <= 5 → mt_Tstart
4. 核序列范围
- 仅保留目标染色体（CHR参数）的比对，无位置限制（原代码有弹性窗口，但最终放宽为全染色体）。

| 列名            | 含义                                                          |
| -------------- | --------------------------------------------------------------|
| **pointGroup** | 断点类型分类（如 `nu_Tend_Bleft`, `mt_Tstart` 等）              |
| **Group**      | 简化分类（`nu` 或 `mt` 前缀 + 方向，如 `nu`, `mtLeft`, `mtRight`）|
| **chr**        | 染色体名称（核为输入的 `CHR` 参数，线粒体为 `MT`）                |
| **Tstart**     | 比对在参考序列上的起始位置（仅右端断点有效，左端断点用 `-1` 填充）  |
| **Tend**       | 比对在参考序列上的终止位置（仅左端断点有效，右端断点用 `-1` 填充）  |
| **strand**     | 比对方向（`+` 或 `-`）                                          |
| **readsCount** | 支持该断点的 reads 数量                                         |
| **sampleID**   | 输入的样本标识符（`SAMPLEID` 参数）                              |
