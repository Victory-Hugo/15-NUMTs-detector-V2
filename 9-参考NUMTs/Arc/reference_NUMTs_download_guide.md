# 人类参考 NUMTs 下载清单

## 说明
“参考 NUMTs”通常指已经嵌入某个**参考核基因组版本**中的固定或参考位点。
它不是唯一的一份清单，数量会随：
1. 参考基因组版本（hg18、hg19、hg38、T2T-CHM13 或文中写作 T2T-CHR13）
2. 检索算法与阈值
3. 是否把相邻 HSP 合并为一个 assembled NUMT
而变化。

因此，真正“全面”的做法不是只给一份表，而是把主流参考目录都列出来。

---

## 1) RHNumtS（2008，hg18，190 个）
- 论文：Lascaro et al. 2008
- 口径：初版 Reference Human NumtS compilation
- 数量：190
- 下载：Additional file 1

```text
https://pmc.ncbi.nlm.nih.gov/articles/PMC2447851/bin/1471-2164-9-267-S1.XLS
```

---

## 2) RHNumtS.2（2011，hg18，585 assembled NUMTs；766 HSP）
- 论文：Simone et al. 2011
- 口径：扩展并验证后的 RHNumtS.2
- 数量：585 assembled NUMTs，对应 766 HSP
- 下载：Additional file 1

```text
https://pmc.ncbi.nlm.nih.gov/articles/PMC3228558/bin/1471-2164-12-517-S1.XLS
```

---

## 3) RHNumtS.2 lifted to hg19（2012，hg19）
- 论文：Calabrese et al. 2012
- 口径：把 RHNumtS.2 从 hg18 remap 到 hg19
- 下载：Additional file 1

```text
https://static-content.springer.com/esm/art%3A10.1186%2F1471-2105-13-S4-S15/MediaObjects/12859_2012_5112_MOESM1_ESM.xls
```

---

## 4) Tao et al. 2023（hg19、hg38、T2T-CHM13/文中写作 T2T-CHR13）
- 论文：Tao et al. 2023, Genes
- 口径：使用 human pan-mitogenome (HPMT) 重新识别 fixed NUMTs
- 合并 complete mtDNA 与 control-region 结果后：
  - hg19：907
  - hg38：908
  - T2T-CHM13/文中写作 T2T-CHR13：958
- 下载：Supplementary File 1 ZIP
- 其中包含：
  - File S2：NUMTs in three reference genomes
  - File S3：Summary of mtDNA-like short segments
  - File S4：Mito-blacklist identified using ATAC-seq
  - File S5：NUMTs-Blacklist

```text
https://www.mdpi.com/article/10.3390/genes14112092/s1
```

---

## 5) UCSC 浏览器轨道
- RHNumtS.2 作者同时把参考 NUMTs 轨道放进了 UCSC Genome Browser
- 适合浏览，不一定比补充表更适合批量下载
- 如果只是做分析，优先下载上面的 xls 或 zip

---

## 6) 目前更“新”的参考集
- Fu et al. 2026 预印本报告在 T2T-CHM13 基础上识别到 918 reference NUMTs，并在更大图谱中给出 fixed 与 polymorphic NUMTs。
- 但在本次检索里，我没有拿到一个可直接核实的公开原始补充表下载地址，所以没有把它列入“可直接下载清单”。

---

## 最实用的建议
如果你是做人类数据分析，优先下载这 4 套：
1. RHNumtS（hg18）
2. RHNumtS.2（hg18）
3. RHNumtS.2 hg19 remap
4. Tao 2023 的 hg19、hg38、T2T-CHM13 全套

这样基本覆盖：
- 早期经典参考集
- UCSC 体系参考集
- 现代完整参考基因组体系
