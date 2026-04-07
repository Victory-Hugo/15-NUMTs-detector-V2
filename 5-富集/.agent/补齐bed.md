# 补齐 Figure 7D 所需 BED 区域清单

## 目标

为了让当前 `5-富集` 的核基因组区域富集分析尽量对齐 `2025.06.21.660892v1.full(1).pdf` 中 Figure 7D，需要补齐文章中用于 nuclear breakpoint enrichment 的全部基因组功能区域 BED。

文章方法段说明：用于 NUMT 插入核基因组序列特征分析的 genome regions 来自 UCSC 和 GENCODE v43，包括 centromere、genomic duplications、simple repeats、functional elements、CpG islands、satellites、LINEs、SINEs、gene 等。Figure 7D 图面进一步给出了完整类别标签。

当前已有 BED：

| 当前文件 | 对应文章类别 | 状态 |
|---|---|---|
| `data/regions/01-Promoter.bed` | 文章 Figure 7D 未单独列 promoter | 作为本项目扩展分析保留并参与运行 |
| `data/regions/02-Protein_coding_region.bed` | 近似 `CDS`，但不能替代 `Gene`、`Exon`、`UTR` | 作为本项目原有自定义区域保留并参与运行 |
| `data/regions/03-Segmental_duplications_custom.bed` | `Genomic_superdups` / genomic duplications | 作为本项目原有自定义区域保留；文章标准版本为 `22-Genomic_superdups.bed` |
| `data/regions/04-STR_without_segmental_duplications.bed` | 近似 `Microsat` 或 `Simple_Repeats` | 作为本项目原有去 superdups 后 STR 区域保留；文章标准版本为 `06-Microsat.bed` 和 `21-Simple_Repeats.bed` |

## 必须补齐的核基因组 BED

以下清单按 Figure 7D 图面 x 轴类别整理。若要复现文章核基因组富集分析，以下类别不要遗漏。

| 建议文件名 | 文章类别 | 推荐来源 | 备注 |
|---|---|---|---|
| `5-CpG_islands.bed` | `CpG_islands` | UCSC hg38 `cpgIslandExt` 或等价 CpG island track | 方法段明确提到 CpG islands |
| `6-Microsat.bed` | `Microsat` | UCSC/RepeatMasker 或 simple repeat 资源中的 microsatellite 注释 | Figure 7D 单独成列，不能被当前 STR 文件完全替代 |
| `7-Start_codon.bed` | `Start_codon` | GENCODE v43 GTF `start_codon` feature | Figure 7D 单独成列 |
| `8-snRNA.bed` | `snRNA` | GENCODE v43，筛选 `gene_type` 或 `gene_biotype` 为 `snRNA` | Figure 7D 单独成列 |
| `9-FuncElems.bed` | `FuncElems` | UCSC functional elements / regulatory elements track | 方法段明确提到 functional elements；需记录具体 UCSC 表名 |
| `10-Retroposon.bed` | `Retroposon` | UCSC RepeatMasker `rmsk` 中 retroposon 相关 `repClass`/`repFamily` | Figure 7D 单独成列 |
| `11-CDS.bed` | `CDS` | GENCODE v43 GTF `CDS` feature | 不能只用当前蛋白编码区域文件替代，建议重新由 GENCODE v43 生成 |
| `12-Stop_codon.bed` | `Stop_codon` | GENCODE v43 GTF `stop_codon` feature | Figure 7D 单独成列 |
| `13-SINE.bed` | `SINE` | UCSC RepeatMasker `rmsk`，筛选 `repClass == SINE` | 方法段和讨论均提到 |
| `14-rmsk-DNA.bed` | `rmsk-DNA` | UCSC RepeatMasker `rmsk`，筛选 `repClass == DNA` | 图注解释为 repetitive DNA |
| `15-Centromeres.bed` | `Centromeres` | UCSC hg38 centromere 注释，或 `gap` table 中 `type=centromere` | 方法段明确提到 centromere |
| `16-Satellite.bed` | `Satellite` | UCSC RepeatMasker `rmsk`，筛选 `repClass == Satellite`，或 UCSC satellite track | 方法段明确提到 satellites |
| `17-LTR.bed` | `LTR` | UCSC RepeatMasker `rmsk`，筛选 `repClass == LTR` | 图注解释为 retrotransposon |
| `18-Exon.bed` | `Exon` | GENCODE v43 GTF `exon` feature | Figure 7D 单独成列 |
| `19-UTR.bed` | `UTR` | GENCODE v43 GTF `UTR` feature，或合并 `five_prime_UTR`/`three_prime_UTR` | 图注解释为 untranslated region |
| `20-LINE.bed` | `LINE` | UCSC RepeatMasker `rmsk`，筛选 `repClass == LINE` | 文章核心结论之一：rare NUMTs 富集于 LINE |
| `21-Simple_Repeats.bed` | `Simple_Repeats` | UCSC hg38 `simpleRepeat` track | 文章正文说 common NUMTs 倾向 simple repeats |
| `22-Genomic_superdups.bed` | `Genomic_superdups` | UCSC hg38 `genomicSuperDups` | 当前 `3-大规模片段重复区域.bed` 基本对应此项 |
| `23-snoRNA_miRNA.bed` | `snoRNA_miRNA` | GENCODE v43，筛选 `snoRNA` 与 `miRNA` gene biotype 后合并 | 文章正文说 common NUMTs 倾向 miRNA 和 snoRNA |
| `24-Gene.bed` | `Gene` | GENCODE v43 GTF `gene` feature | 文章说 ultra-rare / overall breakpoints mainly concentrated in gene regions |

## 正文提到但 Figure 7D 未清楚单独成列的补充项

为避免遗漏文章文字中出现的核基因组区域术语，建议额外生成以下 BED，作为扩展或敏感性分析使用。报告时需注明这些不是 Figure 7D 轴标签的核心类别，除非后续从原图或作者代码确认其确实被纳入。

| 建议文件名 | 文章提到的类别 | 推荐来源 | 备注 |
|---|---|---|---|
| `25-Intron.bed` | `introns` | GENCODE v43：`gene` 区间减去 `exon` 区间 | 讨论中提到 100,000 Genomes Project 中 ultra-rare NUMTs nuclear breakpoints 在 introns 中丰富 |
| `26-Regulatory_elements.bed` | `regulatory elements` | UCSC cCRE/ENCODE regulatory elements 或与 `FuncElems` 同源 track | 讨论中提到 regulatory elements；如果 `FuncElems` 已等价覆盖，可只保留 `FuncElems` 并在说明中注明 |

## 不应混入核基因组 BED 的 mtDNA 区域

文章 Figure 7A/B 还分析 mtDNA breakpoints，但这不是 `data/regions/*.bed` 的核基因组区域富集。若要完整复现 Figure 7A/B，应另建 mtDNA 专用流程，不应和 nuclear BED 混在一起。

需要单独考虑的 mtDNA 注释包括：

| mtDNA 类别 | 文章依据 |
|---|---|
| `CDS` / protein-coding regions | common NUMTs mtDNA breakpoints enriched in coding regions in the middle of mtDNA |
| `D-loop` | rare NUMTs mtDNA breakpoints tended to occur on both sides including D-loop |
| `rRNA`，尤其 `RNR1` | RNR1 region showed higher normalized NUMT mtDNA breakpoint count |
| `tRNA` | mtDNA functional locations 的常规组成 |
| `intergenic` / non-genic regions | mtDNA functional locations 的常规组成 |

## 当前完成状态

已按本文清单生成 22 个文章相关核基因组 BED，并与本项目原有 4 个自定义 BED 统一整理到 `data/regions/`。BED 准备步骤作为长期资源构建流程保留在：

```text
script/prepare_article_regions_bed.sh
```

主富集流程 `pipe/run_enrichment_analysis.sh` 不会自动重新生成这些 BED；后续只有在更新 UCSC/GENCODE 原始注释或需要重建资源时，才需要单独运行该脚本。

当前 `data/regions/` 应包含：

```text
01-Promoter.bed
02-Protein_coding_region.bed
03-Segmental_duplications_custom.bed
04-STR_without_segmental_duplications.bed
05-CpG_islands.bed
06-Microsat.bed
07-Start_codon.bed
08-snRNA.bed
09-FuncElems.bed
10-Retroposon.bed
11-CDS.bed
12-Stop_codon.bed
13-SINE.bed
14-rmsk-DNA.bed
15-Centromeres.bed
16-Satellite.bed
17-LTR.bed
18-Exon.bed
19-UTR.bed
20-LINE.bed
21-Simple_Repeats.bed
22-Genomic_superdups.bed
23-snoRNA_miRNA.bed
24-Gene.bed
25-Intron.bed
26-Regulatory_elements.bed
```

1. 原始 UCSC 表和 GENCODE v43 GTF 已保留在 `data/raw/`；最终参与运行的 BED 统一保留在 `data/regions/`。
2. `script/prepare_article_regions_bed.sh` 已建立，可统一生成上面所有 BED，避免手工下载后不可复现。
3. 所有 BED 已做标准化处理：保留 BED3+ 格式，排序，合并重叠区间，并过滤到标准染色体。
4. GENCODE v43 生成的 `CDS`、`Exon`、`UTR`、`Start_codon`、`Stop_codon`、`Gene`、`snRNA`、`snoRNA_miRNA`、`Intron` 已由 `python/prepare_gencode_article_regions.py` 统一生成。
5. UCSC RepeatMasker 生成的 `LINE`、`SINE`、`LTR`、`rmsk-DNA`、`Satellite`、`Retroposon` 已按 `repClass`/相关类别筛选生成。
6. 当前 `pipe/run_enrichment_analysis.sh` 已切换为扫描 `data/regions/*.bed`，并对所有类别做 `all/common/low-frequency/rare/ultra-rare` 分层的 nuclear breakpoint flank enrichment。

## 下一步工作建议

## 最小优先级

如果先做一版最接近文章主结论的补充，优先生成：

```text
20-LINE.bed
24-Gene.bed
22-Genomic_superdups.bed
21-Simple_Repeats.bed
23-snoRNA_miRNA.bed
13-SINE.bed
17-LTR.bed
18-Exon.bed
19-UTR.bed
11-CDS.bed
```

但最终复现 Figure 7D 时，仍应按“必须补齐的核基因组 BED”全清单执行。
