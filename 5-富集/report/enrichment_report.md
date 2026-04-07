# NUMT nuclear breakpoint flank enrichment report

## 中文摘要

本分析按照已有研究的富集思路，使用 confident nuclear breakpoints 上下游各 100 bp 的 flanks 作为分析区间，并按 `all`、`common`、`low-frequency`、`rare` 和 `ultra-rare` 分层，对目标基因组区域进行全量 permutation 富集分析。

### Breakpoint flank 分层数量

| frequency_class   |   breakpoint_count |   cluster_count |
|:------------------|-------------------:|----------------:|
| all               |               2337 |            1897 |
| common            |                147 |             123 |
| low-frequency     |                151 |             126 |
| rare              |                344 |             259 |
| ultra-rare        |               1695 |            1389 |

### 主要结果

- `03-Segmental_duplications_custom.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 11.60%，经验 P=<0.0001）。
- `03-Segmental_duplications_custom.bed` 中 `common` NUMT breakpoint flanks 显著富集（观测重叠比例 19.05%，经验 P=<0.0001）。
- `03-Segmental_duplications_custom.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 11.34%，经验 P=<0.0001）。
- `03-Segmental_duplications_custom.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 11.15%，经验 P=<0.0001）。
- `04-STR_without_segmental_duplications.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 20.71%，经验 P=<0.0001）。
- `04-STR_without_segmental_duplications.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 19.94%，经验 P=<0.0001）。
- `13-SINE.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 29.35%，经验 P=<0.0001）。
- `13-SINE.bed` 中 `common` NUMT breakpoint flanks 显著富集（观测重叠比例 48.98%，经验 P=<0.0001）。
- `13-SINE.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 50.33%，经验 P=<0.0001）。
- `13-SINE.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 36.63%，经验 P=<0.0001）。
- `21-Simple_Repeats.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 27.90%，经验 P=<0.0001）。
- `21-Simple_Repeats.bed` 中 `common` NUMT breakpoint flanks 显著富集（观测重叠比例 69.39%，经验 P=<0.0001）。
- `21-Simple_Repeats.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 64.24%，经验 P=<0.0001）。
- `21-Simple_Repeats.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 47.09%，经验 P=<0.0001）。
- `21-Simple_Repeats.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 17.17%，经验 P=<0.0001）。
- `22-Genomic_superdups.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 11.60%，经验 P=<0.0001）。
- `22-Genomic_superdups.bed` 中 `common` NUMT breakpoint flanks 显著富集（观测重叠比例 19.05%，经验 P=<0.0001）。
- `22-Genomic_superdups.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 11.34%，经验 P=<0.0001）。
- `22-Genomic_superdups.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 11.15%，经验 P=<0.0001）。
- `04-STR_without_segmental_duplications.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 27.15%，经验 P=0.0001）。
- `10-Retroposon.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 2.65%，经验 P=0.0002）。
- `15-Centromeres.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著减少（观测重叠比例 0.88%，经验 P=0.0005）。
- `10-Retroposon.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 0.51%，经验 P=0.0009）。
- `06-Microsat.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 0.73%，经验 P=0.001）。
- `04-STR_without_segmental_duplications.bed` 中 `common` NUMT breakpoint flanks 显著富集（观测重叠比例 25.85%，经验 P=0.0011）。
- `18-Exon.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 9.09%，经验 P=0.0014）。
- `06-Microsat.bed` 中 `common` NUMT breakpoint flanks 显著富集（观测重叠比例 2.72%，经验 P=0.0019）。
- `15-Centromeres.bed` 中 `all` NUMT breakpoint flanks 显著减少（观测重叠比例 1.20%，经验 P=0.0046）。
- `18-Exon.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 8.47%，经验 P=0.0053）。
- `13-SINE.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 24.31%，经验 P=0.0068）。
- `16-Satellite.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 4.65%，经验 P=0.0124）。
- `03-Segmental_duplications_custom.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 9.93%，经验 P=0.0183）。
- `22-Genomic_superdups.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 9.93%，经验 P=0.0189）。
- `10-Retroposon.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 0.87%，经验 P=0.0208）。
- `24-Gene.bed` 中 `all` NUMT breakpoint flanks 显著富集（观测重叠比例 59.82%，经验 P=0.0222）。
- `06-Microsat.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 1.16%，经验 P=0.0238）。
- `02-Protein_coding_region.bed` 中 `common` NUMT breakpoint flanks 显著减少（观测重叠比例 0.00%，经验 P=0.025）。
- `11-CDS.bed` 中 `common` NUMT breakpoint flanks 显著减少（观测重叠比例 0.00%，经验 P=0.0271）。
- `20-LINE.bed` 中 `all` NUMT breakpoint flanks 显著减少（观测重叠比例 26.44%，经验 P=0.0279）。
- `08-snRNA.bed` 中 `low-frequency` NUMT breakpoint flanks 显著富集（观测重叠比例 0.66%，经验 P=0.0287）。
- `02-Protein_coding_region.bed` 中 `rare` NUMT breakpoint flanks 显著减少（观测重叠比例 0.87%，经验 P=0.0292）。
- `17-LTR.bed` 中 `rare` NUMT breakpoint flanks 显著减少（观测重叠比例 8.72%，经验 P=0.0308）。
- `11-CDS.bed` 中 `rare` NUMT breakpoint flanks 显著减少（观测重叠比例 0.87%，经验 P=0.0323）。
- `20-LINE.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著减少（观测重叠比例 26.31%，经验 P=0.0387）。
- `24-Gene.bed` 中 `ultra-rare` NUMT breakpoint flanks 显著富集（观测重叠比例 59.76%，经验 P=0.0474）。
- `04-STR_without_segmental_duplications.bed` 中 `rare` NUMT breakpoint flanks 显著富集（观测重叠比例 19.48%，经验 P=0.0477）。

简短结论：NUMT 核基因组断点 flanks 在不同功能区域和频率层级中呈现可检测的区域差异；显著结果应优先结合重复序列背景、可比对性和公共 NUMT 注释进一步解释。

## English Summary

Following the nuclear-genome enrichment strategy in Figure 7D of the paper, this analysis used 100 bp flanks on both sides of confident nuclear NUMT breakpoints and performed full permutation-based enrichment testing across target genomic regions, stratified by `all`, `common`, `low-frequency`, `rare`, and `ultra-rare` NUMTs.

### Breakpoint flank strata counts

| frequency_class   |   breakpoint_count |   cluster_count |
|:------------------|-------------------:|----------------:|
| all               |               2337 |            1897 |
| common            |                147 |             123 |
| low-frequency     |                151 |             126 |
| rare              |                344 |             259 |
| ultra-rare        |               1695 |            1389 |

### Main results

- In 03-Segmental_duplications_custom.bed, all NUMT breakpoint flanks were significantly enriched (observed overlap 11.60%, empirical P=<0.0001).
- In 03-Segmental_duplications_custom.bed, common NUMT breakpoint flanks were significantly enriched (observed overlap 19.05%, empirical P=<0.0001).
- In 03-Segmental_duplications_custom.bed, rare NUMT breakpoint flanks were significantly enriched (observed overlap 11.34%, empirical P=<0.0001).
- In 03-Segmental_duplications_custom.bed, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 11.15%, empirical P=<0.0001).
- In 04-STR_without_segmental_duplications.bed, all NUMT breakpoint flanks were significantly enriched (observed overlap 20.71%, empirical P=<0.0001).
- In 04-STR_without_segmental_duplications.bed, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 19.94%, empirical P=<0.0001).
- In SINEs, all NUMT breakpoint flanks were significantly enriched (observed overlap 29.35%, empirical P=<0.0001).
- In SINEs, common NUMT breakpoint flanks were significantly enriched (observed overlap 48.98%, empirical P=<0.0001).
- In SINEs, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 50.33%, empirical P=<0.0001).
- In SINEs, rare NUMT breakpoint flanks were significantly enriched (observed overlap 36.63%, empirical P=<0.0001).
- In simple repeats, all NUMT breakpoint flanks were significantly enriched (observed overlap 27.90%, empirical P=<0.0001).
- In simple repeats, common NUMT breakpoint flanks were significantly enriched (observed overlap 69.39%, empirical P=<0.0001).
- In simple repeats, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 64.24%, empirical P=<0.0001).
- In simple repeats, rare NUMT breakpoint flanks were significantly enriched (observed overlap 47.09%, empirical P=<0.0001).
- In simple repeats, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 17.17%, empirical P=<0.0001).
- In genomic superduplications, all NUMT breakpoint flanks were significantly enriched (observed overlap 11.60%, empirical P=<0.0001).
- In genomic superduplications, common NUMT breakpoint flanks were significantly enriched (observed overlap 19.05%, empirical P=<0.0001).
- In genomic superduplications, rare NUMT breakpoint flanks were significantly enriched (observed overlap 11.34%, empirical P=<0.0001).
- In genomic superduplications, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 11.15%, empirical P=<0.0001).
- In 04-STR_without_segmental_duplications.bed, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 27.15%, empirical P=0.0001).
- In retroposons, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 2.65%, empirical P=0.0002).
- In centromeres, ultra-rare NUMT breakpoint flanks were significantly depleted (observed overlap 0.88%, empirical P=0.0005).
- In retroposons, all NUMT breakpoint flanks were significantly enriched (observed overlap 0.51%, empirical P=0.0009).
- In microsatellites, all NUMT breakpoint flanks were significantly enriched (observed overlap 0.73%, empirical P=0.001).
- In 04-STR_without_segmental_duplications.bed, common NUMT breakpoint flanks were significantly enriched (observed overlap 25.85%, empirical P=0.0011).
- In exons, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 9.09%, empirical P=0.0014).
- In microsatellites, common NUMT breakpoint flanks were significantly enriched (observed overlap 2.72%, empirical P=0.0019).
- In centromeres, all NUMT breakpoint flanks were significantly depleted (observed overlap 1.20%, empirical P=0.0046).
- In exons, all NUMT breakpoint flanks were significantly enriched (observed overlap 8.47%, empirical P=0.0053).
- In SINEs, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 24.31%, empirical P=0.0068).
- In satellites, rare NUMT breakpoint flanks were significantly enriched (observed overlap 4.65%, empirical P=0.0124).
- In 03-Segmental_duplications_custom.bed, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 9.93%, empirical P=0.0183).
- In genomic superduplications, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 9.93%, empirical P=0.0189).
- In retroposons, rare NUMT breakpoint flanks were significantly enriched (observed overlap 0.87%, empirical P=0.0208).
- In gene regions, all NUMT breakpoint flanks were significantly enriched (observed overlap 59.82%, empirical P=0.0222).
- In microsatellites, rare NUMT breakpoint flanks were significantly enriched (observed overlap 1.16%, empirical P=0.0238).
- In 02-Protein_coding_region.bed, common NUMT breakpoint flanks were significantly depleted (observed overlap 0.00%, empirical P=0.025).
- In CDS regions, common NUMT breakpoint flanks were significantly depleted (observed overlap 0.00%, empirical P=0.0271).
- In LINEs, all NUMT breakpoint flanks were significantly depleted (observed overlap 26.44%, empirical P=0.0279).
- In snRNA regions, low-frequency NUMT breakpoint flanks were significantly enriched (observed overlap 0.66%, empirical P=0.0287).
- In 02-Protein_coding_region.bed, rare NUMT breakpoint flanks were significantly depleted (observed overlap 0.87%, empirical P=0.0292).
- In LTR elements, rare NUMT breakpoint flanks were significantly depleted (observed overlap 8.72%, empirical P=0.0308).
- In CDS regions, rare NUMT breakpoint flanks were significantly depleted (observed overlap 0.87%, empirical P=0.0323).
- In LINEs, ultra-rare NUMT breakpoint flanks were significantly depleted (observed overlap 26.31%, empirical P=0.0387).
- In gene regions, ultra-rare NUMT breakpoint flanks were significantly enriched (observed overlap 59.76%, empirical P=0.0474).
- In 04-STR_without_segmental_duplications.bed, rare NUMT breakpoint flanks were significantly enriched (observed overlap 19.48%, empirical P=0.0477).

Concise conclusion: NUMT nuclear breakpoint flanks show region- and frequency-class-specific distributional patterns. Significant enrichment or depletion signals should be interpreted together with repeat context, mappability, and known NUMT annotations before assigning biological mechanism.

## Full result table

| region_name                               | frequency_class   |   breakpoints_total |   breakpoints_target |   observed_percentage |   p_value_less |   p_value_greater |
|:------------------------------------------|:------------------|--------------------:|---------------------:|----------------------:|---------------:|------------------:|
| 01-Promoter.bed                           | all               |                2337 |                   39 |           0.0166881   |         0.583  |            0.4795 |
| 01-Promoter.bed                           | common            |                 147 |                    1 |           0.00680272  |         0.3071 |            0.9115 |
| 01-Promoter.bed                           | low-frequency     |                 151 |                    1 |           0.00662252  |         0.2869 |            0.9162 |
| 01-Promoter.bed                           | rare              |                 344 |                    9 |           0.0261628   |         0.9376 |            0.116  |
| 01-Promoter.bed                           | ultra-rare        |                1695 |                   28 |           0.0165192   |         0.5599 |            0.5112 |
| 02-Protein_coding_region.bed              | all               |                2337 |                   46 |           0.0196834   |         0.0664 |            0.9494 |
| 02-Protein_coding_region.bed              | common            |                 147 |                    0 |           0           |         0.025  |            1      |
| 02-Protein_coding_region.bed              | low-frequency     |                 151 |                    2 |           0.013245    |         0.2814 |            0.8895 |
| 02-Protein_coding_region.bed              | rare              |                 344 |                    3 |           0.00872093  |         0.0292 |            0.9912 |
| 02-Protein_coding_region.bed              | ultra-rare        |                1695 |                   41 |           0.0241888   |         0.4992 |            0.5656 |
| 03-Segmental_duplications_custom.bed      | all               |                2337 |                  271 |           0.115961    |         1      |            0      |
| 03-Segmental_duplications_custom.bed      | common            |                 147 |                   28 |           0.190476    |         1      |            0      |
| 03-Segmental_duplications_custom.bed      | low-frequency     |                 151 |                   15 |           0.0993377   |         0.992  |            0.0183 |
| 03-Segmental_duplications_custom.bed      | rare              |                 344 |                   39 |           0.113372    |         1      |            0      |
| 03-Segmental_duplications_custom.bed      | ultra-rare        |                1695 |                  189 |           0.111504    |         1      |            0      |
| 04-STR_without_segmental_duplications.bed | all               |                2337 |                  484 |           0.207103    |         1      |            0      |
| 04-STR_without_segmental_duplications.bed | common            |                 147 |                   38 |           0.258503    |         0.9995 |            0.0011 |
| 04-STR_without_segmental_duplications.bed | low-frequency     |                 151 |                   41 |           0.271523    |         0.9999 |            0.0001 |
| 04-STR_without_segmental_duplications.bed | rare              |                 344 |                   67 |           0.194767    |         0.9648 |            0.0477 |
| 04-STR_without_segmental_duplications.bed | ultra-rare        |                1695 |                  338 |           0.19941     |         1      |            0      |
| 05-CpG_islands.bed                        | all               |                2337 |                   17 |           0.00727428  |         0.2447 |            0.824  |
| 05-CpG_islands.bed                        | common            |                 147 |                    0 |           0           |         0.274  |            1      |
| 05-CpG_islands.bed                        | low-frequency     |                 151 |                    2 |           0.013245    |         0.8533 |            0.3835 |
| 05-CpG_islands.bed                        | rare              |                 344 |                    5 |           0.0145349   |         0.9166 |            0.194  |
| 05-CpG_islands.bed                        | ultra-rare        |                1695 |                   10 |           0.00589971  |         0.1113 |            0.9355 |
| 06-Microsat.bed                           | all               |                2337 |                   17 |           0.00727428  |         0.9996 |            0.001  |
| 06-Microsat.bed                           | common            |                 147 |                    4 |           0.0272109   |         0.9997 |            0.0019 |
| 06-Microsat.bed                           | low-frequency     |                 151 |                    0 |           0           |         0.6241 |            1      |
| 06-Microsat.bed                           | rare              |                 344 |                    4 |           0.0116279   |         0.995  |            0.0238 |
| 06-Microsat.bed                           | ultra-rare        |                1695 |                    9 |           0.00530973  |         0.9552 |            0.0921 |
| 07-Start_codon.bed                        | all               |                2337 |                    3 |           0.0012837   |         0.3001 |            0.8531 |
| 07-Start_codon.bed                        | common            |                 147 |                    0 |           0           |         0.742  |            1      |
| 07-Start_codon.bed                        | low-frequency     |                 151 |                    0 |           0           |         0.7411 |            1      |
| 07-Start_codon.bed                        | rare              |                 344 |                    1 |           0.00290698  |         0.8501 |            0.4967 |
| 07-Start_codon.bed                        | ultra-rare        |                1695 |                    2 |           0.00117994  |         0.3343 |            0.8614 |
| 08-snRNA.bed                              | all               |                2337 |                    1 |           0.000427899 |         0.9291 |            0.3558 |
| 08-snRNA.bed                              | common            |                 147 |                    0 |           0           |         0.9726 |            1      |
| 08-snRNA.bed                              | low-frequency     |                 151 |                    1 |           0.00662252  |         0.9995 |            0.0287 |
| 08-snRNA.bed                              | rare              |                 344 |                    0 |           0           |         0.9341 |            1      |
| 08-snRNA.bed                              | ultra-rare        |                1695 |                    0 |           0           |         0.7208 |            1      |
| 09-FuncElems.bed                          | all               |                2337 |                    4 |           0.0017116   |         0.7439 |            0.4463 |
| 09-FuncElems.bed                          | common            |                 147 |                    0 |           0           |         0.7984 |            1      |
| 09-FuncElems.bed                          | low-frequency     |                 151 |                    0 |           0           |         0.7961 |            1      |
| 09-FuncElems.bed                          | rare              |                 344 |                    0 |           0           |         0.6077 |            1      |
| 09-FuncElems.bed                          | ultra-rare        |                1695 |                    4 |           0.00235988  |         0.9002 |            0.2358 |
| 10-Retroposon.bed                         | all               |                2337 |                   12 |           0.00513479  |         1      |            0.0009 |
| 10-Retroposon.bed                         | common            |                 147 |                    1 |           0.00680272  |         0.9744 |            0.2169 |
| 10-Retroposon.bed                         | low-frequency     |                 151 |                    4 |           0.0264901   |         1      |            0.0002 |
| 10-Retroposon.bed                         | rare              |                 344 |                    3 |           0.00872093  |         0.9966 |            0.0208 |
| 10-Retroposon.bed                         | ultra-rare        |                1695 |                    4 |           0.00235988  |         0.8456 |            0.3152 |
| 11-CDS.bed                                | all               |                2337 |                   46 |           0.0196834   |         0.0699 |            0.9463 |
| 11-CDS.bed                                | common            |                 147 |                    0 |           0           |         0.0271 |            1      |
| 11-CDS.bed                                | low-frequency     |                 151 |                    2 |           0.013245    |         0.2834 |            0.8854 |
| 11-CDS.bed                                | rare              |                 344 |                    3 |           0.00872093  |         0.0323 |            0.9906 |
| 11-CDS.bed                                | ultra-rare        |                1695 |                   41 |           0.0241888   |         0.5017 |            0.5603 |
| 12-Stop_codon.bed                         | all               |                2337 |                    5 |           0.0021395   |         0.3055 |            0.822  |
| 12-Stop_codon.bed                         | common            |                 147 |                    0 |           0           |         0.6576 |            1      |
| 12-Stop_codon.bed                         | low-frequency     |                 151 |                    0 |           0           |         0.6389 |            1      |
| 12-Stop_codon.bed                         | rare              |                 344 |                    0 |           0           |         0.3616 |            1      |
| 12-Stop_codon.bed                         | ultra-rare        |                1695 |                    5 |           0.00294985  |         0.6119 |            0.5663 |
| 13-SINE.bed                               | all               |                2337 |                  686 |           0.293539    |         1      |            0      |
| 13-SINE.bed                               | common            |                 147 |                   72 |           0.489796    |         1      |            0      |
| 13-SINE.bed                               | low-frequency     |                 151 |                   76 |           0.503311    |         1      |            0      |
| 13-SINE.bed                               | rare              |                 344 |                  126 |           0.366279    |         1      |            0      |
| 13-SINE.bed                               | ultra-rare        |                1695 |                  412 |           0.243068    |         0.9945 |            0.0068 |
| 14-rmsk-DNA.bed                           | all               |                2337 |                  133 |           0.0569106   |         0.2368 |            0.7872 |
| 14-rmsk-DNA.bed                           | common            |                 147 |                   13 |           0.0884354   |         0.9366 |            0.1119 |
| 14-rmsk-DNA.bed                           | low-frequency     |                 151 |                    6 |           0.0397351   |         0.1823 |            0.9026 |
| 14-rmsk-DNA.bed                           | rare              |                 344 |                   17 |           0.0494186   |         0.2257 |            0.8355 |
| 14-rmsk-DNA.bed                           | ultra-rare        |                1695 |                   97 |           0.0572271   |         0.307  |            0.7295 |
| 15-Centromeres.bed                        | all               |                2337 |                   28 |           0.0119812   |         0.0046 |            0.9973 |
| 15-Centromeres.bed                        | common            |                 147 |                    3 |           0.0204082   |         0.6879 |            0.5343 |
| 15-Centromeres.bed                        | low-frequency     |                 151 |                    4 |           0.0264901   |         0.8345 |            0.3315 |
| 15-Centromeres.bed                        | rare              |                 344 |                    6 |           0.0174419   |         0.5068 |            0.6493 |
| 15-Centromeres.bed                        | ultra-rare        |                1695 |                   15 |           0.00884956  |         0.0005 |            0.9999 |
| 16-Satellite.bed                          | all               |                2337 |                   61 |           0.0261018   |         0.683  |            0.3634 |
| 16-Satellite.bed                          | common            |                 147 |                    3 |           0.0204082   |         0.4982 |            0.712  |
| 16-Satellite.bed                          | low-frequency     |                 151 |                    4 |           0.0264901   |         0.6816 |            0.5212 |
| 16-Satellite.bed                          | rare              |                 344 |                   16 |           0.0465116   |         0.9945 |            0.0124 |
| 16-Satellite.bed                          | ultra-rare        |                1695 |                   38 |           0.0224189   |         0.3103 |            0.7465 |
| 17-LTR.bed                                | all               |                2337 |                  283 |           0.121095    |         0.5063 |            0.518  |
| 17-LTR.bed                                | common            |                 147 |                   15 |           0.102041    |         0.2828 |            0.8035 |
| 17-LTR.bed                                | low-frequency     |                 151 |                   21 |           0.139073    |         0.7856 |            0.2903 |
| 17-LTR.bed                                | rare              |                 344 |                   30 |           0.0872093   |         0.0308 |            0.9781 |
| 17-LTR.bed                                | ultra-rare        |                1695 |                  217 |           0.128024    |         0.819  |            0.2004 |
| 18-Exon.bed                               | all               |                2337 |                  198 |           0.084724    |         0.996  |            0.0053 |
| 18-Exon.bed                               | common            |                 147 |                   13 |           0.0884354   |         0.8445 |            0.2401 |
| 18-Exon.bed                               | low-frequency     |                 151 |                   12 |           0.0794702   |         0.7334 |            0.3799 |
| 18-Exon.bed                               | rare              |                 344 |                   19 |           0.0552326   |         0.1556 |            0.8899 |
| 18-Exon.bed                               | ultra-rare        |                1695 |                  154 |           0.0908555   |         0.999  |            0.0014 |
| 19-UTR.bed                                | all               |                2337 |                   65 |           0.0278134   |         0.5398 |            0.5121 |
| 19-UTR.bed                                | common            |                 147 |                    5 |           0.0340136   |         0.7717 |            0.3935 |
| 19-UTR.bed                                | low-frequency     |                 151 |                    2 |           0.013245    |         0.2058 |            0.9273 |
| 19-UTR.bed                                | rare              |                 344 |                   10 |           0.0290698   |         0.6429 |            0.4882 |
| 19-UTR.bed                                | ultra-rare        |                1695 |                   48 |           0.0283186   |         0.5922 |            0.4644 |
| 20-LINE.bed                               | all               |                2337 |                  618 |           0.264442    |         0.0279 |            0.9755 |
| 20-LINE.bed                               | common            |                 147 |                   46 |           0.312925    |         0.8216 |            0.2264 |
| 20-LINE.bed                               | low-frequency     |                 151 |                   35 |           0.231788    |         0.0909 |            0.9346 |
| 20-LINE.bed                               | rare              |                 344 |                   91 |           0.264535    |         0.2495 |            0.7868 |
| 20-LINE.bed                               | ultra-rare        |                1695 |                  446 |           0.263127    |         0.0387 |            0.9661 |
| 21-Simple_Repeats.bed                     | all               |                2337 |                  652 |           0.27899     |         1      |            0      |
| 21-Simple_Repeats.bed                     | common            |                 147 |                  102 |           0.693878    |         1      |            0      |
| 21-Simple_Repeats.bed                     | low-frequency     |                 151 |                   97 |           0.642384    |         1      |            0      |
| 21-Simple_Repeats.bed                     | rare              |                 344 |                  162 |           0.47093     |         1      |            0      |
| 21-Simple_Repeats.bed                     | ultra-rare        |                1695 |                  291 |           0.171681    |         1      |            0      |
| 22-Genomic_superdups.bed                  | all               |                2337 |                  271 |           0.115961    |         1      |            0      |
| 22-Genomic_superdups.bed                  | common            |                 147 |                   28 |           0.190476    |         1      |            0      |
| 22-Genomic_superdups.bed                  | low-frequency     |                 151 |                   15 |           0.0993377   |         0.9917 |            0.0189 |
| 22-Genomic_superdups.bed                  | rare              |                 344 |                   39 |           0.113372    |         1      |            0      |
| 22-Genomic_superdups.bed                  | ultra-rare        |                1695 |                  189 |           0.111504    |         1      |            0      |
| 23-snoRNA_miRNA.bed                       | all               |                2337 |                    0 |           0           |         0.5438 |            1      |
| 23-snoRNA_miRNA.bed                       | common            |                 147 |                    0 |           0           |         0.9612 |            1      |
| 23-snoRNA_miRNA.bed                       | low-frequency     |                 151 |                    0 |           0           |         0.9594 |            1      |
| 23-snoRNA_miRNA.bed                       | rare              |                 344 |                    0 |           0           |         0.914  |            1      |
| 23-snoRNA_miRNA.bed                       | ultra-rare        |                1695 |                    0 |           0           |         0.6502 |            1      |
| 24-Gene.bed                               | all               |                2337 |                 1398 |           0.598203    |         0.9808 |            0.0222 |
| 24-Gene.bed                               | common            |                 147 |                   95 |           0.646259    |         0.9618 |            0.0546 |
| 24-Gene.bed                               | low-frequency     |                 151 |                   92 |           0.609272    |         0.8079 |            0.2429 |
| 24-Gene.bed                               | rare              |                 344 |                  198 |           0.575581    |         0.4885 |            0.5554 |
| 24-Gene.bed                               | ultra-rare        |                1695 |                 1013 |           0.59764     |         0.9568 |            0.0474 |
| 25-Intron.bed                             | all               |                2337 |                 1275 |           0.545571    |         0.7497 |            0.2649 |
| 25-Intron.bed                             | common            |                 147 |                   84 |           0.571429    |         0.8127 |            0.2362 |
| 25-Intron.bed                             | low-frequency     |                 151 |                   86 |           0.569536    |         0.7953 |            0.2495 |
| 25-Intron.bed                             | rare              |                 344 |                  186 |           0.540698    |         0.5484 |            0.4928 |
| 25-Intron.bed                             | ultra-rare        |                1695 |                  919 |           0.542183    |         0.6144 |            0.402  |
| 26-Regulatory_elements.bed                | all               |                2337 |                    4 |           0.0017116   |         0.7446 |            0.4375 |
| 26-Regulatory_elements.bed                | common            |                 147 |                    0 |           0           |         0.8012 |            1      |
| 26-Regulatory_elements.bed                | low-frequency     |                 151 |                    0 |           0           |         0.8054 |            1      |
| 26-Regulatory_elements.bed                | rare              |                 344 |                    0 |           0           |         0.6003 |            1      |
| 26-Regulatory_elements.bed                | ultra-rare        |                1695 |                    4 |           0.00235988  |         0.8966 |            0.2409 |
