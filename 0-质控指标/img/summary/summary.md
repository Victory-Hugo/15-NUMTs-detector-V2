# QC三表归纳总结

- 总表ID集合大小: 8546
- 规则: 三个输入TSV均先与总表ID取交集，再进入归纳与可视化。

## Ancestry.tsv
- rows_before_intersection: 34184
- rows_after_intersection: 34184
- unique_id_after_intersection: 8546

## CopyNumber.tsv
- rows_before_intersection: 8754
- rows_after_intersection: 8754
- unique_id_after_intersection: 8546
- numeric_invalid_rows: 19

## VerifyBam.tsv
- rows_before_intersection: 8542
- rows_after_intersection: 8542
- unique_id_after_intersection: 8542
- fail_with_freemix_le_0.02: 0

详细统计见 summary/*.csv；每张图的绘图配套数据见 data/*.csv。
