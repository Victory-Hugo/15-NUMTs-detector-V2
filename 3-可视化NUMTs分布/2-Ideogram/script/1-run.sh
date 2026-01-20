#!/bin/bash

INPUT_TXT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/input/circos_SLE.txt"
OUTPUT_SVG="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_gene_density.svg"
OUTPUT_NUMT_HEATMAP="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_heatmap.svg"
KARYOTYPE_TXT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/conf/human_karyotype.txt"
R_SRC="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/script/2-NUMTs分布.R"
#? 使用R包Ideogram绘制

Rscript $R_SRC \
    --input "$INPUT_TXT" \
    --karyotype "$KARYOTYPE_TXT" \
    --gene_color "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b" \
    --numt_color "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b" \
    --numt_shape "triangle" \
    --bin_size 1000000 \
    --output "$OUTPUT_SVG" \
    --output_numt_heatmap "$OUTPUT_NUMT_HEATMAP"

# numt_shape 可选：box | circle | triangle
# bin_size: 聚合窗口大小（bp），用于缓解NUMT密集问题
# 输出两个图形：
#   1. NUMTs_gene_density.svg - 基因热图 + NUMT标记
#   2. NUMTs_heatmap.svg - NUMT密度热图（仅热图，无外围标记）
