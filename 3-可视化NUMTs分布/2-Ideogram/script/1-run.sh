#!/bin/bash

INPUT_TXT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/input/circos_10K.txt"
OUTPUT_SVG="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_gene_density.svg"
KARYOTYPE_TXT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/conf/human_karyotype.txt"
R_SRC="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/script/2-NUMTs分布.R"
#? 使用R包Ideogram绘制

Rscript $R_SRC \
    --input "$INPUT_TXT" \
    --karyotype "$KARYOTYPE_TXT" \
    --gene_color "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b" \
    --numt_color "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b" \
    --numt_shape "triangle" \
    --output "$OUTPUT_SVG"

# numt_shape 可选：box | circle | triangle
