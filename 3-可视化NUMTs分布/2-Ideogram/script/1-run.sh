#!/bin/bash

INPUT_TXT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/input/circos_SLE.txt"
OUTPUT_SVG="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_gene_density.svg"
OUTPUT_NUMT_HEATMAP="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_heatmap.svg"
OUTPUT_PNG="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_gene_density.png"
OUTPUT_NUMT_HEATMAP_PNG="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/output/NUMTs_heatmap.png"
PNG_DPI=600
KARYOTYPE_TXT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/conf/human_karyotype.txt"
R_SRC="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-可视化NUMTs分布/2-Ideogram/script/2-NUMTs分布.R"
#? 使用R包Ideogram绘制
# "quantile" 或改成 "log" ,"linear"
# bin_size: 聚合窗口大小（bp），设置越大，运算速度越快
# Rscript $R_SRC \
#     --input "$INPUT_TXT" \
#     --karyotype "$KARYOTYPE_TXT" \
#     --gene_color "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b" \
#     --numt_color "#1f3a5f,#3f6fa6,#6fa6c9,#d9b1a3,#c46a6a,#8b2e3c" \
#     --numt_color_mode "quantile" \
#     --numt_shape "triangle" \
#     --bin_size 10000 \
#     --output "$OUTPUT_SVG" \
#     --output_numt_heatmap "$OUTPUT_NUMT_HEATMAP"

# Convert SVG to PNG (prefer rsvg-convert, fallback to inkscape)
if command -v rsvg-convert >/dev/null 2>&1; then
    rsvg-convert -a --dpi-x "$PNG_DPI" --dpi-y "$PNG_DPI" -o "$OUTPUT_PNG" "$OUTPUT_SVG"
    rsvg-convert -a --dpi-x "$PNG_DPI" --dpi-y "$PNG_DPI" -o "$OUTPUT_NUMT_HEATMAP_PNG" "$OUTPUT_NUMT_HEATMAP"
elif command -v inkscape >/dev/null 2>&1; then
    inkscape "$OUTPUT_SVG" --export-type=png --export-filename="$OUTPUT_PNG" --export-dpi="$PNG_DPI" >/dev/null 2>&1
    inkscape "$OUTPUT_NUMT_HEATMAP" --export-type=png --export-filename="$OUTPUT_NUMT_HEATMAP_PNG" --export-dpi="$PNG_DPI" >/dev/null 2>&1
else
    echo "SVG->PNG failed: install rsvg-convert or inkscape." >&2
fi

# numt_shape 可选：box | circle | triangle
# bin_size: 聚合窗口大小（bp），用于缓解NUMT密集问题
# 输出两个图形：
#   1. NUMTs_gene_density.svg - 基因热图 + NUMT标记
#   2. NUMTs_heatmap.svg - NUMT密度热图（仅热图，无外围标记）
