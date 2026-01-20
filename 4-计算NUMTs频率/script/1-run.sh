#!/bin/bash

INPUT_DIR="/mnt/l/13_SLE_NUMT/1-所有的NUMTs/data"
INPUT_TSV_1="${INPUT_DIR}/all_individuals_ConfidentBreakpoints.tsv"
INPUT_TSV_2="${INPUT_DIR}/all_individuals_mt.disc.sam.cluster.tsv"

PY_SRC="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-计算NUMTs频率/script/3-统计频率.py"
PAIR_SRC="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-计算NUMTs频率/script/2-断点配对_长度聚类.py"
R_SRC="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-计算NUMTs频率/script/4.R"
DATA_DIR="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-计算NUMTs频率/data"
PAIRED_TSV="${DATA_DIR}/all_individuals_ConfidentBreakpoints.paired_len1000.tsv"
OUT_DIR_BASE="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-计算NUMTs频率/output"
OUT_DIR_1="${OUT_DIR_BASE}/frequency-断点"
OUT_DIR_2="${OUT_DIR_BASE}/frequency-聚类"

set -e

echo "Running NUMT frequency plots..."
python3 "$PAIR_SRC" \
  --input "$INPUT_TSV_1" \
  --output "$PAIRED_TSV" \
  --cluster-gap 1000 \
  --max-length 1000

python3 "$PY_SRC" \
  --confident "$INPUT_TSV_1" \
  --cluster "$PAIRED_TSV" \
  --outdir "$OUT_DIR_BASE" \
  --outdir-breakpoint "$OUT_DIR_1" \
  --outdir-cluster "$OUT_DIR_2" \
  --window 1000

echo "Running tidyplots rendering..."
Rscript "$R_SRC" "$OUT_DIR_2"

echo "Done."
