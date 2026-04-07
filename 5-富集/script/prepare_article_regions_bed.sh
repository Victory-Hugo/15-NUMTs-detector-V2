#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  prepare_article_regions_bed.sh --output-dir DIR --source-dir DIR --tmp-dir DIR --bedtools-bin bedtools [--raw-dir DIR] [--force]

Generate hg38 BED files for the nuclear genomic features used in the paper's
Figure 7D-style NUMT breakpoint enrichment analysis.
USAGE
}

OUTPUT_DIR=""
SOURCE_DIR=""
TMP_DIR=""
RAW_DIR=""
BEDTOOLS_BIN="bedtools"
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --source-dir) SOURCE_DIR="$2"; shift 2 ;;
        --tmp-dir) TMP_DIR="$2"; shift 2 ;;
        --raw-dir) RAW_DIR="$2"; shift 2 ;;
        --bedtools-bin) BEDTOOLS_BIN="$2"; shift 2 ;;
        --force) FORCE="true"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" || -z "$SOURCE_DIR" || -z "$TMP_DIR" ]]; then
    usage
    exit 1
fi

mkdir -p "$OUTPUT_DIR" "$SOURCE_DIR" "$TMP_DIR"
if [[ -n "$RAW_DIR" ]]; then
    mkdir -p "$RAW_DIR"
fi

download_if_needed() {
    local url="$1"
    local dest="$2"
    if [[ "$FORCE" == "true" || ! -s "$dest" ]]; then
        wget -O "$dest" "$url"
    fi
}

first_existing() {
    local candidate
    for candidate in "$@"; do
        if [[ -n "$candidate" && -s "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

read_table() {
    local path="$1"
    case "$path" in
        *.gz) zcat "$path" ;;
        *) cat "$path" ;;
    esac
}

std_bed() {
    local input="$1"
    local output="$2"
    awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {print $1,$2,$3}' "$input" \
        | sort -k1,1 -k2,2n \
        | "$BEDTOOLS_BIN" merge -i - > "$output"
}

write_awk() {
    local script_path="$1"
    local script_body="$2"
    printf '%s\n' "$script_body" > "$script_path"
}

UCSC_BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz"
CCRE_URL="https://screen.encodeproject.org/media/files/GRCh38-cCREs.bed.gz"
CCRE_BB_URL="https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"

CPG_FILE=$(first_existing "${RAW_DIR:+$RAW_DIR/5-CpG岛.tsv}" "$SOURCE_DIR/cpgIslandExt.txt.gz" || true)
SIMPLE_REPEAT_FILE=$(first_existing "${RAW_DIR:+$RAW_DIR/7-简单重复.tsv}" "$SOURCE_DIR/simpleRepeat.txt.gz" || true)
MICROSAT_FILE=$(first_existing "${RAW_DIR:+$RAW_DIR/6-微卫星.tsv}" || true)
RMSK_FILE=$(first_existing "${RAW_DIR:+$RAW_DIR/11-全基因组重复序列注释表.txt}" "$SOURCE_DIR/rmsk.txt.gz" || true)
GENCODE_FILE=$(first_existing "${RAW_DIR:+$RAW_DIR/10-gencode.v43.annotation.gtf}" "$SOURCE_DIR/gencode.v43.annotation.gtf.gz" || true)
CCRE_BB_FILE=$(first_existing "${RAW_DIR:+$RAW_DIR/8-非基因型的功能元件.bb}" "$SOURCE_DIR/encodeCcreCombined.bb" || true)

if [[ -z "$CPG_FILE" ]]; then
    download_if_needed "$UCSC_BASE/cpgIslandExt.txt.gz" "$SOURCE_DIR/cpgIslandExt.txt.gz"
    CPG_FILE="$SOURCE_DIR/cpgIslandExt.txt.gz"
fi
if [[ -z "$SIMPLE_REPEAT_FILE" ]]; then
    download_if_needed "$UCSC_BASE/simpleRepeat.txt.gz" "$SOURCE_DIR/simpleRepeat.txt.gz"
    SIMPLE_REPEAT_FILE="$SOURCE_DIR/simpleRepeat.txt.gz"
fi
if [[ -z "$MICROSAT_FILE" ]]; then
    MICROSAT_FILE="$SIMPLE_REPEAT_FILE"
fi
if [[ -z "$RMSK_FILE" ]]; then
    download_if_needed "$UCSC_BASE/rmsk.txt.gz" "$SOURCE_DIR/rmsk.txt.gz"
    RMSK_FILE="$SOURCE_DIR/rmsk.txt.gz"
fi
if [[ -z "$GENCODE_FILE" ]]; then
    download_if_needed "$GENCODE_URL" "$SOURCE_DIR/gencode.v43.annotation.gtf.gz"
    GENCODE_FILE="$SOURCE_DIR/gencode.v43.annotation.gtf.gz"
fi
if [[ -z "$CCRE_BB_FILE" ]]; then
    download_if_needed "$CCRE_BB_URL" "$SOURCE_DIR/encodeCcreCombined.bb"
    CCRE_BB_FILE="$SOURCE_DIR/encodeCcreCombined.bb"
fi
download_if_needed "$UCSC_BASE/genomicSuperDups.txt.gz" "$SOURCE_DIR/genomicSuperDups.txt.gz"
download_if_needed "$UCSC_BASE/gap.txt.gz" "$SOURCE_DIR/gap.txt.gz"
download_if_needed "$UCSC_BASE/centromeres.txt.gz" "$SOURCE_DIR/centromeres.txt.gz"

read_table "$CPG_FILE" | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > "$TMP_DIR/CpG_islands.raw.bed"
std_bed "$TMP_DIR/CpG_islands.raw.bed" "$OUTPUT_DIR/05-CpG_islands.bed"

read_table "$SIMPLE_REPEAT_FILE" | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > "$TMP_DIR/Simple_Repeats.raw.bed"
std_bed "$TMP_DIR/Simple_Repeats.raw.bed" "$OUTPUT_DIR/21-Simple_Repeats.bed"

if [[ "$(basename "$MICROSAT_FILE")" == "6-微卫星.tsv" ]]; then
    read_table "$MICROSAT_FILE" | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > "$TMP_DIR/Microsat.raw.bed"
else
    read_table "$MICROSAT_FILE" | awk 'BEGIN{OFS="\t"} $6 >= 2 && $6 <= 6 {print $2,$3,$4}' > "$TMP_DIR/Microsat.raw.bed"
fi
std_bed "$TMP_DIR/Microsat.raw.bed" "$OUTPUT_DIR/06-Microsat.bed"

zcat "$SOURCE_DIR/genomicSuperDups.txt.gz" | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > "$TMP_DIR/Genomic_superdups.raw.bed"
std_bed "$TMP_DIR/Genomic_superdups.raw.bed" "$OUTPUT_DIR/22-Genomic_superdups.bed"

if zcat "$SOURCE_DIR/centromeres.txt.gz" | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > "$TMP_DIR/Centromeres.raw.bed" && [[ -s "$TMP_DIR/Centromeres.raw.bed" ]]; then
    std_bed "$TMP_DIR/Centromeres.raw.bed" "$OUTPUT_DIR/15-Centromeres.bed"
else
    zcat "$SOURCE_DIR/gap.txt.gz" | awk 'BEGIN{OFS="\t"} $8=="centromere"{print $2,$3,$4}' > "$TMP_DIR/Centromeres.raw.bed"
    std_bed "$TMP_DIR/Centromeres.raw.bed" "$OUTPUT_DIR/15-Centromeres.bed"
fi

read_table "$RMSK_FILE" | awk 'BEGIN{OFS="\t"} $12=="LINE"{print $6,$7,$8}' > "$TMP_DIR/LINE.raw.bed"
read_table "$RMSK_FILE" | awk 'BEGIN{OFS="\t"} $12=="SINE"{print $6,$7,$8}' > "$TMP_DIR/SINE.raw.bed"
read_table "$RMSK_FILE" | awk 'BEGIN{OFS="\t"} $12=="LTR"{print $6,$7,$8}' > "$TMP_DIR/LTR.raw.bed"
read_table "$RMSK_FILE" | awk 'BEGIN{OFS="\t"} $12=="DNA"{print $6,$7,$8}' > "$TMP_DIR/rmsk-DNA.raw.bed"
read_table "$RMSK_FILE" | awk 'BEGIN{OFS="\t"} $12=="Satellite"{print $6,$7,$8}' > "$TMP_DIR/Satellite.raw.bed"
read_table "$RMSK_FILE" | awk 'BEGIN{OFS="\t"} $12 ~ /Retroposon/ || $13 ~ /Retroposon/ || $11 ~ /Retroposon/ {print $6,$7,$8}' > "$TMP_DIR/Retroposon.raw.bed"
std_bed "$TMP_DIR/LINE.raw.bed" "$OUTPUT_DIR/20-LINE.bed"
std_bed "$TMP_DIR/SINE.raw.bed" "$OUTPUT_DIR/13-SINE.bed"
std_bed "$TMP_DIR/LTR.raw.bed" "$OUTPUT_DIR/17-LTR.bed"
std_bed "$TMP_DIR/rmsk-DNA.raw.bed" "$OUTPUT_DIR/14-rmsk-DNA.bed"
std_bed "$TMP_DIR/Satellite.raw.bed" "$OUTPUT_DIR/16-Satellite.bed"
std_bed "$TMP_DIR/Retroposon.raw.bed" "$OUTPUT_DIR/10-Retroposon.bed"

bigBedToBed "$CCRE_BB_FILE" "$TMP_DIR/FuncElems.full.bed"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$TMP_DIR/FuncElems.full.bed" > "$TMP_DIR/FuncElems.raw.bed"
cp "$TMP_DIR/FuncElems.raw.bed" "$TMP_DIR/Regulatory_elements.raw.bed"
std_bed "$TMP_DIR/FuncElems.raw.bed" "$OUTPUT_DIR/09-FuncElems.bed"
std_bed "$TMP_DIR/Regulatory_elements.raw.bed" "$OUTPUT_DIR/26-Regulatory_elements.bed"

"$(dirname "$0")/../python/prepare_gencode_article_regions.py" \
    --gencode-gtf-gz "$GENCODE_FILE" \
    --tmp-dir "$TMP_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --bedtools-bin "$BEDTOOLS_BIN"

cat > "$OUTPUT_DIR/README.article_regions.md" <<'EOF'
# Article Figure 7D BED sources

These BED files were generated for paper-aligned NUMT nuclear breakpoint enrichment.

- UCSC hg38 database: cpgIslandExt, simpleRepeat, genomicSuperDups, rmsk, gap/centromeres.
- ENCODE SCREEN GRCh38 cCREs: used as the practical source for FuncElems and Regulatory_elements.
- GENCODE v43 GTF: used for gene, exon, CDS, UTR, start_codon, stop_codon, snRNA, snoRNA_miRNA and intron.

All BED files were restricted to chr1-22, chrX and chrY, sorted and merged.
EOF

printf 'Generated article region BED files in %s\n' "$OUTPUT_DIR"
