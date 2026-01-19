suppressPackageStartupMessages({
  library(RIdeogram)
})

get_arg_value <- function(args, flag, default = NA_character_) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

split_colors <- function(color_str) {
  if (is.na(color_str) || !nzchar(color_str)) {
    return(character(0))
  }
  trimws(unlist(strsplit(color_str, ",")))
}

args <- commandArgs(trailingOnly = TRUE)

input_file <- get_arg_value(args, "--input", NA_character_)
output_svg <- get_arg_value(args, "--output", NA_character_)
conf_karyotype <- get_arg_value(args, "--karyotype", NA_character_)
gene_color_str <- get_arg_value(
  args,
  "--gene_color",
  "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b"
)
numt_color_str <- get_arg_value(
  args,
  "--numt_color",
  "#273e55,#6aa18a,#b5c3b1,#d7bbaf,#b26966,#782b55,#591c4b"
)
numt_shape <- get_arg_value(args, "--numt_shape", "box")
bin_size <- as.numeric(get_arg_value(args, "--bin_size", "100000"))
min_len <- as.numeric(get_arg_value(args, "--min_len", "10000"))

gene_palette <- split_colors(gene_color_str)
numt_palette <- split_colors(numt_color_str)

if (is.na(input_file) || !nzchar(input_file)) {
  stop("Missing required --input")
}
if (is.na(output_svg) || !nzchar(output_svg)) {
  stop("Missing required --output")
}
if (is.na(conf_karyotype) || !nzchar(conf_karyotype)) {
  inferred_base <- dirname(dirname(input_file))
  conf_karyotype <- file.path(inferred_base, "conf", "human_karyotype.txt")
}

if (length(gene_palette) == 0) {
  stop("gene_color is empty")
}
if (length(numt_palette) == 0) {
  stop("numt_color is empty")
}

output_dir <- dirname(output_svg)


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Avoid creating Rplots.pdf
pdf(NULL)
on.exit(dev.off(), add = TRUE)

# karyotype
karyotype <- read.table(
  conf_karyotype,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("Chr", "Start", "End")
)

if (!file.exists(conf_karyotype)) {
  stop("karyotype file not found: ", conf_karyotype)
}

# use built-in gene density for stronger color contrast
data(gene_density, package = "RIdeogram")
gene_density$Chr <- as.character(gene_density$Chr)
gene_density$Chr <- ifelse(
  grepl("^chr", gene_density$Chr, ignore.case = TRUE),
  gene_density$Chr,
  paste0("chr", gene_density$Chr)
)

# NUMTs data
numt_raw <- read.table(
  input_file,
  sep = "",
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  col.names = c("Chr", "Start", "End", "MtChr", "MtStart", "MtEnd")
)

numt_raw$Chr <- sub("^hs", "chr", numt_raw$Chr, ignore.case = TRUE)
numt_raw$Start <- as.numeric(numt_raw$Start)
numt_raw$End <- as.numeric(numt_raw$End)

midpoint <- floor((numt_raw$Start + numt_raw$End) / 2)
half_min <- floor(min_len / 2)
numt_raw$Start <- pmax(0, midpoint - half_min)
numt_raw$End <- midpoint + half_min

numt_bins <- data.frame(
  Chr = numt_raw$Chr,
  Start = floor((midpoint - 1) / bin_size) * bin_size,
  End = floor((midpoint - 1) / bin_size) * bin_size + bin_size,
  Value = 1,
  stringsAsFactors = FALSE
)

numt_density <- aggregate(Value ~ Chr + Start + End, data = numt_bins, FUN = sum)

karyotype_end <- setNames(karyotype$End, karyotype$Chr)
numt_density <- numt_density[numt_density$Chr %in% karyotype$Chr, ]
numt_density$End <- ifelse(
  numt_density$End > karyotype_end[numt_density$Chr],
  karyotype_end[numt_density$Chr],
  numt_density$End
)
numt_density <- numt_density[numt_density$Start < numt_density$End, ]

numt_density$Color <- cut(
  numt_density$Value,
  breaks = length(numt_palette),
  include.lowest = TRUE,
  labels = numt_palette
)

numt_labels <- data.frame(
  Type = "NUMT",
  Shape = numt_shape,
  Chr = numt_density$Chr,
  Start = numt_density$Start,
  End = numt_density$End,
  color = gsub("#", "", as.character(numt_density$Color)),
  stringsAsFactors = FALSE
)

ideogram(
  karyotype = karyotype,
  overlaid = gene_density,
  label = numt_labels,
  label_type = "marker",
  colorset1 = gene_palette,
  output = output_svg
)
