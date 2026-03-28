#!/usr/bin/env Rscript

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

args <- commandArgs(trailingOnly = TRUE)

input_file <- get_arg_value(args, "--input", NA_character_)
output_marker <- get_arg_value(args, "--output-marker", NA_character_)
output_heatmap <- get_arg_value(args, "--output-heatmap", NA_character_)
conf_karyotype <- get_arg_value(args, "--karyotype", NA_character_)
bin_size <- as.numeric(get_arg_value(args, "--bin-size", "1000000"))

gene_palette <- c("#273e55", "#6aa18a", "#b5c3b1", "#d7bbaf", "#b26966", "#782b55", "#591c4b")
numt_palette <- c(
  "#0b0405", "#30203e", "#3e4d93", "#366b9f", "#3488a6",
  "#36a4ab", "#49c1ad", "#60ceac", "#84d8b0", "#c4e9cf", "#def5e5"
)

if (is.na(input_file) || !nzchar(input_file)) {
  stop("Missing required --input")
}
if (is.na(output_marker) || !nzchar(output_marker)) {
  stop("Missing required --output-marker")
}
if (is.na(output_heatmap) || !nzchar(output_heatmap)) {
  stop("Missing required --output-heatmap")
}
if (is.na(conf_karyotype) || !nzchar(conf_karyotype)) {
  stop("Missing required --karyotype")
}
if (!file.exists(conf_karyotype)) {
  stop("karyotype file not found: ", conf_karyotype)
}

marker_dir <- dirname(output_marker)
heatmap_dir <- dirname(output_heatmap)
if (!dir.exists(marker_dir)) {
  dir.create(marker_dir, recursive = TRUE)
}
if (!dir.exists(heatmap_dir)) {
  dir.create(heatmap_dir, recursive = TRUE)
}

pdf(NULL)
on.exit(dev.off(), add = TRUE)

karyotype <- read.table(
  conf_karyotype,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("Chr", "Start", "End")
)

input_df <- read.table(
  input_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

required_cols <- c("Chr", "Start", "End")
missing_cols <- setdiff(required_cols, names(input_df))
if (length(missing_cols) > 0) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "))
}

input_df$Chr <- as.character(input_df$Chr)
input_df$Start <- as.numeric(input_df$Start)
input_df$End <- as.numeric(input_df$End)
input_df <- input_df[!is.na(input_df$Start) & !is.na(input_df$End), ]
input_df <- input_df[input_df$Chr %in% karyotype$Chr, ]

if (nrow(input_df) == 0) {
  stop("No valid NUMT events after cleaning")
}

data(gene_density, package = "RIdeogram")
gene_density$Chr <- as.character(gene_density$Chr)
gene_density$Chr <- ifelse(
  grepl("^chr", gene_density$Chr, ignore.case = TRUE),
  gene_density$Chr,
  paste0("chr", gene_density$Chr)
)

midpoint <- floor((input_df$Start + input_df$End) / 2)

numt_bins <- data.frame(
  Chr = input_df$Chr,
  Start = floor((midpoint - 1) / bin_size) * bin_size,
  End = floor((midpoint - 1) / bin_size) * bin_size + bin_size,
  Value = 1,
  stringsAsFactors = FALSE
)

numt_density <- aggregate(Value ~ Chr + Start + End, data = numt_bins, FUN = sum)
karyotype_end <- setNames(karyotype$End, karyotype$Chr)
numt_density$End <- ifelse(
  numt_density$End > karyotype_end[numt_density$Chr],
  karyotype_end[numt_density$Chr],
  numt_density$End
)
numt_density <- numt_density[numt_density$Start < numt_density$End, ]

if (nrow(numt_density) == 0) {
  stop("No valid density data after binning")
}

all_bins_list <- lapply(seq_len(nrow(karyotype)), function(i) {
  chr <- karyotype$Chr[i]
  chr_end <- karyotype$End[i]
  starts <- seq(0, chr_end - 1, by = bin_size)
  data.frame(
    Chr = chr,
    Start = starts,
    End = pmin(starts + bin_size, chr_end),
    stringsAsFactors = FALSE
  )
})
numt_full_bins <- do.call(rbind, all_bins_list)
numt_full_bins <- numt_full_bins[numt_full_bins$Start < numt_full_bins$End, ]
numt_full_bins$Value <- 0

numt_density_full <- merge(
  numt_full_bins,
  numt_density,
  by = c("Chr", "Start", "End"),
  all.x = TRUE,
  suffixes = c(".base", ".hit")
)
numt_density_full$Value <- ifelse(
  is.na(numt_density_full$Value.hit),
  numt_density_full$Value.base,
  numt_density_full$Value.hit
)
numt_density_full <- numt_density_full[, c("Chr", "Start", "End", "Value")]

value_range <- range(numt_density_full$Value)
if (diff(value_range) == 0) {
  breaks <- seq(value_range[1] - 0.5, value_range[2] + 0.5, length.out = length(numt_palette) + 1)
} else {
  breaks <- seq(value_range[1] - 0.1, value_range[2] + 0.1, length.out = length(numt_palette) + 1)
}

numt_density_full$Color <- cut(
  numt_density_full$Value,
  breaks = breaks,
  include.lowest = TRUE,
  labels = numt_palette
)

full_key <- paste(numt_density_full$Chr, numt_density_full$Start, numt_density_full$End, sep = ":")
hit_key <- paste(numt_density$Chr, numt_density$Start, numt_density$End, sep = ":")
numt_density$Color <- numt_density_full$Color[match(hit_key, full_key)]

numt_labels <- data.frame(
  Type = "NUMT",
  Shape = "box",
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
  output = output_marker
)

numt_heatmap_data <- numt_density_full[, c("Chr", "Start", "End", "Value")]
ideogram(
  karyotype = karyotype,
  overlaid = numt_heatmap_data,
  colorset1 = numt_palette,
  output = output_heatmap
)
