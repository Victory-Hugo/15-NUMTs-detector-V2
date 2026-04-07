#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pheatmap)
})

get_arg_value <- function(args, key, default = NA_character_) {
  match_idx <- match(key, args)
  if (is.na(match_idx) || match_idx == length(args)) {
    return(default)
  }
  args[[match_idx + 1]]
}

parse_csv <- function(value) {
  value <- trimws(value)
  if (!nzchar(value)) {
    return(character())
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

sanitize_region_label <- function(region_name) {
  label <- sub("\\.bed$", "", region_name)
  label <- sub("^[0-9]+-", "", label)
  label <- gsub("_", " ", label, fixed = TRUE)
  label
}

args <- commandArgs(trailingOnly = TRUE)
summary_tsv <- get_arg_value(args, "--summary-tsv")
output_dir <- get_arg_value(args, "--output-dir")
frequency_classes <- parse_csv(get_arg_value(args, "--frequency-classes", "all,ultra-rare,rare,low-frequency,common"))

if (is.na(summary_tsv) || !nzchar(summary_tsv)) {
  stop("Missing required --summary-tsv")
}
if (is.na(output_dir) || !nzchar(output_dir)) {
  stop("Missing required --output-dir")
}

summary_df <- read.delim(summary_tsv, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("region_name", "frequency_class", "p_value_greater")
missing_cols <- setdiff(required_cols, names(summary_df))
if (length(missing_cols) > 0) {
  stop("Summary TSV missing required columns: ", paste(missing_cols, collapse = ", "))
}

summary_df$p_value_greater <- as.numeric(summary_df$p_value_greater)
summary_df$region_label <- vapply(summary_df$region_name, sanitize_region_label, character(1))

region_order <- unique(summary_df$region_label[order(summary_df$region_name)])
frequency_order <- frequency_classes[frequency_classes %in% unique(summary_df$frequency_class)]
if (length(frequency_order) == 0) {
  frequency_order <- unique(summary_df$frequency_class)
}

p_matrix <- matrix(
  NA_real_,
  nrow = length(frequency_order),
  ncol = length(region_order),
  dimnames = list(frequency_order, region_order)
)

for (idx in seq_len(nrow(summary_df))) {
  row_name <- summary_df$frequency_class[[idx]]
  col_name <- summary_df$region_label[[idx]]
  if (row_name %in% rownames(p_matrix) && col_name %in% colnames(p_matrix)) {
    p_matrix[row_name, col_name] <- summary_df$p_value_greater[[idx]]
  }
}

triangle_matrix <- ifelse(!is.na(p_matrix) & p_matrix <= 0.05, "\u25b2", "")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_pdf <- file.path(output_dir, "enrichment_pvalue_heatmap.pdf")
output_png <- file.path(output_dir, "enrichment_pvalue_heatmap.png")
output_svg <- file.path(output_dir, "enrichment_pvalue_heatmap.svg")
matrix_out <- file.path(output_dir, "enrichment_pvalue_heatmap.matrix.tsv")
write.table(
  data.frame(frequency_class = rownames(p_matrix), p_matrix, check.names = FALSE),
  file = matrix_out,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

breaks <- seq(0, 1, length.out = 101)
palette <- colorRampPalette(c("#630000","#cf6057","#fab85b","#fbd488","#dcd3b8", "#7fb5d6", "#184c88", "#0d3670", "#022259"))(100)

draw_heatmap <- function(filename, device = c("pdf", "png", "svg")) {
  device <- match.arg(device)
  width <- max(9, 0.25 * ncol(p_matrix) + 3) # 调整宽度以适应列数
  height <- max(3.8, 0.20 * nrow(p_matrix) + 2.2) # 调整高度以适应行数
  if (device == "pdf") {
    cairo_pdf(filename, width = width, height = height, family = "sans")
  } else if (device == "png") {
    png(filename, width = width, height = height, units = "in", res = 300, type = "cairo")
  } else {
    svg(filename, width = width, height = height)
  }
  pheatmap(
    p_matrix,
    color = palette,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    border_color = "grey85",
    display_numbers = triangle_matrix,
    number_color = "black",
    fontsize_number = 9,
    fontsize = 8,
    fontsize_row = 8,
    fontsize_col = 7,
    angle_col = 45,
    legend = TRUE,
    na_col = "grey90",
    main = "NUMT nuclear breakpoint flank enrichment P values"
  )
  dev.off()
}

draw_heatmap(output_pdf, "pdf")
draw_heatmap(output_png, "png")
draw_heatmap(output_svg, "svg")

message("Wrote pheatmap P-value heatmap to ", output_pdf)
