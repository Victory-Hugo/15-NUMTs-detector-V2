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
marker_input_file <- get_arg_value(args, "--marker-input", NA_character_)
output_marker <- get_arg_value(args, "--output-marker", NA_character_)
output_heatmap <- get_arg_value(args, "--output-heatmap", NA_character_)
output_heatmap_1kb <- get_arg_value(args, "--output-heatmap-1kb", NA_character_)
conf_karyotype <- get_arg_value(args, "--karyotype", NA_character_)
bin_size <- as.numeric(get_arg_value(args, "--bin-size", "1000000"))
frequency_denominator <- as.numeric(get_arg_value(args, "--frequency-denominator"))

gene_palette <- c("#273e55", "#6aa18a", "#b5c3b1", "#d7bbaf", "#b26966", "#782b55", "#591c4b")
marker_class_palette <- c(
  "ultra-rare" = "#0072b2",
  "rare" = "#56b4e9",
  "low-frequency" = "#009e73",
  "common" = "#f0e442"
)
heatmap_class_palette <- c(
  "singleton" = "#0072b2",
  "rare" = "#56b4e9",
  "low-frequency" = "#009e73",
  "common" = "#f0e442"
)
heatmap_palette <- c(
  "#2d465a",
  "#97b5a1",
  "#a75e62",
  "#5f1f4d"
)

annotate_svg_with_frequency_legend <- function(svg_path, mode = c("marker", "heatmap"), bin_label = NULL, denominator) {
  mode <- match.arg(mode)
  if (!file.exists(svg_path)) {
    return(invisible(NULL))
  }

  xml_escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }

  svg_text <- paste(readLines(svg_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  svg_text <- gsub(">NUMT</text>", ">Distinct NUMT cluster</text>", svg_text, fixed = TRUE)
  svg_text <- gsub(">Low</text>", ">singleton</text>", svg_text, fixed = TRUE)
  svg_text <- gsub(">High</text>", ">common</text>", svg_text, fixed = TRUE)

  box_x <- if (mode == "marker") 500 else 455
  box_y <- if (mode == "marker") 168 else 145
  box_width <- if (mode == "marker") 215 else 260
  box_height <- if (mode == "marker") 110 else 132
  text_x <- box_x + 8
  text_y <- box_y + 16
  line_step <- 12

  header_line <- if (mode == "marker") {
    "Marker color = frequency class of each distinct NUMT cluster"
  } else {
    paste0("Heatmap color = highest cluster frequency class within each ", bin_label, " bin")
  }

  if (mode == "marker") {
    common_min <- ceiling(0.05 * denominator)
    low_min <- ceiling(0.01 * denominator)
    rare_min <- ceiling(0.001 * denominator)
    legend_lines <- c(
      header_line,
      paste0("Frequency = sample_count / ", denominator),
      paste0("ultra-rare: sample_count < ", rare_min, " (<0.1%)"),
      paste0("rare: sample_count = ", rare_min, "-", low_min - 1, " (0.1% to <1%)"),
      paste0("low-frequency: sample_count = ", low_min, "-", common_min - 1, " (1% to <5%)"),
      paste0("common: sample_count >= ", common_min, " (>=5%)")
    )
  } else {
    legend_lines <- c(
      header_line,
      paste0("Frequency = sample_count / ", denominator),
      "singleton: sample_count = 1",
      "rare: sample_count = 2-83",
      "low-frequency: sample_count = 84-418 (1% to <5%)",
      "common: sample_count >= 419 (>=5%)"
    )
  }

  legend_text <- paste(
    vapply(seq_along(legend_lines), function(i) {
      sprintf(
        "<text x=\"%s\" y=\"%s\" font-size=\"9\" font-family=\"Arial\" fill=\"black\">%s</text>",
        text_x,
        text_y + (i - 1) * line_step,
        xml_escape(legend_lines[[i]])
      )
    }, character(1)),
    collapse = ""
  )

  annotation_block <- sprintf(
    "<rect x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\" style=\"fill:white;fill-opacity:0.92;stroke:#666666;stroke-width:0.6\"/>%s",
    box_x, box_y, box_width, box_height, legend_text
  )

  svg_text <- sub("</svg>$", paste0(annotation_block, "</svg>"), svg_text)
  writeLines(svg_text, svg_path, useBytes = TRUE)
  invisible(NULL)
}

if (is.na(input_file) || !nzchar(input_file)) {
  stop("Missing required --input")
}
if (is.na(marker_input_file) || !nzchar(marker_input_file)) {
  stop("Missing required --marker-input")
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
if (!file.exists(marker_input_file)) {
  stop("marker input file not found: ", marker_input_file)
}

marker_dir <- dirname(output_marker)
heatmap_dir <- dirname(output_heatmap)
if (!dir.exists(marker_dir)) {
  dir.create(marker_dir, recursive = TRUE)
}
if (!dir.exists(heatmap_dir)) {
  dir.create(heatmap_dir, recursive = TRUE)
}
if (!is.na(output_heatmap_1kb) && nzchar(output_heatmap_1kb) && !dir.exists(dirname(output_heatmap_1kb))) {
  dir.create(dirname(output_heatmap_1kb), recursive = TRUE)
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
marker_df <- read.table(
  marker_input_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

required_cols <- c("Chr", "Start", "End")
if (all(required_cols %in% names(input_df))) {
  input_df <- input_df[, required_cols]
} else if (all(c("sampleID", "region_chr", "region_start", "region_end") %in% names(input_df))) {
  input_df <- data.frame(
    Chr = as.character(input_df$region_chr),
    Start = as.numeric(input_df$region_start),
    End = as.numeric(input_df$region_end),
    stringsAsFactors = FALSE
  )
} else {
  stop(
    "Input must contain either Chr/Start/End or sampleID/region_chr/region_start/region_end"
  )
}

input_df <- input_df[!is.na(input_df$Start) & !is.na(input_df$End), ]
input_df <- input_df[input_df$Chr %in% karyotype$Chr, ]

if (nrow(input_df) == 0) {
  stop("No valid NUMT events after cleaning")
}

required_marker_cols <- c(
  "chr", "cluster_min_midpoint", "cluster_max_midpoint", "cluster_midpoint_mean",
  "frequency", "sample_count", "frequency_class"
)
missing_marker_cols <- setdiff(required_marker_cols, names(marker_df))
if (length(missing_marker_cols) > 0) {
  stop("Marker input is missing columns: ", paste(missing_marker_cols, collapse = ", "))
}

marker_df <- data.frame(
  Chr = as.character(marker_df$chr),
  Start = as.numeric(marker_df$cluster_min_midpoint),
  End = as.numeric(marker_df$cluster_max_midpoint),
  Midpoint = as.numeric(marker_df$cluster_midpoint_mean),
  Value = as.numeric(marker_df$frequency),
  MarkerClass = as.character(marker_df$frequency_class),
  HeatmapClass = ifelse(
    as.numeric(marker_df$sample_count) == 1,
    "singleton",
    ifelse(
      as.character(marker_df$frequency_class) %in% c("rare", "ultra-rare"),
      "rare",
      as.character(marker_df$frequency_class)
    )
  ),
  stringsAsFactors = FALSE
)
marker_df <- marker_df[!is.na(marker_df$Start) & !is.na(marker_df$End) & !is.na(marker_df$Value), ]
marker_df <- marker_df[marker_df$Chr %in% karyotype$Chr, ]

data(gene_density, package = "RIdeogram")
gene_density$Chr <- as.character(gene_density$Chr)
gene_density$Chr <- ifelse(
  grepl("^chr", gene_density$Chr, ignore.case = TRUE),
  gene_density$Chr,
  paste0("chr", gene_density$Chr)
)

karyotype_end <- setNames(karyotype$End, karyotype$Chr)
marker_df$End <- ifelse(
  marker_df$End > karyotype_end[marker_df$Chr],
  karyotype_end[marker_df$Chr],
  marker_df$End
)
marker_df <- marker_df[marker_df$Start < marker_df$End, ]

if (nrow(marker_df) == 0) {
  stop("No valid marker-frequency data after binning")
}

marker_df$Class <- factor(
  marker_df$HeatmapClass,
  levels = c("singleton", "rare", "low-frequency", "common")
)
marker_df$ClassRank <- as.numeric(marker_df$Class)

marker_df$HeatmapColor <- heatmap_class_palette[as.character(marker_df$Class)]
marker_df$MarkerClass <- factor(
  marker_df$MarkerClass,
  levels = c("ultra-rare", "rare", "low-frequency", "common")
)
marker_df$MarkerColor <- marker_class_palette[as.character(marker_df$MarkerClass)]

build_heatmap_data <- function(marker_data, karyotype_data, karyotype_end_map, current_bin_size) {
  heatmap_bins <- data.frame(
    Chr = marker_data$Chr,
    Start = floor((marker_data$Midpoint - 1) / current_bin_size) * current_bin_size,
    End = floor((marker_data$Midpoint - 1) / current_bin_size) * current_bin_size + current_bin_size,
    Value = marker_data$ClassRank,
    stringsAsFactors = FALSE
  )
  heatmap_bins$End <- ifelse(
    heatmap_bins$End > karyotype_end_map[heatmap_bins$Chr],
    karyotype_end_map[heatmap_bins$Chr],
    heatmap_bins$End
  )
  heatmap_bins <- heatmap_bins[heatmap_bins$Start < heatmap_bins$End, ]
  aggregate(Value ~ Chr + Start + End, data = heatmap_bins, FUN = max)
}

numt_labels <- data.frame(
  Type = "NUMT",
  Shape = "box",
  Chr = marker_df$Chr,
  Start = marker_df$Start,
  End = marker_df$End,
  color = gsub("#", "", as.character(marker_df$MarkerColor)),
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
annotate_svg_with_frequency_legend(
  svg_path = output_marker,
  mode = "marker",
  denominator = frequency_denominator
)

ideogram(
  karyotype = karyotype,
  overlaid = build_heatmap_data(marker_df, karyotype, karyotype_end, bin_size),
  colorset1 = heatmap_palette,
  output = output_heatmap
)
annotate_svg_with_frequency_legend(
  svg_path = output_heatmap,
  mode = "heatmap",
  bin_label = paste0(format(bin_size, scientific = FALSE, trim = TRUE), " bp"),
  denominator = frequency_denominator
)

if (!is.na(output_heatmap_1kb) && nzchar(output_heatmap_1kb)) {
  ideogram(
    karyotype = karyotype,
    overlaid = build_heatmap_data(marker_df, karyotype, karyotype_end, 1000),
    colorset1 = heatmap_palette,
    output = output_heatmap_1kb
  )
  annotate_svg_with_frequency_legend(
    svg_path = output_heatmap_1kb,
    mode = "heatmap",
    bin_label = "1 kb",
    denominator = frequency_denominator
  )
}
