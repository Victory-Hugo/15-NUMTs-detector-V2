#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

get_arg_value <- function(args, flag, default = NA_character_) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

chrom_order_value <- function(chr_vec) {
  suffix <- sub("^chr", "", chr_vec, ignore.case = TRUE)
  out <- suppressWarnings(as.integer(suffix))
  out[is.na(out) & suffix == "X"] <- 23L
  out[is.na(out) & suffix == "Y"] <- 24L
  out[is.na(out)] <- 999L
  out
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- get_arg_value(args, "--input", NA_character_)
output_pdf <- get_arg_value(args, "--output-pdf", NA_character_)
output_png <- get_arg_value(args, "--output-png", NA_character_)
output_svg <- get_arg_value(args, "--output-svg", NA_character_)
frequency_denominator <- as.numeric(get_arg_value(args, "--frequency-denominator"))

if (is.na(input_file) || !nzchar(input_file)) {
  stop("Missing required --input")
}
if (is.na(output_pdf) || !nzchar(output_pdf)) {
  stop("Missing required --output-pdf")
}
if (is.na(output_png) || !nzchar(output_png)) {
  stop("Missing required --output-png")
}
if (is.na(output_svg) || !nzchar(output_svg)) {
  stop("Missing required --output-svg")
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_svg), recursive = TRUE, showWarnings = FALSE)

cluster_df <- read_tsv(input_file, show_col_types = FALSE)
required_cols <- c("merged_cluster_id", "chr", "cluster_min_midpoint", "cluster_max_midpoint", "cluster_midpoint_mean", "POS", "sample_count")
missing_cols <- setdiff(required_cols, names(cluster_df))
if (length(missing_cols) > 0) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "))
}

chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
chromosome_colors <- c(
  "#d43e3e", "#dd6a3b", "#e58f38", "#eab236", "#e8cc36", "#d9df3a",
  "#b8e446", "#8fe051", "#62d76a", "#3dcd8b", "#37c3b3", "#3fb6d4",
  "#4ca1e8", "#5a86ee", "#6769ec", "#7a4fe5", "#9440d9", "#b139c7",
  "#cb3ab0", "#d94795", "#dd5d7e", "#d66a6a", "#424242", "#1b1a1a"
)
names(chromosome_colors) <- chr_levels

plot_df <- cluster_df %>%
  filter(chr %in% chr_levels) %>%
  mutate(
    CHR = factor(chr, levels = chr_levels, ordered = TRUE),
    chr_order = chrom_order_value(chr)
  ) %>%
  arrange(chr_order, cluster_midpoint_mean)

chr_sizes <- plot_df %>%
  group_by(CHR, chr_order) %>%
  summarise(chr_max = max(cluster_max_midpoint), .groups = "drop") %>%
  arrange(chr_order) %>%
  mutate(offset = lag(cumsum(chr_max), default = 0))

plot_df <- plot_df %>%
  left_join(chr_sizes, by = c("CHR", "chr_order")) %>%
  mutate(x = offset + cluster_midpoint_mean)

chr_labels <- plot_df %>%
  group_by(CHR) %>%
  summarise(center = mean(x), .groups = "drop") %>%
  arrange(match(CHR, chr_levels))

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold")
  )

p1 <- ggplot(plot_df, aes(x = x, y = POS, color = CHR)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = chromosome_colors, drop = FALSE) +
  scale_x_continuous(breaks = chr_labels$center, labels = chr_labels$CHR) +
  coord_cartesian(ylim = c(0, 1000)) +
  labs(
    x = "Chromosome",
    y = "Number of unique breakpoint midpoints per clustered NUMT locus",
    title = "Genome-wide distribution of clustered NUMT loci by cluster size"
  ) +
  base_theme

p2 <- ggplot(plot_df, aes(x = x, y = sample_count, color = CHR)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = chromosome_colors, drop = FALSE) +
  scale_x_continuous(breaks = chr_labels$center, labels = chr_labels$CHR) +
  coord_cartesian(ylim = c(0, 1000)) +
  labs(
    x = "Chromosome",
    y = "Number of samples sharing each clustered NUMT locus",
    title = "Genome-wide distribution of clustered NUMT loci by sample recurrence"
  ) +
  base_theme

combined_plot <- wrap_plots(p1, p2, ncol = 2) + plot_annotation(tag_levels = "A")
combined_plot <- combined_plot + plot_annotation(
  subtitle = paste0(
    "Frequency definition: sample_count / ", frequency_denominator,
    ". Class thresholds used in ideogram plots: singleton = 1 sample; rare = 2-83; ",
    "low-frequency = 84-418 (1% to <5%); common >= 419 (>=5%)."
  )
)

ggsave(filename = output_pdf, plot = combined_plot, width = 15, height = 7, units = "in")
ggsave(filename = output_png, plot = combined_plot, width = 15, height = 7, units = "in", dpi = 300)
ggsave(filename = output_svg, plot = combined_plot, width = 15, height = 7, units = "in", device = svglite)
