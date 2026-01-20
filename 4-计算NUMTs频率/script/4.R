library(tidyverse)
library(tidyplots)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 4.R <output_dir>")
}
out_dir <- args[1]
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

length_tsv <- file.path(out_dir, "1-numt_length_hist.tsv")
length_short_tsv <- file.path(out_dir, "2-numt_length_hist_short_lt400.tsv")
freq_pie_tsv <- file.path(out_dir, "3-numt_frequency_pie.tsv")
sample_prop_tsv <- file.path(out_dir, "4-sample_proportion_by_freq.tsv")

length_df <- read_tsv(length_tsv, show_col_types = FALSE)
length_short_df <- read_tsv(length_short_tsv, show_col_types = FALSE)
freq_pie_df <- read_tsv(freq_pie_tsv, show_col_types = FALSE) |>
  mutate(
    pct = count / sum(count),
    label = scales::percent(pct, accuracy = 0.1),
    y_pos = cumsum(count) - count / 2,
    x_pos = 1
  )
sample_prop_df <- read_tsv(sample_prop_tsv, show_col_types = FALSE) |>
  mutate(
    freq_class = factor(freq_class, levels = c("common", "rare", "ultra-rare", "private")),
    label = scales::percent(proportion, accuracy = 0.1)
  )

length_df |>
  tidyplot(x = length_bp) |>
  add_histogram(bins = 1000, fill = "#2980b6") |>
  adjust_x_axis_title("Length (bp)") |>
  adjust_y_axis_title("Count") |>
  add_title("NUMT length distribution") |>
  save_plot(file.path(out_dir, "1-numt_length_hist.pdf"))



length_short_df |>
  tidyplot(x = length_bp) |>
  # add_annotation_text(
  #   x = 200,
  #   y = 1000,
  #   text = sprintf("Mean = %.1f bp SE = %.1f bp", mean_len, se_len),
  #   color = "#d43e3e"
  # ) |>
  add_histogram(
    data = length_short_df |> filter(is_short),
    bins = 50,
    alpha = 0.85,
    fill = "#d43e3e"
  ) |>
  adjust_x_axis(limits = c(0, 400)) |>
  adjust_x_axis_title("Length (bp)") |>
  adjust_y_axis_title("Count") |>
  add_title("NUMT length distribution (highlight < 400 bp)") |>
  save_plot(file.path(out_dir, "2-numt_length_hist_short_lt400.pdf"))



freq_pie_df |>
  tidyplot(color = freq_class, y = count) |>
  add_pie() |>
  add_annotation_text(
    x = freq_pie_df$x_pos,
    y = freq_pie_df$y_pos,
    text = freq_pie_df$label,
    fontsize = 8
  ) |>
  add_title("NUMT population frequency categories") |>
  save_plot(file.path(out_dir, "3-numt_frequency_pie.pdf"))



sample_prop_df |>
  tidyplot(x = freq_class, y = proportion) |>
  add_sum_bar(width = 0.6) |>
  add_data_labels(label = label, label_position = "above") |>
  # adjust_y_axis(limits = c(0, 1)) |>
  adjust_y_axis_title("Proportion of samples") |>
  adjust_title("Samples carrying >= 1 NUMT per frequency class") |>
  save_plot(file.path(out_dir, "4-sample_proportion_by_freq.pdf"))
