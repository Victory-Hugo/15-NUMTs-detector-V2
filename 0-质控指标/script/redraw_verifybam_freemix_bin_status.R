library(readr)
library(dplyr)
library(ggplot2)
library(scales)
setwd('/mnt/l/20-NUMTs/2-质控指标')
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  "."
}
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_dir <- "/mnt/l/20-NUMTs/2-质控指标"

input_path <- file.path(project_dir, "META", "FREEMIX.tsv")
output_data_path <- file.path(project_dir, "img", "data", "verifybam_freemix_bin_status.csv")
output_plot_path <- file.path(project_dir, "img", "plots", "verifybam_freemix_bin_status.pdf")

dir.create(dirname(output_data_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_plot_path), recursive = TRUE, showWarnings = FALSE)

bin_levels <- c("<=0.001", "0.001-0.005", "0.005-0.01", "0.01-0.02", "0.02-0.05", ">0.05")
status_levels <- c("PASS", "FAIL")
status_colors <- c(
  "PASS" = "#0072b2",
  "FAIL" = "#56b4e9"
)

raw_df <- read_tsv(
  input_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) |>
  transmute(
    ID = as.character(ID),
    FREEMIX = suppressWarnings(as.numeric(FREEMIX))
  ) |>
  filter(!is.na(FREEMIX), FREEMIX >= 0)

plot_df <- raw_df |>
  mutate(
    status = if_else(FREEMIX <= 0.02, "PASS", "FAIL"),
    freemix_bin = cut(
      FREEMIX,
      breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.05, Inf),
      labels = bin_levels,
      include.lowest = TRUE,
      right = TRUE
    ),
    freemix_bin = factor(as.character(freemix_bin), levels = bin_levels),
    status = factor(status, levels = status_levels)
  )

summary_df <- plot_df |>
  count(freemix_bin, status, name = "n", .drop = FALSE) |>
  group_by(freemix_bin) |>
  mutate(
    proportion = if (sum(n) > 0) n / sum(n) else rep(0, dplyr::n())
  ) |>
  ungroup()

write_csv(summary_df, output_data_path)

label_df <- summary_df |>
  filter(n > 0)

p <- ggplot(summary_df, aes(x = freemix_bin, y = n, fill = status)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.35) +
  geom_text(
    data = label_df,
    aes(label = comma(n)),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.06)),
    labels = comma
  ) +
  labs(
    title = "FREEMIX Bin by Status",
    x = "FREEMIX Bin",
    y = "Count",
    fill = "Status"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = output_plot_path,
  plot = p,
  width = 8.5,
  height = 6,
  units = "in",
  device = cairo_pdf
)
