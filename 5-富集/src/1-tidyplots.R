library(tidyverse)
library(tidyplots)

setwd('/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/5-富集/script')

df <- read_tsv('../output/3-summary/enrichment_summary.tsv')
# 与 03_plot_enrichment_pheatmap.R 保持同一套 X 轴标签清洗规则。
df <- df |>
  mutate(
    region_label = region_name |>
      sub("\\.bed$", "", x = _) |>
      sub("^[0-9]+-", "", x = _) |>
      str_replace_all("_", " ")
  )
df |>
  tidyplot(
    x = region_label,
    y = observed_percentage,
    color = frequency_class
  ) |>
  add_mean_bar() |>
  # add_data_labels(label = round(observed_percentage, 2), label_position = "above", fontsize = 4) |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_size(width = 200, height = 100)|>
  sort_x_axis_labels(.reverse = TRUE) |>
  adjust_colors(
    new_colors = c(
      "all" = "#032a2c",
      "common" = "#274753",
      "low-frequency" = "#299d8f",
      "rare" = "#f5c710",
      "ultra-rare" = "#d55e00"
    )) |>
  save_plot('../output/4-figures/enrichment_summary-tidyplots.pdf')
