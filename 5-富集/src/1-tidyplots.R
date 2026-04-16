library(tidyverse)
library(tidyplots)

Autosome_file   <- ('/mnt/l/20-NUMTs/7-NUMTs断点富集/2-严格阈值/2-新样本/3-summary/enrichment_summary.tsv')
MtDNA_file      <- ('/mnt/l/20-NUMTs/7-NUMTs断点富集/2-严格阈值/2-新样本/3-mt-summary-fine/mt_enrichment_summary.fine.tsv')
Centromere_file  <- ('/mnt/l/20-NUMTs/7-NUMTs断点富集/2-严格阈值/2-新样本/4-figures/breakpoint_centromere_telomere_distance.tsv')

Autosome_plot   <- ('/mnt/l/20-NUMTs/7-NUMTs断点富集/2-严格阈值/2-新样本/1-tidyplots.pdf')
MtDNA_plot      <- ('/mnt/l/20-NUMTs/7-NUMTs断点富集/2-严格阈值/2-新样本/2-tidyplots_mt.pdf')
Centromere_plot  <- ('/mnt/l/20-NUMTs/7-NUMTs断点富集/2-严格阈值/2-新样本/3-tidyplots_Centromere.pdf')

df_Autosome     <- read_tsv(Autosome_file)
df_MtDNA        <- read_tsv(MtDNA_file)
df_Centromere   <- read_tsv(Centromere_file)

df_Autosome |>
  mutate(
    region_label = region_name |>
      sub("\\.bed$", "", x = _) |>
      sub("^[0-9]+-", "", x = _) |>
      str_replace_all("_", " ")
  ) |>
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
  save_plot(Autosome_plot)

df_MtDNA |>
  mutate(
    region_label = region_name |>
      sub("\\.bed$", "", x = _) |>
      sub("^[0-9]+-", "", x = _) |>
      str_replace_all("_", " ")
  ) |>
  tidyplot(
    x = region_label,
    y = observed_percentage,
    color = frequency_class
  ) |>
  add_mean_bar() |>
  # add_data_labels(label = round(observed_percentage, 2), label_position = "above", fontsize = 4) |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_size(width = 100, height = 50)|>
  sort_x_axis_labels(.reverse = TRUE) |>
  adjust_colors(
    new_colors = c(
      "all" = "#032a2c",
      "common" = "#274753",
      "low-frequency" = "#299d8f",
      "rare" = "#f5c710",
      "ultra-rare" = "#d55e00"
    )) |>
  save_plot(MtDNA_plot)

# 补充 "all" 分组（复制全部断点，赋予 all 标签）并固定分组顺序
df_Centromere <- bind_rows(
  df_Centromere |> mutate(frequency_class = "all"),  #* 全部断点复制一份作为 all 组
  df_Centromere
) |>
  mutate(frequency_class = factor(
    frequency_class,
    levels = c("all", "common", "low-frequency", "rare", "ultra-rare")  #* x 轴从左到右的顺序
  ))

df_Centromere |>
  tidyplot(
    x = frequency_class,
    y = log10_dist,
    color = frequency_class
  ) |>
  add_data_points_beeswarm(
    alpha = 0.1,
    dodge_width = 0.1
  )|>
  add_violin() |>
  adjust_colors(
    new_colors = c(
      "all" = "#032a2c",
      "common" = "#274753",
      "low-frequency" = "#299d8f",
      "rare" = "#f5c710",
      "ultra-rare" = "#d55e00"
    )) |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_x_axis_title("Frequency class") |>
  adjust_y_axis_title("Log10 (distance)") |>
  adjust_title("Distance of NUMT nuclear breakpoints to centromeres or telomeres") |>
  remove_legend() |>
  save_plot(Centromere_plot)
