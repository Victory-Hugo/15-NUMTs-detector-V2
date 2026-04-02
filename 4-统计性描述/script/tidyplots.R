library(tidyverse)
library(tidyplots)

#! 设置输入输出文件路径
BASE_INPUT_DIR <- "/mnt/l/20-NUMTs/6-NUMTs频率分布描述/2-严格阈值/output/排除参考NUMTs严格阈值/"
df_1_file <- paste0(BASE_INPUT_DIR, "02-tables/2-numt-length-by-region.tsv")
df_2_file <- paste0(BASE_INPUT_DIR, "02-tables/4-numt-support-summary.tsv")
df_3_file <- paste0(BASE_INPUT_DIR, "02-tables/4-numt-frequency-class-summary.min-support-1.tsv")
df_4_file <- paste0(BASE_INPUT_DIR, "02-tables/5-numt-size-vs-frequency-scatter.tsv")
df_5_file <- paste0(BASE_INPUT_DIR, "02-tables/5-numt-relative-frequency-percentage.tsv")


df_1_file_out <- paste0(BASE_INPUT_DIR, "03-figures/r/tidyplots_1.pdf")
df_2_file_out <- paste0(BASE_INPUT_DIR, "03-figures/r/tidyplots_2.pdf")
df_3_file_out <- paste0(BASE_INPUT_DIR, "03-figures/r/tidyplots_3.pdf")
df_4_file_out <- paste0(BASE_INPUT_DIR, "03-figures/r/tidyplots_4.pdf")
df_5_file_out <- paste0(BASE_INPUT_DIR, "03-figures/r/tidyplots_5.pdf")
df_6_file_out <- paste0(BASE_INPUT_DIR, "03-figures/r/tidyplots_6.pdf")



# 第一部分

df_length <- read_tsv(df_1_file)

p1 <- df_length |>
  tidyplot(x = chosen_length) |>
  add(
    ggplot2::geom_histogram(
      ggplot2::aes(y = after_stat(count / sum(count))),
      binwidth = 1
    )
  ) |>
  adjust_y_axis(labels = scales::label_percent(accuracy = 1)) |>
  adjust_size(width = 100, height = 100) |>
  adjust_title("NUMT length distribution") |>
  adjust_y_axis_title("Frequency") |>
  adjust_x_axis_title("Length (bp)") |>
  adjust_colors(new_colors ="#0e4c70") |>
  save_plot(df_1_file_out)

p2 <- df_length |> filter(chosen_length <= 750) |>
  tidyplot(x = chosen_length) |>  
  add(
    ggplot2::geom_histogram(
      ggplot2::aes(y = after_stat(count / sum(count))),
      binwidth = 1
    )
  ) |>
  adjust_y_axis(labels = scales::label_percent(accuracy = 1)) |>
  adjust_size(width = 100, height = 100) |>
  adjust_title("NUMT length distribution (<= 750 bp)") |>
  adjust_y_axis_title("Frequency") |>
  adjust_x_axis_title("Length (bp)") |>
  adjust_colors(new_colors ="#0e4c70") |>
  save_plot(df_2_file_out)



df_Frequency_spectrum <- read_tsv(df_3_file)

# p3 <- df_Frequency_spectrum |> 
#   mutate(percentage = cluster_count / sum(cluster_count) * 100) |> # 计算每个类别的%占比
#   tidyplot(
#     y = percentage,
#     color = category
#   ) |>
#   add_donut() |>
#   adjust_colors(
#     new_colors = c(
#     "common" = "#274753",
#     "low-frequency" = "#299d8f",
#     "rare" = "#f5c710",
#     "ultra-rare" = "#d55e00")
#   ) |>
#   adjust_title("Distinct NUMTs frequency spectrum") |>
#   adjust_caption("Common: >= 0.05%; \n
#                             Low-frequency: 0.05% > x >= 0.01%; \n
#                             Rare: 0.01% > x >= 0.001%; \n
#                             Ultra-rare: < 0.001%") |>
#   save_plot(df_3_file_out)


p3 <- df_Frequency_spectrum |>
  mutate(
    percentage = cluster_count / sum(cluster_count) * 100,
    label = sprintf("%.1f%%", percentage)
  ) |>
  tidyplot(
    y = percentage,
    color = category
  ) |>
  add_donut() |>
  adjust_colors(
    new_colors = c(
      "common" = "#274753",
      "low-frequency" = "#299d8f",
      "rare" = "#f5c710",
      "ultra-rare" = "#d55e00"
    )
  ) |>
  adjust_title("Distinct NUMTs frequency spectrum") |>
  adjust_caption(
    "Common: >= 0.05%;\nLow-frequency: 0.05% > x >= 0.01%;\nRare: 0.01% > x >= 0.001%;\nUltra-rare: < 0.001%"
  )

label_df <- ggplot_build(p3)$data[[1]] |>
  mutate(
    x = (xmin + xmax) / 2,
    y = (ymin + ymax) / 2,
    label = df_Frequency_spectrum |>
      mutate(percentage = cluster_count / sum(cluster_count) * 100,
             label = sprintf("%.1f%%", percentage)) |>
      pull(label)
  )

p3_labeled <- p3 +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  )

ggsave(
  filename = df_3_file_out,
  plot = p3_labeled
)

p4 <- df_Frequency_spectrum |>
  tidyplot(
    x = category,
    y = cluster_count,
    color = category
  ) |>
  add_mean_bar() |>
  adjust_colors(
    new_colors = c(
      "common" = "#274753",
      "low-frequency" = "#299d8f",
      "rare" = "#f5c710",
      "ultra-rare" = "#d55e00")
  ) |>
  add_data_labels(label = cluster_count, label_position = "above") |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_title("Distinct NUMTs frequency spectrum") |>
  adjust_x_axis_title("Frequency category") |>
  sort_x_axis_labels(.reverse = TRUE) |>
  adjust_y_axis_title("Count") |>
  save_plot(df_4_file_out)

df_Frequency <- read_tsv(df_2_file)
p5 <- df_Frequency |>
  filter(min_sample_support == 1 ) |> 
  pivot_longer(
    cols = -min_sample_support, 
    names_to = "Sample", 
    values_to = "Frequency") |>
  filter(Sample == "singleton_clusters"| 
           Sample == "doubleton_clusters"|
           Sample == "tripleton_clusters")  |>
  tidyplot(
    x = Sample,
    y = Frequency,
    color = Sample
  ) |>
  add_mean_bar() |>
  sort_x_axis_labels(.reverse = TRUE) |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_colors(
    new_colors = c(
      "singleton_clusters" = "#01275c",
      "doubleton_clusters" = "#1c6e97",
      "tripleton_clusters" = "#1c9389")
  ) |>
  add_data_labels(label = Frequency, label_position = "above") |>
  adjust_title("Distinct NUMTs frequency spectrum (Ultra-rare)") |>
  adjust_x_axis_title("Frequency category") |>
  adjust_y_axis_title("Count") |>
  save_plot(df_5_file_out)

df_Freq_Size <- read_tsv(df_4_file)




df_Freq_Relative <- read_tsv(df_5_file, show_col_types = FALSE) |>
  mutate(
    frequency_category = factor(
      frequency_category,
      levels = c("common", "low-frequency", "rare", "ultra-rare", "private")
    ),
    frequency_category_label = factor(
      frequency_category_label,
      levels = c("Common", "Low-frequency", "Rare", "Ultra-rare", "Private")
    ),
    frequency_category_label_rev = fct_rev(frequency_category_label),
    label_percent = sprintf("%.2f%%", relative_percentage_vs_common)
  )

xmin_positive <- min(
  df_Freq_Relative$relative_percentage_vs_common[
    df_Freq_Relative$relative_percentage_vs_common > 0
  ])
xmin_plot <- xmin_positive / 1.6

df_Freq_Relative <- df_Freq_Relative |>
  mutate(
    text_x = case_when(
      relative_percentage_vs_common >= 85 ~ relative_percentage_vs_common / 1.15,
      TRUE ~ pmin(relative_percentage_vs_common * 1.15, 132)
    ),
    text_hjust = case_when(
      relative_percentage_vs_common >= 85 ~ 1,
      TRUE ~ 0
    ))

p7 <- df_Freq_Relative |>
  tidyplot(
    y = frequency_category_label_rev,
    x = relative_percentage_vs_common,
    color = frequency_category
  ) |>
  add(
    ggplot2::geom_segment(
      aes(
        x = xmin_plot,
        xend = relative_percentage_vs_common,
        yend = frequency_category_label_rev
      ),
      linewidth = 6,
      alpha = 0.7
    )) |>
  add(
    ggplot2::geom_point(size = 2, alpha = 1)
  ) |>
  add(
    ggplot2::geom_text(
      aes(x = text_x, label = label_percent, hjust = text_hjust),
      size = 3,
      alpha = 1
    )) |>
  adjust_x_axis(
    transform = "log10",
    limits = c(xmin_plot, 140)) |>
  adjust_colors(
    new_colors = c(
      "common" = "#274753",
      "low-frequency" = "#299d8f",
      "rare" = "#f5c710",
      "ultra-rare" = "#d55e00",
      "private" = "#1c6e97"
    )) |>
  save_plot(df_6_file_out)
