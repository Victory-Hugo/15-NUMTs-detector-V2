library(tidyverse)
library(tidyplots)

setwd('/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-汇总NUMTs分布/script')
df_allCluster_sum <- read_tsv('../output/merge.tsv.allCluster.manhattan.points.tsv')

chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
df_allCluster_sum <- df_allCluster_sum |>
  mutate(CHR = factor(CHR, levels = chr_levels, ordered = TRUE))

chr_labels <-
  df_allCluster_sum |>
  summarise(center = mean(x), .by = CHR) |>
  arrange(CHR)

p1 <- df_allCluster_sum |>
  tidyplot(x = x, y = POS, color = CHR) |>
  add_data_points(rasterize = TRUE, alpha = 0.5, size = 0.5) |>
  adjust_size(150, 70) |>
  adjust_colors(new_colors = c("#d43e3e", "#dd6a3b", "#e58f38", "#eab236",
                            "#e8cc36", "#d9df3a", "#b8e446", "#8fe051",
                            "#62d76a", "#3dcd8b", "#37c3b3", "#3fb6d4",
                            "#4ca1e8", "#5a86ee", "#6769ec", "#7a4fe5",
                            "#9440d9", "#b139c7", "#cb3ab0", "#d94795",
                            "#dd5d7e", "#d66a6a","#424242","#1b1a1a"))  |>
  adjust_padding(bottom = -0, left = 0) |>
  adjust_x_axis(rotate_labels = 90,label = chr_labels$CHR, breaks= chr_labels$center) |>
  adjust_x_axis_title("Chromosome") |>
  adjust_y_axis(limits = c(0, 1000)) |>
  adjust_y_axis_title("Cluster size (unique positions)") |>
  adjust_title("NUMTs cluster distribution") |>
  remove_legend() 



df_allCluster_sample_count <- read_tsv('../output/merge.tsv.allCluster.manhattan.samplecount.points.tsv')
df_allCluster_sample_count <- df_allCluster_sample_count |>
  mutate(CHR = factor(CHR, levels = chr_levels, ordered = TRUE))

p2 <- df_allCluster_sample_count |>
  tidyplot(x = x, y = sample_count, color = CHR) |>
  add_data_points(rasterize = TRUE,alpha = 0.5, size = 0.5) |>
  adjust_size(150, 70) |>
  adjust_colors(new_colors = c("#d43e3e", "#dd6a3b", "#e58f38", "#eab236",
                            "#e8cc36", "#d9df3a", "#b8e446", "#8fe051",
                            "#62d76a", "#3dcd8b", "#37c3b3", "#3fb6d4",
                            "#4ca1e8", "#5a86ee", "#6769ec", "#7a4fe5",
                            "#9440d9", "#b139c7", "#cb3ab0", "#d94795",
                            "#dd5d7e", "#d66a6a","#424242","#1b1a1a"))  |>
  adjust_padding(bottom = -0, left = 0) |>
  adjust_x_axis(rotate_labels = 90,label = chr_labels$CHR, breaks= chr_labels$center) |>
  adjust_x_axis_title("Chromosome") |>
  adjust_y_axis(limits = c(0, 1000)) |>
  adjust_y_axis_title("Sample count") |>
  adjust_title("NUMTs cluster sample count distribution") |>
  remove_legend() 




#! 使用patchwork v1.3.0组合tidyplot对象需要使用如下代码
patchwork::wrap_plots(p1, p2) + 
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A") -> p_combined

ggsave(filename = "../output/NUMTs_cluster_distribution.pdf")
