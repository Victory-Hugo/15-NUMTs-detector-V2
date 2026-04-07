library(tidyverse)
library(tidyplots)

setwd("/mnt/l/20-NUMTs/2-质控指标/script/")

df_PC <- read_csv("../img/data/ancestry_contam_vs_intended_scatter.csv")

p1 <- df_PC |> 
	tidyplot(
		x = ContaminatingSample,
		y = IntendedSample,
		color = IntendedSample
	) |>
	add_data_points(rasterize = TRUE, rasterize_dpi = 600) |>
	adjust_title("")  |>
	adjust_caption("ancestry_contam_vs_intended_scatter") |>
	adjust_colors(colors_continuous_plasma) |>
	adjust_legend_position("none") |>
	split_plot(by = PC) 

df_freemix <- read_csv("../img/data/verifybam_depth_vs_freemix.csv")
p2 <- df_freemix |>
	tidyplot(
		x = AVG_DP,
		y = FREEMIX,
		color = status
	) |>
	add_data_points(alpha = 0.5, rasterize = TRUE, rasterize_dpi = 600) |>
	adjust_y_axis(transform = "log10") |>
	adjust_colors(new_colors = c(
    "FAIL" = "#d55e00",
    "PASS" = "#0072b2")) |>
	adjust_title("The relationship between FREEMIX and average depth") 

df_copynumber_depth <- read_csv("../img/data/copynumber_depth_scatter.csv")
p3 <- df_copynumber_depth |>
	tidyplot(
		x = mean_autosomal_depth,
		y = mean_mtDNA_depth,
		color = 质控状态
	) |>
	add_data_points(alpha = 0.5, rasterize = TRUE, rasterize_dpi = 600) |>
	adjust_colors(new_colors = c(
		"FAIL" = "#d55e00",
		"PASS" = "#0072b2")) |>
	adjust_title("The mean autosomal depth and mean mtDNA depth")

df_copynumber_distribution <- read_csv("../img/data/copynumber_distribution_by_qc.csv")
p4 <- df_copynumber_distribution |>
	tidyplot(
		x = 质控状态,
		y = mtDNA_copy_number,
		color = 质控状态
	) |>
	add_data_points_beeswarm( 
		alpha = 0.010,
  	jitter_width = 0.001,
		rasterize = TRUE,
		rasterize_dpi = 600) |>
	adjust_colors(new_colors = c(
		"FAIL" = "#d55e00",
		"PASS" = "#0072b2")) |>
	add_boxplot() |>
	adjust_title("The distribution of mtDNA copy number by QC status")

df_freemix_bin <- read_csv("../img/data/verifybam_freemix_bin_status.csv") |>
	mutate(
		freemix_bin = factor(
			freemix_bin,
			levels = c("<=0.001", "0.001-0.005", "0.005-0.01", "0.01-0.02", "0.02-0.05", ">0.05")
		)
	)
p5.1 <- df_freemix_bin |>
	tidyplot(
		x = freemix_bin,
		y = n,
		color = freemix_bin
	) |>
	add_mean_bar() |>
	adjust_x_axis(rotate_labels = 90 ) 

p5.2 <- df_freemix_bin |>
	tidyplot(
		x = freemix_bin,
		y = n,
		color = freemix_bin
	) |>
	add_mean_bar() |>
	add_data_labels(label = n) |>
	adjust_x_axis(rotate_labels = 90 ) 

library(patchwork)
#! 使用patchwork v1.3.0组合tidyplot对象需要使用如下代码
p_all <- patchwork::wrap_plots(p1, p2, p3, p4, p5.1, p5.2) + 
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave(filename = "../img/plots/qc_metrics_scatter.pdf", plot = p_all, width = 40, height = 50, units = "cm")
