#!/usr/bin/env Rscript
# 04_plot_centromere_telomere_distance.R
# 计算 NUMT nuclear breakpoint 到最近着丝粒/端粒的距离，按频率分组绘制小提琴图

library(ggplot2)

get_arg_value <- function(args, key, default = NA_character_) {
  match_idx <- match(key, args)
  if (is.na(match_idx) || match_idx == length(args)) return(default)
  args[[match_idx + 1]]
}

parse_csv <- function(value) {
  trimws(strsplit(trimws(value), ",", fixed = TRUE)[[1]])
}

args              <- commandArgs(trailingOnly = TRUE)
breakpoint_bed    <- get_arg_value(args, "--breakpoint-bed")
centromere_bed    <- get_arg_value(args, "--centromere-bed")
telomere_bed      <- get_arg_value(args, "--telomere-bed")
output_dir        <- get_arg_value(args, "--output-dir")
frequency_classes <- parse_csv(get_arg_value(args, "--frequency-classes",
                                             "all,common,low-frequency,rare,ultra-rare"))

if (is.na(breakpoint_bed) || !nzchar(breakpoint_bed)) stop("缺少参数 --breakpoint-bed")
if (is.na(centromere_bed) || !nzchar(centromere_bed)) stop("缺少参数 --centromere-bed")
if (is.na(telomere_bed)   || !nzchar(telomere_bed))   stop("缺少参数 --telomere-bed")
if (is.na(output_dir)     || !nzchar(output_dir))     stop("缺少参数 --output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#* =====读取数据=====

# 读取断点 BED 文件
dat1 <- read.delim(
  breakpoint_bed,
  header           = FALSE,
  col.names        = c("chr", "start", "end",                  #* 染色体坐标
                       "breakpoint_id", "merged_cluster_id",   #* 断点 ID
                       "breakpoint_side", "frequency_class"),  #* 方向与频率类别
  stringsAsFactors = FALSE
)
dat1$pos <- as.integer((dat1$start + dat1$end) / 2L)  #* 断点坐标 = 区间中点

# 读取着丝粒区域 BED
dat2        <- read.delim(centromere_bed, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
names(dat2) <- c("chr", "start", "end")

# 读取端粒区域 BED
dat3        <- read.delim(telomere_bed, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
names(dat3) <- c("chr", "start", "end")

# 合并着丝粒和端粒为地标区域
dat_landmarks       <- rbind(dat2, dat3)
dat_landmarks$start <- as.integer(dat_landmarks$start)
dat_landmarks$end   <- as.integer(dat_landmarks$end)

#* =====计算最近距离=====

# 按染色体预分组，加速逐行查找
landmark_by_chr <- split(dat_landmarks, dat_landmarks$chr)

# 逐行计算断点到最近地标区域的距离
# 点 p 到区间 [a, b] 的距离 = max(0, a-p, p-b)
min_dist_vec <- numeric(nrow(dat1))
for (i in seq_len(nrow(dat1))) {
  lm_i <- landmark_by_chr[[dat1$chr[[i]]]]
  if (is.null(lm_i) || nrow(lm_i) == 0L) {
    min_dist_vec[[i]] <- NA_real_
  } else {
    p_i              <- dat1$pos[[i]]
    min_dist_vec[[i]] <- min(pmax(0L, pmax(lm_i$start - p_i, p_i - lm_i$end)))
  }
}
dat1$min_dist   <- min_dist_vec
dat1$log10_dist <- log10(dat1$min_dist + 1)  #* log10(距离 + 1)，避免 log(0)

# 输出距离明细表
write.table(
  dat1[, c("chr", "pos", "breakpoint_id", "merged_cluster_id",
           "breakpoint_side", "frequency_class", "min_dist", "log10_dist")],
  file      = file.path(output_dir, "breakpoint_centromere_telomere_distance.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

#* =====绘图准备=====

# 定义频率分组的显示标签
class_labels <- c(
  "all"           = "All",            #* 全部断点
  "common"        = "Common",         #* 常见变异
  "low-frequency" = "Low frequency",  #* 低频变异
  "rare"          = "Rare",           #* 罕见变异
  "ultra-rare"    = "Ultra rare"      #* 极罕见变异
)
class_colors <- c(
  "All"           = "#032a2c",
  "Common"        = "#274753",
  "Low frequency" = "#299d8f",
  "Rare"          = "#f5c710",
  "Ultra rare"    = "#d55e00"
)

# 确定分组顺序
# "all" 为虚拟分组（合并所有断点），不在 all_records.bed 的 frequency_class 列中
real_classes <- frequency_classes[frequency_classes != "all" &
                                  frequency_classes %in% unique(dat1$frequency_class)]
has_all      <- "all" %in% frequency_classes
freq_order   <- frequency_classes[frequency_classes %in% c(if (has_all) "all", real_classes)]
label_order  <- class_labels[freq_order]
label_order  <- label_order[!is.na(label_order)]

# 为 "all" 虚拟分组复制一份数据并赋标签
dat_valid <- dat1[!is.na(dat1$log10_dist), ]
if (has_all) {
  dat_all             <- dat_valid
  dat_all$frequency_class <- "all"
} else {
  dat_all <- dat_valid[0, ]
}
dat_real <- dat_valid[dat_valid$frequency_class %in% real_classes, ]
dat_plot  <- rbind(dat_all, dat_real)

dat_plot$group_label <- class_labels[dat_plot$frequency_class]
dat_plot$group_label[is.na(dat_plot$group_label)] <-
  dat_plot$frequency_class[is.na(dat_plot$group_label)]
dat_plot$group_label <- factor(dat_plot$group_label, levels = label_order)

# 输出分组汇总统计表
dat_summary <- do.call(rbind, lapply(split(dat_plot, dat_plot$group_label), function(g) {
  data.frame(
    frequency_class   = g$frequency_class[[1]],
    group_label       = as.character(g$group_label[[1]]),
    n                 = nrow(g),
    median_log10_dist = median(g$log10_dist,               na.rm = TRUE),
    q25_log10_dist    = quantile(g$log10_dist, 0.25,       na.rm = TRUE),
    q75_log10_dist    = quantile(g$log10_dist, 0.75,       na.rm = TRUE),
    median_dist_bp    = median(g$min_dist,                 na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
}))
write.table(
  dat_summary,
  file      = file.path(output_dir, "breakpoint_centromere_telomere_distance_summary.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

#* =====绘制小提琴图=====

# 基础映射
p <- ggplot(
  dat_plot,
  aes(
    x    = group_label,  #* x 轴为频率分组
    y    = log10_dist,   #* y 轴为 log10(距离+1)
    fill = group_label   #* 按分组填色
  )
)

# 叠加小提琴图层
p2 <- p +
  geom_violin(
    trim  = FALSE,  #* 保留分布完整尾部
    alpha = 0.7,    #* 半透明
    color = NA      #* 不显示边框
  )

# 叠加内嵌箱线图
p3 <- p2 +
  geom_boxplot(
    width        = 0.12,    #* 箱体宽度
    outlier.size = 0.6,     #* 离群点大小
    alpha        = 0.9,
    color        = "grey30",
    fill         = "white"  #* 白色箱体以凸显中位数线
  )

# 设置颜色与主题
p4 <- p3 +
  scale_fill_manual(
    values   = class_colors[levels(dat_plot$group_label)],  #* 按顺序取预设颜色
    na.value = "grey70"
  ) +
  labs(
    x     = NULL,
    y     = expression(log[10](distance + 1)),  #* y 轴数学表达式
    title = "Distance of NUMT nuclear breakpoints to\ncentromeres and telomeres"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 11, face = "bold"),
    axis.text.x     = element_text(size = 10)
  )

# 保存 PDF
cairo_pdf(
  file.path(output_dir, "breakpoint_centromere_telomere_distance.pdf"),
  width  = 6,       #* 宽度（英寸）
  height = 5,       #* 高度（英寸）
  family = "sans"
)
print(p4)
dev.off()

# 保存 PNG
png(
  file.path(output_dir, "breakpoint_centromere_telomere_distance.png"),
  width  = 6,
  height = 5,
  units  = "in",
  res    = 300,     #* 分辨率 300 dpi
  type   = "cairo"
)
print(p4)
dev.off()

# 保存 SVG
svg(
  file.path(output_dir, "breakpoint_centromere_telomere_distance.svg"),
  width  = 6,
  height = 5
)
print(p4)
dev.off()

message("完成。输出目录：", output_dir)
