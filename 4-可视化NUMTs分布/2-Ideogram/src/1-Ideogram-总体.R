library(RIdeogram)
library(dplyr)
#? 安装RIdeogram
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# if (!requireNamespace("RIdeogram", quietly = TRUE)) {
#   devtools::install_github("TickingClock1992/RIdeogram")
# }

base_dir <- "/mnt/c/Users/Administrator/Desktop"
input_dir <- file.path(base_dir, "input")
output_dir <- file.path(base_dir, "output")
gff_path <- "/mnt/c/Users/Administrator/Desktop/conf/gencode.v32.annotation.gff3.gz"
karyotype_path <- "/mnt/c/Users/Administrator/Desktop/conf/human_karyotype.txt"

grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

# karyotype
data(human_karyotype, package = "RIdeogram")
head(human_karyotype)

# gene_density
data(gene_density, package = "RIdeogram")
head(gene_density)


gff_gene_density <- tryCatch(
  GFFex(
    input = gff_path,
    karyotype = karyotype_path,
    feature = "gene",
    window = 1000000
  ),
  error = function(e) {
    message("GFFex failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(gff_gene_density)) {
  gene_density <- gff_gene_density
}

# Random_RNAs_500

data(Random_RNAs_500, package = "RIdeogram")
head(Random_RNAs_500)

# RIdeogram output

ideogram(
  karyotype = human_karyotype,
  output = file.path(output_dir, "chromosome_1.svg")
)

ideogram(
  karyotype = human_karyotype[, 1:3],
  output = file.path(output_dir, "chromosome_2.svg")
)


ideogram(
  karyotype = human_karyotype,
  overlaid = gene_density,
  output = file.path(output_dir, "chromosome_genedensity_3.svg")
)


ideogram(
  karyotype = human_karyotype,
  label = Random_RNAs_500,
  label_type = "marker",
  output = file.path(output_dir, "chromosome_label_4.svg")
)


ideogram(
  karyotype = human_karyotype,
  overlaid = gene_density,
  label = Random_RNAs_500,
  label_type = "marker",
  output = file.path(output_dir, "chromosome_5.svg")
)


ideogram(
  karyotype = human_karyotype,
  overlaid = gene_density,
  label = Random_RNAs_500,
  label_type = "marker",
  colorset1 = c("#fc8d59", "#ffffbf", "#91bfdb"),
  output = file.path(output_dir, "chromosome_6.svg")
)




filter_human_karyotype <- human_karyotype %>% filter(Chr %in% 1:10)
write.csv(filter_human_karyotype, file.path(input_dir, "filter_human_karyotype.csv"), row.names = FALSE)

ideogram(
  karyotype = filter_human_karyotype,
  overlaid = gene_density,
  label = Random_RNAs_500,
  label_type = "marker",
  output = file.path(output_dir, "chromosome_7.svg")
)


ideogram(
  karyotype = filter_human_karyotype,
  overlaid = gene_density,
  label = Random_RNAs_500,
  label_type = "marker",
  width = 100,
  output = file.path(output_dir, "chromosome_8.svg")
)


ideogram(
  karyotype = human_karyotype,
  overlaid = gene_density,
  label = Random_RNAs_500,
  label_type = "marker",
  width = 100,
  Lx = 80,
  Ly = 25,
  output = file.path(output_dir, "chromosome_9.svg")
)

# LTR_density is used as a heatmap label

data(LTR_density, package = "RIdeogram")

ideogram(
  karyotype = human_karyotype,
  overlaid = gene_density,
  label = LTR_density,
  label_type = "heatmap",
  colorset1 = c("#f7f7f7", "#e34a33"),
  colorset2 = c("#f7f7f7", "#2c7fb8"),
  output = file.path(output_dir, "chromosome_10.svg")
)

# line charts

data(liriodendron_karyotype, package = "RIdeogram")
head(liriodendron_karyotype)

data(Fst_between_CE_and_CW, package = "RIdeogram")
head(Fst_between_CE_and_CW)

data(Pi_for_CE, package = "RIdeogram")
head(Pi_for_CE)

ideogram(
  karyotype = liriodendron_karyotype,
  overlaid = Fst_between_CE_and_CW,
  label = Pi_for_CE,
  label_type = "line",
  colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"),
  output = file.path(output_dir, "chromosome_11.svg")
)

# two lines
data(Pi_for_CE_and_CW, package = "RIdeogram")

ideogram(
  karyotype = liriodendron_karyotype,
  overlaid = Fst_between_CE_and_CW,
  label = Pi_for_CE_and_CW,
  label_type = "line",
  colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"),
  output = file.path(output_dir, "chromosome_12.svg")
)


# polygon fill

ideogram(
  karyotype = liriodendron_karyotype,
  overlaid = Fst_between_CE_and_CW,
  label = Pi_for_CE_and_CW,
  label_type = "polygon",
  colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"),
  output = file.path(output_dir, "chromosome_13.svg")
)

# dual comparison (sankey)

data(karyotype_dual_comparison, package = "RIdeogram")
head(karyotype_dual_comparison)

table(karyotype_dual_comparison$species)

data(synteny_dual_comparison, package = "RIdeogram")
head(synteny_dual_comparison)

ideogram(
  karyotype = karyotype_dual_comparison,
  synteny = synteny_dual_comparison,
  output = file.path(output_dir, "sankey_chromosome.svg")
)

# ternary comparison

data(karyotype_ternary_comparison, package = "RIdeogram")
head(karyotype_ternary_comparison)
table(karyotype_ternary_comparison$species)

data(synteny_ternary_comparison, package = "RIdeogram")
head(synteny_ternary_comparison)
tail(synteny_ternary_comparison, n = 20)

ideogram(
  karyotype = karyotype_ternary_comparison,
  synteny = synteny_ternary_comparison,
  output = file.path(output_dir, "tri_chromosome.svg")
)
# ternary gradient

data(synteny_ternary_comparison_graident, package = "RIdeogram")

ideogram(
  karyotype = karyotype_ternary_comparison,
  synteny = synteny_ternary_comparison_graident,
  output = file.path(output_dir, "tri_chromosome_graident.svg")
)
