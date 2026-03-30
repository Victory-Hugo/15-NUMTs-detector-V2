library(tidyverse)
library(tidyplots)

df_length <- read_tsv("/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-统计性描述/output/2-numt-length-by-region.tsv") 

# 统计长度分布
df_length |> select(chosen_length) |> summary()

