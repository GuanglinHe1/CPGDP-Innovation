# ---------------------------------------------------------------------------
# 4-plot-heterozygosity.R
#
# Boxplots of the per-sample level of genetic diversity (heterozygosity),
# grouped by population and coloured by language phylum.
#
# Two figures are produced:
#   1) all populations that passed the heterozygosity filter
#   2) the same data set after removing the populations that are represented
#      by too few samples to give a stable median
#
# Input : <prefix>.het.language.txt written by
#         1-variant-summary-statistics.sh, a headerless table with the columns
#           V1 : sample id
#           V2 : heterozygosity
#           V3 : population label
#           V4 : language phylum
# Output: heterozygosity_all_populations.pdf
#         heterozygosity_selected_populations.pdf
# ---------------------------------------------------------------------------

library(ggplot2)

input_file <- "CPGDP.het.language.txt"

# Upper heterozygosity bound; samples above it are contamination outliers.
max_heterozygosity <- 0.08

# Populations dropped from the second figure.
excluded_populations <- c(
  "Miao_Guizhou", "Tanka_Fujian", "Dong_Guizhou", "Gaoshan_Fujian",
  "Han_Taiwan", "Lahu_Yunnan", "Mongolian_Fujian", "Cajia_Guizhou",
  "Manchu_Fujian", "Miao_Hunan"
)

het <- read.table(input_file, header = FALSE)

# Plot heterozygosity by population, ordered by the population median.
plot_heterozygosity <- function(dat) {
  dat$V3 <- reorder(dat$V3, dat$V2, FUN = median)
  ggplot(dat, aes(x = V3, y = V2, fill = V4)) +
    geom_boxplot() +
    labs(x = "", y = "Level of genetic diversity", fill = "Language Phyla") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, size = 15),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 16),
      axis.text.y  = element_text(size = 15),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
}

het_filtered <- het[het$V2 < max_heterozygosity, ]
ggsave(
  "heterozygosity_all_populations.pdf",
  plot_heterozygosity(het_filtered),
  width = 14, height = 8
)

het_selected <- het_filtered[!het_filtered$V3 %in% excluded_populations, ]
ggsave(
  "heterozygosity_selected_populations.pdf",
  plot_heterozygosity(het_selected),
  width = 14, height = 8
)
