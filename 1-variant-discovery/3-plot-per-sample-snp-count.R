# ---------------------------------------------------------------------------
# 3-plot-per-sample-snp-count.R
#
# Boxplot of the number of SNPs carried by every sample, grouped by population
# and coloured by language phylum.
#
# Input : 27pops_each_sm_snps.txt, a headerless table with the columns
#           V1 : population label
#           V2 : number of SNPs of the sample
#           V3 : language phylum
# Output: per_sample_snp_count.pdf
# ---------------------------------------------------------------------------

library(ggplot2)
library(scales)

input_file  <- "27pops_each_sm_snps.txt"
output_file <- "per_sample_snp_count.pdf"

# Minimum number of SNPs a sample must carry to enter the plot; samples below
# this threshold are low-coverage outliers.
min_snp_count <- 1750000

snp_count <- read.table(input_file, header = FALSE)

# Order the populations by their median SNP count.
snp_count$V1 <- reorder(snp_count$V1, snp_count$V2, FUN = median)

snp_count_filtered <- snp_count[snp_count$V2 > min_snp_count, ]

p <- ggplot(snp_count_filtered, aes(x = V1, y = V2, fill = V3)) +
  geom_boxplot() +
  labs(x = "", y = "Number of SNPs", fill = "Language Phyla") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, size = 15),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  scale_y_continuous(
    label  = scales::comma,
    breaks = pretty(snp_count_filtered$V2)
  )

ggsave(output_file, p, width = 12, height = 8)
