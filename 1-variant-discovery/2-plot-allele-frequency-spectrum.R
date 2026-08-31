# ---------------------------------------------------------------------------
# 2-plot-allele-frequency-spectrum.R
#
# Stacked bar plot of the number of SNPs per non-reference allele frequency
# bin, split by variant type (novel / known).
#
# Input : frq.txt, a header table with the columns
#           group : frequency bin code, one of A-F (see names_map below)
#           all   : number of SNPs in the bin
#           type  : variant category used as the fill of the stacked bars
#         The table is derived from the vcftools --freq output and from the
#         singleton / doubleton lists produced by
#         1-variant-summary-statistics.sh.
# Output: allele_frequency_spectrum.pdf
# ---------------------------------------------------------------------------

library(ggplot2)
library(scales)

input_file  <- "frq.txt"
output_file <- "allele_frequency_spectrum.pdf"

frq <- read.table(input_file, header = TRUE)

# Translate the short bin codes into readable axis labels. The leading letter
# is kept so that the bins stay in the intended order on the x axis.
names_map <- c(
  "A" = "AF >= 0.05",
  "B" = "B0.01 <= AF < 0.05",
  "C" = "C0.001 <= AF < 0.01",
  "D" = "DAF < 0.001",
  "E" = "EAC = 1",
  "F" = "FAC = 2"
)
frq$group <- names_map[frq$group]

# Singletons and doubletons are drawn semi-transparent to set them apart from
# the frequency bins that are based on the allele frequency itself.
p <- ggplot(frq, aes(x = group, y = all, fill = type)) +
  geom_bar(
    stat  = "identity",
    alpha = ifelse(frq$group %in% c("EAC = 1", "FAC = 2"), 0.5, 1)
  ) +
  labs(x = "Non-reference allele frequency", y = "Number of SNPs") +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    legend.position = c(0.1, 0.9),
    legend.text  = element_text(size = 16),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  scale_y_continuous(label = scales::comma, breaks = pretty(frq$all)) +
  geom_text(
    aes(label = scales::comma(all)),
    size     = 6,
    colour   = "black",
    position = position_stack()
  )

ggsave(output_file, p, width = 10, height = 8)
