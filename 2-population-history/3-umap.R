# ---------------------------------------------------------------------------
# 3-umap.R
#
# UMAP embedding of the genetic structure. The embedding is computed on the
# leading principal components of the LD-pruned marker set (2-ld-pruning.sh),
# which keeps the input low-dimensional and removes the noise of the singleton
# markers.
#
# Input : eigenvector table with one row per sample and the columns
#           sample id, PC1 ... PCn, population label
#         (smartpca .evec output of 1-pca-smartpca.sh, or the .eigenvec file
#          of `plink --pca` run on pruned_data.prune.in)
# Output: umap.txt  two-dimensional embedding of every sample
#         umap.pdf  scatter plot of the embedding coloured by population
# ---------------------------------------------------------------------------

library(data.table)
library(uwot)
library(ggplot2)

input_file    <- "pruned_data.eigenvec"
output_table  <- "umap.txt"
output_figure <- "umap.pdf"

# Number of leading principal components fed into UMAP.
n_pc <- 10
# UMAP hyper-parameters: neighbourhood size and minimum embedding distance.
n_neighbors <- 15
min_dist    <- 0.1

set.seed(1234)

evec <- fread(input_file, header = FALSE, data.table = FALSE)

# First column: sample id; last column: population label; in between the PCs.
sample_id  <- evec[[1]]
population <- evec[[ncol(evec)]]
pcs        <- as.matrix(evec[, 2:(1 + n_pc)])

embedding <- umap(
  pcs,
  n_neighbors = n_neighbors,
  min_dist    = min_dist,
  n_components = 2
)

umap_df <- data.frame(
  sample     = sample_id,
  population = population,
  UMAP1      = embedding[, 1],
  UMAP2      = embedding[, 2]
)

write.table(umap_df, output_table, sep = "\t", row.names = FALSE, quote = FALSE)

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = population)) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(x = "UMAP1", y = "UMAP2", colour = "Population") +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    legend.text  = element_text(size = 10),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 14)
  )

ggsave(output_figure, p, width = 10, height = 8)
