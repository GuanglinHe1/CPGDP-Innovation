# ---------------------------------------------------------------------------
# 3-cms-composite-score.R
#
# Composite of multiple signals (CMS) selection score. The four component
# statistics are turned into empirical p values and combined into a single
# score per variant:
#
#   iHS    two-sided empirical p value, extreme values in both directions are
#          evidence of selection
#   XP-EHH one-sided (right tail) empirical p value
#   dDAF   one-sided rank based p value
#   PBS    one-sided rank based p value
#
#   CMS = -log10(p_ihs * p_pbs * p_xpehh * p_ddaf)
#
# Ties in the dDAF and PBS statistics are broken by a tiny amount of jitter so
# that the ranks are unique; the random seeds keep this reproducible. A p
# value of exactly zero is replaced by 1/n, the smallest value the empirical
# distribution can resolve.
#
# Input : <group>.stat.ihs.pbs.xpehh.ddaf.raw.txt, a headerless table with
#           chrpos, ihs, pbs, xpehh, ddaf
#         built from the normalised outputs of 1-ihs.sh and
#         2-pbs-ddaf-xpehh.sh.
# Output: <group>.cms.txt            chrpos, CMS score and CMS p value
#         <group>.cms.top0001.txt    the top 0.1 % of the variants
#
# The <group>.cms.txt files are the input of 4-cms-clumping.sh.
# ---------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggrastr)

# Population groups to process; one input file per group.
groups <- c("South", "island", "inland")

# Empirical p value from the cumulative distribution of a statistic.
calc_empirical_p <- function(x, tail = "two") {
  ec <- ecdf(x)
  if (tail == "two") {
    p <- 2 * pmin(ec(x), 1 - ec(x))
  } else if (tail == "right") {
    p <- 1 - ec(x)
  }
  return(p)
}

for (group in groups) {
  cat("Processing:", group, "\n")

  infile <- paste0(group, ".stat.ihs.pbs.xpehh.ddaf.raw.txt")
  data <- read.table(infile, header = FALSE)
  names(data) <- c("chrpos", "ihs", "pbs", "xpehh", "ddaf")

  data$p_ihs   <- calc_empirical_p(data$ihs, "two")
  data$p_xpehh <- calc_empirical_p(data$xpehh, "right")

  set.seed(1234)
  ddaf_jitter <- data$ddaf + runif(length(data$ddaf), 0, 1e-6)
  data$p_ddaf <- 1 - rank(ddaf_jitter) / length(ddaf_jitter)

  set.seed(123)
  pbs_jitter <- data$pbs + runif(length(data$pbs), 0, 1e-6)
  data$p_pbs <- 1 - rank(pbs_jitter) / length(pbs_jitter)

  # Bound the p values away from zero before they are multiplied.
  data$p_ihs[data$p_ihs == 0]     <- 1 / nrow(data)
  data$p_xpehh[data$p_xpehh == 0] <- 1 / nrow(data)
  data$p_pbs[data$p_pbs == 0]     <- 1 / nrow(data)
  data$p_ddaf[data$p_ddaf == 0]   <- 1 / nrow(data)

  data$cms  <- -log10(data$p_ihs * data$p_pbs * data$p_xpehh * data$p_ddaf)
  data$cmsP <- data$p_ihs * data$p_pbs * data$p_xpehh * data$p_ddaf

  # Columns 1, 10 and 11 are chrpos, cms and cmsP.
  write.table(
    data[, c(1, 10, 11)],
    file      = paste0(group, ".cms.txt"),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  # Split the "<chromosome>-<position>" key into its two fields.
  tmp <- data.frame(do.call(rbind, strsplit(data$chrpos, "-")))
  colnames(tmp) <- c("chr", "pos")

  new_data <- data.frame(
    SNP        = data$chrpos,
    Chromosome = tmp$chr,
    Position   = tmp$pos,
    CMS        = data$cms
  )

  sss <- new_data[order(new_data$CMS, decreasing = TRUE), ]
  top_n <- floor(nrow(sss) / 1000)
  write.table(
    sss[1:top_n, ],
    file      = paste0(group, ".cms.top0001.txt"),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  cat("Finished:", group, "| Top 0.1%:", top_n, "\n")
}
