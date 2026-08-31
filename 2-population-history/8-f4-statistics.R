# ---------------------------------------------------------------------------
# 8-f4-statistics.R
#
# f4(W, X; Y, Z) statistics testing for a differential allele-sharing signal
# between the two reference populations W and X and the studied populations Y,
# with a deep African outgroup as Z. A Z score outside +/- 3 is taken as
# evidence against the null hypothesis of a symmetric relationship.
#
# Input : EIGENSTRAT genotype data loaded as `snps`, plus the vectors of
#         reference and studied population labels.
# Output: f4_results.csv
# ---------------------------------------------------------------------------

library(admixtools)

# Population labels; replace by the labels of the data set under study.
reference1   <- "Reference1"
reference2   <- "Reference2"
studied_pops <- c("StudiedPop1", "StudiedPop2")
# Deep outgroup, held fixed in every comparison.
outgroup     <- "Mbuti.DG"

output_file <- "f4_results.csv"

# EIGENSTRAT data set prefix (<prefix>.geno/.snp/.ind).
snps <- "data"

f4_result <- f4(
  W    = reference1,
  X    = reference2,
  Y    = studied_pops,
  Z    = outgroup,
  data = snps
)

write.table(f4_result, file = output_file, sep = ",", row.names = FALSE)
