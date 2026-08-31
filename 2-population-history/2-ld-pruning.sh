#!/usr/bin/env bash
#
# 2-ld-pruning.sh
#
# Linkage-disequilibrium pruning of the merged PLINK data set. The pruned
# marker set is the input of the UMAP embedding (3-umap.R) and of the
# ADMIXTURE analysis (4-admixture.sh), which both assume independent markers.
#
# Requirements: plink2, plink
#
# Usage:
#   bash 2-ld-pruning.sh [BFILE_PREFIX]

set -euo pipefail

# PLINK binary file set prefix (<prefix>.bed/.bim/.fam).
BFILE="${1:-${BFILE:-poplist}}"

# --- Marker set used for the UMAP embedding -------------------------------
# Window of 50 markers, step of 5 markers, r^2 threshold of 0.2.
plink2 --bfile "${BFILE}" \
  --indep-pairwise 50 5 0.2 \
  --out pruned_data

# --- Marker set used for ADMIXTURE ----------------------------------------
# Window of 200 markers, step of 25 markers, r^2 threshold of 0.4.
plink --bfile "${BFILE}" \
  --indep-pairwise 200 25 0.4 \
  --allow-no-sex \
  --out plink

plink --bfile "${BFILE}" \
  --extract plink.prune.in \
  --make-bed \
  --allow-no-sex \
  --out "${BFILE}_prunned"
