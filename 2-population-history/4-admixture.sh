#!/usr/bin/env bash
#
# 4-admixture.sh
#
# Unsupervised model-based ancestry estimation with ADMIXTURE on the LD-pruned
# marker set produced by 2-ld-pruning.sh. Every ancestral cluster number K is
# run in the background with 100 bootstrap replicates and 10-fold
# cross-validation; the CV error reported in output<K> is used to pick the
# best-fitting K.
#
# Requirements: ADMIXTURE
#
# Usage:
#   bash 4-admixture.sh [PRUNED_BFILE_PREFIX] [K_MIN] [K_MAX]

set -euo pipefail

# LD-pruned PLINK binary file set prefix.
PRUNED_BFILE="${1:-${PRUNED_BFILE:-poplist_prunned}}"
# Range of ancestral cluster numbers to evaluate.
K_MIN="${2:-2}"
K_MAX="${3:-20}"
# Number of CPU threads given to a single ADMIXTURE run.
THREADS="${THREADS:-15}"

for k in $(seq "${K_MIN}" "${K_MAX}"); do
  nohup admixture -B100 --cv=10 -s time -j"${THREADS}" \
    "${PRUNED_BFILE}.bed" "${k}" > "output${k}" &
done
wait
