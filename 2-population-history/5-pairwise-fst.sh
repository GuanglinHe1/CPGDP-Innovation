#!/usr/bin/env bash
#
# 5-pairwise-fst.sh
#
# Pairwise Fst between all population pairs of the merged data set. The helper
# script pairwise.perl wraps smartpca in fst-only mode and writes the full
# pairwise Fst matrix together with its standard errors.
#
# Requirements: perl, EIGENSOFT (smartpca), pairwise.perl
#
# Usage:
#   bash 5-pairwise-fst.sh [BFILE_PREFIX]

set -euo pipefail

# EIGENSTRAT/PLINK data set prefix holding the populations to compare.
BFILE="${1:-${BFILE:-poplist}}"
# Helper script shipped with this repository; override to give a full path.
PAIRWISE_PERL="${PAIRWISE_PERL:-pairwise.perl}"

# The two trailing arguments are the first and the last population index to
# process; 0 0 runs all population pairs.
perl "${PAIRWISE_PERL}" "${BFILE}" 0 0
