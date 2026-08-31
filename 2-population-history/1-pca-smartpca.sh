#!/usr/bin/env bash
#
# 1-pca-smartpca.sh
#
# Principal component analysis of the merged genotype data set with smartpca
# (EIGENSOFT). Modern reference populations are used to build the components
# and the studied populations are projected onto them, following the settings
# given in the smartpca parameter file.
#
# Requirements: EIGENSOFT (smartpca)
#
# Usage:
#   bash 1-pca-smartpca.sh [PAR_FILE]

set -euo pipefail

# smartpca parameter file listing the EIGENSTRAT input triplet, the output
# eigenvector/eigenvalue files and the populations used for the projection.
PAR_FILE="${1:-${PAR_FILE:-parSmartPCA.HOmerge.poplist}}"

smartpca -p "${PAR_FILE}"
