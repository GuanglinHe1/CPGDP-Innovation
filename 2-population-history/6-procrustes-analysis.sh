#!/usr/bin/env bash
#
# 6-procrustes-analysis.sh
#
# Procrustes analysis comparing the genetic PCA space with the geographic
# coordinates of the sampling locations. Calculate_T_stat.py rotates, scales
# and reflects the PCA coordinates onto the geographic ones and reports the
# Procrustes similarity statistic t0 together with its empirical p value.
#
# Requirements: python3, Calculate_T_stat.py
#
# Usage:
#   bash 6-procrustes-analysis.sh

set -euo pipefail

# Helper script shipped with this repository; override to give a full path.
PROCRUSTES_SCRIPT="${PROCRUSTES_SCRIPT:-Calculate_T_stat.py}"

python3 "${PROCRUSTES_SCRIPT}"
