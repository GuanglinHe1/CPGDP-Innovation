#!/usr/bin/env bash
#
# 15-msmc2.sh
#
# Cross-coalescence rate between two populations with MSMC2. Four haplotypes
# per population are used: the first two runs estimate the coalescence rate
# within each population, the third run estimates the cross-population rate.
# The relative cross-coalescence rate derived from the three runs dates the
# separation of the two populations.
#
# Requirements: MSMC2, multihetsep input files
#
# Usage:
#   bash 15-msmc2.sh <MULTIHETSEP_GLOB>

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Per-chromosome multihetsep files of the two populations, e.g. "input/*.txt".
MULTIHETSEP="${1:?usage: bash 15-msmc2.sh <MULTIHETSEP_GLOB>}"
# Number of CPU threads per MSMC2 run.
THREADS="${THREADS:-11}"
# Haplotype indices: 0-3 belong to population A, 4-7 to population B.
POP_A="0,1,2,3"
POP_B="4,5,6,7"
# All cross-population haplotype pairs.
POP_AB="0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7"

# Within population A.
msmc2 -t "${THREADS}" -s -I "${POP_A}" -o AB ${MULTIHETSEP} &
# Within population B.
msmc2 -t "${THREADS}" -s -I "${POP_B}" -o BA ${MULTIHETSEP} &
# Between populations A and B.
msmc2 -t "${THREADS}" -s -I "${POP_AB}" -o A_B ${MULTIHETSEP} &
wait
