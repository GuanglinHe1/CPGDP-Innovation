#!/usr/bin/env bash
#
# 7-outgroup-f3.sh
#
# Outgroup-f3 statistics f3(Outgroup; Pop1, Pop2) measuring the amount of
# shared genetic drift between all pairs of the listed populations.
#
# Requirements: R, admixr/admixtools, allF3pairs.R
#
# Usage:
#   bash 7-outgroup-f3.sh [DATA_PREFIX] [POP_LIST] [OUTGROUP_LIST]

set -euo pipefail

# EIGENSTRAT data set prefix.
DATA_PREFIX="${1:-${DATA_PREFIX:-data}}"
# One population per line; f3 is computed for every pair of these populations.
POP_LIST="${2:-${POP_LIST:-OutgroupPops.txt}}"
# One outgroup per line, used as the reference population of the f3 statistic.
OUTGROUP_LIST="${3:-${OUTGROUP_LIST:-Outgroups.txt}}"
# Helper script shipped with this repository; override to give a full path.
F3_SCRIPT="${F3_SCRIPT:-allF3pairs.R}"

# Arguments of allF3pairs.R:
#   <data prefix> <population list> <outgroup list> <inbreed flag> <mode>
Rscript "${F3_SCRIPT}" "${DATA_PREFIX}" "${POP_LIST}" "${OUTGROUP_LIST}" yes list
