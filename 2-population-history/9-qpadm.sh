#!/usr/bin/env bash
#
# 9-qpadm.sh
#
# qpWave / qpAdm admixture modelling. For a given parameter file the script
# tests whether the target population can be modelled as a mixture of the
# listed source populations relative to a set of right (outgroup) populations,
# and reports the admixture proportions with their standard errors.
#
# Requirements: AdmixTools (qpAdm)
#
# Usage:
#   bash 9-qpadm.sh <PAR_FILE>

set -euo pipefail

# qpAdm parameter file: genotype/snp/ind files plus the left and right
# population lists of the model to test.
PAR_FILE="${1:?usage: bash 9-qpadm.sh <PAR_FILE>}"

qpAdm -p "${PAR_FILE}" > "qpAdm.${PAR_FILE}"
