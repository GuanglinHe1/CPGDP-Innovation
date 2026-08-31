#!/usr/bin/env bash
#
# 13-sourcefind.sh
#
# SOURCEFIND ancestry decomposition. For every target individual the copying
# vector obtained from ChromoPainter is modelled as a mixture of the copying
# vectors of the surrogate populations, using the parameter file template with
# the placeholder AAAAA replaced by the target name.
#
# Requirements: R, sourcefindv2.R, ChromoPainter chunklength output
#
# Usage:
#   bash 13-sourcefind.sh <SAMPLE_LIST>

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# SOURCEFIND v2 R script; override to point at the local installation.
SOURCEFIND_R="${SOURCEFIND_R:-sourcefindv2.R}"
# Chunklength matrix produced by chromocombine over all chromosomes.
CHUNKLENGTHS="${CHUNKLENGTHS:-chromocombineALLfiles.chunklengths.out}"
# Parameter file template; the string AAAAA marks the target name.
PARAM_TEMPLATE="${PARAM_TEMPLATE:-SourcefindParamfile.txt}"
# ChromoPainter id file: <sample id> <population label> <inclusion flag>.
ID_FILE="${ID_FILE:-pops.ids}"
# File with one target individual or population per line.
SAMPLE_LIST="${1:?usage: bash 13-sourcefind.sh <SAMPLE_LIST>}"

while read -r target; do
  [[ -z "${target}" ]] && continue

  sed "s/AAAAA/${target}/g" "${PARAM_TEMPLATE}" > "SourcefindParamfile_${target}.txt"

  nohup Rscript "${SOURCEFIND_R}" \
    --chunklengths "${CHUNKLENGTHS}" \
    --parameters "SourcefindParamfile_${target}.txt" \
    --target "${target}" \
    --output "Sourcefind__${target}" \
    --idfile "${ID_FILE}" \
    > "Sourcefind__${target}.log" 2>&1 &
done < "${SAMPLE_LIST}"
wait
