#!/usr/bin/env bash
#
# 14-smcpp.sh
#
# Effective population size history of single populations with SMC++.
#
#   vcf2smc   convert the VCF of one population and one chromosome into the
#             SMC++ input format, with one distinguished individual per file
#   estimate  fit the size history under a cubic spline and the given
#             per-generation mutation rate
#   plot      draw the size histories on a linear time axis
#
# Requirements: SMC++
#
# Usage:
#   bash 14-smcpp.sh <VCF> <POP_LIST>

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Bgzipped and indexed VCF holding all populations.
VCF="${1:?usage: bash 14-smcpp.sh <VCF> <POP_LIST>}"
# File with one population label per line.
POP_LIST="${2:?usage: bash 14-smcpp.sh <VCF> <POP_LIST>}"
# Directory the SMC++ input files are written to.
OUT_DIR="${OUT_DIR:-./smcpp_input}"
# Per-generation per-site mutation rate and generation time in years.
MUTATION_RATE="${MUTATION_RATE:-1.25e-8}"
GENERATION_TIME="${GENERATION_TIME:-29}"

mkdir -p "${OUT_DIR}"

# Comma-separated list of the samples of every population, in the
# "<population>:<sample1>,<sample2>,..." form expected by vcf2smc.
while read -r pop_spec; do
  [[ -z "${pop_spec}" ]] && continue
  pop="${pop_spec%%:*}"

  for chr in {1..22}; do
    smc++ vcf2smc -d "${pop}" "${pop}" \
      "${VCF}" \
      "${OUT_DIR}/${pop}_chr${chr}.txt" \
      "${chr}" \
      "${pop_spec}"
  done

  smc++ estimate -o "./${pop}" --base "${pop}" --spline cubic \
    "${MUTATION_RATE}" "${OUT_DIR}/${pop}"_chr*.txt

  smc++ plot "${pop}.pdf" -g "${GENERATION_TIME}" -c --linear "./${pop}"/*.final.json
done < "${POP_LIST}"
