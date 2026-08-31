#!/usr/bin/env bash
#
# 2-ibdmix.sh
#
# Detection of Neanderthal introgressed segments with IBDmix, which calls
# segments that are identical by descent between a modern sample and the
# archaic reference genome without relying on a modern reference panel.
#
#   generate_gt  merge the archaic and the modern genotypes into the IBDmix
#                input format
#   ibdmix       call the introgressed segments per sample
#   summary.sh   filter the calls by segment length and LOD score and
#                summarise them per sample
#
# Requirements: IBDmix (generate_gt, ibdmix), summary.sh
#
# Usage:
#   bash 2-ibdmix.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Directory with the archaic reference VCFs and the exclusion masks.
ARCHAIC_REF_DIR="${ARCHAIC_REF_DIR:-./resources/archaic_public_ref/Neanderthal}"
# Helper script shipped with IBDmix; override to give a full path.
SUMMARY_SH="${SUMMARY_SH:-summary.sh}"
# Sample list, one modern sample id per line.
SAMPLE_LIST="${SAMPLE_LIST:-sample.list}"
# Minimum segment length in bp and minimum LOD score of a reported segment.
MIN_LENGTH="${MIN_LENGTH:-50000}"
MIN_LOD="${MIN_LOD:-4}"

for chr in {1..22}; do
  generate_gt \
    -a "${ARCHAIC_REF_DIR}/output_se.hg38.chr${chr}.vcf" \
    -m "chr${chr}.cl.snp.vcf" \
    -o "neaderthal.${chr}"

  ibdmix \
    -g "neaderthal.${chr}" \
    -r "${ARCHAIC_REF_DIR}/GRCh38_Neanderthal_exclude.chr${chr}.chr.bed" \
    -w \
    -s "${SAMPLE_LIST}" \
    -o "AltaiNea.${chr}.out"

  bash "${SUMMARY_SH}" "${MIN_LENGTH}" "${MIN_LOD}" all \
    "AltaiNea.${chr}.out" "all.AltaiNea.${chr}.summary.txt"
done
