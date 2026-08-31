#!/usr/bin/env bash
#
# 10-phasing-and-ibd.sh
#
# Chromosome-wise statistical phasing with SHAPEIT and detection of identity-
# by-descent segments with Refined IBD. The phased haplotypes are also the
# input of the ChromoPainter based analyses (11 to 13).
#
# Requirements: SHAPEIT, java, refined-ibd
#
# Usage:
#   bash 10-phasing-and-ibd.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Prefix of the per-chromosome input files.
PREFIX="${PREFIX:-polist_ibd}"
# Directory holding the genetic maps used by SHAPEIT
# (genetic_map_adjust_chr<N>.txt) and by Refined IBD
# (plink.chr<N>.GRCh37.map).
MAP_DIR="${MAP_DIR:-./maps}"
# Refined IBD jar file.
REFINED_IBD_JAR="${REFINED_IBD_JAR:-./refined-ibd.17Jan20.102.jar}"
# Resources given to a single job.
THREADS="${THREADS:-32}"
JAVA_MEM="${JAVA_MEM:-60g}"
# Minimum length in cM of a reported IBD segment.
MIN_IBD_LENGTH="${MIN_IBD_LENGTH:-0.1}"

for chr in {1..22}; do
  # --- Statistical phasing ------------------------------------------------
  nohup shapeit \
    --input-bed "${PREFIX}.chr${chr}.bed" "${PREFIX}.chr${chr}.bim" "${PREFIX}.chr${chr}.fam" \
    --input-map "${MAP_DIR}/genetic_map_adjust_chr${chr}.txt" \
    --output-max "${PREFIX}.chr${chr}.phased.haps" "${PREFIX}.chr${chr}.phased.sample" \
    --output-log "${PREFIX}.chr${chr}.log" \
    --force --burn 10 --prune 10 --main 30 --thread "${THREADS}" &

  # --- IBD segment detection ----------------------------------------------
  java -Xmx"${JAVA_MEM}" -jar "${REFINED_IBD_JAR}" \
    gt="${PREFIX}.chr${chr}.vcf" \
    map="${MAP_DIR}/plink.chr${chr}.GRCh37.map" \
    out="${PREFIX}.chr${chr}" \
    length="${MIN_IBD_LENGTH}"
done
wait
