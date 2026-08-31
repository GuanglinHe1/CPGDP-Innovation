#!/usr/bin/env bash
#
# 1-sprime.sh
#
# Reference-free detection of archaic introgression with SPrime, followed by
# the assignment of the detected archaic haplotypes to the Neanderthal and the
# Denisovan genome.
#
#   sprime     score putatively introgressed alleles per chromosome, using an
#              African population as the non-introgressed outgroup
#   map_arch   annotate every scored allele with its match state in the
#              high-coverage Vindija Neanderthal and Altai Denisovan genomes,
#              restricted to the callable regions of the archaic genomes
#   R scripts  summarise the match rates and draw the match-rate contour plot
#              that separates Neanderthal-like from Denisovan-like segments
#
# Requirements: java, sprime.jar, map_arch, Rscript
#
# Usage:
#   bash 1-sprime.sh <POPULATION>

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Target population; the input VCF is <POPULATION>.YRI.vcf.gz, i.e. the target
# population merged with the YRI outgroup.
POPULATION="${1:?usage: bash 1-sprime.sh <POPULATION>}"
# SPrime jar file.
SPRIME_JAR="${SPRIME_JAR:-./sprime.jar}"
# Directory with the GRCh38 genetic maps (plink.chr<N>.GRCh38.chr.map) and the
# YRI outgroup sample list (YRI.ids).
MAP_DIR="${MAP_DIR:-./resources/maphg38}"
# Directory with the archaic reference data: per-chromosome VCFs and the
# callability mask of each archaic genome.
ARCHAIC_DIR="${ARCHAIC_DIR:-./resources/archaic/hg38}"
# Memory given to the SPrime run.
JAVA_MEM="${JAVA_MEM:-80G}"
# Helper R scripts shipped with this repository.
SCORE_SUMMARY_R="${SCORE_SUMMARY_R:-score_summary.r}"
PLOT_CONTOUR_R="${PLOT_CONTOUR_R:-plot_contour.r}"

for chr in {1..22}; do
  # --- Score the putatively introgressed alleles --------------------------
  java -Xmx"${JAVA_MEM}" -jar "${SPRIME_JAR}" \
    gt="${POPULATION}.YRI.vcf.gz" \
    map="${MAP_DIR}/plink.chr${chr}.GRCh38.chr.map" \
    outgroup="${MAP_DIR}/YRI.ids" \
    chrom="chr${chr}" \
    out="chr${chr}"

  # map_arch expects the chromosome column without the "chr" prefix.
  sed -i 's/chr//g' "chr${chr}.score"

  # --- Match the scored alleles against the Neanderthal genome ------------
  map_arch --kp --sep '\t' \
    --tag "AltaiNean" \
    --mskbed "${ARCHAIC_DIR}/vindija3319/hg38.chr${chr}_mask.bed" \
    --vcf "${ARCHAIC_DIR}/vindija3319/chr${chr}.Vindija3319.vcf.gz" \
    --score "chr${chr}.score" \
    > "tmp1.${chr}.mscore"

  # --- Match the scored alleles against the Denisovan genome --------------
  map_arch --kp --sep '\t' \
    --tag "AltaiDeni" \
    --mskbed "${ARCHAIC_DIR}/denisova/hg38.chr${chr}_mask.bed" \
    --vcf "${ARCHAIC_DIR}/denisova/chr${chr}.denisova.vcf.gz" \
    --score "tmp1.${chr}.mscore" \
    > "tmp2.${chr}.mscore"

  mv "tmp2.${chr}.mscore" "out.chr${chr}.mscore"
  rm -f "tmp1.${chr}.mscore"
done

# --- Summarise the match rates and draw the contour plot -------------------
Rscript "${SCORE_SUMMARY_R}" ./ match.summary.txt
Rscript "${PLOT_CONTOUR_R}" match.summary.txt "${POPULATION}.contour"
