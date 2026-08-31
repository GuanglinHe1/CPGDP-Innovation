#!/usr/bin/env bash
#
# 1-variant-summary-statistics.sh
#
# Summary statistics of the discovered variant call set:
#   1) genome-wide non-reference allele frequency spectrum
#   2) extraction of biallelic singletons (AC = 1) and doubletons (AC = 2)
#   3) per-sample heterozygosity, joined with the population/language metadata
#
# Requirements: vcftools, bcftools, awk, sort, join
#
# Usage:
#   bash 1-variant-summary-statistics.sh [VCF] [OUT_PREFIX] [POP_METADATA]
#
# Downstream plotting is done by:
#   2-plot-allele-frequency-spectrum.R
#   4-plot-heterozygosity.R

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration (override from the command line or the environment)
# ---------------------------------------------------------------------------
# Joint-called, quality-filtered SNP call set.
VCF="${1:-${VCF:-CPGDP.filtersnp.vcf.gz}}"
# Prefix used for every output file produced by this script.
OUT_PREFIX="${2:-${OUT_PREFIX:-CPGDP}}"
# Two-column table: <sample id> <population/language label>.
POP_METADATA="${3:-${POP_METADATA:-pop2.cov}}"

# ---------------------------------------------------------------------------
# 1. Non-reference allele frequency of every variant
# ---------------------------------------------------------------------------
vcftools --gzvcf "${VCF}" --freq --out "${OUT_PREFIX}"

# ---------------------------------------------------------------------------
# 2. Biallelic singletons and doubletons
#    Keep only records with a single ALT allele and equal REF/ALT length,
#    i.e. biallelic SNPs, and report CHROM POS REF ALT.
# ---------------------------------------------------------------------------
bcftools view "${VCF}" --include 'AC=1' \
  | awk '$0!~/^#/ && $5!~/,/ && length($4)==length($5) {print $1,$2,$4,$5}' \
  > "${OUT_PREFIX}.singleton.bia.snp"

bcftools view "${VCF}" --include 'AC=2' \
  | awk '$0!~/^#/ && $5!~/,/ && length($4)==length($5) {print $1,$2,$4,$5}' \
  > "${OUT_PREFIX}.doubleton.bia.snp"

# ---------------------------------------------------------------------------
# 3. Per-sample heterozygosity
#    vcftools reports O(HOM), E(HOM), N_SITES and F; the level of genetic
#    diversity is computed as (N_SITES - O(HOM)) / N_SITES.
# ---------------------------------------------------------------------------
vcftools --gzvcf "${VCF}" --het --out "${OUT_PREFIX}.hetraw"

awk 'NR>1 {print $1, ($4-$2)/$4}' "${OUT_PREFIX}.hetraw.het" > "${OUT_PREFIX}.het"

# Attach the population / language-family annotation of every sample.
sort -k1,1 "${OUT_PREFIX}.het" \
  | join -1 1 -2 1 - <(sort -k1,1 "${POP_METADATA}") \
  > "${OUT_PREFIX}.het.language.txt"
