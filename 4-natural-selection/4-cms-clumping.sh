#!/usr/bin/env bash
#
# 4-cms-clumping.sh
#
# Post-processing of the composite selection score: the per-chromosome VCFs
# are concatenated into a genome-wide PLINK data set, and the CMS signals are
# clumped by linkage disequilibrium so that every candidate region is
# represented by its lead variant and annotated with the overlapping genes.
#
# Requirements: bcftools, plink
#
# Usage:
#   bash 4-cms-clumping.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Population label; by default the name of the working directory.
POP="${POP:-${PWD##*/}}"
# BED file with the gene annotation used to label the clumped regions.
GENE_BED="${GENE_BED:-./resources/hg38.gtf.gene.name.bed}"
# Clumping thresholds. CLUMP_P1 is the empirical CMS p value an index variant
# must reach, CLUMP_P2 the threshold for the variants clumped around it; both
# are study specific and have to be set explicitly.
CLUMP_P1="${CLUMP_P1:?set CLUMP_P1 to the index variant p value threshold}"
CLUMP_P2="${CLUMP_P2:?set CLUMP_P2 to the secondary variant p value threshold}"
CLUMP_R2="${CLUMP_R2:-0.2}"
CLUMP_KB="${CLUMP_KB:-500}"

# Index the per-chromosome VCFs in parallel, then concatenate them.
for chr in {1..22}; do
  bcftools index "chr${chr}.cl.maf005.final.vcf.gz" &
done
wait

bcftools concat chr{1..22}.cl.maf005.final.vcf.gz -o "${POP}.allchr.vcf.gz"

# The deduplicated call set is the LD reference panel of the clumping step.
plink --vcf "${POP}.allchr.dedup.vcf.gz" --make-bed --out "${POP}" --double-id

plink --bfile "${POP}" \
  --clump "${POP}.cms.txt" \
  --clump-p1 "${CLUMP_P1}" \
  --clump-p2 "${CLUMP_P2}" \
  --clump-r2 "${CLUMP_R2}" \
  --clump-kb "${CLUMP_KB}" \
  --clump-snp-field chrpos \
  --clump-field cmsP \
  --clump-range "${GENE_BED}" \
  --out "CMS/${POP}.cms"
