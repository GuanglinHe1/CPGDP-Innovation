#!/usr/bin/env bash
#
# 3-hmmix.sh
#
# Detection of archaic introgression with hmmix, a hidden Markov model that
# classifies the genome of a single individual into archaic and non-archaic
# states from the density of variants that are absent from an African
# outgroup.
#
#   create_ingroup  collect the outgroup-private variants of the individual
#   train           fit the HMM transition and emission parameters
#   decode          call the archaic segments and assign them to an archaic
#                   source population
#
# Requirements: hmmix, awk
#
# Usage:
#   bash 3-hmmix.sh <SAMPLE_ID>

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Individual to analyse.
SAMPLE="${1:?usage: bash 3-hmmix.sh <SAMPLE_ID>}"
# Directory with the hmmix GRCh38 resources: strict callability mask,
# outgroup sample list, mutation rate file and ancestral allele FASTA.
RESOURCE_DIR="${RESOURCE_DIR:-./resources/hmmix/hg38}"
# Directory with the archaic genotypes used to label the decoded segments.
ARCHAIC_DIR="${ARCHAIC_DIR:-./resources/archaic/hg38archaicvcf}"
# Joint-called genotype VCF of the study samples.
VCF="${VCF:-snp.gt.vcf.gz}"

CALLABILITY_MASK="${RESOURCE_DIR}/hg38_strick_callability_mask.bed"
MUTATION_RATE="${RESOURCE_DIR}/hg38_mutationrate.bed"
OUTGROUP="${RESOURCE_DIR}/hg38_Outgroup_1000g_HGDP.txt"
ANCESTRAL="${RESOURCE_DIR}/hg38_ancestral"

hmmix create_ingroup \
  -ind="${SAMPLE}" \
  -vcf="${VCF}" \
  -weights="${CALLABILITY_MASK}" \
  -out=obs \
  -outgroup="${OUTGROUP}" \
  -ancestral="${ANCESTRAL}"/*fa

hmmix train \
  -obs="obs.${SAMPLE}.txt" \
  -weights="${CALLABILITY_MASK}" \
  -mutrates="${MUTATION_RATE}" \
  -out="trained.${SAMPLE}.json"

# Keep the decoded segments of the autosomes and sex chromosomes only.
hmmix decode \
  -obs="obs.${SAMPLE}.txt" \
  -weights="${CALLABILITY_MASK}" \
  -mutrates="${MUTATION_RATE}" \
  -param="trained.${SAMPLE}.json" \
  -admixpop="${ARCHAIC_DIR}"/*bcf \
  | awk '$1~/chr/' > "output.${SAMPLE}.txt"
