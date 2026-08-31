#!/usr/bin/env bash
#
# 11-chromopainter-finestructure.sh
#
# Fine-scale population structure with ChromoPainter and fineSTRUCTURE, run in
# the HPC mode: fs writes a command file at every stage, the commands are
# executed in parallel, and fs is called again to move on to the next stage.
#
#   stage 1  estimation of the ChromoPainter parameters Ne and mu
#   stage 2  painting of all individuals, producing the coancestry matrix
#   stage 3  fineSTRUCTURE MCMC clustering
#   stage 4  tree building on the inferred clusters
#
# Requirements: finestructure (fs), GNU parallel, phased haplotypes from
#               10-phasing-and-ibd.sh
#
# Usage:
#   bash 11-chromopainter-finestructure.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Name of the fs project; fs creates <PROJECT>/ next to <PROJECT>.cp.
PROJECT="${PROJECT:-poplist_ibdHPC}"
# ChromoPainter id file: <sample id> <population label> <inclusion flag>.
ID_FILE="${ID_FILE:-poplist_ibd.ids}"
# Prefix of the per-chromosome ChromoPainter phase files.
PHASE_PREFIX="${PHASE_PREFIX:-poplist_ibd.chr}"
# Directory holding the per-chromosome recombination rate files.
MAP_DIR="${MAP_DIR:-./maps}"
# Number of MCMC iterations of the fineSTRUCTURE stages.
S3_ITERS="${S3_ITERS:-100000}"
S4_ITERS="${S4_ITERS:-50000}"

# Set up the project and generate the stage 1 command file.
fs "${PROJECT}.cp" -hpc 1 \
  -idfile "${ID_FILE}" \
  -phasefiles "${PHASE_PREFIX}"{1..22}.phase \
  -recombfiles "${MAP_DIR}"/genetic_map_Finalversion_GRCh37_chr{1..22}.recombfile \
  -s3iters "${S3_ITERS}" \
  -s4iters "${S4_ITERS}" \
  -s1minsnps 1000 \
  -s1indfrac 0.1 \
  -go

# Run the four stages: execute the command file, then let fs prepare the next.
for stage in 1 2 3 4; do
  parallel < "${PROJECT}/commandfiles/commandfile${stage}.txt"
  fs "${PROJECT}.cp" -go
done
