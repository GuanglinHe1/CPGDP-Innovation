#!/usr/bin/env bash
#
# 12-chromopainter-fastglobetrotter.sh
#
# ChromoPainter painting runs that provide the input of fastGLOBETROTTER,
# which dates and describes the admixture events of the target populations.
#
#   step 1  expectation-maximisation estimation of Ne and mu
#   step 2  donor versus donor painting, used to build the null expectation
#   step 3  donor versus target painting, used for the admixture inference
#
# Requirements: ChromoPainterv2, phased haplotypes from 10-phasing-and-ibd.sh
#
# Usage:
#   bash 12-chromopainter-fastglobetrotter.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Prefix of the per-chromosome phased ChromoPainter phase files.
PHASE_PREFIX="${PHASE_PREFIX:-poplist.chr}"
# Directory holding the per-chromosome recombination rate files.
MAP_DIR="${MAP_DIR:-./maps}"
# ChromoPainter id file: <sample id> <population label> <inclusion flag>.
ID_FILE="${ID_FILE:-pops.ids}"
# Donor/recipient definition files of the three painting steps.
POPFILE_EM="${POPFILE_EM:-popfileSurr.txt}"
POPFILE_DONOR="${POPFILE_DONOR:-popfileSurr2.txt}"
POPFILE_TARGET="${POPFILE_TARGET:-popfileSurr3.txt}"
# Ne and mu averaged over the chromosomes of the EM step below; they have to
# be filled in before the donor and target paintings are launched.
NE="${NE:-}"
MU="${MU:-}"

# --- Step 1: EM estimation of Ne and mu ------------------------------------
for chr in {1..22}; do
  ChromoPainterv2 \
    -g "${PHASE_PREFIX}${chr}.phased.phase" \
    -r "${MAP_DIR}/genetic_map_Finalversion_GRCh37_chr${chr}.recombfile" \
    -t "${ID_FILE}" \
    -f "${POPFILE_EM}" 1 10 \
    -s 0 -i 10 -in -iM \
    -o "output_estimateEM_Chr${chr}"
done

if [[ -z "${NE}" || -z "${MU}" ]]; then
  echo "Set NE and MU to the values averaged over the EM output before" >&2
  echo "running the donor and target paintings." >&2
  exit 1
fi

# --- Step 2: donor versus donor painting -----------------------------------
for chr in {1..22}; do
  ChromoPainterv2 \
    -g "${PHASE_PREFIX}${chr}.phased.phase" \
    -r "${MAP_DIR}/genetic_map_Finalversion_GRCh37_chr${chr}.recombfile" \
    -t "${ID_FILE}" \
    -f "${POPFILE_DONOR}" 0 0 \
    -s 0 -n "${NE}" -M "${MU}" \
    -o "${PHASE_PREFIX}${chr}_Donor_v_Donor" &
done
wait

# --- Step 3: donor versus target painting ----------------------------------
# Ten painting samples per haplotype are drawn for the admixture inference.
for chr in {1..22}; do
  ChromoPainterv2 \
    -g "${PHASE_PREFIX}${chr}.phased.phase" \
    -r "${MAP_DIR}/genetic_map_Finalversion_GRCh37_chr${chr}.recombfile" \
    -t "${ID_FILE}" \
    -f "${POPFILE_TARGET}" 0 0 \
    -s 10 -n "${NE}" -M "${MU}" \
    -o "${PHASE_PREFIX}${chr}_Donor_v_Target" &
done
wait
