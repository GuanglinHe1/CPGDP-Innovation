#!/usr/bin/env bash
#
# 0-build-reference-panel.sh
#
# Prepare a chromosome-wise reference panel and assess imputation accuracy.
# All input files, executables, and resource settings are supplied at runtime;
# no institutional paths, sample identifiers, or scheduler configuration are
# stored in this script.
#
# Requirements: bcftools, SHAPEIT5, Minimac4, aggRSquare
#
# Usage:
#   bash 0-build-reference-panel.sh clean <INPUT_VCF_TEMPLATE> <OUTPUT_DIR>
#   bash 0-build-reference-panel.sh phase-impute <SAMPLE_LIST> <OUTPUT_DIR>
#   bash 0-build-reference-panel.sh r2 <TRUTH_VCF> <IMPUTED_VCF> <AF_FILE> <OUTPUT_PREFIX>
#
# The clean command expects INPUT_VCF_TEMPLATE to contain {chr}, for example
# "input/cohort.chr{chr}.vcf.gz".  The phase-impute sample list has two
# whitespace-separated columns: input_vcf and sample_id.
#
# Configuration for phase-impute (environment variables):
#   SHAPEIT5_BIN       SHAPEIT5 executable (default: SHAPEIT5_phase_common)
#   MINIMAC4_BIN       Minimac4 executable (default: minimac4)
#   GENETIC_MAP        Genetic map for the target chromosome (required)
#   PHASING_REFERENCE  Phased reference VCF/BCF (required)
#   MINIMAC_REFERENCE  Minimac reference panel in m3vcf/msav format (required)
#   REGION             Target genomic region (default: chr2)
#   THREADS            Threads per phasing/imputation task (default: 8)

set -euo pipefail

die() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

require_file() {
  [[ -f "$1" ]] || die "Required file does not exist: $1"
}

clean_panel() {
  local input_template="$1"
  local output_dir="$2"
  local chr input_vcf normalized_bcf output_vcf

  [[ "${input_template}" == *'{chr}'* ]] || die 'INPUT_VCF_TEMPLATE must contain {chr}.'
  mkdir -p "${output_dir}"

  for chr in {1..22}; do
    input_vcf="${input_template//\{chr\}/${chr}}"
    normalized_bcf="${output_dir}/reference.norm.chr${chr}.bcf"
    output_vcf="${output_dir}/reference.filter.chr${chr}.vcf.gz"
    require_file "${input_vcf}"

    printf 'Cleaning chromosome %s at %s\n' "${chr}" "$(date --iso-8601=seconds)"
    bcftools norm -m - -a "${input_vcf}" \
      | bcftools view -v snps,indels \
        -e 'strlen(REF) - strlen(ALT) > 50 || strlen(ALT) - strlen(REF) > 50' \
        -Ob -o "${normalized_bcf}"
    bcftools filter --SnpGap 3 --IndelGap 5 -e 'FORMAT/GQ[*] < 20' -S . "${normalized_bcf}" \
      | bcftools view -i 'INFO/AN > 0' -Oz -o "${output_vcf}"
    bcftools index --force "${output_vcf}"
  done
}

phase_and_impute() {
  local sample_list="$1"
  local output_dir="$2"
  local shapeit5_bin="${SHAPEIT5_BIN:-SHAPEIT5_phase_common}"
  local minimac4_bin="${MINIMAC4_BIN:-minimac4}"
  local genetic_map="${GENETIC_MAP:?Set GENETIC_MAP to a genetic-map file.}"
  local phasing_reference="${PHASING_REFERENCE:?Set PHASING_REFERENCE to a phased reference VCF/BCF.}"
  local minimac_reference="${MINIMAC_REFERENCE:?Set MINIMAC_REFERENCE to a Minimac reference panel.}"
  local region="${REGION:-chr2}"
  local threads="${THREADS:-8}"
  local input_vcf sample_id clean_dir phase_dir impute_dir clean_vcf phased_bcf imputed_vcf

  require_file "${sample_list}"
  require_file "${genetic_map}"
  require_file "${phasing_reference}"
  require_file "${minimac_reference}"
  mkdir -p "${output_dir}"
  clean_dir="${output_dir}/clean"
  phase_dir="${output_dir}/phased"
  impute_dir="${output_dir}/imputed"
  mkdir -p "${clean_dir}" "${phase_dir}" "${impute_dir}"

  while read -r input_vcf sample_id extra; do
    [[ -z "${input_vcf}" || "${input_vcf}" == \#* ]] && continue
    [[ -z "${sample_id:-}" || -n "${extra:-}" ]] && die "Expected exactly two columns in ${sample_list}."
    require_file "${input_vcf}"
    clean_vcf="${clean_dir}/${sample_id}.vcf.gz"
    phased_bcf="${phase_dir}/${sample_id}.bcf"
    imputed_vcf="${impute_dir}/${sample_id}.vcf.gz"

    bcftools index --force "${input_vcf}"
    bcftools norm -m -both "${input_vcf}" \
      | bcftools view -v snps,indels \
      | bcftools norm -a -d none -Oz -o "${clean_vcf}"
    bcftools index --force "${clean_vcf}"
    "${shapeit5_bin}" --input "${clean_vcf}" --region "${region}" \
      --map "${genetic_map}" --output "${phased_bcf}" \
      --reference "${phasing_reference}" --thread "${threads}"
    "${minimac4_bin}" --region "${region}" --threads "${threads}" \
      --temp-prefix "${impute_dir}/${sample_id}.tmp" --format GT,DS,GP \
      -O vcf.gz -o "${imputed_vcf}" "${minimac_reference}" "${phased_bcf}"
    bcftools index --force "${imputed_vcf}"
  done < "${sample_list}"
}

evaluate_r2() {
  local truth_vcf="$1"
  local imputed_vcf="$2"
  local af_file="$3"
  local output_prefix="$4"

  require_file "${truth_vcf}"
  require_file "${imputed_vcf}"
  require_file "${af_file}"
  aggRSquare -v "${truth_vcf}" -i "${imputed_vcf}" -o "${output_prefix}" --AF "${af_file}" --detail
}

[[ $# -ge 1 ]] || die 'Usage: bash 0-build-reference-panel.sh {clean|phase-impute|r2} ...'
case "$1" in
  clean)
    [[ $# -eq 3 ]] || die 'Usage: bash 0-build-reference-panel.sh clean <INPUT_VCF_TEMPLATE> <OUTPUT_DIR>'
    clean_panel "$2" "$3"
    ;;
  phase-impute)
    [[ $# -eq 3 ]] || die 'Usage: bash 0-build-reference-panel.sh phase-impute <SAMPLE_LIST> <OUTPUT_DIR>'
    phase_and_impute "$2" "$3"
    ;;
  r2)
    [[ $# -eq 5 ]] || die 'Usage: bash 0-build-reference-panel.sh r2 <TRUTH_VCF> <IMPUTED_VCF> <AF_FILE> <OUTPUT_PREFIX>'
    evaluate_r2 "$2" "$3" "$4" "$5"
    ;;
  *)
    die "Unknown command: $1"
    ;;
esac
