#!/usr/bin/env bash
#
# 1-ihs.sh
#
# Within-population selection scan based on the integrated haplotype score
# (iHS).
#
#   vcftools  keep common, well-genotyped variants in Hardy-Weinberg
#             equilibrium
#   beagle    statistical phasing, required by the haplotype-based statistic
#   plink     interpolate the genetic position of every marker; markers whose
#             genetic position cannot be interpolated are dropped
#   selscan   compute the unstandardised iHS and normalise it in frequency
#             bins across all chromosomes
#
# Requirements: vcftools, bgzip, java, beagle, bcftools, plink, selscan
#
# Usage:
#   bash 1-ihs.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Population label; by default the name of the working directory, so that the
# same script can be launched inside one directory per population.
POP="${POP:-${PWD##*/}}"
# Beagle jar file and the directory with its GRCh38 genetic maps.
BEAGLE_JAR="${BEAGLE_JAR:-./software/beagle/beagle.22Jul22.46e.jar}"
BEAGLE_MAP_DIR="${BEAGLE_MAP_DIR:-./software/beagle/plink.GRCh38.map}"
# Directory with the GRCh38 recombination maps used by plink (--cm-map).
GMAP_DIR="${GMAP_DIR:-./resources/hg38gmap}"
# Resources.
THREADS="${THREADS:-8}"
JAVA_MEM="${JAVA_MEM:-120G}"
# Variant filters.
MAF="${MAF:-0.05}"
MAX_MISSING="${MAX_MISSING:-0.8}"
HWE="${HWE:-1e-10}"

for chr in {1..22}; do
  # --- Variant filtering ---------------------------------------------------
  vcftools --gzvcf "chr${chr}.${POP}.vcf.gz" \
    --maf "${MAF}" --max-missing "${MAX_MISSING}" --hwe "${HWE}" \
    --recode --recode-INFO-all --stdout \
    | bgzip -c > "chr${chr}.${POP}.maf005.vcf.gz"

  # --- Phasing -------------------------------------------------------------
  java -Xmx"${JAVA_MEM}" -jar "${BEAGLE_JAR}" \
    gt="chr${chr}.${POP}.maf005.vcf.gz" \
    map="${BEAGLE_MAP_DIR}/plink.chr${chr}.GRCh38.chr.map" \
    out="chr${chr}.beagle" \
    nthreads="${THREADS}"
  bcftools index "chr${chr}.beagle.vcf.gz"

  # --- Genetic positions ---------------------------------------------------
  plink --vcf "chr${chr}.beagle.vcf.gz" --double-id \
    --cm-map "${GMAP_DIR}/chr${chr}.b38.gmap" "${chr}" \
    --make-just-bim --out "tmp_chr${chr}.recom" --threads "${THREADS}"

  # selscan map file: chromosome, marker id, genetic position, physical
  # position, for the markers with a valid genetic position.
  awk '{if ($3 >=0) print "chr"$1,$1"_"$4,$3,$4}' "tmp_chr${chr}.recom.bim" \
    > "chr${chr}.map"
  # Markers with a negative (non-interpolatable) genetic position.
  awk '{if ($3 <0) print "chr"$1 "\t" $4}' "tmp_chr${chr}.recom.bim" \
    > "list_chr${chr}"

  vcftools --gzvcf "chr${chr}.beagle.vcf.gz" \
    --exclude-positions "list_chr${chr}" \
    --recode --recode-INFO-all --stdout \
    | bgzip -c > "chr${chr}.cl.maf005.final.vcf.gz"

  # --- iHS -----------------------------------------------------------------
  selscan --ihs \
    --vcf "chr${chr}.cl.maf005.final.vcf.gz" \
    --map "chr${chr}.map" \
    --out "chr${chr}" \
    --threads "${THREADS}"

  rm -f "tmp_chr${chr}".*
done

# Frequency-bin normalisation over all chromosomes at once.
selscan norm --ihs --files chr{1..22}.ihs.out
