#!/usr/bin/env bash
#
# 2-pbs-ddaf-xpehh.sh
#
# Between-population selection scans against the CHB reference population,
# with CEU as the third population of the branch statistic.
#
#   PBS      population branch statistic from the pairwise Fst of the trio
#            target / CHB / CEU
#   dDAF     difference in derived allele frequency between the target and the
#            CHB population, polarised with the ancestral allele calls
#   XP-EHH   cross-population extended haplotype homozygosity between the
#            target and the CHB population
#
# The three statistics, together with the iHS of 1-ihs.sh, are the components
# of the composite score computed by 3-cms-composite-score.R.
#
# Requirements: java, beagle, bcftools, plink, vcftools, bgzip, selscan,
#               python3, RunPBS.py, derta.DAF.py
#
# Usage:
#   bash 2-pbs-ddaf-xpehh.sh

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Population label; by default the name of the working directory.
POP="${POP:-${PWD##*/}}"
# Reference populations of the scans.
REF_POP="${REF_POP:-CHB}"
OUT_POP="${OUT_POP:-CEU}"
# Beagle jar file and the directory with its GRCh38 genetic maps.
BEAGLE_JAR="${BEAGLE_JAR:-./software/beagle/beagle.22Jul22.46e.jar}"
BEAGLE_MAP_DIR="${BEAGLE_MAP_DIR:-./software/beagle/plink.GRCh38.map}"
# Directory with the GRCh38 recombination maps used by plink (--cm-map).
GMAP_DIR="${GMAP_DIR:-./resources/hg38gmap}"
# Directory with the per-chromosome ancestral allele files.
AA_DIR="${AA_DIR:-./resources/DAF_AA}"
# Helper scripts shipped with this repository.
RUN_PBS_PY="${RUN_PBS_PY:-RunPBS.py}"
DELTA_DAF_PY="${DELTA_DAF_PY:-derta.DAF.py}"
# Resources.
THREADS="${THREADS:-8}"
JAVA_MEM="${JAVA_MEM:-120G}"

mkdir -p pbs DAF xpehh CMS

for chr in {1..22}; do
  # --- PBS -----------------------------------------------------------------
  java -Xmx"${JAVA_MEM}" -jar "${BEAGLE_JAR}" \
    gt="pbs/${POP}.chr${chr}.${REF_POP}.${OUT_POP}.cl.vcf.gz" \
    map="${BEAGLE_MAP_DIR}/plink.chr${chr}.GRCh38.chr.map" \
    out="pbs/${POP}.chr${chr}.${REF_POP}.${OUT_POP}.beagle" \
    nthreads="${THREADS}"
  bcftools index "pbs/${POP}.chr${chr}.${REF_POP}.${OUT_POP}.beagle.vcf.gz"

  python3 "${RUN_PBS_PY}" \
    --popfile pbs/pbs.pop \
    --pbspop "${POP}:${REF_POP}:${OUT_POP}" \
    --dropna "pbs/${POP}.chr${chr}.${REF_POP}.${OUT_POP}.beagle.vcf.gz" \
    > "pbs/${POP}_chr${chr}.pbs"

  # --- Split the phased trio VCF into the target and reference panels ------
  bcftools view -S "${POP}.list" \
    "pbs/${POP}.chr${chr}.${REF_POP}.${OUT_POP}.beagle.vcf.gz" \
    -O z -o "DAF/${POP}.beagle.chr${chr}.vcf.gz"
  bcftools view -S "${REF_POP}.list" \
    "pbs/${POP}.chr${chr}.${REF_POP}.${OUT_POP}.beagle.vcf.gz" \
    -O z -o "DAF/${REF_POP}.beagle.chr${chr}.vcf.gz"

  # --- Genetic positions ---------------------------------------------------
  plink --vcf "DAF/${POP}.beagle.chr${chr}.vcf.gz" --double-id \
    --cm-map "${GMAP_DIR}/chr${chr}.b38.gmap" "${chr}" \
    --make-just-bim --out "DAF/tmp_chr${chr}.recom" --threads 1

  awk '{if ($3 >=0) print "chr"$1,$1"_"$4,$3,$4}' "DAF/tmp_chr${chr}.recom.bim" \
    > "DAF/chr${chr}.map"
  awk '{if ($3 <0) print "chr"$1 "\t" $4}' "DAF/tmp_chr${chr}.recom.bim" \
    > "DAF/list_chr${chr}"

  vcftools --gzvcf "DAF/${POP}.beagle.chr${chr}.vcf.gz" \
    --exclude-positions "DAF/list_chr${chr}" \
    --recode --recode-INFO-all --stdout \
    | bgzip -c > "DAF/${POP}.phased.chr${chr}.vcf.gz"
  vcftools --gzvcf "DAF/${REF_POP}.beagle.chr${chr}.vcf.gz" \
    --exclude-positions "DAF/list_chr${chr}" \
    --recode --recode-INFO-all --stdout \
    | bgzip -c > "DAF/${REF_POP}.phased.chr${chr}.vcf.gz"

  # --- dDAF ----------------------------------------------------------------
  vcftools --gzvcf "DAF/${POP}.phased.chr${chr}.vcf.gz" \
    --freq --out "DAF/${POP}.chr${chr}"
  vcftools --gzvcf "DAF/${REF_POP}.phased.chr${chr}.vcf.gz" \
    --freq --out "DAF/${REF_POP}.chr${chr}"

  python3 "${DELTA_DAF_PY}" \
    --target "DAF/${POP}.chr${chr}.frq" \
    --other "DAF/${REF_POP}.chr${chr}.frq" \
    --AA "${AA_DIR}/ancestral_sites.hg38.chr${chr}.txt" \
    --out "DAF/chr${chr}.DAF"

  # --- XP-EHH --------------------------------------------------------------
  selscan --xpehh \
    --vcf "DAF/${POP}.phased.chr${chr}.vcf.gz" \
    --vcf-ref "DAF/${REF_POP}.phased.chr${chr}.vcf.gz" \
    --out "xpehh/${POP}_${REF_POP}_chr${chr}.xpehh" \
    --map "DAF/chr${chr}.map" \
    --threads 6
done

# Normalisation over all chromosomes at once.
selscan norm --xpehh --files xpehh/"${POP}_${REF_POP}"_chr{1..22}.xpehh.xpehh.out
