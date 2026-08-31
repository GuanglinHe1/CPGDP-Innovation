# Population Genomic Analysis Pipeline

Analysis code used in this study. 

```bash
BEAGLE_JAR=/opt/beagle/beagle.22Jul22.46e.jar bash 1-ihs.sh
```

## Repository layout

### 1-variant-discovery

Summary of the discovered variant call set.

| Script | Description |
| --- | --- |
| `1-variant-summary-statistics.sh` | Allele frequency spectrum, singleton/doubleton extraction and per-sample heterozygosity |
| `2-plot-allele-frequency-spectrum.R` | Number of SNPs per allele frequency bin |
| `3-plot-per-sample-snp-count.R` | Per-sample SNP counts by population and language phylum |
| `4-plot-heterozygosity.R` | Level of genetic diversity by population and language phylum |
| `5-snp-set-venn.py` | Venn diagram and private variants across variant databases |

### 2-population-history

Population structure, admixture modelling and demographic history.

| Script | Description |
| --- | --- |
| `1-pca-smartpca.sh` | Principal component analysis with smartpca |
| `2-ld-pruning.sh` | LD pruning for the UMAP and ADMIXTURE analyses |
| `3-umap.R` | UMAP embedding of the leading principal components |
| `4-admixture.sh` | Model-based ancestry estimation for K = 2 to 20 |
| `5-pairwise-fst.sh` | Pairwise Fst between all population pairs |
| `6-procrustes-analysis.sh` | Procrustes comparison of the genetic and geographic space |
| `7-outgroup-f3.sh` | Outgroup-f3 shared drift statistics |
| `8-f4-statistics.R` | f4 statistics with a deep African outgroup |
| `9-qpadm.sh` | qpWave / qpAdm admixture modelling |
| `10-phasing-and-ibd.sh` | SHAPEIT phasing and Refined IBD segment detection |
| `11-chromopainter-finestructure.sh` | ChromoPainter / fineSTRUCTURE fine-scale structure |
| `12-chromopainter-fastglobetrotter.sh` | ChromoPainter paintings for fastGLOBETROTTER |
| `13-sourcefind.sh` | SOURCEFIND ancestry decomposition |
| `14-smcpp.sh` | SMC++ effective population size history |
| `15-msmc2.sh` | MSMC2 cross-coalescence rate between population pairs |

### 3-archaic-introgression

Detection of Neanderthal and Denisovan ancestry.

| Script | Description |
| --- | --- |
| `1-sprime.sh` | SPrime scan and matching against the archaic genomes |
| `2-ibdmix.sh` | IBDmix Neanderthal segment calling |
| `3-hmmix.sh` | hmmix per-individual archaic segment decoding |

### 4-natural-selection

Genome-wide scans for positive selection and their composite score.

| Script | Description |
| --- | --- |
| `1-ihs.sh` | Filtering, phasing and the iHS within-population scan |
| `2-pbs-ddaf-xpehh.sh` | PBS, dDAF and XP-EHH scans against a reference population |
| `3-cms-composite-score.R` | Composite of multiple signals (CMS) score |
| `4-cms-clumping.sh` | LD clumping and gene annotation of the CMS signals |

## Software

The pipeline calls the following external tools, which have to be available on
`PATH` or through the variables declared in the scripts:

- Variant handling: `bcftools`, `vcftools`, `bgzip`, `plink`, `plink2`
- Phasing and IBD: `SHAPEIT`, `Beagle`, `Refined IBD`
- Structure and history: `EIGENSOFT (smartpca)`, `ADMIXTURE`, `AdmixTools
  (qpAdm)`, `ChromoPainterv2`, `fineSTRUCTURE (fs)`, `SOURCEFIND`, `SMC++`,
  `MSMC2`, `GNU parallel`
- Archaic introgression: `SPrime`, `map_arch`, `IBDmix`, `hmmix`
- Selection: `selscan`
- `R` with `ggplot2`, `ggrepel`, `ggrastr`, `dplyr`, `data.table`, `uwot`,
  `scales`, `admixtools`
- `Python 3` with `matplotlib`, `geneview`
