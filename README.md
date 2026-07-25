# Rare-variant association analysis

This repository contains the real-data analysis code for testing PC-only
ancestry matching followed by conditional logistic regression (CLR).

The current workflow compares:

- plain full-sample logistic regression;
- full-sample logistic regression adjusted for PC1-PC5;
- PC1-PC5-only matching at 1:1, 1:4, and 1:6 ratios followed by CLR;
- asymptotic and within-stratum permutation calibration for rare variants;
- LCT-region burden, SKAT, SKAT-O, and Firth sensitivity analyses.

Sex is not used in the matching or regression adjustment. Genotypes are not
used to construct matched sets.

## Analysis sets

- LCT region and TV: negative-control analysis.
- LCT region and lactose intolerance: positive-control analysis.
- Random chromosome 2 rare variants and TV: independent negative-control
  analysis.
- rs4988235-A: European lactase-persistence sentinel.

Rare variants are defined as MAF < 0.01. Common variants have MAF >= 0.01.

The TV matching analysis uses a fixed pool of 15,000 cases for all three
matching ratios. The lactose-intolerance analysis retains all available cases.

## Main scripts

The scripts are numbered in the order in which the analysis was developed.
The two entry points below run the final workflow:

```bash
Rscript analysis/00_CLR/001_run_clogit.R
Rscript analysis/00_CLR/002_plot.R
```

The main component scripts are:

- `003_run_full_sample_logistic_gwas.R`: full-sample logistic GWAS.
- `004_run_pc_matched_clr_comparison.R`: PC-only matching and CLR helpers.
- `010_pc_matching_helpers.R`: shared matching, CLR, and QQ functions.
- `016_lct_lactose_positive_control_1to1_1to4_1to6.R`: positive control.
- `018_lct_tv_15000_case_ratio6.R` and
  `019_lct_tv_15000_common_pool_1to1_1to4.R`: final TV matching analyses.
- `020`-`025`: LCT and chromosome 2 rare-variant calibration, including
  10,000 within-stratum permutations.
- `027`-`029`: burden, SKAT, SKAT-O, Firth, and matched burden sensitivity
  analyses.
- `036_build_paper_package_recommended.R`: manuscript figures and tables.

## Required input files

Individual-level UK Biobank data are not included in this repository. The
scripts expect the following files under `data/`:

```text
UKBB_pca.eigenvec
UKBB_pca.eigenval
LCT_region_geno.raw
chr2_random100.raw
plink_LCT.afreq
plink_chr2.afreq
processed/TV.RDS
processed/lactose_intolerance.RDS
processed/qualification.RDS
```

Phenotype RDS files contain a list with `case` and `ctrl` participant IDs.
PLINK raw genotype columns follow the `<variant>_<counted allele>` naming
convention.

## R packages

The real-data workflow uses `data.table`, `dplyr`, `tidyr`, `ggplot2`,
`MatchIt`, `survival`, `patchwork`, `SKAT`, and `logistf`.

Generated results, plots, temporary files, and individual-level data are
excluded from version control.
