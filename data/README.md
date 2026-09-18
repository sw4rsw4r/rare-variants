# UK Biobank input files

Individual-level files are not distributed. Place the following files under `data/` after obtaining access under UK Biobank application 98032 or another approved application.

```text
UKBB_pca.eigenvec
LCT_region_geno.raw
chr2_random100.raw
plink_LCT.afreq
ukb672224.tab
```

Then create phenotype ID lists used by the analysis scripts:

```r
source("R/utils.R")
load_TV()
load_lactose_intolerance()
```

This writes `data/processed/TV.RDS` and `data/processed/lactose_intolerance.RDS`. Those files contain participant IDs and must not be committed.

Assessment-centre balance additionally uses UK Biobank field 54 from the phenotype table. Matching scripts write matched-set files with participant IDs under `results/`. Keep those files local.
