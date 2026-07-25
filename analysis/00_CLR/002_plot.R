#!/usr/bin/env Rscript

# Rare-variant calibration, sensitivity analyses, and manuscript figures.

scripts <- c(
    "020_lct_tv_15000_rare_variant_comparison.R",
    "021_lct_tv_15000_rare_permutation_matched.R",
    "022_lct_tv_15000_rare_permutation_vs_chisq.R",
    "024_chr2_random100_tv_15000_rare_permutation.R",
    "025_chr2_random100_tv_15000_rare_perm_vs_chisq.R",
    "027_lct_region_lactose_burden_skat_skato.R",
    "028_lct_lactose_sensitivity_and_variant_plots.R",
    "029_lct_lactose_matched_burden_permutation.R",
    "034_build_firth_sensitivity_supplement.R",
    "036_build_paper_package_recommended.R"
)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) == 1L) {
    script_dir <- dirname(normalizePath(
        sub("^--file=", "", script_arg),
        winslash = "/",
        mustWork = TRUE
    ))
} else {
    script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

for (script in scripts) {
    message("Running ", script)
    source(file.path(script_dir, script), chdir = TRUE)
}
