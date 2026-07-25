#!/usr/bin/env Rscript

# Final real-data workflow used in the manuscript.
# Matching and regression adjustment use PC1-PC5 only.

scripts <- c(
    "003_run_full_sample_logistic_gwas.R",
    "004_run_pc_matched_clr_comparison.R",
    "005_lct_full_sample_two_logistic_and_clr.R",
    "007_lct_1to1_clr_comparison.R",
    "016_lct_lactose_positive_control_1to1_1to4_1to6.R",
    "018_lct_tv_15000_case_ratio6.R",
    "019_lct_tv_15000_common_pool_1to1_1to4.R",
    "023_chr2_random100_tv_15000_rare_chisq.R"
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
