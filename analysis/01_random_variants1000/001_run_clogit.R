#!/usr/bin/env Rscript

# Chromosome 2 random rare-variant negative-control analysis.
# The current analysis uses the 100-variant panel and the common 15,000-case
# PC-only matched sets from the TV analysis.

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

clr_dir <- normalizePath(
    file.path(script_dir, "..", "00_CLR"),
    winslash = "/",
    mustWork = TRUE
)

source(
    file.path(clr_dir, "023_chr2_random100_tv_15000_rare_chisq.R"),
    chdir = TRUE
)
source(
    file.path(clr_dir, "024_chr2_random100_tv_15000_rare_permutation.R"),
    chdir = TRUE
)
