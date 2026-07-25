#!/usr/bin/env Rscript

# Extend the existing LCT/lactose positive-control comparison with a PC-only
# 1:6 matched CLR. Existing 1:1 and 1:4 outputs are reused unchanged. The
# phenotype has 386 available cases, so all cases are retained; the candidate
# control pool is the same 15x pool used by the existing 1:1/1:4 analyses.
# Sex and permutation are not used.

options(stringsAsFactors = FALSE)
if (!requireNamespace("MatchIt", quietly = TRUE)) stop("The MatchIt package is required")
if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required")
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggplot2))

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

data_dir <- file.path(project_root, "data")
comparison_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison")
result_dir <- file.path(comparison_root, "lactose_intolerance")
matching_root <- file.path(comparison_root, "matching_1to6")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
num_controls <- 6L
max_cases <- 1000L
control_pool_multiplier <- 15L
region <- "LCT_region"
phenotype <- "lactose_intolerance"
phenotype_seeds <- c(lactose_intolerance = 1103L)

load_function_definitions <- function(path, names_to_load) {
  parsed <- parse(path)
  for (expr in parsed) {
    if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
      name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (length(name) == 1L && name %in% names_to_load) eval(expr, envir = parent.frame())
    }
  }
}
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"), c("read_header", "read_selected_tsv", "read_genotypes", "load_phenotype_ids", "make_binary_phenotype"))
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "004_run_pc_matched_clr_comparison.R"), c("compute_balance", "add_mahalanobis_distances", "match_phenotype", "run_clogit_gwas"))
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "010_pc_matching_helpers.R"), c("draw_original_style"))
required <- c("read_header", "read_selected_tsv", "read_genotypes", "load_phenotype_ids", "make_binary_phenotype", "compute_balance", "add_mahalanobis_distances", "match_phenotype", "run_clogit_gwas", "draw_original_style")
if (!all(vapply(required, exists, logical(1L), envir = environment(), inherits = FALSE))) stop("Failed to load required helper functions")

pc_data <- read_selected_tsv(file.path(data_dir, "UKBB_pca.eigenvec"), c("IID", pc_names), c("character", rep("numeric", num_pcs)))
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
if (anyDuplicated(pc_data$IID)) stop("Duplicate IID in PCA data")
genotype_frame <- read_genotypes(file.path(data_dir, "LCT_region_geno.raw"))
variant_names <- setdiff(names(genotype_frame), "IID")
genotype_index <- match(pc_data$IID, genotype_frame$IID)
if (anyNA(genotype_index)) stop("Some PC IIDs are missing from LCT genotype data")
genotype_all <- as.matrix(genotype_frame[genotype_index, variant_names, drop = FALSE])
storage.mode(genotype_all) <- "double"
stopifnot(identical(pc_data$IID, genotype_frame$IID[genotype_index]))

message("Matching lactose_intolerance on PC1-PC5 only (1:6)")
match_result <- match_phenotype(pc_data, phenotype)
matched_data <- match_result$data
if (any(grepl("sex", names(matched_data), ignore.case = TRUE))) stop("Sex information entered the matching data")
set_counts <- table(matched_data$STRATUM, matched_data$case)
if (!all(c("0", "1") %in% colnames(set_counts)) || any(set_counts[, "1"] != 1L) || any(set_counts[, "0"] != 6L)) stop("Incomplete 1:6 matching")
if (anyDuplicated(matched_data$IID)) stop("Duplicate IID in matched sets")
write.table(match_result$summary, file.path(matching_root, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(match_result$summary, file.path(matching_root, phenotype, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

genotype_index <- match(matched_data$IID, pc_data$IID)
if (anyNA(genotype_index) || any(as.character(matched_data$IID) != as.character(pc_data$IID[genotype_index]))) stop("Matched IID alignment failed")
genotype_matched <- genotype_all[genotype_index, , drop = FALSE]
message("Running lactose_intolerance 1:6 matched CLR")
clr_1to6 <- run_clogit_gwas(genotype_matched, matched_data)
if ("N_STRATA_COMPLETE_1TO4" %in% names(clr_1to6)) names(clr_1to6)[names(clr_1to6) == "N_STRATA_COMPLETE_1TO4"] <- "N_STRATA_COMPLETE_1TO6"
clr_1to6$TEST <- "MATCHED_CLR_1TO6_WALD_CHISQ1"
write.table(clr_1to6, file.path(result_dir, "matched_clr_1to6_gwas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plain <- read.delim(file.path(result_dir, "logistic_plain.tsv"), check.names = FALSE)
pc <- read.delim(file.path(result_dir, "logistic_pc_adjusted.tsv"), check.names = FALSE)
if (any(grepl("sex", names(pc), ignore.case = TRUE)) || any(grepl("sex", names(plain), ignore.case = TRUE))) stop("Sex information entered full-sample models")
clr4 <- read.delim(file.path(project_root, "results", "00_CLR", "pc_adjusted_vs_matched_clr", region, phenotype, "matched_clr_gwas.tsv"), check.names = FALSE)
clr1 <- read.delim(file.path(result_dir, "matched_clr_1to1_gwas.tsv"), check.names = FALSE)
plain_ok <- plain[plain$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(plain_ok) <- c("SNP", "CHISQ_PLAIN", "P_PLAIN", "N_PLAIN")
pc_ok <- pc[pc$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
clr4_ok <- clr4[clr4$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]; names(clr4_ok) <- c("SNP", "CHISQ_1TO4", "P_1TO4", "N_1TO4", "N_1TO4_BASE")
clr1_ok <- clr1[clr1$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]; names(clr1_ok) <- c("SNP", "CHISQ_1TO1", "P_1TO1", "N_1TO1", "N_1TO1_BASE")
clr6_ok <- clr_1to6[clr_1to6$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]; names(clr6_ok) <- c("SNP", "CHISQ_1TO6", "P_1TO6", "N_1TO6", "N_1TO6_BASE")
comparison <- Reduce(function(x, y) merge(x, y, by = "SNP", sort = FALSE), list(plain_ok, pc_ok, clr1_ok, clr4_ok, clr6_ok))
if (nrow(comparison) == 0L) stop("No common SNPs across positive-control models")
comparison <- comparison[order(match(comparison$SNP, plain_ok$SNP)), , drop = FALSE]
write.table(comparison, file.path(result_dir, "positive_control_comparison_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_columns <- c("CHISQ_1TO1", "CHISQ_1TO4", "CHISQ_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
plot_labels <- c("Matched 1:1", "Matched 1:4", "Matched 1:6", "Full sample + PCs", "Full sample (No PCs)")
full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
plot_full <- draw_original_style(comparison, plot_columns, plot_labels, full_max)
plot_zoom <- draw_original_style(comparison, plot_columns, plot_labels, 10)
ggsave(file.path(result_dir, "qq_positive_control_1to1_1to4_1to6_original_style_full.pdf"), plot_full, width = 19, height = 4.2, units = "in", device = "pdf")
ggsave(file.path(result_dir, "qq_positive_control_1to1_1to4_1to6_original_style_zoom10.pdf"), plot_zoom, width = 19, height = 4.2, units = "in", device = "pdf")
ggsave(file.path(result_dir, "qq_positive_control_1to1_1to4_1to6_original_style_zoom10.png"), plot_zoom, width = 19, height = 4.2, units = "in", dpi = 220)

sum1 <- read.delim(file.path(comparison_root, "matching_1to1", phenotype, "matching_summary.tsv"), check.names = FALSE)
sum4 <- read.delim(file.path(project_root, "results", "00_CLR", "pc_adjusted_vs_matched_clr", "matching", phenotype, "matching_summary.tsv"), check.names = FALSE)
summary <- data.frame(
  REGION = region, PHENOTYPE = phenotype, N_COMMON_MODELS = nrow(comparison),
  N_PLAIN = sum(plain$STATUS == "ok"), N_PC_ADJUSTED = sum(pc$STATUS == "ok"),
  N_CASE_1TO1 = sum1$N_CASE_MATCHED, N_CONTROL_1TO1 = sum1$N_CONTROL_MATCHED,
  N_CASE_1TO4 = sum4$N_CASE_MATCHED, N_CONTROL_1TO4 = sum4$N_CONTROL_MATCHED,
  N_CASE_1TO6 = match_result$summary$N_CASE_MATCHED, N_CONTROL_1TO6 = match_result$summary$N_CONTROL_MATCHED,
  N_1TO1_OK = sum(clr1$STATUS == "ok"), N_1TO4_OK = sum(clr4$STATUS == "ok"), N_1TO6_OK = sum(clr_1to6$STATUS == "ok"),
  MAX_ABS_SMD_POST_1TO1 = sum1$MAX_ABS_SMD_POST, MAX_ABS_SMD_POST_1TO4 = sum4$MAX_ABS_SMD_POST, MAX_ABS_SMD_POST_1TO6 = match_result$summary$MAX_ABS_SMD_POST,
  LAMBDA_1TO1 = median(comparison$CHISQ_1TO1) / qchisq(.5, 1), LAMBDA_1TO4 = median(comparison$CHISQ_1TO4) / qchisq(.5, 1), LAMBDA_1TO6 = median(comparison$CHISQ_1TO6) / qchisq(.5, 1),
  LAMBDA_PC_ADJUSTED = median(comparison$CHISQ_PC_ADJUSTED) / qchisq(.5, 1), LAMBDA_PLAIN = median(comparison$CHISQ_PLAIN) / qchisq(.5, 1),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(result_dir, "analysis_summary_lactose_positive_control_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "LCT_region / lactose_intolerance positive-control comparison",
  "Models: full-sample plain logistic; full-sample PC-adjusted logistic; PC-only matched 1:1, 1:4, and 1:6 CLR",
  "Matching: PC1-PC5 only, Mahalanobis nearest-neighbor, without replacement; no sex",
  "All 386 available lactose_intolerance cases were retained; candidate controls were capped at 15x cases (5,790)",
  "Matched CLR: exact case ~ SNP + strata(STRATUM); squared SNP Wald z referenced to chi-square(1); no permutation",
  "QQ plots use the original ggplot2/facet_wrap/theme_minimal style and fixed 0-10 zoom plus full-scale output",
  "Existing 1:1 and 1:4 matched results were reused; only 1:6 matching/CLR was newly run."
), file.path(result_dir, "METHODS_lactose_positive_control_1to1_1to4_1to6.txt"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo_lactose_positive_control_1to1_1to4_1to6.txt"))
message("Completed lactose positive-control 1:1/1:4/1:6 comparison: ", result_dir)
