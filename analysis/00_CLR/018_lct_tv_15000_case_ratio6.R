#!/usr/bin/env Rscript

# LCT_region / TV: 15,000-case PC-only 1:6 matched CLR analysis.
# The 15,000 cases are fixed for this analysis; all available TV controls
# are candidates, with 90,000 controls selected without replacement.
# No sex covariate and no permutation are used.

options(stringsAsFactors = FALSE)
if (!requireNamespace("MatchIt", quietly = TRUE)) stop("The MatchIt package is required")
if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required")
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggplot2))
setTimeLimit(elapsed = 10 * 60 * 60, transient = FALSE)

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}
data_dir <- file.path(project_root, "data")
output_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_case_ratio6")
matching_root <- file.path(output_root, "matching", "ratio_1to6")
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
phenotype <- "TV"
region <- "LCT_region"
case_seed <- 3101L
target_cases <- 15000L
ratio <- 6L

message_time <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
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
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "004_run_pc_matched_clr_comparison.R"), c("compute_balance"))
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "010_pc_matching_helpers.R"), c("match_ratio", "run_clr", "draw_original_style"))

pc_data <- read_selected_tsv(file.path(data_dir, "UKBB_pca.eigenvec"), c("IID", pc_names), c("character", rep("numeric", num_pcs)))
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
if (anyDuplicated(pc_data$IID)) stop("Duplicate IID in PCA data")
y <- make_binary_phenotype(pc_data$IID, phenotype)
case_available <- which(y == 1L)
control_available <- which(y == 0L)
if (length(case_available) < target_cases || length(control_available) < target_cases * ratio) stop("Insufficient TV cases or controls")
set.seed(case_seed)
selected_case_pool <- sample(case_available, target_cases, replace = FALSE)
selected_case_iids <- pc_data$IID[selected_case_pool]
write.table(data.frame(IID = selected_case_iids), file.path(output_root, "selected_case_pool.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message_time("Selected ", target_cases, " cases; ", length(control_available), " TV controls available; target controls = ", target_cases * ratio)

start_matching <- Sys.time()
match_ratio_result <- match_ratio(ratio, selected_case_pool)
if (any(grepl("sex", names(match_ratio_result$data), ignore.case = TRUE))) stop("Sex entered matched data")
matched_summary <- match_ratio_result$summary
write.table(matched_summary, file.path(output_root, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(match_ratio_result$balance, file.path(output_root, "pc_balance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message_time("Matching finished in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), start_matching, units = "mins"))), " minutes; post-match max |SMD| = ", signif(matched_summary$MAX_ABS_SMD_POST, 6))

genotype_frame <- read_genotypes(file.path(data_dir, "LCT_region_geno.raw"))
variant_names <- setdiff(names(genotype_frame), "IID")
genotype_index <- match(pc_data$IID, genotype_frame$IID)
if (anyNA(genotype_index)) stop("Some PC IIDs are missing from LCT genotype data")
genotype_all <- as.matrix(genotype_frame[genotype_index, variant_names, drop = FALSE])
storage.mode(genotype_all) <- "double"
matched <- match_ratio_result$data
idx <- match(matched$IID, pc_data$IID)
if (anyNA(idx) || any(as.character(matched$IID) != as.character(pc_data$IID[idx]))) stop("Matched IID alignment failed")
start_clr <- Sys.time()
clr <- run_clr(genotype_all[idx, , drop = FALSE], matched, ratio)
write.table(clr, file.path(output_root, "matched_clr_gwas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message_time("CLR finished in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), start_clr, units = "mins"))), " minutes; successful SNPs = ", sum(clr$STATUS == "ok"))

plain <- read.delim(file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison", "TV", "logistic_plain.tsv"), check.names = FALSE)
pc <- read.delim(file.path(project_root, "results", "00_CLR", "full_sample_logistic_gwas", region, phenotype, "logistic_gwas.tsv"), check.names = FALSE)
plain_ok <- plain[plain$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(plain_ok) <- c("SNP", "CHISQ_PLAIN", "P_PLAIN", "N_PLAIN")
pc_ok <- pc[pc$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
comparison <- merge(merge(plain_ok, pc_ok, by = "SNP", sort = FALSE), clr[clr$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE")], by = "SNP", sort = FALSE)
names(comparison)[(ncol(comparison)-3):ncol(comparison)] <- c("CHISQ_1TO6", "P_1TO6", "N_1TO6", "N_1TO6_BASE")
write.table(comparison, file.path(output_root, "comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
plot_columns <- c("CHISQ_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
plot_labels <- c("15,000-case matched 1:6", "Full sample + PCs", "Full sample (No PCs)")
plot_full <- draw_original_style(comparison, plot_columns, plot_labels, max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE))
plot_zoom <- draw_original_style(comparison, plot_columns, plot_labels, 10)
ggsave(file.path(output_root, "qq_15000_case_1to6_original_style_full.pdf"), plot_full, width = 13, height = 4.2, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_15000_case_1to6_original_style_zoom10.pdf"), plot_zoom, width = 13, height = 4.2, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_15000_case_1to6_original_style_zoom10.png"), plot_zoom, width = 13, height = 4.2, units = "in", dpi = 220)

summary <- data.frame(REGION = region, PHENOTYPE = phenotype, N_CASE_SELECTED = target_cases, N_CONTROL_AVAILABLE = length(control_available), N_CONTROL_MATCHED = matched_summary$N_CONTROL_MATCHED, N_STRATA = matched_summary$N_STRATA, MAX_ABS_SMD_PRE = matched_summary$MAX_ABS_SMD_PRE, MAX_ABS_SMD_POST = matched_summary$MAX_ABS_SMD_POST, N_CLR_OK = sum(clr$STATUS == "ok"), LAMBDA_1TO6 = median(comparison$CHISQ_1TO6, na.rm = TRUE) / qchisq(.5, 1), LAMBDA_PC_ADJUSTED = median(comparison$CHISQ_PC_ADJUSTED, na.rm = TRUE) / qchisq(.5, 1), LAMBDA_PLAIN = median(comparison$CHISQ_PLAIN, na.rm = TRUE) / qchisq(.5, 1), stringsAsFactors = FALSE)
write.table(summary, file.path(output_root, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("LCT_region / TV 15,000-case 1:6 matched CLR", "PC1-PC5-only matching; exact without-replacement 1:6; no sex; no permutation", paste0("Selected cases = ", target_cases, "; matched controls = ", matched_summary$N_CONTROL_MATCHED, "; unused available controls = ", length(control_available) - matched_summary$N_CONTROL_MATCHED), "SMD is computed for PC1-PC5 before and after matching; MAX_ABS_SMD is the largest absolute PC SMD.", "Matched CLR uses survival::clogit(method='exact') and reports Wald chi-square(1) p-values."), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message_time("Completed 15,000-case TV 1:6 analysis: ", output_root)
