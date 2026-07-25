#!/usr/bin/env Rscript

# LCT_region / TV: common 15,000-case PC-only matched CLR comparison.
# This script runs the missing 1:1 and 1:4 ratios and reuses the already
# completed 15,000-case 1:6 result, so all three ratios share the same cases.
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
output_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool")
matching_root <- file.path(output_root, "matching")
old_1to6_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_case_ratio6")
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
phenotype <- "TV"
region <- "LCT_region"
case_seed <- 3101L
target_cases <- 15000L
ratios <- c(1L, 4L)

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
if (length(case_available) < target_cases || length(control_available) < target_cases * 6L) stop("Insufficient TV cases or controls")
set.seed(case_seed)
selected_case_pool <- sample(case_available, target_cases, replace = FALSE)
selected_case_iids <- pc_data$IID[selected_case_pool]
write.table(data.frame(IID = selected_case_iids), file.path(output_root, "selected_case_pool.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message_time("Selected ", target_cases, " common cases; TV controls available = ", length(control_available))

matching_summaries <- list()
for (ratio in ratios) {
  start <- Sys.time()
  message_time("Starting TV ratio 1:", ratio, " matching")
  result <- match_ratio(ratio, selected_case_pool)
  if (any(grepl("sex", names(result$data), ignore.case = TRUE))) stop("Sex entered matched data")
  matching_summaries[[paste0("1to", ratio)]] <- result$summary
  message_time("Finished TV ratio 1:", ratio, " matching in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), start, units = "mins"))), " minutes; post max |SMD| = ", signif(result$summary$MAX_ABS_SMD_POST, 6))
  rm(result)
  gc(verbose = FALSE)
}

# Reuse the completed 15,000-case 1:6 result after verifying that its case pool
# is identical to the seeded common pool used above.
old_pool <- read.delim(file.path(old_1to6_root, "selected_case_pool.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!setequal(old_pool$IID, selected_case_iids)) stop("Existing 15,000-case 1:6 pool does not match the new common pool")
old6_summary <- read.delim(file.path(old_1to6_root, "matching_summary.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(old6_summary) != 1L || old6_summary$RATIO != 6L) stop("Existing 1:6 summary is not the expected 15,000-case result")
matching_summaries[["1to6"]] <- old6_summary
dir.create(file.path(matching_root, "ratio_1to6"), recursive = TRUE, showWarnings = FALSE)
old6_match_dir <- file.path(old_1to6_root, "matching", "ratio_1to6", "ratio_1to6")
for (f in c("matched_sets.tsv", "pc_balance.tsv", "matching_summary.tsv")) {
  if (!file.copy(file.path(old6_match_dir, f), file.path(matching_root, "ratio_1to6", f), overwrite = TRUE)) stop("Failed to copy existing 1:6 matching file: ", f)
}

matching_summary <- do.call(rbind, matching_summaries)
write.table(matching_summary, file.path(output_root, "matching_summary_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

genotype_frame <- read_genotypes(file.path(data_dir, "LCT_region_geno.raw"))
variant_names <- setdiff(names(genotype_frame), "IID")
genotype_index <- match(pc_data$IID, genotype_frame$IID)
if (anyNA(genotype_index)) stop("Some PC IIDs are missing from LCT genotype data")
genotype_all <- as.matrix(genotype_frame[genotype_index, variant_names, drop = FALSE])
storage.mode(genotype_all) <- "double"

plain <- read.delim(file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison", "TV", "logistic_plain.tsv"), check.names = FALSE)
pc <- read.delim(file.path(project_root, "results", "00_CLR", "full_sample_logistic_gwas", region, phenotype, "logistic_gwas.tsv"), check.names = FALSE)
plain_ok <- plain[plain$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(plain_ok) <- c("SNP", "CHISQ_PLAIN", "P_PLAIN", "N_PLAIN")
pc_ok <- pc[pc$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
comparison <- merge(plain_ok, pc_ok, by = "SNP", sort = FALSE)
clr_tables <- list()
for (ratio in c(1L, 4L)) {
  matched <- read.delim(file.path(matching_root, paste0("ratio_1to", ratio), "matched_sets.tsv"), check.names = FALSE)
  counts <- table(matched$STRATUM, matched$case)
  if (!all(c("0", "1") %in% colnames(counts)) || any(counts[, "1"] != 1L) || any(counts[, "0"] != ratio) || anyDuplicated(matched$IID)) stop("Invalid matched sets for ratio 1:", ratio)
  idx <- match(matched$IID, pc_data$IID)
  if (anyNA(idx) || any(as.character(matched$IID) != as.character(pc_data$IID[idx]))) stop("Matched IID alignment failed for ratio 1:", ratio)
  message_time("Running TV ratio 1:", ratio, " CLR")
  clr <- run_clr(genotype_all[idx, , drop = FALSE], matched, ratio)
  ratio_dir <- file.path(output_root, paste0("ratio_1to", ratio))
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(clr, file.path(ratio_dir, "matched_clr_gwas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  clr_tables[[paste0("1to", ratio)]] <- clr
  ok <- clr[clr$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]
  names(ok) <- c("SNP", paste0("CHISQ_1TO", ratio), paste0("P_1TO", ratio), paste0("N_1TO", ratio), paste0("N_1TO", ratio, "_BASE"))
  comparison <- merge(comparison, ok, by = "SNP", sort = FALSE)
  rm(matched, idx, clr, ok)
  gc(verbose = FALSE)
}

old6_clr <- read.delim(file.path(old_1to6_root, "matched_clr_gwas.tsv"), check.names = FALSE)
dir.create(file.path(output_root, "ratio_1to6"), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(file.path(old_1to6_root, "matched_clr_gwas.tsv"), file.path(output_root, "ratio_1to6", "matched_clr_gwas.tsv"), overwrite = TRUE)) stop("Failed to copy existing 1:6 CLR output")
ok6 <- old6_clr[old6_clr$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]
names(ok6) <- c("SNP", "CHISQ_1TO6", "P_1TO6", "N_1TO6", "N_1TO6_BASE")
comparison <- merge(comparison, ok6, by = "SNP", sort = FALSE)
write.table(comparison, file.path(output_root, "comparison_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_columns <- c("CHISQ_1TO1", "CHISQ_1TO4", "CHISQ_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
plot_labels <- c("15,000-case matched 1:1", "15,000-case matched 1:4", "15,000-case matched 1:6", "Full sample + PCs", "Full sample (No PCs)")
full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
plot_full <- draw_original_style(comparison, plot_columns, plot_labels, full_max)
plot_zoom <- draw_original_style(comparison, plot_columns, plot_labels, 10)
ggsave(file.path(output_root, "qq_15000_common_pool_1to1_1to4_1to6_original_style_full.pdf"), plot_full, width = 19, height = 4.2, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_15000_common_pool_1to1_1to4_1to6_original_style_zoom10.pdf"), plot_zoom, width = 19, height = 4.2, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_15000_common_pool_1to1_1to4_1to6_original_style_zoom10.png"), plot_zoom, width = 19, height = 4.2, units = "in", dpi = 220)

get_summary_value <- function(ratio, col) matching_summary[[col]][matching_summary$RATIO == ratio][1L]
summary <- data.frame(REGION = region, PHENOTYPE = phenotype, COMMON_CASE_POOL = target_cases, N_CONTROL_AVAILABLE = length(control_available), N_COMMON_MODELS = nrow(comparison), N_CASE_1TO1 = get_summary_value(1, "N_CASE_MATCHED"), N_CONTROL_1TO1 = get_summary_value(1, "N_CONTROL_MATCHED"), N_CASE_1TO4 = get_summary_value(4, "N_CASE_MATCHED"), N_CONTROL_1TO4 = get_summary_value(4, "N_CONTROL_MATCHED"), N_CASE_1TO6 = get_summary_value(6, "N_CASE_MATCHED"), N_CONTROL_1TO6 = get_summary_value(6, "N_CONTROL_MATCHED"), N_1TO1_OK = sum(read.delim(file.path(output_root, "ratio_1to1", "matched_clr_gwas.tsv"), check.names = FALSE)$STATUS == "ok"), N_1TO4_OK = sum(read.delim(file.path(output_root, "ratio_1to4", "matched_clr_gwas.tsv"), check.names = FALSE)$STATUS == "ok"), N_1TO6_OK = sum(old6_clr$STATUS == "ok"), MAX_ABS_SMD_POST_1TO1 = get_summary_value(1, "MAX_ABS_SMD_POST"), MAX_ABS_SMD_POST_1TO4 = get_summary_value(4, "MAX_ABS_SMD_POST"), MAX_ABS_SMD_POST_1TO6 = get_summary_value(6, "MAX_ABS_SMD_POST"), LAMBDA_1TO1 = median(comparison$CHISQ_1TO1) / qchisq(.5, 1), LAMBDA_1TO4 = median(comparison$CHISQ_1TO4) / qchisq(.5, 1), LAMBDA_1TO6 = median(comparison$CHISQ_1TO6) / qchisq(.5, 1), LAMBDA_PC_ADJUSTED = median(comparison$CHISQ_PC_ADJUSTED) / qchisq(.5, 1), LAMBDA_PLAIN = median(comparison$CHISQ_PLAIN) / qchisq(.5, 1), stringsAsFactors = FALSE)
write.table(summary, file.path(output_root, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("LCT_region / TV 15,000-case common-pool matched CLR comparison", "PC1-PC5-only matching; exact without-replacement 1:1, 1:4, and 1:6; no sex; no permutation", paste0("The same seeded 15,000-case pool was used for all ratios; 1:6 reuses the verified completed 15,000-case result."), "Matched CLR uses survival::clogit(method='exact') and reports Wald chi-square(1) p-values.", "The prior 16,564-case maximum common-pool output is intentionally not included here."), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message_time("Completed 15,000-case TV common-pool comparison: ", output_root)
