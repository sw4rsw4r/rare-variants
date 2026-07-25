#!/usr/bin/env Rscript

# Side-by-side display of the same 15,000-case rare-variant matched analyses:
# permutation-calibrated matched score tests versus asymptotic matched CLR
# Wald chi-square results. Full-sample logistic results are references.

options(stringsAsFactors = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required")
suppressPackageStartupMessages(library(ggplot2))

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

perm_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool_rare_permutation")
chisq_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool_rare")
output_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool_rare_perm_vs_chisq")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
load_function_definitions <- function(path, names_to_load) {
  parsed <- parse(path)
  for (expr in parsed) {
    if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
      name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (length(name) == 1L && name %in% names_to_load) eval(expr, envir = parent.frame())
    }
  }
}
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "021_lct_tv_15000_rare_permutation_matched.R"), c("p_lambda", "chisq_lambda", "draw_with_lambdas"))

perm <- read_tab(file.path(perm_root, "comparison_rare_permutation_1to1_1to4_1to6.tsv"))
chisq <- read_tab(file.path(chisq_root, "comparison_rare_1to1_1to4_1to6.tsv"))
merged <- merge(perm, chisq, by = "SNP", suffixes = c("_PERM", "_CHISQ"), sort = FALSE)
if (nrow(merged) == 0L) stop("No common rare variants between permutation and chi-square results")

comparison <- data.frame(SNP = merged$SNP, stringsAsFactors = FALSE)
for (ratio in c(1L, 4L, 6L)) {
  comparison[[paste0("P_PERM_1TO", ratio)]] <- merged[[paste0("P_PERM_1TO", ratio)]]
  comparison[[paste0("CHISQ_PERM_1TO", ratio)]] <- merged[[paste0("CHISQ_PERM_1TO", ratio)]]
  comparison[[paste0("P_WALD_1TO", ratio)]] <- merged[[paste0("P_1TO", ratio)]]
  comparison[[paste0("CHISQ_WALD_1TO", ratio)]] <- merged[[paste0("CHISQ_1TO", ratio)]]
}
comparison$P_PC_ADJUSTED <- merged$P_PC_ADJUSTED_PERM
comparison$CHISQ_PC_ADJUSTED <- merged$CHISQ_PC_ADJUSTED_PERM
comparison$P_PLAIN <- merged$P_PLAIN_PERM
comparison$CHISQ_PLAIN <- merged$CHISQ_PLAIN_PERM
write.table(comparison, file.path(output_root, "comparison_rare_permutation_vs_chisq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_columns <- c("CHISQ_PERM_1TO1", "CHISQ_WALD_1TO1", "CHISQ_PERM_1TO4", "CHISQ_WALD_1TO4", "CHISQ_PERM_1TO6", "CHISQ_WALD_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
plot_labels <- c("1:1 permutation", "1:1 chi-square (Wald)", "1:4 permutation", "1:4 chi-square (Wald)", "1:6 permutation", "1:6 chi-square (Wald)", "Full sample + PCs", "Full sample (No PCs)")
lambda_values <- c(p_lambda(comparison$P_PERM_1TO1), chisq_lambda(comparison$CHISQ_WALD_1TO1), p_lambda(comparison$P_PERM_1TO4), chisq_lambda(comparison$CHISQ_WALD_1TO4), p_lambda(comparison$P_PERM_1TO6), chisq_lambda(comparison$CHISQ_WALD_1TO6), chisq_lambda(comparison$CHISQ_PC_ADJUSTED), chisq_lambda(comparison$CHISQ_PLAIN))
lambda_types <- c("p", "chisq", "p", "chisq", "p", "chisq", "chisq", "chisq")
full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
plot_full <- draw_with_lambdas(comparison, plot_columns, plot_labels, full_max, lambda_values, lambda_types)
plot_zoom <- draw_with_lambdas(comparison, plot_columns, plot_labels, 10, lambda_values, lambda_types)
ggsave(file.path(output_root, "qq_rare_15000_permutation_vs_chisq_full.pdf"), plot_full, width = 30, height = 4.5, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_rare_15000_permutation_vs_chisq_zoom10.pdf"), plot_zoom, width = 30, height = 4.5, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_rare_15000_permutation_vs_chisq_zoom10.png"), plot_zoom, width = 30, height = 4.5, units = "in", dpi = 220)

summary <- data.frame(N_COMMON_VARIANTS = nrow(comparison), LAMBDA_P_PERM_1TO1 = p_lambda(comparison$P_PERM_1TO1), LAMBDA_CHISQ_WALD_1TO1 = chisq_lambda(comparison$CHISQ_WALD_1TO1), LAMBDA_P_PERM_1TO4 = p_lambda(comparison$P_PERM_1TO4), LAMBDA_CHISQ_WALD_1TO4 = chisq_lambda(comparison$CHISQ_WALD_1TO4), LAMBDA_P_PERM_1TO6 = p_lambda(comparison$P_PERM_1TO6), LAMBDA_CHISQ_WALD_1TO6 = chisq_lambda(comparison$CHISQ_WALD_1TO6), LAMBDA_PC_ADJUSTED = chisq_lambda(comparison$CHISQ_PC_ADJUSTED), LAMBDA_PLAIN = chisq_lambda(comparison$CHISQ_PLAIN), stringsAsFactors = FALSE)
write.table(summary, file.path(output_root, "summary_rare_permutation_vs_chisq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("Side-by-side rare-variant comparison using the same 15,000-case matched sets", "Permutation panels use lambda_p = 0.5 / median(P_PERM); chi-square panels use conventional median(CHISQ)/qchisq(0.5,1).", "Full-sample plain and PC-adjusted panels are asymptotic logistic references.", "All panels contain only rare variants successful in every compared model."), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message(sprintf("Created rare permutation-versus-chi-square comparison: %s; variants=%d", output_root, nrow(comparison)))
