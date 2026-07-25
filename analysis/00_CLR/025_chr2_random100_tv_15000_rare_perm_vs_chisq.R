#!/usr/bin/env Rscript

# Compact side-by-side comparison of the same chr2_random100 TV rare variants:
# within-stratum permutation calibration versus asymptotic matched CLR Wald
# chi-square.  Full-sample plain and PC-adjusted logistic panels are included
# as references.

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

perm_root <- file.path(project_root, "results", "00_CLR", "CHR2_random100_TV_15000_rare_permutation")
chisq_root <- file.path(project_root, "results", "00_CLR", "CHR2_random100_TV_15000_rare")
output_root <- file.path(project_root, "results", "00_CLR", "CHR2_random100_TV_15000_rare_perm_vs_chisq")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
p_lambda <- function(p) { p <- p[is.finite(p) & p >= 0 & p <= 1]; if (length(p) == 0L) NA_real_ else 0.5 / median(p) }
chisq_lambda <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else median(x) / qchisq(.5, 1) }

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

draw <- function(max_val) {
  rows <- lapply(seq_along(plot_columns), function(i) { v <- sort(comparison[[plot_columns[i]]][is.finite(comparison[[plot_columns[i]]])]); data.frame(exp_chisq = qchisq(ppoints(length(v)), 1), obs_chisq = v, Model = plot_labels[i], stringsAsFactors = FALSE) })
  qq <- do.call(rbind, rows); qq$Model <- factor(qq$Model, levels = plot_labels)
  anno <- data.frame(Model = factor(plot_labels, levels = plot_labels), label = ifelse(lambda_types == "p", sprintf("lambda_p = %.3f", lambda_values), sprintf("lambda_GC = %.3f", lambda_values)), x = max_val * .1, y = max_val * .9)
  ggplot(qq, aes(exp_chisq, obs_chisq)) + geom_point(size = 1, alpha = .6) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + facet_wrap(~Model, nrow = 1) + geom_text(data = anno, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) + labs(title = "chr2_random100 / TV rare variants: permutation vs chi-square", x = "Expected chi-square(1)", y = "Observed chi-square(1)") + theme_minimal() + scale_x_continuous(limits = c(0, max_val)) + scale_y_continuous(limits = c(0, max_val))
}
full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
ggsave(file.path(output_root, "qq_rare_15000_permutation_vs_chisq_full.pdf"), draw(full_max), width = 30, height = 4.5, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_rare_15000_permutation_vs_chisq_zoom10.pdf"), draw(10), width = 30, height = 4.5, units = "in", device = "pdf")
ggsave(file.path(output_root, "qq_rare_15000_permutation_vs_chisq_zoom10.png"), draw(10), width = 30, height = 4.5, units = "in", dpi = 220)

summary <- data.frame(
  REGION = "chr2_random100", PHENOTYPE = "TV", N_COMMON_VARIANTS = nrow(comparison),
  LAMBDA_P_PERM_1TO1 = p_lambda(comparison$P_PERM_1TO1), LAMBDA_CHISQ_WALD_1TO1 = chisq_lambda(comparison$CHISQ_WALD_1TO1),
  LAMBDA_P_PERM_1TO4 = p_lambda(comparison$P_PERM_1TO4), LAMBDA_CHISQ_WALD_1TO4 = chisq_lambda(comparison$CHISQ_WALD_1TO4),
  LAMBDA_P_PERM_1TO6 = p_lambda(comparison$P_PERM_1TO6), LAMBDA_CHISQ_WALD_1TO6 = chisq_lambda(comparison$CHISQ_WALD_1TO6),
  LAMBDA_PC_ADJUSTED = chisq_lambda(comparison$CHISQ_PC_ADJUSTED), LAMBDA_PLAIN = chisq_lambda(comparison$CHISQ_PLAIN), stringsAsFactors = FALSE
)
write.table(summary, file.path(output_root, "summary_rare_permutation_vs_chisq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "chr2_random100 / TV side-by-side rare-variant comparison",
  "Permutation panels use lambda_p = 0.5 / median(P_PERM); matched chi-square panels use lambda_GC = median(CHISQ)/qchisq(0.5,1).",
  "Full-sample plain and PC-adjusted panels are asymptotic logistic references.",
  "All panels contain rare variants successful in every compared model."
), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message(sprintf("Created chr2_random100 rare permutation-versus-chi-square comparison: %s; variants=%d", output_root, nrow(comparison)))
