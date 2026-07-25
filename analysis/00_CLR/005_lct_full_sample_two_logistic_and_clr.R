#!/usr/bin/env Rscript

# Targeted LCT comparison: full-sample plain logistic, full-sample PC-only
# logistic, and the already-validated PC-only matched CLR. Only TV and
# qualification are processed. Sex and permutation are not used.

options(stringsAsFactors = FALSE)
script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

data_dir <- file.path(project_root, "data")
pc_root <- file.path(project_root, "results", "00_CLR", "full_sample_logistic_gwas")
matched_root <- file.path(project_root, "results", "00_CLR", "pc_adjusted_vs_matched_clr")
output_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
region <- "LCT_region"
phenotypes <- c("TV", "qualification")

message_time <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

# Load only the tested helper functions from 003 without executing its
# all-region top-level analysis.
helper_script <- parse(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"))
helper_names <- c("read_header", "read_selected_tsv", "read_genotypes", "load_phenotype_ids", "make_binary_phenotype", "logistic_wald_gwas")
for (expr in helper_script) {
  if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
    name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
    if (length(name) == 1L && name %in% helper_names) eval(expr, envir = environment())
  }
}
if (!all(vapply(helper_names, exists, logical(1L), envir = environment(), inherits = FALSE))) {
  stop("Could not load all GWAS helper functions")
}

read_pc_adjusted <- function(phenotype) {
  path <- file.path(pc_root, region, phenotype, "logistic_gwas.tsv")
  if (!file.exists(path)) stop("Missing PC-adjusted result: ", path)
  result <- read.delim(path, check.names = FALSE)
  if (!all(result$COVARIATES == paste(pc_names, collapse = "+"))) stop("Unexpected covariates in ", path)
  result
}

read_matched_clr <- function(phenotype) {
  path <- file.path(matched_root, region, phenotype, "matched_clr_gwas.tsv")
  if (!file.exists(path)) stop("Missing matched CLR result: ", path)
  read.delim(path, check.names = FALSE)
}

draw_qq <- function(chisq, model_name, subtitle, axis_max) {
  observed <- sort(chisq[is.finite(chisq)])
  expected <- qchisq(ppoints(length(observed)), df = 1)
  lambda_gc <- median(observed) / qchisq(0.5, df = 1)
  plot(expected, observed, pch = 19, cex = 0.65,
       col = grDevices::adjustcolor("#2166AC", alpha.f = 0.65),
       xlim = c(0, axis_max), ylim = c(0, axis_max),
       xlab = "Expected chi-square(1) quantile",
       ylab = "Observed Wald chi-square", main = model_name)
  abline(a = 0, b = 1, lty = 2, lwd = 1.5, col = "#B2182B")
  legend("topleft", legend = c(sprintf("lambda[GC] = %.3f", lambda_gc),
                                sprintf("variants = %d", length(observed)), subtitle),
         bty = "n", cex = 0.8)
  invisible(lambda_gc)
}

draw_two <- function(comparison, phenotype) {
  axis_max <- max(c(qchisq(ppoints(nrow(comparison)), 1), comparison$CHISQ_UNADJUSTED, comparison$CHISQ_PC_ADJUSTED), na.rm = TRUE)
  old_par <- par(mfrow = c(1, 2), mar = c(4.7, 4.7, 4.2, 1.2), oma = c(0, 0, 2.1, 0))
  l1 <- draw_qq(comparison$CHISQ_UNADJUSTED, "Full sample: plain logistic", sprintf("N = %d", comparison$N_UNADJUSTED[1L]), axis_max)
  l2 <- draw_qq(comparison$CHISQ_PC_ADJUSTED, "Full sample: PC-adjusted logistic", sprintf("N = %d", comparison$N_PC_ADJUSTED[1L]), axis_max)
  mtext(paste(region, phenotype, sep = " / "), outer = TRUE, side = 3, line = 0.3, font = 2, cex = 1.15)
  par(old_par)
  c(unadjusted = l1, pc_adjusted = l2)
}

draw_three <- function(comparison, phenotype) {
  axis_max <- max(c(qchisq(ppoints(nrow(comparison)), 1), comparison$CHISQ_UNADJUSTED, comparison$CHISQ_PC_ADJUSTED, comparison$CHISQ_MATCHED_CLR), na.rm = TRUE)
  old_par <- par(mfrow = c(1, 3), mar = c(4.7, 4.5, 4.2, 1.0), oma = c(0, 0, 2.1, 0))
  l1 <- draw_qq(comparison$CHISQ_UNADJUSTED, "Plain logistic", sprintf("N = %d", comparison$N_UNADJUSTED[1L]), axis_max)
  l2 <- draw_qq(comparison$CHISQ_PC_ADJUSTED, "PC-adjusted logistic", sprintf("N = %d", comparison$N_PC_ADJUSTED[1L]), axis_max)
  l3 <- draw_qq(comparison$CHISQ_MATCHED_CLR, "PC-matched 1:4 CLR", sprintf("baseline N = %d", comparison$N_MATCHED_CLR_BASE[1L]), axis_max)
  mtext(paste(region, phenotype, sep = " / "), outer = TRUE, side = 3, line = 0.3, font = 2, cex = 1.15)
  par(old_par)
  c(unadjusted = l1, pc_adjusted = l2, matched_clr = l3)
}

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
rm(genotype_frame, genotype_index)

all_summary <- list()
all_panels <- list()

for (phenotype in phenotypes) {
  message_time("Running plain and PC-adjusted full-sample GWAS for ", phenotype)
  y_all <- make_binary_phenotype(pc_data$IID, phenotype)
  valid <- !is.na(y_all)
  y <- y_all[valid]
  genotype <- genotype_all[valid, , drop = FALSE]
  pc_covariates <- pc_data[valid, pc_names, drop = FALSE]
  no_covariates <- data.frame(row.names = seq_len(length(y)))

  unadjusted <- logistic_wald_gwas(genotype, no_covariates, y)
  pc_adjusted <- logistic_wald_gwas(genotype, pc_covariates, y)
  for (obj_name in c("unadjusted", "pc_adjusted")) {
    obj <- get(obj_name)
    obj$REGION <- region
    obj$PHENOTYPE <- phenotype
    obj$MODEL <- if (obj_name == "unadjusted") "FULL_SAMPLE_PLAIN_LOGISTIC" else "FULL_SAMPLE_PC_ADJUSTED_LOGISTIC"
    obj <- obj[, c("REGION", "PHENOTYPE", "MODEL", setdiff(names(obj), c("REGION", "PHENOTYPE", "MODEL"))), drop = FALSE]
    assign(obj_name, obj)
  }

  matched <- read_matched_clr(phenotype)
  matched_ok <- matched[matched$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]
  names(matched_ok) <- c("SNP", "CHISQ_MATCHED_CLR", "P_MATCHED_CLR", "N_MATCHED_CLR", "N_MATCHED_CLR_BASE")
  plain_ok <- unadjusted[unadjusted$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]
  names(plain_ok) <- c("SNP", "CHISQ_UNADJUSTED", "P_UNADJUSTED", "N_UNADJUSTED")
  pc_ok <- pc_adjusted[pc_adjusted$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]
  names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
  comparison <- merge(merge(plain_ok, pc_ok, by = "SNP"), matched_ok, by = "SNP")
  if (nrow(comparison) == 0L) stop("No common SNPs for ", phenotype)

  result_dir <- file.path(output_root, phenotype)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(unadjusted, file.path(result_dir, "logistic_plain.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(pc_adjusted, file.path(result_dir, "logistic_pc_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(comparison, file.path(result_dir, "plain_vs_pc_adjusted_vs_matched_clr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  pdf(file.path(result_dir, "qq_full_sample_plain_vs_pc_adjusted.pdf"), width = 11, height = 5.5, useDingbats = FALSE)
  lambdas_two <- draw_two(comparison, phenotype)
  dev.off()
  png(file.path(result_dir, "qq_full_sample_plain_vs_pc_adjusted.png"), width = 2640, height = 1320, res = 220)
  draw_two(comparison, phenotype)
  dev.off()

  pdf(file.path(result_dir, "qq_plain_pc_adjusted_matched_clr.pdf"), width = 16.5, height = 5.5, useDingbats = FALSE)
  lambdas_three <- draw_three(comparison, phenotype)
  dev.off()
  png(file.path(result_dir, "qq_plain_pc_adjusted_matched_clr.png"), width = 3960, height = 1320, res = 220)
  draw_three(comparison, phenotype)
  dev.off()

  all_panels[[phenotype]] <- list(comparison = comparison, phenotype = phenotype)
  all_summary[[phenotype]] <- data.frame(
    REGION = region, PHENOTYPE = phenotype, N_FULL_SAMPLE = length(y),
    N_PLAIN_OK = sum(unadjusted$STATUS == "ok"), N_PC_ADJUSTED_OK = sum(pc_adjusted$STATUS == "ok"),
    N_MATCHED_CLR_OK = sum(matched$STATUS == "ok"), N_COMMON_ALL_THREE = nrow(comparison),
    LAMBDA_GC_PLAIN = unname(lambdas_two["unadjusted"]), LAMBDA_GC_PC_ADJUSTED = unname(lambdas_two["pc_adjusted"]),
    LAMBDA_GC_MATCHED_CLR = unname(lambdas_three["matched_clr"]), MIN_P_PLAIN = min(comparison$P_UNADJUSTED),
    MIN_P_PC_ADJUSTED = min(comparison$P_PC_ADJUSTED), MIN_P_MATCHED_CLR = min(comparison$P_MATCHED_CLR),
    stringsAsFactors = FALSE
  )
  message_time("Saved ", result_dir)
  rm(y_all, valid, y, genotype, pc_covariates, no_covariates, unadjusted, pc_adjusted, matched, comparison)
  gc(verbose = FALSE)
}

analysis_summary <- do.call(rbind, all_summary)
write.table(analysis_summary, file.path(output_root, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(output_root, "qq_full_sample_plain_vs_pc_adjusted_all.pdf"), width = 11, height = 5.5, onefile = TRUE, useDingbats = FALSE)
for (panel in all_panels) draw_two(panel$comparison, panel$phenotype)
dev.off()
pdf(file.path(output_root, "qq_plain_pc_adjusted_matched_clr_all.pdf"), width = 16.5, height = 5.5, onefile = TRUE, useDingbats = FALSE)
for (panel in all_panels) draw_three(panel$comparison, panel$phenotype)
dev.off()

writeLines(c(
  "LCT_region TV and qualification comparison",
  "Full-sample plain logistic: case ~ SNP",
  "Full-sample PC-adjusted logistic: case ~ PC1 + PC2 + PC3 + PC4 + PC5 + SNP",
  "Matched CLR: PC-only 1:4 matching result from 004_run_pc_matched_clr_comparison.R",
  "Sex: not read, matched, or modeled",
  "Permutation: not used",
  "QQ reference: theoretical chi-square(1). The two-panel PDF contains the two full-sample models; the three-panel PDF adds matched CLR."
), file.path(output_root, "METHODS.txt"))
message_time("Completed selected LCT comparisons: ", output_root)
