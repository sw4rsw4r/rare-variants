#!/usr/bin/env Rscript

# Targeted LCT comparison with an additional PC-only 1:1 matched CLR.
#
# Existing scripts 001/002 are intentionally untouched.  This script only
# processes LCT_region / TV and qualification, and compares:
#   (1) full-sample plain logistic GWAS,
#   (2) full-sample PC-adjusted logistic GWAS,
#   (3) PC-only 1:4 matched CLR already computed by script 004, and
#   (4) a newly computed PC-only 1:1 matched CLR.
# No sex term and no permutation are used.

options(stringsAsFactors = FALSE)

if (!requireNamespace("MatchIt", quietly = TRUE)) stop("The MatchIt package is required")
if (!requireNamespace("survival", quietly = TRUE)) stop("The survival package is required")
suppressPackageStartupMessages(library(survival))

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
matching_root <- file.path(output_root, "matching_1to1")
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
num_controls <- 1L
max_cases <- 1000L
control_pool_multiplier <- 15L
region <- "LCT_region"
phenotypes <- c("TV", "qualification")
# Use the same seeded case/control candidate samples as script 004 so that
# the comparison isolates the matching ratio (1:4 versus 1:1).
phenotype_seeds <- c(TV = 1101L, qualification = 1102L)

message_time <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

# Read only helper definitions from 003; its all-region top-level analysis is
# not executed here.
helper_script <- parse(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"))
helper_names <- c("read_header", "read_selected_tsv", "read_genotypes", "load_phenotype_ids", "make_binary_phenotype")
for (expr in helper_script) {
  if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
    name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
    if (length(name) == 1L && name %in% helper_names) eval(expr, envir = environment())
  }
}
if (!all(vapply(helper_names, exists, logical(1L), envir = environment(), inherits = FALSE))) {
  stop("Could not load all genotype/phenotype helper functions")
}

compute_balance <- function(pre_data, matched_data) {
  rows <- lapply(pc_names, function(variable) {
    pre_case <- pre_data[pre_data$case == 1L, variable]
    pre_control <- pre_data[pre_data$case == 0L, variable]
    post_case <- matched_data[matched_data$case == 1L, variable]
    post_control <- matched_data[matched_data$case == 0L, variable]
    denominator <- sqrt((stats::var(pre_case) + stats::var(pre_control)) / 2)
    data.frame(
      VARIABLE = variable,
      PRE_MEAN_CASE = mean(pre_case),
      PRE_MEAN_CONTROL = mean(pre_control),
      PRE_SMD = (mean(pre_case) - mean(pre_control)) / denominator,
      POST_MEAN_CASE = mean(post_case),
      POST_MEAN_CONTROL = mean(post_control),
      POST_SMD = (mean(post_case) - mean(post_control)) / denominator,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

match_phenotype_1to1 <- function(pc_data, phenotype) {
  y <- make_binary_phenotype(pc_data$IID, phenotype)
  case_available <- which(y == 1L)
  control_available <- which(y == 0L)
  if (length(case_available) == 0L || length(control_available) < num_controls) {
    stop("Insufficient cases or controls for ", phenotype)
  }

  set.seed(unname(phenotype_seeds[phenotype]))
  selected_cases <- if (length(case_available) > max_cases) sample(case_available, max_cases, replace = FALSE) else case_available
  candidate_control_n <- min(
    length(control_available),
    max(num_controls * length(selected_cases), control_pool_multiplier * length(selected_cases))
  )
  selected_controls <- if (length(control_available) > candidate_control_n) {
    sample(control_available, candidate_control_n, replace = FALSE)
  } else control_available

  candidate_rows <- c(selected_cases, selected_controls)
  candidate_data <- pc_data[candidate_rows, c("IID", pc_names), drop = FALSE]
  candidate_data$case <- c(rep.int(1L, length(selected_cases)), rep.int(0L, length(selected_controls)))
  candidate_data <- candidate_data[, c("IID", "case", pc_names), drop = FALSE]
  rownames(candidate_data) <- candidate_data$IID

  formula_match <- stats::reformulate(pc_names, response = "case")
  match_object <- MatchIt::matchit(
    formula_match, data = candidate_data, method = "nearest", distance = "mahalanobis",
    estimand = "ATT", replace = FALSE, m.order = "farthest", ratio = num_controls
  )
  matched_data <- MatchIt::match.data(match_object, data = candidate_data, drop.unmatched = TRUE)
  matched_data$STRATUM <- match(as.character(matched_data$subclass), unique(as.character(matched_data$subclass)))
  matched_data <- matched_data[, c("IID", "case", pc_names, "STRATUM", "weights"), drop = FALSE]
  matched_data <- matched_data[order(matched_data$STRATUM, -matched_data$case), , drop = FALSE]
  rownames(matched_data) <- NULL

  set_counts <- table(matched_data$STRATUM, matched_data$case)
  if (!all(c("0", "1") %in% colnames(set_counts)) || any(set_counts[, "1"] != 1L) || any(set_counts[, "0"] != 1L)) {
    stop("Matching did not produce complete 1:1 sets for ", phenotype)
  }
  if (anyDuplicated(matched_data$IID)) stop("A participant was reused in matching for ", phenotype)

  balance <- compute_balance(candidate_data, matched_data)
  summary <- data.frame(
    PHENOTYPE = phenotype,
    MATCH_VARIABLES = paste(pc_names, collapse = "+"),
    RATIO = num_controls,
    MAX_CASES = max_cases,
    CONTROL_POOL_MULTIPLIER = control_pool_multiplier,
    N_CASE_AVAILABLE = length(case_available),
    N_CASE_SELECTED = length(selected_cases),
    N_CONTROL_AVAILABLE = length(control_available),
    N_CONTROL_CANDIDATE = length(selected_controls),
    N_CASE_MATCHED = sum(matched_data$case == 1L),
    N_CONTROL_MATCHED = sum(matched_data$case == 0L),
    N_STRATA = length(unique(matched_data$STRATUM)),
    MAX_ABS_SMD_PRE = max(abs(balance$PRE_SMD)),
    MAX_ABS_SMD_POST = max(abs(balance$POST_SMD)),
    stringsAsFactors = FALSE
  )
  phenotype_dir <- file.path(matching_root, phenotype)
  dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(matched_data, file.path(phenotype_dir, "matched_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(balance, file.path(phenotype_dir, "pc_balance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary, file.path(phenotype_dir, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(data = matched_data, balance = balance, summary = summary)
}

run_clogit_gwas_1to1 <- function(genotype, matched_data) {
  variant_names <- colnames(genotype)
  n_variants <- ncol(genotype)
  result <- data.frame(
    SNP = variant_names,
    N_MATCHED_BASE = nrow(matched_data), N_STRATA_BASE = length(unique(matched_data$STRATUM)),
    N = integer(n_variants), N_CASE = integer(n_variants), N_CONTROL = integer(n_variants),
    N_STRATA = integer(n_variants), N_STRATA_COMPLETE_1TO1 = integer(n_variants),
    N_STRATA_REDUCED = integer(n_variants), N_INFORMATIVE_STRATA = integer(n_variants),
    N_GENOTYPE_MISSING = integer(n_variants), EAF = rep.int(NA_real_, n_variants),
    MAF = rep.int(NA_real_, n_variants), MAC = rep.int(NA_real_, n_variants),
    BETA = rep.int(NA_real_, n_variants), SE = rep.int(NA_real_, n_variants),
    Z = rep.int(NA_real_, n_variants), CHISQ = rep.int(NA_real_, n_variants),
    P = rep.int(NA_real_, n_variants), ITERATIONS = integer(n_variants),
    STATUS = rep.int("not_tested", n_variants), WARNING = rep.int("", n_variants),
    TEST = rep.int("MATCHED_CLR_1TO1_WALD_CHISQ1", n_variants),
    COVARIATES = rep.int("NONE_PC_MATCHED", n_variants), stringsAsFactors = FALSE
  )

  for (j in seq_len(n_variants)) {
    genotype_j <- genotype[, j]
    nonmissing <- is.finite(genotype_j)
    result$N_GENOTYPE_MISSING[j] <- sum(!nonmissing)
    case_present <- tapply(nonmissing & matched_data$case == 1L, matched_data$STRATUM, any)
    control_nonmissing <- tapply(nonmissing & matched_data$case == 0L, matched_data$STRATUM, sum)
    valid_strata <- as.integer(names(case_present)[case_present & control_nonmissing >= 1L])
    use <- nonmissing & matched_data$STRATUM %in% valid_strata
    if (!any(use)) { result$STATUS[j] <- "no_valid_strata"; next }

    model_data <- data.frame(case = matched_data$case[use], SNP = genotype_j[use], STRATUM = matched_data$STRATUM[use])
    controls_per_stratum <- table(model_data$STRATUM[model_data$case == 0L])
    informative <- tapply(model_data$SNP, model_data$STRATUM, function(value) length(unique(value)) > 1L)
    result$N[j] <- nrow(model_data); result$N_CASE[j] <- sum(model_data$case == 1L); result$N_CONTROL[j] <- sum(model_data$case == 0L)
    result$N_STRATA[j] <- length(unique(model_data$STRATUM)); result$N_STRATA_COMPLETE_1TO1[j] <- sum(controls_per_stratum == 1L)
    result$N_STRATA_REDUCED[j] <- sum(controls_per_stratum < 1L); result$N_INFORMATIVE_STRATA[j] <- sum(informative)
    allele_sum <- sum(model_data$SNP); result$EAF[j] <- allele_sum / (2 * nrow(model_data)); result$MAF[j] <- min(result$EAF[j], 1 - result$EAF[j]); result$MAC[j] <- min(allele_sum, 2 * nrow(model_data) - allele_sum)
    if (result$N_INFORMATIVE_STRATA[j] == 0L) { result$STATUS[j] <- "no_within_stratum_variation"; next }

    warning_text <- character()
    fit <- tryCatch(withCallingHandlers(
      survival::clogit(case ~ SNP + strata(STRATUM), data = model_data, method = "exact", control = survival::coxph.control(iter.max = 50L, eps = 1e-9)),
      warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }
    ), error = function(e) e)
    if (inherits(fit, "error")) { result$STATUS[j] <- "fit_error"; result$WARNING[j] <- conditionMessage(fit); next }
    result$ITERATIONS[j] <- if (length(fit$iter) > 0L) max(fit$iter) else NA_integer_; result$WARNING[j] <- paste(unique(warning_text), collapse = " | ")
    if (any(grepl("infinite|converged before", warning_text, ignore.case = TRUE))) { result$STATUS[j] <- "possible_separation"; next }
    coefficient_table <- summary(fit)$coefficients
    if (!("SNP" %in% rownames(coefficient_table))) { result$STATUS[j] <- "snp_coefficient_missing"; next }
    beta <- coefficient_table["SNP", "coef"]; se <- coefficient_table["SNP", "se(coef)"]
    if (!is.finite(beta) || !is.finite(se) || se <= 0) { result$STATUS[j] <- "invalid_coefficient"; next }
    z_value <- beta / se; result$BETA[j] <- beta; result$SE[j] <- se; result$Z[j] <- z_value; result$CHISQ[j] <- z_value^2; result$P[j] <- pchisq(result$CHISQ[j], df = 1, lower.tail = FALSE); result$STATUS[j] <- "ok"
  }
  result
}

draw_qq_zoom <- function(chisq, model_name, subtitle, zoom_max = 10) {
  observed <- sort(chisq[is.finite(chisq)])
  expected <- qchisq(ppoints(length(observed)), df = 1)
  lambda_gc <- median(observed) / qchisq(0.5, df = 1)
  plot(expected, observed, pch = 19, cex = 0.65, col = grDevices::adjustcolor("#2166AC", alpha.f = 0.7),
       xlim = c(0, zoom_max), ylim = c(0, zoom_max), xlab = "Expected chi-square(1) quantile (zoom 0-10)",
       ylab = "Observed Wald chi-square (zoom 0-10)", main = model_name)
  abline(a = 0, b = 1, lty = 2, lwd = 1.5, col = "#B2182B")
  legend("topleft", legend = c(sprintf("lambda[GC] = %.3f", lambda_gc), sprintf("variants = %d", length(observed)), subtitle, ">10 values clipped"), bty = "n", cex = 0.76)
  invisible(lambda_gc)
}

draw_four_zoom <- function(comparison, phenotype) {
  old_par <- par(mfrow = c(1, 4), mar = c(4.9, 4.35, 4.2, 0.65), oma = c(0, 0, 2.1, 0))
  l1 <- draw_qq_zoom(comparison$CHISQ_UNADJUSTED, "Plain logistic", sprintf("N = %d", comparison$N_UNADJUSTED[1L]))
  l2 <- draw_qq_zoom(comparison$CHISQ_PC_ADJUSTED, "PC-adjusted logistic", sprintf("N = %d", comparison$N_PC_ADJUSTED[1L]))
  l3 <- draw_qq_zoom(comparison$CHISQ_MATCHED_CLR, "PC-matched 1:4 CLR", sprintf("baseline N = %d", comparison$N_MATCHED_CLR_BASE[1L]))
  l4 <- draw_qq_zoom(comparison$CHISQ_MATCHED_CLR_1TO1, "PC-matched 1:1 CLR", sprintf("baseline N = %d", comparison$N_MATCHED_CLR_1TO1_BASE[1L]))
  mtext(paste("LCT_region /", phenotype, "/ fixed zoom 0-10"), outer = TRUE, side = 3, line = 0.3, font = 2, cex = 1.15)
  par(old_par)
  c(unadjusted = l1, pc_adjusted = l2, matched_clr_1to4 = l3, matched_clr_1to1 = l4)
}

pc_data <- read_selected_tsv(file.path(data_dir, "UKBB_pca.eigenvec"), c("IID", pc_names), c("character", rep("numeric", num_pcs)))
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
if (anyDuplicated(pc_data$IID)) stop("Duplicate IID in PCA data")
genotype_frame <- read_genotypes(file.path(data_dir, "LCT_region_geno.raw"))
variant_names <- setdiff(names(genotype_frame), "IID")
genotype_index <- match(pc_data$IID, genotype_frame$IID)
if (anyNA(genotype_index)) stop("Some PC IIDs are missing from LCT genotype data")
genotype_all <- as.matrix(genotype_frame[genotype_index, variant_names, drop = FALSE]); storage.mode(genotype_all) <- "double"
stopifnot(identical(pc_data$IID, genotype_frame$IID[genotype_index]))
rm(genotype_frame, genotype_index)

match_results <- setNames(vector("list", length(phenotypes)), phenotypes)
for (phenotype in phenotypes) {
  message_time("Matching ", phenotype, " on PC1-PC5 only (1:1)")
  match_results[[phenotype]] <- match_phenotype_1to1(pc_data, phenotype)
  message_time("Matched ", match_results[[phenotype]]$summary$N_CASE_MATCHED, " cases; max post-match |SMD| = ", sprintf("%.4f", match_results[[phenotype]]$summary$MAX_ABS_SMD_POST))
}
matching_summary <- do.call(rbind, lapply(match_results, `[[`, "summary"))
write.table(matching_summary, file.path(matching_root, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

all_panels <- list(); all_summary <- list()
for (phenotype in phenotypes) {
  matched_data <- match_results[[phenotype]]$data
  genotype_index <- match(matched_data$IID, pc_data$IID)
  if (anyNA(genotype_index)) stop("Matched IID missing from PCA data")
  genotype <- genotype_all[genotype_index, , drop = FALSE]
  stopifnot(identical(matched_data$IID, pc_data$IID[genotype_index]))
  message_time("Running 1:1 matched CLR for ", phenotype)
  clr_1to1 <- run_clogit_gwas_1to1(genotype, matched_data)
  result_dir <- file.path(output_root, phenotype)
  write.table(clr_1to1, file.path(result_dir, "matched_clr_1to1_gwas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  old <- read.delim(file.path(result_dir, "plain_vs_pc_adjusted_vs_matched_clr.tsv"), check.names = FALSE)
  clr_ok <- clr_1to1[clr_1to1$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]
  names(clr_ok) <- c("SNP", "CHISQ_MATCHED_CLR_1TO1", "P_MATCHED_CLR_1TO1", "N_MATCHED_CLR_1TO1", "N_MATCHED_CLR_1TO1_BASE")
  comparison <- merge(old, clr_ok, by = "SNP", sort = FALSE)
  if (nrow(comparison) == 0L) stop("No common SNPs across four models for ", phenotype)
  comparison <- comparison[order(match(comparison$SNP, old$SNP)), , drop = FALSE]
  write.table(comparison, file.path(result_dir, "plain_pc_adjusted_clr1to4_clr1to1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  pdf(file.path(result_dir, "qq_plain_pc_adjusted_clr1to4_clr1to1_zoom10.pdf"), width = 22, height = 5.5, useDingbats = FALSE)
  lambdas <- draw_four_zoom(comparison, phenotype); dev.off()
  png(file.path(result_dir, "qq_plain_pc_adjusted_clr1to4_clr1to1_zoom10.png"), width = 5280, height = 1320, res = 220)
  draw_four_zoom(comparison, phenotype); dev.off()

  all_panels[[phenotype]] <- list(comparison = comparison, phenotype = phenotype)
  all_summary[[phenotype]] <- data.frame(
    REGION = region, PHENOTYPE = phenotype, N_COMMON_FOUR_MODELS = nrow(comparison),
    N_PLAIN = comparison$N_UNADJUSTED[1L], N_PC_ADJUSTED = comparison$N_PC_ADJUSTED[1L],
    N_MATCHED_1TO4_BASE = comparison$N_MATCHED_CLR_BASE[1L], N_MATCHED_1TO1_BASE = comparison$N_MATCHED_CLR_1TO1_BASE[1L],
    N_1TO1_OK = sum(clr_1to1$STATUS == "ok"), N_1TO1_FAILED = sum(clr_1to1$STATUS != "ok"),
    LAMBDA_GC_PLAIN = unname(median(comparison$CHISQ_UNADJUSTED) / qchisq(0.5, 1)),
    LAMBDA_GC_PC_ADJUSTED = unname(median(comparison$CHISQ_PC_ADJUSTED) / qchisq(0.5, 1)),
    LAMBDA_GC_MATCHED_1TO4 = unname(median(comparison$CHISQ_MATCHED_CLR) / qchisq(0.5, 1)),
    LAMBDA_GC_MATCHED_1TO1 = unname(lambdas["matched_clr_1to1"]), stringsAsFactors = FALSE
  )
  rm(matched_data, genotype, clr_1to1, old, clr_ok, comparison); gc(verbose = FALSE)
}

all_summary <- do.call(rbind, all_summary)
write.table(all_summary, file.path(output_root, "analysis_summary_1to1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pdf(file.path(output_root, "qq_plain_pc_adjusted_clr1to4_clr1to1_zoom10_all.pdf"), width = 22, height = 5.5, onefile = TRUE, useDingbats = FALSE)
for (panel in all_panels) draw_four_zoom(panel$comparison, panel$phenotype)
dev.off()

writeLines(c(
  "LCT_region TV and qualification: four-model comparison",
  "Full-sample plain logistic: case ~ SNP",
  "Full-sample PC-adjusted logistic: case ~ PC1 + PC2 + PC3 + PC4 + PC5 + SNP",
  "Matched CLR: PC1-PC5-only nearest-neighbor Mahalanobis matching, without replacement",
  "Matching ratios compared: existing 1:4 result and new 1:1 result; both use the same 1,000-case cap and up to 15x candidate-control pool",
  "Matched CLR model: case ~ SNP + strata(STRATUM), exact conditional logistic regression",
  "Test reference: chi-square(1) from the squared SNP Wald statistic; no permutation",
  "Sex: not read, matched, or modeled",
  "Zoomed QQ plots: fixed x/y limits 0-10; values above 10 are clipped for visual comparability"
), file.path(output_root, "METHODS_1to1.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo_1to1.txt"))
message_time("Completed 1:1 comparison: ", output_root)
