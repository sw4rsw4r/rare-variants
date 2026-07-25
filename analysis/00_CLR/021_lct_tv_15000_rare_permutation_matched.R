#!/usr/bin/env Rscript

# Permutation-calibrated rare-variant analysis for the validated 15,000-case
# TV matched design. Labels are permuted within each matched stratum, keeping
# one case and the exact 1:1, 1:4, or 1:6 control ratio. This is appropriate for
# the conditional matched null; no sex and no global phenotype permutation.

options(stringsAsFactors = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required")
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
source_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool")
rare_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool_rare")
output_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool_rare_permutation")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

n_perm <- 10000L
perm_seed <- 7201L
chunk_size <- 200L

load_function_definitions <- function(path, names_to_load) {
  parsed <- parse(path)
  for (expr in parsed) {
    if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
      name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (length(name) == 1L && name %in% names_to_load) eval(expr, envir = parent.frame())
    }
  }
}
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"), c("read_header", "read_selected_tsv"))
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "010_pc_matching_helpers.R"), c("draw_original_style"))

# A p-value based diagnostic for permutation-calibrated results. Under a
# continuous null, median(P) = 0.5; values below 0.5 indicate an excess of
# small p-values. This is not the conventional chi-square lambda-GC.
p_lambda <- function(p) {
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  if (length(p) == 0L) return(NA_real_)
  0.5 / median(p)
}
chisq_lambda <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x) / qchisq(.5, 1)
}

draw_with_lambdas <- function(comparison, columns, labels, max_val, lambda_values, lambda_types) {
  rows <- lapply(seq_along(columns), function(i) {
    values <- sort(comparison[[columns[i]]][is.finite(comparison[[columns[i]]])])
    data.frame(exp_chisq = qchisq(ppoints(length(values)), 1), obs_chisq = values, Model = labels[i], stringsAsFactors = FALSE)
  })
  qq <- do.call(rbind, rows)
  qq$Model <- factor(qq$Model, levels = labels)
  anno <- data.frame(Model = factor(labels, levels = labels), label = ifelse(lambda_types == "p", sprintf("lambda_p = %.3f", lambda_values), sprintf("lambda_GC = %.3f", lambda_values)), x = max_val * 0.1, y = max_val * 0.9)
  ggplot(qq, aes(exp_chisq, obs_chisq)) + geom_point(size = 1, alpha = .6) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + facet_wrap(~Model, nrow = 1) + geom_text(data = anno, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) + labs(title = "QQ-plot of rare-variant test statistics", x = "Expected chi-square(1)", y = "Observed chi-square(1)") + theme_minimal() + scale_x_continuous(limits = c(0, max_val)) + scale_y_continuous(limits = c(0, max_val))
}

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
rare_list <- read_tab(file.path(rare_root, "rare_variant_list_maf_lt_0.01.tsv"))
rare_ids <- as.character(rare_list$SNP)
if (length(rare_ids) == 0L) stop("Rare variant list is empty")

geno <- read_selected_tsv(file.path(data_dir, "LCT_region_geno.raw"), c("IID", rare_ids), c("character", rep("numeric", length(rare_ids))))
if (anyDuplicated(geno$IID)) stop("Duplicate IID in genotype data")

# Score statistic under the matched conditional null. In each stratum one
# individual is selected uniformly as the permuted case. U is the sum of
# (selected genotype - stratum mean); V is its exact null variance.
permute_one_variant <- function(snp, ratio, matched, genotype, n_perm, seed, chunk_size) {
  idx <- match(matched$IID, genotype$IID)
  if (anyNA(idx)) stop("Matched IID missing from genotype data")
  x <- as.numeric(genotype[[snp]][idx])
  strata_ids <- unique(matched$STRATUM)
  strata <- lapply(strata_ids, function(s) which(matched$STRATUM == s))
  # Match the existing CLR missingness rule: retain a stratum when its
  # observed case genotype is present and at least one control genotype is
  # present. Missing controls are omitted, so the permutation ratio can be
  # reduced within that stratum while the observed case remains defined.
  valid <- vapply(strata, function(ix) {
    case_ix <- ix[matched$case[ix] == 1L]
    control_ix <- ix[matched$case[ix] == 0L]
    length(case_ix) == 1L && is.finite(x[case_ix]) && sum(is.finite(x[control_ix])) >= 1L
  }, logical(1L))
  strata <- lapply(strata[valid], function(ix) ix[is.finite(x[ix])])
  if (length(strata) == 0L) {
    return(data.frame(SNP = snp, RATIO = ratio, STATUS = "no_valid_strata", N_STRATA_VALID = 0L, N_INFORMATIVE_STRATA = 0L, N_PERM = n_perm, OBS_SCORE_CHISQ = NA_real_, P_SCORE_ASYM = NA_real_, P_PERM = NA_real_, CHISQ_PERM = NA_real_, stringsAsFactors = FALSE))
  }
  xlist <- lapply(strata, function(ix) x[ix])
  means <- vapply(xlist, mean, numeric(1L))
  variances <- vapply(xlist, function(v) mean(v^2) - mean(v)^2, numeric(1L))
  informative <- variances > 0 & is.finite(variances)
  if (!any(informative)) {
    return(data.frame(SNP = snp, RATIO = ratio, STATUS = "no_within_stratum_variation", N_STRATA_VALID = length(xlist), N_INFORMATIVE_STRATA = 0L, N_PERM = n_perm, OBS_SCORE_CHISQ = NA_real_, P_SCORE_ASYM = NA_real_, P_PERM = NA_real_, CHISQ_PERM = NA_real_, stringsAsFactors = FALSE))
  }
  xlist <- xlist[informative]
  means <- means[informative]
  variances <- variances[informative]
  strata_informative <- strata[informative]
  observed_case_values <- vapply(strata_informative, function(ix) {
    case_ix <- ix[matched$case[ix] == 1L]
    x[case_ix]
  }, numeric(1L))
  u_obs <- sum(observed_case_values - means)
  v_obs <- sum(variances)
  obs_score <- (u_obs^2) / v_obs
  count_ge <- 0L
  done <- 0L
  set.seed(seed)
  n_strata <- length(xlist)
  # Group strata by the number of nonmissing genotype-bearing individuals so
  # that each permutation chunk can still be generated with vectorized draws.
  groups <- split(seq_len(n_strata), lengths(xlist))
  while (done < n_perm) {
    b <- min(chunk_size, n_perm - done)
    u_perm <- numeric(b)
    for (group_indices in groups) {
      m <- length(xlist[[group_indices[1L]]])
      xmat_group <- do.call(rbind, xlist[group_indices])
      draw_index <- matrix(sample.int(m, length(group_indices) * b, replace = TRUE), nrow = length(group_indices), ncol = b)
      selected <- xmat_group[cbind(rep(seq_len(length(group_indices)), times = b), as.vector(draw_index))]
      u_perm <- u_perm + colSums(matrix(selected, nrow = length(group_indices), ncol = b)) - sum(means[group_indices])
    }
    chisq_perm <- (u_perm^2) / v_obs
    count_ge <- count_ge + sum(chisq_perm >= obs_score - 1e-12)
    done <- done + b
  }
  p_perm <- (1 + count_ge) / (n_perm + 1)
  data.frame(SNP = snp, RATIO = ratio, STATUS = "ok", N_STRATA_VALID = sum(valid), N_INFORMATIVE_STRATA = n_strata, N_PERM = n_perm, OBS_SCORE_CHISQ = obs_score, P_SCORE_ASYM = pchisq(obs_score, 1, lower.tail = FALSE), P_PERM = p_perm, CHISQ_PERM = qchisq(1 - p_perm, 1), stringsAsFactors = FALSE)
}

perm_tables <- list()
for (ratio in c(1L, 4L, 6L)) {
  matched <- read_tab(file.path(source_root, "matching", paste0("ratio_1to", ratio), "matched_sets.tsv"))
  if (any(grepl("sex", names(matched), ignore.case = TRUE))) stop("Sex entered matched permutation data")
  counts <- table(matched$STRATUM, matched$case)
  if (!all(c("0", "1") %in% colnames(counts)) || any(counts[, "1"] != 1L) || any(counts[, "0"] != ratio) || anyDuplicated(matched$IID)) stop("Invalid matched sets for ratio 1:", ratio)
  results <- lapply(seq_along(rare_ids), function(j) {
    message(sprintf("Permutation ratio 1:%d, variant %d/%d: %s", ratio, j, length(rare_ids), rare_ids[j]))
    permute_one_variant(rare_ids[j], ratio, matched, geno, n_perm, perm_seed + ratio * 1000L + j, chunk_size)
  })
  tab <- do.call(rbind, results)
  perm_tables[[paste0("1to", ratio)]] <- tab
  ratio_dir <- file.path(output_root, paste0("ratio_1to", ratio))
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(tab, file.path(ratio_dir, "rare_permutation_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
permutation_all <- do.call(rbind, perm_tables)
write.table(permutation_all, file.path(output_root, "rare_permutation_results_all_ratios.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Join permutation-calibrated matched results to the already computed full-
# sample asymptotic results for a transparent side-by-side comparison.
plain <- read_tab(file.path(rare_root, "logistic_plain_rare.tsv"))
pc <- read_tab(file.path(rare_root, "logistic_pc_adjusted_rare.tsv"))
plain_ok <- plain[plain$STATUS == "ok", c("SNP", "CHISQ", "P"), drop = FALSE]; names(plain_ok) <- c("SNP", "CHISQ_PLAIN", "P_PLAIN")
pc_ok <- pc[pc$STATUS == "ok", c("SNP", "CHISQ", "P"), drop = FALSE]; names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED")
comparison <- merge(plain_ok, pc_ok, by = "SNP", sort = FALSE)
for (ratio in c(1L, 4L, 6L)) {
  tab <- perm_tables[[paste0("1to", ratio)]]
  ok <- tab[tab$STATUS == "ok", c("SNP", "P_PERM", "CHISQ_PERM", "OBS_SCORE_CHISQ", "N_STRATA_VALID", "N_INFORMATIVE_STRATA"), drop = FALSE]
  names(ok) <- c("SNP", paste0("P_PERM_1TO", ratio), paste0("CHISQ_PERM_1TO", ratio), paste0("OBS_SCORE_CHISQ_1TO", ratio), paste0("N_STRATA_VALID_1TO", ratio), paste0("N_INFORMATIVE_STRATA_1TO", ratio))
  comparison <- merge(comparison, ok, by = "SNP", sort = FALSE)
}
write.table(comparison, file.path(output_root, "comparison_rare_permutation_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_columns <- c("CHISQ_PERM_1TO1", "CHISQ_PERM_1TO4", "CHISQ_PERM_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
plot_labels <- c("Rare matched 1:1 (permutation)", "Rare matched 1:4 (permutation)", "Rare matched 1:6 (permutation)", "Rare full sample + PCs (asymptotic)", "Rare full sample (No PCs; asymptotic)")
if (nrow(comparison) > 0L) {
  full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
  plot_lambdas <- c(p_lambda(comparison$P_PERM_1TO1), p_lambda(comparison$P_PERM_1TO4), p_lambda(comparison$P_PERM_1TO6), chisq_lambda(comparison$CHISQ_PC_ADJUSTED), chisq_lambda(comparison$CHISQ_PLAIN))
  plot_lambda_types <- c("p", "p", "p", "chisq", "chisq")
  plot_full <- draw_with_lambdas(comparison, plot_columns, plot_labels, full_max, plot_lambdas, plot_lambda_types)
  plot_zoom <- draw_with_lambdas(comparison, plot_columns, plot_labels, 10, plot_lambdas, plot_lambda_types)
  ggsave(file.path(output_root, "qq_rare_permutation_15000_1to1_1to4_1to6_full.pdf"), plot_full, width = 19, height = 4.2, units = "in", device = "pdf")
  ggsave(file.path(output_root, "qq_rare_permutation_15000_1to1_1to4_1to6_zoom10.pdf"), plot_zoom, width = 19, height = 4.2, units = "in", device = "pdf")
  ggsave(file.path(output_root, "qq_rare_permutation_15000_1to1_1to4_1to6_zoom10.png"), plot_zoom, width = 19, height = 4.2, units = "in", dpi = 220)
}

lambda <- function(x) if (length(x) == 0L || all(!is.finite(x))) NA_real_ else median(x[is.finite(x)]) / qchisq(.5, 1)
balance_max <- function(ratio) max(abs(read_tab(file.path(source_root, "matching", paste0("ratio_1to", ratio), "pc_balance.tsv"))$POST_SMD))
summary <- data.frame(
  REGION = "LCT_region", PHENOTYPE = "TV", MAF_THRESHOLD = 0.01, N_PERM = n_perm, PERMUTATION_UNIT = "within matched stratum; one case retained",
  N_RARE_RS_VARIANTS = length(rare_ids), N_COMMON_MODELS = nrow(comparison), N_PERM_OK_1TO1 = sum(perm_tables$`1to1`$STATUS == "ok"), N_PERM_OK_1TO4 = sum(perm_tables$`1to4`$STATUS == "ok"), N_PERM_OK_1TO6 = sum(perm_tables$`1to6`$STATUS == "ok"),
  MAX_ABS_SMD_POST_1TO1 = balance_max(1), MAX_ABS_SMD_POST_1TO4 = balance_max(4), MAX_ABS_SMD_POST_1TO6 = balance_max(6),
  LAMBDA_P_PERM_1TO1 = p_lambda(comparison$P_PERM_1TO1), LAMBDA_P_PERM_1TO4 = p_lambda(comparison$P_PERM_1TO4), LAMBDA_P_PERM_1TO6 = p_lambda(comparison$P_PERM_1TO6), LAMBDA_CHISQ_TRANSFORMED_PERM_1TO1 = lambda(comparison$CHISQ_PERM_1TO1), LAMBDA_CHISQ_TRANSFORMED_PERM_1TO4 = lambda(comparison$CHISQ_PERM_1TO4), LAMBDA_CHISQ_TRANSFORMED_PERM_1TO6 = lambda(comparison$CHISQ_PERM_1TO6), LAMBDA_PC_ADJUSTED_ASYMPTOTIC = lambda(comparison$CHISQ_PC_ADJUSTED), LAMBDA_PLAIN_ASYMPTOTIC = lambda(comparison$CHISQ_PLAIN),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(output_root, "analysis_summary_rare_permutation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c(
  "LCT_region / TV rare-variant permutation-calibrated matched analysis",
  "Rare rs variants: MAF < 0.01 using data/plink_LCT.afreq",
  paste0("Permutation count per variant and ratio = ", n_perm, "; seed base = ", perm_seed),
  "For each matched stratum, one individual was selected uniformly as the permuted case while preserving one case and the exact control ratio.",
  "The permutation test uses the matched conditional score statistic; CHISQ_PERM is qchisq(1 - empirical permutation p, 1) for QQ visualization.",
  "For permutation panels, lambda_p = 0.5 / median(P_PERM) is reported; this is a p-value diagnostic, not conventional lambda-GC.",
  "Full-sample plain and PC-adjusted values are the existing asymptotic logistic results and are labelled accordingly; no global phenotype permutation was applied.",
  "A stratum is retained when the observed case genotype and at least one control genotype are present; missing controls are omitted and the reduced stratum is permuted. Variants with no within-stratum variation are not testable."
), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message(sprintf("Completed rare matched permutation analysis: %s; common models=%d", output_root, nrow(comparison)))
