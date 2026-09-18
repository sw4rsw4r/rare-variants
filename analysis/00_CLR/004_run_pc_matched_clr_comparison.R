#!/usr/bin/env Rscript

# PC-only matching followed by conditional logistic regression.
#
# The manuscript matching analyses are run from 016, 018, and 019. This file
# is kept because those scripts load shared helper functions from it.
#
#   * all joins and genotype alignments are performed explicitly by IID;
#   * matching uses PC1-PC5 only (no sex or genotype information);
#   * 1:4 Mahalanobis nearest-neighbor matching is without replacement;
#   * matched sets are created once per phenotype and reused across regions;
#   * SNP-specific missingness may reduce a stratum from 1:4 to 1:k; a stratum
#     is removed only when its case is missing or no nonmissing control remains;
#   * the CLR model is case ~ SNP + strata only;
#   * no permutation is performed;
#   * QQ comparisons use the same successfully tested SNPs in both models.

options(stringsAsFactors = FALSE)

if (!requireNamespace("MatchIt", quietly = TRUE)) {
  stop("The MatchIt package is required")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival package is required")
}
suppressPackageStartupMessages(library(survival))

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(
    file.path(dirname(script_file), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

data_dir <- file.path(project_root, "data")
pc_adjusted_root <- file.path(project_root, "results", "00_CLR", "full_sample_logistic_gwas")
output_root <- file.path(project_root, "results", "00_CLR", "pc_adjusted_vs_matched_clr")
matching_root <- file.path(output_root, "matching")
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
num_controls <- 4L
max_cases <- 1000L
control_pool_multiplier <- 15L
pc_names <- paste0("PC", seq_len(num_pcs))

region_phenotypes <- list(
  chr2_random100 = c("TV"),
  LCT_region = c("TV", "lactose_intolerance")
)
region_files <- c(
  chr2_random100 = "chr2_random100.raw",
  LCT_region = "LCT_region_geno.raw"
)
phenotype_seeds <- c(
  TV = 1101L,
  lactose_intolerance = 1103L
)

message_time <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

read_header <- function(path) {
  strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
}

read_selected_tsv <- function(path, selected, classes) {
  header <- read_header(path)
  missing_columns <- setdiff(selected, header)
  if (length(missing_columns) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing_columns, collapse = ", "))
  }
  col_classes <- rep("NULL", length(header))
  names(col_classes) <- header
  col_classes[selected] <- classes
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    colClasses = unname(col_classes),
    check.names = FALSE,
    comment.char = "",
    quote = "",
    na.strings = "NA"
  )
}

read_genotypes <- function(path) {
  header <- read_header(path)
  variant_names <- header[grepl("^rs", header)]
  selected <- c("IID", variant_names)
  classes <- c("character", rep("numeric", length(variant_names)))
  read_selected_tsv(path, selected, classes)
}

load_phenotype_ids <- function(phenotype) {
  path <- file.path(data_dir, "processed", paste0(phenotype, ".RDS"))
  value <- readRDS(path)
  if (!is.list(value) || !all(c("case", "ctrl") %in% names(value))) {
    stop("Phenotype RDS must contain named case and ctrl vectors: ", path)
  }
  list(case = as.character(value$case), ctrl = as.character(value$ctrl))
}

make_binary_phenotype <- function(iids, phenotype) {
  ids <- load_phenotype_ids(phenotype)
  y <- rep.int(NA_integer_, length(iids))
  control_index <- match(ids$ctrl, iids, nomatch = 0L)
  control_index <- control_index[control_index > 0L]
  y[control_index] <- 0L
  case_index <- match(ids$case, iids, nomatch = 0L)
  case_index <- case_index[case_index > 0L]
  y[case_index] <- 1L
  y
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

add_mahalanobis_distances <- function(matched_data, candidate_data) {
  covariance <- stats::cov(candidate_data[, pc_names, drop = FALSE])
  inverse_covariance <- tryCatch(
    solve(covariance),
    error = function(e) solve(covariance + diag(1e-10, nrow(covariance)))
  )
  matched_data$MATCH_DISTANCE <- 0
  for (stratum in unique(matched_data$STRATUM)) {
    rows <- which(matched_data$STRATUM == stratum)
    case_row <- rows[matched_data$case[rows] == 1L]
    control_rows <- rows[matched_data$case[rows] == 0L]
    case_pc <- as.numeric(matched_data[case_row, pc_names, drop = TRUE])
    for (control_row in control_rows) {
      difference <- as.numeric(matched_data[control_row, pc_names, drop = TRUE]) - case_pc
      matched_data$MATCH_DISTANCE[control_row] <- sqrt(
        as.numeric(crossprod(difference, inverse_covariance %*% difference))
      )
    }
  }
  matched_data
}

match_phenotype <- function(pc_data, phenotype) {
  y <- make_binary_phenotype(pc_data$IID, phenotype)
  case_available <- which(y == 1L)
  control_available <- which(y == 0L)
  if (length(case_available) == 0L || length(control_available) < num_controls) {
    stop("Insufficient cases or controls for ", phenotype)
  }

  set.seed(unname(phenotype_seeds[phenotype]))
  selected_cases <- if (length(case_available) > max_cases) {
    sample(case_available, max_cases, replace = FALSE)
  } else {
    case_available
  }
  candidate_control_n <- min(
    length(control_available),
    max(num_controls * length(selected_cases), control_pool_multiplier * length(selected_cases))
  )
  selected_controls <- if (length(control_available) > candidate_control_n) {
    sample(control_available, candidate_control_n, replace = FALSE)
  } else {
    control_available
  }

  candidate_rows <- c(selected_cases, selected_controls)
  candidate_data <- pc_data[candidate_rows, c("IID", pc_names), drop = FALSE]
  candidate_data$case <- c(rep.int(1L, length(selected_cases)), rep.int(0L, length(selected_controls)))
  candidate_data <- candidate_data[, c("IID", "case", pc_names), drop = FALSE]
  rownames(candidate_data) <- candidate_data$IID

  formula_match <- stats::reformulate(pc_names, response = "case")
  match_object <- MatchIt::matchit(
    formula_match,
    data = candidate_data,
    method = "nearest",
    distance = "mahalanobis",
    estimand = "ATT",
    replace = FALSE,
    m.order = "farthest",
    ratio = num_controls
  )
  matched_data <- MatchIt::match.data(match_object, data = candidate_data, drop.unmatched = TRUE)
  matched_data$STRATUM <- match(
    as.character(matched_data$subclass),
    unique(as.character(matched_data$subclass))
  )
  matched_data <- matched_data[, c("IID", "case", pc_names, "STRATUM", "weights"), drop = FALSE]
  matched_data <- matched_data[order(matched_data$STRATUM, -matched_data$case), , drop = FALSE]
  rownames(matched_data) <- NULL

  set_counts <- table(matched_data$STRATUM, matched_data$case)
  if (!all(c("0", "1") %in% colnames(set_counts)) ||
      any(set_counts[, "1"] != 1L) || any(set_counts[, "0"] != num_controls)) {
    stop("Matching did not produce complete 1:", num_controls, " sets for ", phenotype)
  }
  if (anyDuplicated(matched_data$IID)) {
    stop("A participant was reused in matching for ", phenotype)
  }

  balance <- compute_balance(candidate_data, matched_data)
  matched_data <- add_mahalanobis_distances(matched_data, candidate_data)
  control_distances <- matched_data$MATCH_DISTANCE[matched_data$case == 0L]

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
    MEDIAN_MAHALANOBIS_DISTANCE = stats::median(control_distances),
    MAX_MAHALANOBIS_DISTANCE = max(control_distances),
    stringsAsFactors = FALSE
  )

  phenotype_dir <- file.path(matching_root, phenotype)
  dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(matched_data, file.path(phenotype_dir, "matched_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(balance, file.path(phenotype_dir, "pc_balance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary, file.path(phenotype_dir, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  list(data = matched_data, balance = balance, summary = summary)
}

run_clogit_gwas <- function(genotype, matched_data) {
  variant_names <- colnames(genotype)
  n_variants <- ncol(genotype)
  result <- data.frame(
    SNP = variant_names,
    N_MATCHED_BASE = nrow(matched_data),
    N_STRATA_BASE = length(unique(matched_data$STRATUM)),
    N = integer(n_variants),
    N_CASE = integer(n_variants),
    N_CONTROL = integer(n_variants),
    N_STRATA = integer(n_variants),
    N_STRATA_COMPLETE_1TO4 = integer(n_variants),
    N_STRATA_REDUCED = integer(n_variants),
    N_INFORMATIVE_STRATA = integer(n_variants),
    N_GENOTYPE_MISSING = integer(n_variants),
    EAF = rep.int(NA_real_, n_variants),
    MAF = rep.int(NA_real_, n_variants),
    MAC = rep.int(NA_real_, n_variants),
    BETA = rep.int(NA_real_, n_variants),
    SE = rep.int(NA_real_, n_variants),
    Z = rep.int(NA_real_, n_variants),
    CHISQ = rep.int(NA_real_, n_variants),
    P = rep.int(NA_real_, n_variants),
    ITERATIONS = integer(n_variants),
    STATUS = rep.int("not_tested", n_variants),
    WARNING = rep.int("", n_variants),
    TEST = rep.int("MATCHED_CLR_WALD_CHISQ1", n_variants),
    COVARIATES = rep.int("NONE_PC_MATCHED", n_variants),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(n_variants)) {
    genotype_j <- genotype[, j]
    nonmissing <- is.finite(genotype_j)
    result$N_GENOTYPE_MISSING[j] <- sum(!nonmissing)

    case_present <- tapply(nonmissing & matched_data$case == 1L, matched_data$STRATUM, any)
    control_nonmissing <- tapply(nonmissing & matched_data$case == 0L, matched_data$STRATUM, sum)
    valid_strata <- as.integer(names(case_present)[case_present & control_nonmissing >= 1L])
    use <- nonmissing & matched_data$STRATUM %in% valid_strata
    if (!any(use)) {
      result$STATUS[j] <- "no_valid_strata"
      next
    }

    model_data <- data.frame(
      case = matched_data$case[use],
      SNP = genotype_j[use],
      STRATUM = matched_data$STRATUM[use]
    )
    controls_per_stratum <- table(model_data$STRATUM[model_data$case == 0L])
    informative <- tapply(model_data$SNP, model_data$STRATUM, function(value) length(unique(value)) > 1L)

    result$N[j] <- nrow(model_data)
    result$N_CASE[j] <- sum(model_data$case == 1L)
    result$N_CONTROL[j] <- sum(model_data$case == 0L)
    result$N_STRATA[j] <- length(unique(model_data$STRATUM))
    result$N_STRATA_COMPLETE_1TO4[j] <- sum(controls_per_stratum == num_controls)
    result$N_STRATA_REDUCED[j] <- sum(controls_per_stratum < num_controls)
    result$N_INFORMATIVE_STRATA[j] <- sum(informative)
    allele_sum <- sum(model_data$SNP)
    result$EAF[j] <- allele_sum / (2 * nrow(model_data))
    result$MAF[j] <- min(result$EAF[j], 1 - result$EAF[j])
    result$MAC[j] <- min(allele_sum, 2 * nrow(model_data) - allele_sum)

    if (result$N_INFORMATIVE_STRATA[j] == 0L) {
      result$STATUS[j] <- "no_within_stratum_variation"
      next
    }

    warning_text <- character()
    fit <- tryCatch(
      withCallingHandlers(
        survival::clogit(
          case ~ SNP + strata(STRATUM),
          data = model_data,
          method = "exact",
          control = survival::coxph.control(iter.max = 50L, eps = 1e-9)
        ),
        warning = function(w) {
          warning_text <<- c(warning_text, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      result$STATUS[j] <- "fit_error"
      result$WARNING[j] <- conditionMessage(fit)
      next
    }
    result$ITERATIONS[j] <- if (length(fit$iter) > 0L) max(fit$iter) else NA_integer_
    result$WARNING[j] <- paste(unique(warning_text), collapse = " | ")
    if (any(grepl("infinite|converged before", warning_text, ignore.case = TRUE))) {
      result$STATUS[j] <- "possible_separation"
      next
    }

    coefficient_table <- summary(fit)$coefficients
    if (!("SNP" %in% rownames(coefficient_table))) {
      result$STATUS[j] <- "snp_coefficient_missing"
      next
    }
    beta <- coefficient_table["SNP", "coef"]
    se <- coefficient_table["SNP", "se(coef)"]
    if (!is.finite(beta) || !is.finite(se) || se <= 0) {
      result$STATUS[j] <- "invalid_coefficient"
      next
    }
    z_value <- beta / se
    result$BETA[j] <- beta
    result$SE[j] <- se
    result$Z[j] <- z_value
    result$CHISQ[j] <- z_value^2
    result$P[j] <- pchisq(result$CHISQ[j], df = 1, lower.tail = FALSE)
    result$STATUS[j] <- "ok"
  }
  result
}

draw_qq_panel <- function(chisq, model_name, subtitle, axis_max) {
  observed <- sort(chisq[is.finite(chisq)])
  expected <- qchisq(ppoints(length(observed)), df = 1)
  lambda_gc <- stats::median(observed) / qchisq(0.5, df = 1)
  plot(
    expected,
    observed,
    pch = 19,
    cex = 0.65,
    col = grDevices::adjustcolor("#2166AC", alpha.f = 0.65),
    xlim = c(0, axis_max),
    ylim = c(0, axis_max),
    xlab = "Expected chi-square(1) quantile",
    ylab = "Observed Wald chi-square",
    main = model_name
  )
  abline(a = 0, b = 1, lty = 2, lwd = 1.5, col = "#B2182B")
  legend(
    "topleft",
    legend = c(sprintf("lambda[GC] = %.3f", lambda_gc), sprintf("variants = %d", length(observed)), subtitle),
    bty = "n",
    cex = 0.82
  )
  invisible(lambda_gc)
}

draw_comparison <- function(comparison, region, phenotype, pc_subtitle, clr_subtitle) {
  expected <- qchisq(ppoints(nrow(comparison)), df = 1)
  axis_max <- max(
    expected,
    comparison$CHISQ_PC_ADJUSTED,
    comparison$CHISQ_MATCHED_CLR,
    na.rm = TRUE
  )
  old_par <- par(mfrow = c(1, 2), mar = c(4.7, 4.7, 4.3, 1.3), oma = c(0, 0, 2.2, 0))
  lambda_pc <- draw_qq_panel(comparison$CHISQ_PC_ADJUSTED, "Full sample: PC-adjusted logistic", pc_subtitle, axis_max)
  lambda_clr <- draw_qq_panel(comparison$CHISQ_MATCHED_CLR, "PC-matched 1:4 CLR", clr_subtitle, axis_max)
  mtext(paste(region, phenotype, sep = " / "), outer = TRUE, side = 3, line = 0.3, font = 2, cex = 1.15)
  par(old_par)
  c(pc_adjusted = lambda_pc, matched_clr = lambda_clr)
}

message_time("Reading PCs")
pc_data <- read_selected_tsv(
  file.path(data_dir, "UKBB_pca.eigenvec"),
  c("IID", pc_names),
  c("character", rep("numeric", num_pcs))
)
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
if (anyDuplicated(pc_data$IID)) stop("Duplicate IID values found in PCA data")

phenotypes <- unique(unlist(region_phenotypes, use.names = FALSE))
match_results <- setNames(vector("list", length(phenotypes)), phenotypes)
for (phenotype in phenotypes) {
  message_time("Matching ", phenotype, " on PC1-PC5 only")
  match_results[[phenotype]] <- match_phenotype(pc_data, phenotype)
  message_time(
    "Matched ", match_results[[phenotype]]$summary$N_CASE_MATCHED,
    " cases; max post-match |SMD| = ",
    sprintf("%.4f", match_results[[phenotype]]$summary$MAX_ABS_SMD_POST)
  )
}

matching_summary <- do.call(rbind, lapply(match_results, `[[`, "summary"))
matching_balance <- do.call(rbind, lapply(names(match_results), function(phenotype) {
  balance <- match_results[[phenotype]]$balance
  balance$PHENOTYPE <- phenotype
  balance[, c("PHENOTYPE", setdiff(names(balance), "PHENOTYPE")), drop = FALSE]
}))
write.table(matching_summary, file.path(output_root, "matching_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(matching_balance, file.path(output_root, "matching_balance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

comparison_panels <- list()
analysis_summary <- list()
panel_index <- 0L

for (region in names(region_phenotypes)) {
  message_time("Reading genotypes for ", region)
  genotype_frame <- read_genotypes(file.path(data_dir, unname(region_files[region])))
  variant_names <- setdiff(names(genotype_frame), "IID")

  for (phenotype in region_phenotypes[[region]]) {
    matched_data <- match_results[[phenotype]]$data
    genotype_index <- match(matched_data$IID, genotype_frame$IID)
    if (anyNA(genotype_index)) {
      stop("Matched IID missing from ", region, " genotype data")
    }
    genotype <- as.matrix(genotype_frame[genotype_index, variant_names, drop = FALSE])
    storage.mode(genotype) <- "double"
    stopifnot(identical(matched_data$IID, genotype_frame$IID[genotype_index]))

    message_time("Running matched CLR for ", region, " / ", phenotype)
    clr_result <- run_clogit_gwas(genotype, matched_data)
    result_dir <- file.path(output_root, region, phenotype)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    write.table(clr_result, file.path(result_dir, "matched_clr_gwas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    pc_path <- file.path(pc_adjusted_root, region, phenotype, "logistic_gwas.tsv")
    if (!file.exists(pc_path)) {
      stop("PC-adjusted result is missing: ", pc_path)
    }
    pc_result <- read.delim(pc_path, check.names = FALSE)
    if (!all(pc_result$COVARIATES == paste(pc_names, collapse = "+"))) {
      stop("PC-adjusted results contain unexpected covariates: ", pc_path)
    }

    pc_ok <- pc_result[pc_result$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]
    names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
    clr_ok <- clr_result[clr_result$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_STRATA", "N_STRATA_REDUCED"), drop = FALSE]
    names(clr_ok) <- c("SNP", "CHISQ_MATCHED_CLR", "P_MATCHED_CLR", "N_MATCHED_CLR", "N_STRATA_MATCHED_CLR", "N_REDUCED_STRATA_MATCHED_CLR")
    comparison <- merge(pc_ok, clr_ok, by = "SNP", sort = FALSE)
    if (nrow(comparison) == 0L) stop("No common successfully tested SNPs for ", region, " / ", phenotype)
    write.table(comparison, file.path(result_dir, "pc_adjusted_vs_matched_clr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    pc_subtitle <- sprintf("N = %d", unique(pc_result$N_SAMPLE)[1L])
    clr_subtitle <- sprintf(
      "baseline N = %d; median SNP N = %d",
      nrow(matched_data),
      as.integer(stats::median(comparison$N_MATCHED_CLR))
    )

    pdf(file.path(result_dir, "qq_comparison_chisq1.pdf"), width = 11, height = 5.5, useDingbats = FALSE)
    lambdas <- draw_comparison(comparison, region, phenotype, pc_subtitle, clr_subtitle)
    dev.off()
    png(file.path(result_dir, "qq_comparison_chisq1.png"), width = 2640, height = 1320, res = 220)
    draw_comparison(comparison, region, phenotype, pc_subtitle, clr_subtitle)
    dev.off()

    panel_index <- panel_index + 1L
    comparison_panels[[panel_index]] <- list(
      comparison = comparison,
      region = region,
      phenotype = phenotype,
      pc_subtitle = pc_subtitle,
      clr_subtitle = clr_subtitle
    )
    analysis_summary[[panel_index]] <- data.frame(
      REGION = region,
      PHENOTYPE = phenotype,
      N_COMMON_SNPS = nrow(comparison),
      N_PC_ADJUSTED = unique(pc_result$N_SAMPLE)[1L],
      N_MATCHED_BASE = nrow(matched_data),
      MEDIAN_N_MATCHED_CLR = stats::median(comparison$N_MATCHED_CLR),
      MIN_N_MATCHED_CLR = min(comparison$N_MATCHED_CLR),
      MEDIAN_STRATA_MATCHED_CLR = stats::median(comparison$N_STRATA_MATCHED_CLR),
      MIN_STRATA_MATCHED_CLR = min(comparison$N_STRATA_MATCHED_CLR),
      N_CLR_FAILED = sum(clr_result$STATUS != "ok"),
      LAMBDA_GC_PC_ADJUSTED = unname(lambdas["pc_adjusted"]),
      LAMBDA_GC_MATCHED_CLR = unname(lambdas["matched_clr"]),
      MIN_P_PC_ADJUSTED = min(comparison$P_PC_ADJUSTED),
      MIN_P_MATCHED_CLR = min(comparison$P_MATCHED_CLR),
      stringsAsFactors = FALSE
    )
    message_time("Saved comparison for ", region, " / ", phenotype)
    rm(genotype, clr_result, pc_result, pc_ok, clr_ok, comparison)
    gc(verbose = FALSE)
  }
  rm(genotype_frame)
  gc(verbose = FALSE)
}

analysis_summary <- do.call(rbind, analysis_summary)
write.table(analysis_summary, file.path(output_root, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(output_root, "qq_pc_adjusted_vs_matched_clr_all.pdf"), width = 11, height = 5.5, onefile = TRUE, useDingbats = FALSE)
for (panel in comparison_panels) {
  draw_comparison(panel$comparison, panel$region, panel$phenotype, panel$pc_subtitle, panel$clr_subtitle)
}
dev.off()

writeLines(
  c(
    "PC-adjusted logistic GWAS versus PC-matched CLR",
    "",
    "Full-sample model: case ~ PC1 + PC2 + PC3 + PC4 + PC5 + SNP",
    "Matching: PC1-PC5 only; 1:4 nearest-neighbor Mahalanobis distance; no replacement",
    sprintf("Sampling: all cases when <= %d; otherwise a seeded random sample of %d cases", max_cases, max_cases),
    sprintf("Candidate controls: up to %dx the selected case count", control_pool_multiplier),
    "Matched CLR model: case ~ SNP + strata(STRATUM)",
    "SNP missingness: drop a stratum only if its case is missing or no nonmissing control remains",
    "Test: squared SNP Wald z statistic referenced to chi-square(1)",
    "QQ comparison: theoretical chi-square(1), restricted to SNPs successfully tested in both models",
    "Sex: not read, matched, or modeled",
    "Permutation: not used",
    sprintf("The %d-case cap is reported explicitly; matched estimates target the selected case sample.", max_cases)
  ),
  file.path(output_root, "AUDIT_AND_METHODS.txt")
)

writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message_time("All matched CLR comparisons complete: ", output_root)
