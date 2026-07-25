#!/usr/bin/env Rscript

# Full-sample, pre-matching logistic GWAS.
#
# This workflow is deliberately separate from 001_run_clogit.R and 002_plot.R:
#   * no case/control subsampling or matching
#   * no conditional logistic regression
#   * no phenotype permutation
#   * one PC1-PC5 adjusted logistic regression per rs variant
#   * theoretical chi-square(1) QQ plots
#
# Missing genotype dosages are excluded SNP by SNP. Phenotype or covariate
# PC missingness is handled by complete-case exclusion.

options(stringsAsFactors = FALSE)

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
output_root <- file.path(project_root, "results", "00_CLR", "full_sample_logistic_gwas")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
region_phenotypes <- list(
  chr2_random100 = c("lactose_intolerance"),
  LCT_region = c("TV", "qualification", "lactose_intolerance"),
  PCSK9_region = c("TV", "qualification", "CHD")
)
region_files <- c(
  chr2_random100 = "chr2_random100.raw",
  LCT_region = "LCT_region_geno.raw",
  PCSK9_region = "PCSK9_region_geno.raw"
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
  if (length(variant_names) == 0L) {
    stop("No rs variants found in ", path)
  }
  selected <- c("IID", variant_names)
  classes <- c("character", rep("numeric", length(variant_names)))
  read_selected_tsv(path, selected, classes)
}

load_phenotype_ids <- function(phenotype) {
  path <- file.path(data_dir, "processed", paste0(phenotype, ".RDS"))
  if (!file.exists(path)) {
    stop("Processed phenotype file not found: ", path)
  }
  value <- readRDS(path)
  if (!is.list(value) || !all(c("case", "ctrl") %in% names(value))) {
    stop("Phenotype RDS must contain named case and ctrl vectors: ", path)
  }
  list(case = as.character(value$case), ctrl = as.character(value$ctrl))
}

make_binary_phenotype <- function(iids, phenotype) {
  ids <- load_phenotype_ids(phenotype)
  y <- rep.int(NA_integer_, length(iids))

  ctrl_index <- match(ids$ctrl, iids, nomatch = 0L)
  ctrl_index <- ctrl_index[ctrl_index > 0L]
  y[ctrl_index] <- 0L

  # Cases take precedence if malformed input contains an overlap.
  case_index <- match(ids$case, iids, nomatch = 0L)
  case_index <- case_index[case_index > 0L]
  y[case_index] <- 1L
  y
}

logistic_wald_gwas <- function(genotype, covariates, y) {
  variant_names <- colnames(genotype)
  n_variants <- ncol(genotype)
  covariate_design <- cbind("(Intercept)" = 1, as.matrix(covariates))
  family_binomial <- binomial()
  fit_control <- glm.control(epsilon = 1e-8, maxit = 50L, trace = FALSE)

  null_fit <- glm.fit(
    x = covariate_design,
    y = y,
    family = family_binomial,
    control = fit_control
  )
  if (!isTRUE(null_fit$converged) || any(!is.finite(null_fit$coefficients))) {
    stop("Phenotype-specific covariate-only logistic model did not converge")
  }
  start_coefficients <- c(null_fit$coefficients, SNP = 0)

  n <- n_case <- n_control <- n_missing <- integer(n_variants)
  eaf <- maf <- mac <- beta <- se <- z_value <- chisq <- p_value <-
    rep.int(NA_real_, n_variants)
  iterations <- integer(n_variants)
  status <- rep.int("not_tested", n_variants)

  for (j in seq_len(n_variants)) {
    genotype_j <- genotype[, j]
    observed <- is.finite(genotype_j)
    n[j] <- sum(observed)
    n_missing[j] <- length(y) - n[j]
    if (n[j] == 0L) {
      status[j] <- "all_genotypes_missing"
      next
    }

    y_j <- y[observed]
    genotype_j <- genotype_j[observed]
    n_case[j] <- sum(y_j == 1L)
    n_control[j] <- sum(y_j == 0L)
    allele_sum <- sum(genotype_j)
    eaf[j] <- allele_sum / (2 * n[j])
    maf[j] <- min(eaf[j], 1 - eaf[j])
    mac[j] <- min(allele_sum, 2 * n[j] - allele_sum)

    if (n_case[j] == 0L || n_control[j] == 0L) {
      status[j] <- "no_case_control_variation"
      next
    }
    if (length(unique(genotype_j)) < 2L) {
      status[j] <- "monomorphic"
      next
    }

    design_j <- cbind(covariate_design[observed, , drop = FALSE], SNP = genotype_j)
    fit <- tryCatch(
      suppressWarnings(
        glm.fit(
          x = design_j,
          y = y_j,
          family = family_binomial,
          control = fit_control,
          start = start_coefficients
        )
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      status[j] <- "fit_error"
      next
    }
    iterations[j] <- fit$iter
    if (!isTRUE(fit$converged)) {
      status[j] <- "nonconverged"
      next
    }
    if (fit$rank < ncol(design_j) || !is.finite(fit$coefficients[ncol(design_j)])) {
      status[j] <- "singular"
      next
    }

    pivot <- fit$qr$pivot[seq_len(fit$rank)]
    snp_pivot_position <- match(ncol(design_j), pivot)
    if (is.na(snp_pivot_position)) {
      status[j] <- "singular"
      next
    }
    r_matrix <- qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
    covariance_pivot <- tryCatch(
      chol2inv(r_matrix),
      error = function(e) NULL
    )
    if (is.null(covariance_pivot)) {
      status[j] <- "covariance_error"
      next
    }

    beta[j] <- fit$coefficients[ncol(design_j)]
    se[j] <- sqrt(covariance_pivot[snp_pivot_position, snp_pivot_position])
    if (!is.finite(se[j]) || se[j] <= 0) {
      status[j] <- "invalid_se"
      next
    }
    z_value[j] <- beta[j] / se[j]
    chisq[j] <- z_value[j]^2
    p_value[j] <- pchisq(chisq[j], df = 1, lower.tail = FALSE)
    status[j] <- "ok"
  }

  data.frame(
    SNP = variant_names,
    N = n,
    N_CASE = n_case,
    N_CONTROL = n_control,
    N_GENOTYPE_MISSING = n_missing,
    EAF = eaf,
    MAF = maf,
    MAC = mac,
    BETA = beta,
    SE = se,
    Z = z_value,
    CHISQ = chisq,
    P = p_value,
    ITERATIONS = iterations,
    STATUS = status,
    stringsAsFactors = FALSE
  )
}

draw_chisq_qq <- function(result, region, phenotype, sample_label = NULL) {
  observed <- sort(result$CHISQ[is.finite(result$CHISQ)])
  n_variants <- length(observed)
  if (n_variants == 0L) {
    plot.new()
    title(main = paste(region, phenotype, sep = " / "))
    text(0.5, 0.5, "No finite test statistics")
    return(invisible(NULL))
  }

  expected <- qchisq(ppoints(n_variants), df = 1)
  lambda_gc <- median(observed) / qchisq(0.5, df = 1)
  axis_max <- max(c(expected, observed), na.rm = TRUE)

  plot(
    expected,
    observed,
    pch = 19,
    cex = 0.65,
    col = grDevices::adjustcolor("#2166AC", alpha.f = 0.65),
    xlim = c(0, axis_max),
    ylim = c(0, axis_max),
    xlab = "Expected chi-square(1) quantile",
    ylab = "Observed logistic Wald chi-square",
    main = paste(region, phenotype, sep = " / ")
  )
  abline(a = 0, b = 1, lty = 2, lwd = 1.5, col = "#B2182B")
  legend_text <- c(
    sprintf("lambda[GC] = %.3f", lambda_gc),
    sprintf("variants = %d", n_variants)
  )
  if (!is.null(sample_label)) {
    legend_text <- c(legend_text, sample_label)
  }
  legend("topleft", legend = legend_text, bty = "n", cex = 0.85)
  invisible(lambda_gc)
}

message_time("Reading PCs")
pc_names <- paste0("PC", seq_len(num_pcs))
pcs <- read_selected_tsv(
  file.path(data_dir, "UKBB_pca.eigenvec"),
  c("IID", pc_names),
  c("character", rep("numeric", num_pcs))
)
base_sample <- pcs[, c("IID", pc_names), drop = FALSE]
if (anyDuplicated(base_sample$IID)) {
  stop("Duplicate IID values found in the PCA data")
}
base_complete <- complete.cases(base_sample[, pc_names, drop = FALSE])
base_sample <- base_sample[base_complete, , drop = FALSE]
rownames(base_sample) <- NULL
rm(pcs, base_complete)
gc(verbose = FALSE)
message_time("PC sample N = ", nrow(base_sample))

all_summary <- list()
qq_panels <- list()
summary_index <- 0L

for (region in names(region_phenotypes)) {
  message_time("Reading genotypes for ", region)
  genotype_frame <- read_genotypes(file.path(data_dir, unname(region_files[region])))
  base_index <- match(genotype_frame$IID, base_sample$IID)
  aligned <- !is.na(base_index)
  genotype_frame <- genotype_frame[aligned, , drop = FALSE]
  region_sample <- base_sample[base_index[aligned], , drop = FALSE]
  stopifnot(identical(genotype_frame$IID, region_sample$IID))
  rm(base_index, aligned)

  message_time("Preparing ", region, " genotype matrix (N = ", nrow(genotype_frame), ")")
  variant_names <- setdiff(names(genotype_frame), "IID")
  genotype_matrix <- as.matrix(genotype_frame[, variant_names, drop = FALSE])
  storage.mode(genotype_matrix) <- "double"
  rm(genotype_frame)
  gc(verbose = FALSE)

  for (phenotype in region_phenotypes[[region]]) {
    y_all <- make_binary_phenotype(region_sample$IID, phenotype)
    analysis_rows <- !is.na(y_all)
    y <- y_all[analysis_rows]
    covariates <- region_sample[analysis_rows, pc_names, drop = FALSE]
    genotype <- genotype_matrix[analysis_rows, , drop = FALSE]
    n_case <- sum(y == 1L)
    n_control <- sum(y == 0L)

    if (n_case == 0L || n_control == 0L) {
      stop("Phenotype ", phenotype, " has no cases or no controls in ", region)
    }

    message_time(
      "Testing ", region, " / ", phenotype,
      " (N = ", length(y), ", cases = ", n_case,
      ", controls = ", n_control, ", variants = ", ncol(genotype), ")"
    )
    gwas <- logistic_wald_gwas(genotype, covariates, y)

    result <- cbind(data.frame(
      REGION = region,
      PHENOTYPE = phenotype,
      N_SAMPLE = length(y),
      N_CASE_SAMPLE = n_case,
      N_CONTROL_SAMPLE = n_control,
      stringsAsFactors = FALSE
    ), gwas, data.frame(
      TEST = "LOGISTIC_WALD_CHISQ1",
      GENOTYPE_MISSING_METHOD = "COMPLETE_CASE_PER_SNP",
      COVARIATES = paste(pc_names, collapse = "+"),
      stringsAsFactors = FALSE
    ))

    result_dir <- file.path(output_root, region, phenotype)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(result_dir, "logistic_gwas.tsv")
    write.table(
      result,
      result_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )

    sample_label <- sprintf("N = %d (%d cases / %d controls)", length(y), n_case, n_control)
    pdf(file.path(result_dir, "qq_chisq1.pdf"), width = 6.5, height = 6.5, useDingbats = FALSE)
    lambda_gc <- draw_chisq_qq(result, region, phenotype, sample_label)
    dev.off()
    png(file.path(result_dir, "qq_chisq1.png"), width = 1800, height = 1800, res = 240)
    draw_chisq_qq(result, region, phenotype, sample_label)
    dev.off()

    summary_index <- summary_index + 1L
    all_summary[[summary_index]] <- data.frame(
      REGION = region,
      PHENOTYPE = phenotype,
      N = length(y),
      N_CASE = n_case,
      N_CONTROL = n_control,
      N_VARIANTS_TOTAL = nrow(result),
      N_VARIANTS_TESTED = sum(result$STATUS == "ok"),
      N_NONCONVERGED = sum(result$STATUS == "nonconverged"),
      LAMBDA_GC_CHISQ1 = lambda_gc,
      MIN_P = if (any(is.finite(result$P))) min(result$P, na.rm = TRUE) else NA_real_,
      MEDIAN_MODEL_ITERATIONS = if (any(result$ITERATIONS > 0L)) median(result$ITERATIONS[result$ITERATIONS > 0L]) else NA_real_,
      MAX_MODEL_ITERATIONS = if (any(result$ITERATIONS > 0L)) max(result$ITERATIONS) else NA_integer_,
      stringsAsFactors = FALSE
    )
    qq_panels[[summary_index]] <- list(
      result = result,
      region = region,
      phenotype = phenotype,
      sample_label = sample_label
    )

    message_time("Saved ", result_path)
    rm(y_all, analysis_rows, y, covariates, genotype, gwas, result)
    gc(verbose = FALSE)
  }

  rm(genotype_matrix, variant_names, region_sample)
  gc(verbose = FALSE)
}

summary_table <- do.call(rbind, all_summary)
write.table(
  summary_table,
  file.path(output_root, "analysis_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

pdf(
  file.path(output_root, "qq_chisq1_all.pdf"),
  width = 11,
  height = 8.5,
  onefile = TRUE,
  useDingbats = FALSE
)
old_par <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1.5))
for (panel in qq_panels) {
  draw_chisq_qq(panel$result, panel$region, panel$phenotype, panel$sample_label)
}
par(old_par)
dev.off()

writeLines(
  c(
    "Full-sample pre-matching logistic GWAS",
    "Model: binary phenotype ~ PC1 + PC2 + PC3 + PC4 + PC5 + SNP",
    "Test: SNP-coefficient Wald z squared, referenced to chi-square(1)",
    "Genotype missingness: complete-case exclusion within each SNP model",
    "Phenotype/PC missingness: complete-case exclusion",
    "QQ reference: theoretical chi-square distribution with 1 degree of freedom",
    "Matching: not used",
    "Permutation: not used"
  ),
  file.path(output_root, "METHODS.txt")
)

message_time("All analyses complete: ", output_root)
