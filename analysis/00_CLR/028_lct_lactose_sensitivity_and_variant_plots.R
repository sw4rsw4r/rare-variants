#!/usr/bin/env Rscript

# LCT/lactose follow-up:
#   * MAF-mask and weighted/unweighted-burden sensitivity set tests
#   * variant-level carrier enrichment and PC-adjusted Firth estimates
#   * publication-oriented dot plot and Firth forest plot

options(stringsAsFactors = FALSE)
setTimeLimit(elapsed = 10 * 60 * 60, transient = FALSE)
if (!requireNamespace("SKAT", quietly = TRUE)) stop("The SKAT package is required")
if (!requireNamespace("logistf", quietly = TRUE)) stop("The logistf package is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required")
suppressPackageStartupMessages(library(SKAT))
suppressPackageStartupMessages(library(logistf))
suppressPackageStartupMessages(library(ggplot2))

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}
data_dir <- file.path(project_root, "data")
output_root <- file.path(project_root, "results", "00_CLR", "LCT_region_lactose_region_set_tests")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
pc_names <- paste0("PC", 1:5)
n_resampling <- 2000000L
base_seed <- 28101L
message_time <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))

load_function_definitions <- function(path, names_to_load) {
  parsed <- parse(path)
  for (expr in parsed) {
    if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
      nm <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (length(nm) == 1L && nm %in% names_to_load) eval(expr, envir = parent.frame())
    }
  }
}
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"), c("read_header", "read_selected_tsv", "load_phenotype_ids", "make_binary_phenotype"))

af <- read.delim(file.path(data_dir, "plink_LCT.afreq"), check.names = FALSE, stringsAsFactors = FALSE)
af$ALT_FREQS <- as.numeric(af$ALT_FREQS)
af$MAF <- pmin(af$ALT_FREQS, 1 - af$ALT_FREQS)
af$SNP <- paste0(af$ID, "_", af$REF)
raw_path <- file.path(data_dir, "LCT_region_geno.raw")
header <- read_header(raw_path)
region_ids <- header[grepl("^(rs|Affx-)", header)]
af_region <- af[af$SNP %in% region_ids & is.finite(af$MAF), , drop = FALSE]
af_region <- af_region[match(region_ids[region_ids %in% af_region$SNP], af_region$SNP), , drop = FALSE]

pc_data <- read_selected_tsv(file.path(data_dir, "UKBB_pca.eigenvec"), c("IID", pc_names), c("character", rep("numeric", 5)))
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
geno <- read_selected_tsv(raw_path, c("IID", region_ids), c("character", rep("numeric", length(region_ids))))
idx <- match(pc_data$IID, geno$IID)
if (anyNA(idx)) stop("PC IID missing from LCT genotype data")
G_all <- as.matrix(geno[idx, region_ids, drop = FALSE])
storage.mode(G_all) <- "double"
rm(geno, idx)
y_all <- make_binary_phenotype(pc_data$IID, "lactose_intolerance")
valid <- !is.na(y_all)
y <- as.numeric(y_all[valid])
G <- G_all[valid, , drop = FALSE]
PC <- pc_data[valid, pc_names, drop = FALSE]
analysis_data <- data.frame(y = y, PC, check.names = FALSE)
message_time("Lactose sensitivity sample N=", length(y), "; cases=", sum(y == 1), "; controls=", sum(y == 0))

extract_skat <- function(out, test, set_label, weight_label, ids) {
  getn <- function(n) { z <- out[[n]]; if (length(z) == 0L) NA_real_ else as.numeric(z[1L]) }
  data.frame(MODEL = "PC_ADJUSTED", TEST = test, WEIGHTS = weight_label, VARIANT_SET = set_label, N_VARIANTS = length(ids), N_NONZERO_CARRIERS = getn("m"), MAC = getn("MAC"), P_VALUE = getn("p.value"), METHOD_BIN = if (length(out[["method.bin"]]) == 0L) NA_character_ else as.character(out[["method.bin"]][1L]), STATUS = "ok", MESSAGE = "", stringsAsFactors = FALSE)
}
run_set <- function(ids, set_label, seed_offset) {
  if (length(ids) == 0L) return(data.frame())
  Z <- G[, ids, drop = FALSE]
  null <- SKAT::SKAT_Null_Model(y ~ PC1 + PC2 + PC3 + PC4 + PC5, data = analysis_data, out_type = "D")
  rows <- list()
  for (j in seq_along(c("BURDEN_WEIGHTED", "BURDEN_UNWEIGHTED", "SKAT", "SKAT_O"))) {
    test <- c("BURDEN_WEIGHTED", "BURDEN_UNWEIGHTED", "SKAT", "SKAT_O")[j]
    method <- c("Burden", "Burden", "SKAT", "SKATO")[j]
    weights <- if (test == "BURDEN_UNWEIGHTED") rep(1, ncol(Z)) else NULL
    weight_label <- if (test == "BURDEN_UNWEIGHTED") "equal" else "beta(1,25)"
    message_time("Sensitivity ", set_label, " / ", test)
    out <- tryCatch(SKAT::SKATBinary(Z, null, method = method, method.bin = "Hybrid", weights.beta = c(1,25), weights = weights, impute.method = "bestguess", missing_cutoff = 0.15, N.Resampling = n_resampling, seednum = base_seed + seed_offset + j), error = function(e) e)
    if (inherits(out, "error")) {
      rows[[j]] <- data.frame(MODEL = "PC_ADJUSTED", TEST = test, WEIGHTS = weight_label, VARIANT_SET = set_label, N_VARIANTS = length(ids), N_NONZERO_CARRIERS = NA_real_, MAC = NA_real_, P_VALUE = NA_real_, METHOD_BIN = NA_character_, STATUS = "error", MESSAGE = conditionMessage(out), stringsAsFactors = FALSE)
    } else rows[[j]] <- extract_skat(out, test, set_label, weight_label, ids)
  }
  do.call(rbind, rows)
}

sets <- list(
  MAF_0.01_ALL = af_region$SNP[af_region$MAF < 0.01],
  MAF_0.01_RS = af_region$SNP[af_region$MAF < 0.01 & grepl("^rs", af_region$SNP)],
  MAF_0.001_ALL = af_region$SNP[af_region$MAF <= 0.001],
  MAF_0.001_RS = af_region$SNP[af_region$MAF <= 0.001 & grepl("^rs", af_region$SNP)]
)
sensitivity <- do.call(rbind, lapply(seq_along(sets), function(i) run_set(sets[[i]], names(sets)[i], i * 100L)))
sensitivity$PHENOTYPE <- "lactose_intolerance"
write.table(sensitivity, file.path(output_root, "lactose_sensitivity_set_tests.tsv"), sep = "	", quote = FALSE, row.names = FALSE, na = "NA")

# Variant-level carrier table and PC-adjusted Firth estimates.
rare_ids <- sets$MAF_0.01_ALL
firth_rows <- vector("list", length(rare_ids))
for (j in seq_along(rare_ids)) {
  snp <- rare_ids[j]
  g <- G[, snp]
  observed <- is.finite(g)
  yj <- y[observed]
  gj <- g[observed]
  pcj <- PC[observed, , drop = FALSE]
  carrier <- gj > 0
  n_case <- sum(yj == 1); n_ctrl <- sum(yj == 0)
  case_carrier <- sum(carrier & yj == 1); ctrl_carrier <- sum(carrier & yj == 0)
  fisher_p <- tryCatch(fisher.test(matrix(c(case_carrier, n_case - case_carrier, ctrl_carrier, n_ctrl - ctrl_carrier), nrow = 2, byrow = TRUE))$p.value, error = function(e) NA_real_)
  dat <- data.frame(y = yj, SNP = gj, pcj, check.names = FALSE)
  message_time("Firth variant ", j, "/", length(rare_ids), ": ", snp)
  fit_warnings <- character(0L)
  fit <- tryCatch(withCallingHandlers(logistf::logistf(y ~ SNP + PC1 + PC2 + PC3 + PC4 + PC5, data = dat, pl = FALSE, model = FALSE, control = logistf::logistf.control(maxit = 50, maxstep = 5)), warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) e)
  if (inherits(fit, "error")) {
    firth_rows[[j]] <- data.frame(SNP = snp, AF_MAF = af_region$MAF[match(snp, af_region$SNP)], N = length(yj), N_MISSING = sum(!observed), MAC = min(sum(gj), 2*length(gj)-sum(gj)), CASE_CARRIER = case_carrier, CONTROL_CARRIER = ctrl_carrier, CASE_CARRIER_RATE = case_carrier/n_case, CONTROL_CARRIER_RATE = ctrl_carrier/n_ctrl, FISHER_P = fisher_p, FIRTH_BETA = NA_real_, FIRTH_OR = NA_real_, FIRTH_CI_LOWER = NA_real_, FIRTH_CI_UPPER = NA_real_, FIRTH_P = NA_real_, STATUS = "error", MESSAGE = paste(c(conditionMessage(fit), fit_warnings), collapse = " | "), stringsAsFactors = FALSE)
  } else {
    b <- unname(fit$coefficients["SNP"])
    firth_rows[[j]] <- data.frame(SNP = snp, AF_MAF = af_region$MAF[match(snp, af_region$SNP)], N = length(yj), N_MISSING = sum(!observed), MAC = min(sum(gj), 2*length(gj)-sum(gj)), CASE_CARRIER = case_carrier, CONTROL_CARRIER = ctrl_carrier, CASE_CARRIER_RATE = case_carrier/n_case, CONTROL_CARRIER_RATE = ctrl_carrier/n_ctrl, FISHER_P = fisher_p, FIRTH_BETA = b, FIRTH_OR = exp(b), FIRTH_CI_LOWER = exp(unname(fit$ci.lower["SNP"])), FIRTH_CI_UPPER = exp(unname(fit$ci.upper["SNP"])), FIRTH_P = unname(fit$prob["SNP"]), STATUS = "ok", MESSAGE = paste(fit_warnings, collapse = " | "), stringsAsFactors = FALSE)
  }
}
firth_table <- do.call(rbind, firth_rows)
write.table(firth_table, file.path(output_root, "lactose_variant_carrier_firth_pc_adjusted.tsv"), sep = "	", quote = FALSE, row.names = FALSE, na = "NA")

# Dot plot of the current set-level tests; one regional set does not justify a QQ plot.
primary <- read.delim(file.path(output_root, "lactose_region_burden_skat_skato_results.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
primary <- primary[primary$STATUS == "ok", , drop = FALSE]
primary$LABEL <- paste(primary$MODEL, primary$VARIANT_SET, primary$METHOD, sep = " / ")
primary$NEG_LOG10_P <- -log10(primary$P_VALUE)
dot <- ggplot(primary, aes(x = METHOD, y = NEG_LOG10_P, color = VARIANT_SET, shape = MODEL)) + geom_point(position = position_dodge(width = 0.4), size = 3) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") + labs(title = "LCT-region lactose rare-variant set tests", x = NULL, y = expression(-log[10](p)), caption = "PC-adjusted models are primary; plain models are reference") + theme_minimal(base_size = 12)
ggsave(file.path(output_root, "lactose_set_test_pvalue_dotplot.pdf"), dot, width = 9, height = 5.5, units = "in", device = "pdf", bg = "white")
ggsave(file.path(output_root, "lactose_set_test_pvalue_dotplot.png"), dot, width = 9, height = 5.5, units = "in", dpi = 220, bg = "white")

forest_data <- firth_table[firth_table$STATUS == "ok" & is.finite(firth_table$FIRTH_OR), , drop = FALSE]
forest_data$SNP <- factor(forest_data$SNP, levels = rev(forest_data$SNP))
forest <- ggplot(forest_data, aes(x = FIRTH_OR, y = SNP)) + geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") + geom_errorbarh(aes(xmin = FIRTH_CI_LOWER, xmax = FIRTH_CI_UPPER), height = 0.2, na.rm = TRUE) + geom_point(size = 2.2, na.rm = TRUE) + scale_x_log10() + labs(title = "LCT-region rare variants: PC-adjusted Firth estimates", x = "Odds ratio (log scale)", y = NULL) + theme_minimal(base_size = 11)
ggsave(file.path(output_root, "lactose_variant_firth_forest.pdf"), forest, width = 8.5, height = 7.5, units = "in", device = "pdf", bg = "white")
ggsave(file.path(output_root, "lactose_variant_firth_forest.png"), forest, width = 8.5, height = 7.5, units = "in", dpi = 220, bg = "white")

writeLines(c("LCT/lactose sensitivity and variant-level follow-up", "Set-test sensitivity: MAF < 0.01 and <= 0.001, all region variants versus rs-only, weighted beta(1,25) burden versus equal-weight burden.", "Variant-level table: carrier counts, Fisher exact p-value, and PC1-PC5-adjusted Firth logistic estimate.", "The dot plot summarizes set-level p-values; a QQ plot is not used for one region because there are too few independent set p-values.", "The forest plot uses Firth ORs and confidence intervals; rs4988235 is not in the rare set because it is common."), file.path(output_root, "METHODS_sensitivity_and_variant_plots.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo_sensitivity_and_variant_plots.txt"))
message_time("Completed LCT lactose sensitivity, Firth, and figures: ", output_root)
