#!/usr/bin/env Rscript

# LCT-region rare-variant set tests for lactose_intolerance.
# Primary model: full-sample PC1-PC5-adjusted binary-trait tests. Sex is not
# included. Plain/no-PC results are retained as a reference.

options(stringsAsFactors = FALSE)
setTimeLimit(elapsed = 10 * 60 * 60, transient = FALSE)
if (!requireNamespace("SKAT", quietly = TRUE)) stop("The SKAT package is required")
suppressPackageStartupMessages(library(SKAT))

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
num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
n_resampling <- 2000000L
base_seed <- 27001L

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
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"), c("read_header", "read_selected_tsv", "load_phenotype_ids", "make_binary_phenotype"))
required <- c("read_header", "read_selected_tsv", "load_phenotype_ids", "make_binary_phenotype")
if (!all(vapply(required, exists, logical(1L), envir = environment(), inherits = FALSE))) stop("Could not load phenotype/genotype helpers")

af <- read.delim(file.path(data_dir, "plink_LCT.afreq"), check.names = FALSE, stringsAsFactors = FALSE)
af$ALT_FREQS <- as.numeric(af$ALT_FREQS)
af$MAF <- pmin(af$ALT_FREQS, 1 - af$ALT_FREQS)
af$SNP <- paste0(af$ID, "_", af$REF)
raw_path <- file.path(data_dir, "LCT_region_geno.raw")
geno_header <- read_header(raw_path)
region_variants <- geno_header[grepl("^(rs|Affx-)", geno_header)]
rare_all <- af[is.finite(af$MAF) & af$MAF < 0.01 & af$SNP %in% region_variants, c("SNP", "ID", "REF", "ALT", "MAF"), drop = FALSE]
rare_all <- rare_all[match(region_variants[region_variants %in% rare_all$SNP], rare_all$SNP), , drop = FALSE]
rare_rs <- rare_all[grepl("^rs", rare_all$SNP), , drop = FALSE]
if (nrow(rare_all) == 0L) stop("No rare variants found in LCT_region")
write.table(rare_all, file.path(output_root, "rare_variant_list_all_maf_lt_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rare_rs, file.path(output_root, "rare_variant_list_rs_only_maf_lt_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message_time("Rare variants: all=", nrow(rare_all), "; rs-only=", nrow(rare_rs), "; N.Resampling=", n_resampling)

pc_data <- read_selected_tsv(file.path(data_dir, "UKBB_pca.eigenvec"), c("IID", pc_names), c("character", rep("numeric", num_pcs)))
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
if (anyDuplicated(pc_data$IID)) stop("Duplicate IID in PCA data")
all_ids <- rare_all$SNP
geno <- read_selected_tsv(raw_path, c("IID", all_ids), c("character", rep("numeric", length(all_ids))))
if (anyDuplicated(geno$IID)) stop("Duplicate IID in genotype data")
idx <- match(pc_data$IID, geno$IID)
if (anyNA(idx)) stop("Some PC IIDs are missing from LCT genotype data")
genotype_all <- as.matrix(geno[idx, all_ids, drop = FALSE])
storage.mode(genotype_all) <- "double"
stopifnot(identical(pc_data$IID, geno$IID[idx]))
rm(geno, idx)

extract_result <- function(out, method, model, set_label, n, n_case, n_control, n_variants, set_ids) {
  get_num <- function(name) { value <- out[[name]]; if (length(value) == 0L) NA_real_ else as.numeric(value[1L]) }
  get_chr <- function(name) { value <- out[[name]]; if (length(value) == 0L) NA_character_ else as.character(value[1L]) }
  get_lgl <- function(name) { value <- out[[name]]; if (length(value) == 0L) NA else as.logical(value[1L]) }
  data.frame(MODEL = model, METHOD = method, VARIANT_SET = set_label, N = n, N_CASE = n_case, N_CONTROL = n_control, N_VARIANTS = n_variants, N_NONZERO_CARRIERS = get_num("m"), VARIANT_IDS = paste(set_ids, collapse = ";"), P_VALUE = get_num("p.value"), P_VALUE_STANDARD = get_num("p.value.standard"), P_VALUE_RESAMPLING = get_num("p.value.resampling"), MAC = get_num("MAC"), METHOD_BIN = get_chr("method.bin"), IS_EXACT = get_lgl("Is.ExactP"), IS_ACCURATE = get_lgl("is.accurate"), STATUS = "ok", MESSAGE = "", stringsAsFactors = FALSE)
}

run_one_set <- function(Z, set_label, set_ids, y, covariates, model_label, seed_offset) {
  analysis_data <- data.frame(y = as.numeric(y), covariates, check.names = FALSE)
  formula_null <- if (model_label == "PC_ADJUSTED") stats::reformulate(pc_names, response = "y") else y ~ 1
  null_model <- SKAT::SKAT_Null_Model(formula_null, data = analysis_data, out_type = "D")
  method_values <- c(BURDEN = "Burden", SKAT = "SKAT", SKAT_O = "SKATO")
  out_rows <- vector("list", length(method_values))
  for (j in seq_along(method_values)) {
    method_label <- names(method_values)[j]
    message_time("Running ", model_label, " / ", method_label, " / ", set_label, " (N=", nrow(Z), ")")
    out <- tryCatch(SKAT::SKATBinary(Z = Z, obj = null_model, method = method_values[[j]], method.bin = "Hybrid", weights.beta = c(1, 25), impute.method = "bestguess", is_check_genotype = TRUE, is_dosage = FALSE, missing_cutoff = 0.15, max_maf = 1, estimate_MAF = 1, N.Resampling = n_resampling, seednum = base_seed + seed_offset + j), error = function(e) e)
    if (inherits(out, "error")) {
      out_rows[[j]] <- data.frame(MODEL = model_label, METHOD = method_label, VARIANT_SET = set_label, N = nrow(Z), N_CASE = sum(y == 1L), N_CONTROL = sum(y == 0L), N_VARIANTS = ncol(Z), N_NONZERO_CARRIERS = NA_real_, VARIANT_IDS = paste(set_ids, collapse = ";"), P_VALUE = NA_real_, P_VALUE_STANDARD = NA_real_, P_VALUE_RESAMPLING = NA_real_, MAC = NA_real_, METHOD_BIN = NA_character_, IS_EXACT = NA, IS_ACCURATE = NA, STATUS = "error", MESSAGE = conditionMessage(out), stringsAsFactors = FALSE)
    } else {
      out_rows[[j]] <- extract_result(out, method_label, model_label, set_label, nrow(Z), sum(y == 1L), sum(y == 0L), ncol(Z), set_ids)
    }
  }
  do.call(rbind, out_rows)
}

y_all <- make_binary_phenotype(pc_data$IID, "lactose_intolerance")
valid <- !is.na(y_all)
y <- y_all[valid]
Z_all <- genotype_all[valid, , drop = FALSE]
covariates <- pc_data[valid, pc_names, drop = FALSE]
set_list <- list(ALL_RARE = rare_all$SNP, RS_ONLY = rare_rs$SNP)
all_results <- list()
for (model_label in c("PLAIN", "PC_ADJUSTED")) {
  for (set_label in names(set_list)) {
    set_ids <- set_list[[set_label]]
    Z <- Z_all[, set_ids, drop = FALSE]
    all_results[[paste(model_label, set_label, sep = "__")]] <- run_one_set(Z, set_label, set_ids, y, covariates, model_label, length(all_results) * 100L)
  }
}
results <- do.call(rbind, all_results)
write.table(results, file.path(output_root, "lactose_region_burden_skat_skato_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(all_results, file.path(output_root, "lactose_region_burden_skat_skato_objects.rds"))
summary <- data.frame(PHENOTYPE = "lactose_intolerance", MAF_THRESHOLD = 0.01, N_FULL_SAMPLE = length(y), N_CASE = sum(y == 1L), N_CONTROL = sum(y == 0L), N_RARE_ALL = nrow(rare_all), N_RARE_RS_ONLY = nrow(rare_rs), N_RESAMPLING = n_resampling, stringsAsFactors = FALSE)
write.table(summary, file.path(output_root, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c("LCT_region / lactose_intolerance rare-variant set tests", "Primary set: every region variant with MAF < 0.01 using data/plink_LCT.afreq; rs-only set is a sensitivity analysis.", "Tests: weighted burden (beta(1,25)), SKAT, and SKAT-O via SKATBinary Hybrid calibration.", "Models: full-sample plain/no-PC reference and PC1-PC5-adjusted primary model; sex is not included.", paste0("Binary resampling budget N.Resampling = ", n_resampling, "; missing genotypes use best-guess imputation with missing_cutoff = 0.15."), "The current panel is called LCT_region rather than LCT gene because variant coordinates/functional annotation are not present in the supplied AF file."), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message_time("Completed LCT-region lactose burden/SKAT/SKAT-O tests: ", output_root)
