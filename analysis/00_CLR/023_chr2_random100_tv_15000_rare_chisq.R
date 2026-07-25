#!/usr/bin/env Rscript

# Chromosome 2 random-100 rare-variant analysis for TV.
#
# This is the same validated design used for the LCT rare-variant analysis:
# PC1-PC5 only (no sex), a common 15,000-case matching pool with 1:1, 1:4,
# and 1:6 matching, and asymptotic chi-square(1) Wald/CLR statistics.  The
# permutation-calibrated matched analysis is run separately by script 024.

options(stringsAsFactors = FALSE)
setTimeLimit(elapsed = 10 * 60 * 60, transient = FALSE)
suppressPackageStartupMessages(library(survival))

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

data_dir <- file.path(project_root, "data")
region <- "chr2_random100"
phenotype <- "TV"
raw_path <- file.path(data_dir, "chr2_random100.raw")
source_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool")
output_root <- file.path(project_root, "results", "00_CLR", "CHR2_random100_TV_15000_rare")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

num_pcs <- 5L
pc_names <- paste0("PC", seq_len(num_pcs))
ratios <- c(1L, 4L, 6L)

message_time <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

load_function_definitions <- function(path, names_to_load) {
  parsed <- parse(path)
  for (expr in parsed) {
    if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
      name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (length(name) == 1L && name %in% names_to_load) eval(expr, envir = parent.frame())
    }
  }
  invisible(NULL)
}

load_function_definitions(
  file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"),
  c("read_header", "read_selected_tsv", "read_genotypes", "load_phenotype_ids", "make_binary_phenotype", "logistic_wald_gwas")
)
load_function_definitions(
  file.path(project_root, "analysis", "00_CLR", "010_pc_matching_helpers.R"),
  c("run_clr", "draw_original_style")
)
required <- c("read_header", "read_selected_tsv", "read_genotypes", "load_phenotype_ids", "make_binary_phenotype", "logistic_wald_gwas", "run_clr", "draw_original_style")
if (!all(vapply(required, exists, logical(1L), envir = environment(), inherits = FALSE))) stop("Failed to load GWAS helpers")

# Region-specific allele frequency definition.  All 100 random-region
# variants are rs-labelled in this data, but the explicit intersection keeps
# the rs-only convention of the existing pipeline auditable.
af <- read.delim(file.path(data_dir, "plink_chr2.afreq"), check.names = FALSE, stringsAsFactors = FALSE)
af$ALT_FREQS <- as.numeric(af$ALT_FREQS)
af$MAF <- pmin(af$ALT_FREQS, 1 - af$ALT_FREQS)
af$SNP <- paste0(af$ID, "_", af$REF)
geno_header <- read_header(raw_path)
rs_variants <- geno_header[grepl("^rs", geno_header)]
rare_all <- af[is.finite(af$MAF) & af$MAF < 0.01 & af$SNP %in% rs_variants, c("SNP", "ID", "REF", "ALT", "MAF"), drop = FALSE]
rare_all <- rare_all[match(rs_variants[rs_variants %in% rare_all$SNP], rare_all$SNP), , drop = FALSE]
if (nrow(rare_all) == 0L) stop("No MAF < 0.01 rs variants found in chr2_random100")
rare_ids <- as.character(rare_all$SNP)
write.table(rare_all, file.path(output_root, "rare_variant_list_maf_lt_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# PC-complete full sample and rare genotypes are aligned by IID before any
# phenotype filtering.  This avoids accidental row-order or sample leakage.
pc_data <- read_selected_tsv(file.path(data_dir, "UKBB_pca.eigenvec"), c("IID", pc_names), c("character", rep("numeric", num_pcs)))
pc_data <- pc_data[complete.cases(pc_data[, pc_names, drop = FALSE]), , drop = FALSE]
if (anyDuplicated(pc_data$IID)) stop("Duplicate IID in PCA data")
geno <- read_selected_tsv(raw_path, c("IID", rare_ids), c("character", rep("numeric", length(rare_ids))))
if (anyDuplicated(geno$IID)) stop("Duplicate IID in genotype data")
genotype_index <- match(pc_data$IID, geno$IID)
if (anyNA(genotype_index)) stop("Some PC IIDs are missing from chr2_random100 genotype data")
genotype_all <- as.matrix(geno[genotype_index, rare_ids, drop = FALSE])
storage.mode(genotype_all) <- "double"
stopifnot(identical(pc_data$IID, geno$IID[genotype_index]))
rm(geno, genotype_index)

y_all <- make_binary_phenotype(pc_data$IID, phenotype)
valid <- !is.na(y_all)
y <- y_all[valid]
genotype <- genotype_all[valid, , drop = FALSE]
pc_covariates <- pc_data[valid, pc_names, drop = FALSE]
no_covariates <- data.frame(row.names = seq_len(length(y)))
message_time("Full sample TV N = ", length(y), "; cases = ", sum(y == 1L), "; controls = ", sum(y == 0L), "; rare variants = ", length(rare_ids))

plain <- logistic_wald_gwas(genotype, no_covariates, y)
pc_adjusted <- logistic_wald_gwas(genotype, pc_covariates, y)
plain$REGION <- region; plain$PHENOTYPE <- phenotype; plain$MODEL <- "FULL_SAMPLE_PLAIN_LOGISTIC"
pc_adjusted$REGION <- region; pc_adjusted$PHENOTYPE <- phenotype; pc_adjusted$MODEL <- "FULL_SAMPLE_PC_ADJUSTED_LOGISTIC"
plain <- plain[, c("REGION", "PHENOTYPE", "MODEL", setdiff(names(plain), c("REGION", "PHENOTYPE", "MODEL"))), drop = FALSE]
pc_adjusted <- pc_adjusted[, c("REGION", "PHENOTYPE", "MODEL", setdiff(names(pc_adjusted), c("REGION", "PHENOTYPE", "MODEL"))), drop = FALSE]
write.table(plain, file.path(output_root, "logistic_plain_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(pc_adjusted, file.path(output_root, "logistic_pc_adjusted_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plain_ok <- plain[plain$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]
names(plain_ok) <- c("SNP", "CHISQ_PLAIN", "P_PLAIN", "N_PLAIN")
pc_ok <- pc_adjusted[pc_adjusted$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]
names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
comparison <- merge(plain_ok, pc_ok, by = "SNP", sort = FALSE)

# Reuse the already validated PC-only 15,000-case common pool.  Matching is
# phenotype/PC based and therefore identical for a new genotype region.
clr_tables <- list()
for (ratio in ratios) {
  matched_path <- file.path(source_root, "matching", paste0("ratio_1to", ratio), "matched_sets.tsv")
  matched <- read.delim(matched_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (any(grepl("sex", names(matched), ignore.case = TRUE))) stop("Sex entered matched data")
  counts <- table(matched$STRATUM, matched$case)
  if (!all(c("0", "1") %in% colnames(counts)) || any(counts[, "1"] != 1L) || any(counts[, "0"] != ratio) || anyDuplicated(matched$IID)) stop("Invalid matched sets for ratio 1:", ratio)
  matched_index <- match(matched$IID, pc_data$IID)
  if (anyNA(matched_index)) stop("Matched IID missing from PC/genotype matrix")
  genotype_for_ratio <- genotype_all[matched_index, , drop = FALSE]
  clr <- run_clr(genotype_for_ratio, matched, ratio)
  ratio_dir <- file.path(output_root, paste0("ratio_1to", ratio))
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(clr, file.path(ratio_dir, "matched_clr_gwas_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  clr_tables[[paste0("1to", ratio)]] <- clr
  ok <- clr[clr$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]
  names(ok) <- c("SNP", paste0("CHISQ_1TO", ratio), paste0("P_1TO", ratio), paste0("N_1TO", ratio), paste0("N_1TO", ratio, "_BASE"))
  comparison <- merge(comparison, ok, by = "SNP", sort = FALSE)
  message_time("Matched CLR 1:", ratio, " successful variants = ", sum(clr$STATUS == "ok"))
}
write.table(comparison, file.path(output_root, "comparison_rare_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (nrow(comparison) > 0L && requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  plot_columns <- c("CHISQ_1TO1", "CHISQ_1TO4", "CHISQ_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
  plot_labels <- c("Rare matched 1:1", "Rare matched 1:4", "Rare matched 1:6", "Rare full sample + PCs", "Rare full sample (No PCs)")
  full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
  ggsave(file.path(output_root, "qq_rare_15000_1to1_1to4_1to6_original_style_full.pdf"), draw_original_style(comparison, plot_columns, plot_labels, full_max), width = 19, height = 4.2, units = "in", device = "pdf")
  ggsave(file.path(output_root, "qq_rare_15000_1to1_1to4_1to6_original_style_zoom10.pdf"), draw_original_style(comparison, plot_columns, plot_labels, 10), width = 19, height = 4.2, units = "in", device = "pdf")
  ggsave(file.path(output_root, "qq_rare_15000_1to1_1to4_1to6_original_style_zoom10.png"), draw_original_style(comparison, plot_columns, plot_labels, 10), width = 19, height = 4.2, units = "in", dpi = 220)
}

lambda_or_na <- function(x) if (length(x) == 0L || all(!is.finite(x))) NA_real_ else median(x[is.finite(x)]) / qchisq(.5, 1)
read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
balance_max <- function(ratio) max(abs(read_tab(file.path(source_root, "matching", paste0("ratio_1to", ratio), "pc_balance.tsv"))$POST_SMD))
summary <- data.frame(
  REGION = region, PHENOTYPE = phenotype, MAF_THRESHOLD = 0.01,
  N_RARE_RS_VARIANTS = length(rare_ids), N_COMMON_MODELS = nrow(comparison),
  N_FULL_SAMPLE = length(y), N_FULL_SAMPLE_CASE = sum(y == 1L), N_FULL_SAMPLE_CONTROL = sum(y == 0L),
  N_PLAIN_OK = sum(plain$STATUS == "ok"), N_PC_ADJUSTED_OK = sum(pc_adjusted$STATUS == "ok"),
  N_1TO1_OK = sum(clr_tables$`1to1`$STATUS == "ok"), N_1TO4_OK = sum(clr_tables$`1to4`$STATUS == "ok"), N_1TO6_OK = sum(clr_tables$`1to6`$STATUS == "ok"),
  N_CASE_1TO1 = 15000L, N_CONTROL_1TO1 = 15000L, N_CASE_1TO4 = 15000L, N_CONTROL_1TO4 = 60000L, N_CASE_1TO6 = 15000L, N_CONTROL_1TO6 = 90000L,
  MAX_ABS_SMD_POST_1TO1 = balance_max(1), MAX_ABS_SMD_POST_1TO4 = balance_max(4), MAX_ABS_SMD_POST_1TO6 = balance_max(6),
  LAMBDA_1TO1 = lambda_or_na(comparison$CHISQ_1TO1), LAMBDA_1TO4 = lambda_or_na(comparison$CHISQ_1TO4), LAMBDA_1TO6 = lambda_or_na(comparison$CHISQ_1TO6), LAMBDA_PC_ADJUSTED = lambda_or_na(comparison$CHISQ_PC_ADJUSTED), LAMBDA_PLAIN = lambda_or_na(comparison$CHISQ_PLAIN),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(output_root, "analysis_summary_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c(
  "chr2_random100 / TV rare-variant analysis",
  "Region: data/chr2_random100.raw; rare definition: MAF < 0.01 using data/plink_chr2.afreq",
  "Full sample: plain logistic (case ~ SNP) and PC-adjusted logistic (case ~ PC1-PC5 + SNP), all available phenotype/genotype/PC-complete samples",
  "Matched CLR: reuses the validated 15,000-case PC-only matching sets from LCT_TV_15000_common_pool; ratios 1:1, 1:4, 1:6",
  "Sex was not read, matched, or modeled. No permutation is used in this script; script 024 supplies the within-stratum permutation calibration.",
  "QQ plots use the original ggplot2 facet style and theoretical chi-square(1) reference."
), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message_time("Completed chr2_random100 TV rare chi-square analysis: ", output_root, "; common variants = ", nrow(comparison))
