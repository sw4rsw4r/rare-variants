#!/usr/bin/env Rscript

# LCT_region / TV: rare-variant subset of the validated 15,000-case analysis.
# Rare is defined using the region-specific PLINK allele frequency file as
# MAF < 0.01. The current pipeline analyzes rs-labelled variants;
# Affymetrix probe IDs are reported but excluded for comparability with the
# existing full-sample rs-only outputs. No sex and no permutation.

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

data_dir <- file.path(project_root, "data")
source_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool")
output_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool_rare")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

load_function_definitions <- function(path, names_to_load) {
  parsed <- parse(path)
  for (expr in parsed) {
    if (is.call(expr) && length(expr) >= 3L && as.character(expr[[1L]]) %in% c("<-", "=")) {
      name <- if (is.symbol(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (length(name) == 1L && name %in% names_to_load) eval(expr, envir = parent.frame())
    }
  }
}
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "003_run_full_sample_logistic_gwas.R"), c("read_header"))
load_function_definitions(file.path(project_root, "analysis", "00_CLR", "010_pc_matching_helpers.R"), c("draw_original_style"))

af <- read.delim(file.path(data_dir, "plink_LCT.afreq"), check.names = FALSE, stringsAsFactors = FALSE)
af$MAF <- pmin(as.numeric(af$ALT_FREQS), 1 - as.numeric(af$ALT_FREQS))
af$SNP <- paste0(af$ID, "_", af$REF)
rare_all <- af[is.finite(af$MAF) & af$MAF < 0.01, c("SNP", "ID", "REF", "ALT", "MAF"), drop = FALSE]
geno_header <- read_header(file.path(data_dir, "LCT_region_geno.raw"))
rs_variants <- geno_header[grepl("^rs", geno_header)]
rare_rs <- rare_all[rare_all$SNP %in% rs_variants, , drop = FALSE]
rare_rs <- rare_rs[match(intersect(rare_all$SNP, rare_rs$SNP), rare_rs$SNP), , drop = FALSE]
excluded_non_rs <- rare_all[!rare_all$SNP %in% rs_variants, , drop = FALSE]
if (nrow(rare_rs) == 0L) stop("No rare rs variants found")
write.table(rare_rs, file.path(output_root, "rare_variant_list_maf_lt_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(excluded_non_rs, file.path(output_root, "excluded_non_rs_rare_variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
rare_ids <- rare_rs$SNP
subset_model <- function(tab, id_col = "SNP") tab[tab[[id_col]] %in% rare_ids, , drop = FALSE]
plain <- subset_model(read_tab(file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison", "TV", "logistic_plain.tsv")))
pc <- subset_model(read_tab(file.path(project_root, "results", "00_CLR", "full_sample_logistic_gwas", "LCT_region", "TV", "logistic_gwas.tsv")))
write.table(plain, file.path(output_root, "logistic_plain_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(pc, file.path(output_root, "logistic_pc_adjusted_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plain_ok <- plain[plain$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(plain_ok) <- c("SNP", "CHISQ_PLAIN", "P_PLAIN", "N_PLAIN")
pc_ok <- pc[pc$STATUS == "ok", c("SNP", "CHISQ", "P", "N"), drop = FALSE]; names(pc_ok) <- c("SNP", "CHISQ_PC_ADJUSTED", "P_PC_ADJUSTED", "N_PC_ADJUSTED")
comparison <- merge(plain_ok, pc_ok, by = "SNP", sort = FALSE)
clr_tables <- list()
for (ratio in c(1L, 4L, 6L)) {
  source_path <- file.path(source_root, paste0("ratio_1to", ratio), "matched_clr_gwas.tsv")
  clr_all <- read_tab(source_path)
  clr <- clr_all[clr_all$SNP %in% rare_ids, , drop = FALSE]
  clr <- clr[match(intersect(rare_ids, clr$SNP), clr$SNP), , drop = FALSE]
  ratio_dir <- file.path(output_root, paste0("ratio_1to", ratio))
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(clr, file.path(ratio_dir, "matched_clr_gwas_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  clr_tables[[paste0("1to", ratio)]] <- clr
  ok <- clr[clr$STATUS == "ok", c("SNP", "CHISQ", "P", "N", "N_MATCHED_BASE"), drop = FALSE]
  names(ok) <- c("SNP", paste0("CHISQ_1TO", ratio), paste0("P_1TO", ratio), paste0("N_1TO", ratio), paste0("N_1TO", ratio, "_BASE"))
  comparison <- merge(comparison, ok, by = "SNP", sort = FALSE)
}
write.table(comparison, file.path(output_root, "comparison_rare_1to1_1to4_1to6.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_columns <- c("CHISQ_1TO1", "CHISQ_1TO4", "CHISQ_1TO6", "CHISQ_PC_ADJUSTED", "CHISQ_PLAIN")
plot_labels <- c("Rare matched 1:1", "Rare matched 1:4", "Rare matched 1:6", "Rare full sample + PCs", "Rare full sample (No PCs)")
if (nrow(comparison) > 0L) {
  full_max <- max(c(unlist(comparison[plot_columns]), qchisq(ppoints(nrow(comparison)), 1)), na.rm = TRUE)
  plot_full <- draw_original_style(comparison, plot_columns, plot_labels, full_max)
  plot_zoom <- draw_original_style(comparison, plot_columns, plot_labels, 10)
  ggsave(file.path(output_root, "qq_rare_15000_1to1_1to4_1to6_original_style_full.pdf"), plot_full, width = 19, height = 4.2, units = "in", device = "pdf")
  ggsave(file.path(output_root, "qq_rare_15000_1to1_1to4_1to6_original_style_zoom10.pdf"), plot_zoom, width = 19, height = 4.2, units = "in", device = "pdf")
  ggsave(file.path(output_root, "qq_rare_15000_1to1_1to4_1to6_original_style_zoom10.png"), plot_zoom, width = 19, height = 4.2, units = "in", dpi = 220)
}

lambda_or_na <- function(x) if (length(x) == 0L || all(!is.finite(x))) NA_real_ else median(x[is.finite(x)]) / qchisq(.5, 1)
balance_max <- function(ratio) max(abs(read_tab(file.path(source_root, "matching", paste0("ratio_1to", ratio), "pc_balance.tsv"))$POST_SMD))
summary <- data.frame(
  REGION = "LCT_region", PHENOTYPE = "TV", MAF_THRESHOLD = 0.01,
  N_RARE_REGION_VARIANTS = nrow(rare_all), N_RARE_RS_VARIANTS = nrow(rare_rs), N_EXCLUDED_NON_RS = nrow(excluded_non_rs),
  N_COMMON_MODELS = nrow(comparison), N_PLAIN_OK = sum(plain$STATUS == "ok"), N_PC_ADJUSTED_OK = sum(pc$STATUS == "ok"),
  N_1TO1_OK = sum(clr_tables$`1to1`$STATUS == "ok"), N_1TO4_OK = sum(clr_tables$`1to4`$STATUS == "ok"), N_1TO6_OK = sum(clr_tables$`1to6`$STATUS == "ok"),
  N_CASE_1TO1 = 15000L, N_CONTROL_1TO1 = 15000L, N_CASE_1TO4 = 15000L, N_CONTROL_1TO4 = 60000L, N_CASE_1TO6 = 15000L, N_CONTROL_1TO6 = 90000L,
  MAX_ABS_SMD_POST_1TO1 = balance_max(1), MAX_ABS_SMD_POST_1TO4 = balance_max(4), MAX_ABS_SMD_POST_1TO6 = balance_max(6),
  LAMBDA_1TO1 = lambda_or_na(comparison$CHISQ_1TO1), LAMBDA_1TO4 = lambda_or_na(comparison$CHISQ_1TO4), LAMBDA_1TO6 = lambda_or_na(comparison$CHISQ_1TO6), LAMBDA_PC_ADJUSTED = lambda_or_na(comparison$CHISQ_PC_ADJUSTED), LAMBDA_PLAIN = lambda_or_na(comparison$CHISQ_PLAIN),
  LAMBDA_PC_ADJUSTED_ALL_OK = lambda_or_na(pc$CHISQ[pc$STATUS == "ok"]), LAMBDA_PLAIN_ALL_OK = lambda_or_na(plain$CHISQ[plain$STATUS == "ok"]),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(output_root, "analysis_summary_rare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c(
  "LCT_region / TV rare-variant subset of the validated 15,000-case common-pool analysis",
  "Rare definition: MAF < 0.01 using data/plink_LCT.afreq",
  paste0("Rare region variants = ", nrow(rare_all), "; rs-labelled variants analyzed = ", nrow(rare_rs), "; non-rs rare probes excluded = ", nrow(excluded_non_rs)),
  "The existing no-sex, no-permutation full-sample and PC-only matched model results were subset to the rare rs variants; no new matching was needed.",
  paste0("The QQ comparison uses the ", nrow(comparison), " rare variants successful in every model; all-model lambda values are also reported for the full-sample models."),
  "For matched CLR, variants with no within-stratum genotype variation are reported but excluded from the QQ comparison.",
  "Lambda is descriptive because only a small regional rare-variant set is tested."
), file.path(output_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message(sprintf("Completed rare-variant subset: %s; rare rs=%d; common successful models=%d", output_root, nrow(rare_rs), nrow(comparison)))
