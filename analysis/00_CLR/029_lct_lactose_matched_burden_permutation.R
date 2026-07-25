#!/usr/bin/env Rscript

# Permutation-calibrated rare-variant burden test for the LCT/lactose positive
# control.  The case/control labels are permuted within the existing PC-only
# matched strata, preserving exactly one case and the requested 1:1, 1:4, or
# 1:6 control ratio.  Sex is not read or used.

options(stringsAsFactors = FALSE)
setTimeLimit(elapsed = 10 * 60 * 60, transient = FALSE)

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1L) {
  script_file <- sub("^--file=", "", script_file)
  project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

data_dir <- file.path(project_root, "data")
set_root <- file.path(project_root, "results", "00_CLR", "LCT_region_lactose_region_set_tests")
out_root <- file.path(set_root, "matched_burden_permutation")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

n_perm <- 10000L
n_pcs <- 5L
pc_names <- paste0("PC", seq_len(n_pcs))

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
read_raw_selected <- function(path, columns) {
  hdr <- names(read.delim(path, nrows = 0, check.names = FALSE, stringsAsFactors = FALSE))
  missing <- setdiff(columns, hdr)
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             colClasses = ifelse(hdr %in% columns, ifelse(hdr == "IID", "character", "numeric"), "NULL"))
}

af <- read_tab(file.path(data_dir, "plink_LCT.afreq"))
af$MAF <- pmin(as.numeric(af$ALT_FREQS), 1 - as.numeric(af$ALT_FREQS))
af$ID <- as.character(af$ID)
af$SNP <- paste0(af$ID, "_", af$REF)
raw_header <- names(read.delim(file.path(data_dir, "LCT_region_geno.raw"), nrows = 0, check.names = FALSE, stringsAsFactors = FALSE))
region_variants <- raw_header[grepl("^(rs|Affx-)", raw_header)]
all_rare <- af$SNP[is.finite(af$MAF) & af$MAF < 0.01 & af$SNP %in% region_variants]
all_rare <- region_variants[region_variants %in% all_rare]
rs_rare <- all_rare[grepl("^rs", all_rare)]
if (!length(all_rare) || !length(rs_rare)) stop("No LCT rare variants found")

geno <- read_raw_selected(file.path(data_dir, "LCT_region_geno.raw"), c("IID", all_rare))
if (anyDuplicated(geno$IID)) stop("Duplicate IID in genotype file")
G <- as.matrix(geno[, all_rare, drop = FALSE]); storage.mode(G) <- "double"
# SKAT's best-guess convention: rare missing genotypes are almost always 0;
# use the rounded observed allele-count mean for a deterministic imputation.
for (j in seq_len(ncol(G))) {
  m <- mean(G[, j], na.rm = TRUE)
  if (!is.finite(m)) m <- 0
  G[!is.finite(G[, j]), j] <- round(m)
}
burden_all <- rowSums(G[, all_rare, drop = FALSE])
burden_rs <- rowSums(G[, rs_rare, drop = FALSE])

match_paths <- c(
  `1to1` = file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison", "matching_1to1", "lactose_intolerance", "matched_sets.tsv"),
  `1to4` = file.path(project_root, "results", "00_CLR", "pc_adjusted_vs_matched_clr", "matching", "lactose_intolerance", "matched_sets.tsv"),
  `1to6` = file.path(project_root, "results", "00_CLR", "LCT_TV_qualification_comparison", "matching_1to6", "lactose_intolerance", "matched_sets.tsv")
)

permute_burden <- function(x, ratio, matched, iid, n_perm, seed, set_name) {
  if (any(grepl("sex", names(matched), ignore.case = TRUE))) stop("Sex information entered matched data")
  idx <- match(as.character(matched$IID), as.character(iid))
  if (anyNA(idx)) stop("Matched IID missing in LCT genotype data")
  x <- as.numeric(x[idx])
  strata_ids <- unique(matched$STRATUM)
  strata <- lapply(strata_ids, function(s) which(matched$STRATUM == s))
  valid <- vapply(strata, function(ix) {
    ci <- ix[matched$case[ix] == 1L]; ni <- ix[matched$case[ix] == 0L]
    length(ci) == 1L && length(ni) == ratio && all(is.finite(x[ix]))
  }, logical(1L))
  strata <- strata[valid]
  if (!length(strata)) return(data.frame(SET = set_name, RATIO = ratio, STATUS = "no_valid_strata", N_STRATA_VALID = 0L, N_INFORMATIVE_STRATA = 0L, N_CASE = NA_integer_, N_CONTROL = NA_integer_, N_PERM = n_perm, OBS_SCORE_CHISQ = NA_real_, P_SCORE_ASYM = NA_real_, P_PERM = NA_real_, CHISQ_PERM = NA_real_))
  n_case_valid <- length(strata)
  n_control_valid <- length(strata) * ratio
  xlist <- lapply(strata, function(ix) x[ix])
  means <- vapply(xlist, mean, numeric(1L)); variances <- vapply(xlist, function(v) mean(v^2) - mean(v)^2, numeric(1L))
  informative <- variances > 0 & is.finite(variances)
  if (!any(informative)) return(data.frame(SET = set_name, RATIO = ratio, STATUS = "no_within_stratum_variation", N_STRATA_VALID = length(strata), N_INFORMATIVE_STRATA = 0L, N_CASE = n_case_valid, N_CONTROL = n_control_valid, N_PERM = n_perm, OBS_SCORE_CHISQ = NA_real_, P_SCORE_ASYM = NA_real_, P_PERM = NA_real_, CHISQ_PERM = NA_real_))
  xlist <- xlist[informative]; means <- means[informative]; strata <- strata[informative]
  observed_case <- vapply(strata, function(ix) x[ix[matched$case[ix] == 1L]], numeric(1L))
  u_obs <- sum(observed_case - means); v_obs <- sum(variances); obs <- u_obs^2 / v_obs
  count_ge <- 0L; done <- 0L; set.seed(seed); n_strata <- length(xlist)
  groups <- split(seq_len(n_strata), lengths(xlist))
  chunk <- 200L
  while (done < n_perm) {
    b <- min(chunk, n_perm - done); u_perm <- numeric(b)
    for (gi in groups) {
      m <- length(xlist[[gi[1L]]]); xm <- do.call(rbind, xlist[gi])
      draw <- matrix(sample.int(m, length(gi) * b, replace = TRUE), nrow = length(gi), ncol = b)
      sel <- xm[cbind(rep(seq_len(length(gi)), times = b), as.vector(draw))]
      u_perm <- u_perm + colSums(matrix(sel, nrow = length(gi), ncol = b)) - sum(means[gi])
    }
    count_ge <- count_ge + sum((u_perm^2 / v_obs) >= obs - 1e-12); done <- done + b
  }
  p <- (1 + count_ge) / (n_perm + 1)
  data.frame(SET = set_name, RATIO = ratio, STATUS = "ok", N_STRATA_VALID = sum(valid), N_INFORMATIVE_STRATA = n_strata, N_CASE = n_case_valid, N_CONTROL = n_control_valid, N_PERM = n_perm, OBS_SCORE_CHISQ = obs, P_SCORE_ASYM = pchisq(obs, 1, lower.tail = FALSE), P_PERM = p, CHISQ_PERM = qchisq(1 - p, 1), stringsAsFactors = FALSE)
}

out <- list(); k <- 0L
for (ratio_name in names(match_paths)) {
  ratio <- as.integer(sub("1to", "", ratio_name)); path <- match_paths[[ratio_name]]
  if (!file.exists(path)) stop("Missing matched set file: ", path)
  matched <- read_tab(path)
  counts <- table(matched$STRATUM, matched$case)
  if (!all(c("0", "1") %in% colnames(counts)) || any(counts[, "1"] != 1L) || any(counts[, "0"] != ratio)) stop("Invalid matched sets for ", ratio_name)
  for (set_name in c("MAF_0.01_ALL", "MAF_0.01_RS")) {
    k <- k + 1L; x <- if (set_name == "MAF_0.01_ALL") burden_all else burden_rs
    message(sprintf("LCT/lactose burden permutation %s, %s", ratio_name, set_name))
    out[[k]] <- permute_burden(x, ratio, matched, geno$IID, n_perm, 29000L + ratio * 100L + k, set_name)
  }
}
res <- do.call(rbind, out)
write.table(res, file.path(out_root, "matched_burden_permutation_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

main <- read_tab(file.path(set_root, "lactose_region_burden_skat_skato_results.tsv"))
pc_main <- main[main$MODEL == "PC_ADJUSTED" & main$VARIANT_SET %in% c("ALL_RARE", "RS_ONLY") & main$METHOD == "BURDEN", c("VARIANT_SET", "P_VALUE"), drop = FALSE]
names(pc_main) <- c("SET", "P_PC_ADJUSTED_BURDEN")
pc_main$SET <- ifelse(pc_main$SET == "ALL_RARE", "MAF_0.01_ALL", "MAF_0.01_RS")
comparison <- merge(res, pc_main, by = "SET", all.x = TRUE, sort = FALSE)
write.table(comparison, file.path(out_root, "matched_burden_permutation_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  plot_data <- rbind(data.frame(SET = comparison$SET, RATIO = comparison$RATIO, P = comparison$P_PERM, MODEL = paste0("Matched 1:", comparison$RATIO, " permutation")), data.frame(SET = pc_main$SET, RATIO = NA_integer_, P = pc_main$P_PC_ADJUSTED_BURDEN, MODEL = "Full sample + PCs"))
  plot_data$SET <- factor(plot_data$SET, levels = c("MAF_0.01_ALL", "MAF_0.01_RS"))
  plot_data$MODEL <- factor(plot_data$MODEL, levels = c("Matched 1:1 permutation", "Matched 1:4 permutation", "Matched 1:6 permutation", "Full sample + PCs"))
  p <- ggplot(plot_data, aes(x = MODEL, y = -log10(P), color = SET)) + geom_point(position = position_dodge(width = .5), size = 3) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") + labs(title = "LCT/lactose rare burden: matched permutation vs PC-adjusted", x = NULL, y = expression(-log[10](p)), color = "Variant set") + theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(out_root, "matched_burden_permutation_pvalues.pdf"), p, width = 8.5, height = 5.2, units = "in", device = "pdf", bg = "white")
  ggsave(file.path(out_root, "matched_burden_permutation_pvalues.png"), p, width = 8.5, height = 5.2, units = "in", dpi = 220, bg = "white")
}

writeLines(c(
  "LCT_region / lactose_intolerance rare burden permutation",
  "Primary set: MAF < 0.01; RS_ONLY is included as a sensitivity set.",
  "Burden score is the sum of rare genotype dosages after deterministic best-guess imputation (rounded observed dosage mean; rare missing values are generally 0).",
  "Within each PC-only matched stratum, one member is selected uniformly as case; exact 1:1, 1:4, or 1:6 ratio is preserved.",
  paste0("Permutations per set and ratio: ", n_perm, "; no sex variable; no permutation of full-sample models."),
  "P_PERM uses the +1 correction; CHISQ_PERM = qchisq(1-P_PERM, 1)."
), file.path(out_root, "METHODS.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
message("Completed LCT/lactose matched burden permutation analysis: ", out_root)
