#!/usr/bin/env Rscript

# Recommended final manuscript package.  This script reads completed GWAS,
# PC-adjusted, CLR, matching, and permutation outputs; it does not overwrite
# the existing CLR analysis scripts or refit models.
#
# Rare is defined strictly as MAF < 0.01 (1%).  The wording "AF < 1" is
# interpreted as 1 percent; AF < 1.0 would classify every biallelic variant.

options(stringsAsFactors = FALSE)
for (pkg in c("ggplot2", "patchwork", "gridExtra")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required")
}
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(gridExtra); library(grid)
})

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_file) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_file)), "..", ".."), winslash = "/", mustWork = TRUE)
} else normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
real_root <- file.path(project_root, "results", "00_CLR")
sim_root <- Sys.getenv("SIM_ROOT_OVERRIDE", file.path(dirname(project_root), "R simulation", "calibrated_simulation_pc1to5_results"))
sim_is_calibrated <- grepl("calibrated_simulation", basename(normalizePath(sim_root)), fixed = TRUE)
sim_design_label <- if (sim_is_calibrated) "real-data-aligned calibration sensitivity" else "original simulation design"
out_root <- Sys.getenv("PAPER_PACKAGE_OUT_ROOT", file.path(real_root, "paper_package_recommended"))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
old <- list.files(out_root, pattern = "^(Figure|Table1|Supplementary|README_paper_package|Figure_legends_and_results_draft|Methods|sessionInfo)", full.names = TRUE)
if (length(old)) unlink(old, recursive = FALSE, force = TRUE)

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
fmt_p <- function(x) ifelse(is.finite(x) & x > 0, formatC(x, format = "e", digits = 2), NA_character_)
fmt_num <- function(x, digits = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
write_table_pdf <- function(data, path, title, landscape = TRUE, font_size = 7) {
  pdf(path, width = if (landscape) 13 else 10, height = if (landscape) 8.5 else 11, bg = "white", useDingbats = FALSE)
  tg <- gridExtra::tableGrob(data, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = font_size, padding = unit(c(2.5, 2.5), "mm")))
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  grid.draw(gridExtra::arrangeGrob(tg, top = textGrob(title, gp = gpar(fontface = "bold", fontsize = 12))))
  dev.off()
}

paper_theme <- theme_bw(base_size = 10.5) + theme(
  panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank(),
  plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey35"),
  plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA)
)
method_levels <- c("Plain logistic", "PC-adjusted logistic", "CLR 1:6", "CLR 1:4", "CLR 1:1")
ratio_levels <- c("1:6", "1:4", "1:1")
method_cols <- c("Plain logistic" = "#D55E00", "PC-adjusted logistic" = "#0072B2", "CLR 1:6" = "#E69F00", "CLR 1:4" = "#CC79A7", "CLR 1:1" = "#009E73")
p_map <- c("Plain logistic" = "P_PLAIN", "PC-adjusted logistic" = "P_PC_ADJUSTED", "CLR 1:6" = "P_1TO6", "CLR 1:4" = "P_1TO4", "CLR 1:1" = "P_1TO1")
chisq_map <- c("Plain logistic" = "CHISQ_PLAIN", "PC-adjusted logistic" = "CHISQ_PC_ADJUSTED", "CLR 1:6" = "CHISQ_1TO6", "CLR 1:4" = "CHISQ_1TO4", "CLR 1:1" = "CHISQ_1TO1")

## Allele-frequency annotation: strict MAF < 0.01 (1%) rare definition.
af <- read_tab(file.path(project_root, "data", "plink_LCT.afreq"))
af$MAF <- pmin(as.numeric(af$ALT_FREQS), 1 - as.numeric(af$ALT_FREQS))
af$key <- paste0(af$ID, "_", af$REF)
maf_lookup <- setNames(af$MAF, af$key)
add_af_class <- function(d) {
  d$MAF <- as.numeric(maf_lookup[as.character(d$SNP)])
  d$AF_class <- ifelse(is.finite(d$MAF) & d$MAF < 0.01, "Rare (MAF < 1%)", "Common (MAF >= 1%)")
  d$AF_class[!is.finite(d$MAF)] <- NA_character_
  d
}

## Generic QQ plotting on -log10(p).
qq_long <- function(data, panel, map = p_map) {
  out <- lapply(names(map), function(method) {
    p <- suppressWarnings(as.numeric(data[[map[[method]]]]))
    p <- sort(p[is.finite(p) & p > 0 & p <= 1])
    if (!length(p)) return(NULL)
    data.frame(panel = panel, method = method, expected = -log10(ppoints(length(p))), observed = -log10(p), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  out$method <- factor(out$method, levels = method_levels)
  out
}
qq_plot <- function(data, title, subtitle, lim = 8, xlab = "Expected -log10(p)\nTheoretical Uniform(0,1) / chi-square(1)") {
  ggplot(data, aes(expected, observed, colour = method, group = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = .45, alpha = .8) + geom_point(size = .65, alpha = .55) +
    coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) + scale_colour_manual(values = method_cols, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "Observed -log10(p)") +
    paper_theme + theme(legend.position = "none")
}

## Lambda-GC and paper-matched permutation GIF calculations.
lambda_gc_chisq <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x) / qchisq(.5, 1)
}
lambda_gc_p <- function(p) {
  p <- as.numeric(p); p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  lambda_gc_chisq(qchisq(p, 1, lower.tail = FALSE))
}
gif_qqperm <- function(p_observed, p_expected) {
  keep <- is.finite(p_observed) & is.finite(p_expected) & p_observed > 0 & p_observed <= 1 & p_expected > 0 & p_expected <= 1
  if (sum(keep) < 3L) return(NA_real_)
  y <- sort(qchisq(p_observed[keep], 1, lower.tail = FALSE))
  x <- sort(qchisq(p_expected[keep], 1, lower.tail = FALSE))
  if (!any(x > 0)) return(NA_real_)
  unname(coef(lm(y ~ 0 + x))[1L])
}

## Sources for real-data analyses.
lct_tv <- add_af_class(read_tab(file.path(real_root, "LCT_TV_15000_common_pool", "comparison_1to1_1to4_1to6.tsv")))
lct_tv_summary <- read_tab(file.path(real_root, "LCT_TV_15000_common_pool", "analysis_summary.tsv"))
lct_tv_perm <- read_tab(file.path(real_root, "LCT_TV_15000_common_pool_rare_permutation", "comparison_rare_permutation_1to1_1to4_1to6.tsv"))
lactose <- add_af_class(read_tab(file.path(real_root, "LCT_TV_qualification_comparison", "lactose_intolerance", "positive_control_comparison_1to1_1to4_1to6.tsv")))
chr2_asym <- read_tab(file.path(real_root, "CHR2_random100_TV_15000_rare", "comparison_rare_1to1_1to4_1to6.tsv"))

## Pre-specified positive-control sentinel.  The PLINK dosage column is
## rs4988235_A, so all effects below are per copy of the A allele.
sentinel_paths <- c(
  "Plain logistic" = file.path(real_root, "LCT_TV_qualification_comparison", "lactose_intolerance", "logistic_plain.tsv"),
  "PC-adjusted logistic" = file.path(real_root, "LCT_TV_qualification_comparison", "lactose_intolerance", "logistic_pc_adjusted.tsv"),
  "CLR 1:6" = file.path(real_root, "LCT_TV_qualification_comparison", "lactose_intolerance", "matched_clr_1to6_gwas.tsv"),
  "CLR 1:4" = file.path(real_root, "pc_adjusted_vs_matched_clr", "LCT_region", "lactose_intolerance", "matched_clr_gwas.tsv"),
  "CLR 1:1" = file.path(real_root, "LCT_TV_qualification_comparison", "lactose_intolerance", "matched_clr_1to1_gwas.tsv")
)
sentinel_stats <- do.call(rbind, lapply(names(sentinel_paths), function(method) {
  d <- read_tab(sentinel_paths[[method]])
  z <- d[d$SNP == "rs4988235_A" & d$STATUS == "ok", , drop = FALSE]
  if (nrow(z) != 1L) stop("Expected one valid rs4988235_A row for ", method)
  data.frame(
    Method = method,
    SNP = "rs4988235_A",
    Effect_allele = "A",
    N = as.numeric(z$N),
    N_case = as.numeric(z$N_CASE),
    N_control = as.numeric(z$N_CONTROL),
    N_genotype_missing = as.numeric(z$N_GENOTYPE_MISSING),
    EAF = as.numeric(z$EAF),
    Effect_allele_count = round(2 * as.numeric(z$N) * as.numeric(z$EAF)),
    MAC = as.numeric(z$MAC),
    Beta = as.numeric(z$BETA),
    SE = as.numeric(z$SE),
    P = as.numeric(z$P),
    stringsAsFactors = FALSE
  )
}))
sentinel_stats$OR <- exp(sentinel_stats$Beta)
sentinel_stats$OR_lower <- exp(sentinel_stats$Beta - 1.96 * sentinel_stats$SE)
sentinel_stats$OR_upper <- exp(sentinel_stats$Beta + 1.96 * sentinel_stats$SE)
sentinel_stats$Method <- factor(sentinel_stats$Method, levels = method_levels)
sentinel_stats <- sentinel_stats[order(sentinel_stats$Method), , drop = FALSE]

## Main Figure 1: control analyses (TV negative controls, lactose positive control).
p_lct_tv_rare <- lct_tv[lct_tv$AF_class == "Rare (MAF < 1%)", , drop = FALSE]
lct_rare_perm <- merge(p_lct_tv_rare[, c("SNP", "P_PLAIN", "P_PC_ADJUSTED")], lct_tv_perm[, c("SNP", "P_PERM_1TO1", "P_PERM_1TO4", "P_PERM_1TO6")], by = "SNP", sort = FALSE)
names(lct_rare_perm)[names(lct_rare_perm) == "P_PERM_1TO1"] <- "P_1TO1"; names(lct_rare_perm)[names(lct_rare_perm) == "P_PERM_1TO4"] <- "P_1TO4"; names(lct_rare_perm)[names(lct_rare_perm) == "P_PERM_1TO6"] <- "P_1TO6"
p_lct_tv <- qq_plot(qq_long(lct_tv, "LCT TV all"), "LCT-region TV negative control (all variants)", "52 variants; common + rare; asymptotic matched p-values")
p_chr2 <- qq_plot(qq_long(chr2_asym, "chr2 random rare TV"), "Chromosome-2 random rare TV", "100 rare variants; observed chi-square asymptotic p-values")
p_lactose <- qq_plot(qq_long(lactose, "LCT lactose"), "LCT-region lactose positive control", "40 common + 6 rare variants; zoomed to -log10(p) <= 8")

# For the lambda panel, use the observed chi-square statistics directly for
# every method.  This makes LCT rare and chr2 rare comparable on the same
# chi-square(1) scale.  Permutation-p-value-based GIF calibration remains in
# Supplementary Figure S3.
lambda_from_chisq <- function(d, cols) vapply(cols, function(z) lambda_gc_chisq(d[[z]]), numeric(1))
chisq_cols_lct <- unname(chisq_map[method_levels])
chisq_cols_chr2 <- unname(chisq_map[method_levels])
chr2_lambda <- lambda_from_chisq(chr2_asym, chisq_cols_chr2)
lct_all_lambda <- lambda_from_chisq(lct_tv, chisq_cols_lct)
lambda_summary <- data.frame(Control = c(rep("LCT TV all", 5), rep("chr2 random rare TV", 5)), Method = factor(rep(method_levels, 2), levels = method_levels), lambda_GC = c(lct_all_lambda, chr2_lambda), stringsAsFactors = FALSE)
lambda_summary$Control <- factor(lambda_summary$Control, levels = c("LCT TV all", "chr2 random rare TV"))
p_lambda <- ggplot(lambda_summary, aes(Method, lambda_GC, fill = Control)) + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") + geom_col(position = position_dodge(width = .8), width = .7, show.legend = FALSE) + facet_wrap(~Control, nrow = 1, scales = "free_y") + scale_y_log10(breaks = c(.5, 1, 2, 5, 10, 50, 500), labels = c("0.5", "1", "2", "5", "10", "50", "500")) + labs(title = "Negative-control null calibration", subtitle = "Facet-specific y-axis scales; direct chi-square lambda", x = "Increasing matching stringency ->", y = expression(lambda[GC])) + paper_theme + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

sentinel_plot <- sentinel_stats
sentinel_plot$Method <- factor(as.character(sentinel_plot$Method), levels = rev(method_levels))
sentinel_plot$P_label <- paste0("p = ", formatC(sentinel_plot$P, format = "e", digits = 2))
sentinel_plot$minus_log10_P <- -log10(sentinel_plot$P)
p_lead_effect <- ggplot(sentinel_plot, aes(OR, Method, colour = Method)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = OR_lower, xmax = OR_upper), orientation = "y", width = .16, linewidth = .65) +
  geom_point(size = 2.8) +
  scale_x_log10(breaks = c(.5, .6, .7, .8, 1), labels = c("0.5", "0.6", "0.7", "0.8", "1.0"), limits = c(.5, 1.05)) +
  scale_colour_manual(values = method_cols, drop = FALSE) +
  labs(
    x = "OR per A allele (95% CI)",
    y = NULL
  ) +
  paper_theme +
  theme(plot.margin = margin(5.5, 5.5, 22, 5.5)) +
  guides(colour = "none")
p_lead_p <- ggplot(sentinel_plot, aes(minus_log10_P, Method, fill = Method)) +
  geom_vline(xintercept = -log10(.05), linetype = "dotted", colour = "grey45", linewidth = .55) +
  geom_vline(xintercept = -log10(5e-8), linetype = "dashed", colour = "grey35", linewidth = .55) +
  geom_col(width = .62, show.legend = FALSE) +
  geom_text(aes(label = P_label), x = 11.8, hjust = 1, colour = "grey15", size = 2.45, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 12.2), expand = expansion(mult = c(0, .01))) +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  labs(x = expression(-log[10](italic(p))), y = NULL) +
  paper_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(5.5, 5.5, 22, 5.5)) +
  guides(fill = "none")
p_lead_core <- (p_lead_effect | p_lead_p) +
  plot_layout(widths = c(1.12, 1)) +
  plot_annotation(
    title = "Positive-control effect retention and association strength",
    subtitle = "Top-to-bottom: stronger matching; rs4988235-A sentinel; dotted p=.05, dashed p=5e-8"
  )
p_lead <- wrap_elements(full = p_lead_core)
fig_control <- ((p_lct_tv | p_chr2 | p_lactose) / (p_lambda | p_lead)) + plot_layout(heights = c(1.25, .85), guides = "collect") + plot_annotation(title = "Real-data control analyses: LCT, chromosome 2, and lactose", subtitle = "Method order: plain, PC-adjusted, CLR 1:6, 1:4, 1:1 (increasing matching stringency)") & theme(legend.position = "bottom")
ggsave(file.path(out_root, "Figure1_real_data_controls.pdf"), fig_control, width = 12.6, height = 9.4, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Figure1_real_data_controls.png"), fig_control, width = 12.6, height = 9.4, units = "in", dpi = 300, bg = "white")

## Main Figure 2: simulation operating characteristics (existing R simulation project).
sim <- read.csv(file.path(sim_root, "summary_metrics.csv"), check.names = FALSE)
sim_settings <- read.csv(file.path(sim_root, "simulation_settings.csv"), check.names = FALSE)
sim_n_pcs <- as.integer(sim_settings$value[match("n_pcs", sim_settings$parameter)])
if (!is.finite(sim_n_pcs) || sim_n_pcs < 1L) stop("Simulation settings must report a positive n_pcs")
sim_pc_names <- paste0("PC", seq_len(sim_n_pcs))
sim_pc_label <- paste0("PC1-PC", sim_n_pcs)
sim$scenario <- as.character(sim$scenario); sim$method <- as.character(sim$method); sim$ratio_label <- factor(as.character(sim$ratio_label), levels = ratio_levels)
sim$median_chisq_ratio <- sim$lambda_gc
sim$metric_interpretation <- ifelse(sim$scenario == "Negative control", "lambda_GC calibration", "positive-control signal strength (not calibration)")
sim_method <- c("Oracle GWAS", "Standard GWAS", "PC-adjusted GWAS", "Matching + CLR")
sim$method <- factor(sim$method, levels = sim_method)
sim_colours <- c("Oracle GWAS" = "#000000", "Standard GWAS" = "#D55E00", "PC-adjusted GWAS" = "#0072B2", "Matching + CLR" = "#009E73")
sim_shapes <- c("Oracle GWAS" = 18, "Standard GWAS" = 15, "PC-adjusted GWAS" = 17, "Matching + CLR" = 16)
sim_neg <- sim[sim$scenario == "Negative control", , drop = FALSE]
sim_pos <- sim[sim$scenario == "Positive control", , drop = FALSE]
sim_legend_labels <- c(
  "Oracle GWAS" = "Oracle GWAS",
  "Standard GWAS" = "Standard GWAS\u2020",
  "PC-adjusted GWAS" = "PC-adjusted GWAS\u2020",
  "Matching + CLR" = "Matching + CLR"
)
type1_plot <- ggplot(sim_neg, aes(ratio_label, rejection_rate, ymin = ci_lower, ymax = ci_upper, colour = method, shape = method, group = method)) +
  geom_hline(yintercept = .05, linetype = "dashed", colour = "grey35") +
  geom_line(linewidth = .55) + geom_errorbar(width = .1, linewidth = .45) + geom_point(size = 2.1) +
  scale_colour_manual(values = sim_colours, labels = sim_legend_labels, drop = FALSE) +
  scale_shape_manual(values = sim_shapes, labels = sim_legend_labels, drop = FALSE) +
  scale_y_log10(breaks = c(.025, .05, .1, .2, .5, 1), labels = scales::label_percent(accuracy = .1), limits = c(.025, 1.05), expand = expansion(mult = c(.02, .04))) +
  labs(title = "Negative control", subtitle = "Type I error rate (log scale); dashed = 5%", x = "Controls per case\nIncreasing matching stringency ->", y = "Type I error rate") + paper_theme
power_plot <- ggplot(sim_pos, aes(ratio_label, rejection_rate, ymin = ci_lower, ymax = ci_upper, colour = method, shape = method, group = method)) +
  geom_line(linewidth = .55) + geom_errorbar(width = .1, linewidth = .45) + geom_point(size = 2.1) +
  scale_colour_manual(values = sim_colours, labels = sim_legend_labels, drop = FALSE) +
  scale_shape_manual(values = sim_shapes, labels = sim_legend_labels, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1.02), expand = expansion(mult = c(.01, .04))) +
  labs(title = "Positive control", subtitle = "Statistical power", x = "Controls per case\nIncreasing matching stringency ->", y = "Statistical power") + paper_theme
lambda_plot <- ggplot(sim_neg, aes(ratio_label, lambda_gc, colour = method, shape = method, group = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35") + geom_line(linewidth = .55) + geom_point(size = 2.2) +
  scale_colour_manual(values = sim_colours, labels = sim_legend_labels, drop = FALSE) +
  scale_shape_manual(values = sim_shapes, labels = sim_legend_labels, drop = FALSE) +
  scale_y_log10(breaks = c(.5, 1, 2, 4, 8, 20, 80), labels = c("0.5", "1", "2", "4", "8", "20", "80")) +
  labs(title = "Negative-control calibration", subtitle = "Genomic-control inflation factor; dashed target = 1", x = "Controls per case\nIncreasing matching stringency ->", y = expression(lambda[GC])) +
  paper_theme + guides(colour = "none", shape = "none")
positive_chisq_plot <- ggplot(sim_pos, aes(ratio_label, median_chisq_ratio, colour = method, shape = method, group = method)) +
  geom_line(linewidth = .55) + geom_point(size = 2.2) +
  scale_colour_manual(values = sim_colours, labels = sim_legend_labels, drop = FALSE) +
  scale_shape_manual(values = sim_shapes, labels = sim_legend_labels, drop = FALSE) +
  scale_y_log10(breaks = c(5, 10, 20, 50), labels = c("5", "10", "20", "50")) +
  labs(title = "Positive-control signal strength", subtitle = "Median chi-square / 0.4549 (log scale); descriptive, not calibration", x = "Controls per case\nIncreasing matching stringency ->", y = "Median chi-square ratio") +
  paper_theme + guides(colour = "none", shape = "none")
fig_sim <- ((type1_plot | power_plot) / (lambda_plot | positive_chisq_plot)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Simulation operating characteristics",
    subtitle = paste0("1,000 replicates; ", sim_design_label, "; ", sim_pc_label, " matching/adjustment; rejection threshold p < 0.05"),
    caption = "Rates are proportions with p < 0.05. Positive-control median chi-square ratio is descriptive; it has no calibration target.\n\u2020 Null Type I error is inflated; apparent power should be interpreted cautiously.",
    theme = theme(plot.background = element_rect(fill = "white", colour = NA), plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9.5), plot.caption = element_text(size = 8.5, colour = "grey30"))
  ) & theme(legend.position = "bottom")
ggsave(file.path(out_root, "Figure2_simulation_operating_characteristics.pdf"), fig_sim, width = 12.8, height = 8.8, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Figure2_simulation_operating_characteristics.png"), fig_sim, width = 12.8, height = 8.8, units = "in", dpi = 300, bg = "white")

## Main Figure 3: LCT region split into common and rare variants.
split_qq_long <- function(d, phenotype) {
  out <- lapply(names(p_map), function(method) {
    z <- d[is.finite(d$MAF) & !is.na(d$AF_class), , drop = FALSE]
    p_raw <- as.numeric(z[[p_map[[method]]]])
    ok <- is.finite(p_raw) & p_raw > 0 & p_raw <= 1
    z <- z[ok, , drop = FALSE]; p_raw <- p_raw[ok]
    ord <- order(p_raw); p <- p_raw[ord]; z <- z[ord, , drop = FALSE]
    if (!length(p)) return(NULL)
    data.frame(Phenotype = phenotype, AF_class = z$AF_class, method = method, expected = -log10(ppoints(length(p))), observed = -log10(p), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))]); out$AF_class <- factor(out$AF_class, levels = c("Common (MAF >= 1%)", "Rare (MAF < 1%)")); out$method <- factor(out$method, levels = method_levels); out
}
split_qq <- rbind(split_qq_long(lct_tv, "TV negative control"), split_qq_long(lactose, "Lactose positive control"))
write.table(lct_tv[lct_tv$AF_class == "Common (MAF >= 1%)", , drop = FALSE], file.path(out_root, "LCT_TV_common_MAF_ge_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(lct_tv[lct_tv$AF_class == "Rare (MAF < 1%)", , drop = FALSE], file.path(out_root, "LCT_TV_rare_MAF_lt_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(lactose[lactose$AF_class == "Common (MAF >= 1%)", , drop = FALSE], file.path(out_root, "LCT_lactose_common_MAF_ge_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(lactose[lactose$AF_class == "Rare (MAF < 1%)", , drop = FALSE], file.path(out_root, "LCT_lactose_rare_MAF_lt_0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
p_split <- ggplot(split_qq, aes(expected, observed, colour = method, group = method)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") + geom_line(linewidth = .45, alpha = .8) + geom_point(size = .7, alpha = .55) + facet_grid(Phenotype ~ AF_class) + coord_cartesian(xlim = c(0, 8), ylim = c(0, 8)) + scale_colour_manual(values = method_cols, drop = FALSE) + labs(title = "LCT-region analyses split by allele frequency", subtitle = "Legend order progresses toward stronger matching; rare = MAF < 0.01 (1%); common = MAF >= 0.01", x = "Expected -log10(p)", y = "Observed -log10(p)") + paper_theme
ggsave(file.path(out_root, "Figure3_LCT_common_vs_rare_controls.pdf"), p_split, width = 12.2, height = 7.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Figure3_LCT_common_vs_rare_controls.png"), p_split, width = 12.2, height = 7.2, units = "in", dpi = 300, bg = "white")

## Main Table 1: no burden/SKAT/SKAT-O rows; common/rare counts are explicit.
lambda_for <- function(d, class = NULL) {
  if (!is.null(class)) d <- d[d$AF_class == class, , drop = FALSE]
  lambda_from_chisq(d, unname(chisq_map[method_levels]))
}
make_tv_rows <- function(d, class_label) {
  n_cases <- c("Plain logistic" = NA, "PC-adjusted logistic" = NA, "CLR 1:6" = lct_tv_summary$N_CASE_1TO6, "CLR 1:4" = lct_tv_summary$N_CASE_1TO4, "CLR 1:1" = lct_tv_summary$N_CASE_1TO1)
  n_controls <- c("Plain logistic" = NA, "PC-adjusted logistic" = NA, "CLR 1:6" = lct_tv_summary$N_CONTROL_1TO6, "CLR 1:4" = lct_tv_summary$N_CONTROL_1TO4, "CLR 1:1" = lct_tv_summary$N_CONTROL_1TO1)
  max_smd <- c("Plain logistic" = NA, "PC-adjusted logistic" = NA, "CLR 1:6" = lct_tv_summary$MAX_ABS_SMD_POST_1TO6, "CLR 1:4" = lct_tv_summary$MAX_ABS_SMD_POST_1TO4, "CLR 1:1" = lct_tv_summary$MAX_ABS_SMD_POST_1TO1)
  vals <- lambda_for(d); data.frame(Section = "LCT TV negative control", Variant_class = class_label, Method = method_levels, Calibration = "asymptotic chi-square", `lambda_GC or p` = fmt_num(vals), `N variants` = nrow(d), `N cases` = unname(n_cases[method_levels]), `N controls` = unname(n_controls[method_levels]), `max |SMD|` = unname(max_smd[method_levels]), check.names = FALSE)
}
tv_common <- lct_tv[lct_tv$AF_class == "Common (MAF >= 1%)", , drop = FALSE]; tv_rare <- lct_tv[lct_tv$AF_class == "Rare (MAF < 1%)", , drop = FALSE]
table_tv <- rbind(make_tv_rows(lct_tv, "All"), make_tv_rows(tv_common, "Common"), make_tv_rows(tv_rare, "Rare"))
lead_rows <- data.frame(Section = "LCT lactose positive control", Variant_class = "rs4988235_A (common)", Method = method_levels, Calibration = "association p-value", `lambda_GC or p` = fmt_p(sentinel_stats$P), `N variants` = nrow(lactose), `N cases` = NA, `N controls` = NA, `max |SMD|` = NA, check.names = FALSE)
chr2_cases <- c("Plain logistic" = NA, "PC-adjusted logistic" = NA, "CLR 1:6" = 15000, "CLR 1:4" = 15000, "CLR 1:1" = 15000)
chr2_controls <- c("Plain logistic" = NA, "PC-adjusted logistic" = NA, "CLR 1:6" = 90000, "CLR 1:4" = 60000, "CLR 1:1" = 15000)
chr2_rows <- data.frame(Section = "chr2 random rare TV negative control", Variant_class = "Rare", Method = method_levels, Calibration = "observed chi-square statistic", `lambda_GC or p` = fmt_num(lambda_from_chisq(chr2_asym, chisq_cols_chr2)), `N variants` = nrow(chr2_asym), `N cases` = unname(chr2_cases[method_levels]), `N controls` = unname(chr2_controls[method_levels]), `max |SMD|` = NA, check.names = FALSE)
real_table <- rbind(table_tv, chr2_rows, lead_rows)
write.table(real_table, file.path(out_root, "Table1_primary_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_table_pdf(real_table, file.path(out_root, "Table1_primary_results.pdf"), "Table 1. Primary controls with LCT common/rare split", landscape = TRUE, font_size = 5.8)

## Supplementary S3: permutation versus asymptotic p-values and GIF_QQperm.
build_gif_data <- function(region, perm_path, wald_path) {
  pe <- read_tab(perm_path); wa <- read_tab(wald_path); m <- merge(pe, wa, by = "SNP", suffixes = c("_PERM", "_WALD"), sort = FALSE)
  out <- list(); summary <- list()
  for (ratio in c(6L, 4L, 1L)) {
    # These columns occur in only one input table, so merge() leaves their
    # names unchanged (the suffixes apply only to duplicated names).
    po <- as.numeric(m[[paste0("P_1TO", ratio)]])
    px <- as.numeric(m[[paste0("P_PERM_1TO", ratio)]])
    keep <- is.finite(po) & is.finite(px) & po > 0 & po <= 1 & px > 0 & px <= 1
    if (sum(keep) < 3L) next
    x <- sort(qchisq(px[keep], 1, lower.tail = FALSE)); y <- sort(qchisq(po[keep], 1, lower.tail = FALSE)); gif <- gif_qqperm(po[keep], px[keep])
    out[[length(out) + 1L]] <- data.frame(REGION = region, RATIO = paste0("1:", ratio), expected = -log10(pmax(1e-300, pchisq(x, 1, lower.tail = FALSE))), observed = -log10(pmax(1e-300, pchisq(y, 1, lower.tail = FALSE))), GIF_QQperm = gif, stringsAsFactors = FALSE)
    summary[[length(summary) + 1L]] <- data.frame(REGION = region, RATIO = paste0("1:", ratio), N_VARIANTS = sum(keep), GIF_QQperm = gif, lambda_GC_observed = lambda_gc_p(po[keep]), lambda_GC_permutation = lambda_gc_p(px[keep]), stringsAsFactors = FALSE)
  }
  list(plot = do.call(rbind, out), summary = do.call(rbind, summary))
}
gif_lct <- build_gif_data("LCT TV rare", file.path(real_root, "LCT_TV_15000_common_pool_rare_permutation", "comparison_rare_permutation_1to1_1to4_1to6.tsv"), file.path(real_root, "LCT_TV_15000_common_pool_rare", "comparison_rare_1to1_1to4_1to6.tsv"))
gif_chr2 <- build_gif_data("chr2 random100 TV rare", file.path(real_root, "CHR2_random100_TV_15000_rare_permutation", "comparison_rare_permutation_1to1_1to4_1to6.tsv"), file.path(real_root, "CHR2_random100_TV_15000_rare", "comparison_rare_1to1_1to4_1to6.tsv"))
gif_plot_data <- rbind(gif_lct$plot, gif_chr2$plot); gif_summary <- rbind(gif_lct$summary, gif_chr2$summary)
gif_plot_data$REGION <- factor(gif_plot_data$REGION, levels = c("LCT TV rare", "chr2 random100 TV rare")); gif_plot_data$RATIO <- factor(gif_plot_data$RATIO, levels = ratio_levels)
gif_labels <- unique(gif_plot_data[c("REGION", "RATIO", "GIF_QQperm")]); gif_labels$x <- 0.15; gif_labels$y <- 7.7; gif_labels$label <- sprintf("GIF_QQperm = %.2f", gif_labels$GIF_QQperm)
gif_fig <- ggplot(gif_plot_data, aes(expected, observed)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") + geom_point(size = 1.0, alpha = .65, colour = "#0072B2") + facet_grid(REGION ~ RATIO) + geom_text(data = gif_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.1) + coord_cartesian(xlim = c(0, 8), ylim = c(0, 8)) + labs(title = "Permutation versus asymptotic matched rare-variant calibration", subtitle = "Expected axis: within-stratum permutation p-values; observed axis: asymptotic CLR score/Wald p-values", x = "Permutation expected -log10(p)", y = "Asymptotic observed -log10(p)") + paper_theme
write.table(gif_summary, file.path(out_root, "Supplementary_TableS3_permutation_GIF.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gif_pdf <- file.path(out_root, "Supplementary_FigureS3_permutation_vs_no_permutation_GIF.pdf")
pdf(gif_pdf, width = 13, height = 8.5, bg = "white", useDingbats = FALSE); print(gif_fig)
grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
method_text <- c(
  "Permutation-calibrated GIF method (pre-specified for this package)",
  "1. PC-only matched sets were retained: one case and exactly 1, 4, or 6 controls per stratum; sex was not used.",
  "2. For each rare variant, labels were permuted within matched strata (B = 10,000). The empirical p-value was p_perm = (1 + sum_b I[T_b >= T_obs])/(B + 1).",
  "3. To compare calibration with the asymptotic matched test, p_obs is the CLR score/Wald p-value and p_exp is the matched permutation p-value for the same variants.",
  "4. GIF_QQperm is the through-origin slope: beta_hat = argmin_beta sum_i (X_i - beta*Y_i)^2, X_i = chi-square_1^{-1}(1-p_obs,i), Y_i = chi-square_1^{-1}(1-p_exp,i), after sorting each vector.",
  "   This is the chi-square-scale empirical-null QQ regression used by the QQperm implementation; lambda_GC = median(chi-square_1 statistics)/0.454936 is reported separately.",
  "Interpretation: GIF_QQperm near 1 indicates agreement between asymptotic and permutation calibration. It is not the legacy 0.5/median(p_perm) diagnostic, which is not used for the final GIF.",
  "References:",
  "Mitchell et al. Nature Communications (2025), DOI: https://doi.org/10.1038/s41467-025-56944-1",
  "Nature article: https://www.nature.com/articles/s41467-025-56944-1 (n-1 case-control permutation expected p-values).",
  "Wang et al. Nature (2021), DOI: https://doi.org/10.1038/s41586-021-03855-y",
  "Nature article: https://www.nature.com/articles/s41586-021-03855-y (label permutation preserves genotype structure).",
  "QQperm documentation: https://www.rdocumentation.org/packages/QQperm/versions/1.0.1/topics/estlambda2",
  "QQperm source: https://github.com/cran/QQperm/blob/master/R/QQ.R (empirical-null QQ regression on chi-square scale)."
)
for (i in seq_along(method_text)) {
  grid.text(method_text[i], x = unit(.04, "npc"), y = unit(.90 - (i - 1) * .06, "npc"), just = c("left", "top"), gp = gpar(fontfamily = "sans", fontsize = 8.8))
}
dev.off()
writeLines(method_text, file.path(out_root, "Methods_permutation_GIF.txt"))
file.copy(gif_pdf, file.path(out_root, "Supplementary_FigureS3_permutation_vs_no_permutation_GIF_methods.pdf"), overwrite = TRUE)

## Supplementary S4: set-test sensitivity (burden/SKAT/SKAT-O removed from main).
set_res <- read_tab(file.path(real_root, "LCT_region_lactose_region_set_tests", "lactose_region_burden_skat_skato_results.tsv"))
set_res <- set_res[set_res$STATUS == "ok" & set_res$VARIANT_SET %in% c("ALL_RARE", "RS_ONLY") & set_res$METHOD %in% c("BURDEN", "SKAT", "SKAT_O"), , drop = FALSE]
set_res$SET <- ifelse(set_res$VARIANT_SET == "ALL_RARE", "AF-table rare (18 variants)", "rs-only rare (14 variants)")
set_res$MODEL <- ifelse(set_res$MODEL == "PC_ADJUSTED", "PC-adjusted", "Plain"); set_res$TEST <- c(BURDEN = "Burden", SKAT = "SKAT", SKAT_O = "SKAT-O")[set_res$METHOD]
p_set <- ggplot(set_res, aes(TEST, -log10(P_VALUE), colour = MODEL, shape = MODEL, group = MODEL)) + geom_hline(yintercept = -log10(.05), linetype = "dashed", colour = "grey35") + geom_point(position = position_dodge(width = .45), size = 2.4) + facet_wrap(~SET, nrow = 1) + scale_colour_manual(values = c("Plain" = "#999999", "PC-adjusted" = "#0072B2")) + scale_shape_manual(values = c("Plain" = 15, "PC-adjusted" = 16)) + labs(title = "LCT/lactose rare-variant set-test sensitivity", subtitle = "Burden, SKAT, and SKAT-O are exploratory supplementary analyses; dashed line = p=0.05", x = NULL, y = expression(-log[10](p))) + paper_theme
perm_set <- read_tab(file.path(real_root, "LCT_region_lactose_region_set_tests", "matched_burden_permutation", "matched_burden_permutation_results.tsv")); perm_set$SET <- ifelse(perm_set$SET == "MAF_0.01_ALL", "AF-table rare", "rs-only rare")
perm_set <- perm_set[order(match(as.numeric(perm_set$RATIO), c(6, 4, 1)), perm_set$SET), , drop = FALSE]
p_perm_set <- ggplot(perm_set, aes(factor(paste0("1:", RATIO), levels = ratio_levels), -log10(P_PERM), colour = SET, shape = SET)) + geom_hline(yintercept = -log10(.05), linetype = "dashed", colour = "grey35") + geom_point(size = 2.5) + labs(title = "Matched burden permutation sensitivity", subtitle = "10,000 within-stratum permutations", x = "Controls per case\nIncreasing matching stringency ->", y = expression(-log[10](p))) + paper_theme
fig_s4 <- p_set | p_perm_set
ggsave(file.path(out_root, "Supplementary_FigureS4_burden_SKAT_SKATO.pdf"), fig_s4, width = 12.2, height = 5.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Supplementary_FigureS4_burden_SKAT_SKATO.png"), fig_s4, width = 12.2, height = 5.2, units = "in", dpi = 300, bg = "white")
write.table(set_res, file.path(out_root, "Supplementary_TableS4_burden_SKAT_SKATO.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(perm_set, file.path(out_root, "Supplementary_TableS4_matched_burden_permutation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Supplementary S1/S2/S5 and tables.
sim_rep <- read.csv(file.path(sim_root, "replicate_results.csv"), check.names = FALSE)
sim_rep$ratio_label <- factor(as.character(sim_rep$ratio_label), levels = ratio_levels)
sim_rep$scenario <- factor(as.character(sim_rep$scenario), levels = c("Negative control", "Positive control"))
sim_rep$method <- as.character(sim_rep$method)

sim_qq_input <- sim_rep[sim_rep$method == "Matching + CLR" | (sim_rep$method != "Matching + CLR" & sim_rep$ratio == 1), , drop = FALSE]
sim_qq_input$panel <- ifelse(sim_qq_input$method == "Matching + CLR", paste0("Matching + CLR ", as.character(sim_qq_input$ratio_label)), sim_qq_input$method)
sim_qq_panels <- c("Oracle GWAS", "Standard GWAS", "PC-adjusted GWAS", paste0("Matching + CLR ", ratio_levels))
sim_qq_input$panel <- factor(sim_qq_input$panel, levels = sim_qq_panels)
sim_qq_groups <- split(sim_qq_input, interaction(sim_qq_input$scenario, sim_qq_input$panel, drop = TRUE))
sim_qq <- do.call(rbind, lapply(sim_qq_groups, function(d) {
  p <- as.numeric(d$p_value); p <- sort(p[is.finite(p) & p >= 0 & p <= 1]); p <- pmax(p, .Machine$double.xmin)
  data.frame(scenario = as.character(d$scenario[[1L]]), panel = as.character(d$panel[[1L]]), expected = -log10(ppoints(length(p))), observed = -log10(p), stringsAsFactors = FALSE)
}))
sim_qq$scenario <- factor(sim_qq$scenario, levels = c("Negative control", "Positive control")); sim_qq$panel <- factor(sim_qq$panel, levels = sim_qq_panels)
p_sim_qq <- ggplot(sim_qq, aes(expected, observed, colour = scenario, shape = scenario)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey35", linewidth = .5) +
  geom_point(size = .95, alpha = .72) + facet_grid(scenario ~ panel, scales = "free_y") +
  coord_cartesian(xlim = c(0, max(sim_qq$expected, na.rm = TRUE) * 1.03)) +
  scale_colour_manual(values = c("Negative control" = "#0072B2", "Positive control" = "#D55E00")) +
  scale_shape_manual(values = c("Negative control" = 16, "Positive control" = 17)) +
  labs(title = "QQ plots across simulation methods", subtitle = paste0(sim_pc_label, " matching/adjustment; facets progress from 1:6 to 1:1 (increasing stringency)"), x = expression(Expected~~-log[10](italic(p))), y = expression(Observed~~-log[10](italic(p)))) +
  paper_theme + theme(legend.position = "none", strip.text = element_text(face = "bold", size = 8))
ggsave(file.path(out_root, "Supplementary_FigureS1_simulation_QQ.pdf"), p_sim_qq, width = 15.5, height = 6.8, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Supplementary_FigureS1_simulation_QQ.png"), p_sim_qq, width = 15.5, height = 6.8, units = "in", dpi = 300, bg = "white")

sim_balance <- read.csv(file.path(sim_root, "balance_summary.csv"), check.names = FALSE)
sim_balance <- sim_balance[sim_balance$stage == "After matching", , drop = FALSE]
sim_balance$ratio_label <- factor(as.character(sim_balance$ratio_label), levels = ratio_levels)
sim_balance$scenario <- factor(as.character(sim_balance$scenario), levels = c("Negative control", "Positive control"))
sim_balance$variable <- factor(sim_balance$variable, levels = c("ancestry", sim_pc_names), labels = c("True ancestry", sim_pc_names))
sim_balance_labels <- c("True ancestry", sim_pc_names)
sim_balance_colours <- setNames(c("#CC79A7", "#0072B2", "#009E73", "#E69F00", "#D55E00", "#56B4E9")[seq_along(sim_balance_labels)], sim_balance_labels)
sim_balance_shapes <- setNames(c(15, 17, 16, 18, 8, 3)[seq_along(sim_balance_labels)], sim_balance_labels)
p_sim_balance <- ggplot(sim_balance, aes(ratio_label, median_abs_smd, ymin = q05, ymax = q95, colour = variable, shape = variable, group = variable)) +
  geom_hline(yintercept = .1, linetype = "dashed", colour = "grey35") + geom_line(linewidth = .65) +
  geom_errorbar(width = .12, linewidth = .5) + geom_point(size = 2.4) + facet_wrap(~scenario, nrow = 1) +
  scale_colour_manual(values = sim_balance_colours) +
  scale_shape_manual(values = sim_balance_shapes) +
  labs(title = "Balance after PC-based matching", subtitle = paste0(sim_pc_label, "; points: median absolute SMD; bars: 5th-95th percentiles"), x = "Controls per case\nIncreasing matching stringency ->", y = "Absolute standardized mean difference") + paper_theme
ggsave(file.path(out_root, "Supplementary_FigureS2_simulation_balance.pdf"), p_sim_balance, width = 10, height = 4.8, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Supplementary_FigureS2_simulation_balance.png"), p_sim_balance, width = 10, height = 4.8, units = "in", dpi = 300, bg = "white")

clr_rep <- sim_rep[sim_rep$method == "Matching + CLR", , drop = FALSE]
effect_groups <- split(clr_rep, interaction(clr_rep$scenario, clr_rep$ratio_label, drop = TRUE))
effect_sensitivity <- do.call(rbind, lapply(effect_groups, function(d) {
  b <- as.numeric(d$beta); b <- b[is.finite(b)]; se_mean <- if (length(b) > 1L) sd(b) / sqrt(length(b)) else NA_real_
  data.frame(scenario = as.character(d$scenario[[1L]]), ratio_label = as.character(d$ratio_label[[1L]]), mean_beta = mean(b), lower = mean(b) - 1.96 * se_mean, upper = mean(b) + 1.96 * se_mean)
}))
effect_sensitivity$scenario <- factor(effect_sensitivity$scenario, levels = c("Negative control", "Positive control")); effect_sensitivity$ratio_label <- factor(effect_sensitivity$ratio_label, levels = ratio_levels)
beta_positive <- as.numeric(sim_settings$value[match("beta_positive", sim_settings$parameter)])
effect_refs <- data.frame(scenario = factor(c("Negative control", "Positive control"), levels = levels(effect_sensitivity$scenario)), truth = c(0, beta_positive))
rate_sensitivity <- sim[sim$method == "Matching + CLR", , drop = FALSE]
rate_sensitivity$metric <- factor(ifelse(rate_sensitivity$scenario == "Positive control", "Power", "Type I error"), levels = c("Power", "Type I error"))
p_ratio_rate <- ggplot(rate_sensitivity, aes(ratio_label, rejection_rate, ymin = ci_lower, ymax = ci_upper, group = 1)) +
  geom_hline(data = data.frame(metric = factor("Type I error", levels = levels(rate_sensitivity$metric)), ref = .05), aes(yintercept = ref), linetype = "dashed", colour = "grey35") +
  geom_line(colour = "#009E73", linewidth = .65) + geom_errorbar(colour = "#009E73", width = .12) + geom_point(colour = "#009E73", size = 2.4) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") + labs(title = "Matching + CLR: ratio sensitivity", subtitle = "95% Wilson intervals; dashed line = 5% in the null panel", x = "Controls per case\nIncreasing matching stringency ->", y = "Rejection rate") + paper_theme
p_ratio_effect <- ggplot(effect_sensitivity, aes(ratio_label, mean_beta, ymin = lower, ymax = upper, group = 1)) +
  geom_hline(data = effect_refs, aes(yintercept = truth), linetype = "dashed", colour = "grey35") +
  geom_line(colour = "#009E73", linewidth = .65) + geom_errorbar(colour = "#009E73", width = .12) + geom_point(colour = "#009E73", size = 2.4) +
  facet_wrap(~scenario, nrow = 1, scales = "free_y") + labs(title = "Matching + CLR effect estimates", subtitle = "Dashed lines are the data-generating effects; points are mean beta", x = "Controls per case\nIncreasing matching stringency ->", y = "Mean estimated log(OR)") + paper_theme
info_mean <- aggregate(cbind(n_discordant, matched_mac) ~ scenario + ratio_label, data = clr_rep, FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
info_long <- rbind(data.frame(scenario = info_mean$scenario, ratio_label = info_mean$ratio_label, metric = "Discordant strata", value = info_mean$n_discordant), data.frame(scenario = info_mean$scenario, ratio_label = info_mean$ratio_label, metric = "Matched MAC", value = info_mean$matched_mac))
info_long$scenario <- factor(as.character(info_long$scenario), levels = c("Negative control", "Positive control")); info_long$ratio_label <- factor(as.character(info_long$ratio_label), levels = ratio_levels); info_long$metric <- factor(info_long$metric, levels = c("Discordant strata", "Matched MAC"))
p_ratio_info <- ggplot(info_long, aes(ratio_label, value, group = 1)) + geom_line(colour = "#CC79A7", linewidth = .65) + geom_point(colour = "#CC79A7", size = 2.4) + facet_grid(metric ~ scenario, scales = "free_y") + labs(title = "Information retained after matching", subtitle = "Higher matched MAC and more discordant strata improve rare-variant information", x = "Controls per case\nIncreasing matching stringency ->", y = "Mean value") + paper_theme
fig_ratio_sensitivity <- (p_ratio_rate / p_ratio_effect / p_ratio_info) + plot_layout(heights = c(1, 1, 1.35))
ggsave(file.path(out_root, "Supplementary_FigureS2b_ratio_sensitivity_expanded.pdf"), fig_ratio_sensitivity, width = 12.2, height = 11.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_root, "Supplementary_FigureS2b_ratio_sensitivity_expanded.png"), fig_ratio_sensitivity, width = 12.2, height = 11.2, units = "in", dpi = 300, bg = "white")

sim <- sim[order(factor(sim$scenario, levels = c("Negative control", "Positive control")), sim$ratio_label, sim$method), , drop = FALSE]
write.table(sim, file.path(out_root, "Supplementary_TableS1_simulation_full_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
sim_pdf <- sim[, c("scenario", "ratio_label", "method", "rejection_rate", "metric_interpretation", "median_chisq_ratio", "coverage_95", "n_valid", "n_failed", "mean_matched_cases"), drop = FALSE]
write_table_pdf(sim_pdf, file.path(out_root, "Supplementary_TableS1_simulation_full_metrics.pdf"), "Supplementary Table S1. Full simulation operating characteristics", landscape = TRUE, font_size = 6.3)
sens <- read_tab(file.path(real_root, "LCT_region_lactose_region_set_tests", "lactose_sensitivity_set_tests.tsv")); write.table(sens, file.path(out_root, "Supplementary_TableS2_real_set_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
sentinel_tsv <- sentinel_stats[, c("Method", "SNP", "Effect_allele", "N", "N_case", "N_control", "N_genotype_missing", "Effect_allele_count", "EAF", "MAC", "Beta", "SE", "OR", "OR_lower", "OR_upper", "P")]
sentinel_tsv$Method <- as.character(sentinel_tsv$Method)
write.table(sentinel_tsv, file.path(out_root, "Supplementary_TableS5_rs4988235_positive_control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
sentinel_pdf <- data.frame(
  Method = as.character(sentinel_stats$Method),
  `Effect allele` = sentinel_stats$Effect_allele,
  N = format(round(sentinel_stats$N), big.mark = ",", scientific = FALSE, trim = TRUE),
  Cases = format(round(sentinel_stats$N_case), big.mark = ",", scientific = FALSE, trim = TRUE),
  Controls = format(round(sentinel_stats$N_control), big.mark = ",", scientific = FALSE, trim = TRUE),
  `A-allele count` = format(round(sentinel_stats$Effect_allele_count), big.mark = ",", scientific = FALSE, trim = TRUE),
  EAF = fmt_num(sentinel_stats$EAF, 3),
  MAC = format(round(sentinel_stats$MAC), big.mark = ",", scientific = FALSE, trim = TRUE),
  `Beta (SE)` = paste0(fmt_num(sentinel_stats$Beta, 3), " (", fmt_num(sentinel_stats$SE, 3), ")"),
  `OR (95% CI)` = paste0(fmt_num(sentinel_stats$OR, 3), " (", fmt_num(sentinel_stats$OR_lower, 3), "-", fmt_num(sentinel_stats$OR_upper, 3), ")"),
  P = fmt_p(sentinel_stats$P),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table_pdf(sentinel_pdf, file.path(out_root, "Supplementary_TableS5_rs4988235_positive_control.pdf"), "Supplementary Table S5. rs4988235-A European lactase-persistence sentinel", landscape = TRUE, font_size = 9)
writeLines(c(
  "Effect allele A is the allele encoded by the PLINK dosage column rs4988235_A; beta and odds ratios are per A allele.",
  "N, case count, and control count are variant-complete analyzed samples. Effect-allele count is the summed A-allele dosage; MAC is the minor allele count.",
  "Full-sample models are logistic regressions; CLR models are exact conditional logistic regressions in PC-only matched strata."
), file.path(out_root, "Supplementary_TableS5_rs4988235_positive_control_notes.txt"))
file.copy(file.path(real_root, "LCT_region_lactose_region_set_tests", "lactose_variant_firth_forest.pdf"), file.path(out_root, "Supplementary_FigureS5_variant_firth_forest.pdf"), overwrite = TRUE)
file.copy(file.path(real_root, "LCT_region_lactose_region_set_tests", "lactose_variant_firth_forest.png"), file.path(out_root, "Supplementary_FigureS5_variant_firth_forest.png"), overwrite = TRUE)
file.copy(file.path(real_root, "LCT_region_lactose_region_set_tests", "lactose_variant_carrier_firth_pc_adjusted.tsv"), file.path(out_root, "Supplementary_FigureS5_variant_firth_data.tsv"), overwrite = TRUE)

writeLines(c(
  "RECOMMENDED MANUSCRIPT PACKAGE - common/rare LCT split and permutation GIF",
  paste0("MAIN: Figure 1 real-data negative/positive controls; Figure 2 simulation (", sim_design_label, "); Figure 3 LCT TV and lactose split by MAF < 0.01 (1%). Burden/SKAT/SKAT-O are not in the main figures."),
  "SUPPLEMENTARY: S1 simulation QQ; S2 simulation matching balance and S2b ratio-sensitivity display; S3 permutation versus asymptotic rare-variant QQ and GIF formula/citations; S4 burden/SKAT/SKAT-O sensitivity; S5 Firth sensitivity.",
  "Permutation p-values use 10,000 within-stratum label permutations, preserving one case and the exact 1:1, 1:4, or 1:6 control ratio. Matching uses PCs only; sex is excluded.",
  "Figure 1 QQ panels use asymptotic p-values derived from observed chi-square statistics for both LCT and chromosome 2; permutation-based QQ/GIF comparisons are confined to Supplementary Figure S3.",
  "Figure 1 shows the rs4988235-A European lactase-persistence sentinel as per-A-allele odds ratios with 95% confidence intervals and -log10(p) bars with exact p-value labels; genotype counts are in Supplementary Table S5.",
  paste0("Figure 2 uses ", sim_pc_label, " for both matching and PC-adjusted GWAS. Lambda_GC is used only for negative-control calibration. The positive-control median chi-square ratio is shown separately as a descriptive signal-strength metric and is not interpreted against a target of 1."),
  "Across figures and tables, methods are ordered as plain, PC-adjusted, CLR 1:6, CLR 1:4, and CLR 1:1, progressing toward stronger PC-only matching.",
  "Rare definition: strict MAF < 0.01 (1%) from data/plink_LCT.afreq; no variant had exactly MAF=0.01, so strict and inclusive lists coincide here.",
  "The permutation GIF is the chi-square-scale through-origin QQ slope comparing asymptotic matched p-values with empirical within-stratum permutation p-values. Conventional lambda_GC is reported separately.",
  "This package is built by analysis/00_CLR/036_build_paper_package_recommended.R."
), file.path(out_root, "README_paper_package.txt"))
writeLines(c(
  "Figure 1. Real-data controls: LCT TV negative control, chromosome-2 random rare TV negative control, LCT/lactose regional positive control, and per-A-allele odds ratios with 95% confidence intervals plus -log10(p) bars for rs4988235-A. Methods progress from plain and PC-adjusted models through CLR 1:6, 1:4, and 1:1, representing increasing matching stringency.",
  if (sim_is_calibrated) paste0("Figure 2. Real-data-aligned simulation sensitivity: type I error rate and lambda_GC under the negative control, plus statistical power and the median chi-square ratio under the positive control. The positive-control ratio is descriptive signal strength, not calibration, and has no target value of 1. Rates are estimated as the proportion of 1,000 replicates with p < 0.05. N=50,000, MAF_A=0.01, MAF_B=0.001, moderate positive effect OR=1.5, common 2,500-case pool, ", sim_pc_label, "-only Mahalanobis matching without replacement. Standard and PC-adjusted GWAS power estimates are marked because their null type I error is inflated.") else "Figure 2. Simulation operating characteristics: type I error rate and lambda_GC under the negative control, plus statistical power and descriptive median chi-square ratio under the positive control.",
  "Figure 3. LCT-region QQ plots split into common (MAF >= 1%) and rare (MAF < 1%) variants for TV and lactose. Set tests are supplementary.",
  "Supplementary Figure S3. Permutation versus asymptotic matched rare-variant calibration. GIF_QQperm is the through-origin slope on the chi-square scale; formula and citations are on page 2.",
  "Supplementary Table S5. Exact rs4988235-A effect estimates, p-values, analyzed sample counts, effect-allele counts, EAF, and MAC for all five methods.",
  "Mitchell et al. Nat Commun 2025 doi:10.1038/s41467-025-56944-1; Wang et al. Nature 2021 doi:10.1038/s41586-021-03855-y; QQperm documentation/source: rdocumentation.org/packages/QQperm and github.com/cran/QQperm."
), file.path(out_root, "Figure_legends_and_results_draft.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
message("Built recommended package: ", out_root)
