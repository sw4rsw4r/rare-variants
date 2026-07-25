#!/usr/bin/env Rscript

# Build a compact Firth sensitivity supplement from completed outputs.
# The existing CLR/Firth scripts are not modified and no association model is
# refit here; this figure compares the completed PC-adjusted Firth table with
# the completed PC-adjusted full-sample logistic table for overlapping SNPs.

options(stringsAsFactors = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_file) == 1L) normalizePath(file.path(dirname(sub("^--file=", "", script_file)), "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
real_root <- file.path(project_root, "results", "00_CLR")
out_root <- file.path(real_root, "paper_package")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
firth_path <- file.path(real_root, "LCT_region_lactose_region_set_tests", "lactose_variant_carrier_firth_pc_adjusted.tsv")
logistic_path <- file.path(real_root, "LCT_TV_qualification_comparison", "lactose_intolerance", "logistic_pc_adjusted.tsv")
firth <- read_tab(firth_path)
logistic <- read_tab(logistic_path)
keep_logistic <- logistic[, c("SNP", "BETA", "P"), drop = FALSE]
dat <- merge(firth, keep_logistic, by = "SNP", all = FALSE, sort = FALSE)
num_cols <- c("CASE_CARRIER", "CONTROL_CARRIER", "FIRTH_BETA", "FIRTH_P", "BETA", "P")
for (nm in num_cols) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
dat <- dat[is.finite(dat[["FIRTH_BETA"]]) & is.finite(dat[["FIRTH_P"]]) & dat[["FIRTH_P"]] > 0 & dat[["FIRTH_P"]] <= 1 & is.finite(dat[["BETA"]]) & is.finite(dat[["P"]]) & dat[["P"]] > 0 & dat[["P"]] <= 1, , drop = FALSE]
dat$case_carrier_group <- ifelse(dat[["CASE_CARRIER"]] > 0, ">=1 case carrier", "0 case carriers")
dat$neglog10_firth <- -log10(dat[["FIRTH_P"]])
dat$neglog10_logistic <- -log10(dat[["P"]])
dat$stable <- dat$case_carrier_group == ">=1 case carrier"

stable <- dat[dat$stable, , drop = FALSE]
rho_p_all <- suppressWarnings(cor(dat$neglog10_firth, dat$neglog10_logistic, method = "spearman"))
rho_b_all <- suppressWarnings(cor(dat$FIRTH_BETA, dat$BETA, method = "spearman"))
rho_p_stable <- suppressWarnings(cor(stable$neglog10_firth, stable$neglog10_logistic, method = "spearman"))
rho_b_stable <- suppressWarnings(cor(stable$FIRTH_BETA, stable$BETA, method = "spearman"))
sign_stable <- mean(sign(stable$FIRTH_BETA) == sign(stable$BETA))
p_agree_stable <- mean((stable$FIRTH_P < 0.05) == (stable$P < 0.05))
summary_tab <- data.frame(
  comparison = c("All overlapping variants", "Variants with >=1 case carrier"),
  n_variants = c(nrow(dat), nrow(stable)),
  n_zero_case_carrier = c(sum(!dat$stable), 0L),
  spearman_rho_neglog10_p = c(rho_p_all, rho_p_stable),
  spearman_rho_beta = c(rho_b_all, rho_b_stable),
  direction_concordance = c(mean(sign(dat$FIRTH_BETA) == sign(dat$BETA)), sign_stable),
  p_lt_0.05_concordance = c(mean((dat$FIRTH_P < 0.05) == (dat$P < 0.05)), p_agree_stable),
  stringsAsFactors = FALSE
)
write.table(summary_tab, file.path(out_root, "Supplementary_FigureS3_firth_sensitivity_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(dat[, c("SNP", "CASE_CARRIER", "CONTROL_CARRIER", "case_carrier_group", "FIRTH_BETA", "FIRTH_P", "BETA", "P")], file.path(out_root, "Supplementary_FigureS3_firth_sensitivity_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

paper_theme <- theme_bw(base_size = 10.5) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey35"), plot.caption = element_text(size = 8, colour = "grey30", hjust = 0), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA))
point_cols <- c(">=1 case carrier" = "#0072B2", "0 case carriers" = "#BDBDBD")
point_shapes <- c(">=1 case carrier" = 16, "0 case carriers" = 1)
top <- stable[order(stable$FIRTH_P), , drop = FALSE]
top <- head(top, min(3L, nrow(top)))

p_p <- ggplot(dat, aes(x = neglog10_logistic, y = neglog10_firth, colour = case_carrier_group, shape = case_carrier_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.2, alpha = .9) +
  geom_text(data = top, aes(label = SNP), colour = "#0072B2", size = 2.6, nudge_y = .12, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = point_cols) + scale_shape_manual(values = point_shapes) +
  labs(title = "P-value concordance", subtitle = sprintf("Stable subset: Spearman rho = %.2f; %d/%d p<0.05 calls agree", rho_p_stable, round(p_agree_stable * nrow(stable)), nrow(stable)), x = "PC-adjusted logistic -log10(p)", y = "PC-adjusted Firth -log10(p)") + paper_theme

p_beta <- ggplot(dat, aes(x = BETA, y = FIRTH_BETA, colour = case_carrier_group, shape = case_carrier_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.2, alpha = .9) +
  geom_text(data = top, aes(label = SNP), colour = "#0072B2", size = 2.6, nudge_y = .08, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = point_cols) + scale_shape_manual(values = point_shapes) +
  labs(title = "Effect-size concordance", subtitle = sprintf("Stable subset: Spearman rho = %.2f; direction concordance = %.0f%%", rho_b_stable, 100 * sign_stable), x = "PC-adjusted logistic beta", y = "PC-adjusted Firth beta") + paper_theme

fig <- p_p + p_beta + plot_layout(guides = "collect") + plot_annotation(title = "Firth sensitivity analysis for LCT-region rare variants", subtitle = sprintf("PC-adjusted Firth versus PC-adjusted full-sample logistic; %d overlapping variants", nrow(dat)), caption = sprintf("Stable subset: variants with >=1 case carrier (n=%d). Open grey points: zero case carriers (n=%d), flagged for sparse-data separation.", nrow(stable), sum(!dat$stable))) & theme(legend.position = "bottom")

ggsave(file.path(out_root, "Supplementary_FigureS3_firth_sensitivity_comparison.pdf"), fig, width = 12.2, height = 5.8, units = "in", device = "pdf", bg = "white")
ggsave(file.path(out_root, "Supplementary_FigureS3_firth_sensitivity_comparison.png"), fig, width = 12.2, height = 5.8, units = "in", dpi = 300, bg = "white")
message("Built Firth sensitivity comparison from completed outputs: ", out_root)
