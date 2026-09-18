#!/usr/bin/env Rscript

# Combined PC1-PC10 love plot for the 15,000-case LCT-region/TV analysis.
# PC1-PC5 SMDs are the stored matching-balance values (reproduced by script 056).
# PC6-PC10 SMDs are the previously computed unmatched-PC audit.
# Matching is not rerun.

options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing ggplot2")
if (!requireNamespace("systemfonts", quietly = TRUE)) stop("Missing systemfonts")
suppressPackageStartupMessages(library(ggplot2))

aptos_font_dir <- file.path(
  Sys.getenv("LOCALAPPDATA"), "Microsoft", "FontCache", "4", "CloudFonts", "Aptos"
)
aptos_font_files <- list.files(aptos_font_dir, pattern = "[.]ttf$", full.names = TRUE)
if (!length(aptos_font_files)) stop("Aptos font files were not found")
systemfonts::add_fonts(aptos_font_files)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
workspace_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = TRUE)

source_root <- file.path(project_root, "results", "00_CLR", "LCT_TV_15000_common_pool", "matching")
pc6_path <- file.path(project_root, "results", "00_CLR", "paper_package", "pc6_pc10_balance", "PC6_PC10_balance_detailed.tsv")
output_root <- file.path(project_root, "results", "00_CLR", "paper_package", "pc1_pc10_balance")
package_root <- file.path(workspace_root, "cursor_figure_package", "pc1_pc10_balance")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(package_root, recursive = TRUE, showWarnings = FALSE)

ratio_labels <- c("1:6", "1:4", "1:1")
ratio_dirs <- c("ratio_1to6", "ratio_1to4", "ratio_1to1")
matching_pcs <- paste0("PC", 1:5)
unused_pcs <- paste0("PC", 6:10)
all_pcs <- paste0("PC", 1:10)

pc1_rows <- lapply(seq_along(ratio_labels), function(i) {
  part <- read.delim(file.path(source_root, ratio_dirs[[i]], "pc_balance.tsv"), check.names = FALSE)
  if (!identical(as.character(part$VARIABLE), matching_pcs)) {
    stop("Unexpected PC1-PC5 rows for ", ratio_labels[[i]])
  }
  part$Ratio <- ratio_labels[[i]]
  part$GROUP <- "Used in matching"
  part[, c("Ratio", "VARIABLE", "GROUP", "PRE_SMD", "POST_SMD")]
})
pc1 <- do.call(rbind, pc1_rows)

pc6 <- read.delim(pc6_path, check.names = FALSE)
pc6$Ratio <- as.character(pc6$Ratio)
if (!all(pc6$VARIABLE %in% unused_pcs)) stop("PC6-PC10 file contains unexpected variables")
pc6$GROUP <- "Not used in matching"
pc6 <- pc6[, c("Ratio", "VARIABLE", "GROUP", "PRE_SMD", "POST_SMD")]

plot_data <- rbind(pc1, pc6)
plot_data$Ratio <- factor(plot_data$Ratio, levels = ratio_labels)
plot_data$VARIABLE <- factor(plot_data$VARIABLE, levels = all_pcs)
plot_data$GROUP <- factor(plot_data$GROUP, levels = c("Used in matching", "Not used in matching"))
if (nrow(plot_data) != 30L) stop("Expected 10 PCs x 3 ratios")

max_match <- tapply(abs(pc1$POST_SMD), factor(pc1$Ratio, levels = ratio_labels), max)
max_unused <- tapply(abs(pc6$POST_SMD), factor(pc6$Ratio, levels = ratio_labels), max)

plot_data$FACET_LABEL <- factor(
  as.character(plot_data$Ratio),
  levels = ratio_labels
)

ann <- data.frame(
  FACET_LABEL = factor(ratio_labels[[3L]], levels = ratio_labels),
  VARIABLE = factor(c("PC3", "PC8"), levels = all_pcs),
  txt = c("Used in\nmatching", "Not used in\nmatching"),
  stringsAsFactors = FALSE
)

figure <- ggplot(plot_data, aes(y = VARIABLE)) +
  annotate("rect", xmin = -0.1, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = "#EAF4F1", alpha = 0.82) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "#969696", linewidth = 0.55) +
  geom_vline(xintercept = 0, colour = "#7A7A7A", linewidth = 0.60) +
  geom_hline(yintercept = 5.5, colour = "#A8A8A8", linetype = "dotted", linewidth = 0.55) +
  geom_segment(aes(x = PRE_SMD, xend = POST_SMD, yend = VARIABLE), colour = "#C7C7C7", linewidth = 0.65) +
  geom_point(aes(x = PRE_SMD, colour = "Before matching"), size = 3.0) +
  geom_point(aes(x = POST_SMD, colour = "After matching"), size = 3.0) +
  geom_label(
    data = ann,
    aes(x = 0.135, y = VARIABLE, label = txt),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 3.87,
    family = "Aptos",
    fontface = "bold",
    colour = "#374151",
    fill = "white",
    linewidth = 0,
    label.padding = grid::unit(0.10, "lines"),
    label.r = grid::unit(0.08, "lines"),
    lineheight = 1.0
  ) +
  facet_wrap(~FACET_LABEL, nrow = 1) +
  scale_colour_manual(
    values = c("Before matching" = "#D55E00", "After matching" = "#0072B2"),
    breaks = c("Before matching", "After matching")
  ) +
  scale_x_continuous(
    limits = c(-0.17, 0.17),
    breaks = seq(-0.10, 0.10, by = 0.05),
    labels = sprintf("%.2f", seq(-0.10, 0.10, by = 0.05)),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Standardized mean difference (SMD)",
    y = NULL,
    colour = NULL
  ) +
  theme_bw(base_size = 11, base_family = "Aptos") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#6B6B6B", fill = NA, linewidth = 0.55),
    panel.spacing.x = grid::unit(0.14, "in"),
    strip.background = element_rect(fill = "#F0F0F0", colour = "#A8A8A8", linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 11, margin = margin(7, 0, 7, 0)),
    axis.title.x = element_text(size = 11, margin = margin(8, 0, 0, 0)),
    axis.text = element_text(size = 11, colour = "#4B5563"),
    axis.ticks = element_line(colour = "#595959", linewidth = 0.45),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.spacing.x = grid::unit(0.12, "in"),
    legend.margin = margin(0, 0, 8, 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 14, 8, 10)
  )

pdf_path <- file.path(output_root, "PC1_PC10_matching_balance_LCT_TV_clean_draft.pdf")
png_path <- file.path(output_root, "PC1_PC10_matching_balance_LCT_TV_clean_draft.png")
ggsave(pdf_path, figure, width = 12.5, height = 7.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(png_path, figure, width = 12.5, height = 7.2, units = "in", dpi = 300, bg = "white")

summary_tbl <- data.frame(
  Ratio = ratio_labels,
  MAX_ABS_SMD_POST_PC1_PC5 = as.numeric(max_match[ratio_labels]),
  MAX_ABS_SMD_POST_PC6_PC10 = as.numeric(max_unused[ratio_labels]),
  MAX_ABS_SMD_POST_PC1_PC10 = pmax(as.numeric(max_match[ratio_labels]), as.numeric(max_unused[ratio_labels])),
  stringsAsFactors = FALSE
)
write.table(plot_data[, c("Ratio", "VARIABLE", "GROUP", "PRE_SMD", "POST_SMD")],
            file.path(output_root, "PC1_PC10_matching_balance_LCT_TV.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_tbl, file.path(output_root, "PC1_PC10_max_post_match_smd.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
file.copy(file.path(output_root, "PC1_PC10_matching_balance_LCT_TV.tsv"),
          file.path(package_root, "PC1_PC10_matching_balance_LCT_TV.tsv"), overwrite = TRUE)
file.copy(file.path(output_root, "PC1_PC10_max_post_match_smd.tsv"),
          file.path(package_root, "PC1_PC10_max_post_match_smd.tsv"), overwrite = TRUE)

writeLines(
  c(
    "Combined PC1-PC10 matching balance for the 15,000-case LCT-region / TV analysis.",
    "PC1-PC5 are the stored matching-balance SMDs. PC6-PC10 reuse the unmatched-PC audit.",
    "Matching used PC1-PC5 only. Sex and genotype were not used.",
    "SMD uses the pre-match pooled case/control SD for both stages.",
    "Matching was not rerun."
  ),
  file.path(output_root, "METHODS.txt")
)

message("Saved: ", png_path)
print(summary_tbl)
