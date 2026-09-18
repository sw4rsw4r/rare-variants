#!/usr/bin/env Rscript

# Manuscript-style SMD plot for 22 UKB assessment centres.
# Reads saved centre-level SMDs and does not rerun matching.

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

input_path <- file.path(
  project_root,
  "results",
  "00_CLR",
  "paper_package",
  "recruitment_centre_balance",
  "recruitment_centre_balance_detailed.tsv"
)
output_dirs <- c(
  file.path(project_root, "results", "00_CLR", "paper_package", "geographic_phenotype_smd"),
  file.path(workspace_root, "cursor_figure_package", "recruitment_region_smd"),
  file.path(workspace_root, "manuscript", "Figures & Tables")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

ratio_labels <- c("1:6", "1:4", "1:1")
plot_data <- read.delim(input_path, check.names = FALSE)
plot_data <- plot_data[
  !is.na(plot_data$CENTRE) &
    plot_data$CENTRE != "" &
    plot_data$Ratio %in% ratio_labels,
  ,
  drop = FALSE
]
if (nrow(plot_data) != 66L) {
  stop("Expected 66 centre-ratio rows, got ", nrow(plot_data))
}

pre_order <- plot_data[plot_data$Ratio == "1:6", c("CENTRE", "PRE_SMD")]
pre_order <- pre_order[order(abs(pre_order$PRE_SMD), decreasing = FALSE), , drop = FALSE]
plot_data$Ratio <- factor(as.character(plot_data$Ratio), levels = ratio_labels)
plot_data$CENTRE <- factor(plot_data$CENTRE, levels = pre_order$CENTRE)

x_max <- max(0.18, max(abs(c(plot_data$PRE_SMD, plot_data$POST_SMD))) * 1.08)

figure <- ggplot(plot_data, aes(y = CENTRE)) +
  annotate(
    "rect",
    xmin = -0.1,
    xmax = 0.1,
    ymin = -Inf,
    ymax = Inf,
    fill = "#EAF4F1",
    alpha = 0.82
  ) +
  geom_vline(
    xintercept = c(-0.1, 0.1),
    linetype = "dashed",
    colour = "#969696",
    linewidth = 0.55
  ) +
  geom_vline(xintercept = 0, colour = "#7A7A7A", linewidth = 0.60) +
  geom_segment(
    aes(x = PRE_SMD, xend = POST_SMD, yend = CENTRE),
    colour = "#C7C7C7",
    linewidth = 0.65
  ) +
  geom_point(aes(x = PRE_SMD, colour = "Before matching"), size = 2.4) +
  geom_point(aes(x = POST_SMD, colour = "After matching"), size = 2.4) +
  facet_wrap(~Ratio, nrow = 1) +
  scale_colour_manual(
    values = c("Before matching" = "#D55E00", "After matching" = "#0072B2"),
    breaks = c("Before matching", "After matching")
  ) +
  scale_x_continuous(
    limits = c(-x_max, x_max),
    breaks = seq(-0.15, 0.15, by = 0.05),
    labels = function(x) sprintf("%.2f", x),
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
    panel.spacing.x = grid::unit(0.12, "in"),
    strip.background = element_rect(
      fill = "#F0F0F0",
      colour = "#A8A8A8",
      linewidth = 0.55
    ),
    strip.text = element_text(face = "bold", size = 11, margin = margin(6, 0, 6, 0)),
    axis.title.x = element_text(size = 11, margin = margin(8, 0, 0, 0)),
    axis.text = element_text(size = 10, colour = "#4B5563"),
    axis.ticks = element_line(colour = "#595959", linewidth = 0.45),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.spacing.x = grid::unit(0.12, "in"),
    legend.margin = margin(0, 0, 6, 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(6, 12, 6, 8)
  )

stems <- c(
  "Recruitment_centre_matching_balance_LCT_TV_clean",
  "Supplementary Figure 2"
)
for (dest in output_dirs) {
  for (stem in stems) {
    if (basename(dest) != "Figures & Tables" && stem == "Supplementary Figure 2") {
      next
    }
    if (basename(dest) == "Figures & Tables" && stem != "Supplementary Figure 2") {
      next
    }
    pdf_path <- file.path(dest, paste0(stem, ".pdf"))
    png_path <- file.path(dest, paste0(stem, ".png"))
    ggsave(
      pdf_path,
      figure,
      width = 12.8,
      height = 9.6,
      units = "in",
      device = cairo_pdf,
      bg = "white"
    )
    ggsave(
      png_path,
      figure,
      width = 12.8,
      height = 9.6,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    message("Saved: ", png_path)
  }
}
