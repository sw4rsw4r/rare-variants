#!/usr/bin/env Rscript

# Clean permutation versus asymptotic calibration plot for manuscript review.
# Existing permutation and asymptotic results are reused without refitting models.

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

results_root <- file.path(project_root, "results", "00_CLR")
output_root <- file.path(results_root, "paper_package_recommended")
ratio_labels <- c("1:6", "1:4", "1:1")

lambda_perm <- function(p_asymptotic, p_permutation) {
  keep <- (
    is.finite(p_asymptotic) & is.finite(p_permutation) &
      p_asymptotic > 0 & p_asymptotic <= 1 &
      p_permutation > 0 & p_permutation <= 1
  )
  if (sum(keep) < 3L) return(NA_real_)
  observed <- sort(qchisq(p_asymptotic[keep], 1, lower.tail = FALSE))
  expected <- sort(qchisq(p_permutation[keep], 1, lower.tail = FALSE))
  unname(coef(lm(observed ~ 0 + expected))[1L])
}

build_panel <- function(label, permutation_path, asymptotic_path) {
  permutation <- read.delim(permutation_path, check.names = FALSE)
  asymptotic <- read.delim(asymptotic_path, check.names = FALSE)
  data <- merge(
    permutation,
    asymptotic,
    by = "SNP",
    suffixes = c("_PERMUTATION", "_ASYMPTOTIC"),
    sort = FALSE
  )

  plots <- list()
  summaries <- list()
  for (ratio in c(6L, 4L, 1L)) {
    p_asymptotic <- as.numeric(data[[paste0("P_1TO", ratio)]])
    p_permutation <- as.numeric(data[[paste0("P_PERM_1TO", ratio)]])
    keep <- (
      is.finite(p_asymptotic) & is.finite(p_permutation) &
        p_asymptotic > 0 & p_asymptotic <= 1 &
        p_permutation > 0 & p_permutation <= 1
    )
    p_asymptotic <- sort(p_asymptotic[keep])
    p_permutation <- sort(p_permutation[keep])
    ratio_label <- paste0("1:", ratio)

    plots[[length(plots) + 1L]] <- data.frame(
      Panel = label,
      Ratio = ratio_label,
      expected = -log10(p_permutation),
      observed = -log10(p_asymptotic),
      stringsAsFactors = FALSE
    )
    summaries[[length(summaries) + 1L]] <- data.frame(
      Panel = label,
      Ratio = ratio_label,
      lambda_perm = lambda_perm(p_asymptotic, p_permutation),
      N_variants = length(p_asymptotic),
      stringsAsFactors = FALSE
    )
  }

  list(
    plot = do.call(rbind, plots),
    summary = do.call(rbind, summaries)
  )
}

lct <- build_panel(
  "LCT rare variants",
  file.path(
    results_root,
    "LCT_TV_15000_common_pool_rare_permutation",
    "comparison_rare_permutation_1to1_1to4_1to6.tsv"
  ),
  file.path(
    results_root,
    "LCT_TV_15000_common_pool_rare",
    "comparison_rare_1to1_1to4_1to6.tsv"
  )
)

chr2 <- build_panel(
  "Chromosome 2 rare variants",
  file.path(
    results_root,
    "CHR2_random100_TV_15000_rare_permutation",
    "comparison_rare_permutation_1to1_1to4_1to6.tsv"
  ),
  file.path(
    results_root,
    "CHR2_random100_TV_15000_rare",
    "comparison_rare_1to1_1to4_1to6.tsv"
  )
)

plot_data <- rbind(lct$plot, chr2$plot)
summary_data <- rbind(lct$summary, chr2$summary)
plot_data$Panel <- factor(
  plot_data$Panel,
  levels = c("LCT rare variants", "Chromosome 2 rare variants")
)
plot_data$Ratio <- factor(plot_data$Ratio, levels = ratio_labels)
summary_data$Panel <- factor(
  summary_data$Panel,
  levels = levels(plot_data$Panel)
)
summary_data$Ratio <- factor(summary_data$Ratio, levels = ratio_labels)

axis_max <- max(3, ceiling(max(c(plot_data$expected, plot_data$observed)) * 2) / 2)
summary_data$x <- 0.12
summary_data$y <- axis_max - 0.12
summary_data$label <- sprintf("lambda[perm] == %.2f", summary_data$lambda_perm)

figure <- ggplot(plot_data, aes(expected, observed)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "#6B7280",
    linewidth = 0.55
  ) +
  geom_point(size = 1.45, alpha = 0.75, colour = "#0072B2") +
  geom_text(
    data = summary_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    family = "Aptos",
    size = 3.87,
    colour = "#374151"
  ) +
  facet_grid(Panel ~ Ratio) +
  coord_equal(
    xlim = c(0, axis_max),
    ylim = c(0, axis_max),
    expand = FALSE
  ) +
  scale_x_continuous(breaks = seq(0, axis_max, by = 1)) +
  scale_y_continuous(breaks = seq(0, axis_max, by = 1)) +
  labs(
    x = expression("Permutation expected " * -log[10](P)),
    y = expression("Asymptotic observed " * -log[10](P))
  ) +
  theme_bw(base_size = 11, base_family = "Aptos") +
  theme(
    panel.grid = element_line(colour = "#E5E7EB", linewidth = 0.45),
    panel.border = element_rect(colour = "#6B6B6B", fill = NA, linewidth = 0.55),
    panel.spacing = grid::unit(0.13, "in"),
    strip.background = element_rect(
      fill = "#F0F0F0",
      colour = "#A8A8A8",
      linewidth = 0.55
    ),
    strip.text = element_text(face = "bold", size = 11, margin = margin(6, 6, 6, 6)),
    axis.title = element_text(size = 11),
    axis.title.x = element_text(margin = margin(8, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 8, 0, 0)),
    axis.text = element_text(size = 11, colour = "#4B5563"),
    axis.ticks = element_line(colour = "#595959", linewidth = 0.45),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 16, 8, 10)
  )

pdf_path <- file.path(
  output_root,
  "Permutation_vs_asymptotic_lambda_clean_draft.pdf"
)
png_path <- file.path(
  output_root,
  "Permutation_vs_asymptotic_lambda_clean_draft.png"
)

ggsave(
  pdf_path,
  figure,
  width = 12.8,
  height = 7.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  png_path,
  figure,
  width = 12.8,
  height = 7.4,
  units = "in",
  dpi = 300,
  bg = "white"
)
write.table(
  summary_data[, c("Panel", "Ratio", "lambda_perm", "N_variants")],
  file.path(output_root, "Permutation_vs_asymptotic_lambda_clean_draft.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Saved: ", png_path)
print(summary_data[, c("Panel", "Ratio", "lambda_perm", "N_variants")])
