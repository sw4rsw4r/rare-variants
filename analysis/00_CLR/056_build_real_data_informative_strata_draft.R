#!/usr/bin/env Rscript

# Draft supplementary figure showing how the number of informative matched
# sets changes with the case-to-control matching ratio in the two real-data
# rare-variant negative-control panels.

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2", "patchwork", "ragg", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
script_path <- normalizePath(
  sub("^--file=", "", script_arg),
  winslash = "/",
  mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
result_root <- file.path(project_root, "results", "00_CLR")
output_root <- file.path(
  result_root,
  "paper_package",
  "informative_strata_ratio_draft"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

analysis_sources <- c(
  "LCT/TV rare variants" = "LCT_TV_15000_common_pool_rare",
  "Chromosome 2/TV rare variants" = "CHR2_random100_TV_15000_rare"
)
ratios <- c("1:6" = 6L, "1:4" = 4L, "1:1" = 1L)

information <- do.call(rbind, lapply(names(analysis_sources), function(panel) {
  source_name <- unname(analysis_sources[[panel]])
  do.call(rbind, lapply(names(ratios), function(ratio_label) {
    ratio <- unname(ratios[[ratio_label]])
    path <- file.path(
      result_root,
      source_name,
      paste0("ratio_1to", ratio),
      "matched_clr_gwas_rare.tsv"
    )
    data <- read_tab(path)
    required <- c("SNP", "N_INFORMATIVE_STRATA", "MAC", "STATUS")
    if (!all(required %in% names(data))) {
      stop("Missing required columns in ", path)
    }
    data.frame(
      ANALYSIS = panel,
      SNP = as.character(data$SNP),
      RATIO = ratio_label,
      N_INFORMATIVE_STRATA = as.numeric(data$N_INFORMATIVE_STRATA),
      MATCHED_MAC = as.numeric(data$MAC),
      STATUS = as.character(data$STATUS),
      stringsAsFactors = FALSE
    )
  }))
}))

information$ANALYSIS <- factor(
  information$ANALYSIS,
  levels = names(analysis_sources)
)
information$RATIO <- factor(information$RATIO, levels = names(ratios))

expected_counts <- c(
  "LCT/TV rare variants" = 14L,
  "Chromosome 2/TV rare variants" = 100L
)
observed_counts <- tapply(information$SNP, information$ANALYSIS, function(x) {
  length(unique(x))
})
if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
  stop("Unexpected variant counts")
}
if (any(!is.finite(information$N_INFORMATIVE_STRATA)) ||
    any(information$N_INFORMATIVE_STRATA < 0)) {
  stop("Invalid informative-strata counts")
}

reference <- information[information$RATIO == "1:6", c(
  "ANALYSIS", "SNP", "N_INFORMATIVE_STRATA"
)]
names(reference)[3L] <- "N_INFORMATIVE_1TO6"
information <- merge(
  information,
  reference,
  by = c("ANALYSIS", "SNP"),
  all.x = TRUE,
  sort = FALSE
)
information$PERCENT_OF_1TO6 <- ifelse(
  information$N_INFORMATIVE_1TO6 > 0,
  100 * information$N_INFORMATIVE_STRATA /
    information$N_INFORMATIVE_1TO6,
  NA_real_
)

summary_table <- do.call(rbind, lapply(
  split(information, interaction(information$ANALYSIS, information$RATIO)),
  function(data) {
    data.frame(
      ANALYSIS = as.character(data$ANALYSIS[1L]),
      RATIO = as.character(data$RATIO[1L]),
      N_VARIANTS = nrow(data),
      N_WITH_INFORMATIVE_STRATA = sum(data$N_INFORMATIVE_STRATA > 0),
      MEDIAN_INFORMATIVE_STRATA = median(data$N_INFORMATIVE_STRATA),
      Q1_INFORMATIVE_STRATA = unname(
        quantile(data$N_INFORMATIVE_STRATA, 0.25)
      ),
      Q3_INFORMATIVE_STRATA = unname(
        quantile(data$N_INFORMATIVE_STRATA, 0.75)
      ),
      MEDIAN_PERCENT_OF_1TO6 = median(
        data$PERCENT_OF_1TO6,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
))
summary_table$RATIO <- factor(summary_table$RATIO, levels = names(ratios))
summary_table$ANALYSIS <- factor(
  summary_table$ANALYSIS,
  levels = names(analysis_sources)
)
summary_table <- summary_table[
  order(summary_table$ANALYSIS, summary_table$RATIO),
  ,
  drop = FALSE
]
write.table(
  information,
  file.path(output_root, "informative_strata_by_variant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  summary_table,
  file.path(output_root, "informative_strata_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

paper_theme <- theme_bw(base_family = "Aptos", base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.4),
    panel.border = element_blank(),
    axis.line = element_line(colour = "#4B5563", linewidth = 0.45),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5, colour = "#374151"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(
      face = "bold",
      size = 11.5,
      margin = margin(b = 8)
    ),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 9, 8, 8)
  )

distribution_layers <- list(
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    fill = "white",
    colour = "#4B5563",
    linewidth = 0.55
  ),
  geom_jitter(
    width = 0.10,
    height = 0,
    shape = 21,
    size = 1.7,
    stroke = 0.45,
    colour = "#0072B2",
    fill = "white",
    alpha = 0.55
  )
)

p_a <- ggplot(
  information,
  aes(RATIO, N_INFORMATIVE_STRATA)
) +
  distribution_layers +
  facet_wrap(~ANALYSIS, nrow = 1, scales = "free_y") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
    labels = scales::label_comma()
  ) +
  labs(
    x = "Matching ratio",
    y = "Informative matched sets per variant"
  ) +
  paper_theme

figure <- p_a &
  theme(text = element_text(family = "Aptos"))

pdf_path <- file.path(
  output_root,
  "Informative_matched_sets_by_ratio_clean_draft.pdf"
)
png_path <- file.path(
  output_root,
  "Informative_matched_sets_by_ratio_clean_draft.png"
)
ggsave(
  pdf_path,
  figure,
  width = 12.2,
  height = 4.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  png_path,
  figure,
  width = 12.2,
  height = 4.4,
  units = "in",
  dpi = 300,
  device = ragg::agg_png,
  bg = "white"
)

message("Saved: ", pdf_path)
message("Saved: ", png_path)
