#!/usr/bin/env Rscript

# Evaluate recruitment-centre balance after matching on PC1-PC5 only. Saved
# matched sets are reused unchanged; no matching or association model is run.

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

data_root <- file.path(project_root, "data")
result_root <- file.path(project_root, "results", "00_CLR")
source_root <- file.path(result_root, "LCT_TV_15000_common_pool")
analysis_root <- file.path(result_root, "paper_package", "recruitment_centre_balance")
figure_root <- file.path(project_root, "output", "pdf", "recruitment_centre_balance")
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

read_header <- function(path) {
  strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
}

read_selected_tsv <- function(path, selected, classes) {
  header <- read_header(path)
  missing_columns <- setdiff(selected, header)
  if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  col_classes <- rep("NULL", length(header))
  names(col_classes) <- header
  col_classes[selected] <- classes
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    colClasses = unname(col_classes),
    check.names = FALSE,
    comment.char = "",
    quote = "",
    na.strings = c("NA", "")
  )
}

write_tab <- function(data, path) {
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

centre_labels <- c(
  "10003" = "Stockport (pilot)",
  "11001" = "Manchester",
  "11002" = "Oxford",
  "11003" = "Cardiff",
  "11004" = "Glasgow",
  "11005" = "Edinburgh",
  "11006" = "Stoke",
  "11007" = "Reading",
  "11008" = "Bury",
  "11009" = "Newcastle",
  "11010" = "Leeds",
  "11011" = "Bristol",
  "11012" = "Barts",
  "11013" = "Nottingham",
  "11014" = "Sheffield",
  "11016" = "Liverpool",
  "11017" = "Middlesbrough",
  "11018" = "Hounslow",
  "11020" = "Croydon",
  "11021" = "Birmingham",
  "11022" = "Swansea",
  "11023" = "Wrexham"
)

ratios <- c(6L, 4L, 1L)
ratio_labels <- paste0("1:", ratios)

covariate_path <- file.path(data_root, "UKBB_covariate.txt")
covariates <- read_selected_tsv(
  covariate_path,
  c("IID", "assessment_center"),
  c("character", "character")
)
if (anyDuplicated(covariates$IID)) stop("Duplicate IID in UKBB_covariate.txt")
if (anyNA(covariates$assessment_center)) stop("Missing assessment_center values in UKBB_covariate.txt")
unknown_codes <- setdiff(unique(covariates$assessment_center), names(centre_labels))
if (length(unknown_codes)) stop("Unknown assessment-centre codes: ", paste(unknown_codes, collapse = ", "))

phenotype <- readRDS(file.path(data_root, "processed", "TV.RDS"))
if (!is.list(phenotype) || !all(c("case", "ctrl") %in% names(phenotype))) {
  stop("TV phenotype file must contain case and ctrl vectors")
}
raw_case_ids <- as.character(phenotype$case)
raw_control_ids <- as.character(phenotype$ctrl)
if (anyDuplicated(raw_case_ids) || anyDuplicated(raw_control_ids)) stop("Duplicate phenotype IID")
if (length(intersect(raw_case_ids, raw_control_ids))) stop("Case/control overlap in TV phenotype")

selected_pool <- read.delim(
  file.path(source_root, "selected_case_pool.tsv"),
  colClasses = "character",
  check.names = FALSE
)
selected_case_ids <- as.character(selected_pool$IID)
if (length(selected_case_ids) != 15000L || anyDuplicated(selected_case_ids)) {
  stop("The saved common case pool is not 15,000 unique IIDs")
}
if (!all(selected_case_ids %in% raw_case_ids)) stop("Selected case pool is inconsistent with TV phenotype")

selected_case_index <- match(selected_case_ids, covariates$IID)
control_index <- match(raw_control_ids, covariates$IID)
case_join_missing <- sum(is.na(selected_case_index))
control_join_missing <- sum(is.na(control_index))
if (case_join_missing) stop("Selected cases missing from covariate file: ", case_join_missing)
pre_case <- covariates[selected_case_index, , drop = FALSE]
pre_control <- covariates[control_index[!is.na(control_index)], , drop = FALSE]
if (nrow(pre_control) == 0L) stop("No TV controls joined to covariate file")

calculate_level_balance <- function(level, post_case, post_control) {
  pre_case_indicator <- pre_case$assessment_center == level
  pre_control_indicator <- pre_control$assessment_center == level
  post_case_indicator <- post_case$assessment_center == level
  post_control_indicator <- post_control$assessment_center == level

  pre_case_p <- mean(pre_case_indicator)
  pre_control_p <- mean(pre_control_indicator)
  post_case_p <- mean(post_case_indicator)
  post_control_p <- mean(post_control_indicator)
  denominator <- sqrt(
    (pre_case_p * (1 - pre_case_p) + pre_control_p * (1 - pre_control_p)) / 2
  )
  if (!is.finite(denominator) || denominator <= 0) {
    stop("Invalid pre-match pooled Bernoulli SD for centre code ", level)
  }

  data.frame(
    CENTRE_CODE = level,
    CENTRE = unname(centre_labels[[level]]),
    PRE_N_CASE = sum(pre_case_indicator),
    PRE_N_CONTROL = sum(pre_control_indicator),
    POST_N_CASE = sum(post_case_indicator),
    POST_N_CONTROL = sum(post_control_indicator),
    PRE_PROP_CASE = pre_case_p,
    PRE_PROP_CONTROL = pre_control_p,
    POST_PROP_CASE = post_case_p,
    POST_PROP_CONTROL = post_control_p,
    PRE_POOLED_SD = denominator,
    PRE_SMD = (pre_case_p - pre_control_p) / denominator,
    POST_SMD = (post_case_p - post_control_p) / denominator,
    stringsAsFactors = FALSE
  )
}

balance_rows <- list()
summary_rows <- list()
provenance_rows <- list()

for (ratio in ratios) {
  ratio_label <- paste0("1:", ratio)
  match_path <- file.path(source_root, "matching", paste0("ratio_1to", ratio), "matched_sets.tsv")
  matched <- read.delim(match_path, colClasses = c(IID = "character"), check.names = FALSE)
  required <- c("IID", "case", "STRATUM")
  if (!all(required %in% names(matched))) stop("Invalid matched-set file for ", ratio_label)
  if (anyDuplicated(matched$IID)) stop("Duplicate matched IID for ", ratio_label)
  stratum_counts <- table(matched$STRATUM, matched$case)
  if (!all(c("0", "1") %in% colnames(stratum_counts)) ||
      any(stratum_counts[, "1"] != 1L) || any(stratum_counts[, "0"] != ratio)) {
    stop("Incomplete matched strata for ", ratio_label)
  }

  matched_index <- match(as.character(matched$IID), covariates$IID)
  if (anyNA(matched_index)) stop("Matched IID missing from covariate file for ", ratio_label)
  matched <- cbind(
    matched[, c("IID", "case", "STRATUM"), drop = FALSE],
    covariates[matched_index, "assessment_center", drop = FALSE]
  )
  post_case <- matched[matched$case == 1L, , drop = FALSE]
  post_control <- matched[matched$case == 0L, , drop = FALSE]
  if (!setequal(post_case$IID, selected_case_ids)) {
    stop("Matched cases differ from the common 15,000-case pool for ", ratio_label)
  }

  balance <- do.call(
    rbind,
    lapply(names(centre_labels), calculate_level_balance, post_case = post_case, post_control = post_control)
  )
  balance$Ratio <- ratio_label
  balance_rows[[ratio_label]] <- balance

  pre_case_props <- balance$PRE_PROP_CASE
  pre_control_props <- balance$PRE_PROP_CONTROL
  post_case_props <- balance$POST_PROP_CASE
  post_control_props <- balance$POST_PROP_CONTROL
  max_pre_index <- which.max(abs(balance$PRE_SMD))
  max_post_index <- which.max(abs(balance$POST_SMD))
  summary_rows[[ratio_label]] <- data.frame(
    Ratio = ratio_label,
    N_CASE_PRE = nrow(pre_case),
    N_CONTROL_PRE = nrow(pre_control),
    N_CASE_POST = nrow(post_case),
    N_CONTROL_POST = nrow(post_control),
    N_CENTRES = nrow(balance),
    MAX_ABS_SMD_PRE = max(abs(balance$PRE_SMD)),
    MAX_ABS_SMD_POST = max(abs(balance$POST_SMD)),
    MAX_PRE_CENTRE = balance$CENTRE[[max_pre_index]],
    MAX_POST_CENTRE = balance$CENTRE[[max_post_index]],
    N_LEVELS_ABOVE_0_1_PRE = sum(abs(balance$PRE_SMD) >= 0.1),
    N_LEVELS_ABOVE_0_1_POST = sum(abs(balance$POST_SMD) >= 0.1),
    TOTAL_VARIATION_PRE = 0.5 * sum(abs(pre_case_props - pre_control_props)),
    TOTAL_VARIATION_POST = 0.5 * sum(abs(post_case_props - post_control_props)),
    stringsAsFactors = FALSE
  )

  provenance_rows[[ratio_label]] <- data.frame(
    Ratio = ratio_label,
    MATCHED_SET_FILE = normalizePath(match_path, winslash = "/", mustWork = TRUE),
    MD5 = unname(tools::md5sum(match_path)),
    N_CASE = nrow(post_case),
    N_CONTROL = nrow(post_control),
    N_STRATA = length(unique(matched$STRATUM)),
    stringsAsFactors = FALSE
  )
}

detailed <- do.call(rbind, balance_rows)
summary_table <- do.call(rbind, summary_rows)
provenance <- do.call(rbind, provenance_rows)

detailed$Ratio <- factor(detailed$Ratio, levels = ratio_labels)
summary_table$Ratio <- factor(summary_table$Ratio, levels = ratio_labels)
provenance$Ratio <- factor(provenance$Ratio, levels = ratio_labels)
summary_table <- summary_table[order(summary_table$Ratio), , drop = FALSE]
provenance <- provenance[order(provenance$Ratio), , drop = FALSE]

pre_order <- detailed[detailed$Ratio == "1:6", c("CENTRE", "PRE_SMD")]
pre_order <- pre_order[order(abs(pre_order$PRE_SMD), decreasing = FALSE), , drop = FALSE]
detailed$CENTRE <- factor(detailed$CENTRE, levels = pre_order$CENTRE)

max_lookup <- setNames(summary_table$MAX_ABS_SMD_POST, as.character(summary_table$Ratio))
facet_levels <- vapply(
  ratio_labels,
  function(ratio) sprintf("%s\nmax post |SMD| = %.4f", ratio, max_lookup[[ratio]]),
  character(1L)
)
detailed$FACET_LABEL <- factor(
  vapply(
    as.character(detailed$Ratio),
    function(ratio) sprintf("%s\nmax post |SMD| = %.4f", ratio, max_lookup[[ratio]]),
    character(1L)
  ),
  levels = facet_levels
)

plot_limit <- max(0.12, max(abs(c(detailed$PRE_SMD, detailed$POST_SMD))) * 1.12)
plot_breaks <- pretty(c(-plot_limit, plot_limit), n = 7)
plot_breaks <- plot_breaks[plot_breaks >= -plot_limit & plot_breaks <= plot_limit]

detailed_figure <- ggplot(detailed, aes(y = CENTRE)) +
  annotate(
    "rect",
    xmin = -0.1,
    xmax = 0.1,
    ymin = -Inf,
    ymax = Inf,
    fill = "#EAF4F1",
    alpha = 0.82
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "#969696", linewidth = 0.55) +
  geom_vline(xintercept = 0, colour = "#7A7A7A", linewidth = 0.60) +
  geom_segment(
    aes(x = PRE_SMD, xend = POST_SMD, yend = CENTRE),
    colour = "#C7C7C7",
    linewidth = 0.65
  ) +
  geom_point(
    aes(x = PRE_SMD, colour = "Before matching", shape = "Before matching"),
    size = 2.8,
    stroke = 1.1
  ) +
  geom_point(
    aes(x = POST_SMD, colour = "After matching", shape = "After matching"),
    size = 2.8,
    stroke = 1.0
  ) +
  facet_wrap(~FACET_LABEL, nrow = 1) +
  scale_colour_manual(
    values = c("Before matching" = "#D55E00", "After matching" = "#0072B2"),
    breaks = c("Before matching", "After matching")
  ) +
  scale_shape_manual(
    values = c("Before matching" = 1, "After matching" = 16),
    breaks = c("Before matching", "After matching")
  ) +
  scale_x_continuous(
    limits = c(-plot_limit, plot_limit),
    breaks = plot_breaks,
    labels = sprintf("%.2f", plot_breaks),
    expand = expansion(mult = 0)
  ) +
  labs(
    title = "PC-only matching balance for recruitment centre in the LCT-region/TV analysis",
    subtitle = "Centre-specific binary indicators; the same 15,000 cases are used for 1:6, 1:4, and 1:1 matching",
    x = "Standardized mean difference",
    y = NULL,
    colour = NULL,
    shape = NULL,
    caption = "Matching used PC1-PC5 only. The shaded interval denotes |SMD| < 0.1; the pre-match pooled Bernoulli SD is used for both stages."
  ) +
  theme_bw(base_size = 11.2) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#6B6B6B", fill = NA, linewidth = 0.55),
    panel.spacing.x = grid::unit(0.08, "in"),
    strip.background = element_rect(fill = "#F0F0F0", colour = "#A8A8A8", linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 10.2, lineheight = 1.05, margin = margin(6, 0, 5, 0)),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 9.4, colour = "#595959"),
    axis.ticks = element_line(colour = "#595959", linewidth = 0.45),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.margin = margin(-2, 0, 7, 0),
    plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(0, 0, 3, 0)),
    plot.subtitle = element_text(size = 11.1, colour = "#6B7280", hjust = 0, margin = margin(0, 0, 8, 0)),
    plot.caption = element_text(size = 9.1, colour = "#6B7280", hjust = 0, margin = margin(5, 0, 0, 0)),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5, 10, 3, 10)
  )

summary_long <- rbind(
  data.frame(Ratio = summary_table$Ratio, Stage = "Before matching", MAX_ABS_SMD = summary_table$MAX_ABS_SMD_PRE),
  data.frame(Ratio = summary_table$Ratio, Stage = "After matching", MAX_ABS_SMD = summary_table$MAX_ABS_SMD_POST)
)
summary_long$Ratio <- factor(summary_long$Ratio, levels = rev(ratio_labels))
summary_long$Stage <- factor(summary_long$Stage, levels = c("Before matching", "After matching"))

summary_figure <- ggplot(summary_long, aes(y = Ratio)) +
  annotate("rect", xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = "#EAF4F1", alpha = 0.82) +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "#969696", linewidth = 0.55) +
  geom_line(aes(x = MAX_ABS_SMD, group = Ratio), colour = "#C7C7C7", linewidth = 0.75) +
  geom_point(
    aes(x = MAX_ABS_SMD, colour = Stage, shape = Stage),
    size = 3.7,
    stroke = 1.2
  ) +
  geom_text(
    aes(x = MAX_ABS_SMD, label = sprintf("%.4f", MAX_ABS_SMD), colour = Stage),
    hjust = -0.28,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Before matching" = "#D55E00", "After matching" = "#0072B2")) +
  scale_shape_manual(values = c("Before matching" = 1, "After matching" = 16)) +
  scale_x_continuous(
    limits = c(0, max(0.115, max(summary_long$MAX_ABS_SMD) * 1.25)),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Maximum recruitment-centre imbalance",
    subtitle = "Maximum absolute SMD across 22 centre-specific indicators",
    x = "Maximum absolute SMD",
    y = "Matching ratio",
    colour = NULL,
    shape = NULL,
    caption = "The shaded interval denotes maximum |SMD| < 0.1."
  ) +
  theme_bw(base_size = 11.5) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#6B6B6B", fill = NA, linewidth = 0.55),
    axis.text = element_text(colour = "#595959"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11.1, colour = "#6B7280"),
    plot.caption = element_text(size = 9.1, colour = "#6B7280", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 15, 5, 12)
  )

detailed_pdf <- file.path(figure_root, "Recruitment_centre_balance_LCT_TV_all_centres.pdf")
detailed_png <- file.path(figure_root, "Recruitment_centre_balance_LCT_TV_all_centres.png")
summary_pdf <- file.path(figure_root, "Recruitment_centre_balance_LCT_TV_summary.pdf")
summary_png <- file.path(figure_root, "Recruitment_centre_balance_LCT_TV_summary.png")
ggsave(detailed_pdf, detailed_figure, width = 13.5, height = 10.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(detailed_png, detailed_figure, width = 13.5, height = 10.2, units = "in", dpi = 300, bg = "white")
ggsave(summary_pdf, summary_figure, width = 8.2, height = 5.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(summary_png, summary_figure, width = 8.2, height = 5.2, units = "in", dpi = 300, bg = "white")

write_tab(detailed, file.path(analysis_root, "recruitment_centre_balance_detailed.tsv"))
write_tab(summary_table, file.path(analysis_root, "recruitment_centre_balance_summary.tsv"))
write_tab(provenance, file.path(analysis_root, "matched_set_provenance.tsv"))

data_quality <- data.frame(
  CHECK = c(
    "covariate_iid_unique",
    "assessment_centre_complete",
    "assessment_centre_allowed_codes",
    "selected_case_pool_join",
    "tv_control_join",
    "same_common_case_pool",
    "exact_matching_ratios",
    "fixed_pre_match_denominator",
    "no_matching_refit"
  ),
  STATUS = c("PASS", "PASS", "PASS", "PASS", if (control_join_missing == 0L) "PASS" else "NOTE", "PASS", "PASS", "PASS", "PASS"),
  DETAIL = c(
    sprintf("%s unique IIDs in UKBB_covariate.txt", format(nrow(covariates), big.mark = ",")),
    "No missing assessment_center values",
    sprintf("All observed values belong to the expected %d baseline-centre codes", length(centre_labels)),
    sprintf("All %s selected TV cases joined", format(nrow(pre_case), big.mark = ",")),
    sprintf("%s of %s TV controls joined; %s did not join", format(nrow(pre_control), big.mark = ","), format(length(raw_control_ids), big.mark = ","), format(control_join_missing, big.mark = ",")),
    "The same saved 15,000 cases are present in every ratio",
    "1:6, 1:4, and 1:1 contain 90,000, 60,000, and 15,000 matched controls",
    "Pre- and post-match SMDs use the same pre-match pooled Bernoulli SD for each centre",
    "Saved matched IIDs were reused; no matching or association model was rerun"
  ),
  stringsAsFactors = FALSE
)
write_tab(data_quality, file.path(analysis_root, "data_quality_report.tsv"))

writeLines(
  c(
    "Recruitment-centre balance after PC1-PC5-only matching in the LCT-region / TV analysis",
    "The original 15,000-case common pool and saved 1:6, 1:4, and 1:1 matched sets were reused unchanged.",
    "Each of the 22 baseline recruitment centres was represented as a separate binary indicator.",
    "For centre k, SMD_k = (p_case,k - p_control,k) / sqrt((p_case,k*(1-p_case,k) + p_control,k*(1-p_control,k))/2).",
    "The denominator is estimated in the pre-match case/control cohort and reused for the post-match SMD.",
    "The summary maximum absolute SMD is the largest absolute centre-specific SMD.",
    "Total variation distance is 0.5 times the sum of absolute case/control centre-proportion differences.",
    "Sex, assessment centre, and genotype were not used for matching; matching used PC1-PC5 only.",
    "This script evaluates balance and plots saved results only."
  ),
  file.path(analysis_root, "METHODS.txt")
)
capture.output(sessionInfo(), file = file.path(analysis_root, "sessionInfo.txt"))

message("Saved figures to: ", figure_root)
print(summary_table)
