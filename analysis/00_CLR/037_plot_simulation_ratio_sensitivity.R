#!/usr/bin/env Rscript

# Ratio-sensitivity display for the real-data-aligned simulation.  This is a
# visualization of completed simulation output; it does not refit models.
options(stringsAsFactors = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_file) == 1L) normalizePath(file.path(dirname(sub("^--file=", "", script_file)), "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
sim_root <- Sys.getenv("SIM_ROOT_OVERRIDE", file.path(dirname(project_root), "R simulation", "calibrated_simulation_results"))
out_root <- Sys.getenv("SIM_RATIO_OUT_ROOT", sim_root)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

summary <- read.csv(file.path(sim_root, "summary_metrics.csv"), check.names = FALSE)
raw <- read.csv(file.path(sim_root, "replicate_results.csv"), check.names = FALSE)
summary$ratio_label <- factor(as.character(summary$ratio_label), levels = c("1:1", "1:4", "1:6"))
raw$ratio_label <- factor(paste0("1:", raw$ratio), levels = c("1:1", "1:4", "1:6"))

matching <- summary[summary$method == "Matching + CLR", , drop = FALSE]
matching$metric <- ifelse(matching$scenario == "Negative control", "Type-I error", "Power")
rate_plot <- ggplot(matching, aes(ratio_label, rejection_rate, group = 1)) +
  geom_hline(data = data.frame(metric = "Type-I error", y = .05), aes(yintercept = y), linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = .08, colour = "#009E73") +
  geom_line(colour = "#009E73", linewidth = .7) + geom_point(colour = "#009E73", size = 2.4) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(title = "Matching + CLR: ratio sensitivity", subtitle = "Error bars are 95% Wilson intervals; dashed line = 5% in the null panel", x = "Controls per case", y = "Rejection rate") +
  theme_bw(base_size = 10.5) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey35"))

fit <- raw[raw$method == "Matching + CLR" & is.finite(raw$beta), , drop = FALSE]
beta <- do.call(rbind, lapply(split(fit, list(fit$scenario, fit$ratio_label), drop = TRUE), function(d) {
  data.frame(scenario = as.character(d$scenario[[1L]]), ratio_label = d$ratio_label[[1L]], mean_beta = mean(d$beta), se_beta = sd(d$beta) / sqrt(nrow(d)))
}))
beta$truth <- ifelse(beta$scenario == "Negative control", 0, log(1.5))
beta_plot <- ggplot(beta, aes(ratio_label, mean_beta, group = 1)) +
  geom_hline(aes(yintercept = truth), linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = mean_beta - 1.96 * se_beta, ymax = mean_beta + 1.96 * se_beta), width = .08, colour = "#009E73") +
  geom_line(colour = "#009E73", linewidth = .7) + geom_point(colour = "#009E73", size = 2.4) +
  facet_wrap(~scenario, scales = "free_y", nrow = 1) +
  labs(title = "Matching + CLR effect estimates", subtitle = "Dashed lines are the data-generating effects; points are mean beta", x = "Controls per case", y = "Mean estimated log(OR)") +
  theme_bw(base_size = 10.5) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey35"))

mac <- matching[, c("scenario", "ratio_label", "mean_matched_mac", "mean_discordant_strata"), drop = FALSE]
mac_long <- rbind(data.frame(scenario = mac$scenario, ratio_label = mac$ratio_label, value = mac$mean_matched_mac, quantity = "Matched MAC"), data.frame(scenario = mac$scenario, ratio_label = mac$ratio_label, value = mac$mean_discordant_strata, quantity = "Discordant strata"))
mac_plot <- ggplot(mac_long, aes(ratio_label, value, group = 1)) + geom_line(colour = "#CC79A7", linewidth = .7) + geom_point(colour = "#CC79A7", size = 2.4) + facet_grid(quantity ~ scenario, scales = "free_y") + labs(title = "Information retained after matching", subtitle = "Higher matched MAC and more discordant strata improve rare-variant information", x = "Controls per case", y = "Mean value") + theme_bw(base_size = 10.5) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey35"))

fig <- (rate_plot / beta_plot / mac_plot) + plot_layout(heights = c(1, 1, 1.05))
pdf(file.path(out_root, "ratio_sensitivity_expanded.pdf"), width = 11.5, height = 11, bg = "white", useDingbats = FALSE); print(fig); dev.off()
ggsave(file.path(out_root, "ratio_sensitivity_expanded.png"), fig, width = 11.5, height = 11, units = "in", dpi = 240, bg = "white")
writeLines(c("Expanded ratio-sensitivity display for calibrated simulation", "QQ plots are intentionally similar under a well-calibrated null.", "This figure displays rejection rate, mean effect estimates, matched MAC, and discordant strata to reveal ratio-dependent information and power.", "No model refitting was performed."), file.path(out_root, "ratio_sensitivity_expanded_METHODS.txt"))
message("Created ratio sensitivity display: ", out_root)
