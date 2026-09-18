## =====================================================================
## Rare-variant simulation: oracle and matching comparison
## Methods: oracle GWAS, standard GWAS, PC-adjusted GWAS, matching + CLR
## Matching ratios: 1:1, 1:4, and 1:6
## The calibrated real-data-aligned design generates and uses PC1-PC5.
## =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(MatchIt)
  library(patchwork)
})

settings <- list(
  seed = 20260715L,
  N = 20000L,
  n_iter = 1000L,
  prop_A = 0.5,
  maf_A = 0.01,
  maf_B = 0.005,
  alpha = 0.05,
  intercept = -3,
  ancestry_log_odds = log(2),
  beta_positive = log(3),
  n_pcs = 2L,
  pc_separation = c(2.0, 1.0),
  pc_sd = 0.5,
  match_ratios = c(1L, 4L, 6L),
  ps_caliper_sd = 0.2,
  max_matched_cases = 1500L
)

args <- commandArgs(trailingOnly = TRUE)
plot_only <- "--plot-only" %in% args
calibrated_design <- "--calibrated" %in% args
get_arg <- function(prefix, default) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1L]])
}

settings$n_iter <- as.integer(get_arg("--n-iter", settings$n_iter))
requested_cores <- as.integer(get_arg(
  "--cores",
  min(4L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
))
n_cores <- min(settings$n_iter, requested_cores)

## Optional real-data-aligned sensitivity design.  The original design is
## preserved unchanged; this branch uses the same PC-only Mahalanobis,
## farthest-order, without-replacement matching logic as the real-data code,
## a common case pool across ratios, and a moderate positive effect so power
## is not saturated at 100%.
if (calibrated_design) {
  settings$N <- as.integer(get_arg("--N", 50000L))
  settings$maf_A <- 0.01
  settings$maf_B <- 0.001
  settings$beta_positive <- log(1.5)
  settings$max_matched_cases <- as.integer(get_arg("--max-matched-cases", 2500L))
  settings$matching_design <- "real_aligned_pc_mahalanobis"
}

## Keep the original two-PC design reproducible, but make the calibrated
## real-data-aligned analysis use PC1-PC5 by default.  Later PCs carry
## progressively weaker ancestry signal, reflecting the ordered PC structure
## in real data while avoiding five identical copies of the same latent axis.
default_n_pcs <- if (calibrated_design) 5L else 2L
settings$n_pcs <- as.integer(get_arg("--n-pcs", default_n_pcs))
pc_separation_template <- c(2.0, 1.0, 0.5, 0.25, 0.1)
if (!is.finite(settings$n_pcs) || settings$n_pcs < 1L || settings$n_pcs > length(pc_separation_template)) {
  stop("--n-pcs must be an integer between 1 and ", length(pc_separation_template))
}
settings$pc_separation <- pc_separation_template[seq_len(settings$n_pcs)]
settings$pc_names <- paste0("PC", seq_len(settings$n_pcs))

out_dir <- if (calibrated_design) {
  if (settings$n_pcs == 2L) file.path(getwd(), "calibrated_simulation_results") else file.path(getwd(), paste0("calibrated_simulation_pc1to", settings$n_pcs, "_results"))
} else if (settings$n_pcs == 2L) {
  file.path(getwd(), "revised_simulation_results")
} else file.path(getwd(), paste0("revised_simulation_pc1to", settings$n_pcs, "_results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c(
  "Oracle GWAS",
  "Standard GWAS",
  "PC-adjusted GWAS",
  "Matching + CLR"
)

safe_gwas_fit <- function(formula, data) {
  fit <- tryCatch(
    suppressWarnings(glm(formula, family = binomial(), data = data)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(c(beta = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  coef_table <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(coef_table) || !("G" %in% rownames(coef_table))) {
    return(c(beta = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  values <- coef_table["G", ]
  p_col <- grep("Pr\\(", names(values), value = TRUE)
  if (length(p_col) == 0L) {
    return(c(beta = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  c(
    beta = unname(values[["Estimate"]]),
    se = unname(values[["Std. Error"]]),
    p_value = unname(values[[p_col[[1L]]]])
  )
}

weighted_mean <- function(x, w) {
  sum(x * w) / sum(w)
}

pooled_sd <- function(x, group) {
  x1 <- x[group == 1L]
  x0 <- x[group == 0L]
  sqrt((stats::var(x1) + stats::var(x0)) / 2)
}

standardized_difference <- function(x, group, weights = NULL, denominator = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  if (is.null(denominator)) denominator <- pooled_sd(x, group)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)

  mean1 <- weighted_mean(x[group == 1L], weights[group == 1L])
  mean0 <- weighted_mean(x[group == 0L], weights[group == 0L])
  (mean1 - mean0) / denominator
}

## Conditional logistic regression for one case per stratum.
## This evaluates the same conditional likelihood used by clogit. The score-test
## p-value at beta = 0 is used because it is stable for sparse rare variants.
fit_conditional_logistic <- function(genotype_matrix) {
  n_strata <- nrow(genotype_matrix)
  if (n_strata == 0L) {
    return(c(
      beta = NA_real_, se = NA_real_, p_value = NA_real_,
      n_discordant = 0, matched_mac = 0
    ))
  }

  case_g <- genotype_matrix[, 1L]
  mean_g0 <- rowMeans(genotype_matrix)
  var_g0 <- rowMeans(genotype_matrix^2) - mean_g0^2
  information0 <- sum(var_g0)
  score0 <- sum(case_g - mean_g0)

  if (!is.finite(information0) || information0 <= 0) {
    return(c(
      beta = NA_real_, se = NA_real_, p_value = 1,
      n_discordant = 0,
      matched_mac = sum(genotype_matrix)
    ))
  }

  p_value <- 2 * pnorm(-abs(score0 / sqrt(information0)))

  beta <- 0
  information <- information0
  for (step in seq_len(50L)) {
    eta <- beta * genotype_matrix
    eta <- eta - apply(eta, 1L, max)
    exp_eta <- exp(eta)
    denominator <- rowSums(exp_eta)
    mean_g <- rowSums(exp_eta * genotype_matrix) / denominator
    mean_g2 <- rowSums(exp_eta * genotype_matrix^2) / denominator
    variance_g <- pmax(0, mean_g2 - mean_g^2)
    score <- sum(case_g - mean_g)
    information <- sum(variance_g)

    if (!is.finite(information) || information < 1e-10) break
    beta_new <- max(-15, min(15, beta + score / information))
    if (abs(beta_new - beta) < 1e-8) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }

  se <- if (is.finite(information) && information > 0) {
    1 / sqrt(information)
  } else {
    NA_real_
  }

  discordant <- apply(genotype_matrix, 1L, function(x) length(unique(x)) > 1L)
  c(
    beta = beta,
    se = se,
    p_value = p_value,
    n_discordant = sum(discordant),
    matched_mac = sum(genotype_matrix)
  )
}

run_matching <- function(data, ratio, settings, match_seed, fixed_case_indices = NULL) {
  cases <- which(data$Y == 1L)
  if (!is.null(fixed_case_indices)) {
    selected_cases <- intersect(as.integer(fixed_case_indices), cases)
    data <- data[data$Y == 0L | seq_len(nrow(data)) %in% selected_cases, ]
    row.names(data) <- NULL
    cases <- which(data$Y == 1L)
  } else if (length(cases) > settings$max_matched_cases) {
    set.seed(match_seed + 991L)
    selected_cases <- sample(cases, settings$max_matched_cases)
    data <- data[data$Y == 0L | seq_len(nrow(data)) %in% selected_cases, ]
    row.names(data) <- NULL
  }

  pc_formula <- stats::reformulate(settings$pc_names, response = "Y")
  pc_terms_formula <- stats::reformulate(settings$pc_names)
  ps_fit <- tryCatch(
    suppressWarnings(glm(pc_formula, family = binomial(), data = data)),
    error = function(e) NULL
  )
  if (is.null(ps_fit)) return(NULL)
  logit_ps <- as.numeric(predict(ps_fit, type = "link"))

  set.seed(match_seed)
  match_object <- tryCatch(
    suppressWarnings(if (identical(settings$matching_design, "real_aligned_pc_mahalanobis")) {
      MatchIt::matchit(
        pc_formula,
        data = data,
        method = "nearest",
        distance = "mahalanobis",
        ratio = ratio,
        replace = FALSE,
        estimand = "ATT",
        m.order = "farthest"
      )
    } else {
      MatchIt::matchit(
        pc_formula,
        data = data,
        method = "nearest",
        distance = logit_ps,
        mahvars = pc_terms_formula,
        caliper = settings$ps_caliper_sd,
        std.caliper = TRUE,
        ratio = ratio,
        replace = FALSE,
        estimand = "ATT",
        m.order = "random"
      )
    }),
    error = function(e) NULL
  )
  if (is.null(match_object)) return(NULL)

  matched <- MatchIt::match.data(match_object, data = data, drop.unmatched = TRUE)
  if (nrow(matched) == 0L || !"subclass" %in% names(matched)) return(NULL)

  stratum_size <- table(matched$subclass)
  valid_subclasses <- names(stratum_size[stratum_size == ratio + 1L])
  matched <- matched[as.character(matched$subclass) %in% valid_subclasses, ]
  if (nrow(matched) == 0L) return(NULL)

  matched$subclass <- droplevels(factor(matched$subclass))
  matched <- matched[order(as.integer(matched$subclass), -matched$Y), ]
  genotype_matrix <- matrix(
    matched$G,
    ncol = ratio + 1L,
    byrow = TRUE
  )
  clr <- fit_conditional_logistic(genotype_matrix)

  balance_variables <- c("ancestry", settings$pc_names)
  denominators <- vapply(
    balance_variables,
    function(variable) pooled_sd(data[[variable]], data$Y),
    numeric(1L)
  )
  before_smd <- vapply(
    balance_variables,
    function(variable) standardized_difference(
      data[[variable]], data$Y, denominator = denominators[[variable]]
    ),
    numeric(1L)
  )
  after_smd <- vapply(
    balance_variables,
    function(variable) standardized_difference(
      matched[[variable]], matched$Y, weights = matched$weights,
      denominator = denominators[[variable]]
    ),
    numeric(1L)
  )

  list(
    fit = clr,
    n_matched_cases = nrow(genotype_matrix),
    before_smd = before_smd,
    after_smd = after_smd
  )
}

result_row <- function(
    iteration,
    scenario,
    ratio,
    method,
    fit,
    n_matched_cases = NA_integer_,
    n_discordant = NA_integer_,
    matched_mac = NA_real_
) {
  data.frame(
    iteration = iteration,
    scenario = scenario,
    ratio = ratio,
    method = method,
    beta = unname(fit[["beta"]]),
    se = unname(fit[["se"]]),
    p_value = unname(fit[["p_value"]]),
    n_matched_cases = n_matched_cases,
    n_discordant = n_discordant,
    matched_mac = matched_mac
  )
}

run_scenario <- function(iteration, positive, base_data, settings) {
  scenario <- if (positive) "Positive control" else "Negative control"
  truth <- if (positive) settings$beta_positive else 0
  phenotype_seed_offset <- if (positive) 1000003L else 5000003L
  set.seed(settings$seed + iteration * 1009L + phenotype_seed_offset)

  linear_predictor <-
    settings$intercept +
    settings$ancestry_log_odds * base_data$ancestry +
    truth * base_data$G
  Y <- rbinom(settings$N, size = 1L, prob = plogis(linear_predictor))

  data <- data.frame(
    Y = Y,
    G = base_data$G,
    ancestry = base_data$ancestry
  )
  data <- cbind(data, as.data.frame(base_data[settings$pc_names]))

  standard_fit <- safe_gwas_fit(Y ~ G, data)
  pc_fit <- safe_gwas_fit(
    stats::reformulate(c("G", settings$pc_names), response = "Y"),
    data
  )
  oracle_fit <- safe_gwas_fit(Y ~ G + ancestry, data)

  result_rows <- list()
  balance_rows <- list()
  result_position <- 1L
  balance_position <- 1L

  fixed_case_indices <- NULL
  if (identical(settings$matching_design, "real_aligned_pc_mahalanobis")) {
    cases <- which(Y == 1L)
    set.seed(settings$seed + iteration * 7919L + as.integer(positive) * 100003L + 17L)
    fixed_case_indices <- sample(cases, min(length(cases), settings$max_matched_cases), replace = FALSE)
  }

  for (ratio in settings$match_ratios) {
    result_rows[[result_position]] <- result_row(
      iteration, scenario, ratio, "Oracle GWAS", oracle_fit
    )
    result_position <- result_position + 1L
    result_rows[[result_position]] <- result_row(
      iteration, scenario, ratio, "Standard GWAS", standard_fit
    )
    result_position <- result_position + 1L
    result_rows[[result_position]] <- result_row(
      iteration, scenario, ratio, "PC-adjusted GWAS", pc_fit
    )
    result_position <- result_position + 1L

    matching <- run_matching(
      data = data,
      ratio = ratio,
      settings = settings,
      match_seed = settings$seed + iteration * 7919L +
        as.integer(positive) * 100003L + ratio * 101L,
      fixed_case_indices = fixed_case_indices
    )

    if (is.null(matching)) {
      matching_fit <- c(beta = NA_real_, se = NA_real_, p_value = NA_real_)
      result_rows[[result_position]] <- result_row(
        iteration,
        scenario,
        ratio,
        "Matching + CLR",
        matching_fit
      )
    } else {
      result_rows[[result_position]] <- result_row(
        iteration,
        scenario,
        ratio,
        "Matching + CLR",
        matching$fit,
        n_matched_cases = matching$n_matched_cases,
        n_discordant = as.integer(matching$fit[["n_discordant"]]),
        matched_mac = matching$fit[["matched_mac"]]
      )

      for (variable in c("ancestry", settings$pc_names)) {
        balance_rows[[balance_position]] <- data.frame(
          iteration = iteration,
          scenario = scenario,
          ratio = ratio,
          variable = variable,
          before_smd = matching$before_smd[[variable]],
          after_smd = matching$after_smd[[variable]]
        )
        balance_position <- balance_position + 1L
      }
    }
    result_position <- result_position + 1L
  }

  list(
    results = do.call(rbind, result_rows),
    balance = if (length(balance_rows) > 0L) {
      do.call(rbind, balance_rows)
    } else {
      NULL
    }
  )
}

run_replicate <- function(iteration, settings) {
  set.seed(settings$seed + iteration * 1009L)
  ancestry <- rbinom(settings$N, size = 1L, prob = settings$prop_A)
  pc_matrix <- sapply(seq_len(settings$n_pcs), function(index) {
    rnorm(
      settings$N,
      mean = settings$pc_separation[[index]] * ancestry,
      sd = settings$pc_sd
    )
  })
  if (settings$n_pcs == 1L) pc_matrix <- matrix(pc_matrix, ncol = 1L)
  colnames(pc_matrix) <- settings$pc_names
  maf <- ifelse(ancestry == 1L, settings$maf_A, settings$maf_B)
  ## Proper diploid genotype dosage under Hardy-Weinberg equilibrium.
  G <- rbinom(settings$N, size = 2L, prob = maf)

  base_data <- c(list(ancestry = ancestry, G = G), as.list(as.data.frame(pc_matrix)))
  negative <- run_scenario(iteration, FALSE, base_data, settings)
  positive <- run_scenario(iteration, TRUE, base_data, settings)

  list(
    results = rbind(negative$results, positive$results),
    balance = rbind(negative$balance, positive$balance)
  )
}

if (plot_only) {
  results <- read.csv(file.path(out_dir, "replicate_results.csv"))
  balance <- read.csv(file.path(out_dir, "balance_by_replicate.csv"))
} else {
  message(
    "Running ", settings$n_iter, " replicates on ", n_cores,
    " core(s): N=", settings$N, ", ratios=",
    paste0("1:", settings$match_ratios, collapse = ", ")
  )

  if (n_cores > 1L) {
    cluster <- parallel::makeCluster(
      n_cores,
      outfile = file.path(out_dir, "parallel_workers.log")
    )
    parallel::clusterEvalQ(cluster, library(MatchIt))
    parallel::clusterExport(
      cluster,
      varlist = c(
        "safe_gwas_fit", "weighted_mean", "pooled_sd",
        "standardized_difference", "fit_conditional_logistic",
        "run_matching", "result_row", "run_scenario", "run_replicate"
      ),
      envir = environment()
    )
    replicate_results <- parallel::parLapply(
      cluster,
      seq_len(settings$n_iter),
      run_replicate,
      settings = settings
    )
    parallel::stopCluster(cluster)
  } else {
    replicate_results <- lapply(
      seq_len(settings$n_iter),
      run_replicate,
      settings = settings
    )
  }

  results <- do.call(rbind, lapply(replicate_results, `[[`, "results"))
  balance <- do.call(rbind, lapply(replicate_results, `[[`, "balance"))
}

results$ratio_label <- factor(
  paste0("1:", results$ratio),
  levels = paste0("1:", settings$match_ratios)
)
results$method <- factor(results$method, levels = method_levels)
results$scenario <- factor(
  results$scenario,
  levels = c("Negative control", "Positive control")
)

balance$ratio_label <- factor(
  paste0("1:", balance$ratio),
  levels = paste0("1:", settings$match_ratios)
)
balance$scenario <- factor(
  balance$scenario,
  levels = c("Negative control", "Positive control")
)

write.csv(
  results,
  file.path(out_dir, "replicate_results.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  balance,
  file.path(out_dir, "balance_by_replicate.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    parameter = names(settings),
    value = vapply(settings, function(x) paste(x, collapse = ","), character(1L))
  ),
  file.path(out_dir, "simulation_settings.csv"),
  row.names = FALSE
)

wilson_interval <- function(x, n, level = 0.95) {
  if (n == 0L) return(c(lower = NA_real_, upper = NA_real_))
  z <- qnorm(1 - (1 - level) / 2)
  proportion <- x / n
  denominator <- 1 + z^2 / n
  center <- (proportion + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(
    proportion * (1 - proportion) / n + z^2 / (4 * n^2)
  ) / denominator
  c(
    lower = max(0, center - half_width),
    upper = min(1, center + half_width)
  )
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  if (length(p) == 0L) return(NA_real_)
  p <- pmax(p, .Machine$double.xmin)
  median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

summarise_result_group <- function(data) {
  valid <- is.finite(data$p_value)
  p <- data$p_value[valid]
  beta <- data$beta[valid]
  se <- data$se[valid]
  scenario <- as.character(data$scenario[[1L]])
  truth <- if (scenario == "Positive control") settings$beta_positive else 0
  significant <- sum(p < settings$alpha)
  n_valid <- length(p)
  ci <- wilson_interval(significant, n_valid)
  coverage_valid <- is.finite(beta) & is.finite(se) & se > 0
  coverage <- if (any(coverage_valid)) {
    mean(
      beta[coverage_valid] - 1.96 * se[coverage_valid] <= truth &
        beta[coverage_valid] + 1.96 * se[coverage_valid] >= truth
    )
  } else {
    NA_real_
  }
  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }

  data.frame(
    scenario = scenario,
    ratio = data$ratio[[1L]],
    ratio_label = as.character(data$ratio_label[[1L]]),
    method = as.character(data$method[[1L]]),
    rejection_rate = if (n_valid > 0L) significant / n_valid else NA_real_,
    ci_lower = ci[["lower"]],
    ci_upper = ci[["upper"]],
    lambda_gc = lambda_gc(p),
    mean_beta = safe_mean(beta),
    bias = safe_mean(beta - truth),
    rmse = sqrt(safe_mean((beta - truth)^2)),
    coverage_95 = coverage,
    n_valid = n_valid,
    n_failed = nrow(data) - n_valid,
    mean_matched_cases = safe_mean(data$n_matched_cases),
    mean_discordant_strata = safe_mean(data$n_discordant),
    mean_matched_mac = safe_mean(data$matched_mac)
  )
}

result_groups <- split(
  results,
  interaction(results$scenario, results$ratio, results$method, drop = TRUE)
)
summary_metrics <- do.call(
  rbind,
  lapply(result_groups, summarise_result_group)
)
summary_metrics$ratio_label <- factor(
  summary_metrics$ratio_label,
  levels = paste0("1:", settings$match_ratios)
)
summary_metrics$method <- factor(summary_metrics$method, levels = method_levels)
summary_metrics$scenario <- factor(
  summary_metrics$scenario,
  levels = c("Negative control", "Positive control")
)
summary_metrics <- summary_metrics[
  order(summary_metrics$scenario, summary_metrics$ratio, summary_metrics$method),
]
row.names(summary_metrics) <- NULL

write.csv(
  summary_metrics,
  file.path(out_dir, "summary_metrics.csv"),
  row.names = FALSE,
  na = ""
)

balance_long <- rbind(
  data.frame(
    balance[, c("iteration", "scenario", "ratio", "ratio_label", "variable")],
    stage = "Before matching",
    smd = balance$before_smd
  ),
  data.frame(
    balance[, c("iteration", "scenario", "ratio", "ratio_label", "variable")],
    stage = "After matching",
    smd = balance$after_smd
  )
)
balance_long$abs_smd <- abs(balance_long$smd)
balance_groups <- split(
  balance_long,
  interaction(
    balance_long$scenario,
    balance_long$ratio,
    balance_long$variable,
    balance_long$stage,
    drop = TRUE
  )
)
balance_summary <- do.call(rbind, lapply(balance_groups, function(data) {
  values <- data$abs_smd[is.finite(data$abs_smd)]
  data.frame(
    scenario = as.character(data$scenario[[1L]]),
    ratio = data$ratio[[1L]],
    ratio_label = as.character(data$ratio_label[[1L]]),
    variable = data$variable[[1L]],
    stage = data$stage[[1L]],
    mean_abs_smd = mean(values),
    q05 = unname(quantile(values, 0.05)),
    median_abs_smd = median(values),
    q95 = unname(quantile(values, 0.95))
  )
}))
balance_summary$ratio_label <- factor(
  balance_summary$ratio_label,
  levels = paste0("1:", settings$match_ratios)
)

write.csv(
  balance_summary,
  file.path(out_dir, "balance_summary.csv"),
  row.names = FALSE,
  na = ""
)

## =====================================================================
## Publication plots
## =====================================================================

method_colours <- c(
  "Oracle GWAS" = "#000000",
  "Standard GWAS" = "#D55E00",
  "PC-adjusted GWAS" = "#0072B2",
  "Matching + CLR" = "#009E73"
)
method_shapes <- c(
  "Oracle GWAS" = 18,
  "Standard GWAS" = 15,
  "PC-adjusted GWAS" = 17,
  "Matching + CLR" = 16
)

plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey35"),
    plot.caption = element_text(size = 8, colour = "grey35", hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

rate_plot <- ggplot(
  summary_metrics,
  aes(
    x = ratio_label,
    y = rejection_rate,
    ymin = ci_lower,
    ymax = ci_upper,
    colour = method,
    shape = method,
    group = method
  )
) +
  geom_hline(
    data = data.frame(
      scenario = factor("Negative control", levels = levels(summary_metrics$scenario)),
      y = settings$alpha
    ),
    aes(yintercept = y),
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_line(linewidth = 0.65) +
  geom_errorbar(width = 0.12, linewidth = 0.5) +
  geom_point(size = 2.5) +
  facet_wrap(~scenario, nrow = 1) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = method_shapes) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  labs(
    title = "Calibration and alternative-scenario rejection rate",
    subtitle = "Error bars are 95% Wilson intervals",
    x = "Controls per case",
    y = "Significant results (%)"
  ) +
  plot_theme

ggsave(
  file.path(out_dir, "performance_summary.png"),
  rate_plot,
  width = 10,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "performance_summary.pdf"),
  rate_plot,
  width = 10,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

balance_plot_data <- balance_summary[
  balance_summary$stage == "After matching",
]
balance_plot_data$scenario <- factor(
  balance_plot_data$scenario,
  levels = c("Negative control", "Positive control")
)
balance_plot_data$variable <- factor(
  balance_plot_data$variable,
  levels = c("ancestry", settings$pc_names),
  labels = c("True ancestry", settings$pc_names)
)
balance_labels <- c("True ancestry", settings$pc_names)
balance_colours <- setNames(
  c("#CC79A7", "#0072B2", "#009E73", "#E69F00", "#D55E00", "#56B4E9")[seq_along(balance_labels)],
  balance_labels
)
balance_shapes <- setNames(c(15, 17, 16, 18, 8, 3)[seq_along(balance_labels)], balance_labels)

balance_plot <- ggplot(
  balance_plot_data,
  aes(
    x = ratio_label,
    y = median_abs_smd,
    ymin = q05,
    ymax = q95,
    colour = variable,
    shape = variable,
    group = variable
  )
) +
  geom_hline(
    yintercept = 0.1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_line(linewidth = 0.65) +
  geom_errorbar(width = 0.12, linewidth = 0.5) +
  geom_point(size = 2.4) +
  facet_wrap(~scenario, nrow = 1) +
  scale_colour_manual(values = balance_colours) +
  scale_shape_manual(values = balance_shapes) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  labs(
    title = "Balance after PC-based matching",
    subtitle = "Points: median absolute SMD; bars: 5th–95th percentiles",
    x = "Controls per case",
    y = "Absolute standardized mean difference"
  ) +
  plot_theme +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.box.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  file.path(out_dir, "matching_balance.png"),
  balance_plot,
  width = 10,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "matching_balance.pdf"),
  balance_plot,
  width = 10,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

## QQ panels: oracle, standard, PC-adjusted, and matching at each ratio.
qq_input <- results[
  results$method == "Matching + CLR" |
    (results$method != "Matching + CLR" & results$ratio == settings$match_ratios[[1L]]),
]
qq_input$panel <- ifelse(
  qq_input$method == "Matching + CLR",
  paste0("Matching + CLR ", qq_input$ratio_label),
  as.character(qq_input$method)
)
qq_panel_levels <- c(
  "Oracle GWAS",
  "Standard GWAS",
  "PC-adjusted GWAS",
  paste0("Matching + CLR 1:", settings$match_ratios)
)
qq_input$panel <- factor(qq_input$panel, levels = qq_panel_levels)

qq_groups <- split(
  qq_input,
  interaction(qq_input$scenario, qq_input$panel, drop = TRUE)
)
qq_data <- do.call(rbind, lapply(qq_groups, function(data) {
  p <- data$p_value[is.finite(data$p_value) & data$p_value >= 0 & data$p_value <= 1]
  p <- pmax(p, .Machine$double.xmin)
  data.frame(
    scenario = as.character(data$scenario[[1L]]),
    panel = as.character(data$panel[[1L]]),
    expected = -log10(ppoints(length(p))),
    observed = -log10(sort(p))
  )
}))
qq_data$scenario <- factor(
  qq_data$scenario,
  levels = c("Negative control", "Positive control")
)
qq_data$panel <- factor(qq_data$panel, levels = qq_panel_levels)
qq_expected_limit <- max(qq_data$expected, na.rm = TRUE) * 1.03

qq_plot <- ggplot(
  qq_data,
  aes(expected, observed, colour = scenario, shape = scenario)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.5
  ) +
  geom_point(size = 0.95, alpha = 0.72) +
  facet_grid(rows = vars(scenario), cols = vars(panel), scales = "free_y") +
  coord_cartesian(xlim = c(0, qq_expected_limit)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_colour_manual(values = c(
    "Negative control" = "#0072B2",
    "Positive control" = "#D55E00"
  )) +
  scale_shape_manual(values = c("Negative control" = 16, "Positive control" = 17)) +
  labs(
    title = "QQ plots across analysis methods",
    subtitle = paste0(settings$n_iter, " simulation replicates per panel"),
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p))),
    caption = "p = 1 values are retained. Positive-control departures should be interpreted against the corresponding null calibration."
  ) +
  plot_theme +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = grid::unit(0.65, "lines")
  )

ggsave(
  file.path(out_dir, "qq_controls_all_methods.png"),
  qq_plot,
  width = 15.5,
  height = 6.8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "qq_controls_all_methods.pdf"),
  qq_plot,
  width = 15.5,
  height = 6.8,
  units = "in",
  device = cairo_pdf
)

positive_metrics <- summary_metrics[
  summary_metrics$scenario == "Positive control",
]
effect_plot <- ggplot(
  positive_metrics,
  aes(
    x = ratio_label,
    y = mean_beta,
    colour = method,
    shape = method,
    group = method
  )
) +
  geom_hline(
    yintercept = settings$beta_positive,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = method_shapes) +
  labs(
    title = "Positive-control effect estimates",
    subtitle = paste0("True log(OR) = ", round(settings$beta_positive, 3)),
    x = "Controls per case",
    y = "Mean estimated log(OR)"
  ) +
  plot_theme

ggsave(
  file.path(out_dir, "effect_estimates.png"),
  effect_plot,
  width = 7.5,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "effect_estimates.pdf"),
  effect_plot,
  width = 7.5,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

message("Completed. Outputs written to: ", out_dir)
print(summary_metrics, row.names = FALSE)
