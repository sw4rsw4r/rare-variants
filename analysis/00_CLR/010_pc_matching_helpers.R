match_ratio <- function(ratio, selected_rows) {
    selected_controls <- control_available
    candidate_rows <- c(selected_rows, selected_controls)
    candidate_data <- pc_data[
        candidate_rows,
        c("IID", pc_names),
        drop = FALSE
    ]
    candidate_data$case <- c(
        rep.int(1L, length(selected_rows)),
        rep.int(0L, length(selected_controls))
    )
    candidate_data <- candidate_data[
        ,
        c("IID", "case", pc_names),
        drop = FALSE
    ]
    rownames(candidate_data) <- candidate_data$IID

    match_object <- MatchIt::matchit(
        stats::reformulate(pc_names, response = "case"),
        data = candidate_data,
        method = "nearest",
        distance = "mahalanobis",
        estimand = "ATT",
        replace = FALSE,
        m.order = "farthest",
        ratio = ratio
    )
    matched_data <- MatchIt::match.data(
        match_object,
        data = candidate_data,
        drop.unmatched = TRUE
    )
    matched_data$STRATUM <- match(
        as.character(matched_data$subclass),
        unique(as.character(matched_data$subclass))
    )
    matched_data <- matched_data[
        ,
        c("IID", "case", pc_names, "STRATUM", "weights"),
        drop = FALSE
    ]
    matched_data <- matched_data[
        order(matched_data$STRATUM, -matched_data$case),
        ,
        drop = FALSE
    ]
    rownames(matched_data) <- NULL

    counts <- table(matched_data$STRATUM, matched_data$case)
    valid_sets <- all(c("0", "1") %in% colnames(counts)) &&
        all(counts[, "1"] == 1L) &&
        all(counts[, "0"] == ratio)
    if (!valid_sets) {
        stop("Incomplete 1:", ratio, " matching")
    }
    if (anyDuplicated(matched_data$IID)) {
        stop("Duplicate matched IID")
    }

    balance <- compute_balance(candidate_data, matched_data)
    summary <- data.frame(
        PHENOTYPE = phenotype,
        RATIO = ratio,
        N_CASE_AVAILABLE = length(case_available),
        N_CASE_SELECTED = length(selected_rows),
        N_CONTROL_AVAILABLE = length(control_available),
        N_CONTROL_CANDIDATE = length(selected_controls),
        N_CASE_MATCHED = sum(matched_data$case == 1L),
        N_CONTROL_MATCHED = sum(matched_data$case == 0L),
        N_STRATA = length(unique(matched_data$STRATUM)),
        MAX_ABS_SMD_PRE = max(abs(balance$PRE_SMD)),
        MAX_ABS_SMD_POST = max(abs(balance$POST_SMD)),
        stringsAsFactors = FALSE
    )

    ratio_dir <- file.path(
        matching_root,
        paste0("ratio_1to", ratio)
    )
    dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
    write.table(
        matched_data,
        file.path(ratio_dir, "matched_sets.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    write.table(
        balance,
        file.path(ratio_dir, "pc_balance.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    write.table(
        summary,
        file.path(ratio_dir, "matching_summary.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    list(
        data = matched_data,
        balance = balance,
        summary = summary
    )
}

run_clr <- function(genotype, matched_data, ratio) {
    n_variants <- ncol(genotype)
    result <- data.frame(
        SNP = colnames(genotype),
        N_MATCHED_BASE = nrow(matched_data),
        N_STRATA_BASE = length(unique(matched_data$STRATUM)),
        N = integer(n_variants),
        N_CASE = integer(n_variants),
        N_CONTROL = integer(n_variants),
        N_STRATA = integer(n_variants),
        N_STRATA_COMPLETE = integer(n_variants),
        N_STRATA_REDUCED = integer(n_variants),
        N_INFORMATIVE_STRATA = integer(n_variants),
        N_GENOTYPE_MISSING = integer(n_variants),
        EAF = rep.int(NA_real_, n_variants),
        MAF = rep.int(NA_real_, n_variants),
        MAC = rep.int(NA_real_, n_variants),
        BETA = rep.int(NA_real_, n_variants),
        SE = rep.int(NA_real_, n_variants),
        Z = rep.int(NA_real_, n_variants),
        CHISQ = rep.int(NA_real_, n_variants),
        P = rep.int(NA_real_, n_variants),
        ITERATIONS = integer(n_variants),
        STATUS = rep.int("not_tested", n_variants),
        WARNING = rep.int("", n_variants),
        TEST = rep.int(
            paste0("MATCHED_CLR_1TO", ratio, "_WALD_CHISQ1"),
            n_variants
        ),
        COVARIATES = rep.int("NONE_PC_MATCHED", n_variants),
        stringsAsFactors = FALSE
    )

    for (j in seq_len(n_variants)) {
        g <- genotype[, j]
        nonmissing <- is.finite(g)
        result$N_GENOTYPE_MISSING[j] <- sum(!nonmissing)

        case_present <- tapply(
            nonmissing & matched_data$case == 1L,
            matched_data$STRATUM,
            any
        )
        control_nonmissing <- tapply(
            nonmissing & matched_data$case == 0L,
            matched_data$STRATUM,
            sum
        )
        valid_strata <- as.integer(names(case_present)[
            case_present & control_nonmissing >= 1L
        ])
        use <- nonmissing & matched_data$STRATUM %in% valid_strata
        if (!any(use)) {
            result$STATUS[j] <- "no_valid_strata"
            next
        }

        model_data <- data.frame(
            case = matched_data$case[use],
            SNP = g[use],
            STRATUM = matched_data$STRATUM[use]
        )
        controls_per_stratum <- table(
            model_data$STRATUM[model_data$case == 0L]
        )
        informative <- tapply(
            model_data$SNP,
            model_data$STRATUM,
            function(v) length(unique(v)) > 1L
        )

        result$N[j] <- nrow(model_data)
        result$N_CASE[j] <- sum(model_data$case == 1L)
        result$N_CONTROL[j] <- sum(model_data$case == 0L)
        result$N_STRATA[j] <- length(unique(model_data$STRATUM))
        result$N_STRATA_COMPLETE[j] <- sum(
            controls_per_stratum == ratio
        )
        result$N_STRATA_REDUCED[j] <- sum(
            controls_per_stratum < ratio
        )
        result$N_INFORMATIVE_STRATA[j] <- sum(informative)

        allele_sum <- sum(model_data$SNP)
        result$EAF[j] <- allele_sum / (2 * nrow(model_data))
        result$MAF[j] <- min(result$EAF[j], 1 - result$EAF[j])
        result$MAC[j] <- min(
            allele_sum,
            2 * nrow(model_data) - allele_sum
        )
        if (!any(informative)) {
            result$STATUS[j] <- "no_within_stratum_variation"
            next
        }

        warning_text <- character()
        fit <- tryCatch(
            withCallingHandlers(
                survival::clogit(
                    case ~ SNP + strata(STRATUM),
                    data = model_data,
                    method = "exact",
                    control = survival::coxph.control(
                        iter.max = 50L,
                        eps = 1e-9
                    )
                ),
                warning = function(w) {
                    warning_text <<- c(
                        warning_text,
                        conditionMessage(w)
                    )
                    invokeRestart("muffleWarning")
                }
            ),
            error = function(e) e
        )
        if (inherits(fit, "error")) {
            result$STATUS[j] <- "fit_error"
            result$WARNING[j] <- conditionMessage(fit)
            next
        }

        result$ITERATIONS[j] <- if (length(fit$iter) > 0L) {
            max(fit$iter)
        } else {
            NA_integer_
        }
        result$WARNING[j] <- paste(
            unique(warning_text),
            collapse = " | "
        )
        if (any(grepl(
            "infinite|converged before",
            warning_text,
            ignore.case = TRUE
        ))) {
            result$STATUS[j] <- "possible_separation"
            next
        }

        coefficient_table <- summary(fit)$coefficients
        if (!("SNP" %in% rownames(coefficient_table))) {
            result$STATUS[j] <- "snp_coefficient_missing"
            next
        }
        beta <- coefficient_table["SNP", "coef"]
        se <- coefficient_table["SNP", "se(coef)"]
        if (!is.finite(beta) || !is.finite(se) || se <= 0) {
            result$STATUS[j] <- "invalid_coefficient"
            next
        }

        z <- beta / se
        result$BETA[j] <- beta
        result$SE[j] <- se
        result$Z[j] <- z
        result$CHISQ[j] <- z^2
        result$P[j] <- pchisq(z^2, 1, lower.tail = FALSE)
        result$STATUS[j] <- "ok"
    }

    result
}

draw_original_style <- function(comparison, columns, labels, max_val) {
    rows <- lapply(seq_along(columns), function(i) {
        values <- sort(
            comparison[[columns[i]]][
                is.finite(comparison[[columns[i]]])
            ]
        )
        data.frame(
            exp_chisq = qchisq(ppoints(length(values)), 1),
            obs_chisq = values,
            Model = labels[i],
            stringsAsFactors = FALSE
        )
    })
    qq <- do.call(rbind, rows)
    qq$Model <- factor(qq$Model, levels = labels)

    stats <- do.call(rbind, lapply(labels, function(label) {
        values <- qq$obs_chisq[qq$Model == label]
        data.frame(
            Model = label,
            lambda = median(values) / qchisq(.5, 1),
            stringsAsFactors = FALSE
        )
    }))
    stats$Model <- factor(stats$Model, levels = labels)
    anno <- data.frame(
        Model = stats$Model,
        label = sprintf("lambda = %.3f", stats$lambda),
        x = max_val * .1,
        y = max_val * .9
    )

    ggplot(qq, aes(exp_chisq, obs_chisq)) +
        geom_point(size = 1, alpha = .6) +
        geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed",
            color = "red"
        ) +
        facet_wrap(~Model, nrow = 1) +
        geom_text(
            data = anno,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 1,
            size = 3.5
        ) +
        labs(
            title = "QQ-plot of SNP test statistics by model",
            x = "Expected chi-square(1)",
            y = "Observed chi-square(1)"
        ) +
        theme_minimal() +
        scale_x_continuous(limits = c(0, max_val)) +
        scale_y_continuous(limits = c(0, max_val))
}
