source("../../R/utils.R")
library(dplyr)
library(data.table)
library(ggplot2)

num_pcs_used <- 5
xCtrl <- 4
all_qq_stats_uni <- NULL

# 가독성을 위한 모델명 치환 및 ggplot 정렬 순서 고정 함수
rename_models <- function(dt) {
    dt[, Model := factor(Model, 
                         levels = c("clogit", "glm_pc_adjusted", "glm_standard"),
                         labels = c("Matched (CLogit)", "Unmatched + PCs", "Unmatched (No PCs)"))]
    return(dt)
}

# 꼬임 없는 안전한 주석 매핑 함수
get_anno <- function(qq_stats, max_val) {
    data.table(
        Model = qq_stats$Model,
        label = sprintf("lambda = %.3f", qq_stats$lambda_GC),
        x = max_val * 0.1,  # 텍스트가 잘리지 않도록 좌측 상단 배치
        y = max_val * 0.9
    )
}

# 완벽하게 규격화된 QQ-plot 출력 함수
plot_qq <- function(qq_dt_frame, stats_frame, max_val, x_lab) {
    ggplot(qq_dt_frame, aes(x = exp_log, y = obs_log)) +
        geom_point(size = 1, alpha = 0.6) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        facet_wrap(~Model) +
        geom_text(
            data = get_anno(stats_frame, max_val),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5
        ) +
        labs(title = "QQ-plot of SNP p-values by model", x = x_lab, y = "Observed") +
        theme_minimal() +
        scale_x_continuous(limits = c(0, max_val)) +
        scale_y_continuous(limits = c(0, max_val))
}

# --- 메인 루프 가동 ---
for (region in c("chr2_random100", "LCT_region", "PCSK9_region")) {
    if (region == "LCT_region") list_phenotype <- c("TV", "qualification", "lactose_intolerance")
    if (region == "PCSK9_region") list_phenotype <- c("TV", "qualification", "CHD")
    if (region == "chr2_random100") list_phenotype <- c("lactose_intolerance")
    
    for (phenotype in list_phenotype) {
        for (isRare in c(FALSE, TRUE)) {
            
            model_dir <- here::here(paste0("results/00_CLR/", region, "/", phenotype, "/", num_pcs_used, "PCs/", xCtrl, "xCtrl"))
            files <- list.files(model_dir, pattern = "model_", full.names = TRUE)
            if (length(files) == 0) next
            
            # [버그 수정 1] .RDS 확장자를 지운 순수한 rs ID 추출
            clean_names <- gsub(".RDS$", "", basename(files))
            list_snps <- sapply(strsplit(clean_names, "_"), function(x) x[2])

            # [버그 수정 2] 뒤쪽 알릴 부호(_A, _T 등)를 완전히 도려내어 MAF ID와 포맷 통일
            list_snps <- gsub("_[A-Za-z]+$", "", list_snps)

            df_MAF <- load_MAF(region)
            list_MAF <- df_MAF[match(list_snps, df_MAF$ID), "MAF"]
            
            if (isRare) {
                files <- files[list_MAF <= .01]
            } else {
                files <- files[list_MAF > .01]
            }
            files <- files[!is.na(files)]
            if (length(files) == 0) next

            this_seed <- 1
            list_models <- c("clogit", "glm_standard", "glm_pc_adjusted")
            
            # 개별 SNP 모델 결과 취합
            res_list <- lapply(files, function(f) {
                obj <- readRDS(f)
                SNP <- obj$SNP
                d1 <- data.table(
                    SNP   = SNP,
                    Model = list_models,
                    obs_P = sapply(obj[list_models], extract_pval, snp_name = SNP),
                    obs_Z = sapply(obj[list_models], extract_z, snp_name = SNP)
                )
                d2 <- data.table(
                    exp_P = sapply(obj[paste0("perm_", list_models)], extract_pval, snp_name = SNP),
                    exp_Z = sapply(obj[paste0("perm_", list_models)], extract_z, snp_name = SNP)
                )
                cbind(d1, d2)
            })

            pval_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
            pval_dt[, obs_P := clean_list_col(obs_P)]
            pval_dt[, exp_P := clean_list_col(exp_P)]
            pval_dt <- pval_dt[is.finite(obs_P) & is.finite(exp_P)]
            if (nrow(pval_dt) == 0) next

            # [통계학 피드백 1] 모델명 매핑 전 깨끗한 상태에서 올바른 카이제곱 기반 Lambda GC 계산
            qq_stats_uni <- pval_dt[, {
                chisq_obs <- qchisq(1 - obs_P, df = 1)
                lambda_gc <- median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
                .(lambda_GC = lambda_gc, positive_rate = mean(obs_P < .05))
            }, by = Model]

            # [통계학 피드백 2] -log10 나눗셈 버그를 잡은 정통 Permutation Lambda 계산 
            qq_stats_perm <- pval_dt[, {
                chisq_obs <- qchisq(1 - obs_P, df = 1)
                chisq_perm <- qchisq(1 - exp_P, df = 1)
                lambda_perm <- median(chisq_obs, na.rm = TRUE) / median(chisq_perm, na.rm = TRUE)
                .(lambda_GC = lambda_perm, positive_rate = mean(obs_P < .05))
            }, by = Model]

            # 통계량 테이블 가독성 모델명 적용
            qq_stats_uni <- rename_models(qq_stats_uni)
            qq_stats_perm <- rename_models(qq_stats_perm)

            # [버그 수정 3] 대각선 밑으로 처지게 만들던 축 매핑 뒤틀림 버그 원천 차단
            # 각 모델 내부에서 독자적으로 순위를 정렬하여 좌표를 결합합니다.
            qq_dt_uni <- pval_dt[, {
                n <- .N
                data.table(
                    obs_log = -log10(sort(obs_P)),
                    exp_log = -log10(ppoints(n))
                )
            }, by = Model]
            qq_dt_uni <- rename_models(qq_dt_uni)

            qq_dt_perm <- pval_dt[, {
                data.table(
                    obs_log = -log10(sort(obs_P)),
                    exp_log = -log10(sort(exp_P))
                )
            }, by = Model]
            qq_dt_perm <- rename_models(qq_dt_perm)

            # --- 시각화 PDF 파일 저장 ---
            pdf(here::here(paste0("results/plots/QQ_plot_", region, "_", phenotype, "_isRare_", isRare, "_perm", this_seed, "_xCtrl", xCtrl, ".pdf")), width = 7, height = 3)
            
            # 1. Uniform 기댓값 기준 QQ Plot (전체 스케일 & 고정 스케일)
            max_val_uni <- with(qq_dt_uni, max(c(exp_log, obs_log), na.rm = TRUE))
            plot(plot_qq(qq_dt_uni, qq_stats_uni, max_val_uni, "Expected (Uniform)"))
            plot(plot_qq(qq_dt_uni, qq_stats_uni, 3, "Expected (Uniform)"))
            
            # 2. Permutation 기댓값 기준 QQ Plot (전체 스케일 & 고정 스케일)
            max_val_perm <- with(qq_dt_perm, max(c(exp_log, obs_log), na.rm = TRUE))
            plot(plot_qq(qq_dt_perm, qq_stats_perm, max_val_perm, "Expected (Permuted)"))
            plot(plot_qq(qq_dt_perm, qq_stats_perm, 3, "Expected (Permuted)"))
            
            dev.off()

            # 전체 요약 통계량 축적
            all_qq_stats_uni <- rbind(all_qq_stats_uni, data.frame(region, phenotype, isRare = isRare, qq_stats_uni))
        }
    }
}

# 최종 요약 스태츠 텍스트 파일 저장
write.table(all_qq_stats_uni, here::here(paste0("results/plots/all_qq_stats_uni_xCtrl", xCtrl, ".txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")