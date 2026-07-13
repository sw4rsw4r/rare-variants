source("../../R/utils.R")
xCtrl <- 4
num_pcs_used <- 5
input_pcs <- load_pcs(num_pcs_used)
eigen_vals <- load_eigen()
input_cov <- load_covariate()
set.seed(100)
for (region in c("chr2_random100", "LCT_region", "PCSK9_region")) {
    if (region == "LCT_region") list_phenotype <- c("TV", "qualification", "lactose_intolerance")
    if (region == "PCSK9_region") list_phenotype <- c("TV", "qualification", "CHD")
    if (region == "chr2_random100") list_phenotype <- c("lactose_intolerance")

    for (phenotype in list_phenotype) {
        df_variants <- load_variants(region)
        case_eid <- load_phenotype(phenotype)

        input_cov$case <- ifelse(input_cov$IID %in% case_eid$case, 1, ifelse(input_cov$IID %in% case_eid$ctrl, 0, NA))
        list_variants <- colnames(df_variants)[-c(1:6)]

        path_res <- here::here(paste0("results/00_CLR/", region, "/", phenotype, "/", num_pcs_used, "PCs/", xCtrl, "xCtrl/"))
        check_directory(path_res)

        # 변이 데이터에 결측치가 없는 샘플 교집합 (일단 IID와 페노타입 기준으로만)
        base_intersect <- Reduce(
            "intersect",
            list(
                input_pcs$IID,
                subset(input_cov, !is.na(case))$IID,
                df_variants$IID
            )
        )
        base_pcs <- input_pcs[input_pcs$IID %in% base_intersect, ]
        base_cov <- input_cov[input_cov$IID %in% base_intersect, ]
        case_idx <- which(base_cov$case == 1)
        ctrl_idx <- which(base_cov$case == 0)

        n_case_selected <- min(length(case_idx), 1000)
        case_sub_idx <- sample(case_idx, n_case_selected)
        selected_case_iids <- base_cov$IID[case_sub_idx]

        # 제외 처리된 순수 컨트롤 후보군
        pure_ctrl_idx <- ctrl_idx[!(base_cov$IID[ctrl_idx] %in% selected_case_iids)]
        if (length(pure_ctrl_idx) > (n_case_selected * 15)) {
            pure_ctrl_idx <- sample(pure_ctrl_idx, n_case_selected * 15)
        }
        selected_idx <- c(case_sub_idx, pure_ctrl_idx)
        input_pcs_match_ready <- base_pcs[selected_idx, ]
        input_cov_match_ready <- base_cov[selected_idx, ]

        # 공통 매칭 결과 파일명 설정 (SNP 이름 대신 'common' 사용)
        fname_pcamatch_common <- here::here(paste0(path_res, "/PCAmatch_common_1pct.RDS"))

        # if (!file.exists(fname_pcamatch_common)) {
        if(T){
            match_res <- match_maker(
                PC = input_pcs_match_ready,
                eigen_value = eigen_vals[1:num_pcs_used],
                data = input_cov_match_ready,
                ids = c("IID"),
                case_control = c("case"),
                num_controls = xCtrl,
                eigen_sum = sum(eigen_vals)
            )
            saveRDS(match_res, file = fname_pcamatch_common)
        } else {
            match_res <- readRDS(fname_pcamatch_common)
        }

        # 최종 매칭에 선발된 고정 샘플셋 (Case + 4배수 매칭 Ctrl만 남음)
        df_reg_fixed <- match_res$PC_matches
        # ----------------------------------------------------------------------

        # 이제 가벼워진 상태로 SNP 루프 진입
        for (idx in 1:length(list_variants)) {
            this_snp <- list_variants[idx]
            if (!grepl("^rs", this_snp)) next
            if (region == "LCT_region" && this_snp %in% c("rs17653534_C")) next
            if (region == "PCSK9_region" && this_snp %in% c("rs115196099_A", "rs75050571_T")) next
            if (region == "chr2_random100" && this_snp %in% c("rs74374620_G", "rs12615203_C", "rs4668740_T")) next

            fname_model <- here::here(paste0(path_res, "/model_", this_snp, ".RDS"))
            # if (file.exists(fname_model)) next # 이미 계산된 SNP는 스킵

            # 현재 SNP에 결측치가 없는 샘플만 고정 매칭셋에서 필터링 (안전장치)
            valid_iids <- df_variants$IID[!is.na(df_variants[[this_snp]])]
            df_reg <- df_reg_fixed[df_reg_fixed$IID %in% valid_iids, ]

            # 단단하게 묶인 strata 구조 유지를 위해, 탈락한 사람으로 인해 짝이 깨진 stratum 방어
            # (보통 거의 없지만 완벽한 계산을 위해 탈락 샘플 처리)
            tab_strata <- table(df_reg$match_final)
            valid_strata <- names(tab_strata)[tab_strata == (xCtrl + 1)]
            df_reg <- df_reg[df_reg$match_final %in% valid_strata, ]

            # clogit용 데이터 생성
            input_reg <- merge(df_reg, df_variants[, c("IID", this_snp), with = FALSE], by = "IID")
            setDT(df_reg)
            # GLM 분석용 unmatched 데이터 (동일 샘플 매핑)
            input_unmatched <- merge(
                df_reg[, .(IID, case, sex)],
                df_variants[, .(IID, snp_val = get(this_snp))],
                by = "IID"
            )

            setnames(input_unmatched, "snp_val", this_snp)

            input_unmatched <- merge(
                input_unmatched,
                input_pcs,
                by = "IID"
            )

            # ---- model fitting ----
            formula_clogit <- as.formula(paste("case ~ sex +", this_snp, "+ strata(match_final)"))
            clogit_fit <- survival::clogit(formula_clogit, data = input_reg)

            formula_glm1 <- as.formula(paste("case ~ sex +", this_snp))
            glm_fit1 <- stats::glm(formula_glm1, data = input_unmatched, family = binomial)

            pcs <- paste0("PC", 1:num_pcs_used)
            formula_glm2 <- as.formula(paste("case ~ sex +", this_snp, "+", paste(pcs, collapse = " + ")))
            glm_fit2 <- stats::glm(formula_glm2, data = input_unmatched, family = binomial)

            # ---- permutation ----
            setorder(input_reg, match_final) 
            y <- input_reg$case
            stratum <- input_reg$match_final

            input_reg$y_perm <- permute_y(y, stratum)
            input_unmatched$y_perm2 <- sample(input_unmatched$case)

            formula_clogit_perm <- as.formula(paste("y_perm ~ sex +", this_snp, "+ strata(match_final)"))
            perm_clogit_fit <- survival::clogit(formula_clogit_perm, data = input_reg)

            formula_glm1_perm <- as.formula(paste("y_perm2 ~ sex +", this_snp))
            perm_glm_fit1 <- stats::glm(formula_glm1_perm, data = input_unmatched, family = binomial)

            formula_glm2_perm <- as.formula(paste("y_perm2 ~ sex +", this_snp, "+", paste(pcs, collapse = " + ")))
            perm_glm_fit2 <- stats::glm(formula_glm2_perm, data = input_unmatched, family = binomial)

            # 결과 저장
            model_list <- list(
                SNP = this_snp,
                clogit = clogit_fit,
                glm_standard = glm_fit1,
                glm_pc_adjusted = glm_fit2,
                perm_clogit = perm_clogit_fit,
                perm_glm_standard = perm_glm_fit1,
                perm_glm_pc_adjusted = perm_glm_fit2
            )
            saveRDS(model_list, fname_model)
        }
    }
}
