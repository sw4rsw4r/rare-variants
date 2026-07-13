options("optmatch_max_problem_size" = Inf)
library(data.table)
library(here)
library(PCAmatchR)
library(optmatch)
library(dplyr)
library(tidyr)
library(survival)
library(ggplot2)

load_covariate <- function() {
  path_covariate <- here::here("data/UKBB_covariate.txt")
  df_cov <- data.table::fread(path_covariate)
  input_cov <- as.data.frame(df_cov[, c("IID", "sex")])
  return(input_cov)
}

load_phenotype0 <- function(ID) {
  path_pheno_head <- here::here("data/UKBB_phenotype_head100.txt")
  path_pheno <- here::here("data/ukb672224.tab")

  bd_head <- data.table::fread(path_pheno_head, header = TRUE, sep = "\t")
  cols_needed <- grep(ID, colnames(bd_head), value = TRUE)

  bd <- data.table::fread(path_pheno, select = c("f.eid", cols_needed), header = T, sep = "\t")
  bd_long <- bd %>%
    tidyr::pivot_longer(
      cols = all_of(cols_needed),
      names_to = "index",
      values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value))
  return(bd_long)
}

load_qualification <- function() {
  path_phenotype <- here::here("data/processed/qualification.RDS")
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.6138\\.")
    # > table(df_pheno$value)

    #     -7     -3      1      2      3      4      5      6
    #  91531   5833 203979 160835 269170  75534 108572 169252
    df_pheno <- df_pheno %>%
      dplyr::filter(value != -3)

    df_pheno_binary <- df_pheno %>%
      dplyr::group_by(f.eid) %>%
      dplyr::summarise(
        case = ifelse(any(value == 1), 1, 0)
      )
    # > table(df_pheno_binary$case)

    #      0      1
    # 328733 164796
    # > prop.table(table(df_pheno_binary$case))

    #         0         1
    # 0.6660865 0.3339135
    ids_case <- subset(df_pheno_binary, case == 1)$f.eid
    ids_ctrl <- subset(df_pheno_binary, case == 0)$f.eid
    eid_ctrl <- list(case = ids_case, ctrl = ids_ctrl)
    saveRDS(eid_ctrl, file = path_phenotype)
  } else {
    eid_ctrl <- readRDS(path_phenotype)
  }
  return(eid_ctrl)
}

load_qualification <- function() {
  path_phenotype <- here::here(
    "data/processed/qualification.RDS"
  )
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.6138\\.")
    # Remove "Prefer not to answer"
    df_pheno <- df_pheno %>%
      dplyr::filter(value != -3)
    # Case: has College/University degree (code 1)
    df_pheno_binary <- df_pheno %>%
      dplyr::group_by(f.eid) %>%
      dplyr::summarise(
        case = as.integer(any(value == 1)),
        .groups = "drop"
      )
    ids_case <- unique(
      subset(df_pheno_binary, case == 1)$f.eid
    )
    ids_ctrl <- unique(
      subset(df_pheno_binary, case == 0)$f.eid
    )
    case_eid <- list(
      case = ids_case,
      ctrl = ids_ctrl
    )
    saveRDS(case_eid, file = path_phenotype)
  } else {
    case_eid <- readRDS(path_phenotype)
  }
  return(case_eid)
}

load_lactose_intolerance <- function() {
  path_phenotype <- here::here(
    "data/processed/lactose_intolerance.RDS"
  )
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.41270\\.")
    # Lactose intolerance (E73*) excluding E731
    eid_lactose_intolerance <- unique(
      subset(
        df_pheno,
        grepl("^E73", value) & value != "E731"
      )$f.eid
    )
    all_eid <- unique(df_pheno$f.eid)
    eid_ctrl <- setdiff(
      all_eid,
      eid_lactose_intolerance
    )
    case_eid <- list(
      case = eid_lactose_intolerance,
      ctrl = eid_ctrl
    )
    saveRDS(case_eid, file = path_phenotype)
  } else {
    case_eid <- readRDS(path_phenotype)
  }
  return(case_eid)
}

load_CHD <- function() {
  path_phenotype <- here::here("data/processed/CHD.RDS")
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.41270\\.")
    eid_chd <- unique(
      subset(df_pheno, grepl("^I2[0-5]", value))$f.eid
    )
    all_eid <- unique(df_pheno$f.eid)
    eid_ctrl <- setdiff(all_eid, eid_chd)
    case_eid <- list(
      case = eid_chd,
      ctrl = eid_ctrl
    )
    saveRDS(case_eid, file = path_phenotype)
  } else {
    case_eid <- readRDS(path_phenotype)
  }

  return(case_eid)
}

load_TV <- function() {
  path_phenotype <- here::here(
    "data/processed/TV.RDS"
  )
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.1070\\.")
    special_codes <- c(-1, -3)
    df_pheno_binary <- df_pheno %>%
      mutate(
        value2 = case_when(
          value == -10 ~ 0,
          value %in% special_codes ~ NA_real_,
          TRUE ~ as.numeric(value)
        )
      ) %>%
      group_by(f.eid) %>%
      summarise(
        TV_mean = mean(value2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!is.na(TV_mean)) %>%
      mutate(
        case = ifelse(TV_mean >= 2, 1, 0)
      )
    ids_case <- unique(
      subset(df_pheno_binary, case == 1)$f.eid
    )
    ids_ctrl <- unique(
      subset(df_pheno_binary, case == 0)$f.eid
    )
    case_eid <- list(
      case = ids_case,
      ctrl = ids_ctrl
    )
    saveRDS(case_eid, file = path_phenotype)
  } else {
    case_eid <- readRDS(path_phenotype)
  }
  return(case_eid)
}

load_phenotype <- function(phenotype) {
  if (phenotype == "qualification") {
    case_eid <- load_qualification()
  }
  if (phenotype == "lactose_intolerance") {
    case_eid <- load_lactose_intolerance()
  }
  if (phenotype == "TV") {
    case_eid <- load_TV()
  }
  if (phenotype == "CHD") {
    case_eid <- load_CHD()
  }
  return(case_eid)
}

load_eigen <- function() {
  path_eigen <- here::here("data/UKBB_pca.eigenval")
  eigen_vals <- scan(path_eigen)
  return(eigen_vals)
}

load_pcs <- function(num_pcs_used) {
  path_pcs <- here::here("data/UKBB_pca.eigenvec") # plink2 was used
  df_pcs <- read.delim(path_pcs)
  df_pcs_subset <- df_pcs[, c("IID", paste0("PC", 1:num_pcs_used))]
  return(df_pcs_subset)
}

load_variants <- function(Type) {
  if (Type == "chr2_random1000") {
    path_variants <- here::here("data/chr2_random1000.raw")
  } else if (Type == "chr2_random100") {
    path_variants <- here::here("data/chr2_random100.raw")
  } else if (Type == "LCT_region") {
    path_variants <- here::here("data/LCT_region_geno.raw")
  } else if (Type == "PCSK9_region") {
    path_variants <- here::here("data/PCSK9_region_geno.raw")
  } else {
    stop("Type is wrong.")
  }
  df_variants <- fread(path_variants)
  return(df_variants)
}

load_MAF <- function(region) {
  if (region == "PCSK9_region") {
    path_AF <- here::here("data/plink_PCSK9.afreq")
  } else {
    path_AF <- here::here("data/plink_chr2.afreq")
  }

  df_AF <- read.delim(path_AF)
  df_AF$MAF <- with(df_AF, ifelse(ALT_FREQS < .5, ALT_FREQS, 1 - ALT_FREQS))
  df_AF <- df_AF[order(df_AF$MAF), ]

  df_AF$snp <- with(df_AF, paste0(ID, "_", REF))
  return(df_AF)
}

check_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  return(path)
}

permute_y <- function(y, stratum) {
  # 1. 인덱스, y, stratum을 매핑한 임시 data.table 생성
  dt <- data.table(
    idx = seq_along(y),
    y = y,
    stratum = stratum
  )
  
  # 2. 각 stratum 별로 그룹화(.BY)하여 y 값만 무작위 셔플(sample)
  # 이렇게 하면 각 매칭 그룹 내의 Case/Ctrl 비율은 유지되면서 행만 섞입니다.
  dt[, y_perm := sample(y), by = stratum]
  
  # 3. 원래 데이터 순서(idx)대로 정렬 후 섞인 y 값만 반환
  setorder(dt, idx)
  return(dt$y_perm)
}

clean_list_col <- function(col) {
  sapply(col, function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) {
      return(NA_real_) # NULL이나 NA는 숫자형 NA로 통일
    } else {
      return(as.numeric(x)) # 일반 숫자는 그대로 추출
    }
  })
}

extract_pval <- function(fit, snp_name) {
  s <- summary(fit)

  # clogit / glm 공통 처리
  if ("coefficients" %in% names(s)) {
    coef_tab <- s$coefficients
  } else {
    stop("Unknown model type")
  }
  is_firth <- !(snp_name %in% rownames(coef_tab))
  if (is_firth) {
    fit$prob[snp_name]
  } else {
    coef_tab[snp_name, "Pr(>|z|)"]
  }
}
extract_z <- function(fit, snp_name) {
  s <- summary(fit)

  # clogit / glm 공통 처리
  if ("coefficients" %in% names(s)) {
    coef_tab <- s$coefficients
  } else {
    stop("Unknown model type")
  }

  if (!(snp_name %in% rownames(coef_tab))) {
    return(NA_real_)
  }

  if (any(colnames(coef_tab) == "z")) {
    coef_tab[snp_name, "z"]
  } else {
    coef_tab[snp_name, "z value"]
  }
}
