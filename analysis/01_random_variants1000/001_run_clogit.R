source("../../R/utils.R")

df_variants <- load_variants("chr2_random1000")
phenotype <- "qualification"
# phenotype = "LCT"
# phenotype <- "TV"
# phenotype <- "milk"
# phenotype <- "lactose_intolerance"

xCtrl <- 1
num_pcs_used <- 15
input_cov <- load_covariate()
df_pheno_binary <- load_phenotype(phenotype)

ids_case <- subset(df_pheno_binary, case == 1)$f.eid
ids_ctrl <- subset(df_pheno_binary, case == 0)$f.eid
input_cov$case <- ifelse(input_cov$IID %in% ids_case, 1, ifelse(input_cov$IID %in% ids_ctrl, 0, NA))


# if (pheno == "lactose_intolerance") {
#   path_phenotype <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/ICD10.txt"
#   df_pheno <- fread(path_phenotype)
#   df_pheno$index <- NULL
#   df_pheno <- df_pheno[!duplicated(df_pheno), ]
#   eid_lactose_intolerance <- subset(df_pheno, grepl("^E73", value) & value != "E731")$f.eid
#   eid_ctrl <- subset(df_pheno, !(grepl("^E73", value) & value != "E731"))$f.eid

#   input_cov$case <- ifelse(input_cov$IID %in% eid_lactose_intolerance, 1,
#     ifelse(input_cov$IID %in% eid_ctrl, 0, NA)
#   )
# } else if (pheno == "TV") {
#   path_phenotype <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/TV.txt"
#   if (!file.exists(path_phenotype)) {
#     path_pheno_head <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/ukbb_phenotype_head100.txt"
#     path_pheno <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/ukb672224.tab"

#     bd_head <- fread(path_pheno_head, header = TRUE, sep = "\t")
#     eid_col <- "f.eid"
#     ID <- "f.1070"
#     cols_needed <- grep(ID, colnames(bd_head), value = TRUE)
#     bd <- fread(path_pheno, select = c(eid_col, cols_needed), header = T, sep = "\t")

#     special_codes <- c(-1, -3)
#     bd_long <- bd %>%
#       pivot_longer(
#         cols = all_of(cols_needed),
#         names_to = "index",
#         values_to = "value"
#       )

#     bd_mean <- bd_long %>%
#       mutate(
#         value2 = ifelse(value == -10, 0,
#           ifelse(value %in% c(-1, -3), NA, value)
#         )
#       ) %>%
#       group_by(f.eid) %>%
#       summarise(
#         TV_mean = mean(value2, na.rm = TRUE)
#       ) %>%
#       filter(!is.na(TV_mean)) %>%
#       mutate(
#         TV_gt2h = ifelse(TV_mean < 2, 0, 1)
#       )

#     write.table(bd_mean, path_phenotype,
#       quote = F, row.names = F, col.names = T, sep = "\t"
#     )
#   } else {
#     bd_binary <- read.delim(path_phenotype)
#     input_cov$case <- ifelse(input_cov$IID %in% subset(bd_binary, TV_gt2h == 1)$f.eid, 1,
#       ifelse(input_cov$IID %in% subset(bd_binary, TV_gt2h == 0)$f.eid, 0, NA)
#     )
#   }
# } else if (pheno == "milk") {
#   path_phenotype <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/milk_type_used.txt"
#   if (!file.exists(path_phenotype)) {
#     path_pheno_head <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/ukbb_phenotype_head100.txt"
#     path_pheno <- "/home/seongwonhwang/Desktop/projects/rare_variant/pheno/ukb672224.tab"

#     bd_head <- fread(path_pheno_head, header = TRUE, sep = "\t")
#     eid_col <- "f.eid"
#     ID <- "f.1418"
#     cols_needed <- grep(ID, colnames(bd_head), value = TRUE)
#     bd <- fread(path_pheno, select = c(eid_col, cols_needed), header = T, sep = "\t")

#     bd_long <- bd %>%
#       pivot_longer(
#         cols = all_of(cols_needed),
#         names_to = "index",
#         values_to = "value"
#       )

#     bd_long <- bd_long %>%
#       filter(!is.na(value), !value %in% c(-1, -3))

#     bd_binary <- bd_long %>%
#       group_by(f.eid) %>%
#       summarise(
#         case = ifelse(any(value %in% c(1, 2, 3)), 1, 0)
#       )

#     write.table(bd_binary, path_phenotype,
#       quote = F, row.names = F, col.names = T, sep = "\t"
#     )
#   } else {
#     bd_binary <- read.delim(path_phenotype)
#     input_cov$case <- ifelse(input_cov$IID %in% subset(bd_binary, case == 1)$f.eid, 1,
#       ifelse(input_cov$IID %in% subset(bd_binary, case == 0)$f.eid, 0, NA)
#     )
#   }
# }

input_pcs <- load_pcs(num_pcs_used)
eigen_vals <- load_eigen()
df_AF <- load_MAF()

list_variants <- colnames(df_variants)[-c(1:6)]

date()
for (idx in 1:length(list_variants)) {
  print(idx)
  this_snp <- list_variants[idx]

  MAF <- subset(df_AF, snp == this_snp)$MAF

  # ---- PCAmatch 등 기존 코드 ----
  samp_intersect <- Reduce(
    "intersect",
    list(
      input_pcs$IID,
      subset(input_cov, !is.na(case))$IID,
      subset(df_variants, !is.na(df_variants[[this_snp]]))$IID
    )
  )
  input_pcs_sorted <- input_pcs[match(samp_intersect, input_pcs$IID), ]
  input_cov_sorted <- input_cov[match(samp_intersect, input_cov$IID), ]

  set.seed(1)

  if (phenotype %in% c("qualification", "TV", "milk")) {
    fname_pcamatch <- here::here(paste0("results/PCAmatch/", phenotype, "/match_random_1pct_", num_pcs_used, "PCs_", xCtrl, "xCtrl_", this_snp, ".RDS"))
    n_selected <- round(0.01 * length(samp_intersect))
    selected_idx <- sample(1:length(samp_intersect), n_selected)
  } else if (phenotype == "lactose_intolerance") {
    # # if (this_snp == "rs17653534_C") next
    # fname_pcamatch <- paste0(path_PCAmatch, "/match_res_", phenotype, "_random_10pct_", num_pcs_used, "PCs_1xCtrl_", this_snp, ".RDS")
    # case_idx <- which(input_cov_sorted$case == 1)
    # ctrl_idx <- which(input_cov_sorted$case == 0)
    # n_ctrl_selected <- round(0.1 * length(ctrl_idx))
    # ctrl_sub_idx <- sample(ctrl_idx, n_ctrl_selected)
    # selected_idx <- c(case_idx, ctrl_sub_idx)
  }
  input_pcs_sorted2 <- input_pcs_sorted[selected_idx, ]
  input_cov_sorted2 <- input_cov_sorted[selected_idx, ]

  if (!file.exists(fname_pcamatch)) {
    match_res <- match_maker(
      PC = input_pcs_sorted2,
      eigen_value = eigen_vals[1:num_pcs_used],
      data = input_cov_sorted2,
      ids = c("IID"),
      case_control = c("case"),
      num_controls = 1,
      eigen_sum = sum(eigen_vals)
    )
    saveRDS(match_res, file = fname_pcamatch)
  } else {
    match_res <- readRDS(fname_pcamatch)
  }

  df_reg <- match_res$PC_matches
  input_reg <- merge(df_reg, df_variants, by = "IID")

  input_unmatched <- merge(
    input_cov_sorted2, # IID + sex + case
    df_variants[, .(IID, snp_val = get(this_snp))], # IID + snp_val
    by = "IID"
  )
  colnames(input_unmatched)[colnames(input_unmatched) == "snp_val"] <- this_snp
  input_unmatched <- merge(input_unmatched, input_pcs_sorted2, by = "IID")

  # ---- 모델 fitting ----
  formula_clogit <- as.formula(paste("case ~ sex +", this_snp, "+ strata(match_final)"))
  clogit_fit <- survival::clogit(formula_clogit, data = input_reg)

  formula_glm1 <- as.formula(paste("case ~ sex +", this_snp))
  glm_fit1 <- stats::glm(formula_glm1, data = input_unmatched, family = binomial)

  pcs <- paste0("PC", 1:num_pcs_used)
  formula_glm2 <- as.formula(paste("case ~ sex +", this_snp, "+", paste(pcs, collapse = " + ")))
  glm_fit2 <- stats::glm(formula_glm2, data = input_unmatched, family = binomial)

  model_list <- list(
    SNP = this_snp,
    clogit = clogit_fit,
    glm_standard = glm_fit1,
    glm_pc_adjusted = glm_fit2
  )

  fname_model <- here::here(paste0("results/models/", phenotype, "/", this_snp, "_models.RDS"))
  saveRDS(model_list, fname_model)
}
date()
