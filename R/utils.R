options("optmatch_max_problem_size" = Inf)
library(data.table)
library(here)
library(PCAmatchR)
library(optmatch)
library(dplyr)
library(tidyr)
library(survival)
# library(ggplot2)

load_covariate <- function() {
  path_covariate <- here::here("data/UKBB_covariate.txt")
  df_cov <- data.table::fread(path_covariate)
  input_cov <- as.data.frame(df_cov[, c("IID", "sex")])
  return(input_cov)
}

load_phenotype <- function(ID) {
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
  path_phenotype <- here::here("data/processed/qualification.txt")
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype("f\\.6138\\.")
    # > table(bd_long$value)

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
    write.table(df_pheno_binary, path_phenotype,
      quote = F, row.names = F, col.names = T, sep = "\t"
    )
  } else {
    df_pheno_binary <- read.delim(path_phenotype)
  }
  return(df_pheno_binary)
}

load_phenotype <- function(phenotype) {
  if (phenotype == "qualification") {
    df_pheno_binary <- load_qualification()
  }
  return(df_pheno_binary)
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
  } else {
    stop("Type is wrong.")
  }
  df_variants <- fread(path_variants)
  return(df_variants)
}

load_MAF <- function() {
  path_AF_chr2 <- here::here("data/plink_chr2.afreq")
  
  df_AF <- read.delim(path_AF_chr2)
  df_AF$MAF <- with(df_AF, ifelse(ALT_FREQS < .5, ALT_FREQS, 1 - ALT_FREQS))
  df_AF <- df_AF[order(df_AF$MAF), ]

  df_AF$snp <- with(df_AF, paste0(ID, "_", REF))
  return(df_AF)
}


# fname_AF <- file.path(path_git, "data/plink_chr2.afreq")
# check_directory <- function(path) {
#   if (!dir.exists(path)){
#     dir.create(path, recursive = T)
#   }
#   return(path)
# }
