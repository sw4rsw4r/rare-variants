load_phenotype0 <- function(id_pattern, phenotype_tab = "data/ukb672224.tab") {
  header <- data.table::fread(phenotype_tab, nrows = 0, header = TRUE, sep = "\t")
  cols_needed <- grep(id_pattern, names(header), value = TRUE)
  bd <- data.table::fread(
    phenotype_tab,
    select = c("f.eid", cols_needed),
    header = TRUE,
    sep = "\t"
  )
  bd_long <- tidyr::pivot_longer(
    as.data.frame(bd),
    cols = tidyr::all_of(cols_needed),
    names_to = "index",
    values_to = "value"
  )
  dplyr::filter(bd_long, !is.na(.data$value))
}

load_lactose_intolerance <- function(phenotype_tab = "data/ukb672224.tab") {
  path_phenotype <- "data/processed/lactose_intolerance.RDS"
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.41270\\.", phenotype_tab)
    eid_case <- unique(subset(df_pheno, grepl("^E73", value) & value != "E731")$f.eid)
    eid_ctrl <- setdiff(unique(df_pheno$f.eid), eid_case)
    dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(case = eid_case, ctrl = eid_ctrl), path_phenotype)
  }
  readRDS(path_phenotype)
}

load_TV <- function(phenotype_tab = "data/ukb672224.tab") {
  path_phenotype <- "data/processed/TV.RDS"
  if (!file.exists(path_phenotype)) {
    df_pheno <- load_phenotype0("f\\.1070\\.", phenotype_tab)
    special_codes <- c(-1, -3)
    df_pheno_binary <- df_pheno
    df_pheno_binary$value2 <- ifelse(
      df_pheno_binary$value == -10,
      0,
      ifelse(df_pheno_binary$value %in% special_codes, NA_real_, as.numeric(df_pheno_binary$value))
    )
    df_pheno_binary <- dplyr::summarise(
      dplyr::group_by(df_pheno_binary, .data$f.eid),
      TV_mean = mean(.data$value2, na.rm = TRUE),
      .groups = "drop"
    )
    df_pheno_binary <- dplyr::filter(df_pheno_binary, !is.na(.data$TV_mean))
    df_pheno_binary$case <- ifelse(df_pheno_binary$TV_mean >= 2, 1, 0)
    dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      list(
        case = unique(subset(df_pheno_binary, case == 1)$f.eid),
        ctrl = unique(subset(df_pheno_binary, case == 0)$f.eid)
      ),
      path_phenotype
    )
  }
  readRDS(path_phenotype)
}

check_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  path
}
