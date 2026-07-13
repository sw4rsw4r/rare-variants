
library(data.table)

# read .freq
freq <- fread("../../data/plink_LCT.afreq")

# read header only from .raw
header <- names(fread("../../data/LCT_region_geno.raw", nrows = 0))

# SNP columns start after PHENOTYPE
snp_header <- header[-(1:6)]

# extract rsID and allele from header
header_dt <- data.table(
  raw_name = snp_header,
  ID = sub("_(.*)$", "", snp_header),
  header_allele = sub("^.*_", "", snp_header)
)

# merge with freq
check_dt <- merge(
  header_dt,
  freq[, .(ID, REF, ALT)],
  by = "ID",
  all.x = TRUE
)

# check whether header allele matches REF allele
check_dt[, match_REF := header_allele == REF]

table(check_dt$match_REF)
