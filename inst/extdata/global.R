# global.R or top of app.R
ref_freqs <- readRDS("reference_pop_freqs.rds")   # Loads large list/dataframe of pop allele freqs by window
ref_map   <- readRDS("ref_snp_map.rds")           # SNP map: maps SNPs to windows (and possibly alleles, positions)
# (Optionally convert to data.table for fast subset)
library(data.table)
ref_map_dt <- as.data.table(ref_map)              # if ref_map is a data frame
setkey(ref_map_dt, rsid)                          # index by SNP ID for quick joins

