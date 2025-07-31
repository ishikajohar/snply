# 0_preprocess_reference.R  ──────────────────────────────────────────────
library(data.table)

# (a) read the two *.rds that ship with snply v0.1.0
ref_freq  <- readRDS("reference_pop_freqs.rds")   # long data-frame  (~ Mb × pops)
ref_map   <- readRDS("ref_snp_map.rds")           # map: rsid → window_index, alleles

# (b) split the big freq df into a *named list* keyed by window_id.
#     each element is a plain matrix (rows = SNPs   cols = populations)
freq_by_window <- split(ref_freq, ref_freq$window_index)

freq_by_window <- lapply(freq_by_window, function(df) {
  # keep only rsid + pop columns
  df <- as.data.table(df)[, c("rsid", setdiff(names(df), c("rsid","window_index"))), with = FALSE]
  m  <- as.matrix(df[, -1])          # keep only pop columns
  rownames(m) <- df$rsid             # rownames = rsid
  m
})

# (c) save the light-weight list – MUCH faster to load in Shiny
saveRDS(freq_by_window,  "ref_freq_by_window.rds",  compress = "xz")
saveRDS(ref_map,          "ref_snp_map_light.rds", compress = "xz")
cat("✔  Reference split into", length(freq_by_window), "windows\n")
