#' Load AISNP reference allele frequencies from CSV
#'
#' Reads the package-bundled CSV of ancestry-informative SNPs and ensures all population
#' frequency columns are numeric, clipped to (ε, 1-ε), and free of NA rows.
#'
#' @return A data.table with columns: SNP, CHR, A1, CEU, CHB, GIH, JPT, YRI, FST
#' @export
load_aisnp_ref <- function() {
  # Locate the CSV file in the package's extdata
  ref_path <- system.file("extdata", "aisnp_freqs.csv", package = "snply")
  if (!nzchar(ref_path) || !file.exists(ref_path)) {
    stop("Reference frequency file 'aisnp_freqs.csv' not found in package 'snply'.")
  }

  # Read into data.table
  library(data.table)
  ref_dt <- fread(ref_path)

  # Ensure chromosome is character
  if ("CHR" %in% names(ref_dt)) {
    ref_dt[, CHR := as.character(CHR)]
  }

  # Define the population frequency columns explicitly
  pop_cols <- c("CEU", "CHB", "GIH", "JPT", "YRI")

  # Coerce population columns to numeric
  ref_dt[, (pop_cols) := lapply(.SD, function(x) as.numeric(x)), .SDcols = pop_cols]

  # Drop any rows with NA in these freq columns
  na_rows <- !complete.cases(ref_dt[, ..pop_cols])
  if (any(na_rows)) {
    warning("Dropping ", sum(na_rows), " rows with missing frequency values.")
    ref_dt <- ref_dt[!na_rows]
  }

  # Clip frequencies to (epsilon, 1 - epsilon) to avoid log(0)
  eps <- 1e-9
  ref_dt[, (pop_cols) := lapply(.SD, function(x) pmin(pmax(x, eps), 1 - eps)), .SDcols = pop_cols]

  return(ref_dt)
}
