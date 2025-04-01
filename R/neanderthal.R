#' Calculate Neanderthal-derived Allele Percentage and Count
#'
#' Reads a user 23andMe genotype file and a reference Excel from a Neanderthal phenotype study,
#' and calculates how many Neanderthal-derived alleles the user carries at the studied SNPs.
#'
#' The reference Excel is expected to contain columns for chromosome, position, the reference and alternative alleles
#' (with an asterisk marking the Neanderthal-derived allele), and various phenotype columns indicating traits
#' associated with those SNPs (e.g., Hair colour, Skin colour, Weight, Standing height, Ease of skin tanning).
#' The user file is a text file in 23andMe raw data format (tab-separated columns: rsid, chromosome, position, genotype).
#'
#' The function identifies the Neanderthal-derived allele for each SNP (the allele marked with an asterisk in the reference),
#' then checks the user's genotype at those positions. It counts how many copies of the Neanderthal allele the user has (0, 1, or 2 per SNP),
#' and calculates the percentage of Neanderthal-derived alleles among the matched SNPs.
#'
#' @param user_file Path to the user’s SNP data file (23andMe format). Can be a relative or absolute path, using forward slashes or Windows-style backslashes.
#' @param ref_file Path to the reference Excel file from the Neanderthal phenotype study.
#' The Excel file should have column names in row 4 (with the first 3 rows being metadata or headers to skip).
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{count}}{Total number of Neanderthal-derived allele copies the user has at the matched SNPs.}
#'   \item{\code{percentage}}{Percentage of the user's alleles (at those SNPs) that are Neanderthal-derived. This is computed as \code{100 * (total archaic allele copies) / (2 * number of matched SNPs)}.}
#'   \item{\code{phenotypes}}{A data frame with one row per matched SNP, containing the SNP identifier and selected phenotype annotations.
#'         In this version, only the following columns are included: rsid, chromosome, position, Hair colour, Skin colour, Weight, Standing height, and Ease of skin tanning.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- calculate_neanderthal_alleles("path/to/user_snps.txt", "path/to/neanderthal_refs.xlsx")
#' print(result$count)
#' print(result$percentage)
#' head(result$phenotypes)
#' }
#'
#' @export
calculate_neanderthal_alleles <- function(user_file, ref_file) {
  # Normalize file paths for cross-platform compatibility
  user_file <- normalizePath(user_file, winslash = "/", mustWork = TRUE)
  ref_file  <- normalizePath(ref_file, winslash = "/", mustWork = TRUE)

  # Read the reference Excel file, skipping the first 3 metadata rows
  ref_df <- readxl::read_excel(ref_file, skip = 3)
  # Standardize column names: replace spaces with dots
  colnames(ref_df) <- sub(" ", ".", colnames(ref_df))

  # Identify key columns by name (case-insensitive)
  chr_col        <- grep("^chromosome$", tolower(colnames(ref_df)), value = TRUE)
  pos_col        <- grep("^position$",   tolower(colnames(ref_df)), value = TRUE)
  ref_allele_col <- grep("hg19_reference",   tolower(colnames(ref_df)), value = TRUE)
  alt_allele_col <- grep("hg19_alternative", tolower(colnames(ref_df)), value = TRUE)

  # Instead of using all remaining columns as phenotypes, select only the most important ones:
  selected_phenotypes <- c("Hair.colour", "Skin.colour", "Weight", "Standing.height", "Ease.of.skin.tanning")
  phenotype_cols <- intersect(colnames(ref_df), selected_phenotypes)

  # Subset the reference data to only the necessary columns
  ref_df <- ref_df[, c(chr_col, pos_col, ref_allele_col, alt_allele_col, phenotype_cols)]
  # Remove rows with missing chromosome or position
  ref_df <- ref_df[!is.na(ref_df[[chr_col]]) & !is.na(ref_df[[pos_col]]), ]

  # Identify the Neanderthal-derived allele for each SNP (marked with '*')
  ref_df$archaic_allele <- NA_character_
  ref_df$archaic_allele[grepl("\\*", ref_df[[ref_allele_col]])] <-
    gsub("\\*", "", ref_df[[ref_allele_col]][grepl("\\*", ref_df[[ref_allele_col]])])
  ref_df$archaic_allele[grepl("\\*", ref_df[[alt_allele_col]])] <-
    gsub("\\*", "", ref_df[[alt_allele_col]][grepl("\\*", ref_df[[alt_allele_col]])])
  # Clean allele columns by removing any '*' characters
  ref_df[[ref_allele_col]] <- gsub("\\*", "", ref_df[[ref_allele_col]])
  ref_df[[alt_allele_col]] <- gsub("\\*", "", ref_df[[alt_allele_col]])

  # Read the user SNP file (23andMe format) - expecting columns: rsid, chromosome, position, genotype
  user_df <- read.table(user_file, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  colnames(user_df) <- c("rsid", "chromosome", "position", "genotype")

  # Convert chromosome and position columns to character for consistent merging
  user_df$chromosome <- as.character(user_df$chromosome)
  ref_df[[chr_col]]  <- as.character(ref_df[[chr_col]])
  user_df$position   <- as.character(user_df$position)
  ref_df[[pos_col]]  <- as.character(ref_df[[pos_col]])

  # Merge the reference and user data on chromosome and position (inner join)
  merged_df <- merge(ref_df, user_df, by.x = c(chr_col, pos_col), by.y = c("chromosome", "position"))

  # If no SNPs match, return zeros
  if (nrow(merged_df) == 0) {
    return(list(count = 0, percentage = 0, phenotypes = data.frame()))
  }

  # Count Neanderthal-derived alleles in the user's genotype for each matched SNP
  archaic_counts <- mapply(function(geno, arch) {
    alleles <- strsplit(geno, "")[[1]]
    sum(alleles == arch)
  }, merged_df$genotype, merged_df$archaic_allele)

  total_count <- sum(archaic_counts)
  matched_snps <- nrow(merged_df)
  percentage   <- (total_count / (2 * matched_snps)) * 100

  # Prepare phenotype output: include rsid, chromosome, position, and selected phenotype columns
  output_cols <- c("rsid", chr_col, pos_col, phenotype_cols)
  output_cols <- output_cols[output_cols %in% colnames(merged_df)]
  phenotypes_df <- merged_df[, output_cols, drop = FALSE]

  list(
    count      = as.integer(total_count),
    percentage = as.numeric(percentage),
    phenotypes = phenotypes_df
  )
}
