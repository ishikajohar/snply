
#' @importFrom dplyr mutate filter select inner_join
#' @importFrom readr read_csv read_tsv
#' @importFrom stringr str_sub
#' @import ggplot2
#' @import shiny
#'
NULL

# Load required packages
library(readr)
library(dplyr)
library(stringr)
#'
#' Read Reference Data (Neanderthal SNP list)
#'
#' Reads the reference CSV file containing archaic SNP information.
#' It detects which allele is marked with an asterisk (indicating the Neanderthal-derived allele),
#' strips the asterisk, and stores that allele in a new column \code{archaic_allele}.
#'
#' The function attempts to find the CSV in the installed package's \code{extdata} directory,
#' and falls back to a relative path if not found.
#'
#' @return A tibble of reference SNP data with an added `archaic_allele` column.
#' @export
read_reference_data <- function() {
  # Determine the path to the reference CSV (installed or fallback)
  ref_path <- system.file("extdata", "neander_snps.csv", package = "snply")
  if (!file.exists(ref_path)) {
    stop("Reference data not found in package installation")
  }

  # Read the CSV
  ref_df <- read_csv(ref_path, show_col_types = FALSE)

  # Standardize column names: replace spaces with dots
  colnames(ref_df) <- gsub(" ", ".", colnames(ref_df))

  # Convert key columns to character for consistent merging
  ref_df <- ref_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )

  # Check required columns (we no longer require 'rsid' here)
  required_cols <- c("chromosome", "position", "hg19_reference", "hg19_alternative")
  missing_cols <- setdiff(required_cols, colnames(ref_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in reference data: ", paste(missing_cols, collapse = ", "))
  }

  # Identify and extract the archaic allele (marked with '*')
  ref_df <- ref_df %>%
    mutate(
      archaic_allele = case_when(
        grepl("\\*$", hg19_reference)   ~ gsub("\\*", "", hg19_reference),
        grepl("\\*$", hg19_alternative) ~ gsub("\\*", "", hg19_alternative),
        TRUE ~ NA_character_
      ),
      hg19_reference   = gsub("\\*", "", hg19_reference),
      hg19_alternative = gsub("\\*", "", hg19_alternative)
    )

  return(ref_df)
}

#' Read User Data (23andMe genotype file)
#'
#' Reads the user's 23andMe genotype file. The file is expected to be whitespace-delimited
#' with columns: rsid, chromosome, position, genotype. Lines starting with '#' are ignored.
#' (Since we're merging on chromosome and position only, the rsid column is dropped.)
#'
#' @param filepath Path to the user's 23andMe genotype file.
#' @return A tibble with columns: chromosome, position, and genotype.
#' @export
read_user_data <- function(filepath) {
  # Read the file using read_table() which splits on any whitespace.
  user_data <- readr::read_table(
    filepath,
    comment = "#",
    col_names = c("rsid", "chromosome", "position", "genotype"),
    col_types = "cccc"
  )

  # Check if the file has at least 4 columns.
  if(ncol(user_data) < 4) {
    stop("User data file does not have at least 4 columns. Please check the file format.")
  }

  # Restrict to only the first 4 columns (if more exist) and then drop the 'rsid' column.
  user_data <- user_data[, 1:4]
  colnames(user_data) <- c("rsid", "chromosome", "position", "genotype")
  user_data <- user_data %>% select(-rsid)

  # Ensure key columns are characters.
  user_data <- user_data %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )

  return(user_data)
}

#' Merge Reference and User SNP Data
#'
#' Merges the reference archaic SNP dataset with the user's SNP dataset on chromosome and position.
#'
#' @param ref_df Reference data frame from read_reference_data().
#' @param user_df User data frame from read_user_data().
#' @return A merged data frame containing only the SNPs present in both datasets.
#' @export
merge_snp_data <- function(ref_df, user_df) {
  # Convert key columns to character in both datasets
  ref_df <- ref_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )
  user_df <- user_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )
  # Merge by chromosome and position
  merged_df <- inner_join(ref_df, user_df, by = c("chromosome", "position"))
  return(merged_df)
}

#' Compute Neanderthal Allele Statistics
#'
#' For each SNP in the merged data, counts the number of Neanderthal-derived alleles in the user's genotype (0, 1, or 2)
#' and computes summary statistics.
#'
#' @param merged_df Merged data frame from merge_snp_data().
#' @return A list with a summary string, total allele copies, percentage of allele copies, and a phenotype table.
#' @export
compute_allele_statistics <- function(merged_df) {
  merged_df <- merged_df %>%
    mutate(
      allele1 = str_sub(genotype, 1, 1),
      allele2 = str_sub(genotype, 2, 2),
      archaic_count = (allele1 == archaic_allele) + (allele2 == archaic_allele)
    )

  total_tested <- nrow(merged_df)
  total_allele_copies <- sum(merged_df$archaic_count, na.rm = TRUE)
  variant_count <- sum(merged_df$archaic_count > 0, na.rm = TRUE)
  heterozygous_count <- sum(merged_df$archaic_count == 1, na.rm = TRUE)
  homozygous_count <- sum(merged_df$archaic_count == 2, na.rm = TRUE)
  percentage_alleles <- if (total_tested > 0) (total_allele_copies / (2 * total_tested)) * 100 else 0
  percentage_sites <- if (total_tested > 0) (variant_count / total_tested) * 100 else 0

  summary_text <- paste0(
    "Out of ", total_tested, " archaic SNP sites tested, you have ",
    homozygous_count, " sites with two Neanderthal variants and ",
    heterozygous_count, " sites with one Neanderthal variant, totaling ",
    variant_count, " sites with Neanderthal variants (", round(percentage_sites, 1), "%)."
  )

  # Define phenotype columns; adjust these names based on your CSV
  phenotype_cols <- c("Hair.colour", "Skin.colour", "Weight", "Standing.height", "Ease.of.skin.tanning")
  available_pheno_cols <- intersect(phenotype_cols, colnames(merged_df))
  # Since we don't have rsid in the reference, we include only chromosome and position plus phenotype info
  select_cols <- c("chromosome", "position", available_pheno_cols)

  phenotypes_df <- merged_df %>%
    filter(archaic_count > 0) %>%
    select(all_of(select_cols))

  list(
    summary = summary_text,
    count = total_allele_copies,
    percentage = percentage_alleles,
    phenotypes = phenotypes_df
  )
}

#' Calculate Neanderthal Alleles Summary
#'
#' Orchestrates the reading, merging, and computation of Neanderthal allele statistics.
#'
#' @param user_file Path to the user's 23andMe genotype file.
#' @return A list with a summary string, count of allele copies, percentage of allele copies, and a phenotype table.
#' @export
calculate_neanderthal_alleles <- function(user_file) {
  ref_df <- read_reference_data()
  user_df <- read_user_data(user_file)
  merged_df <- merge_snp_data(ref_df, user_df)
  stats <- compute_allele_statistics(merged_df)
  return(stats)
}




