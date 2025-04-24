#' @importFrom dplyr mutate filter select inner_join %>% all_of case_when
#' @importFrom readr read_csv read_tsv read_lines read_table cols
#' @importFrom stringr str_sub str_trim str_split_1
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom shiny runApp fluidPage
#'
NULL

# Load required packages - handled by @importFrom or library() calls in app

#' Read Reference Data (Neanderthal SNP list)
#'
#' Reads the reference CSV file containing archaic SNP information.
#' It detects which allele is marked with an asterisk (indicating the Neanderthal-derived allele),
#' strips the asterisk, and stores that allele in a new column `archaic_allele`.
#'
#' The function attempts to find the CSV in the installed package's `extdata` directory.
#'
#' @return A tibble of reference SNP data with an added `archaic_allele` column.
#' @export
read_reference_data <- function() {
  # Determine the path to the reference CSV (installed or fallback)
  ref_path <- system.file("extdata", "neander_snps.csv", package = "snply")
  if (!file.exists(ref_path)) {
    stop("Reference data 'neander_snps.csv' not found in package 'snply' extdata directory.")
  }

  # Read the CSV
  ref_df <- readr::read_csv(ref_path, show_col_types = FALSE)

  # Standardize column names: replace spaces with dots
  colnames(ref_df) <- gsub(" ", ".", colnames(ref_df))

  # Convert key columns to character for consistent merging
  ref_df <- ref_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )

  # Check required columns (reference and alternative needed for archaic allele ID)
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
        TRUE ~ NA_character_ # Handle cases where neither is marked (though unlikely for this dataset)
      ),
      # Clean the original columns
      hg19_reference   = gsub("\\*", "", hg19_reference),
      hg19_alternative = gsub("\\*", "", hg19_alternative)
    ) %>%
    # Ensure archaic_allele is character type
    mutate(archaic_allele = as.character(archaic_allele))

  return(ref_df)
}

#' Read User Data (e.g., 23andMe genotype file)
#'
#' Reads the user's genotype file, accommodating headers that may
#' or may not be commented out with '#'. The file is expected to contain columns
#' including rsid, chromosome, position, and genotype, typically whitespace-delimited.
#'
#' @param filepath Path to the user's genotype file.
#' @return A tibble with columns: chromosome, position, and genotype.
#' @export
read_user_data <- function(filepath) {
  # Read the first few lines to find the header
  max_lines_to_check <- 50 # Adjust if headers can be further down
  lines <- tryCatch({
    readr::read_lines(filepath, n_max = max_lines_to_check)
  }, error = function(e) {
    stop("Error reading initial lines of the file: ", e$message)
  })


  header_line_num <- 0
  header_content <- NULL
  actual_header_names <- NULL

  # Find the first line containing the expected column names
  for (i in seq_along(lines)) {
    line <- stringr::str_trim(lines[i])
    # Check if it looks like the header (case-insensitive check might be safer)
    # Using fixed=TRUE can be slightly faster if case sensitivity is not needed
    if (grepl("rsid", line, ignore.case = TRUE) &&
        grepl("chromosome", line, ignore.case = TRUE) &&
        grepl("position", line, ignore.case = TRUE) &&
        grepl("genotype", line, ignore.case = TRUE)) {

      header_line_num <- i
      # Remove potential leading '#' and trim whitespace
      header_content <- stringr::str_trim(sub("^#", "", line))
      # Split by whitespace to get actual names
      actual_header_names <- stringr::str_split_1(header_content, "\\s+") # Splits on one or more spaces/tabs
      break
    }
  }

  if (header_line_num == 0) {
    stop("Could not find the header row (containing rsid, chromosome, position, genotype) within the first ", max_lines_to_check, " lines. Please ensure the file format is correct.")
  }

  # Read the data using read_table, skipping lines *before* the header
  # We tell it *not* to use the first data row as header and provide the names we found
  # We also disable the 'comment' character here, as we've already skipped past comments before the header.
  tryCatch({
    user_data <- readr::read_table(
      filepath,
      skip = header_line_num, # Skip lines up to and including the header
      col_names = actual_header_names, # Use the names we found
      col_types = readr::cols(.default = "c"), # Read all found columns as character initially
      comment = "", # Disable comment processing here

    )
  }, error = function(e) {
    stop("Error parsing data table after finding header. Check file structure below header. Original error: ", e$message)
  })


  # Define the standard column names we absolutely need
  standard_names <- c("rsid", "chromosome", "position", "genotype")

  # Check if the necessary standard columns were found in the header
  missing_standard_cols <- setdiff(standard_names, colnames(user_data))
  if (length(missing_standard_cols) > 0) {
    stop("The file header was found, but it is missing required columns: ", paste(missing_standard_cols, collapse = ", "))
  }

  # Select only the standard columns we need, discarding others that might exist
  user_data <- user_data %>%
    select(all_of(standard_names)) # Use all_of for safety

  # Ensure key columns are characters (partly redundant due to col_types="c", but safe)
  # and drop the 'rsid' column as per original logic
  user_data <- user_data %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position),
      genotype = as.character(genotype) # Ensure genotype is character
    ) %>%
    select(-rsid) # Drop rsid after ensuring it was present

  # Final check for expected columns after dropping rsid
  if(!all(c("chromosome", "position", "genotype") %in% colnames(user_data))) {
    stop("Processing error: Required columns (chromosome, position, genotype) are missing after selection.")
  }

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
  # Ensure key columns are characters in both datasets before merging
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

  # Perform the inner join
  merged_df <- tryCatch({
    inner_join(ref_df, user_df, by = c("chromosome", "position"))
  }, error = function(e) {
    stop("Error during merging of reference and user data: ", e$message)
  })

  if(nrow(merged_df) == 0) {
    warning("No matching SNPs found between the reference data and the user file based on chromosome and position.")
  }

  return(merged_df)
}


#' Compute Neanderthal Allele Statistics
#'
#' For each SNP in the merged data, counts the number of Neanderthal-derived alleles
#' in the user's genotype (0, 1, or 2) and computes summary statistics.
#'
#' @param merged_df Merged data frame from merge_snp_data().
#' @return A list with a summary string, total allele copies, percentage of allele copies, and a phenotype table.
#' @export
compute_allele_statistics <- function(merged_df) {
  # Handle empty input gracefully
  if (nrow(merged_df) == 0) {
    return(list(
      summary = "No matching archaic SNP sites found or tested.",
      count = 0,
      percentage = 0,
      phenotypes = data.frame(chromosome=character(), position=character()) # Empty frame with correct cols
    ))
  }

  # Proceed with calculation if data exists
  merged_df <- merged_df %>%
    mutate(
      # Ensure genotype has exactly two characters before splitting
      genotype = ifelse(nchar(genotype) == 2, genotype, "--"), # Replace invalid genotypes
      allele1 = str_sub(genotype, 1, 1),
      allele2 = str_sub(genotype, 2, 2),
      # Ensure archaic_allele is character for comparison
      archaic_allele = as.character(archaic_allele),
      # Calculate count, handling NAs and non-standard genotypes
      archaic_count = case_when(
        is.na(archaic_allele) | archaic_allele == "" ~ NA_integer_, # No archaic allele defined
        genotype == "--" | !grepl("^[ACGTDI]{2}$", genotype) ~ NA_integer_, # Invalid user genotype ('DI' for indel?)
        TRUE ~ (allele1 == archaic_allele) + (allele2 == archaic_allele)
      )
    )

  # Check if archaic_count column was successfully created
  if (!"archaic_count" %in% colnames(merged_df)) {
    stop("Internal error: 'archaic_count' column could not be generated.")
  }

  total_tested <- nrow(merged_df) # Total overlapping sites
  valid_comparisons <- sum(!is.na(merged_df$archaic_count)) # Sites where comparison was possible

  if (valid_comparisons == 0) {
    return(list(
      summary = paste0("Found ", total_tested, " overlapping sites, but none could be compared (check archaic allele definitions and genotype formats)."),
      count = 0,
      percentage = 0,
      phenotypes = data.frame(chromosome=character(), position=character())
    ))
  }

  # Calculate statistics based on valid comparisons
  total_allele_copies <- sum(merged_df$archaic_count, na.rm = TRUE)
  variant_count <- sum(merged_df$archaic_count > 0, na.rm = TRUE)
  heterozygous_count <- sum(merged_df$archaic_count == 1, na.rm = TRUE)
  homozygous_count <- sum(merged_df$archaic_count == 2, na.rm = TRUE)

  # Base percentage calculation on sites where comparison was possible
  percentage_alleles <- (total_allele_copies / (2 * valid_comparisons)) * 100
  percentage_sites <- (variant_count / valid_comparisons) * 100

  summary_text <- paste0(
    "Out of ", valid_comparisons, " comparable archaic SNP sites (out of ", total_tested, " total overlapping), you have ",
    homozygous_count, " sites with two Neanderthal variants and ",
    heterozygous_count, " sites with one Neanderthal variant. This totals ",
    variant_count, " sites carrying Neanderthal variants (", round(percentage_sites, 1), "% of comparable sites)."
  )

  # --- Phenotype Table ---
  # Define phenotype columns; adjust these names based on your reference CSV
  phenotype_cols <- c("Hair.colour", "Skin.colour", "Weight", "Standing.height", "Ease.of.skin.tanning")
  # Find which of these columns actually exist in the merged data
  available_pheno_cols <- intersect(phenotype_cols, colnames(merged_df))
  # Columns to select: always include chromosome and position
  select_cols <- c("chromosome", "position", "archaic_count", available_pheno_cols) # Add archaic_count here

  # Filter for sites with archaic alleles and select relevant columns
  phenotypes_df <- merged_df %>%
    filter(!is.na(archaic_count) & archaic_count > 0) %>%
    select(all_of(select_cols)) # Use all_of for safety

  # If no phenotype columns were found, create a basic table
  if (length(available_pheno_cols) == 0) {
    phenotypes_df <- phenotypes_df %>% select(chromosome, position, archaic_count)
    warning("No standard phenotype columns (e.g., 'Hair.colour') found in merged data.")
  }

  list(
    summary = summary_text,
    count = total_allele_copies, # Total count of archaic alleles found
    percentage = percentage_alleles, # Percentage relative to comparable alleles
    phenotypes = phenotypes_df
  )
}


#' Calculate Neanderthal Alleles Summary
#'
#' Orchestrates the reading, merging, and computation of Neanderthal allele statistics.
#' Includes error handling to return informative messages to the UI.
#'
#' @param user_file Path to the user's genotype file.
#' @return A list containing: 'summary' (string), 'count' (numeric),
#'         'percentage' (numeric), and 'phenotypes' (data.frame).
#'         On error, the list contains an error message in 'summary' and
#'         default/empty values for other elements.
#' @export
calculate_neanderthal_alleles <- function(user_file) {
  # Use tryCatch to gracefully handle errors in any step
  tryCatch({
    ref_df <- read_reference_data()
    user_df <- read_user_data(user_file)
    merged_df <- merge_snp_data(ref_df, user_df)
    stats <- compute_allele_statistics(merged_df)
    return(stats)

  }, error = function(e) {
    # Log the error (optional)
    # message("Error in calculate_neanderthal_alleles: ", e$message)

    # Return a structured list indicating an error
    error_message <- paste("Error processing file:", e$message)
    # Ensure a data frame is returned for phenotypes, even if empty, for DT consistency
    error_phenotypes <- data.frame(Error = error_message, stringsAsFactors = FALSE)

    return(list(
      summary = error_message,
      count = NA_integer_,
      percentage = NA_real_,
      phenotypes = error_phenotypes # Return the specific error data frame
    ))
  })
}
