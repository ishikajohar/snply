# -------- R/neanderthal.R --------

#' @importFrom dplyr mutate filter select inner_join %>% all_of case_when
#' @importFrom readr read_csv read_tsv read_lines read_table cols
#' @importFrom stringr str_sub str_trim str_split_1
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom shiny runApp fluidPage
#'
NULL

#' Read Reference Data (Neanderthal SNP list)
#' @export
read_reference_data <- function() {
  ref_path <- system.file("extdata", "neander_snps.csv", package = "snply")
  if (!file.exists(ref_path)) {
    stop("Reference data 'neander_snps.csv' not found in package 'snply' extdata directory.")
  }
  ref_df <- readr::read_csv(ref_path, show_col_types = FALSE)
  colnames(ref_df) <- gsub(" ", ".", colnames(ref_df))
  ref_df <- ref_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )
  required_cols <- c("chromosome", "position", "hg19_reference", "hg19_alternative")
  missing_cols <- setdiff(required_cols, colnames(ref_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in reference data: ", paste(missing_cols, collapse = ", "))
  }
  ref_df <- ref_df %>%
    mutate(
      archaic_allele = case_when(
        grepl("\\*$", hg19_reference)   ~ gsub("\\*", "", hg19_reference),
        grepl("\\*$", hg19_alternative) ~ gsub("\\*", "", hg19_alternative),
        TRUE ~ NA_character_
      ),
      hg19_reference   = gsub("\\*", "", hg19_reference),
      hg19_alternative = gsub("\\*", "", hg19_alternative)
    ) %>%
    mutate(archaic_allele = as.character(archaic_allele))
  return(ref_df)
}


#' Read User Data (e.g., 23andMe genotype file)
#' Reads the user's genotype file, accommodating headers that may
#' or may not be commented out with '#'. Keeps rsid, chromosome, position, genotype.
#' @param filepath Path to the user's genotype file.
#' @return A tibble with columns: rsid, chromosome, position, and genotype.
#' @export
read_user_data <- function(filepath) {
  max_lines_to_check <- 50
  lines <- tryCatch({
    readr::read_lines(filepath, n_max = max_lines_to_check)
  }, error = function(e) {
    stop("Error reading initial lines of the file: ", e$message)
  })

  header_line_num <- 0
  header_content <- NULL
  actual_header_names <- NULL

  for (i in seq_along(lines)) {
    line <- stringr::str_trim(lines[i])
    if (grepl("rsid", line, ignore.case = TRUE) &&
        grepl("chromosome", line, ignore.case = TRUE) &&
        grepl("position", line, ignore.case = TRUE) &&
        grepl("genotype", line, ignore.case = TRUE)) {
      header_line_num <- i
      header_content <- stringr::str_trim(sub("^#", "", line))
      actual_header_names <- stringr::str_split_1(header_content, "\\s+")
      break
    }
  }

  if (header_line_num == 0) {
    stop("Could not find the header row (containing rsid, chromosome, position, genotype) within the first ", max_lines_to_check, " lines. Please ensure the file format is correct.")
  }

  tryCatch({
    user_data <- readr::read_table(
      filepath,
      skip = header_line_num,
      col_names = actual_header_names,
      col_types = readr::cols(.default = "c"), # Read all as character initially
      comment = ""
    )
  }, error = function(e) {
    stop("Error parsing data table after finding header. Check file structure below header. Original error: ", e$message)
  })

  standard_names <- c("rsid", "chromosome", "position", "genotype")
  missing_standard_cols <- setdiff(standard_names, colnames(user_data))
  if (length(missing_standard_cols) > 0) {
    stop("The file header was found, but it is missing required columns: ", paste(missing_standard_cols, collapse = ", "))
  }

  # Select only the standard columns, KEEPING rsid now
  user_data <- user_data %>%
    select(all_of(standard_names))

  #  key columns are characters
  user_data <- user_data %>%
    mutate(
      rsid = as.character(rsid), # Ensure rsid is character
      chromosome = as.character(chromosome),
      position = as.character(position),
      genotype = as.character(genotype)
    )

  if(!all(c("rsid", "chromosome", "position", "genotype") %in% colnames(user_data))) {
    stop("Processing error: Required columns (rsid, chromosome, position, genotype) are missing after selection.")
  }

  return(user_data)
}


#' Merge Reference and User SNP Data
#' @export
merge_snp_data <- function(ref_df, user_df) {
  ref_df <- ref_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    )
  # user_df now has rsid, but we only merge ref_df (no rsid) with user_df on chr/pos
  user_df_to_merge <- user_df %>%
    mutate(
      chromosome = as.character(chromosome),
      position = as.character(position)
    ) %>%
    select(chromosome, position, genotype) # Select only needed cols for merge + genotype

  merged_df <- tryCatch({
    # Merge ref data with selected user data (chromosome, position, genotype)
    inner_join(ref_df, user_df_to_merge, by = c("chromosome", "position"))
  }, error = function(e) {
    stop("Error during merging of reference and user data: ", e$message)
  })

  if(nrow(merged_df) == 0) {
    warning("No matching SNPs found between the reference data and the user file based on chromosome and position.")
  }
  return(merged_df)
}


#' Compute Neanderthal Allele Statistics
#' @export
compute_allele_statistics <- function(merged_df) {
  if (nrow(merged_df) == 0) {
    return(list(
      summary = "No matching archaic SNP sites found or tested.",
      count = 0,
      percentage = 0,
      # Return phenotype frame with expected columns, even if empty
      phenotypes = data.frame(chromosome=character(), position=character(), archaic_count=integer())
    ))
  }
  merged_df <- merged_df %>%
    mutate(
      genotype = ifelse(nchar(genotype) == 2, genotype, "--"),
      allele1 = str_sub(genotype, 1, 1),
      allele2 = str_sub(genotype, 2, 2),
      archaic_allele = as.character(archaic_allele),
      archaic_count = case_when(
        is.na(archaic_allele) | archaic_allele == "" ~ NA_integer_,
        genotype == "--" | !grepl("^[ACGTDI]{2}$", genotype) ~ NA_integer_,
        TRUE ~ (allele1 == archaic_allele) + (allele2 == archaic_allele)
      )
    )
  if (!"archaic_count" %in% colnames(merged_df)) {
    stop("Internal error: 'archaic_count' column could not be generated.")
  }
  total_tested <- nrow(merged_df)
  valid_comparisons <- sum(!is.na(merged_df$archaic_count))
  if (valid_comparisons == 0) {
    return(list(
      summary = paste0("Found ", total_tested, " overlapping sites, but none could be compared (check archaic allele definitions and genotype formats)."),
      count = 0,
      percentage = 0,
      phenotypes = data.frame(chromosome=character(), position=character(), archaic_count=integer())
    ))
  }
  total_allele_copies <- sum(merged_df$archaic_count, na.rm = TRUE)
  variant_count <- sum(merged_df$archaic_count > 0, na.rm = TRUE)
  heterozygous_count <- sum(merged_df$archaic_count == 1, na.rm = TRUE)
  homozygous_count <- sum(merged_df$archaic_count == 2, na.rm = TRUE)
  percentage_alleles <- (total_allele_copies / (2 * valid_comparisons)) * 100
  percentage_sites <- (variant_count / valid_comparisons) * 100
  summary_text <- paste0(
    "Out of ", valid_comparisons, " comparable archaic SNP sites (out of ", total_tested, " total overlapping), you have ",
    homozygous_count, " sites with two Neanderthal variants and ",
    heterozygous_count, " sites with one Neanderthal variant. This totals ",
    variant_count, " sites carrying Neanderthal variants (", round(percentage_sites, 1), "% of comparable sites)."
  )
  phenotype_cols <- c("Hair.colour", "Skin.colour", "Weight", "Standing.height", "Ease.of.skin.tanning")
  available_pheno_cols <- intersect(phenotype_cols, colnames(merged_df))
  # Define base columns for phenotype df
  select_cols_base <- c("chromosome", "position", "archaic_count")
  # Add available phenotype columns
  select_cols <- c(select_cols_base, available_pheno_cols)
  phenotypes_df <- merged_df %>%
    filter(!is.na(archaic_count) & archaic_count > 0) %>%
    select(all_of(select_cols))
  # Ensure base columns exist even if no phenotypes found/no rows
  if (nrow(phenotypes_df) == 0 && length(available_pheno_cols) == 0) {
    phenotypes_df <- data.frame(chromosome=character(), position=character(), archaic_count=integer())
  } else if (length(available_pheno_cols) == 0 && nrow(phenotypes_df) > 0) {
    phenotypes_df <- phenotypes_df %>% select(all_of(select_cols_base))
    warning("No standard phenotype columns (e.g., 'Hair.colour') found in merged data.")
  } else if (nrow(phenotypes_df) == 0 && length(available_pheno_cols) > 0) {
    # Create empty dataframe with all expected columns if no rows match filter
    empty_df_cols <- setNames(vector("list", length(select_cols)), select_cols)
    # Try to infer types or default to character
    empty_df_cols$chromosome <- character()
    empty_df_cols$position <- character()
    empty_df_cols$archaic_count <- integer()
    # Keep others as default list elements which data.frame handles
    phenotypes_df <- as.data.frame(empty_df_cols)
  }

  list(
    summary = summary_text,
    count = total_allele_copies,
    percentage = percentage_alleles,
    phenotypes = phenotypes_df
  )
}


#' Calculate Neanderthal Alleles Summary
#' @export
calculate_neanderthal_alleles <- function(user_file) {
  tryCatch({
    ref_df <- read_reference_data()
    # read_user_data now returns rsid, chromosome, position, genotype
    user_df_full <- read_user_data(user_file)
    # We need only chr, pos, genotype for merging with Neanderthal ref data
    # Create a version for merging
    user_df_for_merge <- user_df_full %>% select(chromosome, position, genotype)

    # Perform merge using the subsetted user data
    merged_df <- merge_snp_data(ref_df, user_df_for_merge) # Pass the subsetted data

    stats <- compute_allele_statistics(merged_df)
    return(stats)

  }, error = function(e) {
    error_message <- paste("Error processing file:", e$message)
    error_phenotypes <- data.frame(Error = error_message, stringsAsFactors = FALSE)
    return(list(
      summary = error_message,
      count = NA_integer_,
      percentage = NA_real_,
      phenotypes = error_phenotypes
    ))
  })
}

# -------- End of R/neanderthal.R --------
