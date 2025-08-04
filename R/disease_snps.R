# -------- R/disease_snps.R --------

#' @importFrom dplyr filter select inner_join mutate relocate all_of %>%
#' @importFrom readr read_csv cols col_character
#' @importFrom tidyr pivot_longer pivot_wider unite separate replace_na
#' @importFrom rlang .data
#'
NULL

#' Read Disease SNP Reference Data
#'
#' Reads the reference CSV file containing common disease-associated SNP information.
#'
#' @param ref_path Path to the reference CSV file. Defaults to the version
#'   included in the package.
#' @return A tibble of reference disease SNP data.
#' @noRd
read_disease_snp_reference <- function(ref_path = NULL) {
  if (is.null(ref_path)) {
    ref_path <- system.file("extdata", "common_disease_snps.csv", package = "snply")
    if (!file.exists(ref_path)) {
      warning("Default disease SNP reference data not found in package 'snply' extdata directory. Returning empty frame.")
      # Return empty tibble with expected columns
      return(tibble::tibble(
        rsid = character(),
        chromosome = character(),
        position = character(),
        gene = character(),
        disease = character(),
        risk_allele = character(),
        notes = character()
      ))
    }
  }

  # Define column types, especially making sure position is character for joining
  col_types <- readr::cols(
    rsid = readr::col_character(),
    chromosome = readr::col_character(),
    position = readr::col_character(), # Read position as character
    gene = readr::col_character(),
    disease = readr::col_character(),
    risk_allele = readr::col_character(),
    notes = readr::col_character()
  )

  tryCatch({
    ref_df <- readr::read_csv(ref_path, show_col_types = FALSE, col_types = col_types)
    required_cols <- c("rsid", "chromosome", "position", "gene", "disease", "risk_allele")
    missing_cols <- setdiff(required_cols, colnames(ref_df))
    if (length(missing_cols) > 0) {
      stop("Missing required columns in disease SNP reference data: ", paste(missing_cols, collapse = ", "))
    }
    return(ref_df)
  }, error = function(e) {
    warning("Error reading disease SNP reference file '", ref_path, "': ", e$message, ". Returning empty frame.")
    return(tibble::tibble(
      rsid = character(),
      chromosome = character(),
      position = character(),
      gene = character(),
      disease = character(),
      risk_allele = character(),
      notes = character()
    ))
  })
}


#' Check User Genotypes for Common Disease-Associated SNPs
#'
#' Merges user's SNP data with a reference list of disease-associated SNPs
#' and reports the user's genotype at those locations.
#'
#' @param user_data A data frame or tibble containing user SNP data, MUST include
#'   columns 'rsid', 'chromosome', 'position', and 'genotype'. Typically the
#'   output of `read_user_data()`.
#' @param disease_ref A data frame or tibble containing the reference disease
#'   SNP information, typically from `read_disease_snp_reference()`.
#' @return A tibble summarizing the user's genotype for the found disease SNPs,
#'   including risk allele count, or an empty tibble if no matches or errors occur.
#' @export
#' @examples
#' \dontrun{
#'   #  user_snps_df is loaded from read_user_data()
#'   # and includes the rsid column
#'   disease_results <- check_disease_snps(user_snps_df)
#'   print(disease_results)
#' }
check_disease_snps <- function(user_data, disease_ref = NULL) {

  # Input validation
  required_user_cols <- c("rsid", "chromosome", "position", "genotype")
  if (!is.data.frame(user_data) || !all(required_user_cols %in% colnames(user_data))) {
    stop("Invalid 'user_data'. Must be a data frame with columns: ", paste(required_user_cols, collapse = ", "))
  }
  # Ensure key user columns are character for join
  user_data_for_join <- user_data %>%
    mutate(
      chromosome = as.character(.data$chromosome),
      position = as.character(.data$position)
    )

  # Load default reference if not provided
  if (is.null(disease_ref)) {
    disease_ref <- read_disease_snp_reference()
  }
  required_ref_cols <- c("rsid", "chromosome", "position", "gene", "disease", "risk_allele")
  if (!is.data.frame(disease_ref) || !all(required_ref_cols %in% colnames(disease_ref))) {
    warning("Invalid 'disease_ref' data frame provided or loaded. Cannot perform check.")
    return(tibble::tibble()) # Return empty tibble
  }
  # Ensure key ref columns are character for join
  disease_ref_for_join <- disease_ref %>%
    mutate(
      chromosome = as.character(.data$chromosome),
      position = as.character(.data$position),
      rsid = as.character(.data$rsid)
    )

  # Attempt merge - prioritize matching by RSID if available, fallback to Chrom/Pos
  # We need RSID from user_data here.
  merged_data <- tryCatch({
    # Prioritize join by RSID first as it's more stable
    # Ensure RSID is character in both
    user_data_rsid <- user_data_for_join %>% mutate(rsid = as.character(.data$rsid))
    disease_ref_rsid <- disease_ref_for_join %>% mutate(rsid = as.character(.data$rsid))

    # Join by RSID
    joined_rsid <- inner_join(
      disease_ref_rsid,
      user_data_rsid %>% select(all_of(required_user_cols)), # Keep all user cols
      by = "rsid",
      suffix = c(".ref", ".user") # Distinguish cols if needed
    )


    joined_rsid # Return the RSID-matched data

  }, error = function(e) {
    warning("Error during merging of disease SNP reference and user data: ", e$message)
    return(tibble::tibble()) # Return empty tibble on error
  })

  if (nrow(merged_data) == 0) {
    # Return empty tibble with expected output columns if no matches
    return(tibble::tibble(
      rsid = character(),
      chromosome = character(),
      position = character(),
      gene = character(),
      disease = character(),
      your_genotype = character(),
      risk_allele = character(),
      risk_allele_count = integer(),
      notes = character()
    ))
  }

  # Calculate risk allele count for matched SNPs
  results <- merged_data %>%
    # Use the user's genotype column (might be genotype.user if suffix added, check join)
    mutate(
      your_genotype = as.character(.data$genotype), # Use the user's genotype
      genotype_clean = ifelse(nchar(.data$your_genotype) == 2 & grepl("^[ACGTDI]{2}$", .data$your_genotype), .data$your_genotype, "--"),
      allele1 = str_sub(.data$genotype_clean, 1, 1),
      allele2 = str_sub(.data$genotype_clean, 2, 2),
      risk_allele = as.character(.data$risk_allele),
      risk_allele_count = case_when(
        is.na(.data$risk_allele) | .data$risk_allele == "" ~ NA_integer_,
        .data$genotype_clean == "--" ~ NA_integer_,
        TRUE ~ as.integer((.data$allele1 == .data$risk_allele) + (.data$allele2 == .data$risk_allele))
      )
    ) %>%
    # Select and arrange final columns
    select(
      .data$rsid,
      chromosome = .data$chromosome.ref, # Use chromosome/position from reference
      position = .data$position.ref,
      .data$gene,
      .data$disease,
      .data$your_genotype,
      .data$risk_allele,
      .data$risk_allele_count,
      .data$notes
    ) %>%
    # Arrange results for clarity
    arrange(.data$chromosome, as.numeric(.data$position))

  return(results)
}


# -------- End of R/disease_snps.R --------
