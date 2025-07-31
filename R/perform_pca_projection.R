#' Project a 23andMe Genotype File onto Global PCA Space
#'
#' Projects a user's genotype file onto a precomputed PCA space (e.g. from 1000 Genomes).
#' Returns the PCA plot and a subpopulation table grouped by region.
#'
#' @param user_file Path to user-uploaded SNP file (23andMe-style format).
#' @param data_dir Path to reference PCA directory (default = snply extdata).
#' @param scores_file Name of PCA scores file (default = "pca_scores_labeled_codedregions.tsv.gz").
#'
#' @return A list with:
#'   \item{plot}{ggplot2 PCA plot with user projection}
#'   \item{subpop_table}{DataFrame of region_group and region}
#'   \item{closest_region}{User's closest major population group}
#' @export
perform_pca_projection <- function(user_file,
                                   data_dir = NULL,
                                   scores_file = "pca_scores_labeled_codedregions.tsv.gz") {
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)

  # ==== 1. Locate reference files ====
  if (is.null(data_dir) || data_dir == "" || !file.exists(file.path(data_dir, scores_file))) {
    data_dir <- system.file("extdata", package = "snply")
    if (data_dir == "" || !file.exists(file.path(data_dir, scores_file))) {
      stop("Could not find PCA reference files.")
    }
  }

  # ==== 2. Load reference PCA data ====
  scores <- fread(file.path(data_dir, scores_file))
  loadings <- fread(file.path(data_dir, "pca_loadings.tsv.gz"))
  means <- scan(file.path(data_dir, "pca_means.txt"))
  snps_used_raw <- readLines(file.path(data_dir, "pca_snps_used.tsv"))
  snps_df <- data.table(
    snp_raw = snps_used_raw,
    rsid = sub("_[ACGT]+$", "", snps_used_raw),
    alt = sub("^.*_([ACGT]+)$", "\\1", snps_used_raw)
  )

  # ==== 3. Load and clean user file ====
  user <- fread(user_file, skip = 1, header = FALSE,
                col.names = c("rsid", "chrom", "pos", "genotype"))
  user <- user[!is.na(genotype) & genotype != "--"]
  user[, genotype := toupper(genotype)]

  # ==== 4. Merge and convert to 0/1/2 based on ALT ====
  merged <- merge(snps_df, user, by = "rsid")
  call_012 <- function(gt, alt) {
    if (nchar(gt) != 2 || grepl("N", gt)) return(NA)
    a1 <- substr(gt, 1, 1)
    a2 <- substr(gt, 2, 2)
    if (a1 == alt && a2 == alt) return(2)
    if (a1 == alt || a2 == alt) return(1)
    if (a1 != alt && a2 != alt) return(0)
    return(NA)
  }
  merged[, dosage := mapply(call_012, genotype, alt)]

  # ==== 5. Align and impute ====
  geno_vec <- merged[match(snps_used_raw, snp_raw), dosage]
  names(geno_vec) <- snps_used_raw
  geno_vec[is.na(geno_vec)] <- means[is.na(geno_vec)]
  user_centered <- as.numeric(geno_vec) - means
  user_coords <- matrix(user_centered, nrow = 1) %*% as.matrix(loadings[, -1, with = FALSE])

  # ==== 6. Closest region ====
  if (!"region_group" %in% names(scores)) stop("Scores file must include 'region_group'")
  region_centroids <- scores %>%
    group_by(region_group) %>%
    summarise(centroid_PC1 = mean(PC1), centroid_PC2 = mean(PC2), .groups = "drop") %>%
    mutate(distance = sqrt((centroid_PC1 - user_coords[1])^2 + (centroid_PC2 - user_coords[2])^2))
  closest_region <- region_centroids$region_group[which.min(region_centroids$distance)]

  # ==== 7. Colors ====
  region_colors <- c(
    "Africa"     = "#e41a1c",
    "America"    = "#ff7f00",
    "East Asia"  = "#377eb8",
    "South Asia" = "#4daf4a",
    "Europe"     = "#984ea3"
  )

  # ==== 8. Plot ====
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = region_group)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_point(data = data.frame(PC1 = user_coords[1], PC2 = user_coords[2]),
               aes(x = PC1, y = PC2), inherit.aes = FALSE,
               color = "black", size = 4, shape = 8) +
    annotate("text", x = user_coords[1], y = user_coords[2], label = "You",
             size = 5, fontface = "bold", vjust = -1) +
    scale_color_manual(values = region_colors, drop = FALSE) +
    labs(title = "PCA Projection of User",
         subtitle = paste("Closest Region:", closest_region),
         x = "PC1", y = "PC2", color = "Region") +
    theme_minimal(base_size = 14)

  # ==== 9. Subpop table ====
  subpop_table <- scores %>% count(region_group, region, name = "sample_count") %>% arrange(region_group)

  # ==== 10. Return ====
  return(list(plot = p, subpop_table = subpop_table, closest_region = closest_region))
}
