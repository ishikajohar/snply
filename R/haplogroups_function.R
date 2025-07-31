#' Detect Maternal Haplogroup
#'
#' @param user_mt Data frame of mitochondrial SNPs
#' @param mt_ref Reference data for mtDNA haplogroups
#' @return A list with best match haplogroup and match count
#' @export
snply_detect_maternal_haplogroup <- function(user_mt, mt_ref) {
  best_match <- NULL
  max_hits <- 0
  user_mt$position <- as.integer(user_mt$position)

  for (i in seq_len(nrow(mt_ref))) {
    hap <- mt_ref$haplogroup[i]
    muts <- unlist(strsplit(mt_ref$defining_mutations[i], "\\s+|,"))
    muts <- gsub("[()!]", "", muts)

    hits <- sum(sapply(muts, function(m) {
      if (grepl("^[ACGT][0-9]+[ACGTacgt]$", m)) {
        alt <- substr(m, nchar(m), 1)
        pos <- as.numeric(gsub("[^0-9]", "", m))
        match <- user_mt[user_mt$position == pos, ]
        return(nrow(match) > 0 && grepl(alt, match$genotype[1], ignore.case = TRUE))
      }
      return(FALSE)
    }))

    if (hits > max_hits) {
      max_hits <- hits
      best_match <- hap
    }
  }

  return(list(haplogroup = best_match, match_count = max_hits))
}

#' Detect Paternal Haplogroup by Y Position
#'
#' @param user_y Data frame of Y chromosome SNPs
#' @param y_ref Cleaned reference of Y-SNPs
#' @param top_n Number of top matches to return
#' @return A list with haplogroup, match count, and cleaned display name
#' @export
snply_detect_paternal_haplogroup <- function(user_y, y_ref, top_n = 1) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  match_counts <- list()
  user_y$position <- as.integer(user_y$position)

  for (i in seq_len(nrow(y_ref))) {
    pos <- y_ref$position[i]
    hap <- y_ref$haplogroup[i]
    mut <- y_ref$mutation[i]

    if (is.na(pos) || !grepl("->", mut)) next

    derived <- sub(".*->", "", mut)
    match_row <- user_y[user_y$position == pos, ]

    if (nrow(match_row) > 0) {
      geno <- match_row$genotype[1]
      if (geno != "--" && grepl(derived, geno, ignore.case = TRUE)) {
        match_counts[[hap]] <- match_counts[[hap]] %||% 0
        match_counts[[hap]] <- match_counts[[hap]] + 1
      } else if (geno == "--") {
        match_counts[[hap]] <- match_counts[[hap]] %||% 0
        match_counts[[hap]] <- match_counts[[hap]] + 1
      }
    }
  }

  sorted <- sort(unlist(match_counts), decreasing = TRUE)
  best_hap <- names(head(sorted, top_n))
  hits <- unname(head(sorted, top_n))

  if (length(best_hap) == 0 || length(hits) == 0) {
    return(list(haplogroup = NA, match_count = 0, display_haplogroup = "N/A"))
  }

  cleaned <- gsub("\\(.*\\)", "", best_hap[1])
  cleaned <- gsub(" or .*", "", cleaned)
  cleaned <- trimws(cleaned)
  display <- substr(cleaned, 1, 4)

  return(list(haplogroup = best_hap[1], match_count = hits[1], display_haplogroup = display))
}

#' Plot Haplogroup Dispersion Map
#'
#' @param haplogroup Haplogroup string
#' @param haplo_map Reference map of haplogroup locations
#' @return A leaflet map object
#' @export
plot_haplogroup_location <- function(haplogroup, haplo_map) {
  library(leaflet)
  clean <- function(x) gsub("[^A-Za-z0-9]", "", x)

  if (is.null(haplogroup) || is.na(haplogroup) || haplogroup == "") {
    return(leaflet() %>% addTiles() %>% addPopups(0, 0, "No haplogroup detected"))
  }

  haplogroup_clean <- clean(haplogroup)
  haplo_map$cleaned <- clean(haplo_map$haplogroup)

  hg_match <- haplo_map[startsWith(haplo_map$cleaned, haplogroup_clean) |
                          startsWith(haplogroup_clean, haplo_map$cleaned), ]

  if (nrow(hg_match) == 0 && nchar(haplogroup_clean) >= 1) {
    fallback <- substr(haplogroup_clean, 1, 1)
    hg_match <- haplo_map[startsWith(haplo_map$cleaned, fallback), ]
  }

  if (nrow(hg_match) == 0) {
    return(leaflet() %>% addTiles() %>% addPopups(0, 0, paste("Location unknown for:", haplogroup)))
  }

  hg_ordered <- hg_match[order(hg_match$lon), ]
  start <- hg_ordered[1, ]
  end <- hg_ordered[nrow(hg_ordered), ]
  middle <- hg_ordered[2:(nrow(hg_ordered) - 1), ]

  leaflet() %>%
    addProviderTiles("CartoDB.Positron") %>%
    addCircleMarkers(lng = start$lon, lat = start$lat, radius = 6, color = "green",
                     popup = paste("Origin:", start$label)) %>%
    addCircleMarkers(lng = end$lon, lat = end$lat, radius = 6, color = "red",
                     popup = paste("Destination:", end$label)) %>%
    addCircleMarkers(data = middle, ~lon, ~lat, radius = 5,
                     color = "#2c3e50", fillOpacity = 0.7, popup = ~label) %>%
    addPolylines(lng = hg_ordered$lon, lat = hg_ordered$lat,
                 color = "gray", opacity = 0.6, weight = 2, dashArray = "6,6") %>%
    addLegend(position = "bottomright",
              colors = c("green", "red", "#2c3e50"),
              labels = c("Start", "End", "Intermediate"), title = "Migration Key")
}

#' Confidence Level from Match Count
#'
#' @param n SNP match count
#' @return "High", "Medium", "Low", "None", or "Unknown"
#' @export
confidence_level <- function(n) {
  if (is.null(n) || is.na(n)) return("Unknown")
  if (n >= 10) return("High")
  if (n >= 5) return("Medium")
  if (n >= 1) return("Low")
  return("None")
}

#' Load Haplogroup Reference Data
#'
#' @return A list with mt_ref, y_clean, haplo_map
#' @export
load_haplogroup_references <- function() {
  # Attempt to use installed package paths
  mt_ref_path <- system.file("extdata", "mtDNA_haplogroup_reference.csv", package = "snply", mustWork = FALSE)
  y_snp_path  <- system.file("extdata", "Expanded_SNP_Index_Human.csv", package = "snply", mustWork = FALSE)
  map_path    <- system.file("extdata", "haplogroup_dispersions_locations.rds", package = "snply", mustWork = FALSE)

  # Fallback to development file paths if the above fail
  if (!file.exists(mt_ref_path)) mt_ref_path <- "../../extdata/mtDNA_haplogroup_reference.csv"
  if (!file.exists(y_snp_path))  y_snp_path  <- "../../extdata/Expanded_SNP_Index_Human.csv"
  if (!file.exists(map_path))    map_path    <- "../../extdata/haplogroup_dispersions_locations.rds"

  # Final check to fail loudly if files are still missing
  if (!file.exists(mt_ref_path)) stop("Missing file: mtDNA_haplogroup_reference.csv")
  if (!file.exists(y_snp_path))  stop("Missing file: Expanded_SNP_Index_Human.csv")
  if (!file.exists(map_path))    stop("Missing file: haplogroup_dispersions_locations.rds")


  # Load and preprocess files
  mt_ref <- read.csv(mt_ref_path, stringsAsFactors = FALSE)

  y_raw <- read.csv(y_snp_path, stringsAsFactors = FALSE)
  y_clean <- y_raw %>%
    dplyr::filter(
      grepl("^rs", snp_id) | grepl("^[0-9]+$", snp_id),
      !is.na(mutation),
      !is.na(snp_id)
    ) %>%
    dplyr::mutate(position = as.integer(position))

  haplo_map <- readRDS(map_path)

  return(list(mt_ref = mt_ref, y_clean = y_clean, haplo_map = haplo_map))
}
