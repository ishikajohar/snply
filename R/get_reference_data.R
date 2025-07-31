#' Get Reference Data for Chromosome Painting
#'
#' Checks for local cached copies of the pre-calculated reference allele
#' frequencies and SNP map. If not found, downloads them from specified URLs
#' (e.g., from a GitHub Release) and saves them to the cache directory.
#' Loads and returns the data.
#'
#' @param freq_rds_url The direct URL for the reference frequency RDS file.
#'   Default points at v0.1.0 release.
#' @param map_rds_url The direct URL for the reference SNP map RDS file.
#'   Default points at v0.1.0 release.
#' @param package_name Name of the package (used for cache directory). Default: "snply".
#' @param force_download Logical; if TRUE, re-download even if cached.
#'   Default: FALSE.
#' @return A list with elements \code{ref_freq_data} and \code{ref_snp_map}, or
#'   NULL if loading fails.
#' @export
#' @importFrom tools R_user_dir
#' @importFrom utils download.file
get_reference_data <- function(
    freq_rds_url = "https://github.com/ishikajohar/snply/releases/download/v0.1.0/reference_pop_freqs.rds",
    map_rds_url  = "https://github.com/ishikajohar/snply/releases/download/v0.1.0/ref_snp_map.rds",
    package_name = "snply",
    force_download = FALSE
) {
  # validate URLs
  if (!grepl("^https?://", freq_rds_url) || !grepl("^https?://", map_rds_url)) {
    stop("Invalid URL(s); must start with http:// or https://", call. = FALSE)
  }

  # cache directory
  cache_dir <- tools::R_user_dir(package_name, which = "cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # local paths
  local_freq  <- file.path(cache_dir, "reference_pop_freqs.rds")
  local_map   <- file.path(cache_dir, "ref_snp_map.rds")

  # helper to download if missing or forced
  download_if_needed <- function(url, dest, desc) {
    if (!file.exists(dest) || force_download) {
      message("Downloading ", desc, "...")
      tryCatch({
        utils::download.file(url, dest, mode = "wb", timeout = 600)
        if (!file.exists(dest) || file.info(dest)$size < 1000) {
          stop("Downloaded file looks too small or is missing.")
        }
        message(desc, " downloaded.")
        TRUE
      }, error = function(e) {
        warning("Failed to download ", desc, ": ", e$message, call. = FALSE)
        if (file.exists(dest)) file.remove(dest)
        FALSE
      })
    } else {
      message(desc, " found in cache.")
      TRUE
    }
  }

  ok1 <- download_if_needed(freq_rds_url, local_freq, "Reference frequencies")
  ok2 <- download_if_needed(map_rds_url,  local_map,  "Reference SNP map")

  # load
  ref_freq <- NULL
  ref_map  <- NULL
  if (ok1) {
    tryCatch({ ref_freq <- readRDS(local_freq) },
             error = function(e) warning("Failed to read freq RDS:", e$message, call. = FALSE))
  }
  if (ok2) {
    tryCatch({ ref_map <- readRDS(local_map) },
             error = function(e) warning("Failed to read map RDS:", e$message, call. = FALSE))
  }

  if (is.null(ref_freq) || is.null(ref_map)) {
    warning("Could not load one or both reference data files.", call. = FALSE)
    return(NULL)
  }

  message("Reference data loaded successfully.")
  list(ref_freq_data = ref_freq, ref_snp_map = ref_map)
}
