# R/browning_loader.R  (inside your package)
#'@export
read_reference_browning <- function() {
  path <- system.file("extdata", "Browning2018_introgressed.tsv",
                      package = "snply", mustWork = FALSE)
  if (path == "" || !file.exists(path)) {
    path <- file.path("..","inst","extdata","Browning2018_introgressed.tsv")
  }
  readr::read_tsv(
    path,
    col_names = c("chromosome","position","rsid",
                  "archaic_allele","reference_allele"),
    col_types = readr::cols(
      chromosome        = readr::col_character(),
      position          = readr::col_character(),
      rsid              = readr::col_character(),
      archaic_allele    = readr::col_character(),
      reference_allele  = readr::col_character()
    )
  )
}
