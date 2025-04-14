# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(snply)


# Install your package
devtools::install("snply")

# Run analysis
library(snply)
result <- calculate_neanderthal_alleles("their_file.txt")

shiny::runApp("path/to/your/shiny/app.R")

system.file("extdata", "neander_snps.csv", package = "snply")
