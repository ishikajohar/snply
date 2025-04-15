Developing a new R package (provisionally called snply) that provides a series of genomic analyses through a Shiny web interface.
The main purpose is to let users upload their 23andMe SNP files and then run several analyses such as calculating heterozygosity,
estimating the percentage of Neanderthal alleles, plotting PCA positions against global reference data, and “painting” chromosomes by ancestry.
These functions are inspired by academic projects like the Interpretome, which demonstrates how genomic data can be analyzed in a private, customizable way.

# Snply: Neanderthal Allele Analysis  

## Installation  
```r
# installation
devtools::install_github(
  "ishikajohar/snply",
  ref = "feature/neanderthal",
)
```

## Usage  
### R Functions  
```r
library(snply)
result <- calculate_neanderthal_alleles("path/to/23andme.txt")
print(result$summary)
```

### Shiny App  
```r
shiny::runApp(system.file("shinyapps", "app.R", package = "snply"))
```

## Data Sources  
- Reference SNPs: `inst/extdata/neander_snps.csv`
