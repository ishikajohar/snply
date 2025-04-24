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



### Shiny App  
```Shiny App

library(snply)

# Launch the Shiny app
launchSnplyApp()
```

## Data Sources  
- Reference SNPs: `inst/extdata/neander_snps.csv`
