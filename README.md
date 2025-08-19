# Snply <img src="man/figures/logo.png" align="right" width="120"/>

Snply is an R package (in development) that provides a series of genomic analyses through a **Shiny web interface**.  
The main purpose is to let users upload their **23andMe SNP files** and then run several analyses, such as:

- Calculating heterozygosity  
- Estimating the percentage of Neanderthal alleles  
- Plotting PCA positions against global reference data  
- “Painting” chromosomes by ancestry  

These functions are inspired by academic projects like **Interpretome**, which demonstrated how genomic data can be analysed in a **private, customisable, and educational** way.

---

## Features

- **Heterozygosity analysis** – summary statistics and stacked bar plots of homozygous vs heterozygous SNPs  
- **Neanderthal allele detection** – estimate the percentage of Neanderthal-derived SNPs present in user data  
- **PCA projection** – project 23andMe genotypes into a precomputed PCA space of global reference populations  
- **Chromosome painting** – ancestry assignment at the chromosome level  
- **Interactive Shiny app** – a clean web-based interface to run all analyses without requiring coding expertise  

---

##   Data Sources

Snply includes reference datasets packaged in `inst/extdata/` for reproducibility:

- **Reference SNPs for Neanderthal analysis**  
  - `neander_snps.csv`

- **PCA reference files**  
  - `pca_loadings.tsv.gz`  
  - `pca_scores_labeled_codedregions.tsv.gz`  
  - `pca_means.txt`  
  - `pca_snps_used.tsv`

- **Haplogroup lookup tables**  
  - `mtDNA_haplogroup_reference.csv`  
  - `Expanded_SNP_Index_Human.csv`

- **Map data**  
  - `haplogroup_dispersions_locations.rds`

All files are accessed internally using `system.file()`, ensuring reproducibility and independence from the user’s working directory.

---

##  Validation

- Tested with 23andMe v4 and v5 genotype files  
- PCA projection correctly captures broad continental population structure  
- **Limitations**: reduced SNP marker set (~3,992 SNPs) may lower fine-scale resolution, occasionally leading to misclassification  

---

## License

This project is released under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.

## Installation and Running the Shiny App

Snply can be installed directly from GitHub and then launched:

```r
# install devtools if not already installed
install.packages("devtools")

# install Snply
devtools::install_github(
  "ishikajohar/snply",
  ref = "master"
)

# load the package
library(snply)

# launch the app
snply::launchSnplyApp()

# --- If developing locally instead ---

# load the package
library(snply)

# build documentation and reinstall
devtools::document()
devtools::install()

# launch the app
snply::launchSnplyApp()
