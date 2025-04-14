
# Path to your .bim file
bim_file <- "C:/Users/30694/Downloads/snply/plink_data/1000G_phase3_autosomes_23andme_pruned.bim"

# Read the .bim file
bim <- read.table(bim_file, stringsAsFactors = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")

# Replace missing SNP IDs (.) with "chr:pos"
bim$snp[bim$snp == "."] <- paste0("chr", bim$chr[bim$snp == "."], ":", bim$pos[bim$snp == "."])

# Check for duplicates (e.g., same chr:pos but different alleles)
duplicates <- duplicated(bim$snp) | duplicated(bim$snp, fromLast = TRUE)
if (sum(duplicates) > 0) {
  # Append "_1", "_2" to duplicates
  bim$snp <- make.unique(bim$snp, sep = "_")
}

# Save the corrected .bim file
write.table(bim, bim_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

















# Install and load snpStats if not already installed
library(snpStats)

# Define the full path prefix to your PLINK files (without extension)
plink_prefix <- "C:/Users/30694/Downloads/snply/plink_data/1000G_phase3_autosomes_23andme_pruned"

# Read the PLINK data using the prefix; read.plink will look for .bed, .bim, and .fam files
plink_data <- read.plink(plink_prefix)

# Convert genotype data to a numeric matrix (0, 1, 2)
geno_matrix <- as(plink_data$genotypes, "numeric")

# Read the .bim file to get SNP positions
bim_file <- "C:/Users/30694/Downloads/snply/plink_data/1000G_phase3_autosomes_23andme_pruned.bim"
bim <- read.table(bim_file, stringsAsFactors = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")
