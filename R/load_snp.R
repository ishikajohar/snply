# Function to load SNP data (supports .txt, .csv, .tsv, and .vcf)
load_snp_data <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))  # Get file extension

  if (ext %in% c("txt", "csv", "tsv")) {
    # Load SNP data from tab-delimited file
    snp_data <- read.table(file_path, comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(snp_data) <- c("rsid", "chromosome", "position", "genotype")

  } else if (ext == "vcf") {
    # Load SNP data from VCF file
    vcf_data <- readLines(file_path)  # Read the file line-by-line
    vcf_data <- vcf_data[!grepl("^##", vcf_data)]  # Remove metadata lines (## lines)

    # Convert VCF data to a dataframe
    snp_data <- read.table(text = vcf_data, comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

    # Extract only relevant SNP columns (CHROM, POS, ID, REF/ALT, and GENOTYPE)
    colnames(snp_data) <- c("chromosome", "position", "rsid", "ref", "alt", "qual", "filter", "info", "format", "genotype")
    snp_data <- snp_data[, c("rsid", "chromosome", "position", "genotype")]

  } else {
    stop("Unsupported file format. Please upload a .txt, .csv, .tsv, or .vcf file.")
  }

  return(snp_data)
}

# TEST
# Define file path
#snp_file <- "file path"

# Load SNP data
#snp_data <- load_snp_data(snp_file)

# View first few rows
#head(snp_data)
