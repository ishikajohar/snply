calculate_heterozygosity_and_plot <- function(snp_data) {
  # Extract genotype data (assuming genotypes are in the "genotype" column)
  genotype_data <- snp_data$genotype

  # Create counters for heterozygous and homozygous SNPs
  heterozygous_snps <- 0
  homozygous_snps <- 0

  # Loop through all genotype values in the genotype column
  for (genotype in genotype_data) {
    if (!is.na(genotype)) {  # Ignore missing values
      # Check if genotype is heterozygous (two different alleles)
      if (substr(genotype, 1, 1) != substr(genotype, 2, 2)) {
        heterozygous_snps <- heterozygous_snps + 1  # Heterozygous (e.g., AG, TC)
      } else {
        homozygous_snps <- homozygous_snps + 1  # Homozygous (e.g., AA, GG)
      }
    }
  }

  # Compute heterozygosity (avoid division by zero)
  total_snps <- heterozygous_snps + homozygous_snps
  heterozygosity <- ifelse(total_snps > 0, heterozygous_snps / total_snps, NA)

  # Create a dataframe for the bar plot
  plot_data <- data.frame(
    category = c("Heterozygous", "Homozygous"),
    count = c(heterozygous_snps, homozygous_snps)
  )

  # Plot using ggplot2
  library(ggplot2)
  ggplot(plot_data, aes(x = category, y = count, fill = category)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste("Heterozygosity vs Homozygosity\n(Heterozygosity = ", round(heterozygosity, 3), ")"),
         y = "Count", x = "Genotype Type") +
    scale_fill_manual(values = c("lightblue", "#FFB6C1"))  # Light blue & Light pink

  return(heterozygosity)
}








