
ref_df <- read_reference_data()
head(ref_df)
user_df <- read_user_data("C:/Users/30694/Downloads/snply/user_snps.txt")
head(user_df)
merged_df <- merge_snp_data(ref_df, user_df)
head(merged_df)
result <- calculate_neanderthal_alleles("C:/Users/30694/Downloads/snply/user_snps.txt")
print(result$summary)
head(result$phenotypes)


library(snply)
user_file <- system.file("extdata", "example_user_snps.txt", package = "snply")
result <- calculate_neanderthal_alleles(user_file)
print(result$summary)
head(result$phenotypes)
