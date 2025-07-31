

# ==== Calculate and Plot Heterozygosity ====
calculate_heterozygosity_and_plot <- function(user_data) {
  # Clean genotype column
  user_data <- user_data[
    !is.na(genotype) &
      nchar(genotype) == 2 &
      genotype != "--" &
      genotype != "NC"
  ]

  genotype_data <- user_data$genotype

  # Count heterozygous and homozygous SNPs
  hetero <- sum(substr(genotype_data, 1, 1) != substr(genotype_data, 2, 2))
  homo   <- sum(substr(genotype_data, 1, 1) == substr(genotype_data, 2, 2))
  total  <- hetero + homo

  if (total == 0) {
    warning("No valid SNPs found.")
    return(list(value = NA, plot = NULL))
  }

  # Calculate percentages and heterozygosity
  percent_hetero <- (hetero / total) * 100
  percent_homo   <- (homo / total) * 100
  heterozygosity <- hetero / total

  # Create data for stacked bar plot
  plot_data <- data.frame(
    Genotype = c("Heterozygous", "Homozygous"),
    Percent = c(percent_hetero, percent_homo),
    Sample = "Your Sample"
  )

  # Create ggplot with percentage labels
  hetero_plot <- ggplot(plot_data, aes(x = Sample, y = Percent, fill = Genotype)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(round(Percent, 1), "%")),
              position = position_stack(vjust = 0.5), size = 5, color = "black") +
    labs(
      title = paste0("Heterozygosity = ", round(heterozygosity * 100, 2), "%"),
      y = "Percentage (%)", x = ""
    ) +
    scale_fill_manual(values = c("lightblue", "#FFB6C1")) +
    theme_minimal() +
    theme(legend.position = "top")

  # Return both value and plot (for Shiny integration)
  return(list(
    value = heterozygosity,
    plot = hetero_plot
  ))
}
