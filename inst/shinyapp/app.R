# Set max upload size to 50MB
options(shiny.maxRequestSize = 50 * 1024^2)

# Load necessary libraries
library(shiny)
library(DT)
library(ggplot2)
library(readr) # Still needed if read_user_data is used directly elsewhere, but snply loads its dependencies
library(dplyr) # Still needed if used directly elsewhere, but snply loads its dependencies
library(stringr) # Still needed if used directly elsewhere, but snply loads its dependencies

# Load your installed package (ensure it's rebuilt with the updated code)
# Make sure 'snply' is installed and contains the updated functions
library(snply)

ui <- fluidPage(
  titlePanel("Welcome to Snply analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("snp_file", "Upload Genotype File",
                accept = c(".txt", ".csv", ".tsv", ".vcf" # Note: .vcf might need specific parsing not included here
                )),
      helpText("Accepts whitespace-delimited formats (like 23andMe) with columns: rsid, chromosome, position, genotype. Header can optionally start with '#'. Max size: 50MB."),
      hr(),
      helpText("Powered by the snply R package.") # Optional branding
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("SNP Data Preview",
                 DTOutput("snp_table")),
        tabPanel("Heterozygosity Plot",
                 plotOutput("heterozygosity_plot")),
        tabPanel("Neanderthal Allele Analysis",
                 h4("Summary"), # Add a header for clarity
                 textOutput("alleleCount"),
                 hr(), # Add separator
                 h4("Phenotype-Associated Variants"), # Add a header
                 helpText("Table shows your genotypes at archaic SNP locations linked to certain phenotypes (if any found carrying the archaic allele)."),
                 DTOutput("allelePhenotypes"))
      )
    )
  )
)

server <- function(input, output, session) { # Added session for potential future use

  # Reactive expression to load user SNP data using the package function.
  snp_data_reactive <- reactive({
    req(input$snp_file) # Ensure a file is uploaded

    # Use tryCatch directly around the read function for UI feedback
    tryCatch({
      # Using the function from the loaded package
      snply::read_user_data(input$snp_file$datapath)
    }, error = function(e) {
      # Show notification to the user
      showNotification(paste("Error reading file:", e$message), type = "error", duration = 10)
      # Return NULL or an empty frame with an error column to prevent downstream crashes
      # Returning NULL might be better as subsequent reactives can use req()
      return(NULL)
      # Alternatively: return(data.frame(Error = paste("File Reading Error:", e$message)))
    })
  })

  # Render the user SNP data table preview.
  output$snp_table <- renderDT({
    user_data <- snp_data_reactive()
    req(user_data) # Require valid data frame

    # Check if it's the error data frame (if you chose that option above)
    # if ("Error" %in% colnames(user_data) && nrow(user_data) == 1) {
    #   datatable(user_data) # Display the error in the table
    # } else {
    datatable(user_data, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    # }
  })

  # Render the heterozygosity plot.
  output$heterozygosity_plot <- renderPlot({
    data <- snp_data_reactive()
    req(data) # Require successful data loading

    # Check if 'genotype' column exists and is not empty
    if ("genotype" %in% colnames(data) && nrow(data) > 0) {
      # Filter for valid 2-character genotypes before counting
      valid_genotypes <- data$genotype[nchar(data$genotype) == 2 & grepl("^[ACGTDI]+$", data$genotype)] # Allow A,C,G,T,D,I

      if(length(valid_genotypes) == 0) {
        # Plot indicating no valid genotypes found
        ggplot() + labs(title = "No valid genotypes (e.g., 'AA', 'CG') found for analysis", x = "", y = "") + theme_void()
      } else {
        heterozygous_snps <- sum(substr(valid_genotypes, 1, 1) != substr(valid_genotypes, 2, 2))
        homozygous_snps <- sum(substr(valid_genotypes, 1, 1) == substr(valid_genotypes, 2, 2))
        total_valid <- length(valid_genotypes) # Use count of valid genotypes

        # Avoid division by zero if no valid SNPs
        heterozygosity_rate <- if (total_valid > 0) heterozygous_snps / total_valid else 0

        plot_data <- data.frame(category = c("Heterozygous", "Homozygous"),
                                count = c(heterozygous_snps, homozygous_snps))

        ggplot(plot_data, aes(x = category, y = count, fill = category)) +
          geom_bar(stat = "identity") +
          geom_text(aes(label = count), vjust = -0.3) + # Add count labels
          theme_minimal(base_size = 14) +
          labs(title = paste("Genotype Summary (considering", total_valid, "valid SNPs)"),
               subtitle = paste("Overall Heterozygosity:", round(heterozygosity_rate, 4)),
               y = "Count", x = "Genotype Type") +
          scale_fill_manual(values = c("Heterozygous" = "lightblue", "Homozygous" = "lightcoral")) +
          theme(legend.position = "none") # Hide redundant legend
      }
    } else {
      # Plot indicating missing data or read error
      ggplot() +
        labs(title = "Genotype data missing or file could not be read", x = "", y = "") +
        theme_void()
    }
  })

  # Reactive expression for the Neanderthal allele calculation using the package function.
  allele_result_reactive <- reactive({
    req(input$snp_file) # Need a file path

    # The calculate_neanderthal_alleles function now has internal tryCatch
    snply::calculate_neanderthal_alleles(input$snp_file$datapath)
  })

  # Render the Neanderthal summary text. Handles error message from allele_result.
  output$alleleCount <- renderText({
    result <- allele_result_reactive()
    req(result) # Require the result list exists

    # The summary element will contain either the results or an error message
    result$summary
  })

  # Render the phenotype table. Handles error case from allele_result.
  output$allelePhenotypes <- renderDT({
    result <- allele_result_reactive()
    req(result) # Require the result list exists

    # Check if the 'phenotypes' element indicates an error
    # (based on the structure returned by the modified calculate_neanderthal_alleles)
    pheno_data <- result$phenotypes
    if (is.data.frame(pheno_data) && "Error" %in% colnames(pheno_data)) {
      # Display the error message in the table area
      datatable(pheno_data, options = list(dom = 't'), rownames = FALSE) # Minimal table display for error
    } else if (is.data.frame(pheno_data) && nrow(pheno_data) > 0) {
      # Display the actual phenotype data
      datatable(pheno_data, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    } else {
      # Display a message if no variants found or other non-error empty case
      datatable(data.frame(Message = "No Neanderthal variants associated with phenotypes found in your data."),
                options = list(dom = 't'), rownames = FALSE)
    }
  })

  # Note: allelePercentage output was removed as percentage is now in the summary text
  # output$allelePercentage <- renderText({ "" }) # Keep if needed for other purposes

}

# Run the application
shinyApp(ui = ui, server = server)
