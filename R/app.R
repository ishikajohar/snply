# Set max upload size to 50MB
options(shiny.maxRequestSize = 50 * 1024^2)

# Load necessary libraries
library(shiny)
library(DT)        # For rendering data tables
library(ggplot2)   # For plotting
library(readxl)    # For reading Excel files (used in calculate_neanderthal_alleles)

# Source helper functions (if they're in separate files)
# source("load_snp.R")
# source("heterozygosity.R")
# source("neanderthal.R")

# Define UI
ui <- fluidPage(
  titlePanel("Welcome to Snply!"),
  sidebarLayout(
    sidebarPanel(
      fileInput("snp_file", "Upload SNP File",
                accept = c(".txt", ".csv", ".tsv", ".vcf")),
      helpText("Accepted formats: .txt, .csv, .tsv, .vcf | Max size: 50MB."),
      hr()
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("SNP Data",
                 DTOutput("snp_table")),
        tabPanel("Heterozygosity Plot",
                 plotOutput("heterozygosity_plot")),
        tabPanel("Neanderthal Allele Calculation",
                 textOutput("alleleCount"),
                 textOutput("allelePercentage"),
                 DTOutput("allelePhenotypes"))
      )
    )
  )
)

# Define Server
server <- function(input, output) {

  # Reactive expression to load SNP data for heterozygosity and table
  snp_data <- reactive({
    req(input$snp_file)
    tryCatch({
      load_snp_data(input$snp_file$datapath)
    }, error = function(e) {
      data.frame(Error = "Invalid file format. Please upload a valid SNP file.")
    })
  })

  # Render the SNP data table
  output$snp_table <- renderDT({
    datatable(snp_data(), options = list(pageLength = 10))
  })

  # Render the heterozygosity plot
  output$heterozygosity_plot <- renderPlot({
    data <- snp_data()
    if ("genotype" %in% colnames(data)) {
      heterozygous_snps <- sum(grepl("[ACGT]{2}", data$genotype) &
                                 substr(data$genotype, 1, 1) != substr(data$genotype, 2, 2))
      homozygous_snps <- sum(grepl("[ACGT]{2}", data$genotype) &
                               substr(data$genotype, 1, 1) == substr(data$genotype, 2, 2))
      plot_data <- data.frame(category = c("Heterozygous", "Homozygous"),
                              count = c(heterozygous_snps, homozygous_snps))
      ggplot(plot_data, aes(x = category, y = count, fill = category)) +
        geom_bar(stat = "identity") + theme_minimal() +
        labs(title = paste("Heterozygosity vs Homozygosity\n(Heterozygosity = ",
                           round(heterozygous_snps / (heterozygous_snps + homozygous_snps), 4), ")"),
             y = "Count", x = "Genotype Type") +
        scale_fill_manual(values = c("lightblue", "lightpink"))
    } else {
      ggplot() + labs(title = "Error: Genotype data missing", x = "", y = "") + theme_void()
    }
  })

  # Reactive expression for the Neanderthal calculation
  allele_result <- reactive({
    req(input$snp_file)  # Make sure the file is uploaded
    calculate_neanderthal_alleles(input$snp_file$datapath,
                                  "C:/Users/30694/Downloads/Data/neander_snps.rds.xlsx")
  })

  # Render the Neanderthal allele count as text
  output$alleleCount <- renderText({
    req(allele_result())
    paste("Total Neanderthal allele copies:", allele_result()$count)
  })

  # Render the Neanderthal allele percentage as text
  output$allelePercentage <- renderText({
    req(allele_result())
    paste("Percentage of Neanderthal-derived alleles:", round(allele_result()$percentage, 2), "%")
  })

  # Render the phenotype table using DT
  output$allelePhenotypes <- renderDT({
    req(allele_result())
    datatable(allele_result()$phenotypes, options = list(pageLength = 10))
  })

}

# Run the app
shinyApp(ui, server)

