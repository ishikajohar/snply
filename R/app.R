# Set max upload size to 50MB
options(shiny.maxRequestSize = 50 * 1024^2)

# Load necessary libraries
library(shiny)
library(DT)  # For rendering data tables
library(ggplot2)  # For plotting

# Load helper functions
source("load_snp.R")        # Load function to read SNP files
source("heterozygosity.R")  # Load function to calculate heterozygosity

# Define UI
ui <- fluidPage(
  titlePanel("Welcome to Snply!"),
  sidebarLayout(
    sidebarPanel(
      fileInput("snp_file", "Upload SNP File", accept = c(".txt", ".csv", ".tsv", ".vcf")),
      helpText("Accepted formats: .txt, .csv, .tsv, .vcf | Max size: 50MB."),
      hr()
    ),
    mainPanel(
      DTOutput("snp_table"),
      hr(), plotOutput("heterozygosity_plot")
    )
  )
)

# Define Server
server <- function(input, output) {
  snp_data <- reactive({
    req(input$snp_file)
    tryCatch({
      load_snp_data(input$snp_file$datapath)
    }, error = function(e) {
      data.frame(Error = "Invalid file format. Please upload a valid SNP file.")
    })
  })

  output$snp_table <- renderDT({
    datatable(snp_data(), options = list(pageLength = 10))
  })

  output$heterozygosity_plot <- renderPlot({
    data <- snp_data()
    if ("genotype" %in% colnames(data)) {
      heterozygous_snps <- sum(grepl("[ACGT]{2}", data$genotype) &
                                 substr(data$genotype, 1, 1) != substr(data$genotype, 2, 2))
      homozygous_snps <- sum(grepl("[ACGT]{2}", data$genotype) &
                               substr(data$genotype, 1, 1) == substr(data$genotype, 2, 2))
      plot_data <- data.frame(category = c("Heterozygous", "Homozygous"), count = c(heterozygous_snps, homozygous_snps))
      ggplot(plot_data, aes(x = category, y = count, fill = category)) +
        geom_bar(stat = "identity") + theme_minimal() +
        labs(title = paste("Heterozygosity vs Homozygosity\n(Heterozygosity = ", round(heterozygous_snps / (heterozygous_snps + homozygous_snps), 4), ")"),
             y = "Count", x = "Genotype Type") +
        scale_fill_manual(values = c("lightblue", "lightpink"))
    } else {
      ggplot() + labs(title = "Error: Genotype data missing", x = "", y = "") + theme_void()
    }
  })
}

# Run the app
shinyApp(ui, server)

