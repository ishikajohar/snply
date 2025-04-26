# -------- inst/shinyapp/app.R --------

# Set max upload size to 50MB
options(shiny.maxRequestSize = 50 * 1024^2)

# Load necessary libraries
library(shiny)
library(DT)
library(ggplot2)
library(dplyr) # Keep dplyr loaded for potential direct use in server logic
library(stringr)
library(snply) # Load our package functions

# Define UI
ui <- fluidPage(
  titlePanel("Welcome to Snply Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("snp_file", "Upload Genotype File",
                accept = c(".txt", ".csv", ".tsv", ".vcf")),
      helpText("Accepts whitespace-delimited formats (like 23andMe) with columns: rsid, chromosome, position, genotype. Header can optionally start with '#'. Max size: 50MB."),
      hr(),
      helpText("Powered by the snply R package.")
    ),
    mainPanel(
      tabsetPanel(
        # Tab 1: SNP Data Preview
        tabPanel("SNP Data Preview",
                 h4("Your Uploaded SNP Data"),
                 p("This table shows the SNPs found in your uploaded file. It includes the RSID (SNP identifier), chromosome, position (relative to the genome build used by your data provider), and your genotype (the two alleles you have at that position)."),
                 helpText("You can search the table using the search box (e.g., by RSID like 'rs334' or position). Click column headers to sort."),
                 hr(),
                 DTOutput("snp_table")
        ),
        # Tab 2: Heterozygosity Plot
        tabPanel("Heterozygosity Plot",
                 h4("Genome-Wide Heterozygosity"),
                 p("This plot summarizes the proportion of your SNPs where the two alleles are different (heterozygous) versus the same (homozygous)."),
                 helpText("Heterozygous means you inherited a different version of the gene/marker from each parent (e.g., genotype 'AG'). Homozygous means you inherited the same version from both parents (e.g., 'AA' or 'GG'). The overall heterozygosity rate is shown in the subtitle."),
                 hr(),
                 plotOutput("heterozygosity_plot")
        ),
        # Tab 3: Neanderthal Allele Analysis
        tabPanel("Neanderthal Allele Analysis",
                 # Existing content...
                 h4("Summary"),
                 textOutput("alleleCount"),
                 hr(),
                 h4("Phenotype-Associated Variants"),
                 helpText("Table shows your genotypes at archaic SNP locations linked to certain phenotypes in the reference data (if any were found carrying the archaic allele)."),
                 DTOutput("allelePhenotypes")
        ),
        # Tab 4: Common Disease SNPs (NEW)
        tabPanel("Common Disease SNPs",
                 h4("Genotypes for Common Disease-Associated SNPs"),
                 p("This table checks your genotype data against a curated list of SNPs commonly associated with specific diseases or traits."),
                 helpText("Results are shown only for SNPs present in both the reference list and your uploaded file. 'Risk Allele Count' indicates how many copies of the reference risk allele you have (0, 1, or 2)."),
                 strong("Disclaimer: This information is for educational purposes only and is not medical advice or diagnosis."),
                 hr(),
                 DTOutput("disease_snp_table") # New table output
        ),
        # Tab 5: Disclaimer (NEW)
        tabPanel("Disclaimer",
                 h4("Important Information"),
                 p("This Snply application and R package are provided for educational and informational purposes only."),
                 p("The analyses performed (including heterozygosity, Neanderthal allele estimation, and common disease SNP checks) are based on publicly available data and standard methods, but have limitations and should not be considered definitive or exhaustive."),
                 strong("This tool is NOT intended for medical diagnosis, health assessment, or decision-making."),
                 p("Results should not be used to replace consultation with qualified healthcare professionals. Genetic associations are complex, and interpretations require expert knowledge."),
                 p("Data Privacy: When you upload a file, it is processed locally within your R session or web browser. Your genetic data is NOT uploaded to or stored on any external server by this application."),
                 p("Accuracy: While efforts are made to use reliable reference data, genetic databases and annotations change. Positional information depends on the genome build (e.g., GRCh37/hg19) used by your data provider. RSIDs (SNP identifiers) are generally more stable."),
                 p("We are not liable for any decisions made based on the information provided by this tool. Use at your own risk.")
        )
      ) # end tabsetPanel
    ) # end mainPanel
  ) # end sidebarLayout
) # end fluidPage

# Define Server Logic
server <- function(input, output, session) {

  # Reactive expression to load user SNP data. NOW INCLUDES RSID.
  snp_data_reactive <- reactive({
    req(input$snp_file)
    tryCatch({
      # read_user_data now returns rsid, chromosome, position, genotype
      snply::read_user_data(input$snp_file$datapath)
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error", duration = 10)
      return(NULL)
    })
  })

  # Render the user SNP data table preview. INCLUDES RSID and optional link.
  output$snp_table <- renderDT({
    user_data_raw <- snp_data_reactive() # Use a temporary name
    req(user_data_raw) # Require valid data frame

    # --- Make RSID clickable ---
    # Directly modify the 'rsid' column to contain the HTML link string
    user_data_linked <- user_data_raw %>%
      mutate(
        rsid = paste0('<a href="https://www.ncbi.nlm.nih.gov/snp/', .data$rsid, '" target="_blank">', .data$rsid, '</a>')
        # Note: Using .data$rsid inside mutate is good practice
      )
    # No select() needed here just for links, mutate handles it.
    # --- End ---

    # Pass the modified data to datatable.
    # escape = FALSE tells DT to render the HTML in the rsid column.
    datatable(
      user_data_linked, # Use the data frame with the modified rsid column
      escape = FALSE,   # Ensure this is active for links to work
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  # Render the heterozygosity plot (remains the same)
  output$heterozygosity_plot <- renderPlot({
    data <- snp_data_reactive()
    req(data)
    if ("genotype" %in% colnames(data) && nrow(data) > 0) {
      valid_genotypes <- data$genotype[nchar(data$genotype) == 2 & grepl("^[ACGTDI]+$", data$genotype)]
      if(length(valid_genotypes) == 0) {
        ggplot() + labs(title = "No valid genotypes (e.g., 'AA', 'CG') found for analysis", x = "", y = "") + theme_void()
      } else {
        heterozygous_snps <- sum(substr(valid_genotypes, 1, 1) != substr(valid_genotypes, 2, 2))
        homozygous_snps <- sum(substr(valid_genotypes, 1, 1) == substr(valid_genotypes, 2, 2))
        total_valid <- length(valid_genotypes)
        heterozygosity_rate <- if (total_valid > 0) heterozygous_snps / total_valid else 0
        plot_data <- data.frame(category = c("Heterozygous", "Homozygous"),
                                count = c(heterozygous_snps, homozygous_snps))
        ggplot(plot_data, aes(x = category, y = count, fill = category)) +
          geom_bar(stat = "identity") +
          geom_text(aes(label = count), vjust = -0.3) +
          theme_minimal(base_size = 14) +
          labs(title = paste("Genotype Summary (considering", total_valid, "valid SNPs)"),
               subtitle = paste("Overall Heterozygosity:", round(heterozygosity_rate, 4)),
               y = "Count", x = "Genotype Type") +
          scale_fill_manual(values = c("Heterozygous" = "lightblue", "Homozygous" = "lightcoral")) +
          theme(legend.position = "none")
      }
    } else {
      ggplot() +
        labs(title = "Genotype data missing or file could not be read", x = "", y = "") +
        theme_void()
    }
  })

  # Neanderthal Allele Calculation (remains the same)
  allele_result_reactive <- reactive({
    req(input$snp_file)
    # Need the raw file path for calculate_neanderthal_alleles
    # (it calls read_user_data internally)
    snply::calculate_neanderthal_alleles(input$snp_file$datapath)
  })
  output$alleleCount <- renderText({
    result <- allele_result_reactive()
    req(result)
    result$summary
  })
  output$allelePhenotypes <- renderDT({
    result <- allele_result_reactive()
    req(result)
    pheno_data <- result$phenotypes
    if (is.data.frame(pheno_data) && "Error" %in% colnames(pheno_data)) {
      datatable(pheno_data, options = list(dom = 't'), rownames = FALSE)
    } else if (is.data.frame(pheno_data) && nrow(pheno_data) > 0) {
      datatable(pheno_data, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    } else {
      datatable(data.frame(Message = "No Neanderthal variants associated with phenotypes found in your data."),
                options = list(dom = 't'), rownames = FALSE)
    }
  })


  # --- NEW: Disease SNP Analysis ---
  # Reactive expression for disease SNP calculation
  disease_snp_results_reactive <- reactive({
    user_data <- snp_data_reactive() # Get the loaded user data (now includes rsid)
    req(user_data) # Requires valid user data

    # Call the checking function from the package
    # It will load the default reference data internally
    snply::check_disease_snps(user_data)
  })

  # Render the disease SNP results table
  output$disease_snp_table <- renderDT({
    results <- disease_snp_results_reactive()
    req(results) # Require results exist

    # Check for empty results (different from error)
    if (is.data.frame(results) && nrow(results) > 0) {
      datatable(results,
                rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE)
      )
    } else {
      # Display a message if no relevant SNPs were found in user data
      datatable(data.frame(Message = "No matching common disease-associated SNPs found in your data or reference file."),
                options = list(dom = 't'), rownames = FALSE)
    }
  })
  # --- End NEW ---

}

# Run the application (This line should NOT be here if launching via launchSnplyApp)
shinyApp(ui = ui, server = server) # Remove or comment out this line

# -------- End of inst/shinyapp/app.R --------
