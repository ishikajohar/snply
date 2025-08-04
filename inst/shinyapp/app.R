options(shiny.maxRequestSize = 50 * 1024^2) # Sets the maximum file upload size to 50 MB

library(stringr)
library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
#library(snply) # Assuming 'snply' is installed and available, as it's used extensively
library(readr)
library(plotly)
library(data.table)
library(leaflet)

suppressPackageStartupMessages(library(snply))



# global setup
# at the top, replace your pop_cols with a list:
pop_cols <- list(
  CEU = "#4363d8",
  YRI = "#e41a1c",
  CHB = "#4daf4a",
  JPT = "#984ea3",
  GIH = "#ff7f00"
)



# --- Load Browning 2018 Introgressed Data ---
browning_path <- system.file("extdata", "Browning2018_introgressed.tsv", package = "snply", mustWork = FALSE)
if (browning_path == "" || !file.exists(browning_path)) {
  browning_path <- "Browning2018_EUR_full.tsv"
}
browning_df <- read_tsv(browning_path, col_types = cols(
  chromosome = col_character(),
  position = col_double(),
  rsid = col_character(),
  archaic_allele = col_character(),
  reference_allele = col_character()
))

# --- Load Haplogroup References ---
# Assuming these are available from the snply package
refs <- snply::load_haplogroup_references()
mt_ref <- refs$mt_ref
y_clean <- refs$y_clean
haplo_map <- refs$haplo_map

aisnp_ref <- snply::load_aisnp_ref()

# --- UI Definition ---
# --- UI Definition ---
ui <- fluidPage(
  titlePanel("Welcome to Snply Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("snp_file", "Upload Genotype File", accept = c(".txt", ".tsv", ".csv")),
      helpText("Accepts 23andMe style raw data (rsid, chromosome, position, genotype)."),
      hr(),
      conditionalPanel(
        condition = "input.tabs == 'Haplogroup Map'",
        uiOutput("haplogroup_summary")
      ),
      hr(),
      uiOutput("status_ui"),
      helpText("Powered by the snply R package.")
    ),

    mainPanel(
      tabsetPanel(id = "tabs",

                  tabPanel("SNP Data Preview",
                           h4("Your Uploaded SNP Data"),
                           p("This table shows the SNPs found in your uploaded file. It includes the RSID (SNP identifier), chromosome, position (relative to the genome build used by your data provider), and your genotype (the two alleles you have at that position)."),
                           helpText("You can search the table using the search box (e.g., by RSID like 'rs334' or position). Click column headers to sort."),
                           hr(),
                           DTOutput("snp_table")
                  ),

                  tabPanel("Heterozygosity Plot",
                           h4("Genome-Wide Heterozygosity"),
                           p("This plot summarises the proportion of your SNPs where the two alleles are different (heterozygous) versus the same (homozygous)."),
                           helpText("Heterozygous means you inherited a different version of the gene/marker from each parent. Homozygous means you inherited the same version from both parents. The overall heterozygosity rate is shown in the subtitle."),
                           hr(),
                           plotOutput("het_plot")
                  ),

                  tabPanel("Neanderthal Allele Analysis",
                           h4("Phenotypic archaic variants"),
                           p("This section highlights the most statistically significant Neanderthal-introgressed variants found in your genome that have been linked to modern human traits.

Genome-wide significant archaic variants (p < 1.0×10⁻⁸) influencing human traits"),
                           textOutput("pheno_summary"),
                           DTOutput("pheno_table"),
                           h4("Genotypic Archaic Variants"),
                           textOutput("intro_summary"),
                           helpText("The table shows your genotypes at SNP locations associated with archaic alleles in modern humans. These include loci associated with pigmentation, immunity, and metabolic traits."),
                           DTOutput("intro_table")
                  ),



                  tabPanel("Chromosome Painting (Neanderthal)",
                           h4("Chromosome Painting of Neanderthal Alleles"),
                           helpText("This visualisation shows the distribution of Neanderthal alleles across your chromosomes."),
                           helpText("Red blocks indicate regions where both alleles match the archaic version (homozygous), and green blocks indicate heterozygous sites (one archaic allele)."),
                           plotlyOutput("paint_plot", height = "900px")
                  ),

                  tabPanel("Chromosome Painting",
                           sidebarLayout(
                             sidebarPanel(
                               selectInput("paint_chr",  "Chromosome", choices = c(1:22, "X"), selected = "1"),
                               sliderInput("paint_thresh","Mono-ancestry threshold", min = .5, max = 1, value = .8, step = .05),
                               numericInput("paint_block","SNPs per block", value = 50, min = 10, max = 200)
                             ),
                             mainPanel(
                               # ← NEW DESCRIPTION BLOCK
                               h4("Genome-wide Monte Carlo Ancestry Painting"),
                               p("This plot shows your genome as a series of chromosomes, coloured by their most likely ancestry in each SNP block."),
                               tags$ul(
                                 tags$li(strong("Blocks of color:"), " Most-likely ancestry per block (CEU, YRI, CHB, JPT, GIH)."),
                                 tags$li(strong("Black segments:"), " Centromeres (regions of no SNP data)."),
                                 tags$li(strong("Grey segments:"), " Heterochromatic regions (typically excluded from analysis due to low SNP density)."),
                                 tags$li(strong("Interpretation:"), " Long continuous blocks indicate extended ancestry tracts; frequent breaks may indicate recombination or mixed ancestry.")
                               ),
                               hr(),  # optional separator
                               # ← EXISTING PLOTS
                               fluidRow(
                                 column(6, plotlyOutput("paint_ancestry", height = "300px")),
                                 column(6, plotOutput  ("single_chr_gg",  height = "300px"))
                               ),
                               hr(),
                               fluidRow(
                                 column(12, plotOutput("genomePaintingPlot", height = "900px", width = "700px"))
                               )
                             )
                           )
                  ),


                  tabPanel("PCA Projection",
                           h4("Projected PCA Plot"),
                           p("The PCA (Principal Component Analysis) plot shows how your genetic data clusters relative to global reference populations."),
                           p("It visually places you within a space of known population groups from Africa, Europe, East Asia, South Asia, and America. Your closest group is calculated and displayed."),
                           plotOutput("pca_plot"),
                           hr(),
                           h4("Subpopulations by Region"),
                           DTOutput("pca_subpop_table"),
                           hr(),
                           strong(textOutput("closest_region_text"))
                  ),

                  tabPanel("Common Disease SNPs",
                           h4("Genotypes for Common Disease-Associated SNPs"),
                           p("This table checks your genotype data against a curated list of SNPs associated with common diseases or traits based on published genome-wide association studies (GWAS)."),
                           helpText("Only SNPs that match between your uploaded file and the curated reference list are shown."),
                           strong("Disclaimer: This information is for educational purposes only and not intended as medical advice or diagnosis."),
                           hr(),
                           DTOutput("disease_table")
                  ),

                  tabPanel("Haplogroup Map",
                           h4("Maternal and Paternal Haplogroup Origins"),
                           p("Your maternal and paternal haplogroups trace your direct ancestry lines: mitochondria (from mother to all children) and Y chromosome (from father to son)."),
                           p("These lineages reflect ancient human migrations and can show the deep ancestral roots of your family line."),
                           leafletOutput("maternal_map", height = 300),
                           hr(),
                           leafletOutput("paternal_map", height = 300)
                  ),

                  tabPanel("Disclaimer",
                           h4("Important Information"),
                           p("This Snply application and R package are provided for educational and informational purposes only."),
                           p("The analyses performed (including heterozygosity, Neanderthal allele estimation, chromosome painting, and common disease SNP checks) are based on publicly available data and standard methods, but may contain limitations."),
                           strong("This tool is NOT intended for medical diagnosis, health assessment, or clinical decision-making."),
                           p("Data Privacy: When you upload a file, it is processed locally within your R session or browser. No data is uploaded or stored externally.")
                  )
      )
    )
  )
)
# --- Server Logic ---
server <- function(input, output, session) {
  user_df <- reactive({
    req(input$snp_file)
    snply::read_user_data(input$snp_file$datapath)
  })

  # Run paint_mc() *once* for all panels
  blocks <- reactive({
    req(input$snp_file)
    snply::paint_mc(
      user_snps = user_df(),
      ref_freq  = aisnp_ref,
      block_n   = input$paint_block,
      iter      = 150,
      thresh    = input$paint_thresh
    )
  })

  plot_single_chr <- function(blocks, chrom, pop_cols) {
    dt <- melt(
      as.data.table(blocks)[CHR==chrom],
      id.vars       = c("CHR","start","end"),
      measure.vars  = c("hap1","hap2"),
      variable.name = "haplotype",
      value.name    = "ancestry"
    )
    dt[, haplotype := factor(haplotype, levels=c("hap2","hap1"), labels=c("Hap 2","Hap 1"))]

    ggplot(dt, aes(
      x        = start/1e6,
      xend     = end/1e6,
      y        = haplotype,
      yend     = haplotype,
      color    = ancestry
    )) +
      geom_segment(linewidth =5) +
      scale_color_manual(values=pop_cols) +
      scale_y_discrete(expand=expansion(add=.3)) +
      labs(x="Position (Mbp)", y=NULL, title=paste0("Chr ",chrom)) +
      theme_minimal(base_size=14) +
      theme(
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        legend.position="none"
      )
  }




  # Updated plotting function for full ancestry painting across all chromosomes

  plot_all_chromosomes <- function(paint_mc_result, pop_cols, snp_counts){
    library(data.table)
    library(ggplot2)

    # 1) Melt hap1/hap2 into long form (no more melt warnings)
    dt <- as.data.table(paint_mc_result)[, .(CHR, start, end, hap1, hap2)]
    paint_df <- melt(
      dt,
      id.vars       = c("CHR","start","end"),
      measure.vars  = c("hap1","hap2"),     # <- use a character vector
      variable.name = "haplotype",
      value.name    = "ancestry"
    )

    # 2) Chromosome factor & numeric index
    chrom_levels <- c(as.character(1:22),"X","Y")
    paint_df[, chromosome := factor(
      ifelse(CHR=="23","X",
             ifelse(CHR=="24","Y", CHR)),
      levels = chrom_levels
    )]
    paint_df[, chr_idx := as.numeric(chromosome)]

    # 3) Compute rectangle edges for each haplotype
    paint_df[, `:=`(
      x_min    = start/1e6,
      x_max    = end/1e6,
      y_bottom = chr_idx + ifelse(haplotype=="hap1", -0.2, +0.2),
      y_top    = chr_idx + ifelse(haplotype=="hap1",  0.0,  0.4)
    )]
    paint_df_no_XY <- paint_df[!(chromosome %in% c("X","Y"))]


    # 4) Fetch UCSC cytobands (centromeres & het) quietly
    cyto <- suppressMessages(fread(
      "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
      col.names = c("chr","start","end","band","stain"),
      showProgress = FALSE
    ))[, chr := sub("^chr","",chr)]
    cyto <- cyto[chr %in% chrom_levels]
    cyto[, chr_idx := as.numeric(factor(chr, levels=chrom_levels))]

    # ▸ keep 1–22, drop X / Y
    cen <- cyto[stain == "acen"            & !(chr %in% c("X","Y"))]
    het <- cyto[stain %in% c("gvar","stalk") & !(chr %in% c("X","Y"))]

    # 5) Merge in chromosome lengths so we can nudge SNP‐count labels
    chrom_sizes <- data.table(
      chromosome = factor(chrom_levels, levels = chrom_levels),
      length     = c(
        249250621,243199373,198022430,191154276,180915260,171115067,
        159138663,146364022,141213431,135534747,135006516,133851895,
        115169878,107349540,102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566,155270560, 59373566
      )
    )
    snp_counts <- merge(
      snp_counts,
      chrom_sizes[, .(chromosome,length)],
      by = "chromosome",
      all.x = TRUE
    )
    snp_counts[, chr_idx := as.numeric(chromosome)]

    # 6) Build the plot
    ggplot() +
      # heterochromatin background
      geom_rect(data = het,
                aes(xmin = start/1e6, xmax = end/1e6,
                    ymin = chr_idx - 0.2,     ymax = chr_idx + 0.6),
                fill = "grey90", color = NA) +
      # centromere
      geom_rect(data = cen,
                aes(xmin = start/1e6, xmax = end/1e6,
                    ymin = chr_idx - 0.2,     ymax = chr_idx + 0.6),
                fill = "black", color = NA) +
      # ancestry blocks (both haplotypes)
      geom_rect(data = paint_df_no_XY,
                aes(xmin = x_min,   xmax = x_max,
                    ymin = y_bottom, ymax = y_top,
                    fill = ancestry),
                color = NA) +
      # SNP‐count labels, nudged off the right end
      geom_text(data = snp_counts,
                aes(x = (length/1e6) + 5,
                    y = chr_idx + 0.2,
                    label = total_snps),
                size = 3) +
      # keep chr 1 at bottom → Y at top
      scale_y_continuous(
        breaks = seq_along(chrom_levels),
        labels = chrom_levels,
        expand = c(0.02,0)
      ) +
      scale_x_continuous(
        expand = c(0,0)
      ) +
      scale_fill_manual(
        name   = "Ancestry",
        values = unlist(pop_cols)
      ) +
      labs(
        title = "Monte Carlo Chromosome Painting",
        x     = NULL,    # no x‐axis title
        y     = NULL     # no y‐axis title
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.major.y  = element_blank(),
        panel.grid.minor    = element_blank(),
        axis.ticks.y        = element_blank(),
        axis.title.x        = element_blank(),
        axis.title.y        = element_blank(),
        legend.position     = "right"
    )
}






  # 2) static ggplot single-chr
  output$single_chr_gg <- renderPlot({
    req(blocks())
    plot_single_chr(blocks(), input$paint_chr, pop_cols)
  })

  output$genomePaintingPlot <- renderPlot({
    req(blocks())

    # 1) Grab the wide output of paint_mc()
    wide_dt <- data.table::as.data.table(blocks())
    wide_dt[, `:=`(start_bp = start, end_bp = end)]


    # 2) Compute SNP‐counts per chromosome
    #    – .N gives the number of rows (blocks) per CHR
    snp_counts <- wide_dt[, .N, by = .(CHR)]
    #    – rename columns so they match plot_all_chromosomes()
    snp_counts[, chromosome := CHR]
    data.table::setnames(snp_counts, "N", "total_snps")

    # 3) Call your existing plot_all_chromosomes()
    plot_all_chromosomes(
      paint_mc_result = wide_dt,  # wide: CHR, start, end, hap1, hap2
      pop_cols        = pop_cols,  # your color map
      snp_counts      = snp_counts # data.table with chromosome & total_snps
    )
  })




  intro_joined <- reactive({
    u <- user_df(); req(u)
    # Ensure chromosome, position, and rsid are character types for consistent joins
    brow <- browning_df |> mutate(across(c(chromosome, position, rsid), as.character))
    usr <- u |> mutate(across(c(chromosome, position, rsid), as.character)) |> select(chromosome, position, rsid, genotype)
    d <- inner_join(brow, usr, by = c("chromosome", "position", "rsid"))
    if (nrow(d) > 0) {
      d$archaic_count <- mapply(function(geno, arch) {
        if (is.na(geno) || grepl("-", geno)) return(NA_integer_)
        if (nchar(geno) == 1) ifelse(geno == arch, 1L, 0L) else sum(strsplit(geno, "")[[1]] == arch)
      }, d$genotype, d$archaic_allele)
    }
    d
  })

  output$status_ui <- renderUI({
    if (is.null(input$snp_file)) tags$p(strong("Status:"), "Waiting for upload …")
    else tags$p(strong("Status:"), "File loaded ✔")
  })

  output$snp_table <- renderDT({
    df <- user_df(); req(df)
    if ("rsid" %in% names(df)) {
      # Creates clickable links for RSIDs to NCBI dbSNP
      df$rsid <- sprintf('<a href="https://www.ncbi.nlm.nih.gov/snp/%s" target="_blank">%s</a>', df$rsid, df$rsid)
    }
    datatable(df, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$het_plot <- renderPlot({
    df <- user_df(); req(df)
    # Filter for genotypes that are exactly two characters long (e.g., "AA", "AG")
    valid <- df$genotype[nchar(df$genotype) == 2]
    het <- sum(substr(valid, 1, 1) != substr(valid, 2, 2)) # Heterozygous: alleles are different
    hom <- sum(substr(valid, 1, 1) == substr(valid, 2, 2)) # Homozygous: alleles are the same
    total <- het + hom

    df_plot <- data.frame(
      type = c("Heterozygous", "Homozygous"),
      count = c(het, hom),
      pct = round(c(het, hom) / total * 100, 1)
    )

    ggplot(df_plot, aes(x = "Genotype", y = pct, fill = type)) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_text(aes(label = paste0(pct, "%")), position = position_stack(vjust = 0.5), size = 5) +
      scale_fill_manual(values = c("Heterozygous" = "#66c2a5", "Homozygous" = "#fc8d62")) +
      labs(
        title = "Proportion of Heterozygous vs Homozygous SNPs",
        x = NULL, y = "Percentage",
        subtitle = sprintf("Overall heterozygosity = %.2f%%", 100 * het / total)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank()
      )
  })

  calc_pheno <- reactive({
    req(input$snp_file)
    snply::calculate_neanderthal_alleles(input$snp_file$datapath)
  })

  #── build detailed stats for the phenotype‐SNPs ───────────────────────────────
  pheno_stats <- reactive({
    req(input$snp_file)

    # (a) full reference of GW‐significant archaic SNPs
    ref <- snply::read_reference_data() %>% # Assuming read_reference_data from snply
      mutate(across(c(chromosome, position), as.character))

    # (b) your uploaded SNPs
    usr <- user_df() %>%
      mutate(across(c(chromosome, position), as.character))

    # (c) merge and compute counts
    merged <- inner_join(ref, usr, by = c("chromosome","position")) %>%
      mutate(
        genotype = ifelse(nchar(genotype)==2, genotype, NA_character_),
        allele1 = str_sub(genotype, 1, 1),
        allele2 = str_sub(genotype, 2, 2),
        archaic_count = (allele1 == archaic_allele) + (allele2 == archaic_allele)
      )

    hom <- sum(merged$archaic_count == 2, na.rm=TRUE)
    het <- sum(merged$archaic_count == 1, na.rm=TRUE)
    total_hits <- hom + het

    total_ref <- nrow(ref)
    matched <- nrow(merged)
    total_user <- nrow(usr)
    pct_matched <- if(matched > 0) 100 * total_hits / matched else 0
    pct_of_upload <- if(total_user > 0) 100 * total_hits / total_user else 0

    list(
      total_ref = total_ref,
      matched = matched,
      hom = hom,
      het = het,
      total_hits = total_hits,
      pct_matched = pct_matched,
      pct_of_upload = pct_of_upload,
      total_user = total_user
    )
  })

  #── render the clean summary line ─────────────────────────────────────────────
  output$pheno_summary <- renderText({
    s <- pheno_stats()
    sprintf(
      "Out of %d genome‐wide significant archaic SNPs in the reference, %d were found in your data: %d homozygous and %d heterozygous, totalling %d sites with archaic alleles (%.1f%% of matched sites). That represents %.3f%% of all %d SNPs in your file.",
      s$total_ref, s$matched, s$hom, s$het, s$total_hits, s$pct_matched, s$pct_of_upload, s$total_user
    )
  })

  output$neo_summary <- renderText(calc_pheno()$summary)

  output$pheno_table <- renderDT({
    tbl <- calc_pheno()$phenotypes
    keep <- intersect(c("rsid", "chromosome", "position", "genotype", "archaic_allele", "reference_allele",
                        "Hair.colour", "Skin.colour", "Weight", "Standing.height", "Ease.of.skin.tanning"), names(tbl))
    datatable(tbl[, keep, drop = FALSE], options = list(pageLength = 10, scrollX = TRUE))
  })

  output$intro_summary <- renderText({
    d <- intro_joined(); req(d)
    total_ref <- nrow(browning_df)
    comparable <- nrow(d)
    hom <- sum(d$archaic_count == 2, na.rm = TRUE)
    het <- sum(d$archaic_count == 1, na.rm = TRUE)
    total_arch <- hom + het

    # 1) % of matched reference sites
    pct_matched <- if (comparable > 0) 100 * total_arch / comparable else NA_real_

    # 2) % of all user SNPs
    total_user_snps <- nrow(user_df())
    pct_all_snps <- if (total_user_snps > 0) 100 * total_arch / total_user_snps else NA_real_

    sprintf(
      "Out of %d introgressed SNPs in the reference, %d were found in your data: %d homozygous and %d heterozygous, totalling %d sites with archaic alleles (%.1f%% of matched sites).
      That represents %.3f%% of all %d SNPs in your file.
      **Note:** this is *not* your genome-wide Neanderthal ancestry (typically ~1–4%% of the genome).",
      total_ref, comparable, hom, het, total_arch, pct_matched,
      pct_all_snps, total_user_snps
    )
  })

  output$intro_table <- renderDT({
    d <- intro_joined(); req(d)
    if (nrow(d) == 0) return(datatable(data.frame(Message = "No introgressed SNPs found."), options = list(dom = 't')))
    show <- intersect(c("chromosome", "position", "rsid", "genotype", "archaic_allele", "reference_allele", "archaic_count"), names(d))
    datatable(d[, show, drop = FALSE], options = list(pageLength = 10, scrollX = TRUE))
  })

  output$disease_table <- renderDT({
    # Assuming check_disease_snps is a function from snply or a local helper
    datatable(snply::check_disease_snps(user_df()), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$paint_plot <- renderPlotly({
    base <- draw_chromosome_silhouette()
    g    <- paint_snps(base, intro_joined())
    plotly::ggplotly(g, tooltip=c("chr","position","archaic_count"))
  })

  # ==== PCA Result Reactive ====
  pca_result <- reactive({
    req(input$snp_file)
    tryCatch({
      snply::perform_pca_projection(
        user_file = input$snp_file$datapath,
        data_dir = system.file("extdata", package = "snply"),
        scores_file = "pca_scores_labeled_codedregions.tsv.gz"  # must include region_group
      )
    }, error = function(e) {
      showNotification(paste("PCA projection failed:", e$message), type = "error", duration = 10)
      NULL
    })
  })

  # ==== PCA Plot Output ====
  output$pca_plot <- renderPlot({
    req(pca_result())
    pca_result()$plot
  })

  # ==== Subpopulation Table ====
  output$pca_subpop_table <- renderDT({
    req(pca_result())
    datatable(pca_result()$subpop_table, rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE))
  })

  # ==== Optional: Closest Region Text ====
  output$closest_region_text <- renderText({
    req(pca_result())
    paste("Closest major population:", pca_result()$closest_region)
  })
  maternal_info <- reactive({
    df <- user_df()
    user_mt <- df[df$chromosome %in% c("MT", "chrMT", "26", "chr26"), ]
    # Assuming snply_detect_maternal_haplogroup from snply
    snply::snply_detect_maternal_haplogroup(user_mt, mt_ref)
  })

  paternal_info <- reactive({
    df <- user_df()
    user_y <- df[df$chromosome %in% c("Y", "chrY", "24", "chr24"), ]
    # Assuming snply_detect_paternal_haplogroup from snply
    snply::snply_detect_paternal_haplogroup(user_y, y_clean)
  })

  output$maternal_map <- renderLeaflet({
    # Assuming plot_haplogroup_location from snply
    snply::plot_haplogroup_location(maternal_info()$haplogroup, haplo_map)
  })

  output$paternal_map <- renderLeaflet({
    # Assuming plot_haplogroup_location from snply
    snply::plot_haplogroup_location(paternal_info()$display_haplogroup, haplo_map)
  })

  output$haplogroup_summary <- renderUI({
    mat <- maternal_info()
    pat <- paternal_info()
    mat_hap <- mat$haplogroup %||% "N/A"
    mat_hits <- mat$match_count %||% NA
    mat_conf <- snply::confidence_level(mat_hits) # Assuming confidence_level from snply
    pat_hap <- pat$display_haplogroup %||% "N/A"
    pat_hits <- pat$match_count %||% NA
    pat_conf <- snply::confidence_level(pat_hits) # Assuming confidence_level from snply
    confidence_color <- function(level) {
      switch(level, "High" = "green", "Medium" = "orange", "Low" = "red", "None" = "gray", "Unknown" = "gray")
    }
    div(
      style = "border: 1px solid #ccc; border-radius: 10px; padding: 15px; background-color: #f9f9f9;",
      tags$h4("\U0001F50D Haplogroup Summary", style = "margin-top: 0;"),
      tags$p(strong("Maternal Haplogroup: "), mat_hap, br(), strong("SNP Matches: "), mat_hits, br(),
             span(strong("Confidence: "), style = paste0("color:", confidence_color(mat_conf)), mat_conf)),
      tags$hr(),
      tags$p(strong("Paternal Haplogroup: "), pat_hap, br(), strong("SNP Matches: "), pat_hits, br(),
             span(strong("Confidence: "), style = paste0("color:", confidence_color(pat_conf)), pat_conf))
    )
  })

}

shinyApp(ui, server)
