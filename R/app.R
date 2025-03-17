#' Launch the Shiny App
#'
#' @export
run_app <- function() {
  library(shiny)

  # Define UI
  ui <- fluidPage(
    titlePanel("Welcome to snply"),
    sidebarLayout(
      sidebarPanel("Sidebar content here"),
      mainPanel("Main panel content here")
    )
  )

  # Define Server
  server <- function(input, output, session) {}

  # Return the Shiny app object
  shinyApp(ui = ui, server = server)
}


