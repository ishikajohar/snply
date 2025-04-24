#' Launch the Snply Shiny Application
#'
#' This function launches the Shiny application included with the snply package.
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom rlang check_installed
launchSnplyApp <- function() {
  # Check if shiny is installed before trying to load it via ::
  rlang::check_installed("shiny", reason = "to run the Snply Shiny application.")

  appDir <- system.file("shinyapp", package = "snply")
  if (appDir == "") {
    stop(
      "Could not find the 'shinyapp' directory in the snply package.\n",
      "Try re-installing the package.",
      call. = FALSE
    )
  }
  shiny::runApp(appDir, display.mode = "normal")
}
