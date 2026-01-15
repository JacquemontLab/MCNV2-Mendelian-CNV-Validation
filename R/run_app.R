#' @title Launch MCNV2 Shiny app
#' @param bedtools_path Full path to bedtools executable
#' @param results_dir Directory to save results
#' @param python_env Name of the Python virtualenv (default 'r-MCNV2')
#' @export
launch <- function(bedtools_path = NULL, results_dir = NULL, python_env = "r-MCNV2") {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required.")
  }
  library(reticulate)

  # ---- Check & activate Python env ----
  if (!virtualenv_exists(python_env)) {
    stop(
      "Python environment '", python_env,
      "' not found. Run `MCNV2::setup_python_env()` first."
    )
  }
  use_virtualenv(python_env, required = TRUE)

  # ---- Normalize & check bedtools path ----
  if (is.null(bedtools_path) || bedtools_path == "") {
    bedtools_path <- Sys.which("bedtools")
  }

  if (bedtools_path == "" || !file.exists(bedtools_path)) {
    stop(
      "bedtools not found. Please install bedtools and/or provide bedtools_path.\n",
      "Example: MCNV2::launch(bedtools_path = '/usr/local/bin/bedtools')"
    )
  }

  # ---- Normalize & create results directory ----
  if (is.null(results_dir) || results_dir == "") {
    stop("Please provide a results_dir.")
  }

  results_dir <- path.expand(results_dir)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(results_dir)) {
    stop("results_dir does not exist and could not be created: ", results_dir)
  }

  # ---- Locate Shiny app ----
  appDir <- system.file("shiny-app", package = "MCNV2")
  if (appDir == "") stop("Could not find shiny-app directory.")

  # ---- Pass parameters to Shiny ----
  options(MCNV2.params = list(
    bedtools_path = bedtools_path,
    results_dir   = results_dir
  ))

  # ---- Launch app ----
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}

