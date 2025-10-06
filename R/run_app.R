#' @title Launch MCNV2 Shiny app
#' @param bedtools_path Full path to bedtools executable
#' @param results_dir Directory to save results
#' @param python_env Name of the Python virtualenv (default 'r-MCNV2')
#' @export
launch <- function(bedtools_path = NULL, results_dir = NULL, python_env = "r-MCNV2") {
	if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' required.")
	library(reticulate)
	
	# Active le virtualenv
	if (!virtualenv_exists(python_env)) {
		stop("Python environment '", python_env, "' not found. Run `MCNV2::setup_python_env()` first.")
	}
	use_virtualenv(python_env, required = TRUE)
	
	# Chemin de l'app Shiny
	appDir <- system.file("shiny-app", package = "MCNV2")
	if (appDir == "") stop("Could not find shiny-app directory.")
	
	# Passe les paramètres via options()
	options(MCNV2.params = list(
		bedtools_path = bedtools_path,
		results_dir = results_dir
	))
	
	# Lance l'app
	shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
