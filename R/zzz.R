.onLoad <- function(libname, pkgname) {
	if (requireNamespace("reticulate", quietly = TRUE)) {
		library(reticulate)
		envname <- "r-MCNV2"
		if (virtualenv_exists(envname)) {
			use_virtualenv(envname, required = TRUE)
		} else {
			packageStartupMessage(
				"MCNV2: Python virtualenv 'r-MCNV2' not found. Run `MCNV2::setup_python_env()` first."
			)
		}
	}
}
