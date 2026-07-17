# ==============================================================================
# run.R  --  Launch the CRISPR Oligo Designer app
# ==============================================================================
# Usage:
#   Rscript run.R
# or from an R session in this folder:
#   source("run.R")
# ==============================================================================

# Resolve this script's own directory so it works from any working directory.
.args <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args[grep("^--file=", .args)])
app_dir <- if (length(.file_arg) > 0) {
  dirname(normalizePath(.file_arg))
} else if (!is.null(sys.frames()[[1]]$ofile)) {
  dirname(normalizePath(sys.frames()[[1]]$ofile))
} else {
  getwd()
}

setwd(app_dir)
cat("Launching CRISPR Oligo Designer from:", app_dir, "\n")

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The 'shiny' package is not installed. Run  Rscript install.R  first.", call. = FALSE)
}

shiny::runApp(app_dir, launch.browser = TRUE)
