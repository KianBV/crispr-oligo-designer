# ==============================================================================
# packages.R  --  Central package loading for the backend pipeline
# ==============================================================================
#
# The backend assumes its dependencies are already installed (run install.R
# once after cloning). This file only *loads* them and fails with a clear,
# actionable message if something is missing, so non-programmers know exactly
# what to do.
# ==============================================================================

# CRAN packages required by the backend + GUI
.CRAN_PACKAGES <- c("dplyr", "stringr", "data.table", "openxlsx", "yaml")

# Bioconductor packages required by the backend
.BIOC_PACKAGES <- c("Biostrings", "GenomicFeatures", "rtracklayer")

# Packages required only by the Shiny GUI layer
.GUI_PACKAGES <- c("shiny", "DT", "shinyjs", "bslib")

#' Load a set of packages, collecting any that are missing.
#' @return character vector of packages that could not be loaded
.load_packages <- function(pkgs) {
  missing <- character(0)
  for (pkg in pkgs) {
    ok <- suppressWarnings(suppressPackageStartupMessages(
      requireNamespace(pkg, quietly = TRUE)
    ))
    if (ok) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    } else {
      missing <- c(missing, pkg)
    }
  }
  missing
}

#' Load all backend packages (CRAN + Bioconductor).
#' Stops with an informative message if any are missing.
load_backend_packages <- function() {
  missing <- .load_packages(c(.CRAN_PACKAGES, .BIOC_PACKAGES))
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "), ".\n",
      "Please run  Rscript install.R  from the application folder to install them.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Load GUI packages in addition to the backend ones.
load_gui_packages <- function() {
  load_backend_packages()
  missing <- .load_packages(.GUI_PACKAGES)
  if (length(missing) > 0) {
    stop(
      "Missing required GUI package(s): ", paste(missing, collapse = ", "), ".\n",
      "Please run  Rscript install.R  from the application folder to install them.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
