# ==============================================================================
# install.R  --  One-time dependency installation
# ==============================================================================
# Run once after cloning/updating:
#   Rscript install.R
# ==============================================================================

cran_packages <- c("shiny", "DT", "shinyjs", "bslib",
                   "dplyr", "stringr", "data.table", "openxlsx", "yaml")
bioc_packages <- c("Biostrings", "GenomicFeatures", "rtracklayer")

message("== Installing CRAN packages ==")
installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN packages already installed.")
}

message("== Installing Bioconductor packages ==")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
installed <- rownames(installed.packages())
to_install_bioc <- setdiff(bioc_packages, installed)
if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

# Verify
missing <- setdiff(c(cran_packages, bioc_packages), rownames(installed.packages()))
if (length(missing) > 0) {
  stop("These packages failed to install: ", paste(missing, collapse = ", "))
}
message("\nAll dependencies installed. Launch the app with:  Rscript run.R")
