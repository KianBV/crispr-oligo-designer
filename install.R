# ==============================================================================
# install.R  --  One-time dependency installation
# ==============================================================================
# Run once after cloning/updating:
#   Rscript install.R
# ==============================================================================

cran_packages <- c("shiny", "DT", "shinyjs", "bslib",
                   "dplyr", "stringr", "data.table", "openxlsx", "yaml")
bioc_packages <- c("Biostrings", "GenomicFeatures", "rtracklayer")
# Optional: Gviz powers the design inspector's genomic map (the app works
# without it — the map degrades to a message; the sequence view, components and
# off-targets still work). It is heavy, so it is installed best-effort below.
bioc_optional <- c("Gviz")

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
message("\nAll core R dependencies installed. Launch the app with:  Rscript run.R")

# Optional: Gviz (genomic map in the design inspector), best-effort.
to_install_opt <- setdiff(bioc_optional, rownames(installed.packages()))
if (length(to_install_opt) > 0) {
  message("\nInstalling optional Gviz (genomic map; this can take a few minutes)...")
  tryCatch(BiocManager::install(to_install_opt, ask = FALSE, update = FALSE),
           error = function(e) message("  -> Gviz install skipped/failed: ", conditionMessage(e),
                                       "\n     The app still runs; the inspector map will show a hint."))
}
message(if ("Gviz" %in% rownames(installed.packages()))
        "  -> Gviz detected: genomic map is ON." else
        "  -> Gviz not installed; the inspector map will be unavailable (everything else works).")

# Optional: accurate thermodynamic QC (hairpin / dimer Delta-G) via primer3-py.
# The app works without it (falls back to built-in heuristics), but for true
# Delta-G install the Python package:
message("\nOptional, for accurate oligo QC (true hairpin/dimer Delta-G):")
message("    pip install primer3-py      (needs Python on PATH)")
py <- Sys.which("python")
if (nzchar(py)) {
  ok <- tryCatch(system2(py, c("-c", "import primer3"), stdout = FALSE, stderr = FALSE) == 0,
                 error = function(e) FALSE)
  message(if (isTRUE(ok)) "  -> primer3-py detected: accurate QC is ON."
          else "  -> primer3-py not found yet (QC will use the built-in heuristic).")
} else {
  message("  -> Python not on PATH; QC will use the built-in heuristic.")
}
