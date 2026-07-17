# ==============================================================================
# OUTPUT: Protospacer oligos with overhangs + Homology arm with PTO
# Fills the Microsynth order form (.xlsx)
# ==============================================================================
#
# For each design in results_top, this generates 3 oligos:
#   1. Protospacer_Fwd:  5'-TAGT<protospacer>-3'
#   2. Protospacer_Rev:  5'-AAAC<revcomp of protospacer>-3'
#   3. Homology_arm:     5'-X*X*X*XREST...-3'  (PTO on first 4 bases)
#
# When the two protospacer oligos are annealed, they form:
#   5'-TAGT<protospacer>-3'
#        3'-<revcomp>CAAA-5'
# i.e. dsDNA with 5'-TAGT overhang on top and 5'-AAAC overhang on bottom.
#
# ==============================================================================

# --- Load required package ---
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)
library(dplyr)

# --- Helper: reverse complement ---
reverse_complement <- function(seq) {
  complement <- chartr("ACGT", "TGCA", toupper(seq))
  paste(rev(strsplit(complement, "")[[1]]), collapse = "")
}

# --- Helper: add PTO bonds between the first 4 bases ---
add_pto_5prime <- function(seq, n_pto_bases = 4) {
  chars <- strsplit(toupper(seq), "")[[1]]
  pto_part <- paste(chars[1:n_pto_bases], collapse = "*")
  rest <- paste(chars[(n_pto_bases + 1):length(chars)], collapse = "")
  paste0(pto_part, rest)
}

# ==============================================================================
# Generate oligo table from results_top
# ==============================================================================

generate_oligo_orders <- function(results_top) {
  
  oligo_rows <- list()
  
  for (i in seq_len(nrow(results_top))) {
    row <- results_top[i, ]
    mutation_label <- row$mutation
    protospacer <- toupper(row$protospacer)
    homology_arm <- toupper(row$homology_arm)
    
    # Protospacer forward: 5'-TAGT + protospacer
    proto_fwd_seq <- paste0("TAGT", protospacer)
    
    # Protospacer reverse: 5'-AAAC + reverse complement of protospacer
    proto_rev_seq <- paste0("AAAC", reverse_complement(protospacer))
    
    # Homology arm with PTO on first 4 bases
    ha_seq_pto <- add_pto_5prime(homology_arm, n_pto_bases = 4)
    
    # Build oligo entries
    oligo_rows[[length(oligo_rows) + 1]] <- data.frame(
      oligo_name  = paste0(mutation_label, "_gRNA_F"),
      sequence    = proto_fwd_seq,
      stringsAsFactors = FALSE
    )
    oligo_rows[[length(oligo_rows) + 1]] <- data.frame(
      oligo_name  = paste0(mutation_label, "_gRNA_R"),
      sequence    = proto_rev_seq,
      stringsAsFactors = FALSE
    )
    oligo_rows[[length(oligo_rows) + 1]] <- data.frame(
      oligo_name  = paste0(mutation_label, "_HA"),
      sequence    = ha_seq_pto,
      stringsAsFactors = FALSE
    )
  }
  
  oligo_table <- do.call(rbind, oligo_rows)
  
  # Deduplicate: remove oligos with identical sequences, keep first occurrence
  n_before <- nrow(oligo_table)
  oligo_table <- oligo_table[!duplicated(oligo_table$sequence), ]
  n_after <- nrow(oligo_table)
  
  if (n_before > n_after) {
    cat("Removed", n_before - n_after, "duplicate oligos (by sequence).",
        n_after, "unique oligos remaining.\n")
  }
  
  oligo_table
}

# ==============================================================================
# Write into Microsynth order form
# ==============================================================================

fill_microsynth_form <- function(oligo_table,
                                 template_path,
                                 output_path = "MicrosynthOrder_filled.xlsx",
                                 purification = "DES",
                                 scale = "0.04",
                                 condition = "Dried",
                                 datasheet = "Standard",
                                 aliquots = "No") {
  
  wb <- loadWorkbook(template_path)
  sheet <- "DNA Order"
  
  for (i in seq_len(nrow(oligo_table))) {
    row_idx <- i + 1  # row 1 is header, data starts at row 2
    
    writeData(wb, sheet, oligo_table$oligo_name[i], startRow = row_idx, startCol = 1)  # A: Name
    writeData(wb, sheet, oligo_table$sequence[i],   startRow = row_idx, startCol = 2)  # B: Sequence
    # C (Length) is a formula — leave it
    writeData(wb, sheet, purification, startRow = row_idx, startCol = 4)   # D: Purification
    writeData(wb, sheet, scale,        startRow = row_idx, startCol = 5)   # E: Scale
    # F-K (modifications) — leave empty
    writeData(wb, sheet, condition,  startRow = row_idx, startCol = 12)  # L: Physical Condition
    writeData(wb, sheet, datasheet, startRow = row_idx, startCol = 13)  # M: Datasheet
    writeData(wb, sheet, aliquots,  startRow = row_idx, startCol = 14)  # N: Aliquots
  }
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  cat("Microsynth order form saved to:", output_path, "\n")
  cat("Total oligos:", nrow(oligo_table), "\n")
}

# ==============================================================================
# USAGE — paste after your existing pipeline
# ==============================================================================
#
# # Assuming results_top is already generated:
# # results_top <- results_all %>% filter(rank == 1)
#
# oligo_table <- generate_oligo_orders(results_top)
#
# # Preview
# print(oligo_table)
#
# # Fill the Microsynth form
# fill_microsynth_form(
#   oligo_table,
#   template_path = "MicrosynthUploadFormDNA__2_.xlsx",
#   output_path   = "MicrosynthOrder_filled.xlsx"
# )
#
