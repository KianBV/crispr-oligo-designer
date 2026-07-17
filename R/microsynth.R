# ==============================================================================
# microsynth.R  --  Oligo order generation + Microsynth form filling
# ==============================================================================
#
# For each top-ranked design this generates 3 oligos:
#   1. <mutation>_pEcgRNA_fw_gibson : 5'-<protospacer><fwd tail>-3'
#   2. <mutation>_pEcgRNA_rv_gibson : 5'-<revcomp(protospacer)><rev tail>-3'
#   3. <mutation>_HA                : 5'-X*X*X*X<rest>-3'  homology arm, PTO on first N bases
#
# The two gRNA primers are Gibson-assembly primers: the protospacer (or its
# reverse complement) carries a 24-nt 3' homology tail matching the pEcgRNA
# vector. PTO (phosphorothioate) bonds are added ONLY to the homology arm.
# Oligos are deduplicated by sequence (first occurrence name kept).
# ==============================================================================

# --- Helper: add PTO bonds between the first N bases ---------------------------
add_pto_5prime <- function(seq, n_pto_bases = PTO_N_BASES) {
  chars <- strsplit(toupper(seq), "")[[1]]
  if (n_pto_bases < 1 || n_pto_bases >= length(chars)) {
    # Degenerate cases: no protection, or protect nothing beyond the whole oligo
    return(toupper(seq))
  }
  pto_part <- paste(chars[1:n_pto_bases], collapse = "*")
  rest <- paste(chars[(n_pto_bases + 1):length(chars)], collapse = "")
  paste0(pto_part, rest)
}

# ==============================================================================
# Generate oligo table from top-ranked designs
# ==============================================================================

#' @param results_top data.frame of top-ranked designs (one row per mutation)
#' @param fwd_tail 3' homology tail appended to the forward protospacer
#' @param rev_tail 3' homology tail appended to the reverse-complement protospacer
#' @param n_pto_bases number of 5' bases to PTO-protect on the homology arm
#' @param fwd_name_suffix/rev_name_suffix/ha_name_suffix oligo name suffixes
#' @return list(oligos = data.frame(oligo_name, sequence),
#'              n_before, n_after, n_duplicates)
generate_oligo_orders <- function(results_top,
                                   fwd_tail = GRNA_FWD_TAIL,
                                   rev_tail = GRNA_REV_TAIL,
                                   n_pto_bases = PTO_N_BASES,
                                   fwd_name_suffix = GRNA_FWD_SUFFIX,
                                   rev_name_suffix = GRNA_REV_SUFFIX,
                                   ha_name_suffix = HA_SUFFIX) {

  if (nrow(results_top) == 0) {
    return(list(oligos = data.frame(oligo_name = character(0),
                                    sequence = character(0),
                                    stringsAsFactors = FALSE),
                n_before = 0, n_after = 0, n_duplicates = 0))
  }

  oligo_rows <- list()

  for (i in seq_len(nrow(results_top))) {
    row <- results_top[i, ]
    # Sanitise the label for oligo names: identifier delimiters (: | /) and
    # other punctuation are replaced with '_' so the Microsynth Name is clean.
    mutation_label <- gsub("[^A-Za-z0-9_]", "_", row$mutation)
    protospacer <- toupper(row$protospacer)
    homology_arm <- toupper(row$homology_arm)

    # Forward gRNA primer: protospacer + 3' homology tail
    proto_fwd_seq <- paste0(protospacer, fwd_tail)

    # Reverse gRNA primer: reverse complement of protospacer + 3' homology tail
    proto_rev_seq <- paste0(reverse_complement(protospacer), rev_tail)

    # Homology arm with PTO on first N bases
    ha_seq_pto <- add_pto_5prime(homology_arm, n_pto_bases = n_pto_bases)

    oligo_rows[[length(oligo_rows) + 1]] <- data.frame(
      oligo_name = paste0(mutation_label, fwd_name_suffix), sequence = proto_fwd_seq,
      stringsAsFactors = FALSE)
    oligo_rows[[length(oligo_rows) + 1]] <- data.frame(
      oligo_name = paste0(mutation_label, rev_name_suffix), sequence = proto_rev_seq,
      stringsAsFactors = FALSE)
    oligo_rows[[length(oligo_rows) + 1]] <- data.frame(
      oligo_name = paste0(mutation_label, ha_name_suffix), sequence = ha_seq_pto,
      stringsAsFactors = FALSE)
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

  list(oligos = oligo_table, n_before = n_before, n_after = n_after,
       n_duplicates = n_before - n_after)
}

# ==============================================================================
# Write into the Microsynth order form template
# ==============================================================================

#' Fill the Microsynth DNA order form (.xlsx) using the provided template.
#'
#' @param oligo_table data.frame with columns oligo_name, sequence
#' @param template_path path to the official Microsynth .xlsx template
#' @param output_path where to write the filled form
#' @param sheet worksheet name in the template
#' @param purification/scale/condition/datasheet/aliquots default column values
fill_microsynth_form <- function(oligo_table,
                                 template_path,
                                 output_path = "MicrosynthOrder_filled.xlsx",
                                 sheet = "DNA Order",
                                 purification = "DES",
                                 scale = "0.04",
                                 condition = "Dried",
                                 datasheet = "Standard",
                                 aliquots = "No") {

  if (is.null(template_path) || !file.exists(template_path)) {
    stop(
      "Microsynth template not found at: ",
      if (is.null(template_path)) "<not set>" else template_path, "\n",
      "Please place the official 'MicrosynthUploadFormDNA__2_.xlsx' file in the ",
      "templates/ folder (or point config.yml -> microsynth.template_path at it).",
      call. = FALSE
    )
  }

  wb <- openxlsx::loadWorkbook(template_path)

  sheets <- names(wb)
  if (!(sheet %in% sheets)) {
    stop("Sheet '", sheet, "' not found in template. Available sheets: ",
         paste(sheets, collapse = ", "), call. = FALSE)
  }

  for (i in seq_len(nrow(oligo_table))) {
    row_idx <- i + 1  # row 1 is header, data starts at row 2

    openxlsx::writeData(wb, sheet, oligo_table$oligo_name[i], startRow = row_idx, startCol = 1)  # A: Name
    openxlsx::writeData(wb, sheet, oligo_table$sequence[i],   startRow = row_idx, startCol = 2)  # B: Sequence
    # C (Length) is a formula in the template — leave it
    openxlsx::writeData(wb, sheet, purification, startRow = row_idx, startCol = 4)   # D: Purification
    openxlsx::writeData(wb, sheet, scale,        startRow = row_idx, startCol = 5)   # E: Scale
    # F-K (modifications) — leave empty
    openxlsx::writeData(wb, sheet, condition,  startRow = row_idx, startCol = 12)  # L: Physical Condition
    openxlsx::writeData(wb, sheet, datasheet,  startRow = row_idx, startCol = 13)  # M: Datasheet
    openxlsx::writeData(wb, sheet, aliquots,   startRow = row_idx, startCol = 14)  # N: Aliquots
  }

  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  cat("Microsynth order form saved to:", output_path, "\n")
  cat("Total oligos:", nrow(oligo_table), "\n")
  invisible(output_path)
}
