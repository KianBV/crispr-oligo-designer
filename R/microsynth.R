# ==============================================================================
# microsynth.R  --  Oligo order generation + Microsynth form filling
# ==============================================================================
#
# Per top-ranked design, the spacer can be delivered two ways (config:
# microsynth.cloning_method = "gibson" | "golden_gate" | "both"):
#
#   Gibson (whole-plasmid PCR + Gibson/overlap assembly):
#     <mutation>_pEcgRNA_fw_gibson : 5'-<protospacer><fwd tail>-3'
#     <mutation>_pEcgRNA_rv_gibson : 5'-<revcomp(protospacer)><rev tail>-3'
#
#   Golden Gate (annealed oligo pair, ligated into a BsaI-linearised vector):
#     <mutation>_gg_bsaI_top : 5'-<top overhang><protospacer>-3'      (default TAGT)
#     <mutation>_gg_bsaI_bot : 5'-<bot overhang><revcomp(protospacer)>-3' (default AAAC)
#
#   Repair donor (always):
#     <mutation>_HA : 5'-X*X*X*X<rest>-3'  homology arm, PTO on first N bases
#
# PTO (phosphorothioate) bonds are added ONLY to the homology arm. The Golden
# Gate oligos can optionally carry a 5'-phosphate for ligation. Oligos are
# deduplicated by sequence (first occurrence name kept).
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
#' @param method "gibson", "golden_gate" or "both" (which gRNA oligos to emit)
#' @param fwd_tail/rev_tail 3' homology tails for the Gibson primers
#' @param gg_top_overhang/gg_bot_overhang 4-nt BsaI fusion overhangs (Golden Gate)
#' @param n_pto_bases number of 5' bases to PTO-protect on the homology arm
#' @param *_name_suffix oligo name suffixes
#' @param gg_phosphorylate add a 5'-phosphate modification to the Golden Gate oligos
#' @param gg_phos_label the Microsynth 5'-modification term to write for that phosphate
#' @return list(oligos = data.frame(oligo_name, sequence, mod5),
#'              n_before, n_after, n_duplicates)
generate_oligo_orders <- function(results_top,
                                   method = "gibson",
                                   fwd_tail = GRNA_FWD_TAIL,
                                   rev_tail = GRNA_REV_TAIL,
                                   gg_top_overhang = GG_TOP_OVERHANG,
                                   gg_bot_overhang = GG_BOT_OVERHANG,
                                   n_pto_bases = PTO_N_BASES,
                                   fwd_name_suffix = GRNA_FWD_SUFFIX,
                                   rev_name_suffix = GRNA_REV_SUFFIX,
                                   gg_top_name_suffix = GG_TOP_SUFFIX,
                                   gg_bot_name_suffix = GG_BOT_SUFFIX,
                                   ha_name_suffix = HA_SUFFIX,
                                   gg_phosphorylate = FALSE,
                                   gg_phos_label = "Phosphate") {

  method <- match.arg(method, c("gibson", "golden_gate", "both"))
  do_gibson <- method %in% c("gibson", "both")
  do_gg     <- method %in% c("golden_gate", "both")

  empty <- data.frame(oligo_name = character(0), sequence = character(0),
                      mod5 = character(0), stringsAsFactors = FALSE)
  if (nrow(results_top) == 0) {
    return(list(oligos = empty, n_before = 0, n_after = 0, n_duplicates = 0))
  }

  oligo_rows <- list()
  addrow <- function(name, seq, mod5 = "") {
    oligo_rows[[length(oligo_rows) + 1]] <<- data.frame(
      oligo_name = name, sequence = seq, mod5 = mod5, stringsAsFactors = FALSE)
  }

  for (i in seq_len(nrow(results_top))) {
    row <- results_top[i, ]
    # Sanitise the label for oligo names: identifier delimiters (: | /) and
    # other punctuation are replaced with '_' so the Microsynth Name is clean.
    label <- gsub("[^A-Za-z0-9_]", "_", row$mutation)
    protospacer <- toupper(row$protospacer)
    homology_arm <- toupper(row$homology_arm)
    rc <- reverse_complement(protospacer)

    # --- Gibson primers: protospacer/revcomp + 3' homology tail ---
    if (do_gibson) {
      addrow(paste0(label, fwd_name_suffix), paste0(protospacer, fwd_tail))
      addrow(paste0(label, rev_name_suffix), paste0(rc, rev_tail))
    }

    # --- Golden Gate oligos: 4-nt fusion overhang + protospacer/revcomp ---
    if (do_gg) {
      mod <- if (isTRUE(gg_phosphorylate)) gg_phos_label else ""
      addrow(paste0(label, gg_top_name_suffix), paste0(gg_top_overhang, protospacer), mod)
      addrow(paste0(label, gg_bot_name_suffix), paste0(gg_bot_overhang, rc), mod)
    }

    # --- Repair-donor homology arm with PTO on the first N bases (always) ---
    addrow(paste0(label, ha_name_suffix), add_pto_5prime(homology_arm, n_pto_bases = n_pto_bases))
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
    # F: 5' modification (e.g. phosphate for Golden Gate oligos) — only if set
    if (!is.null(oligo_table$mod5) && nzchar(oligo_table$mod5[i])) {
      openxlsx::writeData(wb, sheet, oligo_table$mod5[i], startRow = row_idx, startCol = 6)
    }
    # G-K (modifications) — leave empty
    openxlsx::writeData(wb, sheet, condition,  startRow = row_idx, startCol = 12)  # L: Physical Condition
    openxlsx::writeData(wb, sheet, datasheet,  startRow = row_idx, startCol = 13)  # M: Datasheet
    openxlsx::writeData(wb, sheet, aliquots,   startRow = row_idx, startCol = 14)  # N: Aliquots
  }

  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  cat("Microsynth order form saved to:", output_path, "\n")
  cat("Total oligos:", nrow(oligo_table), "\n")
  invisible(output_path)
}
