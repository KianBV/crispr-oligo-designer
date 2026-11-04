# ==============================================================================
# pipeline.R  --  High-level orchestration bridging backend <-> GUI
# ==============================================================================
# Input reading & validation, running the design, selecting top-ranked results,
# building summaries, and producing the Microsynth order form. This layer keeps
# the GUI thin and free of biology logic.
# ==============================================================================

# ------------------------------------------------------------------------------
# Input reading
# ------------------------------------------------------------------------------

#' Read a mutation list from an uploaded file.
#' Supports .csv/.tsv/.txt (a column or a plain list) and .xlsx.
#' Picks a column named "mutation" (case-insensitive) if present, else the first
#' column. Returns a character vector of trimmed, non-empty mutation strings.
read_mutations_file <- function(path, sheet = 1) {
  ext <- tolower(tools::file_ext(path))

  extract_from_df <- function(df) {
    if (ncol(df) == 0) return(character(0))
    nm <- tolower(names(df))
    col <- if ("mutation" %in% nm) which(nm == "mutation")[1] else 1
    as.character(df[[col]])
  }

  raw <- switch(
    ext,
    "xlsx" = extract_from_df(openxlsx::read.xlsx(path, sheet = sheet)),
    "csv"  = extract_from_df(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)),
    "tsv"  = extract_from_df(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)),
    "txt"  = {
      lines <- readLines(path, warn = FALSE)
      unlist(strsplit(lines, "[,;\\s]+", perl = TRUE))
    },
    # default: try to read as a delimited table, fall back to lines
    {
      out <- tryCatch(extract_from_df(utils::read.csv(path, stringsAsFactors = FALSE,
                                                      check.names = FALSE)),
                      error = function(e) NULL)
      if (is.null(out)) unlist(strsplit(readLines(path, warn = FALSE), "[,;\\s]+", perl = TRUE)) else out
    }
  )

  raw <- trimws(raw)
  raw <- raw[nzchar(raw)]
  raw
}

# ------------------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------------------

#' Validate a vector of mutation strings.
#' @param mutations character vector
#' @param system optional initialised system; if given, gene existence is checked
#' @param id_type identifier type to resolve against ("gene", "locus_tag", …,
#'   or "auto"). Requires system$id_maps.
#' @return data.frame(input, gene, wt_aa, position, mut_aa, valid_format,
#'                    gene_found, status, message)
validate_mutations <- function(mutations, system = NULL, id_type = "gene") {
  mutations <- trimws(as.character(mutations))

  id_maps <- if (!is.null(system)) system$id_maps else NULL
  can_check <- !is.null(id_maps) && length(id_maps) > 0

  rows <- lapply(mutations, function(mut) {
    parsed <- parse_mutation(mut)
    if (is.null(parsed)) {
      return(data.frame(
        input = mut, gene = NA_character_, wt_aa = NA_character_,
        position = NA_integer_, mut_aa = NA_character_,
        valid_format = FALSE, gene_found = NA,
        status = "Invalid format",
        message = "Expected <gene><WT><pos><MUT> (e.g. dpsS100E) or <id>:<WT><pos><MUT> (e.g. b0014:S100E)",
        stringsAsFactors = FALSE))
    }

    gene_found <- NA
    status <- "OK"
    message <- ""
    if (can_check) {
      gene_found <- !is.null(resolve_gene_id(parsed$gene, id_maps, id_type))
      if (!gene_found) {
        status <- "Gene not found"
        message <- paste0("'", parsed$gene, "' not found as ",
                          if (identical(id_type, "auto")) "any identifier" else paste0(id_type),
                          " in this genome")
      }
    }

    data.frame(
      input = mut, gene = parsed$gene, wt_aa = parsed$wt_aa,
      position = parsed$position, mut_aa = parsed$mut_aa,
      valid_format = TRUE, gene_found = gene_found,
      status = status, message = message, stringsAsFactors = FALSE)
  })

  df <- do.call(rbind, rows)

  # Flag duplicates (after the first occurrence)
  dup <- duplicated(df$input)
  df$status[dup & df$status == "OK"] <- "Duplicate"
  df$message[dup] <- "Duplicate of an earlier entry"

  df
}

#' Return the subset of validated mutations that are safe to run
#' (valid format, gene found when known, not a duplicate).
runnable_mutations <- function(validation_df) {
  ok <- validation_df$valid_format &
    (is.na(validation_df$gene_found) | validation_df$gene_found) &
    validation_df$status != "Duplicate"
  unique(validation_df$input[ok])
}

# ------------------------------------------------------------------------------
# Running the design
# ------------------------------------------------------------------------------

#' Run the full mutation-library design.
run_design <- function(mutations, system, cfg, progress_callback = NULL,
                       id_type = "gene") {
  top_n <- cfg$design$top_n %||% 3
  n_cores <- cfg$design$n_cores %||% 1
  design_mutation_library_fast(
    mutations, system, top_n = top_n, n_cores = n_cores,
    progress_callback = progress_callback, id_type = id_type
  )
}

#' Select the top-ranked design per mutation (rank == 1), with optional filters.
#' @param min_score keep only designs with score >= min_score (NULL = no filter)
#' @param exclude_mut_aa character vector of mutant amino acids to drop (e.g. "D")
select_top_ranked <- function(results_all, min_score = NULL, exclude_mut_aa = NULL) {
  if (is.null(results_all) || nrow(results_all) == 0) return(results_all)
  top <- results_all[results_all$rank == 1, , drop = FALSE]
  if (!is.null(min_score)) top <- top[top$score >= min_score, , drop = FALSE]
  if (!is.null(exclude_mut_aa) && length(exclude_mut_aa) > 0) {
    top <- top[!(top$mut_aa %in% exclude_mut_aa), , drop = FALSE]
  }
  top
}

# ------------------------------------------------------------------------------
# Microsynth output
# ------------------------------------------------------------------------------

#' Build the oligo order table for the given top-ranked designs (no file I/O).
#' Centralises the config -> generate_oligo_orders mapping so the GUI and the
#' form writer stay in sync.
build_oligo_orders <- function(results_top, cfg) {
  m <- cfg$microsynth
  generate_oligo_orders(
    results_top,
    method          = m$cloning_method  %||% "gibson",
    fwd_tail        = m$fwd_tail        %||% GRNA_FWD_TAIL,
    rev_tail        = m$rev_tail        %||% GRNA_REV_TAIL,
    gg_top_overhang = m$gg_top_overhang %||% GG_TOP_OVERHANG,
    gg_bot_overhang = m$gg_bot_overhang %||% GG_BOT_OVERHANG,
    n_pto_bases     = m$pto_n_bases     %||% PTO_N_BASES,
    fwd_name_suffix = m$fwd_name_suffix %||% GRNA_FWD_SUFFIX,
    rev_name_suffix = m$rev_name_suffix %||% GRNA_REV_SUFFIX,
    gg_top_name_suffix = m$gg_top_name_suffix %||% GG_TOP_SUFFIX,
    gg_bot_name_suffix = m$gg_bot_name_suffix %||% GG_BOT_SUFFIX,
    ha_name_suffix  = m$ha_name_suffix  %||% HA_SUFFIX,
    gg_phosphorylate = isTRUE(m$gg_phosphorylate),
    gg_phos_label   = m$gg_phos_label   %||% "Phosphate"
  )
}

#' QC the oligos produced for the given designs. Returns a QC data.frame.
qc_oligo_orders <- function(results_top, cfg) {
  orders <- build_oligo_orders(results_top, cfg)
  qc_table(orders$oligos, cfg$qc %||% list(),
           ha_suffix = cfg$microsynth$ha_name_suffix %||% "_HA")
}

#' Generate the Microsynth order form for the given top-ranked designs.
#' @return list with output_path and the dedup counts from generate_oligo_orders
generate_microsynth_output <- function(results_top, cfg, output_path) {
  m <- cfg$microsynth
  orders <- build_oligo_orders(results_top, cfg)

  fill_microsynth_form(
    orders$oligos,
    template_path = m$template_path,
    output_path   = output_path,
    sheet         = m$sheet %||% "DNA Order",
    purification  = m$purification %||% "DES",
    scale         = m$scale %||% "0.04",
    condition     = m$condition %||% "Dried",
    datasheet     = m$datasheet %||% "Standard",
    aliquots      = m$aliquots %||% "No"
  )

  list(output_path = output_path,
       n_oligos = orders$n_after,
       n_duplicates = orders$n_duplicates,
       oligos = orders$oligos)
}

# ------------------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------------------

#' Build a human-readable run summary (character vector of lines).
build_run_summary <- function(validation_df, results_all, results_top, cfg,
                              n_input, elapsed_secs = NULL, genome_label = NULL) {
  L <- c()
  add <- function(...) L[[length(L) + 1]] <<- paste0(...)

  add("CRISPR Oligo Designer — Run Summary")
  add("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  if (!is.null(genome_label)) add("Genome: ", genome_label)
  add("")
  add("--- Input ---")
  add("Mutations submitted:      ", n_input)
  add("Valid format:             ", sum(validation_df$valid_format))
  if (any(!is.na(validation_df$gene_found))) {
    add("Gene found in genome:     ", sum(validation_df$gene_found %in% TRUE))
    add("Gene NOT found:           ", sum(validation_df$gene_found %in% FALSE))
  }
  add("Invalid format:           ", sum(!validation_df$valid_format))
  add("Duplicates:               ", sum(validation_df$status == "Duplicate"))
  add("")
  add("--- Design parameters ---")
  add("Oligo size (homology arm): ", cfg$design$oligo_size, " bp")
  add("Max PAM distance:          ", cfg$design$max_pam_distance, " bp")
  add("Min score threshold:       ", cfg$design$min_score_threshold)
  add("Designs per mutation:      ", cfg$design$top_n)
  add("Cloning method:            ", cfg$microsynth$cloning_method %||% "gibson")
  add("")
  add("--- Results ---")
  n_designs <- if (!is.null(results_all)) nrow(results_all) else 0
  add("Total designs generated:   ", n_designs)
  if (!is.null(results_all) && nrow(results_all) > 0) {
    designed_muts <- length(unique(results_all$mutation))
    add("Mutations with a design:   ", designed_muts, " / ", n_input)
    add("Mutations with NO design:  ", n_input - designed_muts)
  }
  if (!is.null(results_top) && nrow(results_top) > 0) {
    add("Top-ranked (rank 1) kept:  ", nrow(results_top))
    psites <- length(unique(paste0(results_top$gene, results_top$position)))
    add("Unique phosphosites:       ", psites)
    add("")
    add("Score: min ", round(min(results_top$score), 2),
        " / mean ", round(mean(results_top$score), 2),
        " / max ", round(max(results_top$score), 2))
    add("Off-target class distribution (top-ranked):")
    for (cls in 1:4) {
      n <- sum(results_top$off_target_class == cls)
      pct <- if (nrow(results_top) > 0) round(100 * n / nrow(results_top), 1) else 0
      add("  Class ", cls, ": ", n, " (", pct, "%)")
    }
  }
  if (!is.null(elapsed_secs)) {
    add("")
    add("Run time: ", round(elapsed_secs, 1), " s")
  }
  unlist(L)
}
