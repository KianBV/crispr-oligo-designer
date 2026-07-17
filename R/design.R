# ==============================================================================
# design.R  --  gRNA design, homology-arm & oligo assembly, library design
# ==============================================================================
# Core logic preserved verbatim from the original (fixed) pipeline. The only
# behavioural change: the homology-arm length is now a parameter (default
# HOMOLOGY_ARM_LENGTH = 100 bp) instead of a hard-coded 85, so it is tunable.
# ==============================================================================

# --- Small helpers ------------------------------------------------------------

translate_codon <- function(codon) {
  codon <- toupper(codon)
  if (codon %in% names(REVERSE_CODON_TABLE)) unname(REVERSE_CODON_TABLE[codon]) else NA
}

get_codons_for_aa <- function(aa) {
  aa <- toupper(aa)
  if (aa %in% names(CODON_TABLE)) CODON_TABLE[[aa]] else NULL
}

#' Parse a mutation entry into gene identifier + amino-acid substitution.
#'
#' Two accepted forms:
#'   * Compact (gene symbol only, alphabetic): dpsS100E
#'   * Delimited (any identifier, incl. digits): <id>:<WT><pos><MUT>, where the
#'     delimiter is ':', '|' or '/'. e.g. BW25113_0001:S100E, P0A7B8:S100E
#' The gene identifier is returned verbatim in `$gene`; resolve it to a gene with
#' resolve_gene_id() using the chosen identifier type.
parse_mutation <- function(mutation_string) {
  s <- trimws(mutation_string)

  # Delimited form: identifier can contain any non-delimiter character
  m <- regmatches(s, regexec("^(.+)[:|/]([A-Za-z])([0-9]+)([A-Za-z])$", s))[[1]]
  if (length(m) == 5) {
    return(list(gene = m[2], wt_aa = toupper(m[3]), position = as.integer(m[4]),
                mut_aa = toupper(m[5]), original_string = s))
  }

  # Compact form: gene symbol must be purely alphabetic
  m <- regmatches(s, regexec("^([A-Za-z]+)([A-Z])([0-9]+)([A-Z])$", s))[[1]]
  if (length(m) == 5) {
    return(list(gene = m[2], wt_aa = m[3], position = as.integer(m[4]),
                mut_aa = m[5], original_string = s))
  }

  NULL
}

calculate_cassette_score <- function(distance, off_target_class, num_substitutions = 5) {
  base_scores <- c(5.0, 4.41, 3.9, 2.0)
  distance_steps <- floor(min(distance, 30) / 3)
  score <- base_scores[off_target_class] - (distance_steps * 0.2)

  max(score, 0)
}

# --- GC content helper (used by silent-mutation utilities) --------------------
calculate_gc_content <- function(seq) {
  seq <- toupper(seq)
  chars <- strsplit(seq, "")[[1]]
  gc_count <- sum(chars %in% c("G", "C"))
  gc_count / length(chars)
}

generate_silent_pam_mutations <- function(gene_sequence, pam_position, strand, gene_start = 1,
                                           codon_freq = ECOLI_CODON_FREQ) {
  seq_chars <- strsplit(toupper(gene_sequence), "")[[1]]
  mutations <- list()

  if (strand == "+") {
    gg_positions <- c(pam_position + 1, pam_position + 2)
    target_nts <- c("A", "C", "T")
  } else {
    gg_positions <- c(pam_position, pam_position + 1)
    target_nts <- c("A", "G", "T")
  }

  for (pos in gg_positions) {
    if (pos < 1 || pos > length(seq_chars)) next

    codon_num <- ((pos - gene_start) %/% 3) + 1
    codon_start <- (codon_num - 1) * 3 + gene_start
    codon_pos <- ((pos - gene_start) %% 3) + 1

    if (codon_start < 1 || codon_start + 2 > length(seq_chars)) next

    original_codon <- paste(seq_chars[codon_start:(codon_start + 2)], collapse = "")
    original_aa <- translate_codon(original_codon)
    if (is.na(original_aa)) next

    for (new_nt in target_nts) {
      if (new_nt == seq_chars[pos]) next

      new_codon_chars <- seq_chars[codon_start:(codon_start + 2)]
      new_codon_chars[codon_pos] <- new_nt
      new_codon <- paste(new_codon_chars, collapse = "")
      new_aa <- translate_codon(new_codon)

      if (!is.na(new_aa) && new_aa == original_aa && new_aa != "*") {

        original_freq <- as.numeric(codon_freq[original_codon])
        new_freq <- as.numeric(codon_freq[new_codon])
        freq_diff <- abs(new_freq - original_freq)

        mutations[[length(mutations) + 1]] <- list(
          position = pos, original_nt = seq_chars[pos], mutant_nt = new_nt,
          original_codon = original_codon, mutant_codon = new_codon,
          codon_change = paste(original_codon, "->", new_codon), amino_acid = original_aa,
          original_freq = original_freq, new_freq = new_freq, freq_diff = freq_diff
        )
      }
    }
  }
  # Sort by codon usage similarity (smallest freq_diff first)
  if (length(mutations) > 1) {
    freq_diffs <- sapply(mutations, `[[`, "freq_diff")
    mutations <- mutations[order(freq_diffs)]
  }
  mutations
}

# ==============================================================================
# COORDINATE HELPER
# ==============================================================================

#' Convert a gene-local position (within the coding sequence) to a genomic
#' position. For plus-strand genes: genomic_start + local_pos - 1. For
#' minus-strand genes the coding sequence is the reverse complement of the
#' genomic region, so position 1 maps to genomic_end.
gene_local_to_genomic <- function(local_pos, gene_genomic_start, gene_genomic_end, gene_strand) {
  if (gene_strand == "-") {
    gene_genomic_end - local_pos - 2 + 1
  } else {
    gene_genomic_start + local_pos - 1
  }
}

# ==============================================================================
# gRNA DESIGN
# ==============================================================================

#' Design gRNA candidates for a mutation
design_grna_for_mutation_fast <- function(gene_sequence, mutation, gene_genomic_start = 1,
                                          pam_data = NULL, gene_genomic_end = NULL,
                                          gene_strand = "+", codon_freq = ECOLI_CODON_FREQ) {
  seq_upper <- toupper(gene_sequence)
  gene_length_bp <- nchar(seq_upper)
  gene_length_aa <- gene_length_bp %/% 3

  # If gene_genomic_end not provided, estimate it (works correctly for + strand)
  if (is.null(gene_genomic_end)) {
    gene_genomic_end <- gene_genomic_start + gene_length_bp - 1
  }

  codon_start <- (mutation$position - 1) * 3 + 1
  codon_end <- codon_start + 2

  if (mutation$position > gene_length_aa) {
    warning("Position exceeds gene length")
    return(list())
  }

  original_codon <- substr(seq_upper, codon_start, codon_end)
  mut_codons <- get_codons_for_aa(mutation$mut_aa)
  if (is.null(mut_codons)) return(list())

  search_start <- max(1, codon_start - MAX_PAM_DISTANCE - 23)
  search_end <- min(gene_length_bp, codon_end + MAX_PAM_DISTANCE + 3)
  search_region <- substr(seq_upper, search_start, search_end)

  candidates <- list()

  # Forward PAMs (NGG in gene-local coordinates)
  fwd_matches <- gregexpr("(?=.GG)", search_region, perl = TRUE)[[1]]
  if (fwd_matches[1] != -1) {
    for (rel_pos in fwd_matches) {
      abs_pos <- search_start + rel_pos - 1
      distance <- abs(abs_pos - codon_start)

      if (distance <= MAX_PAM_DISTANCE) {
        proto_start <- abs_pos - 20
        proto_end <- abs_pos - 1

        if (proto_start >= 1 && proto_end <= gene_length_bp) {
          protospacer <- substr(seq_upper, proto_start, proto_end)

          if (nchar(protospacer) == 20) {
            silent_pam_muts <- generate_silent_pam_mutations(seq_upper, abs_pos, "+", 1,
                                                              codon_freq = codon_freq)

            if (length(silent_pam_muts) > 0) {
              if (!is.null(pam_data)) {
                target_genomic_pos <- gene_local_to_genomic(
                  abs_pos, gene_genomic_start, gene_genomic_end, gene_strand
                )
                ot_result <- count_off_targets_fast(protospacer, pam_data, target_genomic_pos)
              } else {
                ot_result <- list(off_target_class = 1, seed_matches = 0, off_target_count = 0)
              }

              score <- calculate_cassette_score(distance, ot_result$off_target_class)

              if (score >= MIN_SCORE_THRESHOLD) {
                candidates[[length(candidates) + 1]] <- list(
                  protospacer = protospacer, pam_position = abs_pos, strand = "+",
                  distance = distance, off_target_class = ot_result$off_target_class,
                  seed_matches = ot_result$seed_matches, off_target_count = ot_result$off_target_count,
                  score = score, original_codon = original_codon,
                  mutant_codon = mut_codons[1], pam_mutation = silent_pam_muts[[1]],
                  mutation_position = mutation$position, codon_start = codon_start
                )
              }
            }
          }
        }
      }
    }
  }

  # Reverse PAMs (CCN in gene-local coordinates)
  rev_matches <- gregexpr("(?=CC.)", search_region, perl = TRUE)[[1]]
  if (rev_matches[1] != -1) {
    for (rel_pos in rev_matches) {
      abs_pos <- search_start + rel_pos - 1
      distance <- abs(abs_pos - codon_start)

      if (distance <= MAX_PAM_DISTANCE) {
        proto_start <- abs_pos + 3
        proto_end <- abs_pos + 22

        if (proto_start >= 1 && proto_end <= gene_length_bp) {
          protospacer <- reverse_complement(substr(seq_upper, proto_start, proto_end))

          if (nchar(protospacer) == 20) {
            silent_pam_muts <- generate_silent_pam_mutations(seq_upper, abs_pos, "-", 1,
                                                              codon_freq = codon_freq)

            if (length(silent_pam_muts) > 0) {
              if (!is.null(pam_data)) {
                target_genomic_pos <- gene_local_to_genomic(
                  abs_pos, gene_genomic_start, gene_genomic_end, gene_strand
                )
                ot_result <- count_off_targets_fast(protospacer, pam_data, target_genomic_pos)
              } else {
                ot_result <- list(off_target_class = 1, seed_matches = 0, off_target_count = 0)
              }

              score <- calculate_cassette_score(distance, ot_result$off_target_class)

              if (score >= MIN_SCORE_THRESHOLD) {
                candidates[[length(candidates) + 1]] <- list(
                  protospacer = protospacer, pam_position = abs_pos, strand = "-",
                  distance = distance, off_target_class = ot_result$off_target_class,
                  seed_matches = ot_result$seed_matches, off_target_count = ot_result$off_target_count,
                  score = score, original_codon = original_codon,
                  mutant_codon = mut_codons[1], pam_mutation = silent_pam_muts[[1]],
                  mutation_position = mutation$position, codon_start = codon_start
                )
              }
            }
          }
        }
      }
    }
  }

  # Sort by score
  if (length(candidates) > 0) {
    scores <- sapply(candidates, `[[`, "score")
    candidates <- candidates[order(scores, decreasing = TRUE)]
  }

  candidates
}

# ==============================================================================
# OLIGO GENERATION
# ==============================================================================

#' Build the single-stranded homology-arm (repair template) oligo.
#' @param homology_arm_length Length of the repair oligo in bp. Defaults to the
#'   configurable global HOMOLOGY_ARM_LENGTH (the "Microsynth oligo size").
generate_homology_arm <- function(gene_sequence, design, genome_seq = NULL,
                                   gene_genomic_start = NULL, gene_strand = NULL,
                                   homology_arm_length = HOMOLOGY_ARM_LENGTH) {
  seq_chars <- strsplit(toupper(gene_sequence), "")[[1]]
  gene_length <- length(seq_chars)

  mut_seq_chars <- seq_chars
  codon_start <- design$codon_start
  codon_chars <- strsplit(design$mutant_codon, "")[[1]]
  mut_seq_chars[codon_start:(codon_start + 2)] <- codon_chars

  pam_mut <- design$pam_mutation
  mut_seq_chars[pam_mut$position] <- pam_mut$mutant_nt

  mut_positions <- c(codon_start, codon_start + 1, codon_start + 2, pam_mut$position)
  edit_center <- round(mean(mut_positions))

  half_arm <- homology_arm_length %/% 2
  arm_start <- edit_center - half_arm
  arm_end <- arm_start + homology_arm_length - 1

  if (arm_start < 1) { arm_start <- 1; arm_end <- homology_arm_length }
  if (arm_end > gene_length) { arm_end <- gene_length; arm_start <- max(1, arm_end - homology_arm_length + 1) }

  arm_seq <- paste(mut_seq_chars[arm_start:arm_end], collapse = "")

  list(sequence = arm_seq, start = arm_start, end = arm_end,
       length = nchar(arm_seq), edit_center = edit_center)
}

generate_200bp_oligo <- function(design, homology_arm) {
  oligo_sequence <- paste0(SPACER1, homology_arm$sequence, SPACER2,
                           J23119_PROMOTER, design$protospacer, SGRNA_HANDLE)
  list(oligo_sequence = oligo_sequence, length = nchar(oligo_sequence))
}

# ==============================================================================
# MAIN LIBRARY DESIGN
# ==============================================================================

#' Design mutation library
#' @param mutations Vector of mutation strings
#' @param system Initialised system (from initialize_system_fast)
#' @param top_n Number of designs per mutation
#' @param n_cores Number of cores (NULL/1 = sequential)
#' @param progress_callback Optional function(fraction, message) for GUI progress
#' @param id_type Identifier type for gene lookup ("gene", "locus_tag", …, or
#'   "auto" to search all). Uses system$id_maps.
design_mutation_library_fast <- function(mutations, system, top_n = 3, n_cores = NULL,
                                          progress_callback = NULL, id_type = "gene") {

  gene_sequences <- system$gene_sequences
  gene_metadata <- system$gene_metadata
  id_maps <- system$id_maps
  pam_data <- system$pam_data
  genome_seq <- system$genome_seq
  codon_freq <- system$codon_freq %||% ECOLI_CODON_FREQ

  n_mutations <- length(mutations)
  cat("\nProcessing", n_mutations, "mutations\n")

  # Worker function
  process_one <- function(mut_str) {
    mutation <- parse_mutation(mut_str)
    if (is.null(mutation)) return(NULL)

    key <- resolve_gene_id(mutation$gene, id_maps, id_type)
    if (is.null(key)) return(NULL)

    gene_seq <- gene_sequences[[key]]
    if (is.null(gene_seq)) return(NULL)

    gene_meta <- gene_metadata[[key]]
    gene_genomic_start <- if (!is.null(gene_meta)) gene_meta$genomic_start else 1
    gene_genomic_end   <- if (!is.null(gene_meta)) gene_meta$genomic_end else NULL
    gene_strand        <- if (!is.null(gene_meta)) gene_meta$strand else "+"

    candidates <- tryCatch(
      design_grna_for_mutation_fast(gene_seq, mutation, gene_genomic_start, pam_data,
                                     gene_genomic_end = gene_genomic_end,
                                     gene_strand = gene_strand,
                                     codon_freq = codon_freq),
      error = function(e) list()
    )

    if (length(candidates) == 0) return(NULL)

    top_candidates <- head(candidates, top_n)

    results <- lapply(seq_along(top_candidates), function(i) {
      design <- top_candidates[[i]]
      homology_arm <- generate_homology_arm(gene_seq, design, genome_seq,
                                             gene_meta$genomic_start, gene_meta$strand)
      oligo <- generate_200bp_oligo(design, homology_arm)

      data.frame(
        mutation = mut_str, gene = mutation$gene, position = mutation$position,
        wt_aa = mutation$wt_aa, mut_aa = mutation$mut_aa, rank = i,
        score = round(design$score, 3), off_target_class = design$off_target_class,
        seed_matches = design$seed_matches, off_target_count = design$off_target_count,
        pam_distance = design$distance, strand = design$strand,
        protospacer = design$protospacer, original_codon = design$original_codon,
        mutant_codon = design$mutant_codon, pam_mutation = design$pam_mutation$codon_change,
        pam_mutation_freq_diff = design$pam_mutation$freq_diff,
        homology_arm = homology_arm$sequence, oligo_200bp = oligo$oligo_sequence,
        oligo_length = oligo$length, stringsAsFactors = FALSE
      )
    })

    do.call(rbind, results)
  }

  start_time <- Sys.time()

  if (is.null(n_cores) || n_cores == 1) {
    # Sequential with progress
    results <- vector("list", n_mutations)
    for (i in seq_len(n_mutations)) {
      results[[i]] <- process_one(mutations[i])
      if (!is.null(progress_callback) && (i %% 25 == 0 || i == n_mutations)) {
        progress_callback(i / n_mutations,
                          sprintf("Designing: %d / %d mutations", i, n_mutations))
      }
      if (i %% 100 == 0 || i == n_mutations) {
        cat(sprintf("\r  Progress: %d/%d (%.1f%%)    ", i, n_mutations, 100 * i / n_mutations))
      }
    }
    cat("\n")
  } else {
    # Parallel using mclapply (Unix) or parLapply (Windows)
    cat("Using", n_cores, "cores\n")

    if (.Platform$OS.type == "unix") {
      results <- parallel::mclapply(mutations, process_one, mc.cores = n_cores)
    } else {
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, c("gene_sequences", "gene_metadata", "id_maps",
                                    "id_type", "pam_data",
                                    "genome_seq", "top_n", "codon_freq"), envir = environment())
      parallel::clusterExport(cl, c("parse_mutation", "resolve_gene_id",
                                    "design_grna_for_mutation_fast",
                                    "generate_homology_arm", "generate_200bp_oligo",
                                    "count_off_targets_fast", "generate_silent_pam_mutations",
                                    "calculate_cassette_score", "get_codons_for_aa",
                                    "translate_codon", "reverse_complement",
                                    "gene_local_to_genomic",
                                    "CODON_TABLE", "REVERSE_CODON_TABLE", "ECOLI_CODON_FREQ",
                                    "MAX_PAM_DISTANCE", "HOMOLOGY_ARM_LENGTH",
                                    "MIN_SCORE_THRESHOLD", "MAX_MISMATCHES",
                                    "SPACER1", "SPACER2", "J23119_PROMOTER", "SGRNA_HANDLE", "%||%"),
                              envir = .GlobalEnv)
      results <- parallel::parLapply(cl, mutations, process_one)
      parallel::stopCluster(cl)
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  valid_results <- results[!sapply(results, is.null)]

  cat("\n--- Complete ---\n")
  cat("Time:", round(elapsed, 1), "seconds\n")
  cat("Successful:", length(valid_results), "/", n_mutations, "\n")

  if (length(valid_results) > 0) {
    final_df <- do.call(rbind, valid_results)
    cat("Total designs:", nrow(final_df), "\n")
    return(final_df)
  }

  data.frame()
}

# ==============================================================================
# SYSTEM INITIALISATION
# ==============================================================================

#' Initialise the design system from a genome source.
#' @param source a genome source list from strain_source()/custom_source().
#'   For backward compatibility, a bare accession string or NULL (= default
#'   strain) is also accepted.
#' @param force_download re-download cached NCBI genomes/annotations
#' @param progress_callback optional function(fraction, message) for the GUI
initialize_system_fast <- function(source = NULL, force_download = FALSE,
                                    progress_callback = NULL) {
  # Normalise the source argument
  if (is.null(source)) {
    source <- strain_source(DEFAULT_GENOME_ACCESSION, "BW25113 (default)")
  } else if (is.character(source)) {
    source <- strain_source(source, source)
  }

  cat("\n=== Initializing System (Fast Version) ===\n\n")
  cat("Genome source:", source$name, "\n")
  if (!is.null(progress_callback)) progress_callback(0.05, "Loading genome...")

  empty_gene_data <- list(sequences = list(), metadata = list(), id_maps = list(),
                          id_types = integer(0))

  if (source$type == "custom") {
    if (is.null(source$genome_file) || !file.exists(source$genome_file)) {
      stop("Custom genome FASTA not found.", call. = FALSE)
    }
    genome <- Biostrings::readDNAStringSet(source$genome_file)
    # Custom uploads: annotation from the user's GFF/GTF (via rtracklayer)
    gene_data_loader <- function() {
      ann <- load_annotations_file(source$gff_file, source$gff_format)
      if (is.null(ann)) empty_gene_data else extract_gene_sequences(genome, ann)
    }
  } else {
    genome <- download_genome_fasta(source$accession, force_download)
    # Strains: annotation from the NCBI feature table (rettype=ft)
    gene_data_loader <- function() {
      ft <- download_genome_features(source$accession, force_download)
      if (is.null(ft)) empty_gene_data else extract_gene_sequences_ft(genome, ft)
    }
  }

  if (length(genome) == 0) stop("Genome file contains no sequences.", call. = FALSE)
  if (length(genome) > 1) {
    cat("Note: genome has", length(genome),
        "sequences; using the first record as the chromosome.\n")
  }

  genome_seq <- as.character(genome[[1]])
  cat("Genome length:", format(nchar(genome_seq), big.mark = ","), "bp\n")

  # Build fast PAM index
  if (!is.null(progress_callback)) progress_callback(0.30, "Building PAM index...")
  pam_data <- build_pam_index_fast(genome_seq)

  # Load annotations
  if (!is.null(progress_callback)) progress_callback(0.70, "Loading gene annotations...")
  gene_data <- tryCatch(gene_data_loader(), error = function(e) {
    warning("Annotation loading failed: ", conditionMessage(e), call. = FALSE)
    empty_gene_data
  })

  cat("Genes loaded:", length(gene_data$sequences), "\n")
  if (length(gene_data$id_types) > 0) {
    cat("Identifier types available:",
        paste(sprintf("%s (%d)", names(gene_data$id_types), gene_data$id_types),
              collapse = ", "), "\n")
  }

  # Compute genome-specific codon usage
  if (!is.null(progress_callback)) progress_callback(0.90, "Computing codon usage...")
  codon_freq <- if (length(gene_data$sequences) > 0) {
    compute_codon_usage(gene_data$sequences)
  } else {
    ECOLI_CODON_FREQ
  }

  cat("\nSystem ready!\n")
  if (!is.null(progress_callback)) progress_callback(1.0, "Genome ready.")

  list(
    genome = genome,
    genome_seq = genome_seq,
    genome_name = source$name,
    pam_data = pam_data,
    gene_sequences = gene_data$sequences,
    gene_metadata = gene_data$metadata,
    id_maps = gene_data$id_maps,
    id_types = gene_data$id_types,
    codon_freq = codon_freq
  )
}
