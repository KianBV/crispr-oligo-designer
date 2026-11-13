# ==============================================================================
# offtarget.R  --  PAM index construction and fast off-target counting
# ==============================================================================
# Uses Biostrings (optimised C) for approximate matching. Logic preserved.
# ==============================================================================

# --- Shared helper: reverse complement ----------------------------------------
reverse_complement <- function(seq) {
  complement <- chartr("ACGT", "TGCA", toupper(seq))
  paste(rev(strsplit(complement, "")[[1]]), collapse = "")
}

#' Build optimised PAM index with DNAStringSet for fast matching
build_pam_index_fast <- function(genome_sequence) {
  cat("Building PAM index...\n")

  genome_seq <- toupper(genome_sequence)
  genome_length <- nchar(genome_seq)

  # Find PAM sites
  fwd_matches <- gregexpr("(?=.GG)", genome_seq, perl = TRUE)[[1]]
  rev_matches <- gregexpr("(?=CC.)", genome_seq, perl = TRUE)[[1]]

  cat("  Forward PAMs:", length(fwd_matches), "\n")
  cat("  Reverse PAMs:", length(rev_matches), "\n")

  # Extract protospacers - forward strand
  pam_data <- list()

  if (fwd_matches[1] != -1) {
    valid_fwd <- fwd_matches[fwd_matches > 20]
    if (length(valid_fwd) > 0) {
      protos <- substring(genome_seq, valid_fwd - 20, valid_fwd - 1)
      seeds <- substring(protos, 10, 20)
      pam_data[[1]] <- data.table::data.table(
        pam_position = valid_fwd, strand = "+", protospacer = protos, seed = seeds
      )
    }
  }

  # Reverse strand
  if (rev_matches[1] != -1) {
    valid_rev <- rev_matches[rev_matches + 22 <= genome_length]
    if (length(valid_rev) > 0) {
      protos_fwd <- substring(genome_seq, valid_rev + 3, valid_rev + 22)
      protos <- sapply(protos_fwd, reverse_complement, USE.NAMES = FALSE)
      seeds <- substring(protos, 10, 20)
      pam_data[[2]] <- data.table::data.table(
        pam_position = valid_rev, strand = "-", protospacer = protos, seed = seeds
      )
    }
  }

  pam_index <- data.table::rbindlist(pam_data)

  # Create DNAStringSet for FAST matching with Biostrings
  cat("  Creating DNAStringSet for fast matching...\n")
  proto_dna <- Biostrings::DNAStringSet(pam_index$protospacer)

  # Build seed hash for quick lookups
  cat("  Building seed index...\n")
  seed_index <- new.env(hash = TRUE)
  for (i in seq_len(nrow(pam_index))) {
    seed <- pam_index$seed[i]
    pos <- pam_index$pam_position[i]
    if (exists(seed, envir = seed_index)) {
      seed_index[[seed]] <- c(seed_index[[seed]], pos)
    } else {
      seed_index[[seed]] <- pos
    }
  }

  cat("  Total PAM sites:", nrow(pam_index), "\n")

  list(
    index = pam_index,
    proto_dna = proto_dna,
    seed_index = seed_index
  )
}

#' FAST off-target counting using Biostrings vcountPattern (optimised C)
count_off_targets_fast <- function(protospacer, pam_data, target_position,
                                   exclude_window = OT_EXCLUDE_WINDOW) {
  protospacer <- toupper(protospacer)
  seed <- substr(protospacer, 10, 20)

  pam_index <- pam_data$index
  seed_index <- pam_data$seed_index
  proto_dna <- pam_data$proto_dna

  # Quick seed match check (hash lookup - O(1)). Exclude the on-target footprint
  # (|pos - target| <= exclude_window): the on-target and its frameshifted
  # registers share the seed at a shifted position but are not off-targets.
  seed_matches <- 0
  if (exists(seed, envir = seed_index)) {
    matching_positions <- seed_index[[seed]]
    seed_matches <- sum(abs(matching_positions - target_position) > exclude_window)
  }

  # Use Biostrings for fast approximate matching
  query <- Biostrings::DNAString(protospacer)

  off_target_count <- 0

  # Count all with 1-4 mismatches
  matches_le4 <- Biostrings::vcountPattern(
    query, proto_dna,
    max.mismatch = MAX_MISMATCHES,
    min.mismatch = 1,
    with.indels = FALSE
  )

  # Find which ones matched
  match_indices <- which(matches_le4 > 0)

  # Exclude the on-target footprint (on-target + frameshifted self-overlaps),
  # not just the exact on-target position, so the count matches the inspector's
  # enumeration (find_off_targets uses the same window).
  if (length(match_indices) > 0) {
    mpos <- pam_index$pam_position[match_indices]
    off_target_count <- sum(abs(mpos - target_position) > exclude_window)
  }

  # Determine class
  has_seed_matches <- seed_matches > 0
  has_off_targets <- off_target_count > 0

  if (!has_seed_matches && !has_off_targets) {
    off_target_class <- 1
  } else if (!has_seed_matches && has_off_targets) {
    off_target_class <- 2
  } else if (has_seed_matches && !has_off_targets) {
    off_target_class <- 3
  } else {
    off_target_class <- 4
  }

  list(
    seed_matches = seed_matches,
    off_target_count = off_target_count,
    off_target_class = off_target_class
  )
}

#' Enumerate the actual off-target sites for one protospacer (on-demand; for the
#' design inspector, not the hot design loop). Returns a data.frame of hits with
#' genomic position, strand, mismatch count, whether the seed is shared, the
#' off-target protospacer + PAM, and the gene it falls in.
#'
#' Includes both <= max_mismatch full-length hits ("mismatch") and exact
#' PAM-proximal seed matches with more distal differences ("seed-only").
find_off_targets <- function(protospacer, pam_data, target_position,
                             genome_seq = NULL, gene_metadata = NULL,
                             max_mismatch = MAX_MISMATCHES, max_report = 500,
                             exclude_window = OT_EXCLUDE_WINDOW) {
  protospacer <- toupper(protospacer)
  pam_index <- pam_data$index
  proto_dna <- pam_data$proto_dna
  qchars <- strsplit(protospacer, "")[[1]]
  qseed  <- substr(protospacer, 10, 20)

  # Candidate rows = union of (<=max_mismatch full match) and (exact seed match).
  # Exclude anything within the on-target's own footprint (|pos - target| <=
  # exclude_window): those are the on-target itself and its 1-bp-shifted frames,
  # which share the seed at a frameshifted register and are not real off-targets.
  query  <- Biostrings::DNAString(protospacer)
  mm_hit <- Biostrings::vcountPattern(query, proto_dna, max.mismatch = max_mismatch,
                                      min.mismatch = 0, with.indels = FALSE) > 0
  seed_hit <- pam_index$seed == qseed
  keep <- which((mm_hit | seed_hit) &
                abs(pam_index$pam_position - target_position) > exclude_window)
  if (length(keep) == 0) return(data.frame())

  sub <- pam_index[keep, ]
  protos <- sub$protospacer
  mm <- vapply(protos, function(s) sum(strsplit(s, "")[[1]] != qchars), integer(1))
  seedmatch <- substr(protos, 10, 20) == qseed

  df <- data.frame(
    type        = ifelse(mm <= max_mismatch, "mismatch", "seed-only"),
    position    = sub$pam_position,
    strand      = sub$strand,
    mismatches  = mm,
    seed_match  = seedmatch,
    protospacer = protos,
    stringsAsFactors = FALSE
  )

  # PAM as it reads on the guide strand (NGG)
  if (!is.null(genome_seq)) {
    gl <- nchar(genome_seq)
    df$pam <- mapply(function(pos, st) {
      if (pos < 1 || pos + 2 > gl) return(NA_character_)
      tri <- substr(genome_seq, pos, pos + 2)
      if (identical(st, "-")) reverse_complement(tri) else tri
    }, df$position, df$strand, USE.NAMES = FALSE)
  }

  # Gene containing each off-target (genomic ranges)
  if (!is.null(gene_metadata) && length(gene_metadata) > 0) {
    gn <- names(gene_metadata)
    gs <- vapply(gene_metadata, function(m) as.numeric(m$genomic_start), numeric(1))
    ge <- vapply(gene_metadata, function(m) as.numeric(m$genomic_end), numeric(1))
    df$gene <- vapply(df$position, function(p) {
      hit <- which(gs <= p & ge >= p)
      if (length(hit)) gn[hit[1]] else NA_character_
    }, character(1))
  }

  df <- df[order(df$mismatches, !df$seed_match), , drop = FALSE]  # worst first
  if (nrow(df) > max_report) df <- df[seq_len(max_report), , drop = FALSE]
  rownames(df) <- NULL
  df
}
