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
count_off_targets_fast <- function(protospacer, pam_data, target_position) {
  protospacer <- toupper(protospacer)
  seed <- substr(protospacer, 10, 20)

  pam_index <- pam_data$index
  seed_index <- pam_data$seed_index
  proto_dna <- pam_data$proto_dna

  # Quick seed match check (hash lookup - O(1))
  seed_matches <- 0
  if (exists(seed, envir = seed_index)) {
    matching_positions <- seed_index[[seed]]
    seed_matches <- sum(matching_positions != target_position)
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

  # Exclude on-target
  if (length(match_indices) > 0) {
    on_target_idx <- which(pam_index$pam_position == target_position)
    match_indices <- setdiff(match_indices, on_target_idx)
    off_target_count <- length(match_indices)
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
