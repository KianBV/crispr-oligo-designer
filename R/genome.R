# ==============================================================================
# genome.R  --  Genome/annotation download, caching and gene extraction
# ==============================================================================
# Logic preserved verbatim from the original pipeline; only reorganised.
# ==============================================================================

ensure_cache_dir <- function() {
  if (!dir.exists(GENOME_CACHE_DIR)) dir.create(GENOME_CACHE_DIR, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Genome source constructors
# ------------------------------------------------------------------------------
# A "source" describes where the genome comes from. Two kinds:
#   strain_source(accession, name)          -> downloads from NCBI by accession
#   custom_source(genome_file, gff_file, .)  -> user-uploaded FASTA (+ optional GFF)

strain_source <- function(accession = DEFAULT_GENOME_ACCESSION, name = accession) {
  list(type = "strain", accession = accession, name = name)
}

custom_source <- function(genome_file, gff_file = NULL, gff_format = NULL,
                          name = "Custom genome") {
  list(type = "custom", genome_file = genome_file, gff_file = gff_file,
       gff_format = gff_format, name = name)
}

# ------------------------------------------------------------------------------
# Downloading by accession (cached per-accession)
# ------------------------------------------------------------------------------

#' Download a URL to a cache file robustly: writes to a temp file first and only
#' commits to the cache if the download succeeded and the content passes a
#' validity check. This prevents a failed download (e.g. an HTTP 500 error page)
#' from being cached as a valid file and poisoning later loads.
#' @param validator function(path) -> TRUE if the downloaded content looks valid
#' @return the cache file path on success, or NULL on failure
fetch_to_cache <- function(url, cache_file, force_download = FALSE,
                           validator = function(p) file.info(p)$size > 0,
                           label = "file") {
  ensure_cache_dir()

  if (!force_download && file.exists(cache_file) &&
      isTRUE(tryCatch(validator(cache_file), error = function(e) FALSE))) {
    cat("Loading cached", label, "from:", cache_file, "\n")
    return(cache_file)
  }
  # A cached but invalid file (e.g. a previously poisoned 0-byte file) is dropped
  if (file.exists(cache_file) && !isTRUE(tryCatch(validator(cache_file),
                                                  error = function(e) FALSE))) {
    unlink(cache_file)
  }

  cat("Downloading", label, "...\n")
  tmp <- tempfile()
  ok <- tryCatch({
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    file.exists(tmp) && isTRUE(tryCatch(validator(tmp), error = function(e) FALSE))
  }, error = function(e) FALSE, warning = function(w) {
    file.exists(tmp) && isTRUE(tryCatch(validator(tmp), error = function(e) FALSE))
  })

  if (!ok) {
    if (file.exists(tmp)) unlink(tmp)
    warning("Download failed or returned invalid content for ", label,
            " (", url, ")", call. = FALSE)
    return(NULL)
  }
  file.copy(tmp, cache_file, overwrite = TRUE)
  unlink(tmp)
  cache_file
}

# Validators
.valid_fasta <- function(p) {
  if (file.info(p)$size == 0) return(FALSE)
  first <- readLines(p, n = 1, warn = FALSE)
  length(first) == 1 && startsWith(first, ">")
}
.valid_ft <- function(p) {
  if (file.info(p)$size == 0) return(FALSE)
  first <- readLines(p, n = 1, warn = FALSE)
  length(first) == 1 && grepl("^>Feature", first)
}

download_genome_fasta <- function(accession, force_download = FALSE) {
  cache_file <- file.path(GENOME_CACHE_DIR, paste0(accession, ".fasta"))
  path <- fetch_to_cache(ncbi_fasta_url(accession), cache_file, force_download,
                         validator = .valid_fasta,
                         label = paste0("genome ", accession))
  if (is.null(path)) {
    stop("Could not download genome sequence for '", accession,
         "' from NCBI. Check the accession and your internet connection.",
         call. = FALSE)
  }
  Biostrings::readDNAStringSet(path)
}

#' Download the NCBI feature table (annotation) for an accession.
#' @return path to the cached .ft file, or NULL if it could not be fetched.
download_genome_features <- function(accession, force_download = FALSE) {
  cache_file <- file.path(GENOME_CACHE_DIR, paste0(accession, ".ft"))
  fetch_to_cache(ncbi_ft_url(accession), cache_file, force_download,
                 validator = .valid_ft,
                 label = paste0("annotation (feature table) ", accession))
}

# Load annotations from a user-supplied GFF/GTF file. `format` (e.g. "gff3",
# "gtf") is passed to rtracklayer because the uploaded temp file may not keep a
# meaningful extension.
load_annotations_file <- function(gff_file, format = NULL) {
  if (is.null(gff_file) || !file.exists(gff_file)) return(NULL)
  tryCatch(
    if (is.null(format) || !nzchar(format)) {
      rtracklayer::import(gff_file)
    } else {
      rtracklayer::import(gff_file, format = format)
    },
    error = function(e) NULL
  )
}

# Extract UniProt accessions from a set of db_xref values.
.parse_uniprot <- function(dbxref_val) {
  if (is.null(dbxref_val) || length(dbxref_val) == 0) return(character(0))
  x <- as.character(dbxref_val)
  hits <- grep("UniProt", x, ignore.case = TRUE, value = TRUE)
  unique(sub(".*:", "", hits))   # part after the last colon
}

# Extract protein accession(s) from a (possibly pipe-delimited) protein_id value,
# e.g. "gb|AAC73112.1||gnl|b0001|NP_414542" -> c("AAC73112.1", "NP_414542").
.parse_protein_ids <- function(raw) {
  if (is.null(raw) || length(raw) == 0) return(character(0))
  toks <- unlist(strsplit(as.character(raw), "|", fixed = TRUE))
  toks <- trimws(toks); toks <- toks[nzchar(toks)]
  keep <- grepl("^[A-Z]{2,3}[0-9]{5,}(\\.[0-9]+)?$", toks) |   # GenBank e.g. AAC73112.1
          grepl("^[A-Z]{2}_[0-9]{4,}(\\.[0-9]+)?$", toks)      # RefSeq  e.g. NP_414542
  unique(toks[keep])
}

#' Build the gene index (sequences, metadata, identifier maps) from a list of
#' gene records. Shared by the GFF and feature-table extractors.
#' @param records list of list(key, start, end, strand, ids=named list per type)
assemble_gene_index <- function(genome, records) {
  genome_seq <- genome[[1]]
  glen <- length(genome_seq)
  gene_sequences <- list()
  gene_metadata <- list()

  id_types <- c("gene", "locus_tag", "old_locus_tag", "protein_id", "uniprot")
  id_maps <- setNames(lapply(id_types, function(x) list()), id_types)

  for (r in records) {
    key <- r$key
    if (is.null(key) || length(key) == 0 || is.na(key) || !nzchar(key)) next
    lo <- min(r$start, r$end); hi <- max(r$start, r$end)
    if (is.na(lo) || is.na(hi) || lo < 1 || hi > glen) next

    seq <- Biostrings::subseq(genome_seq, lo, hi)
    if (identical(r$strand, "-")) seq <- Biostrings::reverseComplement(seq)

    gene_sequences[[key]] <- as.character(seq)
    gene_metadata[[key]] <- list(genomic_start = lo, genomic_end = hi, strand = r$strand)

    for (t in id_types) {
      for (v in r$ids[[t]]) {
        v <- trimws(v)
        if (!is.na(v) && nzchar(v)) id_maps[[t]][[v]] <- key
      }
    }
  }

  id_maps <- id_maps[vapply(id_maps, length, integer(1)) > 0]
  list(sequences = gene_sequences, metadata = gene_metadata,
       id_maps = id_maps, id_types = vapply(id_maps, length, integer(1)))
}

#' Extract gene sequences + identifier maps from an rtracklayer GFF/GTF import.
#' Used for user-uploaded annotation files.
extract_gene_sequences <- function(genome, annotations) {
  cds <- annotations[annotations$type == "CDS"]
  if (length(cds) == 0) cds <- annotations[annotations$type == "gene"]

  md <- S4Vectors::mcols(cds)
  col <- function(name) if (name %in% names(md)) md[[name]] else NULL
  gene_v   <- col("gene")
  name_v   <- col("Name")
  locus_v  <- col("locus_tag")
  oldloc_v <- col("old_locus_tag")
  prot_v   <- col("protein_id")
  dbx_v    <- col("Dbxref")

  starts  <- BiocGenerics::start(cds)
  ends    <- BiocGenerics::end(cds)
  strands <- as.character(BiocGenerics::strand(cds))

  at <- function(v, i) {                 # safe scalar accessor (handles NULL/NA)
    if (is.null(v)) return(NA_character_)
    x <- as.character(v[[i]] %||% v[i])
    if (length(x) == 0) NA_character_ else x[1]
  }

  records <- vector("list", length(cds))
  for (i in seq_along(cds)) {
    key <- at(gene_v, i)
    if (is.na(key)) key <- at(name_v, i)
    if (is.na(key)) key <- at(locus_v, i)
    if (is.na(key)) next
    records[[i]] <- list(
      key = key, start = starts[i], end = ends[i], strand = strands[i],
      ids = list(
        gene          = at(gene_v, i),
        locus_tag     = at(locus_v, i),
        old_locus_tag = at(oldloc_v, i),
        protein_id    = at(prot_v, i),
        uniprot       = if (!is.null(dbx_v)) .parse_uniprot(dbx_v[[i]]) else character(0)
      )
    )
  }
  records <- records[!vapply(records, is.null, logical(1))]
  assemble_gene_index(genome, records)
}

#' Parse an NCBI feature table (efetch rettype=ft) into gene records.
#' Each CDS is joined with its preceding gene feature (which carries gene/locus).
parse_feature_table <- function(ft_path) {
  lines <- readLines(ft_path, warn = FALSE)
  pc <- function(x) suppressWarnings(as.integer(gsub("[^0-9]", "", x)))  # parse coord

  records <- list()
  cur <- NULL          # feature currently being read
  last_gene <- NULL    # most recently completed 'gene' feature

  finalize <- function() {
    if (is.null(cur)) return(invisible())
    if (identical(cur$type, "gene")) {
      last_gene <<- cur
    } else if (identical(cur$type, "CDS")) {
      q  <- cur$quals
      lg <- if (!is.null(last_gene)) last_gene$quals else list()
      gene   <- q$gene          %||% lg$gene
      locus  <- q$locus_tag     %||% lg$locus_tag
      oldloc <- q$old_locus_tag %||% lg$old_locus_tag
      prot   <- .parse_protein_ids(q$protein_id)
      uni    <- .parse_uniprot(q$db_xref)
      key    <- (gene %||% locus %||% prot)[1]
      records[[length(records) + 1]] <<- list(
        key = key, start = cur$lo, end = cur$hi, strand = cur$strand,
        ids = list(gene = gene, locus_tag = locus, old_locus_tag = oldloc,
                   protein_id = prot, uniprot = uni)
      )
    }
  }

  for (ln in lines) {
    if (!nzchar(ln) || startsWith(ln, ">")) next
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]

    if (nzchar(f[1])) {
      # feature interval line: start, end, [type]
      s <- pc(f[1]); e <- pc(f[2])
      type <- if (length(f) >= 3) f[3] else ""
      if (nzchar(type)) {
        finalize()
        cur <- list(type = type, strand = if (!is.na(s) && !is.na(e) && s > e) "-" else "+",
                    lo = min(s, e), hi = max(s, e), quals = list())
      } else if (!is.null(cur)) {
        # continuation interval of the current feature
        cur$lo <- min(cur$lo, s, e, na.rm = TRUE)
        cur$hi <- max(cur$hi, s, e, na.rm = TRUE)
      }
    } else {
      # qualifier line: "", "", "", key, [value]
      if (!is.null(cur) && length(f) >= 4 && nzchar(f[4])) {
        key <- f[4]; value <- if (length(f) >= 5) f[5] else ""
        cur$quals[[key]] <- c(cur$quals[[key]], value)
      }
    }
  }
  finalize()  # close the last feature
  records
}

#' Extract gene sequences + identifier maps from an NCBI feature-table file.
extract_gene_sequences_ft <- function(genome, ft_path) {
  assemble_gene_index(genome, parse_feature_table(ft_path))
}

# ------------------------------------------------------------------------------
# Identifier resolution
# ------------------------------------------------------------------------------
# Human-readable labels for identifier types.
ID_TYPE_LABELS <- c(
  gene          = "Gene name",
  locus_tag     = "Locus tag",
  old_locus_tag = "Old locus tag",
  protein_id    = "Protein ID",
  uniprot       = "UniProt ID"
)

#' Resolve a user-supplied identifier to a gene primary key.
#' @param id the identifier as typed by the user
#' @param id_maps the system's id_maps
#' @param id_type "auto" (search all, in priority order) or a specific type
#' @return the primary key, or NULL if not found
resolve_gene_id <- function(id, id_maps, id_type = "gene") {
  if (is.null(id_maps) || length(id_maps) == 0) return(NULL)
  id <- trimws(id)
  types <- if (identical(id_type, "auto") || is.null(id_type)) {
    intersect(c("gene", "locus_tag", "old_locus_tag", "protein_id", "uniprot"),
              names(id_maps))
  } else {
    id_type
  }
  for (t in types) {
    m <- id_maps[[t]]
    if (is.null(m)) next
    if (!is.null(m[[id]])) return(m[[id]])
    # case-insensitive fallback (mainly for gene symbols)
    idx <- match(tolower(id), tolower(names(m)))
    if (!is.na(idx)) return(m[[idx]])
  }
  NULL
}

# Codon usage from genome
compute_codon_usage <- function(gene_sequences) {
  cat("Computing codon usage from", length(gene_sequences), "coding sequences...\n")

  all_codons <- character(0)
  for (seq in gene_sequences) {
    seq <- toupper(seq)
    n_codons <- nchar(seq) %/% 3
    if (n_codons < 1) next
    starts <- seq(1, n_codons * 3, by = 3)
    ends <- seq(3, n_codons * 3, by = 3)
    codons <- substring(seq, starts, ends)
    all_codons <- c(all_codons, codons)
  }

  # Raw counts
  codon_counts <- table(factor(all_codons, levels = names(REVERSE_CODON_TABLE)))

  # Convert to relative frequency within each amino acid family
  freq <- numeric(length(codon_counts))
  names(freq) <- names(codon_counts)

  for (aa in names(CODON_TABLE)) {
    codons_for_aa <- CODON_TABLE[[aa]]
    aa_total <- sum(codon_counts[codons_for_aa])
    if (aa_total > 0) {
      for (cdn in codons_for_aa) {
        freq[cdn] <- as.numeric(codon_counts[cdn]) / aa_total
      }
    }
  }

  cat("  Counted", format(length(all_codons), big.mark = ","), "codons\n")
  freq
}
