# ==============================================================================
# qc.R  --  Oligo quality-control metrics (pure R, no external dependencies)
# ==============================================================================
# Lightweight, dependency-free checks for the ordered oligos:
#   * GC content
#   * Nearest-neighbor melting temperature (SantaLucia 1998 unified parameters)
#   * Longest hairpin stem (self-complementary, heuristic)
#   * 3'-end self-dimer (heuristic)
#   * Longest homopolymer run
#   * Internal type-IIS restriction sites (BsaI / BbsI) -- critical for Golden Gate
#
# The hairpin / dimer values are heuristics (longest complementary stem), NOT
# true thermodynamic Delta-G; they flag the obvious offenders. A hook is left
# for an accurate primer3 / ViennaRNA backend (see qc_oligo()).
# ==============================================================================

# Type-IIS recognition sites (both strands checked by searching seq + revcomp).
ENZYME_SITES <- list(
  BsaI = "GGTCTC",   # revcomp GAGACC
  BbsI = "GAAGAC"    # revcomp GTCTTC  (a.k.a. BpiI)
)

# --- GC content (fraction 0..1) ----------------------------------------------
gc_fraction <- function(seq) {
  chars <- strsplit(toupper(seq), "")[[1]]
  if (length(chars) == 0) return(NA_real_)
  sum(chars %in% c("G", "C")) / length(chars)
}

# --- Nearest-neighbor Tm (SantaLucia 1998 unified NN parameters) --------------
# ΔH (kcal/mol) and ΔS (cal/(mol·K)) for the 10 unique NN pairs.
.NN_H <- c(AA=-7.9, AT=-7.2, TA=-7.2, CA=-8.5, GT=-8.4, CT=-7.8,
           GA=-8.2, CG=-10.6, GC=-9.8, GG=-8.0)
.NN_S <- c(AA=-22.2, AT=-20.4, TA=-21.3, CA=-22.7, GT=-22.4, CT=-21.0,
           GA=-22.2, CG=-27.2, GC=-24.4, GG=-19.9)

.nn_lookup <- function(pair) {
  # map a dinucleotide to its canonical representative (NN symmetry)
  if (!is.na(.NN_H[pair])) return(pair)
  rc <- chartr("ACGT", "TGCA", paste(rev(strsplit(pair, "")[[1]]), collapse = ""))
  rc
}

#' Nearest-neighbor melting temperature (°C).
#' @param seq oligo sequence (ACGT; other characters ignored for NN steps)
#' @param na_mM monovalent cation concentration (mM). Fold Mg2+ in here if wanted.
#' @param oligo_uM total strand concentration (µM)
tm_nn <- function(seq, na_mM = 50, oligo_uM = 0.5) {
  s <- gsub("[^ACGT]", "", toupper(seq))   # uppercase FIRST (tails are lowercase)
  n <- nchar(s)
  if (n < 8) return(NA_real_)  # NN model unreliable for very short oligos

  chars <- strsplit(s, "")[[1]]
  dH <- 0; dS <- 0
  for (i in seq_len(n - 1)) {
    step <- .nn_lookup(paste0(chars[i], chars[i + 1]))
    dH <- dH + .NN_H[step]
    dS <- dS + .NN_S[step]
  }
  # Initiation terms (per SantaLucia 1998): terminal A·T vs G·C
  for (end in c(chars[1], chars[n])) {
    if (end %in% c("G", "C")) { dH <- dH + 0.1;  dS <- dS - 2.8 }
    else                      { dH <- dH + 2.3;  dS <- dS + 4.1 }
  }
  dH <- as.numeric(dH); dS <- as.numeric(dS)

  # Salt correction (SantaLucia 1998): ΔS += 0.368*(N-1)*ln[Na+]
  na <- max(na_mM, 1e-6) / 1000  # mol/L
  dS <- dS + 0.368 * (n - 1) * log(na)

  R <- 1.987  # cal/(mol·K)
  ct <- (oligo_uM * 1e-6)         # mol/L total strand
  # non-self-complementary duplex: use CT/4
  tm_k <- (dH * 1000) / (dS + R * log(ct / 4))
  round(tm_k - 273.15, 1)
}

# --- Longest hairpin stem (heuristic) ----------------------------------------
#' Longest self-complementary stem that could close a hairpin with a loop of at
#' least `min_loop` unpaired bases. Returns the stem length in bp.
longest_hairpin_stem <- function(seq, min_loop = 3) {
  s <- strsplit(gsub("[^ACGT]", "", toupper(seq)), "")[[1]]
  n <- length(s); if (n < 2 * 1 + min_loop) return(0L)
  comp <- c(A = "T", T = "A", G = "C", C = "G")
  best <- 0L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (j - i - 1 < min_loop) next          # loop too small
      L <- 0L; a <- i; b <- j                  # extend stem outward from (i,j)
      while (a >= 1 && b <= n) {
        ca <- s[a]; cb <- s[b]
        if (!is.na(comp[ca]) && comp[ca] == cb) { L <- L + 1L; a <- a - 1L; b <- b + 1L }
        else break
      }
      if (L > best) best <- L
    }
  }
  best
}

# --- 3'-end dimer (heuristic) ------------------------------------------------
#' Longest run of consecutive complementary pairs anchored at the 3' end of `a`
#' against `b` (any antiparallel register). For a self-dimer, pass b = a.
three_prime_dimer <- function(a, b = a) {
  av <- strsplit(gsub("[^ACGT]", "", toupper(a)), "")[[1]]
  bv <- strsplit(gsub("[^ACGT]", "", toupper(b)), "")[[1]]
  na <- length(av); nb <- length(bv)
  if (na == 0 || nb == 0) return(0L)
  comp <- c(A = "T", T = "A", G = "C", C = "G")
  best <- 0L
  for (sft in 0:(nb - 1)) {
    L <- 0L
    while (L < na && (sft + L) < nb) {
      ca <- av[na - L]; cb <- bv[sft + 1 + L]
      if (!is.na(comp[ca]) && comp[ca] == cb) L <- L + 1L else break
    }
    if (L > best) best <- L
  }
  best
}

# --- Longest homopolymer run -------------------------------------------------
longest_homopolymer <- function(seq) {
  s <- gsub("[^ACGTN]", "", toupper(seq))
  if (nchar(s) == 0) return(0L)
  runs <- nchar(unlist(regmatches(s, gregexpr("(.)\\1*", s))))
  if (length(runs) == 0) 0L else max(runs)
}

# --- Internal restriction site? ----------------------------------------------
has_enzyme_site <- function(seq, enzyme) {
  site <- ENZYME_SITES[[enzyme]]
  if (is.null(site)) return(NA)
  s <- toupper(seq)
  rc <- chartr("ACGT", "TGCA", paste(rev(strsplit(site, "")[[1]]), collapse = ""))
  grepl(site, s, fixed = TRUE) || grepl(rc, s, fixed = TRUE)
}

# ==============================================================================
# primer3 backend (optional, accurate thermodynamic Delta-G)
# ==============================================================================
# Shells out to tools/primer3_calc.py (primer3-py). Returns NULL on any failure
# so callers transparently fall back to the built-in heuristics.

.primer3_script <- function(qc_cfg = list()) {
  qc_cfg$primer3_script %||% file.path("tools", "primer3_calc.py")
}
.primer3_python <- function(qc_cfg = list()) {
  qc_cfg$python %||% "python"
}

#' Is the primer3 backend usable? Probes by running the helper on one sequence
#' (avoids fragile cross-platform quoting of a `python -c` one-liner).
primer3_available <- function(qc_cfg = list()) {
  backend <- qc_cfg$thermo_backend %||% "auto"
  if (identical(backend, "builtin") || identical(backend, "off")) return(FALSE)
  if (!file.exists(.primer3_script(qc_cfg))) return(FALSE)
  res <- tryCatch(primer3_run("GAACGTTGTGGATGTCCAGC", qc_cfg), error = function(e) NULL)
  !is.null(res) && nrow(res) == 1 && !is.na(res$tm[1])
}

#' Run the primer3 batch over `seqs`. Returns a data.frame with columns
#' tm, hp_dg, hp_tm, dimer_dg (one row per seq), or NULL on failure.
primer3_run <- function(seqs, qc_cfg = list()) {
  script <- .primer3_script(qc_cfg)
  if (length(seqs) == 0 || !file.exists(script)) return(NULL)

  seqfile <- tempfile(fileext = ".txt")
  on.exit(unlink(seqfile), add = TRUE)
  writeLines(gsub("[^ACGT]", "", toupper(seqs)), seqfile)

  args <- c(script, seqfile,
            qc_cfg$na_mM %||% 50, qc_cfg$mg_mM %||% 1.5, qc_cfg$dntp_mM %||% 0.6,
            (qc_cfg$oligo_uM %||% 0.5) * 1000, qc_cfg$thermo_temp %||% 37)
  out <- tryCatch(
    suppressWarnings(system2(.primer3_python(qc_cfg), as.character(args),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out) || length(out) < 2) return(NULL)
  df <- tryCatch(utils::read.csv(text = paste(out, collapse = "\n"),
                                 stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) != length(seqs) ||
      !all(c("tm", "hp_dg", "hp_tm", "dimer_dg") %in% names(df))) return(NULL)
  df
}

# ==============================================================================
# Per-oligo QC row + table
# ==============================================================================

#' QC a single oligo. `seq` may contain PTO asterisks (stripped) but should be
#' the plain nucleotide sequence for meaningful metrics.
#' @param role "primer" (Tm/dimer flags apply) or "donor" (ssDNA repair oligo:
#'   Tm of a ~100-mer is inherently high and not meaningful, so those flags are
#'   skipped; hairpin/GC/homopolymer still apply).
#' @param thermo optional named list(tm, hp_dg, hp_tm, dimer_dg) from primer3. When
#'   supplied, Tm and the hairpin/dimer flags use these accurate values instead of
#'   the built-in heuristics.
qc_oligo <- function(name, seq, qc_cfg = list(), role = "primer", thermo = NULL) {
  g   <- qc_cfg
  na  <- g$na_mM %||% 50
  oli <- g$oligo_uM %||% 0.5
  enzymes <- g$enzymes %||% c("BsaI", "BbsI")
  use_p3 <- !is.null(thermo)

  clean <- gsub("[^ACGT]", "", toupper(seq))   # uppercase FIRST, then strip PTO '*'/mods
  len <- nchar(clean)
  gc  <- gc_fraction(clean)
  hp  <- longest_hairpin_stem(clean, min_loop = g$hairpin_min_loop %||% 3)
  d3  <- three_prime_dimer(clean, clean)
  hpz <- longest_homopolymer(clean)

  tm        <- if (use_p3 && !is.na(thermo$tm)) round(thermo$tm, 1) else tm_nn(clean, na_mM = na, oligo_uM = oli)
  tm_method <- if (use_p3) "primer3" else "NN"
  hp_dg     <- if (use_p3) thermo$hp_dg else NA_real_
  hp_tm     <- if (use_p3) thermo$hp_tm else NA_real_
  dimer_dg  <- if (use_p3) thermo$dimer_dg else NA_real_

  row <- data.frame(
    oligo_name = name, role = role, length = len,
    gc_pct = round(100 * gc, 1), tm = tm, tm_method = tm_method,
    hairpin_stem = hp, hairpin_dG = hp_dg, hairpin_tm = hp_tm,
    self_dimer3 = d3, dimer_dG = dimer_dg, max_run = hpz,
    stringsAsFactors = FALSE
  )
  for (e in enzymes) row[[paste0(e, "_site")]] <- isTRUE(has_enzyme_site(clean, e))

  # --- flags ---
  flags <- character(0)
  gc_min <- g$gc_min %||% 0.30; gc_max <- g$gc_max %||% 0.65
  if (!is.na(gc) && gc < gc_min) flags <- c(flags, sprintf("low GC (%.0f%%)", 100 * gc))
  if (!is.na(gc) && gc > gc_max) flags <- c(flags, sprintf("high GC (%.0f%%)", 100 * gc))
  if (identical(role, "primer")) {   # Tm / dimer only meaningful for primers
    if (!is.na(tm) && !is.na(g$tm_min) && tm < g$tm_min) flags <- c(flags, sprintf("low Tm (%.0f)", tm))
    if (!is.na(tm) && !is.na(g$tm_max) && tm > g$tm_max) flags <- c(flags, sprintf("high Tm (%.0f)", tm))
    if (use_p3) {
      if (!is.na(dimer_dg) && dimer_dg <= (g$dimer_dg %||% -6)) flags <- c(flags, sprintf("self-dimer %.1f kcal/mol", dimer_dg))
    } else if (d3 >= (g$dimer3 %||% 4)) {
      flags <- c(flags, sprintf("3' self-dimer %dbp", d3))
    }
  }
  # hairpin flag: prefer primer3 Delta-G, else heuristic stem length
  if (use_p3) {
    if (!is.na(hp_dg) && hp_dg <= (g$hairpin_dg %||% -3)) flags <- c(flags, sprintf("hairpin %.1f kcal/mol", hp_dg))
  } else if (hp >= (g$hairpin_stem %||% 6)) {
    flags <- c(flags, sprintf("hairpin %dbp", hp))
  }
  if (hpz >= (g$max_homopolymer %||% 5)) flags <- c(flags, sprintf("%dx homopolymer", hpz))
  for (e in enzymes) if (isTRUE(row[[paste0(e, "_site")]])) flags <- c(flags, paste0(e, " site"))

  row$flags <- if (length(flags)) paste(flags, collapse = "; ") else ""
  row
}

#' QC every oligo in an order table (data.frame with oligo_name, sequence).
#' Oligos whose name ends with `ha_suffix` are treated as donors (role = donor).
#' When the primer3 backend is available it is batch-run once over all oligos and
#' its accurate Tm / Delta-G values are used; otherwise the built-in heuristics are.
qc_table <- function(oligo_table, qc_cfg = list(), ha_suffix = "_HA") {
  if (is.null(oligo_table) || nrow(oligo_table) == 0) {
    return(data.frame())
  }

  # Accurate thermodynamics in one batch call (NULL -> heuristic fallback)
  p3 <- NULL
  backend <- qc_cfg$thermo_backend %||% "auto"
  if (!backend %in% c("builtin", "off")) {
    p3 <- tryCatch(primer3_run(oligo_table$sequence, qc_cfg), error = function(e) NULL)
  }

  rows <- lapply(seq_len(nrow(oligo_table)), function(i) {
    nm <- oligo_table$oligo_name[i]
    role <- if (endsWith(nm, ha_suffix)) "donor" else "primer"
    thermo <- if (!is.null(p3)) list(tm = p3$tm[i], hp_dg = p3$hp_dg[i],
                                     hp_tm = p3$hp_tm[i], dimer_dg = p3$dimer_dg[i]) else NULL
    qc_oligo(nm, oligo_table$sequence[i], qc_cfg, role = role, thermo = thermo)
  })
  do.call(rbind, rows)
}
