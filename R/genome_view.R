# ==============================================================================
# genome_view.R  --  SnapGene-like genomic view of a selected design (Gviz)
# ==============================================================================
# Renders, for one selected design, a gene-local track view showing the
# protospacer, PAM, homology arm (donor), the edited codon and the silent
# PAM-blocking mutation, over the gene's sequence and annotation.
#
# Coordinates are GENE-LOCAL (the CDS in reading-frame orientation), which is
# robust for minus-strand genes and reads like a SnapGene feature map. Gviz is
# an optional dependency: render_design_view() errors cleanly if it is missing.
# ==============================================================================

FEATURE_COLORS <- c(
  protospacer = "#2C7FB8",  # blue
  PAM         = "#D95F0E",  # orange
  HA          = "#31A354",  # green
  edit        = "#E34A33",  # red
  silentPAM   = "#756BB1"   # purple
)

#' Is Gviz available for the genomic view?
gviz_available <- function() requireNamespace("Gviz", quietly = TRUE)

#' Render the design inspector genome view. Call inside a plotting device
#' (e.g. Shiny renderPlot). Returns invisibly.
#' @param design_row one row (data.frame) from the results table
#' @param system the initialised system (for the gene sequence)
#' @param pad bp of flanking context around the features
render_design_view <- function(design_row, system, pad = 18) {
  if (!gviz_available()) {
    stop("Gviz is not installed. Run:  BiocManager::install('Gviz')", call. = FALSE)
  }
  old <- options(ucscChromosomeNames = FALSE)
  on.exit(options(old), add = TRUE)

  # Resolve the gene identifier to the sequence key
  key <- tryCatch(resolve_gene_id(design_row$gene, system$id_maps, "auto"),
                  error = function(e) NULL)
  if (is.null(key)) key <- as.character(design_row$gene)
  gene_seq <- system$gene_sequences[[key]]
  if (is.null(gene_seq)) {
    stop("Gene sequence for '", design_row$gene, "' not found in the loaded genome.",
         call. = FALSE)
  }
  chr  <- as.character(key)
  glen <- nchar(gene_seq)

  dss <- Biostrings::DNAStringSet(toupper(gene_seq)); names(dss) <- chr

  cs <- as.integer(design_row$codon_start)
  pamp <- as.integer(design_row$pam_position)
  ps <- as.integer(design_row$proto_start); pe <- as.integer(design_row$proto_end)
  feats <- data.frame(
    start   = c(min(ps, pe), pamp,        as.integer(design_row$ha_start), cs,     as.integer(design_row$silent_mut_pos)),
    end     = c(max(ps, pe), pamp + 2,    as.integer(design_row$ha_end),   cs + 2, as.integer(design_row$silent_mut_pos)),
    feature = c("protospacer", "PAM", "HA", "edit", "silentPAM"),
    label   = c("protospacer", "PAM", "HA (donor)",
                paste0("edit ", design_row$original_codon, "→", design_row$mutant_codon),
                "silent PAM mut"),
    stringsAsFactors = FALSE
  )
  feats$start <- pmax(1L, feats$start)
  feats$end   <- pmin(glen, feats$end)

  strand_vec <- c(design_row$strand, design_row$strand, "*", "*", "*")
  atrack <- Gviz::AnnotationTrack(
    start = feats$start, end = feats$end, chromosome = chr, strand = strand_vec,
    feature = feats$feature, group = feats$label, id = feats$label,
    name = "design", stacking = "full", groupAnnotation = "group",
    just.group = "below", fontcolor.group = "black", cex.group = 0.7
  )
  # colour each feature type by its `feature` value
  Gviz::displayPars(atrack) <- as.list(FEATURE_COLORS)

  strack  <- Gviz::SequenceTrack(dss, chromosome = chr, name = chr, cex = 0.8)
  axtrack <- Gviz::GenomeAxisTrack(labelPos = "below", cex = 0.7)

  strand_lbl <- if (identical(design_row$strand, "+")) "+ strand" else "− strand"
  ttl <- sprintf("%s   %s→%s   (gene %s, gRNA on %s, PAM %s bp from edit)",
                 design_row$mutation, design_row$original_codon, design_row$mutant_codon,
                 chr, strand_lbl, design_row$pam_distance)

  from <- max(1L, min(feats$start) - pad)
  to   <- min(glen, max(feats$end) + pad)
  Gviz::plotTracks(
    list(axtrack, atrack, strack), from = from, to = to, chromosome = chr,
    main = ttl, cex.main = 0.9,
    background.title = "#2C7FB8", col.title = "white", col.axis = "black",
    fontcolor = "black"
  )
  invisible(TRUE)
}

# ==============================================================================
# Component sequence readout (HTML): guaranteed letters, every component shown
# ==============================================================================

.span <- function(txt, color = NULL, bg = NULL, bold = FALSE, underline = NULL) {
  st <- c(
    if (!is.null(color)) paste0("color:", color),
    if (!is.null(bg)) paste0("background:", bg),
    if (bold) "font-weight:700",
    if (!is.null(underline)) paste0("border-bottom:2px solid ", underline)
  )
  sprintf("<span style='%s'>%s</span>", paste(st, collapse = ";"), txt)
}

#' Build an HTML sequence view for a selected design: the local gene sequence
#' (WT and edited rows) with the protospacer, PAM, HA, edited codon and silent
#' PAM mutation coloured, plus a table of each component's full sequence.
#' Returns an HTML string (wrap with shiny::HTML).
design_sequence_html <- function(design_row, system) {
  key <- tryCatch(resolve_gene_id(design_row$gene, system$id_maps, "auto"),
                  error = function(e) NULL)
  if (is.null(key)) key <- as.character(design_row$gene)
  gene <- toupper(system$gene_sequences[[key]])
  if (is.null(gene) || !nzchar(gene)) return("<em>Gene sequence unavailable.</em>")
  glen <- nchar(gene)

  cs  <- as.integer(design_row$codon_start)
  ps  <- min(as.integer(design_row$proto_start), as.integer(design_row$proto_end))
  pe  <- max(as.integer(design_row$proto_start), as.integer(design_row$proto_end))
  pamS <- as.integer(design_row$pam_position); pamE <- pamS + 2L
  smp <- as.integer(design_row$silent_mut_pos)
  haS <- as.integer(design_row$ha_start); haE <- as.integer(design_row$ha_end)
  ha  <- toupper(design_row$homology_arm)

  C_PROTO <- "#2C7FB8"; C_PAM <- "#D95F0E"; C_HA <- "#31A354"
  BG_EDIT <- "#fde0dd"; BG_SILENT <- "#e6ddf3"

  edited_base <- function(p) {
    if (p >= haS && p <= haE) substr(ha, p - haS + 1L, p - haS + 1L) else substr(gene, p, p)
  }

  from <- max(1L, min(ps, cs, pamS) - 4L)
  to   <- min(glen, max(pe, cs + 2L, pamE) + 4L)

  wt_cells <- character(0); ed_cells <- character(0)
  for (p in from:to) {
    wtb <- substr(gene, p, p); edb <- edited_base(p)
    col <- NULL; bold <- FALSE; und <- NULL; bg <- NULL
    if (p >= ps && p <= pe) col <- C_PROTO
    if (p >= pamS && p <= pamE) { col <- C_PAM; bold <- TRUE }
    if (p >= haS && p <= haE) und <- C_HA
    if (p >= cs && p <= cs + 2L) bg <- BG_EDIT
    if (p == smp) bg <- BG_SILENT
    wt_cells <- c(wt_cells, .span(wtb, color = col, bg = bg, bold = bold, underline = und))
    ed_cells <- c(ed_cells, .span(edb, color = col, bg = bg,
                                  bold = bold || (edb != wtb), underline = und))
  }

  legend <- paste(
    .span("protospacer", color = C_PROTO), "&nbsp;",
    .span("PAM", color = C_PAM, bold = TRUE), "&nbsp;",
    .span("HA/donor", underline = C_HA), "&nbsp;",
    .span("edited codon", bg = BG_EDIT), "&nbsp;",
    .span("silent PAM mut", bg = BG_SILENT)
  )

  mono <- function(cells) paste0(
    "<div style='font-family:monospace;font-size:15px;letter-spacing:2px;",
    "white-space:nowrap;overflow-x:auto'>", paste(cells, collapse = ""), "</div>")

  # HA with the edit codon + silent mutation highlighted (full donor)
  ha_cells <- character(nchar(ha))
  for (i in seq_len(nchar(ha))) {
    gp <- haS + i - 1L
    b <- substr(ha, i, i); bg <- NULL; bold <- FALSE
    if (gp >= cs && gp <= cs + 2L) { bg <- BG_EDIT; bold <- TRUE }
    if (gp == smp) { bg <- BG_SILENT; bold <- TRUE }
    ha_cells[i] <- .span(b, bg = bg, bold = bold)
  }

  tbl <- paste0(
    "<table class='table table-sm' style='font-family:monospace;font-size:13px'>",
    "<tr><td><b>Protospacer (20 nt)</b></td><td>", .span(design_row$protospacer, color = C_PROTO), "</td></tr>",
    "<tr><td><b>Codon change</b></td><td>", design_row$original_codon, " &rarr; ", design_row$mutant_codon, "</td></tr>",
    "<tr><td><b>Silent PAM mutation</b></td><td>", design_row$pam_mutation, "</td></tr>",
    "<tr><td><b>HA / repair donor</b></td><td style='white-space:normal;word-break:break-all'>",
    paste(ha_cells, collapse = ""), "</td></tr>",
    "</table>")

  paste0(
    "<div style='margin-bottom:6px'>", legend, "</div>",
    "<div class='text-muted small'>Gene region ", from, "–", to,
    " (reading frame); top = wild type, bottom = after edit:</div>",
    mono(wt_cells), mono(ed_cells),
    "<div class='mt-3'><b>Component sequences</b></div>", tbl
  )
}
