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
# ggplot "alignment map": actual letters, aligned, edits popping out
# ==============================================================================

gg_available <- function() requireNamespace("ggplot2", quietly = TRUE)

#' Lightweight sequence-alignment map (ggplot): the functional region as letters
#' in two aligned rows (wild type / after edit), coloured by component, with the
#' edited codon and silent PAM mutation popping out in red. No Gviz needed.
render_design_ggmap <- function(design_row, system) {
  if (!gg_available()) stop("ggplot2 is not installed.", call. = FALSE)
  key <- tryCatch(resolve_gene_id(design_row$gene, system$id_maps, "auto"),
                  error = function(e) NULL)
  if (is.null(key)) key <- as.character(design_row$gene)
  gene <- toupper(system$gene_sequences[[key]])
  if (is.null(gene) || !nzchar(gene)) stop("Gene sequence unavailable.", call. = FALSE)
  glen <- nchar(gene)

  cs  <- as.integer(design_row$codon_start)
  ps  <- min(as.integer(design_row$proto_start), as.integer(design_row$proto_end))
  pe  <- max(as.integer(design_row$proto_start), as.integer(design_row$proto_end))
  pamS <- as.integer(design_row$pam_position); pamE <- pamS + 2L
  smp <- as.integer(design_row$silent_mut_pos)
  haS <- as.integer(design_row$ha_start); haE <- as.integer(design_row$ha_end)
  ha  <- toupper(design_row$homology_arm)
  edited_base <- function(p) if (p >= haS && p <= haE) substr(ha, p - haS + 1L, p - haS + 1L) else substr(gene, p, p)

  # Window covers the whole HA (donor) plus the protospacer/PAM and some flank
  from <- max(1L, min(ps, cs, pamS, haS) - 3L)
  to   <- min(glen, max(pe, cs + 2L, pamE, haE) + 3L)
  pos  <- from:to

  feat_of <- function(p) {
    if (p >= cs && p <= cs + 2L) return("edit")
    if (p == smp) return("silentPAM")
    if (p >= pamS && p <= pamE) return("PAM")
    if (p >= ps && p <= pe) return("protospacer")
    if (p >= haS && p <= haE) return("HA")
    "none"
  }
  features <- vapply(pos, feat_of, character(1))
  wt <- vapply(pos, function(p) substr(gene, p, p), character(1))
  ed <- vapply(pos, edited_base, character(1))

  fcol <- c(protospacer = "#2C7FB8", PAM = "#D95F0E", HA = "#31A354",
            edit = "#E34A33", silentPAM = "#756BB1", none = "grey35")
  d <- rbind(
    data.frame(pos = pos, y = 2, base = wt, feature = features, changed = FALSE,
               stringsAsFactors = FALSE),
    data.frame(pos = pos, y = 1, base = ed, feature = features, changed = (wt != ed),
               stringsAsFactors = FALSE)
  )
  # SnapGene-style stagger: edited bases that differ from WT are lifted one block
  # up out of the row so they pop out.
  d$yplot <- ifelse(d$changed, d$y + 0.45, d$y)
  d$tcol  <- ifelse(d$changed, "#C0392B", fcol[d$feature])
  d$tsize <- ifelse(d$changed, 5.6, 3.9)
  d$tface <- ifelse(d$changed, "bold", "plain")

  # one label per feature type (at the midpoint of its longest run)
  rr <- rle(features); ends <- cumsum(rr$lengths); starts <- ends - rr$lengths + 1L
  lab_names <- c(protospacer = "protospacer", PAM = "PAM", HA = "HA (donor)",
                 edit = "edit", silentPAM = "silent PAM")
  runs_df <- data.frame(feature = rr$values, len = rr$lengths,
                        mid = pos[round((starts + ends) / 2)], stringsAsFactors = FALSE)
  runs_df <- runs_df[runs_df$feature != "none", , drop = FALSE]
  labs <- do.call(rbind, lapply(split(runs_df, runs_df$feature), function(x) {
    x <- x[which.max(x$len), ]
    data.frame(pos = x$mid, feature = x$feature, label = unname(lab_names[x$feature]),
               stringsAsFactors = FALSE)
  }))

  gg <- ggplot2::ggplot(d) +
    ggplot2::geom_tile(ggplot2::aes(pos, y, fill = feature), width = 1, height = 0.92,
                       alpha = 0.16, colour = NA) +
    ggplot2::geom_tile(data = d[d$changed, , drop = FALSE],
                       ggplot2::aes(pos, yplot), width = 0.92, height = 0.82,
                       fill = "white", colour = "#C0392B", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(pos, yplot, label = base, colour = tcol, size = tsize),
                       fontface = d$tface, family = "mono") +
    ggplot2::scale_fill_manual(values = fcol, guide = "none") +
    ggplot2::scale_colour_identity() + ggplot2::scale_size_identity() +
    ggplot2::scale_y_continuous(breaks = c(1, 2), labels = c("after edit", "wild type"),
                                limits = c(0.3, 3.3)) +
    ggplot2::labs(x = sprintf("%s position (reading frame)", key), y = NULL,
                  title = sprintf("%s   %s→%s   gRNA on %s strand, PAM %s bp from edit",
                                  design_row$mutation, design_row$original_codon,
                                  design_row$mutant_codon,
                                  if (identical(design_row$strand, "+")) "+" else "−",
                                  design_row$pam_distance)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 11, face = "bold"))

  if (!is.null(labs) && nrow(labs) > 0) {
    # lift the single-base silent-PAM label so it clears the adjacent PAM label
    labs$ylab <- ifelse(labs$feature == "silentPAM", 3.12, 2.8)
    gg <- gg + ggplot2::geom_text(data = labs,
             ggplot2::aes(pos, ylab, label = label, colour = fcol[feature]),
             size = 3.4, fontface = "bold")
  }
  # protospacer 5'->3' direction arrow
  ax <- if (identical(design_row$strand, "+")) c(ps, pe) else c(pe, ps)
  gg <- gg + ggplot2::annotate("segment", x = ax[1], xend = ax[2], y = 0.5, yend = 0.5,
             colour = "#2C7FB8", linewidth = 0.9,
             arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed"))
  gg
}

#' Alignment of the guide protospacer against one off-target (HTML, mismatches
#' in red with a x in the match line). `ot_row` is a row from design_off_targets().
offtarget_alignment_html <- function(guide, ot_row) {
  guide <- toupper(guide); ot <- toupper(ot_row$protospacer)
  n <- min(nchar(guide), nchar(ot))
  g <- strsplit(guide, "")[[1]]; o <- strsplit(ot, "")[[1]]
  # mismatched bases are lifted out of their row (guide up, off-target down) so
  # they pop out, SnapGene-style; matching bases stay in line.
  cell <- function(ch, mm, offset) {
    if (!mm) return(sprintf("<span style='color:#555'>%s</span>", ch))
    sprintf("<span style='color:#C0392B;font-weight:700;position:relative;top:%s'>%s</span>",
            offset, ch)
  }
  topc <- character(n); midc <- character(n); botc <- character(n)
  for (i in seq_len(n)) {
    mm <- g[i] != o[i]
    topc[i] <- cell(g[i], mm, "-9px")
    botc[i] <- cell(o[i], mm, "9px")
    midc[i] <- if (mm) "<span>&nbsp;</span>" else "<span style='color:#9fb6c7'>|</span>"
  }

  # Append the PAM (3' of the protospacer) as a bold, boxed codon so it stands
  # out: the off-target's actual PAM on the OT row, the required NGG motif on the
  # guide row. A small gap separates it from the 20-nt spacer, SnapGene-style.
  ot_pam <- toupper(ot_row$pam %||% "NGG")
  pam_box <- function(txt, col) sprintf(
    paste0("<span style='display:inline-block;margin-left:10px;padding:0 4px;",
           "border:1.5px solid %s;border-radius:4px;color:%s;font-weight:800;",
           "letter-spacing:1px'>%s</span>"), col, col, txt)
  topc <- c(topc, pam_box("NGG", "#8a8a8a"))
  botc <- c(botc, pam_box(ot_pam, "#D95F0E"))
  midc <- c(midc, sprintf(
    "<span style='display:inline-block;margin-left:10px;color:#aaa;font-weight:700'>PAM</span>"))
  # Fixed-width label column (flex) so both sequences start at the SAME x,
  # regardless of label text width.
  row <- function(lead, cells) paste0(
    "<div style='display:flex;align-items:center;font-family:monospace;font-size:16px;line-height:2.0'>",
    "<span style='flex:0 0 92px;color:#888;white-space:nowrap'>", lead, "</span>",
    "<span style='letter-spacing:3px;white-space:nowrap'>", paste(cells, collapse = ""), "</span></div>")

  gene_lbl <- if (is.null(ot_row$gene) || is.na(ot_row$gene)) "intergenic" else paste0("<i>", ot_row$gene, "</i>")
  seedtag <- if (isTRUE(ot_row$seed_match)) " &middot; <span style='color:#C0392B'>seed match</span>" else ""
  hdr <- sprintf(
    "<b>%s</b> at %s (%s strand) in %s &middot; %d mismatch%s &middot; PAM %s%s",
    ot_row$type, format(ot_row$position, big.mark = ","), ot_row$strand, gene_lbl,
    ot_row$mismatches, if (identical(ot_row$mismatches, 1L)) "" else "es",
    ot_row$pam %||% "NGG", seedtag)

  paste0("<div class='mb-2 small'>", hdr, "</div>",
         row("guide 5'-", topc), row("", midc), row("OT 5'-", botc))
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
