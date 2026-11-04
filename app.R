# ==============================================================================
# app.R  --  CRISPR Oligo Designer (Shiny GUI)
# ==============================================================================
# Thin GUI layer over the backend pipeline in R/. Run with:
#     Rscript run.R           (recommended)
#   or  shiny::runApp()       from this folder in an R session.
# ==============================================================================

# --- Load backend modules (order matters for top-level definitions) -----------
.backend_files <- c("packages.R", "constants.R", "config.R", "genome.R",
                    "offtarget.R", "design.R", "microsynth.R", "qc.R", "pipeline.R")
for (.f in .backend_files) source(file.path("R", .f))

load_gui_packages()

# --- Base configuration (from config/config.yml) ------------------------------
CFG0 <- load_config("config/config.yml")
apply_design_config(CFG0)

EXAMPLE_MUTATIONS <- "dpsS100E\ndpsS100D\nlrpS125D\nlrpS125E\nthrBF45A\nmetAL60W"

# ==============================================================================
# UI
# ==============================================================================

ui <- bslib::page_sidebar(
  title = "CRISPR Oligo Designer — E. coli BW25113",
  theme = bslib::bs_theme(version = 5, primary = "#2c7fb8"),
  shinyjs::useShinyjs(),

  sidebar = bslib::sidebar(
    width = 360,
    bslib::accordion(
      open = c("Step 1 · Genome", "Step 2 · Mutations", "Step 3 · Run"),

      # ---- Genome ----
      bslib::accordion_panel(
        "Step 1 · Genome", icon = shiny::icon("dna"),
        shiny::radioButtons("genome_source", "Genome source",
                            choices = c("Standard E. coli strain" = "strain",
                                        "Upload my own genome" = "custom"),
                            selected = "strain"),
        shiny::conditionalPanel(
          condition = "input.genome_source == 'strain'",
          shiny::selectInput("strain", "Strain",
                             choices = names(CFG0$genome$strains),
                             selected = CFG0$genome$default_strain),
          shiny::p(shiny::tags$small(
            "Downloaded once from NCBI and cached. First load may take a minute."))
        ),
        shiny::conditionalPanel(
          condition = "input.genome_source == 'custom'",
          shiny::fileInput("genome_fasta", "Genome FASTA (.fasta / .fa / .fna)",
                           accept = c(".fasta", ".fa", ".fna", ".txt")),
          shiny::fileInput("genome_gff", "Annotation GFF/GTF (recommended)",
                           accept = c(".gff", ".gff3", ".gtf")),
          shiny::p(shiny::tags$small(
            "The first FASTA record is used as the chromosome. An annotation file ",
            "is needed to map mutations to genes; without it, gene lookup is skipped."))
        ),
        shiny::checkboxInput("force_download", "Force re-download / rebuild", value = FALSE),
        shiny::actionButton("init_genome", "Load / initialise genome",
                            class = "btn-primary w-100"),
        shiny::br(), shiny::br(),
        shiny::uiOutput("genome_status")
      ),

      # ---- Mutations ----
      bslib::accordion_panel(
        "Step 2 · Mutations", icon = shiny::icon("list"),
        shiny::uiOutput("id_type_ui"),
        shiny::fileInput("mut_file", "Upload mutations (.csv / .txt / .xlsx)",
                         accept = c(".csv", ".txt", ".tsv", ".xlsx")),
        shiny::tags$small("…or paste them below (one per line, e.g. dpsS100E):"),
        shiny::textAreaInput("mut_text", NULL, rows = 6,
                             placeholder = "dpsS100E\ndpsS100D\nlrpS125D"),
        shiny::actionButton("load_example", "Load example",
                            class = "btn-outline-secondary btn-sm"),
        shiny::actionButton("validate_btn", "Validate",
                            class = "btn-secondary btn-sm")
      ),

      # ---- Run ----
      bslib::accordion_panel(
        "Step 3 · Run", icon = shiny::icon("play"),
        shiny::actionButton("run_btn", "Run design pipeline",
                            class = "btn-success w-100"),
        shiny::br(), shiny::br(),
        shiny::uiOutput("run_status")
      ),

      # ---- Settings ----
      bslib::accordion_panel(
        "Settings", icon = shiny::icon("gear"),
        shiny::numericInput("oligo_size", "Microsynth oligo size (bp)",
                            value = CFG0$design$oligo_size, min = 20, max = 300, step = 1),
        shiny::numericInput("max_pam_distance", "Max PAM distance (bp)",
                            value = CFG0$design$max_pam_distance, min = 3, max = 60, step = 1),
        shiny::numericInput("min_score_threshold", "Min score threshold",
                            value = CFG0$design$min_score_threshold, min = 0, max = 5, step = 0.1),
        shiny::numericInput("top_n", "Designs per mutation",
                            value = CFG0$design$top_n, min = 1, max = 10, step = 1),
        shiny::numericInput("n_cores", "CPU cores",
                            value = CFG0$design$n_cores, min = 1, max = 16, step = 1),
        shiny::hr(),
        shiny::strong("Order form (top-ranked)"),
        shiny::numericInput("export_min_score", "Min score for order form",
                            value = CFG0$export$min_score_filter, min = 0, max = 5, step = 0.1),
        shiny::textInput("exclude_mut_aa", "Exclude mutant amino acids (comma-sep)",
                         value = ""),
        shiny::numericInput("pto_n_bases", "PTO-protected 5' bases (homology arm)",
                            value = CFG0$microsynth$pto_n_bases, min = 0, max = 10, step = 1),
        shiny::hr(),
        shiny::strong("gRNA cloning oligos"),
        shiny::radioButtons("cloning_method", "Cloning method",
                            choices = c("Gibson (whole-plasmid PCR)" = "gibson",
                                        "Golden Gate (BsaI, annealed pair)" = "golden_gate",
                                        "Both" = "both"),
                            selected = CFG0$microsynth$cloning_method),
        shiny::conditionalPanel(
          condition = "input.cloning_method == 'gibson' || input.cloning_method == 'both'",
          shiny::textInput("fwd_tail", "Gibson forward 3' tail",
                           value = CFG0$microsynth$fwd_tail),
          shiny::textInput("rev_tail", "Gibson reverse 3' tail",
                           value = CFG0$microsynth$rev_tail)
        ),
        shiny::conditionalPanel(
          condition = "input.cloning_method == 'golden_gate' || input.cloning_method == 'both'",
          shiny::textInput("gg_top_overhang", "Golden Gate top overhang (4 nt)",
                           value = CFG0$microsynth$gg_top_overhang),
          shiny::textInput("gg_bot_overhang", "Golden Gate bottom overhang (4 nt)",
                           value = CFG0$microsynth$gg_bot_overhang),
          shiny::checkboxInput("gg_phosphorylate", "Order Golden Gate oligos 5'-phosphorylated",
                               value = isTRUE(CFG0$microsynth$gg_phosphorylate))
        ),
        shiny::hr(),
        shiny::strong("QC (Tm / GC / hairpins)"),
        shiny::numericInput("qc_na_mM", "Monovalent cation for Tm (mM)",
                            value = CFG0$qc$na_mM, min = 1, max = 1000, step = 1),
        shiny::fileInput("template_file", "Microsynth template (.xlsx) — optional override",
                         accept = ".xlsx")
      )
    )
  ),

  # ---- Main panel ----
  bslib::navset_card_tab(
    id = "main_tabs",

    bslib::nav_panel(
      "Validation", icon = shiny::icon("clipboard-check"),
      shiny::uiOutput("validation_summary"),
      DT::DTOutput("validation_table")
    ),

    bslib::nav_panel(
      "Results", icon = shiny::icon("table"),
      shiny::div(
        class = "d-flex gap-2 align-items-center mb-2 flex-wrap",
        shiny::radioButtons("results_view", NULL,
                            choices = c("Top-ranked (order candidates)" = "top",
                                        "All designs" = "all"),
                            selected = "top", inline = TRUE),
        shiny::downloadButton("dl_top", "Download top-ranked CSV",
                              class = "btn-outline-primary btn-sm"),
        shiny::downloadButton("dl_all", "Download all designs CSV",
                              class = "btn-outline-secondary btn-sm")
      ),
      DT::DTOutput("results_table")
    ),

    bslib::nav_panel(
      "Order form", icon = shiny::icon("file-excel"),
      shiny::uiOutput("order_summary"),
      shiny::div(
        class = "d-flex gap-2 mb-2",
        shiny::downloadButton("dl_microsynth", "Generate & download Microsynth form",
                              class = "btn-success"),
        shiny::downloadButton("dl_oligos", "Download oligo list CSV",
                              class = "btn-outline-secondary")
      ),
      DT::DTOutput("oligo_table")
    ),

    bslib::nav_panel(
      "QC", icon = shiny::icon("flask"),
      shiny::uiOutput("qc_summary"),
      shiny::div(class = "mb-2",
        shiny::downloadButton("dl_qc", "Download QC CSV",
                              class = "btn-outline-secondary btn-sm")),
      shiny::uiOutput("qc_backend_note"),
      DT::DTOutput("qc_table_out")
    ),

    bslib::nav_panel(
      "Summary", icon = shiny::icon("file-lines"),
      shiny::div(class = "mb-2",
        shiny::downloadButton("dl_summary", "Download summary (.txt)",
                              class = "btn-outline-secondary btn-sm")),
      shiny::verbatimTextOutput("summary_text")
    ),

    bslib::nav_panel(
      "Help", icon = shiny::icon("circle-question"),
      shiny::htmlOutput("help_html")
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  rv <- shiny::reactiveValues(
    system = NULL, mutations = character(0), validation = NULL,
    results_all = NULL, results_top = NULL, summary = NULL,
    n_input = 0, elapsed = NULL
  )

  # ---- Effective config from base file + GUI overrides ----------------------
  effective_cfg <- shiny::reactive({
    cfg <- CFG0
    cfg$design$oligo_size          <- input$oligo_size
    cfg$design$max_pam_distance    <- input$max_pam_distance
    cfg$design$min_score_threshold <- input$min_score_threshold
    cfg$design$top_n               <- input$top_n
    cfg$design$n_cores             <- input$n_cores
    cfg$export$min_score_filter    <- input$export_min_score
    cfg$microsynth$pto_n_bases     <- input$pto_n_bases
    cfg$microsynth$cloning_method  <- input$cloning_method %||% CFG0$microsynth$cloning_method
    # Tails kept verbatim (lowercase convention marks the vector-homology region)
    cfg$microsynth$fwd_tail        <- trimws(input$fwd_tail %||% CFG0$microsynth$fwd_tail)
    cfg$microsynth$rev_tail        <- trimws(input$rev_tail %||% CFG0$microsynth$rev_tail)
    cfg$microsynth$gg_top_overhang <- toupper(trimws(input$gg_top_overhang %||% CFG0$microsynth$gg_top_overhang))
    cfg$microsynth$gg_bot_overhang <- toupper(trimws(input$gg_bot_overhang %||% CFG0$microsynth$gg_bot_overhang))
    cfg$microsynth$gg_phosphorylate <- isTRUE(input$gg_phosphorylate)
    cfg$qc$na_mM                   <- input$qc_na_mM %||% CFG0$qc$na_mM
    if (!is.null(input$template_file)) {
      cfg$microsynth$template_path <- input$template_file$datapath
    }
    cfg
  })

  exclude_aa <- shiny::reactive({
    x <- trimws(unlist(strsplit(toupper(input$exclude_mut_aa), "[,;\\s]+", perl = TRUE)))
    x[nzchar(x)]
  })

  # ==== Step 1: Genome ========================================================

  # Build the genome source spec from the current UI selections.
  genome_source_spec <- function() {
    if (input$genome_source == "custom") {
      if (is.null(input$genome_fasta)) {
        return(list(error = "Please choose a genome FASTA file to upload."))
      }
      gff_path <- NULL; gff_fmt <- NULL
      if (!is.null(input$genome_gff)) {
        gff_path <- input$genome_gff$datapath
        ext <- tolower(tools::file_ext(input$genome_gff$name))
        gff_fmt <- if (ext %in% c("gff", "gff3")) "gff3" else if (ext == "gtf") "gtf" else NULL
      }
      custom_source(genome_file = input$genome_fasta$datapath,
                    gff_file = gff_path, gff_format = gff_fmt,
                    name = paste0("Custom: ", input$genome_fasta$name))
    } else {
      acc <- CFG0$genome$strains[[input$strain]]
      strain_source(accession = acc, name = input$strain)
    }
  }

  shiny::observeEvent(input$init_genome, {
    src <- genome_source_spec()
    if (!is.null(src$error)) {
      shiny::showNotification(src$error, type = "warning"); return()
    }

    prog <- shiny::Progress$new(session)
    prog$set(message = "Initialising genome", value = 0)
    on.exit(prog$close())
    cb <- function(frac, msg) prog$set(value = frac, message = msg)

    cfg <- effective_cfg(); apply_design_config(cfg)
    res <- tryCatch(
      initialize_system_fast(src, force_download = isTRUE(input$force_download),
                             progress_callback = cb),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      shiny::showNotification(paste("Genome init failed:", conditionMessage(res)),
                              type = "error", duration = NULL)
      rv$system <- NULL
    } else {
      rv$system <- res
      if (length(res$gene_sequences) == 0) {
        msg <- if (input$genome_source == "strain") {
          paste("Genome sequence loaded, but its annotation could not be fetched",
                "from NCBI (this is sometimes a transient server error). Tick",
                "'Force re-download' and try again, pick another strain, or upload",
                "the genome + a GFF under 'Upload my own genome'.")
        } else {
          paste("Genome loaded, but no gene annotations were found. Gene-based",
                "mutation design will not work — upload a matching GFF/GTF file.")
        }
        shiny::showNotification(msg, type = "warning", duration = NULL)
      } else {
        shiny::showNotification(paste0("Genome ready: ", res$genome_name, "."),
                                type = "message")
      }
      # re-validate now that genes are known
      if (length(rv$mutations) > 0)
        rv$validation <- validate_mutations(rv$mutations, rv$system, current_id_type())
    }
  })

  output$genome_status <- shiny::renderUI({
    if (is.null(rv$system)) {
      shiny::div(class = "text-muted", shiny::icon("circle"), " Not loaded")
    } else {
      n_genes <- length(rv$system$gene_sequences)
      cls <- if (n_genes > 0) "text-success" else "text-warning"
      icn <- if (n_genes > 0) "circle-check" else "triangle-exclamation"
      id_line <- NULL
      if (length(rv$system$id_types) > 0) {
        parts <- sprintf("%s (%s)",
                         vapply(names(rv$system$id_types),
                                function(t) ID_TYPE_LABELS[[t]] %||% t, character(1)),
                         format(rv$system$id_types, big.mark = ","))
        id_line <- shiny::div(class = "text-muted small mt-1",
                              "Identifiers available: ", paste(parts, collapse = " · "))
      }
      shiny::tagList(
        shiny::div(class = cls,
          shiny::icon(icn),
          sprintf(" %s — %s genes, %s bp",
                  rv$system$genome_name %||% "Genome",
                  format(n_genes, big.mark = ","),
                  format(nchar(rv$system$genome_seq), big.mark = ","))),
        id_line
      )
    }
  })

  # Identifier-type selector, populated from the loaded genome.
  id_type_choices <- shiny::reactive({
    types <- if (!is.null(rv$system) && length(rv$system$id_types) > 0) {
      names(rv$system$id_types)
    } else {
      "gene"
    }
    labels <- vapply(types, function(t) {
      lab <- ID_TYPE_LABELS[[t]] %||% t
      if (!is.null(rv$system) && t %in% names(rv$system$id_types)) {
        sprintf("%s (%s)", lab, format(rv$system$id_types[[t]], big.mark = ","))
      } else lab
    }, character(1))
    c(setNames(types, labels), c("Auto — search all types" = "auto"))
  })

  output$id_type_ui <- shiny::renderUI({
    ch <- id_type_choices()
    sel <- if ("gene" %in% ch) "gene" else unname(ch)[1]
    shiny::selectInput("id_type", "Mutation identifier type", choices = ch, selected = sel)
  })

  current_id_type <- shiny::reactive(input$id_type %||% "gene")

  # ==== Step 2: Mutations & validation =======================================
  shiny::observeEvent(input$load_example, {
    shiny::updateTextAreaInput(session, "mut_text", value = EXAMPLE_MUTATIONS)
  })

  gather_mutations <- function() {
    muts <- character(0)
    if (!is.null(input$mut_file)) {
      muts <- tryCatch(read_mutations_file(input$mut_file$datapath),
                       error = function(e) {
                         shiny::showNotification(paste("Could not read file:",
                                                       conditionMessage(e)), type = "error")
                         character(0)
                       })
    }
    if (nzchar(trimws(input$mut_text %||% ""))) {
      pasted <- unlist(strsplit(input$mut_text, "[,;\\s]+", perl = TRUE))
      muts <- c(muts, trimws(pasted))
    }
    muts <- muts[nzchar(muts)]
    muts
  }

  shiny::observeEvent(input$validate_btn, {
    muts <- gather_mutations()
    if (length(muts) == 0) {
      shiny::showNotification("No mutations provided.", type = "warning"); return()
    }
    rv$mutations <- muts
    rv$validation <- validate_mutations(muts, rv$system, current_id_type())
    bslib::nav_select("main_tabs", "Validation")
  })

  output$validation_summary <- shiny::renderUI({
    v <- rv$validation
    if (is.null(v)) return(shiny::p("Upload or paste mutations, then click ",
                                    shiny::strong("Validate"), "."))
    n_run <- length(runnable_mutations(v))
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Submitted", nrow(v), theme = "primary"),
      bslib::value_box("Valid format", sum(v$valid_format), theme = "secondary"),
      bslib::value_box("Ready to run", n_run, theme = "success"),
      bslib::value_box("Problems",
                       sum(!v$valid_format | v$status %in% c("Gene not found", "Duplicate")),
                       theme = "warning")
    )
  })

  output$validation_table <- DT::renderDT({
    shiny::req(rv$validation)
    DT::datatable(rv$validation, filter = "top", rownames = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  # ==== Step 3: Run ==========================================================
  shiny::observeEvent(input$run_btn, {
    if (is.null(rv$system)) {
      shiny::showNotification("Load the genome first (Step 1).", type = "warning"); return()
    }
    muts <- gather_mutations()
    if (length(muts) == 0) {
      shiny::showNotification("No mutations provided.", type = "warning"); return()
    }
    rv$mutations <- muts
    cfg <- effective_cfg(); apply_design_config(cfg)

    rv$validation <- validate_mutations(muts, rv$system, current_id_type())
    to_run <- runnable_mutations(rv$validation)
    rv$n_input <- length(unique(muts))
    if (length(to_run) == 0) {
      shiny::showNotification("No runnable mutations after validation.", type = "error"); return()
    }

    prog <- shiny::Progress$new(session)
    prog$set(message = "Designing oligos", value = 0)
    on.exit(prog$close())
    cb <- function(frac, msg) prog$set(value = frac, message = msg)

    t0 <- Sys.time()
    res <- tryCatch(run_design(to_run, rv$system, cfg, progress_callback = cb,
                               id_type = current_id_type()),
                    error = function(e) e)
    if (inherits(res, "error")) {
      shiny::showNotification(paste("Design failed:", conditionMessage(res)),
                              type = "error", duration = NULL); return()
    }
    rv$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rv$results_all <- res
    rv$results_top <- select_top_ranked(res, min_score = input$export_min_score,
                                        exclude_mut_aa = exclude_aa())
    rv$summary <- build_run_summary(rv$validation, rv$results_all, rv$results_top,
                                    cfg, rv$n_input, rv$elapsed,
                                    genome_label = rv$system$genome_name)
    shiny::showNotification(
      sprintf("Done: %d designs for %d mutations.",
              nrow(res), length(unique(res$mutation))), type = "message")
    bslib::nav_select("main_tabs", "Results")
  })

  # Recompute top-ranked when filters change (without re-running design)
  shiny::observe({
    shiny::req(rv$results_all)
    rv$results_top <- select_top_ranked(rv$results_all, min_score = input$export_min_score,
                                        exclude_mut_aa = exclude_aa())
  })

  output$run_status <- shiny::renderUI({
    if (is.null(rv$results_all)) {
      shiny::div(class = "text-muted", "No run yet.")
    } else {
      shiny::div(class = "text-success",
        sprintf("%d designs · %d top-ranked",
                nrow(rv$results_all), if (is.null(rv$results_top)) 0 else nrow(rv$results_top)))
    }
  })

  # ==== Results ==============================================================
  output$results_table <- DT::renderDT({
    shiny::req(rv$results_all)
    df <- if (input$results_view == "top") rv$results_top else rv$results_all
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No rows for the current filter."))
    DT::datatable(df, filter = "top", rownames = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_all <- shiny::downloadHandler(
    filename = function() paste0("generated_oligos_all_", Sys.Date(), ".csv"),
    content = function(file) utils::write.csv(rv$results_all, file, row.names = FALSE)
  )
  output$dl_top <- shiny::downloadHandler(
    filename = function() paste0("generated_oligos_top_", Sys.Date(), ".csv"),
    content = function(file) utils::write.csv(rv$results_top, file, row.names = FALSE)
  )

  # ==== Order form ===========================================================
  oligo_orders <- shiny::reactive({
    shiny::req(rv$results_top)
    build_oligo_orders(rv$results_top, effective_cfg())
  })

  # QC of the current oligo set
  oligo_qc <- shiny::reactive({
    shiny::req(rv$results_top)
    cfg <- effective_cfg()
    qc_table(oligo_orders()$oligos, cfg$qc,
             ha_suffix = cfg$microsynth$ha_name_suffix %||% "_HA")
  })

  output$order_summary <- shiny::renderUI({
    if (is.null(rv$results_top) || nrow(rv$results_top) == 0)
      return(shiny::p("Run the pipeline first. The order form is built from the ",
                      shiny::strong("top-ranked"), " designs."))
    o <- oligo_orders()
    cfg <- effective_cfg()
    tmpl_ok <- file.exists(cfg$microsynth$template_path)
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Mutations", nrow(rv$results_top), theme = "primary"),
      bslib::value_box("Unique oligos", o$n_after, theme = "success"),
      bslib::value_box("Duplicates removed", o$n_duplicates, theme = "secondary"),
      bslib::value_box("Template",
                       if (tmpl_ok) "Found" else "Missing",
                       theme = if (tmpl_ok) "success" else "danger")
    )
  })

  output$oligo_table <- DT::renderDT({
    o <- oligo_orders()
    DT::datatable(o$oligos, filter = "top", rownames = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_oligos <- shiny::downloadHandler(
    filename = function() paste0("oligo_list_", Sys.Date(), ".csv"),
    content = function(file) utils::write.csv(oligo_orders()$oligos, file, row.names = FALSE)
  )

  output$dl_microsynth <- shiny::downloadHandler(
    filename = function() paste0("MicrosynthOrder_filled_", Sys.Date(), ".xlsx"),
    content = function(file) {
      cfg <- effective_cfg()
      tryCatch(
        generate_microsynth_output(rv$results_top, cfg, output_path = file),
        error = function(e) {
          shiny::showNotification(conditionMessage(e), type = "error", duration = NULL)
          # write a placeholder so the download does not silently fail
          stop(e)
        }
      )
    }
  )

  # ==== QC ===================================================================
  output$qc_summary <- shiny::renderUI({
    if (is.null(rv$results_top) || nrow(rv$results_top) == 0)
      return(shiny::p("Run the pipeline first — QC is computed on the ordered oligos."))
    q <- oligo_qc()
    n_flag <- sum(nzchar(q$flags))
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Oligos", nrow(q), theme = "primary"),
      bslib::value_box("Flagged", n_flag, theme = if (n_flag > 0) "warning" else "success"),
      bslib::value_box("Hairpin ≥ thr",
                       sum(q$hairpin_stem >= (effective_cfg()$qc$hairpin_stem %||% 6)),
                       theme = "secondary"),
      bslib::value_box("BsaI sites",
                       if ("BsaI_site" %in% names(q)) sum(q$BsaI_site) else 0,
                       theme = if (("BsaI_site" %in% names(q)) && any(q$BsaI_site)) "danger" else "success")
    )
  })

  output$qc_backend_note <- shiny::renderUI({
    q <- oligo_qc()
    p3 <- "tm_method" %in% names(q) && any(q$tm_method == "primer3")
    if (p3) {
      shiny::p(shiny::tags$small(
        shiny::icon("circle-check"), " Thermodynamics: ", shiny::strong("primer3"),
        " — Tm, hairpin & self-dimer are true ΔG (kcal/mol)."))
    } else {
      shiny::p(shiny::tags$small(
        shiny::icon("circle"), " Thermodynamics: ", shiny::strong("built-in heuristic"),
        " (Tm = nearest-neighbor; hairpin/3'-dimer are stem lengths in bp, not ΔG). ",
        "Run ", shiny::code("pip install primer3-py"), " for accurate ΔG."))
    }
  })

  output$qc_table_out <- DT::renderDT({
    q <- oligo_qc()
    shiny::validate(shiny::need(nrow(q) > 0, "No oligos to QC."))
    dt <- DT::datatable(q, filter = "top", rownames = FALSE,
                        options = list(pageLength = 15, scrollX = TRUE))
    # highlight flagged rows
    DT::formatStyle(dt, "flags", target = "row",
                    backgroundColor = DT::styleEqual("", "white", default = "#fff3cd"))
  })

  output$dl_qc <- shiny::downloadHandler(
    filename = function() paste0("oligo_qc_", Sys.Date(), ".csv"),
    content = function(file) utils::write.csv(oligo_qc(), file, row.names = FALSE)
  )

  # ==== Summary ==============================================================
  output$summary_text <- shiny::renderText({
    if (is.null(rv$summary)) return("No run yet.")
    paste(rv$summary, collapse = "\n")
  })
  output$dl_summary <- shiny::downloadHandler(
    filename = function() paste0("run_summary_", Sys.Date(), ".txt"),
    content = function(file) writeLines(rv$summary %||% "No run.", file)
  )

  # ==== Help =================================================================
  output$help_html <- shiny::renderUI({
    shiny::HTML(
      "<h4>How to use</h4>
       <ol>
         <li><b>Choose a genome</b> (Step 1): pick a standard <i>E. coli</i> strain
             from the dropdown, <i>or</i> select <b>Upload my own genome</b> and
             provide a FASTA (+ optional GFF/GTF annotation). Then click
             <b>Load / initialise genome</b>. Downloaded strains are cached.</li>
         <li><b>Provide mutations</b> (Step 2): upload a .csv/.txt/.xlsx with a
             <code>mutation</code> column, or paste one mutation per line.</li>
         <li><b>Validate</b> to check the format and that each gene exists.</li>
         <li><b>Run design pipeline</b> (Step 3).</li>
         <li>Review <b>Results</b>, then download CSVs or the <b>Microsynth order form</b>.</li>
       </ol>
       <h4>Genome</h4>
       <p>The strain list (label &rarr; NCBI accession) lives in
          <code>config/config.yml</code> under <code>genome.strains</code> — edit it
          to add or correct strains. A custom genome needs an annotation file to
          map mutations to genes; the first FASTA record is used as the chromosome.</p>
       <h4>Mutation format</h4>
       <p>Compact (gene symbol):
          <code>&lt;gene&gt;&lt;WT-aa&gt;&lt;position&gt;&lt;MUT-aa&gt;</code>, e.g.
          <code>dpsS100E</code> = serine 100 &rarr; glutamate in gene <i>dps</i>.</p>
       <p>For identifiers that contain digits (locus tags, UniProt IDs, protein
          IDs), separate the identifier from the change with a colon:
          <code>&lt;id&gt;:&lt;WT&gt;&lt;pos&gt;&lt;MUT&gt;</code>, e.g.
          <code>b0014:S100E</code> or <code>P0A6Y8:S100E</code>.</p>
       <h4>Identifier type</h4>
       <p>Choose <b>Mutation identifier type</b> (Step 2) to say which kind of
          identifier your list uses. Only the types actually present in the loaded
          genome are offered — the counts next to each (and under the genome status)
          show how many genes carry that identifier. <b>Auto</b> tries every type.</p>
       <h4>Microsynth template</h4>
       <p>Place the official <code>MicrosynthUploadFormDNA__2_.xlsx</code> in the
          <code>templates/</code> folder (path set in <code>config/config.yml</code>),
          or upload it under <b>Settings</b>. The order form is built from the
          top-ranked design of each mutation, with 3 oligos each
          (<code>_pEcgRNA_fw_gibson</code>, <code>_pEcgRNA_rv_gibson</code>,
          <code>_HA</code>), deduplicated by sequence. The two gRNA primers are
          Gibson-assembly primers: protospacer (or its reverse complement) plus a
          24-nt 3' homology tail matching the pEcgRNA vector.</p>
       <h4>Oligo size</h4>
       <p>The <b>Microsynth oligo size</b> (Settings) sets the homology-arm repair
          oligo length; default 100 bp. It is also configurable in
          <code>config/config.yml</code>.</p>
       <h4>Cloning method</h4>
       <p>Under <b>Settings &rarr; gRNA cloning oligos</b> choose how the spacer is
          delivered: <b>Gibson</b> (whole-plasmid PCR primers,
          <code>_fw_gibson</code>/<code>_rv_gibson</code>), <b>Golden Gate</b> (an
          annealed oligo pair <code>_gg_bsaI_top</code>/<code>_bot</code> with 4-nt
          BsaI fusion overhangs, default <code>TAGT</code>/<code>AAAC</code>, optionally
          5'-phosphorylated), or <b>Both</b>. The <code>_HA</code> repair donor is
          always included.</p>
       <h4>QC tab</h4>
       <p>Every ordered oligo is checked for GC%, nearest-neighbor Tm, longest
          hairpin stem, 3'-end self-dimer, homopolymer runs, and internal
          <b>BsaI/BbsI</b> sites (a BsaI site in the spacer breaks Golden Gate
          assembly). Flagged oligos are highlighted; thresholds are in
          <code>config/config.yml &rarr; qc</code>. With <code>primer3-py</code> installed
          (<code>pip install primer3-py</code>) the hairpin/dimer values are true
          thermodynamic &Delta;G; otherwise they are heuristic stem lengths. The QC tab
          shows which backend is active.</p>")
  })
}

shiny::shinyApp(ui, server)
