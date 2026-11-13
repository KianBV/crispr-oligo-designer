# ==============================================================================
# app.R  --  CRISPR Oligo Designer (Shiny GUI)
# ==============================================================================
# Thin GUI layer over the backend pipeline in R/. Run with:
#     Rscript run.R           (recommended)
#   or  shiny::runApp()       from this folder in an R session.
# ==============================================================================

# --- Load backend modules (order matters for top-level definitions) -----------
.backend_files <- c("packages.R", "constants.R", "config.R", "genome.R",
                    "offtarget.R", "design.R", "microsynth.R", "qc.R",
                    "genome_view.R", "pipeline.R", "tracker_db.R")
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
        shiny::strong("Scoring (A/B toggles)"),
        shiny::p(shiny::tags$small("Off by default = original model. See SCORING.md.")),
        shiny::checkboxInput("score_swap23",
                             "Rank exact seed-match (class 3) above class 2",
                             value = isTRUE(CFG0$design$scoring$swap_class_2_3)),
        shiny::numericInput("score_ot_weight", "Off-target count penalty weight (0 = off)",
                            value = CFG0$design$scoring$offtarget_count_weight,
                            min = 0, max = 3, step = 0.1),
        shiny::checkboxInput("score_activity",
                             "On-target activity penalty (GC band + homopolymer)",
                             value = isTRUE(CFG0$design$scoring$activity_term)),
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
      shiny::p(shiny::tags$small(shiny::icon("hand-pointer"),
               " Click a row to inspect it (genomic map, sequence, off-targets).")),
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
      "Inspector", icon = shiny::icon("magnifying-glass-chart"),
      shiny::uiOutput("inspector_header"),
      shiny::conditionalPanel(
        condition = "output.inspector_ready == true",
        shiny::div(class = "mb-2 d-flex gap-2 align-items-center",
          shiny::actionButton("inspector_log_btn", "Log this chosen design to tracker",
                              class = "btn-success btn-sm",
                              icon = shiny::icon("floppy-disk")),
          shiny::tags$small(class = "text-muted",
            "Records this specific design (marked “chosen”), not the whole run.")),
        bslib::card(
          bslib::card_header("Alignment map (wild type vs edit)"),
          shiny::plotOutput("inspector_map", height = "300px")
        ),
        bslib::card(
          bslib::card_header("Sequence & components"),
          shiny::uiOutput("inspector_seq")
        ),
        bslib::card(
          bslib::card_header(shiny::uiOutput("inspector_ot_header")),
          shiny::uiOutput("inspector_ot_align"),
          DT::DTOutput("inspector_offtargets")
        )
      )
    ),

    bslib::nav_panel(
      "Summary", icon = shiny::icon("file-lines"),
      shiny::div(class = "mb-2",
        shiny::downloadButton("dl_summary", "Download summary (.txt)",
                              class = "btn-outline-secondary btn-sm")),
      shiny::verbatimTextOutput("summary_text")
    ),

    bslib::nav_panel(
      "Tracker", icon = shiny::icon("database"),
      shiny::div(
        class = "d-flex gap-2 align-items-end mb-3 flex-wrap",
        shiny::textInput("tracker_path", "Database file",
                         value = tracker_default_path(CFG0), width = "360px"),
        shiny::actionButton("tracker_connect", "Connect / create",
                            class = "btn-primary", icon = shiny::icon("plug")),
        shiny::actionButton("tracker_log_run", "Log all top picks (best per mutation)",
                            class = "btn-outline-success", icon = shiny::icon("layer-group"))
      ),
      shiny::uiOutput("tracker_status"),
      shiny::div(
        class = "d-flex gap-2 align-items-center mb-3 flex-wrap",
        shiny::tags$small(class = "text-muted", "Export CSV:"),
        shiny::downloadButton("tracker_dl_flat", "Edits + outcomes (flat)",
                              class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton("tracker_dl_edits", "Edits", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton("tracker_dl_attempts", "Attempts", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton("tracker_dl_curing", "Curing", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton("tracker_dl_strains", "Strains", class = "btn-outline-secondary btn-sm")
      ),
      bslib::navset_pill(
        bslib::nav_panel(
          "Edits",
          shiny::p(shiny::tags$small(shiny::icon("hand-pointer"),
                   " Click an edit to see its oligos and to record an attempt.")),
          DT::DTOutput("tracker_edits_tbl"),
          shiny::uiOutput("tracker_edit_detail")
        ),
        bslib::nav_panel(
          "Attempts",
          shiny::div(class = "mb-2",
            shiny::actionButton("tracker_add_attempt_btn", "Record attempt",
                                class = "btn-outline-primary btn-sm",
                                icon = shiny::icon("plus"))),
          DT::DTOutput("tracker_attempts_tbl")
        ),
        bslib::nav_panel(
          "Curing",
          shiny::div(class = "mb-2",
            shiny::actionButton("tracker_add_curing_btn", "Record curing outcome",
                                class = "btn-outline-primary btn-sm",
                                icon = shiny::icon("plus"))),
          DT::DTOutput("tracker_curing_tbl"),
          shiny::h6("Curing success by strain", class = "mt-3"),
          DT::DTOutput("tracker_curing_summary")
        ),
        bslib::nav_panel(
          "Strains",
          shiny::div(class = "mb-2",
            shiny::actionButton("tracker_add_strain_btn", "Add / edit strain",
                                class = "btn-outline-primary btn-sm",
                                icon = shiny::icon("plus"))),
          DT::DTOutput("tracker_strains_tbl")
        )
      )
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
    n_input = 0, elapsed = NULL,
    tracker_path = NULL, tracker_ready = FALSE, tracker_tick = 0
  )

  # ---- Effective config from base file + GUI overrides ----------------------
  effective_cfg <- shiny::reactive({
    cfg <- CFG0
    cfg$design$oligo_size          <- input$oligo_size
    cfg$design$max_pam_distance    <- input$max_pam_distance
    cfg$design$min_score_threshold <- input$min_score_threshold
    cfg$design$top_n               <- input$top_n
    cfg$design$n_cores             <- input$n_cores
    cfg$design$scoring <- modifyList(CFG0$design$scoring %||% list(), list(
      swap_class_2_3         = isTRUE(input$score_swap23),
      offtarget_count_weight = input$score_ot_weight %||% 0,
      activity_term          = isTRUE(input$score_activity)
    ))
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
  results_view_df <- shiny::reactive({
    shiny::req(rv$results_all)
    if (input$results_view == "top") rv$results_top else rv$results_all
  })

  output$results_table <- DT::renderDT({
    df <- results_view_df()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No rows for the current filter."))
    DT::datatable(df, filter = "top", rownames = FALSE, selection = "single",
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  # The design row currently selected in the Results table (for the Inspector)
  selected_design <- shiny::reactive({
    sel <- input$results_table_rows_selected
    df <- results_view_df()
    if (is.null(sel) || length(sel) == 0 || is.null(df) || nrow(df) == 0) return(NULL)
    df[sel[1], , drop = FALSE]
  })

  shiny::observeEvent(input$results_table_rows_selected, {
    if (length(input$results_table_rows_selected)) bslib::nav_select("main_tabs", "Inspector")
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

  # ==== Inspector (select a design in Results) ===============================
  output$inspector_ready <- shiny::reactive(!is.null(selected_design()))
  shiny::outputOptions(output, "inspector_ready", suspendWhenHidden = FALSE)

  output$inspector_header <- shiny::renderUI({
    d <- selected_design()
    if (is.null(d)) {
      return(shiny::p("Select a design row in the ", shiny::strong("Results"),
                      " tab to inspect it — genomic map, sequence, and off-targets."))
    }
    shiny::h5(sprintf("%s  —  %s→%s in %s  (rank %s, score %s, gRNA on %s strand)",
                      d$mutation, d$original_codon, d$mutant_codon, d$gene, d$rank,
                      d$score, if (identical(d$strand, "+")) "+" else "−"))
  })

  inspector_offtargets <- shiny::reactive({
    d <- selected_design(); shiny::req(d, rv$system)
    tryCatch(design_off_targets(d, rv$system), error = function(e) data.frame())
  })

  output$inspector_map <- shiny::renderPlot({
    d <- selected_design(); shiny::req(d, rv$system)
    render_design_ggmap(d, rv$system)
  })

  output$inspector_seq <- shiny::renderUI({
    d <- selected_design(); shiny::req(d, rv$system)
    shiny::HTML(design_sequence_html(d, rv$system))
  })

  output$inspector_ot_header <- shiny::renderUI({
    ot <- inspector_offtargets()
    n <- if (is.null(ot)) 0 else nrow(ot)
    paste0("Off-target sites (", n, ")")
  })

  output$inspector_offtargets <- DT::renderDT({
    ot <- inspector_offtargets()
    shiny::validate(shiny::need(!is.null(ot) && nrow(ot) > 0,
                                "No off-target sites (≤4 mismatches or shared seed)."))
    DT::datatable(ot, rownames = FALSE, filter = "top", selection = "single",
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$inspector_ot_align <- shiny::renderUI({
    ot <- inspector_offtargets(); d <- selected_design()
    if (is.null(ot) || nrow(ot) == 0 || is.null(d)) return(NULL)
    sel <- input$inspector_offtargets_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      return(shiny::p(shiny::tags$small(shiny::icon("hand-pointer"),
                      " Click an off-target row to align it to the guide (mismatches in red).")))
    }
    shiny::div(class = "mb-3 p-2", style = "background:#f8f9fa;border-radius:6px",
               shiny::HTML(offtarget_alignment_html(d$protospacer, ot[sel[1], , drop = FALSE])))
  })

  # ==== Tracker (SQLite edit database) =======================================

  # Require a connected DB; notify and return FALSE otherwise.
  tracker_require <- function() {
    if (!isTRUE(rv$tracker_ready) || is.null(rv$tracker_path)) {
      shiny::showNotification("Connect to a tracker database first (Tracker tab).",
                              type = "warning")
      return(FALSE)
    }
    TRUE
  }

  shiny::observeEvent(input$tracker_connect, {
    path <- trimws(input$tracker_path %||% "")
    if (!nzchar(path)) { shiny::showNotification("Enter a database file path.", type = "warning"); return() }
    ok <- tryCatch({ tracker_init(path); TRUE },
                   error = function(e) { shiny::showNotification(
                     paste("Could not open database:", conditionMessage(e)),
                     type = "error", duration = NULL); FALSE })
    if (ok) {
      rv$tracker_path <- path; rv$tracker_ready <- TRUE
      rv$tracker_tick <- rv$tracker_tick + 1
      shiny::showNotification(paste("Tracker ready:", path), type = "message")
    }
  })

  output$tracker_status <- shiny::renderUI({
    rv$tracker_tick
    if (!isTRUE(rv$tracker_ready)) {
      shiny::div(class = "alert alert-secondary py-2",
        shiny::icon("circle-info"),
        " No database connected. Set a path and click ", shiny::strong("Connect / create"),
        ". A new .sqlite file is created if it doesn't exist.")
    } else {
      e <- tryCatch(tracker_edits(rv$tracker_path), error = function(err) data.frame())
      s <- tryCatch(tracker_strains(rv$tracker_path), error = function(err) data.frame())
      shiny::div(class = "alert alert-success py-2",
        shiny::icon("database"), " Connected to ", shiny::code(rv$tracker_path),
        sprintf(" — %d edit(s), %d strain(s).", nrow(e), nrow(s)))
    }
  })

  # --- Bulk-log the current run's best design per mutation --------------------
  shiny::observeEvent(input$tracker_log_run, {
    if (!tracker_require()) return()
    if (is.null(rv$results_top) || nrow(rv$results_top) == 0) {
      shiny::showNotification("No results to log — run the design pipeline first.",
                              type = "warning"); return() }
    cfg <- effective_cfg()
    method <- cfg$microsynth$cloning_method %||% "golden_gate"
    na_mM  <- cfg$qc$na_mM %||% 50
    genome <- rv$system$genome_name %||% NA
    df <- rv$results_top
    # one (best-scoring) design per mutation
    best <- do.call(rbind, lapply(split(seq_len(nrow(df)), df$mutation), function(ix) {
      sub <- df[ix, , drop = FALSE]; sub[which.max(sub$score), , drop = FALSE] }))
    n <- 0
    withCallingHandlers(
      for (i in seq_len(nrow(best))) {
        row <- best[i, , drop = FALSE]
        oligos <- tryCatch(build_oligo_orders(row, cfg)$oligos, error = function(e) NULL)
        tryCatch({ tracker_log_edit(rv$tracker_path, row, oligos, genome = genome,
                                    cloning_method = method, na_mM = na_mM); n <- n + 1 },
                 error = function(e) shiny::showNotification(
                   paste("Skipped", row$mutation, "-", conditionMessage(e)), type = "error"))
      },
      warning = function(w) invokeRestart("muffleWarning"))
    rv$tracker_tick <- rv$tracker_tick + 1
    shiny::showNotification(sprintf("Logged %d edit(s) to the tracker.", n), type = "message")
  })

  # --- Log the single design currently open in the Inspector -----------------
  shiny::observeEvent(input$inspector_log_btn, {
    if (!tracker_require()) return()
    d <- selected_design()
    if (is.null(d)) { shiny::showNotification("No design selected.", type = "warning"); return() }
    cfg <- effective_cfg()
    oligos <- tryCatch(build_oligo_orders(d, cfg)$oligos, error = function(e) NULL)
    ok <- tryCatch({
      tracker_log_edit(rv$tracker_path, d, oligos,
                       genome = rv$system$genome_name %||% NA,
                       cloning_method = cfg$microsynth$cloning_method %||% "golden_gate",
                       status = "chosen", na_mM = cfg$qc$na_mM %||% 50); TRUE },
      error = function(e) { shiny::showNotification(
        paste("Log failed:", conditionMessage(e)), type = "error"); FALSE })
    if (ok) {
      rv$tracker_tick <- rv$tracker_tick + 1
      shiny::showNotification(paste("Logged", d$mutation, "to the tracker."), type = "message")
    }
  })

  # --- Table outputs (re-read on every write via tracker_tick) ---------------
  .tracker_dt <- function(df, empty_msg) {
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, empty_msg))
    DT::datatable(df, rownames = FALSE, selection = "single", filter = "top",
                  options = list(pageLength = 10, scrollX = TRUE))
  }

  output$tracker_edits_tbl <- DT::renderDT({
    rv$tracker_tick; shiny::req(rv$tracker_ready)
    .tracker_dt(tracker_edits(rv$tracker_path), "No edits logged yet.")
  })
  output$tracker_attempts_tbl <- DT::renderDT({
    rv$tracker_tick; shiny::req(rv$tracker_ready)
    .tracker_dt(tracker_attempts(rv$tracker_path), "No attempts recorded yet.")
  })
  output$tracker_curing_tbl <- DT::renderDT({
    rv$tracker_tick; shiny::req(rv$tracker_ready)
    .tracker_dt(tracker_curing(rv$tracker_path), "No curing outcomes recorded yet.")
  })
  output$tracker_curing_summary <- DT::renderDT({
    rv$tracker_tick; shiny::req(rv$tracker_ready)
    df <- tracker_curing_by_strain(rv$tracker_path)
    shiny::validate(shiny::need(nrow(df) > 0, "No curing data yet."))
    DT::datatable(df, rownames = FALSE, options = list(dom = "t"))
  })
  output$tracker_strains_tbl <- DT::renderDT({
    rv$tracker_tick; shiny::req(rv$tracker_ready)
    .tracker_dt(tracker_strains(rv$tracker_path), "No strains added yet.")
  })

  # --- CSV exports -----------------------------------------------------------
  .tracker_dl <- function(tag, getter) shiny::downloadHandler(
    filename = function() sprintf("tracker_%s_%s.csv", tag, Sys.Date()),
    content = function(file) {
      df <- if (isTRUE(rv$tracker_ready))
        tryCatch(getter(rv$tracker_path), error = function(e) data.frame()) else data.frame()
      utils::write.csv(df, file, row.names = FALSE, na = "")
    })
  output$tracker_dl_flat     <- .tracker_dl("edits_outcomes", tracker_export_flat)
  output$tracker_dl_edits    <- .tracker_dl("edits",    tracker_edits)
  output$tracker_dl_attempts <- .tracker_dl("attempts", tracker_attempts)
  output$tracker_dl_curing   <- .tracker_dl("curing",   tracker_curing)
  output$tracker_dl_strains  <- .tracker_dl("strains",  tracker_strains)

  # --- Edit detail: oligos for the selected edit -----------------------------
  output$tracker_edit_detail <- shiny::renderUI({
    rv$tracker_tick; shiny::req(rv$tracker_ready)
    sel <- input$tracker_edits_tbl_rows_selected
    ed <- tracker_edits(rv$tracker_path)
    if (is.null(sel) || length(sel) == 0 || nrow(ed) == 0) return(NULL)
    eid <- ed$edit_id[sel[1]]
    oligos <- tracker_oligos_for_edit(rv$tracker_path, eid)
    bslib::card(class = "mt-3",
      bslib::card_header(sprintf("Oligos for edit #%d (%s)", eid, ed$mutation[sel[1]])),
      shiny::div(class = "mb-2",
        shiny::actionButton("tracker_attempt_for_edit", "Record attempt for this edit",
                            class = "btn-outline-primary btn-sm", icon = shiny::icon("plus"))),
      if (nrow(oligos) == 0) shiny::em("No oligos linked.")
      else shiny::renderTable(oligos, striped = TRUE, spacing = "xs"))
  })

  # --- Modal forms -----------------------------------------------------------
  shiny::observeEvent(input$tracker_add_strain_btn, {
    if (!tracker_require()) return()
    shiny::showModal(shiny::modalDialog(
      title = "Add / edit strain", easyClose = TRUE,
      shiny::textInput("f_strain_name", "Name (unique) *"),
      shiny::textInput("f_strain_parent", "Parent strain"),
      shiny::textInput("f_strain_background", "Background (e.g. MG1655, BW25113)"),
      shiny::textInput("f_strain_genotype", "Genotype"),
      shiny::textAreaInput("f_strain_notes", "Notes"),
      footer = shiny::tagList(shiny::modalButton("Cancel"),
                              shiny::actionButton("f_strain_save", "Save", class = "btn-primary"))))
  })
  shiny::observeEvent(input$f_strain_save, {
    ok <- tryCatch({
      tracker_add_strain(rv$tracker_path, input$f_strain_name, input$f_strain_parent,
                         input$f_strain_background, input$f_strain_genotype, input$f_strain_notes); TRUE },
      error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); FALSE })
    if (ok) { shiny::removeModal(); rv$tracker_tick <- rv$tracker_tick + 1
      shiny::showNotification("Strain saved.", type = "message") }
  })

  attempt_modal <- function(preselect_edit = NULL) {
    shiny::showModal(shiny::modalDialog(
      title = "Record transformation attempt", easyClose = TRUE,
      shiny::selectInput("f_at_edit", "Edit *", choices = tracker_edit_choices(rv$tracker_path),
                         selected = preselect_edit),
      shiny::selectInput("f_at_strain", "Target strain",
                         choices = c("(none)" = "", tracker_strain_choices(rv$tracker_path))),
      shiny::textInput("f_at_date", "Date", value = as.character(Sys.Date())),
      shiny::selectInput("f_at_method", "Method",
                         choices = c("pEcCas/pEcgRNA (lambda-Red)", "pEcCas only", "other")),
      shiny::selectInput("f_at_transformation", "Transformation success",
                         choices = c("pending", "success", "low efficiency", "failed")),
      shiny::numericInput("f_at_colonies", "Colonies", value = NA, min = 0),
      shiny::selectInput("f_at_screen", "Edit screen result",
                         choices = c("pending", "edited", "not edited", "mixed")),
      shiny::selectInput("f_at_seq", "Sequencing confirmed",
                         choices = c("pending", "yes", "no")),
      shiny::textAreaInput("f_at_notes", "Notes"),
      footer = shiny::tagList(shiny::modalButton("Cancel"),
                              shiny::actionButton("f_at_save", "Save", class = "btn-primary"))))
  }
  shiny::observeEvent(input$tracker_add_attempt_btn, { if (tracker_require()) attempt_modal() })
  shiny::observeEvent(input$tracker_attempt_for_edit, {
    if (!tracker_require()) return()
    sel <- input$tracker_edits_tbl_rows_selected
    ed <- tracker_edits(rv$tracker_path)
    eid <- if (!is.null(sel) && length(sel)) ed$edit_id[sel[1]] else NULL
    attempt_modal(preselect_edit = eid)
  })
  shiny::observeEvent(input$f_at_save, {
    if (is.null(input$f_at_edit) || !nzchar(input$f_at_edit)) {
      shiny::showNotification("Pick an edit.", type = "warning"); return() }
    ok <- tryCatch({
      tracker_add_attempt(rv$tracker_path, edit_id = input$f_at_edit,
                          strain_id = input$f_at_strain, date = input$f_at_date,
                          method = input$f_at_method, transformation = input$f_at_transformation,
                          n_colonies = input$f_at_colonies,
                          screen_result = input$f_at_screen, seq_confirmed = input$f_at_seq,
                          notes = input$f_at_notes); TRUE },
      error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); FALSE })
    if (ok) { shiny::removeModal(); rv$tracker_tick <- rv$tracker_tick + 1
      shiny::showNotification("Attempt recorded.", type = "message") }
  })

  shiny::observeEvent(input$tracker_add_curing_btn, {
    if (!tracker_require()) return()
    ch <- tracker_attempt_choices(rv$tracker_path)
    if (length(ch) == 0) { shiny::showNotification(
      "Record an attempt first — curing attaches to an attempt.", type = "warning"); return() }
    shiny::showModal(shiny::modalDialog(
      title = "Record plasmid curing", easyClose = TRUE,
      shiny::selectInput("f_cu_attempt", "Attempt *", choices = ch),
      shiny::selectInput("f_cu_plasmid", "Plasmid",
                         choices = c("pEcgRNA", "pEcCas", "both")),
      shiny::selectInput("f_cu_method", "Method",
                         choices = c("SacB / sucrose counter-selection", "temperature (repA-ts)", "other")),
      shiny::numericInput("f_cu_passages", "Passages", value = NA, min = 0),
      shiny::selectInput("f_cu_cured", "Cured?", choices = c("yes", "no", "partial")),
      shiny::textInput("f_cu_date", "Date", value = as.character(Sys.Date())),
      shiny::textAreaInput("f_cu_notes", "Notes"),
      footer = shiny::tagList(shiny::modalButton("Cancel"),
                              shiny::actionButton("f_cu_save", "Save", class = "btn-primary"))))
  })
  shiny::observeEvent(input$f_cu_save, {
    ok <- tryCatch({
      tracker_add_curing(rv$tracker_path, attempt_id = input$f_cu_attempt,
                         plasmid = input$f_cu_plasmid, method = input$f_cu_method,
                         passages = input$f_cu_passages, cured = input$f_cu_cured,
                         date = input$f_cu_date, notes = input$f_cu_notes); TRUE },
      error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); FALSE })
    if (ok) { shiny::removeModal(); rv$tracker_tick <- rv$tracker_tick + 1
      shiny::showNotification("Curing outcome recorded.", type = "message") }
  })

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
          top-ranked design of each mutation, with 3 oligos each, deduplicated by
          sequence. By default these are the <b>Golden Gate</b> spacer pair
          (<code>_gg_bsaI_top</code>/<code>_bot</code>) plus the <code>_HA</code>
          repair donor; switch to Gibson or Both under
          <b>Settings &rarr; gRNA cloning oligos</b>.</p>
       <h4>Oligo size</h4>
       <p>The <b>Microsynth oligo size</b> (Settings) sets the homology-arm repair
          oligo length; default 100 bp. It is also configurable in
          <code>config/config.yml</code>.</p>
       <h4>Cloning method</h4>
       <p>Under <b>Settings &rarr; gRNA cloning oligos</b> choose how the spacer is
          delivered: <b>Golden Gate</b> (the default — an annealed oligo pair
          <code>_gg_bsaI_top</code>/<code>_bot</code> with 4-nt BsaI fusion
          overhangs, default <code>TAGT</code>/<code>AAAC</code>, optionally
          5'-phosphorylated), <b>Gibson</b> (whole-plasmid PCR primers,
          <code>_fw_gibson</code>/<code>_rv_gibson</code>), or <b>Both</b>. The
          <code>_HA</code> repair donor is always included.</p>
       <h4>Tracker</h4>
       <p>The <b>Tracker</b> tab keeps a local SQLite database of your edits:
          which strain each went into, the oligos ordered, transformation attempts
          and plasmid-curing outcomes. Set a database file and click
          <b>Connect / create</b>.</p>
       <p><b>Adding designs.</b> Two ways: <b>Log all top picks (best per
          mutation)</b> on the Tracker tab logs the whole run at once (status
          <i>designed</i>); or select a design in <b>Results</b> &rarr; inspect it,
          and click <b>Log this chosen design</b> in the Inspector to record just
          that one (status <i>chosen</i>, which is never downgraded by a later bulk
          log). Re-logging the same design never duplicates it.</p>
       <p><b>Recording outcomes.</b> On the <b>Attempts</b> tab (or the
          <i>Record attempt for this edit</i> button under an edit) log each
          transformation: pick the edit and strain, and record
          <b>Transformation success</b> (did the transformation/recombineering
          work) separately from the <b>Edit screen result</b> and sequencing. On the
          <b>Curing</b> tab, attach a plasmid-curing outcome to an attempt; the
          <i>Curing success by strain</i> table summarises it.</p>
       <p><b>Export.</b> The <b>Export CSV</b> buttons download any table, or a
          single flat <i>Edits + outcomes</i> sheet (edits joined to attempts and
          curing) for sharing. The <code>.sqlite</code> file itself is one portable
          file you can back up or open in any SQLite tool.</p>
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
