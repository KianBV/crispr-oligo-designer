# ==============================================================================
# config.R  --  Load, merge and apply user configuration
# ==============================================================================
#
# Precedence (lowest -> highest):
#   1. built-in defaults (default_config())
#   2. config/config.yml
#   3. overrides passed from the GUI at runtime
#
# apply_design_config() writes the tunable design parameters into the global
# environment so the core pipeline functions (which read them as globals) pick
# up the configured values without any change to their logic.
# ==============================================================================

#' Built-in default configuration (mirrors constants.R defaults).
default_config <- function() {
  list(
    design = list(
      oligo_size          = 100,   # homology arm / repair oligo length (bp)
      max_pam_distance    = 30,
      min_score_threshold = 3,
      seed_length         = 11,
      max_mismatches      = 4,
      top_n               = 3,
      n_cores             = 1
    ),
    microsynth = list(
      template_path   = "templates/MicrosynthUploadFormDNA__2_.xlsx",
      sheet           = "DNA Order",
      # Which gRNA cloning oligos to emit: "gibson", "golden_gate" or "both".
      cloning_method  = "gibson",
      # gRNA Gibson primer 3' homology tails (appended after the protospacer)
      fwd_tail        = "gttttagagctagaaatagcaagt",
      rev_tail        = "actagtattatacctaggactgag",
      fwd_name_suffix = "_pEcgRNA_fw_gibson",
      rev_name_suffix = "_pEcgRNA_rv_gibson",
      # Golden Gate (BsaI) 4-nt fusion overhangs + 5'-phosphate option
      gg_top_overhang = "TAGT",
      gg_bot_overhang = "AAAC",
      gg_top_name_suffix = "_gg_bsaI_top",
      gg_bot_name_suffix = "_gg_bsaI_bot",
      gg_phosphorylate   = FALSE,
      gg_phos_label      = "Phosphate",   # exact Microsynth 5'-mod dropdown term
      ha_name_suffix  = "_HA",
      pto_n_bases     = 4,
      purification  = "DES",
      scale         = "0.04",
      condition     = "Dried",
      datasheet     = "Standard",
      aliquots      = "No"
    ),
    qc = list(
      enabled         = TRUE,
      na_mM           = 50,     # monovalent cation conc for Tm (fold Mg2+ in here)
      oligo_uM        = 0.5,    # strand conc for Tm
      gc_min          = 0.30,
      gc_max          = 0.65,
      tm_min          = 50,
      tm_max          = 75,
      max_homopolymer = 5,
      hairpin_stem    = 6,      # flag hairpin stems >= this many bp
      hairpin_min_loop = 3,
      dimer3          = 4,      # flag 3' self-dimer >= this many bp
      enzymes         = c("BsaI", "BbsI")
    ),
    export = list(
      min_score_filter = 3.0       # default score cut-off for the order form
    ),
    genome = list(
      cache_dir      = "",         # "" -> default ~/.ecoli_grna_cache
      force_download = FALSE,
      default_strain = "BW25113 (Keio parent, CP009273.1)",
      strains        = STANDARD_STRAINS_DEFAULT
    )
  )
}

#' Recursively merge `override` into `base` (override wins).
merge_config <- function(base, override) {
  if (is.null(override)) return(base)
  for (key in names(override)) {
    if (is.list(base[[key]]) && is.list(override[[key]])) {
      base[[key]] <- merge_config(base[[key]], override[[key]])
    } else {
      base[[key]] <- override[[key]]
    }
  }
  base
}

#' Load configuration from a YAML file, merged over the defaults.
#' Missing file -> defaults only (no error).
load_config <- function(path = "config/config.yml") {
  cfg <- default_config()
  if (!is.null(path) && file.exists(path)) {
    file_cfg <- tryCatch(yaml::read_yaml(path), error = function(e) {
      warning("Could not parse ", path, ": ", conditionMessage(e),
              ". Using defaults.", call. = FALSE)
      NULL
    })
    cfg <- merge_config(cfg, file_cfg)
  }
  cfg
}

#' Apply the design-related parameters to the global environment so the core
#' pipeline functions pick them up. Returns the config invisibly.
apply_design_config <- function(cfg) {
  d <- cfg$design
  m <- cfg$microsynth
  g <- cfg$genome

  assign("MAX_PAM_DISTANCE",    as.numeric(d$max_pam_distance),    envir = .GlobalEnv)
  assign("MIN_SCORE_THRESHOLD", as.numeric(d$min_score_threshold), envir = .GlobalEnv)
  assign("SEED_LENGTH",         as.numeric(d$seed_length),         envir = .GlobalEnv)
  assign("MAX_MISMATCHES",      as.numeric(d$max_mismatches),      envir = .GlobalEnv)
  assign("HOMOLOGY_ARM_LENGTH", as.numeric(d$oligo_size),          envir = .GlobalEnv)

  assign("GRNA_FWD_TAIL",   as.character(m$fwd_tail),        envir = .GlobalEnv)
  assign("GRNA_REV_TAIL",   as.character(m$rev_tail),        envir = .GlobalEnv)
  assign("GRNA_FWD_SUFFIX", as.character(m$fwd_name_suffix), envir = .GlobalEnv)
  assign("GRNA_REV_SUFFIX", as.character(m$rev_name_suffix), envir = .GlobalEnv)
  assign("GG_TOP_OVERHANG", as.character(m$gg_top_overhang %||% "TAGT"), envir = .GlobalEnv)
  assign("GG_BOT_OVERHANG", as.character(m$gg_bot_overhang %||% "AAAC"), envir = .GlobalEnv)
  assign("GG_TOP_SUFFIX",   as.character(m$gg_top_name_suffix %||% "_gg_bsaI_top"), envir = .GlobalEnv)
  assign("GG_BOT_SUFFIX",   as.character(m$gg_bot_name_suffix %||% "_gg_bsaI_bot"), envir = .GlobalEnv)
  assign("HA_SUFFIX",       as.character(m$ha_name_suffix),  envir = .GlobalEnv)
  assign("PTO_N_BASES",     as.numeric(m$pto_n_bases),       envir = .GlobalEnv)

  # Optional custom genome cache directory
  if (!is.null(g$cache_dir) && nzchar(g$cache_dir)) {
    assign("GENOME_CACHE_DIR", g$cache_dir, envir = .GlobalEnv)
  }

  invisible(cfg)
}
