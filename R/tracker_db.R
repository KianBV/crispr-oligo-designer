# ==============================================================================
# tracker_db.R  --  SQLite edit-tracking database
# ==============================================================================
# A small relational store for the wet-lab history of each edit: which strain it
# went into, which oligos were ordered, transformation attempts and plasmid
# curing outcomes. One self-contained .sqlite file, queried from R via DBI.
#
# Schema (normalised so one edit -> many oligos and one edit -> many attempts):
#   strains      one row per host strain
#   edits        one row per designed edit (mutation + chosen gRNA/HA design)
#   oligos       one row per unique ordered oligo
#   edit_oligos  link table: which oligos belong to which edit (+ role)
#   attempts     one transformation attempt of an edit into a strain
#   curing       plasmid-curing outcome for an attempt
#
# Everything here opens a connection, does its work, and closes it, so there is
# no long-lived handle to leak or lock. SQLite open/close is cheap.
# ==============================================================================

#' Default location of the tracker database file.
tracker_default_path <- function(cfg = NULL) {
  p <- tryCatch(cfg$tracker$db_path, error = function(e) NULL)
  if (is.null(p) || !nzchar(p)) p <- file.path("tracker", "edits.sqlite")
  p
}

#' Open a connection with foreign keys enabled. Caller must dbDisconnect().
.tracker_open <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  DBI::dbExecute(con, "PRAGMA foreign_keys = ON;")
  con
}

#' Run fun(con) with an open connection, always closing it afterwards.
tracker_with <- function(path, fun) {
  con <- .tracker_open(path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  fun(con)
}

# ------------------------------------------------------------------------------
# Schema creation
# ------------------------------------------------------------------------------

#' Create the tracker schema if it does not yet exist. Idempotent.
tracker_init <- function(path) {
  tracker_with(path, function(con) {
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS strains (
        strain_id   INTEGER PRIMARY KEY AUTOINCREMENT,
        name        TEXT UNIQUE NOT NULL,
        parent      TEXT,
        background  TEXT,
        genotype    TEXT,
        notes       TEXT,
        created     TEXT DEFAULT (datetime('now'))
      );")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS edits (
        edit_id       INTEGER PRIMARY KEY AUTOINCREMENT,
        mutation      TEXT NOT NULL,
        gene          TEXT,
        position      INTEGER,
        wt_aa         TEXT,
        mut_aa        TEXT,
        genome        TEXT,
        score         REAL,
        strand        TEXT,
        protospacer   TEXT,
        pam_mutation  TEXT,
        homology_arm  TEXT,
        cloning_method TEXT,
        designer      TEXT,
        date_designed TEXT DEFAULT (datetime('now')),
        status        TEXT DEFAULT 'designed',
        notes         TEXT,
        UNIQUE(mutation, protospacer, genome)
      );")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS oligos (
        oligo_id   INTEGER PRIMARY KEY AUTOINCREMENT,
        name       TEXT,
        sequence   TEXT NOT NULL,
        role       TEXT,
        mod5       TEXT,
        length     INTEGER,
        tm         REAL,
        order_ref  TEXT,
        order_date TEXT,
        created    TEXT DEFAULT (datetime('now')),
        UNIQUE(name, sequence)
      );")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS edit_oligos (
        edit_id  INTEGER NOT NULL REFERENCES edits(edit_id) ON DELETE CASCADE,
        oligo_id INTEGER NOT NULL REFERENCES oligos(oligo_id) ON DELETE CASCADE,
        role     TEXT,
        PRIMARY KEY (edit_id, oligo_id)
      );")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS attempts (
        attempt_id     INTEGER PRIMARY KEY AUTOINCREMENT,
        edit_id        INTEGER NOT NULL REFERENCES edits(edit_id) ON DELETE CASCADE,
        strain_id      INTEGER REFERENCES strains(strain_id) ON DELETE SET NULL,
        date           TEXT,
        method         TEXT,
        transformation TEXT,
        n_colonies     INTEGER,
        screen_result  TEXT,
        seq_confirmed  TEXT,
        notes          TEXT,
        created        TEXT DEFAULT (datetime('now'))
      );")
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS curing (
        curing_id  INTEGER PRIMARY KEY AUTOINCREMENT,
        attempt_id INTEGER NOT NULL REFERENCES attempts(attempt_id) ON DELETE CASCADE,
        plasmid    TEXT,
        method     TEXT,
        passages   INTEGER,
        cured      TEXT,
        date       TEXT,
        notes      TEXT,
        created    TEXT DEFAULT (datetime('now'))
      );")
    .tracker_migrate(con)
    invisible(TRUE)
  })
}

#' Add columns introduced after a database was first created. SQLite has no
#' "ADD COLUMN IF NOT EXISTS", so we check PRAGMA table_info first. Idempotent.
.tracker_migrate <- function(con) {
  has_col <- function(tbl, col) {
    info <- DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s);", tbl))
    col %in% info$name
  }
  if (!has_col("attempts", "transformation"))
    DBI::dbExecute(con, "ALTER TABLE attempts ADD COLUMN transformation TEXT;")
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Writes
# ------------------------------------------------------------------------------

.na_if_empty <- function(x) if (is.null(x) || length(x) == 0 || !nzchar(trimws(as.character(x)))) NA else x

#' Insert a strain (or update its editable fields if the name already exists).
tracker_add_strain <- function(path, name, parent = NA, background = NA,
                               genotype = NA, notes = NA) {
  name <- trimws(name)
  if (!nzchar(name)) stop("Strain name is required.", call. = FALSE)
  tracker_with(path, function(con) {
    DBI::dbExecute(con,
      "INSERT INTO strains (name, parent, background, genotype, notes)
       VALUES (:n, :p, :b, :g, :o)
       ON CONFLICT(name) DO UPDATE SET
         parent=excluded.parent, background=excluded.background,
         genotype=excluded.genotype, notes=excluded.notes;",
      params = list(n = name, p = .na_if_empty(parent), b = .na_if_empty(background),
                    g = .na_if_empty(genotype), o = .na_if_empty(notes)))
    DBI::dbGetQuery(con, "SELECT strain_id FROM strains WHERE name = :n",
                    params = list(n = name))$strain_id[1]
  })
}

#' Derive a short role label for an oligo from its name suffix (after the
#' sanitised mutation label), e.g. "pEcgRNA_fw_gibson", "GG_top", "HA".
.oligo_role <- function(oligo_name, mutation) {
  lbl <- gsub("[^A-Za-z0-9_]", "_", mutation)
  role <- sub(paste0("^", lbl), "", oligo_name)
  role <- sub("^_+", "", role)
  if (!nzchar(role)) role <- NA
  role
}

#' Log one designed edit and its oligos. `design_row` is a single row of the
#' results_top data.frame; `oligos_df` has columns oligo_name, sequence, mod5
#' (the oligos for THIS edit, from generate_oligo_orders on the single row).
#' Returns the edit_id. Re-logging the same (mutation, protospacer, genome) is a
#' no-op on the edit but still (re)links any new oligos.
tracker_log_edit <- function(path, design_row, oligos_df, genome = NA,
                             designer = NA, cloning_method = NA,
                             status = "designed", na_mM = 50) {
  d <- design_row
  mut <- as.character(d$mutation)
  proto <- toupper(as.character(d$protospacer %||% ""))
  gen <- .na_if_empty(genome)

  tracker_with(path, function(con) {
    # --- edit (idempotent on the natural key) ---
    DBI::dbExecute(con,
      "INSERT INTO edits
         (mutation, gene, position, wt_aa, mut_aa, genome, score, strand,
          protospacer, pam_mutation, homology_arm, cloning_method, designer, status)
       VALUES (:mut,:gene,:pos,:wt,:mt,:gen,:sc,:st,:pr,:pm,:ha,:cm,:ds,:stat)
       ON CONFLICT(mutation, protospacer, genome) DO UPDATE SET
         score=excluded.score, pam_mutation=excluded.pam_mutation,
         homology_arm=excluded.homology_arm, cloning_method=excluded.cloning_method,
         -- upgrade to 'chosen' when logged as such; never downgrade back
         status=CASE WHEN excluded.status='chosen' THEN 'chosen' ELSE edits.status END;",
      params = list(
        mut = mut, gene = .na_if_empty(d$gene), pos = suppressWarnings(as.integer(d$position)),
        wt = .na_if_empty(d$wt_aa), mt = .na_if_empty(d$mut_aa), gen = gen,
        sc = suppressWarnings(as.numeric(d$score)), st = .na_if_empty(d$strand),
        pr = .na_if_empty(proto), pm = .na_if_empty(d$pam_mutation),
        ha = .na_if_empty(toupper(as.character(d$homology_arm %||% ""))),
        cm = .na_if_empty(cloning_method), ds = .na_if_empty(designer),
        stat = status))
    eid <- DBI::dbGetQuery(con,
      "SELECT edit_id FROM edits WHERE mutation=:m AND protospacer=:p AND
         ifnull(genome,'')=ifnull(:g,'')",
      params = list(m = mut, p = proto, g = gen))$edit_id[1]

    # --- oligos + links ---
    if (!is.null(oligos_df) && nrow(oligos_df) > 0) {
      for (i in seq_len(nrow(oligos_df))) {
        nm  <- oligos_df$oligo_name[i]
        seq <- toupper(oligos_df$sequence[i])
        md  <- if ("mod5" %in% names(oligos_df)) oligos_df$mod5[i] else ""
        tm  <- tryCatch(round(tm_nn(seq, na_mM = na_mM), 1), error = function(e) NA_real_)
        role <- .oligo_role(nm, mut)
        DBI::dbExecute(con,
          "INSERT INTO oligos (name, sequence, role, mod5, length, tm)
           VALUES (:n,:s,:r,:m,:l,:t)
           ON CONFLICT(name, sequence) DO UPDATE SET role=excluded.role, mod5=excluded.mod5;",
          params = list(n = nm, s = seq, r = .na_if_empty(role),
                        m = .na_if_empty(md), l = nchar(seq), t = tm))
        oid <- DBI::dbGetQuery(con,
          "SELECT oligo_id FROM oligos WHERE name=:n AND sequence=:s",
          params = list(n = nm, s = seq))$oligo_id[1]
        DBI::dbExecute(con,
          "INSERT OR IGNORE INTO edit_oligos (edit_id, oligo_id, role)
           VALUES (:e,:o,:r);",
          params = list(e = eid, o = oid, r = .na_if_empty(role)))
      }
    }
    eid
  })
}

#' Record a transformation attempt of an edit into a strain. `transformation`
#' captures whether the transformation/recombineering itself worked (e.g.
#' "success", "failed", "low efficiency"); `screen_result`/`seq_confirmed`
#' capture whether the intended edit was obtained.
tracker_add_attempt <- function(path, edit_id, strain_id = NA, date = NA,
                                method = NA, transformation = NA, n_colonies = NA,
                                screen_result = NA, seq_confirmed = NA, notes = NA) {
  tracker_with(path, function(con) {
    DBI::dbExecute(con,
      "INSERT INTO attempts
         (edit_id, strain_id, date, method, transformation, n_colonies, screen_result, seq_confirmed, notes)
       VALUES (:e,:s,:d,:m,:tr,:n,:r,:c,:o);",
      params = list(e = as.integer(edit_id), s = suppressWarnings(as.integer(.na_if_empty(strain_id))),
                    d = .na_if_empty(date), m = .na_if_empty(method), tr = .na_if_empty(transformation),
                    n = suppressWarnings(as.integer(.na_if_empty(n_colonies))),
                    r = .na_if_empty(screen_result), c = .na_if_empty(seq_confirmed),
                    o = .na_if_empty(notes)))
    DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id[1]
  })
}

#' Record a plasmid-curing outcome for an attempt.
tracker_add_curing <- function(path, attempt_id, plasmid = NA, method = NA,
                               passages = NA, cured = NA, date = NA, notes = NA) {
  tracker_with(path, function(con) {
    DBI::dbExecute(con,
      "INSERT INTO curing (attempt_id, plasmid, method, passages, cured, date, notes)
       VALUES (:a,:p,:m,:pa,:c,:d,:o);",
      params = list(a = as.integer(attempt_id), p = .na_if_empty(plasmid),
                    m = .na_if_empty(method), pa = suppressWarnings(as.integer(.na_if_empty(passages))),
                    c = .na_if_empty(cured), d = .na_if_empty(date), o = .na_if_empty(notes)))
    DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id[1]
  })
}

# ------------------------------------------------------------------------------
# Reads (return plain data.frames for DT)
# ------------------------------------------------------------------------------

tracker_edits <- function(path) {
  tracker_with(path, function(con) {
    DBI::dbGetQuery(con,
      "SELECT e.edit_id, e.mutation, e.gene, e.genome, e.strand, e.score,
              e.cloning_method, e.status, e.date_designed,
              (SELECT COUNT(*) FROM edit_oligos eo WHERE eo.edit_id=e.edit_id) AS n_oligos,
              (SELECT COUNT(*) FROM attempts a WHERE a.edit_id=e.edit_id)      AS n_attempts
       FROM edits e ORDER BY e.edit_id DESC;")
  })
}

tracker_strains <- function(path) {
  tracker_with(path, function(con)
    DBI::dbGetQuery(con, "SELECT strain_id, name, parent, background, genotype, notes
                          FROM strains ORDER BY name;"))
}

tracker_oligos_for_edit <- function(path, edit_id) {
  tracker_with(path, function(con)
    DBI::dbGetQuery(con,
      "SELECT o.name, o.role, o.sequence, o.mod5, o.length, o.tm, o.order_ref, o.order_date
       FROM oligos o JOIN edit_oligos eo ON eo.oligo_id=o.oligo_id
       WHERE eo.edit_id=:e ORDER BY o.role;",
      params = list(e = as.integer(edit_id))))
}

tracker_attempts <- function(path) {
  tracker_with(path, function(con)
    DBI::dbGetQuery(con,
      "SELECT a.attempt_id, e.mutation, s.name AS strain, a.date, a.method,
              a.transformation, a.n_colonies, a.screen_result, a.seq_confirmed, a.notes
       FROM attempts a
       JOIN edits e ON e.edit_id=a.edit_id
       LEFT JOIN strains s ON s.strain_id=a.strain_id
       ORDER BY a.attempt_id DESC;"))
}

tracker_curing <- function(path) {
  tracker_with(path, function(con)
    DBI::dbGetQuery(con,
      "SELECT c.curing_id, e.mutation, s.name AS strain, c.plasmid, c.method,
              c.passages, c.cured, c.date, c.notes
       FROM curing c
       JOIN attempts a ON a.attempt_id=c.attempt_id
       JOIN edits e ON e.edit_id=a.edit_id
       LEFT JOIN strains s ON s.strain_id=a.strain_id
       ORDER BY c.curing_id DESC;"))
}

#' Curing success rate grouped by strain background (for a quick at-a-glance).
tracker_curing_by_strain <- function(path) {
  tracker_with(path, function(con)
    DBI::dbGetQuery(con,
      "SELECT COALESCE(s.name,'(unknown)') AS strain,
              COUNT(*) AS n,
              SUM(CASE WHEN lower(c.cured) IN ('yes','y','true','1') THEN 1 ELSE 0 END) AS cured,
              ROUND(100.0*SUM(CASE WHEN lower(c.cured) IN ('yes','y','true','1') THEN 1 ELSE 0 END)/COUNT(*),1) AS pct_cured
       FROM curing c
       JOIN attempts a ON a.attempt_id=c.attempt_id
       LEFT JOIN strains s ON s.strain_id=a.strain_id
       GROUP BY strain ORDER BY n DESC;"))
}

#' One wide, flat table joining edits -> attempts -> curing (+ strain), for
#' CSV export / sharing. LEFT JOINs so an edit with no attempt, or an attempt
#' with no curing record, still appears (outcome columns NA).
tracker_export_flat <- function(path) {
  tracker_with(path, function(con)
    DBI::dbGetQuery(con,
      "SELECT e.edit_id, e.mutation, e.gene, e.wt_aa, e.mut_aa, e.genome,
              e.strand, e.score, e.cloning_method, e.status AS design_status,
              e.protospacer, e.pam_mutation, e.date_designed,
              a.attempt_id, s.name AS strain, a.date AS attempt_date, a.method,
              a.transformation, a.n_colonies, a.screen_result, a.seq_confirmed,
              a.notes AS attempt_notes,
              c.plasmid AS cured_plasmid, c.method AS cure_method,
              c.passages, c.cured, c.date AS cure_date
       FROM edits e
       LEFT JOIN attempts a ON a.edit_id=e.edit_id
       LEFT JOIN strains  s ON s.strain_id=a.strain_id
       LEFT JOIN curing   c ON c.attempt_id=a.attempt_id
       ORDER BY e.edit_id, a.attempt_id, c.curing_id;"))
}

#' Named table-getters keyed by a short label, used by the CSV export UI.
tracker_table_exporters <- function() {
  list(edits = tracker_edits, attempts = tracker_attempts, curing = tracker_curing,
       strains = tracker_strains, flat = tracker_export_flat)
}

#' Convenience: (edit_id, mutation) choices for selectInputs.
tracker_edit_choices <- function(path) {
  e <- tracker_edits(path)
  if (nrow(e) == 0) return(stats::setNames(character(0), character(0)))
  stats::setNames(e$edit_id, sprintf("#%d  %s", e$edit_id, e$mutation))
}

tracker_strain_choices <- function(path) {
  s <- tracker_strains(path)
  if (nrow(s) == 0) return(stats::setNames(character(0), character(0)))
  stats::setNames(s$strain_id, s$name)
}

tracker_attempt_choices <- function(path) {
  a <- tracker_attempts(path)
  if (nrow(a) == 0) return(stats::setNames(character(0), character(0)))
  stats::setNames(a$attempt_id,
                  sprintf("#%d  %s -> %s (%s)", a$attempt_id, a$mutation,
                          a$strain %||% "?", a$date %||% ""))
}
