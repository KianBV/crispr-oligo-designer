# ==============================================================================
# constants.R  --  Biological constants and tunable defaults
# ==============================================================================
#
# Values in the "TUNABLE" block below are DEFAULTS. They can be overridden at
# runtime by config/config.yml (see config.R) or by the GUI. The core pipeline
# functions read these as global variables, so apply_config() reassigns them in
# the global environment. This keeps the finalised core logic untouched while
# still making the key parameters configurable.
# ==============================================================================

# --- Null-coalescing operator (handles NULL / empty / NA) ---------------------
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

# ------------------------------------------------------------------------------
# Reference genome. The default is E. coli BW25113 (NCBI CP009273.1), but any
# strain can be selected by accession, or a custom genome uploaded (see genome.R
# and config.yml -> genome.strains).
# ------------------------------------------------------------------------------
NCBI_EFETCH             <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_GENOME_ACCESSION <- "CP009273.1"   # E. coli BW25113 (Keio parent)
GENOME_CACHE_DIR         <- file.path(Sys.getenv("HOME"), ".ecoli_grna_cache")

# Build NCBI E-utilities download URLs from a nuccore accession.
ncbi_fasta_url <- function(accession) {
  paste0(NCBI_EFETCH, "?db=nuccore&id=", accession, "&rettype=fasta&retmode=text")
}
ncbi_gff_url <- function(accession) {
  paste0(NCBI_EFETCH, "?db=nuccore&id=", accession, "&rettype=gff3&retmode=text")
}
# Feature table: the officially supported efetch annotation format (rettype=ft).
# More reliable than gff3 (which returns HTTP 500 for many accessions) and richer
# (includes locus_tag, protein_id and UniProt db_xref).
ncbi_ft_url <- function(accession) {
  paste0(NCBI_EFETCH, "?db=nuccore&id=", accession, "&rettype=ft&retmode=text")
}

# Fallback list of standard E. coli strains (label -> NCBI accession). This is
# overridden by config.yml -> genome.strains, which is the place users should
# edit to add/correct strains.
STANDARD_STRAINS_DEFAULT <- list(
  "BW25113 (Keio parent, CP009273.1)" = "CP009273.1",
  "MG1655 (K-12, U00096.3)"           = "U00096.3",
  "W3110 (K-12, AP009048.1)"          = "AP009048.1",
  "DH5alpha (CP026085.1)"             = "CP026085.1",
  "BL21(DE3) (CP001509.3)"            = "CP001509.3",
  "Nissle 1917 (CP007799.1)"          = "CP007799.1"
)

# ------------------------------------------------------------------------------
# Codon tables (standard genetic code)
# ------------------------------------------------------------------------------
CODON_TABLE <- list(
  "F" = c("TTT", "TTC"), "L" = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
  "I" = c("ATT", "ATC", "ATA"), "M" = c("ATG"),
  "V" = c("GTT", "GTC", "GTA", "GTG"), "S" = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
  "P" = c("CCT", "CCC", "CCA", "CCG"), "T" = c("ACT", "ACC", "ACA", "ACG"),
  "A" = c("GCT", "GCC", "GCA", "GCG"), "Y" = c("TAT", "TAC"),
  "H" = c("CAT", "CAC"), "Q" = c("CAA", "CAG"),
  "N" = c("AAT", "AAC"), "K" = c("AAA", "AAG"),
  "D" = c("GAT", "GAC"), "E" = c("GAA", "GAG"),
  "C" = c("TGT", "TGC"), "W" = c("TGG"),
  "R" = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
  "G" = c("GGT", "GGC", "GGA", "GGG"), "*" = c("TAA", "TAG", "TGA")
)

REVERSE_CODON_TABLE <- unlist(lapply(names(CODON_TABLE), function(aa) {
  setNames(rep(aa, length(CODON_TABLE[[aa]])), CODON_TABLE[[aa]])
}))

# Fallback E. coli codon usage frequencies (used before the genome-derived
# frequencies are computed at initialisation).
ECOLI_CODON_FREQ <- c(
  "TTT" = 0.58, "TTC" = 0.42,
  "TTA" = 0.14, "TTG" = 0.13, "CTT" = 0.12, "CTC" = 0.10, "CTA" = 0.04, "CTG" = 0.47,
  "ATT" = 0.51, "ATC" = 0.42, "ATA" = 0.08,
  "ATG" = 1.00,
  "GTT" = 0.26, "GTC" = 0.22, "GTA" = 0.15, "GTG" = 0.37,
  "TCT" = 0.17, "TCC" = 0.15, "TCA" = 0.14, "TCG" = 0.14, "AGT" = 0.16, "AGC" = 0.25,
  "CCT" = 0.18, "CCC" = 0.13, "CCA" = 0.20, "CCG" = 0.49,
  "ACT" = 0.19, "ACC" = 0.40, "ACA" = 0.13, "ACG" = 0.27,
  "GCT" = 0.18, "GCC" = 0.26, "GCA" = 0.21, "GCG" = 0.36,
  "TAT" = 0.59, "TAC" = 0.41,
  "CAT" = 0.57, "CAC" = 0.43,
  "CAA" = 0.34, "CAG" = 0.66,
  "AAT" = 0.49, "AAC" = 0.51,
  "AAA" = 0.74, "AAG" = 0.26,
  "GAT" = 0.63, "GAC" = 0.37,
  "GAA" = 0.68, "GAG" = 0.32,
  "TGT" = 0.46, "TGC" = 0.54,
  "TGG" = 1.00,
  "CGT" = 0.36, "CGC" = 0.36, "CGA" = 0.07, "CGG" = 0.11, "AGA" = 0.07, "AGG" = 0.04,
  "GGT" = 0.35, "GGC" = 0.37, "GGA" = 0.13, "GGG" = 0.15,
  "TAA" = 0.61, "TAG" = 0.09, "TGA" = 0.30
)

# ------------------------------------------------------------------------------
# Cloning / cassette fixed sequences
# ------------------------------------------------------------------------------
SPACER1         <- "TCCTCTGGCGGAAAGCCT"
SPACER2         <- "GATC"
J23119_PROMOTER <- "TTGACAGCTAGCTCAGTCCTAGGTATAATACTAGT"
SGRNA_HANDLE    <- "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAG"

# ==============================================================================
# TUNABLE DEFAULTS  (overridable via config/config.yml or the GUI)
# ==============================================================================

# gRNA design parameters
MAX_PAM_DISTANCE    <- 30   # max distance (bp) from PAM to target codon
MIN_SCORE_THRESHOLD <- 3    # minimum cassette score to keep a candidate
SEED_LENGTH         <- 11   # seed length for off-target seed matching
MAX_MISMATCHES      <- 4    # max mismatches counted as off-target

# Homology arm / repair oligo length.  This is the "Microsynth oligo size":
# the length (bp) of the single-stranded repair template oligo that carries the
# desired codon change + silent PAM mutation. Default 100 bp.
HOMOLOGY_ARM_LENGTH <- 100

# gRNA primer 3' homology tails for Gibson assembly into the pEcgRNA vector.
# Each primer is  <protospacer or revcomp> + <tail>  (tail appended at the 3' end).
#   fwd tail = first 24 nt of the sgRNA scaffold (SGRNA_HANDLE)
#   rev tail = reverse complement of the last 24 nt of the J23119 promoter
# Tails are kept lowercase by convention so the homology region is visually
# distinct from the (uppercase) protospacer in the ordered primer.
GRNA_FWD_TAIL <- "gttttagagctagaaatagcaagt"
GRNA_REV_TAIL <- "actagtattatacctaggactgag"

# Oligo name suffixes
GRNA_FWD_SUFFIX <- "_pEcgRNA_fw_gibson"
GRNA_REV_SUFFIX <- "_pEcgRNA_rv_gibson"
HA_SUFFIX       <- "_HA"

# Phosphorothioate (PTO) protection on the homology arm oligo
PTO_N_BASES <- 4  # number of 5' bases protected (asterisks between first N bases)
