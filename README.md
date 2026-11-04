# CRISPR Oligo Designer

A self-contained graphical application for designing CRISPR guide RNAs and
single-stranded repair oligos that introduce **point mutations** in
*E. coli* BW25113, and for generating the **Microsynth DNA order form** for the
best design of each mutation.

It wraps a finalised R pipeline (PAM finding, off-target analysis, scoring,
silent-PAM mutation, homology-arm and oligo assembly) in a browser-based
[Shiny](https://shiny.posit.co/) interface aimed at non-programmers.

---

## What it does

1. **Upload** the mutations you want to make (`.csv` / `.txt` / `.xlsx`, or paste them).
2. **Validate** them — format check plus "does this gene exist in the genome?".
3. **Run** the existing design pipeline.
4. **Display & filter** every design in an interactive table, and **export** to CSV.
5. **Generate the Microsynth order form** automatically for the top-ranked design
   of each mutation (3 oligos each: `gRNA_F`, `gRNA_R`, homology arm), deduplicated.
6. Produce a **run summary** so you can see exactly what was generated.

---

## Installation

You need **R ≥ 4.2** (get it from <https://cran.r-project.org/>). Then, from a
terminal in the application folder:

```bash
# 1. Install dependencies (once). This includes Bioconductor packages and may
#    take several minutes the first time.
Rscript install.R

# 2. Launch the app (opens in your web browser).
Rscript run.R
```

If `Rscript` is not on your PATH, use the full path to it, e.g. on Windows:

```powershell
& "C:\Users\<you>\AppData\Local\Programs\R\R-4.4.2\bin\Rscript.exe" install.R
& "C:\Users\<you>\AppData\Local\Programs\R\R-4.4.2\bin\Rscript.exe" run.R
```

You can also launch it from RStudio: open the folder and click **Run App** on
`app.R`.

### The Microsynth template

Place the official Microsynth order-form template at:

```
templates/MicrosynthUploadFormDNA__2_.xlsx
```

(see [`templates/README.md`](templates/README.md)). The app fills a **copy** and
never edits the original. You can also upload a template from the app's
**Settings** panel to override the configured path for a session.

---

## Using the app

The sidebar is a 3-step workflow:

| Step | Action |
|------|--------|
| **1 · Genome** | Choose a **standard *E. coli* strain** from the dropdown, **or** select **Upload my own genome** and provide a FASTA (+ optional GFF/GTF annotation). Click **Load / initialise genome**. Downloaded strains are cached in `~/.ecoli_grna_cache` and reused on later launches. |
| **2 · Mutations** | Upload a file or paste mutations, then **Validate**. |
| **3 · Run** | Click **Run design pipeline**. Progress is shown; results open automatically. |

### Genome source

- **Standard strain** — pick from the dropdown. The list of strains (label → NCBI
  accession) is defined in `config/config.yml` under `genome.strains`; edit it to
  add, remove or correct strains. The sequence and annotation are downloaded from
  NCBI (E-utilities) by accession and cached per-accession. The annotation uses the
  NCBI **feature table** (`rettype=ft`), which is reliable across accessions —
  unlike the GFF3 endpoint, which returns HTTP 500 for many records.
- **Custom genome** — upload a genome **FASTA** (the first record is used as the
  chromosome) plus, ideally, a **GFF/GTF annotation**. The annotation is what lets
  the app map mutations to genes and read the reading frame — without it, gene
  lookup is skipped and gene-based mutation design won't resolve (you'll get a
  warning). Off-target analysis uses the uploaded genome sequence.

Then use the tabs:

- **Validation** — per-mutation status (valid format, gene found, duplicates).
- **Results** — every design, with per-column filtering; toggle *top-ranked* vs
  *all designs*; download either as CSV.
- **Order form** — preview of the oligo list and one-click **Generate & download
  Microsynth form**.
- **Summary** — a downloadable text report of the run.

### Input format

Each mutation is a string. Two accepted forms:

```
<gene><WT-aa><position><MUT-aa>          e.g. dpsS100E   (compact; gene symbol only)
<identifier>:<WT-aa><position><MUT-aa>   e.g. b0014:S100E, P0A6Y8:S100E   (any identifier)
```

`dpsS100E` = serine 100 → glutamate in gene *dps*. Use the **colon form** for
identifiers that contain digits (locus tags, UniProt IDs, protein IDs), since the
compact form can only be parsed for purely alphabetic gene symbols. The delimiter
may be `:`, `|` or `/`.

- **CSV / Excel**: a column named `mutation` (case-insensitive) is used if
  present, otherwise the first column.
- **Text**: one mutation per line (or comma/space separated).

See [`examples/`](examples/) for ready-made input files.

### Gene identifier types

Genes can be referenced by different identifiers, and **which ones are available
depends on the genome's annotation**. After loading a genome the app shows the
types it found with counts (e.g. `Gene name (4269) · Locus tag (4313) · Protein ID
(4131)`), and the **Mutation identifier type** selector (Step 2) offers exactly
those — plus **Auto**, which tries every type. Pick the type your list uses.

Recognised identifier types (extracted from the GFF/GTF when present):

| Type | Source attribute | Example |
|------|------------------|---------|
| Gene name | `gene` | `dps`, `lacZ` |
| Locus tag | `locus_tag` | `b0344`, `BW25113_0812` |
| Old locus tag | `old_locus_tag` | `b0014` |
| Protein ID | `protein_id` | `AAC73447.1`, `NP_414878` |
| UniProt ID | `db_xref UniProtKB:…` | `P00722` |

> Which types exist depends on how richly the genome is annotated. Well-curated
> genomes (e.g. MG1655, `U00096.3`) expose gene / locus_tag / protein_id / UniProt;
> sparser ones may offer fewer. The app only lists what the chosen genome actually
> contains, with counts, so you always know what you can use.

### Outputs

| Output | Contents |
|--------|----------|
| `generated_oligos_all_<date>.csv` | Every ranked design (protospacer, homology arm, 200-bp cassette, score, off-target class, …). |
| `generated_oligos_top_<date>.csv` | The top-ranked design per mutation, after the score/amino-acid filters. |
| `oligo_list_<date>.csv` | The flat oligo list (name + sequence) that goes into the order form. |
| `MicrosynthOrder_filled_<date>.xlsx` | The filled Microsynth order form. |
| `oligo_qc_<date>.csv` | Per-oligo QC (GC, Tm, hairpin, dimer, enzyme sites, flags). |
| `run_summary_<date>.txt` | Human-readable summary of the run. |

### gRNA cloning oligos

Each mutation always yields the repair donor **`<mutation>_HA`** — the homology-arm
oligo with phosphorothioate (PTO) bonds on the first 4 bases (`X*X*X*X…`). The
**spacer** is delivered by one of two cloning methods
(`microsynth.cloning_method` = `gibson` | `golden_gate` | `both`):

**Gibson** — whole-plasmid PCR primers for the pEcgRNA vector:

- `<mutation>_pEcgRNA_fw_gibson` — protospacer + 24-nt 3′ tail (`gttttagagctagaaatagcaagt`, sgRNA scaffold start)
- `<mutation>_pEcgRNA_rv_gibson` — revcomp(protospacer) + 24-nt 3′ tail (`actagtattatacctaggactgag`, revcomp of J23119 3′ end)

```
_pEcgRNA_fw_gibson : GAACGTTGTGGATGTCCAGCgttttagagctagaaatagcaagt
_pEcgRNA_rv_gibson : GCTGGACATCCACAACGTTCactagtattatacctaggactgag
```

**Golden Gate (BsaI)** — an annealed oligo pair with 4-nt fusion overhangs,
ligated into a BsaI-linearised vector (defaults `TAGT`/`AAAC` = pEcgRNA
promoter/scaffold junctions; **set these to your vector's overhangs**):

- `<mutation>_gg_bsaI_top` — `TAGT` + protospacer
- `<mutation>_gg_bsaI_bot` — `AAAC` + revcomp(protospacer)

```
_gg_bsaI_top : TAGTGAACGTTGTGGATGTCCAGC
_gg_bsaI_bot : AAACGCTGGACATCCACAACGTTC
```

The Golden Gate oligos can be ordered 5′-phosphorylated (`gg_phosphorylate: true`)
for ligation. A spacer containing an internal **BsaI** site (`GGTCTC`/`GAGACC`)
cannot be cloned by Golden Gate — the QC tab flags this.

Oligos are deduplicated by sequence (first occurrence's name is kept).

### Oligo QC

The **QC tab** (and `oligo_qc_<date>.csv`) reports, per oligo:

| Metric | Notes |
|--------|-------|
| GC % | flags <30% / >65% |
| Tm (°C) | nearest-neighbor (SantaLucia); flags apply to **primers**, not the donor |
| Hairpin stem (bp) | longest self-complementary stem — heuristic, not ΔG |
| 3′ self-dimer (bp) | 3′-anchored self-complementarity (primers) |
| Max homopolymer | longest single-base run |
| BsaI / BbsI site | internal type-IIS site (breaks Golden Gate) |

Flagged oligos are highlighted. All thresholds live in `config/config.yml → qc`.
Values are dependency-free heuristics; for true thermodynamic ΔG, wire in
primer3/ViennaRNA (a hook is noted in `R/qc.R`).

---

## Configuration

Defaults live in [`config/config.yml`](config/config.yml) and every key is also
exposed in the app's **Settings** panel (GUI values override the file for that
session). Key parameters:

| Setting | Default | Meaning |
|---------|---------|---------|
| `design.oligo_size` | **100** | **Microsynth oligo size** — homology-arm (repair) oligo length in bp. Tunable. |
| `design.max_pam_distance` | 30 | Max PAM-to-codon distance (bp). |
| `design.min_score_threshold` | 3 | Minimum cassette score to keep a candidate. |
| `design.top_n` | 3 | Ranked designs generated per mutation. |
| `design.n_cores` | 1 | CPU cores (>1 = parallel). |
| `export.min_score_filter` | 3.0 | Min score when selecting top-ranked designs for the order form. |
| `microsynth.template_path` | `templates/MicrosynthUploadFormDNA__2_.xlsx` | Order-form template. |
| `microsynth.cloning_method` | `gibson` | gRNA cloning oligos: `gibson` / `golden_gate` / `both`. |
| `microsynth.fwd_tail` / `rev_tail` | `gttttagagctagaaatagcaagt` / `actagtattatacctaggactgag` | 24-nt 3′ Gibson homology tails appended to the gRNA primers. |
| `microsynth.gg_top_overhang` / `gg_bot_overhang` | `TAGT` / `AAAC` | Golden Gate BsaI 4-nt fusion overhangs (**vector-specific**). |
| `microsynth.gg_phosphorylate` | `false` | Order Golden Gate oligos 5′-phosphorylated. |
| `qc.na_mM` / `oligo_uM` | 50 / 0.5 | Salt & strand conc for the nearest-neighbor Tm. |
| `qc.gc_min`/`gc_max`, `tm_min`/`tm_max`, `hairpin_stem`, `dimer3`, `max_homopolymer` | — | QC flag thresholds. |
| `qc.enzymes` | `[BsaI, BbsI]` | Internal restriction sites checked per oligo. |
| `microsynth.pto_n_bases` | 4 | 5′ PTO-protected bases on the homology arm. |
| `microsynth.purification/scale/condition/datasheet/aliquots` | `DES` / `0.04` / `Dried` / `Standard` / `No` | Order-form column defaults. |
| `genome.default_strain` | `BW25113 …` | Strain selected by default in the app. |
| `genome.strains` | 6 strains | Selectable strains (label → NCBI accession); edit to add/correct. |

> **Behaviour note:** the original notebook hard-coded the homology arm at 85 bp;
> this app defaults it to **100 bp** and makes it configurable, as requested.
> Changing it also changes the length of the 200-bp cassette column accordingly.

---

## Project structure

```
.
├── app.R                 # Shiny GUI (thin layer over the backend)
├── run.R                 # Launcher
├── install.R             # One-time dependency installer
├── config/config.yml     # User-editable settings
├── R/                    # Backend pipeline (no GUI code)
│   ├── packages.R        #   dependency loading
│   ├── constants.R       #   biological constants + tunable defaults
│   ├── config.R          #   config load/merge/apply
│   ├── genome.R          #   genome/annotation download, caching, gene extraction
│   ├── offtarget.R       #   PAM index + fast off-target counting (Biostrings)
│   ├── design.R          #   gRNA design, homology arm, oligo assembly, library design
│   ├── microsynth.R      #   oligo order generation (Gibson + Golden Gate) + form filling
│   ├── qc.R              #   oligo QC: Tm, GC, hairpin, dimer, enzyme sites
│   └── pipeline.R        #   input reading, validation, orchestration, summaries
├── templates/            # Put the Microsynth .xlsx template here
├── examples/             # Example input files
└── legacy/               # Original .Rmd / output module (reference only)
```

The **backend** (`R/`) contains no Shiny code; the **GUI** (`app.R`) contains no
biology logic — so the pipeline can also be scripted directly (see below).

### Scripting the pipeline without the GUI

```r
for (f in c("packages.R","constants.R","config.R","genome.R","offtarget.R",
            "design.R","microsynth.R","pipeline.R")) source(file.path("R", f))
load_backend_packages()
apply_design_config(load_config("config/config.yml"))

system      <- initialize_system_fast()                          # default: BW25113
# or a different strain:   initialize_system_fast(strain_source("U00096.3", "MG1655"))
# or a custom genome:      initialize_system_fast(custom_source("mygenome.fasta", "mygenome.gff"))
results_all <- design_mutation_library_fast(c("dpsS100E","lrpS125D"), system, top_n = 3)
results_top <- select_top_ranked(results_all, min_score = 3.0)
generate_microsynth_output(results_top, load_config("config/config.yml"),
                           output_path = "MicrosynthOrder_filled.xlsx")
```

---

## Updating from GitHub

The app is version-controlled, so updating is:

```bash
git pull
Rscript install.R   # only needed if dependencies changed
Rscript run.R
```

Your `config/config.yml` edits and the cached genome are preserved across
updates. To publish this folder to GitHub the first time:

```bash
git init
git add .
git commit -m "Initial commit: CRISPR Oligo Designer"
git remote add origin https://github.com/<you>/<repo>.git
git push -u origin main
```

---

## Scope

This application implements the **point-mutation** design workflow. The
notebook's silent-mutation *control* generator
(`design_silent_control_library_fast`) is preserved in `legacy/` but is not yet
wired into the GUI; it can be added as a second tab or a companion app later.

## Troubleshooting

- **"Missing required package(s)"** — run `Rscript install.R`.
- **Genome download fails** — check your internet connection / proxy; the files
  come from NCBI E-utilities. Once cached, no network is needed.
- **"Microsynth template not found"** — put the `.xlsx` in `templates/` or set
  `microsynth.template_path` in `config/config.yml`.
- **"Gene not found"** — the gene symbol must match the selected genome's
  annotation (matching is case-insensitive). Different strains/accessions use
  different gene names; for a custom genome, gene lookup needs an uploaded GFF/GTF.
- **Custom genome has 0 genes** — you loaded a FASTA without an annotation. Upload
  a matching GFF/GTF so mutations can be mapped to genes.
