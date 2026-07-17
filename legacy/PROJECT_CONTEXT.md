# CRISPR gRNA & Oligo Generator — Project Context

## What this project is

An R script (`gRNA_oligo_generator.R`) that designs CRISPR guide RNAs and repair oligos for introducing point mutations in E. coli BW25113. It was written by Kian Bigovic Villi for a phosphosite mutagenesis library.

## What the main script does (DO NOT modify the core logic)

1. **Downloads and caches** the E. coli BW25113 genome (CP009273.1) and GFF annotations from NCBI.
2. **Builds a PAM index** across the whole genome using Biostrings for fast off-target analysis.
3. **For each desired amino acid mutation** (e.g. `dpsS100E`):
   - Finds nearby PAM sites (NGG) within 30 bp of the target codon.
   - Extracts the 20-nt protospacer for each candidate gRNA.
   - Identifies silent mutations to destroy the PAM after editing.
   - Scores candidates by PAM distance and off-target class.
   - Generates a 100-nt homology arm carrying the desired codon change + silent PAM mutation.
   - Assembles a full 200-bp oligo: `SPACER1 + homology_arm + SPACER2 + J23119_promoter + protospacer + sgRNA_handle`.
4. **Outputs a data frame** (`results_all`) with columns including `mutation`, `gene`, `protospacer`, `homology_arm`, `oligo_200bp`, `score`, `off_target_class`, `rank`, etc.
5. The user filters to `results_top <- results_all %>% filter(rank == 1)` to get the best design per mutation.

## What was changed (output section only)

A separate output module (`output_microsynth.R`) was created that takes `results_top` and:

1. **Generates 3 oligos per mutation:**
   - `<mutation>_gRNA_F`: `TAGT` + protospacer (24-nt). The `TAGT` is the 5' overhang for cloning into the guide RNA vector.
   - `<mutation>_gRNA_R`: `AAAC` + reverse complement of protospacer (24-nt). The `AAAC` is the 5' overhang on the complementary strand.
   - `<mutation>_HA`: The 100-nt homology arm with **phosphorothioate (PTO) bonds on the first 4 bases at the 5' end**, indicated as `X*X*X*Xrest...` (asterisks between the first 4 nucleotides).

2. **Deduplicates by sequence** — different mutations can share the same protospacer/gRNA pair; duplicates are removed (first occurrence kept), and the count of removed duplicates is reported.

3. **Fills a Microsynth order form** (`.xlsx` template) with the oligos in the "DNA Order" sheet, populating columns: Name (A), Sequence (B), Purification (D, default "DES"), Scale (E, default "0.04"), Physical Condition (L, default "Dried"), Datasheet (M, default "Standard"), Aliquots (N, default "No").

## Key design decisions

- PTO bonds (`*`) are ONLY on the homology arm oligos, NOT on the gRNA oligos.
- PTO is on the **first 4** (not 5) bases at the 5' end: `X*X*X*Xrest...` (3 asterisks).
- When the two gRNA oligos are annealed, they produce dsDNA with sticky ends compatible with BbsI/BpiI-digested guide RNA vectors:
  ```
  5'-TAGT<protospacer>     -3'
  3'-     <revcomp>CAAA-5'
  ```
- The Microsynth template file is `MicrosynthUploadFormDNA__2_.xlsx`. Data starts at row 2 in the "DNA Order" sheet.

## Files

| File | Description |
|------|-------------|
| `gRNA_oligo_generator.R` | Main R script — core pipeline. **Do not modify.** |
| `output_microsynth.R` | Output module — generates oligo table + fills Microsynth form. |
| `MicrosynthUploadFormDNA__2_.xlsx` | Microsynth order form template. |

## How to run

```r
# 1. Run the main pipeline (produces results_top)
source("gRNA_oligo_generator.R")
system <- initialize_system_fast()
results_all <- design_mutation_library_fast(mutations, system, top_n = 3)
results_top <- results_all %>% filter(rank == 1)

# 2. Generate oligos and fill order form
source("output_microsynth.R")
oligo_table <- generate_oligo_orders(results_top)
fill_microsynth_form(
  oligo_table,
  template_path = "MicrosynthUploadFormDNA__2_.xlsx",
  output_path   = "MicrosynthOrder_filled.xlsx"
)
```
