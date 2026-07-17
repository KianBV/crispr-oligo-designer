# Instructions for Claude

You are working on a CRISPR gRNA and oligo design pipeline for E. coli mutagenesis. Read `PROJECT_CONTEXT.md` for full background.

## Critical constraints

1. **DO NOT modify the core R script** (`gRNA_oligo_generator.R`) unless the user explicitly asks you to. The core logic (PAM finding, off-target analysis, scoring, homology arm generation, 200-bp oligo assembly) is finalized and tested.

2. **Only modify the output module** (`output_microsynth.R`) when asked. This module is responsible for:
   - Generating the 3 oligos per mutation (gRNA_F, gRNA_R, HA)
   - Adding overhangs (TAGT to forward protospacer, AAAC to reverse complement)
   - Adding PTO bonds to the homology arm only (first 4 bases: `X*X*X*Xrest...`)
   - Deduplicating by sequence
   - Filling the Microsynth Excel order form

3. **PTO bonds** are asterisks (`*`) between bases, placed ONLY on the homology arm, ONLY on the first 4 bases at the 5' end. Not on gRNA oligos.

4. **Overhang scheme** for gRNA oligos (BbsI/BpiI cloning):
   - Forward: `5'-TAGT` + protospacer
   - Reverse: `5'-AAAC` + reverse complement of protospacer

5. The **Microsynth template** is `MicrosynthUploadFormDNA__2_.xlsx`. The "DNA Order" sheet has these validated columns:
   - A: Name (max 35 chars)
   - B: Sequence
   - C: Length (formula, leave empty)
   - D: Purification (dropdown, default "DES")
   - E: Scale (dropdown, default "0.04")
   - F: 5' modification (leave empty)
   - G-J: Internal modifications (leave empty)
   - K: 3' modification (leave empty)
   - L: Physical Condition (dropdown, default "Dried")
   - M: Datasheet (dropdown, default "Standard")
   - N: Aliquots (dropdown, default "No")

6. **No duplicate oligos** — deduplication is by sequence; first occurrence name is kept.

## Key data flow

```
User's mutation list (e.g. dpsS100E, dpsS100D, ...)
  → initialize_system_fast()  [loads genome, builds PAM index]
  → design_mutation_library_fast()  [returns results_all with protospacer, homology_arm, etc.]
  → filter(rank == 1)  [results_top: best design per mutation]
  → generate_oligo_orders(results_top)  [builds 3 oligos per mutation, deduplicates]
  → fill_microsynth_form()  [writes into Microsynth .xlsx template]
```
