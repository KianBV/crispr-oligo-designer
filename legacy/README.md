# Legacy source

These are the **original** files the application was refactored from. They are
kept for reference and provenance only — the running application does **not**
use them.

| File | What it was |
|------|-------------|
| `26032026_oligost_debugged.Rmd` | Original R Markdown notebook with the full pipeline (genome loading, off-target, gRNA design, homology arm, oligo assembly, silent controls). Split into `R/*.R`. |
| `output_microsynth.R` | Original Microsynth output module. Refactored into `R/microsynth.R`. |
| `PROJECT_CONTEXT.md` | Background notes on the project and design decisions. |
| `INSTRUCTIONS_FOR_CLAUDE.md` | Notes on the output-module constraints. |

The silent-control pipeline from the notebook
(`design_silent_control_library_fast`) was intentionally left out of the GUI
app for now (see the project README, "Scope").
