# Scoring model

How the app ranks and filters candidate gRNAs. The score is a single number per
candidate used to (a) **discard** weak candidates during design and (b) **rank**
the survivors so the best one per mutation becomes `rank == 1`.

Implemented in [`R/design.R`](R/design.R) (`calculate_cassette_score`); all
parameters live in [`config/config.yml`](config/config.yml) under
`design.scoring` and three toggles are in the app's **Settings → Scoring** panel.

> **The defaults reproduce the original model exactly.** Every optional term
> below is off unless you turn it on, so nothing changes until you opt in.

---

## The formula

```
score = base[class]
        − distance_penalty
        − offtarget_penalty      (optional, default off)
        − activity_penalty       (optional, default off)
score = max(score, 0)
```

A candidate is kept only if `score ≥ design.min_score_threshold` (default **3**).
Survivors are sorted by score; the top `top_n` per mutation are returned, and the
best (`rank == 1`) is what the order form uses.

---

## 1. Off-target class → `base[class]`

Each candidate is assigned one of four classes from two genome-wide checks
([`R/offtarget.R`](R/offtarget.R)):

- **seed match** — another site shares the exact **PAM-proximal 11-nt seed**
  (positions 10–20), regardless of the rest.
- **off-target** — another site matches the **full 20-nt** spacer with **1–4
  mismatches**.

| Class | Seed match? | ≤4-mismatch off-target? | Meaning | `base` (default) |
|------:|:-----------:|:-----------------------:|---------|-----------------:|
| 1 | no | no | cleanest | **5.00** |
| 2 | no | yes | near-match, but mismatched in the seed | **4.41** |
| 3 | yes | no | exact seed match elsewhere, distal differences | **3.90** |
| 4 | yes | yes | both | **2.00** |

**Note on class 2 vs 3.** By default class 2 outranks class 3. Because the seed
drives Cas9 recognition, an *exact seed match* (class 3) is arguably at least as
risky as a seed-*mismatched* near-match (class 2). Set `swap_class_2_3: true`
(or tick the GUI box) to rank class 3 below class 2 instead.

**Note on magnitude.** The class is binary in each check — a guide with **1**
off-target and one with **50** get the same class. The optional off-target
penalty (below) restores that magnitude.

---

## 2. Distance penalty

```
distance_penalty = distance_penalty_per_step × floor( min(distance, distance_cap_bp) / distance_step_bp )
```

Default: **−0.2 per 3 bp** of PAM-to-edit distance, capped at **30 bp**. Editing
efficiency falls as the cut moves away from the edit, so closer PAMs score higher.

### Emergent class ↔ distance budget

Because of the `≥ 3` keep-threshold, the base scores + distance penalty impose an
implicit maximum distance per class (with all defaults):

| Class | base | Max distance that still passes (`score ≥ 3`) |
|------:|-----:|:--------------------------------------------:|
| 1 | 5.00 | ≤ 30 bp (all) |
| 2 | 4.41 | ≤ 23 bp |
| 3 | 3.90 | **≤ 14 bp** |
| 4 | 2.00 | **never** (2.0 < 3) |

So class 4 is always rejected, and each class has a distance budget. This is a
*consequence* of the numbers, not a separate rule — worth knowing if you retune
`base_scores` or `min_score_threshold`.

---

## 3. Off-target magnitude penalty *(optional, default off)*

```
offtarget_penalty = offtarget_count_weight × ln( 1 + off_target_count + seed_matches )
```

Config: `offtarget_count_weight` (default **0.0** = off). Uses the counts the
pipeline already computes, on a log scale so the first few off-targets matter
most. Example at weight 0.5: 10 off-targets → −0.5·ln(11) ≈ **−1.20**.

Turning this on makes the score discriminate *within* a class (few vs many
off-targets), which the class alone cannot.

---

## 4. On-target activity penalty *(optional, default off)*

```
activity_penalty = (GC outside [gc_ideal_min, gc_ideal_max] ? gc_penalty : 0)
                 + (longest homopolymer run > homopolymer_max ? homopolymer_penalty : 0)
```

Config: `activity_term: true` to enable. Defaults: ideal GC **0.40–0.65**
(penalty **0.5**), homopolymer run > **4** (penalty **0.5**). Same GC/homopolymer
logic as the QC tab. Flags spacers likely to express or cut poorly (GC extremes,
poly-N tracts).

This is deliberately coarse. It is **not** a trained efficiency model (Doench /
Azimuth); it is a cheap penalty for obviously bad spacers. For a Pol-III-free
*E. coli* J23119 system the exact efficiency determinants differ from mammalian
Cas9, so treat it as a guardrail, not a predictor.

---

## Worked examples (defaults)

| class | distance | off-targets | score | note |
|:-----:|:--------:|:-----------:|:-----:|------|
| 1 | 0 bp | — | 5.00 | perfect |
| 1 | 30 bp | — | 3.00 | clean but far — just passes |
| 2 | 6 bp | — | 4.01 | |
| 3 | 15 bp | — | **filtered** | 3.9 − 1.0 = 2.9 < 3 |
| 4 | 0 bp | — | **filtered** | 2.0 < 3 |

With `offtarget_count_weight = 0.5`, a class-2 guide at 0 bp with 10 off-targets
drops 4.41 → **3.21**. With `activity_term = true`, a class-1 guide at 0 bp with
80 % GC drops 5.00 → **4.50**.

---

## Two different thresholds (don't confuse them)

| Setting | When | Effect |
|---------|------|--------|
| `design.min_score_threshold` (3) | during design | candidates below are **discarded** (never ranked) |
| `export.min_score_filter` (3.0) | in the app, picking order candidates | filters the already-ranked `rank == 1` set for the order form |

If you enable a penalty term, scores shift down — you may want to **lower
`min_score_threshold`** so good candidates are not over-filtered.

---

## Configuration

`config/config.yml`:

```yaml
design:
  min_score_threshold: 3
  scoring:
    base_scores: [5.0, 4.41, 3.9, 2.0]
    swap_class_2_3: false
    distance_step_bp: 3
    distance_penalty_per_step: 0.2
    distance_cap_bp: 30
    offtarget_count_weight: 0.0     # > 0 to enable the magnitude penalty
    activity_term: false            # true to enable the GC/homopolymer penalty
    gc_ideal_min: 0.40
    gc_ideal_max: 0.65
    gc_penalty: 0.5
    homopolymer_max: 4
    homopolymer_penalty: 0.5
```

The three most useful toggles (`swap_class_2_3`, `offtarget_count_weight`,
`activity_term`) are also in **Settings → Scoring** for quick A/B runs without
editing the file. A partial `scoring:` block is merged onto the defaults, so you
only need to list the keys you change.

---

## What the score does *not* capture

- A trained on-target efficiency model (sequence-context, nucleosome, etc.).
- Mismatch **position** within an off-target (all 1–4-mismatch hits treated alike).
- PAM strength (NGG vs NAG) or alternative PAMs.
- sgRNA secondary structure beyond the fixed scaffold.
- The disruptiveness of the silent PAM-blocking mutation (`pam_mutation_freq_diff`
  is computed and used to *pick* the silent change, but does not enter the score).

It is a transparent, tunable **specificity-first filter and ranker** — strong for
rejecting risky or distal guides, intentionally light on predicting raw cutting
efficiency.
