# Pre-specification addendum 5: validation re-run, batch check, positive/negative controls, H0-only and H1-only tests

**Date:** 2026-10-06
**Status:** written after all step 5 results, addenda 1–4 and deviations D1–D8 (technical report v1, commit `713e6ca`), and **before any of the analyses below was run**. Parts (1)–(4) are run in this order, with a stop after each.
**Earlier documents** are unchanged; **the pre-specified tipping-point decision is unchanged.** All parts are quality checks or secondary analyses.
**Rule:** this file is not changed after part (1) starts. Deviations go into `docs/analysis_deviations.md`.

---

## Part 1: re-run of the reference-code validation

- **What:** run `reference/validate_tda_pipeline.py` (full: Part A deterministic checks A1–A6; Part B calibration B1–B4, 100 replicates) unchanged, in the project environment, from `results/validation/`.
- **Commit:** `validation_report.txt` and the console log.
- **Acceptance:** 10/10 checks pass, as on 2026-10-03 (B1–B3 type I error ≤ 0.05 + 2.5 × MC SE; B4 power ≥ 0.8).

## Part 2: batch check, GSE41258 scan dates vs stage

- **Data:** GEO characteristic "scan date" for the samples of set `main` (53 / 48 / 181 / 39) and, for reference, `main_allarrays`.
- **Parsing:** values of the form dd/mm/yyyy are read as dates; plain integers are read as Excel serial dates (days since 1899-12-30). Unparseable or missing values are counted and reported.
- **Descriptive:** scan year and scan day per stage (cross-tables); median and range of scan date per stage.
- **Measures, fixed in advance:**
  1. **Cramér's V** between scan year and stage (4 stages). Interpretation: < 0.2 weak, 0.2 to < 0.5 moderate, ≥ 0.5 strong association.
  2. For each adjacent transition (N→A, A→C, C→M): the **AUC** of scan date for separating the two stages (Mann–Whitney), plus the share of the two stages' samples that were scanned on days when both stages were scanned ("day overlap").
- **Confounding criterion:** a transition is called **confounded with scan date** if AUC ≥ 0.80 or ≤ 0.20 **or** day overlap < 20%. Cramér's V ≥ 0.5 is reported as strong overall association.
- Descriptive only; no change to any result. If confounding is found, its handling is decided with the PI and recorded as a deviation.

## Part 3: positive- and negative-control tests

**Purpose:** check on **real** preprocessed GSE41258 data that the step 5 test:
- keeps its false-positive rate;
- detects inserted changes of known type;
- in shape-only mode, **ignores pure changes of spread**.

**Base data:** the 181 carcinomas of set `main` (one stage, so no stage biology). In every replicate the carcinomas are split at random, without replacement, into **group A (53)** and **group B (48)**, the GSE41258 Normal/Adenoma sizes. Insertions are made into **group B only**.

**Test:** step 5 test unchanged (`tda_crc_pilot.permutation_test`; VR H0+H1; k = 20; m = 31; aligned) with **99 permutations**. Each condition is run in **shape-only and raw** mode, **50 replicates** each. Seeds: `SeedSequence(42, crc32(unit name))`. Reported: detection rate at p ≤ 0.05 and at p ≤ 0.0167 (Bonferroni, 3 transitions), and median z.

### Level 1: embedding level (insertion into the 10-PC scores; isolates steps 4–5)
Embedding: GSE41258 `main` unadjusted, 2,000 genes, 10 PCs (step 3 file); rows of group B are modified. c_B is B's centroid; MPD_B is B's median pairwise distance; r_B is the median distance of B's points to c_B (all computed before insertion).

| Condition | Insertion | Levels |
|---|---|---|
| N0 | none (negative control) | – |
| SP | spread: x ← c_B + s (x − c_B) | s ∈ {1.05, 1.10, 1.25, 1.50} |
| CL | cluster split: random unit direction u; a random half of B gets + δ u, the other half − δ u | δ ∈ {0.25, 0.5, 1.0} × MPD_B (cluster separation 0.5 / 1 / 2 × MPD_B) |
| LP | loop: random 2-D plane (orthonormal u, v). The plane coordinates of B are replaced by w × annulus + (1 − w) × original, where the annulus has radius R = 1.5 × r_B, angles 2π·(random order)/n_B, and radial noise SD 0.1 R. The other 8 coordinates are unchanged. | w ∈ {0.25, 0.5, 1.0} |

### Level 2: expression level (full pipeline, steps 3–5)
- **Expression:** GSE41258 step 1 log2 expression of the 101 samples in A ∪ B.
- **Genes:** 200 random genes among those passing the 25th-percentile mean filter in A ∪ B. SD_g = each gene's SD in B.
- **Insertions** (group B, those 200 genes only; L1, L2 = independent N(0, 1) gene loadings; a ∈ {0.25, 0.5, 1.0}):
  - *spread:* deviations from B's gene means multiplied by s = 1 + a;
  - *cluster:* ± a · SD_g · L1_g for the two random halves of B;
  - *loop:* a · SD_g · (cos θ_i · L1_g + sin θ_i · L2_g), θ_i = 2π·(random order)/n_B.
- **Pipeline:** `tda_crc_pilot.preprocess` (top 2,000 variable genes, z-score, 10 PCs) is then fitted on A ∪ B pooled, followed by the test as above.
- **Conditions:** N0 + 3 types × 3 amplitudes.

### Acceptance criteria (fixed in advance; 50 replicates)

| # | Criterion | Rule (detection at p ≤ 0.05) |
|---|---|---|
| C1 | Negative control calibrated (level 1 and level 2, both modes) | ≤ 7/50 (14%) |
| C2 | Raw mode detects spread | level 1, s = 1.25: ≥ 40/50 (80%) |
| C3 | **Shape-only ignores pure spread** | level 1, every s: ≤ 7/50 (14%) |
| C4 | Shape-only detects the largest cluster split | level 1, δ = 1.0 × MPD_B: ≥ 40/50 |
| C5 | Shape-only detects the full loop | level 1, w = 1.0: ≥ 40/50 |
| C6 | The full pipeline preserves inserted structure | level 2, a = 1.0 cluster and loop, shape-only: ≥ 40/50 |

The 14% limit gives about a 1% chance of failing C1 or C3 when the true rate is 5% (binomial, n = 50).

**Interpretation of failures:**
- C1 fails → miscalibration on real data.
- C3 fails → scale leaks into the "shape-only" statistic.
- C2, C4 or C5 fail → low sensitivity of the topological test itself.
- C4/C5 pass but C6 fails → information is lost in preprocessing (gene selection, z-scoring, PCA truncation).

Level 2 spread is not required to be ignored by shape-only, because a spread change in 200 genes changes the covariance shape after z-scoring and PCA. It is reported only.

## Part 4: H0-only and H1-only versions of the step 5 test (secondary)

- **What:** the step 5 primary test repeated with dims = (0,) (**H0 only**) and dims = (1,) (**H1 only**) instead of (0, 1). Everything else is unchanged: 999 permutations, k = 20, m as in step 5, aligned, primary embeddings.
- **Scope:** all three cohorts (GSE41258 `main`: 3 transitions; GSE20916 `bulk` and GSE37364 `bulk`: 2 transitions), **shape-only and raw**, **unadjusted and adjusted** (primary, stage-blind).
- **Multiplicity:** Holm across the transitions within each cohort × version × mode × dimension.
- **Role:** **secondary.** It shows which homology dimension drives the step 5 results (technical report, systematic-error row 14). It does not change the tipping-point definition (which uses H0 + H1).
- **Seeds:** `SeedSequence(42, crc32(unit name))`; resumable checkpoints.

---

## Part 2b (appended 2026-10-06, after the part 2 results and before running part 2b; parts 1, 2, 3 and 4 above are unchanged)

**Why:** part 2 (with the D9 date correction) showed that Normal and Adenoma arrays of GSE41258 were scanned on different days (1 shared day, 3% overlap). Part 2b measures how large a topological difference a pure scan-batch contrast produces in this dataset.

**Batch-only comparison (same tissue, different scan periods):** GSE41258 `main` normals with corrected scan date (D9) in **spring 2003** (2003-01-01 to 2003-06-30; n = 21) vs normals scanned in **late 2005–2006** (2005-09-01 to 2006-12-31; n = 30). The two normals scanned outside these windows (2003-11-30, 2004-04-02) are not used.

**Test:** frozen step 5 test (`tda_crc_pilot.permutation_test`; VR H0+H1; k = 20; aligned; **999 permutations**), with the **primary step 3 embeddings** of set `main` (unadjusted and adjusted; the same coordinates as the step 5 Normal → Adenoma test), **raw and shape-only**. The groups are smaller than the frozen subsample size (31), so the subsample size follows the pilot pair rule: m = min(31, ⌊0.8 × 21⌋) = **16**. Four tests; **Holm across the four.** Seeds: `SeedSequence(42, crc32(unit name))`.

**Size-matched Normal → Adenoma reference:** z-scores depend on group size and m. The step 5 Normal → Adenoma test (53 vs 48, m = 31) is therefore also repeated at the batch comparison's sizes: **20 random draws** of 21 normals and 30 adenomas (without replacement), m = 16, 999 permutations, same four mode × version combinations. Median z and the range over the 20 draws are reported next to the step 5 z and the batch z.

**Decision, stated in advance:** if any of the four batch-only tests has **Holm p < 0.05**, batch can mimic a stage difference in this dataset, and the GSE41258 Normal → Adenoma results must be treated as possibly batch-driven. If none is significant, the batch contrast is reported with its z-scores, which are compared with the size-matched Normal → Adenoma z-scores to judge how much of that difference batch could plausibly explain. This is descriptive; it changes no step 5 result.

---

## Part 3b (appended 2026-10-06, after the part 3 results and before running part 3b; earlier parts unchanged)

**Why:** level 1 controls passed (C1–C5) but level 2 failed (C6): no insertion into 200 genes was detected after the full pipeline. Part 3b checks whether the insertion really reached the analysis (bug check) and, if it did, at which preprocessing step it is lost. **Descriptive only: no permutation tests, no change to the frozen pipeline.**

**Data:** exact reconstructions of the level 2 data that were tested in part 3. For each insertion, the **first replicate of each of the 5 saved chunks** (5 repeats), rebuilt from the same per-unit seeds (`SeedSequence(42, crc32(unit name))`, shape-only units), so they are byte-identical to the tested data. Insertions: cluster split a = 1, loop a = 1, spread a = 0.5 (+50%) and spread a = 1 (+100%, the largest).

**Measures per repeat:**
- (a) **Sanity / bug check:** for the 200 altered genes, on the exact matrix passed to `tda_crc_pilot.preprocess` (log2, the 101 samples A ∪ B): mean B − A difference and SD ratio B/A **before and after** insertion (medians over genes); the number of changed matrix entries; and confirmation that only the 200 genes × 48 B samples changed.
- (b) Number of the 200 genes in the top-2,000 variance selection, after insertion and (for reference) before.
- (c) Inserted signal in z-space = z-scored(after) − z-scored(before, using the after-insertion gene means and SDs), restricted to the selected genes. Reported: the fraction of its squared norm lying in the first **10** principal components of the after-insertion data (exact SVD); the PC carrying the largest share, with that PC's share of total variance; the inserted signal's share of total z-space variance; and the fraction of the 200 genes' inserted energy kept by gene selection.
- (d) As (c) for the first **50** components.

**Outcome:** the table is used to decide between (i) the insertion did not reach the analysis (bug) and (ii) a real preprocessing loss, located by (b)–(d).
