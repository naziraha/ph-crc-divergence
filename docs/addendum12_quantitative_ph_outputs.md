# Addendum 12: quantitative persistence outputs (landscapes, barcode summaries, dimension split) and standard PH figures

**Date:** 2026-10-10 (locked at commit, before any computation)
**Status:** written before any computation in this addendum. All step 5 results (four GEO cohorts and TCGA) are known, so this is an **informed, secondary** analysis.
**Relation to earlier documents:**
- **Unchanged:** the main pre-specification and addenda 1–11.
- **The tipping-point decision (step 5, shape-only Wasserstein test) remains the only primary test.** Nothing in this addendum can create, remove or override a tipping point.

**Rule:** this file is not changed after the computation starts. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Purpose
1. **Standard PH outputs:** provide the persistent-homology outputs expected in a TDA paper (persistence diagrams, barcodes, persistence landscapes, separation by homology dimension) for every cohort.
2. **Patient-level tests:** make them **quantitative** with tests at the patient level, so they support the primary result with numbers rather than pictures alone.

## 2. Data (frozen)
- **Cohorts and sets:**
  - GSE41258 `main` (N 53 / A 48 / C 181 / M 39);
  - GSE20916 `bulk` (24 / 45 / 36);
  - GSE37364 `bulk` (38 / 29 / 27);
  - GSE117606 `main` (65 / 59 / 71);
  - TCGA-COAD `main` (stages I–IV, 45 / 110 / 80 / 39).
- **Embeddings:** frozen 10-PC embeddings, unadjusted and stage-blind adjusted.
- **Transitions (as step 5):**
  - GEO: N→A and A→C, plus C→M in GSE41258;
  - TCGA: I→II, II→III, III→IV.
- **Subsampling and diagrams (as step 5):**
  - k = 20 subsamples of m (80% of the smallest stage) per stage;
  - Vietoris–Rips with GUDHI (`tda_crc_pilot.persistence`), H0 and H1;
  - shape-only (distances / subsample median) is **primary**, raw is secondary;
  - groups are aligned before permuting (`tda_crc_pilot.align_groups`), exactly as in step 5.
- **Sample counts** are checked against step 5 before computing.

## 3. Quantitative analyses

### 3.1 Persistence landscapes (Bubenik)
- **Landscape per subsample diagram:** `gudhi.representations.Landscape`, first **5 landscapes**, **100 grid points**. Computed separately for **H0** and **H1**; infinite bars removed (as in step 4).
- **Fixed grid range:** set once per cohort × version × mode × dimension, from the pooled observed diagrams of all stages (0 to the maximum death value), and kept fixed for all permutations.
- **Statistic (per transition s → t):** L² distance between the mean landscape of stage s and that of stage t, each averaged over k = 20 subsamples.
- **Patient-level permutation test:**
  - the aligned stage labels are permuted across the patients of the two stages (group sizes kept), and k subsamples are redrawn per permuted group;
  - **999 permutations**; one-sided p = (1 + #{null ≥ observed}) / (1 + 999);
  - z relative to the null.
  - This mirrors `tda_crc_pilot.permutation_test` with the Wasserstein statistic replaced by the landscape distance.
- **Multiplicity:** Holm across the transitions within each cohort × version × mode × dimension.

### 3.2 Barcode summary statistics
- **Per subsample barcode:**
  - **H0:** total persistence, maximum bar length, persistence entropy. The number of H0 bars is fixed at m − 1, so it is not tested.
  - **H1:** total persistence, maximum bar length, number of bars, persistence entropy.
- **Statistic:** the absolute difference between the two stages' means over k = 20 subsamples.
- **Test:** the same patient-level permutation scheme (999 permutations) as §3.1.
- **Multiplicity:** Holm across transitions within each cohort × version × mode × dimension × summary.
  - Because there are 7 summaries, a further **Holm across the 7 summaries** is applied within each cohort × version × mode × transition.
  - Only results significant after **both** corrections are called significant.

### 3.3 Separation by dimension: complete coverage
- **GSE117606:** the H0-only and H1-only versions of the step 5 Wasserstein test are run here (as addendum 5 part 4: unadjusted and adjusted, shape-only and raw), since it was the only cohort without them.
- **Combined table:** all cohorts' dimension results (addendum 5 part 4, addendum 10 §4.1, and this run) are reported together.

## 4. Figures (descriptive, from the same diagrams)
For each cohort × version (shape-only primary; raw in the supplement):
1. **Persistence diagrams:** H0 and H1 for one representative subsample per stage (the first seeded draw), plus an overlay of all k subsamples per stage.
2. **Persistence barcodes:** H0 and H1 for the same representative subsample per stage.
3. **Mean persistence landscapes:** first landscape (λ₁) for H0 and H1 per stage, with the 2.5–97.5 percentile band across the k subsamples. The band is labelled as describing subsample variability, **not** a confidence interval between stages.
4. **Betti-0 and Betti-1 curves** per stage (mean with IQR band). Step 4 showed Betti-1 only.

## 5. Interpretation (fixed in advance)
1. **Secondary status:** these analyses are **secondary**. They **support or qualify** the primary step 5 result and cannot create or remove a tipping point.
2. **Concordance with step 5 (primary summary):**
   - For each cohort × version × transition, the shape-only landscape test (H0 and H1 separately) is compared with the step 5 shape-only call (significant / not significant, Holm p < 0.05).
   - **Supportive:** the landscape tests agree with step 5 in **≥ 80%** of comparisons.
   - Otherwise, the discordant comparisons are listed and described without changing the primary conclusion.
3. **Dimension:** "predominantly H0" is retained if, across all cohorts, H0-only results reproduce at least 80% of the significant combined step 5 results and H1-only results are never significant where the combined test is not. Otherwise, the description is revised to match the data.
4. **Barcode summaries:** reported as quantitative descriptions of which aspect differs (for example, total H1 persistence higher in stage X). A significant summary is interpreted only together with the step 5 result for the same comparison.
5. **Wording:** never "loops drive progression" unless H1 results are significant in both data versions in ≥ 3 of 4 GEO cohorts. That is not expected, and the rule is stated in advance.

## 6. Validation (before any real data)
The landscape and barcode-summary permutation tests are new, so they are checked first on real data with known insertions, using the addendum 5 level 1 design:
- **Setup:**
  - GSE41258 carcinomas split at random into A 53 and B 48;
  - insertions into B (`addendum5_part3_controls.insert_level1`);
  - k = 20, m = 31, aligned;
  - 99 permutations; 20 replicates per condition.
- **Conditions:** N0 none; SP spread s = 1.25; CL cluster split δ = 1 × median pairwise distance; LP loop w = 1.

**Acceptance:**

| Check | Requirement |
|---|---|
| N0, landscape test, H0 and H1, shape-only and raw | ≤ 3/20 detected |
| SP, landscape test, shape-only (should ignore pure spread) | ≤ 3/20 |
| SP, landscape test, raw | ≥ 16/20 |
| CL, landscape test, H0, shape-only | ≥ 16/20 |
| LP, landscape test, H1, shape-only | ≥ 16/20 |
| N0, barcode summaries (each summary, shape-only) | ≤ 3/20 |

**If any criterion fails:** no real data are analysed; the failure is reported and recorded as a deviation, and the PI decides.

## 7. Settings, outputs, run time
- **Seeds:** `SeedSequence(42, crc32(unit name))` per unit.
- **Packages:** GUDHI 3.13.0 (`gudhi.representations`), numpy, scipy. No new installs.
- **Outputs (`results/addendum12/`):**
  - `validation_summary.txt`;
  - `landscape_tests.csv`;
  - `barcode_summary_tests.csv`;
  - `dimension_split_all_cohorts.csv`;
  - `concordance_with_step5.csv`;
  - `summary.txt`;
  - `figures/<cohort>_<version>_{diagrams,barcodes,landscapes,betti}.png`.
- **Expected run time:** several hours in the background (permutation tests over about 5 cohorts × 2 versions × 2 modes × 2 dimensions × 2–3 transitions, plus validation). Checkpointed per unit.
