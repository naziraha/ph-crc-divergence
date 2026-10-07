# Addendum 7: Step 13 benchmark against conventional multivariate methods (PERMDISP, PERMANOVA, stage separation, minimum spanning tree)

**Date:** 2026-10-07 (locked at commit, before any benchmark computation)
**Status:** written after all step 5 TDA results (four cohorts) are known, and **before any benchmark method has been run on any cohort**. Because the TDA results are known, this is an **informed** comparison; its rules are fixed here so that they cannot be chosen after seeing the conventional results.
**Relation to earlier documents:** the main pre-specification and addenda 1–6 are unchanged. **The tipping-point decision and all TDA results are unchanged**; this addendum only adds conventional comparisons.
**Rule:** this file is not changed after the benchmark computation starts. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Purpose

Plan step 13 asks what topological analysis detects beyond conventional methods. The step 5 TDA results showed two things:
- **Adenoma → Carcinoma:** a difference in spread (TDA raw-distance test, Holm p < 0.05 in both data versions in GSE20916, GSE37364 and GSE117606; unadjusted only in GSE41258).
- **Normal → Adenoma:** a shape difference (TDA shape-only test) that is inconsistent across cohorts and data versions, and is carried by H0 (clustering structure), not loops.

This addendum asks whether conventional methods find the same differences on **exactly the same data**:
- **PERMDISP:** dispersion, the conventional counterpart of "spread".
- **PERMANOVA:** location (centroid difference). The TDA test removes location by aligning groups, so PERMANOVA shows what TDA deliberately ignores.
- **Stage separation (Ward clustering + adjusted Rand index):** whether an unsupervised clustering splits the samples along stage lines. This is a **location-type** measure: it responds mainly to the two stages lying in different places. It is **not** a counterpart of H0.
- **Minimum spanning tree (MST) edge-length test:** the conventional counterpart of H0. For a Vietoris–Rips filtration, the H0 death times are **exactly** the edge lengths of the Euclidean minimum spanning tree. Agreement between this test and the TDA test therefore means that **the H0 signal is reproducible without persistent homology**.

Differential expression and WGCNA (the other parts of plan step 13) are **not** part of this addendum.

## 2. Data (identical to step 5)

| Cohort | Set (n per stage) | Embedding |
|---|---|---|
| GSE41258 | `main` (N 53 / A 48 / C 181 / M 39) | `results/step03/embeddings/GSE41258_main_<version>_g2000_pc10.csv` |
| GSE20916 | `bulk` (24 / 45 / 36) | `GSE20916_bulk_<version>_g2000_pc10.csv` |
| GSE37364 | `bulk` (38 / 29 / 27) | `GSE37364_bulk_<version>_g2000_pc10.csv` |
| GSE117606 | `main` (65 / 59 / 71) | `GSE117606_main_<version>_g2000_pc10.csv` |

- **Data versions:** `unadjusted` and `adjusted` (stage-blind; primary). The same frozen **10-PC coordinates** used by the TDA test, with Euclidean distance. Stage-aware data and sensitivity embeddings (5/20 PCs; 1,000/5,000 genes) are not included.
- **Transitions:** the same as step 5:
  - GSE41258: Normal → Adenoma, Adenoma → Carcinoma, Carcinoma → Metastasis.
  - Other cohorts: Normal → Adenoma, Adenoma → Carcinoma.
- **Comparisons:** 18 per method and mode (GSE41258: 3 transitions × 2 versions; three other cohorts: 2 × 2 each).
- **Samples:**
  - PERMDISP, PERMANOVA and stage separation (primary analysis) use **all** samples of the two stages, because these methods do not need equal group sizes.
  - The MST test and the size-matched analysis (section 3.6) use the step 5 subsampling: k = 20 subsamples of m = 80% of the smallest stage, per stage.
- **Sample counts** are checked against step 5 before computing; any mismatch stops the run.

## 3. Methods (fixed in advance)

All methods are implemented in numpy / scipy (no new packages). Each is first checked on simulated data (section 4).

### 3.1 PERMDISP (dispersion; Anderson 2006)
- For each sample: its Euclidean distance to the **centroid of its own stage** (10-PC space).
- **Statistic:** one-way ANOVA F of these distances between the two stages.
- **Permutation:** the stage labels of the distances are permuted 999 times (group sizes kept).
- **p-value:** one-sided, p = (1 + #{F_perm ≥ F_obs}) / (1 + 999). Large F means a dispersion difference in either direction.
- **Direction:** dispersion ratio D = mean distance (later stage) / mean distance (earlier stage). D > 1 means the later stage is more spread out.

### 3.2 PERMANOVA (location; Anderson 2001)
- One-way PERMANOVA on the Euclidean distance matrix of the two stages.
- **Statistic:** pseudo-F = SS_between / (SS_within / (n − 2)).
- **Permutation:** labels permuted 999 times, group sizes kept.
- **Effect size:** R² = SS_between / SS_total.
- **Known limitation:** with unequal group sizes, PERMANOVA also responds to dispersion differences, so it is always read together with PERMDISP.

### 3.3 Stage separation (Ward clustering + adjusted Rand index; location-type)
- The pooled samples of the two stages are clustered with **Ward linkage** (Euclidean) and the tree is cut into **2 clusters**.
- **Statistic:** adjusted Rand index (ARI) between the 2 clusters and the stage labels.
- **Permutation:** the clustering is fixed and the stage labels are permuted 999 times; one-sided p = (1 + #{ARI_perm ≥ ARI_obs}) / (1 + 999).
- **Reading:** ARI near 0 means the stages are not separated by unsupervised clustering; ARI near 1 means they separate cleanly. It measures separation between the stages, not clustering structure **within** a stage.

### 3.4 Minimum spanning tree edge-length test (conventional counterpart of H0)
- **Basis:** for a finite point set, the finite H0 death times of the Vietoris–Rips filtration are exactly the m − 1 edge lengths of the Euclidean minimum spanning tree. So H0 can be computed without persistent homology.
- **Per subsample:** compute the MST of the m points (scipy.sparse.csgraph.minimum_spanning_tree on the Euclidean distance matrix), keeping its m − 1 edge lengths.
- **Two modes:**
  - **raw**;
  - **shape-only analogue:** edge lengths divided by the subsample's median pairwise distance, as in step 5.
- **Distance between two subsamples:** the 1-Wasserstein distance between their two edge-length distributions (scipy.stats.wasserstein_distance; 1-D).
  - This is not identical to the TDA test's diagram distance. The diagram Wasserstein distance (L∞ ground metric) can also match points to the diagonal. The two are therefore expected to agree closely, but not exactly.
- **Subsampling, statistic, alignment and permutations:** identical to step 5 (`tda_crc_pilot.permutation_test`), with the persistence-diagram step replaced by MST edge lengths:
  - k = 20 subsamples of size m per stage, without replacement;
  - statistic = mean between-stage − mean within-stage distance;
  - groups centred before permuting (and also scaled in the shape-only mode);
  - 999 permutations; one-sided p.
- **Cross-check:** on one simulated and one real subsample, the MST edge lengths must equal the GUDHI RipsComplex H0 death times to numerical precision (GUDHI is already installed and used only for this check).

### 3.5 Shared settings and multiplicity
- **Seeds:** `SeedSequence(42, crc32(unit name))` per unit (cohort × version × transition × method × mode).
- **Multiplicity:** Holm across the transitions within each cohort × data version × method (× mode for the MST test). That is 3 transitions for GSE41258 and 2 for the other cohorts, exactly the family structure of the step 5 TDA test.
- **No further correction** across methods, cohorts or versions. Each method is compared one-to-one with the TDA result for the same family, and the four cohorts are separate replications. This is stated as a limitation.
- **TDA results compared against:** the committed step 5 primary results: `results/step05/permutation_results_all.csv` and `results/step05/<GSE>/permutation_results_all.csv`, Holm-adjusted p.

### 3.6 Secondary: size-matched PERMDISP and PERMANOVA
**Purpose:** to separate power (sample size) from the property measured. The conventional tests in 3.1–3.2 use all samples, while the TDA test uses subsamples of size m.

**Method:**
- PERMDISP and PERMANOVA are recomputed on the same k = 20 subsample pairs as the step 5 subsampling scheme (size m per stage; seeds per unit as in 3.5).
- Each subsample pair gets its own 999-permutation p-value.
- **Reported:** the median statistic, median effect size (D, R²) and median p over the 20 subsamples, with Holm applied to the median p within the same families.
- The median p is a descriptive summary, not an exact p-value; this analysis is **secondary**.

**Reading:**
- If a conventional test is significant on all samples but not at size m, the TDA–conventional discordance is attributed to **power**.
- If it remains significant at size m (while TDA is not), the discordance reflects a **different property**.

## 4. Validation on simulated data (before any real data)

Simulation in 10 dimensions, n = 50 per group, 200 replicates per scenario, 999 permutations. The MST test uses k = 20 subsamples of m = 40 (80% of 50), and m = 38 in V5:

| Scenario | Expected |
|---|---|
| V1: identical Gaussian groups (SD 1) | rejection ≈ 5% for all methods and modes |
| V2: group 2 scaled × 1.5 about the same centroid | PERMDISP power ≥ 80%, D > 1 |
| V3: group 2 shifted by 1 SD along one axis, same spread | PERMANOVA power ≥ 80%; PERMDISP rejection ≈ 5% |
| V4: group 1 = one Gaussian (SD 1); group 2 = two tight clusters (SD 0.5, n/2 each) centred at ±2 along one axis (H0 structure) | MST test power ≥ 80% in **both** raw and shape-only modes. The stage-separation test is **not** expected to detect V4 (no requirement) |
| V5: as V1 but unequal sizes (n = 48 vs 181, as GSE41258 A → C) | rejection ≈ 5% for all methods and modes |
| V6: two groups (SD 1) shifted by 4 SD along one axis | stage-separation (ARI) power ≥ 80% |

**Acceptance:**
- Type I error (V1 and V5 for all methods and modes; PERMDISP in V3) ≤ 0.05 + 2.5 × Monte Carlo SE (≈ 0.089 with 200 replicates).
- The power criteria above are met.
- PERMDISP and PERMANOVA statistics match a direct textbook calculation on one simulated data set.
- The MST edge lengths equal the GUDHI H0 death times (section 3.4).
- Rejection rates of all methods in all scenarios are reported, including those without a requirement.

If a check fails, no real data are analysed; the problem is reported and recorded as a deviation.

## 5. Interpretation (fixed in advance)

### 5.1 Concordance between TDA and conventional tests
"Significant" means Holm p < 0.05 within the families of section 3.5. Three pairings are evaluated, each over the 18 comparisons:

| Pairing | Question |
|---|---|
| **P1:** TDA raw vs PERMDISP | Is the TDA spread signal captured by conventional dispersion? |
| **P2:** TDA raw vs MST raw | Is the TDA raw signal reproducible from H0 / MST edge lengths alone? |
| **P3:** TDA shape-only vs MST shape-only | Is the TDA shape signal reproducible from H0 / MST edge lengths alone? |

For each pairing, three quantities are reported:
- **(a1)** of the TDA-significant comparisons, the proportion also significant by the conventional test;
- **(a2)** of the conventionally significant comparisons, the proportion also significant by TDA;
- **(b)** Cohen's kappa over the 18 comparisons (with the 2 × 2 table).

**Pre-stated conclusions:**
1. **Largely captured.** A TDA signal is described as **largely captured by the conventional method** if (a1) ≥ 75% **and** (a2) ≥ 75%.
   - For P1, that means TDA adds little unique information about spread.
   - For P2 and P3, it means the H0 signal is reproducible without persistent homology.
   - If either proportion is below 75%, the discordant comparisons are described individually.
   - If a denominator is 0 (no significant comparisons for one method), the proportion is reported as not applicable and the "largely captured" conclusion cannot be drawn.
   - Kappa is reported descriptively; no threshold.
2. **Adenoma → Carcinoma direction.** The description "spread increases from adenoma to carcinoma" is supported by conventional analysis in a cohort only if PERMDISP is significant **and D > 1** in that cohort.
   - If PERMDISP is significant with D < 1, this is reported as a contradiction of the spread-increase wording.
   - In that case, the wording in the reports will be changed to "a difference in spread".

**Descriptive extra (no conclusion attached): MST vs H0-only TDA.**
- The MST test is the exact H0 counterpart, so the MST raw and shape-only results are also cross-tabulated against the matching H0-only TDA results of addendum 5 part 4 (commit d09a883; `results/addendum5/part4_dims/dims_results_all.csv`, Holm p, same modes).
- That analysis covered GSE41258, GSE20916 and GSE37364 only, so this cross-table has 14 comparisons per mode; GSE117606 has no H0-only TDA result and is not included.
- Reported: 2 × 2 tables, the proportions (a1)/(a2) and kappa, and z values side by side. No conclusion rule is attached, and the extra does not change the conclusions above.

### 5.2 Shape: TDA shape-only test vs conventional methods
- A TDA shape-only difference (Holm p < 0.05) is described as **"not captured by conventional methods"** only if, for the same comparison, **PERMDISP, the stage-separation test and the MST test (shape-only mode) are all non-significant**.
- Otherwise it is described as **"accompanied by conventional differences"**, with those differences named.
- PERMANOVA is reported alongside but is **not** part of this rule, because location is removed by design in the TDA test.
- This is a deliberately conservative criterion for claiming a TDA-specific signal.

### 5.3 PERMANOVA
PERMANOVA is expected to be significant for most transitions: stages differ in average expression, and the TDA test removes that difference by alignment. A significant PERMANOVA with a non-significant TDA test is therefore **not** read as lower TDA sensitivity. It shows that the two methods ask different questions (location vs distribution shape and spread).

### 5.4 What the benchmark cannot show
- Agreement does not make the methods equivalent, because they measure related but different properties.
- Disagreement can reflect different sensitivity (power) as well as different properties. The documented power limits of the TDA test (D7, D8, addendum 5 part 3) apply.
- The primary PERMDISP, PERMANOVA and stage-separation analyses use all samples, so they have more data than the TDA test in cohorts with unbalanced stages (especially GSE41258). The size-matched analysis (3.6) addresses this for PERMDISP and PERMANOVA.

## 6. Outputs

- `scripts/step13_benchmark.py`: the benchmark script; the simulation check (section 4) is part of it and runs first.
- `results/step13/validation_summary.txt`: section 4.
- `results/step13/benchmark_results.csv`: one row per cohort × version × transition × method × mode, with statistic, effect size (D, R², ARI or MST statistic), p and Holm p, next to the TDA raw and shape-only Holm p.
- `results/step13/size_matched_results.csv`: section 3.6.
- `results/step13/summary.txt`: the concordance tables (P1–P3: proportions, kappa, 2 × 2 tables), the shape rule of 5.2 and the pre-stated conclusions.
- `results/step13/benchmark_overview.png`: D, ARI and MST z per comparison, with TDA significance marked.
