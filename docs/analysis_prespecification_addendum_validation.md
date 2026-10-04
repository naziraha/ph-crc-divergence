# Pre-specification addendum: replication in GSE20916 and GSE37364

**Date:** 2026-10-04
**Status:** written after the GSE41258 primary and sensitivity results (step 5), and **before any data from GSE20916 or GSE37364 has been downloaded or looked at**.
**Relation to the main pre-specification:** `docs/analysis_prespecification.md` (commit 4f1f32f) is unchanged and remains the plan for GSE41258. This addendum adds two independent cohorts.
**Rule:** this file is not changed after step 1 starts for either cohort. Deviations go into `docs/analysis_deviations.md`.

**Why it was added:** in GSE41258 no transition met the tipping-point definition. Normal → Adenoma was significant in shape-only only before composition adjustment, and in raw distances in both versions. This addendum tests in independent data:
- whether these findings replicate;
- whether the shape change is explained by tissue composition, using GSE20916's microdissected samples, which are physically enriched for epithelium.

---

## 1. Cohorts and analysis sets

Each cohort, and each subset of GSE20916, is analysed completely separately: no merging, no batch correction.

| Cohort | Platform | Analysis sets |
|---|---|---|
| GSE20916 | GPL570 | **(a) microdissected samples** and **(b) bulk (macro-dissected or whole-tissue) samples**, analysed as two separate sets |
| GSE37364 | GPL570 | One set |

- **Stages:** Normal, Adenoma, Carcinoma, assigned from the GEO sample characteristics and titles. Other sample types (e.g. cell lines, metastases, non-colorectal tissue) are excluded. The exact label → stage mapping is recorded in step 1 for each cohort, before any expression-based analysis.
- **Tested transitions:** Normal → Adenoma and Adenoma → Carcinoma. Carcinoma → Metastasis cannot be tested in these cohorts.

## 2. Frozen pipeline (steps 1–5, identical to GSE41258)

- **Step 1:**
  - probes mapped to gene symbols with the platform annotation (GPL570); probes with no symbol or several symbols dropped;
  - log2 if not already on log scale (99th percentile > 50);
  - highest-mean probe per gene.
- **Step 2, tissue identity and QC:**
  - The same fixed colon-identity panel as in GSE41258 (CDX1, CDX2, CDH17, VIL1, KRT20, CEACAM5, LGALS4) is used to flag samples whose expression contradicts their label.
  - Exclusions are decided with the PI, on tissue identity only, before steps 3–5. They are recorded with reasons; a sensitivity set puts them back.
  - Where a patient id is available: one array per patient per stage, using the same rule as GSE41258 (prefer the original array; re-hybridisations are dropped).
- **Step 2, composition:**
  - Epithelial score = the **same 18 pan-epithelial genes** frozen in GSE41258 (not re-selected).
  - Stromal and immune scores = the same PanglaoDB marker sets.
  - decoupler ULM on gene-centred log2 data.
  - Composition-adjusted data: per gene, OLS on the epithelial, stromal and immune scores, with stage **not** in the model (as in GSE41258). Adjusted = residual + gene mean.
  - No liver or lung covariate and no organ-gene removal (no metastases in these sets).
- **Step 3:** drop the lowest-expressed 25% of genes; top 2,000 most variable genes; z-score; PCA with 10 components. Fitted on all samples of the set pooled, without stage labels (`tda_crc_pilot.preprocess`).
- **Steps 4–5:**
  - Vietoris–Rips, Euclidean, H0 + H1, maximum dimension 2.
  - Subsamples of 80% of the smallest stage in the set, without replacement.
  - Shape-only (primary) and raw (secondary).
  - Exact 1-Wasserstein; statistic = mean between-stage minus mean within-stage distance, with k = 20 subsamples per group.
  - Groups aligned before permuting; 999 permutations.
  - **Holm correction across the two transitions** within each set × data version × mode.
  - Seed 42; per-unit seeds as in the GSE41258 step 5 script.
- **Minimum group size:** a comparison involving a stage with fewer than 13 samples in a set (subsample < 10) is reported descriptively only (step 4), without a permutation test.

## 3. Definitions

- **Topological tipping point (per set):** shape-only Holm p < 0.05 in **both** the unadjusted and composition-adjusted data. This is the same definition as for GSE41258.
- **Replication of the GSE41258 heterogeneity finding (secondary):** Normal → Adenoma raw-distance Holm p < 0.05 in both data versions.
- **Interpretation of GSE20916:**
  - Shape-only Normal → Adenoma significant in the **microdissected** set (both versions) → evidence that the shape change belongs to the epithelium itself.
  - Significant only in the bulk set → supports the composition explanation found in GSE41258.

## 4. Sensitivity analyses (per set; descriptive)

- 5 and 20 PCs; 1,000 and 5,000 genes (same as GSE41258 S5–S8).
- One sample per patient per comparison (patient-disjoint pairs), where patient ids are available.
- QC exclusions put back, if any are made.

## 5. What is not used

- The GSE41258 results are not used to tune any setting for these cohorts. All settings above are those frozen before GSE41258 step 4.
- The pilot results are not used.
