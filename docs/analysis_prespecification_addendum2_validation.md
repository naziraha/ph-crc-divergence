# Pre-specification addendum 2: stage-aware composition adjustment and GSE20916 microdissected set

**Date:** 2026-10-04
**Status:** written after step 1 for GSE20916 (data preparation and label mapping only), and **before step 2 or any expression-based analysis of GSE20916 or GSE37364**. GSE37364 has been downloaded but not opened.
**Relation to earlier documents:** `docs/analysis_prespecification.md` (commit 4f1f32f) and `docs/analysis_prespecification_addendum_validation.md` (commit 6ea31ea) are unchanged. In particular, the **primary tipping-point definition is unchanged**.
**Rule:** this file is not changed after step 2 starts for either cohort. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Secondary analysis: stage-aware composition adjustment (GSE20916, GSE37364)

**Why:** the primary adjustment regresses each gene on the composition scores **without stage in the model**. Because composition differs between stages, this can remove real stage-related expression together with composition effects (over-adjustment). The stage-aware version estimates composition effects **within stages**, which protects stage differences.

**Model** (per gene, ordinary least squares, all samples of the set):

expression = intercept + stage effects (indicator variables for Adenoma and Carcinoma; Normal = reference) + b_epi × epithelial + b_str × stromal + b_imm × immune + error

- **Stage-aware adjusted value** = expression − [b_epi × (epithelial − mean) + b_str × (stromal − mean) + b_imm × (immune − mean)]. Only the composition part is removed; stage effects and residual variation are kept.
- The composition scores and marker genes are identical to the primary adjustment.

**Analyses:** steps 3–5 are repeated on the stage-aware adjusted data with the frozen settings. Shape-only and raw distances are reported; Holm correction is applied across the two transitions.

**Role:** secondary. It is reported next to the unadjusted and primary-adjusted results. It is **not** part of the tipping-point definition, which remains: shape-only Holm p < 0.05 in both the unadjusted and the primary (stage-blind) adjusted data.

**Caveat, stated in advance:** the stage-aware version uses the stage labels inside the adjustment, so it can under-correct when composition and stage are strongly confounded. A result that is significant only in the stage-aware data is reported as such, not as a tipping point.

## 2. GSE20916 microdissected set: which samples, and whether to pool preparations

The microdissected set has two preparations: laser-captured colonic crypt epithelial cells and microdissected mucosa, each with 5 distant normal, 5 adjacent normal, 5 adenoma and 5 carcinoma samples. By the minimum-size rule of addendum 1, it is analysed **descriptively only** (no permutation test).

- **Normals:** the **distant normal** samples are the normal group for analysis.
  - Adjacent (tumour-adjacent) normals are kept for description only: composition scores, and projection onto the PCA embedding in figures.
  - They are **not** used to fit the adjustment, the gene selection or the PCA.
- **Pooling rule** (evaluated once, at step 3, on the unadjusted data):
  1. Fit the step-3 embedding (2,000 genes, 10 PCs) on the 30 analysis samples (distant normal, adenoma, carcinoma; both preparations).
  2. Compute the pseudo-R² of stage and of preparation: the share of the total sum of squared Euclidean distances (10-PC space) explained by each grouping, as in a one-way PERMANOVA.
  3. **Pool the preparations if R²(stage) > R²(preparation)**; otherwise describe each preparation separately (5 / 5 / 5 per stage).

  The R² values and the decision are recorded in the step 3 summary.
- **Composition adjustment:** fitted on the analysis samples of the set actually analysed (pooled set, or each preparation separately). Both primary and stage-aware versions are produced for the pooled set and for each preparation, so that the pooling decision does not require re-running step 2.

## 3. Identity QC in the validation cohorts

There are no normal-liver or normal-lung reference samples in these cohorts. So the identity screen of addendum 1 is applied **within each set**, using the same fixed panels:
- colon: CDX1, CDX2, CDH17, VIL1, KRT20, CEACAM5, LGALS4;
- hepatocyte: ALB, APOA1, APOA2, APOH, FGA, FGB, FGG, AHSG, TF, HP, SERPINA1, TTR;
- lung: SFTPC, SFTPB, SFTPD, SCGB1A1, AGER;
- Müllerian: PAX8, SOX17, MUC16.

A sample is **flagged** if its colon panel mean is more than 3 robust SDs (1.4826 × MAD) below the set median, or if its hepatocyte, lung or Müllerian panel mean is more than 3 robust SDs above the set median. Flagged samples are shown to the PI, who decides on exclusion (tissue identity only) before step 3.
