# Pre-specification addendum 3: random-effects meta-analysis of Normal → Adenoma (secondary)

**Date:** 2026-10-04
**Status:** written **after** steps 1–5 had been completed for GSE41258 and for GSE20916, so their Normal → Adenoma results were known when this was written. For GSE37364 only steps 1–2 (data preparation, identity screen, composition) had been run; no embedding, topology or test results existed.
**Earlier documents:** `analysis_prespecification.md` (4f1f32f), addendum 1 (6ea31ea) and addendum 2 (d73e5d1) are unchanged. The **primary tipping-point definition is unchanged**, and this meta-analysis cannot redefine it.
**Rule:** this file is not changed after step 3 starts for GSE37364.

---

## 1. Role

**Secondary analysis.** It summarises the strength and consistency of the shape-only Normal → Adenoma difference across the three cohorts. Because two of the three inputs were already known, it is an informed summary, not an independent confirmation, and it is reported as such.

## 2. Inputs

One value per cohort and data version: the **z-score of the shape-only Normal → Adenoma permutation test**, z = (T_obs − mean of the permutation null) / SD of the null, from the primary step 5 analysis at frozen settings (2,000 genes, 10 PCs, 999 permutations).

| Cohort | Set | Unadjusted | Adjusted (primary, stage-blind) | Stage-aware |
|---|---|---|---|---|
| GSE41258 | `main` (53 / 48 / 181 / 39) | step 5 primary | step 5 primary | exploratory analysis D3 (`docs/analysis_deviations.md`), to be run with the frozen step 5 settings |
| GSE20916 | `bulk` (24 / 45 / 36) | step 5 primary | step 5 primary | step 5 secondary |
| GSE37364 | `bulk` (all bulk samples after any PI identity exclusions) | step 5 primary | step 5 primary | step 5 secondary |

The microdissected GSE20916 set has no permutation test (minimum-size rule), so it is not included.

## 3. Model

- Random-effects meta-analysis, **DerSimonian–Laird** estimator of the between-cohort variance τ², run separately for each data version (unadjusted, adjusted, stage-aware).
- Each cohort's z-score is the effect estimate, with **within-cohort variance 1**. The z-scores are standardised against their permutation null, so 1 is their sampling variance under the null. Extra between-cohort spread, including the inflation of z-score variance under the alternative, is absorbed by τ².
- **Reported:** pooled z with 95% CI; two-sided p from the normal distribution; τ²; Cochran's Q with p; I²; and each cohort's z with weight. Implemented with `statsmodels.stats.meta_analysis.combine_effects` (method_re = "dl").
- With only three cohorts, τ² and I² are imprecise; this is stated with the results.

## 4. Interpretation (fixed in advance)

- The pooled z describes the average standardised strength of the Normal → Adenoma shape difference across cohorts.
- A pooled 95% CI excluding 0 in the **adjusted** data is reported as "consistent evidence of a composition-robust shape difference across cohorts". It is **not** called a tipping point, because the pre-specified per-cohort definition governs that term.
- Results for unadjusted, adjusted and stage-aware data are always shown together.
