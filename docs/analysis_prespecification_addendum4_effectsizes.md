# Pre-specification addendum 4: standardised effect sizes with bootstrap CIs (descriptive)

**Date:** 2026-10-06
**Status:** written **after** all step 5 permutation results (p-values and z-scores) for GSE41258, GSE20916 and GSE37364 were known. Before writing it, no effect sizes or bootstrap intervals of the kind defined here had been computed. The analysis is therefore **post hoc and descriptive**.
**Earlier documents:** `analysis_prespecification.md` (4f1f32f) and addenda 1–3 are unchanged. **The pre-specified tipping-point decision is unchanged** and is not affected by anything in this addendum.
**Rule:** this file is not changed after the computation starts. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Purpose

To show, for every step 5 transition, **how large** the stage difference is and **how precisely** it is estimated. This distinguishes a non-significant result that is consistent with *no meaningful effect* from one where the data are *too few to tell*.

## 2. Transitions, data and settings

- **Transitions** (all tested in step 5): GSE41258 `main` — Normal → Adenoma, Adenoma → Carcinoma, Carcinoma → Metastasis. GSE20916 `bulk` and GSE37364 `bulk` — Normal → Adenoma, Adenoma → Carcinoma.
- **Data:** the primary step 5 embeddings (2,000 genes, 10 PCs); **shape-only and raw**; **unadjusted and composition-adjusted** (primary, stage-blind adjustment). The stage-aware data are not included.
- **Settings:** frozen step 4/5 settings — Vietoris–Rips, H0 + H1, exact 1-Wasserstein (reference `diag_dist`); subsample size m = 80% of the smallest stage in the set (31, 19 and 21 respectively); k = 20 subsamples per stage.

## 3. Effect size

Within one computation, k = 20 subsample diagrams are drawn per stage (reference `bootstrap_diagrams`). Then:
- **between** = the 400 between-stage W1 distances (all pairs);
- **within** = the 2 × 190 within-stage W1 distances (all pairs within each stage, pooled);
- **effect size d = (mean between − mean within) / SD(within)**, with SD using ddof = 1.

The numerator is the step 5 test statistic. d expresses it in units of the natural variability between subsample diagrams of the same stage. **Point estimate:** one computation on the original samples.

## 4. Bootstrap confidence interval

- **1,000 resamples.** In each, the samples of each stage are resampled with replacement (stage sizes kept), and d is recomputed exactly as in §3 on the resampled stages.
- **95% CI:** percentile interval (2.5th and 97.5th percentiles). The bootstrap SE is the SD of the 1,000 values.
- **Seed:** 42. As in step 5 (deviation D1), each unit (cohort × version × transition × mode) uses `SeedSequence(42, crc32(unit name))`, so results do not depend on run order and an interrupted run can resume.
- **Known limitation, stated in advance:** a resampled stage contains duplicate samples, so a subsample can contain identical points (zero distances). This shrinks the effective number of distinct points and can bias d. It is reported with the results; the interval is a descriptive guide, not an exact inference.

## 5. Pooling across cohorts

For Normal → Adenoma and Adenoma → Carcinoma (three cohorts each), separately for each mode × version: DerSimonian–Laird random-effects meta-analysis of d with within-cohort variance = bootstrap SE², using `statsmodels.stats.meta_analysis.combine_effects` (method_re = "dl"). Reported: pooled d with 95% CI, τ², I². Carcinoma → Metastasis has one cohort only and is not pooled. With three cohorts, τ² and I² are imprecise.

## 6. Plain-language categories (fixed in advance)

The threshold for a "meaningful" effect is **d = 0.5**, the conventional "medium" effect size. It is a convention, not derived from these data.

| Category | Rule |
|---|---|
| **Difference detected** | Step 5 Holm-adjusted p < 0.05 (pre-specified test) |
| **Likely no meaningful effect** | Holm p ≥ 0.05 **and** upper 95% bound of d < 0.5 |
| **Too few samples to tell** | Holm p ≥ 0.05 **and** the 95% CI of d includes 0.5 |
| **Large effect but test not significant** | Holm p ≥ 0.05 **and** lower 95% bound of d > 0.5 (reported as a discordance, not as a finding) |

The same labels are applied to the pooled estimates, using the pooled CI. Pooled estimates have no pre-specified test, so "difference detected" is not used for them; a pooled CI entirely above 0.5 is labelled "consistent meaningful effect across cohorts".

## 7. Role

**Descriptive.** These effect sizes and categories do not replace or modify the step 5 p-values, the Holm correction, or the pre-specified tipping-point definition.
