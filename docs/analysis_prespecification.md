# Analysis pre-specification: topological tipping point in GSE41258

**Date:** 2026-10-03 (S9 added and document approved by the PI on 2026-10-04, still before step 4)
**Status:** written after steps 1–3 (data preparation, composition scoring, embedding) and **before step 4**. No persistent homology, Wasserstein distance or permutation test has been run on the GSE41258 data in this project.
**Rule:** this file is not changed after step 4 starts. Any later deviation goes into a separate file, `docs/analysis_deviations.md`, giving the date, the reason and the effect on the results.

---

## 1. Primary data set

- **Cohort:** GSE41258 (Affymetrix GPL96), discovery cohort.
- **Primary set:** `main`, with one array per patient per stage.

| Stage | n |
|---|---|
| Normal colon mucosa | 53 |
| Adenoma | 48 |
| Primary carcinoma | 181 |
| Metastasis (liver) | 39 |
| **Total** | **321** |

- **How the set was built** (details in `data/processed/GSE41258_exclusions.csv` and `results/step02/summary.txt`):
  - **Stage labels** come from the GEO "tissue" field, using `labels_from_meta` from `reference/tda_crc_pilot.py`. Lung metastases, normal liver, normal lung and cell lines are excluded.
  - **Tissue-identity exclusions (2):** GSM1012647 (expression is that of normal liver) and GSM1012526 (Müllerian/ovarian profile). These were based on tissue identity only and decided with the PI before any topology analysis.
  - **One array per patient per stage:** 15 additional arrays were removed from `main`, plus 1 lung-metastasis array from `withlung`. Patients are identified by the official GEO "patient id" field. Where a choice existed, the array the original authors flagged "included in analysis: Yes" was kept; otherwise the first original array (not an `_ez`/re-hybridisation) was kept. Decided with the PI before any topology analysis.

## 2. Embedding (step 3; fixed)

- **Expression:** probes mapped to HGNC symbols, highest-mean probe per gene, log2.
- **Gene selection and PCA:** drop the lowest-expressed 25% of genes, keep the top 2,000 most variable genes, z-score each gene, PCA with 10 components (`tda_crc_pilot.preprocess`). Gene selection and PCA are fitted on **all samples of the set pooled, without stage labels**.
- **Two data versions, both analysed:**
  - **Unadjusted:** step-1 expression.
  - **Composition-adjusted:** each gene regressed (OLS) on the step-2 epithelial, stromal and immune scores and the liver fraction; residual + gene mean.
- **Epithelial score:** 18 pan-epithelial structural genes whose stage medians are within 1 log2.
- **Stromal and immune scores:** PanglaoDB canonical markers.

## 3. Liver-gene handling

- **Liver signature:** genes at least 8-fold (3 log2) higher in the 13 reference **normal liver** samples than in **both normal colon and primary carcinoma**, with liver mean ≥ 8 log2. That gives 176 genes. **Metastasis samples were not used to define the signature.**
- **Adjusted data only:** the liver-signature genes are removed **before** the top 2,000 genes are chosen. For the `withlung` set, the lung signature (7 normal lung samples) is also removed.
- **Unadjusted data:** no genes are removed.
- **Known limitation, recorded before step 4:** in the adjusted data, within the 39 metastases, some PCs still correlate with liver fraction (Spearman r up to 0.49). The genes involved are non-signature genes: immunoglobulins, blood, liver-enriched metabolic genes below the signature threshold, and inversely desmoplastic stroma genes. This may partly reflect real tumour–liver interface biology (histological growth pattern) and cannot be fully separated from residual liver tissue. The carcinoma → metastasis result is interpreted with this in mind.

## 4. Primary test (steps 4–5)

- **Persistent homology:**
  - Vietoris–Rips complex on Euclidean distances in the 10-PC space, maximum dimension 2, homology dimensions H0 and H1 (GUDHI).
  - **Subsamples:** equal size for both stages of a comparison, 80% of the smallest stage in the set, drawn without replacement (`main`: 0.8 × 39 → 31 samples).
- **Primary statistic: shape-only.** Pairwise distances in each subsample are divided by that subsample's median pairwise distance.
- **Distance between diagrams:** exact 1-Wasserstein distance (L∞ ground metric), summed over H0 and H1.
- **Test statistic:** mean between-stage distance minus mean within-stage distance. It is computed from 20 subsample draws per group (exhaustive pairs; `tda_crc_pilot.test_statistic`, with the pilot's default `k_draws = 20`).
- **Permutation test:**
  - Each group is first centred, and in shape-only mode scaled to unit median pairwise distance (`tda_crc_pilot.align_groups`).
  - Stage labels are then permuted 999 times with group sizes preserved: p = (1 + #{T_perm ≥ T_obs}) / (1 + 999).
- **Multiple testing:** Holm correction across the three adjacent transitions:
  1. Normal → Adenoma
  2. Adenoma → Carcinoma
  3. Carcinoma → Metastasis
- **Random seed:** 42.
- **Secondary (reported, not used for the primary decision):** the raw-distance version of the same test, and descriptive effect sizes (topological separation index, Betti curves, total persistence).
- **Code:** `reference/tda_crc_pilot.py`, unchanged, validated by `reference/validate_tda_pipeline.py` (10/10 checks passed on 2026-10-03 in the project environment).

## 5. Definition of a topological tipping point

A transition is called a **topological tipping point** if its shape-only permutation test has a **Holm-adjusted p < 0.05 in both the unadjusted and the composition-adjusted data** of the primary set.

- A transition significant in only one data version is reported as such, and is **not** called a tipping point.
- If more than one transition meets the definition, all are reported. This definition does not rank them.

## 6. Sensitivity analyses

Each is run with the same test on both data versions. They are interpreted descriptively, as support for or against the primary result, and do not replace it.

| # | Analysis | Change from primary |
|---|---|---|
| S1 | `withlung` | Lung metastases added to metastasis (53 / 48 / 181 / 58); lung signature also removed in the adjusted data |
| S2 | `main_allarrays` | All arrays kept, no de-duplication (54 / 51 / 184 / 47) |
| S3 | `main_inclqc` | The two tissue-identity exclusions put back (53 / 48 / 183 / 39) |
| S4 | One sample per patient | Each compared pair of stages made patient-disjoint (`tda_crc_pilot.drop_shared_patients`, official GEO patient id; for each shared patient, a random choice of which group keeps them; seed 42) |
| S5 | 5 PCs | 2,000 genes, 5 PCs |
| S6 | 20 PCs | 2,000 genes, 20 PCs |
| S7 | 1,000 genes | 1,000 genes, 10 PCs |
| S8 | 5,000 genes | 5,000 genes, 10 PCs |
| S9 | Carcinoma → metastasis without high-liver metastases | Carcinoma → Metastasis only. The 3 metastases with step-2 liver fraction > 20% are removed: GSM1012485 (47.3%), GSM1012650 (44.2%), GSM1012476 (30.7%), leaving 36 metastases. The same `main` embedding is used (not refitted). Subsample size = 80% of the smaller group (0.8 × 36 → 28). Single comparison: unadjusted p reported (no Holm). |

## 7. Pilot results

The earlier pilot analysis (`reference/`, results not part of this project) is **not used for any decision**: not for settings, sample selection, the tipping-point definition or interpretation. Only its validated code is reused.

## 8. Wording

- Never claim clonal evolution; at most "consistent with".
- The late stage is called "metastasis".
