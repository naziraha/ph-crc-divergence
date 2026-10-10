# Addendum 10: TCGA-COAD carcinoma stages I → IV with the frozen pipeline

**Date:** 2026-10-10 (locked at commit, before any TCGA computation)
**Status:** written before any TCGA analysis in this project. Known so far:
- the results of the 4 GEO cohorts;
- an earlier, separate TCGA pilot (OMICLINK project; Stage I vs Stage IV; bootstrap Mann–Whitney on total H1 persistence; reported p < 0.0001).

This is therefore an **informed** analysis. Its rules are fixed here before any TCGA data are processed in this project.
**Relation to earlier documents:** the main pre-specification and addenda 1–9 are unchanged. In the plan, TCGA is listed for "stage and survival only"; this addendum uses it for stage.
**Rule:** this file is not changed after the computation starts. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Questions
1. **Tipping point across carcinoma stages:** does the frozen step 1–5 pipeline find a topological tipping point between carcinoma stages I → II → III → IV in TCGA?
2. **Same pattern as the GEO cohorts?** Is there a spread or heterogeneity change between stages, and is the signal in H0 rather than H1?
3. **Why did the pilot look strongly significant?** The pilot's test (bootstrap subsamples treated as independent) is compared with a patient-level permutation test on the same statistic.

**Note:** these are carcinoma stages (extent of spread: TNM), not the normal → adenoma → carcinoma sequence. A null result here would not contradict the GEO findings; a positive result would add to them.

## 2. Data
- **Source:** UCSC Xena TCGA-COAD, already downloaded in the OMICLINK project:
  - expression: `01_GeneExpression/HiSeqV2_PANCAN_COAD.xlsx` (RSEM, log2(x + 1), PANCAN-normalised);
  - clinical: `02_Clinical/TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt`;
  - mutations: `04_Mutation/COAD_mc3_gene_level.txt`.

  These files are **copied unchanged** into `data/raw/TCGA/` (with SHA-256 checksums recorded). The originals are not modified.
- **Samples:**
  - primary tumours only (sample code `-01`); one sample per patient (if a patient has several, the first in sorted order);
  - **stage:** `pathologic_stage` collapsed to I, II, III, IV (sub-stages A/B/C merged); samples with missing or unknown stage are excluded;
  - **solid-tissue normals** (`-11`) are used only in the secondary analysis (§4.3).
- **Rectal cancer (READ):** **not analysed** (PI decision 2026-10-10: COAD only for now). Any later READ analysis would need its own addendum.
- **Sample counts** per stage are reported after step 1 and checked at every later step, as for the GEO cohorts.

## 3. Pipeline (frozen settings, unchanged)
1. **Step 1:** genes as given, already log2; genes with missing values removed; lowest-expressed 25% dropped.
2. **Step 2:** composition scores with decoupler ULM, using the same marker sets (18 pan-epithelial genes; PanglaoDB stromal and immune). Stage-blind OLS adjustment gives the **adjusted** version (primary, together with unadjusted).
3. **Step 3:** top 2,000 most variable genes, z-score, PCA 10 components, fitted without labels on all included tumours. Sensitivity: 5 / 20 PCs; 1,000 / 5,000 genes.
4. **Step 4:** Vietoris–Rips persistence with GUDHI, H0 and H1, max dimension 2.
5. **Step 5:** the reference `permutation_test`, unchanged:
   - k = 20 subsamples of m = 80% of the smallest stage;
   - shape-only and raw modes;
   - groups aligned; 999 permutations;
   - **Holm across the three transitions** I → II, II → III, III → IV;
   - per-unit seeds `SeedSequence(42, crc32(unit name))`.

**Tipping point (same definition as the GEO cohorts):** a transition with shape-only Holm p < 0.05 in **both** unadjusted and adjusted data.

## 4. Additional analyses

### 4.1 Spread and H0 vs H1
- **PERMDISP** (bias-adjusted distances, as addendum 9 Amendment 1; 999 permutations; Holm across the 3 transitions) per transition, with the dispersion ratio D.
- **H0-only and H1-only** versions of the step 5 test (as addendum 5 part 4), unadjusted and adjusted.

### 4.2 MSI and subtype
MSI-H tumours have very distinct transcriptomes and are unevenly spread across stages, so they could drive stage differences.
- **MSI status:** from the clinical matrix MSI field if present. Otherwise, hypermutation from the mc3 file (≥ 200 mutated genes = MSI-H-like; the pilot showed a gap between 130 and 248).
- **Sensitivity analysis:** steps 3–5 repeated on **MSS / non-hypermutated tumours only**, with the same tipping-point rule.
- **CMS:** CMScaller (r-cms environment, as addendum 9) on the unadjusted data; `RNAseq = FALSE` with gene-wise centring, because the Xena data are already log2-normalised. Reported: CMS proportions by stage (descriptive), with the template coverage flag (< 80%).

### 4.3 Secondary: solid-tissue normal → Stage I (descriptive)
If ≥ 20 solid-tissue normals (`-11`) are available, the step 5 test is also run for Normal → Stage I. It is reported separately and is **not** part of the tipping-point family or the Holm correction. This is adjacent normal tissue vs carcinoma; there are no adenomas in TCGA.

### 4.4 Method demonstration: the pilot's test vs a patient-level test (Stage I vs Stage IV)
On the unadjusted 10-PC embedding, for the pilot's statistic (total H1 persistence per subsample, k = 25 patients per stage):
- **(a) Pilot method:** 100 bootstrap subsamples per stage, Mann–Whitney test across subsamples. Repeated with 1,000 subsamples to show how p depends on the number of subsamples.
- **(b) Patient-level permutation test:** the difference in mean total H1 persistence (20 subsamples per stage), with stage labels permuted across patients 999 times.
- **(c) Shape-only version** of (b): distances scaled by the median pairwise distance.

Reported side by side; descriptive. The pilot's original numbers are not re-used.

## 5. Interpretation (fixed in advance)
| Result | Reading |
|---|---|
| Tipping point met for a transition (§3 rule) | **TCGA shows a topological tipping point between those carcinoma stages**; checked against the MSS-only analysis |
| Shape-only significant in one data version only | **Inconsistent (as in the GEO cohorts)**; not a tipping point |
| No shape-only significance | **No tipping point between carcinoma stages** |
| Raw or PERMDISP significant | A spread difference between stages, reported with its direction (D) |
| H0-only significant, H1-only not | Signal in clustering and spread, as in the GEO cohorts |
| Pilot method significant, permutation test not | The pilot result reflects pseudo-replication, not a patient-level difference |

**Robust** means the same call in the MSS-only analysis. All wording is "consistent with"; no claim of clonal evolution. The late stage is called "Stage IV" here (metastatic disease at diagnosis).

## 6. Validation
- **No new methods,** except the §4.4 demonstration, which re-uses validated functions. The step 5 test, PERMDISP and CMScaller were validated earlier (reference validation; step 13; addendum 9).
- **Before step 5:** a **positive-control check** on TCGA (addendum 5 level 1 type: a known shift inserted into the embedding of a random half of Stage II tumours, which must be detected). This confirms the pipeline works on this data type.

## 7. Outputs
- **Data:** `data/processed/TCGA_COAD_*` (expression, labels).
- **Results by step:** `results/step01/TCGA/`, `results/step02/TCGA/`, `results/step03/TCGA/` (embeddings in `results/step03/embeddings/TCGA_COAD_*`), `results/step04/TCGA/`, `results/step05/TCGA/`.
- **Addendum results:** `results/addendum10/` (spread, H0/H1, MSI, CMS, normal → Stage I, method demonstration, `summary.txt`).
- **Scripts:** `scripts/addendum10_tcga_prepare.py` (steps 1–2); existing step 3–5 scripts extended with a TCGA configuration; `scripts/addendum10_tcga_extra.py` (§4).
