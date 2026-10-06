# Addendum 6: GSE117606 as fourth cohort (replication of steps 1–5)

**Date:** 2026-10-06 (approved by the PI after review of the draft)
**Status:** written **before any expression data of GSE117606 have been downloaded or opened**. Only GEO series and sample metadata (titles, characteristics) have been read: `data/raw/GSE117606_series_brief.txt`, `GSE117606_samples_brief.txt`.
**Relation to earlier documents:** `analysis_prespecification.md` (4f1f32f) and addenda 1–5 are unchanged. This addendum applies the frozen pipeline of addenda 1 and 2 to one further cohort. **The tipping-point definition is unchanged.**
**Rule:** this file is not changed after step 1 starts for GSE117606. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Purpose

Independent replication of plan steps 1–5 in a fourth cohort, approved by Prof Peng Weng Kung.

GSE117606 is the largest replication cohort so far, and the only one in which most patients contribute normal-adjacent, adenoma and carcinoma tissue. It therefore also allows a within-patient (paired) secondary analysis.

## 2. Samples

**Cohort (from GEO metadata):**
- Platform GPL25373: Affymetrix HT HG-U133+ PM array plates, summarised to genes with the Brainarray custom CDF `HTHGU133PlusPM_Hs_ENTREZG_20`.
- FFPE blocks from Avaden Biosciences (USA); submitted 2018.
- 208 arrays from 75 patients (official GEO field `patient id`): tissue "Adjacent Mucosa" 65, "Adenoma" 59, "SSA" 10, "Tumor" 74.

**Stages:**

| GEO `tissue` | Stage in this analysis |
|---|---|
| Adjacent Mucosa | **Normal (tumour-adjacent)** |
| Adenoma (conventional: tubular, tubulovillous, villous, with or without high-grade dysplasia) | **Adenoma** |
| Tumor | **Carcinoma** |
| SSA (sessile serrated adenoma, 10 arrays) | Excluded from the main analysis; described separately (section 2.4) |

There is no metastasis in this cohort. Only Normal → Adenoma and Adenoma → Carcinoma are tested.

### 2.1 One sample per patient per stage
- **The rule:** one array per patient per stage, using the official GEO patient ID.
- **Duplicates in the metadata:** only three patients (1818, 1889, 2488) contribute two "Tumor" arrays:
  - an array coded **B** ("Tumor", diagnosis adenocarcinoma), the code used for the tumour in all other patients;
  - an array coded **A** (diagnosis "Carcinoma"), the code used for the adenoma lesion in all other patients.
- **Which array is kept:** for these patients, the **B-coded array is the carcinoma** of the main analysis. The A-coded carcinoma is set aside: it is neither a conventional adenoma nor the patient's designated tumour.
- **Set-aside arrays**, listed in the sample table as "set aside: A-coded carcinoma": GSM3304801 (1889A), GSM3304813 (2488A) and GSM3304818 (1818A). These three patients therefore have no conventional adenoma in the main set.
- **Any other duplicate found later** (e.g. a re-hybridisation): keep the array with the lowest GSM number, unless it is flagged by the array-quality rule (section 7.2), in which case keep the other one.

**Expected main-analysis counts** (before identity and quality review): **Normal 65, Adenoma 59, Carcinoma 71** (patients: 65, 59, 71).

### 2.2 Minimum group size
13 per stage, as in addendum 1. Expected sizes are well above this.

### 2.3 Identity screen
As in addendum 2 §3: fixed colon, hepatocyte, lung and Müllerian panels; flag at 3 robust SD within the cohort. Flagged arrays are reviewed by the PI; exclusions are on tissue identity only and are decided before step 3. A sensitivity set puts them back.

### 2.4 Sessile serrated adenomas (descriptive only)
The 10 SSA arrays (from 10 different patients) are scored for composition (step 2) and **projected** onto the main step 3 embedding (fitted without them). Their position and within-group spread are described relative to the three main stages. No test is run (n = 10 < 13), and they are not used in any fit.

### 2.5 Limitations stated in advance
- **Normals are tumour-adjacent mucosa** from patients with cancer, not from patients without neoplasia. Field effects may make them more similar to neoplastic tissue than normals in GSE41258, GSE20916 or GSE37364.
- **FFPE tissue:** RNA degradation may add noise or array-quality artefacts; see section 7.
- **Gene-level platform summary** (Brainarray Entrez-gene CDF) instead of probe sets, unlike the other three cohorts.

## 3. Pipeline (frozen; identical to the other cohorts)

| Step | Setting |
|---|---|
| 1 | Series matrix read with `tda_crc_pilot.read_series_matrix`; log2 if the 99th percentile > 50 (rule unchanged); one row per gene (see below) |
| 2 | Identity screen; composition scores; adjustment (section 4) |
| 3 | Drop lowest-expressed 25% of genes; **top 2,000** most variable; z-score; **10 PCs** (`tda_crc_pilot.preprocess`), fitted on all main-analysis samples pooled, without labels |
| 4 | Vietoris–Rips, Euclidean, H0 + H1, max dimension 2; raw and shape-only; subsamples of m = 80% of the smallest stage (expected ⌊0.8 × 59⌋ = 47), B = 100 (descriptive) |
| 5 | 1-Wasserstein (H0 + H1); statistic = mean between − mean within over k = 20 subsamples; groups aligned; **999 permutations**; **Holm across the 2 transitions** (N → A, A → C) per data version × mode; per-unit seeds `SeedSequence(42, crc32(unit name))` |

**The only permitted change: gene identifiers.**
- GPL25373 rows are Brainarray Entrez-gene probe sets (`<EntrezID>_at`).
- **Mapping rule:**
  1. remove Affymetrix control probe sets (`AFFX…`);
  2. take the gene symbol from the GPL25373 annotation table, or, if the table has no symbol column, map the Entrez ID to the official symbol with the NCBI `Homo_sapiens.gene_info` file (download date and file checksum recorded);
  3. drop IDs without a current symbol;
  4. if several IDs map to the same symbol, keep the one with the highest mean expression (the frozen rule).
- The number of rows dropped at each stage is reported.

**No change because of the part 3b findings.** The number of PCs (10) and the gene selection (top 2,000) are **not** changed for this cohort. Any variant (e.g. 50 PCs, a different gene selection) would be **exploratory only**, labelled as such, and recorded as a deviation before being run.

**Sensitivity analyses:** as in addendum 1 §4: 5 and 20 PCs; 1,000 and 5,000 genes; identity exclusions put back. Plus the array-quality sensitivity of section 7.2.

## 4. Composition adjustment

- **Marker lists frozen from GSE41258** (`results/step02/marker_sets_used.csv`): 18 pan-epithelial genes; PanglaoDB stromal (292) and immune (447) sets.
- **Scoring:** decoupler ULM on gene-centred log2 data.
- **No liver score** (no metastases).
- **Primary adjustment: stage-blind,** exactly as in all cohorts. Per-gene OLS on the epithelial, stromal and immune scores; adjusted = residual + gene mean.
- **Secondary adjustment: stage-aware,** exactly as in addendum 2 §1. Per-gene OLS on stage indicators (Adenoma, Carcinoma; Normal = reference) + the epithelial, stromal and immune scores. Only the composition part is removed: y − Σ b_c (c − mean c).
  - Steps 3–5 are repeated on these data with the frozen settings (shape-only and raw; Holm across the 2 transitions).
  - Secondary; not part of the tipping-point definition.
- **Marker genes missing on this platform** (after the mapping rule) are listed for each compartment, with counts present / expected.

## 5. Tipping-point definition (unchanged)

A transition is a topological tipping point if its **shape-only** test has **Holm-adjusted p < 0.05 in both the unadjusted and the composition-adjusted data**.

## 6. Patient overlap: primary and secondary tests

**Primary (comparability with the other cohorts):** the same **unpaired** test as in all other cohorts, on the one-per-patient-per-stage main set. Many patients contribute to two stages (51 patients have both normal and adenoma; 56 have both adenoma and carcinoma). Each patient's samples are therefore in both compared groups, which the unpaired test treats as independent.

**This patient overlap makes the unpaired test conservative:** shared patient effects make the stages look more similar to each other. A weaker result in GSE117606 is therefore **not evidence against the earlier cohorts**.

**Secondary: within-patient restricted permutation**
- **Paired subsets:**
  - N → A: the 51 patients with both a normal-adjacent and a conventional adenoma array;
  - A → C: the 56 patients with both an adenoma and a carcinoma array.
- **Test:** the step 5 statistic (H0 + H1; raw and shape-only; unadjusted and adjusted), with m = ⌊0.8 × n_pairs⌋ (40 and 44).
- **Pairing is kept in every subsample:** each of the k = 20 subsamples draws **m patients** (without replacement), and **both arrays of each drawn patient** are included. The two stage subsamples therefore always contain the same patients. Because the reference `test_statistic` draws the two groups independently, this needs a small wrapper that draws patient sets. Persistence, Wasserstein distance (`persistence`, `diag_dist`) and the between − within statistic are otherwise the unchanged reference functions.
- **Permutation:** stage labels are swapped **within each patient** (each patient's two arrays exchange labels with probability ½), 999 permutations; groups are aligned as in step 5.
- **Multiplicity:** Holm across the 2 transitions.
- **Role:** this tests stage differences net of patient identity. It is **secondary** and is not part of the tipping-point definition.

## 7. Batch and quality checks (before step 4)

### 7.1 Batch
- **Source:** the GEO metadata contain no scan date or plate. After this addendum is committed (and after PI approval of the download size), the raw CEL archive `GSE117606_RAW.tar` is downloaded **only to read the CEL header fields** (scan date; array/plate barcode if present).
  - It is stored **outside the OneDrive project folder**, in `C:\TDA_data\GSE117606\`.
  - **Intensities are not read.**
  - Only the extracted header table (`data/raw/GSE117606_cel_headers.csv`: GSM, file name, scan date/time, barcode/plate fields) is kept in the project.
- **Measures and criteria** as in addendum 5 part 2: cross-table of plate (or scan day) × stage; Cramér's V; per transition the AUC of scan date and the day/plate overlap. **Confounded** if AUC ≥ 0.80 or ≤ 0.20, or overlap < 20%.
- **Default handling if confounding is found** (stated now):
  1. the frozen analysis proceeds **unchanged**;
  2. the confounded transition is **flagged in all tables**;
  3. if a same-stage split by batch (plate or scan period) exists with **at least 13 arrays on each side**, a **batch-only contrast** is run as in addendum 5 part 2b (frozen step 5 test, 999 permutations, raw and shape-only, unadjusted and adjusted, Holm across the four; m by the pilot pair rule). It is reported next to the size-matched stage contrast.

  Any other handling is a recorded deviation.

### 7.2 Array-quality outlier rule (stated in advance)
On the step 1 gene-level log2 matrix (main set), an array is flagged if **either**:
1. its mean Pearson correlation with all other arrays is below Q1 − 1.5 × IQR of these mean correlations (Tukey rule), **or**
2. its median expression or its interquartile range lies more than 3 robust SD (1.4826 × MAD) from the cohort median.

**Handling:** flagged arrays stay in the **primary** analysis (only identity-based exclusions are made, as in the other cohorts). A **sensitivity analysis** repeats steps 3–5 without them.

## 8. Secondary meta-analysis (informed summary)

- **Normal → Adenoma:** GSE117606 is added to the DerSimonian–Laird random-effects meta-analysis of addendum 3 (shape-only z; within-cohort variance 1), separately for unadjusted and adjusted data. Four cohorts.
- **Adenoma → Carcinoma:** the same model is applied to the shape-only A → C z in all four cohorts (GSE41258, GSE20916 bulk, GSE37364, GSE117606), unadjusted and adjusted.
- **Stage-aware data:** both meta-analyses (N → A and A → C) are also run on the stage-aware z-scores of the four cohorts. GSE41258 enters through exploratory analysis D3 (as in addendum 3).
- **Label:** all are **informed summaries**, because the three earlier cohorts' results are known. None defines a tipping point.

## 9. Expected power

- **What earlier power checks showed:** cohorts of about 50 per group detect only strong effects reliably (D7: raw differences at ≥ 75% of the normal → adenoma benchmark; D8: shape differences of the size observed only with the full sample).
- **Controls:** addendum 5 part 3 showed that localised programmes are removed by preprocessing.
- **What this means here:** GSE117606 is somewhat larger (expected m = 47 vs 19–31), so its power is somewhat higher. Still, **a null result will be inconclusive, not evidence of no change.**
