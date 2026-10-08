# Addendum 9: What drives the greater heterogeneity of carcinomas compared with adenomas?

**Date:** 2026-10-08 (locked at commit, before any computation)
**Status:** written before any analysis in this addendum has been run. Known so far:
- carcinomas are more dispersed than adenomas: step 13 PERMDISP, stage-blind adjusted data, D = 1.07 (n.s.) in GSE41258, and significant in GSE20916 (1.25), GSE37364 (1.37) and GSE117606 (1.19);
- progression is consistent with gradual change (addendum 8).

This is therefore an **informed** analysis. Its rules are fixed here before any gene-level variability test, subtype call or dispersion decomposition is computed.
**Relation to earlier documents:** the main pre-specification and addenda 1–8 are unchanged.
**Context:** Prof Peng has decided that the paper focuses on biology. This addendum asks a biological question with conventional statistics.
**Rule:** this file is not changed after the computation starts. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Question

Carcinomas are transcriptionally more heterogeneous than adenomas. Is this mainly because carcinomas fall into **different molecular subtypes** ("heterogeneity reflects subtype diversity"), or because carcinomas are more variable **even within a subtype** ("heterogeneity within subtypes")? And which genes and pathways become more variable?

Only the Adenoma → Carcinoma comparison is studied. All wording is "consistent with"; no claim of clonal evolution is made.

## 2. Data

| Cohort | Adenoma / Carcinoma | Platform |
|---|---|---|
| GSE41258 | 48 / 181 (set `main`) | GPL96 |
| GSE20916 | 45 / 36 (set `bulk`) | GPL570 |
| GSE37364 | 29 / 27 (set `bulk`) | GPL570 |
| GSE117606 | 59 / 71 (set `main`) | array (as in addendum 6) |

- **Dispersion analyses and the variability test** use the frozen **stage-blind adjusted** data:
  - dispersion: the 10-PC embeddings (`results/step03/embeddings/<cohort>_<set>_adjusted_g2000_pc10.csv`);
  - variability test: adjusted gene expression (`data/processed/<cohort>_<set>_expr_genes_adjusted.csv`).
- **Subtype classification** uses the **unadjusted** step 1 expression (`data/processed/<cohort>_<set>_expr_genes.csv`, log2). The classifier was built and validated on unadjusted tumour expression, and the composition adjustment removes stromal and immune signal that partly defines CMS1 and CMS4. The resulting labels are then attached to the samples in the adjusted embedding. A sensitivity analysis classifies the adjusted expression (§4.1).
- **Sample counts** are checked against the table before computing; any mismatch stops the run.

## 3. Part (a): genes whose variability differs between carcinoma and adenoma

### 3.1 Test
- **Genes:** all genes in the cohort's adjusted expression matrix, excluding its lowest-25% expressed genes (mean over all samples of the set) and genes with missing values.
- **Test:** Brown–Forsythe test (Levene's test on absolute deviations from the group **median**; scipy `levene(center="median")`), carcinoma vs adenoma.
- **Direction:** variability ratio V = MAD(carcinoma) / MAD(adenoma), where MAD is the median absolute deviation. V > 1 means more variable in carcinoma.
- **Also reported per gene:** the mean difference (carcinoma − adenoma, log2) and the SD ratio.
- **Multiplicity:** Benjamini–Hochberg FDR within each cohort.

### 3.2 Replication
- **More variable in carcinoma:** BH q < 0.05 with V > 1 in at least one cohort (discovery), **and** nominal p < 0.05 with V > 1 in **at least 2 other cohorts** in which the gene is testable (present and not among that cohort's lowest-25% expressed genes).
- **Less variable in carcinoma:** the same rule with V < 1 (reported for completeness).

### 3.3 Descriptive follow-ups for replicated genes
1. **Mean-shift check:** in log2 microarray data, variability can change when the mean changes. The number of replicated genes whose median |mean difference| across the cohorts is < 0.5 log2 is reported as **"variability change without a mean change"**.
2. **Link to subtypes (connects to part b):** for each replicated gene, the share of carcinoma variance explained by CMS (R² of one-way ANOVA among carcinomas; median over cohorts) is reported. The overlap with the CMS classifier's template genes is also reported (hypergeometric test).
3. **Pathway over-representation:**
   - hypergeometric test against MSigDB Hallmark and Reactome (decoupler/OmniPath, as in addendum 8), gene sets of 10–500 genes;
   - separate lists for "more variable" and "less variable", only if a list has ≥ 5 genes;
   - background: genes testable in ≥ 3 cohorts;
   - BH within each list; descriptive.

## 4. Part (b): consensus molecular subtypes and dispersion

### 4.1 Subtype classification (CMS1–4)
- **Method:** **CMScaller** (Eide et al., *Sci Rep* 2017; R package, GitHub `Lothelab/CMScaller`). It uses nearest template prediction (Hoshida 2010) with templates of cancer-cell-intrinsic genes. These templates were chosen to be little affected by stromal and immune content, and the method is published for microarray and RNA-seq tumour data.
  - **It needs R.** R is not installed on this laptop (checked 2026-10-08).
  - Proposed: a **separate** conda environment `r-cms` (r-base plus CMScaller and its dependencies), so the `tda` environment is not touched. **This requires PI approval.**
  - The Python pipeline calls the R environment only to produce a table of sample → CMS label → FDR.
- **Settings:**
  - input: unadjusted log2 expression of **carcinomas only**, gene symbols as row names;
  - `RNAseq = FALSE`;
  - gene-wise centring within each cohort's carcinomas (`rowCenter = TRUE`, the CMScaller default);
  - `nPerm = 1000`, `seed = 42`, `FDR = 0.05`, `doPlot = FALSE`.
- **Template coverage:** per cohort, the number and percentage of CMScaller's template genes (`templates.CMS`) present in the cohort's expression matrix after symbol mapping are reported, overall and per CMS template. Cohorts with < 80% coverage are flagged in all tables (§5.4).
- **Primary labels:** each carcinoma's **nearest template** (prediction regardless of FDR), so that every carcinoma has a label.
- **Sensitivity labels:** confident calls only (CMScaller FDR < 0.05); unconfident samples are excluded from the within-subtype analyses.
- **Second sensitivity:** CMS called on the adjusted expression. Agreement with the primary labels is reported (Cohen's kappa).
- **Fallback, only if the PI declines installing R:** a Python re-implementation of nearest template prediction using CMScaller's published template gene list, validated by reproducing CMScaller's own example predictions (agreement ≥ 95%) before use. This is less established and would be stated as such.
- **Adenomas are not classified.** CMS is defined for carcinomas only.

### 4.2 Reported subtype composition
- **Per cohort:** CMS1–4 proportions (primary and confident-only), compared descriptively with the published distribution (Guinney et al. 2015: CMS1 14%, CMS2 37%, CMS3 13%, CMS4 23%, unclassified 13%).
- **Subtype diversity per cohort:** Shannon entropy of the CMS proportions.

### 4.3 Dispersion with and without subtype (stage-blind adjusted 10-PC space, Euclidean)
**Leave-one-out (LOO) distances:** each sample's distance to its own group's centroid, with the centroid computed **without** that sample. This removes the small-group bias that otherwise makes small groups look less dispersed. Subtype groups with fewer than 3 carcinomas are excluded from within-subtype analyses (number reported).

1. **Total dispersion ratio:** D_total = mean LOO distance of carcinomas to the carcinoma centroid / mean LOO distance of adenomas to the adenoma centroid. Tested with PERMDISP: ANOVA F of the LOO distances, 999 label permutations; p and D reported. This mirrors step 13 with LOO centroids.
2. **Variance explained by subtype among carcinomas:** PERMANOVA of the carcinomas by CMS label (Euclidean pseudo-F, 999 permutations); R²_CMS and p.
3. **Within-subtype dispersion (main test):**
   - each carcinoma's LOO distance to the centroid of **its own CMS group**, compared with each adenoma's LOO distance to the adenoma centroid;
   - D_within = mean (carcinoma within-subtype distance) / mean (adenoma distance);
   - test: PERMDISP on these distances (ANOVA F, 999 permutations of the stage labels of the distances), with p and D_within reported.
4. **Share of the excess dispersion explained by subtype:** f = (D_total − D_within) / (D_total − 1), computed only if D_total > 1.
   - **Random-label calibration:** removing any group centroids lowers dispersion a little. So carcinoma CMS labels are permuted 999 times (same group sizes), and the 95th percentile of f under random labels, f_rand95, is reported.
   - A subtype-specific explanation requires f > f_rand95.
5. **Within the largest subtype (secondary, as suggested by the PI):** for each CMS group with ≥ 10 carcinomas, PERMDISP (LOO distances) of that group vs all adenomas; D_k and p, with Holm across the subtype groups tested in the cohort.

- **Seeds:** `SeedSequence(42, crc32(unit name))` per cohort × analysis.
- **Implementation:** numpy/scipy, reusing the validated step 13 PERMDISP/PERMANOVA functions with LOO distances added.

## 5. Interpretation (fixed in advance)

### 5.1 Which cohorts are informative
A cohort is **informative** for this question only if carcinomas are more dispersed than adenomas: D_total > 1 with PERMDISP p < 0.05 (§4.3.1). In a non-informative cohort there is no excess heterogeneity to explain; its results are reported but not counted.

### 5.2 Per informative cohort
| Call | Rule |
|---|---|
| **Subtype diversity** | R²_CMS PERMANOVA p < 0.05 **and** f ≥ 0.5 **and** f > f_rand95 **and** within-subtype dispersion not significantly higher (PERMDISP p ≥ 0.05 or D_within ≤ 1) |
| **Within-subtype heterogeneity** | D_within > 1 with PERMDISP p < 0.05 **and** f < 0.5 |
| **Both** | D_within > 1 with p < 0.05 **and** f ≥ 0.5 **and** f > f_rand95: subtype explains at least half of the excess, but excess heterogeneity remains within subtypes |
| **Inconclusive** | anything else (for example, D_within not significant but f < 0.5, which is too imprecise to tell) |

### 5.3 Overall conclusion
| Conclusion | Rule |
|---|---|
| **Heterogeneity reflects subtype diversity** | "subtype diversity" in ≥ 2 informative cohorts, and "within-subtype heterogeneity" in none |
| **Heterogeneity within subtypes** | "within-subtype heterogeneity" in ≥ 2 informative cohorts, and "subtype diversity" in none |
| **Both contribute** | "both" in ≥ 2 informative cohorts; **or** "subtype diversity" and "within-subtype heterogeneity" each in ≥ 1 cohort with "both" in at least 1 |
| **Inconclusive** | anything else, or fewer than 2 informative cohorts |

The sensitivity analyses (confident-only CMS labels; CMS from adjusted expression) are reported side by side. The conclusion is described as **robust** only if the confident-only labels give the same overall conclusion.

### 5.4 Known limitations, stated in advance
- **Small carcinoma groups:** GSE37364 has 27 carcinomas and GSE20916 has 36. Split into four subtypes, the groups are small, so "inconclusive" calls are likely there. GSE41258 has the most carcinomas but no significant excess dispersion in adjusted data (D = 1.07), so it is probably non-informative.
- **Different data versions for labels and dispersion:** CMS labels come from unadjusted expression, while dispersion is measured in the adjusted space. CMS4 (mesenchymal) and CMS1 (immune) partly reflect composition, which the adjustment reduces. Subtype differences in the adjusted space may therefore be understated. The sensitivity analysis with adjusted-data CMS labels addresses this partly.
- **Mean–variance coupling:** gene variability in log2 data is partly tied to mean expression (§3.3 item 1).
- **Template coverage:** cohorts in which < 80% of CMScaller's template genes are present after mapping have less reliable CMS calls. This is especially possible for GSE41258 (GPL96, fewer genes). These cohorts are flagged ("low template coverage") in every table that uses CMS labels; the rules of section 5 are applied to them unchanged.
- **Unconfident calls:** CMS labels for unconfident samples are uncertain (about 10–15% expected); hence the confident-only sensitivity analysis.
- **Purity:** tumour purity is not measured. Residual purity differences after composition adjustment could add to carcinoma dispersion.

## 6. Validation on simulated data (before any real data)

### 6.1 Dispersion decomposition
10 dimensions; adenomas n = 50, carcinomas n = 70 split into 4 subtype groups (30/20/12/8); 200 replicates per scenario; 999 permutations.

| Scenario | Expected |
|---|---|
| V1: subtype diversity. Each carcinoma subtype is a cluster with adenoma-like spread (SD 1); subtype centres 2 SD apart | D_total > 1 significant in ≥ 80%; call "subtype diversity" in ≥ 80% |
| V2: within-subtype heterogeneity. Carcinomas one cluster with SD 1.4; random subtype labels | D_within > 1 significant in ≥ 80%; f ≤ f_rand95 in ≥ 90%; call "within-subtype heterogeneity" in ≥ 80% |
| V3: null. Adenomas and carcinomas identical (SD 1), random labels | D_within PERMDISP rejection ≤ 5% + 2.5 MC SE; D_total rejection ≤ 5% + 2.5 MC SE |
| V4: LOO bias check. As V3, but with subtype groups of 4/5/6/55 | D_within PERMDISP rejection ≤ 5% + 2.5 MC SE |

### 6.2 Variability test
2,000 null genes with Gaussian or heavy-tailed (t, 4 df) noise, at the real group sizes of each cohort. The Brown–Forsythe nominal rejection rate must be ≤ 5% + 2.5 MC SE.

### 6.3 CMS classifier
- **Option A (R):** CMScaller runs on its packaged example data (`crcTCGAsubset`) and reproduces the labels supplied with the package in ≥ 90% of samples.
- **Option B (fallback):** the agreement requirement of §4.1.

If a check fails, no real data are analysed; the problem is reported and recorded as a deviation.

## 7. Outputs (`results/addendum9/`)
- `validation_summary.txt`;
- `variability_<cohort>.csv` (all genes: V, SD ratio, mean difference, p, q);
- `variability_replicated.csv` (more and less variable, with mean-shift flag and CMS R²);
- `enrichment_variability.csv`;
- `cms_labels.csv` (cohort, sample, nearest CMS, FDR, confident flag, adjusted-data CMS);
- `cms_proportions.csv`;
- `dispersion_decomposition.csv` (D_total, R²_CMS, D_within, f, f_rand95, per-subtype PERMDISP; primary and sensitivity labels);
- `summary.txt` (section 5 calls and conclusion);
- figures: carcinoma and adenoma samples in PC1–PC2 of the adjusted space coloured by CMS; D_total vs D_within per cohort; variability ratios of replicated genes.
- **Scripts:** `scripts/addendum9_heterogeneity.py`, and `scripts/addendum9_cms.R` (Option A).
