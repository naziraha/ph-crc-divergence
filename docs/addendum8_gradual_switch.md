# Addendum 8: Gradual vs switch-like change along normal → adenoma → carcinoma

**Date:** 2026-10-08 (locked at commit, before any computation)
**Status:** written before any analysis in this addendum has been run. The step 5 TDA results and the step 13 benchmark are known (no tipping point; adenoma → carcinoma spread increase; H0/clustering signal reproduced by conventional methods), so this is an **informed** analysis. Its rules are fixed here, before the progression score, gene fits or any summaries have been computed.
**Relation to earlier documents:** the main pre-specification and addenda 1–7 are unchanged. **The tipping-point decision is unchanged**; this addendum adds a biological question.
**Context:** Prof Peng has decided that the paper focuses on biology, not method. This analysis asks a biological question directly, with conventional statistics.
**Rule:** this file is not changed after the computation starts. Deviations go into `docs/analysis_deviations.md`.

---

## 1. Question

Does the colorectal transcriptome change **gradually** (a continuum of intermediate states and genes that change steadily) or in a **switch-like** way (few intermediate samples, and genes that jump within a narrow part of the progression), across normal → adenoma → carcinoma?

The two transitions, Normal → Adenoma (N → A) and Adenoma → Carcinoma (A → C), are judged separately.

- **Positive evidence is required for either answer.** A non-significant test of discontinuity is not by itself evidence of gradual change; section 5 states what counts.
- **Cross-sectional limitation:** each patient is sampled once per stage, so a "progression score" orders samples, not time. A gap between stages may reflect a rapid biological transition **or** intermediate lesions that the cohort did not sample (for example advanced adenomas or early carcinomas). All conclusions are worded "consistent with"; no claim of clonal evolution is made.

## 2. Data

| Cohort | Samples used (N / A / C) | Expression | Score source |
|---|---|---|---|
| GSE41258 | 53 / 48 / 181 (metastases excluded) | `data/processed/GSE41258_main_expr_genes_adjusted.csv` | `results/step03/embeddings/GSE41258_main_adjusted_g2000_pc10.csv` |
| GSE20916 | 24 / 45 / 36 | `GSE20916_bulk_expr_genes_adjusted.csv` | `GSE20916_bulk_adjusted_g2000_pc10.csv` |
| GSE37364 | 38 / 29 / 27 | `GSE37364_bulk_expr_genes_adjusted.csv` | `GSE37364_bulk_adjusted_g2000_pc10.csv` |
| GSE117606 | 65 / 59 / 71 | `GSE117606_main_expr_genes_adjusted.csv` | `GSE117606_main_adjusted_g2000_pc10.csv` |

- **Data:** stage-blind composition-adjusted data only (primary data version), with the frozen step 3 embeddings and gene lists.
- **GSE41258:** PC1 is taken as frozen. It was computed with the metastases included, and the metastases are simply left out here.
- **Sample counts** are checked against the table before computing; any mismatch stops the run.

## 3. Part (a): continuity of samples along a progression score

### 3.1 Progression score
- **Score:** per sample, s_raw = adjusted PC1 (frozen coordinates).
- **Orientation and scale:** s = (s_raw − median_N) / (median_C − median_N). The normal median is at 0 and the carcinoma median at 1, so the score runs N → A → C and is comparable across cohorts.
- **Ordering check (pre-stated):** PC1 counts as a progression axis in a cohort only if the adenoma median a* lies strictly between 0 and 1.
  - If it does not, that cohort's part (a) and part (b) results are reported as "PC1 not ordered" and excluded from the cohort counts in section 5.
- **Also reported:** the variance explained by PC1, and the stage R² of PC1.
- **Secondary score (sensitivity, all cohorts):** projection of each sample onto the line from the normal centroid to the carcinoma centroid in the 10-PC space, scaled the same way. This score uses stage labels, so it is secondary and does not enter the section 5 decisions.

### 3.2 Unimodality: Hartigan's dip test
- **Samples:** dip test on s for (i) all N + A + C samples; (ii) pooled N + A samples; (iii) pooled A + C samples.
- **p-value:** Monte Carlo calibration against the uniform distribution (Hartigan's least-favourable unimodal reference), 9,999 simulated samples of the same size.
  - This is the standard calibration for the dip test. Unimodality is a property of one distribution, so there are no labels to permute and a label-permutation test does not apply.
- **Multiplicity:** Holm across the three tests within each cohort.
- **Known limitation:** the dip test has low power when the two groups differ much in size (e.g., GSE41258 A 48 vs C 181) or overlap moderately. It is therefore combined with the counts below.
- **Implementation:** Python `diptest` package 0.11.0 (port of Hartigan's algorithm); the dip statistic only is used, with the Monte Carlo p-value above.

### 3.3 Samples between stage medians, and valley ratio
For each transition (earlier stage X, later stage Y), with medians m_X < m_Y on s and w = m_Y − m_X:

- **Between-median counts:** number and percentage of X samples, and of Y samples, with m_X < s < m_Y. Also the AUC of s for X vs Y, describing overlap: 0.5 means full overlap and 1 means no overlap.
- **Valley ratio (VR):**
  - **Middle window:** [m_X + w/3, m_X + 2w/3].
  - **Reference windows:** windows of width w/3 centred on m_X and on m_Y.
  - **VR** = (number of X + Y samples in the middle window) / (mean number in the two reference windows).
  - **Reading:** VR ≈ 1 means the space between the stages is as densely occupied as the stage centres (a continuum); VR ≈ 0 means an empty gap.
- **Uncertainty:** 95% percentile bootstrap CI, with 2,000 resamples of samples within stage. Medians, windows and VR are recomputed in each resample.
- **Overlapping stages:** if AUC < 0.70, the stages are not separated on the score. VR is then not computed, and the transition is classed as "overlapping" (section 5.1).

## 4. Part (b): gene trajectories along the progression score

### 4.0 Cross-fitted progression scores (primary for part b)
**Why:** the all-gene PC1 score is built from the same genes that are then fitted against it, which is circular. Part (b) therefore uses cross-fitting, so that no gene is ever fitted against a score it helped to build.

- **Split:** each cohort's frozen top-2,000 gene list is split at random into two halves of 1,000 genes, H1 and H2 (seed `SeedSequence(42, crc32("<cohort>__genesplit"))`).
- **Half scores:** for each half, the score is adjusted PC1 computed from that half's genes only, with the same preprocessing as step 3 (`tda_crc_pilot.preprocess`):
  - per-gene z-score, then PCA with random_state 0, fitted without labels on the same samples as the frozen step 3 embedding (for GSE41258 including the metastases, which are then left out);
  - the step 3 gene filter is **not** re-applied within the half (the genes were already selected);
  - each half score is oriented and scaled exactly as in 3.1 (normal median 0, carcinoma median 1), and the ordering check is applied to it.
- **Fitting:** genes in H2 are fitted and classified against the H1 score, and genes in H1 against the H2 score.
  - The adenoma median a* used for a gene (transition interval, 4.2) is the one of the score it is fitted against.
  - The bootstrap grid (4.3) is computed separately for each half score.
- **Ordering check for part (b):** a cohort counts for part (b) only if **both** half scores pass the ordering check. Otherwise its part (b) results are reported as "not ordered", excluded from the cohort counts in section 5, and the cohort cannot provide replication.
- **Reported:** the Spearman correlation of each half score with the all-gene PC1 score and with the other half score.
- **All part (b) classifications, FDR and replication (4.1–4.4) use the cross-fitted scores.** Part (a) keeps the all-gene PC1 score.
- **Sensitivity (part b only):** the full part (b) analysis is repeated with the original all-gene PC1 score (section 3.1) for every gene. The following are reported side by side with the primary analysis:
  - class proportions;
  - candidate and replicated switch genes;
  - the overlap of the switch gene lists.

  The sensitivity analysis does not enter section 5.

### 4.1 Genes and fits
- **Genes:** each cohort's frozen top-2,000 gene list (`<embedding>_genes.csv`), using adjusted log2 expression of N, A and C samples, with x = the cross-fitted score for that gene (4.0).
- **Progression association:** F-test of a natural cubic regression spline (3 df; patsy `cr`) vs an intercept-only model. Benjamini–Hochberg (BH) FDR within the cohort; q < 0.05 = **progression-associated**. The spline is shape-neutral, so this step favours neither linear nor sigmoid genes.
- **Linear model:** y = α + βs (OLS, 2 parameters).
- **Sigmoid model (4-parameter logistic):** y = a + (b − a) / (1 + exp(−k(s − x₀))). Least squares (scipy `curve_fit`) from multiple starting values, keeping the best fit:
  - x₀ ∈ {0, a*, 0.5, 1};
  - k ∈ {±2, ±5, ±10, ±20, ±50};
  - bounds: |k| ≤ 200, x₀ within the observed range of s.
- **Model comparison:** AIC = n·ln(RSS/n) + 2p (p = 2 linear, 4 sigmoid); ΔAIC = AIC_linear − AIC_sigmoid (positive favours the sigmoid).
- **Derived quantities:**
  - **transition width** W = 2·ln(9)/|k| = 4.39/|k|: the stretch of the score over which the gene completes 10% → 90% of its change, in units of the N → C distance;
  - **amplitude** |b − a| (log2);
  - **midpoint** x₀.
- **Non-monotonic check:** the 3-df spline (4 parameters) is also compared. A gene is **non-monotonic** (e.g., adenoma-transient) if the spline AIC is ≥ 10 below both the linear and the sigmoid AIC **and** the fitted spline changes direction within [0, 1].

### 4.2 Classification of progression-associated genes (pre-stated)

**Switch-like:** all of the following:
1. ΔAIC ≥ 10 (strong preference for the sigmoid);
2. W ≤ 0.25: 80% of the change happens within a quarter of the normal-to-carcinoma distance;
3. 0 ≤ x₀ ≤ 1 (the jump lies between the normal and carcinoma medians);
4. amplitude ≥ 1 log2 (two-fold);
5. **resolvable:** at least 5% of the cohort's samples (and at least 5 samples) lie inside the 10–90% window [x₀ − W/2, x₀ + W/2];
6. not non-monotonic.

**Step-like but unresolvable:** meets 1–4 and 6 but fails 5. The jump falls where there are almost no samples, so switch and gradual change cannot be told apart. These genes are reported separately and are **not** counted as switch-like.

**Gradual:** progression-associated, not non-monotonic, and either ΔAIC < 2 (the sigmoid adds nothing over a line) or W ≥ 1 (the change is spread over at least the whole N → C distance).

**Intermediate:** any other progression-associated gene.

**Transition interval of a switch:** N → A if x₀ ≤ a*; A → C if x₀ > a*.

### 4.3 FDR control for switch-likeness
- **Null hypothesis:** each gene follows its linear trend with Gaussian noise (gradual change).
- **Test:** a parametric bootstrap of ΔAIC that uses the cohort's actual score values, so any gaps between stages are built into the null.
- **Grid shortcut:** under this null, ΔAIC depends only on the standardised slope β/σ. So the null distribution is simulated once per cohort on a grid of 21 values of β/σ (0 to the largest observed value, evenly spaced on the square-root scale), with 1,999 simulations each.
- **Gene p-values:** each gene's p-value is the **larger** of the p-values at the two neighbouring grid points (conservative). The grid approximation is checked against a direct bootstrap (1,999 simulations) for 20 randomly chosen genes; the absolute difference in p must be ≤ 0.02.
- **BH FDR** is applied across all progression-associated genes within the cohort (both halves pooled).
- **Candidate switch gene (in a cohort):** switch-like (4.2) **and** bootstrap BH q < 0.05.

### 4.4 Replication across cohorts
A **switch gene** is a candidate switch gene in at least one cohort that **replicates in at least 2 other cohorts**. Replication in a cohort requires:
- the gene is present in that cohort's adjusted data and not among its lowest-25% expressed genes (it need not be in that cohort's top-2,000);
- **score used:** if the gene is in that cohort's top-2,000, the half score that excludes it; if not, the H1 score of that cohort (the gene contributed to neither half score);
- it meets switch-like criteria 1–6 there, with nominal bootstrap p < 0.05 instead of FDR, and the spline association p < 0.05 nominal;
- the same direction (sign of b − a);
- the same transition interval (N → A or A → C).

Cohorts that fail the ordering check (3.1) cannot provide replication.

### 4.5 Reported per cohort
- **Proportions:** of progression-associated genes, the percentage switch-like, candidate switch (FDR), step-like-unresolvable, gradual, intermediate and non-monotonic; switch-like genes are also split by transition interval and direction.
- **Gene list:** all candidate and replicated switch genes, with x₀, W, amplitude, ΔAIC, p and q per cohort.

### 4.5a Replicated non-monotonic genes (descriptive; not part of section 5)
- **Discovery:** non-monotonic genes (4.1) that are progression-associated (BH q < 0.05) in a cohort, using the cross-fitted scores.
- **Pattern of each gene**, from the fitted spline within [0, 1]:
  - **peak** (rises, then falls) or **trough** (falls, then rises);
  - the location of the turning point relative to the adenoma median a*;
  - the adenoma-minus-normal and adenoma-minus-carcinoma mean differences (log2).
- **Adenoma-transient:** a gene whose adenoma mean differs from **both** the normal and the carcinoma means in the same direction (up in adenoma only, or down in adenoma only), each by ≥ 0.5 log2.
- **Replication:** the same rule as 4.4. The gene must be non-monotonic with spline association p < 0.05 (nominal) in **at least 2 other** cohorts, scored with the same score choice as 4.4, with the **same pattern** (peak vs trough, and the same adenoma-transient direction where applicable).
- **Reported:** `nonmonotonic_genes_replicated.csv`, with pattern, turning point and mean differences per cohort. Descriptive only; included in the pathway enrichment (4.6) as separate lists (peak / trough).

### 4.6 Pathway enrichment of replicated switch genes (descriptive)
- **Test:** over-representation (hypergeometric test, scipy).
- **Gene sets:** MSigDB Hallmark and Reactome, obtained through decoupler/OmniPath (already installed); sets of 10–500 genes.
- **Gene lists:** separate lists by transition interval × direction, only if the list has ≥ 5 genes.
- **Background:** genes present in at least one cohort's top-2,000 list and testable (present and expressed) in at least 3 cohorts.
- **Multiplicity:** BH FDR within each list.
- Descriptive; no conclusion in section 5 depends on it.

## 5. Interpretation (fixed in advance)

### 5.1 Population level (part a), per cohort and transition
| Call | Rule |
|---|---|
| **Continuous** | VR ≥ 0.5 **and** VR lower 95% CI ≥ 0.25 **and** dip test for the pooled pair not significant (Holm p ≥ 0.05) |
| **Overlapping** | AUC < 0.70 (stages not separated on the score); counted together with "continuous" as compatible with gradual change |
| **Discontinuous** | dip test for the pooled pair significant (Holm p < 0.05) **and** VR upper 95% CI < 0.25 |
| **Inconclusive** | anything else |

### 5.2 Gene level (part b), per cohort and transition
Here P = number of progression-associated genes in that cohort.

| Call | Rule |
|---|---|
| **Switch-like** | candidate switch genes (FDR) with midpoint in this interval ≥ 10% of P |
| **Gradual** | gradual genes ≥ 50% of P **and** candidate switch genes in this interval < 2% of P |
| **Inconclusive** | anything else |

### 5.3 Overall conclusion per transition (N → A; A → C)
Only cohorts that pass the ordering check (3.1) count. If fewer than 3 cohorts pass, the conclusion is **inconclusive**.

| Conclusion | Rule |
|---|---|
| **Support for gradual change** | population level continuous or overlapping in ≥ 3 cohorts **and** gene level gradual in ≥ 3 cohorts |
| **Support for switch-like change** | population level discontinuous in ≥ 3 cohorts **and** ≥ 10 replicated switch genes (4.4) for this transition |
| **Gradual population with a switch-like gene programme** | the gradual rule is met **and** ≥ 10 replicated switch genes exist for this transition |
| **Sampling gap without gene-level switch** | population level discontinuous in ≥ 3 cohorts, but < 10 replicated switch genes. Reported as **inconclusive**: consistent with unsampled intermediate lesions rather than a demonstrated switch |
| **Inconclusive** | anything else |

All conclusions are worded "consistent with gradual / switch-like change". The thresholds (VR 0.5/0.25, AUC 0.70, 10%/2%/50%, ≥ 10 genes) are judgement calls fixed here before seeing any result.

### 5.4 Known limitations, stated in advance
- **Cross-sectional data:** the score orders samples, it does not measure time (section 1).
- **Score circularity** is removed for part (b) by cross-fitting (4.0): no gene is fitted against a score built from it. Two limits remain:
  - a half score is built from only 1,000 genes, so it may be a little noisier than the all-gene score;
  - genes in the two halves are correlated (co-expression modules), so cross-fitting removes direct but not module-level circularity.

  The all-gene sensitivity analysis shows how much the choice matters. Part (a) uses the all-gene score, where circularity does not arise because no genes are fitted.
- **Batch:** in GSE41258, normal and adenoma arrays were scanned on different days (addendum 5), so an N → A gap in that cohort may partly reflect batch.
- **One axis:** a single progression axis can miss changes orthogonal to PC1.
- **Unadjusted data and stage-aware data** are not analysed here.

## 6. Validation on simulated data (before any real data)

### 6.1 Simulated cohort
- **Samples and stages:** n = 150 samples with a latent progression t, in two layouts:
  - **continuum:** t uniform on [−0.3, 1.3]; stages assigned by t (N: t < 0.25; A: 0.25–0.75; C: > 0.75);
  - **three clusters:** N/A/C clusters at t = 0, 0.5, 1 with SD 0.08, 50 samples each.
- **Expression matrix:** 2,000 genes:
  - 1,000 background genes linear in t (slope ~ N(0, 1)) plus noise (SD 1); these carry the progression axis;
  - 960 pure-noise genes;
  - 40 test genes (20 of the scenario's type and 20 linear test genes).
- **Pipeline:** the genes are randomly split into halves, and **cross-fitted half scores** are built exactly as in 4.0 (z-score, PCA, PC1, orientation and scaling by the simulated stage medians). The test genes are classified against the half score that excludes them.
- **Population checks (dip, VR)** use the all-gene PC1 score, as part (a) does.

### 6.2 Replicates and computing
- **Replicates:** 200 per scenario for the population checks; 50 per scenario for the gene-classifier checks.
- **Computational simplification in the validation only:**
  - the bootstrap grid uses 499 simulations per grid point;
  - the sigmoid fits use 16 starting values (k ∈ {±5, ±20} × the four x₀ starts);
  - BH is applied across the 40 test genes of each replicate.

  The real-data analysis uses the full settings of section 4.

### 6.3 Acceptance criteria
| Check | Acceptance |
|---|---|
| Half scores recover t | median Spearman(half score, t) ≥ 0.9 in both layouts |
| Dip test: unimodal (continuum) | rejection ≤ 5% + 2.5 MC SE |
| Dip test: three clusters | rejection ≥ 80% |
| VR: continuum | median VR ≥ 0.5 |
| VR: three clusters | median VR upper CI < 0.25 |
| Linear test genes (amplitude 2 log2 over t ∈ [0, 1], noise SD 0.5), cross-fitted, continuum and clusters | candidate switch (FDR) rate ≤ 5% |
| Sigmoid test genes (W = 0.1, x₀ = 0.25, amplitude 2, SD 0.5), cross-fitted, continuum | candidate switch power ≥ 80% |
| Sigmoid test genes in three clusters (x₀ = 0.25, jump inside the N–A gap), cross-fitted | classified "step-like unresolvable", not switch-like, in ≥ 80% |
| Peak test genes (adenoma-transient: +2 log2 at t = 0.5, Gaussian bump SD 0.12), cross-fitted, continuum | classified non-monotonic in ≥ 80% |
| Grid vs direct bootstrap p (section 4.3), on one real-size simulated replicate | absolute difference ≤ 0.02 |

If a check fails, no real data are analysed; the problem is reported and recorded as a deviation.

## 7. Settings, seeds and outputs
- **Seeds:** `SeedSequence(42, crc32(unit name))` per cohort × analysis.
- **Packages:** numpy, scipy, pandas, patsy, statsmodels, decoupler; `diptest` 0.11.0 (installed 2026-10-08 with PI approval; core packages unchanged).
- **Script:** `scripts/addendum8_gradual_switch.py`; the simulation check (section 6) runs first.
- **Outputs (`results/addendum8/`):**
  - `validation_summary.txt`;
  - `progression_scores.csv`;
  - `continuity_results.csv` (dip, between-median counts, AUC, VR with CI);
  - `gene_split_<cohort>.csv` (half assignment);
  - `half_scores.csv` (with correlations);
  - `gene_fits_<cohort>.csv` (cross-fitted; all fits and classes);
  - `gene_fits_<cohort>_allgene.csv` (sensitivity);
  - `switch_genes_candidates.csv`;
  - `switch_genes_replicated.csv` (primary, plus the sensitivity version);
  - `nonmonotonic_genes_replicated.csv`;
  - `enrichment_replicated.csv`;
  - `summary.txt` (section 5 calls and conclusions);
  - figures: score densities by stage per cohort; example trajectories of replicated switch genes; class proportions per cohort.

---

## Amendment 1 (2026-10-08): fixed-window resolvability rule (reference: deviation D14)

**Status:** decided by the PI and committed **before** the section 6 validation is rerun. No real data have been analysed.

**Reason:** in the first section 6 validation (commit cca77b1), true switch genes on a continuum were detected in only 61.8% of cases (required ≥ 80%). The main cause was the resolvability criterion (§4.2 item 5). It counted samples inside the **fitted** 10–90% window, which becomes very narrow whenever noise makes the fitted sigmoid too steep. As a result, 29% of true switches were labelled "unresolvable" although samples densely covered the switch region.

**Change:** §4.2 item 5 is replaced by:
> 5. **resolvable:** at least 10% of the cohort's samples (and at least 5 samples) lie within the **fixed** window x₀ ± 0.125 (half the maximum switch width W = 0.25) on the score the gene is fitted against.

The rule applies everywhere item 5 is used: the step-like-unresolvable class (§4.2), replication (§4.4) and the validation (§6).

**Unchanged:** all other criteria and thresholds, including W ≤ 0.25, the 80% power target and all section 5 rules.

**Pre-stated consequence:** if any section 6 criterion still fails after the rerun, the rule will **not** be amended again. The classifier will be reported as conservative, and the PI will be asked before any real data are used.

---

## Amendment 2 (2026-10-08): conditional interpretation of switch-gene calls (reference: deviation D15)

**Status:** decided by the PI and committed **before any real data are analysed**.

**Validation result that motivates it (rerun with Amendment 1, commit 1f15bf9):**
- **On a continuum** the classifier works as intended: true switches detected 84.0%, and 0% false switch calls from linear (gradual) genes.
- **In gapped data** (stages forming separate clusters), jumps placed inside a gap were classified "step-like unresolvable" in 79.3% of cases (required ≥ 80%), and **13.5% were called candidate switch genes**. The midpoint of a jump inside a gap is poorly located; when it falls near a cluster edge, the fixed resolvability window reaches into the neighbouring cluster.
- All other checks passed, including all part (a) checks and the non-monotonic check.

**Change:**
- **Where switch calls count:** switch-gene calls, candidate (§4.3) and replicated (§4.4), for a transition are interpreted **only in cohorts where part (a) classes that transition as "continuous" or "overlapping"** (§5.1).
- **Elsewhere:** in cohorts where part (a) classes the transition as "discontinuous" or "inconclusive", switch-gene calls for that transition are reported as **"not interpretable (possible gap artefact)"**. In those cohorts they:
  - do not count as discovery or as replication (§4.4);
  - do not count towards the gene-level switch-like call (§5.2) or the section 5.3 switch rules.
- **Unchanged:** the part (a) rules, the "gradual" gene rules, the non-monotonic analysis (§4.5a) and all thresholds.

**Note (logical consequence, stated in advance):** if a transition is "discontinuous" or "inconclusive" in 2 or more of the 4 cohorts, at most 2 cohorts can contribute switch genes, so a gene cannot replicate in 2 other cohorts. In that case **no replicated switch genes are possible for that transition**, and the §5.3 conclusion "support for switch-like change" cannot be reached. Such a transition will be reported as **inconclusive** about gene-level switches.
