# Archived scripts (not used in the paper)

These scripts were part of the analysis history and are kept for transparency. None of their outputs appears in the paper's figures or tables, and `run_all.sh` does not run them. Their result files, where they exist, are in `results/`.

| Script | What it did | Why it is not used | Record |
|---|---|---|---|
| `effect_sizes_addendum4.py` | Bootstrap confidence intervals for topological effect sizes (addendum 4). | The intervals were found to be invalid: subsamples contained duplicate samples, and within-stage SDs came from overlapping subsamples. They were replaced by the empirical power analysis (`scripts/power_D7.py`, `scripts/power_D8.py`). | Deviation D7 (`docs/analysis_deviations.md`); addendum 4 |
| `meta_analysis_NA.py` | Random-effects (DerSimonian–Laird) meta-analysis of the normal → adenoma step 5 effect across the three cohorts then available (addendum 3, secondary). | Secondary analysis; the paper reports each cohort separately and does not use the pooled estimate. The omission of this secondary analysis from the paper is not recorded as a deviation in the log. | Addendum 3 (`docs/analysis_prespecification_addendum3_meta.md`) |
| `meta_analysis_4cohort.py` | The same meta-analysis with GSE117606 added (addendum 6 §8), for normal → adenoma and adenoma → carcinoma. | As above. | Addenda 3 and 6 |
| `crosscheck_D6_ripser.py` | Software cross-check: the step 4 persistence diagrams of GSE41258 recomputed with Ripser instead of GUDHI. | A one-off verification of the persistence computation; it agreed with GUDHI and is not part of the reported analysis. | Deviation D6 |

**Running them:** each script resolves the repository root from its own location, so it can be run from the repository root, e.g. `python scripts/archive/meta_analysis_NA.py`. `crosscheck_D6_ripser.py` needs the optional package `ripser` (`pip install ripser==0.6.15`); the other three need only `environment.yml`.
