# Archived scripts (not used in the paper)

These scripts were part of the analysis history and are kept for transparency. None of their outputs appears in the paper's figures or tables, and `run_all.sh` does not run them. Their result files, where they exist, are in `results/`.

| Script | What it did | Why it is not used | Record |
|---|---|---|---|
| `effect_sizes_addendum4.py` | Bootstrap confidence intervals for topological effect sizes (addendum 4). | The intervals were found to be invalid: subsamples contained duplicate samples, and within-stage SDs came from overlapping subsamples. They were replaced by the empirical power analysis (`scripts/power_D7.py`, `scripts/power_D8.py`). | Deviation D7 (`docs/analysis_deviations.md`); addendum 4 |
| `crosscheck_D6_ripser.py` | Software cross-check: the step 4 persistence diagrams of GSE41258 recomputed with Ripser instead of GUDHI. | A one-off verification of the persistence computation; it agreed with GUDHI and is not part of the reported analysis. | Deviation D6 |

**Running them:** each script resolves the repository root from its own location, so it can be run from the repository root, e.g. `python scripts/archive/meta_analysis_NA.py`. `crosscheck_D6_ripser.py` needs the optional package `ripser` (`pip install ripser==0.6.15`); `effect_sizes_addendum4.py` needs only `environment.yml`. The secondary meta-analyses (addenda 3 and 6) were archived here earlier but are reported in the paper (Supplementary Table S_meta_analysis; deviation D29), so their scripts are in `scripts/`.
