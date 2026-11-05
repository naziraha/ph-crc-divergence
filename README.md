# Topological analysis of colorectal carcinogenesis in multi-cohort transcriptomics

Analysis code for: Amran, N. & Peng, W. K. *Topological Analysis Reveals Distinct Transcriptomic Dynamics Across Colorectal Carcinogenesis* (submitted to *Nature Computational Science*).

## The study in brief
- **Question.** Does colorectal carcinogenesis (normal mucosa → adenoma → carcinoma → metastasis) pass through a shared tipping point, or do tumours diverge gradually? Persistent homology was used to separate the **scale** of the gene-expression data (how far apart samples lie) from its **shape** (how samples are arranged once spread is removed). A shared tipping point should change shape; gradual divergence should change scale alone.
- **Data.** Four public GEO microarray cohorts (GSE41258, GSE20916, GSE37364, GSE117606; 715 samples) and TCGA-COAD (274 primary tumours, stage I–IV), each analysed separately.
- **Primary test.** A 1-Wasserstein label-permutation test on Vietoris–Rips persistence diagrams, on raw (scale) and median-rescaled (shape-only) distances, with a tipping-point criterion prespecified in the version-controlled analysis plan (`docs/analysis_prespecification.md`).
- **Result.** Progression changed the scale of the transcriptome but not its shape; the criterion was not met in any cohort. Carcinomas diverged at the adenoma–carcinoma step, and this divergence then reached a plateau.
- **Follow-up analyses** (addenda, written before each was run) examined between-sample dispersion, consensus molecular subtypes, genes more variable in carcinoma, chromosome-arm gains and gradual versus switch-like change. Every departure from the plan is listed in `docs/analysis_deviations.md`.

## Repository layout
| Folder / file | Content |
|---|---|
| `docs/` | Analysis plan (`analysis_prespecification.md`), its 12 addenda, the deviation log (`analysis_deviations.md`) and the commit map (`COMMIT_MAP.md`) |
| `src/` | Core library: preprocessing, Vietoris–Rips persistence (GUDHI), exact 1-Wasserstein distance, permutation test, Holm correction (`tda_crc_pilot.py`), and its validation (`validate_tda_pipeline.py`) |
| `scripts/` | One script per analysis step, addendum, figure and table (run order in `run_all.sh`); `scripts/archive/` holds scripts not used in the paper |
| `data/` | `download_data.py` (all public inputs), `PROVENANCE.md` (file dates and SHA-256 checksums), `raw/GSE117606_cel_headers.csv` (array plate and scan dates, derived from the GEO CEL archive) |
| `results/` | Result tables (`.csv`) and run summaries (`summary.txt`) of every step, including the source data behind each figure (`results/source_data/`, `results/hgnc/source_data/`) |
| `outputs/` | Figures and tables written by the figure and table scripts (created when they are run) |
| `verification/` | `Verify_Methods_Formulas.ipynb`: independent, read-only recomputation of the Methods formulas from the saved results (deviation D32) |
| `environment.yml`, `environment-r.yml` | Python and R environments with the exact package versions used |
| `run_all.sh` | Reproduces the analysis |

## Installation
```bash
conda env create -f environment.yml
conda activate tda-crc
```
The CMS classification (addenda 9–11) also needs R:

```bash
conda env create -f environment-r.yml
conda activate r-cms
Rscript scripts/install_r_packages.R
```
Then put that environment's `Rscript` on the PATH, or set `R_CMS_HOME` to the environment folder. The primary test (`bash run_all.sh table2`) does not need R.

## Data
```bash
python data/download_data.py
```
This puts every input into `data/raw/` (about 250 MB) and compares each file with the checksum in `data/PROVENANCE.md`:
- the GEO series matrices and platform tables;
- NCBI `Homo_sapiens.gene_info`;
- the UCSC Xena TCGA-COAD files.

Four files are taken by default from the archived copies used for the paper ([Zenodo DOI]): NCBI gene_info, the GPL570 platform table, and the TCGA clinical matrix and MC3 mutation table. Their public versions change over time (gene_info daily; the platform header lists newly submitted samples), and the two TCGA copies used for the paper had been re-saved in a spreadsheet program without changing the fields used (details in `data/PROVENANCE.md`). Until the deposit is published, unpack it locally and run `python data/download_data.py --archive-dir <folder>` (for `run_all.sh`: `DOWNLOAD_ARGS="--archive-dir <folder>"`). `--from-source` uses the current public versions instead; their known differences are then reported as warnings, and results can differ slightly. Any other checksum difference stops the script (`--allow-mismatch` overrides).

## Reproducing the results
| Command | What it does | Time (8-core laptop) |
|---|---|---|
| `bash run_all.sh table2` | Downloads the data, prepares all cohorts, runs the primary tests and writes Supplementary Tables 1 and 2 to `outputs/` | about 1 h 45 min (tested from a fresh clone and environment on 14 October 2026) |
| `bash run_all.sh figures` | Redraws every figure and table of the paper from the result files in `results/` | a few minutes (Fig. 5 and the enrichment step need internet access for OmniPath) |
| `bash run_all.sh all` | The complete pipeline (all steps, sensitivity analyses, controls, power simulations, addenda), then all figures and tables | several hours |

On Windows, run these from Git Bash or another bash shell with the `tda-crc` environment active. Each analysis step saves checkpoints, so an interrupted run resumes where it stopped.

**Paper item → script** (from the shipped results; outputs in `outputs/main_figures`, `outputs/main_tables` and `outputs/supplementary`). Run the table scripts before Fig. 2 and Supplementary Fig. 1, which read two of the tables.

| Paper item | Command | Output file |
|---|---|---|
| Fig. 1 | `python scripts/fig1_ph_scale_shape.py` | `Fig1_PH_scale_shape` |
| Fig. 2 | `python scripts/fig2_dispersion.py adjusted` | `Fig2_dispersion` |
| Fig. 3 | `python scripts/fig3_cms.py` | `Fig3_CMS` |
| Fig. 4 | `python scripts/fig4_variable_genes.py` | `Fig4_variable_genes` |
| Fig. 5 | `python scripts/fig5_chromosome_arms.py` | `Fig5_chromosome_arms` |
| Table 1 | `python scripts/make_table1_gradual_vs_switch.py` | `Table1_gradual_vs_switch` |
| Supplementary Fig. 1 | `python scripts/supplfig1_workflow.py` | `SupplementaryFig1_workflow` |
| Supplementary Fig. 2 | `python scripts/supplfig2_barcodes_betti_landscapes.py` | `SupplementaryFig2_barcodes_Betti_landscapes` |
| Supplementary Fig. 3 | `python scripts/fig2_dispersion.py unadjusted` | `SupplementaryFig3_dispersion_unadjusted` |
| Supplementary Fig. 4 | `python scripts/supplfig4_cms_pc12.py` | `SupplementaryFig4_CMS_PC1_PC2` |
| Supplementary Tables 1 and 2 | `python scripts/make_main_tables.py` | `Table1_cohorts`, `Table2_validation`, `TableS_H1_landscape` (combined and renumbered in the paper) |
| Supplementary Table 3 | `python scripts/make_table_s_meta_analysis.py` | `TableS_meta_analysis` |
| Supplementary Table 4 | `python scripts/make_table_s_sensitivity.py` | `TableS_sensitivity` |
| Supplementary Table 5 | `python scripts/make_table_s_tcga_stages.py` | `TableS_TCGA_stages` |
| Supplementary Table 6 | `python scripts/make_table_s_power_controls.py` | `TableS_power_controls` |
| Supplementary Table 7 | `python scripts/make_table_s_cohort_differences.py` | `TableS_cohort_differences` |
| Supplementary Table 8 | `python scripts/make_table_s_122_genes.py` | `TableS_122_genes` |
| Supplementary Table 9 | `python scripts/make_table_s_symbol_harmonisation.py` | `TableS_symbol_harmonisation` |
| Supplementary Table 10 | `python scripts/make_table_s_secondary_analyses.py` | `TableS_secondary_analyses` |
| Supplementary Table 11 | `python scripts/make_table_s_four_arms.py` | `TableS_four_arms` |

The table scripts write the values and a draft layout; the published tables were renumbered and typeset for the journal without changing any value. The earlier figure scripts (`fig1_workflow.py`, `fig2_v2_gse41258_diagrams.py`, `fig2_gse41258_diagrams.py`, `fig3_gse41258_barcodes_betti_landscapes.py`, `fig4_heterogeneity.py`, `fig5_cms_heterogeneity.py`, `fig6_variable_genes.py`, `fig7_chromosome_arms.py`, `figS_cms_pc12.py`) draw the same results in the pre-submission layout and are kept for the record.

**Four-cohort meta-analysis (deviation D32).** `scripts/meta_analysis_4cohort_revised.py` computes the DerSimonian–Laird estimate explicitly with the between-cohort variance truncated at zero; its results (`results/meta_analysis_4cohort_revised/`) are those reported in the paper. The original script and outputs are kept.

**Hardware and software used:**
- Hardware: Windows 11 laptop, Intel Core Ultra 5 (8 cores), 32 GB RAM.
- Software: Python 3.11.16; R 4.5.3.
- Parallel jobs: permutation tests run 6 worker processes.

## Reproducibility notes
- **Random seeds:** every analysis unit uses its own random generator, seeded from 42 and the unit name, so results do not depend on run order or interruptions.
- **Cross-cohort gene matching:** gene symbols are mapped to current HGNC symbols before cross-cohort comparisons (variable genes, chromosome arms; `scripts/hgnc_crosscohort.py`). Within-cohort analyses use the symbols of each platform annotation.
- **"Pre-registration" in commit messages:** some earlier commit messages and script comments use "pre-registration" loosely. The analysis plan and its addenda were prespecified before analysis in this version-controlled repository; they were not registered in any external registry.
- **Analysis plan:** the plan with the prespecified tipping-point criterion is commit `66d1c4f` (dated 4 October 2026, 11:35 UTC+8). The first primary results are commit `a149d04` (dated the same day, 11:54). The analysis documents cite commits of the project repository in which the work was done; `docs/COMMIT_MAP.md` maps them to the commits of this repository.

## Citation and licence
Please cite the paper and this code (`CITATION.cff`). Code: MIT licence (`LICENSE`). The input data are public and subject to the terms of GEO, NCBI and TCGA / UCSC Xena.
