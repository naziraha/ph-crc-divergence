# Topological analysis of colorectal carcinogenesis in multi-cohort transcriptomics

Analysis code for: *[paper title, authors, journal; to be added]*.

## The study in brief
- **Question.** Does the shape of the gene-expression data change abruptly at one step of colorectal carcinogenesis (normal mucosa → adenoma → carcinoma → metastasis)? Such a "topological tipping point" was looked for with persistent homology.
- **Data.** Four public GEO microarray cohorts (GSE41258, GSE20916, GSE37364, GSE117606; 715 samples) and TCGA-COAD (274 primary tumours, stage I–IV), each analysed separately.
- **Primary test.** A 1-Wasserstein label-permutation test on Vietoris–Rips persistence diagrams, with a tipping-point criterion prespecified in the version-controlled analysis plan (`docs/analysis_prespecification.md`).
- **Result.** The criterion was not met in any cohort.
- **Follow-up analyses** (addenda, written before each was run) examined between-sample heterogeneity, consensus molecular subtypes, genes more variable in carcinoma and chromosome-arm gains. Every departure from the plan is listed in `docs/analysis_deviations.md`.

## Repository layout
| Folder / file | Content |
|---|---|
| `docs/` | Analysis plan (`analysis_prespecification.md`), its 12 addenda, the deviation log (`analysis_deviations.md`) and the commit map (`COMMIT_MAP.md`) |
| `src/` | Core library: preprocessing, Vietoris–Rips persistence (GUDHI), exact 1-Wasserstein distance, permutation test, Holm correction (`tda_crc_pilot.py`), and its validation (`validate_tda_pipeline.py`) |
| `scripts/` | One script per analysis step, addendum, figure and table (run order in `run_all.sh`); `scripts/archive/` holds scripts not used in the paper |
| `data/` | `download_data.py` (all public inputs), `PROVENANCE.md` (file dates and SHA-256 checksums), `raw/GSE117606_cel_headers.csv` (array plate and scan dates, derived from the GEO CEL archive) |
| `results/` | Result tables (`.csv`) and run summaries (`summary.txt`) of every step, including the source data behind each figure (`results/source_data/`, `results/hgnc/source_data/`) |
| `outputs/` | Figures and tables written by the figure and table scripts (created when they are run) |
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
Then put that environment's `Rscript` on the PATH, or set `R_CMS_HOME` to the environment folder. Table 2 does not need R.

## Data
```bash
python data/download_data.py
```
This downloads every input into `data/raw/` (about 250 MB) and compares each file with the checksum in `data/PROVENANCE.md`:
- the GEO series matrices and platform tables;
- NCBI `Homo_sapiens.gene_info`;
- the UCSC Xena TCGA-COAD files.

Two sources are updated in place, so a newer file is expected: NCBI gene_info (daily) and the Xena clinical matrix (occasionally). For these, a checksum difference is reported as a warning with an explanation, and the run continues; results can then differ slightly. Any other difference stops the script, because the input is not the one used for the paper (`--allow-mismatch` overrides). The exact files used for the paper are archived at [Zenodo DOI].

## Reproducing the results
| Command | What it does | Time (8-core laptop) |
|---|---|---|
| `bash run_all.sh table2` | Downloads the data, prepares all cohorts, runs the primary tests and writes Table 1 and Table 2 to `outputs/main_tables/` | about 1 h |
| `bash run_all.sh figures` | Redraws every figure and table from the result files in `results/` | a few minutes |
| `bash run_all.sh all` | The complete pipeline (all steps, sensitivity analyses, controls, power simulations, addenda), then all figures and tables | several hours |

On Windows, run these from Git Bash or another bash shell with the `tda-crc` environment active. Each analysis step saves checkpoints, so an interrupted run resumes where it stopped.

**One command per figure or table** (from the shipped results; outputs in `outputs/`):

| Item | Command |
|---|---|
| Figure 1 | `python scripts/fig1_workflow.py` |
| Figure 2 | `python scripts/fig2_v2_gse41258_diagrams.py` |
| Figure 3 | `python scripts/fig3_gse41258_barcodes_betti_landscapes.py` |
| Figure 4 | `python scripts/fig4_heterogeneity.py adjusted` |
| Figure 5 | `python scripts/fig5_cms_heterogeneity.py` |
| Figure 6 | `python scripts/fig6_variable_genes.py` |
| Figure 7 | `python scripts/fig7_chromosome_arms.py` |
| Tables 1–2, Supplementary Table S_H1_landscape | `python scripts/make_main_tables.py` |
| Supplementary figures | `python scripts/fig2_gse41258_diagrams.py`; `python scripts/fig4_heterogeneity.py unadjusted`; `python scripts/figS_cms_pc12.py` |
| Supplementary tables | `python scripts/make_table_s_<name>.py` for `cohort_differences`, `122_genes`, `power_controls`, `four_arms`, `tcga_stages`, `sensitivity`, `symbol_harmonisation` |

**Hardware and software used:**
- Hardware: Windows 11 laptop, Intel Core Ultra 5 (8 cores), 32 GB RAM.
- Software: Python 3.11.16; R 4.5.3.
- Parallel jobs: permutation tests run 6 worker processes.

## Reproducibility notes
- **Random seeds:** every analysis unit uses its own random generator, seeded from 42 and the unit name, so results do not depend on run order or interruptions.
- **Cross-cohort gene matching:** gene symbols are mapped to current HGNC symbols before cross-cohort comparisons (variable genes, chromosome arms; `scripts/hgnc_crosscohort.py`). Within-cohort analyses use the symbols of each platform annotation.
- **Analysis plan:** the plan with the prespecified tipping-point criterion is commit `66d1c4f` (dated 4 October 2026, 11:35 UTC+8). The first primary results are commit `a149d04` (dated the same day, 11:54). The analysis documents cite commits of the project repository in which the work was done; `docs/COMMIT_MAP.md` maps them to the commits of this repository.

## Citation and licence
Please cite the paper and this code (`CITATION.cff`). Code: MIT licence (`LICENSE`). The input data are public and subject to the terms of GEO, NCBI and TCGA / UCSC Xena.
