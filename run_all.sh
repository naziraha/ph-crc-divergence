#!/usr/bin/env bash
# Reproduce the analysis from the public data.
#
#   bash run_all.sh table2    # data download + the primary test (Supplementary Table 1; tested from scratch; ~1 h 45 min)
#   bash run_all.sh all       # the complete pipeline, all figures and tables (several hours; needs the R environment)
#   bash run_all.sh figures   # figures and tables only, from the result files shipped in results/
#
# Python: the conda environment in environment.yml (set PYTHON=/path/to/python to override).
# R (CMS classification, addenda 9-11 only): environment-r.yml + scripts/install_r_packages.R; Rscript must be on the
# PATH or R_CMS_HOME must point to the R environment.
# Every analysis step writes resumable checkpoints; re-running a stage skips finished units.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python}"
run() { echo; echo ">>> $PY $*"; "$PY" "$@"; }

download() {           # DOWNLOAD_ARGS, e.g. "--archive-dir /path/to/deposit", "--from-source", "--allow-mismatch"
  run data/download_data.py ${DOWNLOAD_ARGS:-}
}

prepare() {            # steps 1-3: expression tables, tissue composition, embeddings (all cohorts)
  run scripts/step01_prepare.py
  run scripts/step01_prepare_validation.py GSE20916
  run scripts/step01_prepare_validation.py GSE37364
  run scripts/step01_prepare_GSE117606.py
  run scripts/addendum10_tcga_prepare.py
  run scripts/step02_composition.py
  for c in GSE20916 GSE37364 GSE117606 TCGA_COAD; do run scripts/step02_composition_validation.py "$c"; done
  run scripts/step03_embedding.py
  run scripts/step03_embedding_validation.py GSE20916
  run scripts/step03_embedding_validation.py GSE37364
  run scripts/step03_embedding_GSE117606.py
  run scripts/step03_embedding_TCGA.py
}

table2() {             # primary step 5 tests, H0/H1-only and landscape tests, then Supplementary Tables 1 and 2
  download
  prepare
  run scripts/step05_permutation_test.py primary
  for c in GSE20916 GSE37364 GSE117606; do run scripts/step05_permutation_validation.py "$c" primary; done
  run scripts/addendum5_part4_dims.py
  run scripts/addendum10_tcga_extra.py spread
  run scripts/addendum12_ph_outputs.py run
  run scripts/addendum12_ph_outputs.py summary
  run scripts/source_data.py validation
  run scripts/make_main_tables.py
}

figures() {            # every figure and table of the paper (NCS numbering), from the result files in results/
  run scripts/source_data.py
  run scripts/hgnc_crosscohort.py
  run scripts/hgnc_enrichment.py
  # tables first: Fig. 2 and Supplementary Fig. 1 read two of them
  run scripts/make_main_tables.py                      # Supplementary Tables 1 (cohorts, primary test) and 2 (H1 landscape)
  run scripts/make_table1_gradual_vs_switch.py         # Table 1
  for t in meta_analysis sensitivity tcga_stages power_controls cohort_differences 122_genes symbol_harmonisation secondary_analyses four_arms; do
    run "scripts/make_table_s_${t}.py"                 # Supplementary Tables 3-11
  done
  # main figures
  run scripts/fig1_ph_scale_shape.py
  run scripts/fig2_dispersion.py adjusted
  run scripts/fig3_cms.py
  run scripts/fig4_variable_genes.py
  run scripts/fig5_chromosome_arms.py
  # supplementary figures
  run scripts/supplfig1_workflow.py
  run scripts/supplfig2_barcodes_betti_landscapes.py
  run scripts/fig2_dispersion.py unadjusted          # Supplementary Fig. 3
  run scripts/supplfig4_cms_pc12.py
}

all() {
  download
  run scripts/fetch_geo_metadata.py GSE117606
  prepare
  run scripts/exploratory_D3_stageaware_GSE41258.py
  # step 4 (descriptive) and step 5 (primary and sensitivity)
  run scripts/step04_persistence.py
  for c in GSE20916 GSE37364 GSE117606 TCGA_COAD; do run scripts/step04_persistence_validation.py "$c"; done
  run scripts/step05_permutation_test.py primary
  run scripts/step05_permutation_test.py sensitivity
  for c in GSE20916 GSE37364; do
    run scripts/step05_permutation_validation.py "$c" primary
    run scripts/step05_permutation_validation.py "$c" sensitivity
  done
  run scripts/step05_permutation_validation.py GSE117606 primary
  run scripts/step05_secondary_GSE117606.py run
  run scripts/step05_permutation_TCGA.py run
  run scripts/gse117606_batch_check.py
  # addendum 5: batch checks, controls, H0/H1-only tests
  run scripts/addendum5_part2_scandates.py as_specified
  run scripts/addendum5_part2_scandates.py corrected
  run scripts/addendum5_part2b_batch.py
  run scripts/addendum5_part3_controls.py
  run scripts/addendum5_part3b_diagnostic.py
  run scripts/addendum5_part4_dims.py
  # power (D7, D8), benchmark (step 13), addenda 8-12
  run scripts/power_D7.py run
  run scripts/power_D8.py run
  run scripts/step13_benchmark.py validate
  run scripts/step13_benchmark.py run
  run scripts/addendum8_gradual_switch.py validate
  run scripts/addendum8_gradual_switch.py run
  run scripts/addendum9_heterogeneity.py validate
  run scripts/addendum9_heterogeneity.py run
  run scripts/addendum10_tcga_extra.py spread
  run scripts/addendum10_tcga_extra.py cms
  run scripts/addendum10_tcga_extra.py normal
  run scripts/addendum10_tcga_controls.py
  run scripts/addendum11_chr20.py
  run scripts/addendum11_differentiation.py validate
  run scripts/addendum11_differentiation.py run
  run scripts/addendum12_ph_outputs.py validate
  run scripts/addendum12_ph_outputs.py run
  run scripts/addendum12_ph_outputs.py summary
  run scripts/meta_analysis_NA.py
  run scripts/meta_analysis_4cohort.py
  run scripts/meta_analysis_4cohort_revised.py     # D32: DerSimonian-Laird with tau^2 truncated at zero
  run scripts/check_methods_items.py all
  figures
}

case "${1:-}" in
  download|prepare|table2|figures|all) "$1" ;;
  *) echo "usage: bash run_all.sh table2 | all | figures | download | prepare"; exit 2 ;;
esac
