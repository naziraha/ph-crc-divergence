#!/usr/bin/env bash
# Reproduce the analysis from the public data.
#
#   bash run_all.sh table2    # data download + the steps behind Table 2 (tested from scratch; ~1 h on 6 cores)
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

download() {           # DOWNLOAD_ARGS=--allow-mismatch continues if a downloaded file differs from data/PROVENANCE.md
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

table2() {             # primary step 5 tests, H0/H1-only and landscape tests, then Table 1 and Table 2
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

figures() {
  run scripts/source_data.py
  run scripts/hgnc_crosscohort.py
  run scripts/fig1_workflow.py
  run scripts/fig2_v2_gse41258_diagrams.py
  run scripts/fig2_gse41258_diagrams.py
  run scripts/fig3_gse41258_barcodes_betti_landscapes.py
  run scripts/fig4_heterogeneity.py adjusted
  run scripts/fig4_heterogeneity.py unadjusted
  run scripts/fig4_neighbour_check.py
  run scripts/fig4_composition_check.py
  run scripts/fig5_cms_heterogeneity.py
  run scripts/figS_cms_pc12.py
  run scripts/fig6_variable_genes.py
  run scripts/fig7_chromosome_arms.py
  run scripts/make_main_tables.py
  for t in cohort_differences 122_genes power_controls four_arms tcga_stages sensitivity symbol_harmonisation; do
    run "scripts/make_table_s_${t}.py"
  done
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
  run scripts/check_methods_items.py all
  figures
}

case "${1:-}" in
  download|prepare|table2|figures|all) "$1" ;;
  *) echo "usage: bash run_all.sh table2 | all | figures | download | prepare"; exit 2 ;;
esac
