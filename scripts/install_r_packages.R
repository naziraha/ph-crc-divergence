# Install the R packages used for CMS classification, at the versions used for the paper.
#   Biobase 2.70.0 (Bioconductor 3.22) and CMScaller 2.0.1 (GitHub Lothelab/CMScaller, commit f1541730).
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("Biobase", update = FALSE, ask = FALSE)
remotes::install_github("Lothelab/CMScaller@f1541730b56c2f37c123037476eed0e6de0ee6f0", upgrade = "never")
cat("Biobase", as.character(packageVersion("Biobase")), "| CMScaller", as.character(packageVersion("CMScaller")), "\n")
