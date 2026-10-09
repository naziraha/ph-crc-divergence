# Addendum 9 (5a9b03a; D16 9590f65) - CMS classification with CMScaller (r-cms environment).
# Usage:
#   Rscript scripts/addendum9_cms.R example <out.csv>                 - section 6.3 validation (package example data)
#   Rscript scripts/addendum9_cms.R classify <in.csv> <out.csv>       - carcinomas: genes (symbols) x samples, log2;
#                                                                       centred per gene here (D16a), RNAseq = FALSE
suppressMessages({library(CMScaller); library(Biobase)})
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

nearest <- function(res) {
  d <- as.matrix(res[, grep("^d\\.CMS", colnames(res))])
  sub("^d\\.", "", colnames(d))[apply(d, 1, which.min)]
}

if (mode == "example") {
  data(crcTCGAsubset)
  res <- CMScaller(exprs(crcTCGAsubset), rowNames = "entrez", RNAseq = TRUE, nPerm = 1000, seed = 42,
                   FDR = 0.05, doPlot = FALSE, verbose = FALSE)
  out <- data.frame(sample = rownames(res), nearest = nearest(res), prediction = as.character(res$prediction),
                    FDR = res$FDR, gold = as.character(crcTCGAsubset$CMS.Syn), stringsAsFactors = FALSE)
  write.csv(out, args[2], row.names = FALSE)
} else if (mode == "classify") {
  x <- as.matrix(read.csv(args[2], row.names = 1, check.names = FALSE))
  x <- x - rowMeans(x)                         # D16a: gene-wise centring over the cohort's carcinomas
  cov <- sapply(levels(templates.CMS$class), function(k)
    mean(templates.CMS$symbol[templates.CMS$class == k] %in% rownames(x)))
  res <- CMScaller(x, rowNames = "symbol", RNAseq = FALSE, nPerm = 1000, seed = 42, FDR = 0.05,
                   doPlot = FALSE, verbose = FALSE)
  out <- data.frame(sample = rownames(res), nearest = nearest(res), prediction = as.character(res$prediction),
                    FDR = res$FDR, stringsAsFactors = FALSE)
  write.csv(out, args[3], row.names = FALSE)
  write.csv(templates.CMS, file.path(dirname(args[3]), "templates.csv"), row.names = FALSE)
  write.csv(data.frame(template = names(cov), coverage = cov,
                       overall = mean(templates.CMS$symbol %in% rownames(x))),
            sub("\\.csv$", "_coverage.csv", args[3]), row.names = FALSE)
}
cat("CMScaller", as.character(packageVersion("CMScaller")), "R", R.version.string, "\n")
