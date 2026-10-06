NOTE (2026-10-06): The bootstrap confidence intervals, the plain-language categories and the pooled
estimates in this folder are INVALID and must not be interpreted.

Why: resampling samples with replacement puts duplicates into each stage. Each subsample then has only
about 70-80% distinct samples, and subsamples overlap heavily, so diagrams of the same stage look
artificially similar. The effect size d divides by the standard deviation of within-stage distances
between these overlapping subsamples, which shrinks further. The bootstrap therefore inflates d, most
strongly for transitions without a real difference (point estimate below the whole interval in 11 of
13 non-significant results). The pooled estimates use the bootstrap standard errors and are invalid for
the same reason.

Valid: the point estimates of d (computed on the original samples), reported as descriptive values only.
d is not on Cohen's scale, so d = 0.5 is not a meaningful threshold here.
Power and the classification of null results: deviation D7 (results/power_D7/).
