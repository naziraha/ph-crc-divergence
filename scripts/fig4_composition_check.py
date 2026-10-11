"""
Descriptive check for the Figure 4 reading guide (no statistical tests): median epithelial, stromal and immune
composition scores (step 2, decoupler ULM scores) of adenomas and carcinomas per GEO cohort, their difference, and a
standardised difference (median difference / pooled MAD-based SD of adenomas + carcinomas in that cohort), because the
raw score scale differs between platforms; also the within-stage spread of each score (robust SD, MAD x 1.4826)
in adenomas and carcinomas and its ratio. Primary sample sets only (samples in the frozen embeddings).
Output: results/source_data/fig4_composition_medians.csv (new file).

Usage: python scripts/fig4_composition_check.py
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
FILES = {"GSE41258": ("main", "step02/composition_scores_main.csv"),
         "GSE20916": ("bulk", "step02/GSE20916/composition_scores_bulk.csv"),
         "GSE37364": ("bulk", "step02/GSE37364/composition_scores_bulk.csv"),
         "GSE117606": ("main", "step02/GSE117606/composition_scores_all.csv")}
COMP = ["Epithelial", "Stromal", "Immune"]


def main():
    rows = []
    for c, (st, f) in FILES.items():
        emb = pd.read_csv(os.path.join(R, "step03", "embeddings", f"{c}_{st}_adjusted_g2000_pc10.csv"), index_col=0)
        sc = pd.read_csv(os.path.join(R, f), index_col=0).loc[emb.index]
        stage = emb.stage
        for k in COMP:
            a, b = sc.loc[stage == "Adenoma", k], sc.loc[stage == "Carcinoma", k]
            both = pd.concat([a, b])
            sd = median_abs_deviation(both, scale="normal")
            sa, sb = median_abs_deviation(a, scale="normal"), median_abs_deviation(b, scale="normal")
            rows.append(dict(cohort=c, score=k, n_adenoma=len(a), n_carcinoma=len(b), median_adenoma=a.median(),
                             median_carcinoma=b.median(), diff_C_minus_A=b.median() - a.median(),
                             robust_sd_AC=sd, std_diff=(b.median() - a.median()) / sd,
                             robust_sd_adenoma=sa, robust_sd_carcinoma=sb, sd_ratio_C_over_A=sb / sa))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(R, "source_data", "fig4_composition_medians.csv"), index=False)
    pd.set_option("display.width", 200)
    print(out.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
