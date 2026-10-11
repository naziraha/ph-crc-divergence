"""
Descriptive check for the Figure 4 reading guide (no statistical tests): per cohort x stage, the median nearest-neighbour
distance and the median distance to the stage centroid, in the frozen 10-PC embeddings.
Nearest-neighbour (NN) distance depends on group size, so it is given (i) on all samples of the stage and (ii) within
the 20 equal-size subsamples used by step 5 / the PH figures (raw-mode members, results/source_data/
ph_subsample_members_GEO.csv; median over all samples of all 20 subsamples). Centroid distance: plain and
bias-adjusted (results/source_data/stage_centroid_distances_GEO.csv).
Output: results/source_data/fig4_neighbour_centre_medians.csv (new file).

Usage: python scripts/fig4_neighbour_check.py
"""
import os

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "results", "source_data")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
SETS = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]


def nn(X):
    D = cdist(X, X)
    np.fill_diagonal(D, np.inf)
    return D.min(1)


def main():
    dist = pd.read_csv(os.path.join(SRC, "stage_centroid_distances_GEO.csv"))
    mem = pd.read_csv(os.path.join(SRC, "ph_subsample_members_GEO.csv"))
    rows = []
    for v in ("adjusted", "unadjusted"):
        for c, st in SETS.items():
            e = pd.read_csv(os.path.join(EMB, f"{c}_{st}_{v}_g2000_pc10.csv"), index_col=0)
            pcs = [x for x in e.columns if x.startswith("PC")]
            for s in [x for x in STAGES if x in set(e.stage)]:
                X = e.loc[e.stage == s, pcs].values
                m = mem[(mem.cohort == c) & (mem.version == v) & (mem["mode"] == "raw") & (mem.stage == s)]
                sub = np.concatenate([nn(e.loc[ids.split(";"), pcs].values) for ids in m.samples])
                d = dist[(dist.cohort == c) & (dist.version == v) & (dist.stage == s)]
                rows.append(dict(version=v, cohort=c, stage=s, n=len(X), m=int(m.m.iloc[0]),
                                 nn_median_all=float(np.median(nn(X))), nn_median_subsamples=float(np.median(sub)),
                                 centre_median=float(d.dist_to_centroid.median()),
                                 centre_median_bias_adj=float(d.dist_to_centroid_bias_adj.median())))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(SRC, "fig4_neighbour_centre_medians.csv"), index=False)
    pd.set_option("display.width", 200)
    print(out.round(2).to_string(index=False))
    # ratios between adjacent stages (later / earlier)
    r = []
    for (v, c), g in out.groupby(["version", "cohort"], sort=False):
        g = g.set_index("stage")
        st = [s for s in STAGES if s in g.index]
        for a, b in zip(st[:-1], st[1:]):
            r.append(dict(version=v, cohort=c, transition=f"{a[0]}→{b[0]}",
                          nn_all=g.loc[b, "nn_median_all"] / g.loc[a, "nn_median_all"],
                          nn_sub=g.loc[b, "nn_median_subsamples"] / g.loc[a, "nn_median_subsamples"],
                          centre=g.loc[b, "centre_median"] / g.loc[a, "centre_median"]))
    print(pd.DataFrame(r).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
