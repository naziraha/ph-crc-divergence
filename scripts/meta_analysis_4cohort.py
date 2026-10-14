"""
Secondary 4-cohort meta-analysis (addendum 6 section 8, commit 93613bd; extends addendum 3).

DerSimonian-Laird random effects (statsmodels combine_effects, method_re = "dl") of the shape-only
permutation z-scores, within-cohort variance 1, for Normal -> Adenoma and Adenoma -> Carcinoma,
separately for unadjusted, adjusted (stage-blind) and stage-aware data. Cohorts: GSE41258 (main;
stage-aware = exploratory D3), GSE20916 (bulk), GSE37364 (bulk), GSE117606 (main; stage-aware from
the secondary run). INFORMED SUMMARY: the three earlier cohorts' results were known; it defines no
tipping point. The 3-cohort analysis (results/meta_analysis/) is left unchanged.
"""
import os
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.meta_analysis import combine_effects

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "meta_analysis_4cohort")
os.makedirs(OUT, exist_ok=True)
S5 = os.path.join(ROOT, "results", "step05")
TRS = ("Normal->Adenoma", "Adenoma->Carcinoma")
VERS = ("unadjusted", "adjusted", "stageaware")


def get(df, tr, version, analysis="primary"):
    q = df[(df.analysis == analysis) & (df.transition == tr) & (df["mode"] == "shape-only") & (df.version == version)]
    assert len(q) == 1, (analysis, tr, version, len(q))
    return float(q.z.iloc[0]), float(q.p_holm.iloc[0])


g41 = pd.read_csv(os.path.join(S5, "permutation_results_all.csv"))
d3 = pd.read_csv(os.path.join(ROOT, "results", "exploratory_D3", "permutation_results_D3.csv")).assign(
    analysis="primary", version="stageaware")
g16 = pd.read_csv(os.path.join(S5, "GSE20916", "permutation_results_all.csv"))
g37 = pd.read_csv(os.path.join(S5, "GSE37364", "permutation_results_all.csv"))
g117 = pd.read_csv(os.path.join(S5, "GSE117606", "secondary_results_all.csv"))
g117 = g117.assign(analysis=g117.analysis.replace({"PRIMARY": "primary"}))
g117.loc[g117.analysis == "stageaware", "analysis"] = "primary"

rows = []
for tr in TRS:
    for v in VERS:
        for c, df in (("GSE41258", d3 if v == "stageaware" else g41), ("GSE20916", g16), ("GSE37364", g37),
                      ("GSE117606", g117)):
            z, p = get(df, tr, v)
            rows.append(dict(transition=tr, version=v, cohort=c, z=z, p_holm=p))
inp = pd.DataFrame(rows)
inp.to_csv(os.path.join(OUT, "meta_inputs.csv"), index=False)

res, L = [], ["SECONDARY 4-COHORT META-ANALYSIS (addendum 6 section 8; informed summary)", "=" * 66,
              "Shape-only permutation z; DerSimonian-Laird random effects; within-cohort variance 1.",
              "GSE41258 stage-aware = exploratory D3. Three of the four cohorts were analysed before this was",
              "specified, so it is an informed summary, not independent confirmation; it defines no tipping point.",
              "statsmodels " + metadata.version("statsmodels"), ""]
for tr in TRS:
    for v in VERS:
        d = inp[(inp.transition == tr) & (inp.version == v)]
        r = combine_effects(d.z.values, np.ones(len(d)), method_re="dl", row_names=list(d.cohort))
        re = r.summary_frame().loc["random effect"]
        p = 2 * norm.sf(abs(re["eff"] / re["sd_eff"]))
        res.append(dict(transition=tr, version=v, pooled_z=re["eff"], ci_low=re["ci_low"], ci_upp=re["ci_upp"],
                        p_two_sided=p, tau2=max(r.tau2, 0.0), I2=max(r.i2, 0.0),
                        Q=r.q, Q_p=float(r.test_homogeneity().pvalue)))
        L.append(f"--- {tr}, {v} ---")
        L += [f"  {x.cohort:<10} z = {x.z:6.2f}  (Holm p {x.p_holm:.4f})" for _, x in d.iterrows()]
        L.append(f"  pooled z = {re['eff']:.2f} [95% CI {re['ci_low']:.2f}, {re['ci_upp']:.2f}], p = {p:.2g}; "
                 f"tau2 = {max(r.tau2, 0):.2f}; I2 = {max(r.i2, 0) * 100:.0f}% (Q p {float(r.test_homogeneity().pvalue):.3f})")
        L.append("")
res = pd.DataFrame(res)
res.to_csv(os.path.join(OUT, "meta_results.csv"), index=False)
L += ["Negative tau2 / I2 shown as 0. Caution: 4 cohorts; heterogeneity estimates are imprecise.",
      "Outputs: results/meta_analysis_4cohort/meta_inputs.csv, meta_results.csv, forest_4cohort.png"]
open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
fig, axes = plt.subplots(2, 3, figsize=(14, 6.5), sharey=True)
for i, tr in enumerate(TRS):
    for j, v in enumerate(VERS):
        ax = axes[i, j]
        d = inp[(inp.transition == tr) & (inp.version == v)].reset_index(drop=True)
        rr = res[(res.transition == tr) & (res.version == v)].iloc[0]
        ys = np.arange(len(d))[::-1] + 1
        ax.errorbar(d.z, ys, xerr=1.96, fmt="s", color="k", capsize=3)
        ax.errorbar([rr.pooled_z], [0], xerr=[[rr.pooled_z - rr.ci_low], [rr.ci_upp - rr.pooled_z]], fmt="D",
                    color="#eb6834", capsize=4, ms=8)
        ax.axvline(0, color="0.5", lw=1); ax.axvline(1.96, color="0.7", ls=":", lw=1)
        ax.set_yticks(list(ys) + [0]); ax.set_yticklabels(list(d.cohort) + ["pooled (RE)"])
        ax.set_title(f"{tr.replace('->', ' → ')}, {v}: {rr.pooled_z:.2f} [{rr.ci_low:.2f}, {rr.ci_upp:.2f}]", fontsize=9)
        ax.set_xlabel("shape-only z")
fig.suptitle("Secondary 4-cohort meta-analysis (DerSimonian–Laird); informed summary", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "forest_4cohort.png"), dpi=150)
print("\n".join(L))
