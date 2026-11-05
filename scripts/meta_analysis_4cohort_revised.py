"""
Secondary 4-cohort meta-analysis (addendum 6 section 8, commit 93613bd; extends addendum 3).

DerSimonian-Laird random effects (statsmodels combine_effects, method_re = "dl") of the shape-only
permutation z-scores, within-cohort variance 1, for Normal -> Adenoma and Adenoma -> Carcinoma,
separately for unadjusted, adjusted (stage-blind) and stage-aware data. Cohorts: GSE41258 (main;
stage-aware = exploratory D3), GSE20916 (bulk), GSE37364 (bulk), GSE117606 (main; stage-aware from
the secondary run). INFORMED SUMMARY: the three earlier cohorts' results were known; it defines no
tipping point. The 3-cohort analysis (results/meta_analysis/) is left unchanged.

REVISED 23 Oct 2026 (deviation D32): statsmodels 0.15.0 combine_effects(method_re="dl") does not truncate a
negative DerSimonian-Laird tau2 at zero before computing the random-effects SE, CI and p. This version computes
the DL model explicitly with tau2 = max(0, (Q - (k - 1)) / C), the textbook estimator. It writes to
results/meta_analysis_4cohort_revised/ so the original outputs stay unchanged. Only rows with Q < k - 1 change.
"""
import os
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import chi2


def dersimonian_laird(y, v):
    """Textbook DerSimonian-Laird random-effects meta-analysis with tau2 truncated at zero."""
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    w = 1 / v
    mu_fe = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - mu_fe) ** 2))
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C)                      # truncation at zero
    ws = 1 / (v + tau2)
    mu = np.sum(ws * y) / np.sum(ws)
    se = np.sqrt(1 / np.sum(ws))
    return dict(eff=mu, sd_eff=se, ci_low=mu - norm.ppf(0.975) * se, ci_upp=mu + norm.ppf(0.975) * se,
                tau2=tau2, Q=Q, Q_p=float(chi2.sf(Q, k - 1)), I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "meta_analysis_4cohort_revised")
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
              "DerSimonian-Laird computed explicitly, tau2 truncated at zero (D32); numpy " + metadata.version("numpy"), ""]
for tr in TRS:
    for v in VERS:
        d = inp[(inp.transition == tr) & (inp.version == v)]
        re = dersimonian_laird(d.z.values, np.ones(len(d)))
        p = 2 * norm.sf(abs(re["eff"] / re["sd_eff"]))
        res.append(dict(transition=tr, version=v, pooled_z=re["eff"], se=re["sd_eff"], ci_low=re["ci_low"],
                        ci_upp=re["ci_upp"], p_two_sided=p, tau2=re["tau2"], I2=re["I2"], Q=re["Q"], Q_p=re["Q_p"]))
        L.append(f"--- {tr}, {v} ---")
        L += [f"  {x.cohort:<10} z = {x.z:6.2f}  (Holm p {x.p_holm:.4f})" for _, x in d.iterrows()]
        L.append(f"  pooled z = {re['eff']:.2f} [95% CI {re['ci_low']:.2f}, {re['ci_upp']:.2f}], p = {p:.2g}; "
                 f"tau2 = {re['tau2']:.2f}; I2 = {re['I2'] * 100:.0f}% (Q p {re['Q_p']:.3f})")
        L.append("")
res = pd.DataFrame(res)
res.to_csv(os.path.join(OUT, "meta_results.csv"), index=False)
L += ["tau2 truncated at zero before pooling (textbook DL). Caution: 4 cohorts; heterogeneity estimates are imprecise.",
      "Outputs: results/meta_analysis_4cohort_revised/meta_inputs.csv, meta_results.csv, forest_4cohort.png"]
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
fig.suptitle("Secondary 4-cohort meta-analysis (DerSimonian–Laird, tau2 truncated at 0); informed summary", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "forest_4cohort.png"), dpi=150)
print("\n".join(L))
