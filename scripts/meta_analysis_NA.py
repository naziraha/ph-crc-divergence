"""
Secondary meta-analysis (addendum 3, commit d38db92): DerSimonian-Laird random-effects
meta-analysis of the shape-only Normal -> Adenoma permutation z-scores across GSE41258,
GSE20916 (bulk) and GSE37364, separately for unadjusted, adjusted and stage-aware data.
Effect = z, within-cohort variance = 1. statsmodels combine_effects(method_re="dl").
Written after GSE41258 and GSE20916 were analysed (stated in addendum 3).
"""
import os
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.meta_analysis import combine_effects

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "meta_analysis")
os.makedirs(OUT, exist_ok=True)
TR = "Normal->Adenoma"


def z_of(path, version, analysis="primary"):
    d = pd.read_csv(path)
    d = d[(d.analysis == analysis) & (d.transition == TR) & (d["mode"] == "shape-only") & (d.version == version)]
    assert len(d) == 1, (path, version)
    return float(d.z.iloc[0]), float(d.p_holm.iloc[0])


def main():
    s5 = os.path.join(ROOT, "results", "step05")
    d3 = pd.read_csv(os.path.join(ROOT, "results", "exploratory_D3", "permutation_results_D3.csv"))
    d3 = d3[(d3.transition == TR) & (d3["mode"] == "shape-only")]
    rows = []
    for version in ("unadjusted", "adjusted", "stageaware"):
        if version == "stageaware":
            z41, p41 = float(d3.z.iloc[0]), float(d3.p_holm.iloc[0])
        else:
            z41, p41 = z_of(os.path.join(s5, "permutation_results_all.csv"), version)
        z16, p16 = z_of(os.path.join(s5, "GSE20916", "permutation_results_all.csv"), version)
        z37, p37 = z_of(os.path.join(s5, "GSE37364", "permutation_results_all.csv"), version)
        for c, z, p in (("GSE41258", z41, p41), ("GSE20916", z16, p16), ("GSE37364", z37, p37)):
            rows.append(dict(version=version, cohort=c, z=z, p_holm=p))
    inp = pd.DataFrame(rows)
    inp.to_csv(os.path.join(OUT, "meta_inputs.csv"), index=False)
    res_rows, L = [], ["SECONDARY META-ANALYSIS (addendum 3, commit d38db92)", "=" * 60,
                       "Shape-only Normal -> Adenoma permutation z-scores; DerSimonian-Laird random effects;",
                       "within-cohort variance 1. Written after GSE41258 and GSE20916 had been analysed:",
                       "an informed summary, not an independent confirmation. Does not define a tipping point.",
                       "GSE41258 stage-aware input = exploratory analysis D3.",
                       "statsmodels " + metadata.version("statsmodels"), ""]
    for version in ("unadjusted", "adjusted", "stageaware"):
        d = inp[inp.version == version]
        r = combine_effects(d.z.values, np.ones(len(d)), method_re="dl", row_names=list(d.cohort))
        sf = r.summary_frame()
        re = sf.loc["random effect"]
        p = 2 * norm.sf(abs(re["eff"] / re["sd_eff"]))
        rec = dict(version=version, pooled_z=re["eff"], se=re["sd_eff"], ci_low=re["ci_low"], ci_upp=re["ci_upp"],
                   p_two_sided=p, tau2=r.tau2, Q=r.q, Q_p=float(r.test_homogeneity().pvalue), I2=r.i2)
        res_rows.append(rec)
        L += [f"--- {version} ---"]
        for _, x in d.iterrows():
            L.append(f"  {x.cohort:<9} z = {x.z:5.2f}  (per-cohort Holm p {x.p_holm:.4f})  "
                     f"RE weight {sf.loc[x.cohort, 'w_re'] * 100:4.1f}%")
        L.append(f"  pooled z (RE) = {rec['pooled_z']:.2f} [95% CI {rec['ci_low']:.2f}, {rec['ci_upp']:.2f}], "
                 f"p = {rec['p_two_sided']:.2g}; tau2 = {rec['tau2']:.2f}; Q = {rec['Q']:.2f} "
                 f"(p {rec['Q_p']:.3f}); I2 = {rec['I2'] * 100:.0f}%")
        L.append("")
    res = pd.DataFrame(res_rows)
    res.to_csv(os.path.join(OUT, "meta_results.csv"), index=False)
    L += ["Caution: 3 cohorts only; tau2 and I2 are imprecise.",
          "Outputs: results/meta_analysis/meta_inputs.csv, meta_results.csv, forest_NA_shape.png"]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    for ax, version in zip(axes, ("unadjusted", "adjusted", "stageaware")):
        d = inp[inp.version == version].reset_index(drop=True)
        rr = res[res.version == version].iloc[0]
        ys = np.arange(len(d))[::-1] + 1
        ax.errorbar(d.z, ys, xerr=1.96, fmt="s", color="k", capsize=3)
        ax.errorbar([rr.pooled_z], [0], xerr=[[rr.pooled_z - rr.ci_low], [rr.ci_upp - rr.pooled_z]], fmt="D",
                    color="#eb6834", capsize=4, ms=8)
        ax.axvline(0, color="0.5", lw=1); ax.axvline(1.96, color="0.7", ls=":", lw=1)
        ax.set_yticks(list(ys) + [0]); ax.set_yticklabels(list(d.cohort) + ["pooled (RE)"])
        ax.set_title(f"{version}: pooled z {rr.pooled_z:.2f} [{rr.ci_low:.2f}, {rr.ci_upp:.2f}]", fontsize=9)
        ax.set_xlabel("z, shape-only Normal → Adenoma")
    fig.suptitle("Secondary meta-analysis (DerSimonian–Laird); informed summary, not independent confirmation", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "forest_NA_shape.png"), dpi=150)
    print("\n".join(L))


if __name__ == "__main__":
    main()
