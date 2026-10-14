"""
Figure 4 (heterogeneity, GEO cohorts): (a) each sample's bias-adjusted distance to its stage centroid
(sqrt(n/(n-1)) * ||x - c||, as addendum 9) per stage and cohort; (b) PERMDISP dispersion ratio D per cohort x transition
with Holm p (step 13). Reads results/source_data/stage_centroid_distances_GEO.csv and results/step13/benchmark_results.csv
only. No statistical tests are run.

Usage: python scripts/fig4_heterogeneity.py [adjusted|unadjusted]
  adjusted   -> outputs/main_figures/Fig4_heterogeneity.{png,pdf}            (main)
  unadjusted -> outputs/supplementary/FigS_heterogeneity_unadjusted.{png,pdf} (candidate supplementary)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COL = dict(zip(STAGES, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))   # Fig 2/3 stage colours
TRL = {"Normal->Adenoma": "N→A", "Adenoma->Carcinoma": "A→C", "Carcinoma->Metastasis": "C→M"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def main(version):
    if version == "adjusted":
        out, name = os.path.join(ROOT, "outputs", "main_figures"), "Fig4_heterogeneity"
    else:
        out, name = os.path.join(ROOT, "outputs", "supplementary"), "FigS_heterogeneity_unadjusted"
    os.makedirs(out, exist_ok=True)
    dist = pd.read_csv(os.path.join(ROOT, "results", "source_data", "stage_centroid_distances_GEO.csv"))
    dist = dist[dist.version == version]
    bench = pd.read_csv(os.path.join(ROOT, "results", "step13", "benchmark_results.csv"))
    bench = bench[bench.version == version]
    rng = np.random.default_rng(42)                                       # jitter only

    fig = plt.figure(figsize=(11, 7.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.95], hspace=0.55, wspace=0.32)
    for j, c in enumerate(GEO):
        ax = fig.add_subplot(gs[0, j])
        d = dist[dist.cohort == c]
        st = [s for s in STAGES if s in set(d.stage)]
        data = [d[d.stage == s].dist_to_centroid_bias_adj.values for s in st]
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False, medianprops=dict(color="black", lw=1))
        for k, (patch, s) in enumerate(zip(bp["boxes"], st)):
            patch.set_facecolor(COL[s]); patch.set_alpha(0.35); patch.set_edgecolor(COL[s])
            ax.scatter(k + 1 + rng.uniform(-0.18, 0.18, len(data[k])), data[k], s=6, color=COL[s], alpha=0.75,
                       edgecolor="none", zorder=3)
        ax.set_xticks(range(1, len(st) + 1), [f"{s} ({len(v)})" for s, v in zip(st, data)], fontsize=7, rotation=30, ha="right")
        ax.set_ylim(0, max(np.concatenate(data)) * 1.06)
        ax.set_title(c, fontweight="bold")
        if j == 0:
            ax.set_ylabel("distance to stage centroid\n(10 PCs, bias-adjusted)")
            ax.text(-0.42, 1.12, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.text(-0.2, 1.12, f"Between-sample spread within each stage, {version} data", transform=ax.transAxes, fontsize=9)

    ax = fig.add_subplot(gs[1, 1:3])
    rows = [(c, t) for c in GEO for t in TRL if len(bench[(bench.cohort == c) & (bench.transition == t)])]
    y = np.arange(len(rows))[::-1]
    lo = min(bench.permdisp_effect.min(), 0.9)
    hi = bench.permdisp_effect.max()
    for yi, (c, t) in zip(y, rows):
        r = bench[(bench.cohort == c) & (bench.transition == t)].iloc[0]
        a, b = t.split("->")
        sig = r.permdisp_p_holm < 0.05
        ax.scatter(r.permdisp_effect, yi, s=46, color=COL[b], edgecolor="black" if sig else COL[b], linewidth=0.8, zorder=3)
        ax.text(hi + 0.06, yi, f"Holm p = {r.permdisp_p_holm:.3f}", va="center", fontsize=7.5,
                fontweight="bold" if sig else "normal")
    ax.axvline(1, color="0.4", lw=0.8, ls="--")
    ax.set_yticks(y, [f"{c}  {TRL[t]}" for c, t in rows], fontsize=7.5)
    for k in range(1, len(GEO)):                                           # separators between cohorts
        n_before = sum(1 for c, _ in rows if GEO.index(c) < k)
        ax.axhline(y[n_before - 1] - 0.5, color="0.85", lw=0.6)
    ax.set_xlim(lo - 0.05, hi + 0.38)
    ax.set_xlabel("PERMDISP dispersion ratio D (later stage / earlier stage)")
    ax.set_title(f"Change in spread between adjacent stages, {version} data (step 13)")
    ax.text(-0.36, 1.08, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.text(0.0, -0.26, "D > 1: later stage more spread out; dot colour = later stage; black outline = Holm p < 0.05.",
            transform=ax.transAxes, fontsize=7, color="0.3")

    fig.suptitle(f"Between-sample heterogeneity across stages in the GEO cohorts ({version} data)", fontsize=11, fontweight="bold", y=0.995)
    fig.text(0.5, 0.005, "Panel a is descriptive (points = samples). p-values in b come from the step 13 PERMDISP label-permutation "
             "test (999 permutations; Holm across transitions within cohort), not from the points in a.",
             ha="center", fontsize=7.5, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    # check: D from the plotted (plain) distances reproduces step 13
    for c, t in rows:
        a, b = t.split("->")
        d = dist[dist.cohort == c]
        D = d[d.stage == b].dist_to_centroid.mean() / d[d.stage == a].dist_to_centroid.mean()
        Dref = bench[(bench.cohort == c) & (bench.transition == t)].permdisp_effect.iloc[0]
        assert abs(D - Dref) < 1e-9, (c, t, D, Dref)
    print("check: D from source-data distances = step 13 PERMDISP D (all rows)")
    print("saved", os.path.join(out, name + ".png"))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "adjusted")
