# Figure script for the Nature Computational Science submission (NCS numbering); symbols italic via _italic_symbols.py.
# NCS VERSION (2 Nov 2026), after Prof Peng's comments: copy of fig4_heterogeneity_revised.py for
# Fig. 2 (adjusted) and Supplementary Fig. 3 (unadjusted) of the Nature Computational Science submission.
# Changes: larger fonts and points; Holm-adjusted PERMDISP p-values drawn as brackets between adjacent
# stages in panel a (same stored values as panel b; no new test). Original scripts untouched.
# REVISED COPY (26 Oct 2026): same data and plotting as the original script; only the wording
# inside the figure is cleaned for journal style (no overall title, no internal terms such as
# "step", "addendum" or "frozen", no footnotes that repeat the legend). Outputs go to *_revised folders
# with "_revised" in the file name, so the original figures are not overwritten. No CSV is written.
# Significance asterisks added after p-values (26 Oct 2026): * p < 0.05, ** p < 0.01, *** p < 0.001.
"""
Figure 4 (heterogeneity, GEO cohorts): (a) each sample's bias-adjusted distance to its stage centroid
(sqrt(n/(n-1)) * ||x - c||, as addendum 9) per stage and cohort; (b) PERMDISP dispersion ratio D per cohort x transition
with Holm p (step 13); main (adjusted) figure only: (c) TCGA-COAD PERMDISP D for adjacent stages I-IV with Holm p
(addendum 10, bias-adjusted distances; the values of Supplementary Table S_TCGA_stages part B). Reads
results/source_data/stage_centroid_distances_GEO.csv, results/step13/benchmark_results.csv and
results/addendum10/spread_H0H1/permdisp_by_stage.csv only. No statistical tests are run.

Usage: python scripts/fig2_dispersion.py [adjusted|unadjusted]
  adjusted   -> outputs/main_figures/Fig4_heterogeneity.{png,pdf}            (main)
  unadjusted -> outputs/supplementary/FigS_heterogeneity_unadjusted.{png,pdf} (candidate supplementary)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _italic_symbols  # noqa: F401  italic statistical symbols, Nature style (Zira, 3 Nov 2026)  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COL = dict(zip(STAGES, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))   # Fig 2/3 stage colours
TRL = {"Normal->Adenoma": "N→A", "Adenoma->Carcinoma": "A→C", "Carcinoma->Metastasis": "C→M"}
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False})



def stars(p):
    """Significance asterisks: * p < 0.05, ** p < 0.01, *** p < 0.001 (added 26 Oct 2026)."""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

def main(version):
    if version == "adjusted":
        out, name = os.path.join(ROOT, "outputs", "main_figures"), "Fig2_dispersion"
    else:
        out, name = os.path.join(ROOT, "outputs", "supplementary"), "SupplementaryFig3_dispersion_unadjusted"
    os.makedirs(out, exist_ok=True)
    dist = pd.read_csv(os.path.join(ROOT, "results", "source_data", "stage_centroid_distances_GEO.csv"))
    dist = dist[dist.version == version]
    bench = pd.read_csv(os.path.join(ROOT, "results", "step13", "benchmark_results.csv"))
    bench = bench[bench.version == version]
    rng = np.random.default_rng(42)                                       # jitter only
    tcga = version == "adjusted"                                          # panel c in the main figure only
    if tcga:
        pd_t = pd.read_csv(os.path.join(ROOT, "results", "addendum10", "spread_H0H1", "permdisp_by_stage.csv"))
        pd_t = pd_t[pd_t.version == "adjusted"].set_index("transition")
        tr_t = ["Stage I->Stage II", "Stage II->Stage III", "Stage III->Stage IV"]
        # check: same values as Supplementary Table S_TCGA_stages part B (adjusted)
        tab = pd.read_csv(os.path.join(ROOT, "outputs", "supplementary", "TableS_TCGA_stages.csv"))
        tab = tab[tab.section == "B. PERMDISP"].set_index("transition")
        for t in tr_t:
            lab = t.replace("Stage ", "").replace("->", " → ")
            assert round(float(pd_t.loc[t, "D"]), 3) == float(tab.loc[lab, "D, adjusted"]), t
            assert float(pd_t.loc[t, "p_holm"]) == float(tab.loc[lab, "Holm p, adjusted"]), t

    fig = plt.figure(figsize=(12, 10.6 if tcga else 8.6))
    gs = fig.add_gridspec(3 if tcga else 2, 4, height_ratios=[1, 0.95, 0.42] if tcga else [1, 0.95],
                          hspace=0.62 if tcga else 0.55, wspace=0.32)
    for j, c in enumerate(GEO):
        ax = fig.add_subplot(gs[0, j])
        d = dist[dist.cohort == c]
        st = [s for s in STAGES if s in set(d.stage)]
        data = [d[d.stage == s].dist_to_centroid_bias_adj.values for s in st]
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False, medianprops=dict(color="black", lw=1))
        for k, (patch, s) in enumerate(zip(bp["boxes"], st)):
            patch.set_facecolor(COL[s]); patch.set_alpha(0.35); patch.set_edgecolor(COL[s])
            ax.scatter(k + 1 + rng.uniform(-0.18, 0.18, len(data[k])), data[k], s=11, color=COL[s], alpha=0.75,
                       edgecolor="none", zorder=3)
        ax.set_xticks(range(1, len(st) + 1), [f"{s} ({len(v)})" for s, v in zip(st, data)], fontsize=9.5, rotation=30, ha="right")
        ymax = max(np.concatenate(data))
        # brackets with Holm-adjusted PERMDISP p between adjacent stages (stored values, same as panel b)
        level = 0
        for k in range(len(st) - 1):
            t = f"{st[k]}->{st[k + 1]}"
            r = bench[(bench.cohort == c) & (bench.transition == t)]
            if not len(r):
                continue
            p = float(r.permdisp_p_holm.iloc[0])
            yb = ymax * (1.06 + 0.14 * (level % 2))
            ax.plot([k + 1, k + 1, k + 2, k + 2], [yb - ymax * 0.02, yb, yb, yb - ymax * 0.02], color="black", lw=0.9)
            ax.text(k + 1.5, yb + ymax * 0.005, stars(p) or "ns", ha="center", va="bottom",
                    fontsize=13 if p < 0.05 else 9.5, fontweight="bold" if p < 0.05 else "normal")  # stars only (Zira, 2 Nov)
            level += 1
        ax.set_ylim(0, ymax * 1.38)
        ax.set_title(c, fontweight="bold", pad=10)
        if j == 0:
            ax.set_ylabel("distance to stage centroid\n(10 PCs, bias-adjusted)")
            ax.text(-0.42, 1.2, "a", transform=ax.transAxes, fontsize=15, fontweight="bold")
            ax.text(-0.2, 1.2, f"Between-sample spread within each stage, {version} data", transform=ax.transAxes, fontsize=11)

    ax = fig.add_subplot(gs[1, 1:4])
    rows = [(c, t) for c in GEO for t in TRL if len(bench[(bench.cohort == c) & (bench.transition == t)])]
    y = np.arange(len(rows))[::-1]
    lo = min(bench.permdisp_effect.min(), 0.9)
    hi = bench.permdisp_effect.max()
    for yi, (c, t) in zip(y, rows):
        r = bench[(bench.cohort == c) & (bench.transition == t)].iloc[0]
        a, b = t.split("->")
        sig = r.permdisp_p_holm < 0.05
        ax.scatter(r.permdisp_effect, yi, s=70, color=COL[b], edgecolor="black" if sig else COL[b], linewidth=0.8, zorder=3)
        ax.text(hi + 0.06, yi, f"Holm p = {r.permdisp_p_holm:.3f}{stars(r.permdisp_p_holm)}", va="center", fontsize=9.5,
                fontweight="bold" if sig else "normal")
    ax.axvline(1, color="0.4", lw=0.8, ls="--")
    ax.set_yticks(y, [f"{c}  {TRL[t]}" for c, t in rows], fontsize=9.5)
    for k in range(1, len(GEO)):                                           # separators between cohorts
        n_before = sum(1 for c, _ in rows if GEO.index(c) < k)
        ax.axhline(y[n_before - 1] - 0.5, color="0.85", lw=0.6)
    ax.set_xlim(lo - 0.05, hi + 0.38)
    ax.set_xlabel("PERMDISP dispersion ratio D (later stage / earlier stage)")
    ax.set_title(f"Change in spread between adjacent stages, {version} data")
    ax.text(-0.42, 1.08, "b", transform=ax.transAxes, fontsize=15, fontweight="bold")
    ax.text(0.0, -0.26, "D > 1: later stage more spread out; dot colour = later stage; black outline = Holm p < 0.05.",
            transform=ax.transAxes, fontsize=9, color="0.3")
    xl = ax.get_xlim()

    if tcga:                                                               # panel c: TCGA-COAD stages
        ax = fig.add_subplot(gs[2, 1:4])
        y = np.arange(len(tr_t))[::-1]
        for yi, t in zip(y, tr_t):
            D, p = float(pd_t.loc[t, "D"]), float(pd_t.loc[t, "p_holm"])
            sig = p < 0.05
            ax.scatter(D, yi, s=70, color=COL["Carcinoma"], edgecolor="black" if sig else COL["Carcinoma"], linewidth=0.8, zorder=3)
            ax.text(hi + 0.06, yi, f"Holm p = {p:.3f}{stars(p)}", va="center", fontsize=9.5, fontweight="bold" if sig else "normal")
        ax.axvline(1, color="0.4", lw=0.8, ls="--")
        ax.set_yticks(y, [f"TCGA-COAD  {t.replace('Stage ', '').replace('->', '→')}" for t in tr_t], fontsize=9.5)
        ax.set_ylim(-0.6, len(tr_t) - 0.4)
        ax.set_xlim(*xl)                                                   # same x-axis as panel b
        ax.set_xlabel("PERMDISP dispersion ratio D (later stage / earlier stage)")
        ax.set_title("TCGA-COAD: change in spread between adjacent stages, adjusted data")
        ax.text(-0.42, 1.22, "c", transform=ax.transAxes, fontsize=15, fontweight="bold")

    title = ("Between-sample heterogeneity across stages in the GEO cohorts and TCGA-COAD" if tcga else
             "Between-sample heterogeneity across stages in the GEO cohorts")
    (lambda *a, **k: None)(f"{title} ({version} data)", fontsize=11, fontweight="bold", y=0.995)
    foot = ("Panel a is descriptive (points = samples). p-values in b come from the step 13 PERMDISP label-permutation "
            "test (999 permutations; Holm across transitions within cohort), not from the points in a.")
    if tcga:
        foot += ("\nPanel c: addendum 10 PERMDISP on bias-adjusted distances (999 permutations; Holm across the three "
                 "transitions); unadjusted values in Supplementary Table S_TCGA_stages.")
    (lambda *a, **k: None)(0.5, 0.005, foot, ha="center", fontsize=9.5, color="0.3", **({"linespacing": 1.4} if tcga else {}))
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
