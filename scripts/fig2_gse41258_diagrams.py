"""
Figure 2 (GSE41258, adjusted data): H1 persistence diagrams per stage (raw, shape-only) and H0 death values per stage.
Reads only results/source_data/ph_diagrams_GEO.csv and validation_summary_primary.csv (both from scripts/source_data.py).
No statistical tests are run; p-values are copied from the step 5 patient-level permutation test.

Usage: python scripts/fig2_gse41258_diagrams.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "results", "source_data")
OUT = os.path.join(ROOT, "paper", "final", "main_figures")
COHORT, VERSION = "GSE41258", "adjusted"
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COL = dict(zip(STAGES, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))   # as addendum 12 figures
MODES = [("raw", "raw distances"), ("shape-only", "shape-only (distances / median pairwise distance)")]
TR = [("Normal->Adenoma", "N→A"), ("Adenoma->Carcinoma", "A→C"), ("Carcinoma->Metastasis", "C→M")]

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(os.path.join(SRC, "ph_diagrams_GEO.csv"))
    d = d[(d.cohort == COHORT) & (d.version == VERSION)]
    val = pd.read_csv(os.path.join(SRC, "validation_summary_primary.csv"))
    val = val[val.cohort == COHORT].set_index("transition")

    fig = plt.figure(figsize=(11, 9.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 0.9], hspace=0.62, wspace=0.28)
    for row, (mode, mlab) in enumerate(MODES):
        h1 = d[(d["mode"] == mode) & (d.dim == "H1")]
        hi = h1.death.max() * 1.05
        col_p = "raw_p_holm_adjus" if mode == "raw" else "shape_p_holm_adjus"
        for j, s in enumerate(STAGES):
            ax = fig.add_subplot(gs[row, j])
            x = h1[h1.stage == s]
            ax.scatter(x.birth, x.death, s=7, alpha=0.45, color=COL[s], edgecolor="none")
            ax.plot([0, hi], [0, hi], color="0.5", lw=0.8)
            ax.set(xlim=(0, hi), ylim=(0, hi), xlabel="birth", ylabel="death" if j == 0 else "")
            ax.set_aspect("equal")
            ax.set_title(f"{s}  (n = {len(x)} loops)", color=COL[s], fontweight="bold")
            if j == 0:
                ax.text(-0.42, 1.18, "ab"[row], transform=ax.transAxes, fontsize=13, fontweight="bold")
                ax.text(-0.22, 1.18, f"H1 persistence diagrams, adjusted data, {mlab}; 20 subsamples pooled per stage",
                        transform=ax.transAxes, fontsize=9)
        ptxt = "   ".join(f"{lab} Holm p = {val.loc[t, col_p]:.3f}" for t, lab in TR)
        pos = fig.add_subplot(gs[row, 0]).get_position()
        fig.axes[-1].remove()
        fig.text(0.5, pos.y0 - 0.058, f"Step 5 patient-level permutation test (H0 + H1, adjusted, {mode}):   {ptxt}",
                 ha="center", fontsize=8, style="italic")

    for k, (mode, mlab) in enumerate(MODES):
        ax = fig.add_subplot(gs[2, 2 * k:2 * k + 2])
        h0 = d[(d["mode"] == mode) & (d.dim == "H0")]
        data = [h0[h0.stage == s].death.values for s in STAGES]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        flierprops=dict(marker="o", markersize=2, alpha=0.4, markeredgecolor="none"),
                        medianprops=dict(color="black", lw=1))
        for patch, fl, s in zip(bp["boxes"], bp["fliers"], STAGES):
            patch.set_facecolor(COL[s]); patch.set_alpha(0.75); fl.set_markerfacecolor(COL[s])
        ax.set_xticks(range(1, 5), [f"{s}\n(n = {len(v)})" for s, v in zip(STAGES, data)])
        ax.set(ylabel="H0 death (merge scale)", title=f"H0 death values, adjusted data, {mode}")
        if k == 0:
            ax.text(-0.13, 1.1, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

    fig.suptitle("GSE41258: persistence diagrams by stage (adjusted data; descriptive)", fontsize=11, fontweight="bold", y=0.995)
    fig.text(0.5, 0.005, "Points and boxes pool all features from 20 subsamples per stage (m = 31) and are descriptive. "
             "p-values come from the step 5 patient-level permutation test, not from these points.",
             ha="center", fontsize=7.5, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"Fig2_GSE41258_diagrams.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", os.path.join(OUT, "Fig2_GSE41258_diagrams.png"))


if __name__ == "__main__":
    main()
