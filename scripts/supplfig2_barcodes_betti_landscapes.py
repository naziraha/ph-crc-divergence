# Figure script for the Nature Computational Science submission (NCS numbering); symbols italic via _italic_symbols.py.
# REVISED COPY (26 Oct 2026): same data and plotting as the original script; only the wording
# inside the figure is cleaned for journal style (no overall title, no internal terms such as
# "step", "addendum" or "frozen", no footnotes that repeat the legend). Outputs go to *_revised folders
# with "_revised" in the file name, so the original figures are not overwritten. No CSV is written.
# Significance asterisks added after p-values (26 Oct 2026): * p < 0.05, ** p < 0.01, *** p < 0.001.
"""
Figure 3 (GSE41258, adjusted data): (a) barcodes of the Fig2_v2 representative subsample (raw-mode subsample 1;
H0 blue, H1 orange); (b) Betti-0 and Betti-1 curves per stage, raw and shape-only, mean and IQR over the 20 subsamples;
(c) H1 persistence landscapes (lambda_1) per stage, raw and shape-only, mean and 2.5-97.5 percentile band over the 20
subsamples, with the addendum 12 H1 landscape Holm p-values. Conventions as the addendum 12 figures (landscape:
gudhi Landscape, 100 points on [0, max H1 death]; Betti grid: 150 points on [0, max death]), except that Betti-0 here
includes the one infinite H0 component (so it falls to 1, not 0).
Reads results/source_data/ph_diagrams_GEO.csv and validation_summary_primary.csv only. No statistical tests.

Usage: python scripts/supplfig2_barcodes_betti_landscapes.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _italic_symbols  # noqa: F401  italic statistical symbols (Zira, 3 Nov 2026)  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from gudhi.representations import Landscape  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "results", "source_data")
OUT = os.path.join(ROOT, "outputs", "supplementary")
NAME = "SupplementaryFig2_barcodes_Betti_landscapes"
COHORT, VERSION, REP, K = "GSE41258", "adjusted", 1, 20
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COL = dict(zip(STAGES, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))   # stage colours (addendum 12)
DIMCOL = {"H0": "#1f77b4", "H1": "#ff7f0e"}                              # dimension colours (Fig2_v2)
MODES = ["raw", "shape-only"]
TR = [("Normal->Adenoma", "N→A"), ("Adenoma->Carcinoma", "A→C"), ("Carcinoma->Metastasis", "C→M")]
RES_L, RES_B = 100, 150

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})



def stars(p):
    """Significance asterisks: * p < 0.05, ** p < 0.01, *** p < 0.001 (added 26 Oct 2026)."""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

def diagrams(d, mode, s, dim):
    x = d[(d["mode"] == mode) & (d.stage == s) & (d.dim == dim)]
    return [x[x.subsample == k][["birth", "death"]].values for k in range(1, K + 1)]


def betti(intervals, grid):
    if len(intervals) == 0:
        return np.zeros_like(grid)
    return ((intervals[:, 0][:, None] <= grid) & (intervals[:, 1][:, None] > grid)).sum(axis=0)


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(os.path.join(SRC, "ph_diagrams_GEO.csv"))
    d = d[(d.cohort == COHORT) & (d.version == VERSION)]
    val = pd.read_csv(os.path.join(SRC, "validation_summary_primary.csv"))
    val = val[val.cohort == COHORT].set_index("transition")

    fig = plt.figure(figsize=(11, 10.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.05, 0.9, 0.9], hspace=0.6, wspace=0.3)

    # ---- row a: barcodes of the representative subsample (raw)
    rep = d[(d["mode"] == "raw") & (d.subsample == REP)]
    xmax = rep.death.max() * 1.08
    for j, s in enumerate(STAGES):
        ax = fig.add_subplot(gs[0, j])
        h0 = rep[(rep.stage == s) & (rep.dim == "H0")].sort_values("death")
        h1 = rep[(rep.stage == s) & (rep.dim == "H1")].sort_values("birth", ascending=False)
        y = 0
        for b, e in h0[["birth", "death"]].values:
            ax.plot([b, e], [y, y], color=DIMCOL["H0"], lw=1.4, solid_capstyle="butt"); y += 1
        ax.plot([0, xmax], [y, y], color=DIMCOL["H0"], lw=1.4, solid_capstyle="butt")       # infinite H0 bar
        ax.plot([xmax], [y], marker=">", color=DIMCOL["H0"], markersize=4, clip_on=False)
        y += 2
        for b, e in h1[["birth", "death"]].values:
            ax.plot([b, e], [y, y], color=DIMCOL["H1"], lw=1.4, solid_capstyle="butt"); y += 1
        ax.set(xlim=(0, xmax), ylim=(-1, y), yticks=[], xlabel="filtration scale (raw)")
        ax.spines["left"].set_visible(False)
        ax.set_title(s, color=COL[s], fontweight="bold")
        ax.text(0.98, 0.30, f"H1: {len(h1)} bars", transform=ax.transAxes, ha="right", fontsize=7, color=DIMCOL["H1"])
        ax.text(0.98, 0.22, f"H0: {len(h0)} + 1 (∞)", transform=ax.transAxes, ha="right", fontsize=7, color=DIMCOL["H0"])
        if j == 0:
            ax.text(-0.18, 1.15, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.text(-0.02, 1.15, "Barcodes, adjusted data, raw distances; representative subsample (n = 31; same samples "
                    "as Fig. 1)", transform=ax.transAxes, fontsize=9)

    # ---- row b: Betti curves (20 subsamples; mean and IQR)
    for j, (dim, mode) in enumerate([("H0", "raw"), ("H0", "shape-only"), ("H1", "raw"), ("H1", "shape-only")]):
        ax = fig.add_subplot(gs[1, j])
        dm = d[(d["mode"] == mode) & (d.dim == dim)]
        grid = np.linspace(0, dm.death.max(), RES_B)
        for s in STAGES:
            B = np.array([betti(x, grid) for x in diagrams(d, mode, s, dim)]) + (1 if dim == "H0" else 0)
            ax.plot(grid, B.mean(0), color=COL[s], lw=1.6, label=s)
            ax.fill_between(grid, np.percentile(B, 25, 0), np.percentile(B, 75, 0), color=COL[s], alpha=0.12)
        ax.set(xlabel=f"filtration scale ({mode})",
               ylabel=f"Betti-{dim[1]} ({'clusters' if dim == 'H0' else 'loops'})", xlim=(0, grid[-1]))
        ax.set_title(f"Betti-{dim[1]}, {mode}")
        if j == 0:
            ax.legend(fontsize=6.5, loc="upper right")
            ax.text(-0.3, 1.15, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.text(-0.1, 1.15, "Betti curves, adjusted data; mean (line) and interquartile range (band) over 20 subsamples "
                    "per stage", transform=ax.transAxes, fontsize=9)

    # ---- row c: H1 landscapes lambda_1 (20 subsamples; mean and 2.5-97.5 percentile band)
    for k, mode in enumerate(MODES):
        ax = fig.add_subplot(gs[2, 2 * k:2 * k + 2])
        dm = d[(d["mode"] == mode) & (d.dim == "H1")]
        rmax = dm.death.max()
        grid = np.linspace(0, rmax, RES_L)
        L = Landscape(num_landscapes=1, resolution=RES_L, sample_range=[0.0, rmax])
        for s in STAGES:
            V = np.array([np.zeros(RES_L) if len(x) == 0 else L.fit_transform([x])[0] for x in diagrams(d, mode, s, "H1")])
            ax.plot(grid, V.mean(0), color=COL[s], lw=1.8, label=s)
            ax.fill_between(grid, np.percentile(V, 2.5, 0), np.percentile(V, 97.5, 0), color=COL[s], alpha=0.10)
        ax.set(xlabel=f"filtration scale ({mode})", ylabel="H1 landscape λ1", xlim=(0, rmax))
        ax.set_title(f"H1 persistence landscape λ1, {mode}")
        col = "H1landscape_raw_p_holm_adjus" if mode == "raw" else "H1landscape_shape_p_holm_adjus"
        ptxt = "   ".join(f"{lab} Holm p = {val.loc[t, col]:.3f}{stars(val.loc[t, col])}" for t, lab in TR)
        ax.text(0.5, -0.3, f"Holm-adjusted p-values (adjusted, {mode}):\n{ptxt}",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5, style="italic")
        if k == 0:
            ax.legend(fontsize=6.5, loc="upper left")
            ax.text(-0.13, 1.12, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

    (lambda *a, **k: None)("GSE41258: barcodes, Betti curves and H1 landscapes by stage (adjusted data; descriptive)",
                 fontsize=11, fontweight="bold", y=0.995)
    (lambda *a, **k: None)(0.5, -0.01, "Curves and bars are descriptive. Bands show variability across subsamples, not confidence intervals. "
             "p-values come from the label-permutation test (one sample per patient per stage). H0 landscapes are not tested (D26).",
             ha="center", fontsize=7.5, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    # numbers for the caption
    for s in STAGES:
        print(s, "rep H0", len(rep[(rep.stage == s) & (rep.dim == "H0")]), "H1", len(rep[(rep.stage == s) & (rep.dim == "H1")]))
    for mode in MODES:
        print(mode, "landscape p:", [val.loc[t, "H1landscape_raw_p_holm_adjus" if mode == "raw" else "H1landscape_shape_p_holm_adjus"] for t, _ in TR])
    print("saved", os.path.join(OUT, NAME + ".png"))


if __name__ == "__main__":
    main()
