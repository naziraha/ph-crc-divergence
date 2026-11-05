# Figure script for the Nature Computational Science submission (NCS numbering); symbols italic via _italic_symbols.py.
# REVISED COPY (26 Oct 2026): same data and plotting as the original script; only the wording
# inside the figure is cleaned for journal style (no overall title, no internal terms such as
# "step", "addendum" or "frozen", no footnotes that repeat the legend). Outputs go to *_revised folders
# with "_revised" in the file name, so the original figures are not overwritten. No CSV is written.
"""
Supplementary Figure S_CMS_PC1_PC2 (formerly Figure 5 panel a): GEO carcinomas on PC1 vs PC2 of the frozen adjusted
embedding, coloured by primary CMS call (addendum 9); hollow = not confidently classified. Also writes the descriptive
share of carcinoma variation that lines up with CMS (R2 = 1 - within-subtype / total sum of squares) in PC1-2 and in
all 10 PCs (results/source_data/figS_cms_r2_pcs.csv). No statistical tests.

Usage: python scripts/supplfig4_cms_pc12.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _italic_symbols  # noqa: F401  italic statistical symbols (Zira, 3 Nov 2026)  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A9 = os.path.join(ROOT, "results", "addendum9")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "outputs", "supplementary")
NAME = "SupplementaryFig4_CMS_PC1_PC2"
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
SETS = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
CMS = ["CMS1", "CMS2", "CMS3", "CMS4"]
CCOL = {"CMS1": "#e79e00", "CMS2": "#0071b1", "CMS3": "#ca78a6", "CMS4": "#009f74"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def r2(P, g):
    tot = ((P - P.mean(0)) ** 2).sum()
    w = sum(((P[g == k] - P[g == k].mean(0)) ** 2).sum() for k in np.unique(g))
    return 1 - w / tot


def main():
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_csv(os.path.join(A9, "cms_labels.csv"), index_col=0)
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.4))
    rows = []
    for j, c in enumerate(GEO):
        ax = axes[j]
        e = pd.read_csv(os.path.join(EMB, f"{c}_{SETS[c]}_adjusted_g2000_pc10.csv"), index_col=0)
        car = e.index[e.stage == "Carcinoma"]
        L = lab[lab.cohort == c].reindex(car)
        conf = (L.confident.astype(str).str.lower() == "true").values
        for k in CMS:
            m = (L.nearest == k).values
            ax.scatter(e.loc[car[m & conf], "PC1"], e.loc[car[m & conf], "PC2"], s=14, color=CCOL[k],
                       edgecolor="white", linewidth=0.3)
            ax.scatter(e.loc[car[m & ~conf], "PC1"], e.loc[car[m & ~conf], "PC2"], s=14, facecolor="none",
                       edgecolor=CCOL[k], linewidth=0.9)
        ax.set(xlabel="PC1", ylabel="PC2" if j == 0 else "")
        ax.set_title(f"{c}  (n = {len(car)})", fontweight="bold")
        P = e.loc[car, [f"PC{i}" for i in range(1, 11)]].values
        g = L.nearest.values
        tot = ((P - P.mean(0)) ** 2).sum()
        rows.append(dict(cohort=c, n_carcinoma=len(car), n_not_confident=int((~conf).sum()),
                         r2_cms_pc1_2=r2(P[:, :2], g), r2_cms_pc1_10=r2(P, g),
                         share_variance_in_pc1_2=float(((P[:, :2] - P[:, :2].mean(0)) ** 2).sum() / tot)))
    handles = [Line2D([], [], marker="o", ls="", color=CCOL[k], label=k) for k in CMS] + \
              [Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor="0.4",
                      label="hollow = not confidently classified")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=5, frameon=False, fontsize=7.5)
    (lambda *a, **k: None)("GEO carcinomas on PC1 vs PC2 of the frozen adjusted embedding, coloured by primary CMS call",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    out = pd.DataFrame(rows)
    (lambda *a, **k: None)(os.path.join(ROOT, "results", "source_data", "figS_cms_r2_pcs.csv"), index=False)
    print(out.round(3).to_string(index=False))
    print("saved", os.path.join(OUT, NAME + ".png"))


if __name__ == "__main__":
    main()
