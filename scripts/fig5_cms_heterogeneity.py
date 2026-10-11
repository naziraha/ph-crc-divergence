"""
Figure 5 (CMS subtypes and carcinoma heterogeneity, GEO cohorts). From existing files only; no statistical tests.
(a) carcinomas on PC1 vs PC2 of the frozen adjusted embedding, coloured by the primary CMS call (nearest template,
    addendum 9); hollow = not confidently classified (unclassified under the confident rule).
(b) CMS proportions per cohort: primary calls (CMS1-4) and confident calls (CMS1-4 + unclassified).
(c) share f of the excess carcinoma dispersion explained by CMS (primary labels), with f_rand95; GSE41258 not applicable.
Inputs: results/addendum9/cms_labels.csv, cms_proportions.csv, dispersion_decomposition.csv;
results/step03/embeddings/<cohort>_<set>_adjusted_g2000_pc10.csv.

Usage: python scripts/fig5_cms_heterogeneity.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A9 = os.path.join(ROOT, "results", "addendum9")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "paper", "final", "main_figures")
NAME = "Fig5_CMS_heterogeneity"
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
SETS = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
CMS = ["CMS1", "CMS2", "CMS3", "CMS4"]
CCOL = {"CMS1": "#e79e00", "CMS2": "#0071b1", "CMS3": "#ca78a6", "CMS4": "#009f74", "unclassified": "#b0b0b0"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def main():
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_csv(os.path.join(A9, "cms_labels.csv"), index_col=0)
    prop = pd.read_csv(os.path.join(A9, "cms_proportions.csv")).set_index("cohort")
    dec = pd.read_csv(os.path.join(A9, "dispersion_decomposition.csv"))
    dec = dec[dec["labels"] == "primary"].set_index("cohort")

    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.95], hspace=0.55, wspace=0.42)
    counts = {}
    for j, c in enumerate(GEO):
        ax = fig.add_subplot(gs[0, j])
        e = pd.read_csv(os.path.join(EMB, f"{c}_{SETS[c]}_adjusted_g2000_pc10.csv"), index_col=0)
        car = e.index[e.stage == "Carcinoma"]
        L = lab[lab.cohort == c].reindex(car)
        if L.nearest.isna().any():
            raise SystemExit(f"{c}: carcinoma without CMS label")
        conf = L.confident.astype(str).str.lower() == "true"
        counts[c] = (len(car), int((~conf).sum()))
        for k in CMS:
            m = (L.nearest == k).values
            ax.scatter(e.loc[car[m & conf.values], "PC1"], e.loc[car[m & conf.values], "PC2"], s=14, color=CCOL[k],
                       edgecolor="white", linewidth=0.3, label=k)
            ax.scatter(e.loc[car[m & ~conf.values], "PC1"], e.loc[car[m & ~conf.values], "PC2"], s=14,
                       facecolor="none", edgecolor=CCOL[k], linewidth=0.9)
        ax.set(xlabel="PC1", ylabel="PC2" if j == 0 else "")
        ax.set_title(f"{c}  (n = {len(car)})", fontweight="bold")
        if j == 0:
            ax.text(-0.3, 1.13, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.text(-0.1, 1.13, "Carcinomas on PC1 vs PC2 (frozen adjusted embedding), coloured by primary CMS call",
                    transform=ax.transAxes, fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", color=CCOL[k], label=k) for k in CMS] + \
              [Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor="0.4",
                      label="hollow = not confidently classified")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.535), ncol=5, frameon=False, fontsize=7.5)

    # panel b: stacked proportions (primary and confident) per cohort
    ax = fig.add_subplot(gs[1, 0:2])
    x, labels = [], []
    pos = 0
    for c in GEO:
        for kind in ("primary", "confident"):
            bottom = 0.0
            keys = CMS + (["unclassified"] if kind == "confident" else [])
            for k in keys:
                v = prop.loc[c, f"{kind}_{k}"]
                ax.bar(pos, v, bottom=bottom, color=CCOL[k], width=0.8, edgecolor="white", linewidth=0.4)
                bottom += v
            x.append(pos)
            labels.append("prim." if kind == "primary" else "conf.")
            pos += 1
        ax.text(pos - 1.5, 1.04, f"{c}\n(n = {int(prop.loc[c, 'n'])})", ha="center", fontsize=7.5, fontweight="bold")
        pos += 0.6
    ax.set_xticks(x, labels, fontsize=7)
    ax.set(ylim=(0, 1.15), ylabel="proportion of carcinomas", yticks=[0, 0.25, 0.5, 0.75, 1.0])
    ax.legend(handles=[Patch(color=CCOL[k], label=k) for k in CMS + ["unclassified"]], fontsize=7, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=5)
    ax.text(-0.13, 1.12, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.set_title("CMS composition: primary vs confident calls", pad=24)

    # panel c: f vs f_rand95
    ax = fig.add_subplot(gs[1, 2:4])
    for i, c in enumerate(GEO):
        r = dec.loc[c]
        if r["call"] == "not informative":
            ax.text(i, 0.30, "not applicable\n(no excess carcinoma\ndispersion after\nadjustment)", ha="center",
                    va="bottom", fontsize=6.8, color="0.35")
            continue
        above = r.f > r.f_rand95
        ax.bar(i, r.f, width=0.6, color="#5b5b5b" if above else "#c8c8c8")
        ax.plot([i - 0.38, i + 0.38], [r.f_rand95, r.f_rand95], color="#c0392b", lw=2)
        ax.text(i, r.f + 0.02, f"f = {r.f:.2f}", ha="center", fontsize=7.5)
    ax.plot([], [], color="#c0392b", lw=2, label="random-label 95th percentile (f_rand95)")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.set_xticks(range(len(GEO)), GEO, fontsize=7.5)
    ax.set(ylim=(0, 1.0), xlim=(-0.6, len(GEO) - 0.4), ylabel="share of excess carcinoma dispersion\nexplained by CMS (f; adjusted data)")
    ax.text(-0.18, 1.12, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.set_title("Excess carcinoma dispersion explained by CMS", pad=24)

    fig.suptitle("CMS subtypes and carcinoma heterogeneity in the GEO cohorts", fontsize=11, fontweight="bold", y=0.995)
    fig.text(0.5, 0.0, "CMS calls from bulk expression partly reflect stromal content (especially CMS4). Panels a-b are "
             "descriptive; f and f_rand95 from addendum 9 (999 random-label permutations).", ha="center", fontsize=7.5,
             color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    for c in GEO:
        print(c, "carcinomas", counts[c][0], "not confident", counts[c][1])
    print("saved", os.path.join(OUT, NAME + ".png"))


if __name__ == "__main__":
    main()
