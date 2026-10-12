"""
Figure 5 (CMS subtypes and carcinoma heterogeneity, GEO cohorts). From existing files only; no statistical tests.
(a) per carcinoma (adjusted embedding, all 10 PCs, primary CMS calls): bias-adjusted distance to the centre of all
    carcinomas and to the centre of its own CMS subtype (as addendum 9 loo_dist); checked against addendum 9 D_within.
    Medians/means written to results/source_data/fig5_carcinoma_centre_distances.csv (new file).
    (The former PC1/PC2 panel is now scripts/figS_cms_pc12.py.)
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


def bias_adj_dist(P, groups):
    """As addendum 9 loo_dist: sqrt(n/(n-1)) * ||x - c|| to own-group centroid; NaN for groups with n < 3."""
    d = np.full(len(P), np.nan)
    for gname in np.unique(groups):
        m = groups == gname
        n = m.sum()
        if n < 3:
            continue
        d[m] = np.sqrt(n / (n - 1)) * np.linalg.norm(P[m] - P[m].mean(0), axis=1)
    return d


def main():
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_csv(os.path.join(A9, "cms_labels.csv"), index_col=0)
    prop = pd.read_csv(os.path.join(A9, "cms_proportions.csv")).set_index("cohort")
    dec = pd.read_csv(os.path.join(A9, "dispersion_decomposition.csv"))
    dec = dec[dec["labels"] == "primary"].set_index("cohort")

    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.95], hspace=0.55, wspace=0.42)
    counts, rows = {}, []
    rng = np.random.default_rng(42)                                         # jitter only
    for j, c in enumerate(GEO):
        ax = fig.add_subplot(gs[0, j])
        e = pd.read_csv(os.path.join(EMB, f"{c}_{SETS[c]}_adjusted_g2000_pc10.csv"), index_col=0)
        car = e.index[e.stage == "Carcinoma"]
        L = lab[lab.cohort == c].reindex(car)
        if L.nearest.isna().any():
            raise SystemExit(f"{c}: carcinoma without CMS label")
        conf = L.confident.astype(str).str.lower() == "true"
        counts[c] = (len(car), int((~conf).sum()))
        P = e.loc[car, [f"PC{i}" for i in range(1, 11)]].values
        g = L.nearest.values
        d_all = bias_adj_dist(P, np.zeros(len(P)))                          # centre of all carcinomas
        d_own = bias_adj_dist(P, g)                                         # centre of own CMS subtype
        # check against addendum 9: D_within = D_total * mean(d_own) / mean(d_all) (same adenoma reference)
        Dw = dec.loc[c, "D_total"] * np.nanmean(d_own) / np.nanmean(d_all)
        if abs(Dw - dec.loc[c, "D_within"]) > 1e-9:
            raise SystemExit(f"{c}: D_within check failed ({Dw} vs {dec.loc[c, 'D_within']})")
        bp = ax.boxplot([d_all, d_own[~np.isnan(d_own)]], patch_artist=True, widths=0.55, showfliers=False,
                        medianprops=dict(color="black", lw=1))
        bp["boxes"][0].set(facecolor="#c0392b", alpha=0.25, edgecolor="#c0392b")
        bp["boxes"][1].set(facecolor="0.8", alpha=0.5, edgecolor="0.4")
        ax.scatter(1 + rng.uniform(-0.17, 0.17, len(d_all)), d_all, s=6, color="#c0392b", alpha=0.7, edgecolor="none",
                   zorder=3)
        jit = rng.uniform(-0.17, 0.17, len(d_own))
        for k in CMS:
            m = g == k
            ax.scatter(2 + jit[m], d_own[m], s=6, color=CCOL[k], alpha=0.85, edgecolor="none", zorder=3)
        ax.set_xticks([1, 2], ["all carcinomas", "own CMS\nsubtype"], fontsize=7.5)
        ax.set_ylim(0, max(np.nanmax(d_all), np.nanmax(d_own)) * 1.06)
        ax.set_title(f"{c}  (n = {len(car)})", fontweight="bold")
        ratio = np.nanmedian(d_own) / np.median(d_all)
        ax.text(0.98, 0.97, f"median ratio {ratio:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=7,
                color="0.3")
        rows.append(dict(cohort=c, n=len(car), median_all=float(np.median(d_all)), median_own=float(np.nanmedian(d_own)),
                         mean_all=float(np.mean(d_all)), mean_own=float(np.nanmean(d_own)), median_ratio=float(ratio)))
        if j == 0:
            ax.set_ylabel("carcinoma distance to centre\n(10 PCs, bias-adjusted)")
            ax.text(-0.42, 1.13, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.text(-0.2, 1.13, "Carcinoma spread around the centre of all carcinomas vs around the centre of their "
                    "own CMS subtype (adjusted, 10 PCs)", transform=ax.transAxes, fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", color="#c0392b", label="carcinomas (all)")] + \
              [Line2D([], [], marker="o", ls="", color=CCOL[k], label=k) for k in CMS]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.535), ncol=5, frameon=False, fontsize=7.5)
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "results", "source_data", "fig5_carcinoma_centre_distances.csv"),
                              index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print("check: D_within recomputed from panel a distances = addendum 9 (all cohorts)")

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
