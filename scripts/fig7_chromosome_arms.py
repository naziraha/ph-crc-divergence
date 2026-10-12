"""
Figure 7 (chromosome arms; addendum 11), from existing files only; no statistical tests.
(a) per arm (genome order) and GEO cohort: median ln V of the arm's genes (V = MAD carcinoma / MAD adenoma, adjusted
    expression; 117 genes excluded; arms >= 20 genes) from arm_variability_ranking.csv (median_logV).
(b) 20q arm score per sample (adjusted, 117 genes excluded), adenoma vs carcinoma, per GEO cohort; adenomas grey,
    carcinomas by primary CMS; Brown-Forsythe p and median difference from arm_scores_tests.csv.
(c) TCGA: 20q expression score (tcga_20q_score.csv) vs mean GISTIC2 copy number of 20q genes per tumour (raw GISTIC2
    file, same computation as addendum 11 A4); Spearman rho checked against tcga_expression_vs_copynumber.csv.
(d) share f of excess carcinoma dispersion explained by 20q, by the four arms (20q+20p+13q+8q) and by CMS, with f_rand95
    (B3_dispersion_decomposition.csv; CMS f also in addendum 9), informative cohorts only.

Usage: python scripts/fig7_chromosome_arms.py
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
R = os.path.join(ROOT, "results")
A11 = os.path.join(R, "addendum11")
OUT = os.path.join(ROOT, "paper", "final", "main_figures")
NAME = "Fig7_chromosome_arms"
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
KCOL = {"GSE41258": "#332288", "GSE20916": "#44aa99", "GSE37364": "#c9a227", "GSE117606": "#cc6677"}   # cohort colours
CCOL = {"CMS1": "#e79e00", "CMS2": "#0071b1", "CMS3": "#ca78a6", "CMS4": "#009f74"}
HIGH = ["20q", "20p", "13q", "8q"]
COMP = ["7p", "7q"]
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def armkey(a):
    m = re.match(r"(\d+)([pq])", a)
    return (int(m.group(1)), m.group(2))


def tcga_cn():
    import addendum11_chr20 as A11m
    arms, _ = A11m.arm_map()
    s = pd.read_csv(os.path.join(A11, "partA_chr", "tcga_20q_score.csv"), index_col=0)
    cn = pd.read_csv(os.path.join(ROOT, "data", "raw", "TCGA", "TCGA.COAD.sampleMap_Gistic2_all_data_by_genes.gz"),
                     sep="\t", index_col=0)
    common = [x for x in s.index if x in cn.columns]
    g20 = [x for x in cn.index if arms.get(x) == "20q"]
    return s.loc[common, "TCGA_20q_score"], cn.loc[g20, common].mean(axis=0), s.loc[common, "CMS"]


def main():
    os.makedirs(OUT, exist_ok=True)
    rk = pd.read_csv(os.path.join(A11, "partA_chr", "arm_variability_ranking.csv"))
    tests = pd.read_csv(os.path.join(A11, "partA_chr", "arm_scores_tests.csv"))
    tests = tests[(tests.version == "adjusted") & (tests.genes == "excl117") & (tests.arm == "20q")].set_index("cohort")
    tcv = pd.read_csv(os.path.join(A11, "partA_chr", "tcga_expression_vs_copynumber.csv")).set_index("arm")
    b3 = pd.read_csv(os.path.join(A11, "partB_differentiation", "B3_dispersion_decomposition.csv")).set_index("cohort")
    a9 = pd.read_csv(os.path.join(R, "addendum9", "dispersion_decomposition.csv"))
    a9 = a9[a9["labels"] == "primary"].set_index("cohort")
    armsc = pd.read_csv(os.path.join(R, "source_data", "arm_scores_GEO.csv"))
    armsc = armsc[(armsc.version == "adjusted") & (armsc.genes == "excl117")]
    cms = pd.read_csv(os.path.join(R, "addendum9", "cms_labels.csv"), index_col=0)["nearest"]

    fig = plt.figure(figsize=(11, 12.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.9, 0.85, 1.0], hspace=0.55, wspace=0.38, top=0.94, bottom=0.04)

    # ---- a: arm variability, genome order
    ax = fig.add_subplot(gs[0, :])
    arms = sorted(rk.arm.unique(), key=armkey)
    xpos = {a: i for i, a in enumerate(arms)}
    for a in HIGH:
        ax.axvspan(xpos[a] - 0.45, xpos[a] + 0.45, color="#f4c7c3", zorder=0)
    for a in COMP:
        ax.axvspan(xpos[a] - 0.45, xpos[a] + 0.45, color="#e3e3e3", zorder=0)
    off = dict(zip(GEO, (-0.24, -0.08, 0.08, 0.24)))
    for c in GEO:
        d = rk[rk.cohort == c]
        ax.scatter([xpos[a] + off[c] for a in d.arm], d.median_logV, s=16, color=KCOL[c], edgecolor="white",
                   linewidth=0.3, label=c, zorder=3)
    ax.axhline(0, color="0.4", lw=0.8, ls="--")
    ax.set_xticks(range(len(arms)), arms, fontsize=6.5, rotation=90)
    for t in ax.get_xticklabels():
        if t.get_text() in HIGH:
            t.set_color("#c0392b"); t.set_fontweight("bold")
        elif t.get_text() in COMP:
            t.set_color("0.35"); t.set_fontweight("bold")
    ax.set_xlim(-0.7, len(arms) - 0.3)
    ax.set_ylabel("median ln V of the arm's genes\n(V = MAD carcinoma ÷ MAD adenoma)")
    ax.legend(fontsize=7, frameon=False, ncol=4, loc="lower right", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Increase in gene variability from adenoma to carcinoma, per chromosome arm\n(adjusted data; 117 selected "
                 "genes excluded)", loc="left", fontsize=8.5)
    ax.text(0.005, 0.95, "red band: 20q, 20p, 13q, 8q; grey band: 7p, 7q (comparators)", transform=ax.transAxes,
            fontsize=7, color="0.3", va="top")
    ax.text(-0.07, 1.06, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")

    # ---- b: 20q arm score per sample
    rng = np.random.default_rng(42)                                   # jitter only
    for j, c in enumerate(GEO):
        ax = fig.add_subplot(gs[1, j])
        d = armsc[armsc.cohort == c].set_index("sample")
        data = [d.loc[d.stage == s, "arm_20q"].values for s in ("Adenoma", "Carcinoma")]
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False, medianprops=dict(color="black", lw=1))
        bp["boxes"][0].set(facecolor="#e0a100", alpha=0.15, edgecolor="#e0a100")
        bp["boxes"][1].set(facecolor="#c0392b", alpha=0.15, edgecolor="#c0392b")
        ax.scatter(1 + rng.uniform(-0.17, 0.17, len(data[0])), data[0], s=5, color="0.55", alpha=0.8, edgecolor="none",
                   zorder=3)
        car = d.index[d.stage == "Carcinoma"]
        ax.scatter(2 + rng.uniform(-0.17, 0.17, len(car)), d.loc[car, "arm_20q"], s=5,
                   c=[CCOL[cms[s]] for s in car], alpha=0.85, edgecolor="none", zorder=3)
        t = tests.loc[c]
        ax.set_xticks([1, 2], [f"Adenoma\n({len(data[0])})", f"Carcinoma\n({len(car)})"], fontsize=7)
        ax.set_title(c, fontweight="bold", loc="left", pad=16)
        ax.text(0.0, 1.02, f"BF p = {t.p_BF:.3g}; median diff = {t.median_diff:+.2f}", transform=ax.transAxes,
                fontsize=6.8, color="0.25")
        if j == 0:
            ax.set_ylabel("20q arm score\n(mean gene z-score, adjusted)")
            ax.text(-0.42, 1.2, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")
    hb = [Line2D([], [], marker="o", ls="", color="0.55", label="adenoma")] + \
         [Line2D([], [], marker="o", ls="", color=CCOL[k], label=f"carcinoma, {k}") for k in CCOL]
    fig.legend(handles=hb, loc="upper center", bbox_to_anchor=(0.5, 0.375), ncol=5, frameon=False, fontsize=7)

    # ---- c: TCGA 20q expression vs copy number
    ax = fig.add_subplot(gs[2, 0:2])
    expr, cnm, tcms = tcga_cn()
    rho = spearmanr(expr, cnm)[0]
    if abs(rho - tcv.loc["20q", "rho"]) > 1e-9 or len(expr) != int(tcv.loc["20q", "n_tumours"]):
        raise SystemExit(f"TCGA rho/n mismatch: {rho} vs {tcv.loc['20q', 'rho']}")
    ax.scatter(cnm, expr, s=9, color="#4d4d4d", alpha=0.55, edgecolor="none")
    ax.set(xlabel="mean 20q copy number per tumour (GISTIC2)", ylabel="20q expression score")
    ax.set_title("TCGA validation (DNA copy number)", fontweight="bold")
    ax.text(0.03, 0.95, f"Spearman ρ = {rho:.2f}\nn = {len(expr)} tumours", transform=ax.transAxes, va="top", fontsize=8)
    ax.text(-0.18, 1.06, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

    # ---- d: shares of excess dispersion
    ax = fig.add_subplot(gs[2, 2:4])
    inf = [c for c in GEO if bool(b3.loc[c, "informative"])]
    items = [("20q", "chr20q_f", "chr20q_f_rand95", "#c0392b"), ("20q+20p+13q+8q", "chr4_f", "chr4_f_rand95", "#e08a7a"),
             ("CMS", "f_CMS", None, "#5b5b5b")]
    w = 0.26
    for i, c in enumerate(inf):
        if abs(b3.loc[c, "f_CMS"] - a9.loc[c, "f"]) > 1e-9:
            raise SystemExit(f"{c}: f_CMS differs between addendum 11 and addendum 9")
        for k, (lab, fc, rc, col) in enumerate(items):
            x = i + (k - 1) * w
            f = b3.loc[c, fc]
            r95 = b3.loc[c, rc] if rc else a9.loc[c, "f_rand95"]
            ax.bar(x, f, width=w * 0.9, color=col, label=lab if i == 0 else None)
            ax.plot([x - w * 0.45, x + w * 0.45], [r95, r95], color="black", lw=1.6)
            ax.text(x, f + (0.03 if f >= 0 else -0.03), f"{f:.2f}", ha="center", va="bottom" if f >= 0 else "top",
                    fontsize=6.3)
    ax.plot([], [], color="black", lw=1.6, label="f_rand95 (random-score 95th pct.)")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(range(len(inf)), inf, fontsize=7.5)
    ax.set(ylim=(-0.95, 1.0), ylabel="share of excess carcinoma dispersion\nexplained (f; adjusted data)")
    ax.legend(fontsize=6.8, frameon=False, loc="upper left", ncol=2)
    ax.set_title("Excess carcinoma dispersion explained by arm scores and CMS")
    ax.text(0.99, 0.02, "GSE41258 not informative (no excess dispersion)", transform=ax.transAxes, ha="right",
            fontsize=6.5, color="0.4")
    ax.text(-0.2, 1.06, "d", transform=ax.transAxes, fontsize=13, fontweight="bold")

    fig.suptitle("Chromosome-arm gains as a source of carcinoma heterogeneity (expression-based arm scores; descriptive "
                 "panels a-c)", fontsize=11, fontweight="bold", y=0.995)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("arms plotted:", len(arms), "| 20q rank per cohort:",
          {c: int(rk[(rk.cohort == c) & (rk.arm == "20q")]["rank"].iloc[0]) for c in GEO})
    print(f"check: TCGA rho {rho:.4f} = file; n {len(expr)}; f_CMS addendum 11 = addendum 9")
    print("saved", os.path.join(OUT, NAME + ".png"))


if __name__ == "__main__":
    main()
