# Figure script for the Nature Computational Science submission (NCS numbering); symbols italic via _italic_symbols.py.
# NCS VERSION (2 Nov 2026), after Prof Peng's comments: copy of fig6_variable_genes_revised.py for Fig. 4.
# Changes: larger fonts; heatmap gene labels enlarged and spread with leader lines so they do not overlap;
# panel b shows the stored Brown-Forsythe p-value (adenoma vs carcinoma, per cohort, from variability_<cohort>.csv)
# as a bracket, as in Prof Peng's figure style; shorter panel b heading. No new test. Original scripts untouched.
# REVISED COPY (26 Oct 2026): same data and plotting as the original script; only the wording
# inside the figure is cleaned for journal style (no overall title, no internal terms such as
# "step", "addendum" or "frozen", no footnotes that repeat the legend). Outputs go to *_revised folders
# with "_revised" in the file name, so the original figures are not overwritten. No CSV is written.
"""
Figure 6 (the 122 replicated 'more variable in carcinoma' genes, addendum 9 rule with current HGNC symbols), GEO cohorts.
Descriptive; no tests. Inputs from results/hgnc/ (scripts/hgnc_crosscohort.py; Supplementary Methods S7).
Expression version: composition-ADJUSTED (as used by addendum 9 to find the genes; per_sample 'adj_' columns).
Groups (addendum 9): 'no mean shift' = median |mean difference carcinoma - adenoma| < 0.5 log2 across cohorts
(no_mean_shift); otherwise 'lower' / 'higher' in carcinoma by the sign of the median signed mean difference
(recomputed here from variability_<cohort>.csv and checked against med_mean_diff in the annotated file).
(a) heatmaps, rows = genes (grouped; ordered by median_V within group), columns = N/A/C per cohort:
    left = mean per stage of per-gene z-scores (z over the cohort's N+A+C samples);
    right = within-stage SD per stage divided by the gene's SD over the cohort's N+A+C samples. Grey = not measured.
(b) four example genes (fixed rule: top-2 median_V in 'lower', top-1 in 'higher', top-1 in 'no mean shift'):
    box plots by stage per cohort; carcinoma points coloured by primary CMS call.
Inputs: results/addendum9/variability_replicated_annotated.csv, variability_<cohort>.csv,
results/source_data/per_sample_<cohort>.csv.

Usage: python scripts/fig4_variable_genes.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _italic_symbols  # noqa: F401  italic statistical symbols, Nature style (Zira, 3 Nov 2026)  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A9 = os.path.join(ROOT, "results", "hgnc", "addendum9")
SD = os.path.join(ROOT, "results", "hgnc", "source_data")
OUT = os.path.join(ROOT, "outputs", "main_figures")
NAME = "Fig4_variable_genes"
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
ST3 = ["Normal", "Adenoma", "Carcinoma"]
SCOL = {"Normal": "#2a78d6", "Adenoma": "#e0a100", "Carcinoma": "#c0392b"}
CCOL = {"CMS1": "#e79e00", "CMS2": "#0071b1", "CMS3": "#ca78a6", "CMS4": "#009f74"}
GROUPS = ["lower in carcinoma", "higher in carcinoma", "no mean shift"]
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False})


def groups():
    a = pd.read_csv(os.path.join(A9, "variability_replicated_annotated.csv")).set_index("gene")
    md = {}
    for g in a.index:
        vals = []
        for c in GEO:
            v = pd.read_csv(os.path.join(A9, f"variability_{c}.csv"), index_col=0) if c not in md else md[c]
            md[c] = v
        md_g = [md[c].loc[g, "mean_diff"] for c in GEO if g in md[c].index]
        a.loc[g, "signed_md"] = float(np.median(md_g))
        a.loc[g, "abs_md"] = float(np.median(np.abs(md_g)))
    if not np.allclose(a.signed_md, a.med_mean_diff, atol=1e-9) or not np.allclose(a.abs_md, a.median_abs_mean_diff, atol=1e-9):
        raise SystemExit("mean-difference recomputation does not match the annotated file")
    a["group"] = np.where(a.no_mean_shift, "no mean shift", np.where(a.signed_md < 0, "lower in carcinoma", "higher in carcinoma"))
    a = a.sort_values("median_V", ascending=False)
    order = [g for grp in GROUPS for g in a.index[a.group == grp]]
    return a, order


def main():
    os.makedirs(OUT, exist_ok=True)
    a, order = groups()
    ex = {c: pd.read_csv(os.path.join(SD, f"per_sample_{c}.csv"), index_col=0) for c in GEO}
    cols = [(c, s) for c in GEO for s in ST3]
    M = np.full((len(order), len(cols)), np.nan)
    S = np.full((len(order), len(cols)), np.nan)
    for j, (c, s) in enumerate(cols):
        d = ex[c][ex[c].stage.isin(ST3)]
        for i, g in enumerate(order):
            col = f"adj_{g}"
            if col not in d or d[col].isna().all():
                continue
            x = d[col]
            z = (x - x.mean()) / x.std(ddof=1)
            M[i, j] = z[d.stage == s].mean()
            S[i, j] = x[d.stage == s].std(ddof=1) / x.std(ddof=1)

    fig = plt.figure(figsize=(13, 18))
    top = fig.add_gridspec(1, 2, left=0.16, right=0.97, top=0.95, bottom=0.50, wspace=0.62)
    cmap_m = plt.get_cmap("RdBu_r").copy(); cmap_m.set_bad("#d9d9d9")
    cmap_s = plt.get_cmap("magma_r").copy(); cmap_s.set_bad("#d9d9d9")
    for k, (mat, cmap, vmin, vmax, ttl, lab) in enumerate((
            (M, cmap_m, -1.5, 1.5, "Mean expression per stage", "mean z-score"),
            (S, cmap_s, 0.3, 1.6, "Within-stage spread per stage", "SD ÷ gene SD in cohort"))):
        ax = fig.add_subplot(top[0, k])
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="none")
        ax.set_xticks(range(len(cols)), [s[0] for _, s in cols], fontsize=9)
        for t, (_, s) in zip(ax.get_xticklabels(), cols):
            t.set_color(SCOL[s])
        for cj, c in enumerate(GEO):
            ax.text(cj * 3 + 1, -1.5, c, ha="center", va="bottom", fontsize=9.5, fontweight="bold", rotation=0)
            if cj:
                ax.axvline(cj * 3 - 0.5, color="white", lw=1.5)
        # fixed label set: all 'higher' genes + top 5 by median V in 'lower' and in 'no mean shift'
        labset = set(a.index[a.group == "higher in carcinoma"]) | set(a.index[a.group == "lower in carcinoma"][:5]) | \
            set(a.index[a.group == "no mean shift"][:5])
        idx = [i for i, g in enumerate(order) if g in labset]
        ax.set_yticks([])
        # labels spread with leader lines (minimum gap in rows) so they stay readable
        gap = 3.0
        ys = [float(i) for i in idx]
        for _ in range(200):
            moved = False
            for q in range(1, len(ys)):
                if ys[q] - ys[q - 1] < gap:
                    d = (gap - (ys[q] - ys[q - 1])) / 2
                    ys[q - 1] -= d; ys[q] += d; moved = True
            if not moved:
                break
        for i, yl in zip(idx, ys):
            ax.annotate(order[i], xy=(-0.5, i), xytext=(-2.3, yl), textcoords="data", ha="right", va="center",
                        fontsize=9, annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4", shrinkA=0, shrinkB=0))
        bounds = np.cumsum([(a.group == g).sum() for g in GROUPS])
        for b in bounds[:-1]:
            ax.axhline(b - 0.5, color="black", lw=0.8)
        start = 0
        for g, b in zip(GROUPS, bounds):
            if k == 0:                                                     # group labels once, left heatmap
                # dashed bracket marking the row range of each group (Zira, 2 Nov)
                y0, y1 = 1 - start / len(order) - 0.004, 1 - b / len(order) + 0.004
                xb, cap = -0.55, 0.025
                ax.plot([xb, xb], [y0, y1], transform=ax.transAxes, clip_on=False,
                        color="0.25", lw=1.0, ls=(0, (4, 2.5)))
                for yy in (y0, y1):
                    ax.plot([xb, xb + cap], [yy, yy], transform=ax.transAxes, clip_on=False, color="0.25", lw=1.0)
                ymid = 1 - (start + b) / 2 / len(order)
                if b - start >= 15:
                    ax.text(-0.62, ymid, f"{g}\n(n = {b - start})",
                            transform=ax.transAxes, ha="center", va="center", rotation=90, fontsize=10)
                else:                                                      # small group: horizontal label fits its range
                    ax.text(-0.57, ymid, g.replace(" in ", " in\n") + f"\n(n = {b - start})",
                            transform=ax.transAxes, ha="right", va="center", fontsize=9.5, linespacing=1.1)
            start = b
        ax.set_title(ttl, pad=22)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label(lab, fontsize=9.5)
        cb.ax.tick_params(labelsize=8.5)
        if k == 0:
            ax.text(-0.75, 1.05, "a", transform=ax.transAxes, fontsize=15, fontweight="bold")

    # panel b: example genes (fixed rule)
    ex_genes = (list(a.index[a.group == "lower in carcinoma"][:2]) + list(a.index[a.group == "higher in carcinoma"][:1])
                + list(a.index[a.group == "no mean shift"][:1]))
    bot = fig.add_gridspec(4, 4, left=0.08, right=0.98, top=0.43, bottom=0.05, hspace=0.7, wspace=0.32)
    rng = np.random.default_rng(42)                                      # jitter only
    for i, g in enumerate(ex_genes):
        for j, c in enumerate(GEO):
            ax = fig.add_subplot(bot[i, j])
            d = ex[c][ex[c].stage.isin(ST3)]
            col = f"adj_{g}"
            if col not in d or d[col].isna().all():
                ax.text(0.5, 0.5, "not measured", ha="center", va="center", transform=ax.transAxes, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
            else:
                data = [d.loc[d.stage == s, col].dropna().values for s in ST3]
                bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False, medianprops=dict(color="black", lw=1))
                for p, s in zip(bp["boxes"], ST3):
                    p.set(facecolor=SCOL[s], alpha=0.15, edgecolor=SCOL[s])
                for k, s in enumerate(ST3):
                    sub = d[d.stage == s]
                    y = sub[col].values
                    xj = k + 1 + rng.uniform(-0.17, 0.17, len(y))
                    if s == "Carcinoma":
                        colr = [CCOL.get(v, "0.5") for v in sub.CMS_nearest]
                    else:
                        colr = "0.55"                                  # normal / adenoma points grey
                    ax.scatter(xj, y, s=7, c=colr, alpha=0.8, edgecolor="none", zorder=3)
                ax.set_xticks([1, 2, 3], ["N", "A", "C"], fontsize=9.5)
                V = a.loc[g, f"V_{c}"]
                p = pd.read_csv(os.path.join(A9, f"variability_{c}.csv"), index_col=0).p.get(g, np.nan)
                ymax = max(np.concatenate(data)); ymin = min(np.concatenate(data)); rngy = ymax - ymin
                if pd.notna(p):
                    yb = ymax + 0.08 * rngy
                    ax.plot([2, 2, 3, 3], [yb - 0.04 * rngy, yb, yb, yb - 0.04 * rngy], color="black", lw=0.9)
                    star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                    ax.text(2.5, yb + 0.005 * rngy, star, ha="center", va="bottom", fontsize=12 if p < 0.05 else 9,
                            fontweight="bold" if p < 0.05 else "normal")  # stars only (Zira, 2 Nov)
                    ax.set_ylim(ymin - 0.05 * rngy, ymax + 0.32 * rngy)
                ax.set_title("V not tested" if pd.isna(V) else f"V = {V:.2f}", fontsize=9, color="0.25", loc="right", pad=2)
                ax.tick_params(axis="y", labelsize=9)
            if i == 0:
                ax.set_title(c, fontweight="bold", loc="left", pad=10, fontsize=11)
            if j == 0:
                ax.set_ylabel(f"{g}\n({a.loc[g, 'group']})\nadjusted log2", fontsize=9.5)
            if i == 0 and j == 0:
                ax.text(-0.55, 1.45, "b", transform=ax.transAxes, fontsize=15, fontweight="bold")
                ax.text(-0.33, 1.45, "Example genes; carcinoma points coloured by CMS", transform=ax.transAxes,
                        fontsize=11)
    handles = [Line2D([], [], marker="o", ls="", color="0.55", label="normal / adenoma sample")] + \
              [Line2D([], [], marker="o", ls="", color=CCOL[k], label=f"carcinoma, {k}") for k in CCOL]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=5, frameon=False, fontsize=10)
    (lambda *a, **k: None)("The 122 genes more variable in carcinoma than in adenoma (addendum 9 rule; current HGNC symbols; adjusted data; descriptive)",
                 fontsize=11, fontweight="bold", y=0.995)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    meas = {c: int(sum(f"adj_{g}" in ex[c] and ex[c][f"adj_{g}"].notna().any() for g in order)) for c in GEO}
    print("groups:", {g: int((a.group == g).sum()) for g in GROUPS})
    print("example genes:", [(g, a.loc[g, "group"], round(a.loc[g, "median_V"], 2)) for g in ex_genes])
    print("genes measured per cohort:", meas)
    print("check: signed and absolute median mean differences recomputed = annotated file")
    print("saved", os.path.join(OUT, NAME + ".png"))


if __name__ == "__main__":
    main()
