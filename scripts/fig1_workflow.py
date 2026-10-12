"""
Figure 1: study design and workflow schematic (no data analysis, no tests). Sample counts are read from
paper/final/main_tables/Table1_cohorts.csv; TCGA tumour counts from results/addendum11/partA_chr/tcga_20q_score.csv and
tcga_expression_vs_copynumber.csv; pipeline settings from results/step05/summary.txt and Table 2 footnotes; the last
deviation number from docs/analysis_deviations.md. Drawn in matplotlib for later refinement (e.g. BioRender).

Usage: python scripts/fig1_workflow.py
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "paper", "final", "main_figures")
NAME = "Fig1_workflow"
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
SCOL = {"Normal": "#2a78d6", "Adenoma": "#e0a100", "Carcinoma": "#c0392b", "Metastasis": "#7b4fb3"}
INK, MUTED, FILL, EDGE = "#1f1f1f", "#555555", "#f4f4f4", "#9a9a9a"


def box(ax, x, y, w, h, text, fc=FILL, ec=EDGE, fs=7.5, bold=False, color=INK, lw=0.9, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012", fc=fc, ec=ec, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color,
            fontweight="bold" if bold else "normal", linespacing=1.25)


def arrow(ax, x0, y0, x1, y1, color=MUTED):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9, color=color, lw=1.0,
                                 shrinkA=0, shrinkB=0))


def numbers():
    t1 = pd.read_csv(os.path.join(ROOT, "paper", "final", "main_tables", "Table1_cohorts.csv")).set_index("Cohort")
    tc = pd.read_csv(os.path.join(ROOT, "results", "addendum11", "partA_chr", "tcga_20q_score.csv"))
    tcv = pd.read_csv(os.path.join(ROOT, "results", "addendum11", "partA_chr", "tcga_expression_vs_copynumber.csv"))
    dev = open(os.path.join(ROOT, "docs", "analysis_deviations.md"), encoding="utf-8").read()
    dmax = max(int(x) for x in re.findall(r"^\| D(\d+) \|", dev, flags=re.M))
    s5 = open(os.path.join(ROOT, "results", "step05", "summary.txt"), encoding="utf-8").read()
    assert "999 permutations" in s5 and "k = 20" in s5 and "m = 80% of smallest stage" in s5
    total = int(t1["Total samples"].astype(int).sum())
    return t1, len(tc), int(tcv.set_index("arm").loc["20q", "n_tumours"]), dmax, total


def main():
    os.makedirs(OUT, exist_ok=True)
    t1, n_tcga, n_tcga_cn, dmax, total = numbers()
    fig = plt.figure(figsize=(12, 9.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---------------- panel a: study design
    ax.text(0.015, 0.975, "a", fontsize=14, fontweight="bold")
    ax.text(0.045, 0.975, f"Study design: four GEO cohorts (n = {total}) analysed separately; TCGA-COAD for copy-number "
            "validation", fontsize=10, va="baseline")
    x0, w, gap, y0, h = 0.03, 0.163, 0.022, 0.735, 0.205
    for i, c in enumerate(t1.index):
        x = x0 + i * (w + gap)
        r = t1.loc[c]
        ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.006,rounding_size=0.012", fc="white", ec=EDGE, lw=0.9))
        ax.text(x + w / 2, y0 + h - 0.022, c, ha="center", fontsize=9, fontweight="bold")
        plat = r.Platform.split(" (")[0]
        ax.text(x + w / 2, y0 + h - 0.045, f"{plat} · set '{r['Sample set']}' · n = {r['Total samples']}", ha="center",
                fontsize=6.6, color=MUTED)
        st = [s for s in STAGES if str(r[s]) not in ("–", "-", "nan")]
        nmax = 181
        for k, s in enumerate(st):
            yy = y0 + h - 0.075 - k * 0.032
            n = int(r[s])
            ax.text(x + 0.008, yy, s, fontsize=7, va="center", color=SCOL[s], fontweight="bold")
            ax.add_patch(plt.Rectangle((x + 0.068, yy - 0.009), 0.072 * n / nmax, 0.018, fc=SCOL[s], ec="none", alpha=0.85))
            ax.text(x + 0.068 + 0.072 * n / nmax + 0.004, yy, str(n), fontsize=7, va="center")
    x = x0 + 4 * (w + gap)
    box(ax, x, y0, 0.205, h, f"TCGA-COAD\nprimary tumours (RNA-seq)\nn = {n_tcga} ({n_tcga_cn} with GISTIC2 copy number)\n\n"
        "validation of 20q copy number only\n(Figure 7c)", fc="#fbfbfb", ec=EDGE, fs=7.6, ls="--")
    ax.text(0.03, 0.712, "Cohorts are never merged and no batch correction is applied; results are compared across cohorts.",
            fontsize=7, color=MUTED, style="italic")

    # ---------------- panel b: pipeline
    ax.text(0.015, 0.665, "b", fontsize=14, fontweight="bold")
    ax.text(0.045, 0.665, "Frozen, pre-registered pipeline (each cohort separately)", fontsize=10, va="baseline")
    steps = ["Quality control\n\none sample per\npatient per stage;\ntissue-identity\nexclusions",
             "Two data versions\n\nunadjusted and\ncomposition-adjusted\n(stage-blind; epithelial,\nstromal, immune)",
             "Gene selection\n\n2,000 most\nvariable genes\n(lowest-expressed\n25% removed)",
             "Embedding\n\nz-score, PCA:\n10 principal\ncomponents",
             "Subsampling\n\n20 subsamples\nper stage of\n80% of the\nsmallest stage",
             "Persistent homology\n\nVietoris–Rips:\nH0 clusters, H1 loops\ndiagrams, barcodes,\nBetti curves, landscapes\nraw and shape-only",
             "Statistical test\n\n1-Wasserstein distance\nlabel-permutation test\n(999 permutations)\nHolm correction"]
    n = len(steps)
    bw, bg, by, bh = 0.116, 0.022, 0.455, 0.175
    for i, s in enumerate(steps):
        bx = 0.03 + i * (bw + bg)
        title, rest = s.split("\n\n", 1)
        ax.add_patch(FancyBboxPatch((bx, by), bw, bh, boxstyle="round,pad=0.006,rounding_size=0.012",
                                    fc="#eef3fa" if i >= 5 else FILL, ec=EDGE, lw=0.9))
        ax.text(bx + bw / 2, by + bh - 0.025, title, ha="center", va="center", fontsize=8, fontweight="bold")
        ax.text(bx + bw / 2, by + (bh - 0.045) / 2 + 0.003, rest, ha="center", va="center", fontsize=7, linespacing=1.3)
        if i < n - 1:
            arrow(ax, bx + bw + 0.002, by + bh / 2, bx + bw + bg - 0.002, by + bh / 2)

    # ---------------- panel c: decision
    ax.text(0.015, 0.365, "c", fontsize=14, fontweight="bold")
    ax.text(0.045, 0.365, "Pre-registered decision and follow-up", fontsize=10, va="baseline")
    box(ax, 0.03, 0.17, 0.2, 0.15, "Pre-registered tipping-point rule\n\nshape-only Holm p < 0.05\nin BOTH unadjusted and\n"
        "composition-adjusted data", fc="#eef3fa", fs=7.8)
    arrow(ax, 0.232, 0.245, 0.262, 0.245)
    box(ax, 0.264, 0.19, 0.12, 0.11, "Not met in\nany cohort\n(Table 2)", fc="#fdecea", ec="#c0392b", fs=8.5, bold=True,
        color="#922b21")
    arrow(ax, 0.386, 0.245, 0.416, 0.245)
    ax.add_patch(FancyBboxPatch((0.418, 0.06), 0.355, 0.28, boxstyle="round,pad=0.006,rounding_size=0.012", fc="white",
                                ec=EDGE, lw=0.9, ls="--"))
    ax.text(0.5955, 0.318, "Follow-up analyses: pre-registered addenda, after the primary result", ha="center",
            fontsize=8, fontweight="bold")
    subs = [("Between-sample heterogeneity", "dispersion by stage (PERMDISP)", "Figure 4"),
            ("CMS subtypes", "share of carcinoma dispersion", "Figure 5"),
            ("Variable genes", "117 genes more variable\nin carcinoma", "Figure 6"),
            ("Chromosome arms", "20q, 20p, 13q, 8q;\nTCGA copy number", "Figure 7")]
    for k, (t, d, f) in enumerate(subs):
        cx = 0.43 + (k % 2) * 0.174
        cy = 0.19 - (k // 2) * 0.118
        box(ax, cx, cy, 0.158, 0.098, f"{t}\n{d}\n({f})", fc=FILL, fs=7.2)
    box(ax, 0.795, 0.06, 0.19, 0.28, "Supporting checks\n\n• power simulations and\n  positive/negative controls\n"
        "  (Supplementary Table\n  S_power_controls)\n• sensitivity analyses\n  (PCs, genes, sample sets)\n"
        f"• deviation log\n  (D1–D{dmax})", fc="#fbfbfb", fs=7.3)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("total GEO", total, "| TCGA", n_tcga, "with CN", n_tcga_cn, "| deviations D1-D%d" % dmax)
    print(t1[["Platform", "Sample set", "Normal", "Adenoma", "Carcinoma", "Metastasis", "Total samples"]].to_string())
    print("saved", os.path.join(OUT, NAME + ".png"))


if __name__ == "__main__":
    main()
