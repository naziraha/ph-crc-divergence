# Figure script for the Nature Computational Science submission (NCS numbering); symbols italic via _italic_symbols.py.
# NCS VERSION 2 (2 Nov 2026), after Prof Peng's comments: copy of fig1_ph_scale_shape_revised.py.
# Changes: less white space between the panel a complexes, larger points and larger fonts throughout.
# Same data and checks; no new test. Original scripts untouched.
# NEW REVISED FIGURE (29 Oct 2026): main Figure 1 of the Nature Computational Science version.
# Built from fig2_v2_gse41258_diagrams_revised.py (unchanged panels a-c) plus a new panel d that summarises the
# primary PH test in all four GEO cohorts (Holm-adjusted p-values already stored in
# results/source_data/validation_summary_primary.csv; no new statistical test is run). The original scripts and
# figures are not changed. Output: outputs/main_figures/Fig1_PH_scale_shape.png/.pdf
# REVISED COPY (26 Oct 2026): same data and plotting as the original script; only the wording
# inside the figure is cleaned for journal style (no overall title, no internal terms such as
# "step", "addendum" or "frozen", no footnotes that repeat the legend). Outputs go to *_revised folders
# with "_revised" in the file name, so the original figures are not overwritten. No CSV is written.
# Significance asterisks added after p-values (26 Oct 2026): * p < 0.05, ** p < 0.01, *** p < 0.001.
"""
Figure 2, version 2 (GSE41258, adjusted data), three panels; a drawing of existing data, no statistical test is run.
  (a) Vietoris-Rips complexes of the representative subsample per stage (subsample 1 of the raw-mode draws = the
      'representative subsample' of the addendum 12 figures) at two scales (14 and 20, raw distances): an edge joins two
      samples whose 10-dimensional Euclidean distance (PC1-PC10 of the frozen embedding) is <= the scale, a triangle is
      shaded when all three of its edges are <= the scale. Samples are drawn on PC1/PC2 only, centred on each stage, with
      the same axis limits in all 8 panels; the number of connected components (clusters) at that scale is printed.
  (b) standard persistence diagrams (H0 + H1) of the same subsamples, raw distances (stored source data), with the
      label-permutation Holm p-values of both modes underneath (added 2026-10-18: panel a; insets removed from b).
  (c) H0 death values over all 20 subsamples per stage, raw and shape-only.
Check: persistence recomputed with gudhi.RipsComplex (max_dimension 2; tda_crc_pilot.persistence) from the same 31
samples per stage reproduces the stored raw diagrams in results/source_data/ph_diagrams_GEO.csv exactly (all 4 stages).
Reads results/source_data/ (ph_diagrams_GEO.csv, ph_subsample_members_GEO.csv, validation_summary_primary.csv) and
results/step03/embeddings/GSE41258_main_adjusted_g2000_pc10.csv.

Usage: python scripts/fig1_ph_scale_shape.py
"""
import itertools
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _italic_symbols  # noqa: F401  italic statistical symbols, Nature style (Zira, 3 Nov 2026)  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402
SRC = os.path.join(ROOT, "results", "source_data")
EMB = os.path.join(ROOT, "results", "step03", "embeddings", "GSE41258_main_adjusted_g2000_pc10.csv")
OUT = os.path.join(ROOT, "outputs", "main_figures")
COHORT, VERSION, REP = "GSE41258", "adjusted", 1
SCALES = (14, 20)
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COL = dict(zip(STAGES, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))   # stage colours (addendum 12)
DIMCOL = {"H0": "#1f77b4", "H1": "#ff7f0e"}                              # dimension colours (standard)
MODES = ["raw", "shape-only"]
TR = [("Normal->Adenoma", "N→A"), ("Adenoma->Carcinoma", "A→C"), ("Carcinoma->Metastasis", "C→M")]

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False})



def stars(p):
    """Significance asterisks: * p < 0.05, ** p < 0.01, *** p < 0.001 (added 26 Oct 2026)."""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(os.path.join(SRC, "ph_diagrams_GEO.csv"))
    d = d[(d.cohort == COHORT) & (d.version == VERSION)]
    mem = pd.read_csv(os.path.join(SRC, "ph_subsample_members_GEO.csv"))
    mem = mem[(mem.cohort == COHORT) & (mem.version == VERSION) & (mem.subsample == REP)]
    val = pd.read_csv(os.path.join(SRC, "validation_summary_primary.csv"))
    val = val[val.cohort == COHORT].set_index("transition")
    emb = pd.read_csv(EMB, index_col=0)
    pcs = [c for c in emb.columns if c.startswith("PC")]

    # same 31 samples per stage in panels a and b: raw-mode representative subsample (subsample 1)
    members = {s: mem[(mem["mode"] == "raw") & (mem.stage == s)].samples.iloc[0].split(";") for s in STAGES}
    rep = d[(d["mode"] == "raw") & (d.subsample == REP)]
    ok = True
    for s in STAGES:                      # check: recomputed persistence = stored raw diagrams
        dg = T.persistence(emb.loc[members[s], pcs].values, False)
        ref = rep[rep.stage == s]
        for k, dim in enumerate(("H0", "H1")):
            r = ref[ref.dim == dim][["birth", "death"]].values
            ok &= r.shape == dg[k].shape and np.allclose(r, dg[k], rtol=0, atol=1e-12)
    if not ok:
        raise SystemExit("recomputed persistence does not reproduce the stored diagrams")
    print("check: persistence recomputed (gudhi.RipsComplex, max_dimension 2) = stored source data (all 4 stages)")

    pts = {}
    for s in STAGES:
        X = emb.loc[members[s], pcs].values
        X = X - X.mean(axis=0)                                            # centred on the stage (distances unchanged)
        pts[s] = (X, squareform(pdist(X)))
    # Common limits fitted to the data range of all stages (not a symmetric square), so the
    # clouds fill the panels; limits stay identical across stages so scale remains comparable (Zira, 3 Nov)
    allp = np.vstack([v[0][:, :2] for v in pts.values()])
    (x0, y0), (x1, y1) = allp.min(axis=0), allp.max(axis=0)
    px, py = (x1 - x0) * 0.04, (y1 - y0) * 0.04
    XL, YL = (x0 - px, x1 + px), (y0 - py, y1 + (y1 - y0) * 0.13)        # top headroom for the cluster label

    fig = plt.figure(figsize=(14, 21.6))
    sfa, sfb, sfc, sfd = fig.subfigures(4, 1, height_ratios=[2.25, 1.55, 1.4, 1.8], hspace=0.01)

    # ---------------- a: Vietoris-Rips complexes
    sfa.text(0.005, 0.985, "a", fontsize=20, fontweight="bold", va="top")
    sfa.suptitle("Vietoris–Rips complexes at two scales (adjusted data, raw distances; one subsample of "
                 f"{len(members[STAGES[0]])} per stage)", x=0.04, ha="left", fontsize=14.5, y=0.985)
    ax = sfa.subplots(2, 4, gridspec_kw=dict(left=0.06, right=0.995, top=0.915, bottom=0.045, wspace=0.025, hspace=0.035))
    counts = {}
    for i, e in enumerate(SCALES):
        for j, s in enumerate(STAGES):
            a = ax[i, j]
            X, dm = pts[s]
            P, n = X[:, :2], len(X)
            tri = [P[[p, q, r]] for p, q, r in itertools.combinations(range(n), 3) if max(dm[p, q], dm[p, r], dm[q, r]) <= e]
            if tri:
                a.add_collection(PolyCollection(tri, facecolor=COL[s], alpha=0.07, edgecolor="none"))
            seg = [P[[p, q]] for p, q in itertools.combinations(range(n), 2) if dm[p, q] <= e]
            if seg:
                a.add_collection(LineCollection(seg, colors="0.35", linewidths=0.5, alpha=0.6))
            a.scatter(P[:, 0], P[:, 1], s=34, color=COL[s], edgecolor="k", linewidth=0.45, zorder=3)
            k = connected_components(dm <= e, directed=False)[0]
            counts[(e, s)] = int(k)
            a.text(0.04, 0.96, f"{k} cluster{'s' if k > 1 else ''}", transform=a.transAxes, va="top", fontsize=13,
                   zorder=2)
            a.set_xlim(*XL); a.set_ylim(*YL)   # headroom so the cluster label clears the points (Zira, 2 Nov)
            a.set_xticks([]); a.set_yticks([])
            for sp in a.spines.values():
                sp.set_visible(True)
            if i == 0:
                a.set_title(s, color=COL[s], fontweight="bold", fontsize=16)
            if j == 0:
                a.set_ylabel(f"scale = {e}\n(PC2)", fontsize=13)
            if i == 1:
                a.set_xlabel("PC1", fontsize=12)
    (lambda *a, **k: None)(0.07, 0.012, "Edges join samples within the scale in the 10-dimensional PCA space; shaded triangles are 2-simplices.\n"
             "Samples are drawn on PC1 and PC2 only (centred on each stage), so some edges look longer than others. "
             "Illustrative; no test is performed on this panel.", fontsize=9, style="italic")

    # ---------------- b: persistence diagrams (stored), no insets
    sfb.text(0.005, 0.985, "b", fontsize=20, fontweight="bold", va="top")
    sfb.suptitle("Persistence diagrams (H0 + H1) of the same subsamples, raw distances", x=0.04, ha="left",
                 fontsize=14.5, y=0.985)
    ax = sfb.subplots(1, 4, gridspec_kw=dict(left=0.07, right=0.99, top=0.82, bottom=0.25, wspace=0.28))
    top = rep.death.max() * 1.08                   # finite features below; infinite H0 drawn at 'inf_y'
    inf_y = top * 1.06
    lim_b = inf_y * 1.04
    for j, s in enumerate(STAGES):
        a = ax[j]
        x = rep[rep.stage == s]
        h0, h1 = x[x.dim == "H0"], x[x.dim == "H1"]
        a.plot([0, lim_b], [0, lim_b], color="0.4", lw=0.8, ls="--")
        a.axhline(inf_y, color="0.4", lw=0.8, ls="-.")
        a.text(lim_b * 0.98, inf_y, "∞", ha="right", va="bottom", fontsize=8, color="0.3")
        a.scatter(h0.birth, h0.death, s=24, alpha=0.8, color=DIMCOL["H0"], edgecolor="none", label=f"H0 ({len(h0)} + ∞)")
        a.scatter([0], [inf_y], s=14, color=DIMCOL["H0"], edgecolor="none", clip_on=False)   # the infinite H0 bar
        a.scatter(h1.birth, h1.death, s=24, alpha=0.8, color=DIMCOL["H1"], edgecolor="none", label=f"H1 ({len(h1)})")
        a.set(xlim=(-0.02 * lim_b, lim_b), ylim=(0, lim_b), xlabel="birth", ylabel="death" if j == 0 else "")
        a.set_aspect("equal")
        a.set_title(s, color=COL[s], fontweight="bold", fontsize=15)
        a.legend(loc="lower right", fontsize=10.5, handletextpad=0.2, borderaxespad=0.3)
    lines = []
    for mode in MODES:
        col_p = "raw_p_holm_adjus" if mode == "raw" else "shape_p_holm_adjus"
        lines.append(f"{mode}:  " + "   ".join(f"{lab} Holm p = {val.loc[t, col_p]:.3f}{stars(val.loc[t, col_p])}" for t, lab in TR))
        print(lines[-1])
    sfb.text(0.5, 0.03, "Holm-adjusted p-values (adjusted data):\n"
             + "\n".join(lines), ha="center", fontsize=12, style="italic", linespacing=1.4)

    # ---------------- c: H0 death values (unchanged content)
    sfc.text(0.005, 0.985, "c", fontsize=20, fontweight="bold", va="top")
    sfc.suptitle("H0 death values (merge scales), all 20 subsamples per stage", x=0.04, ha="left", fontsize=14.5, y=0.985)
    ax = sfc.subplots(1, 2, gridspec_kw=dict(left=0.07, right=0.99, top=0.8, bottom=0.17, wspace=0.18))
    for k, mode in enumerate(MODES):
        a = ax[k]
        h0 = d[(d["mode"] == mode) & (d.dim == "H0")]
        data = [h0[h0.stage == s].death.values for s in STAGES]
        bp = a.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                       flierprops=dict(marker="o", markersize=2, alpha=0.4, markeredgecolor="none"),
                       medianprops=dict(color="black", lw=1))
        for patch, fl, s in zip(bp["boxes"], bp["fliers"], STAGES):
            patch.set_facecolor(COL[s]); patch.set_alpha(0.75); fl.set_markerfacecolor(COL[s])
        a.set_xticks(range(1, 5), [f"{s}\n(n = {len(v)})" for s, v in zip(STAGES, data)])
        a.set(ylabel="H0 death (merge scale)")
        a.set_title(f"H0 death values, all 20 subsamples, adjusted data, {mode}", fontsize=13)

    # ---------------- d (new, 29 Oct 2026): primary PH test in all four GEO cohorts, stored Holm p-values
    sfd.text(0.005, 0.985, "d", fontsize=20, fontweight="bold", va="top")
    sfd.suptitle("Primary PH test in four cohorts (Holm-adjusted p-values)", x=0.04, ha="left", fontsize=14.5, y=0.985)
    allv = pd.read_csv(os.path.join(SRC, "validation_summary_primary.csv"))
    cols = [("raw_p_holm_unadj", "Raw\nunadjusted"), ("raw_p_holm_adjus", "Raw\nadjusted"),
            ("shape_p_holm_unadj", "Shape-only\nunadjusted"), ("shape_p_holm_adjus", "Shape-only\nadjusted")]
    lab = {"Normal->Adenoma": "N→A", "Adenoma->Carcinoma": "A→C", "Carcinoma->Metastasis": "C→M"}
    P = allv[[c for c, _ in cols]].values.astype(float)
    a = sfd.subplots(1, 1, gridspec_kw=dict(left=0.24, right=0.80, top=0.86, bottom=0.04))
    a.imshow(np.clip(-np.log10(P), 0, 3), cmap="Blues", vmin=0, vmax=3.6, aspect="auto")
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            p = P[i, j]
            a.text(j, i, f"{p:.3f}{stars(p)}", ha="center", va="center", fontsize=12.5,
                   color="white" if p < 0.01 else "black", fontweight="bold" if p < 0.05 else "normal")
    a.set_xticks(range(len(cols)), [t for _, t in cols], fontsize=12)
    a.xaxis.tick_top()
    a.set_yticks(range(len(allv)), [f"{c}   {lab[t]}" for c, t in zip(allv.cohort, allv.transition)], fontsize=12)
    a.tick_params(length=0)
    for sp in a.spines.values():
        sp.set_visible(False)
    a.axvline(1.5, color="k", lw=1.2)
    for i in range(1, len(allv)):
        if allv.cohort.iloc[i] != allv.cohort.iloc[i - 1]:
            a.axhline(i - 0.5, color="white", lw=3)
    print("panel d: rule met in", int((allv.rule_met == "yes").sum()), "of", len(allv), "transitions")

    (lambda *a, **k: None)(0.5, 0.002, "Panels a–b: one representative subsample per stage (descriptive). Panel c pools 20 subsamples per "
             "stage. p-values come from the label-permutation test (one sample per patient per stage), not from these "
             "points.", ha="center", fontsize=8, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"Fig1_PH_scale_shape.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("cluster counts (scale: N/A/C/M):", {e: [counts[(e, s)] for s in STAGES] for e in SCALES})
    print("saved", os.path.join(OUT, "Fig1_PH_scale_shape.png"))


if __name__ == "__main__":
    main()
