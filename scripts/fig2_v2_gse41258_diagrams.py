"""
Figure 2, version 2 (GSE41258, adjusted data): standard persistence diagrams (H0 + H1 in one panel) of the
representative subsample per stage (subsample 1 of the raw-mode draws = the 'representative subsample' of the addendum
12 figures), with an inset of that subsample's samples on PC1 vs PC2; raw (row a) and shape-only (row b) for the SAME
samples; H0 death box plots (panel c, as v1). Reads results/source_data/ (ph_diagrams_GEO.csv,
ph_subsample_members_GEO.csv, validation_summary_primary.csv) and the frozen embedding
results/step03/embeddings/GSE41258_main_adjusted_g2000_pc10.csv (insets; row b shape-only diagrams recomputed for the
raw-mode members with tda_crc_pilot.persistence, display only; the same recomputation in raw mode is checked against
the stored row-a diagrams). No statistical tests are run.

Usage: python scripts/fig2_v2_gse41258_diagrams.py
"""
import os

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402
SRC = os.path.join(ROOT, "results", "source_data")
EMB = os.path.join(ROOT, "results", "step03", "embeddings", "GSE41258_main_adjusted_g2000_pc10.csv")
OUT = os.path.join(ROOT, "paper", "final", "main_figures")
COHORT, VERSION, REP = "GSE41258", "adjusted", 1
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COL = dict(zip(STAGES, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))   # stage colours (addendum 12)
DIMCOL = {"H0": "#1f77b4", "H1": "#ff7f0e"}                              # dimension colours (standard)
MODES = [("raw", "raw distances"), ("shape-only", "shape-only (÷ median pairwise distance)")]
TR = [("Normal->Adenoma", "N→A"), ("Adenoma->Carcinoma", "A→C"), ("Carcinoma->Metastasis", "C→M")]

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(os.path.join(SRC, "ph_diagrams_GEO.csv"))
    d = d[(d.cohort == COHORT) & (d.version == VERSION)]
    mem = pd.read_csv(os.path.join(SRC, "ph_subsample_members_GEO.csv"))
    mem = mem[(mem.cohort == COHORT) & (mem.version == VERSION) & (mem.subsample == REP)]
    val = pd.read_csv(os.path.join(SRC, "validation_summary_primary.csv"))
    val = val[val.cohort == COHORT].set_index("transition")
    emb = pd.read_csv(EMB, index_col=0)
    pad = 0.05 * (emb.PC1.max() - emb.PC1.min())
    xl = (emb.PC1.min() - pad, emb.PC1.max() + pad)
    pad = 0.05 * (emb.PC2.max() - emb.PC2.min())
    yl = (emb.PC2.min() - pad, emb.PC2.max() + pad)

    # Rows a and b use the SAME samples: the raw-mode representative subsample (subsample 1). Row a = stored raw
    # diagrams; row b = shape-only diagrams recomputed from the frozen embedding for these members (display only).
    members = {s: mem[(mem["mode"] == "raw") & (mem.stage == s)].samples.iloc[0].split(";") for s in STAGES}
    pcs = [c for c in emb.columns if c.startswith("PC")]
    reps = {"raw": d[(d["mode"] == "raw") & (d.subsample == REP)]}
    rows, ok = [], True
    for s in STAGES:
        X = emb.loc[members[s], pcs].values
        for normalize, store in ((False, None), (True, rows)):
            dg = T.persistence(X, normalize)
            if store is None:                     # check: raw recomputation reproduces the stored row-a diagrams
                ref = reps["raw"][reps["raw"].stage == s]
                for k, dim in enumerate(("H0", "H1")):
                    r = ref[ref.dim == dim][["birth", "death"]].values
                    ok &= r.shape == dg[k].shape and np.allclose(r, dg[k], rtol=0, atol=1e-12)
            else:
                for k, dim in enumerate(("H0", "H1")):
                    store += [dict(stage=s, dim=dim, birth=b, death=dd) for b, dd in dg[k]]
    if not ok:
        raise SystemExit("raw recomputation does not reproduce the stored diagrams")
    print("check: raw diagrams recomputed from the frozen embedding = stored source data (all 4 stages)")
    reps["shape-only"] = pd.DataFrame(rows)
    for mode in reps:
        print(mode, "H1 counts:", reps[mode][reps[mode].dim == "H1"].groupby("stage").size().reindex(STAGES).tolist())

    fig = plt.figure(figsize=(11, 9.8))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 0.9], hspace=0.62, wspace=0.28)
    for row, (mode, mlab) in enumerate(MODES):
        rep = reps[mode]
        top = rep.death.max() * 1.08                   # finite features below; infinite H0 drawn at 'inf_y'
        inf_y = top * 1.06
        lim = inf_y * 1.04
        col_p = "raw_p_holm_adjus" if mode == "raw" else "shape_p_holm_adjus"
        for j, s in enumerate(STAGES):
            ax = fig.add_subplot(gs[row, j])
            x = rep[rep.stage == s]
            for dim in ("H0", "H1"):
                z = x[x.dim == dim]
                ax.scatter(z.birth, z.death, s=12, alpha=0.8, color=DIMCOL[dim], edgecolor="none",
                           label=f"{dim} ({len(z)})")
            ax.scatter([0], [inf_y], s=14, color=DIMCOL["H0"], edgecolor="none", clip_on=False)   # the infinite H0 bar
            ax.plot([0, lim], [0, lim], color="0.4", lw=0.8, ls="--")
            ax.axhline(inf_y, color="0.4", lw=0.8, ls="-.")
            ax.text(lim * 0.98, inf_y, "∞", ha="right", va="bottom", fontsize=8, color="0.3")
            ax.set(xlim=(-0.02 * lim, lim), ylim=(0, lim), xlabel="birth", ylabel="death" if j == 0 else "")
            ax.set_aspect("equal")
            ax.set_title(s, color=COL[s], fontweight="bold")
            ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.9), fontsize=6.5, handletextpad=0.2, borderaxespad=0.2)
            # inset: this subsample's samples on PC1 vs PC2 (same limits in every panel)
            ids = members[s]
            ins = inset_axes(ax, width="38%", height="38%", loc="lower right", borderpad=0.6)
            ins.scatter(emb.loc[ids, "PC1"], emb.loc[ids, "PC2"], s=4, color=COL[s], edgecolor="none")
            ins.set(xlim=xl, ylim=yl, xticks=[], yticks=[])
            ins.set_xlabel("PC1", fontsize=5.5, labelpad=1)
            ins.set_ylabel("PC2", fontsize=5.5, labelpad=1)
            for sp in ins.spines.values():
                sp.set_visible(True); sp.set_color("0.6"); sp.set_linewidth(0.6)
            if j == 0:
                ax.text(-0.42, 1.18, "ab"[row], transform=ax.transAxes, fontsize=13, fontweight="bold")
                ax.text(-0.22, 1.18, f"Persistence diagrams (H0 + H1), adjusted data, {mlab}; representative subsample "
                        f"(n = {len(ids)}; same samples in a and b)", transform=ax.transAxes, fontsize=9)
        ptxt = "   ".join(f"{lab} Holm p = {val.loc[t, col_p]:.3f}" for t, lab in TR)
        pos = ax.get_position()
        fig.text(0.5, pos.y0 - 0.058, f"Step 5 patient-level permutation test (H0 + H1, 20 subsamples, adjusted, {mode}):   {ptxt}",
                 ha="center", fontsize=8, style="italic")

    for k, (mode, mlab) in enumerate(MODES):                 # panel c, as v1
        ax = fig.add_subplot(gs[2, 2 * k:2 * k + 2])
        h0 = d[(d["mode"] == mode) & (d.dim == "H0")]
        data = [h0[h0.stage == s].death.values for s in STAGES]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        flierprops=dict(marker="o", markersize=2, alpha=0.4, markeredgecolor="none"),
                        medianprops=dict(color="black", lw=1))
        for patch, fl, s in zip(bp["boxes"], bp["fliers"], STAGES):
            patch.set_facecolor(COL[s]); patch.set_alpha(0.75); fl.set_markerfacecolor(COL[s])
        ax.set_xticks(range(1, 5), [f"{s}\n(n = {len(v)})" for s, v in zip(STAGES, data)])
        ax.set(ylabel="H0 death (merge scale)", title=f"H0 death values, all 20 subsamples, adjusted data, {mode}")
        if k == 0:
            ax.text(-0.13, 1.1, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

    fig.suptitle("GSE41258: persistence diagrams by stage (adjusted data; descriptive)", fontsize=11, fontweight="bold", y=0.995)
    fig.text(0.5, 0.005, "Rows a-b: the same representative subsample per stage (descriptive). Panel c pools 20 subsamples per stage. "
             "p-values come from the step 5 patient-level permutation test, not from these points.",
             ha="center", fontsize=7.5, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"Fig2_v2_GSE41258_diagrams.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", os.path.join(OUT, "Fig2_v2_GSE41258_diagrams.png"))


if __name__ == "__main__":
    main()
