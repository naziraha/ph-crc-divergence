"""
Step 3 for GSE117606 (addendum 6, commit 93613bd; deviations D10, D11).

Frozen: tda_crc_pilot.preprocess unchanged (drop lowest 25% by mean, top 2,000 variable, z-score,
10 PCs), fitted on all 195 main-set arrays pooled, without labels; versions unadjusted, adjusted
(stage-blind, primary) and stageaware (secondary).
Saved for later steps (not interpreted here): sensitivity embeddings 5 / 20 PCs, 1,000 / 5,000 genes
(addendum 1 section 4); without the 13 quality-flagged arrays (addendum 6 section 7.2); without
GSM3304812 (D11). The sensitivity sets refit preprocess on their own samples.
Descriptive: the 10 SSA and 3 set-aside carcinomas projected onto the unadjusted embedding (z-scored
with the main set's gene means / SDs); D10 = Spearman r of PCs 1-10 with array quality.
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

GSE = "GSE117606"
SEED = 42
np.random.seed(SEED)
PROC = os.path.join(ROOT, "data", "processed")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "results", "step03", GSE)
os.makedirs(EMB, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
EXPECTED = {"Normal": 65, "Adenoma": 59, "Carcinoma": 71}
COMP = ["Epithelial", "Stromal", "Immune"]
QA_COLS = ["mean_corr", "median", "IQR"]
D11_DROP = "GSM3304812"


def selected_genes(expr, n_top):
    X = expr.dropna(axis=0, how="any")
    mean = X.mean(axis=1)
    X = X[mean > mean.quantile(0.25)]
    return list(X.var(axis=1).sort_values(ascending=False).index[:n_top])


def pseudo_r2(E, groups):
    tot = ((E - E.mean(0)) ** 2).sum()
    return float(sum((groups == g).sum() * ((E[groups == g].mean(0) - E.mean(0)) ** 2).sum()
                     for g in np.unique(groups)) / tot)


lab = pd.read_csv(os.path.join(PROC, f"{GSE}_main_samples_qc.csv"), index_col=0)
if lab["stage"].value_counts().to_dict() != EXPECTED:
    raise SystemExit(f"Stage counts changed: {lab['stage'].value_counts().to_dict()}")
expr = {"unadjusted": pd.read_csv(os.path.join(PROC, f"{GSE}_main_expr_genes.csv"), index_col=0)[lab.index],
        "adjusted": pd.read_csv(os.path.join(PROC, f"{GSE}_main_expr_genes_adjusted.csv"), index_col=0)[lab.index],
        "stageaware": pd.read_csv(os.path.join(PROC, f"{GSE}_main_expr_genes_adjusted_stageaware.csv"),
                                  index_col=0)[lab.index]}
qa = pd.read_csv(os.path.join(ROOT, "results", "step01", GSE, "array_quality.csv"), index_col=0).loc[lab.index]
sc = pd.read_csv(os.path.join(ROOT, "results", "step02", GSE, "composition_scores_all.csv"), index_col=0)
other = pd.read_csv(os.path.join(PROC, f"{GSE}_other_expr_genes.csv"), index_col=0)
olab = pd.read_csv(os.path.join(PROC, f"{GSE}_other_samples.csv"), index_col=0)

flagged = list(qa.index[qa.flagged])
sets = {"main": list(lab.index), "main_noQAflag": [s for s in lab.index if s not in flagged],
        "main_noGSM3304812": [s for s in lab.index if s != D11_DROP]}
configs = {"main": [(2000, 10), (2000, 5), (2000, 20), (1000, 10), (5000, 10)],
           "main_noQAflag": [(2000, 10)], "main_noGSM3304812": [(2000, 10)]}
idx_rows, prim = [], {}
for set_name, ids in sets.items():
    for version, g in expr.items():
        gv = g[ids]   # adjusted versions of the sensitivity sets reuse the adjustment fitted on all 195 arrays
        for n_top, n_pcs in configs[set_name]:
            T.LOG2_APPLIED[0] = False
            emb, pca, _ = T.preprocess(gv, n_top, n_pcs)
            if T.LOG2_APPLIED[0]:
                raise SystemExit("unexpected log2")
            e = emb.copy()
            e.columns = [f"PC{k + 1}" for k in range(n_pcs)]
            e.insert(0, "stage", lab.loc[e.index, "stage"])
            tag = f"{GSE}_{set_name}_{version}_g{n_top}_pc{n_pcs}"
            e.to_csv(os.path.join(EMB, tag + ".csv"), index_label="sample", float_format="%.6f")
            pd.Series(selected_genes(gv, n_top), name="gene").to_csv(os.path.join(EMB, tag + "_genes.csv"), index=False)
            idx_rows.append(dict(set=set_name, version=version, n_genes=n_top, n_pcs=n_pcs, n_samples=len(ids),
                                 var_explained=float(pca.explained_variance_ratio_.sum())))
            if set_name == "main" and (n_top, n_pcs) == (2000, 10):
                prim[version] = (e, pca, selected_genes(gv, n_top), gv)
pd.DataFrame(idx_rows).to_csv(os.path.join(OUT, "embedding_index.csv"), index=False)

# ---- projection of SSA / set-aside (unadjusted only)
e_u, pca_u, genes_u, g_u = prim["unadjusted"]
Xs = g_u.loc[genes_u]
mu, sd = Xs.mean(axis=1), Xs.std(axis=1) + 1e-9
Zo = other.loc[genes_u].sub(mu, axis=0).div(sd, axis=0).T.values
P = pd.DataFrame(pca_u.transform(Zo), index=other.columns, columns=[f"PC{k + 1}" for k in range(10)])
P.insert(0, "group", olab.loc[P.index, "set"])
P.to_csv(os.path.join(OUT, "projected_SSA_setaside_unadjusted.csv"), index_label="sample", float_format="%.4f")

# ---- descriptions
L = [f"STEP 3 - Feature selection and PCA embedding ({GSE})", "=" * 55,
     "Pre-specification: docs/addendum6_GSE117606.md (93613bd); deviations D10, D11.",
     "Method: tda_crc_pilot.preprocess unchanged (drop lowest 25% by mean, top 2,000 variable genes, z-score,",
     "  10 PCs; PCA random_state = 0), fitted on all 195 main-set arrays pooled WITHOUT stage labels.",
     f"Seed: {SEED}. Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("scikit-learn", "numpy", "pandas", "scipy")),
     "Samples per stage (verified): " + ", ".join(f"{s} {EXPECTED[s]}" for s in STAGES3), ""]
d10_rows = []
for version in ("unadjusted", "adjusted", "stageaware"):
    e, pca, genes, gv = prim[version]
    pcs = [c for c in e.columns if c.startswith("PC")]
    evr = pca.explained_variance_ratio_ * 100
    st = e["stage"].values
    r2 = pseudo_r2(e[pcs].values, st)
    kw = {p: kruskal(*[e.loc[e.stage == s, p] for s in STAGES3]).pvalue for p in ("PC1", "PC2")}
    medp = e.groupby("stage")[["PC1", "PC2"]].median().reindex(STAGES3)
    comp = pd.DataFrame({c: [spearmanr(e[p], sc.loc[e.index, c])[0] for p in pcs] for c in COMP}, index=pcs)
    qcor = pd.DataFrame({c: [spearmanr(e[p], qa.loc[e.index, c])[0] for p in pcs] for c in QA_COLS}, index=pcs)
    for p in pcs:
        d10_rows.append(dict(version=version, PC=p, var_pct=evr[int(p[2:]) - 1], **{f"r_{c}": qcor.loc[p, c] for c in QA_COLS},
                             **{f"r_{c}": comp.loc[p, c] for c in COMP}))
    dom = [p for p in pcs[:3] if qcor.loc[p].abs().max() >= 0.5]
    strong = [f"{p} (|r| {qcor.loc[p].abs().max():.2f})" for p in pcs if qcor.loc[p].abs().max() >= 0.5]
    L += [f"=== {version} ===",
          f"Variance explained: 10 PCs {evr.sum():.1f}% | per PC: " + ", ".join(f"{x:.1f}" for x in evr),
          f"Stage pseudo-R2 in 10-PC space: {r2:.3f}",
          "PC1 / PC2 stage medians: " + "; ".join(f"{s} {medp.loc[s, 'PC1']:.1f} / {medp.loc[s, 'PC2']:.1f}" for s in STAGES3)
          + f"  (Kruskal-Wallis p: PC1 {kw['PC1']:.1e}, PC2 {kw['PC2']:.1e})",
          "D10 - Spearman r of PCs with array quality (mean inter-array r, median, IQR) and composition:",
          "        " + "".join(f"{c:>10}" for c in QA_COLS + COMP)]
    for p in pcs:
        L.append(f"  {p:<5}" + "".join(f"{qcor.loc[p, c]:>10.2f}" for c in QA_COLS) + "".join(f"{comp.loc[p, c]:>10.2f}" for c in COMP))
    L += [f"  PCs with |r| >= 0.5 to any quality measure: {', '.join(strong) or 'none'}; among PC1-PC3: {', '.join(dom) or 'none'}", ""]
pd.DataFrame(d10_rows).round(3).to_csv(os.path.join(OUT, "D10_pc_vs_quality_and_composition.csv"), index=False)
pm = P.groupby("group")[["PC1", "PC2"]].median()
L += ["Projected descriptive arrays (unadjusted embedding), median PC1 / PC2:",
      "  " + pm.round(1).to_string().replace("\n", "\n  "), "",
      "Sensitivity embeddings saved for later steps (not interpreted): 5/20 PCs, 1,000/5,000 genes (main);",
      f"  main_noQAflag (n = {len(sets['main_noQAflag'])}; 13 flagged arrays removed), main_noGSM3304812 (n = {len(sets['main_noGSM3304812'])}).",
      "  For the adjusted versions of these two sets, the composition adjustment fitted on all 195 arrays is used.", "",
      "Outputs: results/step03/embeddings/GSE117606_<set>_<version>_g<genes>_pc<PCs>.csv (+ _genes.csv);",
      "  results/step03/GSE117606/embedding_index.csv, D10_pc_vs_quality_and_composition.csv,",
      "  projected_SSA_setaside_unadjusted.csv, pca_overview.png"]

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
for ax, version in zip(axes, ("unadjusted", "adjusted", "stageaware")):
    e, pca, _, _ = prim[version]
    for s in STAGES3:
        d = e[e.stage == s]
        ax.scatter(d.PC1, d.PC2, s=14, c=T.COLORS[s], label=f"{s} ({len(d)})", alpha=0.85)
    fl = e.loc[[x for x in flagged if x in e.index]]
    ax.scatter(fl.PC1, fl.PC2, s=45, facecolors="none", edgecolors="k", lw=0.8, label="quality-flagged")
    if version == "unadjusted":
        for grp, mk in (("SSA (descriptive)", "^"), ("set aside (A-coded carcinoma)", "s")):
            q = P[P.group == grp]
            ax.scatter(q.PC1, q.PC2, marker=mk, s=40, c="#7b4fb3", label=f"{grp} (projected)")
    ev = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
    ax.set_title(f"{GSE} main: {version}")
    ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "pca_overview.png"), dpi=140)
open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
