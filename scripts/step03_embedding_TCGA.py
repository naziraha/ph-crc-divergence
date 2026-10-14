"""
Step 3 for TCGA-COAD (addendum 10, commit 296eaed; D18, D19).

Frozen: tda_crc_pilot.preprocess unchanged (drop lowest 25% by mean, top 2,000 variable, z-score, 10 PCs;
PCA random_state 0), fitted without labels; versions unadjusted and adjusted (stage-blind, primary).
Sets (each refits preprocess on its own samples):
  main (274 tumours, stage I-IV)            : 2,000/10 primary; sensitivity 5 / 20 PCs, 1,000 / 5,000 genes
  main_noLowColon (D19b)                    : without the 18 low-colon-identity tumours
  main_MSS (addendum 10 section 4.2)        : MSS / MSI-L / MSS-like tumours only
  withnormal (addendum 10 section 4.3)      : 41 solid-tissue normals + 274 tumours (adjustment fitted on this set)
Adjusted versions of main_noLowColon and main_MSS reuse the adjustment fitted on all 274 tumours (as GSE117606).
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

GSE, SEED = "TCGA_COAD", 42
np.random.seed(SEED)
PROC = os.path.join(ROOT, "data", "processed")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "results", "step03", "TCGA")
os.makedirs(EMB, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
STAGES = ["Stage I", "Stage II", "Stage III", "Stage IV"]
EXPECTED = {"Stage I": 45, "Stage II": 110, "Stage III": 80, "Stage IV": 39}
COMP = ["Epithelial", "Stromal", "Immune"]
COL = {"Normal": "#2a78d6", "Stage I": "#7fb800", "Stage II": "#e0a100", "Stage III": "#e05a00", "Stage IV": "#c0392b"}


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
labw = pd.read_csv(os.path.join(PROC, f"{GSE}_withnormal_samples_qc.csv"), index_col=0)
raw = pd.read_csv(os.path.join(PROC, f"{GSE}_withnormal_expr_genes.csv"), index_col=0)
adj_main = pd.read_csv(os.path.join(PROC, f"{GSE}_main_expr_genes_adjusted.csv"), index_col=0)
adj_wn = pd.read_csv(os.path.join(PROC, f"{GSE}_withnormal_expr_genes_adjusted.csv"), index_col=0)
sc = pd.read_csv(os.path.join(ROOT, "results", "step02", GSE, "composition_scores_withnormal.csv"), index_col=0)
idq = pd.read_csv(os.path.join(ROOT, "results", "step02", GSE, "tissue_identity_qc_withnormal.csv"), index_col=0)

lowcol = [s for s in lab.index if isinstance(idq.loc[s, "flag"], str) and "low colon" in idq.loc[s, "flag"]]
mss = [s for s in lab.index if lab.loc[s, "msi_class"] in ("MSS/MSI-L", "MSS-like")]
sets = {"main": (list(lab.index), adj_main), "main_noLowColon": ([s for s in lab.index if s not in lowcol], adj_main),
        "main_MSS": (mss, adj_main), "withnormal": (list(labw.index), adj_wn)}
configs = {"main": [(2000, 10), (2000, 5), (2000, 20), (1000, 10), (5000, 10)],
           "main_noLowColon": [(2000, 10)], "main_MSS": [(2000, 10)], "withnormal": [(2000, 10)]}
stage_all = pd.concat([labw["stage"]])
rows, prim = [], {}
for set_name, (ids, adj) in sets.items():
    for version, g in (("unadjusted", raw), ("adjusted", adj)):
        gv = g[ids]
        for n_top, n_pcs in configs[set_name]:
            T.LOG2_APPLIED[0] = False
            emb, pca, _ = T.preprocess(gv, n_top, n_pcs)
            if T.LOG2_APPLIED[0]:
                raise SystemExit("unexpected log2")
            e = emb.copy()
            e.columns = [f"PC{k + 1}" for k in range(n_pcs)]
            e.insert(0, "stage", stage_all.loc[e.index])
            tag = f"{GSE}_{set_name}_{version}_g{n_top}_pc{n_pcs}"
            e.to_csv(os.path.join(EMB, tag + ".csv"), index_label="sample", float_format="%.6f")
            genes = selected_genes(gv, n_top)
            pd.Series(genes, name="gene").to_csv(os.path.join(EMB, tag + "_genes.csv"), index=False)
            rows.append(dict(set=set_name, version=version, n_genes=n_top, n_pcs=n_pcs, n_samples=len(ids),
                             **{f"n_{s}": int((e.stage == s).sum()) for s in ["Normal"] + STAGES},
                             var_explained=float(pca.explained_variance_ratio_.sum())))
            if (n_top, n_pcs) == (2000, 10):
                prim[(set_name, version)] = (e, pca)
idx = pd.DataFrame(rows)
idx.to_csv(os.path.join(OUT, "embedding_index.csv"), index=False)

msiH = lab["msi_class"].str.startswith("MSI-H").astype(float)
L = [f"STEP 3 - Feature selection and PCA embedding (TCGA-COAD)", "=" * 55,
     "Pre-specification: docs/addendum10_TCGA_stage.md (296eaed); deviations D18, D19.",
     "Method: tda_crc_pilot.preprocess unchanged (drop lowest 25% by mean, top 2,000 variable genes, z-score, 10 PCs;",
     "  PCA random_state 0), fitted per set WITHOUT stage labels; versions unadjusted and adjusted (stage-blind).",
     f"Seed {SEED}. Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("scikit-learn", "numpy", "pandas", "scipy")),
     "Sample counts per set (verified against step 1/2):"]
for s, (ids, _) in sets.items():
    st = stage_all.loc[ids].value_counts()
    L.append(f"  {s:<16} n {len(ids)}: " + ", ".join(f"{k} {st.get(k, 0)}" for k in ["Normal"] + STAGES if st.get(k, 0)))
L.append(f"  low-colon-identity tumours removed in main_noLowColon: {len(lowcol)}; MSI-H/MSI-H-like/unknown removed in main_MSS: {len(lab) - len(mss)}")
L.append("")
pcrows = []
for (set_name, version), (e, pca) in prim.items():
    if set_name not in ("main", "withnormal"):
        L.append(f"{set_name} {version}: 10 PCs {100 * pca.explained_variance_ratio_.sum():.1f}% variance")
        continue
    pcs = [c for c in e.columns if c.startswith("PC")]
    evr = pca.explained_variance_ratio_ * 100
    st = e.stage.values
    groups = [s for s in ["Normal"] + STAGES if (st == s).any()]
    kw = {p: kruskal(*[e.loc[e.stage == s, p] for s in groups]).pvalue for p in ("PC1", "PC2")}
    medp = e.groupby("stage")[["PC1", "PC2"]].median().reindex(groups)
    comp = pd.DataFrame({c: [spearmanr(e[p], sc.loc[e.index, c])[0] for p in pcs] for c in COMP}, index=pcs)
    L += [f"=== {set_name} / {version} ===",
          f"Variance explained: 10 PCs {evr.sum():.1f}% | per PC: " + ", ".join(f"{x:.1f}" for x in evr),
          f"Stage pseudo-R2 in 10-PC space: {pseudo_r2(e[pcs].values, st):.3f}",
          "PC1 / PC2 stage medians: " + "; ".join(f"{s} {medp.loc[s, 'PC1']:.1f} / {medp.loc[s, 'PC2']:.1f}" for s in groups)
          + f"  (Kruskal-Wallis p: PC1 {kw['PC1']:.1e}, PC2 {kw['PC2']:.1e})",
          "Spearman r of PCs with composition" + (" and MSI-H (tumours):" if set_name == "main" else ":"),
          "        " + "".join(f"{c:>11}" for c in COMP) + ("      MSI-H" if set_name == "main" else "")]
    for p in pcs:
        rm = spearmanr(e[p], msiH.loc[e.index])[0] if set_name == "main" else np.nan
        L.append(f"  {p:<5}" + "".join(f"{comp.loc[p, c]:>11.2f}" for c in COMP) + (f"{rm:>11.2f}" if set_name == "main" else ""))
        pcrows.append(dict(set=set_name, version=version, PC=p, var_pct=evr[int(p[2:]) - 1],
                           **{f"r_{c}": comp.loc[p, c] for c in COMP}, r_MSIH=rm))
    L.append("")
pd.DataFrame(pcrows).round(3).to_csv(os.path.join(OUT, "pc_vs_composition_msi.csv"), index=False)
L += ["Sensitivity embeddings saved for later steps (not interpreted): main 5/20 PCs, 1,000/5,000 genes;",
      "  main_noLowColon, main_MSS (adjusted versions reuse the adjustment fitted on all 274 tumours); withnormal.",
      "Outputs: results/step03/embeddings/TCGA_COAD_<set>_<version>_g<genes>_pc<PCs>.csv (+ _genes.csv);",
      "  results/step03/TCGA/embedding_index.csv, pc_vs_composition_msi.csv, pca_overview.png"]

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
fig, axes = plt.subplots(1, 4, figsize=(20, 4.8))
for ax, key in zip(axes, [("main", "unadjusted"), ("main", "adjusted"), ("withnormal", "unadjusted"), ("withnormal", "adjusted")]):
    e, pca = prim[key]
    for s in ["Normal"] + STAGES:
        d = e[e.stage == s]
        if len(d):
            ax.scatter(d.PC1, d.PC2, s=12, c=COL[s], label=f"{s} ({len(d)})", alpha=0.8)
    if key[0] == "main":
        m = e.loc[[x for x in e.index if msiH.get(x, 0) == 1]]
        ax.scatter(m.PC1, m.PC2, s=40, facecolors="none", edgecolors="k", lw=0.7, label="MSI-H")
    ev = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
    ax.set_title(f"TCGA-COAD {key[0]}: {key[1]}"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "pca_overview.png"), dpi=140)
open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
