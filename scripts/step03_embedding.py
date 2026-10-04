"""
Step 3 - Feature selection and PCA embedding (GSE41258).

Uses the validated reference function tda_crc_pilot.preprocess unchanged:
  drop lowest-expressed 25% genes (by mean) -> top-N most variable genes ->
  z-score per gene -> PCA (k components).
Gene selection and PCA are fitted on ALL samples of a set pooled, without
stage labels. Labels are used only afterwards, to describe the embedding.

Embeddings
  primary:      2,000 genes, 10 PCs
  sensitivity:  2,000 genes, 5 and 20 PCs; 1,000 and 5,000 genes, 10 PCs
  for sets main (primary), withlung, main_allarrays, main_inclqc; unadjusted and adjusted.

Outputs: results/step03/*, results/step03/embeddings/*.csv
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kruskal
from sklearn.decomposition import PCA

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402

SEED = 42
np.random.seed(SEED)
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "step03")
EMB = os.path.join(OUT, "embeddings")
os.makedirs(EMB, exist_ok=True)
EXPECTED = {"main": {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39},
            "withlung": {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 58},
            "main_allarrays": {"Normal": 54, "Adenoma": 51, "Carcinoma": 184, "Metastasis": 47},
            "main_inclqc": {"Normal": 53, "Adenoma": 48, "Carcinoma": 183, "Metastasis": 39}}
PRIMARY = (2000, 10)
SENS = [(2000, 5), (2000, 20), (1000, 10), (5000, 10)]
COMP = ["Epithelial", "Stromal", "Immune", "Liver_fraction"]


def load(set_name, adjusted):
    lab = pd.read_csv(os.path.join(PROC, f"GSE41258_{set_name}_samples_qc.csv"), index_col=0)
    if lab["stage"].value_counts().to_dict() != EXPECTED[set_name]:
        sys.exit(f"Stage counts changed for {set_name}: {lab['stage'].value_counts().to_dict()}")
    if adjusted:
        g = pd.read_csv(os.path.join(PROC, f"GSE41258_{set_name}_expr_genes_adjusted.csv"), index_col=0)
    else:
        base = "withlung" if set_name == "withlung" else "main"
        g = pd.read_csv(os.path.join(PROC, f"GSE41258_{base}_expr_genes.csv"), index_col=0)
    return g[lab.index], lab


def selected_genes(expr, n_top):
    """Same gene selection as T.preprocess (recorded for traceability)."""
    X = expr.dropna(axis=0, how="any")
    mean = X.mean(axis=1)
    X = X[mean > mean.quantile(0.25)]
    var = X.var(axis=1).sort_values(ascending=False)
    return list(var.index[:n_top]), int(X.shape[0])


rows_var, rows_meta, primary = [], [], {}
ORGAN_REMOVED = {}
for set_name in EXPECTED:
    for adjusted in (False, True):
        g, lab = load(set_name, adjusted)
        version = "adjusted" if adjusted else "unadjusted"
        organ_removed = []
        if adjusted:   # adjusted data only: remove host-organ signature genes before selection
            for f in [f"liver_signature_{set_name}.csv"] + \
                     (["lung_signature_withlung.csv"] if set_name == "withlung" else []):
                organ_removed += list(pd.read_csv(os.path.join(ROOT, "results", "step02", f))["gene"])
            organ_removed = sorted(set(organ_removed) & set(g.index))
            g = g.drop(index=organ_removed)
            ORGAN_REMOVED[set_name] = organ_removed
        configs = [PRIMARY] + (SENS if set_name == "main" else [])
        for n_top, n_pcs in configs:
            T.LOG2_APPLIED[0] = False
            emb, pca, n_used = T.preprocess(g, n_top, n_pcs)     # validated reference code
            if T.LOG2_APPLIED[0]:
                sys.exit("preprocess applied log2 to already-log data - stop")
            genes, n_after_filter = selected_genes(g, n_top)
            assert n_used == len(genes)
            # exact-SVD check: reference PCA (random_state=0, auto solver) vs full SVD
            Xs = g.loc[genes]
            Z = Xs.sub(Xs.mean(axis=1), axis=0).div(Xs.std(axis=1) + 1e-9, axis=0).T.values
            exact = PCA(n_components=n_pcs, svd_solver="full").fit(Z)
            d_var = float(np.max(np.abs(exact.explained_variance_ratio_ - pca.explained_variance_ratio_)))
            ex_emb = exact.transform(Z)
            min_abs_r = float(min(abs(np.corrcoef(ex_emb[:, k], emb.values[:, k])[0, 1])
                                  for k in range(n_pcs)))
            tag = f"GSE41258_{set_name}_{version}_g{n_top}_pc{n_pcs}"
            out = emb.copy()
            out.columns = [f"PC{k + 1}" for k in range(n_pcs)]
            out.insert(0, "stage", lab.loc[out.index, "stage"])
            out.to_csv(os.path.join(EMB, tag + ".csv"), index_label="sample", float_format="%.6f")
            pd.Series(genes, name="gene").to_csv(os.path.join(EMB, tag + "_genes.csv"), index=False)
            for k, (v, c) in enumerate(zip(pca.explained_variance_ratio_,
                                           np.cumsum(pca.explained_variance_ratio_))):
                rows_var.append((set_name, version, n_top, n_pcs, f"PC{k + 1}", v, c))
            rows_meta.append(dict(set=set_name, version=version, n_genes=n_top, n_pcs=n_pcs,
                                  n_samples=out.shape[0], genes_after_low_filter=n_after_filter,
                                  var_explained=float(pca.explained_variance_ratio_.sum()),
                                  exact_svd_max_var_diff=d_var, exact_svd_min_abs_r=min_abs_r,
                                  file=os.path.relpath(os.path.join(EMB, tag + ".csv"), ROOT)))
            if (n_top, n_pcs) == PRIMARY:
                primary[(set_name, version)] = (out, genes, pca)

var_df = pd.DataFrame(rows_var, columns=["set", "version", "n_genes", "n_pcs", "PC",
                                         "var_explained", "cumulative"])
var_df.to_csv(os.path.join(OUT, "variance_explained.csv"), index=False)
meta_df = pd.DataFrame(rows_meta)
meta_df.to_csv(os.path.join(OUT, "embedding_index.csv"), index=False)

# ---------------------------------------------------------------- describe (labels used here only)
def describe(set_name, version):
    out, genes, pca = primary[(set_name, version)]
    sc = pd.read_csv(os.path.join(ROOT, "results", "step02", f"composition_scores_{set_name}.csv"),
                     index_col=0).loc[out.index]
    pcs = [c for c in out.columns if c.startswith("PC")]
    corr = pd.DataFrame(index=pcs, columns=COMP, dtype=float)
    pval = corr.copy()
    for p in pcs:
        for c in COMP:
            r, pv = spearmanr(out[p], sc[c])
            corr.loc[p, c], pval.loc[p, c] = r, pv
    stage_p = pd.Series({p: kruskal(*[out.loc[out.stage == s, p] for s in T.STAGES]).pvalue
                         for p in pcs})
    stage_med = out.groupby("stage")[pcs].median().reindex(T.STAGES).T
    return corr, pval, stage_p, stage_med


desc = {k: describe(*k) for k in primary}
for (s, v), (corr, pval, stage_p, stage_med) in desc.items():
    t = corr.round(3)
    t["stage_KW_p"] = stage_p
    t = t.join(stage_med.round(2))
    t.to_csv(os.path.join(OUT, f"pc_vs_composition_{s}_{v}.csv"), index_label="PC")

# overlap of selected genes, unadjusted vs adjusted (main)
gu, ga = set(primary[("main", "unadjusted")][1]), set(primary[("main", "adjusted")][1])

# ---------------------------------------------------------------- figures
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
for row, version in enumerate(("unadjusted", "adjusted")):
    out, genes, pca = primary[("main", version)]
    evr = pca.explained_variance_ratio_ * 100
    ax = axes[row, 0]
    for s in T.STAGES:
        d = out[out.stage == s]
        ax.scatter(d.PC1, d.PC2, s=12, c=T.COLORS[s], label=f"{s} ({len(d)})", alpha=0.85)
    ax.set_xlabel(f"PC1 ({evr[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({evr[1]:.1f}%)")
    ax.set_title(f"main, {version}: PC1 vs PC2"); ax.legend(fontsize=8)
    ax = axes[row, 1]
    for s in T.STAGES:
        d = out[out.stage == s]
        ax.scatter(d.PC3, d.PC4, s=12, c=T.COLORS[s], alpha=0.85)
    ax.set_xlabel(f"PC3 ({evr[2]:.1f}%)"); ax.set_ylabel(f"PC4 ({evr[3]:.1f}%)")
    ax.set_title(f"main, {version}: PC3 vs PC4")
    ax = axes[row, 2]
    corr = desc[("main", version)][0]
    im = ax.imshow(corr.values.astype(float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(COMP))); ax.set_xticklabels([c.replace("_fraction", "") for c in COMP])
    ax.set_yticks(range(len(corr))); ax.set_yticklabels(
        [f"{p} ({e:.1f}%)" for p, e in zip(corr.index, evr)])
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            v = float(corr.iloc[i, j])
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(v) > 0.6 else "black")
    ax.set_title(f"Spearman r: PC vs composition ({version})")
    fig.colorbar(im, ax=ax, fraction=0.04)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "pca_main_overview.png"), dpi=150)

fig, ax = plt.subplots(figsize=(7, 4))
for version, st in (("unadjusted", "-"), ("adjusted", "--")):
    d = var_df[(var_df.set == "main") & (var_df.version == version) & (var_df.n_genes == 2000)
               & (var_df.n_pcs == 20)]
    ax.plot(range(1, 21), d.var_explained * 100, st, marker="o", ms=3, label=version)
ax.axvline(10.5, color="0.6", ls=":")
ax.set_xlabel("principal component"); ax.set_ylabel("% variance explained")
ax.set_title("Scree plot, main set, 2,000 genes (dotted line = 10 PCs used)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT, "scree_main.png"), dpi=150)

# ---------------------------------------------------------------- summary
vers = {p: metadata.version(p) for p in ("scikit-learn", "pandas", "numpy", "scipy", "gudhi")}
lines = ["STEP 3 - Feature selection and PCA embedding (GSE41258)", "=" * 50,
         "Inputs: data/processed/GSE41258_<set>_samples_qc.csv; unadjusted = step-1 expression,",
         "  adjusted = data/processed/GSE41258_<set>_expr_genes_adjusted.csv (step 2)",
         f"Seed: {SEED} (numpy). PCA inside the validated reference function uses random_state=0;",
         "  checked against an exact full-SVD PCA (see 'exact SVD check' below).",
         "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()), "",
         "METHOD",
         "  Reference function tda_crc_pilot.preprocess, unchanged: drop lowest-expressed 25% of genes",
         "  (by mean), keep top-N most variable, z-score per gene, PCA.",
         "  Gene selection and PCA were fitted on ALL samples of each set pooled, WITHOUT stage labels.",
         "  Stage labels were used only afterwards, to describe the embedding (tables and figures below).",
         "  Data were already log2 (no further transform applied).",
         "  Adjusted data only: host-organ signature genes removed BEFORE gene selection (PI decision",
         "  2026-10-03). Liver signature (step 2) = genes >= 8-fold higher in the 13 reference normal-liver",
         "  samples than in BOTH normal colon and primary carcinoma, organ mean >= 8 log2; metastases",
         "  were not used to define it. withlung also removes the lung signature (7 normal lung samples).",
         "  Genes removed: " + "; ".join(f"{k} {len(v)}" for k, v in ORGAN_REMOVED.items()),
         "  Unadjusted data: no genes removed.",
         f"  Primary: {PRIMARY[0]} genes, {PRIMARY[1]} PCs. Sensitivity (main set): "
         + "; ".join(f"{a} genes/{b} PCs" for a, b in SENS),
         "  Sets: main (primary; one array per patient per stage), withlung, main_allarrays, main_inclqc;",
         "  each unadjusted and composition-adjusted.", "",
         "SAMPLES PER STAGE"]
for s in EXPECTED:
    lines.append(f"  {s:<12} " + ", ".join(f"{k} {v}" for k, v in EXPECTED[s].items())
                 + f" (total {sum(EXPECTED[s].values())}) - verified")
lines += ["", "EMBEDDINGS (variance explained by all PCs; exact SVD check)"]
for _, r in meta_df.iterrows():
    lines.append(f"  {r['set']:<12} {r['version']:<10} {r['n_genes']:>4} genes {r['n_pcs']:>2} PCs: "
                 f"{r['var_explained'] * 100:5.1f}%  | max diff in var. ratio {r['exact_svd_max_var_diff']:.1e}, "
                 f"min |r| per PC {r['exact_svd_min_abs_r']:.6f}")
lines.append("  Genes passing the low-expression filter: " + "; ".join(
    f"{r['set']} {r['version']} {r['genes_after_low_filter']}" for _, r in
    meta_df[(meta_df.n_genes == 2000) & (meta_df.n_pcs == 10)].iterrows()))
lines.append(f"  Top-2,000 genes shared by unadjusted and adjusted (main): {len(gu & ga)}/2000")
for version in ("unadjusted", "adjusted"):
    out, genes, pca = primary[("main", version)]
    corr, pval, stage_p, stage_med = desc[("main", version)]
    evr = pca.explained_variance_ratio_ * 100
    lines += ["", f"MAIN SET, {version.upper()} (2,000 genes, 10 PCs)",
              "  Variance explained per PC (%), cumulative:"]
    lines += [f"    PC{k + 1:<2} {e:5.1f}  cum {c:5.1f}" for k, (e, c) in enumerate(zip(evr, np.cumsum(evr)))]
    lines += ["  Spearman correlation of each PC with composition scores (step 2):",
              "         " + "".join(f"{c.replace('_fraction', ''):>12}" for c in COMP)]
    for p in corr.index:
        lines.append(f"    {p:<5}" + "".join(
            f"{corr.loc[p, c]:>9.2f}{'***' if pval.loc[p, c] < 1e-3 else '** ' if pval.loc[p, c] < 0.01 else '*  ' if pval.loc[p, c] < 0.05 else '   '}"
            for c in COMP))
    lines.append("    (* p<0.05, ** p<0.01, *** p<0.001; strongest |r| per PC listed in the CSV)")
    lines += ["  Stage separation per PC (descriptive; Kruskal-Wallis p, stage medians):"]
    for p in corr.index:
        lines.append(f"    {p:<5} p = {stage_p[p]:.1e}   " + "  ".join(
            f"{s[:4]} {stage_med.loc[p, s]:6.2f}" for s in T.STAGES))
lines += ["", "Outputs:",
          "  results/step03/embeddings/GSE41258_<set>_<version>_g<genes>_pc<PCs>.csv   sample coordinates",
          "  results/step03/embeddings/..._genes.csv                                     genes used",
          "  results/step03/embedding_index.csv, variance_explained.csv",
          "  results/step03/pc_vs_composition_<set>_<version>.csv",
          "  results/step03/pca_main_overview.png, scree_main.png"]
with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
