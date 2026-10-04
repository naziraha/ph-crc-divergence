"""
Step 3 (validation cohorts) - feature selection and PCA embedding.

Frozen settings (tda_crc_pilot.preprocess, unchanged): drop lowest-expressed 25% of
genes, top-N most variable, z-score, PCA. Fitted on all analysis samples of a set
pooled, without stage labels. Versions: unadjusted, adjusted (primary, stage-blind),
stageaware (secondary, addendum 2).

GSE20916 micro pooling rule (addendum 2, evaluated once on unadjusted micro_pooled,
2,000 genes / 10 PCs): pseudo-R2 = between-group / total sum of squared Euclidean
distances (one-way PERMANOVA form) for stage and for preparation; pool if
R2(stage) > R2(preparation), else describe micro_crypt and micro_mucosa separately.
Adjacent normals are projected onto the unadjusted embeddings (description only).

Usage: python scripts/step03_embedding_validation.py GSE20916
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kruskal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402

SEED = 42
np.random.seed(SEED)
PROC = os.path.join(ROOT, "data", "processed")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
COMP = ["Epithelial", "Stromal", "Immune"]
PRIMARY, SENS = (2000, 10), [(2000, 5), (2000, 20), (1000, 10), (5000, 10)]
COHORT = {"GSE20916": {"sets": {"bulk": ("bulk", {"Normal": 24, "Adenoma": 45, "Carcinoma": 36}),
                                "micro_pooled": ("micro", {"Normal": 10, "Adenoma": 10, "Carcinoma": 10}),
                                "micro_crypt": ("micro", {"Normal": 5, "Adenoma": 5, "Carcinoma": 5}),
                                "micro_mucosa": ("micro", {"Normal": 5, "Adenoma": 5, "Carcinoma": 5})}},
          "GSE37364": {"sets": {"bulk": ("bulk", {"Normal": 38, "Adenoma": 29, "Carcinoma": 27})}}}


def selected_genes(expr, n_top):
    X = expr.dropna(axis=0, how="any")
    mean = X.mean(axis=1)
    X = X[mean > mean.quantile(0.25)]
    return list(X.var(axis=1).sort_values(ascending=False).index[:n_top])


def pseudo_r2(E, groups):
    tot = ((E - E.mean(0)) ** 2).sum()
    btw = sum((groups == g).sum() * ((E[groups == g].mean(0) - E.mean(0)) ** 2).sum()
              for g in np.unique(groups))
    return float(btw / tot)


def main():
    gse = sys.argv[1]
    out = os.path.join(ROOT, "results", "step03", gse)
    os.makedirs(out, exist_ok=True)
    os.makedirs(EMB, exist_ok=True)
    cfg = COHORT[gse]["sets"]
    sc_all = {raw: pd.read_csv(os.path.join(ROOT, "results", "step02", gse, f"composition_scores_{raw}.csv"),
                               index_col=0) for raw in {v[0] for v in cfg.values()}}
    raw_expr = {raw: pd.read_csv(os.path.join(PROC, f"{gse}_{raw}_expr_genes.csv"), index_col=0)
                for raw in sc_all}
    rows, prim = [], {}
    for aset, (raw, expected) in cfg.items():
        lab = pd.read_csv(os.path.join(PROC, f"{gse}_{aset}_samples_qc.csv"), index_col=0)
        cnt = lab["stage"].value_counts().to_dict()
        if cnt != expected:
            raise SystemExit(f"Stage counts changed for {aset}: {cnt}")
        for version in ("unadjusted", "adjusted", "stageaware"):
            if version == "unadjusted":
                g = raw_expr[raw][lab.index]
            else:
                suffix = "adjusted" if version == "adjusted" else "adjusted_stageaware"
                g = pd.read_csv(os.path.join(PROC, f"{gse}_{aset}_expr_genes_{suffix}.csv"), index_col=0)[lab.index]
            for n_top, n_pcs in [PRIMARY] + SENS:
                n_pcs_eff = min(n_pcs, len(lab) - 1)
                T.LOG2_APPLIED[0] = False
                emb, pca, n_used = T.preprocess(g, n_top, n_pcs_eff)
                if T.LOG2_APPLIED[0]:
                    raise SystemExit("unexpected log2 on log data")
                genes = selected_genes(g, n_top)
                tag = f"{gse}_{aset}_{version}_g{n_top}_pc{n_pcs}"
                e = emb.copy()
                e.columns = [f"PC{k + 1}" for k in range(n_pcs_eff)]
                e.insert(0, "stage", lab.loc[e.index, "stage"])
                e.to_csv(os.path.join(EMB, tag + ".csv"), index_label="sample", float_format="%.6f")
                pd.Series(genes, name="gene").to_csv(os.path.join(EMB, tag + "_genes.csv"), index=False)
                rows.append(dict(set=aset, version=version, n_genes=n_top, n_pcs=n_pcs, n_pcs_used=n_pcs_eff,
                                 n_samples=len(lab), var_explained=float(pca.explained_variance_ratio_.sum())))
                if (n_top, n_pcs) == PRIMARY:
                    prim[(aset, version)] = (e, genes, pca, g)
    idx = pd.DataFrame(rows)
    idx.to_csv(os.path.join(out, "embedding_index.csv"), index=False)

    # ---------------- pooling decision (micro)
    decision = None
    if "micro_pooled" in cfg:
        e, genes, pca, g = prim[("micro_pooled", "unadjusted")]
        lab = pd.read_csv(os.path.join(PROC, f"{gse}_micro_pooled_samples_qc.csv"), index_col=0)
        E = e[[c for c in e.columns if c.startswith("PC")]].values
        r2_stage = pseudo_r2(E, lab.loc[e.index, "stage"].values)
        r2_prep = pseudo_r2(E, lab.loc[e.index, "morphology"].values)
        decision = dict(r2_stage=r2_stage, r2_prep=r2_prep,
                        pool=r2_stage > r2_prep,
                        r2_stage_within_crypt=None, r2_stage_within_mucosa=None)
        for prep, aset in (("colonic crypt epithelial cells", "micro_crypt"), ("mucosa", "micro_mucosa")):
            ee = prim[(aset, "unadjusted")][0]
            decision[f"r2_stage_within_{aset.split('_')[1]}"] = pseudo_r2(
                ee[[c for c in ee.columns if c.startswith("PC")]].values, ee["stage"].values)
        pd.Series(decision).to_csv(os.path.join(out, "micro_pooling_decision.csv"), header=["value"])

    # ---------------- description: PC vs composition, stage separation
    desc_lines = []
    for (aset, version), (e, genes, pca, g) in prim.items():
        raw = cfg[aset][0]
        sc = sc_all[raw].loc[e.index]
        pcs = [c for c in e.columns if c.startswith("PC")]
        corr = pd.DataFrame({c: [spearmanr(e[p], sc[c])[0] for p in pcs] for c in COMP}, index=pcs)
        stage_p = [kruskal(*[e.loc[e.stage == s, p] for s in STAGES3]).pvalue for p in pcs]
        tab = corr.round(2)
        tab["stage_KW_p"] = stage_p
        tab = tab.join(e.groupby("stage")[pcs].median().reindex(STAGES3).T.round(2))
        tab["var_pct"] = (pca.explained_variance_ratio_ * 100).round(1)
        tab.to_csv(os.path.join(out, f"pc_vs_composition_{aset}_{version}.csv"), index_label="PC")
        lab = pd.read_csv(os.path.join(PROC, f"{gse}_{aset}_samples_qc.csv"), index_col=0)
        r2s = pseudo_r2(e[pcs].values, e["stage"].values)
        desc_lines.append((aset, version, pca.explained_variance_ratio_, corr, stage_p, tab, r2s))

    # ---------------- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = [s for s in cfg]
    fig, axes = plt.subplots(len(show), 3, figsize=(15, 4.3 * len(show)), squeeze=False)
    for r, aset in enumerate(show):
        raw = cfg[aset][0]
        for c, version in enumerate(("unadjusted", "adjusted", "stageaware")):
            ax = axes[r, c]
            e, genes, pca, g = prim[(aset, version)]
            lab = pd.read_csv(os.path.join(PROC, f"{gse}_{aset}_samples_qc.csv"), index_col=0)
            for s in STAGES3:
                d = e[e.stage == s]
                mk = None
                if raw == "micro" and aset == "micro_pooled":
                    for prep, m in (("colonic crypt epithelial cells", "o"), ("mucosa", "^")):
                        dd = d[lab.loc[d.index, "morphology"] == prep]
                        ax.scatter(dd.PC1, dd.PC2, c=T.COLORS[s], marker=m, s=30,
                                   label=f"{s} {'crypt' if m == 'o' else 'mucosa'}")
                else:
                    ax.scatter(d.PC1, d.PC2, c=T.COLORS[s], s=18 if raw == "bulk" else 34, label=f"{s} ({len(d)})")
            if version == "unadjusted" and raw == "micro":
                # project adjacent normals of the same preparation(s)
                full = pd.read_csv(os.path.join(PROC, f"{gse}_micro_samples.csv"), index_col=0)
                adjn = full.index[(full.stage == "Normal") & (full.normal_type != "distant")]
                if aset == "micro_crypt":
                    adjn = [s for s in adjn if full.loc[s, "morphology"] == "colonic crypt epithelial cells"]
                elif aset == "micro_mucosa":
                    adjn = [s for s in adjn if full.loc[s, "morphology"] == "mucosa"]
                Xs = g.loc[genes]
                mu, sd = Xs.mean(axis=1), Xs.std(axis=1) + 1e-9
                Z = raw_expr["micro"].loc[genes, adjn].sub(mu, axis=0).div(sd, axis=0).T.values
                P = pca.transform(Z)
                ax.scatter(P[:, 0], P[:, 1], facecolors="none", edgecolors="0.4", s=34,
                           label="adjacent normal (projected)")
            ev = pca.explained_variance_ratio_ * 100
            ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
            ax.set_title(f"{gse} {aset}: {version}")
            if c == 0:
                ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "pca_overview.png"), dpi=140)

    # ---------------- summary
    vers = {p: metadata.version(p) for p in ("scikit-learn", "pandas", "numpy", "scipy")}
    L = [f"STEP 3 - Feature selection and PCA embedding ({gse})", "=" * 55,
         "Pre-specification: addendum 1 (6ea31ea), addendum 2 (d73e5d1)",
         f"Inputs: data/processed/{gse}_<set>_samples_qc.csv; unadjusted = step-1 expression; "
         "adjusted / stageaware = step-2 outputs",
         f"Seed: {SEED}; reference preprocess (PCA random_state=0)",
         "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()), "",
         "METHOD: tda_crc_pilot.preprocess unchanged (drop lowest 25% by mean, top-N variable, z-score, PCA),",
         "  fitted on all analysis samples of each set pooled, WITHOUT stage labels; labels used only to describe.",
         f"  Primary {PRIMARY[0]} genes / {PRIMARY[1]} PCs; sensitivity "
         + "; ".join(f"{a}/{b}" for a, b in SENS)
         + " (PCs capped at n-1 for 15-sample sets).",
         "  Versions: unadjusted; adjusted (primary, stage not in model); stageaware (secondary).", "",
         "SAMPLES PER STAGE (verified): " + "; ".join(
             f"{k} " + "/".join(str(v[1][s]) for s in STAGES3) for k, v in cfg.items()), ""]
    if decision:
        L += ["MICRO POOLING DECISION (pre-specified, unadjusted micro_pooled, 2,000 genes / 10 PCs):",
              f"  pseudo-R2 stage       = {decision['r2_stage']:.3f}",
              f"  pseudo-R2 preparation = {decision['r2_prep']:.3f}",
              "  -> " + ("POOL the two preparations (stage explains more)." if decision["pool"] else
                         "DO NOT POOL: preparation explains more than stage; micro_crypt and micro_mucosa "
                         "are described separately."),
              f"  (descriptive) stage pseudo-R2 within crypt epithelium {decision['r2_stage_within_crypt']:.3f}, "
              f"within mucosa {decision['r2_stage_within_mucosa']:.3f}", ""]
    L.append("VARIANCE EXPLAINED (10 PCs, 2,000 genes) and stage pseudo-R2 in PC space:")
    for aset, version, evr, corr, sp, tab, r2s in desc_lines:
        L.append(f"  {aset:<13} {version:<11} {evr.sum() * 100:5.1f}%  PC1 {evr[0] * 100:4.1f}%  "
                 f"PC2 {evr[1] * 100:4.1f}%   stage R2 {r2s:.3f}")
    for aset, version, evr, corr, sp, tab, r2s in desc_lines:
        if aset == "bulk" or (aset == "micro_crypt" and version == "unadjusted"):
            L += ["", f"{aset.upper()}, {version}: PC vs composition (Spearman r), stage KW p, stage medians",
                  "  " + tab.to_string().replace("\n", "\n  ")]
    L += ["", "NOTE (PI): for the microdissected set, the unadjusted crypt-epithelium data are emphasised,",
          "  because laser capture is itself the physical purification. Adjusted micro versions are fitted",
          "  on 15 (per preparation) or 30 (pooled) samples with 3-5 regression terms and overfit",
          "  (step 2: median 48-79% of each gene's variance removed); read them with caution.",
          "  The microdissected set is descriptive only (addendum 1: < 13 samples per stage).", "",
          "Outputs:",
          f"  results/step03/embeddings/{gse}_<set>_<version>_g<genes>_pc<PCs>.csv (+ _genes.csv)",
          f"  results/step03/{gse}/embedding_index.csv, pc_vs_composition_<set>_<version>.csv,",
          f"  results/step03/{gse}/micro_pooling_decision.csv, pca_overview.png"]
    with open(os.path.join(out, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
