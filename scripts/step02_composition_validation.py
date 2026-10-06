"""
Step 2 (validation cohorts) - tissue identity screen and composition scoring.

Pre-specification: addendum 1 (commit 6ea31ea) and addendum 2 (commit d73e5d1).
  * Identity screen within each raw set: fixed colon / hepatocyte / lung / Mullerian
    panels; flag if colon panel < median - 3 robust SD, or another panel > median + 3 robust SD.
  * Composition: marker genes FROZEN from GSE41258 (results/step02/marker_sets_used.csv:
    18 pan-epithelial genes, PanglaoDB stromal and immune sets), decoupler ULM on
    gene-centred log2 data (centred on the analysis samples of the set).
  * Adjusted data, per gene OLS:
      primary     : composition covariates only (stage NOT in the model)
      stage-aware : stage indicators + composition covariates (secondary; addendum 2)
    In both, only the composition part is removed: adjusted = y - sum b_c (c - mean c).
  * GSE20916 micro: analysis samples = distant normal + adenoma + carcinoma; adjacent
    normals scored for description only. Sets: micro_pooled, micro_crypt, micro_mucosa.

Usage: python scripts/step02_composition_validation.py GSE20916
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402
import decoupler as dc  # noqa: E402

SEED = 42
PROC = os.path.join(ROOT, "data", "processed")
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
COMP = ["Epithelial", "Stromal", "Immune"]
PANELS = {"colon": ["CDX1", "CDX2", "CDH17", "VIL1", "KRT20", "CEACAM5", "LGALS4"],
          "hepatocyte": ["ALB", "APOA1", "APOA2", "APOH", "FGA", "FGB", "FGG", "AHSG", "TF",
                         "HP", "SERPINA1", "TTR"],
          "lung": ["SFTPC", "SFTPB", "SFTPD", "SCGB1A1", "AGER"],
          "mullerian": ["PAX8", "SOX17", "MUC16"]}
SANITY = {"Epithelial": ["EPCAM", "CDH1", "KRT8"], "Stromal": ["COL1A2", "COL3A1", "DCN", "VWF"],
          "Immune": ["PTPRC", "CD3E", "CD74", "LCP1"]}
NET = pd.read_csv(os.path.join(ROOT, "results", "step02", "marker_sets_used.csv"))   # frozen


def raw_sets(gse):
    """raw set name -> (expr, labels, analysis sample ids per analysis set, descriptive-only ids)"""
    out = {}
    if gse == "GSE20916":
        for raw in ("bulk", "micro"):
            g = pd.read_csv(os.path.join(PROC, f"{gse}_{raw}_expr_genes.csv"), index_col=0)
            lab = pd.read_csv(os.path.join(PROC, f"{gse}_{raw}_samples.csv"), index_col=0)
            assert list(g.columns) == list(lab.index)
            if raw == "bulk":
                out[raw] = (g, lab, {"bulk": list(lab.index)}, [])
            else:
                keep = lab.index[(lab.stage != "Normal") | (lab.normal_type == "distant")]
                adj_n = list(lab.index[(lab.stage == "Normal") & (lab.normal_type != "distant")])
                crypt = [s for s in keep if lab.loc[s, "morphology"] == "colonic crypt epithelial cells"]
                muc = [s for s in keep if lab.loc[s, "morphology"] == "mucosa"]
                out[raw] = (g, lab, {"micro_pooled": list(keep), "micro_crypt": crypt,
                                     "micro_mucosa": muc}, adj_n)
    elif gse == "GSE117606":
        # addendum 6: main set (65/59/71) for screen, scores and adjustment; SSA (10) and set-aside
        # A-coded carcinomas (3) are screened and scored for description only (not used in any fit)
        gm = pd.read_csv(os.path.join(PROC, f"{gse}_main_expr_genes.csv"), index_col=0)
        go = pd.read_csv(os.path.join(PROC, f"{gse}_other_expr_genes.csv"), index_col=0)
        lm = pd.read_csv(os.path.join(PROC, f"{gse}_main_samples.csv"), index_col=0)
        lo = pd.read_csv(os.path.join(PROC, f"{gse}_other_samples.csv"), index_col=0)
        assert list(gm.columns) == list(lm.index) and list(go.columns) == list(lo.index)
        g = pd.concat([gm, go], axis=1)
        lab = pd.concat([lm, lo])
        lab["raw_label"] = lab["tissue"] + "; " + lab["title"]
        lab["desc_group"] = np.where(lab["set"] == "main", lab["stage"], lab["set"])
        out["all"] = (g, lab, {"main": list(lm.index)}, list(lo.index))
    elif gse == "GSE37364":
        g = pd.read_csv(os.path.join(PROC, f"{gse}_bulk_expr_genes.csv"), index_col=0)
        lab = pd.read_csv(os.path.join(PROC, f"{gse}_bulk_samples.csv"), index_col=0)
        assert list(g.columns) == list(lab.index)
        out["bulk"] = (g, lab, {"bulk": list(lab.index)}, [])
    else:
        raise SystemExit(f"cohort {gse} not configured")
    return out


def repeat_array_check(g, lab):
    """No patient ids: look for near-identical arrays. Pearson r on log2 genes above the
    25th percentile of mean expression; each sample's nearest neighbour; flag if r > 0.99
    and > median + 4 robust SD of the nearest-neighbour distribution. Shown to the PI."""
    m = g.mean(axis=1)
    x = g[m > m.quantile(0.25)]
    r = np.corrcoef(x.T.values)
    np.fill_diagonal(r, np.nan)
    r = pd.DataFrame(r, index=g.columns, columns=g.columns)
    nn = r.max(axis=1)
    nn_id = r.idxmax(axis=1)
    thr = nn.median() + 4 * 1.4826 * (nn - nn.median()).abs().median()
    tab = pd.DataFrame({"nearest": nn_id, "r": nn.round(4), "stage": lab["stage"],
                        "nearest_stage": lab.loc[nn_id, "stage"].values})
    tab["flag"] = (tab.r > 0.99) & (tab.r > thr)
    return tab, float(thr)


def identity_screen(g):
    pm = pd.DataFrame({p: g.loc[[x for x in gl if x in g.index]].mean() for p, gl in PANELS.items()})
    med = pm.median()
    rsd = 1.4826 * (pm - med).abs().median()
    z = (pm - med) / rsd.replace(0, np.nan)
    flags = pd.Series("", index=pm.index)
    for s in pm.index:
        f = []
        if z.loc[s, "colon"] < -3:
            f.append("low colon identity")
        for p in ("hepatocyte", "lung", "mullerian"):
            if z.loc[s, p] > 3:
                f.append(f"high {p}")
        flags[s] = "; ".join(f)
    return pm, z, flags


def scores(g, centre_ids):
    centred = g.sub(g[centre_ids].mean(axis=1), axis=0)
    net = NET[NET["target"].isin(g.index)][["source", "target", "weight"]]
    res = dc.mt.ulm(centred.T, net, tmin=5)
    es = res[0] if isinstance(res, tuple) else res
    return es[COMP].copy(), net["source"].value_counts().reindex(COMP).astype(int)


def adjust(g, cov, stage=None):
    """Remove only the composition part; stage indicators (if given) protect stage effects."""
    C = (cov - cov.mean()).values
    X = [np.ones(len(cov))]
    if stage is not None:
        X += [(stage.values == s).astype(float) for s in STAGES3[1:] if (stage.values == s).any()]
    X = np.c_[tuple(X) + (C,)]
    Y = g.T.values
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    bc = beta[-C.shape[1]:]
    adj = Y - C @ bc
    comp_var = (C @ bc).var(axis=0) / Y.var(axis=0)
    return pd.DataFrame(adj.T, index=g.index, columns=g.columns), pd.Series(comp_var, index=g.index)


def main():
    gse = sys.argv[1]
    out = os.path.join(ROOT, "results", "step02", gse)
    os.makedirs(out, exist_ok=True)
    lines = [f"STEP 2 - Identity screen and composition scoring ({gse})", "=" * 55,
             "Pre-specification: addendum 1 (6ea31ea) + addendum 2 (d73e5d1)",
             f"Inputs: data/processed/{gse}_<set>_expr_genes.csv, _samples.csv (step 1)",
             "Marker genes: frozen from GSE41258 (results/step02/marker_sets_used.csv)",
             f"Seed: {SEED}",
             "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("decoupler", "pandas", "numpy")),
             ""]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    all_scores = []
    for raw, (g, lab, asets, desc_ids) in raw_sets(gse).items():
        pm, z, flags = identity_screen(g)
        qc = pm.round(3).join(z.round(2).add_suffix("_robust_z")).join(lab[["stage", "raw_label"]])
        qc["flag"] = flags
        qc.to_csv(os.path.join(out, f"tissue_identity_qc_{raw}.csv"), index_label="sample")
        fl = qc[qc.flag != ""]
        lines += [f"=== raw set {raw} ({g.shape[1]} samples) ===",
                  "Identity screen (within set; flag at 3 robust SD):",
                  "  panel medians: " + ", ".join(f"{k} {v:.2f}" for k, v in pm.median().items()),
                  f"  flagged samples: {len(fl)}"]
        for s, r in fl.iterrows():
            lines.append(f"    {s} [{r.stage}; {r.raw_label}] {r.flag}: colon {r.colon:.2f} "
                         f"(z {r.colon_robust_z:.1f}), hepatocyte z {r.hepatocyte_robust_z:.1f}, "
                         f"lung z {r.lung_robust_z:.1f}, mullerian z {r.mullerian_robust_z:.1f}")
        if "patient id" not in lab.columns and \
                lab.get("patient_id", pd.Series("", index=lab.index)).replace("", np.nan).isna().all():
            rep, thr = repeat_array_check(g, lab)
            rep.to_csv(os.path.join(out, f"repeat_array_check_{raw}.csv"), index_label="sample")
            nn = rep.r
            lines += ["Repeat-array check (no patient ids; nearest-neighbour Pearson r, expressed genes):",
                      f"  nearest-neighbour r: median {nn.median():.4f}, IQR {nn.quantile(.25):.4f}-{nn.quantile(.75):.4f}, "
                      f"max {nn.max():.4f}; flag threshold r > 0.99 and > {thr:.4f}",
                      f"  flagged: {int(rep.flag.sum())}"]
            top = rep.sort_values("r", ascending=False).head(5)
            for s, rr in top.iterrows():
                lines.append(f"    {s} [{rr.stage}] <-> {rr.nearest} [{rr.nearest_stage}] r = {rr.r:.4f}"
                             f"{'  FLAG' if rr.flag else ''}")
        centre = asets["micro_pooled"] if raw == "micro" else (asets[raw] if raw in asets else asets["main"])
        sc, sizes = scores(g, centre)
        sc.insert(0, "stage", lab["stage"])
        for c in ("normal_type", "morphology"):
            if c in lab.columns:
                sc[c] = lab[c]
        sc["analysis_sample"] = sc.index.isin(centre)
        sc.to_csv(os.path.join(out, f"composition_scores_{raw}.csv"), index_label="sample")
        all_scores.append((raw, sc))
        lines.append("Marker genes present: " + ", ".join(f"{k} {v}" for k, v in sizes.items())
                     + "  (GSE41258: Epithelial 18, Stromal 292, Immune 447)")
        miss = NET[~NET["target"].isin(g.index)][["source", "target"]]
        miss.to_csv(os.path.join(out, f"marker_genes_missing_{raw}.csv"), index=False)
        lines.append("  missing marker genes: " + ", ".join(
            f"{c} {int((miss.source == c).sum())}" for c in COMP)
            + (f"; missing epithelial: {', '.join(miss[miss.source == 'Epithelial'].target)}"
               if (miss.source == "Epithelial").any() else "; no epithelial genes missing"))
        for comp, gl in SANITY.items():
            gl = [x for x in gl if x in g.index]
            m = g.loc[gl].sub(g.loc[gl].mean(axis=1), axis=0).mean()
            lines.append(f"  sanity {comp:<10} vs mean({','.join(gl)}): Spearman r = "
                         f"{sc.loc[centre, comp].corr(m[centre], method='spearman'):.2f}")
        epi = NET.loc[NET.source == "Epithelial", "target"]
        epi = [x for x in epi if x in g.index]
        flat = g.loc[epi, centre].T.groupby(lab.loc[centre, "stage"]).median().T.reindex(columns=STAGES3)
        flat["range"] = flat.max(axis=1) - flat.min(axis=1)
        flat.round(2).to_csv(os.path.join(out, f"pan_epithelial_medians_{raw}.csv"), index_label="gene")
        lines.append(f"  frozen pan-epithelial genes, stage-median range in this set: median "
                     f"{flat['range'].median():.2f}, max {flat['range'].max():.2f} log2 "
                     f"({int((flat['range'] > 1).sum())} of {len(flat)} genes > 1 log2; descriptive, not re-selected)")
        grp = ["stage"] + (["morphology"] if raw == "micro" else [])
        tab = sc[sc.analysis_sample].groupby(grp)[COMP].median().round(2)
        lines += ["Median composition scores (analysis samples):", "  " + tab.to_string().replace("\n", "\n  ")]
        if desc_ids:
            if "desc_group" in lab.columns:
                t2 = sc.loc[desc_ids, COMP].groupby(lab.loc[desc_ids, "desc_group"]).median().round(2)
                head = "Descriptive-only arrays (SSA, set-aside A-coded carcinomas; not used for fitting):"
            else:
                t2 = sc.loc[desc_ids].groupby("morphology")[COMP].median().round(2)
                head = "Adjacent normals (descriptive only, not used for fitting):"
            lines += [head, "  " + t2.to_string().replace("\n", "\n  ")]
        for aset, ids in asets.items():
            labs = lab.loc[ids].copy()
            labs["flag_identity"] = flags[ids]
            labs.to_csv(os.path.join(PROC, f"{gse}_{aset}_samples_qc.csv"), index_label="sample")
            ga, cov, st = g[ids], sc.loc[ids, COMP], lab.loc[ids, "stage"]
            a1, v1 = adjust(ga, cov)
            a2, v2 = adjust(ga, cov, st)
            for a, suffix in ((a1, "adjusted"), (a2, "adjusted_stageaware")):
                a.index.name = "gene"
                a.to_csv(os.path.join(PROC, f"{gse}_{aset}_expr_genes_{suffix}.csv"), float_format="%.5f")
            pd.DataFrame({"R2_primary": v1, "R2_stageaware": v2}).to_csv(
                os.path.join(out, f"variance_removed_by_composition_{aset}.csv"), index_label="gene")
            cnt = st.value_counts().reindex(STAGES3).fillna(0).astype(int)
            lines.append(f"Analysis set {aset}: " + ", ".join(f"{s} {cnt[s]}" for s in STAGES3)
                         + f" | variance removed per gene, median: primary {v1.median():.3f}, "
                           f"stage-aware {v2.median():.3f}; genes > 0.5: {int((v1 > .5).sum())} / {int((v2 > .5).sum())}")
        lines.append("")

    fig, axes = plt.subplots(len(all_scores), 3, figsize=(13, 3.6 * len(all_scores)), squeeze=False)
    rng = np.random.default_rng(SEED)
    for row, (raw, sc) in enumerate(all_scores):
        d = sc.copy()
        d["group"] = d["stage"]
        if raw != "micro":
            d = d[d.analysis_sample]
        if raw == "micro":
            d["group"] = np.where(d.analysis_sample, d.stage, "Normal (adjacent)")
            d["prep"] = d["morphology"].map({"colonic crypt epithelial cells": "crypt", "mucosa": "mucosa"})
        groups = STAGES3 + (["Normal (adjacent)"] if raw == "micro" else [])
        for col, comp in enumerate(COMP):
            ax = axes[row, col]
            for i, gname in enumerate(groups):
                v = d[d.group == gname]
                colr = T.COLORS.get(gname, "#999999")
                if raw == "micro":
                    for prep, mk in (("crypt", "o"), ("mucosa", "^")):
                        vv = v[v.prep == prep]
                        ax.scatter(i + rng.uniform(-0.15, 0.15, len(vv)), vv[comp], s=18, c=colr, marker=mk,
                                   label=prep if (i == 0 and col == 0) else None)
                else:
                    ax.boxplot(v[comp], positions=[i], widths=0.5, showfliers=False)
                    ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v[comp], s=8, c=colr)
            ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, rotation=20, fontsize=8)
            ax.set_title(f"{gse} {raw}: {comp}")
        if raw == "micro":
            axes[row, 0].legend(fontsize=7, title="o crypt epithelium / ^ mucosa", title_fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "composition_by_stage.png"), dpi=150)
    lines += ["Outputs:",
              f"  results/step02/{gse}/tissue_identity_qc_<raw>.csv, composition_scores_<raw>.csv",
              f"  results/step02/{gse}/pan_epithelial_medians_<raw>.csv, variance_removed_by_composition_<set>.csv",
              f"  results/step02/{gse}/composition_by_stage.png",
              f"  data/processed/{gse}_<set>_samples_qc.csv",
              f"  data/processed/{gse}_<set>_expr_genes_adjusted.csv (primary) and _adjusted_stageaware.csv (secondary)",
              "  (unadjusted expression = step-1 file restricted to the set's samples)"]
    with open(os.path.join(out, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
