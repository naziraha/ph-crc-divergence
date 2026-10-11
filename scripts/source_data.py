"""
Source-data tables for the paper (GEO cohorts only; same sample sets as before). No new statistical tests, no figures,
no existing file changed. Everything is written to results/source_data/.

1. ph_diagrams_GEO.csv / ph_subsample_members_GEO.csv - persistence diagram birth/death values, recomputed from the
   frozen 10-PC embeddings with the same seeded draws as the addendum 12 / paper figures
   (rng 'a12fig__<cohort>__<version>__<mode>', k = 20 subsamples of m per stage).
2. validation_summary_primary.csv - step 5 primary results per cohort x transition with a 'rule_met' column
   (shape-only Holm p < 0.05 in both unadjusted and adjusted), plus raw, H0-only / H1-only (addendum 12) and the
   H1 landscape Holm p (addendum 12), copied from the saved result files.
3. stage_centroid_distances_GEO.csv - each sample's Euclidean distance to its own stage centroid in the 10-PC
   embedding (step 13 PERMDISP definition) and the bias-adjusted distance sqrt(n/(n-1))*||x - c|| (addendum 9).
4. arm_scores_GEO.csv - per-sample arm scores (20q, 20p, 13q, 8q) for adenomas + carcinomas, computed with
   addendum11_chr20.arm_map / arm_scores exactly as in Part A (z over adenomas + carcinomas; the 117 genes excluded
   = primary; full arm = sensitivity; adjusted = primary, unadjusted = sensitivity).
5. per_sample_<cohort>.csv - sample, stage, patient id, CMS labels (carcinomas), differentiation scores, arm scores,
   centroid distances, and the 117 variable genes (unadjusted log2 'expr_' and composition-adjusted 'adj_').
Checks (checks.csv): recomputed values are compared with the saved results.

Usage: python scripts/source_data.py
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import tda_crc_pilot as T  # noqa: E402
import addendum12_ph_outputs as A12  # noqa: E402
import addendum11_chr20 as A11  # noqa: E402

R = os.path.join(ROOT, "results")
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(R, "source_data")
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
SETS = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
VERS = ("unadjusted", "adjusted")
ARMS4 = ["20q", "20p", "13q", "8q"]
STEP5 = {"GSE41258": "step05/permutation_results_all.csv", "GSE20916": "step05/GSE20916/permutation_results_all.csv",
         "GSE37364": "step05/GSE37364/permutation_results_all.csv",
         "GSE117606": "step05/GSE117606/permutation_results_all.csv"}
CHECKS = []


def check(name, ok, detail):
    CHECKS.append(dict(check=name, ok=bool(ok), detail=detail))
    print(f"[{'OK' if ok else 'MISMATCH'}] {name}: {detail}", flush=True)


def load_emb(c, v):
    emb, st, stages, m = A12.load12(c, v)          # raises if stage counts changed
    df = pd.read_csv(os.path.join(A12.EMB, f"{c}_{SETS[c]}_{v}_g2000_pc10.csv"), index_col=0)
    return df, emb, st, stages, m


# ---------------------------------------------------------------- 1. persistence diagrams
def diagrams():
    rows, mem = [], []
    n_same = n_tot = 0
    for c in GEO:
        for v in VERS:
            df, emb, st, stages, m = load_emb(c, v)
            for mode in ("shape", "raw"):
                rng = A12.rng_for(f"a12fig__{c}__{v}__{mode}")
                rng2 = A12.rng_for(f"a12fig__{c}__{v}__{mode}")
                for s in stages:
                    idx = np.where(st == s)[0]
                    ref = A12.draw(emb, idx, m, mode == "shape", rng2)        # the figure draws, for the check
                    for k in range(A12.K):
                        take = rng.choice(idx, size=m, replace=False)          # same sequence as A12.draw
                        dg = T.persistence(emb[take], mode == "shape")
                        mem.append(dict(cohort=c, version=v, mode="shape-only" if mode == "shape" else "raw", stage=s,
                                        subsample=k + 1, m=m, samples=";".join(df.index[take])))
                        for dim in (0, 1):
                            n_tot += 1
                            n_same += int(dg[dim].shape == ref[k][dim].shape and np.array_equal(dg[dim], ref[k][dim]))
                            for b, d in dg[dim]:
                                rows.append(dict(cohort=c, version=v, mode="shape-only" if mode == "shape" else "raw",
                                                 stage=s, subsample=k + 1, dim=f"H{dim}", birth=b, death=d,
                                                 persistence=d - b))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "ph_diagrams_GEO.csv"), index=False)
    pd.DataFrame(mem).to_csv(os.path.join(OUT, "ph_subsample_members_GEO.csv"), index=False)
    check("1 diagrams identical to the seeded figure draws", n_same == n_tot, f"{n_same}/{n_tot} diagrams identical")
    h0 = out[out.dim == "H0"].groupby(["cohort", "version", "mode", "stage", "subsample"]).size()
    mm = pd.DataFrame(mem).set_index(["cohort", "version", "mode", "stage", "subsample"]).m
    check("1 H0 has m-1 finite features per subsample", bool((h0 == mm.loc[h0.index] - 1).all()),
          f"{int((h0 == mm.loc[h0.index] - 1).sum())}/{len(h0)} subsamples")
    return out


# ---------------------------------------------------------------- 2. validation summary
def validation():
    dims = pd.read_csv(os.path.join(R, "addendum12", "dimension_split_all_cohorts.csv"))
    land = pd.read_csv(os.path.join(R, "addendum12", "landscape_tests.csv"))
    rows = []
    for c in GEO:
        d = pd.read_csv(os.path.join(R, STEP5[c]))
        d = d[(d.analysis.str.lower() == "primary") & d.version.isin(VERS)]
        for tr in d.transition.unique():
            g = d[d.transition == tr]
            get = lambda v, mo, col: float(g[(g.version == v) & (g["mode"] == mo)][col].iloc[0])
            r = dict(cohort=c, transition=tr, n_from=int(g.n_from.iloc[0]), n_to=int(g.n_to.iloc[0]),
                     subsample_m=int(g.subsample_size.iloc[0]))
            for mo, tag in (("shape-only", "shape"), ("raw", "raw")):
                for v in VERS:
                    r[f"{tag}_z_{v[:5]}"] = get(v, mo, "z")
                    r[f"{tag}_p_holm_{v[:5]}"] = get(v, mo, "p_holm")
            r["rule_met"] = "yes" if (r["shape_p_holm_unadj"] < 0.05 and r["shape_p_holm_adjus"] < 0.05) else "no"
            for dim in ("H0", "H1"):
                for v in VERS:
                    x = dims[(dims.cohort == c) & (dims.version == v) & (dims["mode"] == "shape-only") &
                             (dims.dim == dim) & (dims.transition == tr)]
                    r[f"shape_{dim}only_p_holm_{v[:5]}"] = float(x.p_holm.iloc[0]) if len(x) else np.nan
            for mo, tag in (("shape-only", "shape"), ("raw", "raw")):
                for v in VERS:
                    x = land[(land.cohort == c) & (land.version == v) & (land["mode"] == mo) & (land.transition == tr)]
                    r[f"H1landscape_{tag}_p_holm_{v[:5]}"] = float(x.p_holm_trans.iloc[0]) if len(x) else np.nan
            rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "validation_summary_primary.csv"), index=False)
    check("2 primary rule met (expected: none)", (out.rule_met == "no").all(),
          f"rule met in {int((out.rule_met == 'yes').sum())} of {len(out)} cohort x transition rows")
    check("2 H0/H1 and landscape values found for every row", not out.filter(like="only_p_holm").isna().any().any()
          and not out.filter(like="H1landscape").isna().any().any(), "no missing values")
    return out


# ---------------------------------------------------------------- 3. distances to stage centroid
def distances():
    rows = []
    for c in GEO:
        for v in VERS:
            df, emb, st, stages, m = load_emb(c, v)
            for s in stages:
                msk = st == s
                n = int(msk.sum())
                d = np.linalg.norm(emb[msk] - emb[msk].mean(0), axis=1)
                for smp, x in zip(df.index[msk], d):
                    rows.append(dict(cohort=c, version=v, sample=smp, stage=s, n_stage=n, dist_to_centroid=x,
                                     dist_to_centroid_bias_adj=np.sqrt(n / (n - 1)) * x))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "stage_centroid_distances_GEO.csv"), index=False)
    b = pd.read_csv(os.path.join(R, "step13", "benchmark_results.csv"))
    diffs = []
    for _, r in b.iterrows():
        a, t = r.transition.split("->")
        x = out[(out.cohort == r.cohort) & (out.version == r.version)]
        D = x[x.stage == t].dist_to_centroid.mean() / x[x.stage == a].dist_to_centroid.mean()
        diffs.append(abs(D - r.permdisp_effect))
    check("3 PERMDISP D reproduces step 13 (all rows)", max(diffs) < 1e-9,
          f"{len(diffs)} rows; max |difference| {max(diffs):.2e}")
    dec = pd.read_csv(os.path.join(R, "addendum9", "dispersion_decomposition.csv"))
    dec = dec.groupby("cohort").D_total.first()
    diffs = []
    for c in GEO:
        x = out[(out.cohort == c) & (out.version == "adjusted")]
        D = x[x.stage == "Carcinoma"].dist_to_centroid_bias_adj.mean() / x[x.stage == "Adenoma"].dist_to_centroid_bias_adj.mean()
        diffs.append(abs(D - dec[c]))
    check("3 bias-adjusted D_total reproduces addendum 9 (adjusted, A->C)", max(diffs) < 1e-9,
          f"4 cohorts; max |difference| {max(diffs):.2e}")
    return out


# ---------------------------------------------------------------- 4. arm scores
def armscores(arms, excl):
    rows, ref = [], pd.read_csv(os.path.join(R, "addendum11", "partA_chr", "arm_scores_tests.csv"))
    diffs, n_cmp = [], 0
    for c in GEO:
        st_name = SETS[c]
        var = pd.read_csv(os.path.join(R, "addendum9", f"variability_{c}.csv"), index_col=0)
        var = var[np.isfinite(np.log(var.V))]
        emb = pd.read_csv(os.path.join(A11.EMB, f"{c}_{st_name}_adjusted_g2000_pc10.csv"), index_col=0)
        ac = emb.index[emb.stage.isin(["Adenoma", "Carcinoma"])]
        st = emb.loc[ac, "stage"].values
        for version, f in (("adjusted", f"{c}_{st_name}_expr_genes_adjusted.csv"), ("unadjusted", f"{c}_{st_name}_expr_genes.csv")):
            ex = pd.read_csv(os.path.join(PROC, f), index_col=0)[ac]
            genes_ok = [x for x in var.index if x in ex.index]
            for setname, exc in (("excl117", excl), ("fullarm", set())):
                if version == "unadjusted" and setname == "fullarm":
                    continue
                S, n = A11.arm_scores(ex, arms, genes_ok, exc)
                S = S[[a for a in ARMS4 if a in S.columns]]
                tmp = S.copy()
                tmp.insert(0, "stage", st)
                tmp.insert(0, "genes", setname)
                tmp.insert(0, "version", version)
                tmp.insert(0, "cohort", c)
                rows.append(tmp.rename(columns={a: f"arm_{a}" for a in ARMS4}))
                for a in S.columns:
                    r = A11.bf(S[a].values, st)
                    x = ref[(ref.cohort == c) & (ref.version == version) & (ref.genes == setname) & (ref.arm == a)]
                    if len(x):
                        n_cmp += 1
                        diffs += [abs(r[k] - float(x[k].iloc[0])) for k in ("V_arm", "p_BF", "median_diff")]
                        diffs.append(abs(int(n[a]) - int(x.n_genes.iloc[0])))
    out = pd.concat(rows).rename_axis("sample").reset_index()
    out.to_csv(os.path.join(OUT, "arm_scores_GEO.csv"), index=False)
    check("4 arm V_arm / BF p / median diff / n genes reproduce addendum 11", n_cmp == 4 * 3 * 4 and max(diffs) < 1e-9,
          f"{n_cmp} cohort x version x arm rows compared (expected 48); max |difference| {max(diffs):.2e}")
    return out


# ---------------------------------------------------------------- 5. per-sample tables
def per_sample(arm_df, dist_df):
    genes = pd.read_csv(os.path.join(R, "addendum9", "variability_replicated_annotated.csv")).gene.tolist()
    cms = pd.read_csv(os.path.join(R, "addendum9", "cms_labels.csv"), index_col=0)
    dif = pd.read_csv(os.path.join(R, "addendum11", "partB_differentiation", "differentiation_scores_GEO.csv"), index_col=0)
    for c in GEO:
        st_name = SETS[c]
        lab = pd.read_csv(os.path.join(PROC, f"{c}_{st_name}_samples_qc.csv"), index_col=0)
        emb = pd.read_csv(os.path.join(A12.EMB, f"{c}_{st_name}_unadjusted_g2000_pc10.csv"), index_col=0)
        t = pd.DataFrame(index=emb.index)
        t.index.name = "sample"
        t["cohort"] = c
        t["stage"] = emb.stage
        if "patient_id" in lab.columns:
            t["patient_id"] = lab.patient_id.reindex(t.index)
        cc = cms[cms.cohort == c]
        for col in ("nearest", "prediction", "confident", "nearest_adjusted"):
            t[f"CMS_{col}"] = cc[col].reindex(t.index)
        for col in ("Differentiation", "Differentiation_fullset", "Goblet", "Enterocyte"):
            t[col] = dif[dif.cohort == c][col].reindex(t.index)
        a = arm_df[(arm_df.cohort == c) & (arm_df.version == "adjusted") & (arm_df.genes == "excl117")].set_index("sample")
        for arm in ARMS4:
            t[f"arm_{arm}"] = a[f"arm_{arm}"].reindex(t.index)
        for v in VERS:
            dd = dist_df[(dist_df.cohort == c) & (dist_df.version == v)].set_index("sample")
            t[f"dist_centroid_{v}"] = dd.dist_to_centroid.reindex(t.index)
        raw = pd.read_csv(os.path.join(PROC, f"{c}_{st_name}_expr_genes.csv"), index_col=0)
        adj = pd.read_csv(os.path.join(PROC, f"{c}_{st_name}_expr_genes_adjusted.csv"), index_col=0)
        raw, adj = raw[~raw.index.duplicated()], adj[~adj.index.duplicated()]
        E = pd.concat([raw.reindex(genes)[t.index].T.add_prefix("expr_"),
                       adj.reindex(genes).reindex(columns=t.index).T.add_prefix("adj_")], axis=1)
        t = pd.concat([t, E], axis=1)
        t.to_csv(os.path.join(OUT, f"per_sample_{c}.csv"))
        # checks
        nC = int((t.stage == "Carcinoma").sum())
        check(f"5 {c} CMS labels joined for all carcinomas", int(t.CMS_nearest.notna().sum()) == nC == int((cms.cohort == c).sum()),
              f"{int(t.CMS_nearest.notna().sum())} labelled / {nC} carcinomas")
        nAC = int(t.stage.isin(["Adenoma", "Carcinoma"]).sum())
        check(f"5 {c} differentiation and arm scores joined for all adenomas + carcinomas",
              int(t.Differentiation.notna().sum()) == nAC == int(t.arm_20q.notna().sum()),
              f"differentiation {int(t.Differentiation.notna().sum())}, arm_20q {int(t.arm_20q.notna().sum())}, A+C {nAC}")
        var = pd.read_csv(os.path.join(R, "addendum9", f"variability_{c}.csv"), index_col=0)
        from scipy.stats import median_abs_deviation as mad
        A, C = t[t.stage == "Adenoma"], t[t.stage == "Carcinoma"]
        dv, n_g = [], 0
        for gname in genes:
            if gname in var.index and f"adj_{gname}" in t and t[f"adj_{gname}"].notna().all():
                n_g += 1
                V = mad(C[f"adj_{gname}"]) / mad(A[f"adj_{gname}"])
                dv.append(abs(V - var.loc[gname, "V"]))
        check(f"5 {c} variability ratio V of the 117 genes reproduces addendum 9", max(dv) < 1e-9,
              f"{n_g} genes compared; max |difference| {max(dv):.2e}; "
              f"{int(t.filter(like='expr_').isna().all().sum())} of 117 genes absent from this cohort's expression")


def main():
    os.makedirs(OUT, exist_ok=True)
    arms, _ = A11.arm_map()
    excl = set(pd.read_csv(os.path.join(R, "addendum9", "variability_replicated_annotated.csv")).gene)
    diagrams()
    val = validation()
    dist = distances()
    arm = armscores(arms, excl)
    per_sample(arm, dist)
    ck = pd.DataFrame(CHECKS)
    ck.to_csv(os.path.join(OUT, "checks.csv"), index=False)
    pk = ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "pandas", "scipy", "decoupler"))
    L = ["SOURCE DATA (GEO cohorts; no new statistical tests; no figures; no existing file changed)", "=" * 72,
         "Script: scripts/source_data.py. Packages: " + pk + ". Seeds: as the source analyses (42-based named streams).",
         "Sample sets: GSE41258 main (53/48/181/39), GSE20916 bulk (24/45/36), GSE37364 bulk (38/29/27), "
         "GSE117606 main (65/59/71); stage counts verified on load.", "",
         "FILES",
         "  ph_diagrams_GEO.csv            birth/death/persistence of every H0 and H1 feature; cohort x version x mode x stage x",
         "                                 subsample (k = 20, m = 80% of the smallest stage); same seeded draws as the PH figures",
         "  ph_subsample_members_GEO.csv   which samples were in each subsample",
         "  validation_summary_primary.csv step 5 primary: shape and raw z / Holm p (unadjusted, adjusted), rule_met, H0-only and",
         "                                 H1-only shape Holm p (addendum 12), H1 landscape Holm p (addendum 12); copied values",
         "  stage_centroid_distances_GEO.csv each sample's distance to its stage centroid (10 PCs), plain and bias-adjusted",
         "  arm_scores_GEO.csv             per-sample 20q/20p/13q/8q scores, adenomas + carcinomas only (score defined as",
         "                                 z over adenomas + carcinomas); adjusted/excl117 = primary; fullarm and unadjusted = sensitivity",
         "  per_sample_<cohort>.csv        one row per sample: stage, patient id, CMS (carcinomas), differentiation, arm scores",
         "                                 (adjusted, excl117), centroid distances, the 117 genes (expr_ = unadjusted log2, adj_ = adjusted)",
         "  checks.csv                     reproduction checks against the saved results", "",
         "RULE MET (primary): " + ", ".join(f"{r.cohort} {r.transition} {r.rule_met}" for r in val.itertuples()), "",
         "CHECKS"] + [f"  [{'OK' if r.ok else 'MISMATCH'}] {r.check}: {r.detail}" for r in ck.itertuples()]
    L += ["", f"All checks passed: {bool(ck.ok.all())}"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
