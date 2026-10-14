"""
Cross-cohort analyses with harmonised gene symbols (PI decision 2026-10-14; deviation recorded in Supplementary
Methods S7). Original results are not changed; everything is written to results/hgnc/ (mirroring the original folder
names) and data/processed_hgnc/ (harmonised copies of the inputs).

Symbols: every symbol is mapped to the current nomenclature-authority (HGNC) symbol with the NCBI
Homo_sapiens.gene_info file already used in step 1 (check_methods_items.hgnc_map / classify): a previous symbol or
alias with exactly one current match is renamed; if the new name already exists in the same table, the renamed row is
dropped (the row that already carried the current symbol is kept); ambiguous and unknown symbols are left unchanged.

Steps
  1. harmonised inputs: GEO unadjusted + adjusted expression (primary sets), TCGA adjusted expression, TCGA GISTIC2.
  2. per-gene variability (addendum9_heterogeneity.variability, unchanged) on harmonised adjusted expression.
  3. replication rule (as addendum9_heterogeneity.py lines 436-471) and the descriptive annotation (median signed
     mean difference; CMS template class). The same code applied to the ORIGINAL inputs must reproduce the original
     variability_replicated.csv and variability_replicated_annotated.csv exactly (checked first).
  4. addendum 11 Part A (arms), Part C (residual composition) and Part B (differentiation / arm dispersion shares),
     functions unchanged, with paths pointed at the harmonised inputs and the new gene list.
  5. source data for Figures 6-7 (source_data.armscores / per_sample, unchanged).
CMS labels and the addendum 9 dispersion decomposition are not symbol-dependent across cohorts and are copied.

Usage: python scripts/hgnc_crosscohort.py
"""
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import check_methods_items as CK  # noqa: E402

PROC = os.path.join(ROOT, "data", "processed")
HP = os.path.join(ROOT, "data", "processed_hgnc")
R = os.path.join(ROOT, "results")
H = os.path.join(R, "hgnc")
EMB = os.path.join(R, "step03", "embeddings")
A9 = os.path.join(R, "addendum9")
GEO = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


CUR, ALIAS = CK.hgnc_map()


def harmonise(df):
    cls = CK.classify(df.index.astype(str), CUR, ALIAS)
    out, ndrop = CK.harmonise_index(df, cls)
    nren = sum(cls[g][0] == "renamed" for g in df.index.astype(str))
    return out, nren, ndrop


# ---------------------------------------------------------------- 1. inputs
def inputs():
    os.makedirs(HP, exist_ok=True)
    rows = []
    files = []
    for c, st in GEO.items():
        files += [f"{c}_{st}_expr_genes.csv", f"{c}_{st}_expr_genes_adjusted.csv"]
        shutil.copy2(os.path.join(PROC, f"{c}_{st}_samples_qc.csv"), HP)
    files += ["TCGA_COAD_main_expr_genes_adjusted.csv"]
    shutil.copy2(os.path.join(PROC, "TCGA_COAD_main_samples_qc.csv"), HP)
    for f in files:
        df = pd.read_csv(os.path.join(PROC, f), index_col=0)
        h, nren, ndrop = harmonise(df)
        h.index.name = df.index.name
        h.to_csv(os.path.join(HP, f), float_format="%.5f")
        rows.append(dict(file=f, genes=len(df), renamed=nren, dropped_collision=ndrop, genes_after=len(h)))
    cn = pd.read_csv(os.path.join(ROOT, "data", "raw", "TCGA", "TCGA.COAD.sampleMap_Gistic2_all_data_by_genes.gz"),
                     sep="\t", index_col=0)
    h, nren, ndrop = harmonise(cn)
    h.index.name = cn.index.name
    h.to_csv(os.path.join(HP, "TCGA_COAD_Gistic2_all_data_by_genes_hgnc.tsv.gz"), sep="\t", compression="gzip")
    rows.append(dict(file="TCGA GISTIC2 gene-level", genes=len(cn), renamed=nren, dropped_collision=ndrop, genes_after=len(h)))
    t = pd.DataFrame(rows)
    os.makedirs(H, exist_ok=True)
    t.to_csv(os.path.join(H, "harmonised_inputs.csv"), index=False)
    say("1. HARMONISED INPUTS\n" + t.to_string(index=False))


# ---------------------------------------------------------------- 2-3. variability, replication, annotation
def replicate(var, adj, labs):
    """addendum9_heterogeneity.py lines 436-471, unchanged logic."""
    COH9 = list(GEO)
    tested = {g: set(var[g].index) for g in COH9}
    rep = []
    for direction, sgn in (("more variable in carcinoma", 1), ("less variable in carcinoma", -1)):
        for g in COH9:
            v = var[g]
            disc = v[(v.q < 0.05) & (np.sign(np.log(v.V)) == sgn)]
            for gene, r in disc.iterrows():
                reps = [o for o in COH9 if o != g and gene in tested[o]
                        and var[o].loc[gene, "p"] < 0.05 and np.sign(np.log(var[o].loc[gene, "V"])) == sgn]
                rep.append(dict(gene=gene, direction=direction, discovery=g, n_replicating=len(reps),
                                replicating=";".join(reps)))
    rep = pd.DataFrame(rep)
    repl = rep[rep.n_replicating >= 2] if len(rep) else rep
    uniq = []
    for (gene, direction), v in (repl.groupby(["gene", "direction"]) if len(repl) else []):
        cs = [o for o in COH9 if gene in tested[o]]
        md = np.median([abs(var[o].loc[gene, "mean_diff"]) for o in cs])
        Vs = {o: var[o].loc[gene, "V"] for o in cs}
        r2s = []
        for o in cs:
            car = labs[o].index
            y = adj[o].loc[gene, car].values
            grp = labs[o].nearest.values
            tot = ((y - y.mean()) ** 2).sum()
            btw = sum((grp == k).sum() * (y[grp == k].mean() - y.mean()) ** 2 for k in np.unique(grp))
            r2s.append(btw / tot if tot > 0 else np.nan)
        uniq.append(dict(gene=gene, direction=direction, discovery_cohorts=";".join(sorted(set(v.discovery))),
                         max_replicating=int(v.n_replicating.max()), n_cohorts_tested=len(cs),
                         median_V=float(np.median(list(Vs.values()))), median_abs_mean_diff=float(md),
                         no_mean_shift=bool(md < 0.5), median_CMS_R2=float(np.nanmedian(r2s)),
                         **{f"V_{o}": Vs.get(o, np.nan) for o in COH9}))
    uniq = pd.DataFrame(uniq)
    if len(uniq):
        uniq = uniq.sort_values(["direction", "median_V"], ascending=[True, False])
    return uniq


def annotate(uniq, var, harmonise_templates):
    tpl = pd.read_csv(os.path.join(A9, "cache", "templates.csv"))
    if harmonise_templates:
        cls = CK.classify(tpl.symbol.astype(str), CUR, ALIAS)
        tpl["symbol"] = [cls[s][1] if cls[s][0] == "renamed" else s for s in tpl.symbol.astype(str)]
    cls_of = tpl.drop_duplicates("symbol").set_index("symbol")["class"]
    a = uniq.copy()
    a["med_mean_diff"] = [float(np.median([var[o].loc[g, "mean_diff"] for o in GEO if g in var[o].index])) for g in a.gene]
    a["CMS_template"] = [cls_of.get(g, np.nan) for g in a.gene]
    return a


def load_labs():
    cms = pd.read_csv(os.path.join(A9, "cms_labels.csv"), index_col=0)
    return {c: cms[cms.cohort == c] for c in GEO}


def check_reproduction(labs):
    var = {c: pd.read_csv(os.path.join(A9, f"variability_{c}.csv"), index_col=0) for c in GEO}
    adj = {c: pd.read_csv(os.path.join(PROC, f"{c}_{st}_expr_genes_adjusted.csv"), index_col=0) for c, st in GEO.items()}
    u = replicate(var, adj, labs)
    a = annotate(u, var, harmonise_templates=False)
    o = pd.read_csv(os.path.join(A9, "variability_replicated.csv"))
    oa = pd.read_csv(os.path.join(A9, "variability_replicated_annotated.csv"))
    same_genes = list(u.gene) == list(o.gene)
    num = [c for c in o.columns if o[c].dtype != object and c in u.columns]
    dmax = max(float(np.nanmax(np.abs(u[c].values.astype(float) - o[c].values.astype(float)))) for c in num)
    dmd = float(np.max(np.abs(a.med_mean_diff.values - oa.med_mean_diff.values)))
    tpl_same = (a.CMS_template.fillna("").values == oa.CMS_template.fillna("").values).all()
    ok = same_genes and dmax < 1e-9 and dmd < 1e-9 and tpl_same
    say(f"3a. REPRODUCTION on original inputs: genes identical {same_genes} ({len(u)}); max numeric difference {dmax:.1e}; "
        f"med_mean_diff {dmd:.1e}; CMS_template identical {bool(tpl_same)} -> {'OK' if ok else 'MISMATCH'}")
    if not ok:
        raise SystemExit("re-implementation does not reproduce the original files")


def variability_and_replication(labs):
    import addendum9_heterogeneity as A9m
    out = os.path.join(H, "addendum9")
    os.makedirs(out, exist_ok=True)
    var, adj = {}, {}
    for c, st in GEO.items():
        emb = pd.read_csv(os.path.join(EMB, f"{c}_{st}_adjusted_g2000_pc10.csv"), index_col=0)
        adj[c] = pd.read_csv(os.path.join(HP, f"{c}_{st}_expr_genes_adjusted.csv"), index_col=0)[emb.index]
        var[c] = A9m.variability(c, emb, adj[c])
        var[c].to_csv(os.path.join(out, f"variability_{c}.csv"), index_label="gene")
    u = replicate(var, adj, labs)
    u.to_csv(os.path.join(out, "variability_replicated.csv"), index=False)
    a = annotate(u, var, harmonise_templates=True)
    a.to_csv(os.path.join(out, "variability_replicated_annotated.csv"), index=False)
    for f in ("cms_labels.csv", "dispersion_decomposition.csv", "cms_proportions.csv"):
        shutil.copy2(os.path.join(A9, f), out)
    o = set(pd.read_csv(os.path.join(A9, "variability_replicated.csv")).query("direction == 'more variable in carcinoma'").gene)
    n = set(u[u.direction == "more variable in carcinoma"].gene)
    o_h = {CK.classify([g], CUR, ALIAS)[g][1] if CK.classify([g], CUR, ALIAS)[g][0] == "renamed" else g for g in o}
    say(f"3b. HARMONISED: genes tested " + ", ".join(f"{c} {len(var[c])}" for c in GEO)
        + f"; more variable in carcinoma {len(n)} (original {len(o)}); less variable "
        + f"{int((u.direction == 'less variable in carcinoma').sum())}")
    say(f"    added: {sorted(n - o_h)}; lost: {sorted(o_h - n)}; renamed originals: "
        + ", ".join(f"{g}->{CK.classify([g], CUR, ALIAS)[g][1]}" for g in sorted(o) if CK.classify([g], CUR, ALIAS)[g][0] == "renamed"))
    return n


# ---------------------------------------------------------------- 4. addendum 11
def addendum11(excl):
    import addendum11_chr20 as A
    import addendum11_differentiation as D
    A.A9, A.PROC = os.path.join(H, "addendum9"), HP
    A.OUTA, A.OUTC = os.path.join(H, "addendum11", "partA_chr"), os.path.join(H, "addendum11", "partC_composition")
    A.CN_FILE = os.path.join(HP, "TCGA_COAD_Gistic2_all_data_by_genes_hgnc.tsv.gz")
    arms, bands = A.arm_map()
    A.part_a(arms, bands, excl)
    A.part_c()
    D.A9, D.PROC, D.OUT = os.path.join(H, "addendum9"), HP, os.path.join(H, "addendum11", "partB_differentiation")
    D.part_b()
    return arms


# ---------------------------------------------------------------- 5. source data
def source_data(arms, excl):
    import source_data as S
    S.R, S.PROC, S.OUT = H, HP, os.path.join(H, "source_data")
    os.makedirs(S.OUT, exist_ok=True)
    dist = pd.read_csv(os.path.join(R, "source_data", "stage_centroid_distances_GEO.csv"))
    arm = S.armscores(arms, excl)
    S.per_sample(arm, dist)
    ck = pd.DataFrame(S.CHECKS)
    ck.to_csv(os.path.join(S.OUT, "checks.csv"), index=False)
    say("5. SOURCE DATA CHECKS")
    for r in ck.itertuples():
        say(f"  [{'OK' if r.ok else 'MISMATCH'}] {r.check}: {r.detail}")


def relabel(nvar):
    """The addendum 11 code labels the primary gene set 'excl117'; here it excludes the harmonised list."""
    lab = f"excl{nvar}"
    for f in (os.path.join(H, "addendum11", "partA_chr", "arm_scores_tests.csv"), os.path.join(H, "source_data", "arm_scores_GEO.csv")):
        d = pd.read_csv(f)
        d["genes"] = d["genes"].replace({"excl117": lab})
        d.to_csv(f, index=False)
    say(f"   gene-set label 'excl117' replaced by '{lab}' in arm_scores_tests.csv and arm_scores_GEO.csv")


if __name__ == "__main__":
    t0 = time.time()
    labs = load_labs()
    check_reproduction(labs)
    inputs()
    excl = variability_and_replication(labs)
    arms = addendum11(excl)
    source_data(arms, excl)
    relabel(len(excl))
    say(f"Run time {(time.time() - t0) / 60:.1f} min")
    with open(os.path.join(H, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("CROSS-COHORT ANALYSES WITH HARMONISED GENE SYMBOLS (scripts/hgnc_crosscohort.py)\n" + "=" * 72 + "\n"
                 + "\n".join(LOG) + "\n")
