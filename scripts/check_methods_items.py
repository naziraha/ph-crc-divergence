"""
Checks requested during the Methods review (2026-10-13). Read-only: no existing result is changed;
everything is written to results/checks_methods/.

  w1        1-Wasserstein: what tda_crc_pilot.wasserstein receives, and agreement with GUDHI's Hera
            implementation (gudhi.hera.wasserstein_distance, order 1, internal_p = inf) on real diagrams.
  permdisp  Fig 4 PERMDISP recomputed with bias-adjusted distances sqrt(n/(n-1))*||x - c|| (addendum 9
            definition); same data, same seeds and permutations as step 13; Holm within cohort x version.
  symbols   gene symbols mapped to current HGNC symbols (NCBI Homo_sapiens.gene_info: nomenclature-authority
            symbol and synonyms); per-cohort counts; 117-gene replication, Fig 6 label set and arm ranking
            recomputed with harmonised symbols.
  repeat    GSE20916 near-duplicate array check (step02_composition_validation.repeat_array_check, unchanged),
            which was added to step 2 after GSE20916 had been processed.

Usage: python scripts/check_methods_items.py [w1|permdisp|symbols|repeat|all]
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import tda_crc_pilot as T  # noqa: E402

OUT = os.path.join(ROOT, "results", "checks_methods")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
PROC = os.path.join(ROOT, "data", "processed")
A9 = os.path.join(ROOT, "results", "addendum9")
os.makedirs(OUT, exist_ok=True)
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ---------------------------------------------------------------- 1. Wasserstein
def check_w1():
    from gudhi.hera import wasserstein_distance as hera_w
    say("== 1. WASSERSTEIN ==")
    df = pd.read_csv(os.path.join(EMB, "GSE41258_main_adjusted_g2000_pc10.csv"), index_col=0)
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    rng = np.random.default_rng(42)
    rows = []
    for rep in range(10):
        for normalize in (False, True):
            a = T.persistence(emb[rng.choice(np.where(df.stage == "Normal")[0], 31, replace=False)], normalize)
            b = T.persistence(emb[rng.choice(np.where(df.stage == "Adenoma")[0], 31, replace=False)], normalize)
            for dim in (0, 1):
                ours = T.wasserstein(a[dim], b[dim])
                hera = hera_w(a[dim], b[dim], order=1, internal_p=np.inf, delta=1e-6)
                rows.append(dict(rep=rep, mode="shape-only" if normalize else "raw", dim=f"H{dim}",
                                 n_points_a=len(a[dim]), n_points_b=len(b[dim]), ref_W1=ours, hera_W1=hera,
                                 rel_diff=abs(ours - hera) / max(hera, 1e-12)))
    if rep == 9:
        say("example input, H1 diagram (first 3 rows of [birth, death]):", np.round(a[1][:3], 4).tolist())
        say("example input, H0 diagram (first 3 rows of [birth, death]):", np.round(a[0][:3], 4).tolist())
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(OUT, "w1_vs_gudhi_hera.csv"), index=False)
    say(f"H0 births all zero: {bool(np.all(a[0][:, 0] == 0))}; H0 points per diagram = m - 1 = {len(a[0])}")
    say(f"comparisons: {len(r)}; max relative difference reference vs GUDHI Hera: {r.rel_diff.max():.2e}")


# ---------------------------------------------------------------- 4. PERMDISP bias-adjusted
def check_permdisp():
    import step13_benchmark as S13
    say("== 4. PERMDISP, BIAS-ADJUSTED (check only) ==")
    bench = pd.read_csv(os.path.join(ROOT, "results", "step13", "benchmark_results.csv"))
    rows = []
    for g in S13.COH:
        for v in S13.VERS:
            emb, st, _ = S13.load(g, v)
            for s, t in S13.COH[g][2]:
                A, B = emb[st == s], emb[st == t]
                n1, n2 = len(A), len(B)
                z = np.r_[np.sqrt(n1 / (n1 - 1)) * np.linalg.norm(A - A.mean(0), axis=1),
                          np.sqrt(n2 / (n2 - 1)) * np.linalg.norm(B - B.mean(0), axis=1)]
                rng = S13.rng_for(S13.bunit_name(g, v, s, t, "permdisp"))   # same seed as step 13
                x0 = np.r_[np.zeros(n1), np.ones(n2)]
                F = float(S13._anova_F(z, x0[None, :], n1, n2)[0])
                null = S13._anova_F(z, S13.perm_indicators(n1, n2, S13.N_PERM, rng), n1, n2)
                rows.append(dict(cohort=g, version=v, transition=f"{s}->{t}", D_bias=float(z[n1:].mean() / z[:n1].mean()),
                                 F_bias=F, p_bias=S13.perm_p(F, null)))
    r = pd.DataFrame(rows)
    r["p_holm_bias"] = np.nan
    for _, gg in r.groupby(["cohort", "version"]):
        r.loc[gg.index, "p_holm_bias"] = T.holm(gg.p_bias.values)
    b = bench[["cohort", "version", "transition", "permdisp_effect", "permdisp_p_holm"]] if "permdisp_effect" in bench else None
    if b is None:
        raise SystemExit("benchmark_results.csv has no permdisp_effect column")
    r = r.merge(b, on=["cohort", "version", "transition"])
    r["call_plain"] = r.permdisp_p_holm < 0.05
    r["call_bias"] = r.p_holm_bias < 0.05
    r["call_changes"] = r.call_plain != r.call_bias
    r.to_csv(os.path.join(OUT, "permdisp_bias_adjusted_check.csv"), index=False)
    for _, x in r.iterrows():
        say(f"  {x.cohort:<9} {x.version:<10} {x.transition:<21} D plain {x.permdisp_effect:.3f} -> bias {x.D_bias:.3f};"
            f" Holm p plain {x.permdisp_p_holm:.3f} -> bias {x.p_holm_bias:.3f}{'  CALL CHANGES' if x.call_changes else ''}")
    say(f"significance calls changed: {int(r.call_changes.sum())} of {len(r)}")


# ---------------------------------------------------------------- 6. gene symbols
def hgnc_map():
    gi = pd.read_csv(os.path.join(ROOT, "data", "raw", "Homo_sapiens.gene_info.gz"), sep="\t", dtype=str,
                     usecols=["GeneID", "Symbol", "Synonyms", "Symbol_from_nomenclature_authority", "type_of_gene"])
    gi["current"] = np.where(gi.Symbol_from_nomenclature_authority != "-", gi.Symbol_from_nomenclature_authority, gi.Symbol)
    current = set(gi.current)
    alias = {}
    for cur, syn, sym in zip(gi.current, gi.Synonyms, gi.Symbol):
        for a in set(([] if syn == "-" else syn.split("|")) + [sym]):
            if a != cur:
                alias.setdefault(a, set()).add(cur)
    return current, alias


def classify(symbols, current, alias):
    out = {}
    for s in symbols:
        if s in current:
            out[s] = ("current", s)
        elif s in alias and len(alias[s]) == 1:
            out[s] = ("renamed", next(iter(alias[s])))
        elif s in alias:
            out[s] = ("ambiguous alias", "")
        else:
            out[s] = ("not found", "")
    return out


def harmonise_index(df, cls):
    """Rename uniquely renamed symbols; drop a renamed row if its new name already exists in the table."""
    new = [cls[g][1] if cls[g][0] == "renamed" else g for g in df.index]
    d = df.copy()
    d.index = new
    keep = ~pd.Series(new).duplicated(keep=False).values | np.array([cls[g][0] != "renamed" for g in df.index])
    return d[keep], int((~keep).sum())


def check_symbols():
    say("== 6. GENE SYMBOLS vs CURRENT HGNC (check only) ==")
    current, alias = hgnc_map()
    say(f"reference: NCBI Homo_sapiens.gene_info.gz (data/raw), {len(current)} current symbols, {len(alias)} aliases/previous")
    files = {"GSE41258": "GSE41258_main_expr_genes.csv", "GSE20916": "GSE20916_bulk_expr_genes.csv",
             "GSE37364": "GSE37364_bulk_expr_genes.csv", "GSE117606": "GSE117606_main_expr_genes.csv",
             "TCGA-COAD": "TCGA_COAD_main_expr_genes.csv"}
    rows, cls_all = [], {}
    for c, f in files.items():
        genes = pd.read_csv(os.path.join(PROC, f), usecols=[0]).iloc[:, 0].astype(str)
        cls = classify(genes, current, alias)
        cls_all[c] = cls
        k = pd.Series([v[0] for v in cls.values()]).value_counts()
        rows.append(dict(cohort=c, genes=len(genes), current=int(k.get("current", 0)), renamed=int(k.get("renamed", 0)),
                         ambiguous_alias=int(k.get("ambiguous alias", 0)), not_found=int(k.get("not found", 0))))
        pd.DataFrame([(g, *v) for g, v in cls.items() if v[0] != "current"], columns=["symbol", "status", "current_symbol"]) \
            .to_csv(os.path.join(OUT, f"symbols_noncurrent_{c}.csv"), index=False)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, "symbols_per_cohort.csv"), index=False)
    say(tab.to_string(index=False))

    # 117-gene replication with harmonised symbols (rule of addendum9_heterogeneity.py lines 438-448)
    COH9 = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
    var = {}
    for g in COH9:
        v = pd.read_csv(os.path.join(A9, f"variability_{g}.csv"), index_col=0)
        var[g], ndrop = harmonise_index(v, classify(v.index.astype(str), current, alias))
        say(f"  {g}: variability table renamed {sum(classify(v.index.astype(str), current, alias)[x][0] == 'renamed' for x in v.index)}, "
            f"dropped as duplicates after renaming {ndrop}")
    tested = {g: set(var[g].index) for g in COH9}
    new117 = set()
    for g in COH9:
        v = var[g]
        for gene in v.index[(v.q < 0.05) & (v.V > 1)]:
            reps = [o for o in COH9 if o != g and gene in tested[o] and var[o].loc[gene, "p"] < 0.05 and var[o].loc[gene, "V"] > 1]
            if len(reps) >= 2:
                new117.add(gene)
    old = pd.read_csv(os.path.join(A9, "variability_replicated.csv"))
    old117 = set(old[old.direction == "more variable in carcinoma"].gene)
    old117_h = {classify([x], current, alias)[x][1] if classify([x], current, alias)[x][0] == "renamed" else x for x in old117}
    renamed_in_117 = sorted(x for x in old117 if classify([x], current, alias)[x][0] != "current")
    say(f"  117 list: {len(old117)} originally; {len(new117)} with harmonised symbols; "
        f"gained {sorted(new117 - old117_h)}; lost {sorted(old117_h - new117)}")
    say(f"  non-current symbols inside the 117: {renamed_in_117} -> "
        + ", ".join(f"{x}: {classify([x], current, alias)[x]}" for x in renamed_in_117))
    pd.DataFrame({"gene": sorted(new117)}).to_csv(os.path.join(OUT, "variable_genes_harmonised.csv"), index=False)

    # Fig 6 label set (fig6_variable_genes.py lines 52, 92-94)
    import fig6_variable_genes as F6
    a, _ = F6.groups()
    lab = list(a.index[a.group == "higher in carcinoma"]) + list(a.index[a.group == "lower in carcinoma"][:5]) + \
        list(a.index[a.group == "no mean shift"][:5])
    say(f"  Fig 6 labelled genes ({len(lab)}): " + ", ".join(f"{x}{'' if classify([x], current, alias)[x][0] == 'current' else ' [' + str(classify([x], current, alias)[x]) + ']'}" for x in lab))
    same_117 = new117 == old117_h
    say(f"  Fig 6 label set unchanged: {same_117 and all(classify([x], current, alias)[x][0] == 'current' for x in lab)}"
        + ("" if same_117 else " (117 list changes; label set would need re-derivation)"))

    # arm ranking with harmonised symbols (addendum11_chr20.py lines 87-94)
    import addendum11_chr20 as A
    arms, _ = A.arm_map()
    rank_old = pd.read_csv(os.path.join(ROOT, "results", "addendum11", "partA_chr", "arm_variability_ranking.csv"))
    rows = []
    for g in COH9:
        v = var[g][np.isfinite(np.log(var[g].V))]
        lv = np.log(v.V)
        d = pd.DataFrame({"logV": lv, "arm": arms.reindex(lv.index)}).dropna()
        d = d[~d.index.isin(old117_h)]
        med = d.groupby("arm").logV.agg(["median", "size"])
        med = med[med["size"] >= A.MIN_GENES].sort_values("median", ascending=False)
        med["rank"] = np.arange(1, len(med) + 1)
        o = rank_old[rank_old.cohort == g].set_index("arm")
        for arm in ("20q", "20p", "13q", "8q", "7p", "7q"):
            rows.append(dict(cohort=g, arm=arm, rank_original=int(o.loc[arm, "rank"]) if arm in o.index else None,
                             rank_harmonised=int(med.loc[arm, "rank"]) if arm in med.index else None,
                             n_genes_original=int(o.loc[arm, "n_genes"]) if arm in o.index else None,
                             n_genes_harmonised=int(med.loc[arm, "size"]) if arm in med.index else None))
        top5_old = list(o.sort_values("rank").index[:5])
        top5_new = list(med.index[:5])
        say(f"  {g}: top-5 arms original {top5_old} | harmonised {top5_new}")
    ar = pd.DataFrame(rows)
    ar.to_csv(os.path.join(OUT, "arm_ranking_harmonised.csv"), index=False)
    say(ar.to_string(index=False))


# ---------------------------------------------------------------- 2a. primary PH analysis vs symbol harmonisation
def check_primary():
    """Symbols enter the primary PH pipeline in two places only: (i) probe -> gene collapse within a cohort (two
    rows whose symbols are an old and a current name of one gene); (ii) the frozen composition marker sets, matched
    to each cohort by symbol. (i): are any colliding genes among the selected genes of the primary embeddings?
    (ii): composition scores recomputed with harmonised markers and expression, correlated with the originals."""
    import decoupler as dc
    say("== 2a. PRIMARY PH ANALYSIS vs SYMBOL HARMONISATION ==")
    current, alias = hgnc_map()
    sets = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
    rows = []
    for c, st in sets.items():
        g = pd.read_csv(os.path.join(PROC, f"{c}_{st}_expr_genes.csv"), usecols=[0]).iloc[:, 0].astype(str)
        cls = classify(g, current, alias)
        present = set(g)
        coll = [(s, cls[s][1]) for s in g if cls[s][0] == "renamed" and cls[s][1] in present]
        cg = {x for pair in coll for x in pair}
        for v in ("unadjusted", "adjusted"):
            sel = set(pd.read_csv(os.path.join(EMB, f"{c}_{st}_{v}_g2000_pc10_genes.csv")).gene)
            rows.append(dict(cohort=c, version=v, collision_pairs=len(coll), colliding_genes_in_top2000=len(cg & sel),
                             genes=";".join(sorted(cg & sel))))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, "primary_symbol_collisions.csv"), index=False)
    say(t.to_string(index=False))

    # (ii) composition markers: rebuild the GSE41258 marker sets on harmonised symbols (step02_composition.py 277-287)
    net0 = pd.read_csv(os.path.join(ROOT, "results", "step02", "marker_sets_used.csv"))
    epi = list(net0[net0.source == "Epithelial"].target)
    pg = dc.op.resource("PanglaoDB", organism="human")
    pg = pg[(pg["human"].astype(bool)) & (pg["canonical_marker"].astype(bool))]
    comps = {"Stromal": ["Fibroblasts", "Myofibroblasts", "Smooth muscle cells", "Endothelial cells", "Pericytes"],
             "Immune": ["T cells", "B cells", "Plasma cells", "NK cells", "Macrophages", "Monocytes", "Dendritic cells",
                        "Mast cells", "Neutrophils"]}

    def hz(symbols):
        cl = classify(symbols, current, alias)
        return [cl[s][1] if cl[s][0] == "renamed" else s for s in symbols]

    g41 = pd.read_csv(os.path.join(PROC, "GSE41258_main_expr_genes.csv"), index_col=0)
    g41h, _ = harmonise_index(g41, classify(g41.index.astype(str), current, alias))
    rows = [("Epithelial", x) for x in hz(epi)]
    for comp, cts in comps.items():
        rows += [(comp, x) for x in set(hz(list(pg.loc[pg["cell_type"].isin(cts), "genesymbol"].astype(str)))) if x not in epi]
    net = pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()
    sh = net["target"].value_counts()
    net = net[net["target"].isin(sh.index[sh == 1]) & net["target"].isin(g41h.index)].assign(weight=1.0)
    n0 = net0.source.value_counts()
    n1 = net.source.value_counts()
    say("  marker sets, original -> harmonised (GSE41258): " + ", ".join(f"{k} {n0.get(k, 0)} -> {n1.get(k, 0)}" for k in ("Epithelial", "Stromal", "Immune")))
    sc_files = {"GSE41258": os.path.join(ROOT, "results", "step02", "composition_scores_main.csv"),
                "GSE20916": os.path.join(ROOT, "results", "step02", "GSE20916", "composition_scores_bulk.csv"),
                "GSE37364": os.path.join(ROOT, "results", "step02", "GSE37364", "composition_scores_bulk.csv"),
                "GSE117606": os.path.join(ROOT, "results", "step02", "GSE117606", "composition_scores_all.csv")}
    rows = []
    for c, st in sets.items():
        g = pd.read_csv(os.path.join(PROC, f"{c}_{st}_expr_genes.csv"), index_col=0)
        lab = pd.read_csv(os.path.join(PROC, f"{c}_{st}_samples_qc.csv"), index_col=0)
        g = g[lab.index]
        gh, _ = harmonise_index(g, classify(g.index.astype(str), current, alias))
        old = pd.read_csv(sc_files[c], index_col=0).loc[lab.index]
        found0 = net0[net0.target.isin(g.index)].source.value_counts()
        found1 = net[net.target.isin(gh.index)].source.value_counts()
        cen = gh.sub(gh.mean(axis=1), axis=0)
        res = dc.mt.ulm(cen.T, net[net.target.isin(gh.index)][["source", "target", "weight"]], tmin=5)
        es = res[0] if isinstance(res, tuple) else res
        for k in ("Epithelial", "Stromal", "Immune"):
            rows.append(dict(cohort=c, score=k, markers_found_original=int(found0.get(k, 0)), markers_found_harmonised=int(found1.get(k, 0)),
                             pearson_r=float(np.corrcoef(old[k], es.loc[lab.index, k])[0, 1])))
    t2 = pd.DataFrame(rows)
    t2.to_csv(os.path.join(OUT, "primary_composition_markers.csv"), index=False)
    say(t2.to_string(index=False))


# ---------------------------------------------------------------- 9. GSE20916 repeat-array check
def check_repeat():
    import step02_composition_validation as S2
    say("== 9. GSE20916 NEAR-DUPLICATE ARRAY CHECK (function unchanged) ==")
    for raw in ("bulk", "micro"):
        g = pd.read_csv(os.path.join(PROC, f"GSE20916_{raw}_expr_genes.csv"), index_col=0)
        lab = pd.read_csv(os.path.join(PROC, f"GSE20916_{raw}_samples.csv"), index_col=0)
        rep, thr = S2.repeat_array_check(g, lab)
        rep.to_csv(os.path.join(OUT, f"repeat_array_check_GSE20916_{raw}.csv"), index_label="sample")
        nn = rep.r
        say(f"  {raw} ({g.shape[1]} arrays): nearest-neighbour r median {nn.median():.4f}, max {nn.max():.4f}; "
            f"threshold r > 0.99 and > {thr:.4f}; flagged {int(rep.flag.sum())}")
        for s, x in rep.sort_values("r", ascending=False).head(3).iterrows():
            say(f"    {s} [{x.stage}] <-> {x.nearest} [{x.nearest_stage}] r = {x.r:.4f}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name, f in (("w1", check_w1), ("permdisp", check_permdisp), ("symbols", check_symbols), ("repeat", check_repeat),
                    ("primary", check_primary)):
        if which in (name, "all"):
            f()
            say("")
    with open(os.path.join(OUT, f"summary_{which}.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
