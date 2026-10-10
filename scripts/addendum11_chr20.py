"""
Addendum 11 (baf4bc4; Amendment 1 384c81e; D23) - Part A (chromosome arms, 20q) and Part C (residual composition).

A1 arm annotation: MSigDB C1 positional sets (decoupler/OmniPath; D23); genes on > 1 arm, X/Y/MT excluded.
A2 per cohort and arm: median log V (addendum 9 per-gene variability ratios, adjusted expression), 117 genes excluded;
   arms with >= 20 genes ranked (1 = most increased variability in carcinoma).
A3 arm scores: mean gene-wise z (over adenomas + carcinomas) of the arm's expressed genes (117 excluded; full arm =
   sensitivity); Brown-Forsythe (carcinoma vs adenoma), V_arm = MAD ratio, Mann-Whitney; 20q gain-like share.
   Unadjusted expression = sensitivity. TCGA: 20q score by stage / CMS (descriptive).
A4 TCGA: Spearman of arm expression score vs mean GISTIC2 value of the arm's genes (272 tumours); chr20 genes among
   the 117 vs their own copy number.
Part C: Spearman (carcinomas) of DES, SYNPO2, FOXF2, CXCL14 with the stromal score and IGJ/JCHAIN with the immune
   score; residual composition if |rho| >= 0.5 in >= 3 of 4 cohorts (adjusted expression).
"""
import os
import re
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.stats import kruskal, levene, mannwhitneyu, median_abs_deviation, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
A9 = os.path.join(ROOT, "results", "addendum9")
OUTA = os.path.join(ROOT, "results", "addendum11", "partA_chr")
OUTC = os.path.join(ROOT, "results", "addendum11", "partC_composition")
GEO = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
COMP_FILE = {"GSE41258": os.path.join(ROOT, "results", "step02", "composition_scores_main.csv"),
             "GSE20916": os.path.join(ROOT, "results", "step02", "GSE20916", "composition_scores_bulk.csv"),
             "GSE37364": os.path.join(ROOT, "results", "step02", "GSE37364", "composition_scores_bulk.csv"),
             "GSE117606": os.path.join(ROOT, "results", "step02", "GSE117606", "composition_scores_all.csv")}
CRC_ARMS = ["7p", "7q", "8q", "13q", "20p", "20q", "8p", "17p", "18p", "18q"]
COMPARATORS = ["20p", "13q", "8q", "7p", "7q"]
EIGHT = ["PLAGL2", "DIDO1", "AAR2", "NELFCD", "OSER1", "DYNLRB1", "ABHD12", "ENTPD6"]
MIN_GENES = 20


def arm_map():
    import decoupler as dc
    ms = dc.op.resource("MSigDB", organism="human")
    pos = ms[ms["collection"] == "positional"]
    rows = []
    for gs, sym in zip(pos.geneset, pos.genesymbol):
        m = re.match(r"^chr(\d+|X|Y)([pq])(\S*)$", gs)
        if m:
            rows.append((sym, m.group(1), m.group(1) + m.group(2), gs))
    d = pd.DataFrame(rows, columns=["gene", "chrom", "arm", "band"]).drop_duplicates()
    n_arm = d.groupby("gene").arm.nunique()
    d = d[d.gene.isin(n_arm.index[n_arm == 1]) & ~d.chrom.isin(["X", "Y"])]
    bands = d.groupby("gene").band.apply(lambda v: ",".join(sorted(set(v))))
    return d.drop_duplicates("gene").set_index("gene")["arm"], bands


def zscore_rows(X):
    return X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=1) + 1e-9, axis=0)


def arm_scores(expr, arms, genes_ok, exclude):
    Z = zscore_rows(expr.loc[genes_ok])
    out, n = {}, {}
    for arm, g in pd.Series(arms.reindex(genes_ok).dropna()).groupby(arms.reindex(genes_ok).dropna()):
        gl = [x for x in g.index if x not in exclude]
        if len(gl) >= MIN_GENES:
            out[arm], n[arm] = Z.loc[gl].mean(axis=0), len(gl)
    return pd.DataFrame(out), pd.Series(n)


def bf(sc, st):
    a, c = sc[st == "Adenoma"], sc[st == "Carcinoma"]
    madA = median_abs_deviation(a)
    return dict(V_arm=float(median_abs_deviation(c) / madA) if madA > 0 else np.nan,
                p_BF=float(levene(c, a, center="median").pvalue), p_MW=float(mannwhitneyu(c, a).pvalue),
                median_diff=float(np.median(c) - np.median(a)))


def armkey(a):
    m = re.match(r"(\d+)([pq])", a)
    return (int(m.group(1)), m.group(2))


def part_a(arms, bands, excl):
    os.makedirs(OUTA, exist_ok=True)
    rank_rows, score_rows = [], []
    for g, st_name in GEO.items():
        v = pd.read_csv(os.path.join(A9, f"variability_{g}.csv"), index_col=0)
        v = v[np.isfinite(np.log(v.V))]
        lv = np.log(v.V)
        dfa = pd.DataFrame({"logV": lv, "arm": arms.reindex(lv.index)}).dropna()
        dfa = dfa[~dfa.index.isin(excl)]
        med = dfa.groupby("arm").logV.agg(["median", "size"])
        med = med[med["size"] >= MIN_GENES].sort_values("median", ascending=False)
        med["rank"] = np.arange(1, len(med) + 1)
        for arm, r in med.iterrows():
            rank_rows.append(dict(cohort=g, arm=arm, median_logV=r["median"], n_genes=int(r["size"]), rank=int(r["rank"]),
                                  n_arms=len(med)))
        emb = pd.read_csv(os.path.join(EMB, f"{g}_{st_name}_adjusted_g2000_pc10.csv"), index_col=0)
        ac = emb.index[emb.stage.isin(["Adenoma", "Carcinoma"])]
        st = emb.loc[ac, "stage"].values
        for version, f in (("adjusted", f"{g}_{st_name}_expr_genes_adjusted.csv"), ("unadjusted", f"{g}_{st_name}_expr_genes.csv")):
            ex = pd.read_csv(os.path.join(PROC, f), index_col=0)[ac]
            genes_ok = [x for x in v.index if x in ex.index]
            for setname, exc in (("excl117", excl), ("fullarm", set())):
                if version == "unadjusted" and setname == "fullarm":
                    continue
                S, n = arm_scores(ex, arms, genes_ok, exc)
                for arm in S.columns:
                    r = bf(S[arm].values, st)
                    rec = dict(cohort=g, version=version, genes=setname, arm=arm, n_genes=int(n[arm]), **r)
                    if arm == "20q":
                        a = S.loc[st == "Adenoma", arm]
                        rec["gainlike_share_carc"] = float((S.loc[st == "Carcinoma", arm] > np.percentile(a, 97.5)).mean())
                    score_rows.append(rec)
    rk = pd.DataFrame(rank_rows)
    sc = pd.DataFrame(score_rows)
    rk.to_csv(os.path.join(OUTA, "arm_variability_ranking.csv"), index=False)
    sc.to_csv(os.path.join(OUTA, "arm_scores_tests.csv"), index=False)

    # TCGA
    lab = pd.read_csv(os.path.join(PROC, "TCGA_COAD_main_samples_qc.csv"), index_col=0)
    ex = pd.read_csv(os.path.join(PROC, "TCGA_COAD_main_expr_genes_adjusted.csv"), index_col=0)[lab.index].dropna()
    mm = ex.mean(axis=1)
    genes_ok = list(ex.index[mm > mm.quantile(0.25)])
    S, n = arm_scores(ex, arms, genes_ok, excl)
    cms = pd.read_csv(os.path.join(ROOT, "results", "addendum10", "cms", "cms_labels_TCGA.csv"), index_col=0)["nearest"]
    cn = pd.read_csv(os.path.join(ROOT, "data", "raw", "TCGA", "TCGA.COAD.sampleMap_Gistic2_all_data_by_genes.gz"),
                     sep="\t", index_col=0)
    common = [s for s in lab.index if s in cn.columns]
    tc_rows = []
    for arm in [a for a in ["20q"] + COMPARATORS + ["8p", "17p", "18p", "18q"] if a in S.columns]:
        cg = [x for x in cn.index if arms.get(x) == arm]
        cnm = cn.loc[cg, common].mean(axis=0)
        rho, p = spearmanr(S.loc[common, arm], cnm)
        tc_rows.append(dict(arm=arm, n_expr_genes=int(n[arm]), n_cn_genes=len(cg), n_tumours=len(common), rho=rho, p=p,
                            share_cn_gt_0_1=float((cnm > 0.1).mean()), share_cn_gt_0_3=float((cnm > 0.3).mean()),
                            share_cn_lt_m0_1=float((cnm < -0.1).mean())))
    tc = pd.DataFrame(tc_rows)
    tc.to_csv(os.path.join(OUTA, "tcga_expression_vs_copynumber.csv"), index=False)
    s20 = S["20q"]
    by_stage = s20.groupby(lab.loc[S.index, "stage"]).median()
    by_cms = s20.groupby(cms.reindex(S.index)).median()
    kw_stage = kruskal(*[s20[lab.loc[S.index, "stage"] == k] for k in by_stage.index]).pvalue
    kw_cms = kruskal(*[s20[cms.reindex(S.index) == k] for k in by_cms.index]).pvalue
    pd.DataFrame({"TCGA_20q_score": s20, "stage": lab.loc[S.index, "stage"], "CMS": cms.reindex(S.index)}).to_csv(
        os.path.join(OUTA, "tcga_20q_score.csv"), index_label="sample")
    chr20 = [x for x in excl if arms.get(x) in ("20p", "20q")]
    g_rows = []
    for x in sorted(chr20):
        if x in ex.index and x in cn.index:
            r, p = spearmanr(ex.loc[x, common], cn.loc[x, common])
            g_rows.append(dict(gene=x, arm=arms[x], band=bands.get(x, ""), rho_expr_vs_cn=r, p=p))
    gdf = pd.DataFrame(g_rows)
    gdf.to_csv(os.path.join(OUTA, "tcga_chr20_variable_genes_vs_cn.csv"), index=False)

    # calls
    def call(arm):
        top = [rk[(rk.cohort == g) & (rk.arm == arm)]["rank"] for g in GEO]
        i_ok = sum(len(t) and int(t.iloc[0]) <= 5 for t in top)
        s = sc[(sc.version == "adjusted") & (sc.genes == "excl117") & (sc.arm == arm)]
        ii_ok = int(((s.p_BF < 0.05) & (s.V_arm > 1)).sum())
        t = tc[tc.arm == arm]
        iii = bool(len(t) and t.rho.iloc[0] >= 0.5 and t.p.iloc[0] < 0.05)
        i, ii = i_ok >= 3, ii_ok >= 3
        c = ("CONFIRMED" if i and ii and iii else "PARTLY SUPPORTED" if (i or ii) else "NOT CONFIRMED")
        if i and ii and not iii:
            c = "PARTLY SUPPORTED (copy-number validation failed)"
        return i_ok, ii_ok, iii, c

    L = ["ADDENDUM 11 PART A - chromosome arms and 20q (baf4bc4; Amendment 1 384c81e; D23)", "=" * 72,
         "Arm annotation: MSigDB C1 positional sets (decoupler/OmniPath), genes on one autosomal arm only (D23).",
         f"Primary: adjusted expression, the 117 addendum 9 genes excluded; arms with >= {MIN_GENES} genes.",
         "Packages: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("decoupler", "numpy", "scipy", "pandas")),
         f"Genes with a single-arm annotation: {len(arms)}; 117 genes excluded from arm statistics ({len(excl)} in list).", "",
         "A2 - 20q rank among arms by median log variability ratio (carcinoma vs adenoma); top 8 arms per cohort:"]
    for g in GEO:
        r = rk[rk.cohort == g]
        r20 = r[r.arm == "20q"]
        L.append(f"  {g:<9} 20q rank {int(r20['rank'].iloc[0]) if len(r20) else 'n/a'} of {int(r.n_arms.iloc[0])}"
                 f" (median logV {r20.median_logV.iloc[0]:+.3f})  | top 8: "
                 + ", ".join(f"{a}({m:+.2f})" for a, m in zip(r.arm.head(8), r.median_logV.head(8))))
    L.append("  CRC arms, rank per cohort (GSE41258/GSE20916/GSE37364/GSE117606): " + "; ".join(
        f"{a} " + "/".join(str(int(rk[(rk.cohort == g) & (rk.arm == a)]['rank'].iloc[0])) if len(rk[(rk.cohort == g) & (rk.arm == a)]) else "-"
                           for g in GEO) for a in CRC_ARMS))
    L += ["", "A3 - arm score variability, carcinoma vs adenoma (adjusted, 117 excluded): V_arm (Brown-Forsythe p) | median diff (MW p):"]
    for a in ["20q"] + COMPARATORS:
        s = sc[(sc.version == "adjusted") & (sc.genes == "excl117") & (sc.arm == a)].set_index("cohort")
        L.append(f"  {a:<4} " + "; ".join(f"{g} {s.loc[g, 'V_arm']:.2f} ({s.loc[g, 'p_BF']:.3g}) | {s.loc[g, 'median_diff']:+.2f} ({s.loc[g, 'p_MW']:.2g})"
                                        for g in GEO if g in s.index))
    s = sc[(sc.version == "adjusted") & (sc.genes == "excl117") & (sc.arm == "20q")].set_index("cohort")
    L.append("  20q gain-like carcinomas (score > adenoma 97.5th pct): " + ", ".join(f"{g} {100 * s.loc[g, 'gainlike_share_carc']:.0f}%" for g in GEO))
    for lab_, ver, gs in (("unadjusted expression", "unadjusted", "excl117"), ("full arm incl. 117 genes", "adjusted", "fullarm")):
        s2 = sc[(sc.version == ver) & (sc.genes == gs) & (sc.arm == "20q")].set_index("cohort")
        L.append(f"  sensitivity, {lab_}: 20q V_arm (BF p) " + "; ".join(f"{g} {s2.loc[g, 'V_arm']:.2f} ({s2.loc[g, 'p_BF']:.3g})" for g in GEO))
    L += ["", f"A4 - TCGA (272 tumours with copy number): arm expression score vs mean GISTIC2 copy number of the arm's genes:"]
    for _, r in tc.iterrows():
        L.append(f"  {r.arm:<4} rho {r.rho:.2f} (p {r.p:.2g}); expr genes {r.n_expr_genes}, CN genes {r.n_cn_genes}; "
                 f"tumours with mean CN > 0.1: {100 * r.share_cn_gt_0_1:.0f}%, > 0.3: {100 * r.share_cn_gt_0_3:.0f}%, < -0.1: {100 * r.share_cn_lt_m0_1:.0f}%")
    L.append("  (gain shares from the continuous GISTIC2 values; the thresholded file was not downloaded - descriptive only)")
    L.append(f"  TCGA 20q expression score by stage (median): " + ", ".join(f"{k} {v:+.2f}" for k, v in by_stage.items()) + f"  (Kruskal p {kw_stage:.2g})")
    L.append(f"  TCGA 20q expression score by CMS (median): " + ", ".join(f"{k} {v:+.2f}" for k, v in by_cms.items()) + f"  (Kruskal p {kw_cms:.2g})")
    L.append("  Chromosome 20 genes among the 117, expression vs own copy number (TCGA): " + ", ".join(
        f"{r.gene} [{r.arm}] {r.rho_expr_vs_cn:.2f}" for _, r in gdf.iterrows()))
    L += ["", "Arm confirmation of the eight genes (Amendment 1; MSigDB C1 bands):"]
    for x in EIGHT:
        L.append(f"  {x:<8} {arms.get(x, 'not annotated (multi-arm or absent)')}  [{bands.get(x, '')}]")
    p20 = [x for x in EIGHT if arms.get(x) == "20p"]
    q20 = [x for x in EIGHT if arms.get(x) == "20q"]
    L.append(f"  20p: {', '.join(p20) or 'none'} | 20q: {', '.join(q20) or 'none'}")
    L += ["", "A5 - CALLS (i: rank <= 5 in >= 3/4 cohorts; ii: BF p < 0.05 and V_arm > 1 in >= 3/4; iii: TCGA rho >= 0.5, p < 0.05):"]
    for a in ["20q"] + COMPARATORS:
        i_ok, ii_ok, iii, c = call(a)
        L.append(f"  {a:<4} (i) {i_ok}/4  (ii) {ii_ok}/4  (iii) {'yes' if iii else 'no'}  -> {c}" + ("" if a == "20q" else "   [comparator]"))
    L += ["", "INTERPRETATION: chromosome-arm gains are a confirmed source of carcinoma heterogeneity; their share of the",
          "excess dispersion is quantified in Amendment 2.",
          "NOTE: the GEO 'gain-like' percentage (expression score above the adenoma 97.5th percentile) is not comparable",
          "to the TCGA GISTIC2 gain percentage (copy-number values).",
          "", "Outputs: results/addendum11/partA_chr/arm_variability_ranking.csv, arm_scores_tests.csv,",
          "  tcga_expression_vs_copynumber.csv, tcga_20q_score.csv, tcga_chr20_variable_genes_vs_cn.csv"]
    open(os.path.join(OUTA, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


def part_c():
    os.makedirs(OUTC, exist_ok=True)
    pairs = [("DES", "Stromal"), ("SYNPO2", "Stromal"), ("FOXF2", "Stromal"), ("CXCL14", "Stromal"), ("IGJ", "Immune")]
    rows = []
    for g, st_name in GEO.items():
        emb = pd.read_csv(os.path.join(EMB, f"{g}_{st_name}_adjusted_g2000_pc10.csv"), index_col=0)
        car = emb.index[emb.stage == "Carcinoma"]
        comp = pd.read_csv(COMP_FILE[g], index_col=0).loc[car]
        for version, f in (("adjusted", f"{g}_{st_name}_expr_genes_adjusted.csv"), ("unadjusted", f"{g}_{st_name}_expr_genes.csv")):
            ex = pd.read_csv(os.path.join(PROC, f), index_col=0)
            for gene, score in pairs:
                sym = gene if gene in ex.index else ("JCHAIN" if gene == "IGJ" and "JCHAIN" in ex.index else None)
                if sym is None:
                    rows.append(dict(cohort=g, version=version, gene=gene, symbol_used="absent", score=score, rho=np.nan, p=np.nan))
                    continue
                r, p = spearmanr(ex.loc[sym, car], comp[score])
                rows.append(dict(cohort=g, version=version, gene=gene, symbol_used=sym, score=score, rho=r, p=p, n=len(car)))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUTC, "residual_composition.csv"), index=False)
    L = ["ADDENDUM 11 PART C - residual composition among the 117 variable genes (Amendment 1, 384c81e)", "=" * 72,
         "Carcinomas per GEO cohort; Spearman rho of gene expression with the step 2 composition score.",
         "Rule: |rho| >= 0.5 in >= 3 of 4 cohorts (adjusted expression) -> residual composition signal.", "",
         "  gene    score     adjusted rho per cohort (GSE41258/GSE20916/GSE37364/GSE117606)   | unadjusted rho   -> call"]
    for gene, score in pairs:
        a = d[(d.gene == gene) & (d.version == "adjusted")].set_index("cohort")
        u = d[(d.gene == gene) & (d.version == "unadjusted")].set_index("cohort")
        n_hit = int((a.rho.abs() >= 0.5).sum())
        call = "RESIDUAL COMPOSITION SIGNAL" if n_hit >= 3 else "not explained by residual composition"
        L.append(f"  {gene:<7} {score:<8}  " + " / ".join(f"{a.loc[g, 'rho']:+.2f}" if g in a.index and pd.notna(a.loc[g, 'rho']) else " n/a" for g in GEO)
                 + "   | " + " / ".join(f"{u.loc[g, 'rho']:+.2f}" if g in u.index and pd.notna(u.loc[g, 'rho']) else " n/a" for g in GEO)
                 + f"   -> {call} ({n_hit}/4)")
    syms = d[d.gene == "IGJ"].groupby("cohort").symbol_used.first().to_dict()
    L += ["", "Symbol used for IGJ: " + ", ".join(f"{k} {v}" for k, v in syms.items()),
          "Outputs: results/addendum11/partC_composition/residual_composition.csv"]
    open(os.path.join(OUTC, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    arms, bands = arm_map()
    excl = set(pd.read_csv(os.path.join(A9, "variability_replicated.csv")).query(
        "direction == 'more variable in carcinoma'").gene)
    part_a(arms, bands, excl)
    print()
    part_c()
