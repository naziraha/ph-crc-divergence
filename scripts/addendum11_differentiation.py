"""
Addendum 11 (baf4bc4) Part B - differentiation axis. This file currently implements only the section 2 simulation
check of the new procedure (removing a score axis within stage and calibrating f); Part B on real data is added
after the PI has reviewed Part A.

Procedure: within each stage group, regress every PC on the score (with intercept) and keep residuals; dispersion
ratio from bias-adjusted PERMDISP distances (addendum 9 Amendment 1); f = (D_total - D_resid) / (D_total - 1);
f_rand95 = 95th percentile of f with scores permuted within stage (999 times).

Usage: python scripts/addendum11_differentiation.py validate
"""
import os
import sys
import time
import zlib
from importlib import metadata
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from addendum9_heterogeneity import loo_dist, permdisp_dist  # noqa: E402  (bias-adjusted distances)

SEED, NPERM, WORKERS = 42, 999, 6
OUT = os.path.join(ROOT, "results", "addendum11", "partB_differentiation")


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def resid_within(P, score, groups):
    R = np.empty_like(P)
    for g in np.unique(groups):
        m = groups == g
        X = np.c_[np.ones(m.sum()), score[m]]
        b, *_ = np.linalg.lstsq(X, P[m], rcond=None)
        R[m] = P[m] - X @ b
    return R


def disp_ratio(PA, PC):
    dA, dC = loo_dist(PA, np.zeros(len(PA))), loo_dist(PC, np.zeros(len(PC)))
    return dA, dC, float(dC.mean() / dA.mean())


def axis_f(PA, PC, sA, sC, rng, nperm=NPERM):
    """D_total, D_resid (with PERMDISP p on residuals), f, f_rand95 for removing the score axis within stage."""
    dA, dC, Dt = disp_ratio(PA, PC)
    _, pt, _ = permdisp_dist(dA, dC, rng, nperm)
    P = np.vstack([PA, PC])
    grp = np.r_[np.zeros(len(PA)), np.ones(len(PC))]
    s = np.r_[sA, sC]
    R = resid_within(P, s, grp)
    rA, rC = R[grp == 0], R[grp == 1]
    dAr, dCr, Dr = disp_ratio(rA, rC)
    _, pr, _ = permdisp_dist(dAr, dCr, rng, nperm)
    f = (Dt - Dr) / (Dt - 1) if Dt > 1 else np.nan
    fr = []
    for _ in range(nperm):
        sp = np.r_[rng.permutation(sA), rng.permutation(sC)]
        Rp = resid_within(P, sp, grp)
        Drp = disp_ratio(Rp[grp == 0], Rp[grp == 1])[2]
        fr.append((Dt - Drp) / (Dt - 1) if Dt > 1 else np.nan)
    return dict(D_total=Dt, p_total=pt, D_resid=Dr, p_resid=pr, f=f,
                f_rand95=float(np.nanpercentile(fr, 95)) if Dt > 1 else np.nan)


def sim_rep(task):
    scen, rep = task
    rng = rng_for(f"a11sim__{scen}__{rep:03d}")
    d, nA, nC = 10, 50, 70
    u = np.zeros(d); u[0] = 1.0
    a = np.sqrt(10 * (1.25 ** 2 - 1))              # excess spread along u giving D_total ~ 1.25
    PA = rng.standard_normal((nA, d))
    tA, tC = rng.standard_normal(nA), rng.standard_normal(nC)
    PC = rng.standard_normal((nC, d))
    if scen in ("W1", "W2"):
        PC = PC + a * np.outer(tC, u)
    if scen == "W1":                                # observed score correlates 0.8 with the true axis
        sA, sC = 0.8 * tA + 0.6 * rng.standard_normal(nA), 0.8 * tC + 0.6 * rng.standard_normal(nC)
    else:                                           # W2: score unrelated to the excess; W3: null
        sA, sC = rng.standard_normal(nA), rng.standard_normal(nC)
    r = axis_f(PA, PC, sA, sC, rng)
    r.update(scenario=scen, rep=rep)
    return r


def validate():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    tasks = [(s, r) for s in ("W1", "W2", "W3") for r in range(200)]
    with Pool(WORKERS) as pool:
        d = pd.DataFrame(pool.map(sim_rep, tasks, chunksize=5))
    d.to_csv(os.path.join(OUT, "validation_replicates.csv"), index=False)
    g = {s: d[d.scenario == s] for s in ("W1", "W2", "W3")}
    thr = 0.05 + 2.5 * np.sqrt(0.05 * 0.95 / 200)
    w1 = float(((g["W1"].f >= 0.25) & (g["W1"].f > g["W1"].f_rand95)).mean())
    w2 = float((g["W2"].f <= g["W2"].f_rand95).mean())
    w3 = float((g["W3"].p_resid < 0.05).mean())
    crit = [("W1 score-driven excess: f >= 0.25 and > f_rand95", ">= 0.80", w1, w1 >= 0.8),
            ("W2 excess unrelated to score: f <= f_rand95", ">= 0.90", w2, w2 >= 0.9),
            ("W3 null: D_resid PERMDISP rejection", f"<= {thr:.3f}", w3, w3 <= thr)]
    ok = all(c[3] for c in crit)
    L = ["ADDENDUM 11 section 2 - simulation check of the Part B axis-removal procedure", "=" * 66,
         "No real data. 10 dims; adenomas 50, carcinomas 70; carcinoma excess spread along one axis (D_total ~ 1.25);",
         "W1 score correlates 0.8 with that axis; W2 score unrelated; W3 null. 200 replicates, 999 permutations.",
         f"Seed {SEED} per replicate. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("numpy", "scipy")), ""]
    L += [f"  {'PASS' if c[3] else 'FAIL'}  {c[0]:<52} {c[1]:<8} {c[2]:.3f}" for c in crit]
    L += ["", "Median D_total / D_resid / f / f_rand95: " + "; ".join(
        f"{s} {v.D_total.median():.2f}/{v.D_resid.median():.2f}/{v.f.median():.2f}/{v.f_rand95.median():.2f}" for s, v in g.items()),
          f"OVERALL: {'ALL CRITERIA PASSED' if ok else 'FAILED - do not analyse real data for Part B'}",
          f"Run time {(time.time() - t0) / 60:.1f} min."]
    open(os.path.join(OUT, "validation_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ================================================================ Part B + Amendment 2 on real data
from scipy.stats import kruskal, levene, mannwhitneyu, median_abs_deviation, spearmanr  # noqa: E402

PROC = os.path.join(ROOT, "data", "processed")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
A9 = os.path.join(ROOT, "results", "addendum9")
GEO = {"GSE41258": ("main", 48, 181), "GSE20916": ("bulk", 45, 36), "GSE37364": ("bulk", 29, 27),
       "GSE117606": ("main", 59, 71)}
ARMS4 = ["20q", "20p", "13q", "8q"]


def resid_multi(P, S, groups, dummies=None):
    """Within each group regress P on [1, S] (+ dummies for that group); return residuals."""
    R = np.empty_like(P)
    for g in np.unique(groups):
        m = groups == g
        cols = [np.ones(m.sum())]
        if S is not None:
            cols.append(S[m].reshape(m.sum(), -1))
        if dummies is not None and g == 1:
            D = pd.get_dummies(dummies[m]).values.astype(float)
            cols.append(D[:, 1:])
        X = np.column_stack(cols)
        b, *_ = np.linalg.lstsq(X, P[m], rcond=None)
        R[m] = P[m] - X @ b
    return R


def disp_from(R, grp, cgroups=None):
    dA = loo_dist(R[grp == 0], np.zeros((grp == 0).sum()))
    dC = loo_dist(R[grp == 1], np.zeros((grp == 1).sum()) if cgroups is None else np.asarray(cgroups)[grp == 1])
    keep = np.isfinite(dC)
    return dA, dC[keep], float(dC[keep].mean() / dA.mean())


def removal(P, grp, S, Dt, rng, cms=None, nperm=NPERM, perm_within=None):
    """f for removing score matrix S (n x k, or None) within stage (+ CMS centroids if cms given)."""
    R = resid_multi(P, S, grp, cms)
    dA, dC, Dr = disp_from(R, grp, cms if cms is not None else None)
    _, pr, _ = permdisp_dist(dA, dC, rng, nperm)
    f = (Dt - Dr) / (Dt - 1) if Dt > 1 else np.nan
    fr = []
    if S is not None:
        blocks = perm_within if perm_within is not None else grp
        for _ in range(nperm):
            idx = np.arange(len(grp))
            for b in np.unique(blocks):
                m = np.where(blocks == b)[0]
                idx[m] = rng.permutation(m)
            Rp = resid_multi(P, S[idx], grp, cms)
            fr.append((Dt - disp_from(Rp, grp, cms)[2]) / (Dt - 1) if Dt > 1 else np.nan)
    return dict(D_resid=Dr, p_resid_descriptive=pr, f=f,
                f_rand95=float(np.nanpercentile(fr, 95)) if fr and Dt > 1 else np.nan, null=fr)


def diff_net(exclude, full=False):
    import decoupler as dc
    pg = dc.op.resource("PanglaoDB", organism="human")
    pg = pg[(pg["human"].astype(bool)) & (pg["canonical_marker"].astype(bool))]
    gob = set(pg.loc[pg["cell_type"] == "Goblet cells", "genesymbol"])
    ent = set(pg.loc[pg["cell_type"] == "Enterocytes", "genesymbol"])
    if not full:
        gob, ent = gob - exclude, ent - exclude
    rows = [("Differentiation", x) for x in gob | ent] + [("Goblet", x) for x in gob] + [("Enterocyte", x) for x in ent]
    net = pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()
    net["weight"] = 1.0
    return net


def ulm_scores(expr, centre_ids, net):
    import decoupler as dc
    X = expr.dropna(axis=0, how="any")
    centred = X.sub(X[centre_ids].mean(axis=1), axis=0)
    n = net[net.target.isin(X.index)]
    res = dc.mt.ulm(centred.T, n, tmin=10)
    es = res[0] if isinstance(res, tuple) else res
    return es, n.source.value_counts()


def r2_by(y, g):
    y, g = np.asarray(y, float), np.asarray(g)
    tot = ((y - y.mean()) ** 2).sum()
    return float(sum((g == k).sum() * (y[g == k].mean() - y.mean()) ** 2 for k in np.unique(g)) / tot)


def part_b():
    import addendum11_chr20 as A
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    excl = set(pd.read_csv(os.path.join(A9, "variability_replicated.csv")).query(
        "direction == 'more variable in carcinoma'").gene)
    net, net_full = diff_net(excl), diff_net(excl, full=True)
    arms, _ = A.arm_map()
    cms_all = pd.read_csv(os.path.join(A9, "cms_labels.csv"), index_col=0)["nearest"]
    b2, b3, b4, setsz = [], [], [], []
    score_tab = []
    for g, (stn, nA, nC) in GEO.items():
        emb = pd.read_csv(os.path.join(EMB, f"{g}_{stn}_adjusted_g2000_pc10.csv"), index_col=0)
        ac = emb.index[emb.stage.isin(["Adenoma", "Carcinoma"])]
        if (emb.loc[ac, "stage"] == "Adenoma").sum() != nA or (emb.loc[ac, "stage"] == "Carcinoma").sum() != nC:
            raise SystemExit(f"Stage counts changed: {g}")
        st = emb.loc[ac, "stage"].values
        grp = (st == "Carcinoma").astype(int)
        P = emb.loc[ac, [c for c in emb.columns if c.startswith("PC")]].values
        adj = pd.read_csv(os.path.join(PROC, f"{g}_{stn}_expr_genes_adjusted.csv"), index_col=0)[ac]
        una = pd.read_csv(os.path.join(PROC, f"{g}_{stn}_expr_genes.csv"), index_col=0)[ac]
        sc, sz = ulm_scores(adj, list(ac), net)
        scf, szf = ulm_scores(adj, list(ac), net_full)
        scu, _ = ulm_scores(una, list(ac), net)
        setsz.append(dict(cohort=g, **{f"n_{k}": int(sz.get(k, 0)) for k in ("Differentiation", "Goblet", "Enterocyte")},
                          **{f"n_full_{k}": int(szf.get(k, 0)) for k in ("Differentiation", "Goblet", "Enterocyte")}))
        for lab_, S_ in (("primary", sc), ("full sets", scf), ("unadjusted", scu)):
            for col in [c for c in ("Differentiation", "Goblet", "Enterocyte") if c in S_.columns]:
                y = S_.loc[ac, col].values
                a, c = y[grp == 0], y[grp == 1]
                b2.append(dict(cohort=g, scores=lab_, score=col, V=float(median_abs_deviation(c) / median_abs_deviation(a)),
                               p_BF=float(levene(c, a, center="median").pvalue), median_diff=float(np.median(c) - np.median(a)),
                               p_MW=float(mannwhitneyu(c, a).pvalue)))
        s = sc.loc[ac, "Differentiation"].values
        sfull = scf.loc[ac, "Differentiation"].values
        score_tab.append(pd.DataFrame({"cohort": g, "stage": st, "Differentiation": s, "Differentiation_fullset": sfull,
                                       "Goblet": sc.loc[ac, "Goblet"].values if "Goblet" in sc else np.nan,
                                       "Enterocyte": sc.loc[ac, "Enterocyte"].values if "Enterocyte" in sc else np.nan}, index=ac))
        # arm scores (Part A definition)
        v = pd.read_csv(os.path.join(A9, f"variability_{g}.csv"), index_col=0)
        genes_ok = [x for x in v.index if x in adj.index]
        AS, _ = A.arm_scores(adj, arms, genes_ok, excl)
        C4 = AS.loc[ac, ARMS4].values
        # dispersion
        dA = loo_dist(P[grp == 0], np.zeros((grp == 0).sum()))
        dC = loo_dist(P[grp == 1], np.zeros((grp == 1).sum()))
        _, pt, Dt = permdisp_dist(dA, dC, rng_for(f"a11B__{g}__total"))
        informative = Dt > 1 and pt < 0.05
        rec = dict(cohort=g, D_total=Dt, p_total=pt, informative=informative)
        cmsC = cms_all.reindex(ac).values
        cm = np.where(grp == 1, cmsC, "A")
        for name, S in (("diff", s[:, None]), ("diff_fullset", sfull[:, None]), ("chr20q", AS.loc[ac, ["20q"]].values),
                        ("chr4", C4), ("diff+chr", np.column_stack([s, C4]))):
            r = removal(P, grp, S, Dt, rng_for(f"a11B__{g}__{name}"))
            rec.update({f"{name}_{k}": r[k] for k in ("D_resid", "p_resid_descriptive", "f", "f_rand95")})
        # CMS (carcinomas with labels); groups of < 3 merged out by loo_dist (NaN) as in addendum 9
        rc = removal(P, grp, None, Dt, rng_for(f"a11B__{g}__cms"), cms=cm)
        rb = removal(P, grp, s[:, None], Dt, rng_for(f"a11B__{g}__cmsdiff"), cms=cm, perm_within=np.where(grp == 1, cmsC, "A"))
        inc = rb["f"] - rc["f"]
        inc95 = float(np.nanpercentile(np.array(rb["null"]) - rc["f"], 95)) if Dt > 1 else np.nan
        rec.update(f_CMS=rc["f"], f_CMS_diff=rb["f"], incremental_f=inc, incremental_f_rand95=inc95)
        dcent = dC
        rho_per, p_per = spearmanr(dcent, np.abs(s[grp == 1] - np.median(s[grp == 1])))
        rec.update(rho_periphery=rho_per, p_periphery=p_per)
        b3.append(rec)
        y = s[grp == 1]
        kw = kruskal(*[y[cmsC[grp == 1] == k] for k in np.unique(cmsC[grp == 1])]).pvalue
        b4.append(dict(cohort=g, R2_score_CMS=r2_by(y, cmsC[grp == 1]), p_KW=kw,
                       **{f"median_{k}": float(np.median(y[cmsC[grp == 1] == k])) for k in ("CMS1", "CMS2", "CMS3", "CMS4")
                          if (cmsC[grp == 1] == k).any()}))
        print(f"[a11B] {g} done [{(time.time() - t0) / 60:.1f} min]", flush=True)
    # TCGA (descriptive: score by CMS and stage)
    lab = pd.read_csv(os.path.join(PROC, "TCGA_COAD_main_samples_qc.csv"), index_col=0)
    tadj = pd.read_csv(os.path.join(PROC, "TCGA_COAD_main_expr_genes_adjusted.csv"), index_col=0)[lab.index]
    ts, tsz = ulm_scores(tadj, list(lab.index), net)
    tcms = pd.read_csv(os.path.join(ROOT, "results", "addendum10", "cms", "cms_labels_TCGA.csv"), index_col=0)["nearest"].reindex(lab.index)
    y = ts.loc[lab.index, "Differentiation"].values
    b4.append(dict(cohort="TCGA_COAD", R2_score_CMS=r2_by(y, tcms.values),
                   p_KW=kruskal(*[y[tcms.values == k] for k in np.unique(tcms.values)]).pvalue,
                   **{f"median_{k}": float(np.median(y[tcms.values == k])) for k in ("CMS1", "CMS2", "CMS3", "CMS4")}))
    tstage = pd.Series(y, index=lab.index).groupby(lab.stage).median()
    tkw = kruskal(*[y[lab.stage.values == k] for k in tstage.index]).pvalue
    B2, B3, B4, SZ = pd.DataFrame(b2), pd.DataFrame(b3), pd.DataFrame(b4), pd.DataFrame(setsz)
    B2.to_csv(os.path.join(OUT, "B2_variability_level.csv"), index=False)
    B3.to_csv(os.path.join(OUT, "B3_dispersion_decomposition.csv"), index=False)
    B4.to_csv(os.path.join(OUT, "B4_cms_relation.csv"), index=False)
    SZ.to_csv(os.path.join(OUT, "marker_set_sizes.csv"), index=False)
    pd.concat(score_tab).to_csv(os.path.join(OUT, "differentiation_scores_GEO.csv"), index_label="sample")

    def calls_f(prefix):
        inf = B3[B3.informative]
        sub = int(((inf[f"{prefix}_f"] >= 0.25) & (inf[f"{prefix}_f"] > inf[f"{prefix}_f_rand95"])).sum())
        minor = int(((inf[f"{prefix}_f"] < 0.25) | (inf[f"{prefix}_f"] <= inf[f"{prefix}_f_rand95"])).sum())
        return ("SUBSTANTIAL" if sub >= 2 else "MINOR" if minor >= 2 else "INCONCLUSIVE") + f" ({sub}/{len(inf)} substantial)"
    pb = B2[(B2.scores == "primary") & (B2.score == "Differentiation")]
    axis_ok = int(((pb.p_BF < 0.05) & (pb.V > 1)).sum())
    pbf = B2[(B2.scores == "full sets") & (B2.score == "Differentiation")]
    axis_full = int(((pbf.p_BF < 0.05) & (pbf.V > 1)).sum())
    inf = B3[B3.informative]
    same = int((B4.R2_score_CMS >= 0.30).sum())
    inc_small = int((inf.incremental_f < 0.10).sum())
    adds = int(((inf.incremental_f >= 0.10) & (inf.incremental_f > inf.incremental_f_rand95)).sum())
    L = ["ADDENDUM 11 PART B + AMENDMENT 2 - differentiation axis and chromosome-arm scores vs carcinoma dispersion", "=" * 78,
         "Scores: decoupler ULM, PanglaoDB human canonical Goblet cells + Enterocytes (union = Differentiation), the 117",
         "addendum 9 genes removed (primary); full sets = sensitivity. Adjusted expression centred over adenomas + carcinomas.",
         "Dispersion: adjusted 10-PC embedding, bias-adjusted PERMDISP distances; axis removal within stage; f = (D_total - D_resid)",
         "/ (D_total - 1); f_rand95 from 999 within-stage permutations of the score rows. D_resid p descriptive only (D24).",
         "Packages: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("decoupler", "numpy", "scipy")), "",
         "Marker set sizes (primary | full): " + "; ".join(
             f"{r.cohort} Diff {r.n_Differentiation}/{r.n_full_Differentiation}, Goblet {r.n_Goblet}/{r.n_full_Goblet}, "
             f"Enterocyte {r.n_Enterocyte}/{r.n_full_Enterocyte}" for _, r in SZ.iterrows()), "",
         "B2 - Differentiation score, carcinoma vs adenoma: V (Brown-Forsythe p) | median difference (Mann-Whitney p)"]
    for lab_ in ("primary", "full sets", "unadjusted"):
        for col in ("Differentiation", "Goblet", "Enterocyte"):
            d = B2[(B2.scores == lab_) & (B2.score == col)].set_index("cohort")
            if len(d):
                L.append(f"  {lab_:<10} {col:<15} " + "; ".join(f"{g} {d.loc[g, 'V']:.2f} ({d.loc[g, 'p_BF']:.2g}) | {d.loc[g, 'median_diff']:+.2f} ({d.loc[g, 'p_MW']:.2g})"
                                                          for g in GEO if g in d.index))
    L += ["", "B3 / Amendment 2 - share of excess carcinoma dispersion explained (f vs within-stage permutation 95th pct):"]
    for _, r in B3.iterrows():
        L.append(f"  {r.cohort:<9} D_total {r.D_total:.2f} (p {r.p_total:.3f}) {'informative' if r.informative else 'NOT informative'}")
        for name, lab_ in (("diff", "differentiation"), ("diff_fullset", "differentiation (full sets)"), ("chr20q", "20q score"),
                           ("chr4", "20q+20p+13q+8q"), ("diff+chr", "differentiation + 4 arms")):
            L.append(f"      {lab_:<28} f {r[name + '_f']:+.2f} vs rand95 {r[name + '_f_rand95']:+.2f}   (D_resid {r[name + '_D_resid']:.2f}, descriptive p {r[name + '_p_resid_descriptive']:.3f})")
        L.append(f"      CMS centroids: f_CMS {r.f_CMS:+.2f}; CMS + differentiation {r.f_CMS_diff:+.2f}; incremental {r.incremental_f:+.2f} vs rand95 {r.incremental_f_rand95:+.2f}")
        L.append(f"      periphery: rho(distance to carcinoma centroid, |score - median|) = {r.rho_periphery:+.2f} (p {r.p_periphery:.2g})")
    L += ["", "B4 - Differentiation score by CMS among carcinomas (R2; Kruskal p; medians CMS1/2/3/4):"]
    for _, r in B4.iterrows():
        L.append(f"  {r.cohort:<10} R2 {r.R2_score_CMS:.2f} (p {r.p_KW:.2g}); " + " / ".join(
            f"{r.get('median_' + k, np.nan):+.2f}" for k in ("CMS1", "CMS2", "CMS3", "CMS4")))
    L.append(f"  TCGA differentiation score by stage (median): " + ", ".join(f"{k} {v:+.2f}" for k, v in tstage.items()) + f" (Kruskal p {tkw:.2g})")
    L += ["", "CALLS (B5 / Amendment 2):",
          f"  Differentiation axis confirmed (BF p < 0.05 and V > 1 in >= 3/4): {axis_ok}/4 -> {'CONFIRMED' if axis_ok >= 3 else 'NOT CONFIRMED'}"
          f"   [full-set sensitivity {axis_full}/4 -> {'same' if (axis_ok >= 3) == (axis_full >= 3) else 'DIFFERENT'}]",
          f"  Differentiation share of excess dispersion: {calls_f('diff')}   [full sets: {calls_f('diff_fullset')}]",
          f"  20q score share: {calls_f('chr20q')}",
          f"  20q+20p+13q+8q share: {calls_f('chr4')}",
          f"  Differentiation + 4 arms share: {calls_f('diff+chr')}",
          f"  Same axis as subtype (R2 >= 0.30 in >= 3/5 cohorts: {same}/5; incremental f < 0.10 in >= 2 informative: {inc_small}/{len(inf)})"
          f" -> {'LARGELY THE SAME AXIS AS SUBTYPE' if same >= 3 and inc_small >= 2 else 'not met'}",
          f"  Adds beyond subtype (incremental f >= 0.10 and > rand95 in >= 2 informative: {adds}/{len(inf)}) -> {'ADDS BEYOND SUBTYPE' if adds >= 2 else 'not met'}",
          "", f"Run time {(time.time() - t0) / 60:.1f} min.",
          "Outputs: results/addendum11/partB_differentiation/B2_variability_level.csv, B3_dispersion_decomposition.csv,",
          "  B4_cms_relation.csv, marker_set_sizes.csv, differentiation_scores_GEO.csv"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "validate":
        validate()
    elif cmd == "run":
        part_b()
    else:
        raise SystemExit("Usage: addendum11_differentiation.py validate | run")
