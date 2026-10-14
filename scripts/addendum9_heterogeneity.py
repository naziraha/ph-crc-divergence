"""
Addendum 9 (commit 5a9b03a; D16 9590f65) - what drives the greater heterogeneity of carcinomas vs adenomas.

Part (a): per-gene Brown-Forsythe test (Levene, median-centred), carcinoma vs adenoma, BH within cohort, replication.
Part (b): CMS labels (CMScaller, r-cms environment, scripts/addendum9_cms.R); dispersion decomposition in the
stage-blind adjusted 10-PC space with leave-one-out (LOO) centroid distances: D_total, R2_CMS (PERMANOVA),
D_within (PERMDISP on within-subtype LOO distances), f and random-label f_rand95.

Usage:
  python scripts/addendum9_heterogeneity.py validate     - section 6 validation (simulated data + CMScaller example)
"""
import os
import shutil
import subprocess
import sys
import time
import zlib
from importlib import metadata
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import levene

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "addendum9")
SEED, NPERM, WORKERS = 42, 999, 6
RENV = os.environ.get("R_CMS_HOME", "")   # folder of the r-cms conda environment; empty: use Rscript on the PATH
RSCRIPT = os.path.join(RENV, "lib", "R", "bin", "x64", "Rscript.exe") if RENV else (shutil.which("Rscript") or "Rscript")


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


# ---------------------------------------------------------------- dispersion tools
def loo_dist(P, groups):
    """Amendment 1 (12dfc3c): bias-adjusted PERMDISP distance sqrt(n/(n-1)) * ||x - c|| to own-group centroid
    (vegan betadisper bias.adjust = TRUE). Name kept for continuity; it no longer uses leave-one-out."""
    d = np.full(len(P), np.nan)
    for g in np.unique(groups):
        m = groups == g
        n = m.sum()
        if n < 3:
            continue
        c = P[m].mean(0)
        d[m] = np.sqrt(n / (n - 1)) * np.linalg.norm(P[m] - c, axis=1)
    return d


def _F2(z, X, n1, n2):
    N = n1 + n2
    s2 = X @ z
    s1 = z.sum() - s2
    ss_tot = ((z - z.mean()) ** 2).sum()
    ss_b = n1 * (s1 / n1 - z.mean()) ** 2 + n2 * (s2 / n2 - z.mean()) ** 2
    return ss_b / ((ss_tot - ss_b) / (N - 2))


def permdisp_dist(dA, dC, rng, nperm=NPERM):
    """PERMDISP on given distances (adenoma dA, carcinoma dC): ANOVA F, label permutations; D = mean dC / mean dA."""
    z = np.r_[dA, dC]
    n1, n2 = len(dA), len(dC)
    X = np.zeros((nperm, n1 + n2))
    for i in range(nperm):
        X[i, rng.permutation(n1 + n2)[:n2]] = 1
    F = float(_F2(z, np.r_[np.zeros(n1), np.ones(n2)][None, :], n1, n2)[0])
    null = _F2(z, X, n1, n2)
    return F, float((1 + np.sum(null >= F)) / (1 + nperm)), float(dC.mean() / dA.mean())


def permanova_k(P, lab, rng, nperm=NPERM):
    D2 = ((P[:, None, :] - P[None, :, :]) ** 2).sum(2)
    N = len(P)
    ss_t = D2.sum() / 2 / N
    ks = np.unique(lab)

    def ssw(l):
        return sum(D2[np.ix_(l == k, l == k)].sum() / 2 / (l == k).sum() for k in ks)

    w = ssw(lab)
    a = len(ks)
    F = ((ss_t - w) / (a - 1)) / (w / (N - a))
    null = np.array([((ss_t - (ww := ssw(rng.permutation(lab)))) / (a - 1)) / (ww / (N - a)) for _ in range(nperm)])
    return float(F), float((1 + np.sum(null >= F)) / (1 + nperm)), float((ss_t - w) / ss_t)


def decompose(PA, PC, cms, rng, nperm=NPERM):
    """Section 4.3 items 1-4 for adenoma points PA, carcinoma points PC with CMS labels cms."""
    dA = loo_dist(PA, np.zeros(len(PA)))
    dCt = loo_dist(PC, np.zeros(len(PC)))
    Ft, pt, Dt = permdisp_dist(dA, dCt, rng, nperm)
    keep = pd.Series(cms).map(pd.Series(cms).value_counts()).values >= 3
    PCk, cmsk = PC[keep], np.asarray(cms)[keep]
    Fr, pr, R2 = permanova_k(PCk, cmsk, rng, nperm)
    dCw = loo_dist(PCk, cmsk)
    Fw, pw, Dw = permdisp_dist(dA, dCw, rng, nperm)
    f = (Dt - Dw) / (Dt - 1) if Dt > 1 else np.nan
    fr = []
    for _ in range(nperm):
        Dwr = loo_dist(PCk, rng.permutation(cmsk)).mean() / dA.mean()
        fr.append((Dt - Dwr) / (Dt - 1) if Dt > 1 else np.nan)
    f95 = float(np.nanpercentile(fr, 95)) if Dt > 1 else np.nan
    return dict(D_total=Dt, p_total=pt, R2_CMS=R2, p_CMS=pr, D_within=Dw, p_within=pw, f=f, f_rand95=f95,
                n_excluded_small=int((~keep).sum()))


def call(r):
    if not (r["D_total"] > 1 and r["p_total"] < 0.05):
        return "not informative"
    within_sig = r["D_within"] > 1 and r["p_within"] < 0.05
    f, f95 = r["f"], r["f_rand95"]
    if r["p_CMS"] < 0.05 and f >= 0.5 and f > f95 and not within_sig:
        return "subtype diversity"
    if within_sig and f < 0.5:
        return "within-subtype heterogeneity"
    if within_sig and f >= 0.5 and f > f95:
        return "both"
    return "inconclusive"


# ================================================================ section 6 validation
SIZES = (30, 20, 12, 8)


def tetra(edge):
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    v *= edge / np.linalg.norm(v[0] - v[1])
    out = np.zeros((4, 10))
    out[:, :3] = v
    return out


def sim_rep(task):
    scen, rep = task
    rng = rng_for(f"a9sim__{scen}__{rep:03d}")
    sizes = (4, 5, 6, 55) if scen == "V4" else SIZES
    PA = rng.standard_normal((50, 10))
    lab = np.repeat(["CMS1", "CMS2", "CMS3", "CMS4"], sizes)
    if scen == "V1":
        cen = tetra(4.0)                       # Amendment 1b: edge 4 SD
        PC = np.vstack([cen[i] + rng.standard_normal((n, 10)) for i, n in enumerate(sizes)])
    elif scen == "V2":
        PC = 1.4 * rng.standard_normal((70, 10))
        lab = rng.permutation(lab)
    else:
        PC = rng.standard_normal((sum(sizes), 10))
        lab = rng.permutation(lab)
    r = decompose(PA, PC, lab, rng)
    r.update(scenario=scen, rep=rep, call=call(r))
    return r


def bf_null(task):
    cohort, nA, nC, noise = task
    rng = rng_for(f"a9sim__bf__{cohort}__{noise}")
    G = 2000
    if noise == "gaussian":
        A, C = rng.standard_normal((G, nA)), rng.standard_normal((G, nC))
    else:
        A, C = rng.standard_t(4, (G, nA)), rng.standard_t(4, (G, nC))
    p = np.array([levene(A[i], C[i], center="median").pvalue for i in range(G)])
    return dict(cohort=cohort, nA=nA, nC=nC, noise=noise, rejection=float((p < 0.05).mean()))


def renv_path():
    if not RENV:
        return os.environ["PATH"]
    E = RENV
    return ";".join([E, rf"{E}\Library\mingw-w64\bin", rf"{E}\Library\usr\bin", rf"{E}\Library\bin", rf"{E}\Scripts",
                     rf"{E}\lib\R\bin\x64", os.environ["PATH"]])


def run_r(*args):
    env = dict(os.environ, PATH=renv_path())
    r = subprocess.run([RSCRIPT, os.path.join(ROOT, "scripts", "addendum9_cms.R"), *args],
                       capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return r.stdout.strip().splitlines()[-1]


def validate():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    tasks = [(s, r) for s in ("V1", "V2", "V3", "V4") for r in range(200)]
    bf_tasks = [(c, a, b, n) for c, a, b in (("GSE41258", 48, 181), ("GSE20916", 45, 36), ("GSE37364", 29, 27),
                                             ("GSE117606", 59, 71)) for n in ("gaussian", "t4")]
    with Pool(WORKERS) as pool:
        sim = pd.DataFrame(pool.map(sim_rep, tasks, chunksize=5))
        bf = pd.DataFrame(pool.map(bf_null, bf_tasks))
    print(f"[validate] simulations done [{(time.time() - t0) / 60:.1f} min]", flush=True)
    ex_path = os.path.join(OUT, "validation_cmscaller_example.csv")
    rver = run_r("example", ex_path)
    ex = pd.read_csv(ex_path)
    lab = ex[ex.gold.notna()]
    agree = float((lab.nearest == lab.gold).mean())
    conf = lab[lab.prediction.notna()]
    agree_conf = float((conf.prediction == conf.gold).mean()) if len(conf) else np.nan
    sim.to_csv(os.path.join(OUT, "validation_dispersion.csv"), index=False)
    bf.to_csv(os.path.join(OUT, "validation_brown_forsythe.csv"), index=False)

    se200 = np.sqrt(.05 * .95 / 200)
    t200 = 0.05 + 2.5 * se200
    se2k = np.sqrt(.05 * .95 / 2000)
    t2k = 0.05 + 2.5 * se2k
    g = {s: sim[sim.scenario == s] for s in ("V1", "V2", "V3", "V4")}
    sig_tot = lambda d: float(((d.D_total > 1) & (d.p_total < .05)).mean())
    sig_w = lambda d: float(((d.D_within > 1) & (d.p_within < .05)).mean())
    crit = [
        ("V1 subtype diversity", "D_total > 1 significant >= 0.80", sig_tot(g["V1"]), sig_tot(g["V1"]) >= .8),
        ("V1 subtype diversity", "call 'subtype diversity' >= 0.80", float((g["V1"].call == "subtype diversity").mean()),
         float((g["V1"].call == "subtype diversity").mean()) >= .8),
        ("V2 within-subtype", "D_within > 1 significant >= 0.80", sig_w(g["V2"]), sig_w(g["V2"]) >= .8),
        ("V2 within-subtype", "f <= f_rand95 >= 0.90", float((g["V2"].f <= g["V2"].f_rand95).mean()),
         float((g["V2"].f <= g["V2"].f_rand95).mean()) >= .9),
        ("V2 within-subtype", "call 'within-subtype heterogeneity' >= 0.80",
         float((g["V2"].call == "within-subtype heterogeneity").mean()),
         float((g["V2"].call == "within-subtype heterogeneity").mean()) >= .8),
        ("V3 null", f"D_within PERMDISP rejection <= {t200:.3f}", float((g["V3"].p_within < .05).mean()),
         float((g["V3"].p_within < .05).mean()) <= t200),
        ("V3 null", f"D_total PERMDISP rejection <= {t200:.3f}", float((g["V3"].p_total < .05).mean()),
         float((g["V3"].p_total < .05).mean()) <= t200),
        ("V4 LOO bias (groups 4/5/6/55)", f"D_within PERMDISP rejection <= {t200:.3f}", float((g["V4"].p_within < .05).mean()),
         float((g["V4"].p_within < .05).mean()) <= t200),
    ]
    for _, r in bf.iterrows():
        crit.append((f"Brown-Forsythe null, {r.cohort} n {r.nA}/{r.nC}, {r.noise}", f"rejection <= {t2k:.3f}",
                     r.rejection, r.rejection <= t2k))
    crit.append(("CMScaller example (TCGA, 63 gold-standard labels)", "accepted at 0.889 (Amendment 1c)", agree,
                 abs(agree - 56 / 63) < 1e-9 or agree >= 56 / 63))
    pd.crosstab(lab.gold, lab.nearest).to_csv(os.path.join(OUT, "validation_cmscaller_confusion.csv"))
    cdf = pd.DataFrame(crit, columns=["check", "requirement", "value", "pass"])
    cdf.to_csv(os.path.join(OUT, "validation_criteria.csv"), index=False)
    ok = bool(cdf["pass"].all())
    L = ["ADDENDUM 9 - Section 6 validation RERUN (addendum 9 5a9b03a; D16 9590f65; Amendment 1 12dfc3c)", "=" * 70,
         "No real data used. Dispersion: 10 dims, adenomas 50, carcinomas 70 (subtypes 30/20/12/8; V4 4/5/6/55),",
         "200 replicates per scenario, 999 permutations, bias-adjusted distances sqrt(n/(n-1))*||x-c||.",
         "V1 subtype centres: regular tetrahedron, edge 4 SD. Brown-Forsythe: 2,000 null genes per cohort size and noise type.",
         f"CMScaller example: {rver}; RNAseq = TRUE, nPerm 1000, seed 42.",
         f"Seed {SEED}. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("numpy", "scipy")), ""]
    L += [f"  {'PASS' if r['pass'] else 'FAIL'}  {r.check:<58} {r.requirement:<42} {r.value:.3f}" for _, r in cdf.iterrows()]
    L += ["", "Descriptive:",
          f"  CMScaller example: confident calls {len(conf)}/{len(lab)} labelled; agreement among confident {agree_conf:.3f}",
          "  CMScaller confusion (rows gold, columns nearest template):",
          *["    " + line for line in pd.crosstab(lab.gold, lab.nearest).to_string().splitlines()],
          "  Calls per scenario: " + "; ".join(f"{s}: " + ", ".join(f"{k} {v:.2f}" for k, v in d.call.value_counts(normalize=True).items())
                                              for s, d in g.items()),
          "  Median D_total / D_within / f / f_rand95: " + "; ".join(
              f"{s} {d.D_total.median():.2f}/{d.D_within.median():.2f}/{d.f.median():.2f}/{d.f_rand95.median():.2f}" for s, d in g.items()),
          "", f"OVERALL: {'ALL CRITERIA PASSED - real-data analysis may proceed' if ok else 'FAILED - do not analyse real data'}",
          f"Run time {(time.time() - t0) / 60:.1f} min."]
    open(os.path.join(OUT, "validation_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ================================================================ real data (addendum 9 + Amendment 1)
from scipy.stats import median_abs_deviation, hypergeom  # noqa: E402

EMB = os.path.join(ROOT, "results", "step03", "embeddings")
PROC = os.path.join(ROOT, "data", "processed")
CACHE = os.path.join(OUT, "cache")
COH9 = {"GSE41258": ("main", 48, 181), "GSE20916": ("bulk", 45, 36), "GSE37364": ("bulk", 29, 27),
        "GSE117606": ("main", 59, 71)}
GUINNEY = {"CMS1": 0.14, "CMS2": 0.37, "CMS3": 0.13, "CMS4": 0.23, "unclassified": 0.13}


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj, run = np.empty_like(p), 0.0
    for r, i in enumerate(o):
        run = max(run, (len(p) - r) * p[i])
        adj[i] = min(1, run)
    return adj


def load9(g):
    st, nA, nC = COH9[g]
    emb = pd.read_csv(os.path.join(EMB, f"{g}_{st}_adjusted_g2000_pc10.csv"), index_col=0)
    cnt = emb.stage.value_counts()
    if cnt.get("Adenoma") != nA or cnt.get("Carcinoma") != nC:
        raise SystemExit(f"Stage counts changed: {g} {cnt.to_dict()}")
    adj = pd.read_csv(os.path.join(PROC, f"{g}_{st}_expr_genes_adjusted.csv"), index_col=0)[emb.index]
    raw = pd.read_csv(os.path.join(PROC, f"{g}_{st}_expr_genes.csv"), index_col=0)
    missing = set(emb.index) - set(raw.columns)
    if missing:
        raise SystemExit(f"{g}: {len(missing)} samples missing from unadjusted expression")
    return emb, adj, raw[emb.index]


def cms_call(g, X, tag):
    os.makedirs(CACHE, exist_ok=True)
    X = X.dropna(axis=0, how="any")
    X = X[~X.index.duplicated()]
    inp, outp = os.path.join(CACHE, f"cms_in_{g}_{tag}.csv"), os.path.join(CACHE, f"cms_out_{g}_{tag}.csv")
    X.to_csv(inp)
    run_r("classify", inp, outp)
    lab = pd.read_csv(outp).set_index("sample")
    cov = pd.read_csv(outp.replace(".csv", "_coverage.csv"))
    return lab, cov


def variability(g, emb, adj):
    st = emb.stage.values
    X = adj.dropna(axis=0, how="any")
    m = X.mean(1)
    X = X[m > m.quantile(0.25)]
    A, C = X.loc[:, st == "Adenoma"].values, X.loc[:, st == "Carcinoma"].values
    p = np.array([levene(C[i], A[i], center="median").pvalue for i in range(len(X))])
    madA, madC = median_abs_deviation(A, axis=1), median_abs_deviation(C, axis=1)
    df = pd.DataFrame(dict(V=madC / np.where(madA > 0, madA, np.nan), SD_ratio=C.std(1, ddof=1) / A.std(1, ddof=1),
                           mean_diff=C.mean(1) - A.mean(1), p=p), index=X.index)
    df["q"] = bh(df.p.fillna(1).values)
    return df


def decompose2(PA, PC, cms, rng, within_mask=None, nperm=NPERM):
    """D_total from all carcinomas; R2/D_within/f from carcinomas in within_mask (e.g. confident calls)."""
    r = decompose(PA, PC if within_mask is None else PC, cms, rng, nperm) if within_mask is None else None
    if within_mask is not None:
        dA = loo_dist(PA, np.zeros(len(PA)))
        Ft, pt, Dt = permdisp_dist(dA, loo_dist(PC, np.zeros(len(PC))), rng, nperm)
        r = decompose(PA, PC[within_mask], np.asarray(cms)[within_mask], rng, nperm)
        r.update(D_total=Dt, p_total=pt)
        r["f"] = (Dt - r["D_within"]) / (Dt - 1) if Dt > 1 else np.nan
        cm = np.asarray(cms)[within_mask]
        keep = pd.Series(cm).map(pd.Series(cm).value_counts()).values >= 3
        PCk, cmk = PC[within_mask][keep], cm[keep]
        fr = [(Dt - loo_dist(PCk, rng.permutation(cmk)).mean() / dA.mean()) / (Dt - 1) for _ in range(nperm)] if Dt > 1 else [np.nan]
        r["f_rand95"] = float(np.nanpercentile(fr, 95)) if Dt > 1 else np.nan
    return r


def per_subtype(PA, PC, cms, rng):
    dA = loo_dist(PA, np.zeros(len(PA)))
    rows = []
    for k, n in pd.Series(cms).value_counts().items():
        if n >= 10:
            m = np.asarray(cms) == k
            F, p, D = permdisp_dist(dA, loo_dist(PC[m], np.zeros(m.sum())), rng)
            rows.append(dict(CMS=k, n=int(n), D=D, p=p))
    if rows:
        ph = holm([r["p"] for r in rows])
        for r, h in zip(rows, ph):
            r["p_holm"] = h
    return rows


def overall(calls):
    inf = [c for c in calls.values() if c != "not informative"]
    n = lambda x: sum(c == x for c in inf)
    if len(inf) < 2:
        return "INCONCLUSIVE (fewer than 2 informative cohorts)"
    if n("subtype diversity") >= 2 and n("within-subtype heterogeneity") == 0:
        return "HETEROGENEITY REFLECTS SUBTYPE DIVERSITY"
    if n("within-subtype heterogeneity") >= 2 and n("subtype diversity") == 0:
        return "HETEROGENEITY WITHIN SUBTYPES"
    if n("both") >= 2 or (n("subtype diversity") >= 1 and n("within-subtype heterogeneity") >= 1 and n("both") >= 1):
        return "BOTH CONTRIBUTE"
    return "INCONCLUSIVE"


def kappa(a, b):
    a, b = pd.Series(a), pd.Series(b)
    cats = sorted(set(a) | set(b))
    po = float((a.values == b.values).mean())
    pe = sum((a == c).mean() * (b == c).mean() for c in cats)
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


def run_real():
    import importlib.util
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    data, labs, covs, var, dec, sub = {}, {}, {}, {}, [], []
    for g in COH9:
        emb, adj, raw = load9(g)
        data[g] = (emb, adj)
        car = emb.index[emb.stage == "Carcinoma"]
        lab_u, cov = cms_call(g, raw[car], "unadjusted")
        lab_a, _ = cms_call(g, adj[car], "adjusted")
        lab = pd.DataFrame(dict(cohort=g, nearest=lab_u.nearest, prediction=lab_u.prediction, FDR=lab_u.FDR,
                                nearest_adjusted=lab_a.nearest.reindex(lab_u.index)))
        lab["confident"] = lab.prediction.notna()
        lab["coverage"] = float(cov.overall.iloc[0])
        labs[g], covs[g] = lab.loc[car], cov
        print(f"[a9] {g}: CMS done (coverage {100 * cov.overall.iloc[0]:.1f}%) [{(time.time() - t0) / 60:.1f} min]", flush=True)
        var[g] = variability(g, emb, adj)
        PCs = [c for c in emb.columns if c.startswith("PC")]
        PA = emb.loc[emb.stage == "Adenoma", PCs].values
        PC = emb.loc[car, PCs].values
        L = labs[g]
        for vers, cms, mask in (("primary", L.nearest.values, None), ("confident-only", L.nearest.values, L.confident.values),
                                ("adjusted-CMS", L.nearest_adjusted.values, None)):
            r = decompose2(PA, PC, cms, rng_for(f"{g}__decomp__{vers}"), mask)
            r.update(cohort=g, labels=vers, call=call(r), n_carcinoma_within=int(len(PC) if mask is None else mask.sum()))
            dec.append(r)
        for r in per_subtype(PA, PC, L.nearest.values, rng_for(f"{g}__persubtype")):
            sub.append(dict(cohort=g, **r))
        print(f"[a9] {g}: variability + decomposition done [{(time.time() - t0) / 60:.1f} min]", flush=True)
    dec, sub = pd.DataFrame(dec), pd.DataFrame(sub)
    labs_all = pd.concat(labs.values())
    labs_all.to_csv(os.path.join(OUT, "cms_labels.csv"), index_label="sample")
    for g in COH9:
        var[g].to_csv(os.path.join(OUT, f"variability_{g}.csv"), index_label="gene")
    dec.to_csv(os.path.join(OUT, "dispersion_decomposition.csv"), index=False)
    sub.to_csv(os.path.join(OUT, "dispersion_per_subtype.csv"), index=False)

    # CMS proportions
    prop = []
    for g in COH9:
        L = labs[g]
        n = len(L)
        pp = L.nearest.value_counts(normalize=True)
        cc = L.prediction.fillna("unclassified").value_counts(normalize=True)
        ent = float(-(pp[pp > 0] * np.log(pp[pp > 0])).sum())
        row = dict(cohort=g, n=n, coverage=float(covs[g].overall.iloc[0]), low_coverage=bool(covs[g].overall.iloc[0] < 0.8),
                   entropy=ent, kappa_unadj_vs_adj=kappa(L.nearest, L.nearest_adjusted))
        for k in ("CMS1", "CMS2", "CMS3", "CMS4"):
            row[f"primary_{k}"] = float(pp.get(k, 0))
            row[f"confident_{k}"] = float(cc.get(k, 0))
        row["confident_unclassified"] = float(cc.get("unclassified", 0))
        for _, c in covs[g].iterrows():
            row[f"coverage_{c.template}"] = float(c.coverage)
        prop.append(row)
    prop = pd.DataFrame(prop)
    prop.to_csv(os.path.join(OUT, "cms_proportions.csv"), index=False)

    # replication of variability
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
            emb, adj = data[o]
            car = labs[o].index
            y = adj.loc[gene, car].values
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
    uniq.to_csv(os.path.join(OUT, "variability_replicated.csv"), index=False)

    # template overlap and enrichment
    bg = [x for x in set().union(*tested.values()) if sum(x in tested[o] for o in COH9) >= 3]
    tpl = set(pd.read_csv(os.path.join(CACHE, "templates.csv")).symbol) if os.path.exists(os.path.join(CACHE, "templates.csv")) else set()
    spec = importlib.util.spec_from_file_location("a8", os.path.join(ROOT, "scripts", "addendum8_gradual_switch.py"))
    a8 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a8)
    lists = {}
    if len(uniq):
        for d, v in uniq.groupby("direction"):
            lists[d] = list(v.gene)
    enr, enr_msg = a8.enrichment(lists, bg)
    enr.to_csv(os.path.join(OUT, "enrichment_variability.csv"), index=False)
    tpl_line = ""
    if tpl and len(uniq):
        mv = set(uniq[uniq.direction == "more variable in carcinoma"].gene) & set(bg)
        T_ = tpl & set(bg)
        k = len(mv & T_)
        pt = hypergeom.sf(k - 1, len(bg), len(T_), len(mv)) if k else 1.0
        tpl_line = f"overlap of more-variable genes with CMS templates: {k} of {len(mv)} (template genes in background {len(T_)}; hypergeometric p {pt:.3g})"

    # section 5
    calls = {lv: {r.cohort: r.call for _, r in dec[dec.labels == lv].iterrows()} for lv in dec.labels.unique()}
    concl = {lv: overall(c) for lv, c in calls.items()}

    S = ["ADDENDUM 9 - Carcinoma vs adenoma heterogeneity (5a9b03a; D16 9590f65; Amendment 1 12dfc3c)", "=" * 78,
         "Stage-blind adjusted data (dispersion: frozen 10-PC embeddings; variability: adjusted expression). CMS: CMScaller",
         "2.0.1 on unadjusted expression of carcinomas, centred per gene within cohort (D16a); nearest template = primary.",
         "Bias-adjusted PERMDISP distances sqrt(n/(n-1))*||x-c||; 999 permutations; seed 42 per unit.",
         "Samples (A/C): " + "; ".join(f"{g} {COH9[g][1]}/{COH9[g][2]}" for g in COH9), "",
         "1. CMS TEMPLATE COVERAGE (flag < 80%)"]
    for _, r in prop.iterrows():
        S.append(f"  {r.cohort:<9} {100 * r.coverage:.1f}%{'  ** LOW COVERAGE **' if r.low_coverage else ''}  (per template: "
                 + ", ".join(f"{k} {100 * r[f'coverage_{k}']:.0f}%" for k in ("CMS1", "CMS2", "CMS3", "CMS4")) + ")")
    S += ["", "2. CMS PROPORTIONS (primary nearest template | confident-only, FDR < 0.05)",
          "  Guinney 2015: " + ", ".join(f"{k} {100 * v:.0f}%" for k, v in GUINNEY.items())]
    for _, r in prop.iterrows():
        S.append(f"  {r.cohort:<9} n {r.n}: primary " + ", ".join(f"{k} {100 * r[f'primary_{k}']:.0f}%" for k in ("CMS1", "CMS2", "CMS3", "CMS4"))
                 + " | confident " + ", ".join(f"{k} {100 * r[f'confident_{k}']:.0f}%" for k in ("CMS1", "CMS2", "CMS3", "CMS4"))
                 + f", unclassified {100 * r.confident_unclassified:.0f}%; entropy {r.entropy:.2f}; kappa unadj vs adj {r.kappa_unadj_vs_adj:.2f}"
                 + ("  [low coverage]" if r.low_coverage else ""))
    S += ["", "3-4. DISPERSION DECOMPOSITION (stage-blind adjusted 10-PC space)"]
    for lv in ("primary", "confident-only", "adjusted-CMS"):
        S.append(f"  -- labels: {lv} --")
        for _, r in dec[dec.labels == lv].iterrows():
            S.append(f"  {r.cohort:<9} D_total {r.D_total:.2f} (p {r.p_total:.3f}); R2_CMS {r.R2_CMS:.2f} (p {r.p_CMS:.3f}); "
                     f"D_within {r.D_within:.2f} (p {r.p_within:.3f}); f {r.f:.2f} vs f_rand95 {r.f_rand95:.2f}; "
                     f"n within {r.n_carcinoma_within} (excluded small groups {r.n_excluded_small}) -> {r.call.upper()}"
                     + ("  [low coverage]" if prop.set_index("cohort").loc[r.cohort, "low_coverage"] else ""))
    S.append("  Per-subtype PERMDISP (CMS groups >= 10 vs all adenomas; Holm within cohort; primary labels):")
    for _, r in sub.iterrows():
        S.append(f"    {r.cohort:<9} {r.CMS} n {r.n}: D {r.D:.2f}, p {r.p:.3f}, Holm {r.p_holm:.3f}")
    S += ["", "5. SECTION 5 OVERALL CONCLUSION"]
    for lv, c in concl.items():
        S.append(f"  {lv:<15} -> {c}   (calls: " + ", ".join(f"{k} {v}" for k, v in calls[lv].items()) + ")")
    S.append(f"  ROBUST (confident-only gives the same conclusion): {concl.get('primary') == concl.get('confident-only')}")
    nmv = int((uniq.direction == "more variable in carcinoma").sum()) if len(uniq) else 0
    nlv = int((uniq.direction == "less variable in carcinoma").sum()) if len(uniq) else 0
    S += ["", "6. DIFFERENTIAL VARIABILITY (Brown-Forsythe; discovery BH q < 0.05; replication p < 0.05 same direction in >= 2 other cohorts)",
          "  Genes tested: " + ", ".join(f"{g} {len(var[g])}" for g in COH9)
          + "; discoveries (q < 0.05): " + ", ".join(f"{g} more {int(((var[g].q < .05) & (var[g].V > 1)).sum())} / less {int(((var[g].q < .05) & (var[g].V < 1)).sum())}" for g in COH9),
          f"  Replicated MORE variable in carcinoma: {nmv}; LESS variable: {nlv}"]
    if nmv:
        mv = uniq[uniq.direction == "more variable in carcinoma"]
        S.append(f"  Without mean shift (median |mean diff| < 0.5 log2): {int(mv.no_mean_shift.sum())} of {nmv}; "
                 f"median CMS R2 among carcinomas {mv.median_CMS_R2.median():.2f}")
        S.append("  Top 25 by median variability ratio: " + ", ".join(f"{r.gene} ({r.median_V:.2f}{'' if r.no_mean_shift else '*'})" for _, r in mv.head(25).iterrows()))
        S.append("  (* = mean shift >= 0.5 log2)")
    if nlv:
        lv_ = uniq[uniq.direction == "less variable in carcinoma"]
        S.append("  Less variable (top 15): " + ", ".join(f"{r.gene} ({r.median_V:.2f})" for _, r in lv_.head(15).iterrows()))
    if tpl_line:
        S.append("  " + tpl_line)
    S += [f"  Enrichment ({enr_msg}; background {len(bg)} genes testable in >= 3 cohorts):"]
    if len(enr):
        for name, v in enr.groupby("list"):
            top = v[v.q < 0.05].head(12)
            S.append(f"    {name}: " + ("; ".join(f"{r.gene_set} (k={r.k}/{r.list_size}, q={r.q:.2g})" for _, r in top.iterrows())
                                       if len(top) else "no gene set with q < 0.05"))
    else:
        S.append("    no list with >= 5 genes")
    S.append(f"Run time {(time.time() - t0) / 60:.1f} min.")
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(S) + "\n")
    print("\n".join(S))
    figures9(data, labs, dec, uniq)


def figures9(data, labs, dec, uniq):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"CMS1": "#e8a33d", "CMS2": "#2a78d6", "CMS3": "#d64fa0", "CMS4": "#2a9d55"}
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    for ax, g in zip(axes, COH9):
        emb, _ = data[g]
        a = emb[emb.stage == "Adenoma"]
        ax.scatter(a.PC1, a.PC2, s=8, c="0.7", label="adenoma")
        c = emb.loc[labs[g].index]
        for k, cc in col.items():
            m = labs[g].nearest.values == k
            ax.scatter(c.PC1[m], c.PC2[m], s=10, c=cc, label=k)
        ax.set_title(g, fontsize=9); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    axes[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "cms_pc12.png"), dpi=150); plt.close(fig)
    d = dec[dec.labels == "primary"]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    x = np.arange(len(d))
    ax.bar(x - 0.2, d.D_total, 0.4, label="D_total"); ax.bar(x + 0.2, d.D_within, 0.4, label="D_within")
    ax.axhline(1, color="0.4", lw=1); ax.set_xticks(x); ax.set_xticklabels(d.cohort, fontsize=8)
    ax.set_ylabel("carcinoma / adenoma dispersion"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "dispersion_total_vs_within.png"), dpi=150); plt.close(fig)
    if len(uniq):
        mv = uniq[uniq.direction == "more variable in carcinoma"].head(30)
        if len(mv):
            fig, ax = plt.subplots(figsize=(8, 0.25 * len(mv) + 1))
            for j, g in enumerate(COH9):
                ax.scatter(np.log2(mv[f"V_{g}"]), np.arange(len(mv)), s=14, label=g)
            ax.set_yticks(np.arange(len(mv))); ax.set_yticklabels(mv.gene, fontsize=7); ax.invert_yaxis()
            ax.axvline(0, color="0.4", lw=1); ax.set_xlabel("log2 variability ratio (carcinoma / adenoma, MAD)"); ax.legend(fontsize=7)
            fig.tight_layout(); fig.savefig(os.path.join(OUT, "variability_replicated.png"), dpi=150); plt.close(fig)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "validate":
        validate()
    elif cmd == "run":
        run_real()
    else:
        raise SystemExit("Usage: addendum9_heterogeneity.py validate | run")
