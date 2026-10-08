"""
Addendum 8 (commit 8e8df5b) - gradual vs switch-like change along normal -> adenoma -> carcinoma.
Implementation deviation D13 (248d295): sigmoid fitted by profile grid over (k > 0, x0) with closed-form
plateaus, then local refinement (curve_fit); same procedure for genes and bootstrap simulations.

Part (a): progression score (adjusted PC1, scaled N median 0 / C median 1), Hartigan dip test with Monte Carlo
uniform calibration, between-median counts, AUC, valley ratio with stratified bootstrap CI.
Part (b): cross-fitted half scores; per gene spline association (BH), linear vs 4-parameter sigmoid (AIC),
classification (switch-like / step-like unresolvable / gradual / intermediate / non-monotonic), parametric
bootstrap of dAIC on a grid of standardised slopes (BH), replication across cohorts.

Usage:
  python scripts/addendum8_gradual_switch.py validate    - section 6 simulation only (no real data)
"""
import os
import sys
import time
import zlib
from importlib import metadata
from multiprocessing import Pool

import numpy as np
import pandas as pd
from patsy import build_design_matrices, dmatrix
from scipy.optimize import curve_fit
from scipy.special import expit
from scipy.stats import f as f_dist, mannwhitneyu, spearmanr
from sklearn.decomposition import PCA

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "addendum8")
SEED, WORKERS = 42, 6
LN81 = 2 * np.log(9)                       # W = 2 ln 9 / k  (10% -> 90% width)
RES_HALF, RES_FRAC = 0.125, 0.10          # Amendment 1 (0275b7a): resolvable if >= 10% (and >= 5) within x0 +/- 0.125
K_GRID = np.unique(np.r_[np.geomspace(1, 200, 40), [2, 5, 10, 20, 50]])


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ---------------------------------------------------------------- progression scores
def pc1(X):
    """X genes x samples (array). Step 3 preprocessing without gene filter: row z-score, PCA (random_state 0)."""
    X = np.asarray(X, float)
    Z = (X - X.mean(1, keepdims=True)) / (X.std(1, ddof=1, keepdims=True) + 1e-9)
    return PCA(n_components=min(10, Z.shape[1] - 1), random_state=0).fit_transform(Z.T)[:, 0]


def scale_score(raw, stage):
    mN, mC = np.median(raw[stage == "Normal"]), np.median(raw[stage == "Carcinoma"])
    s = (raw - mN) / (mC - mN)
    astar = float(np.median(s[stage == "Adenoma"]))
    return s, astar, bool(0 < astar < 1)


# ---------------------------------------------------------------- part (a)
_DIPNULL = {}


def dip_null(n, nsim=9999):
    if n not in _DIPNULL:
        import diptest
        r = rng_for(f"dipnull__n{n}")
        _DIPNULL[n] = np.sort([diptest.dipstat(r.uniform(size=n)) for _ in range(nsim)])
    return _DIPNULL[n]


def dip_test(x):
    import diptest
    d = float(diptest.dipstat(np.asarray(x, float)))
    null = dip_null(len(x))
    return d, float((1 + np.sum(null >= d)) / (1 + len(null)))


def valley_ratio(sx, sy):
    mx, my = np.median(sx), np.median(sy)
    w = my - mx
    if w <= 0:
        return np.nan
    allv = np.r_[sx, sy]
    mid = np.sum((allv >= mx + w / 3) & (allv <= mx + 2 * w / 3))
    rx = np.sum((allv >= mx - w / 6) & (allv <= mx + w / 6))
    ry = np.sum((allv >= my - w / 6) & (allv <= my + w / 6))
    den = (rx + ry) / 2
    return mid / den if den > 0 else np.nan


def continuity(sx, sy, rng, nboot=2000):
    vr = valley_ratio(sx, sy)
    boots = np.array([valley_ratio(rng.choice(sx, len(sx)), rng.choice(sy, len(sy))) for _ in range(nboot)])
    boots = boots[np.isfinite(boots)]
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) else (np.nan, np.nan))
    auc = mannwhitneyu(sy, sx).statistic / (len(sx) * len(sy))
    mx, my = np.median(sx), np.median(sy)
    return dict(VR=vr, VR_lo=lo, VR_hi=hi, AUC=auc, n_x_between=int(np.sum((sx > mx) & (sx < my))),
                n_y_between=int(np.sum((sy > mx) & (sy < my))), n_x=len(sx), n_y=len(sy))


# ---------------------------------------------------------------- part (b) fits
def sig(x, a, b, k, x0):
    return a + (b - a) * expit(k * (x - x0))


class SigFitter:
    """4-parameter sigmoid fitted by profile grid over (k > 0, x0) + local refinement (D13)."""

    def __init__(self, x, astar):
        self.x = np.asarray(x, float)
        x0s = np.unique(np.r_[np.linspace(self.x.min(), self.x.max(), 60),
                              np.clip([0, astar, 0.5, 1], self.x.min(), self.x.max())])
        kk, xx = np.meshgrid(K_GRID, x0s, indexing="ij")
        kk, xx = kk.ravel(), xx.ravel()
        S = expit(kk[:, None] * (self.x[None, :] - xx[:, None]))
        Sc = S - S.mean(1, keepdims=True)
        nrm = np.sqrt((Sc ** 2).sum(1))
        keep = nrm > 1e-6
        self.k, self.x0, self.Sm, self.nrm = kk[keep], xx[keep], S.mean(1)[keep], nrm[keep]
        self.Q = Sc[keep] / nrm[keep, None]
        self.lo, self.hi = self.x.min(), self.x.max()

    def fit(self, Y, refine=True):
        Y = np.asarray(Y, float)
        if Y.ndim == 1:
            Y = Y[:, None]
        yb = Y.mean(0)
        Yc = Y - yb
        P = self.Q @ Yc
        j = np.argmax(P ** 2, axis=0)
        cols = np.arange(Y.shape[1])
        pj = P[j, cols]
        rss = (Yc ** 2).sum(0) - pj ** 2
        beta = pj / self.nrm[j]
        a = yb - beta * self.Sm[j]
        b = a + beta
        k, x0 = self.k[j].copy(), self.x0[j].copy()
        if refine:
            for c in cols:
                try:
                    p, _ = curve_fit(sig, self.x, Y[:, c], p0=[a[c], b[c], k[c], x0[c]],
                                     bounds=([-np.inf, -np.inf, 1e-3, self.lo], [np.inf, np.inf, 200, self.hi]),
                                     maxfev=2000)
                    r = float(((Y[:, c] - sig(self.x, *p)) ** 2).sum())
                    if r < rss[c]:
                        rss[c], a[c], b[c], k[c], x0[c] = r, p[0], p[1], p[2], p[3]
                except (RuntimeError, ValueError):
                    pass
        return np.maximum(rss, 1e-12), a, b, k, x0


def linfit(x, Y):
    X = np.c_[np.ones_like(x), x]
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = ((Y - X @ coef) ** 2).sum(0)
    return np.maximum(rss, 1e-12), coef[1]


def aic(rss, n, p):
    return n * np.log(rss / n) + 2 * p


def daic_null(x, fitter, r_grid, nsim, rng, E=None):
    """Null dAIC distributions (linear trend, Gaussian noise) for each standardised slope in r_grid."""
    n = len(x)
    out = []
    for r in r_grid:
        e = rng.standard_normal((n, nsim)) if E is None else E
        Y = r * x[:, None] + e
        rl, _ = linfit(x, Y)
        rs = fitter.fit(Y)[0]
        out.append(np.sort(aic(rl, n, 2) - aic(rs, n, 4)))
    return out


def p_from_null(null_sorted, obs):
    return (1 + len(null_sorted) - np.searchsorted(null_sorted, obs, side="left")) / (1 + len(null_sorted))


def grid_p(r, d, r_grid, nulls):
    j = int(np.clip(np.searchsorted(r_grid, r, side="right") - 1, 0, len(r_grid) - 2))
    return max(p_from_null(nulls[j], d), p_from_null(nulls[j + 1], d))


def fit_genes(x, Y, astar, n_all=None):
    """All per-gene fits against score x. Y: samples x genes. Returns DataFrame (one row per gene)."""
    n = len(x)
    n_all = n if n_all is None else n_all
    fitter = SigFitter(x, astar)
    rl, slope = linfit(x, Y)
    rs, a, b, k, x0 = fitter.fit(Y)
    dm = dmatrix("cr(x, df=3, constraints='center')", {"x": x}, return_type="matrix")
    B = np.asarray(dm)
    coef, *_ = np.linalg.lstsq(B, Y, rcond=None)
    rsp = np.maximum(((Y - B @ coef) ** 2).sum(0), 1e-12)
    rss0 = ((Y - Y.mean(0)) ** 2).sum(0)
    F = ((rss0 - rsp) / 3) / (rsp / (n - 4))
    p_assoc = f_dist.sf(F, 3, n - 4)
    g = np.linspace(0, 1, 101)
    Bg = np.asarray(build_design_matrices([dm.design_info], {"x": g})[0])
    curves = Bg @ coef
    dif = np.diff(curves, axis=0)
    sgn = np.sign(np.where(np.abs(dif) < 1e-9, 0, dif))
    change = np.array([len(np.unique(s[s != 0])) > 1 for s in sgn.T])
    ext = np.where(np.abs(curves.max(0) - curves[[0, -1]].max(0)) >= np.abs(curves.min(0) - curves[[0, -1]].min(0)),
                   "peak", "trough")
    turn_x = np.where(ext == "peak", g[curves.argmax(0)], g[curves.argmin(0)])
    A_l, A_s, A_p = aic(rl, n, 2), aic(rs, n, 4), aic(rsp, n, 4)
    W = LN81 / k
    win_fit = np.array([np.sum((x >= c - w / 2) & (x <= c + w / 2)) for c, w in zip(x0, W)])  # record only
    win = np.array([np.sum(np.abs(x - c) <= RES_HALF) for c in x0])   # Amendment 1: fixed window x0 +/- 0.125
    sigma = np.sqrt(rl / (n - 2))
    return pd.DataFrame(dict(p_assoc=p_assoc, slope=slope, r_std=np.abs(slope) / sigma, aic_lin=A_l, aic_sig=A_s,
                             aic_spl=A_p, dAIC=A_l - A_s, a=a, b=b, k=k, x0=x0, W=W, amp=np.abs(b - a),
                             direction=np.sign(b - a), n_window=win, n_window_fitted=win_fit,
                             resolvable=win >= max(5, RES_FRAC * n_all), spline_turn=change, turn_type=ext,
                             turn_x=turn_x,
                             astar=astar)), fitter


def classify(df):
    df = df.copy()
    df["q_assoc"] = bh(df.p_assoc.values)
    assoc = df.q_assoc < 0.05
    nonmono = (df.aic_spl <= np.minimum(df.aic_lin, df.aic_sig) - 10) & df.spline_turn
    c14 = (df.dAIC >= 10) & (df.W <= 0.25) & (df.x0 >= 0) & (df.x0 <= 1) & (df.amp >= 1)
    cls = np.where(~assoc, "not associated",
          np.where(nonmono, "non-monotonic",
          np.where(c14 & df.resolvable, "switch-like",
          np.where(c14 & ~df.resolvable, "step-like unresolvable",
          np.where((df.dAIC < 2) | (df.W >= 1), "gradual", "intermediate")))))
    df["cls"] = cls
    df["interval"] = np.where(df.x0 <= df.astar, "N->A", "A->C")
    return df


# ================================================================ section 6 validation
N_SIM = 150


def sim_cohort(layout, scen, rng):
    if layout == "continuum":
        t = rng.uniform(-0.3, 1.3, N_SIM)
        stage = np.where(t < 0.25, "Normal", np.where(t <= 0.75, "Adenoma", "Carcinoma"))
    else:
        t = np.r_[rng.normal(0, .08, 50), rng.normal(.5, .08, 50), rng.normal(1, .08, 50)]
        stage = np.repeat(["Normal", "Adenoma", "Carcinoma"], 50)
    bg = rng.standard_normal(1000)[:, None] * t[None, :] + rng.standard_normal((1000, N_SIM))
    noise = rng.standard_normal((960, N_SIM))
    lin = 2 * t[None, :] + 0.5 * rng.standard_normal((20, N_SIM))
    if scen == "sigmoid":
        test = 2 * expit((LN81 / 0.1) * (t[None, :] - 0.25)) + 0.5 * rng.standard_normal((20, N_SIM))
    elif scen == "peak":
        test = 2 * np.exp(-(t[None, :] - 0.5) ** 2 / (2 * 0.12 ** 2)) + 0.5 * rng.standard_normal((20, N_SIM))
    else:
        test = 2 * t[None, :] + 0.5 * rng.standard_normal((20, N_SIM))
    X = np.vstack([bg, noise, test, lin])
    gtype = np.r_[["bg"] * 1000, ["noise"] * 960, [scen] * 20, ["linear"] * 20]
    return X, t, stage, gtype


def half_scores(X, stage, rng):
    perm = rng.permutation(X.shape[0])
    h = np.zeros(X.shape[0], int)
    h[perm[X.shape[0] // 2:]] = 1
    sc = [scale_score(pc1(X[h == i]), stage) for i in (0, 1)]
    return h, sc


def pop_rep(task):
    layout, rep = task
    rng = rng_for(f"a8sim__pop__{layout}__{rep:03d}")
    X, t, stage, _ = sim_cohort(layout, "linear", rng)
    s, astar, ordered = scale_score(pc1(X), stage)
    d, p = dip_test(s)
    h, sc = half_scores(X, stage, rng)
    out = dict(layout=layout, rep=rep, ordered=ordered, astar=astar, dip=d, dip_p=p,
               rho_h1=abs(spearmanr(sc[0][0], t)[0]), rho_h2=abs(spearmanr(sc[1][0], t)[0]))
    for nm, (x, y) in (("NA", ("Normal", "Adenoma")), ("AC", ("Adenoma", "Carcinoma"))):
        c = continuity(s[stage == x], s[stage == y], rng)
        out.update({f"{nm}_{k}": v for k, v in c.items()})
    return out


def gene_rep(task, nsim=499, refine_null=True):
    layout, scen, rep = task
    rng = rng_for(f"a8sim__gene__{layout}__{scen}__{rep:03d}")
    X, t, stage, gtype = sim_cohort(layout, scen, rng)
    h, sc = half_scores(X, stage, rng)
    test = np.where(np.isin(gtype, [scen, "linear"]))[0]
    rows = []
    for half in (0, 1):
        genes = test[h[test] != half]           # genes NOT in this half are fitted against its score
        x, astar, ordered = sc[half]
        df, fitter = fit_genes(x, X[genes].T, astar)
        rg = (np.linspace(0, np.sqrt(df.r_std.max()), 21)) ** 2
        nulls = daic_null(x, fitter, rg, nsim, rng)
        df["p_boot"] = [grid_p(r, d, rg, nulls) for r, d in zip(df.r_std, df.dAIC)]
        df["gtype"], df["gene"], df["score_half"], df["ordered"] = gtype[genes], genes, half, ordered
        rows.append(df)
    df = classify(pd.concat(rows, ignore_index=True))
    df["q_boot"] = np.nan
    a = df.q_assoc < 0.05
    df.loc[a, "q_boot"] = bh(df.loc[a, "p_boot"].values)
    df["cand_switch"] = (df.cls == "switch-like") & (df.q_boot < 0.05)
    df["layout"], df["scen"], df["rep"] = layout, scen, rep
    return df


def grid_check():
    """Grid vs direct bootstrap p, full settings (1,999 sims), common random numbers (D13b)."""
    rng = rng_for("a8sim__gridcheck")
    X, t, stage, gtype = sim_cohort("continuum", "sigmoid", rng)
    h, sc = half_scores(X, stage, rng)
    test = np.where(np.isin(gtype, ["sigmoid", "linear"]))[0]
    pick = rng.choice(test, 20, replace=False)
    E = rng.standard_normal((N_SIM, 1999))
    res = []
    for half in (0, 1):
        genes = pick[h[pick] != half]
        if len(genes) == 0:
            continue
        x, astar, _ = sc[half]
        allg = test[h[test] != half]
        df_all, fitter = fit_genes(x, X[allg].T, astar)
        rg = (np.linspace(0, np.sqrt(df_all.r_std.max()), 21)) ** 2
        nulls = daic_null(x, fitter, rg, 1999, rng, E=E)
        df, _ = fit_genes(x, X[genes].T, astar)
        for gname, r, d in zip(genes, df.r_std, df.dAIC):
            pg = grid_p(r, d, rg, nulls)
            pd_ = p_from_null(daic_null(x, fitter, [r], 1999, rng, E=E)[0], d)
            res.append(dict(gene=int(gname), type=gtype[gname], r_std=r, dAIC=d, p_grid=pg, p_direct=pd_,
                            abs_diff=abs(pg - pd_)))
    return pd.DataFrame(res)


def validate():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    dip_null(N_SIM)
    pop_tasks = [(l, r) for l in ("continuum", "clusters") for r in range(200)]
    gene_tasks = [(l, s, r) for l, s in (("continuum", "sigmoid"), ("clusters", "sigmoid"), ("continuum", "peak"))
                  for r in range(50)]
    with Pool(WORKERS, initializer=dip_null, initargs=(N_SIM,)) as pool:
        pop = pd.DataFrame(pool.map(pop_rep, pop_tasks, chunksize=4))
        print(f"[validate] population done [{(time.time() - t0) / 60:.1f} min]", flush=True)
        genes = []
        for i, df in enumerate(pool.imap_unordered(gene_rep, gene_tasks), 1):
            genes.append(df)
            if i % 10 == 0:
                print(f"[validate] gene reps {i}/{len(gene_tasks)} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    genes = pd.concat(genes, ignore_index=True)
    gc = grid_check()
    pop.to_csv(os.path.join(OUT, "validation_population.csv"), index=False)
    genes.to_csv(os.path.join(OUT, "validation_genes.csv"), index=False)
    gc.to_csv(os.path.join(OUT, "validation_gridcheck.csv"), index=False)

    se = np.sqrt(0.05 * 0.95 / 200)
    crit = []
    for l in ("continuum", "clusters"):
        p = pop[pop.layout == l]
        rho = float(np.median(np.r_[p.rho_h1, p.rho_h2]))
        crit.append((f"Half scores recover t ({l})", "median Spearman >= 0.90", rho, rho >= 0.9))
    pc, pk = pop[pop.layout == "continuum"], pop[pop.layout == "clusters"]
    rj = float((pc.dip_p < .05).mean())
    crit.append(("Dip test, continuum (unimodal)", f"rejection <= {0.05 + 2.5 * se:.3f}", rj, rj <= 0.05 + 2.5 * se))
    rj = float((pk.dip_p < .05).mean())
    crit.append(("Dip test, three clusters", "rejection >= 0.80", rj, rj >= 0.8))
    for nm in ("NA", "AC"):
        v = float(np.nanmedian(pc[f"{nm}_VR"]))
        crit.append((f"VR continuum, {nm[0]}->{nm[1]}", "median VR >= 0.5", v, v >= 0.5))
    for nm in ("NA", "AC"):
        v = float(np.nanmedian(pk[f"{nm}_VR_hi"]))
        crit.append((f"VR three clusters, {nm[0]}->{nm[1]}", "median VR upper CI < 0.25", v, v < 0.25))
    lc = genes[(genes.gtype == "linear") & (genes.layout == "continuum")]
    lk = genes[(genes.gtype == "linear") & (genes.layout == "clusters")]
    for nm, d in (("continuum", lc), ("clusters", lk)):
        v = float(d.cand_switch.mean())
        crit.append((f"Linear genes, {nm}: candidate switch rate", "<= 0.05", v, v <= 0.05))
    sc_ = genes[(genes.gtype == "sigmoid") & (genes.layout == "continuum")]
    v = float(sc_.cand_switch.mean())
    crit.append(("Sigmoid genes, continuum: candidate switch power", ">= 0.80", v, v >= 0.8))
    sk = genes[(genes.gtype == "sigmoid") & (genes.layout == "clusters")]
    v = float((sk.cls == "step-like unresolvable").mean())
    v2 = float((sk.cls == "switch-like").mean())
    crit.append(("Sigmoid genes, clusters: step-like unresolvable", ">= 0.80 (switch-like reported)", v, v >= 0.8))
    pkg = genes[genes.gtype == "peak"]
    v = float((pkg.cls == "non-monotonic").mean())
    crit.append(("Peak genes, continuum: non-monotonic", ">= 0.80", v, v >= 0.8))
    v = float(gc.abs_diff.max())
    crit.append(("Grid vs direct bootstrap p (20 genes)", "max |diff| <= 0.02", v, v <= 0.02))
    cdf = pd.DataFrame(crit, columns=["check", "requirement", "value", "pass"])
    cdf.to_csv(os.path.join(OUT, "validation_criteria.csv"), index=False)
    ok = bool(cdf["pass"].all())

    L = ["ADDENDUM 8 - Section 6 validation on simulated data (addendum 8 8e8df5b; D13 248d295;",
         "Amendment 1 0275b7a: fixed-window resolvability, x0 +/- 0.125, >= 10% and >= 5 samples) - RERUN", "=" * 72,
         "No real data used. n = 150 per simulated cohort; 2,000 genes (1,000 progression background, 960 noise,",
         "20 scenario test genes, 20 linear test genes). Population checks: 200 replicates per layout, all-gene PC1.",
         "Gene checks: 50 replicates per scenario, cross-fitted half scores, grid 21 x 499 simulations, BH over 40",
         "test genes. Dip p: Monte Carlo vs uniform (9,999). VR bootstrap 2,000. Grid check: 1,999 sims, CRN.",
         f"Seed {SEED} (SeedSequence per replicate). " + ", ".join(
             f"{p} {metadata.version(p)}" for p in ("numpy", "scipy", "scikit-learn", "patsy", "diptest")), ""]
    L += [f"  {'PASS' if r['pass'] else 'FAIL'}  {r.check:<50} {r.requirement:<32} {r.value:.3f}" for _, r in cdf.iterrows()]
    L += ["", "Descriptive:",
          f"  PC1 ordered (a* in (0,1)): continuum {pc.ordered.mean():.2f}, clusters {pk.ordered.mean():.2f}",
          f"  Median VR: continuum N->A {np.nanmedian(pc.NA_VR):.2f}, A->C {np.nanmedian(pc.AC_VR):.2f}; "
          f"clusters N->A {np.nanmedian(pk.NA_VR):.2f}, A->C {np.nanmedian(pk.AC_VR):.2f}",
          f"  Sigmoid genes in clusters classified switch-like: {v2:.3f}"]
    for (l, s_), d in genes[genes.gtype != "bg"].groupby(["layout", "gtype"]):
        L.append(f"  Classes, {l} / {s_} genes: " + ", ".join(f"{k} {v:.2f}" for k, v in d.cls.value_counts(normalize=True).items()))
    L += ["", f"OVERALL: {'ALL CRITERIA PASSED - real-data analysis may proceed' if ok else 'FAILED - do not analyse real data'}",
          f"Run time {(time.time() - t0) / 60:.1f} min."]
    open(os.path.join(OUT, "validation_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ================================================================ real data (addendum 8 + Amendments 1, 2)
import pickle  # noqa: E402

EMB = os.path.join(ROOT, "results", "step03", "embeddings")
PROC = os.path.join(ROOT, "data", "processed")
CACHE = os.path.join(OUT, "cache")
COH8 = {"GSE41258": ("main", {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39}),
        "GSE20916": ("bulk", {"Normal": 24, "Adenoma": 45, "Carcinoma": 36}),
        "GSE37364": ("bulk", {"Normal": 38, "Adenoma": 29, "Carcinoma": 27}),
        "GSE117606": ("main", {"Normal": 65, "Adenoma": 59, "Carcinoma": 71})}
NSIM_REAL = 1999
TRS = (("N->A", "Normal", "Adenoma", "NA"), ("A->C", "Adenoma", "Carcinoma", "AC"))


def load_cohort(g):
    st, exp = COH8[g]
    tag = f"{g}_{st}_adjusted_g2000_pc10"
    emb = pd.read_csv(os.path.join(EMB, tag + ".csv"), index_col=0)
    if emb.stage.value_counts().to_dict() != exp:
        raise SystemExit(f"Stage counts changed: {g} {emb.stage.value_counts().to_dict()}")
    genes = pd.read_csv(os.path.join(EMB, tag + "_genes.csv")).gene.tolist()
    expr = pd.read_csv(os.path.join(PROC, f"{g}_{st}_expr_genes_adjusted.csv"), index_col=0)[emb.index]
    return emb, genes, expr


def pop_call(c, dip_holm):
    if c["AUC"] < 0.70:
        return "overlapping"
    if c["VR"] >= 0.5 and c["VR_lo"] >= 0.25 and dip_holm >= 0.05:
        return "continuous"
    if dip_holm < 0.05 and c["VR_hi"] < 0.25:
        return "discontinuous"
    return "inconclusive"


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj, run = np.empty_like(p), 0.0
    for r, i in enumerate(o):
        run = max(run, (len(p) - r) * p[i])
        adj[i] = min(1, run)
    return adj


def part_a(s, st, rng):
    d = {}
    pairs = {"all": np.ones(len(s), bool), "NA": np.isin(st, ["Normal", "Adenoma"]), "AC": np.isin(st, ["Adenoma", "Carcinoma"])}
    dips = {k: dip_test(s[m]) for k, m in pairs.items()}
    ph = holm([dips[k][1] for k in pairs])
    for (k, (dv, pv)), h in zip(dips.items(), ph):
        d[f"dip_{k}"], d[f"dip_p_{k}"], d[f"dip_holm_{k}"] = dv, pv, h
    for tr, x, y, nm in TRS:
        c = continuity(s[st == x], s[st == y], rng)
        if c["AUC"] < 0.70:
            c["VR"] = c["VR_lo"] = c["VR_hi"] = np.nan
        d.update({f"{nm}_{k}": v for k, v in c.items()})
        d[f"{nm}_call"] = pop_call(c, d[f"dip_holm_{nm}"])
    return d


def score_fits(x, Y, astar, n_all, rng, nsim=NSIM_REAL):
    df, fitter = fit_genes(x, Y.values, astar, n_all)
    df.index = Y.columns
    rg = (np.linspace(0, np.sqrt(df.r_std.max()), 21)) ** 2
    nulls = daic_null(x, fitter, rg, nsim, rng)
    df["p_boot"] = [grid_p(r, d, rg, nulls) for r, d in zip(df.r_std, df.dAIC)]
    return df, (rg, nulls)


def finish(df):
    df = classify(df)
    df["q_boot"] = np.nan
    a = df.q_assoc < 0.05
    df.loc[a, "q_boot"] = bh(df.loc[a, "p_boot"].values)
    df["cand_switch"] = (df.cls == "switch-like") & (df.q_boot < 0.05)
    return df


def stage1(g):
    path = os.path.join(CACHE, f"{g}.pkl")
    if os.path.exists(path):
        return g, "cached"
    t0 = time.time()
    emb, genes, expr = load_cohort(g)
    stage_all = emb.stage.values
    nac = stage_all != "Metastasis"
    st = stage_all[nac]
    res = dict(cohort=g, n=dict(pd.Series(st).value_counts()), samples=list(emb.index[nac]), stage=st)
    # part (a): frozen all-gene PC1
    s, astar, ordered = scale_score(emb.PC1.values[nac], st)
    Xt = expr.loc[genes].values
    Z = (Xt - Xt.mean(1, keepdims=True)) / (Xt.std(1, ddof=1, keepdims=True) + 1e-9)
    pca = PCA(n_components=10, random_state=0).fit(Z.T)
    pcre = pca.transform(Z.T)[:, 0]
    grp = pd.Series(s).groupby(st)
    r2 = float(sum(len(v) * (v.mean() - s.mean()) ** 2 for _, v in grp) / ((s - s.mean()) ** 2).sum())
    pa = part_a(s, st, rng_for(f"{g}__parta"))
    pa.update(cohort=g, astar=astar, ordered=ordered, pc1_var=float(pca.explained_variance_ratio_[0]),
              pc1_stage_R2=r2, pc1_recompute_absr=float(abs(np.corrcoef(pcre, emb.PC1.values)[0, 1])))
    res["part_a"], res["s"] = pa, s
    # secondary: N->C centroid axis (stage-informed; descriptive)
    P = emb[[c for c in emb.columns if c.startswith("PC")]].values[nac]
    cN, cC = P[st == "Normal"].mean(0), P[st == "Carcinoma"].mean(0)
    s2, a2, o2 = scale_score((P - cN) @ (cC - cN), st)
    sec = part_a(s2, st, rng_for(f"{g}__parta_secondary"))
    sec.update(cohort=g, astar=a2, ordered=o2, rho_with_primary=float(spearmanr(s, s2)[0]))
    res["part_a_secondary"] = sec
    # part (b): cross-fitted half scores
    perm = rng_for(f"{g}__genesplit").permutation(len(genes))
    H = [[genes[i] for i in perm[:1000]], [genes[i] for i in perm[1000:]]]
    hs = [scale_score(pc1(expr.loc[h].values)[nac], st) for h in H]
    res["H"], res["half_scores"] = H, hs
    res["half_rho"] = dict(h1_all=float(spearmanr(hs[0][0], s)[0]), h2_all=float(spearmanr(hs[1][0], s)[0]),
                           h1_h2=float(spearmanr(hs[0][0], hs[1][0])[0]))
    Y = expr.loc[:, nac].T
    n_all = int(nac.sum())
    parts, grids = [], {}
    for i in (0, 1):                     # genes of the OTHER half fitted against half i's score
        x, a_, _ = hs[i]
        df, grids[f"H{i + 1}"] = score_fits(x, Y[H[1 - i]], a_, n_all, rng_for(f"{g}__grid__H{i + 1}"))
        df["score"] = f"H{i + 1}"
        parts.append(df)
    res["fits"] = finish(pd.concat(parts))
    res["ordered_b"] = bool(hs[0][2] and hs[1][2])
    df, grids["all"] = score_fits(s, Y[genes], astar, n_all, rng_for(f"{g}__grid__all"))
    df["score"] = "all"
    res["fits_allgene"] = finish(df)
    res["grids"] = grids
    gm = expr.mean(1)
    res["expressed"] = set(gm.index[gm > gm.quantile(0.25)])
    res["genes"] = genes
    res["stage_means"] = expr.loc[:, nac].T.groupby(st).mean().T
    res["Y"] = Y
    res["seconds"] = round(time.time() - t0)
    with open(path + ".tmp", "wb") as fh:
        pickle.dump(res, fh)
    os.replace(path + ".tmp", path)
    return g, f"done in {res['seconds']} s"


def extra_fits(R, genes, mode):
    """Fit genes not in this cohort's top-2,000 against its H1 score (primary) or all-gene score (sensitivity)."""
    genes = [x for x in genes if x in R["expressed"] and x not in set(R["genes"]) and x in R["Y"].columns
             and not R["Y"][x].isna().any()]
    if not genes:
        return pd.DataFrame()
    if mode == "primary":
        x, a_, _ = R["half_scores"][0]
        rg, nulls = R["grids"]["H1"]
    else:
        x, a_ = R["s"], R["part_a"]["astar"]
        rg, nulls = R["grids"]["all"]
    df, _ = fit_genes(x, R["Y"][genes].values, a_, len(x))
    df.index = genes
    df["p_boot"] = [grid_p(r, d, rg, nulls) for r, d in zip(df.r_std, df.dAIC)]
    df["score"] = "H1-extra" if mode == "primary" else "all-extra"
    return df


def meets_switch_nominal(r):
    nonmono = (r.aic_spl <= min(r.aic_lin, r.aic_sig) - 10) and r.spline_turn
    return bool(r.p_assoc < 0.05 and r.dAIC >= 10 and r.W <= 0.25 and 0 <= r.x0 <= 1 and r.amp >= 1
                and r.resolvable and not nonmono and r.p_boot < 0.05)


def nonmono_nominal(r):
    return bool(r.p_assoc < 0.05 and r.aic_spl <= min(r.aic_lin, r.aic_sig) - 10 and r.spline_turn)


def adeno_transient(R, gene):
    m = R["stage_means"]
    if gene not in m.index:
        return "", np.nan, np.nan
    dAN, dAC = m.loc[gene, "Adenoma"] - m.loc[gene, "Normal"], m.loc[gene, "Adenoma"] - m.loc[gene, "Carcinoma"]
    lab = ("up in adenoma only" if dAN >= 0.5 and dAC >= 0.5 else
           "down in adenoma only" if dAN <= -0.5 and dAC <= -0.5 else "")
    return lab, float(dAN), float(dAC)


def interpretable(R, interval):
    nm = "NA" if interval == "N->A" else "AC"
    return bool(R["part_a"]["ordered"] and R["ordered_b"]
                and R["part_a"][f"{nm}_call"] in ("continuous", "overlapping"))


def replication(Rs, mode):
    key = "fits" if mode == "primary" else "fits_allgene"
    coh = list(Rs)
    disc = []
    for c in coh:
        f = Rs[c][key]
        for gname, r in f[f.cand_switch].iterrows():
            disc.append((c, gname, r.interval, r.direction, interpretable(Rs[c], r.interval)))
    need = {g: set() for g in coh}
    nm_disc = []
    for c in coh:
        f = Rs[c][key]
        for gname, r in f[f.cls == "non-monotonic"].iterrows():
            nm_disc.append((c, gname))
    for c, gname, *_ in disc:
        for o in coh:
            if o != c:
                need[o].add(gname)
    for c, gname in nm_disc:
        for o in coh:
            if o != c:
                need[o].add(gname)
    ext = {o: extra_fits(Rs[o], sorted(need[o]), mode) for o in coh}

    def row(o, gname):
        f = Rs[o][key]
        if gname in f.index:
            return f.loc[gname]
        if len(ext[o]) and gname in ext[o].index:
            return ext[o].loc[gname]
        return None

    sw = []
    for c, gname, interval, direction, interp in disc:
        reps, detail = [], {}
        for o in coh:
            if o == c:
                continue
            r = row(o, gname)
            if r is None:
                detail[o] = "not testable"; continue
            if not Rs[o]["ordered_b"]:
                detail[o] = "not ordered"; continue
            iv = "N->A" if r.x0 <= r.astar else "A->C"
            ok = meets_switch_nominal(r) and r.direction == direction and iv == interval
            if ok and not interpretable(Rs[o], interval):
                detail[o] = "meets rule; not interpretable (gap)"; continue
            detail[o] = "replicates" if ok else "no"
            if ok:
                reps.append(o)
        sw.append(dict(gene=gname, discovery=c, interval=interval, direction="up" if direction > 0 else "down",
                       discovery_interpretable=interp, n_replicating=len(reps), replicating=";".join(reps),
                       **{f"rep_{o}": detail.get(o, "discovery") for o in coh}))
    sw = pd.DataFrame(sw)
    if len(sw):
        sw["replicated"] = sw.discovery_interpretable & (sw.n_replicating >= 2)
    nmr = []
    for c, gname in nm_disc:
        r0 = Rs[c][key].loc[gname]
        lab, dAN, dAC = adeno_transient(Rs[c], gname)
        reps, detail = [], {}
        for o in coh:
            if o == c:
                continue
            r = row(o, gname)
            if r is None or not Rs[o]["ordered_b"]:
                detail[o] = "not testable" if r is None else "not ordered"; continue
            lab_o, *_ = adeno_transient(Rs[o], gname)
            ok = nonmono_nominal(r) and r.turn_type == r0.turn_type and (lab == "" or lab_o == lab)
            detail[o] = "replicates" if ok else "no"
            if ok:
                reps.append(o)
        nmr.append(dict(gene=gname, discovery=c, pattern=r0.turn_type, turn_x=r0.turn_x, astar=r0.astar,
                        adenoma_transient=lab, A_minus_N=dAN, A_minus_C=dAC, n_replicating=len(reps),
                        replicating=";".join(reps), **{f"rep_{o}": detail.get(o, "discovery") for o in coh}))
    nmr = pd.DataFrame(nmr)
    if len(nmr):
        nmr["replicated"] = nmr.n_replicating >= 2
    return sw, nmr


def unique_replicated(sw, cols=("gene", "interval", "direction")):
    if not len(sw):
        return sw
    r = sw[sw.replicated]
    if not len(r):
        return r
    return (r.sort_values("n_replicating", ascending=False).groupby(list(cols), as_index=False)
            .agg(discovery_cohorts=("discovery", lambda v: ";".join(sorted(set(v)))),
                 max_replicating=("n_replicating", "max")))


def gene_level_calls(R):
    f = R["fits"]
    P = int((f.q_assoc < 0.05).sum())
    out = {}
    for tr, *_ in TRS:
        nc = int((f.cand_switch & (f.interval == tr)).sum())
        ng = int((f.cls == "gradual").sum())
        if not R["ordered_b"]:
            call = "not ordered"
        elif P > 0 and nc / P >= 0.10:
            call = "switch-like" if interpretable(R, tr) else "not interpretable (possible gap artefact)"
        elif P > 0 and ng / P >= 0.5 and nc / P < 0.02:
            call = "gradual"
        else:
            call = "inconclusive"
        out[tr] = dict(P=P, n_cand_interval=nc, n_gradual=ng, call=call)
    return out


def enrichment(lists, background):
    try:
        import decoupler as dc
        hm = dc.op.hallmark(organism="human")
        sets = {f"HALLMARK:{k}": set(v) for k, v in hm.groupby("source").target}
        try:
            ms = dc.op.resource("MSigDB", organism="human")
            col = [c for c in ms.columns if "collection" in c.lower()]
            if col:
                rx = ms[ms[col[0]].astype(str).str.contains("reactome", case=False)]
                sets.update({f"REACTOME:{k}": set(v) for k, v in rx.groupby("geneset").genesymbol})
        except Exception as e:  # noqa: BLE001
            print(f"[enrich] Reactome not available: {e}", flush=True)
    except Exception as e:  # noqa: BLE001
        return pd.DataFrame(), f"gene sets not available: {e}"
    from scipy.stats import hypergeom
    bg = set(background)
    rows = []
    for name, genes in lists.items():
        L = set(genes) & bg
        if len(L) < 5:
            continue
        res = []
        for sname, S in sets.items():
            S2 = S & bg
            if not 10 <= len(S2) <= 500:
                continue
            k = len(L & S2)
            res.append(dict(list=name, gene_set=sname, k=k, list_size=len(L), set_size=len(S2),
                            p=hypergeom.sf(k - 1, len(bg), len(S2), len(L)) if k else 1.0,
                            genes=";".join(sorted(L & S2))))
        r = pd.DataFrame(res)
        r["q"] = bh(r.p.values)
        rows.append(r.sort_values("p"))
    nsets = sum(1 for _ in sets)
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()), f"{nsets} gene sets (Hallmark + Reactome where available)"


def run_real():
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    with Pool(4) as pool:
        for g, msg in pool.imap_unordered(stage1, list(COH8)):
            print(f"[a8] {g}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summarise_real(t0)


def summarise_real(t0=None):
    t0 = t0 or time.time()
    Rs = {g: pickle.load(open(os.path.join(CACHE, f"{g}.pkl"), "rb")) for g in COH8}
    pa = pd.DataFrame([Rs[g]["part_a"] for g in COH8])
    pa.to_csv(os.path.join(OUT, "continuity_results.csv"), index=False)
    pd.DataFrame([Rs[g]["part_a_secondary"] for g in COH8]).to_csv(os.path.join(OUT, "continuity_results_secondary.csv"), index=False)
    sc = []
    for g in COH8:
        R = Rs[g]
        sc.append(pd.DataFrame(dict(cohort=g, sample=R["samples"], stage=R["stage"], s_allgene=R["s"],
                                    s_H1=R["half_scores"][0][0], s_H2=R["half_scores"][1][0])))
        pd.DataFrame(dict(gene=R["H"][0] + R["H"][1], half=["H1"] * 1000 + ["H2"] * 1000)).to_csv(
            os.path.join(OUT, f"gene_split_{g}.csv"), index=False)
        R["fits"].to_csv(os.path.join(OUT, f"gene_fits_{g}.csv"), index_label="gene")
        R["fits_allgene"].to_csv(os.path.join(OUT, f"gene_fits_{g}_allgene.csv"), index_label="gene")
    pd.concat(sc).to_csv(os.path.join(OUT, "progression_scores.csv"), index=False)
    pd.DataFrame([dict(cohort=g, **Rs[g]["half_rho"], H1_astar=Rs[g]["half_scores"][0][1],
                       H2_astar=Rs[g]["half_scores"][1][1], ordered_b=Rs[g]["ordered_b"]) for g in COH8]).to_csv(
        os.path.join(OUT, "half_scores.csv"), index=False)

    sw, nmr = replication(Rs, "primary")
    sw_s, nmr_s = replication(Rs, "sensitivity")
    cands = []
    for g in COH8:
        f = Rs[g]["fits"]
        c = f[f.cand_switch].copy()
        c["cohort"] = g
        c["interpretable"] = [interpretable(Rs[g], iv) for iv in c.interval]
        cands.append(c)
    pd.concat(cands).to_csv(os.path.join(OUT, "switch_genes_candidates.csv"), index_label="gene")
    sw.to_csv(os.path.join(OUT, "switch_genes_replicated.csv"), index=False)
    sw_s.to_csv(os.path.join(OUT, "switch_genes_replicated_allgene.csv"), index=False)
    nmr.to_csv(os.path.join(OUT, "nonmonotonic_genes_replicated.csv"), index=False)
    usw, usw_s = unique_replicated(sw), unique_replicated(sw_s)
    unm = (nmr[nmr.replicated].groupby(["gene", "pattern", "adenoma_transient"], as_index=False)
           .agg(discovery_cohorts=("discovery", lambda v: ";".join(sorted(set(v)))),
                max_replicating=("n_replicating", "max"))) if len(nmr) and nmr.replicated.any() else pd.DataFrame()

    # enrichment
    allg = set().union(*[set(Rs[g]["genes"]) for g in COH8])
    bg = [x for x in allg if sum(x in Rs[g]["expressed"] for g in COH8) >= 3]
    lists = {}
    if len(usw):
        for (iv, d), v in usw.groupby(["interval", "direction"]):
            lists[f"switch {iv} {d}"] = list(v.gene)
    if len(unm):
        for p_, v in unm.groupby("pattern"):
            lists[f"non-monotonic {p_}"] = list(v.gene)
        at = unm[unm.adenoma_transient != ""]
        for lab, v in at.groupby("adenoma_transient"):
            lists[f"adenoma-transient: {lab}"] = list(v.gene)
    enr, enr_msg = enrichment(lists, bg)
    enr.to_csv(os.path.join(OUT, "enrichment_replicated.csv"), index=False)

    # section 5
    glc = {g: gene_level_calls(Rs[g]) for g in COH8}
    concl = []
    for tr, x, y, nm in TRS:
        oa = [g for g in COH8 if Rs[g]["part_a"]["ordered"]]
        ob = [g for g in COH8 if Rs[g]["ordered_b"]]
        pop_cont = sum(Rs[g]["part_a"][f"{nm}_call"] in ("continuous", "overlapping") for g in oa)
        pop_disc = sum(Rs[g]["part_a"][f"{nm}_call"] == "discontinuous" for g in oa)
        gene_grad = sum(glc[g][tr]["call"] == "gradual" for g in ob)
        nrep = int((usw.interval == tr).sum()) if len(usw) else 0
        if len(oa) < 3 or len(ob) < 3:
            c = "INCONCLUSIVE (fewer than 3 ordered cohorts)"
        else:
            grad = pop_cont >= 3 and gene_grad >= 3
            if grad and nrep >= 10:
                c = "GRADUAL POPULATION WITH A SWITCH-LIKE GENE PROGRAMME"
            elif grad:
                c = "SUPPORT FOR GRADUAL CHANGE"
            elif pop_disc >= 3 and nrep >= 10:
                c = "SUPPORT FOR SWITCH-LIKE CHANGE"
            elif pop_disc >= 3:
                c = "INCONCLUSIVE - sampling gap without gene-level switch (consistent with unsampled intermediates)"
            else:
                c = "INCONCLUSIVE"
        concl.append(dict(transition=tr, ordered_a=len(oa), ordered_b=len(ob), pop_continuous_or_overlapping=pop_cont,
                          pop_discontinuous=pop_disc, gene_gradual=gene_grad, replicated_switch_genes=nrep, conclusion=c))
    concl = pd.DataFrame(concl)
    concl.to_csv(os.path.join(OUT, "section5_conclusions.csv"), index=False)

    # summary text
    L = ["ADDENDUM 8 - Gradual vs switch-like change (8e8df5b; Amendment 1 0275b7a; Amendment 2 a3af4f9; D13)", "=" * 78,
         "Stage-blind adjusted data; N/A/C samples (GSE41258 metastases excluded). Part (a): frozen all-gene PC1.",
         "Part (b): cross-fitted half scores (primary); all-gene score = sensitivity. Bootstrap grid 21 x 1,999.",
         f"Seed {SEED}. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("numpy", "scipy", "scikit-learn", "patsy", "diptest", "decoupler")),
         "Samples (N/A/C): " + "; ".join(f"{g} {Rs[g]['n'].get('Normal')}/{Rs[g]['n'].get('Adenoma')}/{Rs[g]['n'].get('Carcinoma')}" for g in COH8), "",
         "PART (a) - continuity along adjusted PC1 (s: normal median 0, carcinoma median 1)"]
    for _, r in pa.iterrows():
        L.append(f"  {r.cohort:<9} a* = {r.astar:.2f} ({'ordered' if r.ordered else 'NOT ORDERED'}); PC1 var {100 * r.pc1_var:.1f}%, "
                 f"stage R2 {r.pc1_stage_R2:.2f}; dip Holm p: all {r.dip_holm_all:.4f}, N+A {r.dip_holm_NA:.4f}, A+C {r.dip_holm_AC:.4f}")
        for tr, x, y, nm in TRS:
            L.append(f"      {tr}: AUC {r[f'{nm}_AUC']:.2f}; VR {r[f'{nm}_VR']:.2f} [{r[f'{nm}_VR_lo']:.2f}, {r[f'{nm}_VR_hi']:.2f}]; "
                     f"between medians: {r[f'{nm}_n_x_between']}/{r[f'{nm}_n_x']} {x[:4]}, {r[f'{nm}_n_y_between']}/{r[f'{nm}_n_y']} {y[:4]}"
                     f"  -> {r[f'{nm}_call'].upper()}")
    L += ["  (PC1 recomputed from the frozen genes matches the frozen PC1: |r| = "
          + ", ".join(f"{r.cohort} {r.pc1_recompute_absr:.4f}" for _, r in pa.iterrows()) + ")", "",
          "PART (b) - class proportions of progression-associated genes (cross-fitted)"]
    for g in COH8:
        f, R = Rs[g]["fits"], Rs[g]
        P = int((f.q_assoc < 0.05).sum())
        vc = f[f.q_assoc < 0.05].cls.value_counts()
        L.append(f"  {g:<9} half-score rho with all-gene {R['half_rho']['h1_all']:.3f}/{R['half_rho']['h2_all']:.3f}; "
                 f"half a* {R['half_scores'][0][1]:.2f}/{R['half_scores'][1][1]:.2f} ({'ordered' if R['ordered_b'] else 'NOT ORDERED'}); "
                 f"progression-associated {P}/2000")
        L.append("      " + ", ".join(f"{k} {100 * v / max(P, 1):.1f}%" for k, v in vc.items()))
        for tr, *_ in TRS:
            gl = glc[g][tr]
            L.append(f"      {tr}: candidate switch genes (FDR) {gl['n_cand_interval']} "
                     f"({100 * gl['n_cand_interval'] / max(P, 1):.1f}%); interpretable: {interpretable(R, tr)}; gene-level call: {gl['call'].upper()}")
    L += ["", "  Sensitivity (all-gene score): " + "; ".join(
        f"{g} cand {int(Rs[g]['fits_allgene'].cand_switch.sum())} vs cross-fitted {int(Rs[g]['fits'].cand_switch.sum())}, "
        f"gradual {100 * (Rs[g]['fits_allgene'].cls == 'gradual').sum() / max((Rs[g]['fits_allgene'].q_assoc < .05).sum(), 1):.0f}%"
        for g in COH8)]
    L += ["", "REPLICATED SWITCH GENES (interpretable discovery + >= 2 other interpretable cohorts):"]
    if len(usw):
        for _, r in usw.iterrows():
            L.append(f"  {r.gene:<12} {r.interval} {r.direction:<5} discovered in {r.discovery_cohorts}; replicating cohorts {r.max_replicating}")
    else:
        L.append("  none")
    if len(sw):
        L.append(f"  (candidate discoveries: {len(sw)}; interpretable {int(sw.discovery_interpretable.sum())}; "
                 f"not interpretable (gap) {int((~sw.discovery_interpretable).sum())})")
    L.append(f"  Sensitivity (all-gene score): {len(usw_s)} replicated switch genes"
             + (f"; overlap with primary {len(set(usw.gene) & set(usw_s.gene)) if len(usw) and len(usw_s) else 0}" if len(usw_s) else ""))
    L += ["", "REPLICATED NON-MONOTONIC GENES (>= 2 other cohorts, same pattern; descriptive):"]
    if len(unm):
        for p_ in ("peak", "trough"):
            v = unm[unm.pattern == p_]
            L.append(f"  {p_} ({len(v)}): " + ", ".join(v.gene))
        at = unm[unm.adenoma_transient != ""]
        for lab, v in at.groupby("adenoma_transient"):
            L.append(f"  adenoma-transient, {lab} ({len(v)}): " + ", ".join(v.gene))
    else:
        L.append("  none")
    L += ["", f"ENRICHMENT (descriptive; {enr_msg}; background {len(bg)} genes):"]
    if len(enr):
        for name, v in enr.groupby("list"):
            top = v[v.q < 0.05].head(8)
            L.append(f"  {name}: " + ("; ".join(f"{r.gene_set} (k={r.k}, q={r.q:.2g})" for _, r in top.iterrows())
                                     if len(top) else "no gene set with q < 0.05"))
    else:
        L.append("  no list with >= 5 genes" if not lists else "  not computed")
    L += ["", "SECTION 5 - OVERALL CONCLUSION PER TRANSITION"]
    for _, r in concl.iterrows():
        L.append(f"  {r.transition}: population continuous/overlapping {r.pop_continuous_or_overlapping}, discontinuous "
                 f"{r.pop_discontinuous} (of {r.ordered_a} ordered); gene-level gradual {r.gene_gradual} (of {r.ordered_b}); "
                 f"replicated switch genes {r.replicated_switch_genes}")
        L.append(f"      -> {r.conclusion}")
    L += ["", "Secondary score (N->C centroid axis; descriptive): " + "; ".join(
        f"{r['cohort']} " + ", ".join(f"{tr} {r[f'{nm}_call']}" for tr, _, _, nm in TRS)
        for r in [Rs[g]['part_a_secondary'] for g in COH8]),
          f"Run time {(time.time() - t0) / 60:.1f} min (summary stage)."]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))
    figures(Rs, usw, unm)


def figures(Rs, usw, unm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"Normal": "#2a78d6", "Adenoma": "#e0a100", "Carcinoma": "#c0392b"}
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.6))
    for ax, g in zip(axes, COH8):
        R = Rs[g]
        bins = np.linspace(np.min(R["s"]), np.max(R["s"]), 40)
        for stg in col:
            ax.hist(R["s"][R["stage"] == stg], bins=bins, alpha=0.6, color=col[stg], label=stg)
        ax.axvline(0, color="0.4", lw=.8); ax.axvline(1, color="0.4", lw=.8)
        ax.set_title(g, fontsize=9); ax.set_xlabel("progression score s (adjusted PC1)")
    axes[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "score_densities.png"), dpi=150); plt.close(fig)
    classes = ["switch-like", "step-like unresolvable", "gradual", "intermediate", "non-monotonic"]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    bottom = np.zeros(4)
    for c in classes:
        v = np.array([((Rs[g]["fits"].cls == c).sum() / max((Rs[g]["fits"].q_assoc < .05).sum(), 1)) for g in COH8])
        ax.bar(list(COH8), v, bottom=bottom, label=c); bottom += v
    ax.set_ylabel("share of progression-associated genes"); ax.legend(fontsize=7, bbox_to_anchor=(1, 1))
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "class_proportions.png"), dpi=150); plt.close(fig)
    show = list(usw.gene[:6]) if len(usw) else (list(unm.gene[:6]) if len(unm) else [])
    if show:
        fig, axes = plt.subplots(len(show), 4, figsize=(14, 2.4 * len(show)), squeeze=False)
        for i, gname in enumerate(show):
            for j, g in enumerate(COH8):
                ax = axes[i, j]
                R = Rs[g]
                if gname not in R["Y"].columns:
                    ax.axis("off"); continue
                hs = R["half_scores"]
                x = hs[1][0] if gname in R["H"][0] else hs[0][0]
                for stg in col:
                    m = R["stage"] == stg
                    ax.scatter(x[m], R["Y"][gname].values[m], s=6, color=col[stg])
                ax.set_title(f"{gname} - {g}", fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "example_trajectories.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "validate":
        validate()
    elif cmd == "run":
        run_real()
    elif cmd == "summary":
        summarise_real()
    else:
        raise SystemExit("Usage: addendum8_gradual_switch.py validate | run | summary")
