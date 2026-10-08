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
    A_l, A_s, A_p = aic(rl, n, 2), aic(rs, n, 4), aic(rsp, n, 4)
    W = LN81 / k
    win = np.array([np.sum((x >= c - w / 2) & (x <= c + w / 2)) for c, w in zip(x0, W)])
    sigma = np.sqrt(rl / (n - 2))
    return pd.DataFrame(dict(p_assoc=p_assoc, slope=slope, r_std=np.abs(slope) / sigma, aic_lin=A_l, aic_sig=A_s,
                             aic_spl=A_p, dAIC=A_l - A_s, a=a, b=b, k=k, x0=x0, W=W, amp=np.abs(b - a),
                             direction=np.sign(b - a), n_window=win,
                             resolvable=win >= max(5, 0.05 * n_all), spline_turn=change, turn_type=ext,
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

    L = ["ADDENDUM 8 - Section 6 validation on simulated data (addendum 8 8e8df5b; D13 248d295)", "=" * 72,
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


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        validate()
    else:
        raise SystemExit("Only 'validate' is implemented at this stage (addendum 8 section 6).")
