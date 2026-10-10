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


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        validate()
    else:
        raise SystemExit("Only 'validate' is implemented at this stage (Part B awaits PI review of Part A).")
