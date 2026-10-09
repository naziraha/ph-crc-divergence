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
RENV = r"C:\Users\naziraha\AppData\Local\anaconda3\envs\r-cms"


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
    E = RENV
    return ";".join([E, rf"{E}\Library\mingw-w64\bin", rf"{E}\Library\usr\bin", rf"{E}\Library\bin", rf"{E}\Scripts",
                     rf"{E}\lib\R\bin\x64", os.environ["PATH"]])


def run_r(*args):
    env = dict(os.environ, PATH=renv_path())
    r = subprocess.run([rf"{RENV}\lib\R\bin\x64\Rscript.exe", os.path.join(ROOT, "scripts", "addendum9_cms.R"), *args],
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


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        validate()
    else:
        raise SystemExit("Only 'validate' is implemented at this stage (addendum 9 section 6).")
