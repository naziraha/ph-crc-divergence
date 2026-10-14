"""
Step 13 - benchmark of the TDA test against conventional multivariate methods.

Prespecified before analysis: docs/addendum7_benchmark.md (commit c165512).
Methods (numpy/scipy):
  PERMDISP   - ANOVA F of distances to own-group centroid, 999 label permutations; D = mean dist later/earlier
  PERMANOVA  - Euclidean pseudo-F, 999 label permutations; R2
  stage separation - Ward (Euclidean), 2 clusters, ARI vs labels, 999 label permutations (location-type)
  MST test   - conventional counterpart of H0: VR H0 death times = Euclidean MST edge lengths. Same subsampling
               (k = 20, m), statistic (mean between - mean within), group alignment and 999 permutations as the
               step 5 test (tda_crc_pilot.permutation_test), with diagrams replaced by MST edge lengths and the
               diagram distance by the 1-D 1-Wasserstein distance between edge-length distributions.
               Modes: raw, shape-only (edge lengths / subsample median pairwise distance).
  For equal sample sizes the 1-D W1 equals the mean absolute difference of the sorted values (checked against
  scipy.stats.wasserstein_distance in the validation).

Usage:
  python scripts/step13_benchmark.py validate   - section 4 simulation only (no real data)
"""
import os
import sys
import time
import zlib
from importlib import metadata
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway, wasserstein_distance

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
OUT = os.path.join(ROOT, "results", "step13")

SEED, N_PERM, K_DRAWS, WORKERS = 42, 999, 20, 6
N_REP, DIM = 200, 10


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def perm_indicators(n1, n2, n_perm, rng):
    """Rows = permutations; 1 = sample assigned to group 2 (group sizes kept)."""
    N = n1 + n2
    X = np.zeros((n_perm, N))
    for i in range(n_perm):
        X[i, rng.permutation(N)[:n2]] = 1.0
    return X


def perm_p(obs, null):
    return float((1 + np.sum(null >= obs)) / (1 + len(null)))


# ---------------------------------------------------------------- PERMDISP
def _anova_F(z, X, n1, n2):
    """Two-group one-way ANOVA F of values z for each indicator row of X (1 = group 2)."""
    N = n1 + n2
    s2 = X @ z
    s1 = z.sum() - s2
    ss_tot = ((z - z.mean()) ** 2).sum()
    ss_b = n1 * (s1 / n1 - z.mean()) ** 2 + n2 * (s2 / n2 - z.mean()) ** 2
    return ss_b / ((ss_tot - ss_b) / (N - 2))


def permdisp(A, B, rng, n_perm=N_PERM):
    z = np.concatenate([np.linalg.norm(A - A.mean(0), axis=1), np.linalg.norm(B - B.mean(0), axis=1)])
    n1, n2 = len(A), len(B)
    x0 = np.r_[np.zeros(n1), np.ones(n2)]
    F = float(_anova_F(z, x0[None, :], n1, n2)[0])
    null = _anova_F(z, perm_indicators(n1, n2, n_perm, rng), n1, n2)
    D = float(z[n1:].mean() / z[:n1].mean())
    return F, perm_p(F, null), D


# ---------------------------------------------------------------- PERMANOVA
def _permanova_F(D2, X, n1, n2):
    N = n1 + n2
    ss_t = D2.sum() / 2 / N
    w2 = ((X @ D2) * X).sum(1) / 2
    Y = 1 - X
    w1 = ((Y @ D2) * Y).sum(1) / 2
    ss_w = w1 / n1 + w2 / n2
    return (ss_t - ss_w) / (ss_w / (N - 2)), (ss_t - ss_w) / ss_t


def permanova(A, B, rng, n_perm=N_PERM):
    P = np.vstack([A, B])
    D2 = squareform(pdist(P)) ** 2
    n1, n2 = len(A), len(B)
    x0 = np.r_[np.zeros(n1), np.ones(n2)][None, :]
    F, R2 = _permanova_F(D2, x0, n1, n2)
    null, _ = _permanova_F(D2, perm_indicators(n1, n2, n_perm, rng), n1, n2)
    return float(F[0]), perm_p(F[0], null), float(R2[0])


# ---------------------------------------------------------------- stage separation (Ward + ARI)
def _ari_2x2(c, X):
    """ARI between fixed 2-cluster labels c (0/1) and each indicator row of X (0/1)."""
    N = len(c)
    comb = lambda x: x * (x - 1) / 2
    n11 = X @ c
    n10 = X.sum(1) - n11
    n01 = c.sum() - n11
    n00 = N - n11 - n10 - n01
    sij = comb(n11) + comb(n10) + comb(n01) + comb(n00)
    a = comb(c.sum()) + comb(N - c.sum())
    b = comb(X.sum(1)) + comb(N - X.sum(1))
    exp = a * b / comb(N)
    return (sij - exp) / (0.5 * (a + b) - exp)


def stage_separation(A, B, rng, n_perm=N_PERM):
    P = np.vstack([A, B])
    c = (fcluster(linkage(P, method="ward"), 2, criterion="maxclust") == 2).astype(float)
    n1, n2 = len(A), len(B)
    ari = float(_ari_2x2(c, np.r_[np.zeros(n1), np.ones(n2)][None, :])[0])
    null = _ari_2x2(c, perm_indicators(n1, n2, n_perm, rng))
    return ari, perm_p(ari, null), c


# ---------------------------------------------------------------- MST edge-length test (H0 counterpart)
def mst_edges(points, normalize):
    D = squareform(pdist(points))
    if normalize:
        D = D / np.median(D[np.triu_indices_from(D, 1)])
    e = minimum_spanning_tree(D).data
    assert len(e) == len(points) - 1, "MST edge count (duplicate points?)"
    return np.sort(e)


def _mst_stat(emb, idx_s, idx_t, m, k, normalize, rng):
    S = np.array([mst_edges(emb[rng.choice(idx_s, size=m, replace=False)], normalize) for _ in range(k)])
    Tt = np.array([mst_edges(emb[rng.choice(idx_t, size=m, replace=False)], normalize) for _ in range(k)])
    between = np.abs(S[:, None, :] - Tt[None, :, :]).mean(2).ravel()
    iu = np.triu_indices(k, 1)
    within = np.r_[np.abs(S[:, None, :] - S[None, :, :]).mean(2)[iu], np.abs(Tt[:, None, :] - Tt[None, :, :]).mean(2)[iu]]
    return float(between.mean() - within.mean())


def align_groups(emb, idx_s, idx_t, normalize):
    """As tda_crc_pilot.align_groups (centre each group; shape-only: unit median pairwise distance)."""
    out = emb.copy()
    for idx in (idx_s, idx_t):
        out[idx] = emb[idx] - emb[idx].mean(0)
        if normalize:
            out[idx] /= np.median(pdist(out[idx]))
    return out


def mst_test(emb, idx_s, idx_t, m, normalize, rng, k=K_DRAWS, n_perm=N_PERM):
    emb = align_groups(emb, idx_s, idx_t, normalize)
    obs = _mst_stat(emb, idx_s, idx_t, m, k, normalize, rng)
    pool, ns = np.concatenate([idx_s, idx_t]), len(idx_s)
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(pool)
        null[i] = _mst_stat(emb, perm[:ns], perm[ns:], m, k, normalize, rng)
    return obs, perm_p(obs, null), float((obs - null.mean()) / (null.std(ddof=1) + 1e-12))


# ================================================================ section 4 validation
SCEN = {
    "V1": "identical Gaussian groups (SD 1), n 50/50",
    "V2": "group 2 scaled x1.5, same centroid, n 50/50",
    "V3": "group 2 shifted 1 SD along one axis, n 50/50",
    "V4": "group 1 one Gaussian (SD 1); group 2 two clusters (SD 0.5) at +/-2 on one axis, n 50/50",
    "V5": "identical Gaussian groups, unequal n 48/181",
    "V6": "group 2 shifted 4 SD along one axis, n 50/50",
}


def simulate(scen, rng):
    n1, n2 = (48, 181) if scen == "V5" else (50, 50)
    A = rng.standard_normal((n1, DIM))
    B = rng.standard_normal((n2, DIM))
    if scen == "V2":
        B *= 1.5
    elif scen == "V3":
        B[:, 0] += 1.0
    elif scen == "V4":
        B *= 0.5
        B[: n2 // 2, 0] -= 2.0
        B[n2 // 2:, 0] += 2.0
    elif scen == "V6":
        B[:, 0] += 4.0
    return A, B


def run_rep(task):
    scen, rep = task
    base = f"sim__{scen}__rep{rep:03d}"
    A, B = simulate(scen, rng_for(base + "__data"))
    F, p_disp, D = permdisp(A, B, rng_for(base + "__permdisp"))
    Fp, p_perm, R2 = permanova(A, B, rng_for(base + "__permanova"))
    ari, p_ari, _ = stage_separation(A, B, rng_for(base + "__stagesep"))
    emb = np.vstack([A, B])
    idx_s, idx_t = np.arange(len(A)), np.arange(len(A), len(A) + len(B))
    m = int(0.8 * min(len(A), len(B)))
    _, p_mr, z_mr = mst_test(emb, idx_s, idx_t, m, False, rng_for(base + "__mst_raw"))
    _, p_ms, z_ms = mst_test(emb, idx_s, idx_t, m, True, rng_for(base + "__mst_shape"))
    return dict(scenario=scen, rep=rep, permdisp_F=F, permdisp_p=p_disp, D=D, permanova_F=Fp, permanova_p=p_perm,
                R2=R2, ARI=ari, stagesep_p=p_ari, mst_raw_p=p_mr, mst_raw_z=z_mr, mst_shape_p=p_ms, mst_shape_z=z_ms)


def exact_checks():
    """Textbook calculations, W1 shortcut, and MST = GUDHI H0 deaths (simulated data only)."""
    import tda_crc_pilot as T
    rng = rng_for("sim__checks")
    A, B = simulate("V2", rng)
    L = []
    # PERMDISP F vs scipy f_oneway on centroid distances
    da, db = np.linalg.norm(A - A.mean(0), axis=1), np.linalg.norm(B - B.mean(0), axis=1)
    F, _, _ = permdisp(A, B, rng, n_perm=9)
    ok1 = np.isclose(F, f_oneway(da, db).statistic)
    L.append(("PERMDISP F = scipy f_oneway on centroid distances", f"{F:.6f} vs {f_oneway(da, db).statistic:.6f}", ok1))
    # PERMANOVA pseudo-F vs centroid sums of squares (Euclidean)
    P = np.vstack([A, B])
    ssw = ((A - A.mean(0)) ** 2).sum() + ((B - B.mean(0)) ** 2).sum()
    sst = ((P - P.mean(0)) ** 2).sum()
    Ftb = (sst - ssw) / (ssw / (len(P) - 2))
    Fp, _, R2 = permanova(A, B, rng, n_perm=9)
    ok2 = np.isclose(Fp, Ftb) and np.isclose(R2, (sst - ssw) / sst)
    L.append(("PERMANOVA pseudo-F / R2 = centroid sums of squares", f"{Fp:.6f} vs {Ftb:.6f}; R2 {R2:.6f}", ok2))
    # ARI vs scikit-learn
    from sklearn.metrics import adjusted_rand_score
    ari, _, c = stage_separation(A, B, rng, n_perm=9)
    lab = np.r_[np.zeros(len(A)), np.ones(len(B))]
    ok3 = np.isclose(ari, adjusted_rand_score(lab, c))
    L.append(("ARI = sklearn adjusted_rand_score", f"{ari:.6f} vs {adjusted_rand_score(lab, c):.6f}", ok3))
    # 1-D W1 shortcut vs scipy
    e1, e2 = mst_edges(A[:40], False), mst_edges(B[:40], False)
    w_fast, w_sp = np.abs(e1 - e2).mean(), wasserstein_distance(e1, e2)
    ok4 = np.isclose(w_fast, w_sp)
    L.append(("1-D W1 shortcut = scipy wasserstein_distance", f"{w_fast:.8f} vs {w_sp:.8f}", ok4))
    # MST edge lengths = GUDHI Rips H0 death times (raw and shape-only)
    for norm in (False, True):
        for X, nm in ((A[:40], "Gaussian"), (simulate("V4", rng)[1][:40], "two-cluster")):
            e = mst_edges(X, norm)
            h0 = np.sort(T.persistence(X, norm)[0][:, 1])
            ok = len(e) == len(h0) and np.allclose(e, h0, atol=1e-6)
            L.append((f"MST edges = GUDHI H0 deaths ({nm}, {'shape-only' if norm else 'raw'})",
                      f"{len(e)} edges, max |diff| {np.max(np.abs(e - h0)) if len(e) == len(h0) else np.nan:.2e}", ok))
    return L


def validate():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    checks = exact_checks()
    for c in checks:
        print(("PASS " if c[2] else "FAIL ") + c[0] + ": " + c[1], flush=True)
    tasks = [(s, r) for s in SCEN for r in range(N_REP)]
    rows = []
    with Pool(WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(run_rep, tasks, chunksize=2), 1):
            rows.append(r)
            if i % 100 == 0:
                print(f"[validate] {i}/{len(tasks)} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    df = pd.DataFrame(rows).sort_values(["scenario", "rep"])
    df.to_csv(os.path.join(OUT, "validation_replicates.csv"), index=False)

    meths = [("PERMDISP", "permdisp_p"), ("PERMANOVA", "permanova_p"), ("Stage separation (ARI)", "stagesep_p"),
             ("MST raw", "mst_raw_p"), ("MST shape-only", "mst_shape_p")]
    rate = {s: {m: float((df[df.scenario == s][c] < 0.05).mean()) for m, c in meths} for s in SCEN}
    se = np.sqrt(0.05 * 0.95 / N_REP)
    t1 = 0.05 + 2.5 * se
    crit = []  # (scenario, method, requirement, value, pass)
    for s in ("V1", "V5"):
        for m, _ in meths:
            crit.append((s, m, f"type I <= {t1:.3f}", rate[s][m], rate[s][m] <= t1))
    crit.append(("V3", "PERMDISP", f"type I <= {t1:.3f}", rate["V3"]["PERMDISP"], rate["V3"]["PERMDISP"] <= t1))
    crit.append(("V2", "PERMDISP", "power >= 0.80", rate["V2"]["PERMDISP"], rate["V2"]["PERMDISP"] >= 0.8))
    dmed = float(df[df.scenario == "V2"].D.median())
    crit.append(("V2", "PERMDISP", "D > 1 (median)", dmed, dmed > 1))
    crit.append(("V3", "PERMANOVA", "power >= 0.80", rate["V3"]["PERMANOVA"], rate["V3"]["PERMANOVA"] >= 0.8))
    for m in ("MST raw", "MST shape-only"):
        crit.append(("V4", m, "power >= 0.80", rate["V4"][m], rate["V4"][m] >= 0.8))
    crit.append(("V6", "Stage separation (ARI)", "power >= 0.80", rate["V6"]["Stage separation (ARI)"],
                 rate["V6"]["Stage separation (ARI)"] >= 0.8))
    cdf = pd.DataFrame(crit, columns=["scenario", "method", "requirement", "value", "pass"])
    cdf.to_csv(os.path.join(OUT, "validation_criteria.csv"), index=False)
    all_pass = bool(cdf["pass"].all()) and all(c[2] for c in checks)

    L = ["STEP 13 - Section 4 validation on simulated data (addendum 7, commit c165512)", "=" * 70,
         f"No real data used. {DIM} dimensions; {N_REP} replicates per scenario; {N_PERM} permutations; "
         f"MST test k = {K_DRAWS}, m = 80% of smaller group (40; 38 in V5).",
         f"Seed {SEED} (SeedSequence per scenario x replicate x method). "
         + ", ".join(f"{p} {metadata.version(p)}" for p in ("numpy", "scipy", "scikit-learn", "gudhi")),
         f"Monte Carlo SE at 5% = {se:.4f}; type I threshold = {t1:.4f}.", "",
         "Exact checks:"] + [f"  {'PASS' if c[2] else 'FAIL'}  {c[0]}: {c[1]}" for c in checks] + ["",
         "Rejection rates (p < 0.05):", f"  {'scenario':<9}" + "".join(f"{m:>24}" for m, _ in meths)]
    for s in SCEN:
        L.append(f"  {s:<9}" + "".join(f"{rate[s][m]:>24.3f}" for m, _ in meths) + f"   ({SCEN[s]})")
    L += ["", "Acceptance criteria:"]
    L += [f"  {'PASS' if r['pass'] else 'FAIL'}  {r.scenario} {r.method:<24} {r.requirement:<18} {r.value:.3f}"
          for _, r in cdf.iterrows()]
    L += ["", f"OVERALL: {'ALL CRITERIA PASSED - real-data benchmark may proceed' if all_pass else 'FAILED - do not analyse real data'}",
          f"Run time {(time.time() - t0) / 60:.1f} min.",
          "Outputs: validation_summary.txt, validation_criteria.csv, validation_replicates.csv"]
    open(os.path.join(OUT, "validation_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ================================================================ real-data benchmark (addendum 7 sections 2-3, 5)
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
CK = os.path.join(OUT, "checkpoints")
VERS = ("unadjusted", "adjusted")
COH = {  # set, expected counts, transitions (as step 5)
    "GSE41258": ("main", {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39},
                 [("Normal", "Adenoma"), ("Adenoma", "Carcinoma"), ("Carcinoma", "Metastasis")]),
    "GSE20916": ("bulk", {"Normal": 24, "Adenoma": 45, "Carcinoma": 36}, [("Normal", "Adenoma"), ("Adenoma", "Carcinoma")]),
    "GSE37364": ("bulk", {"Normal": 38, "Adenoma": 29, "Carcinoma": 27}, [("Normal", "Adenoma"), ("Adenoma", "Carcinoma")]),
    "GSE117606": ("main", {"Normal": 65, "Adenoma": 59, "Carcinoma": 71}, [("Normal", "Adenoma"), ("Adenoma", "Carcinoma")]),
}
METHODS = ("permdisp", "permanova", "stagesep", "mst_raw", "mst_shape", "sizematched")


def load(gse, v):
    st, exp, _ = COH[gse]
    df = pd.read_csv(os.path.join(EMB, f"{gse}_{st}_{v}_g2000_pc10.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != exp:
        raise SystemExit(f"Stage counts changed: {gse} {v}: {df['stage'].value_counts().to_dict()}")
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    m = max(10, int(0.8 * min(exp.values())))  # as step 5
    return emb, df["stage"].values, m


def bench_units():
    return [(g, v, s, t, meth) for g in COH for v in VERS for s, t in COH[g][2] for meth in METHODS]


def bunit_name(g, v, s, t, meth):
    return f"{g}__{v[:5]}__{s[0]}-{t[0]}__{meth}"


def run_bunit(u):
    import json
    g, v, s, t, meth = u
    name = bunit_name(*u)
    path = os.path.join(CK, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = rng_for(name)
    emb, st, m = load(g, v)
    i_s, i_t = np.where(st == s)[0], np.where(st == t)[0]
    A, B = emb[i_s], emb[i_t]
    rec = dict(unit=name, cohort=g, version=v, transition=f"{s}->{t}", method=meth, n_from=len(i_s), n_to=len(i_t), m=m)
    if meth == "permdisp":
        rec["stat"], rec["p"], rec["effect"] = permdisp(A, B, rng)
    elif meth == "permanova":
        rec["stat"], rec["p"], rec["effect"] = permanova(A, B, rng)
    elif meth == "stagesep":
        ari, p, _ = stage_separation(A, B, rng)
        rec["stat"], rec["p"], rec["effect"] = ari, p, ari
    elif meth in ("mst_raw", "mst_shape"):
        obs, p, z = mst_test(emb, i_s, i_t, m, meth == "mst_shape", rng)
        rec["stat"], rec["p"], rec["effect"] = obs, p, z
    else:  # size-matched PERMDISP and PERMANOVA on the same k subsample pairs
        res = []
        for _ in range(K_DRAWS):
            a, b = emb[rng.choice(i_s, m, replace=False)], emb[rng.choice(i_t, m, replace=False)]
            res.append(permdisp(a, b, rng) + permanova(a, b, rng))
        res = np.array(res)
        med = np.median(res, axis=0)
        rec.update(sm_permdisp_F=med[0], sm_permdisp_p=med[1], sm_D=med[2], sm_permanova_F=med[3],
                   sm_permanova_p=med[4], sm_R2=med[5], sm_permdisp_frac_sig=float((res[:, 1] < .05).mean()),
                   sm_permanova_frac_sig=float((res[:, 4] < .05).mean()))
    rec = {k: (float(x) if isinstance(x, (np.floating, np.integer)) else x) for k, x in rec.items()}
    rec["seconds"] = round(time.time() - t0, 1)
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"done in {rec['seconds']:.0f} s"


def real_crosscheck():
    """MST edge lengths = GUDHI H0 deaths on one real subsample (addendum 7 section 3.4)."""
    import tda_crc_pilot as T
    emb, st, m = load("GSE41258", "unadjusted")
    rng = rng_for("GSE41258__crosscheck")
    X = emb[rng.choice(np.where(st == "Normal")[0], m, replace=False)]
    out = []
    for norm in (False, True):
        e, h0 = mst_edges(X, norm), np.sort(T.persistence(X, norm)[0][:, 1])
        ok = len(e) == len(h0) and np.allclose(e, h0, atol=1e-6)
        out.append((f"GSE41258 unadjusted Normal subsample (m={m}), {'shape-only' if norm else 'raw'}",
                    f"{len(e)} MST edges vs {len(h0)} H0 deaths, max |diff| "
                    f"{np.max(np.abs(e - h0)) if len(e) == len(h0) else float('nan'):.2e}", ok))
    return out


def tda_results():
    rows = []
    S5 = os.path.join(ROOT, "results", "step05")
    for g in COH:
        f = os.path.join(S5, "permutation_results_all.csv") if g == "GSE41258" else os.path.join(S5, g, "permutation_results_all.csv")
        d = pd.read_csv(f)
        d = d[(d.analysis == "primary") & d.version.isin(VERS)]
        for _, r in d.iterrows():
            rows.append(dict(cohort=g, version=r.version, transition=r.transition,
                             mode="raw" if r["mode"] == "raw" else "shape", tda_z=r.z, tda_p_holm=r.p_holm))
    tda = pd.DataFrame(rows).pivot_table(index=["cohort", "version", "transition"], columns="mode",
                                         values=["tda_z", "tda_p_holm"])
    tda.columns = [f"{a}_{b}" for a, b in tda.columns]
    h0 = pd.read_csv(os.path.join(ROOT, "results", "addendum5", "part4_dims", "dims_results_all.csv"))
    h0 = h0[h0.dim == "H0"].assign(mode=lambda x: np.where(x["mode"] == "raw", "raw", "shape"))
    h0 = h0.pivot_table(index=["cohort", "version", "transition"], columns="mode", values=["z", "p_holm"])
    h0.columns = [f"h0_{a}_{b}" for a, b in h0.columns]
    return tda.join(h0, how="left").reset_index()


def concord(x, y):
    """x = TDA significant, y = conventional significant (boolean arrays)."""
    x, y = np.asarray(x, bool), np.asarray(y, bool)
    n = len(x)
    a, b, c, d = (x & y).sum(), (x & ~y).sum(), (~x & y).sum(), (~x & ~y).sum()
    a1 = a / (a + b) if a + b else np.nan
    a2 = a / (a + c) if a + c else np.nan
    po = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n ** 2
    kappa = (po - pe) / (1 - pe) if pe < 1 else np.nan
    return dict(n=n, both=a, tda_only=b, conv_only=c, neither=d, a1=a1, a2=a2, kappa=kappa,
                largely_captured=bool(a1 >= .75 and a2 >= .75) if a1 == a1 and a2 == a2 else None)


def fmt_pct(x):
    return "n/a" if x != x else f"{100 * x:.0f}%"


def bench_summary():
    import json
    import tda_crc_pilot as T
    recs = [json.load(open(os.path.join(CK, f))) for f in sorted(os.listdir(CK)) if f.endswith(".json")]
    df = pd.DataFrame(recs)
    if len(df) < len(bench_units()):
        print(f"[summary] {len(df)}/{len(bench_units())} units finished"); return
    main = df[df.method != "sizematched"].copy()
    main["p_holm"] = np.nan
    for _, g in main.groupby(["cohort", "version", "method"]):
        main.loc[g.index, "p_holm"] = T.holm(g.p.values)
    sm = df[df.method == "sizematched"].copy()
    for c in ("sm_permdisp_p", "sm_permanova_p"):
        sm[c + "_holm"] = np.nan
        for _, g in sm.groupby(["cohort", "version"]):
            sm.loc[g.index, c + "_holm"] = T.holm(g[c].values)
    key = ["cohort", "version", "transition"]
    wide = main.pivot_table(index=key, columns="method", values=["stat", "effect", "p", "p_holm"])
    wide.columns = [f"{m}_{a}" for a, m in wide.columns]
    wide = wide.reset_index().merge(tda_results(), on=key, how="left").merge(
        sm[key + [c for c in sm.columns if c.startswith("sm_")]], on=key, how="left")
    order = {g: i for i, g in enumerate(COH)}
    wide["_o"] = wide.cohort.map(order) * 100 + wide.version.map({"unadjusted": 0, "adjusted": 1}) * 10 + \
        wide.transition.map({"Normal->Adenoma": 0, "Adenoma->Carcinoma": 1, "Carcinoma->Metastasis": 2})
    wide = wide.sort_values("_o").drop(columns="_o").reset_index(drop=True)
    assert len(wide) == 18 and wide.tda_p_holm_raw.notna().all()
    main.to_csv(os.path.join(OUT, "benchmark_results_long.csv"), index=False)
    wide.to_csv(os.path.join(OUT, "benchmark_results.csv"), index=False)
    sm.to_csv(os.path.join(OUT, "size_matched_results.csv"), index=False)

    sig = lambda c: (wide[c] < .05).values
    L = ["STEP 13 - Benchmark against conventional methods (addendum 7, commit c165512)", "=" * 70,
         "Embeddings: frozen 10-PC step 3 embeddings (unadjusted, stage-blind adjusted); transitions as step 5.",
         f"PERMDISP / PERMANOVA / stage separation on all samples; MST test and size-matched analysis with k = {K_DRAWS},",
         "m as step 5 (GSE41258 31, GSE20916 19, GSE37364 21, GSE117606 47). 999 permutations; Holm across transitions",
         "within cohort x version x method (x mode). Significant = Holm p < 0.05.",
         f"Seed {SEED} (per-unit SeedSequence). " + ", ".join(f"{p} {metadata.version(p)}" for p in ("numpy", "scipy", "gudhi")),
         "Samples per stage: " + "; ".join(f"{g} " + "/".join(str(x) for x in COH[g][1].values()) for g in COH), ""]
    L += ["GUDHI cross-check on one real subsample:"] + [f"  {'PASS' if c[2] else 'FAIL'}  {c[0]}: {c[1]}"
                                                         for c in real_crosscheck()] + [""]
    L.append("Per-comparison results (Holm p; * < 0.05):")
    L.append(f"  {'cohort':<9} {'version':<10} {'transition':<21} {'TDAraw':>7} {'TDAshp':>7} {'PDISP':>7} {'D':>5} "
             f"{'PANOVA':>7} {'R2':>5} {'StSep':>7} {'ARI':>5} {'MSTraw':>7} {'MSTshp':>7}")
    st = lambda p: f"{p:.3f}{'*' if p < .05 else ' '}"
    for _, r in wide.iterrows():
        L.append(f"  {r.cohort:<9} {r.version:<10} {r.transition:<21} {st(r.tda_p_holm_raw):>7} {st(r.tda_p_holm_shape):>7} "
                 f"{st(r.permdisp_p_holm):>7} {r.permdisp_effect:5.2f} {st(r.permanova_p_holm):>7} {r.permanova_effect:5.2f} "
                 f"{st(r.stagesep_p_holm):>7} {r.stagesep_effect:5.2f} {st(r.mst_raw_p_holm):>7} {st(r.mst_shape_p_holm):>7}")
    L += ["", "5.1 CONCORDANCE (18 comparisons). a1 = of TDA-significant, % also conventional-significant;",
          "    a2 = of conventional-significant, % also TDA-significant. 'Largely captured' if a1 >= 75% AND a2 >= 75%."]
    conc = []
    for lab, x, y in (("P1 TDA raw vs PERMDISP", "tda_p_holm_raw", "permdisp_p_holm"),
                      ("P2 TDA raw vs MST raw", "tda_p_holm_raw", "mst_raw_p_holm"),
                      ("P3 TDA shape-only vs MST shape-only", "tda_p_holm_shape", "mst_shape_p_holm")):
        c = concord(sig(x), sig(y)); c["pairing"] = lab; conc.append(c)
        concl = ("LARGELY CAPTURED" if c["largely_captured"] else "not applicable (denominator 0)"
                 if c["largely_captured"] is None else "NOT largely captured - discordant comparisons listed below")
        L.append(f"  {lab}: both {c['both']}, TDA only {c['tda_only']}, conventional only {c['conv_only']}, "
                 f"neither {c['neither']}; a1 {fmt_pct(c['a1'])}, a2 {fmt_pct(c['a2'])}, kappa "
                 f"{'n/a' if c['kappa'] != c['kappa'] else round(c['kappa'], 2)} -> {concl}")
        disc = wide[sig(x) != sig(y)]
        for _, r in disc.iterrows():
            L.append(f"      discordant: {r.cohort} {r.version} {r.transition}: TDA {r[x]:.3f}, conventional {r[y]:.3f}")
    L += ["", "  Descriptive extra (no conclusion): MST vs H0-only TDA (addendum 5 part 4; 3 cohorts, 14 comparisons)"]
    h = wide[wide.h0_p_holm_raw.notna()]
    for lab, x, y in (("MST raw vs H0-only raw", "h0_p_holm_raw", "mst_raw_p_holm"),
                      ("MST shape vs H0-only shape", "h0_p_holm_shape", "mst_shape_p_holm")):
        c = concord(h[x] < .05, h[y] < .05); c["pairing"] = lab + " (descriptive)"; conc.append(c)
        L.append(f"  {lab}: n {c['n']}; both {c['both']}, H0-TDA only {c['tda_only']}, MST only {c['conv_only']}, "
                 f"neither {c['neither']}; a1 {fmt_pct(c['a1'])}, a2 {fmt_pct(c['a2'])}, kappa "
                 f"{'n/a' if c['kappa'] != c['kappa'] else round(c['kappa'], 2)}")
    pd.DataFrame(conc).to_csv(os.path.join(OUT, "concordance.csv"), index=False)

    L += ["", "5.1(2) ADENOMA -> CARCINOMA DIRECTION (supported if PERMDISP Holm p < 0.05 AND D > 1):"]
    for _, r in wide[wide.transition == "Adenoma->Carcinoma"].iterrows():
        call = ("SUPPORTED (spread increases)" if r.permdisp_p_holm < .05 and r.permdisp_effect > 1 else
                "CONTRADICTION (significant, D < 1)" if r.permdisp_p_holm < .05 else "not significant (no support)")
        L.append(f"  {r.cohort:<9} {r.version:<10} D = {r.permdisp_effect:.2f}, PERMDISP Holm p {r.permdisp_p_holm:.3f} -> {call}"
                 f"   [TDA raw Holm p {r.tda_p_holm_raw:.3f}]")

    L += ["", "5.2 SHAPE RULE (each TDA shape-only significant comparison; 'not captured' only if PERMDISP,",
          "    stage separation and MST shape-only are ALL non-significant; PERMANOVA reported only):"]
    shp = wide[wide.tda_p_holm_shape < .05]
    for _, r in shp.iterrows():
        hits = [n for n, c in (("PERMDISP", "permdisp_p_holm"), ("stage separation", "stagesep_p_holm"),
                               ("MST shape-only", "mst_shape_p_holm")) if r[c] < .05]
        L.append(f"  {r.cohort} {r.version} {r.transition} (TDA shape Holm p {r.tda_p_holm_shape:.3f}): "
                 + ("NOT CAPTURED by conventional methods" if not hits else "accompanied by " + ", ".join(hits))
                 + f"   [PERMANOVA Holm p {r.permanova_p_holm:.3f}]")
    if shp.empty:
        L.append("  (no TDA shape-only significant comparisons)")

    L += ["", "3.6 SIZE-MATCHED (secondary; medians over 20 subsamples of size m; median p is descriptive):",
          f"  {'cohort':<9} {'version':<10} {'transition':<21} {'PDISP all':>9} {'PDISP m':>8} {'D m':>5} {'%sig':>5} "
          f"{'PANOVA all':>10} {'PANOVA m':>9} {'%sig':>5}"]
    for _, r in wide.iterrows():
        L.append(f"  {r.cohort:<9} {r.version:<10} {r.transition:<21} {st(r.permdisp_p_holm):>9} {st(r.sm_permdisp_p_holm):>8} "
                 f"{r.sm_D:5.2f} {100 * r.sm_permdisp_frac_sig:4.0f}% {st(r.permanova_p_holm):>10} {st(r.sm_permanova_p_holm):>9} "
                 f"{100 * r.sm_permanova_frac_sig:4.0f}%")
    L.append("  Reading for comparisons where the conventional test is significant on all samples but TDA raw is not:")
    any_ = False
    for _, r in wide.iterrows():
        for n, ca, cm in (("PERMDISP", "permdisp_p_holm", "sm_permdisp_p_holm"),
                          ("PERMANOVA", "permanova_p_holm", "sm_permanova_p_holm")):
            if r[ca] < .05 and r.tda_p_holm_raw >= .05:
                any_ = True
                L.append(f"    {r.cohort} {r.version} {r.transition}, {n}: at size m {r[cm]:.3f} -> "
                         + ("still significant: DIFFERENT PROPERTY" if r[cm] < .05 else "not significant: POWER"))
    if not any_:
        L.append("    (none)")
    L += ["", "Outputs: benchmark_results.csv (wide, with TDA), benchmark_results_long.csv, size_matched_results.csv,",
          "concordance.csv, benchmark_overview.png, checkpoints/<unit>.json"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lab = [f"{r.cohort[3:]} {r.version[:3]} {r.transition[0]}→{r.transition.split('>')[1][0]}" for _, r in wide.iterrows()]
    y = np.arange(len(wide))[::-1]
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
    for ax, col, xl, ref in ((axes[0], "permdisp_effect", "PERMDISP dispersion ratio D (later/earlier)", 1),
                             (axes[1], "stagesep_effect", "Stage separation ARI", 0),
                             (axes[2], "mst_shape_effect", "MST z (filled = shape-only, open = raw)", 0)):
        sigc = {"permdisp_effect": "permdisp_p_holm", "stagesep_effect": "stagesep_p_holm",
                "mst_shape_effect": "mst_shape_p_holm"}[col]
        cols = np.where(wide.tda_p_holm_raw < .05, "#eb6834", "0.6")
        ax.scatter(wide[col], y, c=cols, marker="o", s=np.where(wide[sigc] < .05, 70, 25), zorder=3)
        if col == "mst_shape_effect":
            ax.scatter(wide.mst_raw_effect, y, facecolors="none", edgecolors=cols, s=np.where(wide.mst_raw_p_holm < .05, 70, 25))
        ax.axvline(ref, color="0.5", lw=1); ax.set_xlabel(xl, fontsize=9)
    axes[0].set_yticks(y); axes[0].set_yticklabels(lab, fontsize=8)
    fig.suptitle("Step 13 benchmark: orange = TDA raw significant (Holm); large marker = conventional test significant (Holm)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "benchmark_overview.png"), dpi=150); plt.close(fig)


def run_benchmark():
    os.makedirs(CK, exist_ok=True)
    cc = real_crosscheck()
    for c in cc:
        print(("PASS " if c[2] else "FAIL ") + c[0] + ": " + c[1], flush=True)
    if not all(c[2] for c in cc):
        raise SystemExit("GUDHI cross-check on real subsample FAILED - stopping before any benchmark computation.")
    todo = [u for u in bench_units() if not os.path.exists(os.path.join(CK, bunit_name(*u) + ".json"))]
    # longest (MST) units first
    todo.sort(key=lambda u: 0 if u[4].startswith("mst") else 1)
    print(f"[step13] {len(todo)} units to run", flush=True)
    t0 = time.time()
    with Pool(WORKERS) as pool:
        for i, (n, msg) in enumerate(pool.imap_unordered(run_bunit, todo), 1):
            print(f"[step13] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    bench_summary()
    print(f"[step13] finished in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "validate":
        validate()
    elif cmd == "run":
        run_benchmark()
    elif cmd == "summary":
        bench_summary()
    else:
        raise SystemExit("Usage: step13_benchmark.py validate | run | summary")
