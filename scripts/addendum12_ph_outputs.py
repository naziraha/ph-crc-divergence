"""
Addendum 12 (73b4e71; D25 4018ea9) - quantitative persistence outputs.

Test (per transition / per validation replicate): groups aligned (tda_crc_pilot.align_groups); k = 20 subsamples of
size m per group; VR diagrams (tda_crc_pilot.persistence, H0 and H1, infinite bars removed); statistics computed from
the same draws (D25b):
  landscape distance (H0, H1): L2 distance between mean landscapes (gudhi.representations.Landscape, 5 landscapes,
                                100 points, grid 0..max finite death of the pooled observed diagrams; D25c)
  barcode summaries: H0 total / max / entropy; H1 total / max / n features / entropy -> |mean_s - mean_t|
Patient-level permutation of the aligned labels (group sizes kept), k subsamples redrawn per permuted group;
one-sided p = (1 + #null >= obs) / (1 + n_perm).

Usage: python scripts/addendum12_ph_outputs.py validate
"""
import json
import os
import sys
import time
import zlib
from importlib import metadata
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import tda_crc_pilot as T  # noqa: E402

SEED, K, WORKERS = 42, 20, 6
NL, RES = 5, 100
OUT = os.path.join(ROOT, "results", "addendum12")
SUMMARIES = ["H0_total_persistence", "H0_max_persistence", "H0_entropy",
             "H1_total_persistence", "H1_max_persistence", "H1_n_features", "H1_entropy"]
STATS = ["landscape_H0", "landscape_H1"] + SUMMARIES


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def draw(emb, idx, m, normalize, rng):
    return [T.persistence(emb[rng.choice(idx, size=m, replace=False)], normalize) for _ in range(K)]


def landscape_vecs(diags, dim, rmax):
    from gudhi.representations import Landscape
    L = Landscape(num_landscapes=NL, resolution=RES, sample_range=[0.0, rmax])
    out = []
    for d in diags:
        x = d[dim]
        out.append(np.zeros(NL * RES) if len(x) == 0 else L.fit_transform([x])[0])
    return np.array(out)


def group_stats(dA, dB, rmax):
    s = {}
    for dim in (0, 1):
        a, b = landscape_vecs(dA, dim, rmax[dim]), landscape_vecs(dB, dim, rmax[dim])
        s[f"landscape_H{dim}"] = float(np.linalg.norm(a.mean(0) - b.mean(0)))
    sa = pd.DataFrame([T.summaries(d) for d in dA])
    sb = pd.DataFrame([T.summaries(d) for d in dB])
    for c in SUMMARIES:
        s[c] = float(abs(sa[c].mean() - sb[c].mean()))
    return s


def multi_test(emb, idx_s, idx_t, m, normalize, n_perm, rng):
    """All addendum 12 statistics with one shared patient-level permutation scheme."""
    emb = T.align_groups(emb, idx_s, idx_t, normalize)
    dA, dB = draw(emb, idx_s, m, normalize, rng), draw(emb, idx_t, m, normalize, rng)
    rmax = {}
    for dim in (0, 1):
        deaths = [d[dim][:, 1].max() for d in dA + dB if len(d[dim])]
        rmax[dim] = float(max(deaths)) if deaths else 1.0
    obs = group_stats(dA, dB, rmax)
    pool, ns = np.concatenate([idx_s, idx_t]), len(idx_s)
    null = {k: [] for k in STATS}
    for _ in range(n_perm):
        perm = rng.permutation(pool)
        pa, pb = draw(emb, perm[:ns], m, normalize, rng), draw(emb, perm[ns:], m, normalize, rng)
        st = group_stats(pa, pb, rmax)
        for k in STATS:
            null[k].append(st[k])
    res = {}
    for k in STATS:
        nl = np.array(null[k])
        res[k] = dict(obs=obs[k], p=float((1 + np.sum(nl >= obs[k])) / (1 + n_perm)),
                      z=float((obs[k] - nl.mean()) / (nl.std(ddof=1) + 1e-12)))
    return res


# ================================================================ section 6 validation
VAL_COND = [("N0", 0), ("SP", 1.25), ("CL", 1.0), ("LP", 1.0)]
NA_, NB_, M_VAL, NPERM_VAL, R_VAL = 53, 48, 31, 99, 20


def val_unit(u):
    cond, x, rep = u
    name = f"a12val__{cond}{x:g}__r{rep:02d}"
    ck = os.path.join(OUT, "validation_checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    from addendum5_part3_controls import insert_level1
    rng = rng_for(name)
    df = pd.read_csv(os.path.join(ROOT, "results", "step03", "embeddings", "GSE41258_main_unadjusted_g2000_pc10.csv"),
                     index_col=0)
    carc = df.index[df.stage == "Carcinoma"]
    assert len(carc) == 181
    X = df.loc[carc, [c for c in df.columns if c.startswith("PC")]].values.copy()
    perm = rng.permutation(len(carc))
    ia, ib = perm[:NA_], perm[NA_:NA_ + NB_]
    X[ib] = insert_level1(X[ib], cond, x, rng)
    emb = np.r_[X[ia], X[ib]]
    gA, gB = np.arange(NA_), np.arange(NA_, NA_ + NB_)
    rec = dict(condition=cond, rep=rep)
    for mode in ("shape", "raw"):
        r = multi_test(emb, gA, gB, M_VAL, mode == "shape", NPERM_VAL, rng)
        for k, v in r.items():
            rec[f"{mode}__{k}__p"] = v["p"]
            rec[f"{mode}__{k}__z"] = v["z"]
    json.dump(rec, open(ck, "w"))
    return rec


def validate():
    os.makedirs(os.path.join(OUT, "validation_checkpoints"), exist_ok=True)
    t0 = time.time()
    units = [(c, x, r) for c, x in VAL_COND for r in range(R_VAL)]
    with Pool(WORKERS) as pool:
        recs = []
        for i, rec in enumerate(pool.imap_unordered(val_unit, units), 1):
            recs.append(rec)
            if i % 10 == 0:
                print(f"[a12 val] {i}/{len(units)} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    d = pd.DataFrame(recs)
    d.to_csv(os.path.join(OUT, "validation_replicates.csv"), index=False)
    det = lambda cond, col: int((d[d.condition == cond][col + "__p"] <= 0.05).sum())
    crit = []
    for mode in ("shape", "raw"):
        for dim in ("H0", "H1"):
            v = det("N0", f"{mode}__landscape_{dim}")
            crit.append((f"N0 landscape {dim} {mode}", "<= 3/20", v, v <= 3))
    for dim in ("H0", "H1"):
        v = det("SP", f"shape__landscape_{dim}")
        crit.append((f"SP landscape {dim} shape-only (ignore pure spread)", "<= 3/20", v, v <= 3))
    v = det("SP", "raw__landscape_H0")
    crit.append(("SP landscape H0 raw (D25a)", ">= 16/20", v, v >= 16))
    v = det("CL", "shape__landscape_H0")
    crit.append(("CL landscape H0 shape-only", ">= 16/20", v, v >= 16))
    v = det("LP", "shape__landscape_H1")
    crit.append(("LP landscape H1 shape-only", ">= 16/20", v, v >= 16))
    for s in SUMMARIES:
        v = det("N0", f"shape__{s}")
        crit.append((f"N0 barcode summary {s} shape-only", "<= 3/20", v, v <= 3))
    cdf = pd.DataFrame(crit, columns=["check", "requirement", "detected", "pass"])
    cdf.to_csv(os.path.join(OUT, "validation_criteria.csv"), index=False)
    ok = bool(cdf["pass"].all())
    tab = []
    for cond, _ in VAL_COND:
        for mode in ("shape", "raw"):
            tab.append(f"  {cond:<3} {mode:<6} " + "  ".join(f"{k.replace('_persistence', '').replace('landscape_', 'L_')}:{det(cond, f'{mode}__{k}'):>2}"
                                                       for k in STATS))
    L = ["ADDENDUM 12 section 6 - validation of landscape and barcode-summary permutation tests (D25)", "=" * 72,
         "GSE41258 carcinomas split A 53 / B 48; insertion into B (addendum 5 level 1); k = 20, m = 31, aligned;",
         "99 permutations; 20 replicates per condition; detection = p <= 0.05; shared draws for all statistics (D25b).",
         f"Seed {SEED} per replicate. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "pandas")), ""]
    L += [f"  {'PASS' if r['pass'] else 'FAIL'}  {r.check:<52} {r.requirement:<9} {r.detected}/20" for _, r in cdf.iterrows()]
    L += ["", "All detection counts (/20) by condition, mode and statistic:"] + tab
    L += ["", f"OVERALL: {'ALL CRITERIA PASSED - real-data analysis may proceed' if ok else 'FAILED - do not analyse real data'}",
          f"Run time {(time.time() - t0) / 60:.1f} min."]
    open(os.path.join(OUT, "validation_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        validate()
    else:
        raise SystemExit("Only 'validate' is implemented at this stage (addendum 12 section 6).")
