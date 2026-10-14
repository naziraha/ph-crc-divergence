"""
Addendum 5, part 3 (commit d5673a7) - positive and negative controls on real GSE41258 data.

Base: the 181 carcinomas of set main; each replicate splits them at random into A (53) and B (48);
insertions into B only. Frozen step 5 test (k = 20, m = 31, aligned) with 99 permutations;
shape-only and raw; 50 replicates per condition (chunks of 10, resumable).
Level 1 (10-PC scores of GSE41258_main_unadjusted_g2000_pc10):
  N0 none | SP spread s in {1.05,1.10,1.25,1.50} | CL cluster split delta in {.25,.5,1} x MPD_B |
  LP loop w in {.25,.5,1}: random 2-D plane, w*annulus + (1-w)*original, R = 1.5 r_B, noise 0.1 R
Level 2 (step 1 log2 expression of A u B; 200 random genes passing the 25th-percentile filter;
  insertion into B; then tda_crc_pilot.preprocess (2,000 genes, z-score, 10 PCs) on A u B):
  N0 | SP s = 1 + a | CL +/- a SD_g L1_g | LP a SD_g (cos t L1_g + sin t L2_g); a in {.25,.5,1}
Acceptance criteria C1-C6 as in the addendum (detection at p <= 0.05, 50 replicates).
"""
import contextlib
import io
import json
import os
import sys
import time
import zlib
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

SEED, K, NPERM, DIMS, M, WORKERS, R, CHUNK = 42, 20, 99, (0, 1), 31, 6, 50, 10
NA_, NB_ = 53, 48
OUT = os.path.join(ROOT, "results", "addendum5", "part3_controls")
CKPT = os.path.join(OUT, "checkpoints")
EMB_FILE = os.path.join(ROOT, "results", "step03", "embeddings", "GSE41258_main_unadjusted_g2000_pc10.csv")
EXPR_FILE = os.path.join(ROOT, "data", "processed", "GSE41258_main_expr_genes.csv")
L1 = [("N0", 0)] + [("SP", s) for s in (1.05, 1.10, 1.25, 1.50)] + [("CL", d) for d in (0.25, 0.5, 1.0)] + \
     [("LP", w) for w in (0.25, 0.5, 1.0)]
L2 = [("N0", 0)] + [(c, a) for c in ("SP", "CL", "LP") for a in (0.25, 0.5, 1.0)]


def units():
    u = [(1, c, x, mo, ch) for c, x in L1 for mo in ("norm", "raw") for ch in range(R // CHUNK)]
    u += [(2, c, x, mo, ch) for c, x in L2 for mo in ("norm", "raw") for ch in range(R // CHUNK)]
    return u


def unit_name(lv, c, x, mo, ch):
    return f"L{lv}__{c}{x:g}__{'shape' if mo == 'norm' else 'raw'}__c{ch}"


def insert_level1(XB, cond, x, rng):
    XB = XB.copy()
    c = XB.mean(0)
    if cond == "SP":
        return c + x * (XB - c)
    if cond == "CL":
        u = rng.normal(size=XB.shape[1]); u /= np.linalg.norm(u)
        delta = x * np.median(pdist(XB))
        half = rng.permutation(len(XB))[: len(XB) // 2]
        sign = -np.ones(len(XB)); sign[half] = 1.0
        return XB + (sign[:, None] * delta) * u[None, :]
    if cond == "LP":
        Q, _ = np.linalg.qr(rng.normal(size=(XB.shape[1], 2)))
        u, v = Q[:, 0], Q[:, 1]
        D = XB - c
        a, b = D @ u, D @ v
        Rr = 1.5 * np.median(np.linalg.norm(D, axis=1))
        theta = 2 * np.pi * rng.permutation(len(XB)) / len(XB)
        rad = Rr + rng.normal(scale=0.1 * Rr, size=len(XB))
        a2 = x * rad * np.cos(theta) + (1 - x) * a
        b2 = x * rad * np.sin(theta) + (1 - x) * b
        return XB + np.outer(a2 - a, u) + np.outer(b2 - b, v)
    return XB


def insert_level2(E, colsB, cond, a, rng):
    """E: genes x samples (A u B), modified in place for 200 random filtered genes, columns of B."""
    mean = E.mean(axis=1)
    pool = np.where(mean > np.quantile(mean, 0.25))[0]
    genes = rng.choice(pool, size=200, replace=False)
    XB = E[np.ix_(genes, colsB)]
    sd = XB.std(axis=1, ddof=1)
    nB = len(colsB)
    if cond == "SP":
        mB = XB.mean(axis=1, keepdims=True)
        XB = mB + (1 + a) * (XB - mB)
    elif cond == "CL":
        L1g = rng.normal(size=len(genes))
        half = rng.permutation(nB)[: nB // 2]
        sign = -np.ones(nB); sign[half] = 1.0
        XB = XB + a * (sd * L1g)[:, None] * sign[None, :]
    elif cond == "LP":
        L1g, L2g = rng.normal(size=len(genes)), rng.normal(size=len(genes))
        theta = 2 * np.pi * rng.permutation(nB) / nB
        XB = XB + a * sd[:, None] * (np.outer(L1g, np.cos(theta)) + np.outer(L2g, np.sin(theta)))
    E[np.ix_(genes, colsB)] = XB
    return E


def run_unit(u):
    lv, cond, x, mo, ch = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    emb_df = pd.read_csv(EMB_FILE, index_col=0)
    carc = emb_df.index[emb_df.stage == "Carcinoma"]
    assert len(carc) == 181
    if lv == 2:
        expr = pd.read_csv(EXPR_FILE, index_col=0)[carc]
    out = []
    for _ in range(CHUNK):
        perm = rng.permutation(len(carc))
        ia, ib = perm[:NA_], perm[NA_:NA_ + NB_]
        if lv == 1:
            X = emb_df.loc[carc, [c for c in emb_df.columns if c.startswith("PC")]].values.copy()
            X[ib] = insert_level1(X[ib], cond, x, rng)
            emb = np.r_[X[ia], X[ib]]
        else:
            E = expr.iloc[:, np.r_[ia, ib]].values.astype(float).copy()
            E = insert_level2(E, np.arange(NA_, NA_ + NB_), cond, x, rng)
            Edf = pd.DataFrame(E, index=expr.index)
            with contextlib.redirect_stdout(io.StringIO()):
                e, _, _ = T.preprocess(Edf, 2000, 10)
            emb = e.values
        gA, gB = np.arange(NA_), np.arange(NA_, NA_ + NB_)
        _, p, z, _ = T.permutation_test(emb, gA, gB, M, K, DIMS, mo == "norm", NPERM, rng, align=True)
        out.append((float(p), float(z)))
    rec = dict(unit=name, level=lv, condition=cond, level_value=x, mode="shape-only" if mo == "norm" else "raw",
               chunk=ch, p=[o[0] for o in out], z=[o[1] for o in out], seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"{sum(p <= 0.05 for p, _ in out)}/{CHUNK} detected ({rec['seconds']:.0f} s)"


LABEL = {"N0": "negative control", "SP": "spread", "CL": "cluster split", "LP": "loop"}


def summary():
    rows = [json.load(open(os.path.join(CKPT, f))) for f in os.listdir(CKPT) if f.endswith(".json")]
    reps = pd.DataFrame([dict(level=r["level"], condition=r["condition"], value=r["level_value"], mode=r["mode"],
                              p=p, z=z) for r in rows for p, z in zip(r["p"], r["z"])])
    tab = reps.groupby(["level", "condition", "value", "mode"]).agg(
        n=("p", "size"), det05=("p", lambda s: int((s <= 0.05).sum())),
        det_bonf=("p", lambda s: int((s <= 0.05 / 3).sum())), median_z=("z", "median")).reset_index()
    order = {"N0": 0, "SP": 1, "CL": 2, "LP": 3}
    tab = tab.sort_values(["level", "condition", "value", "mode"], key=lambda c: c.map(order) if c.name == "condition" else c)
    tab.to_csv(os.path.join(OUT, "controls_table.csv"), index=False)

    def det(lv, c, v, mode):
        q = tab[(tab.level == lv) & (tab.condition == c) & (np.isclose(tab.value, v)) & (tab["mode"] == mode)]
        return (int(q.det05.iloc[0]), int(q.n.iloc[0])) if len(q) else (None, None)

    crit = []
    c1 = [det(lv, "N0", 0, mo) for lv in (1, 2) for mo in ("shape-only", "raw")]
    crit.append(("C1", "negative control <= 7/50 (levels 1 and 2, both modes)", c1, all(d is not None and d <= 7 for d, _ in c1)))
    c2 = [det(1, "SP", 1.25, "raw")]
    crit.append(("C2", "raw detects spread s = 1.25 >= 40/50", c2, all(d is not None and d >= 40 for d, _ in c2)))
    c3 = [det(1, "SP", s, "shape-only") for s in (1.05, 1.10, 1.25, 1.50)]
    crit.append(("C3", "shape-only ignores pure spread: every s <= 7/50", c3, all(d is not None and d <= 7 for d, _ in c3)))
    c4 = [det(1, "CL", 1.0, "shape-only")]
    crit.append(("C4", "shape-only detects cluster split delta = 1 x MPD >= 40/50", c4, all(d is not None and d >= 40 for d, _ in c4)))
    c5 = [det(1, "LP", 1.0, "shape-only")]
    crit.append(("C5", "shape-only detects full loop w = 1 >= 40/50", c5, all(d is not None and d >= 40 for d, _ in c5)))
    c6 = [det(2, "CL", 1.0, "shape-only"), det(2, "LP", 1.0, "shape-only")]
    crit.append(("C6", "level 2: a = 1 cluster and loop, shape-only, >= 40/50", c6, all(d is not None and d >= 40 for d, _ in c6)))

    L = ["ADDENDUM 5, PART 3 - Positive and negative controls (real GSE41258 carcinomas)", "=" * 70,
         "Pre-specified: docs/analysis_prespecification_addendum5_checks.md part 3 (commit d5673a7).",
         "Random carcinoma split A = 53 / B = 48; insertion into B; step 5 test (k = 20, m = 31, aligned), 99 permutations.",
         f"Replicates finished: {len(reps)}/{len(units()) * CHUNK}", "",
         f"{'lvl':<4}{'condition':<18}{'value':>7}  {'mode':<11}{'det p<=.05':>11}{'det p<=.0167':>13}{'median z':>10}"]
    for _, r in tab.iterrows():
        L.append(f"{r.level:<4}{LABEL[r.condition]:<18}{r.value:>7g}  {r['mode']:<11}{r.det05:>6}/{r.n:<4}"
                 f"{r.det_bonf:>8}/{r.n:<4}{r.median_z:>10.2f}")
    L += ["", "ACCEPTANCE CRITERIA (fixed in advance)"]
    for cid, txt, vals, ok in crit:
        L.append(f"  {cid}: {txt}: " + ", ".join(f"{d}/{n}" for d, n in vals) + f"  -> {'PASS' if ok else 'FAIL'}")
    L += ["", "Outputs: results/addendum5/part3_controls/controls_table.csv, checkpoints/"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
        todo.sort(key=lambda u: u[0])
        print(f"[P3] {len(todo)} chunks to run", flush=True)
        t0 = time.time()
        with Pool(WORKERS) as pool:
            for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
                print(f"[P3] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
