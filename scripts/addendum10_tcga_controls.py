"""
Addendum 10 section 6 positive-control check on TCGA-COAD, as made concrete in D20 (92f9de4).

Base: the 110 Stage II tumours of TCGA_COAD_main_unadjusted_g2000_pc10; each replicate splits them at random into
A (55) and B (55); insertion into B only (addendum5_part3_controls.insert_level1, unchanged). Frozen step 5 test
(k = 20, m = 31, aligned), 99 permutations, shape-only and raw, 20 replicates per condition and mode (2 chunks of 10).
Conditions: N0 none; SP s = 1.25; CL delta = 1 x MPD_B; LP w = 1.
"""
import json
import os
import sys
import time
import zlib
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import tda_crc_pilot as T  # noqa: E402
from addendum5_part3_controls import insert_level1  # noqa: E402

SEED, K, NPERM, DIMS, M, WORKERS, R, CHUNK = 42, 20, 99, (0, 1), 31, 6, 20, 10
OUT = os.path.join(ROOT, "results", "addendum10", "controls")
CKPT = os.path.join(OUT, "checkpoints")
EMB_FILE = os.path.join(ROOT, "results", "step03", "embeddings", "TCGA_COAD_main_unadjusted_g2000_pc10.csv")
CONDS = [("N0", 0), ("SP", 1.25), ("CL", 1.0), ("LP", 1.0)]


def units():
    return [(c, x, mo, ch) for c, x in CONDS for mo in ("norm", "raw") for ch in range(R // CHUNK)]


def unit_name(c, x, mo, ch):
    return f"TCGActl__{c}{x:g}__{'shape' if mo == 'norm' else 'raw'}__c{ch}"


def run_unit(u):
    cond, x, mo, ch = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(EMB_FILE, index_col=0)
    s2 = df.index[df.stage == "Stage II"]
    assert len(s2) == 110
    Xall = df.loc[s2, [c for c in df.columns if c.startswith("PC")]].values
    out = []
    for _ in range(CHUNK):
        perm = rng.permutation(len(s2))
        ia, ib = perm[:55], perm[55:]
        X = Xall.copy()
        X[ib] = insert_level1(X[ib], cond, x, rng)
        emb = np.r_[X[ia], X[ib]]
        _, p, z, _ = T.permutation_test(emb, np.arange(55), np.arange(55, 110), M, K, DIMS, mo == "norm", NPERM, rng,
                                        align=True)
        out.append((float(p), float(z)))
    rec = dict(unit=name, condition=cond, level_value=x, mode="shape-only" if mo == "norm" else "raw", chunk=ch,
               p=[o[0] for o in out], z=[o[1] for o in out], seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"{sum(p <= 0.05 for p, _ in out)}/{CHUNK} detected ({rec['seconds']:.0f} s)"


def summary():
    recs = [json.load(open(os.path.join(CKPT, f))) for f in sorted(os.listdir(CKPT)) if f.endswith(".json")]
    rows = []
    for r in recs:
        for p, z in zip(r["p"], r["z"]):
            rows.append(dict(condition=r["condition"], mode=r["mode"], p=p, z=z))
    d = pd.DataFrame(rows)
    tab = d.groupby(["condition", "mode"]).agg(n=("p", "size"), detected=("p", lambda v: int((v <= 0.05).sum())),
                                               median_z=("z", "median")).reset_index()
    tab.to_csv(os.path.join(OUT, "controls_table.csv"), index=False)
    det = lambda c, m: int(tab[(tab.condition == c) & (tab["mode"] == m)].detected.iloc[0])
    crit = [("N0 shape-only <= 3/20", det("N0", "shape-only"), det("N0", "shape-only") <= 3),
            ("N0 raw <= 3/20", det("N0", "raw"), det("N0", "raw") <= 3),
            ("SP raw >= 16/20", det("SP", "raw"), det("SP", "raw") >= 16),
            ("SP shape-only <= 3/20", det("SP", "shape-only"), det("SP", "shape-only") <= 3),
            ("CL shape-only >= 16/20", det("CL", "shape-only"), det("CL", "shape-only") >= 16),
            ("LP shape-only >= 16/20", det("LP", "shape-only"), det("LP", "shape-only") >= 16)]
    ok = all(c[2] for c in crit)
    L = ["ADDENDUM 10 section 6 - positive-control check on TCGA-COAD (D20, 92f9de4)", "=" * 64,
         "Stage II tumours (110) split at random A 55 / B 55; insertion into B (addendum 5 level 1 functions);",
         "step 5 test k = 20, m = 31, aligned, 99 permutations; 20 replicates per condition and mode; seed 42.", "",
         "  condition   mode         detected/20   median z"]
    L += [f"  {r.condition:<10}  {r['mode']:<11}  {r.detected:>4}/{r.n:<6}   {r.median_z:6.2f}" for _, r in tab.iterrows()]
    L += ["", "Acceptance:"] + [f"  {'PASS' if c[2] else 'FAIL'}  {c[0]:<26} {c[1]}" for c in crit]
    L.append(f"OVERALL: {'ALL CRITERIA PASSED' if ok else 'FAILED - report with limitation; PI decides before interpretation'}")
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    os.makedirs(CKPT, exist_ok=True)
    todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
    t0 = time.time()
    with Pool(WORKERS) as pool:
        for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
            print(f"[ctl] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()
