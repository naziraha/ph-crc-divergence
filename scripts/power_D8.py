"""
Deviation D8 (commit 9d814f9) - sample-size sensitivity of the shape-only test (descriptive).

Real shape-only Normal -> Adenoma differences detected in step 5 (unadjusted): GSE41258 main (53/48)
and GSE20916 bulk (24/45). Groups reduced to f in {1, .75, .5} of their size without replacement;
m = min(cohort m, int(0.8 x smaller reduced group)); step 5 test unchanged (k = 20, aligned,
shape-only) with 99 permutations; 50 replicates per cell (chunks of 10, resumable).
Detection p <= 0.05/K (K = 3 GSE41258, 2 GSE20916); p <= 0.05 also reported.

Usage: python scripts/power_D8.py [run|summary]
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
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

SEED, K_DRAWS, N_PERM, DIMS, WORKERS = 42, 20, 99, (0, 1), 6
FRACS, R, CHUNK = (1.0, 0.75, 0.5), 50, 10
OUT = os.path.join(ROOT, "results", "power_D8")
CKPT = os.path.join(OUT, "checkpoints")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
COH = {"GSE41258": dict(file="GSE41258_main_unadjusted_g2000_pc10.csv", K=3, m=31, step5_z=6.90,
                        expected={"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39}),
       "GSE20916": dict(file="GSE20916_bulk_unadjusted_g2000_pc10.csv", K=2, m=19, step5_z=4.27,
                        expected={"Normal": 24, "Adenoma": 45, "Carcinoma": 36})}


def units():
    return [(c, f, ch) for c in COH for f in FRACS for ch in range(R // CHUNK)]


def unit_name(c, f, ch):
    return f"{c}__shape__N-A__f{int(f * 100)}__c{ch}"


def run_unit(u):
    c, f, ch = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    cfg = COH[c]
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(os.path.join(EMB, cfg["file"]), index_col=0)
    if df["stage"].value_counts().to_dict() != cfg["expected"]:
        raise SystemExit(f"Stage counts changed: {name}")
    emb = df[[x for x in df.columns if x.startswith("PC")]].values
    stv = df["stage"].values
    iN, iA = np.where(stv == "Normal")[0], np.where(stv == "Adenoma")[0]
    nN, nA = int(round(f * len(iN))), int(round(f * len(iA)))
    m = min(cfg["m"], int(0.8 * min(nN, nA)))
    out = []
    for _ in range(CHUNK):
        g1 = rng.choice(iN, size=nN, replace=False)
        g2 = rng.choice(iA, size=nA, replace=False)
        _, p, z, _ = T.permutation_test(emb, g1, g2, m, K_DRAWS, DIMS, True, N_PERM, rng, align=True)
        out.append((float(p), float(z)))
    rec = dict(unit=name, cohort=c, frac=f, n_N=nN, n_A=nA, m=m, p=[o[0] for o in out], z=[o[1] for o in out],
               seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"detected {sum(p <= 0.05 / cfg['K'] for p, _ in out)}/{CHUNK} in {rec['seconds']:.0f} s"


def summary():
    rows = [json.load(open(os.path.join(CKPT, f))) for f in sorted(os.listdir(CKPT)) if f.endswith(".json")]
    reps = pd.DataFrame([dict(cohort=r["cohort"], frac=r["frac"], n_N=r["n_N"], n_A=r["n_A"], m=r["m"], p=p, z=z)
                         for r in rows for p, z in zip(r["p"], r["z"])])
    reps["K"] = reps.cohort.map({c: COH[c]["K"] for c in COH})
    tab = reps.groupby(["cohort", "frac"]).agg(n_N=("n_N", "first"), n_A=("n_A", "first"), m=("m", "first"),
                                               replicates=("p", "size"),
                                               power_bonf=("p", lambda x: np.nan),
                                               median_z=("z", "median")).reset_index()
    tab["power_bonf"] = [float((reps[(reps.cohort == c) & (reps.frac == f)].p <= 0.05 / COH[c]["K"]).mean())
                         for c, f in zip(tab.cohort, tab.frac)]
    tab["power_05"] = [float((reps[(reps.cohort == c) & (reps.frac == f)].p <= 0.05).mean())
                       for c, f in zip(tab.cohort, tab.frac)]
    tab = tab.sort_values(["cohort", "frac"], ascending=[True, False])
    tab.to_csv(os.path.join(OUT, "power_table.csv"), index=False)
    L = ["DEVIATION D8 - Sample-size sensitivity of the shape-only test (DESCRIPTIVE)", "=" * 64,
         "docs/analysis_deviations.md D8 (commit 9d814f9), specified before computation.",
         "Real shape-only Normal -> Adenoma differences (unadjusted) detected in step 5; groups reduced without",
         "replacement; step 5 test unchanged except 99 permutations; detection p <= 0.05/K; 50 replicates per cell.",
         "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy")),
         f"Replicates finished: {len(reps)}/{len(COH) * len(FRACS) * R}", ""]
    for _, r in tab.iterrows():
        L.append(f"  {r.cohort} (step 5 z = {COH[r.cohort]['step5_z']:.2f})  groups {int(r.frac * 100):3d}%: "
                 f"n = {r.n_N}/{r.n_A}, m = {r.m}: detected {r.power_bonf * 100:3.0f}% (strict) / "
                 f"{r.power_05 * 100:3.0f}% (p <= 0.05); median z {r.median_z:.2f}")
    L += ["", "Reading: how reliably a shape difference of the size actually observed is found with fewer samples.",
          "It does not estimate power for smaller shape effects.",
          "Outputs: results/power_D8/power_table.csv, checkpoints/"]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
        todo.sort(key=lambda u: u[0] != "GSE41258")
        print(f"[D8] {len(todo)} chunks to run", flush=True)
        t0 = time.time()
        with Pool(WORKERS) as pool:
            for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
                print(f"[D8] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
