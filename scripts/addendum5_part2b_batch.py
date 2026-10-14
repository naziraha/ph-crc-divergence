"""
Addendum 5, part 2b (commit 0a3a86a) - batch-only comparison in GSE41258.

Normals scanned spring 2003 (2003-01-01..2003-06-30) vs late 2005-2006 (2005-09-01..2006-12-31),
corrected scan dates (D9), primary step 3 embeddings of set main (unadjusted, adjusted), raw and
shape-only; frozen step 5 test (k = 20, aligned, 999 permutations), m = min(31, int(0.8 * 21)) = 16;
Holm across the 4 tests. Size-matched Normal -> Adenoma reference: 20 random draws of 21 normals and
30 adenomas, m = 16, 999 permutations. Resumable checkpoints; per-unit seeds.
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
import tda_crc_pilot as T  # noqa: E402

SEED, K, NPERM, DIMS, WORKERS, NDRAW = 42, 20, 999, (0, 1), 6, 20
OUT = os.path.join(ROOT, "results", "addendum5", "part2b_batch")
CKPT = os.path.join(OUT, "checkpoints")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
DATES = os.path.join(ROOT, "results", "addendum5", "part2_scandates", "scan_dates_main_corrected.csv")
EXPECTED = {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39}
VERSIONS, MODES = ("unadjusted", "adjusted"), ("norm", "raw")


def groups():
    d = pd.read_csv(DATES, index_col=0, parse_dates=["date"])
    n = d[d.stage == "Normal"]
    early = n.index[(n.date >= "2003-01-01") & (n.date <= "2003-06-30")]
    late = n.index[(n.date >= "2005-09-01") & (n.date <= "2006-12-31")]
    return list(early), list(late)


def units():
    u = [("batch", v, mo, 0) for v in VERSIONS for mo in MODES]
    u += [("matchedNA", v, mo, i) for v in VERSIONS for mo in MODES for i in range(NDRAW)]
    return u


def unit_name(kind, v, mo, i):
    return f"{kind}__{v[:5]}__{'shape' if mo == 'norm' else 'raw'}__d{i}"


def run_unit(u):
    kind, v, mo, i = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(os.path.join(EMB, f"GSE41258_main_{v}_g2000_pc10.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != EXPECTED:
        raise SystemExit("Stage counts changed")
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    pos = {s: i_ for i_, s in enumerate(df.index)}
    early, late = groups()
    if kind == "batch":
        g1, g2 = np.array([pos[s] for s in early]), np.array([pos[s] for s in late])
    else:
        st = df["stage"].values
        g1 = rng.choice(np.where(st == "Normal")[0], size=len(early), replace=False)
        g2 = rng.choice(np.where(st == "Adenoma")[0], size=len(late), replace=False)
    m = min(31, int(0.8 * min(len(g1), len(g2))))
    Tobs, p, z, _ = T.permutation_test(emb, g1, g2, m, K, DIMS, mo == "norm", NPERM, rng, align=True)
    rec = dict(unit=name, kind=kind, version=v, mode="shape-only" if mo == "norm" else "raw", draw=i,
               n1=int(len(g1)), n2=int(len(g2)), m=int(m), T_obs=float(Tobs), p=float(p), z=float(z),
               seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"z={z:.2f} p={p:.3f} ({rec['seconds']:.0f} s)"


def summary():
    r = pd.DataFrame([json.load(open(os.path.join(CKPT, f))) for f in os.listdir(CKPT) if f.endswith(".json")])
    b = r[r.kind == "batch"].copy()
    b["p_holm"] = T.holm(b.p.values) if len(b) == 4 else np.nan
    s5 = pd.read_csv(os.path.join(ROOT, "results", "step05", "permutation_results_all.csv"))
    s5 = s5[(s5.analysis == "primary") & (s5.transition == "Normal->Adenoma")]
    early, late = groups()
    L = ["ADDENDUM 5, PART 2b - Batch-only comparison, GSE41258 normals (descriptive)", "=" * 64,
         "Pre-specified: docs/analysis_prespecification_addendum5_checks.md, part 2b (commit 0a3a86a).",
         f"Normals scanned spring 2003: n = {len(early)}; late 2005-2006: n = {len(late)} (corrected dates, D9).",
         "Frozen step 5 test, 999 permutations, k = 20, aligned, m = 16; Holm across the 4 batch tests.",
         "Decision (fixed in advance): any batch Holm p < 0.05 -> batch can mimic a stage difference in this dataset.", "",
         f"{'version':<11}{'mode':<11}| {'BATCH z':>8} {'p':>6} {'Holm p':>7} | {'N->A size-matched z: median [min, max]':>40} "
         f"| {'step 5 N->A z (53 vs 48, m=31)':>30}"]
    rows = []
    for v in VERSIONS:
        for mode in ("shape-only", "raw"):
            bb = b[(b.version == v) & (b["mode"] == mode)]
            mm = r[(r.kind == "matchedNA") & (r.version == v) & (r["mode"] == mode)]
            ss = s5[(s5.version == v) & (s5["mode"] == mode)]
            bz, bp, bh = (float(bb.z.iloc[0]), float(bb.p.iloc[0]), float(bb.p_holm.iloc[0])) if len(bb) else (np.nan,) * 3
            rows.append(dict(version=v, mode=mode, batch_z=bz, batch_p=bp, batch_p_holm=bh,
                             matched_NA_median_z=mm.z.median(), matched_NA_min_z=mm.z.min(),
                             matched_NA_max_z=mm.z.max(), matched_NA_n=len(mm),
                             matched_NA_frac_p05=float((mm.p <= 0.05).mean()) if len(mm) else np.nan,
                             step5_NA_z=float(ss.z.iloc[0]), step5_NA_p_holm=float(ss.p_holm.iloc[0])))
            L.append(f"{v:<11}{mode:<11}| {bz:8.2f} {bp:6.3f} {bh:7.3f} | {mm.z.median():9.2f} [{mm.z.min():5.2f}, {mm.z.max():5.2f}]"
                     f" (p<=0.05 in {int((mm.p <= 0.05).sum())}/{len(mm)}) | {float(ss.z.iloc[0]):12.2f} (Holm p {float(ss.p_holm.iloc[0]):.3f})")
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, "batch_vs_NA.csv"), index=False)
    sig = bool((b.p_holm < 0.05).any())
    L += ["", "DECISION: " + ("at least one batch-only test has Holm p < 0.05 -> BATCH CAN MIMIC A STAGE DIFFERENCE in "
                              "GSE41258; Normal->Adenoma results must be treated as possibly batch-driven."
                              if sig else "no batch-only test has Holm p < 0.05 -> batch did not produce a significant "
                              "topological difference; compare its z with the size-matched N->A z above."),
          "", "Outputs: results/addendum5/part2b_batch/batch_vs_NA.csv, checkpoints/"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    e, l = groups()
    print(f"[2b] early n={len(e)}, late n={len(l)}", flush=True)
    todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
    with Pool(WORKERS) as pool:
        for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
            print(f"[2b] {i}/{len(todo)} {n}: {msg}", flush=True)
    summary()


if __name__ == "__main__":
    main()
