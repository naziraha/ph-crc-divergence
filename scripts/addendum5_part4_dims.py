"""
Addendum 5, part 4 (commit d5673a7) - H0-only and H1-only versions of the step 5 test (SECONDARY).

Step 5 primary test unchanged except dims = (0,) or (1,): 999 permutations, k = 20, m as in step 5,
aligned, primary embeddings (2,000 genes, 10 PCs). Cohorts: GSE41258 main (3 transitions), GSE20916
bulk and GSE37364 bulk (2 transitions); shape-only and raw; unadjusted and adjusted.
Holm across transitions within cohort x version x mode x dimension. Resumable; per-unit seeds.
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

SEED, K, NPERM, WORKERS = 42, 20, 999, 6
OUT = os.path.join(ROOT, "results", "addendum5", "part4_dims")
CKPT = os.path.join(OUT, "checkpoints")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
COH = {"GSE41258": dict(prefix="GSE41258_main", stages=list(T.STAGES),
                        expected={"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39},
                        step5=os.path.join(ROOT, "results", "step05", "permutation_results_all.csv")),
       "GSE20916": dict(prefix="GSE20916_bulk", stages=["Normal", "Adenoma", "Carcinoma"],
                        expected={"Normal": 24, "Adenoma": 45, "Carcinoma": 36},
                        step5=os.path.join(ROOT, "results", "step05", "GSE20916", "permutation_results_all.csv")),
       "GSE37364": dict(prefix="GSE37364_bulk", stages=["Normal", "Adenoma", "Carcinoma"],
                        expected={"Normal": 38, "Adenoma": 29, "Carcinoma": 27},
                        step5=os.path.join(ROOT, "results", "step05", "GSE37364", "permutation_results_all.csv"))}
VERSIONS, MODES, DIMSETS = ("unadjusted", "adjusted"), ("norm", "raw"), ((0,), (1,))


def units():
    out = []
    for c, cfg in COH.items():
        st = cfg["stages"]
        for v in VERSIONS:
            for mo in MODES:
                for dims in DIMSETS:
                    for s, t in zip(st[:-1], st[1:]):
                        out.append((c, v, mo, dims[0], s, t))
    return out


def unit_name(c, v, mo, d, s, t):
    return f"{c}__{v[:5]}__{'shape' if mo == 'norm' else 'raw'}__H{d}__{s[0]}-{t[0]}"


def run_unit(u):
    c, v, mo, d, s, t = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    cfg = COH[c]
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(os.path.join(EMB, f"{cfg['prefix']}_{v}_g2000_pc10.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != cfg["expected"]:
        raise SystemExit(f"Stage counts changed: {name}")
    emb = df[[x for x in df.columns if x.startswith("PC")]].values
    stv = df["stage"].values
    m = max(10, int(0.8 * min((stv == x).sum() for x in cfg["stages"])))
    idx_s, idx_t = np.where(stv == s)[0], np.where(stv == t)[0]
    Tobs, p, z, _ = T.permutation_test(emb, idx_s, idx_t, m, K, (d,), mo == "norm", NPERM, rng, align=True)
    rec = dict(unit=name, cohort=c, version=v, mode="shape-only" if mo == "norm" else "raw", dim=f"H{d}",
               transition=f"{s}->{t}", n_from=int(len(idx_s)), n_to=int(len(idx_t)), m=int(m),
               T_obs=float(Tobs), p=float(p), z=float(z), seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"z={z:.2f} p={p:.3f} ({rec['seconds']:.0f} s)"


def summary():
    r = pd.DataFrame([json.load(open(os.path.join(CKPT, f))) for f in os.listdir(CKPT) if f.endswith(".json")])
    r["p_holm"] = np.nan
    for _, g in r.groupby(["cohort", "version", "mode", "dim"]):
        if len(g) == len(COH[g.cohort.iloc[0]]["stages"]) - 1:
            r.loc[g.index, "p_holm"] = T.holm(g.p.values)
    comb = []
    for c, cfg in COH.items():
        s5 = pd.read_csv(cfg["step5"])
        s5 = s5[(s5.analysis == "primary") & (s5.version.isin(VERSIONS))]
        for _, x in s5.iterrows():
            row = dict(cohort=c, version=x.version, mode=x["mode"], transition=x.transition,
                       z_H0H1=x.z, p_holm_H0H1=x.p_holm)
            for d in ("H0", "H1"):
                q = r[(r.cohort == c) & (r.version == x.version) & (r["mode"] == x["mode"]) & (r.dim == d)
                      & (r.transition == x.transition)]
                row[f"z_{d}"] = float(q.z.iloc[0]) if len(q) else np.nan
                row[f"p_holm_{d}"] = float(q.p_holm.iloc[0]) if len(q) else np.nan
            comb.append(row)
    comb = pd.DataFrame(comb)
    comb.to_csv(os.path.join(OUT, "dims_comparison.csv"), index=False)
    r.to_csv(os.path.join(OUT, "dims_results_all.csv"), index=False)

    def cell(z, p):
        return f"{z:6.2f}{'*' if p < 0.05 else ' '} ({p:.3f})" if not np.isnan(z) else "   n/a"

    L = ["ADDENDUM 5, PART 4 - H0-only and H1-only versions of the step 5 test (SECONDARY)", "=" * 70,
         "Pre-specified: docs/analysis_prespecification_addendum5_checks.md part 4 (commit d5673a7).",
         "Step 5 test unchanged except homology dimensions; 999 permutations; Holm across transitions within",
         "cohort x version x mode x dimension. Values: z (Holm p); * = Holm p < 0.05. H0+H1 = step 5 primary.",
         f"Units finished: {len(r)}/{len(units())}", "",
         f"{'cohort':<9}{'version':<11}{'mode':<11}{'transition':<24}{'H0+H1 (step 5)':>18}{'H0 only':>18}{'H1 only':>18}"]
    for _, x in comb.iterrows():
        L.append(f"{x.cohort:<9}{x.version:<11}{x['mode']:<11}{x.transition:<24}{cell(x.z_H0H1, x.p_holm_H0H1):>18}"
                 f"{cell(x.z_H0, x.p_holm_H0):>18}{cell(x.z_H1, x.p_holm_H1):>18}")
    sig = comb[comb.p_holm_H0H1 < 0.05]
    L += ["", f"Of the {len(sig)} significant step 5 results (H0+H1): significant with H0 only in "
              f"{int((sig.p_holm_H0 < 0.05).sum())}, with H1 only in {int((sig.p_holm_H1 < 0.05).sum())}.",
          f"Results significant with H1 only but not with H0+H1: "
          f"{int(((comb.p_holm_H1 < 0.05) & (comb.p_holm_H0H1 >= 0.05)).sum())}; with H0 only but not H0+H1: "
          f"{int(((comb.p_holm_H0 < 0.05) & (comb.p_holm_H0H1 >= 0.05)).sum())}.",
          f"Correlation of z across all {len(comb)} results: H0+H1 vs H0 r = {comb.z_H0H1.corr(comb.z_H0):.2f}; "
          f"H0+H1 vs H1 r = {comb.z_H0H1.corr(comb.z_H1):.2f}",
          "", "Outputs: results/addendum5/part4_dims/dims_comparison.csv, dims_results_all.csv, checkpoints/"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
        todo.sort(key=lambda u: u[0] != "GSE41258")
        print(f"[P4] {len(todo)} units to run", flush=True)
        t0 = time.time()
        with Pool(WORKERS) as pool:
            for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
                print(f"[P4] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
