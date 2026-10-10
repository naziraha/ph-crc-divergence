"""
Addendum 10 (296eaed) section 4 - additional TCGA-COAD analyses, run one at a time.

  spread  : section 4.1 - PERMDISP per transition (bias-adjusted distances, as addendum 9 Amendment 1; 999
            permutations; Holm across the 3 transitions) with dispersion ratio D; and H0-only / H1-only versions
            of the step 5 test (as addendum 5 part 4), unadjusted and adjusted, shape-only and raw.

Usage: python scripts/addendum10_tcga_extra.py spread
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

SEED, NPERM, K, WORKERS = 42, 999, 20, 6
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "results", "addendum10")
STAGES = ["Stage I", "Stage II", "Stage III", "Stage IV"]
TRANS = list(zip(STAGES[:-1], STAGES[1:]))
ROMAN = {"Stage I": "I", "Stage II": "II", "Stage III": "III", "Stage IV": "IV"}
VERSIONS = ("unadjusted", "adjusted")
EXPECTED = [45, 110, 80, 39]


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def load(version):
    df = pd.read_csv(os.path.join(EMB, f"TCGA_COAD_main_{version}_g2000_pc10.csv"), index_col=0)
    if [int((df.stage == s).sum()) for s in STAGES] != EXPECTED:
        raise SystemExit("Stage counts changed")
    return df[[c for c in df.columns if c.startswith("PC")]].values, df.stage.values


# ---------------------------------------------------------------- 4.1 spread and H0 / H1
def dims_unit(u):
    v, mo, d, s, t = u
    name = f"TCGA41__{v[:5]}__{'shape' if mo == 'norm' else 'raw'}__H{d}__{ROMAN[s]}-{ROMAN[t]}"
    ck = os.path.join(OUT, "spread_H0H1", "checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    emb, st = load(v)
    m = int(0.8 * min(EXPECTED))
    i_s, i_t = np.where(st == s)[0], np.where(st == t)[0]
    Tobs, p, z, _ = T.permutation_test(emb, i_s, i_t, m, K, (d,), mo == "norm", NPERM, rng_for(name), align=True)
    rec = dict(unit=name, version=v, mode="shape-only" if mo == "norm" else "raw", dim=f"H{d}", transition=f"{s}->{t}",
               T_obs=float(Tobs), p=float(p), z=float(z))
    json.dump(rec, open(ck, "w"))
    return rec


def spread():
    from addendum9_heterogeneity import loo_dist, permdisp_dist   # bias-adjusted (addendum 9 Amendment 1)
    out = os.path.join(OUT, "spread_H0H1")
    os.makedirs(os.path.join(out, "checkpoints"), exist_ok=True)
    t0 = time.time()
    pd_rows = []
    for v in VERSIONS:
        emb, st = load(v)
        for s, t in TRANS:
            A, B = emb[st == s], emb[st == t]
            dA, dB = loo_dist(A, np.zeros(len(A))), loo_dist(B, np.zeros(len(B)))
            F, p, D = permdisp_dist(dA, dB, rng_for(f"TCGA41__permdisp__{v}__{ROMAN[s]}-{ROMAN[t]}"))
            pd_rows.append(dict(version=v, transition=f"{s}->{t}", n_from=len(A), n_to=len(B), F=F, p=p, D=D,
                                mean_dist_from=float(dA.mean()), mean_dist_to=float(dB.mean())))
    pdf = pd.DataFrame(pd_rows)
    pdf["p_holm"] = np.nan
    for _, g in pdf.groupby("version"):
        pdf.loc[g.index, "p_holm"] = T.holm(g.p.values)
    pdf.to_csv(os.path.join(out, "permdisp_by_stage.csv"), index=False)
    units = [(v, mo, d, s, t) for v in VERSIONS for mo in ("norm", "raw") for d in (0, 1) for s, t in TRANS]
    with Pool(WORKERS) as pool:
        recs = pool.map(dims_unit, units)
    dd = pd.DataFrame(recs)
    dd["p_holm"] = np.nan
    for _, g in dd.groupby(["version", "mode", "dim"]):
        dd.loc[g.index, "p_holm"] = T.holm(g.p.values)
    s5 = pd.read_csv(os.path.join(ROOT, "results", "step05", "TCGA", "permutation_results_all.csv"))
    s5 = s5[s5.analysis == "primary"].assign(dim="H0+H1")[["version", "mode", "dim", "transition", "z", "p", "p_holm"]]
    allr = pd.concat([s5, dd[["version", "mode", "dim", "transition", "z", "p", "p_holm"]]], ignore_index=True)
    allr.to_csv(os.path.join(out, "dims_results_all.csv"), index=False)
    L = ["ADDENDUM 10 section 4.1 - TCGA-COAD spread (PERMDISP) and H0-only / H1-only tests by stage", "=" * 72,
         "Embeddings: TCGA_COAD_main_<version>_g2000_pc10 (274 tumours; I/II/III/IV 45/110/80/39).",
         "PERMDISP: bias-adjusted distances sqrt(n/(n-1))*||x-c|| (addendum 9 Amendment 1), 999 permutations,",
         "  Holm across the 3 transitions per version; D = mean distance (later) / mean distance (earlier).",
         "H0/H1-only: step 5 test unchanged except dims (k = 20, m = 31, aligned, 999 permutations), Holm across the",
         "  3 transitions per version x mode x dimension (as addendum 5 part 4).",
         f"Seed {SEED} per unit. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "scipy")), "",
         "PERMDISP (spread):"]
    for _, r in pdf.iterrows():
        L.append(f"  {r.version:<11} {r.transition:<22} D = {r.D:.3f}  (mean distance {r.mean_dist_from:.1f} -> {r.mean_dist_to:.1f})"
                 f"  p {r.p:.3f}  Holm {r.p_holm:.3f}{'  *' if r.p_holm < .05 else ''}")
    L += ["", "z (Holm p) by homology dimension:  H0+H1 (step 5) | H0 only | H1 only"]
    for v in VERSIONS:
        for mo in ("shape-only", "raw"):
            for s, t in TRANS:
                tr = f"{s}->{t}"
                cell = []
                for d in ("H0+H1", "H0", "H1"):
                    r = allr[(allr.version == v) & (allr["mode"] == mo) & (allr.dim == d) & (allr.transition == tr)].iloc[0]
                    cell.append(f"{r.z:6.2f} ({r.p_holm:.3f}){'*' if r.p_holm < .05 else ' '}")
                L.append(f"  {v:<11} {mo:<10} {tr:<22} " + " | ".join(cell))
    nsig = int((dd.p_holm < .05).sum())
    L += ["", f"Significant (Holm p < 0.05): PERMDISP {int((pdf.p_holm < .05).sum())}/6; H0-only / H1-only {nsig}/24.",
          f"Run time {(time.time() - t0) / 60:.1f} min.",
          "Outputs: results/addendum10/spread_H0H1/permdisp_by_stage.csv, dims_results_all.csv, checkpoints/"]
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "spread":
        spread()
    else:
        raise SystemExit("Usage: addendum10_tcga_extra.py spread")
