"""
Step 5 for TCGA-COAD (addendum 10, 296eaed; D18-D20) - pre-specified Wasserstein permutation test.

Identical test to the GEO cohorts (tda_crc_pilot.permutation_test unchanged; VR H0+H1; exact W1; k = 20;
m = 80% of the smallest stage of the set; groups aligned; 999 permutations), on carcinoma stages with
Holm across the THREE transitions I->II, II->III, III->IV per analysis x version x mode.
Versions: unadjusted, adjusted (stage-blind). Tipping point = shape-only Holm p < 0.05 in BOTH versions.
Analyses: primary (2,000 genes / 10 PCs); sensitivity 5 / 20 PCs, 1,000 / 5,000 genes; noLowColon (D19b);
MSS only (addendum 10 section 4.2). Resumable; per-unit seeds SeedSequence(42, crc32(unit name)).

Usage: python scripts/step05_permutation_TCGA.py [run|summary]
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

SEED, N_PERM, K_DRAWS, DIMS, WORKERS = 42, 999, 20, (0, 1), 6
GSE = "TCGA_COAD"
STAGES = ["Stage I", "Stage II", "Stage III", "Stage IV"]
TRANS = list(zip(STAGES[:-1], STAGES[1:]))
VERSIONS = ("unadjusted", "adjusted")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "results", "step05", "TCGA")
CK = os.path.join(OUT, "checkpoints")
ANALYSES = {"primary": ("main", "g2000_pc10"), "S5_pc5": ("main", "g2000_pc5"), "S6_pc20": ("main", "g2000_pc20"),
            "S7_genes1000": ("main", "g1000_pc10"), "S8_genes5000": ("main", "g5000_pc10"),
            "noLowColon": ("main_noLowColon", "g2000_pc10"), "MSS": ("main_MSS", "g2000_pc10")}
EXPECTED = {"main": [45, 110, 80, 39], "main_noLowColon": [43, 101, 75, 37], "main_MSS": [34, 79, 72, 37]}
ROMAN = {"Stage I": "I", "Stage II": "II", "Stage III": "III", "Stage IV": "IV"}


def unit_name(a, v, s, t, mode):
    return f"TCGA__{a}__{v[:5]}__{ROMAN[s]}-{ROMAN[t]}__{'shape' if mode == 'norm' else 'raw'}"


def units():
    return [(a, v, s, t, mode) for a in ANALYSES for v in VERSIONS for s, t in TRANS for mode in ("norm", "raw")]


def run_unit(u):
    a, v, s, t, mode = u
    name = unit_name(*u)
    path = os.path.join(CK, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    aset, suffix = ANALYSES[a]
    df = pd.read_csv(os.path.join(EMB, f"{GSE}_{aset}_{v}_{suffix}.csv"), index_col=0)
    cnt = [int((df.stage == x).sum()) for x in STAGES]
    if cnt != EXPECTED[aset] or len(df) != sum(cnt):
        raise SystemExit(f"Stage counts changed: {name} {cnt}")
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    st = df.stage.values
    m = max(10, int(0.8 * min(cnt)))
    i_s, i_t = np.where(st == s)[0], np.where(st == t)[0]
    Tobs, p, z, null = T.permutation_test(emb, i_s, i_t, m, K_DRAWS, DIMS, mode == "norm", N_PERM, rng, align=True)
    rec = dict(unit=name, cohort=GSE, analysis=a, set=aset, version=v, transition=f"{s}->{t}",
               mode="shape-only" if mode == "norm" else "raw", n_from=int(len(i_s)), n_to=int(len(i_t)),
               subsample_size=int(m), T_obs=float(Tobs), p=float(p), z=float(z), null=[float(x) for x in null],
               seconds=round(time.time() - t0, 1), finished=time.strftime("%Y-%m-%dT%H:%M:%S"))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"done in {rec['seconds']:.0f} s (p={p:.4f}, z={z:.1f})"


def summary():
    rows = []
    for f in sorted(os.listdir(CK)):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(CK, f)))
            r.pop("null")
            rows.append(r)
    df = pd.DataFrame(rows)
    df["p_holm"] = np.nan
    for _, g in df.groupby(["analysis", "version", "mode"]):
        if len(g) == 3:
            df.loc[g.index, "p_holm"] = T.holm(g["p"].values)
    order = {f"{s}->{t}": i for i, (s, t) in enumerate(TRANS)}
    df["_a"] = df.analysis.map({a: i for i, a in enumerate(ANALYSES)})
    df["_t"] = df.transition.map(order)
    df = df.sort_values(["_a", "version", "mode", "_t"], ascending=[True, False, False, True]).drop(columns=["_a", "_t"])
    df.to_csv(os.path.join(OUT, "permutation_results_all.csv"), index=False)
    L = ["STEP 5 - Pre-specified permutation test (TCGA-COAD, carcinoma stages)", "=" * 66,
         "Pre-specification: addendum 10 (296eaed); D18, D19, D20.",
         f"Seed {SEED} (per-unit SeedSequence, resumable); " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "scipy")),
         f"Settings: VR H0+H1, exact W1, k = {K_DRAWS}, m = 80% of smallest stage, aligned, {N_PERM} permutations,",
         "  Holm across the 3 transitions per analysis x version x mode.",
         "Samples per stage (I/II/III/IV): " + "; ".join(f"{k} {'/'.join(map(str, v))}" for k, v in EXPECTED.items()),
         f"Units finished: {len(df)}/{len(units())}", ""]
    for a in ANALYSES:
        sub = df[df.analysis == a]
        if sub.empty:
            continue
        L.append(f"--- {a} ({' / '.join(ANALYSES[a])}) ---")
        for _, r in sub.iterrows():
            L.append(f"  {r.version:<11} {r['mode']:<10} {r.transition:<22} n={r.n_from}/{r.n_to} m={r.subsample_size}"
                     f"  T={r.T_obs:9.4f}  z={r.z:6.2f}  p={r.p:.4f}  p_Holm {r.p_holm:.4f}")
        if len(sub) == 12:
            L.append("  TIPPING-POINT DECISION (shape-only Holm p < 0.05 in BOTH unadjusted and adjusted):")
            for s, t in TRANS:
                tr = f"{s}->{t}"
                get = lambda v, mo: sub[(sub.version == v) & (sub["mode"] == mo) & (sub.transition == tr)].p_holm.iloc[0]
                pu, pa = get("unadjusted", "shape-only"), get("adjusted", "shape-only")
                ru, ra = get("unadjusted", "raw"), get("adjusted", "raw")
                call = ("TIPPING POINT" if pu < .05 and pa < .05 else "significant in unadjusted only" if pu < .05
                        else "significant in adjusted only" if pa < .05 else "not significant")
                L.append(f"    {tr:<22} shape-only: unadjusted {pu:.4f}, adjusted {pa:.4f} -> {call}"
                         f"   | raw: unadjusted {ru:.4f}, adjusted {ra:.4f}")
        L.append("")
    L += ["Outputs: results/step05/TCGA/permutation_results_all.csv, checkpoints/<unit>.json, step05_primary.png"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))
    prim = df[df.analysis == "primary"]
    if len(prim) == 12:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
        for ax, mode in zip(axes, ("shape-only", "raw")):
            for k, (v, col) in enumerate(zip(VERSIONS, ("#2a78d6", "#eb6834"))):
                g = prim[(prim.version == v) & (prim["mode"] == mode)].set_index("transition").loc[[f"{s}->{t}" for s, t in TRANS]]
                x = np.arange(3) + (k - 0.5) * 0.35
                ax.bar(x, g.z.values, width=0.33, color=col, label=v)
                for xi, zz, ph in zip(x, g.z.values, g.p_holm.values):
                    ax.text(xi, max(zz, 0) + 0.3, f"p={ph:.3f}", ha="center", fontsize=7)
            ax.axhline(0, color="0.4", lw=1); ax.axhline(1.96, color="0.5", ls="--", lw=1)
            ax.set_xticks(range(3)); ax.set_xticklabels(["I→II", "II→III", "III→IV"])
            ax.set_title(mode); ax.legend(fontsize=8)
        axes[0].set_ylabel("z vs permutation null (labels: Holm p)")
        fig.suptitle("TCGA-COAD carcinoma stages: permutation test (999 permutations)")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "step05_primary.png"), dpi=150); plt.close(fig)


def main():
    os.makedirs(CK, exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "summary":
        summary(); return
    todo = [u for u in units() if not os.path.exists(os.path.join(CK, unit_name(*u) + ".json"))]
    print(f"[step5 TCGA] {len(todo)} units to run", flush=True)
    t0 = time.time()
    with Pool(WORKERS) as pool:
        for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
            print(f"[step5] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
