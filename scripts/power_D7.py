"""
Deviation D7 (commit a446a34) - empirical power of the step 5 test, descriptive.

Benchmark: each cohort's real Normal vs Adenoma samples (primary step 3 embedding, unadjusted and
adjusted). Effect scaled to alpha in {1, .75, .5, .25} by swapping k = round((1-alpha)/(1/nN + 1/nA))
random adenomas and normals between the groups (real samples only, each once, group sizes kept).
Each replicate = tda_crc_pilot.permutation_test unchanged (k = 20, cohort m, aligned) with 99
permutations. Detection: p <= 0.05/K (K = 3 GSE41258, 2 otherwise); p <= 0.05 also reported.
50 replicates per cohort x version x mode x alpha, run in chunks of 10 (resumable checkpoints).
Then every non-significant step 5 primary result is classified:
  'likely no meaningful effect' if Bonferroni-level detection at alpha = 50% >= 80%
  (same cohort, version, mode), else 'too few samples to tell'.

Usage: python scripts/power_D7.py [run|summary]
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
import tda_crc_pilot as T  # noqa: E402

SEED, K_DRAWS, N_PERM, DIMS, WORKERS = 42, 20, 99, (0, 1), 6
ALPHAS, R, CHUNK = (1.0, 0.75, 0.5, 0.25), 50, 10
OUT = os.path.join(ROOT, "results", "power_D7")
CKPT = os.path.join(OUT, "checkpoints")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
COH = {"GSE41258": dict(prefix="GSE41258_main", K=3, stages=list(T.STAGES),
                        expected={"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39},
                        step5=os.path.join(ROOT, "results", "step05", "permutation_results_all.csv")),
       "GSE20916": dict(prefix="GSE20916_bulk", K=2, stages=["Normal", "Adenoma", "Carcinoma"],
                        expected={"Normal": 24, "Adenoma": 45, "Carcinoma": 36},
                        step5=os.path.join(ROOT, "results", "step05", "GSE20916", "permutation_results_all.csv")),
       "GSE37364": dict(prefix="GSE37364_bulk", K=2, stages=["Normal", "Adenoma", "Carcinoma"],
                        expected={"Normal": 38, "Adenoma": 29, "Carcinoma": 27},
                        step5=os.path.join(ROOT, "results", "step05", "GSE37364", "permutation_results_all.csv"))}
VERSIONS, MODES = ("unadjusted", "adjusted"), ("norm", "raw")


def units():
    return [(c, v, mo, a, ch) for c in COH for v in VERSIONS for mo in MODES for a in ALPHAS
            for ch in range(R // CHUNK)]


def unit_name(c, v, mo, a, ch):
    return f"{c}__{v[:5]}__{'shape' if mo == 'norm' else 'raw'}__a{int(a * 100)}__c{ch}"


def run_unit(u):
    c, v, mo, a, ch = u
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
    iN, iA = np.where(stv == "Normal")[0], np.where(stv == "Adenoma")[0]
    k = int(round((1 - a) / (1 / len(iN) + 1 / len(iA))))
    ps = []
    for _ in range(CHUNK):
        sw_n = rng.choice(iN, size=k, replace=False) if k else np.array([], int)
        sw_a = rng.choice(iA, size=k, replace=False) if k else np.array([], int)
        g1 = np.concatenate([np.setdiff1d(iN, sw_n), sw_a])
        g2 = np.concatenate([np.setdiff1d(iA, sw_a), sw_n])
        _, p, z, _ = T.permutation_test(emb, g1, g2, m, K_DRAWS, DIMS, mo == "norm", N_PERM, rng, align=True)
        ps.append((float(p), float(z)))
    rec = dict(unit=name, cohort=c, version=v, mode="shape-only" if mo == "norm" else "raw", alpha=a, chunk=ch,
               n_N=int(len(iN)), n_A=int(len(iA)), swaps=k,
               composition_contrast=float(1 - k / len(iN) - k / len(iA)), m=int(m),
               p=[x[0] for x in ps], z=[x[1] for x in ps], seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"detected {sum(p <= 0.05 / cfg['K'] for p, _ in ps)}/{CHUNK} in {rec['seconds']:.0f} s"


def summary():
    rows = [json.load(open(os.path.join(CKPT, f))) for f in sorted(os.listdir(CKPT)) if f.endswith(".json")]
    reps = []
    for r in rows:
        for p, z in zip(r["p"], r["z"]):
            reps.append(dict(cohort=r["cohort"], version=r["version"], mode=r["mode"], alpha=r["alpha"],
                             swaps=r["swaps"], contrast=r["composition_contrast"], p=p, z=z))
    reps = pd.DataFrame(reps)
    reps["K"] = reps.cohort.map({c: COH[c]["K"] for c in COH})
    reps["det_bonf"] = reps.p <= 0.05 / reps.K
    reps["det_05"] = reps.p <= 0.05
    pw = reps.groupby(["cohort", "version", "mode", "alpha"]).agg(
        n=("p", "size"), swaps=("swaps", "first"), contrast=("contrast", "first"),
        power_bonf=("det_bonf", "mean"), power_05=("det_05", "mean"), median_z=("z", "median")).reset_index()
    pw["mc_se"] = np.sqrt(pw.power_bonf * (1 - pw.power_bonf) / pw.n)
    pw.to_csv(os.path.join(OUT, "power_table.csv"), index=False)

    cls = []
    for c, cfg in COH.items():
        s5 = pd.read_csv(cfg["step5"])
        s5 = s5[(s5.analysis == "primary") & (s5.version.isin(VERSIONS))]
        for _, r in s5.iterrows():
            q = pw[(pw.cohort == c) & (pw.version == r.version) & (pw["mode"] == r["mode"]) & (pw.alpha == 0.5)]
            p50 = float(q.power_bonf.iloc[0]) if len(q) else np.nan
            p100 = pw[(pw.cohort == c) & (pw.version == r.version) & (pw["mode"] == r["mode"]) & (pw.alpha == 1.0)]
            if r.p_holm < 0.05:
                lab = "Difference detected (step 5)"
            elif np.isnan(p50):
                lab = "not yet computed"
            else:
                lab = "Likely no meaningful effect" if p50 >= 0.8 else "Too few samples to tell"
            cls.append(dict(cohort=c, version=r.version, mode=r["mode"], transition=r.transition, z=r.z,
                            p_holm=r.p_holm, power_at_50=p50,
                            power_at_100=float(p100.power_bonf.iloc[0]) if len(p100) else np.nan,
                            approx=("power from Normal->Adenoma group sizes" if r.transition != "Normal->Adenoma"
                                    else "power at actual group sizes"),
                            classification=lab))
    cls = pd.DataFrame(cls)
    cls.to_csv(os.path.join(OUT, "classification.csv"), index=False)

    L = ["DEVIATION D7 - Empirical power of the step 5 test (DESCRIPTIVE)", "=" * 64,
         "docs/analysis_deviations.md D7 (commit a446a34), specified before computation.",
         "Benchmark = each cohort's real raw Normal -> Adenoma difference; scaled by swapping samples between the",
         "groups so their composition differs by alpha of the full contrast. Step 5 test unchanged except 99",
         f"permutations; detection p <= 0.05/K (K=3 GSE41258, 2 others); {R} replicates per cell (MC SE <= 7%).",
         "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy")),
         f"Replicates finished: {len(reps)}/{len(COH) * 2 * 2 * len(ALPHAS) * R}", "",
         "POWER (Bonferroni-level detection rate; in brackets: p <= 0.05 rate)"]
    for (c, v, mo), g in pw.groupby(["cohort", "version", "mode"], sort=False):
        g = g.sort_values("alpha", ascending=False)
        L.append(f"  {c} {v:<10} {mo:<10} n={COH[c]['expected']['Normal']}/{COH[c]['expected']['Adenoma']}: " +
                 "  ".join(f"{int(a * 100)}%: {pb * 100:3.0f}% ({p5 * 100:3.0f}%)"
                           for a, pb, p5 in zip(g.alpha, g.power_bonf, g.power_05)))
    L += ["", "CLASSIFICATION OF STEP 5 PRIMARY RESULTS (rule: power at 50% >= 80% -> likely no meaningful effect)"]
    for _, r in cls.iterrows():
        L.append(f"  {r.cohort} {r.version:<10} {r['mode']:<10} {r.transition:<23} z={r.z:5.2f} Holm p {r.p_holm:.3f} "
                 f"| power@50% {r.power_at_50 * 100:3.0f}%  -> {r.classification}"
                 + ("   [approx.: N->A group sizes]" if r.transition != "Normal->Adenoma" and r.p_holm >= 0.05 else ""))
    L += ["", "Limitation: power for Adenoma->Carcinoma and Carcinoma->Metastasis is approximated by the Normal->Adenoma",
          "group sizes of the same cohort (same subsample size m); larger real groups would give somewhat higher power.",
          "Outputs: results/power_D7/power_table.csv, classification.csv, power_curves.png, checkpoints/"]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for ax, c in zip(axes, COH):
        for (v, mo), g in pw[pw.cohort == c].groupby(["version", "mode"]):
            g = g.sort_values("alpha")
            ax.plot(g.alpha * 100, g.power_bonf * 100, marker="o",
                    ls="-" if mo == "shape-only" else "--",
                    color="#2a78d6" if v == "unadjusted" else "#eb6834", label=f"{mo}, {v}")
        ax.axhline(80, color="0.6", ls=":"); ax.axvline(50, color="0.8", ls=":")
        ax.set_xlabel("effect size, % of the real raw Normal → Adenoma difference")
        ax.set_title(f"{c} (n = {COH[c]['expected']['Normal']}/{COH[c]['expected']['Adenoma']})")
        ax.set_xticks([25, 50, 75, 100])
    axes[0].set_ylabel("% of replicates detected (p ≤ 0.05/K)")
    axes[0].legend(fontsize=8)
    fig.suptitle("D7: empirical power of the step 5 test (descriptive)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "power_curves.png"), dpi=150)
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
        # slow units (GSE41258) first for load balance
        todo.sort(key=lambda u: u[0] != "GSE41258")
        print(f"[D7] {len(todo)} chunks to run", flush=True)
        t0 = time.time()
        with Pool(WORKERS) as pool:
            for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
                print(f"[D7] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
