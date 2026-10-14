"""
Step 5 (validation cohorts) - pre-specified Wasserstein permutation test, bulk sets.

Addendum 1 (6ea31ea) / addendum 2 (d73e5d1): identical test to GSE41258 step 5
(tda_crc_pilot.permutation_test unchanged; VR H0+H1; exact W1; k = 20; m = 80% of the
smallest stage; groups aligned; 999 permutations), Holm across the TWO transitions
(Normal->Adenoma, Adenoma->Carcinoma) per set x version x mode.
Versions: unadjusted, adjusted (primary, stage-blind), stageaware (secondary).
Tipping point = shape-only Holm p < 0.05 in BOTH unadjusted and adjusted.
Secondary replication = Normal->Adenoma raw Holm p < 0.05 in both.
Sensitivity: 5 / 20 PCs, 1,000 / 5,000 genes.
Resumable: one checkpoint file per unit; per-unit seeds from 42 and the unit name.

Usage: python scripts/step05_permutation_validation.py GSE20916 [primary|sensitivity|summary]
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
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
TRANS = [("Normal", "Adenoma"), ("Adenoma", "Carcinoma")]
VERSIONS = ("unadjusted", "adjusted", "stageaware")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
COHORT = {"GSE20916": {"set": "bulk", "expected": {"Normal": 24, "Adenoma": 45, "Carcinoma": 36}},
          "GSE37364": {"set": "bulk", "expected": {"Normal": 38, "Adenoma": 29, "Carcinoma": 27}},
          # addendum 6: this run = PRIMARY data versions only; stage-aware and sensitivities in a separate run
          "GSE117606": {"set": "main", "expected": {"Normal": 65, "Adenoma": 59, "Carcinoma": 71},
                        "primary_versions": ("unadjusted", "adjusted"), "prespec": "addendum 6 (93613bd)"}}
ANALYSES = {"primary": "g2000_pc10", "S5_pc5": "g2000_pc5", "S6_pc20": "g2000_pc20",
            "S7_genes1000": "g1000_pc10", "S8_genes5000": "g5000_pc10"}


def unit_name(gse, a, v, s, t, mode):
    return f"{gse}__{a}__{v[:5]}__{s[0]}-{t[0]}__{'shape' if mode == 'norm' else 'raw'}"


def units(gse, which):
    an = ["primary"] if which == "primary" else [a for a in ANALYSES if a != "primary"]
    vers = COHORT[gse].get("primary_versions", VERSIONS) if which == "primary" else VERSIONS
    return [(gse, a, v, s, t, mode) for a in an for v in vers for s, t in TRANS for mode in ("norm", "raw")]


def run_unit(u):
    gse, a, v, s, t, mode = u
    ck = os.path.join(ROOT, "results", "step05", gse, "checkpoints")
    name = unit_name(*u)
    path = os.path.join(ck, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    c = COHORT[gse]
    df = pd.read_csv(os.path.join(EMB, f"{gse}_{c['set']}_{v}_{ANALYSES[a]}.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != c["expected"]:
        raise SystemExit(f"Stage counts changed: {name}")
    emb = df[[x for x in df.columns if x.startswith("PC")]].values
    st = df["stage"].values
    m = max(10, int(0.8 * min((st == x).sum() for x in STAGES3)))
    idx_s, idx_t = np.where(st == s)[0], np.where(st == t)[0]
    Tobs, p, z, null = T.permutation_test(emb, idx_s, idx_t, m, K_DRAWS, DIMS, mode == "norm", N_PERM, rng,
                                          align=True)
    rec = dict(unit=name, cohort=gse, analysis=a, version=v, transition=f"{s}->{t}",
               mode="shape-only" if mode == "norm" else "raw", n_from=int(len(idx_s)), n_to=int(len(idx_t)),
               subsample_size=int(m), T_obs=float(Tobs), p=float(p), z=float(z), null=[float(x) for x in null],
               seconds=round(time.time() - t0, 1), finished=time.strftime("%Y-%m-%dT%H:%M:%S"))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"done in {rec['seconds']:.0f} s (p={p:.4f}, z={z:.1f})"


def summary(gse):
    out = os.path.join(ROOT, "results", "step05", gse)
    ck = os.path.join(out, "checkpoints")
    rows = []
    for f in sorted(os.listdir(ck)):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(ck, f)))
            r.pop("null")
            rows.append(r)
    if not rows:
        return
    df = pd.DataFrame(rows)
    df["p_holm"] = np.nan
    for _, g in df.groupby(["analysis", "version", "mode"]):
        if len(g) == 2:
            df.loc[g.index, "p_holm"] = T.holm(g["p"].values)
    df["_a"] = df.analysis.map({a: i for i, a in enumerate(ANALYSES)})
    df["_v"] = df.version.map({v: i for i, v in enumerate(VERSIONS)})
    df = df.sort_values(["_a", "_v", "mode", "transition"], ascending=[True, True, False, False]) \
        .drop(columns=["_a", "_v"])
    df.to_csv(os.path.join(out, "permutation_results_all.csv"), index=False)
    L = [f"STEP 5 - Pre-specified permutation test ({gse}, set {COHORT[gse]['set']})", "=" * 60,
         "Pre-specification: " + COHORT[gse].get("prespec", "addendum 1 (6ea31ea), addendum 2 (d73e5d1)"),
         f"Seed {SEED} (per-unit SeedSequence, resumable); "
         + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "scipy")),
         f"Settings: VR H0+H1, exact W1, k={K_DRAWS}, m = 80% of smallest stage, aligned, {N_PERM} permutations, "
         "Holm across 2 transitions per analysis x version x mode.",
         "Samples per stage: " + ", ".join(f"{k} {v}" for k, v in COHORT[gse]["expected"].items()),
         f"Units finished: {len(df)}/{len(units(gse, 'primary')) + len(units(gse, 'sensitivity'))}", ""]
    for a in ANALYSES:
        sub = df[df.analysis == a]
        if sub.empty:
            L.append(f"{a}: not run yet"); continue
        L.append(f"--- {a} ({ANALYSES[a]}) ---")
        for _, r in sub.iterrows():
            L.append(f"  {r.version:<11} {r['mode']:<10} {r.transition:<20} n={r.n_from}/{r.n_to} m={r.subsample_size}"
                     f"  T={r.T_obs:9.4f}  z={r.z:6.2f}  p={r.p:.4f}  p_Holm {r.p_holm:.4f}")
        if a == "primary" and len(sub) == len(units(gse, "primary")):
            L.append("  TIPPING-POINT DECISION (shape-only Holm p < 0.05 in BOTH unadjusted and adjusted):")
            for s, t in TRANS:
                tr = f"{s}->{t}"
                get = lambda v, mo: sub[(sub.version == v) & (sub["mode"] == mo) & (sub.transition == tr)].p_holm.iloc[0]
                pu, pa = get("unadjusted", "shape-only"), get("adjusted", "shape-only")
                ps = get("stageaware", "shape-only") if (sub.version == "stageaware").any() else float("nan")
                call = ("TIPPING POINT" if pu < .05 and pa < .05 else "significant in unadjusted only" if pu < .05
                        else "significant in adjusted only" if pa < .05 else "not significant")
                L.append(f"    {tr:<20} unadjusted {pu:.4f}, adjusted {pa:.4f} -> {call}   "
                         + (f"[secondary stage-aware: {ps:.4f}]" if ps == ps else "[stage-aware: separate run]"))
            ru, ra = (sub[(sub.version == v) & (sub["mode"] == "raw") & (sub.transition == "Normal->Adenoma")]
                      .p_holm.iloc[0] for v in ("unadjusted", "adjusted"))
            L.append(f"  SECONDARY REPLICATION (Normal->Adenoma raw Holm p < 0.05 in both): unadjusted {ru:.4f}, "
                     f"adjusted {ra:.4f} -> {'REPLICATED' if ru < .05 and ra < .05 else 'not replicated'}")
        L.append("")
    L += ["Outputs: results/step05/<GSE>/permutation_results_all.csv, checkpoints/<unit>.json, step05_primary.png",
          "Resume after interruption: re-run the same command."]
    with open(os.path.join(out, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    prim = df[df.analysis == "primary"]
    if len(prim) == len(units(gse, "primary")):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
        for ax, mode in zip(axes, ("shape-only", "raw")):
            for k, (v, col) in enumerate(zip([x for x in VERSIONS if (prim.version == x).any()], ("#2a78d6", "#eb6834", "#7b4fb3"))):
                g = prim[(prim.version == v) & (prim["mode"] == mode)].set_index("transition").loc[
                    [f"{s}->{t}" for s, t in TRANS]]
                x = np.arange(2) + (k - 1) * 0.27
                ax.bar(x, g.z.values, width=0.25, color=col, label=v)
                for xi, zz, ph in zip(x, g.z.values, g.p_holm.values):
                    ax.text(xi, max(zz, 0) + 0.3, f"p={ph:.3f}", ha="center", fontsize=7)
            ax.axhline(0, color="0.4", lw=1); ax.axhline(1.96, color="0.5", ls="--", lw=1)
            ax.set_xticks(range(2)); ax.set_xticklabels(["Norm→Aden", "Aden→Carc"])
            ax.set_title(f"{mode}{' (primary)' if mode == 'shape-only' else ' (secondary)'}"); ax.legend(fontsize=8)
        axes[0].set_ylabel("z vs permutation null (labels: Holm p)")
        fig.suptitle(f"{gse} {COHORT[gse]['set']}: permutation test (999 permutations)")
        fig.tight_layout(); fig.savefig(os.path.join(out, "step05_primary.png"), dpi=150); plt.close(fig)


def main():
    gse, which = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "primary")
    os.makedirs(os.path.join(ROOT, "results", "step05", gse, "checkpoints"), exist_ok=True)
    if which == "summary":
        summary(gse); return
    todo = [u for u in units(gse, which)
            if not os.path.exists(os.path.join(ROOT, "results", "step05", gse, "checkpoints", unit_name(*u) + ".json"))]
    print(f"[step5 {gse}] {which}: {len(todo)} units to run", flush=True)
    t0 = time.time()
    with Pool(WORKERS) as pool:
        for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
            print(f"[step5] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
            summary(gse)
    summary(gse)
    print(f"[step5 {gse}] {which} finished in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
