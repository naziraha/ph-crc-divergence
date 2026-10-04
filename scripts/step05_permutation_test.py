"""
Step 5 - Pre-specified Wasserstein permutation test (GSE41258).

Pre-specification: docs/analysis_prespecification.md (git commit 4f1f32f, 2026-10-04).
Validated reference code used unchanged: tda_crc_pilot.permutation_test
(align_groups -> test_statistic -> label permutations), holm, drop_shared_patients.

Test: Vietoris-Rips H0+H1, exact 1-Wasserstein; statistic = mean between -
mean within over k = 20 subsamples per group of size m (80% of the smallest
stage); groups centred (and scaled, shape-only) before permuting; 999
permutations; Holm across the three adjacent transitions (per mode).
Primary = shape-only; raw reported as secondary.
Tipping point = shape-only Holm p < 0.05 in BOTH unadjusted and adjusted (primary set).

Resumable: every unit (analysis x version x transition x mode) is saved to
results/step05/checkpoints/<unit>.json as soon as it finishes; finished units are
skipped on re-run. Each unit has its own seed derived from 42 and the unit name,
so results do not depend on order, parallelism or interruptions.

Usage
  python scripts/step05_permutation_test.py primary        # primary analysis
  python scripts/step05_permutation_test.py sensitivity    # S1-S9
  python scripts/step05_permutation_test.py summary        # rebuild tables/summary only
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

SEED = 42
N_PERM = 999
K_DRAWS = 20
DIMS = (0, 1)
WORKERS = 6
OUT = os.path.join(ROOT, "results", "step05")
CKPT = os.path.join(OUT, "checkpoints")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
PROC = os.path.join(ROOT, "data", "processed")
TRANSITIONS = list(zip(T.STAGES[:-1], T.STAGES[1:]))
S9_DROP = ["GSM1012485", "GSM1012650", "GSM1012476"]   # metastases with liver fraction > 20%

# analysis -> (sample set, embedding tag, special handling)
ANALYSES = {
    "primary": ("main", "g2000_pc10", None),
    "S1_withlung": ("withlung", "g2000_pc10", None),
    "S2_main_allarrays": ("main_allarrays", "g2000_pc10", None),
    "S3_main_inclqc": ("main_inclqc", "g2000_pc10", None),
    "S4_one_per_patient": ("main", "g2000_pc10", "one_per_patient"),
    "S5_pc5": ("main", "g2000_pc5", None),
    "S6_pc20": ("main", "g2000_pc20", None),
    "S7_genes1000": ("main", "g1000_pc10", None),
    "S8_genes5000": ("main", "g5000_pc10", None),
    "S9_no_high_liver_mets": ("main", "g2000_pc10", "drop_high_liver"),
}
EXPECTED = {"main": {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39},
            "withlung": {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 58},
            "main_allarrays": {"Normal": 54, "Adenoma": 51, "Carcinoma": 184, "Metastasis": 47},
            "main_inclqc": {"Normal": 53, "Adenoma": 48, "Carcinoma": 183, "Metastasis": 39}}


def unit_name(analysis, version, s, t, mode):
    # short names keep Windows paths < 260 characters
    v = "unadj" if version == "unadjusted" else "adj"
    return f"{analysis}__{v}__{s[0]}-{t[0]}__{'shape' if mode == 'norm' else 'raw'}"


def units_for(analyses):
    out = []
    for a in analyses:
        trans = [("Carcinoma", "Metastasis")] if ANALYSES[a][2] == "drop_high_liver" else TRANSITIONS
        for v in ("unadjusted", "adjusted"):
            for s, t in trans:
                for mode in ("norm", "raw"):
                    out.append((a, v, s, t, mode))
    return out


def load(analysis, version):
    set_name, tag, special = ANALYSES[analysis]
    df = pd.read_csv(os.path.join(EMB, f"GSE41258_{set_name}_{version}_{tag}.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != EXPECTED[set_name]:
        raise SystemExit(f"Stage counts changed: {analysis} {version}")
    lab = pd.read_csv(os.path.join(PROC, f"GSE41258_{set_name}_samples_qc.csv"), index_col=0,
                      dtype={"patient_id": str})
    if special == "drop_high_liver":
        assert all(s in df.index for s in S9_DROP)
        df = df.drop(index=S9_DROP)
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    stages = df["stage"].values
    patient = lab.loc[df.index, "patient_id"].astype(object).where(
        lab.loc[df.index, "patient_id"].notna(), None).values
    return df, emb, stages, patient


def run_unit(u):
    analysis, version, s, t, mode = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df, emb, stages, patient = load(analysis, version)
    special = ANALYSES[analysis][2]
    counts = {st: int((stages == st).sum()) for st in T.STAGES}
    m = max(10, int(0.8 * min(counts.values())))          # 80% of smallest stage (pilot rule)
    idx_s, idx_t = np.where(stages == s)[0], np.where(stages == t)[0]
    shared_p = len(set(patient[idx_s]) & set(patient[idx_t]) - {None})
    n_dropped = 0
    if special == "one_per_patient":
        idx_s, idx_t, n_dropped = T.drop_shared_patients(idx_s, idx_t, patient, rng)
    m_pair = min(m, int(0.8 * min(len(idx_s), len(idx_t))))   # pilot rule
    Tobs, p, z, null = T.permutation_test(emb, idx_s, idx_t, m_pair, K_DRAWS, DIMS,
                                          mode == "norm", N_PERM, rng, align=True)
    rec = dict(unit=name, analysis=analysis, version=version, transition=f"{s}->{t}",
               mode="shape-only" if mode == "norm" else "raw", n_from=int(len(idx_s)),
               n_to=int(len(idx_t)), subsample_size=int(m_pair), k_draws=K_DRAWS, n_perm=N_PERM,
               shared_patients=int(shared_p), shared_patients_resolved=int(n_dropped),
               T_obs=float(Tobs), p=float(p), z=float(z), null_mean=float(null.mean()),
               null_sd=float(null.std(ddof=1)), null=[float(x) for x in null],
               seconds=round(time.time() - t0, 1),
               finished=time.strftime("%Y-%m-%dT%H:%M:%S"),
               versions={k: metadata.version(k) for k in ("gudhi", "numpy", "scipy")})
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(rec, fh)
    os.replace(tmp, path)          # atomic: a half-written file never counts as finished
    return name, f"done in {rec['seconds']:.0f} s (p={p:.4f}, z={z:.1f})"


# ---------------------------------------------------------------- results tables
def collect():
    rows = []
    for f in sorted(os.listdir(CKPT)):
        if f.endswith(".json"):
            with open(os.path.join(CKPT, f)) as fh:
                r = json.load(fh)
            r.pop("null")
            rows.append(r)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["p_holm"] = np.nan
    for (a, v, mode), g in df.groupby(["analysis", "version", "mode"]):
        if a.startswith("S9") or len(g) != 3:
            continue
        df.loc[g.index, "p_holm"] = T.holm(g["p"].values)
    order = {f"{s}->{t}": i for i, (s, t) in enumerate(TRANSITIONS)}
    df["_o"] = df["transition"].map(order)
    an_order = {a: i for i, a in enumerate(ANALYSES)}
    df["_a"] = df["analysis"].map(an_order)
    return df.sort_values(["_a", "version", "mode", "_o"]).drop(columns=["_o", "_a"])


def write_summary():
    df = collect()
    if df.empty:
        return
    df.drop(columns=["versions"]).to_csv(os.path.join(OUT, "permutation_results_all.csv"), index=False)
    total = len(units_for(list(ANALYSES)))
    vers = {p: metadata.version(p) for p in ("gudhi", "numpy", "scipy", "pandas")}
    lines = ["STEP 5 - Pre-specified Wasserstein permutation test (GSE41258)", "=" * 60,
             "Pre-specification: docs/analysis_prespecification.md, git commit 4f1f32f "
             "(2026-10-04T11:35:19+08:00)",
             "Inputs: results/step03/embeddings/*.csv; data/processed/GSE41258_<set>_samples_qc.csv",
             f"Seed: {SEED}; each unit seeded by SeedSequence(42, crc32(unit name)) -> resumable, order-free",
             "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()),
             f"Settings: VR H0+H1, exact 1-Wasserstein, k = {K_DRAWS} subsamples per group, m = 80% of "
             f"smallest stage, groups aligned before permuting, {N_PERM} permutations, Holm across the 3 "
             "transitions per analysis x version x mode (S9: single comparison, no Holm).",
             "Reference functions used unchanged: permutation_test, holm, drop_shared_patients.",
             f"Units finished: {len(df)}/{total}", ""]

    def block(a):
        sub = df[df.analysis == a]
        if sub.empty:
            return [f"{a}: not run yet", ""]
        set_name = ANALYSES[a][0]
        out = [f"--- {a} (set {set_name}, embedding {ANALYSES[a][1]}) ---"]
        for v in ("unadjusted", "adjusted"):
            for mode in ("shape-only", "raw"):
                g = sub[(sub.version == v) & (sub["mode"] == mode)]
                for _, r in g.iterrows():
                    holm_txt = f"p_Holm {r.p_holm:.4f}" if pd.notna(r.p_holm) else "p_Holm  n/a "
                    out.append(f"  {v:<10} {mode:<10} {r.transition:<22} n={r.n_from}/{r.n_to} m={r.subsample_size}"
                               f"  T={r.T_obs:8.4f}  z={r.z:6.2f}  p={r.p:.4f}  {holm_txt}"
                               + (f"  shared patients {r.shared_patients}" if a == "primary" else "")
                               + (f"  (patients resolved {r.shared_patients_resolved})"
                                  if a.startswith("S4") else ""))
        return out + [""]

    lines += ["SAMPLES PER STAGE (primary set main): " +
              ", ".join(f"{k} {v}" for k, v in EXPECTED["main"].items())]
    lines += ["", "PRIMARY ANALYSIS"] + block("primary")
    prim = df[(df.analysis == "primary") & (df["mode"] == "shape-only")]
    if len(prim) == 6:
        lines.append("TIPPING-POINT DECISION (pre-specified: shape-only Holm p < 0.05 in BOTH versions):")
        for s, t in TRANSITIONS:
            tr = f"{s}->{t}"
            pu = prim[(prim.version == "unadjusted") & (prim.transition == tr)].p_holm.iloc[0]
            pa = prim[(prim.version == "adjusted") & (prim.transition == tr)].p_holm.iloc[0]
            call = ("TIPPING POINT" if (pu < 0.05 and pa < 0.05) else
                    "significant in unadjusted only" if pu < 0.05 else
                    "significant in adjusted only" if pa < 0.05 else "not significant")
            lines.append(f"  {tr:<24} Holm p unadjusted {pu:.4f}, adjusted {pa:.4f} -> {call}")
        lines.append("")
    lines += ["SENSITIVITY ANALYSES (descriptive support; do not replace the primary result)"]
    for a in ANALYSES:
        if a != "primary":
            lines += block(a)
    lines += ["Outputs:",
              "  results/step05/permutation_results_all.csv   one row per unit (T, z, p, Holm p, n, m)",
              "  results/step05/checkpoints/<unit>.json        per-unit result incl. full null distribution",
              "  results/step05/step05_primary.png             primary result figure",
              "Resume after interruption: re-run the same command; finished units are skipped."]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    # primary figure: z-scores with Holm p, both versions
    if len(df[df.analysis == "primary"]) == 12:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
        for ax, mode in zip(axes, ("shape-only", "raw")):
            for k, (v, col) in enumerate((("unadjusted", "#2a78d6"), ("adjusted", "#eb6834"))):
                g = df[(df.analysis == "primary") & (df.version == v) & (df["mode"] == mode)]
                x = np.arange(3) + (k - 0.5) * 0.36
                ax.bar(x, g.z.values, width=0.34, color=col, label=v)
                for xi, zz, ph in zip(x, g.z.values, g.p_holm.values):
                    ax.text(xi, max(zz, 0) + 0.3, f"p={ph:.3f}", ha="center", fontsize=7)
            ax.axhline(0, color="0.4", lw=1)
            ax.axhline(1.96, color="0.5", ls="--", lw=1)
            ax.set_xticks(range(3)); ax.set_xticklabels([f"{s[:4]}→{t[:4]}" for s, t in TRANSITIONS])
            ax.set_title(f"{mode}{' (primary)' if mode == 'shape-only' else ' (secondary)'}")
            ax.legend(fontsize=8)
        axes[0].set_ylabel("z vs permutation null (labels: Holm p)")
        fig.suptitle("Step 5: primary permutation test, set main (999 permutations)")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "step05_primary.png"), dpi=150)
        plt.close(fig)


def main():
    os.makedirs(CKPT, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "primary"
    if which == "summary":
        write_summary(); print(open(os.path.join(OUT, "summary.txt"), encoding="utf-8").read()); return
    analyses = ["primary"] if which == "primary" else [a for a in ANALYSES if a != "primary"]
    todo = units_for(analyses)
    left = [u for u in todo if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
    print(f"[step5] {which}: {len(todo)} units, {len(todo) - len(left)} already finished, "
          f"{len(left)} to run on {WORKERS} workers", flush=True)
    t0 = time.time()
    with Pool(WORKERS) as pool:
        for i, (name, msg) in enumerate(pool.imap_unordered(run_unit, left), 1):
            print(f"[step5] {i}/{len(left)} {name}: {msg}  [{(time.time() - t0) / 60:.1f} min]", flush=True)
            write_summary()
    write_summary()
    print(f"[step5] {which} finished in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
