"""
Addendum 4 (commit 15cbf68) - standardised effect sizes with bootstrap 95% CIs (DESCRIPTIVE).

d = (mean between-stage W1 - mean within-stage W1) / SD(within-stage W1), from k = 20 subsample
diagrams per stage (reference bootstrap_diagrams / diag_dist, unchanged), m = 80% of the smallest
stage. Point estimate on the original samples; 1,000 bootstrap resamples of samples within each
stage (with replacement); percentile CI; bootstrap SE. Pooled across cohorts with DerSimonian-Laird
(variance = SE^2). Categories at d = 0.5 combined with the step 5 Holm p. Does not change the
pre-specified tipping-point decision. Resumable: one checkpoint per unit.

Usage: python scripts/effect_sizes_addendum4.py [run|summary]
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

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # scripts/archive/ -> repository root
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

SEED, K, N_BOOT, DIMS, WORKERS, THR = 42, 20, 1000, (0, 1), 6, 0.5
OUT = os.path.join(ROOT, "results", "effect_sizes_addendum4")
CKPT = os.path.join(OUT, "checkpoints")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
COHORTS = {
    "GSE41258": dict(prefix="GSE41258_main", stages=T.STAGES,
                     expected={"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39},
                     step5=os.path.join(ROOT, "results", "step05", "permutation_results_all.csv")),
    "GSE20916": dict(prefix="GSE20916_bulk", stages=["Normal", "Adenoma", "Carcinoma"],
                     expected={"Normal": 24, "Adenoma": 45, "Carcinoma": 36},
                     step5=os.path.join(ROOT, "results", "step05", "GSE20916", "permutation_results_all.csv")),
    "GSE37364": dict(prefix="GSE37364_bulk", stages=["Normal", "Adenoma", "Carcinoma"],
                     expected={"Normal": 38, "Adenoma": 29, "Carcinoma": 27},
                     step5=os.path.join(ROOT, "results", "step05", "GSE37364", "permutation_results_all.csv")),
}
VERSIONS, MODES = ("unadjusted", "adjusted"), ("norm", "raw")


def units():
    out = []
    for c, cfg in COHORTS.items():
        st = cfg["stages"]
        for v in VERSIONS:
            for s, t in zip(st[:-1], st[1:]):
                for mo in MODES:
                    out.append((c, v, s, t, mo))
    return out


def unit_name(c, v, s, t, mo):
    return f"{c}__{v[:5]}__{s[0]}-{t[0]}__{'shape' if mo == 'norm' else 'raw'}"


def effect(emb, idx_s, idx_t, m, normalize, rng):
    ds = T.bootstrap_diagrams(emb, idx_s, m, K, normalize, rng)
    dt = T.bootstrap_diagrams(emb, idx_t, m, K, normalize, rng)
    btw = np.array([T.diag_dist(a, b, DIMS) for a in ds for b in dt])
    wth = np.array([T.diag_dist(ds[a], ds[b], DIMS) for a in range(K) for b in range(a + 1, K)] +
                   [T.diag_dist(dt[a], dt[b], DIMS) for a in range(K) for b in range(a + 1, K)])
    return float((btw.mean() - wth.mean()) / wth.std(ddof=1)), float(btw.mean() - wth.mean())


def run_unit(u):
    c, v, s, t, mo = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    cfg = COHORTS[c]
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(os.path.join(EMB, f"{cfg['prefix']}_{v}_g2000_pc10.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != cfg["expected"]:
        raise SystemExit(f"Stage counts changed: {name}")
    emb = df[[x for x in df.columns if x.startswith("PC")]].values
    stv = df["stage"].values
    m = max(10, int(0.8 * min((stv == x).sum() for x in cfg["stages"])))
    idx_s, idx_t = np.where(stv == s)[0], np.where(stv == t)[0]
    d0, num0 = effect(emb, idx_s, idx_t, m, mo == "norm", rng)
    boots = []
    for _ in range(N_BOOT):
        bs = rng.choice(idx_s, size=len(idx_s), replace=True)
        bt = rng.choice(idx_t, size=len(idx_t), replace=True)
        boots.append(effect(emb, bs, bt, m, mo == "norm", rng)[0])
    boots = np.array(boots)
    rec = dict(unit=name, cohort=c, version=v, transition=f"{s}->{t}",
               mode="shape-only" if mo == "norm" else "raw", n_from=int(len(idx_s)), n_to=int(len(idx_t)),
               subsample_size=int(m), d=d0, numerator=num0, ci_low=float(np.percentile(boots, 2.5)),
               ci_upp=float(np.percentile(boots, 97.5)), se=float(boots.std(ddof=1)),
               boot_median=float(np.median(boots)), boots=[float(x) for x in boots],
               seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"d={d0:.2f} [{rec['ci_low']:.2f}, {rec['ci_upp']:.2f}] in {rec['seconds']:.0f} s"


def category(p_holm, lo, hi, pooled=False):
    if not pooled and p_holm < 0.05:
        return "Difference detected"
    if pooled and lo > THR:
        return "Consistent meaningful effect across cohorts"
    if hi < THR:
        return "Likely no meaningful effect"
    if lo > THR:
        return "Large effect but test not significant"
    return "Too few samples to tell"


def summary():
    rows = []
    for f in sorted(os.listdir(CKPT)):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(CKPT, f)))
            r.pop("boots")
            rows.append(r)
    df = pd.DataFrame(rows)
    ph = []
    for _, r in df.iterrows():
        s5 = pd.read_csv(COHORTS[r.cohort]["step5"])
        s5 = s5[(s5.analysis == "primary") & (s5.version == r.version) & (s5.transition == r.transition)
                & (s5["mode"] == r["mode"])]
        ph.append((float(s5.z.iloc[0]), float(s5.p_holm.iloc[0])))
    df["step5_z"], df["step5_p_holm"] = zip(*ph)
    df["category"] = [category(p, lo, hi) for p, lo, hi in zip(df.step5_p_holm, df.ci_low, df.ci_upp)]
    order = {u: i for i, u in enumerate(unit_name(*u) for u in units())}
    df = df.sort_values("unit", key=lambda x: x.map(order))
    df.to_csv(os.path.join(OUT, "effect_sizes_per_cohort.csv"), index=False)

    from statsmodels.stats.meta_analysis import combine_effects
    pooled = []
    for tr in ("Normal->Adenoma", "Adenoma->Carcinoma"):
        for mode in ("shape-only", "raw"):
            for v in VERSIONS:
                d = df[(df.transition == tr) & (df["mode"] == mode) & (df.version == v)]
                if len(d) != 3:
                    continue
                r = combine_effects(d.d.values, d.se.values ** 2, method_re="dl", row_names=list(d.cohort))
                re = r.summary_frame().loc["random effect"]
                pooled.append(dict(transition=tr, mode=mode, version=v, pooled_d=re["eff"], ci_low=re["ci_low"],
                                   ci_upp=re["ci_upp"], tau2=max(r.tau2, 0.0), I2=max(r.i2, 0.0),   # negative -> 0
                                   category=category(1.0, re["ci_low"], re["ci_upp"], pooled=True)))
    pooled = pd.DataFrame(pooled)
    pooled.to_csv(os.path.join(OUT, "effect_sizes_pooled.csv"), index=False)

    L = ["ADDENDUM 4 - Standardised effect sizes with bootstrap 95% CIs (DESCRIPTIVE)", "=" * 70,
         "*** NOTE (2026-10-06, PI review): THE BOOTSTRAP INTERVALS, THE CATEGORIES AND THE POOLED ESTIMATES ARE",
         "*** INVALID AND MUST NOT BE INTERPRETED. Reasons: (1) resampling with replacement puts duplicate samples",
         "*** into each stage; subsamples then contain only ~70-80% distinct samples and overlap heavily, so",
         "*** same-stage diagrams look artificially similar; (2) d divides by the SD of within-stage distances",
         "*** between overlapping subsamples, which this shrinks further. Net effect: bootstrap d is inflated,",
         "*** most for null transitions (point estimate below the whole interval in 11 of 13 non-significant",
         "*** results). The pooled estimates use bootstrap SEs and are invalid for the same reason.",
         "*** ONLY THE POINT ESTIMATES OF d (original samples) ARE REPORTED, AS DESCRIPTIVE VALUES. d is not on",
         "*** Cohen's scale (the within-stage SD comes from 80%-overlapping subsamples), so d = 0.5 is not a",
         "*** meaningful threshold here. Power and 'no effect vs too few samples' -> deviation D7.",
         "*** Negative I2 / tau2 are shown as 0.", "",
         "Pre-specification: docs/analysis_prespecification_addendum4_effectsizes.md (commit 15cbf68), written after",
         "all step 5 p-values were known. Does NOT change the pre-specified tipping-point decision.",
         f"d = (mean between W1 - mean within W1) / SD(within W1); k = {K} subsamples per stage; m = 80% of smallest stage;",
         f"CI = percentile of {N_BOOT} bootstrap resamples of samples within stage; seed {SEED} (per-unit SeedSequence).",
         f"Categories at d = {THR} (conventional 'medium'); combined with the step 5 Holm p.",
         "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "statsmodels")),
         f"Units finished: {len(df)}/{len(units())}", "",
         "PER COHORT"]
    for _, r in df.iterrows():
        L.append(f"  {r.cohort} {r.version:<10} {r['mode']:<10} {r.transition:<23} n={r.n_from}/{r.n_to} m={r.subsample_size} "
                 f"d = {r.d:5.2f} [{r.ci_low:5.2f}, {r.ci_upp:5.2f}] (boot median {r.boot_median:5.2f})  "
                 f"step5 Holm p {r.step5_p_holm:.3f}  -> {r.category}")
    L += ["", "POOLED (DerSimonian-Laird; variance = bootstrap SE^2)"]
    for _, r in pooled.iterrows():
        L.append(f"  {r.transition:<20} {r['mode']:<10} {r.version:<10} pooled d = {r.pooled_d:5.2f} "
                 f"[{r.ci_low:5.2f}, {r.ci_upp:5.2f}]  tau2 {r.tau2:.2f}  I2 {r.I2 * 100:3.0f}%  -> {r.category}")
    L += ["", "Limitation (stated in addendum 4): bootstrap resamples contain duplicate samples, which can bias d;",
          "the interval is a descriptive guide. Compare point estimate with bootstrap median as a bias indicator.",
          "Outputs: results/effect_sizes_addendum4/effect_sizes_per_cohort.csv, effect_sizes_pooled.csv,",
          "  checkpoints/<unit>.json (incl. all bootstrap values), effect_sizes_forest.png"]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    combos = [(m, v) for m in ("shape-only", "raw") for v in VERSIONS]
    fig, axes = plt.subplots(1, 4, figsize=(17, 6), sharey=True)
    labels = []
    for ax, (mode, v) in zip(axes, combos):
        sub = df[(df["mode"] == mode) & (df.version == v)]
        y, labels = 0, []
        for tr in ("Normal->Adenoma", "Adenoma->Carcinoma", "Carcinoma->Metastasis"):
            for _, r in sub[sub.transition == tr].iterrows():
                col = "#2a78d6" if r.step5_p_holm < 0.05 else "0.35"
                # CI drawn separately: the point estimate can lie outside the bootstrap interval (bias)
                ax.hlines(-y, r.ci_low, r.ci_upp, color=col, lw=1.5)
                ax.plot(r.d, -y, "o", color=col)
                ax.plot(r.boot_median, -y, "|", color=col, ms=10)
                labels.append(f"{tr.replace('->', '→')} {r.cohort}")
                y += 1
            p = pooled[(pooled.transition == tr) & (pooled["mode"] == mode) & (pooled.version == v)]
            if len(p):
                p = p.iloc[0]
                ax.hlines(-y, p.ci_low, p.ci_upp, color="#eb6834", lw=2)
                ax.plot(p.pooled_d, -y, "D", color="#eb6834", ms=7)
                labels.append(f"{tr.replace('->', '→')} pooled")
                y += 1
            y += 0.5
            labels.append("")
        ax.axvline(0, color="0.5", lw=1); ax.axvline(THR, color="0.6", ls=":", lw=1)
        ax.set_title(f"{mode}, {v}", fontsize=10); ax.set_xlabel("standardised effect size d (95% bootstrap CI)")
    ticks, lab, y = [], [], 0
    for l in labels:
        if l == "":
            y += 0.5
            continue
        ticks.append(-y); lab.append(l); y += 1
    axes[0].set_yticks(ticks); axes[0].set_yticklabels(lab, fontsize=8)
    fig.suptitle("Addendum 4 (descriptive): dot = point estimate, | = bootstrap median, line = 95% bootstrap CI; "
                 "blue = step 5 Holm p < 0.05, grey = n.s., orange = pooled; dotted d = 0.5", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "effect_sizes_forest.png"), dpi=150)
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
        print(f"[A4] {len(todo)} units to run", flush=True)
        t0 = time.time()
        with Pool(WORKERS) as pool:
            for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
                print(f"[A4] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
