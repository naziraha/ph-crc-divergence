"""
Step 4 - Persistent homology profiles per stage (GSE41258, primary set `main`).

Descriptive only (the inferential test is step 5). Pre-specification:
docs/analysis_prespecification.md (git commit 4f1f32f, 2026-10-04).

For each stage: B = 100 random subsamples of m = 80% of the smallest stage
(without replacement); Vietoris-Rips persistence (Euclidean, max dimension 2,
H0 and H1) on each subsample, raw and shape-only (distances / median pairwise
distance). Validated reference functions are used unchanged:
tda_crc_pilot.bootstrap_diagrams, persistence, summaries, betti_curve, tsi,
make_figure.

Per stage: total / maximum persistence, persistence entropy and number of
features for H0 and H1. Between stages: topological separation index
(TSI = mean between-stage W1 - mean within-stage W1, H0+H1) with bootstrap
95% CI, ratio between/within, centroid shift. No p-values here.

Inputs : results/step03/embeddings/GSE41258_main_{unadjusted,adjusted}_g2000_pc10.csv
Outputs: results/step04/*
"""
import os
import sys
import time
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

SEED = 42
B = 100                 # bootstrap subsamples per stage (pilot default; descriptive)
DIMS = (0, 1)
N_PAIRS_ADJ, N_PAIRS_ALL, N_CI = 400, 200, 500   # pilot defaults
EXPECTED = {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39}
OUT = os.path.join(ROOT, "results", "step04")
os.makedirs(OUT, exist_ok=True)
EMB = os.path.join(ROOT, "results", "step03", "embeddings")

t0 = time.time()
res = {}
for version in ("unadjusted", "adjusted"):
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(os.path.join(EMB, f"GSE41258_main_{version}_g2000_pc10.csv"), index_col=0)
    cnt = df["stage"].value_counts().to_dict()
    if cnt != EXPECTED:
        sys.exit(f"Stage counts changed ({version}): {cnt}")
    pcs = [c for c in df.columns if c.startswith("PC")]
    emb = df[pcs].values
    stages = df["stage"]
    pos = {s: np.where(stages.values == s)[0] for s in T.STAGES}
    m = max(10, int(0.8 * min(len(v) for v in pos.values())))   # as in the pilot

    diags = {"raw": {}, "norm": {}}
    rows = []
    for mode in ("raw", "norm"):
        for s in T.STAGES:
            diags[mode][s] = T.bootstrap_diagrams(emb, pos[s], m, B, mode == "norm", rng)
            for b, d in enumerate(diags[mode][s]):
                rows.append({"stage": s, "mode": mode, "draw": b, **T.summaries(d)})
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(OUT, f"stage_topology_bootstrap_{version}.csv"), index=False)

    adj = []
    for s, t in zip(T.STAGES[:-1], T.STAGES[1:]):
        rec = {"from": s, "to": t, "n_from": int(len(pos[s])), "n_to": int(len(pos[t])),
               "subsample_size": m}
        for mode in ("raw", "norm"):
            obs, btw, wth = T.tsi(diags[mode][s], diags[mode][t], DIMS, rng, N_PAIRS_ADJ)
            boots = [btw[rng.choice(len(btw), len(btw))].mean() - wth[rng.choice(len(wth), len(wth))].mean()
                     for _ in range(N_CI)]
            rec[f"tsi_{mode}"] = float(obs)
            rec[f"tsi_{mode}_ci"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
            rec[f"ratio_{mode}"] = float(btw.mean() / wth.mean())
        cs, ct = emb[pos[s]].mean(0), emb[pos[t]].mean(0)
        pooled = np.sqrt((emb[pos[s]].var(0).sum() + emb[pos[t]].var(0).sum()) / 2)
        rec["centroid_shift_sd"] = float(np.linalg.norm(cs - ct) / pooled)
        adj.append(rec)

    pairmat = pd.DataFrame(np.nan, index=T.STAGES, columns=T.STAGES)
    for i, s in enumerate(T.STAGES):
        for t in T.STAGES[i + 1:]:
            v = T.tsi(diags["norm"][s], diags["norm"][t], DIMS, rng, N_PAIRS_ALL)[0]
            pairmat.loc[s, t] = pairmat.loc[t, s] = v
    pairmat.to_csv(os.path.join(OUT, f"all_pairs_tsi_shape_only_{version}.csv"))

    # spread of each stage (median pairwise distance in PC space, full stage)
    from scipy.spatial.distance import pdist
    spread = {s: float(np.median(pdist(emb[pos[s]]))) for s in T.STAGES}

    # figure via the validated pilot function (columns renamed to 0..k-1 as it expects)
    allmax = max(d[1][:, 1].max() for s in T.STAGES for d in diags["norm"][s] if len(d[1]))
    grid = np.linspace(0, allmax * 1.05, 200)
    emb_df = pd.DataFrame(emb, index=df.index)
    T.make_figure(emb_df, stages, diags, adj, pairmat, summ, grid, OUT)
    for ext in ("png", "pdf"):
        dst = os.path.join(OUT, f"step04_topology_{version}.{ext}")
        if os.path.exists(dst):
            os.remove(dst)
        os.rename(os.path.join(OUT, f"tda_pilot_figure.{ext}"), dst)

    # Betti-1 curve summary (shape-only): peak height and location per stage
    betti = {}
    for s in T.STAGES:
        curves = np.array([T.betti_curve(d[1], grid) for d in diags["norm"][s]])
        mu = curves.mean(0)
        betti[s] = (float(mu.max()), float(grid[mu.argmax()]))
    res[version] = dict(summ=summ, adj=adj, pairmat=pairmat, m=m, spread=spread, betti=betti)
    pd.DataFrame(adj).to_json(os.path.join(OUT, f"adjacent_transitions_descriptive_{version}.json"),
                              orient="records", indent=2)
    print(f"[{version}] done, m = {m}, {time.time() - t0:.0f} s")

# shape-only TSI per adjacent transition, both versions (pilot panel D shares an axis with raw)
import matplotlib.pyplot as plt  # noqa: E402
fig, ax = plt.subplots(figsize=(7, 4))
labels = [f"{a['from'][:4]}→{a['to'][:4]}" for a in res["unadjusted"]["adj"]]
for k, (version, col) in enumerate((("unadjusted", "#2a78d6"), ("adjusted", "#eb6834"))):
    v = np.array([a["tsi_norm"] for a in res[version]["adj"]])
    ci = np.array([a["tsi_norm_ci"] for a in res[version]["adj"]])
    x = np.arange(3) + (k - 0.5) * 0.36
    ax.bar(x, v, width=0.34, color=col, label=version)
    ax.errorbar(x, v, yerr=[v - ci[:, 0], ci[:, 1] - v], fmt="none", ecolor="k", capsize=3, lw=1)
ax.set_xticks(range(3)); ax.set_xticklabels(labels)
ax.set_ylabel("TSI, shape-only (between − within W1)")
ax.set_title("Step 4: adjacent-stage topological separation (descriptive, 95% CI of pairs)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT, "step04_tsi_shape_only.png"), dpi=150)


# ---------------------------------------------------------------- summary
def q(v):
    return f"{np.median(v):.3f} [{np.percentile(v, 25):.3f}-{np.percentile(v, 75):.3f}]"


vers = {p: metadata.version(p) for p in ("gudhi", "numpy", "pandas", "scipy", "matplotlib")}
lines = ["STEP 4 - Persistent homology profiles per stage (GSE41258, set main)", "=" * 60,
         "Pre-specification: docs/analysis_prespecification.md, locked in git commit 4f1f32f "
         "(2026-10-04T11:35:19+08:00) before this step.",
         "Inputs: results/step03/embeddings/GSE41258_main_{unadjusted,adjusted}_g2000_pc10.csv",
         "  (adjusted embedding: liver-signature genes removed before gene selection)",
         f"Seed: {SEED} (numpy default_rng, same seed for both versions -> identical subsample draws)",
         "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()), "",
         "SETTINGS",
         "  Vietoris-Rips (GUDHI RipsComplex), Euclidean distance in 10-PC space, max dimension 2, H0 and H1",
         f"  {B} subsamples per stage of m = 80% of the smallest stage, without replacement",
         "  raw distances and shape-only (distances / subsample median pairwise distance)",
         f"  TSI = mean between-stage minus mean within-stage 1-Wasserstein distance, H0+H1;",
         f"  {N_PAIRS_ADJ} random pairs (adjacent stages; {N_PAIRS_ALL} for all pairs), 95% CI from "
         f"{N_CI} bootstrap resamples of the pairs",
         "  Descriptive only - no p-values. Inference is step 5 (pre-specified permutation test).", "",
         "SAMPLES PER STAGE: " + ", ".join(f"{k} {v}" for k, v in EXPECTED.items()) + " (verified)"]
for version in ("unadjusted", "adjusted"):
    r = res[version]
    s_ = r["summ"]
    lines += ["", f"=== {version.upper()} (subsample size m = {r['m']}) ===",
              "Stage spread (median pairwise distance in PC space, all samples of the stage):",
              "  " + ", ".join(f"{k} {v:.1f}" for k, v in r["spread"].items()),
              "Per-stage topology, median [IQR] over 100 subsamples:"]
    for mode, name in (("norm", "shape-only"), ("raw", "raw")):
        lines.append(f"  {name}:")
        lines.append(f"    {'stage':<11}{'H0 total':>24}{'H1 total':>24}{'H1 max':>24}"
                     f"{'H1 n loops':>16}{'H1 entropy':>24}")
        for st in T.STAGES:
            d = s_[(s_.stage == st) & (s_["mode"] == mode)]
            lines.append(f"    {st:<11}{q(d.H0_total_persistence):>24}{q(d.H1_total_persistence):>24}"
                         f"{q(d.H1_max_persistence):>24}{np.median(d.H1_n_features):>16.0f}"
                         f"{q(d.H1_entropy):>24}")
    lines.append("  Betti-1 curve (shape-only), peak mean number of loops @ scale:")
    lines.append("    " + ", ".join(f"{k} {v[0]:.1f} @ {v[1]:.2f}" for k, v in r["betti"].items()))
    lines.append("Adjacent transitions (descriptive):")
    for a in r["adj"]:
        lines.append(f"  {a['from']:>10} -> {a['to']:<10} TSI shape-only {a['tsi_norm']:.3f} "
                     f"[{a['tsi_norm_ci'][0]:.3f}, {a['tsi_norm_ci'][1]:.3f}] (between/within {a['ratio_norm']:.2f}) | "
                     f"TSI raw {a['tsi_raw']:.3f} [{a['tsi_raw_ci'][0]:.3f}, {a['tsi_raw_ci'][1]:.3f}] | "
                     f"centroid shift {a['centroid_shift_sd']:.2f} SD")
    lines.append("All stage pairs, TSI shape-only:")
    lines.append("  " + r["pairmat"].round(3).to_string().replace("\n", "\n  "))
lines += ["", "Outputs:",
          "  results/step04/stage_topology_bootstrap_<version>.csv   per-subsample topology summaries",
          "  results/step04/adjacent_transitions_descriptive_<version>.json",
          "  results/step04/all_pairs_tsi_shape_only_<version>.csv",
          "  results/step04/step04_topology_<version>.png/.pdf       6-panel figure (pilot layout;",
          "    panel D puts raw and shape-only TSI on one axis, so shape-only bars are not visible)",
          "  results/step04/step04_tsi_shape_only.png                 shape-only TSI, both versions",
          "Note: the TSI confidence intervals resample pairs of the same 100 subsamples; they show",
          "  the precision of the estimate, not sampling uncertainty, and are not a significance test.",
          f"Runtime: {time.time() - t0:.0f} s"]
with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
