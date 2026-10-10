"""
Step 4 (validation cohorts) - persistent homology profiles per stage (descriptive).

Bulk sets: frozen rule, identical to GSE41258 step 4 (B = 100 subsamples of m = 80% of the
smallest stage; raw and shape-only; H0+H1; TSI between adjacent stages with 95% CI of pairs).
Microdissected sets (deviation D4): one diagram per stage on all samples; W1 (H0+H1) between
stage diagrams; within-stage spread for micro_crypt unadjusted.
Validated reference functions used unchanged: bootstrap_diagrams, persistence, summaries,
betti_curve, tsi, diag_dist, make_figure.

Usage: python scripts/step04_persistence_validation.py GSE20916
"""
import json
import os
import sys
import time
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402

SEED = 42
B, DIMS = 100, (0, 1)
N_PAIRS_ADJ, N_PAIRS_ALL, N_CI = 400, 200, 500
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
VERSIONS = ("unadjusted", "adjusted", "stageaware")
COHORT = {"GSE20916": {"bulk": ["bulk"], "micro": ["micro_pooled", "micro_crypt", "micro_mucosa"],
                       "expected": {"bulk": {"Normal": 24, "Adenoma": 45, "Carcinoma": 36},
                                    "micro_pooled": {"Normal": 10, "Adenoma": 10, "Carcinoma": 10},
                                    "micro_crypt": {"Normal": 5, "Adenoma": 5, "Carcinoma": 5},
                                    "micro_mucosa": {"Normal": 5, "Adenoma": 5, "Carcinoma": 5}}},
          "GSE37364": {"bulk": ["bulk"], "micro": [],
                       "expected": {"bulk": {"Normal": 38, "Adenoma": 29, "Carcinoma": 27}}},
          "GSE117606": {"bulk": ["main"], "micro": [], "versions": ("unadjusted", "adjusted"),
                        "expected": {"main": {"Normal": 65, "Adenoma": 59, "Carcinoma": 71}}},
          # addendum 10: carcinoma stages; same frozen rule (m = 80% of smallest stage); 3 adjacent transitions
          "TCGA_COAD": {"bulk": ["main"], "micro": [], "versions": ("unadjusted", "adjusted"),
                        "stages": ["Stage I", "Stage II", "Stage III", "Stage IV"],
                        "colors": {"Stage I": "#7fb800", "Stage II": "#e0a100", "Stage III": "#e05a00",
                                   "Stage IV": "#c0392b"},
                        "expected": {"main": {"Stage I": 45, "Stage II": 110, "Stage III": 80, "Stage IV": 39}}}}


def load(gse, aset, version, expected):
    df = pd.read_csv(os.path.join(EMB, f"{gse}_{aset}_{version}_g2000_pc10.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != expected:
        raise SystemExit(f"Stage counts changed: {aset} {version}")
    pcs = [c for c in df.columns if c.startswith("PC")]
    return df, df[pcs].values, df["stage"]


def q(v):
    return f"{np.median(v):.3f} [{np.percentile(v, 25):.3f}-{np.percentile(v, 75):.3f}]"


def main():
    global STAGES3
    gse = sys.argv[1]
    cfg = COHORT[gse]
    versions = cfg.get("versions", VERSIONS)
    if "stages" in cfg:            # cohort-specific stage order; figure labels/colours set on the reference module
        STAGES3 = cfg["stages"]
        T.STAGES = list(cfg["stages"])
        T.COLORS.update(cfg["colors"])
    out = os.path.join(ROOT, "results", "step04", gse)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    L = [f"STEP 4 - Persistent homology profiles per stage ({gse}; descriptive)", "=" * 60,
         "Pre-specification: addendum 1 (6ea31ea), addendum 2 (d73e5d1); deviation D4 (ac2a591)",
         "Inputs: results/step03/embeddings/<GSE>_<set>_<version>_g2000_pc10.csv",
         f"Seed: {SEED}",
         "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "pandas", "scipy")),
         "Settings: Vietoris-Rips (GUDHI), Euclidean in 10-PC space, max dimension 2, H0 and H1; raw and shape-only.",
         "  Bulk: 100 subsamples per stage of m = 80% of the smallest stage (frozen rule). TSI = mean between - mean",
         f"  within W1 (H0+H1), {N_PAIRS_ADJ} pairs, 95% CI from {N_CI} pair resamples (precision, not a test).",
         "  Micro (D4): one diagram per stage on all samples; W1 between stage diagrams. No test.", ""]

    # ---------------------------------------------------------------- bulk (frozen rule)
    for aset in cfg["bulk"]:
        res = {}
        for version in versions:
            rng = np.random.default_rng(SEED)
            df, emb, stages = load(gse, aset, version, cfg["expected"][aset])
            pos = {s: np.where(stages.values == s)[0] for s in STAGES3}
            m = max(10, int(0.8 * min(len(v) for v in pos.values())))
            diags = {"raw": {}, "norm": {}}
            rows = []
            for mode in ("raw", "norm"):
                for s in STAGES3:
                    diags[mode][s] = T.bootstrap_diagrams(emb, pos[s], m, B, mode == "norm", rng)
                    rows += [{"stage": s, "mode": mode, "draw": b, **T.summaries(d)}
                             for b, d in enumerate(diags[mode][s])]
            summ = pd.DataFrame(rows)
            summ.to_csv(os.path.join(out, f"stage_topology_bootstrap_{aset}_{version}.csv"), index=False)
            adj = []
            for s, t in zip(STAGES3[:-1], STAGES3[1:]):
                rec = {"from": s, "to": t, "n_from": int(len(pos[s])), "n_to": int(len(pos[t])), "subsample_size": m}
                for mode in ("raw", "norm"):
                    obs, btw, wth = T.tsi(diags[mode][s], diags[mode][t], DIMS, rng, N_PAIRS_ADJ)
                    boots = [btw[rng.choice(len(btw), len(btw))].mean() - wth[rng.choice(len(wth), len(wth))].mean()
                             for _ in range(N_CI)]
                    rec[f"tsi_{mode}"] = float(obs)
                    rec[f"tsi_{mode}_ci"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
                cs, ct = emb[pos[s]].mean(0), emb[pos[t]].mean(0)
                pooled = np.sqrt((emb[pos[s]].var(0).sum() + emb[pos[t]].var(0).sum()) / 2)
                rec["centroid_shift_sd"] = float(np.linalg.norm(cs - ct) / pooled)
                adj.append(rec)
            pairmat = pd.DataFrame(np.nan, index=STAGES3, columns=STAGES3)
            for i, s in enumerate(STAGES3):
                for t in STAGES3[i + 1:]:
                    pairmat.loc[s, t] = pairmat.loc[t, s] = T.tsi(diags["norm"][s], diags["norm"][t], DIMS, rng,
                                                                  N_PAIRS_ALL)[0]
            allmax = max(d[1][:, 1].max() for s in STAGES3 for d in diags["norm"][s] if len(d[1]))
            grid = np.linspace(0, allmax * 1.05, 200)
            T.make_figure(pd.DataFrame(emb, index=df.index), stages, diags, adj, pairmat, summ, grid, out)
            for ext in ("png", "pdf"):
                dst = os.path.join(out, f"step04_topology_{aset}_{version}.{ext}")
                if os.path.exists(dst):
                    os.remove(dst)
                os.rename(os.path.join(out, f"tda_pilot_figure.{ext}"), dst)
            with open(os.path.join(out, f"adjacent_transitions_descriptive_{aset}_{version}.json"), "w") as fh:
                json.dump(adj, fh, indent=2)
            spread = {s: float(np.median(pdist(emb[pos[s]]))) for s in STAGES3}
            res[version] = (summ, adj, pairmat, m, spread)
        L += [f"=== {aset} (bulk; " + ", ".join(f"{k} {v}" for k, v in cfg["expected"][aset].items()) + ") ==="]
        for version in versions:
            summ, adj, pairmat, m, spread = res[version]
            L += [f"--- {version} (m = {m}) ---",
                  "  stage spread (median pairwise distance): " + ", ".join(f"{k} {v:.1f}" for k, v in spread.items()),
                  "  shape-only per stage, median [IQR]:  H1 total | H1 n loops"]
            for st in STAGES3:
                d = summ[(summ.stage == st) & (summ["mode"] == "norm")]
                L.append(f"    {st:<10} {q(d.H1_total_persistence)} | {np.median(d.H1_n_features):.0f}")
            for a in adj:
                L.append(f"  {a['from']:>9} -> {a['to']:<9} TSI shape-only {a['tsi_norm']:.3f} "
                         f"[{a['tsi_norm_ci'][0]:.3f}, {a['tsi_norm_ci'][1]:.3f}] | TSI raw {a['tsi_raw']:.2f} "
                         f"[{a['tsi_raw_ci'][0]:.2f}, {a['tsi_raw_ci'][1]:.2f}] | centroid shift {a['centroid_shift_sd']:.2f} SD")
        # shape-only TSI figure across versions
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        for k, (version, col) in enumerate(zip(versions, ("#2a78d6", "#eb6834", "#7b4fb3"))):
            adj = res[version][1]
            v = np.array([a["tsi_norm"] for a in adj]); ci = np.array([a["tsi_norm_ci"] for a in adj])
            x = np.arange(len(adj)) + (k - 1) * 0.27
            ax.bar(x, v, width=0.25, color=col, label=version)
            ax.errorbar(x, v, yerr=[v - ci[:, 0], ci[:, 1] - v], fmt="none", ecolor="k", capsize=3, lw=1)
        ax.set_xticks(range(len(adj)))
        ax.set_xticklabels([f"{a['from'][:4]}→{a['to'][:4]}" if "stages" not in cfg
                            else f"{a['from'].replace('Stage ', '')}→{a['to'].replace('Stage ', '')}" for a in adj])
        ax.set_ylabel("TSI, shape-only"); ax.legend(fontsize=8)
        ax.set_title(f"{gse} {aset}: adjacent-stage separation (descriptive)")
        fig.tight_layout(); fig.savefig(os.path.join(out, f"step04_tsi_shape_only_{aset}.png"), dpi=150)
        plt.close(fig)
        L.append("")

    if not cfg["micro"]:      # bulk-only cohort
        L += ["Outputs:",
              f"  results/step04/{gse}/step04_topology_{cfg['bulk'][0]}_<version>.png, step04_tsi_shape_only_{cfg['bulk'][0]}.png",
              f"  results/step04/{gse}/stage_topology_bootstrap_{cfg['bulk'][0]}_<version>.csv, adjacent_transitions_descriptive_*.json",
              f"Runtime: {time.time() - t0:.0f} s"]
        with open(os.path.join(out, "summary.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(L) + "\n")
        print("\n".join(L))
        return

    # ---------------------------------------------------------------- micro (D4)
    L += ["=== Microdissected sets (D4: one diagram per stage, all samples; descriptive) ==="]
    rows = []
    for aset in cfg["micro"]:
        for version in versions:
            df, emb, stages = load(gse, aset, version, cfg["expected"][aset])
            dg = {mode: {s: T.persistence(emb[stages.values == s], mode == "norm") for s in STAGES3}
                  for mode in ("raw", "norm")}
            rec = {"set": aset, "version": version}
            for mode in ("raw", "norm"):
                for s in STAGES3:
                    sm = T.summaries(dg[mode][s])
                    rec[f"{mode}_{s}_H0_total"] = sm["H0_total_persistence"]
                    rec[f"{mode}_{s}_H1_total"] = sm["H1_total_persistence"]
                    rec[f"{mode}_{s}_H1_n"] = sm["H1_n_features"]
                for s, t in [("Normal", "Adenoma"), ("Adenoma", "Carcinoma"), ("Normal", "Carcinoma")]:
                    rec[f"{mode}_W1_{s[0]}{t[0]}"] = T.diag_dist(dg[mode][s], dg[mode][t], DIMS)
            for s in STAGES3:
                rec[f"spread_{s}"] = float(np.median(pdist(emb[stages.values == s])))
            rows.append(rec)
    micro = pd.DataFrame(rows)
    micro.to_csv(os.path.join(out, "micro_stage_diagrams_summary.csv"), index=False)
    for aset in cfg["micro"]:
        L.append(f"--- {aset} ---")
        for _, r in micro[micro.set == aset].iterrows():
            L.append(f"  {r.version:<11} shape-only W1 between stage diagrams: N-A {r.norm_W1_NA:.3f}, "
                     f"A-C {r.norm_W1_AC:.3f}, N-C {r.norm_W1_NC:.3f} | H1 loops N/A/C "
                     f"{r.norm_Normal_H1_n:.0f}/{r.norm_Adenoma_H1_n:.0f}/{r.norm_Carcinoma_H1_n:.0f} | "
                     f"spread N/A/C {r.spread_Normal:.1f}/{r.spread_Adenoma:.1f}/{r.spread_Carcinoma:.1f}")
    L += ["  Caution: 5-10 points per diagram; H1 loops need >= 4 points, so loop counts are near-uninformative.",
          "  Adjusted micro versions overfit (step 2); unadjusted crypt epithelium is the emphasised view.", ""]

    # ---------------------------------------------------------------- crypt spread (D4b)
    df, emb, stages = load(gse, "micro_crypt", "unadjusted", cfg["expected"]["micro_crypt"])
    sp = []
    for s in STAGES3:
        E = emb[stages.values == s]
        cen = np.linalg.norm(E - E.mean(0), axis=1)
        for sid, dcen in zip(df.index[stages.values == s], cen):
            sp.append({"sample": sid, "stage": s, "distance_to_centroid": float(dcen)})
    sp = pd.DataFrame(sp)
    sp.to_csv(os.path.join(out, "crypt_unadjusted_within_stage_spread.csv"), index=False)
    pw = {s: pdist(emb[stages.values == s]) for s in STAGES3}
    L += ["=== Within-stage spread, crypt epithelium, unadjusted (D4b; 2,000 genes, 10 PCs; descriptive) ==="]
    for s in STAGES3:
        c = sp[sp.stage == s].distance_to_centroid
        L.append(f"  {s:<10} median pairwise distance {np.median(pw[s]):6.2f} (range {pw[s].min():.1f}-{pw[s].max():.1f}); "
                 f"distance to centroid median {c.median():.2f} (range {c.min():.1f}-{c.max():.1f})")
    L.append(f"  Ratio adenoma/normal (median pairwise): {np.median(pw['Adenoma']) / np.median(pw['Normal']):.2f}; "
             f"carcinoma/normal: {np.median(pw['Carcinoma']) / np.median(pw['Normal']):.2f}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for i, s in enumerate(STAGES3):
        v = sp[sp.stage == s].distance_to_centroid
        ax.scatter(np.full(len(v), i) + np.linspace(-0.1, 0.1, len(v)), v, c=T.COLORS[s], s=40)
        ax.hlines(v.median(), i - 0.25, i + 0.25, color="k")
    ax.set_xticks(range(3)); ax.set_xticklabels(STAGES3)
    ax.set_ylabel("distance to stage centroid (10-PC space)")
    ax.set_title(f"{gse} crypt epithelium, unadjusted: within-stage spread")
    fig.tight_layout(); fig.savefig(os.path.join(out, "crypt_within_stage_spread.png"), dpi=150)
    L += ["", "Outputs:",
          f"  results/step04/{gse}/step04_topology_{cfg['bulk'][0]}_<version>.png, step04_tsi_shape_only_{cfg['bulk'][0]}.png",
          f"  results/step04/{gse}/stage_topology_bootstrap_{cfg['bulk'][0]}_<version>.csv, adjacent_transitions_descriptive_*.json",
          f"  results/step04/{gse}/micro_stage_diagrams_summary.csv",
          f"  results/step04/{gse}/crypt_unadjusted_within_stage_spread.csv, crypt_within_stage_spread.png",
          f"Runtime: {time.time() - t0:.0f} s"]
    with open(os.path.join(out, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
