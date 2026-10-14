"""
Deviation D6 - software cross-check of the GSE41258 step 4 persistence diagrams: GUDHI vs Ripser.

Cross-check ONLY (no result, test or decision uses Ripser). Set main, unadjusted and adjusted
10-PC embeddings. Replays the step 4 random draws exactly (default_rng(42) per version; order:
mode raw then shape-only, stages Normal..Metastasis, 100 draws of m = 31 each, as in
tda_crc_pilot.bootstrap_diagrams), and checks that the GUDHI summaries reproduce the saved step 4
file before comparing.
Per subsample: GUDHI diagram = tda_crc_pilot.persistence (unchanged); Ripser diagram =
ripser(distance matrix, maxdim=1), same distance matrix (shape-only: divided by the median
pairwise distance), infinite H0 bar removed. Compared with the exact 1-Wasserstein distance
(tda_crc_pilot.wasserstein), feature counts (bars with persistence > 1e-12), total persistence.
"""
import os
import sys
import time
from importlib import metadata

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from ripser import ripser

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # scripts/archive/ -> repository root
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

SEED, B = 42, 100
OUT = os.path.join(ROOT, "results", "step04", "ripser_crosscheck")
os.makedirs(OUT, exist_ok=True)
EXPECTED = {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39}
EPS = 1e-12


def ripser_diagram(points, normalize):
    D = squareform(pdist(points))
    if normalize:
        D = D / np.median(D[np.triu_indices_from(D, 1)])
    dg = ripser(D, distance_matrix=True, maxdim=1)["dgms"]
    out = {}
    for dim in (0, 1):
        d = np.asarray(dg[dim], dtype=float)
        out[dim] = d[np.isfinite(d[:, 1])] if len(d) else np.zeros((0, 2))
    return out


def nfeat(d):
    return int(((d[:, 1] - d[:, 0]) > EPS).sum()) if len(d) else 0


def main():
    t0 = time.time()
    rows = []
    for version in ("unadjusted", "adjusted"):
        df = pd.read_csv(os.path.join(ROOT, "results", "step03", "embeddings",
                                      f"GSE41258_main_{version}_g2000_pc10.csv"), index_col=0)
        if df["stage"].value_counts().to_dict() != EXPECTED:
            raise SystemExit("Stage counts changed")
        emb = df[[c for c in df.columns if c.startswith("PC")]].values
        st = df["stage"].values
        pos = {s: np.where(st == s)[0] for s in T.STAGES}
        m = max(10, int(0.8 * min(len(v) for v in pos.values())))
        rng = np.random.default_rng(SEED)
        for mode in ("raw", "norm"):
            for s in T.STAGES:
                for b in range(B):
                    take = rng.choice(pos[s], size=m, replace=False)      # identical to bootstrap_diagrams
                    g = T.persistence(emb[take], mode == "norm")
                    r = ripser_diagram(emb[take], mode == "norm")
                    rec = {"version": version, "mode": "shape-only" if mode == "norm" else "raw",
                           "stage": s, "draw": b}
                    sg = T.summaries(g)
                    rec["gudhi_H0_total"], rec["gudhi_H1_total"] = sg["H0_total_persistence"], sg["H1_total_persistence"]
                    for dim in (0, 1):
                        rec[f"W1_H{dim}"] = T.wasserstein(g[dim], r[dim])
                        rec[f"n_gudhi_H{dim}"], rec[f"n_ripser_H{dim}"] = nfeat(g[dim]), nfeat(r[dim])
                        lg = (g[dim][:, 1] - g[dim][:, 0]).sum() if len(g[dim]) else 0.0
                        lr = (r[dim][:, 1] - r[dim][:, 0]).sum() if len(r[dim]) else 0.0
                        rec[f"tot_gudhi_H{dim}"], rec[f"tot_ripser_H{dim}"] = float(lg), float(lr)
                        rec[f"maxabs_bar_H{dim}"] = (float(np.max(np.abs(np.sort(g[dim], axis=0) - np.sort(r[dim], axis=0))))
                                                     if len(g[dim]) == len(r[dim]) and len(g[dim]) else np.nan)
                    rows.append(rec)
        print(f"[D6] {version} done ({time.time() - t0:.0f} s)", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "gudhi_vs_ripser_per_subsample.csv"), index=False)

    # reproduction check against the saved step 4 summaries (same draws?)
    repro = {}
    for version in ("unadjusted", "adjusted"):
        s4 = pd.read_csv(os.path.join(ROOT, "results", "step04", f"stage_topology_bootstrap_{version}.csv"))
        s4["mode"] = s4["mode"].map({"raw": "raw", "norm": "shape-only"})
        mrg = res[res.version == version].merge(s4, on=["stage", "mode", "draw"])
        repro[version] = (len(mrg), float(np.max(np.abs(mrg.gudhi_H0_total - mrg.H0_total_persistence))),
                          float(np.max(np.abs(mrg.gudhi_H1_total - mrg.H1_total_persistence))))

    # typical scale for context: median within-stage GUDHI-GUDHI distance is not needed; use diagram size
    res["rel_tot_H0"] = (res.tot_ripser_H0 - res.tot_gudhi_H0).abs() / res.tot_gudhi_H0.clip(lower=EPS)
    res["rel_tot_H1"] = np.where(res.tot_gudhi_H1 > EPS,
                                 (res.tot_ripser_H1 - res.tot_gudhi_H1).abs() / res.tot_gudhi_H1.clip(lower=EPS), 0.0)
    g = res.groupby(["version", "mode"])
    tab = pd.DataFrame({
        "subsamples": g.size(),
        "max_W1_H0": g.W1_H0.max(), "max_W1_H1": g.W1_H1.max(),
        "median_W1_H0+H1": g.apply(lambda x: (x.W1_H0 + x.W1_H1).median()),
        "identical_counts_H0_%": g.apply(lambda x: 100 * (x.n_gudhi_H0 == x.n_ripser_H0).mean()),
        "identical_counts_H1_%": g.apply(lambda x: 100 * (x.n_gudhi_H1 == x.n_ripser_H1).mean()),
        "max_rel_diff_total_H0": g.rel_tot_H0.max(), "max_rel_diff_total_H1": g.rel_tot_H1.max(),
        "max_abs_bar_diff": g.apply(lambda x: np.nanmax(x[["maxabs_bar_H0", "maxabs_bar_H1"]].values)),
    })
    tab.to_csv(os.path.join(OUT, "gudhi_vs_ripser_summary.csv"))
    # context: typical within-stage GUDHI W1 (H0+H1) from step 4 TSI would be larger by orders of magnitude
    worst = res.assign(W=res.W1_H0 + res.W1_H1).sort_values("W", ascending=False).head(5)
    L = ["DEVIATION D6 - Software cross-check: GUDHI vs Ripser (GSE41258 set main, step 4 subsamples)",
         "=" * 70, "CROSS-CHECK ONLY: no result, test or decision uses Ripser output (CLAUDE.md: GUDHI for all TDA).",
         "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "ripser", "numpy", "scipy")),
         f"Seed {SEED}; draws replayed exactly as step 4 (100 per stage x 4 stages x raw/shape-only x 2 versions"
         f" = {len(res)} subsample pairs of 31 samples)", "",
         "Reproduction of saved step 4 GUDHI summaries (same draws):"]
    for v, (n, d0, d1) in repro.items():
        L.append(f"  {v:<10} {n} subsamples matched; max |H0 total diff| {d0:.2e}, max |H1 total diff| {d1:.2e}")
    L += ["", "GUDHI vs Ripser, per (version, mode):", tab.to_string(float_format=lambda x: f"{x:.3g}"), "",
          "Five largest GUDHI-Ripser distances (W1 H0 + H1):"]
    for _, r in worst.iterrows():
        L.append(f"  {r.version} {r['mode']} {r.stage} draw {r.draw}: W1 H0 {r.W1_H0:.2e}, H1 {r.W1_H1:.2e}; "
                 f"H1 bars gudhi/ripser {r.n_gudhi_H1}/{r.n_ripser_H1}")
    L += ["", "Note: Ripser stores filtration values in single precision (float32), GUDHI in double precision;",
          "differences of order 1e-6 relative to the distance scale are expected from this alone.",
          f"Runtime: {time.time() - t0:.0f} s",
          "Outputs: results/step04/ripser_crosscheck/gudhi_vs_ripser_per_subsample.csv, gudhi_vs_ripser_summary.csv"]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
