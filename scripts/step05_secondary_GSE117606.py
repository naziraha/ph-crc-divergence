"""
GSE117606 step 5 SECONDARY analyses (addendum 6 commit 93613bd; deviations D11, D12).

All unpaired tests: tda_crc_pilot.permutation_test unchanged (VR H0+H1, exact W1, k = 20,
m = max(10, 80% of the smallest stage of the set), aligned, 999 permutations), Holm across the
2 transitions per analysis x version x mode.
  stageaware       stage-aware adjustment (addendum 2 section 1), main set
  noQAflag         without the 13 quality-flagged arrays (addendum 6 section 7.2), unadjusted / adjusted
  noGSM3304812     without GSM3304812 (D11), unadjusted / adjusted
  adjqual          stage-blind adjustment with composition + 3 quality covariates (D12), then frozen step 3
  S5_pc5, S6_pc20, S7_genes1000, S8_genes5000   unadjusted / adjusted / stage-aware (addendum 1 section 4)
Paired within-patient test (addendum 6 section 6): patients with both arrays (N->A 51, A->C 56);
each of k = 20 subsamples draws m = int(0.8 x n_pairs) patients and includes both arrays of each;
groups aligned once as in permutation_test; 999 permutations swapping the two labels within patients
(probability 1/2 each); persistence and diag_dist are the unchanged reference functions.
Resumable (one checkpoint per unit); per-unit seeds SeedSequence(42, crc32(unit name)).

Usage: python scripts/step05_secondary_GSE117606.py [run|summary]
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

GSE, SEED, K, NPERM, DIMS, WORKERS = "GSE117606", 42, 20, 999, (0, 1), 6
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
TRANS = [("Normal", "Adenoma"), ("Adenoma", "Carcinoma")]
PROC = os.path.join(ROOT, "data", "processed")
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "results", "step05", GSE)
CKPT = os.path.join(OUT, "checkpoints_secondary")
EXP = {"main": {"Normal": 65, "Adenoma": 59, "Carcinoma": 71},
       "main_noQAflag": {"Normal": 61, "Adenoma": 53, "Carcinoma": 68},
       "main_noGSM3304812": {"Normal": 65, "Adenoma": 59, "Carcinoma": 70}}
# analysis -> list of (embedding file stem, expected-count key, version label)
ANALYSES = {
    "stageaware": [("main_stageaware_g2000_pc10", "main", "stageaware")],
    "noQAflag": [(f"main_noQAflag_{v}_g2000_pc10", "main_noQAflag", v) for v in ("unadjusted", "adjusted")],
    "noGSM3304812": [(f"main_noGSM3304812_{v}_g2000_pc10", "main_noGSM3304812", v) for v in ("unadjusted", "adjusted")],
    "adjqual": [("main_adjqual_g2000_pc10", "main", "adjusted+quality")],
    **{a: [(f"main_{v}_{tag}", "main", v) for v in ("unadjusted", "adjusted", "stageaware")]
       for a, tag in (("S5_pc5", "g2000_pc5"), ("S6_pc20", "g2000_pc20"), ("S7_genes1000", "g1000_pc10"),
                      ("S8_genes5000", "g5000_pc10"))},
    "paired": [(f"main_{v}_g2000_pc10", "main", v) for v in ("unadjusted", "adjusted")],
}


def build_adjqual():
    """D12: stage-blind adjustment with composition + quality covariates, then frozen step 3."""
    path = os.path.join(EMB, f"{GSE}_main_adjqual_g2000_pc10.csv")
    if os.path.exists(path):
        return
    lab = pd.read_csv(os.path.join(PROC, f"{GSE}_main_samples_qc.csv"), index_col=0)
    g = pd.read_csv(os.path.join(PROC, f"{GSE}_main_expr_genes.csv"), index_col=0)[lab.index]
    sc = pd.read_csv(os.path.join(ROOT, "results", "step02", GSE, "composition_scores_all.csv"), index_col=0)
    qa = pd.read_csv(os.path.join(ROOT, "results", "step01", GSE, "array_quality.csv"), index_col=0)
    cov = pd.concat([sc.loc[lab.index, ["Epithelial", "Stromal", "Immune"]],
                     qa.loc[lab.index, ["mean_corr", "median", "IQR"]]], axis=1)
    C = (cov - cov.mean()).values
    X = np.c_[np.ones(len(cov)), C]
    Y = g.T.values
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    adj = pd.DataFrame((Y - C @ beta[1:]).T, index=g.index, columns=g.columns)
    adj.index.name = "gene"
    adj.to_csv(os.path.join(PROC, f"{GSE}_main_expr_genes_adjusted_quality.csv"), float_format="%.5f")
    emb, pca, _ = T.preprocess(adj, 2000, 10)
    e = emb.copy()
    e.columns = [f"PC{k + 1}" for k in range(10)]
    e.insert(0, "stage", lab.loc[e.index, "stage"])
    e.to_csv(path, index_label="sample", float_format="%.6f")
    json.dump({"var_explained_10PC": float(pca.explained_variance_ratio_.sum()),
               "covariates": list(cov.columns)}, open(os.path.join(OUT, "D12_adjqual_info.json"), "w"), indent=2)


def units():
    return [(a, stem, key, v, s, t, mo) for a, lst in ANALYSES.items() for stem, key, v in lst
            for s, t in TRANS for mo in ("norm", "raw")]


def unit_name(a, stem, key, v, s, t, mo):
    return f"{a}__{v.replace('+', '')[:9]}__{s[0]}-{t[0]}__{'shape' if mo == 'norm' else 'raw'}"


def paired_test(emb, ps, pt, m, normalize, rng):
    """ps[i], pt[i] = row indices of patient i's two arrays (stage s, stage t)."""
    emb = T.align_groups(emb, ps, pt, normalize)
    n = len(ps)

    def stat(a_idx, b_idx):
        ds, dt = [], []
        for _ in range(K):
            sel = rng.choice(n, size=m, replace=False)          # same m patients in both groups
            ds.append(T.persistence(emb[a_idx[sel]], normalize))
            dt.append(T.persistence(emb[b_idx[sel]], normalize))
        btw = [T.diag_dist(x, y, DIMS) for x in ds for y in dt]
        wth = [T.diag_dist(ds[i], ds[j], DIMS) for i in range(K) for j in range(i + 1, K)]
        wth += [T.diag_dist(dt[i], dt[j], DIMS) for i in range(K) for j in range(i + 1, K)]
        return float(np.mean(btw) - np.mean(wth))

    obs = stat(ps, pt)
    null = np.empty(NPERM)
    for i in range(NPERM):
        swap = rng.random(n) < 0.5
        null[i] = stat(np.where(swap, pt, ps), np.where(swap, ps, pt))
    p = (1 + np.sum(null >= obs)) / (1 + NPERM)
    z = (obs - null.mean()) / (null.std(ddof=1) + 1e-12)
    return obs, float(p), float(z)


def run_unit(u):
    a, stem, key, v, s, t, mo = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(os.path.join(EMB, f"{GSE}_{stem}.csv"), index_col=0)
    if df["stage"].value_counts().to_dict() != EXP[key]:
        raise SystemExit(f"Stage counts changed: {name} {df['stage'].value_counts().to_dict()}")
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    st = df["stage"].values
    if a == "paired":
        lab = pd.read_csv(os.path.join(PROC, f"{GSE}_main_samples_qc.csv"), index_col=0).loc[df.index]
        pos = {sid: i for i, sid in enumerate(df.index)}
        ps, pt = [], []
        for pid, grp in lab.groupby("patient id"):
            if (grp.stage == s).sum() == 1 and (grp.stage == t).sum() == 1:
                ps.append(pos[grp.index[grp.stage == s][0]]); pt.append(pos[grp.index[grp.stage == t][0]])
        ps, pt = np.array(ps), np.array(pt)
        m = int(0.8 * len(ps))
        Tobs, p, z = paired_test(emb, ps, pt, m, mo == "norm", rng)
        n1 = n2 = len(ps)
    else:
        m = max(10, int(0.8 * min((st == x).sum() for x in STAGES3)))
        idx_s, idx_t = np.where(st == s)[0], np.where(st == t)[0]
        Tobs, p, z, _ = T.permutation_test(emb, idx_s, idx_t, m, K, DIMS, mo == "norm", NPERM, rng, align=True)
        n1, n2 = len(idx_s), len(idx_t)
    rec = dict(unit=name, analysis=a, version=v, transition=f"{s}->{t}", mode="shape-only" if mo == "norm" else "raw",
               n_from=int(n1), n_to=int(n2), subsample_size=int(m), T_obs=float(Tobs), p=float(p), z=float(z),
               seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"z={z:.2f} p={p:.3f} ({rec['seconds']:.0f} s)"


def summary():
    r = pd.DataFrame([json.load(open(os.path.join(CKPT, f))) for f in os.listdir(CKPT) if f.endswith(".json")])
    prim = pd.read_csv(os.path.join(OUT, "permutation_results_all.csv"))
    prim = prim[prim.analysis == "primary"].assign(analysis="PRIMARY")
    r = pd.concat([prim[r.columns.intersection(prim.columns)], r], ignore_index=True)
    r["p_holm"] = np.nan
    for _, g in r.groupby(["analysis", "version", "mode"]):
        if len(g) == 2:
            r.loc[g.index, "p_holm"] = T.holm(g.p.values)
    r.to_csv(os.path.join(OUT, "secondary_results_all.csv"), index=False)
    order = ["PRIMARY"] + list(ANALYSES)
    L = [f"STEP 5 SECONDARY ANALYSES ({GSE})", "=" * 60,
         "Pre-specified: docs/addendum6_GSE117606.md (93613bd); deviations D11, D12 (docs/analysis_deviations.md).",
         f"999 permutations, k = {K}, aligned; Holm across the 2 transitions per analysis x version x mode.",
         f"Units finished: {len(r) - len(prim)}/{len(units())}",
         "Values: z (Holm p); * = Holm p < 0.05.", "",
         f"{'analysis':<15}{'version':<17}{'n (N/A/C or pairs)':<20}{'N->A shape':>16}{'N->A raw':>16}{'A->C shape':>16}{'A->C raw':>16}"]

    def cell(g, tr, mode):
        q = g[(g.transition == tr) & (g["mode"] == mode)]
        if q.empty:
            return "–"
        z, ph = float(q.z.iloc[0]), float(q.p_holm.iloc[0])
        return f"{z:5.2f}{'*' if ph < 0.05 else ' '} ({ph:.3f})"

    for a in order:
        for v in [x for x in ("unadjusted", "adjusted", "stageaware", "adjusted+quality") if ((r.analysis == a) & (r.version == x)).any()]:
            g = r[(r.analysis == a) & (r.version == v)]
            na = g[g.transition == "Normal->Adenoma"].iloc[0] if (g.transition == "Normal->Adenoma").any() else None
            ac = g[g.transition == "Adenoma->Carcinoma"].iloc[0] if (g.transition == "Adenoma->Carcinoma").any() else None
            nstr = (f"{na.n_from}/{na.n_to} | {ac.n_from}/{ac.n_to}" if na is not None and ac is not None else "")
            L.append(f"{a:<15}{v:<17}{nstr:<20}{cell(g, 'Normal->Adenoma', 'shape-only'):>16}{cell(g, 'Normal->Adenoma', 'raw'):>16}"
                     f"{cell(g, 'Adenoma->Carcinoma', 'shape-only'):>16}{cell(g, 'Adenoma->Carcinoma', 'raw'):>16}")
    L += ["", "Outputs: results/step05/GSE117606/secondary_results_all.csv, checkpoints_secondary/"]
    open(os.path.join(OUT, "summary_secondary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    os.makedirs(CKPT, exist_ok=True)
    if (sys.argv[1] if len(sys.argv) > 1 else "run") == "run":
        build_adjqual()
        todo = [u for u in units() if not os.path.exists(os.path.join(CKPT, unit_name(*u) + ".json"))]
        print(f"[sec] {len(todo)} units to run", flush=True)
        t0 = time.time()
        with Pool(WORKERS) as pool:
            for i, (n, msg) in enumerate(pool.imap_unordered(run_unit, todo), 1):
                print(f"[sec] {i}/{len(todo)} {n}: {msg} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summary()


if __name__ == "__main__":
    main()
