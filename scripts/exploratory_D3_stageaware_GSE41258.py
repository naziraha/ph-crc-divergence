"""
EXPLORATORY analysis D3 (docs/analysis_deviations.md) - GSE41258 stage-aware composition adjustment.

Not pre-specified; cannot change the primary result or the tipping-point definition.
Same model as addendum 2 section 1, with the liver fraction added:
  per gene OLS: expression ~ 1 + stage indicators (Adenoma, Carcinoma, Metastasis; Normal ref)
                + Epithelial + Stromal + Immune + Liver_fraction   (step-2 scores, set main)
  stage-aware adjusted = expression - sum b_c (c - mean c)   (only the composition part removed)
Then the frozen steps 3 and 5 for set `main`: liver-signature genes removed before gene
selection (as in the primary adjusted analysis), 2,000 genes, 10 PCs, tda_crc_pilot.preprocess;
pre-specified permutation test (k = 20, m = 31, aligned, 999 permutations), Holm across the
three transitions per mode. Resumable: one checkpoint per unit; per-unit seeds from 42.

Usage: python scripts/exploratory_D3_stageaware_GSE41258.py
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
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "exploratory_D3")
CKPT = os.path.join(OUT, "checkpoints")
EMB_FILE = os.path.join(ROOT, "results", "step03", "embeddings", "GSE41258_main_stageaware_g2000_pc10.csv")
EXPECTED = {"Normal": 53, "Adenoma": 48, "Carcinoma": 181, "Metastasis": 39}
COV = ["Epithelial", "Stromal", "Immune", "Liver_fraction"]
TRANS = list(zip(T.STAGES[:-1], T.STAGES[1:]))


def build_embedding():
    lab = pd.read_csv(os.path.join(PROC, "GSE41258_main_samples_qc.csv"), index_col=0)
    if lab["stage"].value_counts().to_dict() != EXPECTED:
        raise SystemExit("Stage counts changed")
    g = pd.read_csv(os.path.join(PROC, "GSE41258_main_expr_genes.csv"), index_col=0)[lab.index]
    sc = pd.read_csv(os.path.join(ROOT, "results", "step02", "composition_scores_main.csv"), index_col=0).loc[lab.index]
    C = (sc[COV] - sc[COV].mean()).values
    st = lab["stage"].values
    X = np.c_[np.ones(len(st)), (st == "Adenoma"), (st == "Carcinoma"), (st == "Metastasis"), C].astype(float)
    Y = g.T.values
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    comp = C @ beta[-len(COV):]
    adj = pd.DataFrame((Y - comp).T, index=g.index, columns=g.columns)
    removed = pd.Series(comp.var(axis=0) / Y.var(axis=0), index=g.index)
    adj.index.name = "gene"
    adj.to_csv(os.path.join(PROC, "GSE41258_main_expr_genes_adjusted_stageaware.csv"), float_format="%.5f")
    liver = set(pd.read_csv(os.path.join(ROOT, "results", "step02", "liver_signature_main.csv"))["gene"])
    adj = adj.drop(index=[x for x in adj.index if x in liver])
    emb, pca, n = T.preprocess(adj, 2000, 10)
    e = emb.copy()
    e.columns = [f"PC{k + 1}" for k in range(10)]
    e.insert(0, "stage", lab.loc[e.index, "stage"])
    e.to_csv(EMB_FILE, index_label="sample", float_format="%.6f")
    return dict(var_removed_median=float(removed.median()), genes_var_removed_gt_half=int((removed > .5).sum()),
                pca_var=float(pca.explained_variance_ratio_.sum()), liver_genes_removed=len(liver & set(g.index)))


def unit_name(s, t, mode):
    return f"D3_stageaware__{s[0]}-{t[0]}__{'shape' if mode == 'norm' else 'raw'}"


def run_unit(u):
    s, t, mode = u
    name = unit_name(*u)
    path = os.path.join(CKPT, name + ".json")
    if os.path.exists(path):
        return name, "skipped"
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))
    df = pd.read_csv(EMB_FILE, index_col=0)
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    stv = df["stage"].values
    m = max(10, int(0.8 * min((stv == x).sum() for x in T.STAGES)))
    idx_s, idx_t = np.where(stv == s)[0], np.where(stv == t)[0]
    Tobs, p, z, null = T.permutation_test(emb, idx_s, idx_t, m, K_DRAWS, DIMS, mode == "norm", N_PERM, rng, align=True)
    rec = dict(unit=name, version="stageaware (exploratory D3)", transition=f"{s}->{t}",
               mode="shape-only" if mode == "norm" else "raw", n_from=int(len(idx_s)), n_to=int(len(idx_t)),
               subsample_size=int(m), T_obs=float(Tobs), p=float(p), z=float(z), null=[float(x) for x in null],
               seconds=round(time.time() - t0, 1))
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return name, f"done in {rec['seconds']:.0f} s (p={p:.4f}, z={z:.1f})"


def main():
    os.makedirs(CKPT, exist_ok=True)
    info = build_embedding() if not os.path.exists(EMB_FILE) else None
    if info:
        with open(os.path.join(OUT, "adjustment_info.json"), "w") as fh:
            json.dump(info, fh, indent=2)
    info = json.load(open(os.path.join(OUT, "adjustment_info.json")))
    todo = [(s, t, mo) for s, t in TRANS for mo in ("norm", "raw")
            if not os.path.exists(os.path.join(CKPT, unit_name(s, t, mo) + ".json"))]
    print(f"[D3] {len(todo)} units to run", flush=True)
    with Pool(WORKERS) as pool:
        for n, msg in pool.imap_unordered(run_unit, todo):
            print(f"[D3] {n}: {msg}", flush=True)
    rows = []
    for f in sorted(os.listdir(CKPT)):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(CKPT, f)))
            r.pop("null")
            rows.append(r)
    df = pd.DataFrame(rows)
    df["p_holm"] = np.nan
    for _, gg in df.groupby("mode"):
        df.loc[gg.index, "p_holm"] = T.holm(gg["p"].values)
    df["_o"] = df.transition.map({f"{s}->{t}": i for i, (s, t) in enumerate(TRANS)})
    df = df.sort_values(["mode", "_o"], ascending=[False, True]).drop(columns="_o")
    df.to_csv(os.path.join(OUT, "permutation_results_D3.csv"), index=False)
    prim = pd.read_csv(os.path.join(ROOT, "results", "step05", "permutation_results_all.csv"))
    prim = prim[prim.analysis == "primary"]
    L = ["EXPLORATORY ANALYSIS D3 - GSE41258 stage-aware composition adjustment", "=" * 60,
         "NOT pre-specified (docs/analysis_deviations.md, D3). Cannot change the primary result or the",
         "tipping-point definition. Set main: " + ", ".join(f"{k} {v}" for k, v in EXPECTED.items()),
         "Model: per gene OLS ~ stage indicators + Epithelial + Stromal + Immune + Liver_fraction;",
         "  only the composition part removed. Liver-signature genes removed before gene selection.",
         f"  composition variance removed per gene: median {info['var_removed_median']:.3f}; genes > 0.5: "
         f"{info['genes_var_removed_gt_half']}; liver genes removed {info['liver_genes_removed']}; "
         f"10 PCs explain {info['pca_var'] * 100:.1f}%",
         f"Test: frozen step 5 (k={K_DRAWS}, m=31, aligned, {N_PERM} permutations, Holm across 3 transitions).",
         f"Seed {SEED}; " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy")), "",
         "Results (with the pre-specified primary results for comparison):"]
    for _, r in df.iterrows():
        ref = prim[(prim.transition == r.transition) & (prim["mode"] == r["mode"])]
        pu = ref[ref.version == "unadjusted"].iloc[0]
        pa = ref[ref.version == "adjusted"].iloc[0]
        L.append(f"  {r['mode']:<10} {r.transition:<24} stage-aware z={r.z:5.2f} p_Holm {r.p_holm:.4f} | "
                 f"unadjusted z={pu.z:5.2f} p_Holm {pu.p_holm:.4f} | adjusted z={pa.z:5.2f} p_Holm {pa.p_holm:.4f}")
    L += ["", "Outputs: results/exploratory_D3/permutation_results_D3.csv, checkpoints/;",
          "  results/step03/embeddings/GSE41258_main_stageaware_g2000_pc10.csv;",
          "  data/processed/GSE41258_main_expr_genes_adjusted_stageaware.csv"]
    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
