"""
Addendum 12 (73b4e71; D25 4018ea9) - quantitative persistence outputs.

Test (per transition / per validation replicate): groups aligned (tda_crc_pilot.align_groups); k = 20 subsamples of
size m per group; VR diagrams (tda_crc_pilot.persistence, H0 and H1, infinite bars removed); statistics computed from
the same draws (D25b):
  landscape distance (H0, H1): L2 distance between mean landscapes (gudhi.representations.Landscape, 5 landscapes,
                                100 points, grid 0..max finite death of the pooled observed diagrams; D25c)
  barcode summaries: H0 total / max / entropy; H1 total / max / n features / entropy -> |mean_s - mean_t|
Patient-level permutation of the aligned labels (group sizes kept), k subsamples redrawn per permuted group;
one-sided p = (1 + #null >= obs) / (1 + n_perm).

Usage: python scripts/addendum12_ph_outputs.py validate
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
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import tda_crc_pilot as T  # noqa: E402

SEED, K, WORKERS = 42, 20, 6
NL, RES = 5, 100
OUT = os.path.join(ROOT, "results", "addendum12")
SUMMARIES = ["H0_total_persistence", "H0_max_persistence", "H0_entropy",
             "H1_total_persistence", "H1_max_persistence", "H1_n_features", "H1_entropy"]
STATS = ["landscape_H0", "landscape_H1"] + SUMMARIES


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def draw(emb, idx, m, normalize, rng):
    return [T.persistence(emb[rng.choice(idx, size=m, replace=False)], normalize) for _ in range(K)]


def landscape_vecs(diags, dim, rmax):
    from gudhi.representations import Landscape
    L = Landscape(num_landscapes=NL, resolution=RES, sample_range=[0.0, rmax])
    out = []
    for d in diags:
        x = d[dim]
        out.append(np.zeros(NL * RES) if len(x) == 0 else L.fit_transform([x])[0])
    return np.array(out)


def group_stats(dA, dB, rmax, land_dims=(0, 1)):
    s = {}
    for dim in land_dims:
        a, b = landscape_vecs(dA, dim, rmax[dim]), landscape_vecs(dB, dim, rmax[dim])
        s[f"landscape_H{dim}"] = float(np.linalg.norm(a.mean(0) - b.mean(0)))
    sa = pd.DataFrame([T.summaries(d) for d in dA])
    sb = pd.DataFrame([T.summaries(d) for d in dB])
    for c in SUMMARIES:
        s[c] = float(abs(sa[c].mean() - sb[c].mean()))
    return s


def multi_test(emb, idx_s, idx_t, m, normalize, n_perm, rng, land_dims=(0, 1)):
    """All addendum 12 statistics with one shared patient-level permutation scheme."""
    stats = [f"landscape_H{d}" for d in land_dims] + SUMMARIES
    emb = T.align_groups(emb, idx_s, idx_t, normalize)
    dA, dB = draw(emb, idx_s, m, normalize, rng), draw(emb, idx_t, m, normalize, rng)
    rmax = {}
    for dim in (0, 1):
        deaths = [d[dim][:, 1].max() for d in dA + dB if len(d[dim])]
        rmax[dim] = float(max(deaths)) if deaths else 1.0
    obs = group_stats(dA, dB, rmax, land_dims)
    pool, ns = np.concatenate([idx_s, idx_t]), len(idx_s)
    null = {k: [] for k in stats}
    for _ in range(n_perm):
        perm = rng.permutation(pool)
        pa, pb = draw(emb, perm[:ns], m, normalize, rng), draw(emb, perm[ns:], m, normalize, rng)
        st = group_stats(pa, pb, rmax, land_dims)
        for k in stats:
            null[k].append(st[k])
    res = {}
    for k in stats:
        nl = np.array(null[k])
        res[k] = dict(obs=obs[k], p=float((1 + np.sum(nl >= obs[k])) / (1 + n_perm)),
                      z=float((obs[k] - nl.mean()) / (nl.std(ddof=1) + 1e-12)))
    return res


# ================================================================ section 6 validation
VAL_COND = [("N0", 0), ("SP", 1.25), ("CL", 1.0), ("LP", 1.0)]
NA_, NB_, M_VAL, NPERM_VAL, R_VAL = 53, 48, 31, 99, 20


def val_unit(u):
    cond, x, rep = u
    name = f"a12val__{cond}{x:g}__r{rep:02d}"
    ck = os.path.join(OUT, "validation_checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    from addendum5_part3_controls import insert_level1
    rng = rng_for(name)
    df = pd.read_csv(os.path.join(ROOT, "results", "step03", "embeddings", "GSE41258_main_unadjusted_g2000_pc10.csv"),
                     index_col=0)
    carc = df.index[df.stage == "Carcinoma"]
    assert len(carc) == 181
    X = df.loc[carc, [c for c in df.columns if c.startswith("PC")]].values.copy()
    perm = rng.permutation(len(carc))
    ia, ib = perm[:NA_], perm[NA_:NA_ + NB_]
    X[ib] = insert_level1(X[ib], cond, x, rng)
    emb = np.r_[X[ia], X[ib]]
    gA, gB = np.arange(NA_), np.arange(NA_, NA_ + NB_)
    rec = dict(condition=cond, rep=rep)
    for mode in ("shape", "raw"):
        r = multi_test(emb, gA, gB, M_VAL, mode == "shape", NPERM_VAL, rng)
        for k, v in r.items():
            rec[f"{mode}__{k}__p"] = v["p"]
            rec[f"{mode}__{k}__z"] = v["z"]
    json.dump(rec, open(ck, "w"))
    return rec


def validate():
    os.makedirs(os.path.join(OUT, "validation_checkpoints"), exist_ok=True)
    t0 = time.time()
    units = [(c, x, r) for c, x in VAL_COND for r in range(R_VAL)]
    with Pool(WORKERS) as pool:
        recs = []
        for i, rec in enumerate(pool.imap_unordered(val_unit, units), 1):
            recs.append(rec)
            if i % 10 == 0:
                print(f"[a12 val] {i}/{len(units)} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    d = pd.DataFrame(recs)
    d.to_csv(os.path.join(OUT, "validation_replicates.csv"), index=False)
    det = lambda cond, col: int((d[d.condition == cond][col + "__p"] <= 0.05).sum())
    crit = []
    for mode in ("shape", "raw"):
        for dim in ("H0", "H1"):
            v = det("N0", f"{mode}__landscape_{dim}")
            crit.append((f"N0 landscape {dim} {mode}", "<= 3/20", v, v <= 3))
    for dim in ("H0", "H1"):
        v = det("SP", f"shape__landscape_{dim}")
        crit.append((f"SP landscape {dim} shape-only (ignore pure spread)", "<= 3/20", v, v <= 3))
    v = det("SP", "raw__landscape_H0")
    crit.append(("SP landscape H0 raw (D25a)", ">= 16/20", v, v >= 16))
    v = det("CL", "shape__landscape_H0")
    crit.append(("CL landscape H0 shape-only", ">= 16/20", v, v >= 16))
    v = det("LP", "shape__landscape_H1")
    crit.append(("LP landscape H1 shape-only", ">= 16/20", v, v >= 16))
    for s in SUMMARIES:
        v = det("N0", f"shape__{s}")
        crit.append((f"N0 barcode summary {s} shape-only", "<= 3/20", v, v <= 3))
    cdf = pd.DataFrame(crit, columns=["check", "requirement", "detected", "pass"])
    cdf.to_csv(os.path.join(OUT, "validation_criteria.csv"), index=False)
    ok = bool(cdf["pass"].all())
    tab = []
    for cond, _ in VAL_COND:
        for mode in ("shape", "raw"):
            tab.append(f"  {cond:<3} {mode:<6} " + "  ".join(f"{k.replace('_persistence', '').replace('landscape_', 'L_')}:{det(cond, f'{mode}__{k}'):>2}"
                                                       for k in STATS))
    L = ["ADDENDUM 12 section 6 - validation of landscape and barcode-summary permutation tests (D25)", "=" * 72,
         "GSE41258 carcinomas split A 53 / B 48; insertion into B (addendum 5 level 1); k = 20, m = 31, aligned;",
         "99 permutations; 20 replicates per condition; detection = p <= 0.05; shared draws for all statistics (D25b).",
         f"Seed {SEED} per replicate. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "pandas")), ""]
    L += [f"  {'PASS' if r['pass'] else 'FAIL'}  {r.check:<52} {r.requirement:<9} {r.detected}/20" for _, r in cdf.iterrows()]
    L += ["", "All detection counts (/20) by condition, mode and statistic:"] + tab
    L += ["", f"OVERALL: {'ALL CRITERIA PASSED - real-data analysis may proceed' if ok else 'FAILED - do not analyse real data'}",
          f"Run time {(time.time() - t0) / 60:.1f} min."]
    open(os.path.join(OUT, "validation_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ================================================================ real data (addendum 12 + Amendment 1)
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
GEO_ST = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
TC_ST = ["Stage I", "Stage II", "Stage III", "Stage IV"]
COH12 = {"GSE41258": ("GSE41258_main", GEO_ST, [53, 48, 181, 39]),
         "GSE20916": ("GSE20916_bulk", GEO_ST[:3], [24, 45, 36]),
         "GSE37364": ("GSE37364_bulk", GEO_ST[:3], [38, 29, 27]),
         "GSE117606": ("GSE117606_main", GEO_ST[:3], [65, 59, 71]),
         "TCGA_COAD": ("TCGA_COAD_main", TC_ST, [45, 110, 80, 39])}
VERS12 = ("unadjusted", "adjusted")
NPERM_REAL = 999


def load12(c, v):
    stem, stages, exp = COH12[c]
    df = pd.read_csv(os.path.join(EMB, f"{stem}_{v}_g2000_pc10.csv"), index_col=0)
    if [int((df.stage == s).sum()) for s in stages] != exp:
        raise SystemExit(f"Stage counts changed: {c} {v}")
    return df[[x for x in df.columns if x.startswith("PC")]].values, df.stage.values, stages, max(10, int(0.8 * min(exp)))


def short(s):
    return s.replace("Stage ", "") if s.startswith("Stage") else s[0]


def real_unit(u):
    c, v, mode, s, t = u
    name = f"a12__{c}__{v[:5]}__{mode}__{short(s)}-{short(t)}"
    ck = os.path.join(OUT, "checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    t0 = time.time()
    emb, st, stages, m = load12(c, v)
    r = multi_test(emb, np.where(st == s)[0], np.where(st == t)[0], m, mode == "shape", NPERM_REAL, rng_for(name),
                   land_dims=(1,))                                  # Amendment 1: landscape test H1 only
    rec = dict(cohort=c, version=v, mode=mode, transition=f"{s}->{t}", m=m, seconds=round(time.time() - t0))
    for k, x in r.items():
        rec[f"{k}__obs"], rec[f"{k}__p"], rec[f"{k}__z"] = x["obs"], x["p"], x["z"]
    json.dump(rec, open(ck, "w"))
    return rec


def dims_unit(u):
    v, mode, d, s, t = u
    name = f"a12dims__GSE117606__{v[:5]}__{mode}__H{d}__{s[0]}-{t[0]}"
    ck = os.path.join(OUT, "checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    emb, st, _, m = load12("GSE117606", v)
    Tobs, p, z, _ = T.permutation_test(emb, np.where(st == s)[0], np.where(st == t)[0], m, K, (d,), mode == "shape",
                                       NPERM_REAL, rng_for(name), align=True)
    rec = dict(cohort="GSE117606", version=v, mode="shape-only" if mode == "shape" else "raw", dim=f"H{d}",
               transition=f"{s}->{t}", z=float(z), p=float(p))
    json.dump(rec, open(ck, "w"))
    return rec


def real_units():
    u = []
    for c, (_, stages, _) in COH12.items():
        for v in VERS12:
            for mode in ("shape", "raw"):
                for s, t in zip(stages[:-1], stages[1:]):
                    u.append((c, v, mode, s, t))
    return u


def dims_units():
    return [(v, mode, d, s, t) for v in VERS12 for mode in ("shape", "raw") for d in (0, 1)
            for s, t in (("Normal", "Adenoma"), ("Adenoma", "Carcinoma"))]


def figures(c, v):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gudhi.representations import Landscape
    os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
    emb, st, stages, m = load12(c, v)
    colors = dict(zip(stages, ["#2a78d6", "#e0a100", "#c0392b", "#7b4fb3"]))
    for mode in ("shape", "raw"):
        rng = rng_for(f"a12fig__{c}__{v}__{mode}")
        D = {s: draw(emb, np.where(st == s)[0], m, mode == "shape", rng) for s in stages}
        tag = f"{c}_{v}_{mode}"
        # diagrams: representative + overlay
        fig, axes = plt.subplots(2, 2, figsize=(10, 9))
        for j, dim in enumerate((0, 1)):
            hi = max([d[dim][:, 1].max() for s in stages for d in D[s] if len(d[dim])] + [1e-9]) * 1.05
            for i, (lab, rep) in enumerate((("representative subsample", True), ("all 20 subsamples", False))):
                ax = axes[i, j]
                for s in stages:
                    dl = D[s][:1] if rep else D[s]
                    pts = np.vstack([d[dim] for d in dl if len(d[dim])] or [np.zeros((0, 2))])
                    ax.scatter(pts[:, 0], pts[:, 1], s=10 if rep else 4, alpha=0.9 if rep else 0.3, color=colors[s], label=s)
                ax.plot([0, hi], [0, hi], color="0.5", lw=0.8)
                ax.set(xlim=(0, hi), ylim=(0, hi), xlabel="birth", ylabel="death", title=f"H{dim} diagram, {lab}")
            axes[0, j].legend(fontsize=7)
        fig.suptitle(f"{c} {v} {'shape-only' if mode == 'shape' else 'raw'}: persistence diagrams (descriptive)")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "figures", f"{tag}_diagrams.png"), dpi=130); plt.close(fig)
        # barcodes (representative subsample)
        fig, axes = plt.subplots(len(stages), 2, figsize=(10, 2.4 * len(stages)), squeeze=False)
        for i, s in enumerate(stages):
            for j, dim in enumerate((0, 1)):
                ax = axes[i, j]
                bars = D[s][0][dim]
                bars = bars[np.argsort(bars[:, 1] - bars[:, 0])] if len(bars) else bars
                for k, (b, e) in enumerate(bars):
                    ax.plot([b, e], [k, k], color=colors[s], lw=1.5)
                ax.set(title=f"{s}: H{dim} barcode ({len(bars)} bars)", xlabel="filtration scale", yticks=[])
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "figures", f"{tag}_barcodes.png"), dpi=130); plt.close(fig)
        # landscapes (lambda_1) and Betti curves
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        for j, dim in enumerate((0, 1)):
            rmax = max([d[dim][:, 1].max() for s in stages for d in D[s] if len(d[dim])] + [1e-9])
            grid = np.linspace(0, rmax, RES)
            L = Landscape(num_landscapes=1, resolution=RES, sample_range=[0.0, rmax])
            for s in stages:
                V = np.array([np.zeros(RES) if len(d[dim]) == 0 else L.fit_transform([d[dim]])[0] for d in D[s]])
                axes[0, j].plot(grid, V.mean(0), color=colors[s], lw=2, label=s)
                axes[0, j].fill_between(grid, np.percentile(V, 2.5, 0), np.percentile(V, 97.5, 0), color=colors[s], alpha=0.12)
                B = np.array([T.betti_curve(d[dim], grid) for d in D[s]])
                axes[1, j].plot(grid, B.mean(0), color=colors[s], lw=2, label=s)
                axes[1, j].fill_between(grid, np.percentile(B, 25, 0), np.percentile(B, 75, 0), color=colors[s], alpha=0.12)
            axes[0, j].set(title=f"H{dim} landscape λ1" + (" (figure only)" if dim == 0 else "") + "\nmean; band = subsample variability",
                           xlabel="filtration scale")
            axes[1, j].set(title=f"Betti-{dim} curve (mean, IQR)", xlabel="filtration scale")
            axes[0, j].legend(fontsize=7)
        fig.suptitle(f"{c} {v} {'shape-only' if mode == 'shape' else 'raw'}: landscapes and Betti curves (descriptive)")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "figures", f"{tag}_landscapes_betti.png"), dpi=130); plt.close(fig)


def holm(p):
    return T.holm(np.asarray(p, float))


def summarise():
    recs = [json.load(open(os.path.join(OUT, "checkpoints", f))) for f in os.listdir(os.path.join(OUT, "checkpoints"))
            if f.startswith("a12__")]
    d = pd.DataFrame(recs)
    rows = []
    for _, r in d.iterrows():
        for k in ["landscape_H1"] + SUMMARIES:
            rows.append(dict(cohort=r.cohort, version=r.version, mode="shape-only" if r["mode"] == "shape" else "raw",
                             transition=r.transition, statistic=k, obs=r[f"{k}__obs"], z=r[f"{k}__z"], p=r[f"{k}__p"]))
    L_ = pd.DataFrame(rows)
    L_["p_holm_trans"] = np.nan
    for _, g in L_.groupby(["cohort", "version", "mode", "statistic"]):
        L_.loc[g.index, "p_holm_trans"] = holm(g.p.values)
    land = L_[L_.statistic == "landscape_H1"].copy()
    summ = L_[L_.statistic != "landscape_H1"].copy()
    summ["p_holm_both"] = np.nan
    for _, g in summ.groupby(["cohort", "version", "mode", "transition"]):
        summ.loc[g.index, "p_holm_both"] = holm(g.p_holm_trans.values)
    land.to_csv(os.path.join(OUT, "landscape_tests.csv"), index=False)
    summ.to_csv(os.path.join(OUT, "barcode_summary_tests.csv"), index=False)
    # dimension split, all cohorts
    dd = pd.DataFrame([json.load(open(os.path.join(OUT, "checkpoints", f))) for f in os.listdir(os.path.join(OUT, "checkpoints"))
                       if f.startswith("a12dims__")])
    dd["p_holm"] = np.nan
    for _, g in dd.groupby(["version", "mode", "dim"]):
        dd.loc[g.index, "p_holm"] = holm(g.p.values)
    a5 = pd.read_csv(os.path.join(ROOT, "results", "addendum5", "part4_dims", "dims_results_all.csv"))
    a5 = a5[["cohort", "version", "mode", "dim", "transition", "z", "p_holm"]]
    tc = pd.read_csv(os.path.join(ROOT, "results", "addendum10", "spread_H0H1", "dims_results_all.csv"))
    tc = tc[tc.dim != "H0+H1"].assign(cohort="TCGA_COAD")[["cohort", "version", "mode", "dim", "transition", "z", "p_holm"]]
    dims_all = pd.concat([a5, tc, dd[["cohort", "version", "mode", "dim", "transition", "z", "p_holm"]]], ignore_index=True)
    # step 5 combined (primary)
    s5 = []
    for c in COH12:
        f = {"GSE41258": "permutation_results_all.csv", "TCGA_COAD": os.path.join("TCGA", "permutation_results_all.csv")}.get(
            c, os.path.join(c, "permutation_results_all.csv"))
        x = pd.read_csv(os.path.join(ROOT, "results", "step05", f))
        x = x[(x.analysis == "primary") & x.version.isin(VERS12)].assign(cohort=c, dim="H0+H1")
        s5.append(x[["cohort", "version", "mode", "dim", "transition", "z", "p_holm"]])
    s5 = pd.concat(s5, ignore_index=True)
    dims_all = pd.concat([s5, dims_all], ignore_index=True)
    dims_all.to_csv(os.path.join(OUT, "dimension_split_all_cohorts.csv"), index=False)
    # concordance (shape-only): step 5 vs H1 landscape
    con = land[land["mode"] == "shape-only"].merge(
        s5[s5["mode"] == "shape-only"][["cohort", "version", "transition", "p_holm"]].rename(columns={"p_holm": "step5_p_holm"}),
        on=["cohort", "version", "transition"], how="left")
    con["step5_sig"] = con.step5_p_holm < 0.05
    con["land_sig"] = con.p_holm_trans < 0.05
    con["agree"] = con.step5_sig == con.land_sig
    con.to_csv(os.path.join(OUT, "concordance_with_step5.csv"), index=False)
    agree = float(con.agree.mean())
    # dimension rules
    piv = dims_all.pivot_table(index=["cohort", "version", "mode", "transition"], columns="dim", values="p_holm").dropna()
    sig_comb = piv[piv["H0+H1"] < 0.05]
    h0_rep = float((sig_comb["H0"] < 0.05).mean()) if len(sig_comb) else np.nan
    h1_alone = int(((piv["H1"] < 0.05) & (piv["H0+H1"] >= 0.05)).sum())
    pred_h0 = (h0_rep >= 0.8) and h1_alone == 0
    geo = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
    h1_both = []
    for c in geo:
        for tr in piv.loc[c].index.get_level_values("transition").unique() if c in piv.index.get_level_values(0) else []:
            ok = all(((piv.loc[(c, v, "shape-only", tr), "H1"] < 0.05) if (c, v, "shape-only", tr) in piv.index else False) for v in VERS12)
            if ok:
                h1_both.append(c); break
    lands_both = []
    for c in geo:
        x = land[(land.cohort == c) & (land["mode"] == "shape-only")]
        for tr in x.transition.unique():
            if all(((x.version == v) & (x.transition == tr) & (x.p_holm_trans < 0.05)).any() for v in VERS12):
                lands_both.append(c); break
    loops = len(set(h1_both)) >= 3 or len(set(lands_both)) >= 3
    # summary text
    S = ["ADDENDUM 12 - quantitative persistence outputs (73b4e71; Amendment 1 d5fe9e6; D25, D26)", "=" * 72,
         "Secondary analyses; the step 5 tipping-point decision is unchanged. Patient-level permutation (999), aligned groups,",
         "k = 20 subsamples of m per stage (as step 5); shared draws for all statistics (D25b).",
         "Landscape test: H1 only (Amendment 1); H0 landscapes are figures only (validation: H0 landscape underpowered).",
         f"Seed {SEED} per unit. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "scipy")), "",
         "1. H1 LANDSCAPE TEST (L2 distance between mean landscapes; Holm across transitions):"]
    for _, r in land.sort_values(["cohort", "version", "mode", "transition"]).iterrows():
        S.append(f"  {r.cohort:<10} {r.version:<10} {r['mode']:<10} {r.transition:<22} z {r.z:6.2f}  p {r.p:.3f}  Holm {r.p_holm_trans:.3f}{'  *' if r.p_holm_trans < .05 else ''}")
    sig_s = summ[summ.p_holm_both < 0.05]
    S += ["", f"2. BARCODE SUMMARIES (Holm across transitions, then across 7 summaries): {len(sig_s)} significant of {len(summ)} tests"]
    for _, r in sig_s.sort_values(["cohort", "version", "mode", "transition"]).iterrows():
        S.append(f"  {r.cohort:<10} {r.version:<10} {r['mode']:<10} {r.transition:<22} {r.statistic:<22} z {r.z:6.2f}  final p {r.p_holm_both:.3f}")
    S += ["", "3. CONCORDANCE, shape-only: step 5 (H0+H1 Wasserstein) vs H1 landscape test",
          f"  agreement {100 * agree:.0f}% of {len(con)} comparisons -> {'SUPPORTIVE' if agree >= 0.8 else 'NOT >= 80% - discordances listed (expected lower; step 5 signal predominantly H0, Amendment 1)'}"]
    for _, r in con[~con.agree].iterrows():
        S.append(f"    discordant: {r.cohort} {r.version} {r.transition}: step 5 Holm {r.step5_p_holm:.3f}, H1 landscape Holm {r.p_holm_trans:.3f}")
    S += ["", "4. DIMENSION SPLIT, all cohorts (Wasserstein; H0+H1 = step 5)",
          f"  significant combined results: {len(sig_comb)}; reproduced by H0-only: {100 * h0_rep:.0f}%; H1-only significant where combined not: {h1_alone}",
          f"  -> 'predominantly H0' {'RETAINED' if pred_h0 else 'REVISED - see table'}",
          f"  'loops drive progression' rule (H1 significant in both versions in >= 3 of 4 GEO cohorts): {'MET' if loops else 'NOT MET'}"
          f" (H1-only Wasserstein: {sorted(set(h1_both)) or 'none'}; H1 landscape: {sorted(set(lands_both)) or 'none'})",
          "  GSE117606 H0/H1 split (new):"]
    for _, r in dd.sort_values(["version", "mode", "dim", "transition"]).iterrows():
        S.append(f"    {r.version:<10} {r['mode']:<10} {r.dim} {r.transition:<20} z {r.z:6.2f}  Holm {r.p_holm:.3f}{'  *' if r.p_holm < .05 else ''}")
    S += ["", "5. VALIDATION (section 6; full table incl. the two failures, for the supplement): results/addendum12/validation_summary.txt",
          "", "Figures: results/addendum12/figures/<cohort>_<version>_<mode>_{diagrams,barcodes,landscapes_betti}.png",
          "Outputs: landscape_tests.csv, barcode_summary_tests.csv, concordance_with_step5.csv, dimension_split_all_cohorts.csv"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(S) + "\n")
    print("\n".join(S))


def run_real():
    os.makedirs(os.path.join(OUT, "checkpoints"), exist_ok=True)
    t0 = time.time()
    units = [u for u in real_units()]
    du = dims_units()
    with Pool(WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(dims_unit, du), 1):
            print(f"[a12 dims] {i}/{len(du)} [{(time.time() - t0) / 60:.1f} min]", flush=True)
        for i, r in enumerate(pool.imap_unordered(real_unit, units), 1):
            print(f"[a12] {i}/{len(units)} {r['cohort']} {r['version']} {r['mode']} {r['transition']} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    for c in COH12:
        for v in VERS12:
            figures(c, v)
    print(f"[a12] figures done [{(time.time() - t0) / 60:.1f} min]", flush=True)
    summarise()


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "validate":
        validate()
    elif cmd == "run":
        run_real()
    elif cmd == "summary":
        summarise()
    else:
        raise SystemExit("Usage: addendum12_ph_outputs.py validate | run | summary")
