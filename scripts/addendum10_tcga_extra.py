"""
Addendum 10 (296eaed) section 4 - additional TCGA-COAD analyses, run one at a time.

  spread  : section 4.1 - PERMDISP per transition (bias-adjusted distances, as addendum 9 Amendment 1; 999
            permutations; Holm across the 3 transitions) with dispersion ratio D; and H0-only / H1-only versions
            of the step 5 test (as addendum 5 part 4), unadjusted and adjusted, shape-only and raw.

Usage: python scripts/addendum10_tcga_extra.py spread
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

SEED, NPERM, K, WORKERS = 42, 999, 20, 6
EMB = os.path.join(ROOT, "results", "step03", "embeddings")
OUT = os.path.join(ROOT, "results", "addendum10")
STAGES = ["Stage I", "Stage II", "Stage III", "Stage IV"]
TRANS = list(zip(STAGES[:-1], STAGES[1:]))
ROMAN = {"Stage I": "I", "Stage II": "II", "Stage III": "III", "Stage IV": "IV"}
VERSIONS = ("unadjusted", "adjusted")
EXPECTED = [45, 110, 80, 39]


def rng_for(name):
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(zlib.crc32(name.encode()),)))


def load(version):
    df = pd.read_csv(os.path.join(EMB, f"TCGA_COAD_main_{version}_g2000_pc10.csv"), index_col=0)
    if [int((df.stage == s).sum()) for s in STAGES] != EXPECTED:
        raise SystemExit("Stage counts changed")
    return df[[c for c in df.columns if c.startswith("PC")]].values, df.stage.values


# ---------------------------------------------------------------- 4.1 spread and H0 / H1
def dims_unit(u):
    v, mo, d, s, t = u
    name = f"TCGA41__{v[:5]}__{'shape' if mo == 'norm' else 'raw'}__H{d}__{ROMAN[s]}-{ROMAN[t]}"
    ck = os.path.join(OUT, "spread_H0H1", "checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    emb, st = load(v)
    m = int(0.8 * min(EXPECTED))
    i_s, i_t = np.where(st == s)[0], np.where(st == t)[0]
    Tobs, p, z, _ = T.permutation_test(emb, i_s, i_t, m, K, (d,), mo == "norm", NPERM, rng_for(name), align=True)
    rec = dict(unit=name, version=v, mode="shape-only" if mo == "norm" else "raw", dim=f"H{d}", transition=f"{s}->{t}",
               T_obs=float(Tobs), p=float(p), z=float(z))
    json.dump(rec, open(ck, "w"))
    return rec


def spread():
    from addendum9_heterogeneity import loo_dist, permdisp_dist   # bias-adjusted (addendum 9 Amendment 1)
    out = os.path.join(OUT, "spread_H0H1")
    os.makedirs(os.path.join(out, "checkpoints"), exist_ok=True)
    t0 = time.time()
    pd_rows = []
    for v in VERSIONS:
        emb, st = load(v)
        for s, t in TRANS:
            A, B = emb[st == s], emb[st == t]
            dA, dB = loo_dist(A, np.zeros(len(A))), loo_dist(B, np.zeros(len(B)))
            F, p, D = permdisp_dist(dA, dB, rng_for(f"TCGA41__permdisp__{v}__{ROMAN[s]}-{ROMAN[t]}"))
            pd_rows.append(dict(version=v, transition=f"{s}->{t}", n_from=len(A), n_to=len(B), F=F, p=p, D=D,
                                mean_dist_from=float(dA.mean()), mean_dist_to=float(dB.mean())))
    pdf = pd.DataFrame(pd_rows)
    pdf["p_holm"] = np.nan
    for _, g in pdf.groupby("version"):
        pdf.loc[g.index, "p_holm"] = T.holm(g.p.values)
    pdf.to_csv(os.path.join(out, "permdisp_by_stage.csv"), index=False)
    units = [(v, mo, d, s, t) for v in VERSIONS for mo in ("norm", "raw") for d in (0, 1) for s, t in TRANS]
    with Pool(WORKERS) as pool:
        recs = pool.map(dims_unit, units)
    dd = pd.DataFrame(recs)
    dd["p_holm"] = np.nan
    for _, g in dd.groupby(["version", "mode", "dim"]):
        dd.loc[g.index, "p_holm"] = T.holm(g.p.values)
    s5 = pd.read_csv(os.path.join(ROOT, "results", "step05", "TCGA", "permutation_results_all.csv"))
    s5 = s5[s5.analysis == "primary"].assign(dim="H0+H1")[["version", "mode", "dim", "transition", "z", "p", "p_holm"]]
    allr = pd.concat([s5, dd[["version", "mode", "dim", "transition", "z", "p", "p_holm"]]], ignore_index=True)
    allr.to_csv(os.path.join(out, "dims_results_all.csv"), index=False)
    L = ["ADDENDUM 10 section 4.1 - TCGA-COAD spread (PERMDISP) and H0-only / H1-only tests by stage", "=" * 72,
         "Embeddings: TCGA_COAD_main_<version>_g2000_pc10 (274 tumours; I/II/III/IV 45/110/80/39).",
         "PERMDISP: bias-adjusted distances sqrt(n/(n-1))*||x-c|| (addendum 9 Amendment 1), 999 permutations,",
         "  Holm across the 3 transitions per version; D = mean distance (later) / mean distance (earlier).",
         "H0/H1-only: step 5 test unchanged except dims (k = 20, m = 31, aligned, 999 permutations), Holm across the",
         "  3 transitions per version x mode x dimension (as addendum 5 part 4).",
         f"Seed {SEED} per unit. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "scipy")), "",
         "PERMDISP (spread):"]
    for _, r in pdf.iterrows():
        L.append(f"  {r.version:<11} {r.transition:<22} D = {r.D:.3f}  (mean distance {r.mean_dist_from:.1f} -> {r.mean_dist_to:.1f})"
                 f"  p {r.p:.3f}  Holm {r.p_holm:.3f}{'  *' if r.p_holm < .05 else ''}")
    L += ["", "z (Holm p) by homology dimension:  H0+H1 (step 5) | H0 only | H1 only"]
    for v in VERSIONS:
        for mo in ("shape-only", "raw"):
            for s, t in TRANS:
                tr = f"{s}->{t}"
                cell = []
                for d in ("H0+H1", "H0", "H1"):
                    r = allr[(allr.version == v) & (allr["mode"] == mo) & (allr.dim == d) & (allr.transition == tr)].iloc[0]
                    cell.append(f"{r.z:6.2f} ({r.p_holm:.3f}){'*' if r.p_holm < .05 else ' '}")
                L.append(f"  {v:<11} {mo:<10} {tr:<22} " + " | ".join(cell))
    nsig = int((dd.p_holm < .05).sum())
    L += ["", f"Significant (Holm p < 0.05): PERMDISP {int((pdf.p_holm < .05).sum())}/6; H0-only / H1-only {nsig}/24.",
          f"Run time {(time.time() - t0) / 60:.1f} min.",
          "Outputs: results/addendum10/spread_H0H1/permdisp_by_stage.csv, dims_results_all.csv, checkpoints/"]
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ---------------------------------------------------------------- 4.2 CMS by stage (descriptive)
def cms():
    """CMScaller (r-cms; as addendum 9 / D16a): unadjusted log2 expression of the 274 tumours, gene-wise centred,
    RNAseq = FALSE (Xena data already log2-normalised), nPerm 1000, seed 42. Proportions by stage (descriptive)."""
    from addendum9_heterogeneity import run_r
    out = os.path.join(OUT, "cms")
    cache = os.path.join(out, "cache")
    os.makedirs(cache, exist_ok=True)
    lab = pd.read_csv(os.path.join(ROOT, "data", "processed", "TCGA_COAD_main_samples_qc.csv"), index_col=0)
    if [int((lab.stage == s).sum()) for s in STAGES] != EXPECTED:
        raise SystemExit("Stage counts changed")
    X = pd.read_csv(os.path.join(ROOT, "data", "processed", "TCGA_COAD_main_expr_genes.csv"), index_col=0)[lab.index]
    X = X.dropna(axis=0, how="any")
    X = X[X.var(axis=1) > 0]                       # constant genes carry no information for template correlation
    inp, outp = os.path.join(cache, "cms_in_TCGA_COAD.csv"), os.path.join(cache, "cms_out_TCGA_COAD.csv")
    X.to_csv(inp)
    rver = run_r("classify", inp, outp)
    res = pd.read_csv(outp).set_index("sample").loc[lab.index]
    cov = pd.read_csv(outp.replace(".csv", "_coverage.csv"))
    tab = lab[["stage", "msi_class"]].join(res)
    tab["confident"] = tab.prediction.notna()
    tab["coverage"] = float(cov.overall.iloc[0])
    tab.to_csv(os.path.join(out, "cms_labels_TCGA.csv"), index_label="sample")
    cms4 = ["CMS1", "CMS2", "CMS3", "CMS4"]
    prim = pd.crosstab(tab.stage, tab.nearest).reindex(index=STAGES, columns=cms4, fill_value=0)
    conf = pd.crosstab(tab.stage, tab.prediction.fillna("unclassified")).reindex(index=STAGES, columns=cms4 + ["unclassified"],
                                                                                 fill_value=0)
    msi = pd.crosstab(tab.nearest, tab.msi_class.str.startswith("MSI-H").map({True: "MSI-H", False: "MSS/other"}))
    prim.to_csv(os.path.join(out, "cms_by_stage_primary.csv"))
    conf.to_csv(os.path.join(out, "cms_by_stage_confident.csv"))
    pct = lambda t: (100 * t.div(t.sum(axis=1), axis=0)).round(0).astype(int)
    covall = float(cov.overall.iloc[0])
    L = ["ADDENDUM 10 section 4.2 - CMS subtypes by stage, TCGA-COAD (descriptive)", "=" * 66,
         f"CMScaller via r-cms: {rver}; unadjusted log2(RSEM+1) of 274 tumours, gene-wise centred (D16a), RNAseq = FALSE,",
         "nPerm 1000, seed 42, FDR 0.05. Primary = nearest template; confident = FDR < 0.05.",
         f"Genes passed to CMScaller: {X.shape[0]} (missing values and zero-variance genes removed).",
         f"Template coverage: {100 * covall:.1f}%" + ("  ** LOW COVERAGE (< 80%) **" if covall < .8 else "")
         + " (per template: " + ", ".join(f"{r.template} {100 * r.coverage:.0f}%" for _, r in cov.iterrows()) + ")",
         "", "Primary labels, counts:", "  " + prim.to_string().replace("\n", "\n  "),
         "Primary labels, % of stage:", "  " + pct(prim).to_string().replace("\n", "\n  "),
         "", "Confident labels (FDR < 0.05), % of stage:", "  " + pct(conf).to_string().replace("\n", "\n  "),
         "", "All tumours: " + ", ".join(f"{k} {100 * v:.0f}%" for k, v in tab.nearest.value_counts(normalize=True).reindex(cms4).items())
         + f"; unclassified (confident) {100 * (~tab.confident).mean():.0f}%",
         "Guinney 2015: CMS1 14%, CMS2 37%, CMS3 13%, CMS4 23%, unclassified 13%",
         "", "CMS (primary) by MSI status (counts):", "  " + msi.to_string().replace("\n", "\n  "),
         "", "Caution: CMS4 calls on unadjusted expression partly reflect stromal content of the sample, so the rise",
         "of CMS4 with stage may partly reflect more stroma in advanced tumours rather than tumour-cell biology alone.",
         "", "Outputs: results/addendum10/cms/cms_labels_TCGA.csv, cms_by_stage_primary.csv, cms_by_stage_confident.csv"]
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ---------------------------------------------------------------- 4.3 normal -> Stage I (D22)
def load_wn(version):
    df = pd.read_csv(os.path.join(EMB, f"TCGA_COAD_withnormal_{version}_g2000_pc10.csv"), index_col=0)
    cnt = [int((df.stage == s).sum()) for s in ["Normal"] + STAGES]
    if cnt != [41] + EXPECTED:
        raise SystemExit(f"Stage counts changed: {cnt}")
    return df


def normal_unit(u):
    kind, v, mo = u
    name = f"TCGA43__{kind}__{v[:5]}__{'shape' if mo == 'norm' else 'raw'}"
    ck = os.path.join(OUT, "normal_stageI", "checkpoints", name + ".json")
    if os.path.exists(ck):
        return json.load(open(ck))
    rng = rng_for(name)
    df = load_wn(v)
    emb = df[[c for c in df.columns if c.startswith("PC")]].values
    st = df.stage.values
    patient = np.array([s[:12] for s in df.index], dtype=object)
    if kind == "main":
        i_s, i_t = np.where(st == "Normal")[0], np.where(st == "Stage I")[0]
        shared = len(set(patient[i_s]) & set(patient[i_t]))
        i_s, i_t, nres = T.drop_shared_patients(i_s, i_t, patient, rng)
        m = min(int(0.8 * min(EXPECTED)), int(0.8 * min(len(i_s), len(i_t))))
        Tobs, p, z, _ = T.permutation_test(emb, i_s, i_t, m, K, (0, 1), mo == "norm", NPERM, rng, align=True)
        rec = dict(n_from=int(len(i_s)), n_to=int(len(i_t)), shared_patients=shared, resolved=int(nres))
    else:
        import step05_secondary_GSE117606 as S2
        S2.K, S2.NPERM, S2.DIMS = K, NPERM, (0, 1)
        nidx = np.where(st == "Normal")[0]
        tum = {patient[i]: i for i in np.where(st != "Normal")[0]}
        pairs = [(i, tum[patient[i]]) for i in nidx if patient[i] in tum]
        ps, pt = np.array([a for a, _ in pairs]), np.array([b for _, b in pairs])
        m = int(0.8 * len(pairs))
        Tobs, p, z = S2.paired_test(emb, ps, pt, m, mo == "norm", rng)
        rec = dict(n_pairs=len(pairs), tumour_stages=pd.Series(st[pt]).value_counts().to_dict())
    rec.update(unit=name, analysis=kind, version=v, mode="shape-only" if mo == "norm" else "raw", subsample_size=int(m),
               T_obs=float(Tobs), p=float(p), z=float(z))
    json.dump(rec, open(ck, "w"))
    return rec


def normal():
    out = os.path.join(OUT, "normal_stageI")
    os.makedirs(os.path.join(out, "checkpoints"), exist_ok=True)
    t0 = time.time()
    units = [(k, v, mo) for k in ("main", "paired") for v in VERSIONS for mo in ("norm", "raw")]
    with Pool(WORKERS) as pool:
        recs = pool.map(normal_unit, units)
    d = pd.DataFrame(recs)
    d.to_csv(os.path.join(out, "normal_stageI_results.csv"), index=False)
    m0 = d[d.analysis == "main"].iloc[0]
    p0 = d[d.analysis == "paired"].iloc[0]
    L = ["ADDENDUM 10 section 4.3 - TCGA-COAD solid-tissue normal -> Stage I (secondary, descriptive; D22)", "=" * 72,
         "Embeddings: TCGA_COAD_withnormal_<version>_g2000_pc10 (41 normals + 274 tumours; adjusted = adjustment fitted",
         "on all 315). Step 5 test unchanged (k = 20, aligned, 999 permutations); single comparison, no Holm.",
         "Caveat: TCGA normals are tumour-adjacent tissue (possible field effect).",
         f"Seed {SEED} per unit. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy")), "",
         f"(a) MAIN: Normal vs Stage I, one sample per patient ({int(m0.shared_patients)} shared patients assigned at random to one",
         f"    group): n = {int(m0.n_from)} normals / {int(m0.n_to)} Stage I, subsample m = {int(m0.subsample_size)}"]
    for _, r in d[d.analysis == "main"].iterrows():
        L.append(f"    {r.version:<11} {r['mode']:<10} T = {r.T_obs:9.4f}  z = {r.z:6.2f}  p = {r.p:.4f}{'  *' if r.p < .05 else ''}")
    L += ["", f"(b) PAIRED SENSITIVITY: normal vs the same patient's own tumour (any stage; {int(p0.n_pairs)} pairs; tumour stages "
          + ", ".join(f"{k} {v}" for k, v in sorted(p0.tumour_stages.items())) + f"), subsample m = {int(p0.subsample_size)} patients,",
          "    labels swapped within patient pairs (GSE117606 paired test). Related but different question (normal -> carcinoma",
          "    within patient); does not replace (a)."]
    for _, r in d[d.analysis == "paired"].iterrows():
        L.append(f"    {r.version:<11} {r['mode']:<10} T = {r.T_obs:9.4f}  z = {r.z:6.2f}  p = {r.p:.4f}{'  *' if r.p < .05 else ''}")
    L += ["", f"Run time {(time.time() - t0) / 60:.1f} min.", "Outputs: results/addendum10/normal_stageI/normal_stageI_results.csv, checkpoints/"]
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


# ---------------------------------------------------------------- 4.4 pilot method vs label-permutation test
def h1_total(points, normalize):
    d = T.persistence(points, normalize)[1]
    return float((d[:, 1] - d[:, 0]).sum()) if len(d) else 0.0


def subsample_stats(emb, idx, k, n, normalize, rng):
    return np.array([h1_total(emb[rng.choice(idx, k, replace=False)], normalize) for _ in range(n)])


def perm_unit(u):
    normalize, chunk = u
    rng = rng_for(f"TCGA44__perm__{'shape' if normalize else 'raw'}__c{chunk}")
    emb, st = load("unadjusted")
    pool = np.r_[np.where(st == "Stage I")[0], np.where(st == "Stage IV")[0]]
    n1 = int((st == "Stage I").sum())
    out = []
    for _ in range(111):                      # 9 chunks x 111 = 999 permutations
        perm = rng.permutation(pool)
        a, b = subsample_stats(emb, perm[:n1], 25, 20, normalize, rng), subsample_stats(emb, perm[n1:], 25, 20, normalize, rng)
        out.append(abs(a.mean() - b.mean()))
    return out


def pilot():
    from scipy.stats import mannwhitneyu
    out = os.path.join(OUT, "pilot_vs_patient")
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    emb, st = load("unadjusted")
    i1, i4 = np.where(st == "Stage I")[0], np.where(st == "Stage IV")[0]
    rows = []
    # (a) pilot method: bootstrap subsamples (k = 25, without replacement within a subsample), Mann-Whitney across subsamples
    rng = rng_for("TCGA44__pilot")
    big1, big4 = subsample_stats(emb, i1, 25, 1000, False, rng), subsample_stats(emb, i4, 25, 1000, False, rng)
    for nb in (25, 50, 100, 250, 500, 1000):
        u, p = mannwhitneyu(big1[:nb], big4[:nb], alternative="two-sided")
        rows.append(dict(method="(a) pilot: Mann-Whitney across bootstrap subsamples", n_subsamples=nb,
                         mean_I=big1[:nb].mean(), mean_IV=big4[:nb].mean(), p=float(p)))
    # (b), (c) label permutation of stage labels; statistic |mean_I - mean_IV| over 20 subsamples per stage
    with Pool(WORKERS) as pool:
        res = {}
        for normalize in (False, True):
            r0 = rng_for(f"TCGA44__obs__{'shape' if normalize else 'raw'}")
            a, b = subsample_stats(emb, i1, 25, 20, normalize, r0), subsample_stats(emb, i4, 25, 20, normalize, r0)
            obs = abs(a.mean() - b.mean())
            null = np.concatenate(pool.map(perm_unit, [(normalize, c) for c in range(9)]))
            p = (1 + np.sum(null >= obs)) / (1 + len(null))
            res[normalize] = (a.mean(), b.mean(), obs, p, (obs - null.mean()) / null.std(ddof=1))
            rows.append(dict(method=f"({'c' if normalize else 'b'}) label permutation, {'shape-only' if normalize else 'raw'}",
                             n_subsamples=20, mean_I=a.mean(), mean_IV=b.mean(), p=float(p)))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(out, "pilot_vs_patient_results.csv"), index=False)
    L = ["ADDENDUM 10 section 4.4 - pilot method vs label-permutation test (TCGA-COAD Stage I 45 vs Stage IV 39)", "=" * 72,
         "Embedding: TCGA_COAD_main_unadjusted_g2000_pc10 (GUDHI Vietoris-Rips). Statistic: total H1 persistence per subsample",
         "of k = 25 patients per stage (the pilot's statistic). The pilot's original numbers are not re-used.",
         f"Seed {SEED} per unit. " + ", ".join(f"{p} {metadata.version(p)}" for p in ("gudhi", "numpy", "scipy")), "",
         "(a) Pilot method - subsamples drawn from the same 45 / 39 patients, Mann-Whitney test across subsamples (raw):"]
    for _, r in d[d.method.str.startswith("(a)")].iterrows():
        L.append(f"    {int(r.n_subsamples):>5} subsamples per stage: mean H1 total I {r.mean_I:.2f} vs IV {r.mean_IV:.2f}  ->  p = {r.p:.2g}")
    L += ["", "(b)/(c) Label-permutation test - stage labels permuted across the 84 patients, 999 permutations, 20 subsamples per stage:"]
    for normalize in (False, True):
        mI, mIV, obs, p, z = res[normalize]
        L.append(f"    {'shape-only' if normalize else 'raw       '}: mean H1 total I {mI:.3f} vs IV {mIV:.3f}; |difference| {obs:.3f}; z {z:.2f}; p = {p:.3f}")
    L += ["", "Reading: in (a) the p-value depends on how many subsamples are drawn from the same patients (pseudo-replication);",
          "the label-permutation test (b)/(c) treats patients as the unit and gives the valid p-value.",
          f"Run time {(time.time() - t0) / 60:.1f} min.", "Outputs: results/addendum10/pilot_vs_patient/pilot_vs_patient_results.csv"]
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "pilot":
        pilot()
        sys.exit()
    if cmd == "spread":
        spread()
    elif cmd == "cms":
        cms()
    elif cmd == "normal":
        normal()
    else:
        raise SystemExit("Usage: addendum10_tcga_extra.py spread | cms | normal")
