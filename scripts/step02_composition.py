"""
Step 2 - Tissue-composition scoring for GSE41258.

  0. Tissue-identity screen of every sample (fixed textbook panels: colon
     epithelium, hepatocyte, lung) against normal colon / liver / lung.
     Samples whose expression contradicts their label are listed; those in
     EXCLUDE (decided with the PI) are removed, with the reason recorded.
  A. Epithelial / stromal / immune scores per sample, decoupler ULM on
     gene-centred log2 data.
       Epithelial: pre-specified pan-epithelial structural genes, kept only if
         their stage medians differ by <= 1 log2 (flat across progression), so the
         score measures epithelial CONTENT, not differentiation state.
       Stromal / immune: PanglaoDB human canonical markers (decoupler resource);
         genes shared between compartments removed.
  B. Host-organ content: organ signature = genes high in the excluded normal
     liver (or lung) samples and low in normal colon and primary carcinoma
     (metastases not used). Organ fraction per sample = two-component linear
     mixture (carcinoma median + organ median, linear scale), per-gene fraction,
     median across signature genes.
  C. Composition-adjusted data: each gene regressed (OLS) on all composition
     covariates; adjusted = residual + gene mean.

Analysis sets
  main           liver metastases only, EXCLUDE removed, one array per patient per stage (primary)
  withlung       + lung metastases, same rules                                         (sensitivity)
  main_allarrays main without de-duplication (all repeat arrays kept)                  (sensitivity)
  main_inclqc    main with the EXCLUDE samples put back                                 (sensitivity)

Outputs: results/step02/*, data/processed/GSE41258_exclusions.csv,
         data/processed/GSE41258_<set>_samples_qc.csv,
         data/processed/GSE41258_<set>_expr_genes_adjusted.csv
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402
import decoupler as dc  # noqa: E402

SEED = 42
np.random.seed(SEED)
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "step02")
os.makedirs(OUT, exist_ok=True)
MATRIX = os.path.join(ROOT, "data", "raw", "GSE41258_series_matrix.txt.gz")

EXCLUDE = {
    "GSM1012647": "Labelled 'Primary Tumor' but expression is normal liver: no colon identity "
                  "(CDX2 0.0, KRT20 4.7, EPCAM 5.6 log2) and full hepatocyte program "
                  "(ALB 13.7, APOA1 15.3, FGB 15.3 log2; liver fraction 100%); no patient id. "
                  "Excluded after PI review, 2026-10-03; re-included in set main_inclqc.",
    "GSM1012526": "Labelled 'Primary Tumor' but no colorectal identity (CDX2 2.2, KRT20 4.1, "
                  "SATB2 4.5 log2, lowest in cohort) and a Mullerian/gynaecological program "
                  "(PAX8 9.2, SOX17 9.7, MUC16 10.4 vs carcinoma medians 7.7, 4.0, 2.8); "
                  "more consistent with an ovarian/Mullerian carcinoma than CRC. "
                  "Excluded after PI review, 2026-10-03; re-included in set main_inclqc.",
}
EXCLUSION_STATEMENT = ("All exclusions were based on tissue identity only (expression of fixed "
                       "colon, hepatocyte, lung and lineage marker panels) and were decided with "
                       "the PI on 2026-10-03, before any topology analysis (steps 3-6) was run.")
# flagged by the identity screen, reviewed with the PI and KEPT
KEPT_NOTES = {
    "GSM1012429": "Carcinoma; low CDX2/KRT20/SATB2, low MLH1 (5.8 vs 7.2), immune-rich: "
                  "consistent with MSI-high / medullary CRC. Kept.",
    "GSM1012488": "Carcinoma; low KRT20, highest immune score: consistent with MSI / medullary CRC. Kept.",
    "GSM1012326": "Carcinoma; low CDX2 but colonic KRT20/SATB2, squamous genes (KRT5, KRT14) on: "
                  "consistent with squamous differentiation within CRC. Kept.",
    "GSM1012396": "Carcinoma; moderately low colon markers, consistent with CDX2-low CRC. Kept.",
    "GSM1012525": "Carcinoma; moderately low colon markers, consistent with CDX2-low CRC. Kept.",
    "GSM1012470": "Liver metastasis; CDX2 near absent (0.4), KRT20 reduced: possible CDX2-negative "
                  "CRC metastasis. Kept, noted.",
    "GSM1012347": "Normal; little epithelium (EPCAM 7.8), muscle-rich (DES 12.7): deep sampling of "
                  "colon wall. Kept; handled by composition adjustment.",
    "GSM1012640": "Normal; colon identity normal, ~2% liver signal; same patient (C0415) as a liver "
                  "metastasis, likely handling contamination. Kept.",
}
STEP1_COUNTS = {"main": {"Normal": 54, "Adenoma": 51, "Carcinoma": 186, "Metastasis": 47},
                "withlung": {"Normal": 54, "Adenoma": 51, "Carcinoma": 186, "Metastasis": 67}}

# fixed identity panels (textbook markers; not derived from these data)
PANELS = {"colon": ["CDX1", "CDX2", "CDH17", "VIL1", "KRT20", "CEACAM5", "LGALS4"],
          "hepatocyte": ["ALB", "APOA1", "APOA2", "APOH", "FGA", "FGB", "FGG", "AHSG", "TF",
                         "HP", "SERPINA1", "TTR"],
          "lung": ["SFTPC", "SFTPB", "SFTPD", "SCGB1A1", "AGER"]}

PAN_EPI_CANDIDATES = ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4", "CLDN7",
                      "ELF3", "ESRP1", "SPINT1", "SPINT2", "RAB25", "GRHL2", "OVOL1", "PRSS8",
                      "ST14", "CDS1", "TJP3"]
COMPARISON_GENES = ["KRT20", "TACSTD2"]   # shown for contrast: differentiation / tumour genes
FLAT_MAX = 1.0                            # max - min of stage medians, log2
FLAT_STRICT = 0.75                        # sensitivity: stricter epithelial set

COMPARTMENTS = {
    "Stromal": ["Fibroblasts", "Myofibroblasts", "Smooth muscle cells", "Endothelial cells",
                "Pericytes"],
    "Immune": ["T cells", "B cells", "Plasma cells", "NK cells", "Macrophages", "Monocytes",
               "Dendritic cells", "Mast cells", "Neutrophils"],
}
SANITY = {"Stromal": ["COL1A2", "COL3A1", "DCN", "VWF"],
          "Immune": ["PTPRC", "CD3E", "CD74", "LCP1"]}
ORGAN_LFC, ORGAN_MIN = 3.0, 8.0   # organ signature: >= 8-fold above colon/carcinoma, organ mean >= 8

# ---------------------------------------------------------------- load step-1 data
step1 = {}
for tag in ("main", "withlung"):
    g = pd.read_csv(os.path.join(PROC, f"GSE41258_{tag}_expr_genes.csv"), index_col=0)
    lab = pd.read_csv(os.path.join(PROC, f"GSE41258_{tag}_samples.csv"), index_col=0)
    assert list(g.columns) == list(lab.index)
    if lab["stage"].value_counts().to_dict() != STEP1_COUNTS[tag]:
        sys.exit(f"Step-1 stage counts changed for {tag}")
    step1[tag] = (g, lab)

expr_all, meta_all = T.read_series_matrix(MATRIX)
expr_all = np.log2(expr_all.clip(lower=1))          # step 1 applied log2
all_lab = pd.read_csv(os.path.join(ROOT, "results", "step01", "all_samples_labels.csv"),
                      index_col=0)


def organ_samples(tag, pattern):
    pick = pd.read_csv(os.path.join(ROOT, "results", "step01", f"probe_to_gene_{tag}.csv"))
    ids = all_lab.index[all_lab["raw_label"].str.lower().str.contains(pattern)]
    x = expr_all.loc[pick["probe"], ids]
    x.index = pick["gene"].values
    return x.sort_index()


liver_ref = organ_samples("main", "normal liver")
lung_ref = organ_samples("withlung", "normal lung")

# ---------------------------------------------------------------- 0. identity screen
g_all, lab_all = step1["withlung"]          # superset: 338 main + 20 lung metastases
grp = lab_all["stage"].where(~lab_all["raw_label"].str.lower().str.contains("lung"), "LungMet")
pm = pd.DataFrame({p: g_all.loc[[x for x in gl if x in g_all.index]].mean() for p, gl in PANELS.items()})
ref_l = pd.Series({p: liver_ref.loc[[x for x in gl if x in liver_ref.index]].mean().median()
                   for p, gl in PANELS.items()})
ref_u = pd.Series({p: lung_ref.loc[[x for x in gl if x in lung_ref.index]].mean().median()
                   for p, gl in PANELS.items()})
ref_c = pm.loc[grp == "Normal"].median()
mid = {"colon_vs_liver": (ref_c["colon"] + ref_l["colon"]) / 2,
       "colon_vs_lung": (ref_c["colon"] + ref_u["colon"]) / 2,
       "hepatocyte": (ref_c["hepatocyte"] + ref_l["hepatocyte"]) / 2,
       "lung": (ref_c["lung"] + ref_u["lung"]) / 2}
idq = pm.copy()
idq.insert(0, "group", grp)
idq.insert(1, "raw_label", lab_all["raw_label"])
idq["low_colon"] = idq["colon"] < min(mid["colon_vs_liver"], mid["colon_vs_lung"])
idq["high_liver"] = idq["hepatocyte"] > mid["hepatocyte"]
idq["high_lung"] = idq["lung"] > mid["lung"]


def call(r):
    if r.low_colon and r.high_liver:
        return "MISMATCH: liver-like"
    if r.low_colon and r.high_lung:
        return "MISMATCH: lung-like"
    if r.low_colon:
        return "check: low colon identity only"
    if r.high_liver or r.high_lung:
        expected = (r.group == "Metastasis" and r.high_liver) or (r.group == "LungMet" and r.high_lung)
        return "host organ present (expected for site)" if expected else "check: other-organ signal"
    return "ok"


idq["call"] = idq.apply(call, axis=1)
idq["excluded"] = idq.index.isin(list(EXCLUDE))
idq["review_note"] = [EXCLUDE.get(s, KEPT_NOTES.get(s, "")) for s in idq.index]
idq.round(3).to_csv(os.path.join(OUT, "tissue_identity_qc.csv"), index_label="sample")
flagged = idq[idq["call"] != "ok"].sort_values(["call", "colon"])
print("[identity] calls:\n" + idq["call"].value_counts().to_string())

# ---------------------------------------------------------------- repeat arrays
# Official GEO fields: patient id, authors' 'included in analysis' flag, array code.
def ch_field(name):
    out = pd.Series("", index=list(expr_all.columns))
    for row in meta_all.get("!Sample_characteristics_ch1", []):
        for i, v in enumerate(row):
            if v.lower().startswith(name):
                out.iloc[i] = v.split(":", 1)[1].strip()
    return out


geo = pd.DataFrame({"patient_id": ch_field("patient id"),
                    "authors_included": ch_field("included in analysis"),
                    "array_code": pd.Series(meta_all["!Sample_title"][0],
                                            index=list(expr_all.columns)).str.split().str[-1]})
geo["repeat_suffix"] = geo["array_code"].str.contains("_")


def dedup(lab):
    """One array per patient per stage: prefer authors' 'Yes'; then an original array
    (no '_ez'/'_2'/... suffix); then the lowest array code. Returns kept ids + records."""
    d = lab[["stage", "raw_label"]].join(geo)
    d["_rank"] = list(zip(d["authors_included"] != "Yes", d["repeat_suffix"], d["array_code"]))
    keep, rec = [], []
    for (pid, stage), grp in d.groupby(["patient_id", "stage"], sort=False):
        grp = grp.sort_values("_rank")
        k0 = grp.index[0]
        keep.append(k0)
        for s in grp.index[1:]:
            base_s, base_k = (grp.loc[x, "array_code"].split("_")[0] for x in (s, k0))
            # array code = patient + tissue letter + 'R' + piece number (+ _suffix for re-hybridisation)
            if base_s == base_k:
                kind = "re-hybridisation of the same RNA"
            elif base_s[:-3] == base_k[:-3] and base_s[-3] == base_k[-3]:
                kind = "another piece of the same tissue (different R number)"
            else:
                kind = "possibly a separate lesion (different tissue letter)"
            st = "Metastasis (lung)" if "lung" in grp.loc[s, "raw_label"].lower() else stage
            rec.append({"sample": s, "step1_stage": st, "category": "same patient and stage",
                        "kind": kind,
                        "reason": f"Additional array of patient {pid}, {st}: {kind} (code {grp.loc[s, 'array_code']}, "
                                  f"authors' flag {grp.loc[s, 'authors_included'] or 'NA'}); kept "
                                  f"{grp.index[0]} (code {grp.loc[grp.index[0], 'array_code']}, flag "
                                  f"{grp.loc[grp.index[0], 'authors_included'] or 'NA'}). One array per "
                                  "patient per stage, decided with the PI 2026-10-03 before topology; "
                                  "all arrays kept in set main_allarrays."})
    return [s for s in lab.index if s in set(keep)], rec


DEDUP_STATEMENT = ("Where one patient had several arrays of the same stage (official "
                   "GEO 'patient id' field), one was kept: the array flagged 'included in analysis: "
                   "Yes' by the original authors, otherwise the first original array (not an '_ez' "
                   "or other re-hybridisation). Applies to re-hybridisations, other pieces of the same tissue and, in one case, a possible second lesion. Decided with the PI on 2026-10-03, before any "
                   "topology analysis. Set main_allarrays keeps all arrays (sensitivity).")

# ---------------------------------------------------------------- analysis sets
SET_DEFS = [  # name, step-1 table, drop identity exclusions, de-duplicate
    ("main", "main", True, True),
    ("withlung", "withlung", True, True),
    ("main_allarrays", "main", True, False),
    ("main_inclqc", "main", False, True),
]
sets, dedup_rec = {}, {}
for name, tag, drop, dd in SET_DEFS:
    g, lab = step1[tag]
    keep = [s for s in lab.index if not (drop and s in EXCLUDE)]
    lab = lab.loc[keep].copy()
    if dd:
        keep, rec = dedup(lab)
        dedup_rec[name] = rec
        lab = lab.loc[keep]
    lab = lab.rename(columns={"patient": "patient_pilot_regex"}).join(geo)
    sets[name] = (g[lab.index], lab)
    lab.to_csv(os.path.join(PROC, f"GSE41258_{name}_samples_qc.csv"))

excl_rows = [{"sample": s, "step1_stage": lab_all.loc[s, "stage"], "category": "tissue identity",
              "reason": r} for s, r in EXCLUDE.items()]
seen = set()
for name in ("main", "withlung"):
    for r in dedup_rec[name]:
        if r["sample"] not in seen:
            excl_rows.append(r); seen.add(r["sample"])
excl_df = pd.DataFrame(excl_rows)[["sample", "step1_stage", "category", "kind", "reason"]]
excl_df["kind"] = excl_df["kind"].fillna("")
excl_df.to_csv(os.path.join(PROC, "GSE41258_exclusions.csv"), index=False)
print("[sets] " + "; ".join(f"{n}: " + ", ".join(f"{s} {v}" for s, v in
                                                sets[n][1]["stage"].value_counts().reindex(T.STAGES).items())
                            for n in sets))

# ---------------------------------------------------------------- A. marker sets
g_main, lab_main = sets["main"]
flat = g_main.loc[[x for x in PAN_EPI_CANDIDATES + COMPARISON_GENES if x in g_main.index]] \
    .T.groupby(lab_main["stage"]).median().T[T.STAGES]
flat["range"] = flat.max(axis=1) - flat.min(axis=1)
flat["role"] = ["comparison (not used)" if x in COMPARISON_GENES else "candidate" for x in flat.index]
flat["used"] = (flat["role"] == "candidate") & (flat["range"] <= FLAT_MAX)
flat.round(2).to_csv(os.path.join(OUT, "pan_epithelial_marker_flatness.csv"), index_label="gene")
epi_genes = list(flat.index[flat["used"]])
epi_strict = list(flat.index[flat["used"] & (flat["range"] <= FLAT_STRICT)])
missing_epi = [x for x in PAN_EPI_CANDIDATES if x not in g_main.index]

pg = dc.op.resource("PanglaoDB", organism="human")
pg = pg[(pg["human"].astype(bool)) & (pg["canonical_marker"].astype(bool))]
rows = [("Epithelial", x) for x in epi_genes]
for comp, cts in COMPARTMENTS.items():
    rows += [(comp, x) for x in set(pg.loc[pg["cell_type"].isin(cts), "genesymbol"])
             if x not in epi_genes]
net = pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()
shared = net["target"].value_counts()
net = net[net["target"].isin(shared.index[shared == 1])].copy()
net["weight"] = 1.0
net = net[net["target"].isin(g_main.index)]
COMP_NAMES = ["Epithelial", "Stromal", "Immune"]
set_sizes = net["source"].value_counts().reindex(COMP_NAMES).astype(int)
net.sort_values(["source", "target"]).to_csv(os.path.join(OUT, "marker_sets_used.csv"), index=False)
print("[markers] genes per compartment:\n" + set_sizes.to_string())


def ulm_scores(g):
    centred = g.sub(g.mean(axis=1), axis=0)
    res = dc.mt.ulm(centred.T, net[["source", "target", "weight"]], tmin=5)
    es = res[0] if isinstance(res, tuple) else res
    return es[COMP_NAMES].copy()


def epi_strict_score(g):
    centred = g.sub(g.mean(axis=1), axis=0)
    n2 = pd.DataFrame({"source": "Epithelial_strict", "target": epi_strict, "weight": 1.0})
    res = dc.mt.ulm(centred.T, n2, tmin=5)
    es = res[0] if isinstance(res, tuple) else res
    return es["Epithelial_strict"]


# ---------------------------------------------------------------- B. organ fraction
def organ_signature(g, lab, ref):
    m_ref = ref.mean(axis=1)
    m_norm = g[lab.index[lab["stage"] == "Normal"]].mean(axis=1)
    m_carc = g[lab.index[lab["stage"] == "Carcinoma"]].mean(axis=1)
    diff = m_ref - np.maximum(m_norm, m_carc)
    sig = diff[(diff >= ORGAN_LFC) & (m_ref >= ORGAN_MIN)].sort_values(ascending=False)
    return pd.DataFrame({"organ_mean": m_ref[sig.index], "normal_colon_mean": m_norm[sig.index],
                         "carcinoma_mean": m_carc[sig.index], "log2_diff": sig})


def organ_fraction(g, lab, ref, sig_genes):
    """x = (1-p)*C + p*L on linear scale; p per signature gene, median across genes.
    C = carcinoma median, L = organ median (robust to single contaminated samples)."""
    C = 2 ** g.loc[sig_genes, lab.index[lab["stage"] == "Carcinoma"]].median(axis=1)
    L = 2 ** ref.loc[sig_genes].median(axis=1)
    X = 2 ** g.loc[sig_genes]
    p = (X.sub(C, axis=0)).div(L - C, axis=0).median(axis=0)
    return pd.Series(np.clip(p.values, 0, 1), index=g.columns)


# ---------------------------------------------------------------- C. adjustment
def adjust(g, cov):
    Xc = np.c_[np.ones(len(cov)), (cov - cov.mean()).values]
    Y = g.T.values
    beta, *_ = np.linalg.lstsq(Xc, Y, rcond=None)
    resid = Y - Xc @ beta
    adj = pd.DataFrame(resid.T + g.mean(axis=1).values[:, None], index=g.index, columns=g.columns)
    r2 = 1 - resid.var(axis=0) / Y.var(axis=0)
    return adj, pd.Series(r2, index=g.index)


results = {}
for name, (g, lab) in sets.items():
    sc = ulm_scores(g)
    liv_tab = organ_signature(g, lab, liver_ref)
    liv_tab.to_csv(os.path.join(OUT, f"liver_signature_{name}.csv"), index_label="gene")
    sc["Liver_fraction"] = organ_fraction(g, lab, liver_ref, liv_tab.index)
    lung_tab = None
    if name == "withlung":
        lung_tab = organ_signature(g, lab, lung_ref)
        lung_tab.to_csv(os.path.join(OUT, "lung_signature_withlung.csv"), index_label="gene")
        sc["Lung_fraction"] = organ_fraction(g, lab, lung_ref, lung_tab.index)
    sc.insert(0, "stage", lab["stage"])
    sc.insert(1, "raw_label", lab["raw_label"])
    adj, r2 = adjust(g, sc.drop(columns=["stage", "raw_label"]))
    strict = epi_strict_score(g)          # sensitivity only, not used for adjustment
    sc_out = sc.copy()
    sc_out["Epithelial_strict15_sensitivity"] = strict
    sc_out.to_csv(os.path.join(OUT, f"composition_scores_{name}.csv"), index_label="sample")
    epi_corr = (sc["Epithelial"].corr(strict, method="pearson"),
                sc["Epithelial"].corr(strict, method="spearman"))
    adj.index.name = "gene"
    adj.to_csv(os.path.join(PROC, f"GSE41258_{name}_expr_genes_adjusted.csv"), float_format="%.5f")
    r2.rename("R2_composition").to_csv(
        os.path.join(OUT, f"variance_explained_by_composition_{name}.csv"), index_label="gene")
    results[name] = (sc, liv_tab, lung_tab, r2, g, epi_corr)

# ---------------------------------------------------------------- figures
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# pan-epithelial flatness
fig, ax = plt.subplots(figsize=(8, 6))
label_y = flat["Metastasis"].astype(float).sort_values()
placed, last = {}, -np.inf
for gname, y in label_y.items():          # spread labels so they do not overlap
    last = max(y, last + 0.17)
    placed[gname] = last
for gname, r in flat.iterrows():
    style = dict(color="0.6", lw=1, alpha=0.8) if r["used"] else \
        dict(color="tab:red", lw=1.8, ls="--")
    y = r[T.STAGES].values.astype(float)
    ax.plot(range(4), y, marker="o", ms=3, **style)
    ax.plot([3, 3.12], [y[-1], placed[gname]], color="0.8", lw=0.5)
    ax.text(3.15, placed[gname], gname, fontsize=7, va="center",
            color="tab:red" if not r["used"] else "0.25")
ax.set_xticks(range(4)); ax.set_xticklabels(T.STAGES)
ax.set_xlim(-0.2, 3.7)
ax.set_ylabel("median log2 expression")
ax.set_title(f"Pan-epithelial candidates (grey = used, range <= {FLAT_MAX:g} log2; red dashed = not used)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "pan_epithelial_flatness.png"), dpi=150)

# identity screen
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
cmap = {**T.COLORS, "LungMet": "#7b4fb3"}
for ax, organ, refx, midkey in [(axes[0], "hepatocyte", ref_l, "colon_vs_liver"),
                                (axes[1], "lung", ref_u, "colon_vs_lung")]:
    for gr, c in cmap.items():
        d = idq[idq.group == gr]
        ax.scatter(d[organ], d["colon"], s=10, c=c, label=gr, alpha=0.8)
    ex = idq[idq.excluded]
    ax.scatter(ex[organ], ex["colon"], s=80, facecolors="none", edgecolors="k", label="excluded")
    ax.scatter([refx[organ]], [refx["colon"]], marker="*", s=200, c="k",
               label=f"normal {'liver' if organ == 'hepatocyte' else 'lung'} (median)")
    ax.axhline(mid[midkey], color="0.5", ls=":")
    ax.axvline(mid[organ], color="0.5", ls=":")
    ax.set_xlabel(f"{organ} panel (mean log2)")
    ax.set_ylabel("colon epithelial panel (mean log2)")
    ax.legend(fontsize=7, loc="lower left")
axes[0].set_title("Tissue identity: colon vs liver")
axes[1].set_title("Tissue identity: colon vs lung")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "tissue_identity_screen.png"), dpi=150)

sc = results["main"][0]
cols = ["Epithelial", "Stromal", "Immune", "Liver_fraction"]
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
rng = np.random.default_rng(SEED)
for ax, c in zip(axes, cols):
    vals = [sc.loc[sc.stage == s, c].values for s in T.STAGES]
    bp = ax.boxplot(vals, patch_artist=True, widths=0.6, showfliers=False)
    for patch, s in zip(bp["boxes"], T.STAGES):
        patch.set_facecolor(T.COLORS[s]); patch.set_alpha(0.5)
    for i, v in enumerate(vals):
        ax.scatter(i + 1 + rng.uniform(-0.18, 0.18, len(v)), v, s=6, c=T.COLORS[T.STAGES[i]])
    ax.set_xticks(range(1, 5)); ax.set_xticklabels(T.STAGES, rotation=30)
    ax.set_title(c.replace("_", " ") + ("" if c == "Liver_fraction" else " score (ULM)"))
fig.tight_layout()
fig.savefig(os.path.join(OUT, "composition_by_stage_main.png"), dpi=150)

met = sc[sc.stage == "Metastasis"].sort_values("Liver_fraction")
fig, ax = plt.subplots(figsize=(10, 3.5))
ax.bar(range(len(met)), met["Liver_fraction"] * 100, color=T.COLORS["Metastasis"])
ax.set_xticks([]); ax.set_xlabel("liver metastasis samples (sorted)")
ax.set_ylabel("estimated liver tissue (%)")
ax.set_title("Estimated normal-liver content of each liver metastasis")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "liver_fraction_metastases.png"), dpi=150)

# ---------------------------------------------------------------- summary
from scipy.stats import kruskal  # noqa: E402


def sanity(sc, g):
    out = []
    for comp, gl in SANITY.items():
        gl = [x for x in gl if x in g.index]
        m = g.loc[gl].sub(g.loc[gl].mean(axis=1), axis=0).mean(axis=0)
        out.append(f"  {comp:<10} vs mean of {','.join(gl)}: Spearman r = "
                   f"{sc[comp].corr(m, method='spearman'):.2f}")
    return out


vers = {p: metadata.version(p) for p in ("decoupler", "pandas", "numpy", "scipy")}
fl = flat.copy()
fl[T.STAGES + ["range"]] = fl[T.STAGES + ["range"]].astype(float).round(2)
lines = ["STEP 2 - Tissue-composition scoring (GSE41258)", "=" * 50,
         "Inputs: data/processed/GSE41258_{main,withlung}_expr_genes.csv (+ samples.csv, step 1);",
         "  reference organ samples from data/raw/GSE41258_series_matrix.txt.gz",
         f"  ({liver_ref.shape[1]} normal liver, {lung_ref.shape[1]} normal lung; same probe per gene as step 1)",
         f"Seed: {SEED}", "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()), "",
         "TISSUE-IDENTITY SCREEN (all 358 samples: 338 main + 20 lung metastases)",
         "  Panels: " + "; ".join(f"{k}: {', '.join(x for x in v if x in g_all.index)}"
                                  for k, v in PANELS.items()),
         f"  Normal colon median: colon {ref_c['colon']:.2f}, hepatocyte {ref_c['hepatocyte']:.2f}, "
         f"lung {ref_c['lung']:.2f}",
         f"  Normal liver median: colon {ref_l['colon']:.2f}, hepatocyte {ref_l['hepatocyte']:.2f}; "
         f"normal lung median: colon {ref_u['colon']:.2f}, lung {ref_u['lung']:.2f}",
         "  Rule: low colon = below midpoint between normal colon and normal liver/lung (log2);",
         "        high liver/lung = above midpoint between normal colon and that organ.",
         "        MISMATCH = low colon AND high other-organ (sample closer to liver/lung than colon).",
         "  Calls:", *[f"    {k:<42} {v}" for k, v in idq["call"].value_counts().items()],
         "  Samples not called 'ok' or 'host organ present':"]
odd = idq[~idq["call"].isin(["ok", "host organ present (expected for site)"])]
for s, r in odd.iterrows():
    lines.append(f"    {s} [{r.group}] {r.call}: colon {r.colon:.2f}, hepatocyte {r.hepatocyte:.2f}, "
                 f"lung {r.lung:.2f}{'  -> EXCLUDED' if r.excluded else ''}")
lines += ["", "EXCLUSIONS (data/processed/GSE41258_exclusions.csv):",
          "  " + EXCLUSION_STATEMENT,
          *[f"  {k}: {v}" for k, v in EXCLUDE.items()],
          "FLAGGED BY THE SCREEN, REVIEWED AND KEPT:",
          *[f"  {k}: {v}" for k, v in KEPT_NOTES.items()], "",
          "ONE ARRAY PER PATIENT PER STAGE (de-duplication):",
          "  " + DEDUP_STATEMENT,
          f"  Patient ids from the official GEO 'patient id' field: "
          f"{(sets['main_allarrays'][1]['patient_id'] != '').sum()}/{len(sets['main_allarrays'][1])} samples, "
          f"{sets['main_allarrays'][1]['patient_id'].nunique()} patients (main_allarrays); the pilot's "
          f"title regex found {sets['main_allarrays'][1]['patient_pilot_regex'].notna().sum()}.",
          *[f"  removed {r['sample']} [{r['step1_stage']}]: " + r["reason"].split("). One array")[0] + ")"
            for r in excl_df[excl_df.category == "same patient and stage"].to_dict("records")],
          "  Removed repeat arrays per stage: " + ", ".join(
              f"{s} {v}" for s, v in excl_df[excl_df.category == "same patient and stage"]["step1_stage"]
              .value_counts().reindex(T.STAGES + ["Metastasis (lung)"]).dropna().astype(int).items()),
          "  Authors' 'included in analysis' flag in the primary set (Yes/No/NA): " + ", ".join(
              f"{s} {'/'.join(str(int(((sets['main'][1].stage == s) & (sets['main'][1].authors_included == f)).sum())) for f in ('Yes', 'No', ''))}"
              for s in T.STAGES), "",
          "SETTINGS",
          f"  Epithelial: pan-epithelial structural candidates, used if max-min of stage medians <= {FLAT_MAX:g} log2",
          f"    candidates missing from platform: {', '.join(missing_epi) or 'none'}",
          "  Stromal / immune: PanglaoDB human canonical markers (decoupler dc.op.resource);",
          *[f"    {k}: {', '.join(v)}" for k, v in COMPARTMENTS.items()],
          "    genes in >1 compartment removed",
          "  Genes per compartment: " + ", ".join(f"{k} {v}" for k, v in set_sizes.items()),
          f"  Epithelial genes used: {', '.join(epi_genes)}",
          f"  Strict epithelial set (sensitivity, range <= {FLAT_STRICT:g}): {len(epi_strict)} genes, "
          f"drops {', '.join(x for x in epi_genes if x not in epi_strict)}",
          "  Scoring: decoupler ULM (tmin=5) on gene-centred log2 expression",
          f"  Organ signature: organ mean - max(normal colon, carcinoma mean) >= {ORGAN_LFC:g} log2, "
          f"organ mean >= {ORGAN_MIN:g}; metastases not used",
          "  Organ fraction: carcinoma median + organ median mixture (linear scale), per-gene fraction,",
          "    median across signature genes, clipped 0-1",
          "  Adjustment: per gene OLS on all composition covariates; adjusted = residual + gene mean", "",
          "PAN-EPITHELIAL MARKER MEDIANS PER STAGE (main set, log2):",
          fl[T.STAGES + ["range", "used", "role"]].to_string(), ""]
for name in sets:
    sc, liv_tab, lung_tab, r2, g, epi_corr = results[name]
    cov = [c for c in sc.columns if c not in ("stage", "raw_label")]
    cnt = sc["stage"].value_counts().reindex(T.STAGES)
    lines += [f"--- set {name} ---",
              "Samples per stage: " + ", ".join(f"{s} {cnt[s]}" for s in T.STAGES) + f" (total {cnt.sum()})",
              f"Liver signature genes: {len(liv_tab)} (top: {', '.join(liv_tab.index[:10])})"]
    if lung_tab is not None:
        lines.append(f"Lung signature genes: {len(lung_tab)} (top: {', '.join(lung_tab.index[:10])})")
    lines += ["Median score per stage:", sc.groupby("stage")[cov].median().reindex(T.STAGES)
              .round(3).to_string()]
    for c in cov:
        p = kruskal(*[sc.loc[sc.stage == s, c] for s in T.STAGES]).pvalue
        lines.append(f"  Kruskal-Wallis across stages, {c}: p = {p:.1e}")
    lines += ["Sanity check (score vs classic markers):", *sanity(sc, g)]
    m = sc[sc.stage == "Metastasis"]
    lines.append("Liver fraction, metastases: median {:.1%}, IQR {:.1%}-{:.1%}, max {:.1%}, n > 20%: {}; "
                 "other stages max {:.1%}".format(m.Liver_fraction.median(), m.Liver_fraction.quantile(.25),
                                                  m.Liver_fraction.quantile(.75), m.Liver_fraction.max(),
                                                  int((m.Liver_fraction > .2).sum()),
                                                  sc.loc[sc.stage != "Metastasis", "Liver_fraction"].max()))
    if lung_tab is not None:
        lm = sc[sc.raw_label.str.lower().str.contains("lung")]
        vm = sc[(sc.stage == "Metastasis") & ~sc.raw_label.str.lower().str.contains("lung")]
        lines.append("Lung fraction: lung metastases median {:.1%} (max {:.1%}); liver metastases median {:.1%}"
                     .format(lm.Lung_fraction.median(), lm.Lung_fraction.max(), vm.Lung_fraction.median()))
    lines += [f"Epithelial score, {len(epi_genes)}-gene (used) vs {len(epi_strict)}-gene strict "
              f"(range <= {FLAT_STRICT:g}): Pearson r = {epi_corr[0]:.3f}, Spearman r = {epi_corr[1]:.3f}",
              f"Variance explained by composition (R2 per gene): median {r2.median():.2f}, "
              f"genes with R2 > 0.5: {int((r2 > 0.5).sum())}/{len(r2)}", ""]
lines += ["Outputs:",
          "  results/step02/tissue_identity_qc.csv, tissue_identity_screen.png",
          "  results/step02/pan_epithelial_marker_flatness.csv, pan_epithelial_flatness.png",
          "  results/step02/composition_scores_<set>.csv        per-sample scores",
          "  results/step02/marker_sets_used.csv                 marker genes per compartment",
          "  results/step02/liver_signature_<set>.csv, lung_signature_withlung.csv",
          "  results/step02/variance_explained_by_composition_<set>.csv",
          "  results/step02/composition_by_stage_main.png, liver_fraction_metastases.png",
          "  data/processed/GSE41258_exclusions.csv",
          "  data/processed/GSE41258_<set>_samples_qc.csv       samples per analysis set",
          "  data/processed/GSE41258_<set>_expr_genes_adjusted.csv  composition-adjusted expression",
          "  (unadjusted expression for each set = step-1 file restricted to that set's samples)"]
with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))

