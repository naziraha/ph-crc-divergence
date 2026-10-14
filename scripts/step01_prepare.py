"""
Step 1 - Data preparation for GSE41258 (GPL96).

  1. Read the series matrix and assign stages with the validated pilot code
     (src/tda_crc_pilot.py: read_series_matrix, labels_from_meta).
     Main version: metastasis = liver metastases; lung metastases, normal liver /
     lung and cell lines excluded. Sensitivity version: lung metastases added.
  2. Map probes to gene symbols with the GPL96 annotation (GEOparse). Probes with
     no symbol or with several symbols ("A /// B") are dropped.
  3. log2 if not already on log scale (same rule as the pilot).
  4. One probe per gene: the probe with the highest mean expression.

Outputs: data/processed/GSE41258_*.csv, results/step01/*
"""
import os
import sys
import re
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402
import GEOparse  # noqa: E402

SEED = 42
np.random.seed(SEED)
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "step01")
os.makedirs(PROC, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
MATRIX = os.path.join(RAW, "GSE41258_series_matrix.txt.gz")
EXPECTED = {"Normal": 54, "Adenoma": 51, "Carcinoma": 186, "Metastasis": 47}


def excluded_type(label):
    t = label.lower()
    for pattern, name in [(r"cell line", "cell line"), (r"normal liver", "normal liver"),
                          (r"normal lung", "normal lung"), (r"lung metastasis", "lung metastasis")]:
        if re.search(pattern, t):
            return name
    return "other / unmatched: " + label


def collapse_to_genes(expr, probe2gene):
    """expr: probes x samples (log scale). Keep highest-mean probe per gene."""
    x = expr.loc[expr.index.intersection(probe2gene.index)]
    mean = x.mean(axis=1)
    pick = (pd.DataFrame({"gene": probe2gene.loc[x.index], "mean": mean})
            .sort_values("mean", ascending=False)
            .groupby("gene", sort=True).head(1))
    genes = x.loc[pick.index]
    genes.index = pick["gene"].values
    return genes.sort_index(), pick.rename_axis("probe").reset_index()


# ---------------------------------------------------------------- 1. load + labels
expr, meta = T.read_series_matrix(MATRIX)
samples = list(expr.columns)
print(f"[load] {expr.shape[0]} probes x {expr.shape[1]} samples")
stages_main, raw_labels = T.labels_from_meta(meta, samples, keep_micro=True, keep_lung=False)
stages_lung, _ = T.labels_from_meta(meta, samples, keep_micro=True, keep_lung=True)
patient = T.patients_from_meta(meta, samples)
label_field = [k for k in ("!Sample_source_name_ch1", "!Sample_characteristics_ch1",
                           "!Sample_title", "!Sample_description")
               if any(list(r) == list(raw_labels.values) for r in meta.get(k, []))][0]

counts_main = stages_main.value_counts().reindex(T.STAGES).fillna(0).astype(int)
counts_lung = stages_lung.value_counts().reindex(T.STAGES).fillna(0).astype(int)
excl = raw_labels[stages_main.isna()].map(excluded_type)
excl_counts = excl.value_counts()
micro_n = int(raw_labels[stages_main == "Adenoma"].str.lower().str.contains("microadenoma").sum())
print("[labels] main:\n" + counts_main.to_string())
print("[labels] excluded:\n" + excl_counts.to_string())

# ---------------------------------------------------------------- 2. platform annotation
gpl = GEOparse.get_GEO(geo="GPL96", destdir=RAW, silent=True)
ann = gpl.table[["ID", "Gene Symbol"]].copy()
ann["ID"] = ann["ID"].astype(str)
ann["Gene Symbol"] = ann["Gene Symbol"].fillna("").astype(str).str.strip()
n_probes_platform = len(ann)
no_symbol = ann["Gene Symbol"] == ""
multi = ann["Gene Symbol"].str.contains("///")
probe2gene = ann.loc[~no_symbol & ~multi].set_index("ID")["Gene Symbol"]

# ---------------------------------------------------------------- 3. log2 check
n_na_probes = int(expr.isna().any(axis=1).sum())
expr = expr.dropna(axis=0, how="any")
p99_before = float(np.nanpercentile(expr.values, 99))
log2_applied = p99_before > 50          # same rule as the pilot (preprocess)
if log2_applied:
    expr = np.log2(expr.clip(lower=1))
print(f"[log2] 99th percentile before = {p99_before:.1f}; log2 applied: {log2_applied}")


# ---------------------------------------------------------------- 4. collapse + save
def build(stages, tag):
    keep = stages.dropna()
    keep = keep[keep.isin(T.STAGES)]
    e = expr[keep.index]
    genes, pick = collapse_to_genes(e, probe2gene)
    genes.index.name = "gene"
    genes.to_csv(os.path.join(PROC, f"GSE41258_{tag}_expr_genes.csv"), float_format="%.5f")
    lab = pd.DataFrame({"sample": keep.index, "stage": keep.values,
                        "raw_label": raw_labels.loc[keep.index].values,
                        "patient": patient.reindex(keep.index).values if patient is not None else None})
    lab.to_csv(os.path.join(PROC, f"GSE41258_{tag}_samples.csv"), index=False)
    pick.to_csv(os.path.join(OUT, f"probe_to_gene_{tag}.csv"), index=False)
    return genes, pick, lab


genes_main, pick_main, lab_main = build(stages_main, "main")
genes_lung, pick_lung, lab_lung = build(stages_lung, "withlung")
same_probe = pick_main.set_index("gene")["probe"].eq(
    pick_lung.set_index("gene")["probe"].reindex(pick_main["gene"])).sum()

# all samples incl. excluded, for traceability
pd.DataFrame({"sample": samples, "raw_label": raw_labels.values,
              "stage_main": stages_main.values, "stage_withlung": stages_lung.values,
              "patient": patient.values if patient is not None else None}) \
    .to_csv(os.path.join(OUT, "all_samples_labels.csv"), index=False)
excl_counts.rename_axis("excluded_type").reset_index(name="n") \
    .to_csv(os.path.join(OUT, "excluded_samples.csv"), index=False)

# QC figure: per-sample expression distribution by stage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
order = lab_main.sort_values("stage", key=lambda s: s.map(T.STAGES.index))["sample"]
fig, ax = plt.subplots(figsize=(14, 4))
q = genes_main[order].quantile([0.05, 0.25, 0.5, 0.75, 0.95])
xs = np.arange(len(order))
cols = [T.COLORS[s] for s in lab_main.set_index("sample").loc[order, "stage"]]
ax.vlines(xs, q.loc[0.05], q.loc[0.95], color=cols, lw=0.6, alpha=0.5)
ax.vlines(xs, q.loc[0.25], q.loc[0.75], color=cols, lw=2)
ax.scatter(xs, q.loc[0.5], s=3, c="k")
for s in T.STAGES:
    ax.plot([], [], color=T.COLORS[s], lw=3, label=f"{s} (n={counts_main[s]})")
ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.15), frameon=False)
ax.set_xlabel("samples (ordered by stage)")
ax.set_ylabel("log2 expression")
ax.set_xticks([])
ax.set_title("GSE41258 per-sample expression distribution (median, IQR, 5-95%)", pad=28)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "qc_sample_distributions.png"), dpi=150)

# ---------------------------------------------------------------- summary
med = genes_main.median(axis=0)
outl = med[(med - med.median()).abs() > 3 * (med - med.median()).abs().median() * 1.4826]
vers = {p: metadata.version(p) for p in ("GEOparse", "pandas", "numpy", "matplotlib")}
ok = all(counts_main[s] == EXPECTED[s] for s in EXPECTED)
lines = [
    "STEP 1 - Data preparation (GSE41258, GPL96)",
    "=" * 50,
    f"Input: {os.path.relpath(MATRIX, ROOT)}; platform annotation GPL96 (GEOparse, data/raw/GPL96.txt)",
    "Labels: src/tda_crc_pilot.py labels_from_meta (keep_micro=True)",
    f"  metadata field used: {label_field}",
    f"Seed: {SEED}",
    "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()) + f", python {sys.version.split()[0]}",
    "",
    f"Total samples in series matrix: {len(samples)}",
    "Samples per stage, MAIN (metastasis = liver):",
    *[f"  {s:<11} {counts_main[s]:>4}   (expected {EXPECTED[s]})" for s in T.STAGES],
    f"  total       {counts_main.sum():>4}",
    f"  of which adenomas labelled microadenoma: {micro_n}",
    f"  matches expected counts: {ok}",
    "Excluded samples (main version):",
    *[f"  {k:<20} {v:>4}" for k, v in excl_counts.items()],
    f"  total excluded {excl_counts.sum()}",
    "Samples per stage, WITH-LUNG sensitivity version:",
    *[f"  {s:<11} {counts_lung[s]:>4}" for s in T.STAGES],
    f"Patient ids parsed: {lab_main['patient'].notna().sum()}/{len(lab_main)} samples, "
    f"{lab_main['patient'].nunique()} patients (main version)",
    "",
    "Probe -> gene mapping:",
    f"  probes on platform: {n_probes_platform}; in matrix: {expr.shape[0] + n_na_probes}",
    f"  probes with missing values dropped: {n_na_probes}",
    f"  probes without gene symbol (incl. AFFX controls): {int(no_symbol.sum())}",
    f"  probes mapping to >1 gene ('///'), dropped: {int((multi & ~no_symbol).sum())}",
    f"  single-gene probes used: {len(probe2gene)}",
    f"  genes after keeping highest-mean probe: main {genes_main.shape[0]}, with-lung {genes_lung.shape[0]}",
    f"  genes with the same chosen probe in both versions: {same_probe}/{len(pick_main)}",
    f"Scale: 99th percentile before transform = {p99_before:.1f}; log2 applied: {log2_applied}",
    f"After processing: median log2 expression {float(np.median(genes_main.values)):.2f}, "
    f"range {float(genes_main.values.min()):.2f} to {float(genes_main.values.max()):.2f}",
    f"Per-sample medians: IQR {med.quantile(0.25):.2f}-{med.quantile(0.75):.2f}; "
    f"largest deviation from cohort median {float((med - med.median()).abs().max()):.2f} log2 units",
    f"  samples >3 robust SD from cohort median: {len(outl)}"
    + (f" ({', '.join(outl.index)}; deviations of <0.1 log2 units, kept)" if len(outl) else ""),
    "Note: gene symbols are those in the GPL96 annotation; harmonise with GPL570 /",
    "  RNA-seq symbols before cross-cohort comparison (step 11).",
    "",
    "Outputs:",
    "  data/processed/GSE41258_main_expr_genes.csv      genes x samples, log2 (main analysis)",
    "  data/processed/GSE41258_main_samples.csv         sample, stage, raw label, patient",
    "  data/processed/GSE41258_withlung_expr_genes.csv  same, lung metastases added (sensitivity)",
    "  data/processed/GSE41258_withlung_samples.csv",
    "  results/step01/all_samples_labels.csv           every sample incl. excluded",
    "  results/step01/excluded_samples.csv",
    "  results/step01/probe_to_gene_*.csv              probe chosen per gene",
    "  results/step01/qc_sample_distributions.png",
]
with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
if not ok:
    print("\nWARNING: stage counts differ from expected - check before continuing.")
