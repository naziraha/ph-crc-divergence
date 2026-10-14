"""
Step 1 for GSE117606, exactly as addendum 6 (commit 93613bd).

* Series matrix read with tda_crc_pilot.read_series_matrix; log2(max(x,1)) if the 99th percentile > 50.
* Gene IDs (GPL25373 = Brainarray ENTREZG CDF, rows '<Entrez>_at'; the GPL table has no symbol column):
  remove AFFX controls -> Entrez -> official symbol via NCBI Homo_sapiens.gene_info (download date and
  md5 recorded) -> drop unmapped -> several IDs per symbol: keep the highest mean (main set).
* Main set: Adjacent Mucosa -> Normal, Adenoma -> Adenoma, Tumor -> Carcinoma; one array per patient per
  stage; A-coded carcinomas GSM3304801 / GSM3304813 / GSM3304818 set aside; 10 SSA kept separately.
* Section 7.2 array-quality rule (main set, all mapped genes): flag if mean Pearson correlation with all
  other arrays < Q1 - 1.5 IQR, or median or IQR more than 3 robust SD (1.4826 MAD) from the cohort median.
  Flagged arrays stay in the primary analysis.
"""
import datetime as dt
import gzip
import hashlib
import io
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402

RAW, PROC = os.path.join(ROOT, "data", "raw"), os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "step01", "GSE117606")
os.makedirs(OUT, exist_ok=True)
SET_ASIDE = {"GSM3304801": "1889A", "GSM3304813": "2488A", "GSM3304818": "1818A"}
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]
EXPECTED = {"Normal": 65, "Adenoma": 59, "Carcinoma": 71}


def fetch(url, path):
    if not os.path.exists(path):
        r = requests.get(url, timeout=900)
        r.raise_for_status()
        open(path, "wb").write(r.content)
    b = open(path, "rb").read()
    return len(b), hashlib.md5(b).hexdigest(), dt.datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M")


sm_info = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117606/matrix/GSE117606_series_matrix.txt.gz",
                os.path.join(RAW, "GSE117606_series_matrix.txt.gz"))
gi_info = fetch("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                os.path.join(RAW, "Homo_sapiens.gene_info.gz"))

expr, meta = T.read_series_matrix(os.path.join(RAW, "GSE117606_series_matrix.txt.gz"))
samples = list(expr.columns)
n_rows0 = expr.shape[0]
n_na = int(expr.isna().any(axis=1).sum())
p99 = float(np.nanpercentile(expr.values, 99))
log2_applied = p99 > 50
if log2_applied:
    expr = np.log2(expr.clip(lower=1))

# ---- sample table from the series-matrix metadata
ch = {}
for row in meta["!Sample_characteristics_ch1"]:
    for i, v in enumerate(row):
        k, _, val = v.partition(":")
        ch.setdefault(k.strip().lower(), [""] * len(samples))[i] = val.strip()
lab = pd.DataFrame(ch, index=samples)
lab.insert(0, "title", meta["!Sample_title"][0])
lab["stage"] = lab["tissue"].map({"Adjacent Mucosa": "Normal", "Adenoma": "Adenoma", "Tumor": "Carcinoma", "SSA": "SSA"})
lab["set"] = np.where(lab.index.isin(list(SET_ASIDE)), "set aside (A-coded carcinoma)",
                      np.where(lab.stage == "SSA", "SSA (descriptive)", "main"))
hdr = pd.read_csv(os.path.join(RAW, "GSE117606_cel_headers.csv")).set_index("gsm")
lab = lab.join(hdr[["affymetrix-ArrayPlateGUID", "affymetrix-scan-date"]].rename(
    columns={"affymetrix-ArrayPlateGUID": "plate", "affymetrix-scan-date": "scan_date"}))
main_ids = list(lab.index[lab.set == "main"])
dups = lab.loc[main_ids].groupby(["patient id", "stage"]).size()
cnt = lab.loc[main_ids, "stage"].value_counts().to_dict()
if cnt != EXPECTED or (dups > 1).any():
    raise SystemExit(f"Unexpected main-set counts {cnt} or duplicates {int((dups > 1).sum())}")

# ---- gene-ID mapping rule
ids = pd.Series(expr.index.astype(str), index=expr.index)
is_affx = ids.str.startswith("AFFX")
entrez = ids[~is_affx].str.replace("_at", "", regex=False)
gi = pd.read_csv(os.path.join(RAW, "Homo_sapiens.gene_info.gz"), sep="\t", usecols=["GeneID", "Symbol"], dtype=str)
sym_map = dict(zip(gi.GeneID, gi.Symbol))
symbols = entrez.map(sym_map)
unmapped = symbols.isna()
mapped = symbols[~unmapped]
e_main = expr.loc[mapped.index, main_ids]
mean = e_main.mean(axis=1)
pick = (pd.DataFrame({"symbol": mapped, "mean": mean}).sort_values("mean", ascending=False)
        .groupby("symbol", sort=True).head(1))
n_dup_ids = int(len(mapped) - len(pick))
genes_all = expr.loc[pick.index]
genes_all.index = pick["symbol"].values
genes_all = genes_all.sort_index()
genes_all.index.name = "gene"
pick.rename_axis("probe_id").reset_index().assign(entrez=lambda d: d.probe_id.str.replace("_at", "", regex=False)) \
    .to_csv(os.path.join(OUT, "probe_to_gene.csv"), index=False)

main = genes_all[main_ids]
main.to_csv(os.path.join(PROC, "GSE117606_main_expr_genes.csv"), float_format="%.5f")
other_ids = list(lab.index[lab.set != "main"])
genes_all[other_ids].to_csv(os.path.join(PROC, "GSE117606_other_expr_genes.csv"), float_format="%.5f")
cols = ["title", "patient id", "tissue", "stage", "set", "adenomadiagnosis", "tnm", "uicc", "age", "gender",
        "plate", "scan_date"]
lab.loc[main_ids, cols].to_csv(os.path.join(PROC, "GSE117606_main_samples.csv"), index_label="sample")
lab.loc[other_ids, cols].to_csv(os.path.join(PROC, "GSE117606_other_samples.csv"), index_label="sample")
lab[cols].to_csv(os.path.join(OUT, "all_samples_labels.csv"), index_label="sample")

# ---- section 7.2 array-quality rule (main set)
C = np.corrcoef(main.T.values)
np.fill_diagonal(C, np.nan)
mcorr = pd.Series(np.nanmean(C, axis=1), index=main_ids)
q1, q3 = mcorr.quantile(0.25), mcorr.quantile(0.75)
thr_corr = q1 - 1.5 * (q3 - q1)
med = main.median(axis=0)
iqr = main.quantile(0.75, axis=0) - main.quantile(0.25, axis=0)


def rz(x):
    return (x - x.median()) / (1.4826 * (x - x.median()).abs().median())


qa = pd.DataFrame({"stage": lab.loc[main_ids, "stage"], "plate": lab.loc[main_ids, "plate"],
                   "mean_corr": mcorr.round(4), "median": med.round(3), "IQR": iqr.round(3),
                   "median_robust_z": rz(med).round(2), "IQR_robust_z": rz(iqr).round(2)})
qa["flag_corr"] = qa.mean_corr < thr_corr
qa["flag_median"] = qa.median_robust_z.abs() > 3
qa["flag_IQR"] = qa.IQR_robust_z.abs() > 3
qa["flagged"] = qa.flag_corr | qa.flag_median | qa.flag_IQR
qa.to_csv(os.path.join(OUT, "array_quality.csv"), index_label="sample")
fl = qa[qa.flagged]

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
fig, axes = plt.subplots(1, 2, figsize=(13, 4))
order = qa.sort_values(["stage", "mean_corr"]).index
colors = {"Normal": T.COLORS["Normal"], "Adenoma": T.COLORS["Adenoma"], "Carcinoma": T.COLORS["Carcinoma"]}
axes[0].scatter(range(len(order)), qa.loc[order, "mean_corr"], c=[colors[s] for s in qa.loc[order, "stage"]], s=12)
axes[0].axhline(thr_corr, color="k", ls=":", lw=1)
axes[0].set_title("Mean correlation with all other arrays (dotted: Q1 - 1.5 IQR)")
axes[0].set_xticks([])
axes[1].scatter(qa["median"], qa["IQR"], c=[colors[s] for s in qa.stage], s=12)
for s, r in fl.iterrows():
    axes[1].annotate(s, (r["median"], r["IQR"]), fontsize=7)
axes[1].set_xlabel("array median (log2)"); axes[1].set_ylabel("array IQR (log2)")
axes[1].set_title("Array median vs IQR (flagged arrays labelled)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "array_quality.png"), dpi=140)

vers = {p: metadata.version(p) for p in ("pandas", "numpy", "requests")}
L = ["STEP 1 - Data preparation (GSE117606, GPL25373)", "=" * 55,
     "Pre-specification: docs/addendum6_GSE117606.md (commit 93613bd)",
     f"Series matrix: data/raw/GSE117606_series_matrix.txt.gz ({sm_info[0] / 1e6:.1f} MB, md5 {sm_info[1]}, obtained {sm_info[2]})",
     f"Gene symbols: NCBI Homo_sapiens.gene_info.gz ({gi_info[0] / 1e6:.1f} MB, md5 {gi_info[1]}, downloaded {gi_info[2]})",
     "Versions: " + ", ".join(f"{k} {v}" for k, v in vers.items()) + f", python {sys.version.split()[0]}; seed 42 (not used)", "",
     f"Arrays in series matrix: {len(samples)}",
     "Sample sets:", "  " + pd.crosstab(lab.tissue, lab.set).to_string().replace("\n", "\n  "),
     "Main set per stage: " + ", ".join(f"{s} {cnt[s]}" for s in STAGES3) + " (expected 65/59/71; verified); "
     "patient-stage duplicates: 0",
     "Set aside (A-coded carcinomas): " + ", ".join(f"{k} ({v})" for k, v in SET_ASIDE.items()),
     f"SSA (descriptive): {int((lab.stage == 'SSA').sum())} arrays from {lab.loc[lab.stage == 'SSA', 'patient id'].nunique()} patients", "",
     "GENE-ID MAPPING (rows dropped at each stage)",
     f"  rows in series matrix: {n_rows0} (rows with missing values: {n_na})",
     f"  1. AFFX control probe sets removed: {int(is_affx.sum())}  -> {int((~is_affx).sum())}",
     f"  2. Entrez IDs without a current NCBI symbol (dropped): {int(unmapped.sum())}  -> {len(mapped)}",
     f"  3. several IDs mapping to the same symbol, lower-mean IDs dropped: {n_dup_ids}  -> {len(pick)} genes",
     f"Scale: 99th percentile = {p99:.2f}; log2 applied: {log2_applied}",
     f"Per-array median log2 (main set): IQR {med.quantile(.25):.2f}-{med.quantile(.75):.2f}, range {med.min():.2f}-{med.max():.2f}", "",
     "ARRAY-QUALITY RULE (section 7.2; main set, all mapped genes)",
     f"  mean correlation with other arrays: median {mcorr.median():.3f}, range {mcorr.min():.3f}-{mcorr.max():.3f}; "
     f"threshold Q1 - 1.5 IQR = {thr_corr:.3f}",
     f"  flagged arrays: {len(fl)} (correlation {int(qa.flag_corr.sum())}, median {int(qa.flag_median.sum())}, IQR {int(qa.flag_IQR.sum())})"]
for s, r in fl.iterrows():
    reasons = [n for n, f in (("low correlation", r.flag_corr), ("median", r.flag_median), ("IQR", r.flag_IQR)) if f]
    L.append(f"    {s} [{r.stage}; {lab.loc[s, 'title']}] {', '.join(reasons)}: mean r {r.mean_corr:.3f}, "
             f"median z {r.median_robust_z:.1f}, IQR z {r.IQR_robust_z:.1f}")
L += [f"  flagged per stage: {fl.stage.value_counts().reindex(STAGES3).fillna(0).astype(int).to_dict()}",
      "  Flagged arrays stay in the primary analysis; sensitivity analysis without them (section 7.2).", "",
      "Outputs:",
      "  data/processed/GSE117606_main_expr_genes.csv / _main_samples.csv   main set (genes x 195 arrays)",
      "  data/processed/GSE117606_other_expr_genes.csv / _other_samples.csv SSA (10) + set-aside (3), descriptive",
      "  results/step01/GSE117606/probe_to_gene.csv, all_samples_labels.csv, array_quality.csv, array_quality.png"]
open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
