"""
Addendum 6 section 7.1 (commit 93613bd) - GSE117606 batch check: plate and scan day vs stage.

Stages from GEO metadata (tissue): Adjacent Mucosa -> Normal, Adenoma -> Adenoma, Tumor -> Carcinoma;
main set = one array per patient per stage (A-coded carcinomas of patients 1818, 1889, 2488 set aside;
SSA excluded from the main set, shown for description). Batch units: plate (ArrayPlateGUID, plate
barcode) and scan day (CEL header). Measures (as addendum 5 part 2): Cramer's V (batch x stage);
per transition AUC of scan date (Mann-Whitney) and overlap = share of the two stages' arrays on
plates (days) containing both stages. Confounded if AUC >= 0.80 or <= 0.20, or overlap < 20%.
"""
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "addendum6", "batch_check")
os.makedirs(OUT, exist_ok=True)
STAGES = ["Normal", "Adenoma", "Carcinoma"]
SET_ASIDE = {"GSM3304801", "GSM3304813", "GSM3304818"}

# --- GEO sample metadata (titles / characteristics only)
meta, cur = [], None
for line in open(os.path.join(ROOT, "data", "raw", "GSE117606_samples_brief.txt"), encoding="utf-8"):
    if line.startswith("^SAMPLE"):
        cur = {"gsm": line.split("=")[1].strip()}
        meta.append(cur)
    elif cur is not None and line.startswith("!Sample_title"):
        cur["title"] = line.partition(" = ")[2].strip()
    elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
        k, _, v = line.partition(" = ")[2].partition(":")
        cur[k.strip().lower()] = v.strip()
meta = pd.DataFrame(meta).set_index("gsm")
stage_map = {"Adjacent Mucosa": "Normal", "Adenoma": "Adenoma", "Tumor": "Carcinoma", "SSA": "SSA"}
meta["stage"] = meta["tissue"].map(stage_map)
meta["set"] = np.where(meta.index.isin(SET_ASIDE), "set aside (A-coded carcinoma)",
                       np.where(meta.stage == "SSA", "SSA (descriptive)", "main"))

hdr = pd.read_csv(os.path.join(ROOT, "data", "raw", "GSE117606_cel_headers.csv")).set_index("gsm")
d = meta.join(hdr, how="left")
d["scan_dt"] = pd.to_datetime(d["affymetrix-scan-date"], utc=True)
d["scan_day"] = d.scan_dt.dt.strftime("%Y-%m-%d")
d["plate"] = d["affymetrix-ArrayPlateGUID"]
d[["title", "patient id", "tissue", "stage", "set", "scan_day", "plate", "affymetrix-plate-barcode",
   "affymetrix-plate-peg-wellposition"]].to_csv(os.path.join(OUT, "sample_batch_table.csv"), index_label="gsm")

main = d[d.set == "main"]
# duplicates check (one per patient per stage)
dup = main.groupby(["patient id", "stage"]).size()
L = ["ADDENDUM 6 section 7.1 - GSE117606 batch check (descriptive)", "=" * 60,
     "Pre-specified: docs/addendum6_GSE117606.md section 7.1 (commit 93613bd). CEL headers only (no intensities).",
     "Confounded if AUC >= 0.80 or <= 0.20, or overlap < 20%; Cramer's V >= 0.5 = strong association.", "",
     f"CEL headers matched to GEO samples: {int(d['affymetrix-scan-date'].notna().sum())}/{len(d)}",
     "Main set per stage: " + ", ".join(f"{s} {int((main.stage == s).sum())}" for s in STAGES)
     + f"; patient-stage duplicates in main set: {int((dup > 1).sum())}",
     f"Set aside: {', '.join(sorted(SET_ASIDE))}; SSA (descriptive): {int((d.stage == 'SSA').sum())}",
     f"Plates (ArrayPlateGUID): {d.plate.nunique()} (barcodes: {d['affymetrix-plate-barcode'].nunique()}); "
     f"scan days: {d.scan_day.nunique()}; scan period {d.scan_dt.min():%Y-%m-%d} to {d.scan_dt.max():%Y-%m-%d}", ""]
rows = []
for unit, col in (("plate", "plate"), ("scan day", "scan_day")):
    tab = pd.crosstab(main[col], main.stage).reindex(columns=STAGES).fillna(0).astype(int)
    chi2, p, dof, _ = chi2_contingency(tab.values)
    V = float(np.sqrt(chi2 / (tab.values.sum() * (min(tab.shape) - 1))))
    tab_all = pd.crosstab(d[col], d.stage.where(d.set != "set aside (A-coded carcinoma)", "set aside")) \
        .reindex(columns=STAGES + ["SSA", "set aside"]).fillna(0).astype(int)
    tab_all.to_csv(os.path.join(OUT, f"{unit.replace(' ', '_')}_x_stage.csv"))
    L += [f"--- {unit} x stage (main set; SSA and set-aside shown for description) ---",
          "  " + tab_all.to_string().replace("\n", "\n  "),
          f"  Cramer's V ({unit} x stage, main set) = {V:.2f} (chi2 = {chi2:.1f}, df = {dof}, p = {p:.2g}) -> "
          + ("strong" if V >= 0.5 else "moderate" if V >= 0.2 else "weak") + " association"]
    for s, t in (("Normal", "Adenoma"), ("Adenoma", "Carcinoma")):
        a, b = main[main.stage == s], main[main.stage == t]
        both = set(a[col]) & set(b[col])
        ov = float((a[col].isin(both).sum() + b[col].isin(both).sum()) / (len(a) + len(b)))
        xa = a.scan_dt.map(lambda x: x.timestamp()).values
        xb = b.scan_dt.map(lambda x: x.timestamp()).values
        auc = float(mannwhitneyu(xb, xa, alternative="two-sided").statistic / (len(xa) * len(xb)))
        conf = auc >= 0.80 or auc <= 0.20 or ov < 0.20
        rows.append(dict(unit=unit, transition=f"{s}->{t}", n_from=len(a), n_to=len(b), shared=len(both),
                         overlap=ov, auc_scan_date=auc, cramers_v=V, confounded=conf))
        L.append(f"  {s:>9} -> {t:<9} shared {unit}s {len(both):>2}; overlap {ov * 100:5.1f}%; AUC (scan date) {auc:.2f}"
                 f" -> {'CONFOUNDED' if conf else 'not confounded'}")
    L.append("")
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "batch_confounding.csv"), index=False)
# same-stage batch split availability (default handling, section 7.1 point 3)
L.append("Same-stage splits with >= 13 arrays on each side (for a batch-only contrast, if needed):")
for s in STAGES:
    x = main[main.stage == s].plate.value_counts()
    L.append(f"  {s}: arrays per plate {x.to_dict()}")
L += ["", "Outputs: results/addendum6/batch_check/sample_batch_table.csv, plate_x_stage.csv, scan_day_x_stage.csv,",
      "  batch_confounding.csv"]
open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
