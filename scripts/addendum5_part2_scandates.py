"""
Addendum 5, part 2 (commit d5673a7) - GSE41258 scan date vs stage (batch check, descriptive).

Scan date from the GEO characteristic "scan date": dd/mm/yyyy, or an integer = Excel serial date
(days since 1899-12-30). Sets main (primary) and main_allarrays (reference).
Measures fixed in advance: Cramer's V (scan year x stage); per adjacent transition the AUC of scan
date (Mann-Whitney) and the day overlap (share of the two stages' samples scanned on days when both
stages were scanned). Confounded if AUC >= 0.80 or <= 0.20, or day overlap < 20%.
"""
import os
import sys
from datetime import date, timedelta

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import tda_crc_pilot as T  # noqa: E402

OUT = os.path.join(ROOT, "results", "addendum5", "part2_scandates")
os.makedirs(OUT, exist_ok=True)


def parse(v):
    v = v.strip()
    if not v:
        return pd.NaT
    if v.isdigit():
        return pd.Timestamp(date(1899, 12, 30) + timedelta(days=int(v)))
    try:
        return pd.Timestamp(pd.to_datetime(v, format="%d/%m/%Y"))
    except ValueError:
        return pd.NaT


_, meta = T.read_series_matrix(os.path.join(ROOT, "data", "raw", "GSE41258_series_matrix.txt.gz"))
gsm = meta["!Sample_geo_accession"][0]
raw = pd.Series("", index=gsm)
for row in meta["!Sample_characteristics_ch1"]:
    for i, v in enumerate(row):
        if v.lower().startswith("scan date"):
            raw.iloc[i] = v.split(":", 1)[1].strip()
fmt = raw.map(lambda v: "excel serial" if v.isdigit() else ("dd/mm/yyyy" if "/" in v else ("missing" if not v else "other")))


def parse_corrected(v):
    """Deviation D9: integer = Excel serial produced by a day/month -> month/day locale error; swap back."""
    x = parse(v)
    if v.strip().isdigit() and not pd.isna(x):
        return pd.Timestamp(year=x.year, month=x.day, day=x.month)
    return x


MODE = sys.argv[1] if len(sys.argv) > 1 else "as_specified"
dates = raw.map(parse_corrected if MODE == "corrected" else parse)
SUF = "" if MODE == "as_specified" else "_corrected"

L = ["ADDENDUM 5, PART 2 - GSE41258 scan date vs stage (batch check, descriptive)", "=" * 64,
     "Pre-specified in docs/analysis_prespecification_addendum5_checks.md (commit d5673a7).",
     f"PARSING: {MODE.upper()}" + ("  (deviation D9: Excel serial dates with day and month swapped back)"
                                  if MODE == "corrected" else "  (integer = Excel serial, as in addendum 5; "
                                  "found erroneous, see D9)"),
     f"Check: text dates with day > 12: {int(raw[fmt == 'dd/mm/yyyy'].str.split('/').str[0].astype(int).gt(12).sum())}"
     f"/{int((fmt == 'dd/mm/yyyy').sum())}; integer dates falling on day <= 12 as read: "
     f"{int(raw[fmt == 'excel serial'].map(parse).dt.day.le(12).sum())}/{int((fmt == 'excel serial').sum())}",
     "Confounded if AUC >= 0.80 or <= 0.20, or day overlap < 20%; Cramer's V >= 0.5 = strong association.", ""]
rows = []
for set_name in ("main", "main_allarrays"):
    lab = pd.read_csv(os.path.join(ROOT, "data", "processed", f"GSE41258_{set_name}_samples_qc.csv"), index_col=0)
    d = pd.DataFrame({"stage": lab["stage"], "raw": raw.reindex(lab.index), "format": fmt.reindex(lab.index),
                      "date": dates.reindex(lab.index)})
    d["year"] = d.date.dt.year.astype("Int64")
    d.to_csv(os.path.join(OUT, f"scan_dates_{set_name}{SUF}.csv"), index_label="sample")
    n_missing = int(d.date.isna().sum())
    ok = d.dropna(subset=["date"])
    yt = pd.crosstab(ok.year, ok.stage).reindex(columns=T.STAGES).fillna(0).astype(int)
    chi2, p_chi, dof, _ = chi2_contingency(yt.values)
    V = float(np.sqrt(chi2 / (yt.values.sum() * (min(yt.shape) - 1))))
    L += [f"=== set {set_name} (" + ", ".join(f"{s} {int((lab.stage == s).sum())}" for s in T.STAGES) + ") ===",
          f"Date formats: {d.format.value_counts().to_dict()}; unparseable/missing: {n_missing}",
          "Scan date range per stage (median [min - max]):"]
    for s in T.STAGES:
        x = ok.date[ok.stage == s]
        L.append(f"  {s:<11} n={len(x):<4} median {x.median():%Y-%m-%d}  [{x.min():%Y-%m-%d} - {x.max():%Y-%m-%d}]  "
                 f"distinct scan days {x.dt.normalize().nunique()}")
    L += ["Scan year x stage:", "  " + yt.to_string().replace("\n", "\n  "),
          f"Cramer's V (year x stage) = {V:.2f} (chi2 = {chi2:.1f}, df = {dof}, p = {p_chi:.1e}) -> "
          + ("strong" if V >= 0.5 else "moderate" if V >= 0.2 else "weak") + " association", "",
          "Adjacent transitions:"]
    for s, t in zip(T.STAGES[:-1], T.STAGES[1:]):
        a, b = ok[ok.stage == s], ok[ok.stage == t]
        xa, xb = a.date.map(pd.Timestamp.toordinal).values, b.date.map(pd.Timestamp.toordinal).values
        U = mannwhitneyu(xb, xa, alternative="two-sided").statistic
        auc = float(U / (len(xa) * len(xb)))            # P(later stage scanned later)
        both = set(a.date.dt.normalize()) & set(b.date.dt.normalize())
        ov = float((a.date.dt.normalize().isin(both).sum() + b.date.dt.normalize().isin(both).sum())
                   / (len(a) + len(b)))
        conf = auc >= 0.80 or auc <= 0.20 or ov < 0.20
        rows.append(dict(set=set_name, transition=f"{s}->{t}", n_from=len(a), n_to=len(b), auc=auc,
                         day_overlap=ov, shared_days=len(both), confounded=conf, cramers_v=V))
        L.append(f"  {s:>10} -> {t:<10} AUC = {auc:.2f}  day overlap = {ov * 100:5.1f}% ({len(both)} shared scan days)"
                 f"  -> {'CONFOUNDED' if conf else 'not confounded'}")
    L.append("")
pd.DataFrame(rows).to_csv(os.path.join(OUT, f"scan_date_confounding{SUF}.csv"), index=False)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
d = pd.read_csv(os.path.join(OUT, f"scan_dates_main{SUF}.csv"), index_col=0, parse_dates=["date"])
fig, ax = plt.subplots(figsize=(10, 3.8))
rng = np.random.default_rng(42)
for i, s in enumerate(T.STAGES):
    x = d.date[d.stage == s].dropna()
    ax.scatter(x, i + rng.uniform(-0.25, 0.25, len(x)), s=10, color=T.COLORS[s], alpha=0.7)
ax.set_yticks(range(4)); ax.set_yticklabels([f"{s} (n={(d.stage == s).sum()})" for s in T.STAGES])
ax.set_xlabel("array scan date"); ax.set_title(f"GSE41258 set main: scan date by stage ({MODE.replace(chr(95), chr(32))})")
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"scan_dates_by_stage{SUF}.png"), dpi=150)
L += ["Outputs: results/addendum5/part2_scandates/scan_dates_<set>.csv, scan_date_confounding.csv,",
      "  scan_dates_by_stage.png"]
open(os.path.join(OUT, f"summary{SUF}.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
