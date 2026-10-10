"""
Addendum 10 (296eaed; D18 438d686) - TCGA-COAD step 1: data preparation.

Inputs (data/raw/TCGA/, never edited):
  TCGA.COAD.sampleMap_HiSeqV2.gz              - gene-level log2(RSEM normalised count + 1), genes x samples (Xena)
  TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt - clinical (pathologic_stage, sample_type, _PATIENT, MSI CDE_ID_3226963)
  COAD_mc3_gene_level.txt                     - mc3 gene-level mutation calls (0/1), genes x samples
Rules: primary tumours (code -01), one per patient (first in sorted order); stage = pathologic_stage collapsed to
I-IV (sub-stages merged), unknown excluded; solid-tissue normals (-11, one per patient) kept for the secondary set.
MSI: clinical CDE_ID_3226963 (MSI-H / MSI-L / MSS); if missing, mc3 hypermutation (>= 200 mutated genes = MSI-H-like).
Genes: as given (already log2); genes with missing values removed. The lowest-25% expression filter is applied in
step 3 (D18b), as for the GEO cohorts.

Outputs: data/processed/TCGA_COAD_{main,withnormal}_expr_genes.csv and _samples.csv; results/step01/TCGA/.
"""
import hashlib
import os
import re
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "TCGA")
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "step01", "TCGA")
os.makedirs(OUT, exist_ok=True)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stage_of(s):
    if pd.isna(s):
        return None
    m = re.match(r"\s*Stage\s+(IV|III|II|I)[ABC]?\s*$", str(s), flags=re.I)
    return {"I": "Stage I", "II": "Stage II", "III": "Stage III", "IV": "Stage IV"}[m.group(1).upper()] if m else None


def main():
    files = {k: os.path.join(RAW, v) for k, v in (("expr", "TCGA.COAD.sampleMap_HiSeqV2.gz"),
                                                    ("clin", "TCGA.COAD.sampleMap_COAD_clinicalMatrix.txt"),
                                                    ("mut", "COAD_mc3_gene_level.txt"))}
    expr = pd.read_csv(files["expr"], sep="\t", index_col=0)
    clin = pd.read_csv(files["clin"], sep="\t", index_col=0)
    mut = pd.read_csv(files["mut"], sep="\t", index_col=0)
    n_genes_raw = expr.shape[0]
    p99 = float(np.nanpercentile(expr.values, 99))
    if p99 > 50 or expr.values.min() < -0.01:
        raise SystemExit(f"Expression not on log2(x+1) scale as expected (p99 {p99:.1f}, min {expr.values.min():.2f})")
    expr = expr[~expr.index.duplicated()].dropna(axis=0, how="any")

    # sample table
    rows = []
    for s in expr.columns:
        code = s.split("-")[-1][:2] if s.count("-") >= 3 else ""
        c = clin.loc[s] if s in clin.index else None
        rows.append(dict(sample=s, patient_id=s[:12], code=code,
                         sample_type=None if c is None else c.get("sample_type"),
                         pathologic_stage=None if c is None else c.get("pathologic_stage"),
                         msi_clinical=None if c is None else c.get("CDE_ID_3226963")))
    tab = pd.DataFrame(rows).set_index("sample")
    tab["stage"] = tab.pathologic_stage.map(stage_of)
    # mutations: number of mutated genes per sample (mc3 columns are sample barcodes)
    nmut = mut.sum(axis=0)
    nmut.index = [c[:15] for c in nmut.index]
    tab["n_mutated_genes"] = [nmut.get(s[:15], np.nan) for s in tab.index]
    msi = tab.msi_clinical.astype(str).str.upper()
    tab["msi_source"] = np.where(msi.isin(["MSI-H", "MSI-L", "MSS"]), "clinical",
                                 np.where(tab.n_mutated_genes.notna(), "mc3 hypermutation", "unavailable"))
    tab["msi_class"] = np.where(msi == "MSI-H", "MSI-H",
                       np.where(msi.isin(["MSI-L", "MSS"]), "MSS/MSI-L",
                       np.where(tab.n_mutated_genes >= 200, "MSI-H-like (hypermutated)",
                       np.where(tab.n_mutated_genes.notna(), "MSS-like", "unknown"))))

    # selection with reasons
    tab["included_main"] = False
    tab["included_normal"] = False
    tab["reason"] = ""
    tum = tab[tab.code == "01"].sort_index()
    first_t = tum[~tum.patient_id.duplicated()].index
    tab.loc[tum.index.difference(first_t), "reason"] = "duplicate primary tumour of patient"
    nost = [s for s in first_t if tab.loc[s, "stage"] is None]
    tab.loc[nost, "reason"] = "stage missing / not I-IV"
    keep = [s for s in first_t if tab.loc[s, "stage"] is not None]
    tab.loc[keep, "included_main"] = True
    nor = tab[tab.code == "11"].sort_index()
    first_n = nor[~nor.patient_id.duplicated()].index
    tab.loc[nor.index.difference(first_n), "reason"] = "duplicate normal of patient"
    tab.loc[first_n, "included_normal"] = True
    tab.loc[~tab.code.isin(["01", "11"]), "reason"] = "other sample type (code " + tab.code + ")"
    tab.loc[first_n, "stage"] = "Normal"
    tab.to_csv(os.path.join(OUT, "all_samples_labels.csv"), index_label="sample")

    main_ids = list(tab.index[tab.included_main])
    wn_ids = list(tab.index[tab.included_normal]) + main_ids
    cols = ["patient_id", "stage", "pathologic_stage", "msi_class", "msi_source", "n_mutated_genes", "sample_type"]
    for name, ids in (("main", main_ids), ("withnormal", wn_ids)):
        expr[ids].to_csv(os.path.join(PROC, f"TCGA_COAD_{name}_expr_genes.csv"), float_format="%.5f")
        tab.loc[ids, cols].to_csv(os.path.join(PROC, f"TCGA_COAD_{name}_samples.csv"), index_label="sample")

    # QC figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"]
    colr = dict(zip(order, ["#2a78d6", "#7fb800", "#e0a100", "#e05a00", "#c0392b"]))
    ids = sorted(wn_ids, key=lambda s: (order.index(tab.loc[s, "stage"]), s))
    q = expr[ids].quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    fig, ax = plt.subplots(figsize=(14, 3.6))
    xs = np.arange(len(ids))
    cs = [colr[tab.loc[s, "stage"]] for s in ids]
    ax.vlines(xs, q.loc[0.05], q.loc[0.95], color=cs, lw=0.4)
    ax.vlines(xs, q.loc[0.25], q.loc[0.75], color=cs, lw=1.5)
    ax.scatter(xs, q.loc[0.5], s=2, c="k")
    ax.set_ylabel("log2(RSEM + 1)"); ax.set_xticks([])
    ax.set_title("TCGA-COAD per-sample expression distribution (5-95%, IQR, median), ordered by stage")
    for k in order:
        ax.plot([], [], color=colr[k], lw=3, label=k)
    ax.legend(fontsize=7, ncol=5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "qc_sample_distributions.png"), dpi=150); plt.close(fig)

    m = tab.loc[main_ids]
    med = expr[wn_ids].median()
    L = ["STEP 1 - TCGA-COAD data preparation (addendum 10, 296eaed; D18 438d686)", "=" * 66,
         "Inputs (data/raw/TCGA/, SHA-256):"] + [f"  {os.path.basename(p)}  {sha(p)}" for p in files.values()] + [
         f"Expression: {n_genes_raw} genes x {expr.shape[1]} samples; log2(RSEM+1) check: 99th pct {p99:.2f} (no log2 applied)",
         f"Genes after removing duplicates / missing values: {expr.shape[0]} (lowest-25% filter applied in step 3, D18b)",
         f"Samples by code: " + ", ".join(f"-{k} {v}" for k, v in tab.code.value_counts().items()),
         "", "Exclusions:"] + [f"  {k}: {v}" for k, v in tab.reason[tab.reason != ""].value_counts().items()] + [
         "", f"PRIMARY SET 'main' (primary tumours, one per patient, stage I-IV): {len(main_ids)}"]
    L += [f"  {k}: {v}" for k, v in m.stage.value_counts().reindex(order[1:]).items()]
    L += [f"  smallest stage {m.stage.value_counts().min()} -> step 5 subsample m = {int(0.8 * m.stage.value_counts().min())}",
          "", "MSI class by stage (main):"]
    ct = pd.crosstab(m.msi_class, m.stage).reindex(columns=order[1:], fill_value=0)
    L += ["  " + line for line in ct.to_string().splitlines()]
    L += [f"  MSI source: " + ", ".join(f"{k} {v}" for k, v in m.msi_source.value_counts().items()),
          "", f"SECONDARY SET 'withnormal': solid-tissue normals (one per patient) {int(tab.included_normal.sum())} + main",
          f"  normals from patients also in main: {len(set(tab.loc[tab.included_normal, 'patient_id']) & set(m.patient_id))}",
          "", f"Per-sample median log2 expression: IQR {med.quantile(.25):.2f}-{med.quantile(.75):.2f}, "
              f"range {med.min():.2f}-{med.max():.2f}",
          "Packages: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("pandas", "numpy")),
          "Outputs: data/processed/TCGA_COAD_{main,withnormal}_expr_genes.csv / _samples.csv; results/step01/TCGA/"]
    open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
