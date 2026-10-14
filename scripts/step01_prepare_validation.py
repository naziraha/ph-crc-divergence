"""
Step 1 (validation cohorts) - data preparation for GSE20916 / GSE37364 (GPL570).

Frozen rules (docs/analysis_prespecification_addendum_validation.md, commit 6ea31ea):
probes -> gene symbols with the GPL570 annotation (no-symbol and multi-symbol probes
dropped); log2 if 99th percentile > 50; highest-mean probe per gene. Each analysis
set is processed separately (probe filtering and scaling differ between GSE20916
sub-studies). Series matrices are read with the validated reference reader.

Usage: python scripts/step01_prepare_validation.py GSE20916
Outputs: data/processed/<GSE>_<set>_expr_genes.csv, <GSE>_<set>_samples.csv,
         results/step01/<GSE>/*
"""
import os
import sys
from importlib import metadata

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import tda_crc_pilot as T  # noqa: E402
import GEOparse  # noqa: E402

SEED = 42
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
STAGES3 = ["Normal", "Adenoma", "Carcinoma"]


def ch_field(meta, samples, name):
    out = pd.Series("", index=samples)
    for row in meta.get("!Sample_characteristics_ch1", []):
        for i, v in enumerate(row):
            if v.lower().startswith(name):
                out.iloc[i] = v.split(":", 1)[1].strip()
    return out


def labels_gse20916(meta, samples):
    """Microdissected sub-study (title tag '(micro)'): stage from 'malignacy' field;
    bulk sub-study: stage from 'tissue' field. All mappings recorded."""
    title = pd.Series(meta["!Sample_title"][0], index=samples)
    tissue = ch_field(meta, samples, "tissue:")
    malig = ch_field(meta, samples, "malignacy")
    morph = ch_field(meta, samples, "tissue morphology")
    micro = title.str.contains(r"\(micro\)")
    stage, sub, label = pd.Series(None, index=samples, dtype=object), pd.Series("", index=samples), \
        pd.Series("", index=samples)
    map_micro = {"normal colon": "Normal", "distant normal colon": "Normal",
                 "adenoma": "Adenoma", "carcinoma": "Carcinoma"}
    map_bulk = {"normal colon": "Normal", "adenoma": "Adenoma", "adenocarcinoma": "Carcinoma"}
    for s in samples:
        if micro[s]:
            stage[s] = map_micro.get(malig[s])
            sub[s] = "micro"
            label[s] = f"malignacy: {malig[s]}; morphology: {morph[s]}"
        else:
            stage[s] = map_bulk.get(tissue[s])
            sub[s] = "bulk"
            label[s] = f"tissue: {tissue[s]}"
    extra = pd.DataFrame({
        "normal_type": np.where(malig == "distant normal colon", "distant",
                                np.where((malig == "normal colon") | (micro & (malig == "normal colon")),
                                         "adjacent to tumour", np.where(stage == "Normal", "not stated", ""))),
        "morphology": morph.replace("", np.nan),
        "protocol": ch_field(meta, samples, "protocol").replace("", np.nan),
        "title": title})
    sets = {"micro": samples_where(sub, "micro"), "bulk": samples_where(sub, "bulk")}
    return stage, label, sets, extra


def samples_where(series, value):
    return list(series.index[series == value])


def labels_gse37364(meta, samples):
    """One bulk set. Stage from 'disease state'; the 13 Dukes C/D carcinomas carry only
    'Stage: Dukes C/D' in that field and are mapped to Carcinoma. Grade/Dukes recorded."""
    title = pd.Series(meta["!Sample_title"][0], index=samples)
    c0 = pd.Series(meta["!Sample_characteristics_ch1"][0], index=samples)
    c1 = pd.Series(meta["!Sample_characteristics_ch1"][1], index=samples)
    mapping = {"disease state: Normal colonic mucosa": "Normal", "disease state: Adenoma": "Adenoma",
               "disease state: Colorectal adenocarcinoma": "Carcinoma", "Stage: Dukes C/D": "Carcinoma"}
    stage = c0.map(mapping)
    grade = np.where(c0 == "Stage: Dukes C/D", "Dukes C/D", c1.str.replace("Stage: ", "", regex=False))
    label = c0 + np.where(c1 != "", "; " + c1, "")
    extra = pd.DataFrame({"normal_type": np.where(stage == "Normal", "patient without neoplasia", ""),
                          "morphology": np.nan, "protocol": "biopsy", "title": title,
                          "grade_or_dukes": grade}, index=samples)
    return stage, label, {"bulk": samples}, extra


COHORTS = {"GSE20916": labels_gse20916, "GSE37364": labels_gse37364}


def collapse_to_genes(expr, probe2gene):
    x = expr.loc[expr.index.intersection(probe2gene.index)]
    mean = x.mean(axis=1)
    pick = (pd.DataFrame({"gene": probe2gene.loc[x.index], "mean": mean})
            .sort_values("mean", ascending=False).groupby("gene", sort=True).head(1))
    genes = x.loc[pick.index]
    genes.index = pick["gene"].values
    return genes.sort_index(), pick.rename_axis("probe").reset_index()


def main():
    gse = sys.argv[1]
    out = os.path.join(ROOT, "results", "step01", gse)
    os.makedirs(out, exist_ok=True)
    expr, meta = T.read_series_matrix(os.path.join(RAW, f"{gse}_series_matrix.txt.gz"))
    samples = list(expr.columns)
    stage, label, sets, extra = COHORTS[gse](meta, samples)

    gpl = GEOparse.get_GEO(geo="GPL570", destdir=RAW, silent=True)
    ann = gpl.table[["ID", "Gene Symbol"]].copy()
    ann["ID"] = ann["ID"].astype(str)
    ann["Gene Symbol"] = ann["Gene Symbol"].fillna("").astype(str).str.strip()
    no_sym = ann["Gene Symbol"] == ""
    multi = ann["Gene Symbol"].str.contains("///")
    probe2gene = ann.loc[~no_sym & ~multi].set_index("ID")["Gene Symbol"]

    all_tab = pd.DataFrame({"raw_label": label, "stage": stage}).join(extra)
    all_tab["set"] = ""
    for k, v in sets.items():
        all_tab.loc[v, "set"] = k
    all_tab.to_csv(os.path.join(out, "all_samples_labels.csv"), index_label="sample")

    lines = [f"STEP 1 - Data preparation ({gse}, GPL570)", "=" * 50,
             f"Input: data/raw/{gse}_series_matrix.txt.gz; GPL570 annotation (data/raw/GPL570.txt)",
             "Pre-specification: docs/analysis_prespecification_addendum_validation.md (commit 6ea31ea)",
             f"Seed: {SEED}",
             "Versions: " + ", ".join(f"{p} {metadata.version(p)}" for p in ("GEOparse", "pandas", "numpy"))
             + f", python {sys.version.split()[0]}", "",
             f"Samples in series matrix: {len(samples)}",
             "Label -> stage mapping (all samples):",
             "  " + pd.crosstab([all_tab["set"], all_tab["raw_label"]], all_tab["stage"].fillna("excluded"))
             .to_string().replace("\n", "\n  "), "",
             f"GPL570 probes: {len(ann)}; without symbol {int(no_sym.sum())}; multi-symbol dropped "
             f"{int((multi & ~no_sym).sum())}; single-gene probes {len(probe2gene)}",
             f"Probes in series matrix: {expr.shape[0]}", ""]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(sets), 1, figsize=(12, 3.2 * len(sets)))
    axes = np.atleast_1d(axes)
    for ax, (name, ids) in zip(axes, sets.items()):
        ids = [s for s in ids if stage[s] in STAGES3]
        e = expr[ids]
        n_na = int(e.isna().any(axis=1).sum())
        e = e.dropna(axis=0, how="any")
        p99 = float(np.nanpercentile(e.values, 99))
        log2_applied = p99 > 50
        if log2_applied:
            e = np.log2(e.clip(lower=1))
        genes, pick = collapse_to_genes(e, probe2gene)
        genes.index.name = "gene"
        keep_cols = ["stage", "raw_label", "normal_type", "morphology", "protocol", "title"] + \
            [c for c in ("grade_or_dukes",) if c in all_tab.columns]
        tab = all_tab.loc[ids, keep_cols]
        tab["patient_id"] = ""            # no patient identifier in GEO metadata
        genes.to_csv(os.path.join(PROC, f"{gse}_{name}_expr_genes.csv"), float_format="%.5f")
        tab.to_csv(os.path.join(PROC, f"{gse}_{name}_samples.csv"), index_label="sample")
        pick.to_csv(os.path.join(out, f"probe_to_gene_{name}.csv"), index=False)
        med = genes.median(axis=0)
        cnt = tab["stage"].value_counts().reindex(STAGES3).fillna(0).astype(int)
        lines += [f"--- set {name} ---",
                  "Samples per stage: " + ", ".join(f"{s} {cnt[s]}" for s in STAGES3) + f" (total {cnt.sum()})"]
        if name == "micro":
            lines.append("  by morphology x stage:\n    " + pd.crosstab(tab["morphology"], tab["stage"])
                         .reindex(columns=STAGES3).to_string().replace("\n", "\n    "))
            lines.append("  normal type: " + ", ".join(f"{k} {v}" for k, v in
                                                       tab.loc[tab.stage == "Normal", "normal_type"]
                                                       .value_counts().items()))
        small = [s for s in STAGES3 if cnt[s] < 13]
        lines += [f"  stages with < 13 samples (descriptive only per addendum): {', '.join(small) or 'none'}",
                  f"  probes with missing values in this set (dropped): {n_na}",
                  f"  99th percentile = {p99:.2f}; log2 applied: {log2_applied}",
                  f"  genes after highest-mean probe per gene: {genes.shape[0]}",
                  f"  per-sample median log2: IQR {med.quantile(.25):.2f}-{med.quantile(.75):.2f}, "
                  f"range {med.min():.2f}-{med.max():.2f}",
                  "  patient ids: not available in GEO metadata (de-duplication / one-per-patient not possible)",
                  ""]
        order = tab.sort_values("stage", key=lambda x: x.map(STAGES3.index)).index
        q = genes[order].quantile([0.05, 0.25, 0.5, 0.75, 0.95])
        xs = np.arange(len(order))
        cols = [T.COLORS[s] for s in tab.loc[order, "stage"]]
        ax.vlines(xs, q.loc[0.05], q.loc[0.95], color=cols, lw=0.6, alpha=0.5)
        ax.vlines(xs, q.loc[0.25], q.loc[0.75], color=cols, lw=2)
        ax.scatter(xs, q.loc[0.5], s=3, c="k")
        ax.set_xticks([]); ax.set_ylabel("log2 expression")
        ax.set_title(f"{gse} {name}: " + ", ".join(f"{s} {cnt[s]}" for s in STAGES3))
    fig.tight_layout()
    fig.savefig(os.path.join(out, "qc_sample_distributions.png"), dpi=150)
    excl = all_tab[all_tab["stage"].isna()]
    lines += [f"Samples not mapped to a stage: {len(excl)}"
              + (" -> " + ", ".join(excl["raw_label"].unique()) if len(excl) else ""), "",
              "Outputs:",
              f"  data/processed/{gse}_<set>_expr_genes.csv   genes x samples, log2",
              f"  data/processed/{gse}_<set>_samples.csv      stage, original label, normal type, morphology",
              f"  results/step01/{gse}/all_samples_labels.csv, probe_to_gene_<set>.csv, qc_sample_distributions.png"]
    with open(os.path.join(out, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
