"""
Supplementary Table S_sensitivity: step 5 sensitivity analyses, copied from existing results; no tests.
(A) PC and gene-number settings, GEO cohorts (primary 10 PCs / 2,000 genes; 5 PCs; 20 PCs; 1,000 genes; 5,000 genes);
(B) alternative sample sets and adjustments, GEO cohorts (GSE41258 S1-S4, S9, D3 stage-aware; GSE117606 noGSM3304812,
    noQAflag, paired, quality-adjusted; stage-aware adjustment in GSE20916, GSE37364, GSE117606 for every setting);
(C) TCGA-COAD (PC/gene settings, no-low-colon-identity, MSS-only).
Rule: shape-only Holm p < 0.05 in both unadjusted and adjusted data. Raw-mode Holm p in the .csv only.
Output: paper/final/supplementary/TableS_sensitivity.{csv,docx}

Usage: python scripts/make_table_s_sensitivity.py
"""
import os

import pandas as pd
from docx import Document
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S5 = os.path.join(ROOT, "results", "step05")
OUT = os.path.join(ROOT, "paper", "final", "supplementary")
SET = {"primary": "primary (10 PCs, 2,000 genes)", "S5_pc5": "5 PCs", "S6_pc20": "20 PCs",
       "S7_genes1000": "1,000 genes", "S8_genes5000": "5,000 genes"}
TR_GEO = {"Normal->Adenoma": "N → A", "Adenoma->Carcinoma": "A → C", "Carcinoma->Metastasis": "C → M"}
TR_TCGA = {"Stage I->Stage II": "I → II", "Stage II->Stage III": "II → III", "Stage III->Stage IV": "III → IV"}


def load():
    d = {"GSE41258": pd.read_csv(os.path.join(S5, "permutation_results_all.csv")),
         "GSE20916": pd.read_csv(os.path.join(S5, "GSE20916", "permutation_results_all.csv")),
         "GSE37364": pd.read_csv(os.path.join(S5, "GSE37364", "permutation_results_all.csv")),
         "GSE117606": pd.read_csv(os.path.join(S5, "GSE117606", "secondary_results_all.csv")),
         "TCGA-COAD": pd.read_csv(os.path.join(S5, "TCGA", "permutation_results_all.csv"))}
    d["GSE117606"]["analysis"] = d["GSE117606"]["analysis"].replace({"PRIMARY": "primary"})
    d3 = pd.read_csv(os.path.join(ROOT, "results", "exploratory_D3", "permutation_results_D3.csv"))
    d3["analysis"], d3["version"] = "D3_stageaware", "stageaware"
    d["GSE41258"] = pd.concat([d["GSE41258"], d3], ignore_index=True)
    return d


def get(df, analysis, version, tr, mode, col="p_holm"):
    g = df[(df.analysis == analysis) & (df.version == version) & (df.transition == tr) & (df["mode"] == mode)]
    if len(g) != 1:
        raise ValueError(f"{analysis} {version} {tr} {mode}: {len(g)} rows")
    return g.iloc[0] if col is None else float(g.iloc[0][col])


def row(part, cohort, df, tr, trlab, label, a1, v1, a2, v2, pair, col="p_holm"):
    """a1/v1 = first (unadjusted) analysis/version; a2/v2 = second (adjusted-type) analysis/version."""
    r1, r2 = get(df, a1, v1, tr, "shape-only", None), get(df, a2, v2, tr, "shape-only", None)
    su, sa = float(r1[col]), float(r2[col])
    return {"part": part, "cohort": cohort, "transition": trlab, "setting / sample set": label,
            "versions compared": pair, "n (from / to)": f"{int(r2.n_from)} / {int(r2.n_to)}", "m": int(r2.subsample_size),
            "shape-only Holm p, unadjusted": su, "shape-only Holm p, adjusted": sa,
            "rule met": ("yes" if pair.startswith(UA) else "analogue only") if (su < 0.05 and sa < 0.05) else "no",
            "raw Holm p, unadjusted": get(df, a1, v1, tr, "raw", col), "raw Holm p, adjusted": get(df, a2, v2, tr, "raw", col)}


UA = "unadjusted + stage-blind adjusted"


def part_a(d):
    rows = []
    for c in ("GSE41258", "GSE20916", "GSE37364", "GSE117606"):
        trs = list(TR_GEO) if c == "GSE41258" else list(TR_GEO)[:2]
        for tr in trs:
            for s, lab in SET.items():
                rows.append(row("A", c, d[c], tr, TR_GEO[tr], lab, s, "unadjusted", s, "adjusted", UA))
    return pd.DataFrame(rows)


def part_b(d):
    rows, g41 = [], d["GSE41258"]
    sets41 = {"S1_withlung": "S1: lung metastases added (withlung)",
              "S2_main_allarrays": "S2: all arrays, no de-duplication (some patients counted more than once)",
              "S3_main_inclqc": "S3: tissue-identity exclusions put back",
              "S4_one_per_patient": "S4: each compared pair made patient-disjoint"}
    for s, lab in sets41.items():
        for tr in TR_GEO:
            rows.append(row("B", "GSE41258", g41, tr, TR_GEO[tr], lab, s, "unadjusted", s, "adjusted", UA))
    r = row("B", "GSE41258", g41, "Carcinoma->Metastasis", "C → M", "S9: 3 liver-rich metastases removed (single comparison)",
            "S9_no_high_liver_mets", "unadjusted", "S9_no_high_liver_mets", "adjusted", UA + " (p, no Holm)", col="p")
    rows.append(r)
    for tr in TR_GEO:
        rows.append(row("B", "GSE41258", g41, tr, TR_GEO[tr], "stage-aware adjustment (exploratory, D3)",
                        "primary", "unadjusted", "D3_stageaware", "stageaware", "unadjusted + stage-aware"))
    g117 = d["GSE117606"]
    sets117 = {"noGSM3304812": "without GSM3304812 (hepatocyte signal; D11)",
               "noQAflag": "without arrays flagged by the quality rule",
               "paired": "paired subsets (patients with both stages)"}
    for s, lab in sets117.items():
        for tr in list(TR_GEO)[:2]:
            rows.append(row("B", "GSE117606", g117, tr, TR_GEO[tr], lab, s, "unadjusted", s, "adjusted", UA))
    for tr in list(TR_GEO)[:2]:
        rows.append(row("B", "GSE117606", g117, tr, TR_GEO[tr], "quality-adjusted (composition + array quality; D12)",
                        "primary", "unadjusted", "adjqual", "adjusted+quality", "unadjusted + quality-adjusted"))
    for c in ("GSE20916", "GSE37364", "GSE117606"):
        for s, lab in SET.items():
            a2 = "stageaware" if (c == "GSE117606" and s == "primary") else s
            for tr in list(TR_GEO)[:2]:
                rows.append(row("B", c, d[c], tr, TR_GEO[tr], f"stage-aware adjustment, {lab}",
                                s, "unadjusted", a2, "stageaware", "unadjusted + stage-aware"))
    return pd.DataFrame(rows)


def part_c(d):
    rows, t = [], d["TCGA-COAD"]
    sets = dict(SET, noLowColon="without 18 low-colon-identity tumours (D19)", MSS="MSS / non-hypermutated only")
    for s, lab in sets.items():
        for tr in TR_TCGA:
            rows.append(row("C", "TCGA-COAD", t, tr, TR_TCGA[tr], lab, s, "unadjusted", s, "adjusted", UA))
    return pd.DataFrame(rows)


def add(doc, df, cols, title):
    doc.add_paragraph().add_run(title).bold = True
    t = doc.add_table(rows=1, cols=len(cols))
    t.style = "Table Grid"
    for j, c in enumerate(cols):
        r = t.rows[0].cells[j].paragraphs[0].add_run(c)
        r.bold, r.font.size = True, Pt(7)
    for _, rw in df.iterrows():
        cells = t.add_row().cells
        for j, c in enumerate(cols):
            v = rw[c]
            run = cells[j].paragraphs[0].add_run(f"{v:.3f}" if isinstance(v, float) else str(v))
            run.font.size = Pt(7)
            run.bold = rw["rule met"] == "yes"


def main():
    os.makedirs(OUT, exist_ok=True)
    d = load()
    A, B, C = part_a(d), part_b(d), part_c(d)
    allr = pd.concat([A, B, C], ignore_index=True)
    allr.to_csv(os.path.join(OUT, "TableS_sensitivity.csv"), index=False)

    nonprim = A[A["setting / sample set"] != SET["primary"]]
    met = nonprim[nonprim["rule met"] == "yes"]
    n_met, n_tot = len(met), len(nonprim)
    n_prim_met = int((A[A["setting / sample set"] == SET["primary"]]["rule met"] == "yes").sum())
    # replication check: same transition and setting met in more than one cohort?
    assert met.groupby(["transition", "setting / sample set"]).cohort.nunique().max() == 1
    summ = (f"Rule met in {n_met} of {n_tot} PC/gene combinations (Part A, excluding the primary rows; "
            f"{n_prim_met} of 9 primary rows met), each in a single cohort: GSE41258 N → A with 5 PCs; GSE20916 N → A "
            "with 1,000 genes; GSE117606 A → C with 5 PCs and with 5,000 genes. None was replicated in another cohort "
            "with the same setting, and none was met in the primary setting.")

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    sec = doc.sections[0]
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(1.3))
    doc.add_paragraph().add_run("Supplementary Table S_sensitivity. Sensitivity analyses of the step 5 test "
                                "(shape-only Holm p; rows where the tipping-point rule is met in bold).").bold = True
    cols = ["cohort", "transition", "setting / sample set", "n (from / to)", "m",
            "shape-only Holm p, unadjusted", "shape-only Holm p, adjusted", "rule met"]
    add(doc, A, cols, "A. PC and gene-number settings (GEO cohorts; unadjusted and stage-blind adjusted data)")
    add(doc, B, ["cohort", "transition", "setting / sample set", "versions compared", "n (from / to)", "m",
                 "shape-only Holm p, unadjusted", "shape-only Holm p, adjusted", "rule met"],
        "B. Alternative sample sets and adjustments (GEO cohorts)")
    add(doc, C, cols, "C. TCGA-COAD (stage I–IV; unadjusted and stage-blind adjusted data)")
    notes = [
        "Primary rule (prespecified before analysis): a transition is a topological tipping point only if the shape-only "
        "label-permutation test of the difference in persistence-diagram shape between stages (mean 1-Wasserstein "
        "distance between diagrams of different stages minus the mean within stages; H0 and H1 combined) (one sample per patient per stage where patient identifiers exist; k = 20 subsamples of m "
        "per group, 80% of the smallest stage; 999 permutations) gives Holm p < 0.05 in both unadjusted and stage-blind "
        "composition-adjusted data. Holm correction is across the transitions of each cohort (3 for GSE41258 and TCGA-COAD, "
        "2 for the other cohorts), within each version and setting. The rule was not met in any primary analysis.",
        "Sensitivity results are descriptive. They do not replace the primary analysis and cannot create a tipping point.",
        summ,
        "Part B, 'versions compared': for the stage-aware and quality-adjusted rows, the 'adjusted' column holds the "
        "stage-aware or quality-adjusted p-value and the 'unadjusted' column the unadjusted p-value of the same setting, "
        "so 'analogue only' there means an analogue of the rule is met, not the pre-specified rule; these rows are not "
        "bold (3 rows: GSE20916 N → A with 5 PCs and with 5,000 genes; GSE117606 A → C with 5,000 genes). The stage-aware adjustment uses the "
        "stage labels and can under-correct when composition and stage are confounded; a result significant only with it "
        "is not reported as a tipping point (addendum 2). GSE41258 stage-aware (D3) is exploratory, not pre-specified.",
        "GSE41258 S2 (all arrays, no de-duplication) counts some patients more than once within a stage, so samples within "
        "a group are not independent; the primary set keeps one array per patient per stage. S4 additionally makes each "
        "compared pair of stages patient-disjoint (subsample size per D2). S9 is a single comparison (p, not Holm-adjusted).",
        "GSE117606: many patients contribute to two stages; the unpaired test is conservative (addendum 6). The paired "
        "rows use only patients with both stages. GSE20916 microdissected sets had no step 5 test (too few samples; "
        "descriptive only, D4) and are not shown.",
        "Raw-mode (spread and shape) Holm p-values for every row are in the .csv file only. No new tests were run; all "
        "values are copied from existing result files.",
    ]
    for n in notes:
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(os.path.join(OUT, "TableS_sensitivity.docx"))

    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    show = ["part", "cohort", "transition", "setting / sample set", "versions compared",
            "shape-only Holm p, unadjusted", "shape-only Holm p, adjusted"]
    print(f"rows: A {len(A)}, B {len(B)}, C {len(C)}")
    print(summ)
    print("RULE MET ROWS:")
    print(allr[allr["rule met"] != "no"][show + ["rule met"]].to_string(index=False))


if __name__ == "__main__":
    main()
