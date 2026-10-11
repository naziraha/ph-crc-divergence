"""
Main Tables 1 (cohorts) and 2 (primary validation results) for the GEO cohorts. All numbers are copied from existing
files; no statistical tests are run. Outputs (paper/final/main_tables/): Table1_cohorts.csv/.docx,
Table2_validation.csv/.docx (editable Word tables; landscape page for Table 2).

Usage: python scripts/make_main_tables.py
"""
import os

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "paper", "final", "main_tables")
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
SETS = {"GSE41258": "main", "GSE20916": "bulk", "GSE37364": "bulk", "GSE117606": "main"}
PLATFORM = {"GSE41258": "GPL96 (Affymetrix Human Genome U133A Array)",      # !Platform_title in data/raw/GPL*.txt
            "GSE20916": "GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)",
            "GSE37364": "GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)",
            "GSE117606": "GPL25373 (Affymetrix HT HG-U133+ PM Array Plate)"}
PATIENT_COL = {"GSE41258": "patient_id", "GSE20916": "patient_id", "GSE37364": "patient_id", "GSE117606": "patient id"}
# Source publications: PubMed IDs from each series matrix (!Series_pubmed_id); citation details verified in PubMed
# (title and sample numbers match the GEO series). GSE117606 lists no PubMed ID.
PUB = {"GSE41258": "Sheffer M et al. Proc Natl Acad Sci U S A 2009;106:7131-6 (PMID 19359472; doi:10.1073/pnas.0902232106)",
       "GSE20916": "Skrzypczak M et al. PLoS One 2010;5:e13091 (PMID 20957034; doi:10.1371/journal.pone.0013091)",
       "GSE37364": "Galamb O et al. PLoS One 2012;7:e48547 (PMID 23155391; doi:10.1371/journal.pone.0048547)",
       "GSE117606": "to verify (no PubMed ID in the GEO record)"}
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
TRL = {"Normal->Adenoma": "N→A", "Adenoma->Carcinoma": "A→C", "Carcinoma->Metastasis": "C→M"}


# ---------------------------------------------------------------- Table 1
def table1():
    rows = []
    for c in GEO:
        st = SETS[c]
        emb = pd.read_csv(os.path.join(R, "step03", "embeddings", f"{c}_{st}_adjusted_g2000_pc10.csv"), index_col=0)
        lab = pd.read_csv(os.path.join(PROC, f"{c}_{st}_samples_qc.csv")).set_index("sample").loc[emb.index]
        pid = lab[PATIENT_COL[c]]
        n_pat = str(pid.nunique()) if pid.notna().all() else "not stated in GEO"
        cnt = emb.stage.value_counts()
        r = {"Cohort": c, "Platform": PLATFORM[c], "Sample set": st, "Patients": n_pat}
        for s in STAGES:
            r[s] = str(int(cnt[s])) if s in cnt else "–"
        r["Total samples"] = str(len(emb))
        r["Source publication"] = PUB[c]
        rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- Table 2
def dim_call(dims, c, tr):
    sig = []
    for (v, mo), g in dims[(dims.cohort == c) & (dims.transition == tr)].groupby(["version", "mode"]):
        p = g.set_index("dim").p_holm
        if p["H0+H1"] < 0.05:
            sig.append((p["H0"] < 0.05, p["H1"] < 0.05))
    if not sig:
        return "–", "no significant combined result"
    n, h0, h1 = len(sig), sum(a for a, _ in sig), sum(b for _, b in sig)
    main = "H0" if h0 > h1 else "H1" if h1 > h0 else ("H0 and H1" if h0 > 0 else "neither alone")
    return main, f"H0 {h0}/{n}, H1 {h1}/{n}"


def table2():
    v = pd.read_csv(os.path.join(R, "source_data", "validation_summary_primary.csv"))
    dims = pd.read_csv(os.path.join(R, "addendum12", "dimension_split_all_cohorts.csv"))
    rows = []
    for _, r in v.iterrows():
        main, detail = dim_call(dims, r.cohort, r.transition)
        rows.append({"Cohort": r.cohort, "Transition": TRL[r.transition], "n (from / to)": f"{r.n_from} / {r.n_to}",
                     "Raw, unadjusted": r.raw_p_holm_unadj, "Raw, adjusted": r.raw_p_holm_adjus,
                     "Shape-only, unadjusted": r.shape_p_holm_unadj, "Shape-only, adjusted": r.shape_p_holm_adjus,
                     "Tipping-point rule met": r.rule_met,
                     "Main dimension": main, "Dimension detail": detail,
                     "H1 landscape raw, unadjusted": r.H1landscape_raw_p_holm_unadj,
                     "H1 landscape raw, adjusted": r.H1landscape_raw_p_holm_adjus,
                     "H1 landscape shape-only, unadjusted": r.H1landscape_shape_p_holm_unadj,
                     "H1 landscape shape-only, adjusted": r.H1landscape_shape_p_holm_adjus})
    return pd.DataFrame(rows)


PCOLS = ["Raw, unadjusted", "Raw, adjusted", "Shape-only, unadjusted", "Shape-only, adjusted",
         "H1 landscape raw, unadjusted", "H1 landscape raw, adjusted", "H1 landscape shape-only, unadjusted",
         "H1 landscape shape-only, adjusted"]


# ---------------------------------------------------------------- Word output
def word_table(df, path, title, notes, landscape=False, pcols=(), fs=8):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(10)
    sec = doc.sections[0]
    if landscape:
        sec.orientation = WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(1.5))
    p = doc.add_paragraph()
    p.add_run(title).bold = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    for j, col in enumerate(df.columns):
        cell = t.rows[0].cells[j].paragraphs[0]
        run = cell.add_run(col)
        run.bold, run.font.size = True, Pt(fs)
    for _, r in df.iterrows():
        cells = t.add_row().cells
        for j, col in enumerate(df.columns):
            val = r[col]
            par = cells[j].paragraphs[0]
            if col in pcols:
                run = par.add_run(f"{val:.3f}")
                run.bold = bool(val < 0.05)
            else:
                run = par.add_run(str(val))
            run.font.size = Pt(fs)
    for n in notes:
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        r = q.add_run(n)
        r.font.size = Pt(8)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    t1 = table1()
    t1.to_csv(os.path.join(OUT, "Table1_cohorts.csv"), index=False)
    word_table(t1, os.path.join(OUT, "Table1_cohorts.docx"),
               "Table 1. GEO cohorts and samples analysed (primary sample sets).",
               ["Samples per stage are those in the primary analysis, after the pre-specified tissue-identity exclusions and, "
                "for GSE41258, one array per patient per stage. – = stage not present in the cohort.",
                "Patients: number of distinct patient identifiers in the GEO sample annotation; GSE20916 and GSE37364 "
                "do not provide patient identifiers.",
                "Source publication: PubMed ID listed in the GEO series record, verified in PubMed; GSE117606 lists none "
                "(to verify)."], landscape=True, fs=8)
    t2 = table2()
    t2.to_csv(os.path.join(OUT, "Table2_validation.csv"), index=False)
    word_table(t2, os.path.join(OUT, "Table2_validation.docx"),
               "Table 2. Pre-registered step 5 permutation test and secondary persistent-homology tests in the GEO cohorts "
               "(primary analysis; Holm-adjusted p-values).",
               ["Values are Holm-adjusted p-values; p < 0.05 in bold. Step 5 test: 1-Wasserstein distance between stage "
                "persistence diagrams (H0 + H1), 20 subsamples of 80% of the smallest stage per group, patient-level "
                "permutation test with 999 permutations, Holm correction across the transitions of each cohort "
                "(three in GSE41258, two elsewhere).",
                "Tipping-point rule (pre-registered): shape-only Holm p < 0.05 in both unadjusted and adjusted data.",
                "Unadjusted = expression as measured; adjusted = stage-blind composition adjustment (per-gene regression on "
                "epithelial, stromal, immune and, where relevant, liver scores, with stage not in the model). Raw = Euclidean "
                "distances; shape-only = distances divided by each subsample's median pairwise distance.",
                "Main dimension: among the version × mode combinations in which the combined (H0 + H1) test is significant, "
                "the homology dimension whose separate test (H0 only or H1 only; addendum 12) reproduces more of them; "
                "detail gives the counts. – = no significant combined result.",
                "H1 landscape: patient-level permutation test (999 permutations) of the L2 distance between mean H1 "
                "persistence landscapes, Holm across transitions (addendum 12). H0 landscapes were not tested (D26).",
                "Sensitivity analyses (5 and 20 PCs; 1,000 and 5,000 genes; alternative sample sets) are reported in a "
                "supplementary table."], landscape=True, pcols=PCOLS, fs=7)
    print(t1.to_string())
    print(t2.to_string())


if __name__ == "__main__":
    main()
