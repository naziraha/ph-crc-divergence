"""
Main Tables 1 (cohorts) and 2 (primary validation results) for the GEO cohorts. All numbers are copied from existing
files; no statistical tests are run. Outputs (outputs/main_tables/): Table1_cohorts.csv/.docx,
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
OUT = os.path.join(ROOT, "outputs", "main_tables")
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
       "GSE117606": "to verify (no PubMed ID in the GEO record; the related adenoma-only series GSE117607 lists none either)"}
SUPP = os.path.join(ROOT, "outputs", "supplementary")
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


def shared_patients(c, tr):
    """Patients with a sample in both stages of the comparison (primary sample set); None if no patient IDs."""
    st = SETS[c]
    emb = pd.read_csv(os.path.join(R, "step03", "embeddings", f"{c}_{st}_adjusted_g2000_pc10.csv"), index_col=0)
    lab = pd.read_csv(os.path.join(PROC, f"{c}_{st}_samples_qc.csv"), dtype={PATIENT_COL[c]: str}).set_index("sample")
    pid = lab.loc[emb.index, PATIENT_COL[c]]
    if pid.isna().all():
        return None
    a, b = tr.split("->")
    return len(set(pid[emb.stage == a].dropna()) & set(pid[emb.stage == b].dropna()))


def table2():
    v = pd.read_csv(os.path.join(R, "source_data", "validation_summary_primary.csv"))
    dims = pd.read_csv(os.path.join(R, "addendum12", "dimension_split_all_cohorts.csv"))
    s5 = pd.read_csv(os.path.join(R, "step05", "permutation_results_all.csv"))
    s5 = s5[(s5.analysis == "primary")].groupby("transition").shared_patients.first()   # recorded by step 5 (GSE41258)
    rows = []
    for _, r in v.iterrows():
        main, detail = dim_call(dims, r.cohort, r.transition)
        sh = shared_patients(r.cohort, r.transition)
        if r.cohort == "GSE41258" and sh != int(s5[r.transition]):
            raise SystemExit(f"shared-patient count differs from step 5 record: {r.transition} {sh} vs {s5[r.transition]}")
        rows.append({"Cohort": r.cohort, "Transition": TRL[r.transition], "n (from / to)": f"{r.n_from} / {r.n_to}",
                     "Patients in both groups": "not known" if sh is None else str(sh),
                     "Raw, unadjusted": r.raw_p_holm_unadj, "Raw, adjusted": r.raw_p_holm_adjus,
                     "Shape-only, unadjusted": r.shape_p_holm_unadj, "Shape-only, adjusted": r.shape_p_holm_adjus,
                     "Tipping-point rule met": r.rule_met,
                     "Main dimension": main, "Dimension detail": detail,
                     "H1 landscape shape-only, unadjusted": r.H1landscape_shape_p_holm_unadj,
                     "H1 landscape shape-only, adjusted": r.H1landscape_shape_p_holm_adjus})
    return pd.DataFrame(rows)


def table_s_landscape():
    v = pd.read_csv(os.path.join(R, "source_data", "validation_summary_primary.csv"))
    return pd.DataFrame([{"Cohort": r.cohort, "Transition": TRL[r.transition], "n (from / to)": f"{r.n_from} / {r.n_to}",
                          "Raw, unadjusted": r.H1landscape_raw_p_holm_unadj,
                          "Raw, adjusted": r.H1landscape_raw_p_holm_adjus,
                          "Shape-only, unadjusted": r.H1landscape_shape_p_holm_unadj,
                          "Shape-only, adjusted": r.H1landscape_shape_p_holm_adjus} for _, r in v.iterrows()])


PCOLS = ["Raw, unadjusted", "Raw, adjusted", "Shape-only, unadjusted", "Shape-only, adjusted",
         "H1 landscape shape-only, unadjusted", "H1 landscape shape-only, adjusted"]


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
                "Sample sets. GSE41258 'main': of 390 arrays, normal colon, adenoma (including 2 microadenomas), primary "
                "carcinoma and liver metastasis (338); excluded were lung metastases (20; a with-lung sensitivity set adds "
                "them back), normal liver (13), normal lung (7) and cell lines (12); then 2 arrays excluded for tissue "
                "identity (expression of normal liver; Müllerian profile) and 15 additional arrays of the same patient and "
                "stage removed (re-hybridisations or further pieces of the same tissue), leaving 321. GSE20916 'bulk': the "
                "105 whole-tissue sections (normal colon, adenoma, adenocarcinoma); the 40 laser-microdissected samples "
                "(crypt epithelium and mucosa) form a separate 'micro' set analysed descriptively only. GSE37364 'bulk': all "
                "94 biopsies (normal mucosa; low- and high-grade adenoma; Dukes A/B and C/D carcinoma). GSE117606 'main': of "
                "208 arrays, adjacent mucosa (normal), adenoma and tumour (carcinoma), one array per patient per stage (195); "
                "10 sessile serrated adenomas (described separately) and 3 additional 'A-coded' carcinomas of patients whose "
                "designated tumour is another array were set aside. Apart from the two GSE41258 identity exclusions, arrays "
                "flagged by the identity screen or (GSE117606) by array quality were kept in the primary sets, as "
                "pre-specified; sensitivity analyses remove or add them.",
                "Patients: number of distinct patient identifiers in the GEO sample annotation; GSE20916 and GSE37364 "
                "do not provide patient identifiers. GSE117606: the 195 analysed arrays (and all 208 arrays of the series) "
                "carry 75 distinct patient identifiers (sample characteristic 'patient id', matching the sample titles), "
                "whereas the free-text series summary states 70 patients; the difference cannot be resolved from the GEO "
                "record.",
                "Source publication: PubMed ID listed in the GEO series record, verified in PubMed; GSE117606 lists none, "
                "and neither does the related adenoma-only series GSE117607 (to verify)."], landscape=True, fs=8)
    t2 = table2()
    t2.to_csv(os.path.join(OUT, "Table2_validation.csv"), index=False)
    word_table(t2.drop(columns=["Dimension detail"]), os.path.join(OUT, "Table2_validation.docx"),   # detail kept in CSV
               "Table 2. Step 5 permutation test (prespecified before analysis) and secondary persistent-homology tests in the GEO cohorts "
               "(primary analysis; Holm-adjusted p-values).",
               ["Values are Holm-adjusted p-values; p < 0.05 in bold. Step 5 test: the difference in persistence-diagram shape "
                "between stages, measured as the mean 1-Wasserstein distance between diagrams of different stages minus the "
                "mean distance within stages (H0 and H1 combined; 20 subsamples per stage, each of 80% of the smallest "
                "stage), label-permutation "
                "test with 999 permutations (groups aligned before permuting), Holm correction across the transitions of "
                "each cohort (three in GSE41258, two elsewhere).",
                "Patients in both groups: patients contributing a sample to both stages of the comparison (GSE41258: as "
                "recorded by step 5; GSE117606: from the GEO 'patient id'; GSE20916 and GSE37364: not known, no patient "
                "identifiers). Within each stage every sample comes from a different patient (GSE41258 after one-array-per-"
                "patient de-duplication). In the primary step 5 test, stage labels were permuted between individual samples, "
                "and the samples of shared patients were kept in both groups; the permutation did not keep a patient's "
                "samples together. Shared patients make the two groups more alike, which the analysis plan, prespecified before analysis, expected to "
                "make the test conservative. Two pre-specified analyses addressed this: GSE41258 'one per patient' "
                "(each shared patient assigned at random to one group only) and GSE117606 'paired' (labels swapped within "
                "patients); neither met the tipping-point rule (supplementary sensitivity table).",
                "Tipping-point rule (prespecified before analysis): shape-only Holm p < 0.05 in both unadjusted and adjusted data.",
                "Unadjusted = expression as measured; adjusted = stage-blind composition adjustment (per-gene regression on "
                "epithelial, stromal, immune and, where relevant, liver scores, with stage not in the model). Raw = Euclidean "
                "distances; shape-only = distances divided by each subsample's median pairwise distance.",
                "Main dimension: among the version × mode combinations in which the combined (H0 + H1) test is significant, "
                "the homology dimension whose separate test (H0 only or H1 only; addendum 12) reproduces more of them "
                "(counts per row in the accompanying CSV file). – = no significant combined result.",
                "H1 landscape (shape-only): label-permutation test (999 permutations) of the L2 distance between mean H1 "
                "persistence landscapes, Holm across transitions (addendum 12); raw-distance results are in Supplementary "
                "Table S_H1_landscape. H0 landscapes were not tested (D26).",
                "Sensitivity analyses (5 and 20 PCs; 1,000 and 5,000 genes; alternative sample sets) are reported in a "
                "supplementary table."], landscape=True, pcols=PCOLS, fs=7)
    ts = table_s_landscape()
    os.makedirs(SUPP, exist_ok=True)
    ts.to_csv(os.path.join(SUPP, "TableS_H1_landscape.csv"), index=False)
    word_table(ts, os.path.join(SUPP, "TableS_H1_landscape.docx"),
               "Supplementary Table S_H1_landscape. H1 persistence landscape tests in the GEO cohorts, raw and shape-only, "
               "unadjusted and adjusted (primary analysis; Holm-adjusted p-values).",
               ["Values are Holm-adjusted p-values; p < 0.05 in bold. Test (addendum 12, secondary): L2 distance between "
                "the mean H1 persistence landscapes (first five landscapes) of the two stages, 20 subsamples of 80% of the "
                "smallest stage per group, label-permutation test with 999 permutations (groups aligned before permuting), "
                "Holm correction across the transitions of each cohort.",
                "Raw = Euclidean distances; shape-only = distances divided by each subsample's median pairwise distance. "
                "Unadjusted = expression as measured; adjusted = stage-blind composition adjustment.",
                "H0 landscapes were not tested because the test was underpowered in validation (D26). The two shape-only "
                "columns are also shown in main Table 2."], landscape=False, pcols=["Raw, unadjusted", "Raw, adjusted",
                                                                                "Shape-only, unadjusted", "Shape-only, adjusted"], fs=8)
    print(t1.to_string())
    print(t2.to_string())
    print(ts.to_string())


if __name__ == "__main__":
    main()
