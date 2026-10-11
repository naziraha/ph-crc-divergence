"""
Supplementary Table S_cohort_differences: four GEO cohorts side by side, to explain why GSE41258 behaves differently in
Figure 4. All numbers are copied from existing files; no statistical tests are run.
Output: paper/final/supplementary/TableS_cohort_differences.{csv,docx}

Usage: python scripts/make_table_s_cohort_differences.py
"""
import os

import pandas as pd
from docx import Document
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
FIN = os.path.join(ROOT, "paper", "final")
OUT = os.path.join(FIN, "supplementary")
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
PLATFORM_SHORT = {"GSE41258": "GPL96 (HG-U133A)", "GSE20916": "GPL570 (HG-U133 Plus 2.0)",
                  "GSE37364": "GPL570 (HG-U133 Plus 2.0)", "GSE117606": "GPL25373 (HT HG-U133+ PM; FFPE)"}
# batch / quality issues: copied from results/addendum5/part2_scandates/summary_corrected.txt (GSE41258),
# results/addendum6/batch_check/batch_confounding.csv and results/step01/GSE117606/summary.txt, results/step03/GSE117606
# (GSE117606); no scan-date check was possible for GSE20916 and GSE37364 (CEL headers not downloaded).
BATCH = {"GSE41258": "Scan date confounded with stage at N→A (day overlap 3.0%, AUC 0.40); not at A→C (overlap 52.0%, "
                     "AUC 0.53) (addendum 5, D9)",
         "GSE20916": "None recorded (scan dates not checked; raw CEL files not downloaded)",
         "GSE37364": "None recorded (scan dates not checked; raw CEL files not downloaded)",
         "GSE117606": "Plate and scan day not confounded with stage (overlap ≥ 99%); 13 arrays flagged for low quality "
                      "(N 4, A 6, C 3; kept, sensitivity analysis); quality loads on unadjusted PC3 (D10)"}


def fmt_p(p):
    return f"{p:.3f}"


def build():
    t1 = pd.read_csv(os.path.join(FIN, "main_tables", "Table1_cohorts.csv")).set_index("Cohort")
    t2 = pd.read_csv(os.path.join(FIN, "main_tables", "Table2_validation.csv"))
    t2 = t2[t2.Transition == "A→C"].set_index("Cohort")
    b = pd.read_csv(os.path.join(R, "step13", "benchmark_results.csv"))
    b = b[b.transition == "Adenoma->Carcinoma"].set_index(["cohort", "version"])
    comp = pd.read_csv(os.path.join(R, "source_data", "fig4_composition_medians.csv")).set_index(["cohort", "score"])
    dec = pd.read_csv(os.path.join(R, "addendum9", "dispersion_decomposition.csv"))
    dec = dec[dec["labels"] == "primary"].set_index("cohort")
    rows = {}
    rows["Platform"] = {c: PLATFORM_SHORT[c] for c in GEO}
    rows["Adenoma / carcinoma n"] = {c: f"{t1.loc[c, 'Adenoma']} / {t1.loc[c, 'Carcinoma']}" for c in GEO}
    rows["Patients shared between adenoma and carcinoma"] = {c: str(t2.loc[c, "Patients in both groups"]) for c in GEO}
    rows["Known batch / quality issue"] = BATCH
    for v, lab in (("unadjusted", "unadjusted"), ("adjusted", "adjusted")):
        rows[f"PERMDISP A→C, {lab}: D (Holm p)"] = {
            c: f"{b.loc[(c, v), 'permdisp_effect']:.2f} ({fmt_p(b.loc[(c, v), 'permdisp_p_holm'])})" for c in GEO}
    for k in ("Epithelial", "Stromal", "Immune"):
        rows[f"Standardised A→C shift, {k.lower()} score"] = {c: f"{comp.loc[(c, k), 'std_diff']:+.2f}" for c in GEO}
    for k in ("Epithelial", "Stromal", "Immune"):
        rows[f"Within-stage variability ratio C ÷ A, {k.lower()} score"] = {
            c: f"{comp.loc[(c, k), 'sd_ratio_C_over_A']:.2f}" for c in GEO}
    rows["Excess carcinoma dispersion explained by CMS: f (f_rand95)"] = {
        c: (f"not informative (D_total {dec.loc[c, 'D_total']:.2f}, p {dec.loc[c, 'p_total']:.3f})"
            if dec.loc[c, "call"] == "not informative"
            else f"{dec.loc[c, 'f']:.2f} ({dec.loc[c, 'f_rand95']:.2f})") for c in GEO}
    return pd.DataFrame(rows).T.rename_axis("Row").reset_index()


NOTES = [
    "A→C = adenoma → carcinoma. PERMDISP: dispersion ratio D (carcinoma ÷ adenoma mean distance to centroid, 10 PCs) and "
    "Holm-adjusted p from the step 13 label-permutation test (one sample per patient per stage; 999 permutations); p < "
    "0.05 in bold. Adjusted = stage-blind composition adjustment.",
    "Standardised shift = (median carcinoma − median adenoma) ÷ robust SD of adenomas and carcinomas combined, per cohort "
    "(step 2 decoupler scores; positive = higher in carcinoma). Variability ratio = robust SD (MAD × 1.4826) in "
    "carcinomas ÷ in adenomas. Descriptive only.",
    "CMS share f: proportion of the excess carcinoma dispersion (bias-adjusted D_total, adjusted data) removed when "
    "dispersion is measured within consensus molecular subtypes; f_rand95 = 95th percentile under random subtype labels "
    "(addendum 9, primary CMS labels). 'Not informative' = carcinomas not significantly more dispersed than adenomas.",
    "What stands out for GSE41258 (numbers only): its adenoma → carcinoma shift in composition is not larger than in the "
    "other cohorts (stromal +1.16, immune +0.52, epithelial −0.81), but composition varies 1.5–1.7 times more among its "
    "carcinomas than among its adenomas for all three scores, which no other cohort shows. It is also the only cohort in "
    "which the adenoma → carcinoma increase in spread is lost after composition adjustment (D 1.19 → 1.07), so the CMS "
    "decomposition is not informative there. Its recorded batch issue concerns normal → adenoma, not adenoma → "
    "carcinoma. Patients shared between adenoma and carcinoma are fewer in GSE41258 (12) than in GSE117606 (56). It is "
    "also the only cohort on the older GPL96 platform and the only one with many more carcinomas than adenomas "
    "(181 vs 48).",
]


def word(df, path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(10)
    sec = doc.sections[0]
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(1.5))
    p = doc.add_paragraph()
    p.add_run("Supplementary Table S_cohort_differences. Characteristics of the four GEO cohorts relevant to the "
              "adenoma → carcinoma change in between-sample heterogeneity (Figure 4).").bold = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    for j, col in enumerate(["", *GEO]):
        r = t.rows[0].cells[j].paragraphs[0].add_run(col)
        r.bold, r.font.size = True, Pt(8)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, col in enumerate(df.columns):
            val = str(row[col])
            run = cells[j].paragraphs[0].add_run(val)
            run.font.size = Pt(8)
            if j == 0:
                run.bold = True
            elif row["Row"].startswith("PERMDISP") and "(" in val and float(val.split("(")[1].rstrip(")")) < 0.05:
                run.bold = True
    for n in NOTES:
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = build()
    df.to_csv(os.path.join(OUT, "TableS_cohort_differences.csv"), index=False)
    word(df, os.path.join(OUT, "TableS_cohort_differences.docx"))
    pd.set_option("display.width", 250, "display.max_colwidth", 60)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
