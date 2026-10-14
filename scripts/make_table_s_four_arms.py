"""
Supplementary Table S_four_arms: share of excess carcinoma dispersion explained by removing the four confirmed arm scores
together (20q, 20p, 13q, 8q), with 20q alone for reference (addendum 11 Amendment 2). Copied from
results/addendum11/partB_differentiation/B3_dispersion_decomposition.csv; no statistical tests.
Output: outputs/supplementary/TableS_four_arms.{csv,docx}

Usage: python scripts/make_table_s_four_arms.py
"""
import os

import pandas as pd
from docx import Document
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B3 = os.path.join(ROOT, "results", "hgnc", "addendum11", "partB_differentiation", "B3_dispersion_decomposition.csv")
OUT = os.path.join(ROOT, "outputs", "supplementary")
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]


def build():
    b = pd.read_csv(B3).set_index("cohort")
    rows = []
    for c in GEO:
        r = b.loc[c]
        inf = bool(r.informative)
        rows.append({"Cohort": c, "Informative": "yes" if inf else "no",
                     "D_total (p)": f"{r.D_total:.2f} ({r.p_total:.3f})",
                     "Four arms: D_resid (p, descriptive)": f"{r.chr4_D_resid:.2f} ({r.chr4_p_resid_descriptive:.3f})",
                     "Four arms: f": f"{r.chr4_f:.2f}" if inf else "—",
                     "Four arms: f_rand95": f"{r.chr4_f_rand95:.2f}" if inf else "—",
                     "Four arms: f > f_rand95": ("yes" if r.chr4_f > r.chr4_f_rand95 else "no") if inf else "—",
                     "20q alone: f (f_rand95)": f"{r.chr20q_f:.2f} ({r.chr20q_f_rand95:.2f})" if inf else "—"})
    return pd.DataFrame(rows)


NOTES = [
    "Four arms = the arm scores of 20q, 20p, 13q and 8q (adjusted expression, the 122 variable genes excluded; current HGNC symbols), removed "
    "together: within each stage, every principal component is regressed on the four scores (with intercept) and the "
    "residuals are kept. D_total = bias-adjusted carcinoma/adenoma dispersion ratio; D_resid = the same after removal; "
    "f = (D_total − D_resid) ÷ (D_total − 1); f_rand95 = 95th percentile of f when the scores are permuted within stage "
    "(999 times). 'Substantial' (pre-set) = f ≥ 0.25 and f > f_rand95 in ≥ 2 informative cohorts.",
    "Negative f: D_resid > D_total, i.e. removing the arm scores shrinks adenoma dispersion proportionally more than "
    "carcinoma dispersion, so the carcinoma excess grows. It means the arm scores do not account for the excess there; it "
    "is not a 'negative explanation'.",
    "Caution (D24): the residual-dispersion p-values are descriptive only, because the procedure slightly exceeded its "
    "false-positive limit in simulation (18/200 = 0.090; limit 0.089). The joint four-arm removal was not separately "
    "simulated, so its f values should be interpreted cautiously.",
    "Not informative: GSE41258 carcinomas are not significantly more dispersed than adenomas after adjustment, so f is "
    "not defined (D_total p from addendum 11's own permutations; addendum 9 gives p = 0.309).",
]


def word(df, path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(10)
    sec = doc.sections[0]
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(1.5))
    doc.add_paragraph().add_run("Supplementary Table S_four_arms. Excess carcinoma dispersion explained by the four "
                                "confirmed chromosome-arm scores together (20q, 20p, 13q, 8q).").bold = True
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    for j, c in enumerate(df.columns):
        r = t.rows[0].cells[j].paragraphs[0].add_run(c)
        r.bold, r.font.size = True, Pt(8)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, c in enumerate(df.columns):
            cells[j].paragraphs[0].add_run(str(row[c])).font.size = Pt(8)
    for n in NOTES:
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = build()
    df.to_csv(os.path.join(OUT, "TableS_four_arms.csv"), index=False)
    word(df, os.path.join(OUT, "TableS_four_arms.docx"))
    pd.set_option("display.width", 220)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
