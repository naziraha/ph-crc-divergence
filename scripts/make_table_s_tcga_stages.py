"""
Supplementary Table S_TCGA_stages (addendum 10): TCGA-COAD stage results, copied from existing files; no tests.
(A) step 5 test per stage transition (raw and shape-only Holm p, unadjusted and adjusted, tipping-point rule);
(B) PERMDISP dispersion ratio D and Holm p per transition (bias-adjusted, addendum 10 section 4.1);
(C) CMS proportions by stage (primary and confident calls; descriptive; CMS4 stromal caution);
(D) secondary, descriptive: solid-tissue normal vs Stage I (single comparison, no Holm) and the paired sensitivity.
Output: paper/final/supplementary/TableS_TCGA_stages.{csv,docx}

Usage: python scripts/make_table_s_tcga_stages.py
"""
import os

import pandas as pd
from docx import Document
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
A10 = os.path.join(R, "addendum10")
OUT = os.path.join(ROOT, "paper", "final", "supplementary")
TR = ["Stage I->Stage II", "Stage II->Stage III", "Stage III->Stage IV"]
STG = ["Stage I", "Stage II", "Stage III", "Stage IV"]


def lab(t):
    a, b = t.split("->")
    return f"{a.replace('Stage ', '')} → {b.replace('Stage ', '')}"


def part_a():
    d = pd.read_csv(os.path.join(R, "step05", "TCGA", "permutation_results_all.csv"))
    d = d[d.analysis == "primary"]
    rows = []
    for t in TR:
        g = d[d.transition == t]
        get = lambda v, m: float(g[(g.version == v) & (g["mode"] == m)].p_holm.iloc[0])
        su, sa = get("unadjusted", "shape-only"), get("adjusted", "shape-only")
        rows.append({"section": "A. Step 5 test", "transition": lab(t), "n (from / to)": f"{int(g.n_from.iloc[0])} / {int(g.n_to.iloc[0])}",
                     "m": int(g.subsample_size.iloc[0]),
                     "raw Holm p, unadjusted": get("unadjusted", "raw"), "raw Holm p, adjusted": get("adjusted", "raw"),
                     "shape-only Holm p, unadjusted": su, "shape-only Holm p, adjusted": sa,
                     "rule met": "yes" if (su < 0.05 and sa < 0.05) else "no"})
    return pd.DataFrame(rows)


def part_b():
    d = pd.read_csv(os.path.join(A10, "spread_H0H1", "permdisp_by_stage.csv"))
    rows = []
    for t in TR:
        g = d[d.transition == t].set_index("version")
        rows.append({"section": "B. PERMDISP", "transition": lab(t),
                     "D, unadjusted": round(float(g.loc["unadjusted", "D"]), 3), "Holm p, unadjusted": float(g.loc["unadjusted", "p_holm"]),
                     "D, adjusted": round(float(g.loc["adjusted", "D"]), 3), "Holm p, adjusted": float(g.loc["adjusted", "p_holm"])})
    return pd.DataFrame(rows)


def part_c():
    p = pd.read_csv(os.path.join(A10, "cms", "cms_by_stage_primary.csv")).set_index("stage")
    c = pd.read_csv(os.path.join(A10, "cms", "cms_by_stage_confident.csv")).set_index("stage")
    rows = []
    for s in STG:
        n = int(p.loc[s].sum())
        rec = {"section": "C. CMS by stage", "stage": s, "n": n}
        for k in ("CMS1", "CMS2", "CMS3", "CMS4"):
            rec[f"{k} primary, n (%)"] = f"{int(p.loc[s, k])} ({100 * p.loc[s, k] / n:.0f}%)"
        rec["unclassified (confident calls), n (%)"] = f"{int(c.loc[s, 'unclassified'])} ({100 * c.loc[s, 'unclassified'] / n:.0f}%)"
        rows.append(rec)
    return pd.DataFrame(rows)


def part_d():
    d = pd.read_csv(os.path.join(A10, "normal_stageI", "normal_stageI_results.csv"))
    rows = []
    for _, r in d.iterrows():
        main = r.analysis == "main"
        rows.append({"section": "D. Normal vs Stage I (secondary, descriptive)",
                     "comparison": "normal vs Stage I, one sample per patient" if main else "normal vs own tumour, paired (any stage)",
                     "n": f"{int(r.n_from)} normal / {int(r.n_to)} Stage I" if main else f"{int(r.n_pairs)} pairs",
                     "m": int(r.subsample_size), "version": r.version, "mode": r["mode"], "z": round(float(r.z), 2),
                     "p (single comparison, no Holm)": float(r.p)})
    return pd.DataFrame(rows)


def fmt(v):
    return f"{v:.3f}" if isinstance(v, float) else str(v)


def add(doc, df, cols, title, pcols=()):
    doc.add_paragraph().add_run(title).bold = True
    t = doc.add_table(rows=1, cols=len(cols))
    t.style = "Table Grid"
    for j, c in enumerate(cols):
        r = t.rows[0].cells[j].paragraphs[0].add_run(c)
        r.bold, r.font.size = True, Pt(7.5)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, c in enumerate(cols):
            run = cells[j].paragraphs[0].add_run(fmt(row[c]))
            run.font.size = Pt(7.5)
            if c in pcols and isinstance(row[c], float) and row[c] < 0.05:
                run.bold = True


def main():
    os.makedirs(OUT, exist_ok=True)
    A, B, C, D = part_a(), part_b(), part_c(), part_d()
    pd.concat([A, B, C, D], ignore_index=True).to_csv(os.path.join(OUT, "TableS_TCGA_stages.csv"), index=False)
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    for m in ("left_margin", "right_margin"):
        setattr(doc.sections[0], m, Cm(1.5))
    doc.add_paragraph().add_run("Supplementary Table S_TCGA_stages. TCGA-COAD primary tumours by stage (addendum 10).").bold = True
    pa = ["raw Holm p, unadjusted", "raw Holm p, adjusted", "shape-only Holm p, unadjusted", "shape-only Holm p, adjusted"]
    add(doc, A, ["transition", "n (from / to)", "m"] + pa + ["rule met"],
        "A. Step 5 test between adjacent stages (primary analysis; Holm across the three transitions)", pa)
    pb = ["Holm p, unadjusted", "Holm p, adjusted"]
    add(doc, B, ["transition", "D, unadjusted", "Holm p, unadjusted", "D, adjusted", "Holm p, adjusted"],
        "B. Between-sample spread (PERMDISP dispersion ratio D = later ÷ earlier stage)", pb)
    add(doc, C, ["stage", "n", "CMS1 primary, n (%)", "CMS2 primary, n (%)", "CMS3 primary, n (%)", "CMS4 primary, n (%)",
                 "unclassified (confident calls), n (%)"], "C. Consensus molecular subtypes by stage (descriptive)")
    add(doc, D, ["comparison", "n", "m", "version", "mode", "z", "p (single comparison, no Holm)"],
        "D. Secondary, descriptive: solid-tissue normal vs Stage I", ["p (single comparison, no Holm)"])
    for n in ("A: 1-Wasserstein distance between stage persistence diagrams (H0 + H1), k = 20 subsamples of m = 31 per group "
              "(80% of the smallest stage, 39 Stage IV), label-permutation test (one sample per patient per stage; 999 "
              "permutations), Holm across the three transitions per version and mode. Tipping-point rule: shape-only Holm "
              "p < 0.05 in both unadjusted and adjusted data. Bold = p < 0.05.",
              "B: bias-adjusted PERMDISP distances √(n/(n−1))·‖x − c‖ in the 10-PC embedding, 999 permutations, Holm across "
              "the three transitions per version; D > 1 = later stage more spread out.",
              "C: CMScaller on unadjusted expression; primary = nearest template (every tumour assigned); unclassified = no "
              "significant template match under the confident rule. Caution: CMS4 calls on bulk expression partly reflect the "
              "stromal content of the sample, so the rise of CMS4 with stage may partly reflect more stroma in advanced "
              "tumours rather than tumour-cell biology alone. Descriptive; no test.",
              "D: secondary and descriptive (pre-specified in addendum 10 section 4.3, D22). 'Main' = 41 adjacent normals vs "
              "45 Stage I tumours, one sample per patient: the 3 patients with both were assigned at random to one group in "
              "each test, so n is 39–41 normals and 42–44 Stage I depending on the row. 'Paired' = 25 patients' normal vs own tumour (Stage I 3, II 13, III 4, IV 5), "
              "labels swapped within patients. Single comparisons, no Holm. TCGA normals are tumour-adjacent tissue (possible "
              "field effect). TCGA has no adenomas, so this cannot locate a change along normal → adenoma → carcinoma."):
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(os.path.join(OUT, "TableS_TCGA_stages.docx"))
    pd.set_option("display.width", 220)
    for df in (A, B, C, D):
        print(df.drop(columns="section").to_string(index=False))


if __name__ == "__main__":
    main()
