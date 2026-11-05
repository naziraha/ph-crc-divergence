"""
Supplementary Table S_secondary_analyses: secondary and sensitivity results that were computed under the prespecified
addenda but not shown in any other item. Every value is COPIED from an existing result file; nothing is recomputed.
A value that no file stores is left empty (and listed in the caption).

  A  differentiation score (addendum 11 Part B + Amendment 2), harmonised symbols:
       results/hgnc/addendum11/partB_differentiation/{B2_variability_level,B3_dispersion_decomposition,B4_cms_relation}.csv;
       B5 calls copied as text from summary.txt
  B  barcode summaries (addendum 12): results/addendum12/barcode_summary_tests.csv (all 336 tests)
  C  CMS sensitivity (addendum 9): results/hgnc/addendum9/dispersion_decomposition.csv (identical to results/addendum9/),
       results/addendum9/dispersion_per_subtype.csv; overall calls copied as text from results/addendum9/summary.txt
  D  chromosome-arm sensitivity (addendum 11 Part A), harmonised symbols:
       results/hgnc/addendum11/partA_chr/arm_scores_tests.csv (unadjusted / full-arm rows; gain-like share)
  E  addendum 8 sensitivity: results/addendum8/summary.txt (all-gene score line, copied as text) and
       results/addendum8/continuity_results_secondary.csv
Only string handling is done (e.g. H0/H1 split from the summary name, reading numbers out of summary.txt lines) and
comparisons with 0.05 for the significance flag.
Output: outputs/supplementary/TableS_secondary_analyses.{csv,docx}

Usage: python scripts/make_table_s_secondary_analyses.py
"""
import os
import re

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "outputs", "supplementary")
NAME = "TableS_secondary_analyses"
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
COLS = ["section", "block", "cohort", "data", "labels_or_set", "mode", "transition", "dimension", "measure", "value",
        "ci_low", "ci_high", "calibration", "z", "n", "p", "p_adjusted", "call", "significant", "source"]
TR = {"Normal->Adenoma": "N → A", "Adenoma->Carcinoma": "A → C", "Carcinoma->Metastasis": "C → M"}
rows = []


def add(**k):
    rows.append({c: k.get(c, np.nan) for c in COLS})


def lines_after(path, header, n):
    L = open(path, encoding="utf-8").read().split("\n")
    i = next(k for k, l in enumerate(L) if l.startswith(header))
    return [l.strip() for l in L[i + 1:i + 1 + n]]


# ---------------- A. differentiation score (addendum 11 Part B), harmonised symbols
PB = os.path.join(R, "hgnc", "addendum11", "partB_differentiation")
src = "results/hgnc/addendum11/partB_differentiation/"
b2 = pd.read_csv(os.path.join(PB, "B2_variability_level.csv"))
for _, r in b2.iterrows():
    base = dict(section="A", block="A1 variability and level, carcinoma vs adenoma", cohort=r.cohort,
                labels_or_set=f"{r.score} score, {r.scores} marker set" if r.scores != "unadjusted"
                else f"{r.score} score, unadjusted expression", source=src + "B2_variability_level.csv")
    add(**base, measure="V (MAD carcinoma ÷ MAD adenoma); Brown–Forsythe p", value=r.V, p=r.p_BF)
    add(**base, measure="median difference (carcinoma − adenoma); Mann–Whitney p", value=r.median_diff, p=r.p_MW)
b3 = pd.read_csv(os.path.join(PB, "B3_dispersion_decomposition.csv"))
for _, r in b3.iterrows():
    base = dict(section="A", block="A2 share of excess carcinoma dispersion", cohort=r.cohort,
                call="informative" if r.informative else "not informative (no excess dispersion)",
                source=src + "B3_dispersion_decomposition.csv")
    add(**base, measure="D_total (carcinoma ÷ adenoma); permutation p", value=r.D_total, p=r.p_total)
    add(**base, measure="f, differentiation (primary marker set) vs f_rand95", value=r.diff_f, calibration=r.diff_f_rand95)
    add(**base, measure="f, differentiation (full marker set) vs f_rand95", value=r.diff_fullset_f,
        calibration=r.diff_fullset_f_rand95)
    add(**base, measure="f, CMS centroids", value=r.f_CMS)
    add(**base, measure="f, CMS + differentiation", value=r.f_CMS_diff)
    add(**base, measure="incremental f beyond CMS vs its calibration (95th pct)", value=r.incremental_f,
        calibration=r.incremental_f_rand95)
b4 = pd.read_csv(os.path.join(PB, "B4_cms_relation.csv"))
for _, r in b4.iterrows():
    add(section="A", block="A3 differentiation score by CMS among carcinomas", cohort=r.cohort,
        measure="R² (score | CMS); Kruskal–Wallis p", value=r.R2_score_CMS, p=r.p_KW, source=src + "B4_cms_relation.csv")
for l in lines_after(os.path.join(PB, "summary.txt"), "CALLS (B5", 7):
    what, _, call = l.partition("->") if "->" in l else (l, "", "")
    if not call:                                       # share lines have the form "<what>: MINOR (...)"
        what, _, call = l.partition(": ")
    add(section="A", block="A4 prespecified calls (B5 / Amendment 2)", measure=what.strip(), call=call.strip(),
        source=src + "summary.txt (CALLS)")

# ---------------- B. barcode summaries (addendum 12)
bc = pd.read_csv(os.path.join(R, "addendum12", "barcode_summary_tests.csv"))
for _, r in bc.iterrows():
    dim, _, summ = r.statistic.partition("_")
    add(section="B", block="B barcode summaries", cohort=r.cohort, data=r.version, mode=r["mode"],
        transition=TR.get(r.transition, r.transition), dimension=dim, measure=summ.replace("_", " "), value=r.obs, z=r.z,
        p=r.p, p_adjusted=r.p_holm_both, significant="yes" if r.p_holm_both < 0.05 else "no",
        source="results/addendum12/barcode_summary_tests.csv")

# ---------------- C. CMS sensitivity (addendum 9)
dd = pd.read_csv(os.path.join(R, "hgnc", "addendum9", "dispersion_decomposition.csv"))
for _, r in dd.iterrows():
    add(section="C", block="C1 PERMANOVA of CMS among carcinomas and per-cohort call", cohort=r.cohort,
        labels_or_set=r["labels"], measure="PERMANOVA R² (CMS); permutation p", value=r.R2_CMS, p=r.p_CMS, call=r.call,
        n=r.n_carcinoma_within, source="results/hgnc/addendum9/dispersion_decomposition.csv")
for l in lines_after(os.path.join(R, "addendum9", "summary.txt"), "5. SECTION 5 OVERALL CONCLUSION", 3):
    m = re.match(r"(\S+)\s+->\s+(.+?)\s+\(calls:", l)
    add(section="C", block="C2 overall call (addendum 9 section 5)", labels_or_set=m.group(1), call=m.group(2),
        measure="overall conclusion", source="results/addendum9/summary.txt (section 5)")
ps = pd.read_csv(os.path.join(R, "addendum9", "dispersion_per_subtype.csv"))
for _, r in ps.iterrows():
    add(section="C", block="C3 per-subtype PERMDISP (CMS group ≥ 10 carcinomas vs all adenomas; primary labels)",
        cohort=r.cohort, labels_or_set=r.CMS, n=r.n, measure="D (subtype ÷ adenoma); permutation p, Holm within cohort",
        value=r.D, p=r.p, p_adjusted=r.p_holm, source="results/addendum9/dispersion_per_subtype.csv")

# ---------------- D. chromosome-arm sensitivity (addendum 11 Part A), harmonised symbols
at = pd.read_csv(os.path.join(R, "hgnc", "addendum11", "partA_chr", "arm_scores_tests.csv"))
ARMS = ["20q", "20p", "13q", "8q", "7p", "7q"]
for lab, ver, genes in (("unadjusted expression", "unadjusted", "excl122"), ("full arm (122 genes kept)", "adjusted", "fullarm")):
    for arm in ARMS:
        for c in GEO:
            r = at[(at.version == ver) & (at.genes == genes) & (at.arm == arm) & (at.cohort == c)]
            assert len(r) == 1, (ver, genes, arm, c)
            r = r.iloc[0]
            add(section="D", block=f"D1 arm score sensitivity: {lab}", cohort=c, data=ver, labels_or_set=f"{arm} ({genes})",
                measure="V_arm; Brown–Forsythe p (rank: not stored for this variant)", value=r.V_arm, p=r.p_BF,
                n=r.n_genes, source="results/hgnc/addendum11/partA_chr/arm_scores_tests.csv")
for c in GEO:
    r = at[(at.version == "adjusted") & (at.genes == "excl122") & (at.arm == "20q") & (at.cohort == c)].iloc[0]
    add(section="D", block="D2 20q gain-like share (primary: adjusted, 122 genes excluded)", cohort=c, data="adjusted",
        labels_or_set="20q (excl122)", measure="share of carcinomas with 20q score > adenoma 97.5th percentile",
        value=r.gainlike_share_carc, source="results/hgnc/addendum11/partA_chr/arm_scores_tests.csv")

# ---------------- E. addendum 8 sensitivity
s8 = open(os.path.join(R, "addendum8", "summary.txt"), encoding="utf-8").read()
line = re.search(r"Sensitivity \(all-gene score\): (.+)", s8).group(1)
for c in GEO:
    m = re.search(rf"{c} cand (\d+) vs cross-fitted (\d+), gradual (\d+)%", line)
    base = dict(section="E", block="E1 part (b) with the all-gene PC1 score (sensitivity)", cohort=c,
                source="results/addendum8/summary.txt (Sensitivity (all-gene score) line)")
    add(**base, measure="candidate switch genes (FDR), all-gene score [cross-fitted primary in brackets]",
        value=int(m.group(1)), call=f"[primary {m.group(2)}]")
    add(**base, measure="gradual genes (% of progression-associated genes, as printed, rounded)", value=int(m.group(3)))
    for cls in ("intermediate", "switch-like", "step-like unresolvable", "non-monotonic"):
        add(**base, measure=f"{cls} genes (% of progression-associated genes): not stored")
m = re.search(r"Sensitivity \(all-gene score\): (\d+) replicated switch genes", s8)
add(section="E", block="E1 part (b) with the all-gene PC1 score (sensitivity)", measure="replicated switch genes, all-gene score",
    value=int(m.group(1)), source="results/addendum8/summary.txt")
cs = pd.read_csv(os.path.join(R, "addendum8", "continuity_results_secondary.csv")).set_index("cohort")
for c in GEO:
    r = cs.loc[c]
    for key, lab in (("NA", "N → A"), ("AC", "A → C")):
        base = dict(section="E", block="E2 part (a) with the secondary N→C centroid-line score (descriptive)", cohort=c,
                    transition=lab, source="results/addendum8/continuity_results_secondary.csv")
        add(**base, measure="valley ratio (95% CI); dip test Holm p (pooled pair); continuity call", value=r[f"{key}_VR"],
            ci_low=r[f"{key}_VR_lo"], ci_high=r[f"{key}_VR_hi"], p_adjusted=r[f"dip_holm_{key}"], call=r[f"{key}_call"])
        add(**base, measure="AUC of the score, earlier vs later stage", value=r[f"{key}_AUC"])

T = pd.DataFrame(rows, columns=COLS)


def fmt(v, kind="num"):
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if kind == "p":
        return f"{v:.3g}"
    if float(v).is_integer() and abs(v) < 1e6:
        return f"{int(v)}"
    return f"{v:.3f}"


def word(path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(8)
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for mm in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, mm, Cm(1.2))
    doc.add_paragraph().add_run("Supplementary Table S_secondary_analyses. Secondary and sensitivity results of the "
                                "prespecified addenda not shown in other items (values copied from existing result files).").bold = True
    show = [("cohort", "Cohort"), ("data", "Data"), ("labels_or_set", "Labels / set"), ("mode", "Mode"),
            ("transition", "Transition"), ("dimension", "Dim."), ("measure", "Measure"), ("value", "Value"),
            ("ci_low", "CI low"), ("ci_high", "CI high"), ("calibration", "Calibration"), ("z", "z"), ("n", "n"),
            ("p", "p"), ("p_adjusted", "Adjusted p"), ("call", "Call")]
    titles = {"A": "A. Differentiation score (addendum 11 Part B; harmonised symbols)",
              "B": "B. Barcode summaries (addendum 12); bold = adjusted p < 0.05",
              "C": "C. CMS sensitivity (addendum 9)",
              "D": "D. Chromosome-arm sensitivity (addendum 11 Part A; harmonised symbols)",
              "E": "E. Addendum 8 sensitivity analyses"}
    for s in "ABCDE":
        sub = T[T.section == s]
        doc.add_paragraph().add_run(titles[s]).bold = True
        for blk, g in sub.groupby("block", sort=False):
            doc.add_paragraph().add_run(blk).italic = True
            cols = [(c, h) for c, h in show if g[c].notna().any() and not (g[c].astype(str) == "").all()]
            tab = doc.add_table(rows=1, cols=len(cols))
            tab.style = "Table Grid"
            for j, (_, h) in enumerate(cols):
                run = tab.rows[0].cells[j].paragraphs[0].add_run(h)
                run.bold, run.font.size = True, Pt(7)
            for _, r in g.iterrows():
                cells = tab.add_row().cells
                for j, (c, _) in enumerate(cols):
                    run = cells[j].paragraphs[0].add_run(fmt(r[c], "p" if c in ("p", "p_adjusted") else "num"))
                    run.font.size = Pt(7)
                    run.bold = r["significant"] == "yes"
    for n in ("All values are copied from existing result files (sources in the CSV column 'source' and in the caption); "
              "no analysis was rerun. Empty cells: the value is not stored in any result file.",
              "Barcode-summary adjusted p: Holm across the transitions of the cohort, then across the 7 summaries "
              "(addendum 12). Part A and C p-values are as stored: Part B D_resid p-values are descriptive only (D24)."):
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(7.5)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    T.to_csv(os.path.join(OUT, f"{NAME}.csv"), index=False)
    word(os.path.join(OUT, f"{NAME}.docx"))
    print("rows per section:", T.groupby("section").size().to_dict())
    print("rows per block:")
    print(T.groupby(["section", "block"], sort=False).size().to_string())
    sig = T[T.significant == "yes"]
    # check: significant barcode rows = the 23 listed in results/addendum12/summary.txt
    assert len(sig) == 23, len(sig)
    print(f"\nsignificant barcode-summary rows ({len(sig)}):")
    print(sig[["cohort", "data", "mode", "transition", "dimension", "measure", "z", "p_adjusted"]].to_string(index=False))


if __name__ == "__main__":
    main()
