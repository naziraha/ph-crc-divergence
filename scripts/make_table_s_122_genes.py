"""
Supplementary Table S_122_genes: the 122 replicated 'more variable in carcinoma' genes (addendum 9 rule, current HGNC symbols;
scripts/hgnc_crosscohort.py), in the Figure 6
order (group, then median V), with per-cohort V, Brown-Forsythe p and BH q, and the number of cohorts meeting the
addendum 9 criteria. All values copied from existing files; no statistical tests.
Addendum 9 rule (scripts/addendum9_heterogeneity.py, 'replication of variability'): discovery in a cohort if BH q < 0.05
and V > 1; replication in another cohort (where the gene passed the expression filter) if p < 0.05 and V > 1; a gene is
kept if at least one discovery cohort has >= 2 replicating cohorts.
Output: outputs/supplementary/TableS_122_genes.{csv,docx}

Usage: python scripts/make_table_s_122_genes.py
"""
import os

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A9 = os.path.join(ROOT, "results", "hgnc", "addendum9")      # harmonised symbols (S7)
A9_ORIG = os.path.join(ROOT, "results", "addendum9")
OUT = os.path.join(ROOT, "outputs", "supplementary")
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
GROUPS = ["lower in carcinoma", "higher in carcinoma", "no mean shift"]


def symbol_notes(genes):
    """For each harmonised gene: was it in the original 117 (and under which symbol), or added by harmonisation
    (older symbol(s) in the platform annotations)."""
    import sys
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import check_methods_items as CK
    cur, al = CK.hgnc_map()
    orig = pd.read_csv(os.path.join(A9_ORIG, "variability_replicated.csv"))
    orig = list(orig[orig.direction == "more variable in carcinoma"].gene)
    to_cur = {}
    for c in GEO:
        for s in pd.read_csv(os.path.join(A9_ORIG, f"variability_{c}.csv"), usecols=[0]).iloc[:, 0].astype(str):
            k = CK.classify([s], cur, al)[s]
            if k[0] == "renamed":
                to_cur.setdefault(k[1], set()).add(s)
    orig_cur = {}
    for s in orig:
        k = CK.classify([s], cur, al)[s]
        orig_cur[k[1] if k[0] == "renamed" else s] = s
    out = {}
    for g in genes:
        if g in orig_cur:
            out[g] = "in original 117" + ("" if orig_cur[g] == g else f" (as {orig_cur[g]})")
        else:
            out[g] = "added by symbol harmonisation" + (f" (older symbol {', '.join(sorted(to_cur[g]))})" if g in to_cur else "")
    missing = set(orig_cur) - set(genes)
    if missing:
        raise SystemExit(f"original genes missing from the harmonised list: {missing}")
    return out


def build():
    a = pd.read_csv(os.path.join(A9, "variability_replicated_annotated.csv")).set_index("gene")
    var = {c: pd.read_csv(os.path.join(A9, f"variability_{c}.csv"), index_col=0) for c in GEO}
    notes = symbol_notes(list(a.index))
    rows = []
    for g, r in a.iterrows():
        md = [var[c].loc[g, "mean_diff"] for c in GEO if g in var[c].index]
        grp = "no mean shift" if r.no_mean_shift else ("lower in carcinoma" if np.median(md) < 0 else "higher in carcinoma")
        rec = dict(gene=g, group=grp, median_V=r.median_V, median_mean_diff=float(np.median(md)))
        n_tested = n_rep = n_disc = 0
        for c in GEO:
            if g in var[c].index:
                v = var[c].loc[g]
                n_tested += 1
                n_rep += int(v.p < 0.05 and v.V > 1)
                n_disc += int(v.q < 0.05 and v.V > 1)
                rec[f"V_{c}"], rec[f"p_{c}"], rec[f"q_{c}"] = v.V, v.p, v.q
            else:
                rec[f"V_{c}"] = rec[f"p_{c}"] = rec[f"q_{c}"] = np.nan
        rec.update(n_cohorts_tested=n_tested, n_cohorts_p05_V_gt1=n_rep, n_cohorts_q05_V_gt1=n_disc,
                   discovery_cohorts=r.discovery_cohorts, CMS_template=r.CMS_template, symbol_note=notes[g])
        if n_disc < 1 or n_rep < 3:
            raise SystemExit(f"{g}: does not meet the addendum 9 rule as recomputed")
        if not np.isclose(rec["median_V"], np.median([rec[f"V_{c}"] for c in GEO if pd.notna(rec[f"V_{c}"])])):
            raise SystemExit(f"{g}: median V mismatch")
        rows.append(rec)
    t = pd.DataFrame(rows)
    t["group"] = pd.Categorical(t.group, GROUPS, ordered=True)
    t = t.sort_values(["group", "median_V"], ascending=[True, False]).reset_index(drop=True)
    t.insert(0, "order", np.arange(1, len(t) + 1))
    return t


def word(t, path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Cm(1.3))
    doc.add_paragraph().add_run("Supplementary Table S_122_genes. The 122 genes more variable in carcinoma than in "
                                "adenoma (addendum 9 rule; current HGNC symbols), in Figure 6 order.").bold = True
    cols = ["order", "gene", "group", "median_V"] + [f"V_{c}" for c in GEO] + ["n_cohorts_tested", "n_cohorts_p05_V_gt1",
                                                                             "n_cohorts_q05_V_gt1", "symbol_note"]
    heads = ["#", "Gene", "Group", "Median V"] + [f"V {c}" for c in GEO] + ["Cohorts tested", "Cohorts p < 0.05, V > 1",
                                                                          "Cohorts q < 0.05, V > 1", "Symbol note"]
    tab = doc.add_table(rows=1, cols=len(cols))
    tab.style = "Table Grid"
    for j, h in enumerate(heads):
        run = tab.rows[0].cells[j].paragraphs[0].add_run(h)
        run.bold, run.font.size = True, Pt(7.5)
    for _, r in t.iterrows():
        cells = tab.add_row().cells
        for j, c in enumerate(cols):
            v = r[c]
            if c.startswith("V_") or c == "median_V":
                s = "–" if pd.isna(v) else f"{v:.2f}"
            else:
                s = str(v)
            run = cells[j].paragraphs[0].add_run(s)
            run.font.size = Pt(7)
            if c.startswith("V_") and pd.notna(v):
                cohort = c[2:]
                run.bold = bool(r[f"p_{cohort}"] < 0.05 and v > 1)
    for n in ("V = variability ratio, median absolute deviation in carcinoma ÷ in adenoma, composition-adjusted "
              "expression; bold = Brown–Forsythe p < 0.05 with V > 1 in that cohort; – = gene not tested in that cohort "
              "(not on the platform, or below the cohort's expression filter).",
              "Addendum 9 rule: discovery in a cohort if BH q < 0.05 and V > 1; replication in another cohort if p < 0.05 "
              "and V > 1; a gene is kept if at least one discovery cohort has ≥ 2 replicating cohorts, i.e. the criterion "
              "is met in at least 3 cohorts, not necessarily all 4.",
              "Groups: 'no mean shift' = median |mean difference (carcinoma − adenoma)| < 0.5 log2 across tested cohorts; "
              "otherwise 'lower' or 'higher' in carcinoma by the sign of the median mean difference. Order within group: "
              "median V, descending (as Figure 6). The CSV adds per-cohort p and q, discovery cohorts, median mean "
              "difference and CMS-template membership.",
              "Gene symbols were mapped to current HGNC symbols before the cross-cohort comparison (Supplementary Methods "
              "S7). All 117 genes of the original analysis are retained (IGJ is now JCHAIN, SEPP1 is now SELENOP); five "
              "genes were added because their symbol had differed between platforms (Symbol note). The original 117-gene "
              "results are given in Supplementary Table S_symbol_harmonisation."):
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(7.5)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    t = build()
    t.to_csv(os.path.join(OUT, "TableS_122_genes.csv"), index=False)
    word(t, os.path.join(OUT, "TableS_122_genes.docx"))
    print(t.group.value_counts().to_dict())
    print("cohorts meeting p<0.05 & V>1:", t.n_cohorts_p05_V_gt1.value_counts().sort_index().to_dict())
    print("cohorts tested:", t.n_cohorts_tested.value_counts().sort_index().to_dict())
    ex = ["FCGBP", "AGR3", "EREG", "MOSPD1"]
    print(t.set_index("gene").loc[ex, [f"V_{c}" for c in GEO] + [f"p_{c}" for c in GEO] +
                                  ["n_cohorts_tested", "n_cohorts_p05_V_gt1", "n_cohorts_q05_V_gt1"]].round(3).to_string())


if __name__ == "__main__":
    main()
