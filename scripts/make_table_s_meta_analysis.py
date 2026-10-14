"""
Supplementary Table S_meta_analysis: every prespecified secondary meta-analysis (addendum 3; addendum 6 section 8),
copied from the saved results; no statistical test is run here.
  Model (both addenda): DerSimonian-Laird random effects of the per-cohort SHAPE-ONLY permutation z-scores (within-
  cohort variance 1), separately for unadjusted, stage-blind adjusted and stage-aware data. Raw-mode z-scores were
  not part of the prespecified meta-analysis and were not pooled.
  (A) addendum 3: normal -> adenoma, three cohorts (GSE41258, GSE20916, GSE37364)
  (B) addendum 6 section 8: normal -> adenoma and adenoma -> carcinoma, four cohorts (+ GSE117606)
Output: outputs/supplementary/TableS_meta_analysis.{csv,docx}

Usage: python scripts/make_table_s_meta_analysis.py
"""
import os

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "outputs", "supplementary")
COH = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
VLAB = {"unadjusted": "unadjusted", "adjusted": "adjusted (stage-blind)", "stageaware": "stage-aware"}
TLAB = {"Normal->Adenoma": "N → A", "Adenoma->Carcinoma": "A → C"}


def rows(part, analysis, res, inp, cohorts):
    out = []
    for _, r in res.iterrows():
        tr = r.get("transition", "Normal->Adenoma")
        x = inp[(inp.version == r.version) & ((inp.transition == tr) if "transition" in inp else True)].set_index("cohort")
        rec = dict(part=part, analysis=analysis, transition=TLAB[tr], data=VLAB[r.version], mode="shape-only",
                   n_cohorts=len(x), pooled_z=r.pooled_z, ci_low=r.ci_low, ci_high=r.ci_upp, p=r.p_two_sided,
                   tau2=r.tau2, Q=r.Q, Q_p=r.Q_p, I2_percent=100 * r.I2,
                   ci_excludes_0="yes" if (r.ci_low > 0 or r.ci_upp < 0) else "no")
        for c in cohorts:
            rec[f"z_{c}"] = float(x.loc[c, "z"]) if c in x.index else np.nan
        out.append(rec)
    return out


def word(t, path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Cm(1.3))
    doc.add_paragraph().add_run("Supplementary Table S_meta_analysis. Prespecified secondary random-effects meta-analyses "
                                "of the shape-only step 5 z-scores (informed summaries).").bold = True
    cols = ["transition", "data", "n_cohorts", "pooled_z", "ci", "p", "tau2", "Q (p)", "I2_percent"] + [f"z_{c}" for c in COH]
    heads = ["Transition", "Data", "Cohorts", "Pooled z", "95% CI", "p", "τ²", "Q (p)", "I² (%)"] + [f"z {c}" for c in COH]
    for part, title in (("A", "A. Addendum 3: normal → adenoma, three cohorts"),
                        ("B", "B. Addendum 6 §8: normal → adenoma and adenoma → carcinoma, four cohorts")):
        doc.add_paragraph().add_run(title).bold = True
        sub = t[t.part == part]
        tab = doc.add_table(rows=1, cols=len(cols))
        tab.style = "Table Grid"
        for j, h in enumerate(heads):
            run = tab.rows[0].cells[j].paragraphs[0].add_run(h)
            run.bold, run.font.size = True, Pt(7.5)
        for _, r in sub.iterrows():
            cells = tab.add_row().cells
            vals = {"transition": r.transition, "data": r.data, "n_cohorts": str(r.n_cohorts), "pooled_z": f"{r.pooled_z:.2f}",
                    "ci": f"{r.ci_low:.2f} to {r.ci_high:.2f}", "p": f"{r.p:.3g}", "tau2": f"{r.tau2:.2f}",
                    "Q (p)": f"{r.Q:.2f} ({r.Q_p:.2g})", "I2_percent": f"{r.I2_percent:.0f}",
                    **{f"z_{c}": ("–" if pd.isna(r[f"z_{c}"]) else f"{r[f'z_{c}']:.2f}") for c in COH}}
            for j, c in enumerate(cols):
                run = cells[j].paragraphs[0].add_run(vals[c])
                run.font.size = Pt(7.5)
                run.bold = c in ("pooled_z", "ci", "p") and r.ci_excludes_0 == "yes"
    for n in ("Model (prespecified in addenda 3 and 6): DerSimonian–Laird random effects of the per-cohort shape-only "
              "step 5 permutation z-scores, z = (T_obs − mean of the null) / SD of the null, each with within-cohort "
              "variance 1; two-sided p from the normal distribution; τ², Cochran's Q and I² describe between-cohort "
              "heterogeneity. Raw-mode z-scores were not part of the prespecified meta-analysis and were not pooled.",
              "Informed summaries: most input cohorts had been analysed before these meta-analyses were specified, so they "
              "are not independent confirmation, and none defines a tipping point (the tipping-point criterion is "
              "applied per cohort; Table 2). GSE41258 stage-aware z = exploratory analysis D3. With three or four "
              "cohorts, τ² and I² are imprecise.",
              "Bold = 95% CI excluding 0. Interpretation fixed in addendum 3: a pooled CI excluding 0 in the adjusted data "
              "is reported as consistent evidence of a composition-robust shape difference across cohorts, not as a "
              "tipping point.",
              "No new tests were run; values are copied from results/meta_analysis/ and results/meta_analysis_4cohort/."):
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    a_res = pd.read_csv(os.path.join(R, "meta_analysis", "meta_results.csv"))
    a_inp = pd.read_csv(os.path.join(R, "meta_analysis", "meta_inputs.csv"))
    b_res = pd.read_csv(os.path.join(R, "meta_analysis_4cohort", "meta_results.csv"))
    b_inp = pd.read_csv(os.path.join(R, "meta_analysis_4cohort", "meta_inputs.csv"))
    t = pd.DataFrame(rows("A", "addendum 3 (three cohorts)", a_res, a_inp, COH[:3]) +
                     rows("B", "addendum 6 section 8 (four cohorts)", b_res, b_inp, COH))
    order = {"unadjusted": 0, "adjusted (stage-blind)": 1, "stage-aware": 2}
    t["_o"] = t.data.map(order)
    t = t.sort_values(["part", "transition", "_o"], ascending=[True, False, True]).drop(columns="_o").reset_index(drop=True)
    for c in COH:
        if f"z_{c}" not in t:
            t[f"z_{c}"] = np.nan
    t.to_csv(os.path.join(OUT, "TableS_meta_analysis.csv"), index=False)
    word(t, os.path.join(OUT, "TableS_meta_analysis.docx"))
    pd.set_option("display.width", 220)
    print(t[["part", "transition", "data", "n_cohorts", "pooled_z", "ci_low", "ci_high", "p", "I2_percent", "ci_excludes_0"]]
          .round(4).to_string(index=False))


if __name__ == "__main__":
    main()
