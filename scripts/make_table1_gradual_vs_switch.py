"""
Table 1 (main; Table 3 in earlier drafts): addendum 8 (gradual vs switch-like change along normal -> adenoma -> carcinoma;
Amendments 1 and 2), copied from the saved results in results/addendum8/; no statistical test is run here.
  Part A, per cohort x transition: stage-order check (adenoma median a* on the all-gene score, part a, and on the two
  cross-fitted half scores, part b), sample continuity (AUC, valley ratio with 95% CI, Holm dip p, call), gene classes
  (share of progression-associated genes gradual and switch-like), candidate switch genes (FDR) in the interval,
  interpretability (Amendment 2), number replicated in >= 2 other cohorts, gene-level call.
  Part B, per transition: the section 5.3 counts and the prespecified conclusion.
  Counts are recounted from the saved gene_fits_<cohort>.csv (cls, interval, cand_switch, q_assoc) and checked against
  results/addendum8/summary.txt; the gene-level calls (stored only in summary.txt) are read from it.
Output: outputs/main_tables/Table1_gradual_vs_switch.{csv,docx}

Usage: python scripts/make_table1_gradual_vs_switch.py
"""
import os
import re

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "addendum8")
OUT = os.path.join(ROOT, "outputs", "main_tables")
COH = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
TR = [("NA", "N->A", "N → A", "Normal", "Adenoma"), ("AC", "A->C", "A → C", "Adenoma", "Carcinoma")]


def summary_calls():
    """Gene-level call, candidate count and interpretability per cohort x transition, and class percentages, from summary.txt."""
    calls, pct, coh = {}, {}, None
    for line in open(os.path.join(R, "summary.txt"), encoding="utf-8"):
        m = re.match(r"  (GSE\d+)\s+half-score", line)
        if m:
            coh = m.group(1)
            continue
        if coh and re.match(r"\s+gradual [\d.]+%", line):
            pct[coh] = {k.strip(): float(v) for k, v in re.findall(r"([a-z -]+) ([\d.]+)%", line)}
            continue
        m = re.match(r"\s+(N->A|A->C): candidate switch genes \(FDR\) (\d+) \(([\d.]+)%\); interpretable: (True|False); "
                     r"gene-level call: (\w+)", line)
        if m and coh:
            calls[(coh, m.group(1))] = dict(n=int(m.group(2)), interpretable=m.group(4) == "True", call=m.group(5).lower())
        if line.startswith("REPLICATED SWITCH GENES"):
            coh = None
    return calls, pct


def word(t, c, path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Cm(1.2))
    doc.add_paragraph().add_run("Table 1. Gradual versus switch-like change along normal → "
                                "adenoma → carcinoma (addendum 8; stage-blind adjusted data; informed analysis).").bold = True
    heads = ["Cohort", "Trans-ition", "n", "Stage order: a* (part a); half scores (part b)", "AUC", "Valley ratio (95% CI)",
             "Dip Holm p (pair)", "Sample continuity", "Progression-assoc. genes P", "Gradual (% of P)",
             "Switch-like (% of P)", "Candidate switch genes, this interval (% of P)", "Switch calls interpretable",
             "Replicated in ≥ 2 other cohorts", "Gene-level call"]
    doc.add_paragraph().add_run("A. Per cohort and transition").bold = True
    tab = doc.add_table(rows=1, cols=len(heads))
    tab.style = "Table Grid"
    for j, h in enumerate(heads):
        run = tab.rows[0].cells[j].paragraphs[0].add_run(h)
        run.bold, run.font.size = True, Pt(7)
    for _, r in t.iterrows():
        vr = "not computed (AUC < 0.70)" if pd.isna(r.valley_ratio) else f"{r.valley_ratio:.2f} ({r.vr_ci_low:.2f}–{r.vr_ci_high:.2f})"
        vals = [r.cohort, r.transition, f"{r.n_from} / {r.n_to}",
                f"{r.astar_part_a:.2f}; {r.astar_half1:.2f} / {r.astar_half2:.2f} ({'ordered' if r.ordered_a and r.ordered_b else 'not ordered'})",
                f"{r.auc:.2f}", vr, f"{r.dip_holm_p_pair:.2f}", r.continuity_call, str(r.P), f"{r.gradual_pct:.1f}",
                f"{r.switch_like_pct:.1f}", f"{r.candidate_switch_interval} ({r.candidate_switch_pct_of_P:.1f})",
                "yes" if r.switch_calls_interpretable else "no (possible gap artefact)",
                str(r.replicated_switch_genes), r.gene_level_call]
        cells = tab.add_row().cells
        for j, v in enumerate(vals):
            cells[j].paragraphs[0].add_run(v).font.size = Pt(7)
    doc.add_paragraph().add_run("B. Prespecified conclusion per transition (section 5.3; ordered cohorts only)").bold = True
    h2 = ["Transition", "Cohorts ordered (part a / part b)", "Continuous or overlapping", "Discontinuous",
          "Gene-level gradual", "Replicated switch genes", "Conclusion"]
    tab = doc.add_table(rows=1, cols=len(h2))
    tab.style = "Table Grid"
    for j, h in enumerate(h2):
        run = tab.rows[0].cells[j].paragraphs[0].add_run(h)
        run.bold, run.font.size = True, Pt(8)
    for _, r in c.iterrows():
        vals = [r.transition, f"{r.ordered_a} / {r.ordered_b} of 4", f"{r.pop_continuous_or_overlapping} of 4",
                f"{r.pop_discontinuous} of 4", f"{r.gene_gradual} of 4", str(r.replicated_switch_genes), r.conclusion]
        cells = tab.add_row().cells
        for j, v in enumerate(vals):
            run = cells[j].paragraphs[0].add_run(v)
            run.font.size = Pt(8)
            run.bold = j == 6
    notes = [
        "Progression score: adjusted PC1 scaled so that the normal median is 0 and the carcinoma median 1; a* = adenoma "
        "median, which must lie strictly between 0 and 1 (ordering check). Part (a) uses the all-gene PC1; part (b) uses two "
        "cross-fitted half scores (each gene fitted against the score built without it).",
        "Sample continuity (section 5.1): continuous = VR ≥ 0.5 and VR lower 95% CI ≥ 0.25 and dip test for the pooled pair "
        "not significant (Holm p ≥ 0.05); overlapping = AUC < 0.70; discontinuous = dip test for the pooled pair significant "
        "(Holm p < 0.05) and VR upper 95% CI < 0.25; inconclusive = anything else.",
        "Gene level (section 5.2; P = progression-associated genes, spline F-test BH q < 0.05): switch-like = candidate "
        "switch genes (FDR) with midpoint in this interval ≥ 10% of P; gradual = gradual genes ≥ 50% of P and candidate switch "
        "genes in this interval < 2% of P; inconclusive = anything else. Gene classes are per cohort (the same in both rows).",
        "Switch-like gene (section 4.2, Amendment 1): ΔAIC ≥ 10, W ≤ 0.25, 0 ≤ x₀ ≤ 1, amplitude ≥ 1 log2, ≥ 10% of samples "
        "(and ≥ 5) within x₀ ± 0.125, not non-monotonic. Candidate = switch-like and bootstrap BH q < 0.05. Gradual gene: "
        "progression-associated, not non-monotonic, and ΔAIC < 2 or W ≥ 1.",
        "Amendment 2: switch-gene calls for a transition count only in cohorts where that transition is continuous or "
        "overlapping; elsewhere they are 'not interpretable (possible gap artefact)'. Replication (section 4.4): the same "
        "switch criteria at nominal p < 0.05, same direction and interval, in ≥ 2 other cohorts.",
        "For A→C, replication of switch genes was not possible by design: under Amendment 2, a transition that is "
        "inconclusive or discontinuous in ≥ 2 of 4 cohorts cannot yield replicated switch genes. The A→C count of 0 is "
        "therefore not independent evidence against a switch.",
        "Conclusion (section 5.3): support for gradual change = population level continuous or overlapping in ≥ 3 cohorts "
        "and gene level gradual in ≥ 3 cohorts; support for switch-like change = discontinuous in ≥ 3 cohorts and ≥ 10 "
        "replicated switch genes; otherwise inconclusive (or the two mixed categories, which require ≥ 10 replicated switch "
        "genes). All conclusions are worded 'consistent with'.",
        "Informed analysis: written after the step 5 and step 13 results were known, before any addendum 8 computation. "
        "Cross-sectional data: the score orders samples, not time. No new tests were run; values are copied from "
        "results/addendum8/.",
    ]
    for n in notes:
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(7.5)
    doc.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    cont = pd.read_csv(os.path.join(R, "continuity_results.csv")).set_index("cohort")
    half = pd.read_csv(os.path.join(R, "half_scores.csv")).set_index("cohort")
    rep = pd.read_csv(os.path.join(R, "switch_genes_replicated.csv"))
    concl = pd.read_csv(os.path.join(R, "section5_conclusions.csv"))
    calls, pct = summary_calls()
    rows = []
    for g in COH:
        f = pd.read_csv(os.path.join(R, f"gene_fits_{g}.csv"))
        assoc = f[f.cls != "not associated"]
        P = len(assoc)
        assert P == int((f.q_assoc < 0.05).sum())
        share = (assoc.cls.value_counts() / P * 100)
        for key, tr, lab, x, y in TR:
            ca = cont.loc[g]
            n_cand = int((f.cand_switch & (f.interval == tr)).sum())
            s = calls[(g, tr)]
            assert s["n"] == n_cand, (g, tr, s["n"], n_cand)
            for cls, v in share.items():
                if cls in pct[g]:
                    assert abs(round(v, 1) - pct[g][cls]) < 1e-9, (g, cls, v, pct[g][cls])
            reps = rep[(rep.discovery == g) & (rep.interval == tr)]
            rows.append(dict(
                cohort=g, transition=lab, n_from=int(ca[f"{key}_n_x"]), n_to=int(ca[f"{key}_n_y"]),
                astar_part_a=ca.astar, ordered_a=bool(ca.ordered), astar_half1=half.loc[g, "H1_astar"],
                astar_half2=half.loc[g, "H2_astar"], ordered_b=bool(half.loc[g, "ordered_b"]),
                auc=ca[f"{key}_AUC"], valley_ratio=ca[f"{key}_VR"], vr_ci_low=ca[f"{key}_VR_lo"], vr_ci_high=ca[f"{key}_VR_hi"],
                dip_holm_p_pair=ca[f"dip_holm_{key}"], continuity_call=ca[f"{key}_call"],
                n_between_medians_from=int(ca[f"{key}_n_x_between"]), n_between_medians_to=int(ca[f"{key}_n_y_between"]),
                P=P, gradual_pct=share.get("gradual", 0.0), intermediate_pct=share.get("intermediate", 0.0),
                switch_like_pct=share.get("switch-like", 0.0), step_like_unresolvable_pct=share.get("step-like unresolvable", 0.0),
                non_monotonic_pct=share.get("non-monotonic", 0.0),
                switch_like_interval=int(((f.cls == "switch-like") & (f.interval == tr)).sum()),
                candidate_switch_interval=n_cand, candidate_switch_pct_of_P=100 * n_cand / P,
                switch_calls_interpretable=s["interpretable"],
                replicated_switch_genes=int(reps.replicated.sum()) if len(reps) else 0,
                max_other_cohorts_replicating=int(reps.n_replicating.max()) if len(reps) else 0,
                gene_level_call=s["call"],
                conclusion_transition=concl.set_index("transition").loc[tr, "conclusion"]))
    t = pd.DataFrame(rows)
    t["_o"] = t.transition.map({"N → A": 0, "A → C": 1})
    t = t.sort_values(["_o"], kind="stable").drop(columns="_o").reset_index(drop=True)
    t.to_csv(os.path.join(OUT, "Table1_gradual_vs_switch.csv"), index=False)
    c = concl.copy()
    c["transition"] = c.transition.map({"N->A": "N → A", "A->C": "A → C"})
    word(t, c, os.path.join(OUT, "Table1_gradual_vs_switch.docx"))
    pd.set_option("display.width", 250)
    print(t[["cohort", "transition", "astar_part_a", "auc", "valley_ratio", "dip_holm_p_pair", "continuity_call", "P",
             "gradual_pct", "switch_like_pct", "candidate_switch_interval", "switch_calls_interpretable",
             "replicated_switch_genes", "gene_level_call"]].round(3).to_string(index=False))
    print(c.to_string(index=False))


if __name__ == "__main__":
    main()
