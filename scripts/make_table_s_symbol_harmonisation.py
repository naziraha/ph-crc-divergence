"""
Supplementary Table S_symbol_harmonisation: original (platform symbols, 117 genes) versus harmonised (current HGNC
symbols, 122 genes) cross-cohort results, copied from results/addendum9, results/addendum11 (original) and
results/hgnc (harmonised; scripts/hgnc_crosscohort.py). No statistical tests.
  A. variable genes: every gene in either list, with symbol, inclusion, group and median V in both analyses
  B. arm variability index (median ln V) ranks for 20q, 20p, 13q, 8q and comparators 7p, 7q
  C. per-sample arm scores: Brown-Forsythe p, V_arm and median difference (GEO); 20q gain-like share; TCGA rho
  D. share of excess carcinoma dispersion explained (f, f_rand95): 20q alone and four arms
  E. Figure 6 labelled genes
Output: outputs/supplementary/TableS_symbol_harmonisation.{csv,docx}

Usage: python scripts/make_table_s_symbol_harmonisation.py
"""
import os
import sys

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import check_methods_items as CK  # noqa: E402

R = os.path.join(ROOT, "results")
O9, N9 = os.path.join(R, "addendum9"), os.path.join(R, "hgnc", "addendum9")
O11, N11 = os.path.join(R, "addendum11"), os.path.join(R, "hgnc", "addendum11")
OUT = os.path.join(ROOT, "outputs", "supplementary")
GEO = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
ARMS = ["20q", "20p", "13q", "8q", "7p", "7q"]
CUR, AL = CK.hgnc_map()


def cur(s):
    k = CK.classify([s], CUR, AL)[s]
    return k[1] if k[0] == "renamed" else s


def groups(ann):
    return pd.Series(np.where(ann.no_mean_shift, "no mean shift",
                              np.where(ann.med_mean_diff < 0, "lower in carcinoma", "higher in carcinoma")), index=ann.gene)


def labset(ann):
    g = groups(ann)
    a = ann.assign(group=g.values).sort_values("median_V", ascending=False)
    return (list(a[a.group == "higher in carcinoma"].gene) + list(a[a.group == "lower in carcinoma"].gene[:5])
            + list(a[a.group == "no mean shift"].gene[:5]))


def part_a():
    o = pd.read_csv(os.path.join(O9, "variability_replicated_annotated.csv"))
    n = pd.read_csv(os.path.join(N9, "variability_replicated_annotated.csv"))
    go, gn = groups(o), groups(n)
    o["current"] = [cur(x) for x in o.gene]
    om = o.set_index("current")
    nm = n.set_index("gene")
    rows = []
    for g in sorted(set(om.index) | set(nm.index)):
        rows.append(dict(section="A. Variable genes", gene_current=g,
                         symbol_original=om.loc[g, "gene"] if g in om.index else "–",
                         in_original_117="yes" if g in om.index else "no", in_harmonised_122="yes" if g in nm.index else "no",
                         group_original=go[om.loc[g, "gene"]] if g in om.index else "–", group_harmonised=gn[g] if g in nm.index else "–",
                         median_V_original=float(om.loc[g, "median_V"]) if g in om.index else np.nan,
                         median_V_harmonised=float(nm.loc[g, "median_V"]) if g in nm.index else np.nan))
    return pd.DataFrame(rows), o, n


def part_b():
    o = pd.read_csv(os.path.join(O11, "partA_chr", "arm_variability_ranking.csv"))
    n = pd.read_csv(os.path.join(N11, "partA_chr", "arm_variability_ranking.csv"))
    rows = []
    for a in ARMS:
        rec = dict(section="B. Arm variability index rank", arm=a)
        for lab, d in (("original", o), ("harmonised", n)):
            ranks = [int(d[(d.cohort == g) & (d.arm == a)]["rank"].iloc[0]) for g in GEO]
            for g, r in zip(GEO, ranks):
                rec[f"rank_{g}_{lab}"] = r
            rec[f"cohorts_rank_le5_{lab}"] = int(sum(r <= 5 for r in ranks))
        rows.append(rec)
    return pd.DataFrame(rows)


def part_c():
    o = pd.read_csv(os.path.join(O11, "partA_chr", "arm_scores_tests.csv"))
    n = pd.read_csv(os.path.join(N11, "partA_chr", "arm_scores_tests.csv"))
    o = o[(o.version == "adjusted") & (o.genes == "excl117")]
    n = n[(n.version == "adjusted") & (n.genes == "excl122")]
    rows = []
    for a in ["20q", "20p", "13q", "8q"]:
        for g in GEO:
            x, y = o[(o.cohort == g) & (o.arm == a)].iloc[0], n[(n.cohort == g) & (n.arm == a)].iloc[0]
            rows.append(dict(section="C. Per-sample arm score, carcinoma vs adenoma", arm=a, cohort=g,
                             n_genes_original=int(x.n_genes), n_genes_harmonised=int(y.n_genes),
                             V_arm_original=x.V_arm, V_arm_harmonised=y.V_arm, BF_p_original=x.p_BF, BF_p_harmonised=y.p_BF,
                             median_diff_original=x.median_diff, median_diff_harmonised=y.median_diff,
                             gainlike_share_original=x.get("gainlike_share_carc", np.nan),
                             gainlike_share_harmonised=y.get("gainlike_share_carc", np.nan)))
    to = pd.read_csv(os.path.join(O11, "partA_chr", "tcga_expression_vs_copynumber.csv")).set_index("arm")
    tn = pd.read_csv(os.path.join(N11, "partA_chr", "tcga_expression_vs_copynumber.csv")).set_index("arm")
    for a in ["20q", "20p", "13q", "8q"]:
        rows.append(dict(section="C. TCGA copy-number check", arm=a, cohort="TCGA-COAD",
                         n_genes_original=int(to.loc[a, "n_expr_genes"]), n_genes_harmonised=int(tn.loc[a, "n_expr_genes"]),
                         rho_original=to.loc[a, "rho"], rho_harmonised=tn.loc[a, "rho"],
                         n_tumours_original=int(to.loc[a, "n_tumours"]), n_tumours_harmonised=int(tn.loc[a, "n_tumours"])))
    return pd.DataFrame(rows)


def part_d():
    o = pd.read_csv(os.path.join(O11, "partB_differentiation", "B3_dispersion_decomposition.csv")).set_index("cohort")
    n = pd.read_csv(os.path.join(N11, "partB_differentiation", "B3_dispersion_decomposition.csv")).set_index("cohort")
    rows = []
    for name, lab in (("chr20q", "20q alone"), ("chr4", "four arms (20q, 20p, 13q, 8q)")):
        for g in GEO:
            rows.append(dict(section="D. Dispersion share f", score=lab, cohort=g, informative="yes" if o.loc[g, "informative"] else "no",
                             f_original=o.loc[g, f"{name}_f"], f_rand95_original=o.loc[g, f"{name}_f_rand95"],
                             f_harmonised=n.loc[g, f"{name}_f"], f_rand95_harmonised=n.loc[g, f"{name}_f_rand95"],
                             f_gt_rand95_original="yes" if o.loc[g, f"{name}_f"] > o.loc[g, f"{name}_f_rand95"] else "no",
                             f_gt_rand95_harmonised="yes" if n.loc[g, f"{name}_f"] > n.loc[g, f"{name}_f_rand95"] else "no"))
    return pd.DataFrame(rows)


def fmt(v, c):
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return "–"
        return f"{v:.3g}" if "p" in c.split("_") else f"{v:.2f}"
    return str(v)


def add(doc, df, cols, heads, title, bold_if=None):
    doc.add_paragraph().add_run(title).bold = True
    t = doc.add_table(rows=1, cols=len(cols))
    t.style = "Table Grid"
    for j, h in enumerate(heads):
        r = t.rows[0].cells[j].paragraphs[0].add_run(h)
        r.bold, r.font.size = True, Pt(7)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, c in enumerate(cols):
            run = cells[j].paragraphs[0].add_run(fmt(row[c], c))
            run.font.size = Pt(7)
            if bold_if is not None and bold_if(row):
                run.bold = True


def main():
    A, o_ann, n_ann = part_a()
    B, C, D = part_b(), part_c(), part_d()
    lo, ln = labset(o_ann), labset(n_ann)
    E = pd.DataFrame([dict(section="E. Figure 6 labelled genes", analysis="original (117)", genes=", ".join(lo)),
                      dict(section="E. Figure 6 labelled genes", analysis="harmonised (122)", genes=", ".join(ln))])
    allr = pd.concat([A, B, C, D, E], ignore_index=True)
    allr.to_csv(os.path.join(OUT, "TableS_symbol_harmonisation.csv"), index=False)

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Cm(1.3))
    doc.add_paragraph().add_run("Supplementary Table S_symbol_harmonisation. Cross-cohort results with the original platform "
                                "gene symbols (117 variable genes) and with current HGNC symbols (122 variable genes).").bold = True
    ch = A[(A.in_original_117 == "no") | (A.symbol_original != A.gene_current)]
    add(doc, ch, ["gene_current", "symbol_original", "in_original_117", "in_harmonised_122", "group_harmonised", "median_V_harmonised"],
        ["Gene (current)", "Original symbol", "In original 117", "In harmonised 122", "Group", "Median V"],
        f"A. Variable genes that differ between the analyses ({len(ch)} of {len(A)}; the other {len(A) - len(ch)} are identical, "
        "with unchanged median V; full list in the .csv)")
    bc = ["arm"] + [f"rank_{g}_{l}" for l in ("original", "harmonised") for g in GEO] + ["cohorts_rank_le5_original", "cohorts_rank_le5_harmonised"]
    bh = ["Arm"] + [f"{g} {l[:4]}." for l in ("original", "harmonised") for g in GEO] + ["Rank ≤ 5, orig.", "Rank ≤ 5, harm."]
    add(doc, B, bc, bh, "B. Arm variability index (median ln V) rank among 39 arms (1 = largest increase)",
        bold_if=lambda r: r.cohorts_rank_le5_original != r.cohorts_rank_le5_harmonised)
    cg = C[C.cohort != "TCGA-COAD"]
    add(doc, cg, ["arm", "cohort", "n_genes_original", "n_genes_harmonised", "V_arm_original", "V_arm_harmonised", "BF_p_original",
                  "BF_p_harmonised", "median_diff_original", "median_diff_harmonised"],
        ["Arm", "Cohort", "Genes orig.", "Genes harm.", "V_arm orig.", "V_arm harm.", "BF p orig.", "BF p harm.", "Median diff orig.", "Median diff harm."],
        "C1. Per-sample arm score, carcinoma vs adenoma (adjusted data; bold = significance call at p < 0.05 changes)",
        bold_if=lambda r: (r.BF_p_original < 0.05) != (r.BF_p_harmonised < 0.05))
    g20 = C[(C.arm == "20q") & (C.cohort != "TCGA-COAD")]
    add(doc, g20, ["cohort", "gainlike_share_original", "gainlike_share_harmonised"], ["Cohort", "20q gain-like, orig.", "20q gain-like, harm."],
        "C2. Share of carcinomas with a 20q gain-like score (> adenoma 97.5th percentile)")
    tc = C[C.cohort == "TCGA-COAD"]
    add(doc, tc, ["arm", "n_genes_original", "n_genes_harmonised", "n_tumours_harmonised", "rho_original", "rho_harmonised"],
        ["Arm", "Expressed genes orig.", "Expressed genes harm.", "Tumours", "Spearman ρ orig.", "Spearman ρ harm."],
        "C3. TCGA-COAD: per-sample arm score vs mean GISTIC2 copy number")
    add(doc, D, ["score", "cohort", "informative", "f_original", "f_rand95_original", "f_harmonised", "f_rand95_harmonised",
                 "f_gt_rand95_original", "f_gt_rand95_harmonised"],
        ["Score", "Cohort", "Informative", "f orig.", "f_rand95 orig.", "f harm.", "f_rand95 harm.", "f > rand95 orig.", "f > rand95 harm."],
        "D. Share f of excess carcinoma dispersion explained (adjusted 10-PC space, bias-adjusted distances)",
        bold_if=lambda r: r.f_gt_rand95_original != r.f_gt_rand95_harmonised)
    add(doc, E, ["analysis", "genes"], ["Analysis", "Labelled genes (Figure 6)"], "E. Figure 6 labelled genes")
    for n in ("Original = gene symbols as given by each platform annotation (GPL96, GPL570, NCBI gene_info for GSE117606, "
              "Xena for TCGA). Harmonised = every symbol mapped to the current HGNC symbol (NCBI gene_info, nomenclature "
              "authority; previous symbols and aliases with one current match), before the cross-cohort comparison "
              "(Supplementary Methods S7). The primary persistent-homology analysis (Table 2, Figures 2–4, Supplementary "
              "Table S_sensitivity) was not repeated.",
              "No statistical test was run for this table; all values are copied from the two sets of result files. In "
              "GSE41258 f is not interpretable (no significant excess dispersion; not informative)."):
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(os.path.join(OUT, "TableS_symbol_harmonisation.docx"))
    pd.set_option("display.width", 250)
    print(ch.drop(columns="section").to_string(index=False))
    print(B.drop(columns="section").to_string(index=False))
    print(E.to_string(index=False))
    for lab, ann in (("original", o_ann), ("harmonised", n_ann)):
        gr = groups(ann).value_counts().to_dict()
        print(lab, "groups", gr)


if __name__ == "__main__":
    main()
