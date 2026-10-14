"""
Supplementary Table S_power_controls: (A) addendum 5 part 3 positive/negative controls (GSE41258 carcinomas),
(B) addendum 10 TCGA-COAD positive control (D20/D21), (C) step 5 subsample sizes per cohort and transition.
All numbers copied from existing result files; no statistical tests.
Outputs: outputs/supplementary/TableS_power_controls.{csv,docx}

Usage: python scripts/make_table_s_power_controls.py
"""
import os

import pandas as pd
from docx import Document
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "outputs", "supplementary")
COND = {"N0": "negative control (no change)", "SP": "spread", "CL": "cluster split", "LP": "loop"}
L1_VALUE = {"SP": "s = {v} (spread factor)", "CL": "δ = {v} × MPD", "LP": "w = {v} (loop weight)"}
L2_VALUE = "a = {v}"


def part_a():
    t = pd.read_csv(os.path.join(R, "addendum5", "part3_controls", "controls_table.csv"))
    rows = []
    for _, r in t.iterrows():
        lvl = int(r.level)
        if r.condition == "N0":
            size = "–"
        elif lvl == 1:
            size = L1_VALUE[r.condition].format(v=f"{r.value:g}")
        else:
            size = L2_VALUE.format(v=f"{r.value:g}") + (" (spread × (1 + a))" if r.condition == "SP" else " × gene SD")
        rows.append({"section": "A. Addendum 5 controls", "cohort": "GSE41258 (181 carcinomas; random split 53 / 48)",
                     "level": "1: embedding (10 PCs)" if lvl == 1 else "2: expression (full pipeline)",
                     "inserted change": COND[r.condition], "size": size,
                     "genes changed": "– (PC scores)" if lvl == 1 else ("0" if r.condition == "N0" else "200"),
                     "mode": r["mode"], "detected (p ≤ 0.05)": f"{int(r.det05)}/{int(r.n)}",
                     "detected (p ≤ 0.0167)": f"{int(r.det_bonf)}/{int(r.n)}", "median z": round(float(r.median_z), 2)})
    return pd.DataFrame(rows)


def part_b():
    t = pd.read_csv(os.path.join(R, "addendum10", "controls", "controls_table.csv"))
    return t


def part_c():
    files = {"GSE41258": "step05/permutation_results_all.csv", "GSE20916": "step05/GSE20916/permutation_results_all.csv",
             "GSE37364": "step05/GSE37364/permutation_results_all.csv",
             "GSE117606": "step05/GSE117606/permutation_results_all.csv", "TCGA-COAD": "step05/TCGA/permutation_results_all.csv"}
    rows = []
    for c, f in files.items():
        d = pd.read_csv(os.path.join(R, f))
        d = d[(d.analysis.str.lower() == "primary") & (d.version == "unadjusted") & (d["mode"] == "raw")]
        for _, r in d.iterrows():
            a, b = r.transition.split("->")
            rows.append({"section": "C. Step 5 subsample sizes", "cohort": c, "transition": f"{a} → {b}",
                         "n (from / to)": f"{int(r.n_from)} / {int(r.n_to)}", "subsample size m per group":
                         int(r.subsample_size), "subsamples k per group": 20, "permutations": 999})
    return pd.DataFrame(rows)


def part_d():
    """D7 (raw mode only; shape-only part withdrawn in D7) and D8 empirical power, copied from results/power_D7, power_D8."""
    d7 = pd.read_csv(os.path.join(R, "power_D7", "power_table.csv"))
    d7 = d7[d7["mode"] == "raw"].sort_values(["cohort", "version", "alpha"], ascending=[True, False, False])
    K = {"GSE41258": 3, "GSE20916": 2, "GSE37364": 2}
    D7 = pd.DataFrame([{"section": "D1. D7 power (raw mode)", "cohort": r.cohort, "version": r.version,
                        "effect level (share of real N→A difference)": f"{int(round(r.alpha * 100))}%",
                        "samples swapped": int(r.swaps), "detected (p ≤ 0.05/K)": f"{int(round(r.power_bonf * r.n))}/{int(r.n)}",
                        "detected (p ≤ 0.05)": f"{int(round(r.power_05 * r.n))}/{int(r.n)}", "K": K[r.cohort],
                        "median z": round(float(r.median_z), 2)} for _, r in d7.iterrows()])
    d8 = pd.read_csv(os.path.join(R, "power_D8", "power_table.csv")).sort_values(["cohort", "frac"], ascending=[False, False])
    K8 = {"GSE41258": 3, "GSE20916": 2}
    D8 = pd.DataFrame([{"section": "D2. D8 sample-size sensitivity (shape-only, unadjusted N→A)", "cohort": r.cohort,
                        "share of samples kept": f"{int(round(r.frac * 100))}%", "n normal / adenoma": f"{int(r.n_N)} / {int(r.n_A)}",
                        "subsample size m": int(r.m), "detected (p ≤ 0.05/K)": f"{int(round(r.power_bonf * r.replicates))}/{int(r.replicates)}",
                        "detected (p ≤ 0.05)": f"{int(round(r.power_05 * r.replicates))}/{int(r.replicates)}", "K": K8[r.cohort],
                        "median z": round(float(r.median_z), 2)} for _, r in d8.iterrows()])
    return D7, D8


def add_table(doc, df, cols, title, fs=7.5):
    doc.add_paragraph().add_run(title).bold = True
    t = doc.add_table(rows=1, cols=len(cols))
    t.style = "Table Grid"
    for j, c in enumerate(cols):
        r = t.rows[0].cells[j].paragraphs[0].add_run(c)
        r.bold, r.font.size = True, Pt(fs)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, c in enumerate(cols):
            cells[j].paragraphs[0].add_run(str(row[c])).font.size = Pt(fs)


def main():
    os.makedirs(OUT, exist_ok=True)
    A, Braw, C = part_a(), part_b(), part_c()
    B = pd.DataFrame([{"section": "B. Addendum 10 TCGA control", "cohort": "TCGA-COAD (110 Stage II tumours; random split 55 / 55)",
                       "level": "1: embedding (10 PCs)", "inserted change": COND[r.condition],
                       "size": {"N0": "–", "SP": "s = 1.25", "CL": "δ = 1 × MPD", "LP": "w = 1"}[r.condition],
                       "genes changed": "– (PC scores)", "mode": r["mode"],
                       "detected (p ≤ 0.05)": f"{int(r.detected)}/{int(r.n)}", "detected (p ≤ 0.0167)": "",
                       "median z": round(float(r.median_z), 2)} for _, r in Braw.iterrows()])
    D7, D8 = part_d()
    allrows = pd.concat([A, B, C, D7, D8], ignore_index=True)
    allrows.to_csv(os.path.join(OUT, "TableS_power_controls.csv"), index=False)

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = "Calibri", Pt(9)
    sec = doc.sections[0]
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(1.5))
    doc.add_paragraph().add_run("Supplementary Table S_power_controls. Positive and negative controls of the step 5 "
                                "test, and subsample sizes.").bold = True
    cA = ["level", "inserted change", "size", "genes changed", "mode", "detected (p ≤ 0.05)", "detected (p ≤ 0.0167)",
          "median z"]
    add_table(doc, A, cA, "A. Addendum 5 controls, GSE41258 carcinomas (random split 53 / 48; 50 replicates per row; "
                          "99 permutations; k = 20, m = 31)")
    add_table(doc, B, ["inserted change", "size", "mode", "detected (p ≤ 0.05)", "median z"],
              "B. Addendum 10 positive control, TCGA-COAD Stage II tumours (random split 55 / 55; 20 replicates per row; "
              "99 permutations; k = 20, m = 31)")
    add_table(doc, C, ["cohort", "transition", "n (from / to)", "subsample size m per group", "subsamples k per group",
                       "permutations"], "C. Step 5 subsample sizes (primary analysis; same m in all modes and versions)")
    add_table(doc, D7, ["cohort", "version", "effect level (share of real N→A difference)", "samples swapped",
                        "detected (p ≤ 0.05/K)", "detected (p ≤ 0.05)", "median z"],
              "D1. Empirical power, deviation D7 (raw mode; real normal → adenoma difference diluted by swapping samples "
              "between groups; 50 replicates; 99 permutations). The shape-only part of D7 was withdrawn (see note).")
    add_table(doc, D8, ["cohort", "share of samples kept", "n normal / adenoma", "subsample size m",
                        "detected (p ≤ 0.05/K)", "detected (p ≤ 0.05)", "median z"],
              "D2. Sample-size sensitivity, deviation D8 (shape-only; real unadjusted normal → adenoma difference with "
              "fewer samples; 50 replicates; 99 permutations)")
    for n in ("Detection = permutation p ≤ 0.05 in a replicate (A also at p ≤ 0.0167, Bonferroni for 3 transitions). "
              "Insertions are made into group B only; group A is unchanged. MPD = median pairwise distance of group B; "
              "s = spread factor about B's centroid; δ = half-separation of the two clusters along a random direction; "
              "w = weight of an annulus (radius 1.5 × median distance to centroid) in a random 2-D plane; level 2: "
              "a = amplitude in units of each gene's SD in B, in 200 random genes passing the expression filter, followed "
              "by the full preprocessing (top 2,000 genes, z-score, 10 PCs).",
              "Pre-set acceptance (A): C1 negative control ≤ 7/50 (both levels, both modes); C2 raw detects spread s = 1.25 "
              "≥ 40/50; C3 shape-only ignores every pure spread ≤ 7/50; C4 shape-only detects cluster split δ = 1.0 ≥ 40/50; "
              "C5 shape-only detects loop w = 1.0 ≥ 40/50; C6 level 2 a = 1.0 cluster and loop, shape-only ≥ 40/50. "
              "C1–C5 passed; C6 failed (3/50 and 5/50 detected): changes carried by 200 genes are largely lost in gene "
              "selection and the 10-PC embedding (part 3b).",
              "Pre-set acceptance (B; D20): N0 ≤ 3/20 both modes; SP raw ≥ 16/20; SP shape-only ≤ 3/20; CL and LP shape-only "
              "≥ 16/20. Five of six passed; shape-only flagged pure spread in 4/20 (limit 3/20; D21).",
              "C: m = 80% of the smallest stage of the cohort (GSE41258: 39 metastases → 31; TCGA: 39 Stage IV → 31), drawn "
              "without replacement k = 20 times per group; how patients contributing to both groups are handled is described "
              "in Table 2.",
              "Both control sets used the unadjusted 10-PC embeddings (A: GSE41258 `main`; B: TCGA-COAD `main`).",
              "D1 (D7): for effect level α, k = round((1 − α) ÷ (1/n_N + 1/n_A)) normals and adenomas swap groups, so the "
              "group contrast is α of the real one; group sizes unchanged, each sample used once; same m as step 5 "
              "(31 / 19 / 21); detection at p ≤ 0.05/K, K = number of transitions (Bonferroni stand-in for Holm). D7 also "
              "withdrew the addendum 4 bootstrap effect-size intervals (invalid: duplicate samples). Shape-only rows are not "
              "shown: swapping creates mixed two-cluster groups, an artificial shape difference, so they are not "
              "interpretable as power (results/power_D7/README_NOTE.txt). All non-significant raw step 5 results were "
              "classified 'too few samples to tell' (detection < 80% at α = 50%). GSE117606 was not part of D7.",
              "D2 (D8): each group randomly reduced to the stated share without replacement; m = min(cohort m, ⌊0.8 × "
              "smaller reduced group⌋); only cohorts with a significant unadjusted shape-only normal → adenoma difference "
              "(GSE41258, GSE20916). It shows how reliably a shape difference of the observed size is detected with fewer "
              "samples, not power for smaller effects."):
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(2)
        q.add_run(n).font.size = Pt(8)
    doc.save(os.path.join(OUT, "TableS_power_controls.docx"))
    pd.set_option("display.width", 220, "display.max_colwidth", 60)
    print(A[["level", "inserted change", "size", "mode", "detected (p ≤ 0.05)", "detected (p ≤ 0.0167)", "median z"]].to_string(index=False))
    print(B[["inserted change", "size", "mode", "detected (p ≤ 0.05)", "median z"]].to_string(index=False))
    print(C[["cohort", "transition", "n (from / to)", "subsample size m per group"]].to_string(index=False))


if __name__ == "__main__":
    main()
