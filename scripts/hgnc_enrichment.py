"""
Addendum 9 pathway enrichment and CMS-template overlap, rerun on the harmonised 122-gene list (PI decision 2026-10-15;
part of the symbol-harmonisation deviation, Supplementary Methods S7). The original 117-gene results in
results/addendum9/ are not changed.

Method exactly as addendum9_heterogeneity.py lines 473-491:
  - background = genes tested (present in the variability analysis) in >= 3 of the 4 GEO cohorts;
  - enrichment = addendum8_gradual_switch.enrichment, unchanged: MSigDB Hallmark + Reactome via decoupler/OmniPath,
    each set intersected with the background and tested only if it then has 10-500 genes; hypergeometric
    over-representation; BH FDR across the tested sets;
  - CMS-template overlap = hypergeometric test of (list n background) vs (CMS template genes n background).
Inputs
  original (117):   results/addendum9/variability_<cohort>.csv, variability_replicated.csv; templates as stored
  harmonised (122): results/hgnc/addendum9/variability_<cohort>.csv, variability_replicated.csv; template symbols
                    harmonised exactly as hgnc_crosscohort.annotate (renamed where the symbol has one current match)
Check first: the original run must reproduce results/addendum9/enrichment_variability.csv (set names, set sizes, k, p,
q) and the summary line "30 of 117 (template genes in background 380; p 8.09e-18)"; gene sets are downloaded from
OmniPath at run time, so this also checks that the resource has not changed since addendum 9 was run.
Outputs: results/hgnc/addendum9/enrichment_variability.csv, enrichment_comparison.csv, enrichment_summary.txt

Usage: python scripts/hgnc_enrichment.py
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import check_methods_items as CK  # noqa: E402

A9 = os.path.join(ROOT, "results", "addendum9")
H9 = os.path.join(ROOT, "results", "hgnc", "addendum9")
COH9 = ["GSE41258", "GSE20916", "GSE37364", "GSE117606"]
LIST = "more variable in carcinoma"


def load_a8():
    spec = importlib.util.spec_from_file_location("a8", os.path.join(ROOT, "scripts", "addendum8_gradual_switch.py"))
    a8 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a8)
    return a8


def run(folder, tpl, a8):
    tested = {c: set(pd.read_csv(os.path.join(folder, f"variability_{c}.csv"), index_col=0).index) for c in COH9}
    bg = [x for x in set().union(*tested.values()) if sum(x in tested[o] for o in COH9) >= 3]
    uniq = pd.read_csv(os.path.join(folder, "variability_replicated.csv"))
    lists = {d: list(v.gene) for d, v in uniq.groupby("direction")}
    enr, msg = a8.enrichment(lists, bg)
    mv = set(uniq[uniq.direction == LIST].gene) & set(bg)
    T_ = tpl & set(bg)
    k = len(mv & T_)
    pt = hypergeom.sf(k - 1, len(bg), len(T_), len(mv)) if k else 1.0
    e = enr[enr.list == LIST]
    return enr, dict(genes_in_list=int((uniq.direction == LIST).sum()), list_in_background=len(mv), background=len(bg),
                     gene_sets_loaded=int(msg.split()[0]), gene_sets_tested=len(e), smallest_p=float(e.p.min()),
                     smallest_q=float(e.q.min()), sets_q_below_0_05=int((e.q < 0.05).sum()),
                     sets_p_below_0_05=int((e.p < 0.05).sum()), top_set=e.sort_values("p").gene_set.iloc[0],
                     cms_overlap_k=k, cms_template_genes_in_background=len(T_), cms_overlap_p=float(pt))


def main():
    a8 = load_a8()
    tpl_tab = pd.read_csv(os.path.join(A9, "cache", "templates.csv"))
    tpl_orig = set(tpl_tab.symbol)
    CUR, ALIAS = CK.hgnc_map()
    cls = CK.classify(tpl_tab.symbol.astype(str), CUR, ALIAS)
    tpl_h = {cls[s][1] if cls[s][0] == "renamed" else s for s in tpl_tab.symbol.astype(str)}

    # check: original inputs reproduce the stored addendum 9 enrichment and template overlap
    enr_o, so = run(A9, tpl_orig, a8)
    ref = pd.read_csv(os.path.join(A9, "enrichment_variability.csv"))
    key = ["list", "gene_set"]
    m = ref.merge(enr_o, on=key, how="outer", suffixes=("_ref", "_new"), indicator=True)
    same_sets = bool((m._merge == "both").all())
    same_vals = same_sets and all(np.allclose(m[f"{c}_ref"], m[f"{c}_new"], rtol=0, atol=1e-12)
                                  for c in ("k", "list_size", "set_size", "p", "q"))
    tpl_ok = (so["cms_overlap_k"], so["list_in_background"], so["cms_template_genes_in_background"]) == (30, 117, 380) \
        and f"{so['cms_overlap_p']:.3g}" == "8.09e-18" and so["background"] == 10467
    print(f"check original: sets identical {same_sets}; values identical {same_vals}; template line and background "
          f"identical {tpl_ok}")
    if not (same_vals and tpl_ok):
        print("only in stored:", m[m._merge == "left_only"].gene_set.tolist()[:10])
        print("only in rerun:", m[m._merge == "right_only"].gene_set.tolist()[:10])
        raise SystemExit("original enrichment not reproduced (gene-set resource or code differs); stopping")

    enr_h, sh = run(H9, tpl_h, a8)
    enr_h.to_csv(os.path.join(H9, "enrichment_variability.csv"), index=False)
    comp = pd.DataFrame([dict(analysis="original symbols (addendum 9)", **so),
                         dict(analysis="harmonised symbols (S7)", **sh)])
    comp.to_csv(os.path.join(H9, "enrichment_comparison.csv"), index=False)
    top = enr_h[enr_h.list == LIST].sort_values("p").head(8)
    show = comp.astype(object)
    for c in ("smallest_p", "cms_overlap_p"):
        show[c] = [f"{v:.3g}" for v in comp[c]]
    show["smallest_q"] = [f"{v:.3f}" for v in comp.smallest_q]
    L = ["ADDENDUM 9 ENRICHMENT AND CMS-TEMPLATE OVERLAP, HARMONISED SYMBOLS (scripts/hgnc_enrichment.py)", "=" * 78,
         "Method unchanged (addendum9_heterogeneity.py lines 473-491; addendum8_gradual_switch.enrichment).",
         f"Check: original inputs reproduce results/addendum9/enrichment_variability.csv exactly ({so['gene_sets_tested']} "
         "rows) and the template line 30 of 117 / 380 / p 8.09e-18 -> gene-set resource unchanged",
         "", show.T.to_string(header=False), "",
         "Harmonised, 8 smallest p:", top[["gene_set", "k", "set_size", "p", "q", "genes"]].to_string(index=False),
         "",
         f"Gene sets: {sh['gene_sets_loaded']} loaded (Hallmark + Reactome); a set is tested only if 10-500 of its genes are "
         f"in the background, so {sh['gene_sets_tested']} (harmonised) / {so['gene_sets_tested']} (original) were tested, "
         "and BH is applied across the tested sets. The CSV lists every tested set (including k = 0, p = 1)."]
    open(os.path.join(H9, "enrichment_summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
