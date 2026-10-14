"""
Addendum 5, part 3b (commit 5c5791a) - level 2 diagnostic with bug check (descriptive; no tests).

Exact reconstruction of the first replicate of each of the 5 saved level-2 shape-only chunks
(same per-unit seeds as part 3), for cluster a = 1, loop a = 1, spread a = 0.5 and a = 1.
(a) insertion present in the matrix passed to preprocess? B-A mean difference and SD ratio before/after;
    changed entries confined to the intended 200 genes x 48 B samples (intended genes reconstructed
    independently from a copy of the random generator);
(b) altered genes in the top-2,000 variance selection (after; before for reference);
(c,d) share of the inserted z-space signal inside the first 10 / 50 PCs (exact SVD of the after data).
"""
import copy
import os
import sys
import zlib

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import addendum5_part3_controls as P  # noqa: E402

OUT = os.path.join(ROOT, "results", "addendum5", "part3b_diagnostic")
os.makedirs(OUT, exist_ok=True)
CONDS = [("CL", 1.0, "cluster split a=1"), ("LP", 1.0, "loop a=1"), ("SP", 0.5, "spread +50% (a=0.5)"),
         ("SP", 1.0, "spread +100% (a=1, largest)")]


def selection(E, n_top=2000):
    """Same gene selection as tda_crc_pilot.preprocess (genes x samples numpy array)."""
    mean = E.mean(axis=1)
    keep = np.where(mean > np.quantile(mean, 0.25))[0]          # pandas quantile = linear, same as numpy
    var = E[keep].var(axis=1, ddof=1)
    return keep[np.argsort(-var, kind="stable")[:n_top]]


emb_df = pd.read_csv(P.EMB_FILE, index_col=0)
carc = emb_df.index[emb_df.stage == "Carcinoma"]
expr = pd.read_csv(P.EXPR_FILE, index_col=0)[carc]
rows = []
for cond, a, label in CONDS:
    for ch in range(5):
        name = P.unit_name(2, cond, a, "norm", ch)
        rng = np.random.default_rng(np.random.SeedSequence(P.SEED, spawn_key=(zlib.crc32(name.encode()),)))
        perm = rng.permutation(len(carc))
        ia, ib = perm[:P.NA_], perm[P.NA_:P.NA_ + P.NB_]
        Eb = expr.iloc[:, np.r_[ia, ib]].values.astype(float).copy()
        colsB = np.arange(P.NA_, P.NA_ + P.NB_)
        rcopy = copy.deepcopy(rng)                                   # independent reconstruction of genes
        mean0 = Eb.mean(axis=1)
        intended = np.sort(rcopy.choice(np.where(mean0 > np.quantile(mean0, 0.25))[0], size=200, replace=False))
        Ea = P.insert_level2(Eb.copy(), colsB, cond, a, rng)        # the matrix passed to preprocess
        diff = Ea - Eb
        ch_rows, ch_cols = np.nonzero(np.abs(diff) > 0)
        rows_ok = np.array_equal(np.unique(ch_rows), intended)
        cols_ok = set(np.unique(ch_cols)) <= set(colsB)
        g = intended
        mdiff_b = np.median(Eb[g][:, colsB].mean(1) - Eb[g][:, :P.NA_].mean(1))
        mdiff_a = np.median(Ea[g][:, colsB].mean(1) - Ea[g][:, :P.NA_].mean(1))
        sdr_b = np.median(Eb[g][:, colsB].std(1, ddof=1) / Eb[g][:, :P.NA_].std(1, ddof=1))
        sdr_a = np.median(Ea[g][:, colsB].std(1, ddof=1) / Ea[g][:, :P.NA_].std(1, ddof=1))
        sel_a, sel_b = selection(Ea), selection(Eb)
        in_a, in_b = np.isin(g, sel_a).sum(), np.isin(g, sel_b).sum()
        # z-space (after-insertion means / SDs), selected genes
        Xa = Ea[sel_a]
        mu, sd = Xa.mean(1, keepdims=True), Xa.std(1, ddof=1, keepdims=True) + 1e-9
        Za = ((Xa - mu) / sd).T                                     # samples x genes
        S = diff[sel_a]
        Sz = ((S - S.mean(1, keepdims=True)) / sd).T
        U, sv, Vt = np.linalg.svd(Za, full_matrices=False)
        proj = Sz @ Vt.T                                             # signal coordinates on each PC
        e_pc = (proj ** 2).sum(0)
        e_tot = (Sz ** 2).sum()
        f10, f50 = e_pc[:10].sum() / e_tot, e_pc[:50].sum() / e_tot
        top = int(np.argmax(e_pc)) + 1
        var_share = sv ** 2 / (sv ** 2).sum()
        # energy of the inserted signal over all 200 genes (z-scaled by each gene's after-insertion SD)
        sd_all = Ea[g].std(1, ddof=1, keepdims=True) + 1e-9
        Sg = (diff[g] - diff[g].mean(1, keepdims=True)) / sd_all
        e200 = (Sg ** 2).sum(1)
        kept = e200[np.isin(g, sel_a)].sum() / e200.sum()
        rows.append(dict(insertion=label, repeat=ch + 1, unit=name, changed_entries=int(len(ch_rows)),
                         changed_rows_equal_intended_200=bool(rows_ok), changed_cols_only_B=bool(cols_ok),
                         meandiff_BA_before=mdiff_b, meandiff_BA_after=mdiff_a, sdratio_BA_before=sdr_b,
                         sdratio_BA_after=sdr_a, genes_selected_after=int(in_a), genes_selected_before=int(in_b),
                         energy_kept_by_selection=kept, signal_share_of_total_var=e_tot / (Za ** 2).sum(),
                         frac_in_10PCs=f10, frac_in_50PCs=f50, top_PC=top, top_PC_signal_share=e_pc[top - 1] / e_tot,
                         top_PC_var_share=var_share[top - 1], retained_overall_10PCs=kept * f10,
                         retained_overall_50PCs=kept * f50))
        print(f"{label:<28} rep {ch + 1}: entries {len(ch_rows)}, rows ok {rows_ok}, cols ok {cols_ok}, "
              f"SD ratio {sdr_b:.2f}->{sdr_a:.2f}, selected {in_a}/200 (before {in_b}), "
              f"10PC {f10:.2f}, 50PC {f50:.2f}, top PC {top}", flush=True)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "diagnostic_per_repeat.csv"), index=False)


def med(x):
    return f"{np.median(x):.2f} [{np.min(x):.2f}-{np.max(x):.2f}]"


L = ["ADDENDUM 5, PART 3b - Level 2 diagnostic with bug check (descriptive; no permutation tests)", "=" * 72,
     "Pre-specified: addendum 5 part 3b (commit 5c5791a). Exact reconstructions of the tested level-2 data",
     "(first replicate of each of the 5 saved shape-only chunks, same seeds). Values: median [min-max] over 5 repeats.", ""]
for label in df.insertion.unique():
    d = df[df.insertion == label]
    L += [f"--- {label} ---",
          f"  (a) changed entries {med(d.changed_entries)}; changed rows = intended 200 genes: "
          f"{int(d.changed_rows_equal_intended_200.sum())}/5; changes only in B samples: {int(d.changed_cols_only_B.sum())}/5",
          f"      B-A mean difference (median over 200 genes): before {med(d.meandiff_BA_before)}, after {med(d.meandiff_BA_after)}",
          f"      B/A SD ratio (median over 200 genes):        before {med(d.sdratio_BA_before)}, after {med(d.sdratio_BA_after)}",
          f"  (b) altered genes in top 2,000: after {med(d.genes_selected_after)} (before {med(d.genes_selected_before)})"
          f"; share of inserted energy kept by selection {med(d.energy_kept_by_selection)}",
          f"  (c) inserted signal = {med(d.signal_share_of_total_var * 100)}% of total z-space variance; "
          f"share inside 10 PCs {med(d.frac_in_10PCs)}; top PC {', '.join(map(str, d.top_PC))} "
          f"(signal share {med(d.top_PC_signal_share)}, PC variance share {med(d.top_PC_var_share)})",
          f"  (d) share inside 50 PCs {med(d.frac_in_50PCs)}",
          f"      overall retained (selection x PCs): 10 PCs {med(d.retained_overall_10PCs)}, 50 PCs {med(d.retained_overall_50PCs)}", ""]
L += ["Outputs: results/addendum5/part3b_diagnostic/diagnostic_per_repeat.csv"]
open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
