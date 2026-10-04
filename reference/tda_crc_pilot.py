#!/usr/bin/env python3
"""
tda_crc_pilot.py  (v2, 25 Sep 2026)
Does the TOPOLOGY of the colorectal transcriptome change across the ordered
progression Normal -> Adenoma -> Carcinoma -> Metastasis, and at which step?

Method (Vietoris-Rips persistent homology, GUDHI backend)
  1. Load GSE41258 series matrix (GPL96), map samples to 4 ordered stages.
  2. log2 if needed, drop the lowest-expressed 25% of probes, keep the top-N
     variable probes, z-score per probe, ONE shared PCA embedding (all samples).
  3. Each stage is represented by repeated equal-size subsamples (without
     replacement), so unequal group sizes cannot create a topological difference.
  4. Vietoris-Rips persistence (H0, H1; 2-skeleton) on each subsample, twice:
       raw        : Euclidean distances in PC space
       shape-only : distances divided by the subsample's median pairwise distance
                    (removes global spread; keeps shape)
     VR persistence depends only on pairwise distances, so it is invariant to
     translation and rotation: a shift of the stage centroid alone cannot
     produce a topological difference.
  5. Test statistic for stages s, t (identical for observed and permuted labels):
       T = mean_{a,b} W1(D_s^a, D_t^b) - mean within-stage W1 over all pairs
     W1 = exact 1-Wasserstein distance between persistence diagrams (L-inf
     ground metric, diagonal matching), summed over the chosen homology dims.
     p-value from label permutations (group sizes preserved) after centring
     each group (and, shape-only, scaling it to unit median distance) so that
     H0 = 'same shape up to translation/scale' is exchangeable;
     p = (1 + #{T_perm >= T_obs}) / (1 + n_perm); Holm correction across the
     adjacent transitions. z = (T_obs - mean T_perm) / sd T_perm.
  6. Descriptive effect sizes (TSI, Betti curves, total persistence) come from a
     separate larger bootstrap and are used only for figures.

Changes from v1: observed and null statistics now use the identical procedure;
groups aligned before permuting (see align_groups);
exhaustive within/between pairs; Holm correction; default 199 permutations;
optional patient-disjoint comparisons (--one-per-patient); versions and log2
status recorded in summary.txt.

Usage
  python tda_crc_pilot.py --matrix GSE41258_series_matrix.txt.gz --out results_gse41258
  python tda_crc_pilot.py --expr expr.csv --labels labels.csv --out results_x
  python tda_crc_pilot.py --demo --out results_demo      # planted transition
  Final paper settings: add --n-perm 999

Requirements
  pip install gudhi numpy scipy pandas scikit-learn matplotlib
"""

import argparse
import gzip
import io
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

try:
    import gudhi
except ImportError:
    sys.exit("GUDHI is required:  pip install gudhi")

LOG2_APPLIED = [False]
STAGES = ["Normal", "Adenoma", "Carcinoma", "Metastasis"]
COLORS = {"Normal": "#2a78d6", "Adenoma": "#eb6834",
          "Carcinoma": "#1baf7a", "Metastasis": "#eda100"}

# GSE41258 tissue labels -> progression stage (None = excluded)
STAGE_RULES = [
    (r"cell line", None),
    (r"normal liver|normal lung", None),
    (r"normal colon|normal mucosa", "Normal"),
    (r"microadenoma", "Adenoma"),            # can be excluded with --no-microadenoma
    (r"polyp|adenoma", "Adenoma"),
    (r"primary tumou?r", "Carcinoma"),
    (r"liver metastasis", "Metastasis"),
    (r"lung metastasis", "Metastasis_lung"),  # kept out unless --include-lung-met
]


# --------------------------------------------------------------------------
# 1. Data loading
# --------------------------------------------------------------------------
def read_series_matrix(path):
    """Parse a GEO series matrix into (expr probes x samples, metadata dict)."""
    opener = gzip.open if path.endswith(".gz") else open
    meta, table_lines, in_table = {}, [], False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key = parts[0]
                vals = [v.strip('"') for v in parts[1:]]
                meta.setdefault(key, []).append(vals)
    expr = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    expr.columns = [c.strip('"') for c in expr.columns]
    return expr.apply(pd.to_numeric, errors="coerce"), meta


def stage_from_text(text, keep_micro, keep_lung):
    t = text.lower()
    if not keep_micro and "microadenoma" in t:
        return None
    for pattern, stage in STAGE_RULES:
        if re.search(pattern, t):
            if stage == "Metastasis_lung":
                return "Metastasis" if keep_lung else None
            return stage
    return None


def labels_from_meta(meta, samples, keep_micro, keep_lung):
    """Find the metadata row that best explains tissue type, map to stages."""
    candidates = []
    for key in ("!Sample_source_name_ch1", "!Sample_characteristics_ch1",
                "!Sample_title", "!Sample_description"):
        for row in meta.get(key, []):
            if len(row) == len(samples):
                candidates.append((key, row))
    best, best_hits = None, -1
    for key, row in candidates:
        hits = sum(stage_from_text(v, True, True) is not None for v in row)
        if hits > best_hits:
            best, best_hits = (key, row), hits
    if best is None or best_hits == 0:
        sys.exit("Could not find tissue labels in the series matrix. "
                 "Use --labels with a sample,stage CSV instead.")
    key, row = best
    print(f"[labels] using metadata field {key} ({best_hits}/{len(samples)} matched)")
    raw = pd.Series(row, index=samples)
    stages = raw.map(lambda v: stage_from_text(v, keep_micro, keep_lung))
    return stages, raw


def make_demo(seed=0, n_genes=2000):
    """Synthetic cohort, GSE41258-like sizes, planted change at Adenoma->Carcinoma.
    Normal & Adenoma: compact Gaussian clouds (no loops).
    Carcinoma & Metastasis: samples spread around a ring in a latent 2D plane
    (heterogeneous subtypes -> an H1 loop) plus noise."""
    rng = np.random.default_rng(seed)
    sizes = {"Normal": 54, "Adenoma": 49, "Carcinoma": 186, "Metastasis": 47}
    loadings = rng.normal(size=(6, n_genes))
    blocks, labels = [], []
    for s, n in sizes.items():
        z = rng.normal(scale=0.6, size=(n, 6))
        if s in ("Carcinoma", "Metastasis"):
            theta = rng.uniform(0, 2 * np.pi, n)
            r = 3.0 if s == "Carcinoma" else 3.3
            z[:, 0] = r * np.cos(theta) + rng.normal(scale=0.35, size=n)
            z[:, 1] = r * np.sin(theta) + rng.normal(scale=0.35, size=n)
            z[:, 2] += 2.0
        elif s == "Adenoma":
            z[:, 2] += 0.8
        x = z @ loadings + rng.normal(scale=2.0, size=(n, n_genes)) + 8
        blocks.append(x)
        labels += [s] * n
    X = np.vstack(blocks).T
    samples = [f"S{i:03d}" for i in range(X.shape[1])]
    expr = pd.DataFrame(X, index=[f"g{i}" for i in range(n_genes)], columns=samples)
    return expr, pd.Series(labels, index=samples)


# --------------------------------------------------------------------------
# 2. Preprocessing and embedding
# --------------------------------------------------------------------------
def preprocess(expr, n_top, n_pcs):
    X = expr.dropna(axis=0, how="any")
    if np.nanpercentile(X.values, 99) > 50:          # not yet log scale
        X = np.log2(X.clip(lower=1))
        LOG2_APPLIED[0] = True
        print("[prep] applied log2 transform")
    mean = X.mean(axis=1)
    X = X[mean > mean.quantile(0.25)]                  # drop lowest-expressed 25%
    var = X.var(axis=1).sort_values(ascending=False)
    X = X.loc[var.index[:n_top]]
    Z = (X.sub(X.mean(axis=1), axis=0)).div(X.std(axis=1) + 1e-9, axis=0)
    pca = PCA(n_components=n_pcs, random_state=0)
    emb = pca.fit_transform(Z.T.values)
    print(f"[prep] {X.shape[0]} probes, {n_pcs} PCs explain "
          f"{pca.explained_variance_ratio_.sum() * 100:.1f}% variance")
    return pd.DataFrame(emb, index=Z.columns), pca, X.shape[0]


# --------------------------------------------------------------------------
# 3. Persistent homology and diagram distances
# --------------------------------------------------------------------------
def persistence(points, normalize):
    D = squareform(pdist(points))
    if normalize:
        D = D / np.median(D[np.triu_indices_from(D, 1)])
    rc = gudhi.RipsComplex(distance_matrix=D)
    st = rc.create_simplex_tree(max_dimension=2)
    st.compute_persistence()
    out = {}
    for dim in (0, 1):
        d = st.persistence_intervals_in_dimension(dim)
        d = d[np.isfinite(d[:, 1])] if len(d) else np.zeros((0, 2))
        out[dim] = d
    return out


def wasserstein(A, B):
    """Exact 1-Wasserstein distance between persistence diagrams (L-inf ground)."""
    n, m = len(A), len(B)
    if n == 0 and m == 0:
        return 0.0
    C = np.zeros((n + m, n + m))
    if n and m:
        C[:n, :m] = np.maximum(np.abs(A[:, None, 0] - B[None, :, 0]),
                               np.abs(A[:, None, 1] - B[None, :, 1]))
    diagA = (A[:, 1] - A[:, 0]) / 2 if n else np.zeros(0)
    diagB = (B[:, 1] - B[:, 0]) / 2 if m else np.zeros(0)
    C[:n, m:] = diagA[:, None]
    C[n:, :m] = diagB[None, :]
    r, c = linear_sum_assignment(C)
    return float(C[r, c].sum())


def summaries(dg):
    out = {}
    for dim in (0, 1):
        life = dg[dim][:, 1] - dg[dim][:, 0] if len(dg[dim]) else np.zeros(0)
        tot = life.sum()
        p = life / tot if tot > 0 else np.zeros(0)
        out[f"H{dim}_total_persistence"] = float(tot)
        out[f"H{dim}_max_persistence"] = float(life.max()) if len(life) else 0.0
        out[f"H{dim}_entropy"] = float(-(p[p > 0] * np.log(p[p > 0])).sum())
        out[f"H{dim}_n_features"] = int(len(life))
    return out


def betti_curve(d, grid):
    if len(d) == 0:
        return np.zeros_like(grid)
    return ((d[:, 0][:, None] <= grid) & (d[:, 1][:, None] > grid)).sum(axis=0)


def bootstrap_diagrams(emb, idx, m, B, normalize, rng):
    diags = []
    for _ in range(B):
        take = rng.choice(idx, size=m, replace=False)
        diags.append(persistence(emb[take], normalize))
    return diags


def diag_dist(d1, d2, dims):
    return sum(wasserstein(d1[k], d2[k]) for k in dims)


def tsi(diags_s, diags_t, dims, rng, n_pairs):
    """Topological separation index = between - within (mean W distances)."""
    B = min(len(diags_s), len(diags_t))
    between, within = [], []
    for _ in range(n_pairs):
        i, j = rng.integers(B, size=2)
        between.append(diag_dist(diags_s[i], diags_t[j], dims))
        a, b = rng.choice(B, size=2, replace=False)
        within.append(diag_dist(diags_s[a], diags_s[b], dims))
        a, b = rng.choice(B, size=2, replace=False)
        within.append(diag_dist(diags_t[a], diags_t[b], dims))
    between, within = np.array(between), np.array(within)
    return between.mean() - within.mean(), between, within


def test_statistic(emb, idx_s, idx_t, m, k, dims, normalize, rng):
    """Test statistic T used for BOTH the observed labels and every permutation.

    Draw k subsamples of size m (without replacement) from each group, compute
    their persistence diagrams, then
        T = mean over all k*k between-group pairs of W(D_s^a, D_t^b)
          - mean over all within-group pairs (k(k-1)/2 per group) of W(D^a, D^b).
    Using exhaustive pairs (not random pairs) removes one layer of Monte Carlo
    noise; using the identical procedure for observed and null labels keeps the
    permutation test exact up to the Monte Carlo randomisation."""
    ds = bootstrap_diagrams(emb, idx_s, m, k, normalize, rng)
    dt = bootstrap_diagrams(emb, idx_t, m, k, normalize, rng)
    between = [diag_dist(a, b, dims) for a in ds for b in dt]
    within = [diag_dist(ds[a], ds[b], dims) for a in range(k) for b in range(a + 1, k)]
    within += [diag_dist(dt[a], dt[b], dims) for a in range(k) for b in range(a + 1, k)]
    return float(np.mean(between) - np.mean(within))


def align_groups(emb, idx_s, idx_t, normalize):
    """Remove the nuisance differences the test must NOT detect, before pooling.

    VR persistence ignores translation, and the shape-only statistic ignores
    global scale. But a permutation test assumes the two groups are
    exchangeable under H0; two clouds that differ only by a shift (or, for
    shape-only, a scale) are NOT exchangeable, because a permuted group becomes
    a two-cloud mixture with a different shape. So each group is centred on its
    own centroid (and, in shape-only mode, divided by its own median pairwise
    distance). This leaves every observed diagram unchanged, and makes H0
    'same shape up to translation (and scale)' exchangeable."""
    out = emb.copy()
    for idx in (idx_s, idx_t):
        out[idx] = emb[idx] - emb[idx].mean(0)
        if normalize:
            out[idx] /= np.median(pdist(out[idx]))
    return out


def permutation_test(emb, idx_s, idx_t, m, k, dims, normalize, n_perm, rng, align=True):
    """Label-permutation test of H0: the two stages share one distribution of
    point clouds up to translation (and scale in shape-only mode).
    Group sizes are preserved."""
    if align:
        emb = align_groups(emb, idx_s, idx_t, normalize)
    observed = test_statistic(emb, idx_s, idx_t, m, k, dims, normalize, rng)
    pool = np.concatenate([idx_s, idx_t])
    ns = len(idx_s)
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(pool)
        null[i] = test_statistic(emb, perm[:ns], perm[ns:], m, k, dims, normalize, rng)
    p = (1 + np.sum(null >= observed)) / (1 + n_perm)
    z = (observed - null.mean()) / (null.std(ddof=1) + 1e-12)
    return observed, float(p), float(z), null


def holm(pvals):
    """Holm step-down adjustment (controls family-wise error rate)."""
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(p) - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


# --------------------------------------------------------------------------
# Patient identity (GSE41258 contains several samples from the same patient)
# --------------------------------------------------------------------------
PATIENT_RE = re.compile(r"\b([A-Z]?\d{4,5})[A-Z]{1,2}R\d+\b")


def patients_from_meta(meta, samples):
    """Return a Series sample -> patient id, parsed from sample codes such as
    '02679BR1' (patient 02679). Returns None if the codes are not found."""
    best, best_hits = None, 0
    for key in ("!Sample_title", "!Sample_source_name_ch1",
                "!Sample_characteristics_ch1", "!Sample_description"):
        for row in meta.get(key, []):
            if len(row) != len(samples):
                continue
            ids = [PATIENT_RE.search(v) for v in row]
            hits = sum(x is not None for x in ids)
            if hits > best_hits:
                best_hits = hits
                best = [x.group(1) if x else None for x in ids]
    if best is None or best_hits < 0.5 * len(samples):
        return None
    return pd.Series(best, index=samples)


def drop_shared_patients(idx_s, idx_t, patient, rng):
    """Make two groups patient-disjoint: for every patient present in both,
    keep that patient's samples in one randomly chosen group only."""
    ps, pt = patient[idx_s], patient[idx_t]
    shared = set(p for p in ps if p is not None) & set(p for p in pt if p is not None)
    keep_s = np.ones(len(idx_s), bool)
    keep_t = np.ones(len(idx_t), bool)
    for p in shared:
        if rng.random() < 0.5:
            keep_s &= ps != p
        else:
            keep_t &= pt != p
    return idx_s[keep_s], idx_t[keep_t], len(shared)


# --------------------------------------------------------------------------
# 4. Figure
# --------------------------------------------------------------------------
def make_figure(emb_df, stages, diags, adj, pairmat, stage_summ, grid, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid_c = "#1f1f1e", "#6b6a63", "#e6e5df"
    plt.rcParams.update({
        "font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
        "xtick.color": muted, "ytick.color": muted, "axes.titleweight": "bold",
        "axes.titlesize": 10, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": grid_c, "grid.linewidth": 0.6,
        "axes.axisbelow": True, "legend.frameon": False,
    })
    present = [s for s in STAGES if s in diags["norm"]]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))

    # A. shared PCA embedding
    ax = axes[0, 0]
    for s in present:
        pts = emb_df.loc[stages[stages == s].index]
        ax.scatter(pts[0], pts[1], s=16, color=COLORS[s], edgecolor="white",
                   linewidth=0.6, label=f"{s} (n={len(pts)})")
    ax.set(xlabel="PC1", ylabel="PC2", title="A  Shared PCA embedding")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=8)

    # B. representative H1 diagrams (shape-only)
    ax = axes[0, 1]
    hi = 0
    for s in present:
        d = diags["norm"][s][0][1]
        if len(d):
            ax.scatter(d[:, 0], d[:, 1], s=18, color=COLORS[s], edgecolor="white",
                       linewidth=0.6, label=s)
            hi = max(hi, d.max())
    hi = hi * 1.05 if hi > 0 else 1
    ax.plot([0, hi], [0, hi], color=muted, lw=1)
    ax.set(xlim=(0, hi), ylim=(0, hi), xlabel="Birth", ylabel="Death",
           title="B  H1 persistence diagrams (one draw, shape-only)")
    ax.legend(loc="lower right", fontsize=8)

    # C. Betti-1 curves with 95% band
    ax = axes[0, 2]
    for s in present:
        curves = np.array([betti_curve(d[1], grid) for d in diags["norm"][s]])
        mu, lo, up = curves.mean(0), np.percentile(curves, 25, 0), np.percentile(curves, 75, 0)
        ax.fill_between(grid, lo, up, color=COLORS[s], alpha=0.12, linewidth=0)
        ax.plot(grid, mu, color=COLORS[s], lw=2, label=s)
    ax.set(xlabel="Filtration scale (x median distance)", ylabel="Betti-1 (loops)",
           title="C  Betti-1 curves, mean and IQR band")
    ax.legend(fontsize=8)

    # D. headline: adjacent transitions as z-scores vs permutation null
    ax = axes[1, 0]
    labels = [f"{r['from'][:4]}→{r['to'][:4]}" for r in adj]
    x = np.arange(len(adj))
    use_z = all("z_norm" in r for r in adj)
    for mode, off, col, hatch, name in [("raw", -0.2, "#2a78d6", None, "raw distances"),
                                         ("norm", 0.2, "#86b6ef", "//", "shape-only (scale-normalized)")]:
        key = f"z_{mode}" if use_z else f"tsi_{mode}"
        vals = [r[key] for r in adj]
        ax.bar(x + off, vals, width=0.36, color=col, edgecolor="white", linewidth=1.5,
               hatch=hatch, label=name)
        for xi, r, v in zip(x, adj, vals):
            p = r.get(f"p_{mode}_holm", r.get(f"p_{mode}"))
            if p is not None:
                ax.text(xi + off, max(v, 0), f"p={p:.3f}", ha="center", va="bottom",
                        fontsize=7, color=ink)
    ax.axhline(0, color=muted, lw=1)
    if use_z:
        ax.axhline(1.96, color=muted, lw=1, ls="--")
        ax.text(-0.45, 1.96, "z = 1.96", fontsize=7, color=muted, va="bottom", ha="left")
    ax.set_xticks(x, labels)
    ax.set(ylabel="z-score vs permutation null\n(labels: Holm-adjusted p)" if use_z else "TSI (between − within)",
           title="D  Stage-to-stage topological jump")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.08)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)

    # E. all-pairs TSI heatmap (shape-only)
    ax = axes[1, 1]
    M = pairmat.loc[present, present].values.astype(float)
    vmax = np.nanmax(np.abs(M)) or 1
    im = ax.imshow(M, cmap="Blues", vmin=0, vmax=vmax)
    ax.set_xticks(range(len(present)), present, rotation=30)
    ax.set_yticks(range(len(present)), present)
    ax.grid(False)
    for i in range(len(present)):
        for j in range(len(present)):
            if i != j:
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if M[i, j] > 0.6 * vmax else ink)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="TSI (shape-only)")
    ax.set_title("E  All stage pairs")

    # F. per-stage total H1 persistence
    ax = axes[1, 2]
    for i, s in enumerate(present):
        v = stage_summ[(stage_summ.stage == s) & (stage_summ["mode"] == "norm")]["H1_total_persistence"]
        ax.boxplot(v, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor=COLORS[s], edgecolor=ink, linewidth=0.8),
                   medianprops=dict(color=ink, lw=1.5),
                   whiskerprops=dict(color=ink), capprops=dict(color=ink))
    ax.set_xticks(range(len(present)), present)
    ax.set(ylabel="Total H1 persistence", title="F  Loop content per stage (shape-only)")

    fig.tight_layout()
    fig.savefig(os.path.join(out, "tda_pilot_figure.png"), dpi=200)
    fig.savefig(os.path.join(out, "tda_pilot_figure.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------------
# 5. Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--matrix", help="GEO series matrix (.txt or .txt.gz)")
    src.add_argument("--expr", help="genes x samples CSV (first column = gene id)")
    src.add_argument("--demo", action="store_true", help="synthetic sanity check")
    ap.add_argument("--labels", help="CSV with columns sample,stage (for --expr)")
    ap.add_argument("--out", default="tda_pilot_results")
    ap.add_argument("--n-top", type=int, default=2000, help="top variable probes")
    ap.add_argument("--n-pcs", type=int, default=10)
    ap.add_argument("--boot", type=int, default=100, help="bootstrap draws per stage")
    ap.add_argument("--subsample", type=int, default=0,
                    help="points per draw (default 80%% of smallest stage)")
    ap.add_argument("--n-perm", type=int, default=199,
                    help="permutations per transition (use 999 for the final paper)")
    ap.add_argument("--k-draws", "--perm-boot", dest="k_draws", type=int, default=20,
                    help="subsample draws per group inside the test statistic")
    ap.add_argument("--dims", default="0,1", help="homology dims to compare, e.g. 1 or 0,1")
    ap.add_argument("--no-microadenoma", action="store_true")
    ap.add_argument("--include-lung-met", action="store_true")
    ap.add_argument("--no-align", action="store_true",
                    help="do not centre/scale groups before permuting (not recommended)")
    ap.add_argument("--one-per-patient", action="store_true",
                    help="make each compared pair of stages patient-disjoint")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    t0 = time.time()
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)
    dims = tuple(int(d) for d in args.dims.split(","))

    # ---- load
    patient = None
    if args.demo:
        expr, stages = make_demo(args.seed)
        raw_labels = stages.copy()
    elif args.matrix:
        expr, meta = read_series_matrix(args.matrix)
        stages, raw_labels = labels_from_meta(meta, list(expr.columns),
                                              not args.no_microadenoma,
                                              args.include_lung_met)
        patient = patients_from_meta(meta, list(expr.columns))
    else:
        expr = pd.read_csv(args.expr, index_col=0)
        if not args.labels:
            sys.exit("--expr needs --labels")
        lab = pd.read_csv(args.labels)
        stages = pd.Series(lab.iloc[:, 1].values, index=lab.iloc[:, 0].astype(str))
        raw_labels = stages.copy()
    stages = stages.dropna()
    stages = stages[stages.isin(STAGES)]
    expr = expr[[c for c in expr.columns if c in stages.index]]
    stages = stages.loc[expr.columns]
    counts = stages.value_counts().reindex(STAGES).dropna().astype(int)
    print("[labels] samples per stage:\n" + counts.to_string())
    pd.DataFrame({"raw_label": raw_labels.reindex(stages.index), "stage": stages,
                  "patient": patient.reindex(stages.index) if patient is not None else None}) \
        .to_csv(os.path.join(args.out, "sample_stages.csv"))
    pat_arr = patient.reindex(stages.index).values if patient is not None else None
    if pat_arr is not None:
        print(f"[labels] patient ids parsed for {pd.notna(pat_arr).sum()} samples, "
              f"{pd.Series(pat_arr).nunique()} patients")
    elif args.one_per_patient:
        print("[labels] WARNING: patient ids not found; --one-per-patient ignored")
    present = list(counts.index)
    if len(present) < 2:
        sys.exit("Need at least two stages.")

    # ---- embed
    emb_df, pca, n_probes = preprocess(expr, args.n_top, args.n_pcs)
    emb = emb_df.values
    pos = {s: np.where(stages.values == s)[0] for s in present}
    m = args.subsample or max(10, int(0.8 * counts.min()))
    print(f"[tda] {args.boot} draws of {m} samples per stage, dims={dims}")

    # ---- bootstrap diagrams
    diags = {"raw": {}, "norm": {}}
    rows = []
    for mode in ("raw", "norm"):
        for s in present:
            diags[mode][s] = bootstrap_diagrams(emb, pos[s], m, args.boot, mode == "norm", rng)
            for b, d in enumerate(diags[mode][s]):
                rows.append({"stage": s, "mode": mode, "draw": b, **summaries(d)})
    stage_summ = pd.DataFrame(rows)
    stage_summ.to_csv(os.path.join(args.out, "stage_topology_bootstrap.csv"), index=False)

    # ---- adjacent transitions: descriptive TSI + permutation tests
    adj = []
    for s, t in zip(present[:-1], present[1:]):
        rec = {"from": s, "to": t}
        idx_s, idx_t = pos[s], pos[t]
        if pat_arr is not None:
            shared = (set(x for x in pat_arr[idx_s] if pd.notna(x)) &
                      set(x for x in pat_arr[idx_t] if pd.notna(x)))
            rec["shared_patients"] = len(shared)
            if args.one_per_patient:
                idx_s, idx_t, _ = drop_shared_patients(idx_s, idx_t, pat_arr, rng)
        rec["n_from"], rec["n_to"] = int(len(idx_s)), int(len(idx_t))
        m_pair = min(m, int(0.8 * min(len(idx_s), len(idx_t))))
        rec["subsample_size"] = m_pair
        for mode in ("raw", "norm"):
            # descriptive effect size from the large bootstrap (figures)
            obs_d, btw, wth = tsi(diags[mode][s], diags[mode][t], dims, rng, 400)
            boots = []
            for _ in range(500):
                bi = rng.choice(len(btw), len(btw))
                wi = rng.choice(len(wth), len(wth))
                boots.append(btw[bi].mean() - wth[wi].mean())
            rec[f"tsi_{mode}"] = float(obs_d)
            rec[f"tsi_{mode}_ci"] = [float(np.percentile(boots, 2.5)),
                                     float(np.percentile(boots, 97.5))]
            rec[f"ratio_{mode}"] = float(btw.mean() / wth.mean())
            # inferential test: identical statistic for observed and permuted labels
            if args.n_perm > 0:
                T, p, z, _ = permutation_test(emb, idx_s, idx_t, m_pair, args.k_draws,
                                              dims, mode == "norm", args.n_perm, rng,
                                              align=not args.no_align)
                rec[f"T_{mode}"], rec[f"p_{mode}"], rec[f"z_{mode}"] = T, p, z
        cs, ct = emb[idx_s].mean(0), emb[idx_t].mean(0)
        pooled = np.sqrt((emb[idx_s].var(0).sum() + emb[idx_t].var(0).sum()) / 2)
        rec["centroid_shift_sd"] = float(np.linalg.norm(cs - ct) / pooled)
        adj.append(rec)
        print(f"[tda] {s:>10} -> {t:<10} z_shape={rec.get('z_norm', float('nan')):.1f} "
              f"p_shape={rec.get('p_norm', float('nan')):.3f} "
              f"z_raw={rec.get('z_raw', float('nan')):.1f}")
    if args.n_perm > 0:
        for mode in ("raw", "norm"):
            for r, pa in zip(adj, holm([r[f"p_{mode}"] for r in adj])):
                r[f"p_{mode}_holm"] = float(pa)

    # ---- all pairs (shape-only, no permutation)
    pairmat = pd.DataFrame(np.nan, index=present, columns=present)
    for i, s in enumerate(present):
        for t in present[i + 1:]:
            v = tsi(diags["norm"][s], diags["norm"][t], dims, rng, 200)[0]
            pairmat.loc[s, t] = pairmat.loc[t, s] = v
    pairmat.to_csv(os.path.join(args.out, "all_pairs_tsi_shape_only.csv"))

    # ---- figure
    allmax = max(d[1][:, 1].max() for s in present for d in diags["norm"][s] if len(d[1])) \
        if any(len(d[1]) for s in present for d in diags["norm"][s]) else 1.0
    grid = np.linspace(0, allmax * 1.05, 200)
    make_figure(emb_df, stages, diags, adj, pairmat, stage_summ, grid, args.out)

    # ---- verdict (uses Holm-adjusted shape-only p-values)
    names = [f"{r['from']}->{r['to']}" for r in adj]
    key = "Adenoma->Carcinoma"
    nan = float("nan")
    lines = ["TDA pilot: topological change across CRC progression (script v2)", "=" * 64,
             f"Samples per stage: {counts.to_dict()}",
             f"Embedding: {n_probes} probes, {args.n_pcs} PCs "
             f"({pca.explained_variance_ratio_.sum() * 100:.1f}% variance), log2 applied: {LOG2_APPLIED[0]}",
             f"Homology dims compared: {dims}; descriptive bootstrap: {args.boot} draws x {m}",
             f"Test: {args.n_perm} label permutations, statistic from {args.k_draws} draws per group; "
             f"Holm correction over {len(adj)} transitions; groups aligned before permuting: {not args.no_align}",
             f"Patient-disjoint comparisons: {bool(args.one_per_patient and pat_arr is not None)}",
             f"Versions: gudhi {gudhi.__version__}, numpy {np.__version__}, "
             f"pandas {pd.__version__}; seed {args.seed}", ""]
    for r in adj:
        lines.append(
            f"{r['from']:>10} -> {r['to']:<10} n={r['n_from']}/{r['n_to']}"
            + (f", shared patients={r['shared_patients']}" if "shared_patients" in r else ""))
        lines.append(
            f"    shape-only: z={r.get('z_norm', nan):.1f}, p={r.get('p_norm', nan):.4f}, "
            f"p_Holm={r.get('p_norm_holm', nan):.4f} | TSI {r['tsi_norm']:.3f} "
            f"[{r['tsi_norm_ci'][0]:.3f}, {r['tsi_norm_ci'][1]:.3f}]")
        lines.append(
            f"    raw       : z={r.get('z_raw', nan):.1f}, p={r.get('p_raw', nan):.4f}, "
            f"p_Holm={r.get('p_raw_holm', nan):.4f} | centroid shift {r['centroid_shift_sd']:.2f} SD")
    lines.append("")
    if key in names and args.n_perm > 0:
        r = adj[names.index(key)]
        sig_list = [x for x in adj if x["p_norm_holm"] < 0.05]
        largest = max(adj, key=lambda x: x["z_norm"])
        if r["p_norm_holm"] < 0.05 and largest is r:
            verdict = ("GO: shape-only topology changes significantly at Adenoma->Carcinoma "
                       "and it is the largest adjacent jump.")
        elif r["p_norm_holm"] < 0.05:
            verdict = (f"PARTIAL: Adenoma->Carcinoma is significant but the largest jump is "
                       f"{largest['from']}->{largest['to']}.")
        elif largest in sig_list:
            verdict = (f"REDIRECT: no significant change at Adenoma->Carcinoma; the shape-only "
                       f"topology jumps significantly at {largest['from']}->{largest['to']} "
                       f"(z={largest['z_norm']:.1f}, Holm p={largest['p_norm_holm']:.4f}).")
        elif r["p_raw_holm"] < 0.05:
            verdict = ("CAUTION: only the raw (spread-sensitive) topology differs at "
                       "Adenoma->Carcinoma; the change is dispersion rather than shape.")
        else:
            verdict = "NO-GO: no significant shape-only topological change at any step."
        lines.append("Verdict: " + verdict)
    lines += ["", "Sensitivity checks: --dims 1 | --n-pcs 5 / 20 | --n-top 1000 / 5000 | "
              "--no-microadenoma | --one-per-patient | --seed", "",
              f"Runtime: {time.time() - t0:.0f} s"]
    text = "\n".join(lines)
    print("\n" + text)
    with open(os.path.join(args.out, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    with open(os.path.join(args.out, "transitions.json"), "w") as fh:
        json.dump(adj, fh, indent=2)


if __name__ == "__main__":
    main()
