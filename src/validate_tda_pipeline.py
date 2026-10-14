#!/usr/bin/env python3
"""
validate_tda_pipeline.py
Independent checks of the maths and statistics in tda_crc_pilot.py (v2).

Part A - deterministic correctness checks (seconds)
  A1  Our exact 1-Wasserstein distance equals GUDHI's reference implementation
      (gudhi.wasserstein, POT backend, order=1, internal_p=inf).
  A2  Vietoris-Rips on a noisy circle gives exactly one dominant H1 class.
  A3  Two well-separated clusters give exactly one dominant finite H0 bar.
  A4  Diagrams are invariant to translation and rotation of the point cloud.
  A5  Shape-only diagrams are invariant to global scaling; raw diagrams scale.
  A6  Holm adjustment matches statsmodels.

Part B - Monte Carlo calibration of the permutation test (minutes)
  B1  H0, identical distributions            -> rejection rate should be ~5%
  B2  H0, same shape but shifted (translation) -> ~5% only WITH group alignment
  B3  H0 (shape-only), same shape, 2x scale   -> ~5% only WITH group alignment
  B4  H1, loop (ring) vs blob                -> power should be high

Usage:  python validate_tda_pipeline.py            (A + B, about 10 min on 2 cores)
        python validate_tda_pipeline.py --quick    (A only)
"""
import argparse
import sys
import time
from multiprocessing import Pool

import numpy as np

import tda_crc_pilot as T

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")


# ---------------------------------------------------------------- Part A
def part_a():
    import gudhi.wasserstein as gw
    from statsmodels.stats.multitest import multipletests

    rng = np.random.default_rng(1)

    # A1
    worst = 0.0
    for _ in range(200):
        def rand_diag(n):
            b = rng.uniform(0, 1, n)
            return np.c_[b, b + rng.exponential(0.3, n)]
        A, B = rand_diag(rng.integers(0, 25)), rand_diag(rng.integers(0, 25))
        ours = T.wasserstein(A, B)
        ref = gw.wasserstein_distance(A, B, order=1, internal_p=np.inf)
        worst = max(worst, abs(ours - ref))
    check("A1 Wasserstein matches GUDHI reference (200 random pairs)",
          worst < 1e-9, f"max |diff| = {worst:.2e}")

    # A2
    th = rng.uniform(0, 2 * np.pi, 60)
    circle = np.c_[np.cos(th), np.sin(th)] + rng.normal(0, 0.05, (60, 2))
    h1 = T.persistence(circle, False)[1]
    life = np.sort(h1[:, 1] - h1[:, 0])[::-1]
    second = life[1] if len(life) > 1 else 0.0
    check("A2 noisy circle has one dominant loop", life[0] > 5 * max(second, 1e-9),
          f"top two H1 lifetimes = {life[0]:.3f}, {second:.3f}")

    # A3
    two = np.r_[rng.normal(0, 0.1, (30, 3)), rng.normal(5, 0.1, (30, 3))]
    h0 = T.persistence(two, False)[0]
    life0 = np.sort(h0[:, 1] - h0[:, 0])[::-1]
    check("A3 two clusters give one dominant finite H0 bar",
          life0[0] > 10 * life0[1], f"top two H0 lifetimes = {life0[0]:.3f}, {life0[1]:.3f}")

    # A4
    X = rng.normal(size=(40, 6))
    Q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    Y = X @ Q + 7.5
    dX, dY = T.persistence(X, False), T.persistence(Y, False)
    d = sum(T.wasserstein(dX[k], dY[k]) for k in (0, 1))
    check("A4 invariance to rotation + translation", d < 1e-8, f"W1 = {d:.2e}")

    # A5
    dX_n, dY_n = T.persistence(X, True), T.persistence(3.0 * X, True)
    d_norm = sum(T.wasserstein(dX_n[k], dY_n[k]) for k in (0, 1))
    dX_r, dY_r = T.persistence(X, False), T.persistence(3.0 * X, False)
    d_raw = sum(T.wasserstein(dX_r[k], dY_r[k]) for k in (0, 1))
    check("A5 shape-only is scale invariant, raw is not",
          d_norm < 1e-8 and d_raw > 1e-2, f"shape-only W1 = {d_norm:.2e}, raw W1 = {d_raw:.3f}")

    # A6
    p = rng.uniform(0, 0.2, 7)
    ref = multipletests(p, method="holm")[1]
    check("A6 Holm adjustment matches statsmodels", np.allclose(T.holm(p), ref),
          f"max |diff| = {np.max(np.abs(T.holm(p) - ref)):.1e}")


# ---------------------------------------------------------------- Part B
N, DIM, M, K, NPERM = 40, 6, 30, 8, 39   # P(p <= 0.05) = 2/40 = 0.05 exactly under H0


def ring(n, rng):
    th = rng.uniform(0, 2 * np.pi, n)
    X = rng.normal(0, 0.15, (n, DIM))
    X[:, 0] += 2 * np.cos(th)
    X[:, 1] += 2 * np.sin(th)
    return X


def blob(n, rng):
    return rng.normal(0, 1, (n, DIM))


def one_rep(job):
    scenario, align, normalize, seed = job
    rng = np.random.default_rng(seed)
    if scenario == "B1":
        Xs, Xt = blob(N, rng), blob(N, rng)
    elif scenario == "B2":
        Xs, Xt = blob(N, rng), blob(N, rng) + np.r_[4.0, np.zeros(DIM - 1)]
    elif scenario == "B3":
        Xs, Xt = blob(N, rng), 2.0 * blob(N, rng)
    else:
        Xs, Xt = blob(N, rng), ring(N, rng)
    emb = np.r_[Xs, Xt]
    idx_s, idx_t = np.arange(N), np.arange(N, 2 * N)
    _, p, _, _ = T.permutation_test(emb, idx_s, idx_t, M, K, (0, 1), normalize,
                                    NPERM, rng, align=align)
    return p


def part_b(reps):
    plan = [("B1", "H0 identical", True, True),
            ("B2", "H0 shift only", False, True),
            ("B2", "H0 shift only", True, True),
            ("B3", "H0 scale only", False, True),
            ("B3", "H0 scale only", True, True),
            ("B4", "H1 ring vs blob", True, True)]
    rows = []
    with Pool(2) as pool:
        for code, label, align, normalize in plan:
            t0 = time.time()
            jobs = [(code, align, normalize, 1000 + i) for i in range(reps)]
            ps = np.array(pool.map(one_rep, jobs))
            rate = np.mean(ps <= 0.05)
            se = np.sqrt(0.05 * 0.95 / reps)
            rows.append((code, label, align, rate))
            print(f"  {code} {label:<16} align={str(align):<5} rejection rate "
                  f"{rate:.3f}  ({reps} reps, {time.time() - t0:.0f}s)")
            if code == "B4":
                check(f"{code} power against a loop (shape-only)", rate >= 0.8,
                      f"power = {rate:.2f}")
            elif align:
                check(f"{code} type-I error with alignment (nominal 0.05)",
                      rate <= 0.05 + 2.5 * se, f"rate = {rate:.3f} (MC s.e. {se:.3f})")
            else:
                print(f"         (without alignment; shown for comparison)")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--reps", type=int, default=100)
    args = ap.parse_args()
    print("Part A: deterministic checks")
    part_a()
    rows = []
    if not args.quick:
        print(f"\nPart B: permutation-test calibration (n={N}/group, dim={DIM}, "
              f"subsample {M}, k={K}, {NPERM} permutations)")
        rows = part_b(args.reps)
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    with open("validation_report.txt", "w") as fh:
        fh.write("Validation of tda_crc_pilot.py (v2)\n" + "=" * 40 + "\n")
        for name, ok, detail in RESULTS:
            fh.write(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}\n")
        if rows:
            fh.write("\nCalibration (rejection rate at alpha = 0.05):\n")
            for code, label, align, rate in rows:
                fh.write(f"  {code} {label:<16} aligned={align!s:<5} {rate:.3f}\n")
        fh.write(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed\n")
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed "
          "(report written to validation_report.txt)")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
