"""K2 sensitivity table (WP-144 G1; sheet 144a section 14; plan 5.2 K2).

Runs the AB06 6.1 stress-point localization protocol (S.43) for BOTH oracles over the swept
unknowns of sheet section 14 and prints the full table as markdown:

    pi_i0 in {-60.4, -80, -100}  x  chi in {-3.0, -3.5, -4.0}  x  v_c0 in {1.80, 1.81, 1.82}
    x both crossing criteria ("first": first step with min det <= 0;
                              "interp": linear-interpolated zero crossing of the min-det curve)
    x (rho, rho_bar) in {(0.7, 0.8), (1, 1)}.

Usage (from norsand_oracle/):   python tests/k2_sensitivity_table.py > tests/out/k2_sensitivity.md

This module also holds the shared runner `run_case` and the sweep grids, imported by
test_g1_k2_benchmark.py. The expected bands (ordering, gap, nominal n bands) are NOT here:
they live in the test file, written before any oracle output was looked at.
"""
from __future__ import annotations

import itertools
import os
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conftest import K2_BASE, ORACLES, make_params  # noqa: E402

# Sheet 14 sweep grids (the swept unknowns) and the nominal combination.
PI_I0_GRID = (-60.4, -80.0, -100.0)
CHI_GRID = (-3.0, -3.5, -4.0)
VC0_GRID = (1.80, 1.81, 1.82)
NOMINAL = dict(pi_i0=-60.4, chi=-3.5, v_c0=1.81)
CASES_RHO = {"0.7": (0.7, 0.8), "1.0": (1.0, 1.0)}   # label -> (rho, rho_bar)
N_MAX = 60          # stop the search here; a crossing beyond it is reported as None
V0 = 1.59           # AB06 6.1 initial specific volume (sheet 14)
SIGMA0 = -100.0 * __import__("numpy").eye(3)   # sheet 14 "Inferable": isotropic -100 kPa


def run_case(oracle: str, pi_i0: float, chi: float, v_c0: float, rho: float, rho_bar: float):
    """One K2 localization run. Returns (n_first, n_interp); None where no crossing by N_MAX."""
    mod = ORACLES[oracle]
    kw = dict(K2_BASE)
    kw.update(chi=chi, v_c0=v_c0, rho=rho, rho_bar=rho_bar)
    params = make_params(oracle, **kw)
    st0 = mod.initial_state(params, SIGMA0, V0, pi_i0)
    if oracle == "O1":
        r = mod.k2_path(params, st0, N_MAX, extra_after=0)
    else:
        r = mod.k2_path(params, st0, N_MAX)
    return r["n_first"], r["n_interp"]


def _job(args):
    t = time.time()
    return args, run_case(*args), time.time() - t


def sweep_jobs():
    jobs = []
    for oracle in ORACLES:
        for pi0, chi, vc in itertools.product(PI_I0_GRID, CHI_GRID, VC0_GRID):
            for lab in CASES_RHO:
                rho, rb = CASES_RHO[lab]
                jobs.append((oracle, pi0, chi, vc, rho, rb))
    return jobs


def _fmt(x):
    return "none" if x is None else (f"{x:.2f}" if isinstance(x, float) else str(x))


def main():
    import test_g1_k2_benchmark as g1   # deferred: it imports this module at its top
    jobs = sweep_jobs()
    t0 = time.time()
    with Pool(min(len(jobs), os.cpu_count() or 1)) as pool:
        res = {a: (v, dt) for a, v, dt in pool.imap_unordered(_job, jobs)}
    wall = time.time() - t0
    print("# K2 sensitivity table (AB06 6.1, sheet 144a section 14)\n")
    print(f"Both oracles, N_MAX = {N_MAX}, v0 = {V0}, sigma0 = -100 kPa isotropic; "
          f"wall time {wall:.0f} s on {os.cpu_count()} cores "
          f"(sum of per-run times {sum(dt for _, dt in res.values()):.0f} s).\n")
    print(f"Gate per row and criterion: n_0.7 < n_1.0 and gap in [{g1.GAP_LO}, {g1.GAP_HI}]. "
          f"Paper: n = 22 (rho 0.7 / rho_bar 0.8), 26 (rho = rho_bar = 1). "
          f"Nominal = (pi_i0 {NOMINAL['pi_i0']}, chi {NOMINAL['chi']}, v_c0 {NOMINAL['v_c0']}), "
          f"first-step criterion.\n")
    print("| oracle | pi_i0 | chi | v_c0 | n_0.7 first | n_1.0 first | gap first | ok | "
          "n_0.7 interp | n_1.0 interp | gap interp | ok |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    bad = 0
    for oracle in ORACLES:
        for pi0, chi, vc in itertools.product(PI_I0_GRID, CHI_GRID, VC0_GRID):
            a = res[(oracle, pi0, chi, vc, *CASES_RHO["0.7"])][0]
            b = res[(oracle, pi0, chi, vc, *CASES_RHO["1.0"])][0]
            cells = []
            for k in (0, 1):
                n07, n10 = a[k], b[k]
                ok = g1.ordering_and_gap_ok(n07, n10)
                bad += (not ok)
                gap = None if (n07 is None or n10 is None) else n10 - n07
                cells += [_fmt(n07), _fmt(n10), _fmt(gap), "yes" if ok else "**NO**"]
            print(f"| {oracle} | {pi0} | {chi} | {vc} | " + " | ".join(cells) + " |")
    print(f"\nCells failing ordering + gap: {bad} of {len(res)} "
          f"({len(res) // 2} (oracle, combination) pairs x 2 criteria; {len(res)} runs).")


if __name__ == "__main__":
    main()
