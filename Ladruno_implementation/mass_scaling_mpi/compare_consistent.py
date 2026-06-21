"""Compare CONSISTENT (Olovsson) SMS tip-disp histories: np=1 (reference) vs np=2
(split at the shared node). Stdlib only.

Unlike the lumped path (a direct diagonal solve -> bit-identical), the consistent
path is an iterative matrix-free PCG, so np=2 matches np=1 to the SOLVER tolerance
accumulated over the run, not to round-off. PASS = agreement well inside what a
broken (rank-local, boundary-coupling-dropped) PCG would produce.
"""
import math
import sys

# Physical sanity bound: the unit-load bar's tip disp is O(1e-3). Anything far above
# this (or non-finite) means the run DIVERGED -- a diverged np=1 and np=2 can match
# each other's garbage (identical overflow) and falsely "pass", so guard explicitly.
PHYS_MAX = 1.0


def finite_ok(rows, tag):
    for _, d in rows:
        if not math.isfinite(d) or abs(d) > PHYS_MAX:
            print(f"FAIL: {tag} DIVERGED (non-finite or |disp|>{PHYS_MAX}: {d:.3e}) -- "
                  "the consistent solve did not stabilize")
            return False
    return True


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                p = line.split()
                rows.append((float(p[0]), float(p[-1])))
    return rows


def main():
    ref = load("tipc_np1.out")
    par = load("tipc_np2.out")
    n = min(len(ref), len(par))
    if n == 0:
        print("FAIL: empty output(s)")
        sys.exit(1)
    if len(ref) != len(par):
        print(f"WARN: step count differs ref={len(ref)} par={len(par)}")

    if not finite_ok(ref[:n], "np=1") or not finite_ok(par[:n], "np=2"):
        sys.exit(1)

    max_abs = peak = 0.0
    for i in range(n):
        max_abs = max(max_abs, abs(ref[i][1] - par[i][1]))
        peak = max(peak, abs(ref[i][1]))

    rel = max_abs / peak if peak > 0 else max_abs
    print(f"steps compared : {n}")
    print(f"peak |tipDisp| : {peak:.6e}")
    print(f"final ref/par  : {ref[n-1][1]:.10e}  {par[n-1][1]:.10e}")
    print(f"max |abs diff| : {max_abs:.3e}")
    print(f"max  rel diff  : {rel:.3e}")

    tol_rel = 1.0e-6   # iterative PCG agreement over 150 steps; a dropped-coupling bug is O(1)
    if rel <= tol_rel:
        print(f"PASS: np=2 matches np=1 within {tol_rel:.0e} rel -> distributed consistent PCG correct")
        sys.exit(0)
    print(f"FAIL: np=2 diverges from np=1 (rel {rel:.3e} > {tol_rel:.0e}) -> boundary coupling dropped")
    sys.exit(1)


if __name__ == "__main__":
    main()
