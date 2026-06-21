"""Compare tip-displacement histories from the np=1 (reference) and np=2 (split)
runs of sms_bar_mp.tcl. Stdlib only.

PASS criterion: the np=2 (shared-node, cross-rank DeltaM reduction) history must
match the np=1 (whole-model, all-local injection) history to round-off. A failure
means the shared-node fictitious mass was NOT summed across the partition.
"""
import math
import sys

# Physical sanity bound (unit-load bar tip disp is O(1e-3)). A diverged np=1 and np=2
# can match each other's overflow garbage and falsely "pass" -- reject non-finite/huge.
PHYS_MAX = 1.0


def finite_ok(rows, tag):
    for _, d in rows:
        if not math.isfinite(d) or abs(d) > PHYS_MAX:
            print(f"FAIL: {tag} DIVERGED (non-finite or |disp|>{PHYS_MAX}: {d:.3e})")
            return False
    return True


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            rows.append((float(parts[0]), float(parts[-1])))
    return rows


def main():
    ref = load("tip_np1.out")
    par = load("tip_np2.out")
    n = min(len(ref), len(par))
    if n == 0:
        print("FAIL: empty output(s)")
        sys.exit(1)
    if len(ref) != len(par):
        print(f"WARN: step count differs ref={len(ref)} par={len(par)} (compare first {n})")

    if not finite_ok(ref[:n], "np=1") or not finite_ok(par[:n], "np=2"):
        sys.exit(1)

    max_abs = 0.0
    max_rel = 0.0
    peak = 0.0
    for i in range(n):
        d = abs(ref[i][1] - par[i][1])
        max_abs = max(max_abs, d)
        peak = max(peak, abs(ref[i][1]))
        denom = abs(ref[i][1])
        if denom > 1e-12:
            max_rel = max(max_rel, d / denom)

    print(f"steps compared : {n}")
    print(f"peak |tipDisp| : {peak:.6e}")
    print(f"final ref/par  : {ref[n-1][1]:.10e}  {par[n-1][1]:.10e}")
    print(f"max |abs diff| : {max_abs:.3e}")
    print(f"max  rel diff  : {max_rel:.3e}")

    # Round-off tolerance: identical arithmetic, only the assembly path differs.
    tol_abs = 1e-9 * max(peak, 1.0)
    if max_abs <= tol_abs:
        print(f"PASS: np=2 matches np=1 within {tol_abs:.1e} -> shared-node DeltaM reduced correctly")
        sys.exit(0)
    print(f"FAIL: np=2 diverges from np=1 (tol {tol_abs:.1e}) -> shared-node DeltaM NOT reduced")
    sys.exit(1)


if __name__ == "__main__":
    main()
