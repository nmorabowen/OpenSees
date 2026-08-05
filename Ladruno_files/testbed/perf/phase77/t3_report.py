"""ADR-77 T3 — the G3 verdict: cost per simulated second AT MATCHED ACCURACY.

Reading a winner off any single dt column is the mistake this report exists to
avoid: a scheme with more numerical damping is cheaper per step AND less
accurate at that step, so a same-dt comparison rewards it twice. The G3 number
is therefore obtained by interpolating each scheme's (error, cost) curve in
log-log to a COMMON target error and comparing the costs there.

G3 (from ADR-77 §6) authorizes an implicit Bathe beta1/beta2 implementation only
if TRBDF2 *and* HHT/GeneralizedAlpha each lose >= 1.3x on cost-per-simulated-
second at matched accuracy, or fail outright. This prints exactly that test.

Run with any python3.12:  python3.12 -S t3_report.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGETS = [float(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1
                              else ["5e-3", "1e-2"])]


def cost_at(points, target):
    """Log-log interpolate cost-per-simulated-second at a target relative error.

    Returns (cost, how) or (None, reason). Refuses to extrapolate beyond the
    measured bracket -- an extrapolated G3 verdict would be a guess wearing a
    number's clothes.
    """
    pts = sorted(((p["rel_err"], p["cost_per_sim_s"], p["dt"])
                  for p in points.values() if "rel_err" in p), key=lambda x: x[0])
    if len(pts) < 2:
        return None, "insufficient points"
    lo = [p for p in pts if p[0] <= target]
    hi = [p for p in pts if p[0] >= target]
    if not lo:
        return None, f"target below measured range (best err {pts[0][0]:.2e})"
    if not hi:
        return None, f"target above measured range (worst err {pts[-1][0]:.2e})"
    a, b = lo[-1], hi[0]
    if a[0] == b[0]:
        return a[1], "exact"
    # log-log linear interpolation of cost vs error
    w = (math.log(target) - math.log(a[0])) / (math.log(b[0]) - math.log(a[0]))
    cost = math.exp(math.log(a[1]) + w * (math.log(b[1]) - math.log(a[1])))
    return cost, f"interp dt {a[2]:g}..{b[2]:g}"


def order(points):
    """Observed convergence order from the two finest dt that both converged."""
    pts = sorted(((p["dt"], p["rel_err"]) for p in points.values()
                  if "rel_err" in p), key=lambda x: x[0])
    if len(pts) < 2:
        return None
    (d1, e1), (d2, e2) = pts[0], pts[1]
    if e1 <= 0 or e2 <= 0 or d1 == d2:
        return None
    return math.log(e2 / e1) / math.log(d2 / d1)


def main():
    p = HERE / "t3_scheme_shootout.json"
    if not p.exists():
        print("no t3_scheme_shootout.json")
        return
    d = json.loads(p.read_text())
    print(f"ADR-77 T3 — n={d['n']} load x{d['load']:g} T_end={d['T_end']:g} "
          f"algorithm={d['algorithm']}")
    print(f"reference dt={d['dt_ref']:g}, self-convergence (L2) = "
          f"{d['ref_self_conv']:.2e}")
    floor = d["ref_self_conv"]

    print(f"\n{'scheme':<34} {'rho_inf':>8} {'obs.order':>10} "
          + "".join(f"{'cost@'+f'{t:g}':>14}" for t in TARGETS))
    rows = []
    for key, s in d["schemes"].items():
        o = order(s["points"])
        cells, costs = [], {}
        for t in TARGETS:
            c, how = cost_at(s["points"], t)
            costs[t] = c
            cells.append(f"{c:14.1f}" if c else f"{'--':>14}")
        rho = f"{s['rho_inf']:.3f}" if s.get("rho_inf") is not None else "-"
        print(f"{s['label']:<34} {rho:>8} "
              f"{(f'{o:.2f}' if o else '-'):>10} " + "".join(cells))
        rows.append((key, s["label"], costs))

    print("\n  cost = wall-clock seconds per simulated second, lower is better")
    print(f"  NOTE: errors below the reference floor {floor:.1e} are not resolved;")
    print("        any target at or under that floor is not a real comparison.")

    # ---- the G3 test, stated explicitly -------------------------------------
    print("\n=== G3 test — is an implicit Bathe beta1/beta2 authorized? ===")
    print("G3: authorize ONLY if TRBDF2 *and* HHT/GeneralizedAlpha each lose")
    print("    >=1.3x on cost-per-simulated-second at matched accuracy, or fail.")
    for t in TARGETS:
        if t <= floor:
            print(f"\n  target {t:g}: SKIPPED — at/below the reference floor")
            continue
        avail = [(lbl, c[t]) for _, lbl, c in rows if c.get(t)]
        if not avail:
            print(f"\n  target {t:g}: no scheme bracketed this error")
            continue
        best_lbl, best = min(avail, key=lambda x: x[1])
        print(f"\n  target err {t:g} — best = {best_lbl} at {best:.1f}")
        for lbl, c in sorted(avail, key=lambda x: x[1]):
            print(f"    {lbl:<34} {c:9.1f}  {c/best:5.2f}x")
        trb = dict(avail).get("TRBDF2 (composite Trap+BDF2)")
        damped = [c for lbl, c in avail
                  if lbl.startswith("HHT") or lbl.startswith("GeneralizedAlpha")]
        if trb and damped:
            trb_r, damp_r = trb / best, min(damped) / best
            verdict = ("AUTHORIZED" if trb_r >= 1.3 and damp_r >= 1.3
                       else "NOT authorized")
            print(f"    -> TRBDF2 {trb_r:.2f}x, best damped {damp_r:.2f}x "
                  f"=> G3 {verdict}")


if __name__ == "__main__":
    main()
