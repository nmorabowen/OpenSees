"""WP-106 / ADR-93 II.1 -- roll up the two BVP arms of the strip-footing leg.

    python3.12 bvp_summary.py <pRe0_dir> <pRe1_dir>

Reads the `a2_<tag>_curve.csv` and `a2_<tag>.json` each `sanisand_tau0_band.py`
run writes, and reports what the F13 question asks for:

  * per-step ModifiedEuler substep totals over the WHOLE mesh and over the
    free-surface RING beside the footing edge, plus the per-point max in each
  * wall seconds per step
  * the committed load-settlement delta between the arms over the settlement
    range BOTH reached (interpolated onto the coarser arm's own s/B grid, so the
    delta is a curve difference and not an artefact of different step sizes)
"""
from __future__ import annotations

import csv
import glob
import json
import os
import statistics
import sys


def load(d):
    cur = glob.glob(os.path.join(d, "a2_*_curve.csv"))
    jsn = [p for p in glob.glob(os.path.join(d, "a2_*.json"))]
    assert cur, f"no curve csv under {d}"
    body = [l for l in open(cur[0], encoding="utf-8") if not l.startswith("#")]
    rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(body)]
    meta = json.load(open(jsn[0], encoding="utf-8")) if jsn else {}
    return meta, rows


def interp(rows, key, x):
    """Linear in s/B."""
    xs = [r["s_over_B"] for r in rows]
    ys = [r[key] for r in rows]
    if x <= xs[0]:
        return ys[0]
    for i in range(1, len(xs)):
        if xs[i] >= x:
            t = (x - xs[i - 1]) / max(xs[i] - xs[i - 1], 1e-30)
            return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return ys[-1]


def col(rows, k):
    return [r[k] for r in rows if k in r]


def main():
    dirs = sys.argv[1:]
    assert len(dirs) == 2, __doc__
    arms = [load(d) for d in dirs]

    print("| arm | pRe | mode | steps | s/B reached | q_max [kPa] | wall [s] | "
          "wall/step [s] | substeps total | /step (median) | max @ any point | "
          "RING total | RING /step (median) | RING max |")
    print("|" + "|".join("---" for _ in range(14)) + "|")
    for d, (meta, rows) in zip(dirs, arms):
        st = col(rows, "substeps_tot")
        rg = col(rows, "substeps_ring_tot")
        ws = col(rows, "wall_step_s")
        print("| " + " | ".join([
            os.path.basename(d), f"{meta.get('pre', 0):g}",
            str(meta.get("mode", "?")), str(len(rows)),
            f"{rows[-1]['s_over_B']:.5f}",
            f"{max(r['q_foot_kPa'] for r in rows):.2f}",
            f"{rows[-1]['wall_s']:.0f}",
            f"{statistics.median(ws):.2f}" if ws else "-",
            f"{sum(st):.0f}" if st else "-",
            f"{statistics.median(st):.0f}" if st else "-",
            f"{max(col(rows, 'substeps_max')):.0f}" if st else "-",
            f"{sum(rg):.0f}" if rg else "-",
            f"{statistics.median(rg):.0f}" if rg else "-",
            f"{max(col(rows, 'substeps_ring_max')):.0f}" if rg else "-",
        ]) + " |")

    (m0, r0), (m1, r1) = arms
    print()
    print(f"gravity state: p_min after K0 gravity = "
          f"{m0.get('p_min_grav')} kPa (pRe 0) / {m1.get('p_min_grav')} kPa (pRe 1); "
          f"surcharge {m0.get('surcharge_kpa')} kPa; "
          f"nodes {m0.get('nodes')}, hexes {m0.get('hexes')}")

    # --- the committed load-settlement delta over the COMMON range -----------
    smax = min(r0[-1]["s_over_B"], r1[-1]["s_over_B"])
    smin = max(r0[0]["s_over_B"], r1[0]["s_over_B"])
    n = 60
    ds = []
    for i in range(n + 1):
        x = smin + (smax - smin) * i / n
        a = interp(r0, "q_foot_kPa", x)
        b = interp(r1, "q_foot_kPa", x)
        ds.append((x, a, b, (b - a) / a if abs(a) > 1e-12 else 0.0))
    rel = [abs(t[3]) for t in ds]
    worst = max(ds, key=lambda t: abs(t[3]))
    print(f"\nload-settlement delta over the COMMON range s/B in "
          f"[{smin:.5f}, {smax:.5f}] ({n + 1} samples):")
    print(f"  median |dq/q| = {statistics.median(rel) * 100:.3f} %   "
          f"mean = {sum(rel) / len(rel) * 100:.3f} %   "
          f"max = {max(rel) * 100:.3f} % at s/B = {worst[0]:.5f} "
          f"(q {worst[1]:.3f} -> {worst[2]:.3f} kPa)")
    print(f"  q at the common end s/B = {smax:.5f}: "
          f"{ds[-1][1]:.3f} -> {ds[-1][2]:.3f} kPa "
          f"({(ds[-1][2] - ds[-1][1]) / ds[-1][1] * 100:+.3f} %)")


if __name__ == "__main__":
    main()
