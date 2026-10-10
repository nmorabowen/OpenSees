"""WP-130 / F18(c) rollup of the bearing-deck arms (bvp/<arm>/).

Per arm: committed steps, s/B reached, q at the end, wall; ladder use; Newton
iterations per COMMITTED step (the analyze call that returned 0); the observed
convergence order on those calls; depth at matched wall clock; the mesh-summed
`substepStats` census (CPPM columns). Prints Markdown tables.

Convergence order: on each committed call with >= 3 decreasing unbalance norms
r_0 > r_1 > ... (testNorm), p_k = log(r_{k+1}/r_k) / log(r_k/r_{k-1}) over the
last three; quadratic is p ~ 2, linear p ~ 1. We report the median of the LAST
p of each call, and the share of calls whose last p >= 1.8.
"""
import csv
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LEG = "a2_h1.0_e0.6944"
WALLS = (60, 300, 600, 900, 1200)


def curve(arm):
    p = os.path.join(HERE, "bvp", arm, LEG + "_curve.csv")
    if not os.path.exists(p):
        return []
    with open(p, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines))


def calls(arm):
    p = os.path.join(HERE, "bvp", arm, "analyze_log.jsonl")
    if not os.path.exists(p):
        return []
    return [json.loads(ln) for ln in open(p, encoding="utf-8") if ln.strip()]


def order(norms):
    r = [x for x in norms if x > 0]
    if len(r) < 3:
        return None
    a, b, c = r[-3], r[-2], r[-1]
    if not (a > b > c) or b / a >= 1 or c / b >= 1:
        return None
    return math.log(c / b) / math.log(b / a)


def main(arms):
    print("| arm | steps | s/B reached | q end (kPa) | wall s | analyze calls | failed calls | wall in failed calls (s) | it/committed step (mean, max) | rung use N/NLS/K | last-order median | calls with order>=1.8 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    census = {}
    depth = {}
    for arm in arms:
        rows = curve(arm)
        cs = calls(arm)
        push = cs[10:] if len(cs) > 10 else []   # the first 10 calls are the stage-0 gravity steps
        ok = [c for c in push if c["rc"] == 0]
        bad = [c for c in push if c["rc"] != 0]
        its = [c["iters"] for c in ok]
        rung = {"Newton": 0, "NewtonLineSearch": 0, "KrylovNewton": 0}
        for c in ok:
            rung[c["algo"]] = rung.get(c["algo"], 0) + 1
        ords = [o for o in (order(c["norms"]) for c in ok) if o is not None]
        sB = float(rows[-1]["s_over_B"]) if rows else 0.0
        q = float(rows[-1]["q_foot_kPa"]) if rows else float("nan")
        wall = float(rows[-1]["wall_s"]) if rows else 0.0
        print(f"| {arm} | {len(rows)} | {sB:.5f} | {q:.2f} | {wall:.0f} | {len(push)} | {len(bad)} | "
              f"{sum(c['wall'] for c in bad):.1f} | "
              f"{(statistics.mean(its) if its else float('nan')):.2f}, {max(its) if its else 0} | "
              f"{rung['Newton']}/{rung['NewtonLineSearch']}/{rung['KrylovNewton']} | "
              f"{(statistics.median(ords) if ords else float('nan')):.2f} | "
              f"{sum(1 for o in ords if o >= 1.8)}/{len(ords)} |")
        depth[arm] = [max([float(r["s_over_B"]) for r in rows if float(r["wall_s"]) <= w] or [0.0]) for w in WALLS]
        cp = os.path.join(HERE, "bvp", arm, "census.json")
        if os.path.exists(cp):
            census[arm] = json.load(open(cp)).get("census") or {}
    print()
    print("| arm | " + " | ".join(f"s/B at {w} s" for w in WALLS) + " |")
    print("|---|" + "---|" * len(WALLS))
    for arm in arms:
        print(f"| {arm} | " + " | ".join(f"{x:.5f}" for x in depth[arm]) + " |")
    keys = ["meCalls", "substeps", "cppmCalls", "cppmNewtonFail", "cppmHalvings", "cppmExplicitFail",
            "cppmExplicitLowP", "cppmRefusals", "cppmGuessTries", "cppmGuessOk"]
    print()
    print("| arm | " + " | ".join(keys) + " |")
    print("|---|" + "---|" * len(keys))
    for arm in arms:
        c = census.get(arm, {})
        print(f"| {arm} | " + " | ".join(f"{int(c.get(k, 0))}" for k in keys) + " |")


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(os.listdir(os.path.join(HERE, "bvp"))))
