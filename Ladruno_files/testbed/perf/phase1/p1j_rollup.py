"""ADR-75 P1j — roll the per-size phase-split JSONs into the trend table.

Usage (needs h5py-free input: run p1h_parse_split.py first to make the JSONs):
    python3.12 p1j_rollup.py

Reads p1j_<tag>_n<N>_r<R>.json (written by p1h_parse_split.py) and prints the trend.

THE HEADLINE METRIC IS fac/(fac+tri), NOT fac/solver.

`soe.symbolic` runs once per sparsity pattern while factor/trisolve run once per
solve, so symbolic's share of the solver depends on how many solves amortize it --
i.e. on the STEP COUNT, not on the model size. Including it in the denominator would
let a step-count choice leak into a size trend. The sweep holds steps constant so
either denominator would be defensible, but fac/(fac+tri) is the one that cannot be
gamed, and it is also the ratio the reuse levers actually bound: they skip
factorization and always pay the triangular solve.

Also fits a log-log exponent for each phase against nDOF, because the QUESTION is
whether factorization is superlinear relative to substitution -- a slope comparison
answers it directly, where two endpoint percentages could be a fluke of the box.
"""
import json
import glob
import math
import os
import re
import statistics
import collections

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    rows = collections.defaultdict(list)   # (tag, n) -> [split dicts]
    for f in sorted(glob.glob(os.path.join(HERE, "p1j_*_n*_r*.json"))):
        m = re.match(r"p1j_(.+)_n(\d+)_r(\d+)\.json$", os.path.basename(f))
        if not m:
            continue
        tag, n = m.group(1), int(m.group(2))
        rows[(tag, n)].append(json.load(open(f)))
    return rows


def med(rs, *keys):
    """Explicit key tuple: the JSON has keys that CONTAIN dots (soe.factor), so a
    dotted-path string accessor silently returns 0 for every one of them."""
    vals = []
    for r in rs:
        d = r
        for k in keys:
            d = d.get(k) if isinstance(d, dict) else None
            if d is None:
                break
        vals.append(d if isinstance(d, (int, float)) else 0.0)
    return statistics.median(vals) if vals else 0.0


def loglog_slope(xs, ys):
    """Least-squares slope of log(y) vs log(x) -- the scaling exponent."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(pts) < 2:
        return float("nan")
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    num = sum((p[0] - mx) * (p[1] - my) for p in pts)
    den = sum((p[0] - mx) ** 2 for p in pts)
    return num / den if den else float("nan")


def check_not_accumulated(tol=0.25):
    """Cross-check every profile's `step` total against the wall clock the sweep
    measured around the SAME loop, and REFUSE to print a table if they disagree.

    Why this guard exists: the Profiler is a process-global singleton and
    `ops.wipe()` does not touch it, so a multi-run script that forgets
    `profiler reset` gets a report containing the sum of every run so far. The first
    P1j sweep did exactly that -- soe.factor call counts ran 44, 88, 132 ... 936 --
    and it did NOT look like an error. It looked like a clean, monotone, entirely
    plausible size trend, with every ratio wrong. Accumulation is invisible to any
    check that only reads the profiler; the only way to catch it is to compare
    against a clock the profiler cannot influence.
    """
    src = os.path.join(HERE, "p1j_size_trend.json")
    if not os.path.exists(src):
        print("!! p1j_size_trend.json missing -- cannot verify accumulation; refusing")
        return False
    walls = json.load(open(src)).get("wall_by_run") or {}
    if not walls:
        print("!! p1j_size_trend.json has no wall_by_run (pre-guard sweep) -- refusing.\n"
              "   Re-run p1j_size_trend.py; a table from an unverified sweep is not "
              "publishable.")
        return False
    bad = []
    for stem, wall in sorted(walls.items()):
        f = os.path.join(HERE, stem + ".json")
        if not os.path.exists(f):
            bad.append((stem, wall, None, "no parsed json"))
            continue
        step_s = json.load(open(f)).get("step_ms", 0.0) / 1000.0
        rel = abs(step_s - wall) / wall if wall else 1.0
        if rel > tol:
            bad.append((stem, wall, step_s, f"{rel*100:.0f}% off"))
    if bad:
        print("!! PROFILER/WALL MISMATCH -- almost certainly missing `profiler reset`:")
        for stem, wall, step_s, why in bad:
            got = f"{step_s:8.2f}" if step_s is not None else "     n/a"
            print(f"     {stem:<28} wall={wall:8.2f}s  profiler step={got}s  ({why})")
        print("   Refusing to print a trend table from these files.")
        return False
    print(f"accumulation guard: OK -- all {len(walls)} profiles agree with the wall "
          f"clock to within {int(tol*100)}%")
    return True


def main():
    rows = load()
    if not rows:
        print("no p1j_*_n*_r*.json found -- run p1h_parse_split.py over the h5 first")
        return
    if not check_not_accumulated():
        return
    tags = sorted({t for (t, _) in rows})
    series = {}

    for tag in tags:
        ns = sorted(n for (t, n) in rows if t == tag)
        print(f"\n=== {tag} ===")
        h = (f"{'n':>4} {'nDOF':>8} {'step s':>8} {'solver s':>9} {'%step':>6} | "
             f"{'symb ms':>8} {'factor ms':>10} {'trisol ms':>10} | "
             f"{'FAC/(FAC+TRI)':>14} {'fac/tri':>8}")
        print(h)
        print("-" * len(h))
        cols = {"ndof": [], "fac": [], "tri": [], "sym": [], "step": [], "solver": [],
                "ratio": []}
        for n in ns:
            rs = rows[(tag, n)]
            ndof = med(rs, "meta", "nDOF")
            step = med(rs, "step_ms")
            solver = med(rs, "split", "solver_total_ms")
            fac = med(rs, "solver", "soe.factor", "ms")
            tri = med(rs, "solver", "soe.trisolve", "ms")
            sym = med(rs, "solver", "soe.symbolic", "ms")
            ratio = 100.0 * fac / (fac + tri) if (fac + tri) else 0.0
            print(f"{n:>4} {ndof:>8.0f} {step/1000:>8.2f} {solver/1000:>9.2f} "
                  f"{100.0*solver/step:>6.2f} | {sym:>8.1f} {fac:>10.1f} {tri:>10.1f} | "
                  f"{ratio:>13.2f}% {fac/tri if tri else 0:>8.2f}")
            for k, v in (("ndof", ndof), ("fac", fac), ("tri", tri), ("sym", sym),
                         ("step", step), ("solver", solver), ("ratio", ratio)):
                cols[k].append(v)
        series[tag] = cols

        print(f"  scaling exponents vs nDOF (log-log slope, n={len(ns)} points):")
        for k, label in (("fac", "soe.factor"), ("tri", "soe.trisolve"),
                         ("sym", "soe.symbolic"), ("step", "step"),
                         ("solver", "solver total")):
            print(f"    {label:<14} ~ N^{loglog_slope(cols['ndof'], cols[k]):.3f}")
        d = cols["ratio"][-1] - cols["ratio"][0]
        print(f"  fac/(fac+tri): {cols['ratio'][0]:.2f}% -> {cols['ratio'][-1]:.2f}% "
              f"({d:+.2f} pp across {cols['ndof'][0]:.0f} -> {cols['ndof'][-1]:.0f} DOF)")

    out = os.path.join(HERE, "p1j_trend_rollup.json")
    with open(out, "w") as f:
        json.dump(series, f, indent=2)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
