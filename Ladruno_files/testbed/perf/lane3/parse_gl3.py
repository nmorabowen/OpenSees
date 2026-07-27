#!/usr/bin/env python3
"""ADR-75b G-L3 — reduce the coarse profiler rollup to the per-LOOP fractions the
gate is defined on.

The gate is on a SINGLE loop's share of step, NOT the aggregate element-kernel
fraction (ADR-75b §12.4: Lane D reads kernel 85.30% but gate FAIL at loop A
38.95%). Loops, per §2's taxonomy:

    loop A = Domain::update   -> rollup node `update`
    loop B = element tangent  -> rollup node `formTangent`
    loop C = element residual -> rollup node `formUnbalance`

Usage: parse_gl3.py <file.h5> [more.h5 ...]
"""
import sys
import h5py

PHASES = [("A", "update"), ("B", "formTangent"), ("C", "formUnbalance")]
OTHER = ["linearSolve", "convTest", "commit", "newStep"]


def walk(node, out, prefix=""):
    for k in node:
        o = node[k]
        if isinstance(o, h5py.Group):
            out[k] = int(o.attrs.get("wall_ns", 0))
            walk(o, out, prefix + k + "/")


def main(paths):
    print(f"{'file':<26} {'step_ms':>9} {'A upd%':>8} {'B tan%':>8} "
          f"{'C unb%':>8} {'solve%':>8} {'maxloop':>8} {'gate':>6}")
    for p in paths:
        with h5py.File(p, "r") as f:
            try:
                root = f["runs/run0/rollup/root"]
            except KeyError:
                print(f"{p:<26} NO ROLLUP")
                continue
            vals = {}
            walk(root, vals)
            step = vals.get("step", 0)
            if step <= 0:
                print(f"{p:<26} step=0")
                continue
            fr = {tag: 100.0 * vals.get(name, 0) / step for tag, name in PHASES}
            solve = 100.0 * vals.get("linearSolve", 0) / step
            mx = max(fr.values())
            gate = "PASS" if mx >= 40.0 else "FAIL"
            name = p.split("/")[-1]
            print(f"{name:<26} {step/1e6:>9.1f} {fr['A']:>8.2f} {fr['B']:>8.2f} "
                  f"{fr['C']:>8.2f} {solve:>8.2f} {mx:>8.2f} {gate:>6}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
