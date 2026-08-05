"""ADR-77 T0 — turn the profiler .h5 rollups into the five-loop step-anatomy table.

Reads the per-arm profiles written by t0_step_anatomy.py and prints, per solver,
each loop's share of `step`. The two rows that did not exist before this ADR are
`dof.tangent` and `dof.residual` (the DOF_Group loops).

The T2 headroom line is the sum of the *constant-matrix* work: element inertia
(getMass, re-formed every Newton iteration in both formTangent and formUnbalance)
plus the DOF_Group loops. That number is the ceiling of ADR-77's G2 cache, and G2
requires >= 10% of step to authorize a design.

Run (needs h5py — use the opensees_env interpreter, NOT python3.12 -S):
  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe t0_report.py [n]
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
N = sys.argv[1] if len(sys.argv) > 1 else "15"

# loop -> h5 path under rollup/root (leaf name is matched, path is for clarity)
LOOPS = [
    ("newStep",                 "step/newStep"),
    ("  elem.update (newStep)", "step/newStep/elem.update"),
    ("formTangent",             "step/solveCurrentStep/formTangent"),
    ("  dof.tangent  *NEW*",    "step/solveCurrentStep/formTangent/dof.tangent"),
    ("  elem.tangent",          "step/solveCurrentStep/formTangent/elem.tangent"),
    ("    brick.inertia",       "step/solveCurrentStep/formTangent/elem.tangent/brick.inertia"),
    ("formUnbalance",           "step/solveCurrentStep/formUnbalance"),
    ("  dof.residual *NEW*",    "step/solveCurrentStep/formUnbalance/dof.residual"),
    ("  elem.residual",         "step/solveCurrentStep/formUnbalance/elem.residual"),
    ("    brick.inertia",       "step/solveCurrentStep/formUnbalance/elem.residual/brick.inertia"),
    ("linearSolve",             "step/solveCurrentStep/linearSolve"),
    ("update",                  "step/solveCurrentStep/update"),
    ("  elem.update",           "step/solveCurrentStep/update/elem.update"),
    ("convTest",                "step/solveCurrentStep/convTest"),
    ("commit",                  "step/commit"),
    ("domainChanged",           "step/domainChanged"),
]


def wall(root, path):
    node = root
    for part in path.split("/"):
        if part not in node:
            return None
        node = node[part]
    return int(node.attrs.get("wall_ns", 0))


def report(arm, n):
    p = HERE / f"t0_{arm}_n{n}.h5"
    if not p.exists():
        print(f"  (missing {p.name})")
        return None
    with h5py.File(p, "r") as f:
        root = f[f"runs/{arm}/rollup/root"]
        step = int(root["step"].attrs.get("wall_ns", 0))
        print(f"\n=== {arm.upper()}  (step total = {step/1e6:.1f} ms) ===")
        print(f"{'loop':<26} {'ms':>10} {'% step':>8}")
        vals = {}
        for label, path in LOOPS:
            w = wall(root, path)
            if w is None:
                continue
            vals[label.strip().split()[0]] = w
            print(f"{label:<26} {w/1e6:>10.2f} {100.0*w/step:>7.2f}%")
        # G2 headroom: constant-matrix work re-done every iteration
        inertia = (wall(root, "step/solveCurrentStep/formTangent/elem.tangent/brick.inertia") or 0) \
                + (wall(root, "step/solveCurrentStep/formUnbalance/elem.residual/brick.inertia") or 0)
        dofl = (wall(root, "step/solveCurrentStep/formTangent/dof.tangent") or 0) \
             + (wall(root, "step/solveCurrentStep/formUnbalance/dof.residual") or 0)
        head = inertia + dofl
        print(f"{'-'*46}")
        print(f"{'G2 headroom (inertia+dof)':<26} {head/1e6:>10.2f} {100.0*head/step:>7.2f}%"
              f"   [gate: >=10% to authorize]")
        return {"step": step, "headroom": head,
                "elem_tangent": wall(root, "step/solveCurrentStep/formTangent/elem.tangent") or 0,
                "linearSolve": wall(root, "step/solveCurrentStep/linearSolve") or 0}


def main():
    print(f"ADR-77 T0 five-loop step anatomy | n={N}")
    res = {}
    for arm in ("umfpack", "pardiso"):
        r = report(arm, N)
        if r:
            res[arm] = r
    if len(res) == 2:
        print("\n=== L3-0 TRANSPLANT CHECK (G0) ===")
        for arm, r in res.items():
            et = 100.0 * r["elem_tangent"] / r["step"]
            ls = 100.0 * r["linearSolve"] / r["step"]
            print(f"  {arm:<8} elem.tangent={et:5.1f}%  linearSolve={ls:5.1f}%")
        print("  (ADR-75b L3-0 measured element 35.8% -> 74.9% on the STATIC lane;"
              " G0 forbids assuming the same flip here)")


if __name__ == "__main__":
    main()
