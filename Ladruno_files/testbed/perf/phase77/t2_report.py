"""ADR-77 T2 — the G2 verdict on a nodal-mass deck.

Prints the headroom a constant-M/C tangent cache could delete, per arm, against
G2's >=10%-of-step gate. The point of the nodal arm is that T0's deck routed all
its inertia through the ELEMENT loop, so the DOF_Group loop it would have to pay
for measured 0.27% -- which is not evidence about decks that lump mass on nodes.

Headroom components, kept separate because they are deleted by different means:
  dof.tangent / dof.residual   nodal M/C assembly (transient-only loops)
  brick.inertia (x2)           element mass re-formed inside elem.tangent and
                               elem.residual every iteration

Run with the h5py interpreter:
  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe t2_report.py
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
GATE = 10.0
ARMS = [("material", "A: material mass (= T0 deck)"),
        ("nodal", "B: nodal mass"),
        ("nodal_rayleigh", "C: nodal + rayleigh betaK0 (C const)"),
        ("nodal_rayleigh_kt", "D: nodal + rayleigh betaKcurr")]


def wall(root, path):
    node = root
    for part in path.split("/"):
        if part not in node:
            return 0
        node = node[part]
    return int(node.attrs.get("wall_ns", 0))


def read(key, n):
    p = HERE / f"t2_{key}_n{n}.h5"
    if not p.exists():
        return None
    with h5py.File(p, "r") as f:
        r = f[f"runs/{key}/rollup/root"]
        step = int(r["step"].attrs.get("wall_ns", 0))
        if not step:
            return None
        scs = "step/solveCurrentStep"
        d = {
            "step_ms": step / 1e6,
            "dof_tan": wall(r, f"{scs}/formTangent/dof.tangent"),
            "dof_res": wall(r, f"{scs}/formUnbalance/dof.residual"),
            "inertia_tan": wall(r, f"{scs}/formTangent/elem.tangent/brick.inertia"),
            "inertia_res": wall(r, f"{scs}/formUnbalance/elem.residual/brick.inertia"),
            "elem_tan": wall(r, f"{scs}/formTangent/elem.tangent"),
            "solve": wall(r, f"{scs}/linearSolve"),
        }
        d["headroom"] = (d["dof_tan"] + d["dof_res"]
                         + d["inertia_tan"] + d["inertia_res"])
        d["step_ns"] = step
        return d


def main():
    jp = HERE / "t2_nodal_mass.json"
    meta = json.loads(jp.read_text()) if jp.exists() else {}
    n = meta.get("n", 15)
    print(f"ADR-77 T2 — G2 headroom vs where the mass lives "
          f"(n={n}, PARDISO, {meta.get('threads','?')} threads)")
    print(f"G2 gate: headroom >= {GATE:.0f}% of step wall\n")
    print(f"{'arm':<30} {'step_ms':>9} {'dof.tan':>9} {'dof.res':>9} "
          f"{'inertia':>9} {'HEADROOM':>10} {'verdict':>9}")
    rows = []
    for key, label in ARMS:
        d = read(key, n)
        if not d:
            continue
        pct = lambda v: 100.0 * v / d["step_ns"]
        hp = pct(d["headroom"])
        rows.append((key, label, d, hp))
        print(f"{label:<30} {d['step_ms']:>9.1f} {pct(d['dof_tan']):>8.2f}% "
              f"{pct(d['dof_res']):>8.2f}% "
              f"{pct(d['inertia_tan']+d['inertia_res']):>8.2f}% "
              f"{hp:>9.2f}% {'PASS' if hp >= GATE else 'fail':>9}")

    if not rows:
        print("no profiles found")
        return

    print("\n=== does moving mass to the nodes reopen G2? ===")
    base = next((r for r in rows if r[0] == "material"), None)
    for key, label, d, hp in rows:
        if key == "material":
            continue
        if base:
            dt_ratio = (d["dof_tan"] / base[2]["dof_tan"]) if base[2]["dof_tan"] else float("inf")
            print(f"  {label}: dof.tangent {100.0*d['dof_tan']/d['step_ns']:.2f}% "
                  f"vs {100.0*base[2]['dof_tan']/base[2]['step_ns']:.2f}% "
                  f"({dt_ratio:.1f}x the absolute time)")
        print(f"     headroom {hp:.2f}% vs gate {GATE:.0f}%  -> "
              f"{'G2 REOPENS' if hp >= GATE else 'G2 STAYS CLOSED'}")

    worst = max(r[3] for r in rows)
    print(f"\n  best headroom across all arms: {worst:.2f}%  "
          f"({'>=' if worst >= GATE else '<'} {GATE:.0f}% gate)")
    print("  => G2 " + ("AUTHORIZED — design the cache"
                        if worst >= GATE else
                        "CLOSED for good — the nodal-mass caveat is settled"))


if __name__ == "__main__":
    main()
