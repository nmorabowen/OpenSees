"""ADR-77 T1 — cross-nonlinearity view of the reuse levers.

Aggregates t1_reuse_levers_L{1,4,10}.json into the deliverable T1 actually owes:
not "which algorithm is fastest" (that depends on the deck) but "which lever at
what nonlinearity", plus the per-arm assembly/solve split that says WHY.

The load sweep exists because the T0 deck is only 1% past yield (peak von Mises
252.5 vs s0 250), where K_t ~= K_i and the initial-stiffness arms flatter
themselves. x4 = 2.4x yield, x10 = 7x yield.

Run with the h5py interpreter:
  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe t1_report.py
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
LOADS = [1, 4, 10]
ORDER = ["newton", "newton_init", "modnewton", "modnewton_init",
         "modnewton_init_fo", "krylov"]
VM = {1: "252 MPa (1% past yield)", 4: "612 MPa (2.4x yield)",
      10: "1773 MPa (7x yield)"}


def wall(root, path):
    node = root
    for part in path.split("/"):
        if part not in node:
            return 0
        node = node[part]
    return int(node.attrs.get("wall_ns", 0))


def split(key, load, n=15):
    """(assembly%, solve%, derived?) or None.

    INSTRUMENTATION GAP, found by T1 (2026-07-26): only THREE EquiSolnAlgo
    implementations carry P3 profile scopes -- Linear, ModifiedNewton,
    NewtonRaphson. Eight others (KrylovNewton, BFGS, Broyden, NewtonLineSearch,
    AcceleratedNewton, ExpressNewton, NewtonHallM, PeriodicNewton) have ZERO, so
    their steps have no formTangent/formUnbalance/linearSolve wrappers and the
    solve time silently lands in the unattributed remainder of solveCurrentStep.
    Reporting that as "0.0% solve" would be a fabrication.

    For the uninstrumented arms the leaf scopes (elem.tangent, elem.residual,
    dof.*) ARE still emitted -- they live in the integrator, not the algorithm --
    so assembly is measured directly and solve is DERIVED BY DIFFERENCE as
    solveCurrentStep minus its accounted children. Flagged so the two kinds of
    number are never silently mixed.
    """
    p = HERE / f"t1_{key}_n{n}_L{load}.h5"
    if not p.exists():
        return None
    with h5py.File(p, "r") as f:
        if f"runs/{key}/rollup/root" not in f:
            return None
        r = f[f"runs/{key}/rollup/root"]
        step = int(r["step"].attrs.get("wall_ns", 0))
        if not step:
            return None
        scs = "step/solveCurrentStep"
        if wall(r, f"{scs}/formTangent") or wall(r, f"{scs}/linearSolve"):
            asm = (wall(r, f"{scs}/formTangent/elem.tangent")
                   + wall(r, f"{scs}/formUnbalance/elem.residual")
                   + wall(r, f"{scs}/formTangent/dof.tangent")
                   + wall(r, f"{scs}/formUnbalance/dof.residual"))
            return (100.0 * asm / step, 100.0 * wall(r, f"{scs}/linearSolve") / step, False)
        # uninstrumented algorithm: flat leaves, solve by difference
        asm = (wall(r, f"{scs}/elem.tangent") + wall(r, f"{scs}/elem.residual")
               + wall(r, f"{scs}/dof.tangent") + wall(r, f"{scs}/dof.residual"))
        total = wall(r, scs)
        accounted = asm + wall(r, f"{scs}/elem.update")
        return (100.0 * asm / step, 100.0 * max(0, total - accounted) / step, True)


def main():
    data = {}
    for L in LOADS:
        p = HERE / f"t1_reuse_levers_L{L}.json"
        if p.exists():
            data[L] = json.loads(p.read_text())
    if not data:
        print("no t1_reuse_levers_L*.json found")
        return

    print("ADR-77 T1 — reuse levers vs nonlinearity (n=15, PARDISO, 4 threads)")
    for L, d in data.items():
        print(f"\n=== load x{L} — peak von Mises {VM.get(L, '?')} ===")
        base = d["arms"].get("newton", {}).get("wall_s")
        print(f"{'algorithm':<38} {'wall_s':>8} {'vs Newton':>10} {'it/step':>8} "
              f"{'ms/iter':>8} {'asm%':>6} {'solve%':>7}  answer")
        for key in ORDER:
            a = d["arms"].get(key)
            if not a:
                continue
            if a.get("failed"):
                # "failed" here means it hit maxIter, which for the
                # initial-stiffness family is SLOW CONVERGENCE, not divergence
                # -- see the maxIter probes in the results doc. Do not read
                # this row as "cannot converge".
                print(f"{a['label']:<38} {'--':>8} {'hit maxIter':>10}")
                continue
            sp = f"{base / a['wall_s']:.2f}x" if base else "-"
            sl = split(key, L)
            asm = f"{sl[0]:.1f}" if sl else "-"
            slv = (f"{sl[1]:.1f}~" if sl and sl[2] else f"{sl[1]:.1f}") if sl else "-"
            print(f"{a['label']:<38} {a['wall_s']:>8.3f} {sp:>10} "
                  f"{a['iters_per_step']:>8.2f} {a['ms_per_iter']:>8.1f} "
                  f"{asm:>6} {slv:>7}  {a['answer']}")

    # the cross-load story: does the winner change with nonlinearity?
    print("\n=== winner by load (wall clock) ===")
    for L, d in data.items():
        ok = [(k, a) for k, a in d["arms"].items() if a.get("wall_s")]
        if not ok:
            continue
        k, a = min(ok, key=lambda kv: kv[1]["wall_s"])
        base = d["arms"].get("newton", {}).get("wall_s")
        sp = f"{base / a['wall_s']:.2f}x vs Newton" if base else ""
        print(f"  load x{L:<3} -> {a['label']:<38} {sp}")

    print("\n=== iteration cost vs iteration count (the whole tradeoff) ===")
    print("if ms/iter falls faster than it/step rises, the cheap-iteration arm wins")
    for L, d in data.items():
        nw = d["arms"].get("newton")
        if not nw or not nw.get("wall_s"):
            continue
        for key in ORDER[1:]:
            a = d["arms"].get(key)
            if not a or not a.get("wall_s"):
                continue
            ir = a["iters_per_step"] / nw["iters_per_step"]
            cr = a["ms_per_iter"] / nw["ms_per_iter"]
            print(f"  x{L:<3} {a['label']:<36} iters {ir:5.2f}x  cost/iter "
                  f"{cr:5.2f}x  net {ir * cr:5.2f}x")


if __name__ == "__main__":
    main()
