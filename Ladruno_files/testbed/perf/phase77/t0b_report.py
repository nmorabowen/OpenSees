"""ADR-77 T0b — the step-anatomy-vs-threads table.

Answers the question T0 left open: at 1 thread `linearSolve` was still the
largest loop on the transient step (45.2% vs `elem.tangent` 40.8%), unlike the
L3-0 static lane where the element loop reached 74.9% under PARDISO. Is that a
real property of the transient step, or a 1-thread artifact?

Reads the per-thread profiles written by t0b_thread_sweep.sh and prints, for
each arm, how the split moves with MKL_NUM_THREADS. The crossover row -- where
elem.tangent overtakes linearSolve -- is what decides whether ADR-77's
assembly-side premise is restored.

UmfPack is the control: a serial factorization. If its split also moves with
threads, the mover is BLAS inside the element kernels, not the solver.

Run with the h5py interpreter:
  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe t0b_report.py [n]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
N = sys.argv[1] if len(sys.argv) > 1 else "15"
LADDER = [int(t) for t in (sys.argv[2].split(",") if len(sys.argv) > 2
                           else ["1", "2", "4", "8", "16"])]


def wall(root, path):
    node = root
    for part in path.split("/"):
        if part not in node:
            return None
        node = node[part]
    return int(node.attrs.get("wall_ns", 0))


def read(arm, n, t):
    p = HERE / f"t0_{arm}_n{n}_t{t}.h5"
    if not p.exists():
        return None
    with h5py.File(p, "r") as f:
        root = f[f"runs/{arm}/rollup/root"]
        step = int(root["step"].attrs.get("wall_ns", 0))
        if not step:
            return None
        g = lambda pth: (wall(root, pth) or 0)
        return {
            "step_ms": step / 1e6,
            # absolute ms matters as much as the share: element kernels work on
            # 24x24 matrices, so MKL should run them serially and elem.tangent_ms
            # should be ~thread-INVARIANT. If it is not, threaded BLAS inside the
            # kernels is confounding the comparison and the shares cannot be read
            # as "the solver got faster".
            "elem_tangent_ms": (wall(root, "step/solveCurrentStep/formTangent/elem.tangent") or 0) / 1e6,
            "linear_solve_ms": (wall(root, "step/solveCurrentStep/linearSolve") or 0) / 1e6,
            "elem_tangent": 100.0 * g("step/solveCurrentStep/formTangent/elem.tangent") / step,
            "linear_solve": 100.0 * g("step/solveCurrentStep/linearSolve") / step,
            "elem_residual": 100.0 * g("step/solveCurrentStep/formUnbalance/elem.residual") / step,
            "dof_loops": 100.0 * (g("step/solveCurrentStep/formTangent/dof.tangent")
                                  + g("step/solveCurrentStep/formUnbalance/dof.residual")) / step,
            "update": 100.0 * g("step/solveCurrentStep/update") / step,
            "assembly_total": 100.0 * (
                g("step/solveCurrentStep/formTangent/elem.tangent")
                + g("step/solveCurrentStep/formUnbalance/elem.residual")
                + g("step/solveCurrentStep/formTangent/dof.tangent")
                + g("step/solveCurrentStep/formUnbalance/dof.residual")) / step,
        }


def wall_seconds(n, t, arm_key):
    p = HERE / f"t0_step_anatomy_t{t}.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())["arms"][arm_key]["wall_s"]
    except (KeyError, ValueError):
        return None


def main():
    print(f"ADR-77 T0b — step anatomy vs threads | n={N}")
    print("(L3-0 static lane, for contrast: element 35.8% -> 74.9% under PARDISO)")
    for arm, key in (("pardiso", "Pardiso"), ("umfpack", "UmfPack")):
        rows = [(t, read(arm, N, t)) for t in LADDER]
        rows = [(t, r) for t, r in rows if r]
        if not rows:
            continue
        ctl = " (CONTROL: serial factorization)" if arm == "umfpack" else ""
        print(f"\n=== {arm.upper()}{ctl} ===")
        print(f"{'thr':>4} {'wall_s':>8} {'step_ms':>9} {'elem.tan':>9} {'(ms)':>8} "
              f"{'lin.solve':>10} {'(ms)':>8} {'dof':>6} {'ASSEMBLY':>9}  verdict")
        for t, r in rows:
            ws = wall_seconds(N, t, key)
            v = "elem>solve" if r["elem_tangent"] > r["linear_solve"] else "solve>elem"
            print(f"{t:>4} {(('%8.3f' % ws) if ws else '       -')} {r['step_ms']:>9.1f} "
                  f"{r['elem_tangent']:>8.1f}% {r['elem_tangent_ms']:>8.0f} "
                  f"{r['linear_solve']:>9.1f}% {r['linear_solve_ms']:>8.0f} "
                  f"{r['dof_loops']:>5.2f}% {r['assembly_total']:>8.1f}%  {v}")
        # Thread-invariance check on the assembly side, reported over two regions.
        # MEASURED 2026-07-26: elem.tangent ms is flat to 4 threads and then RISES
        # at 8/16 -- in BOTH arms identically, which rules out a solver-specific
        # cause. That is MKL worker threads (spin-wait by default) contending with
        # the serial element loop on a 8P+8E core part, not threaded BLAS speeding
        # the kernels up. Consequence: the shares in the <=4 rows are clean, the
        # 8/16 rows are contaminated and must not be read as "assembly got worse".
        def spread(sub):
            v = [r["elem_tangent_ms"] for t, r in rows if t in sub]
            return (max(v) - min(v)) / min(v) * 100.0 if len(v) > 1 and min(v) > 0 else None
        lo, hi = spread({1, 2, 4}), spread({1, 2, 4, 8, 16})
        if lo is not None:
            flag = "OK -- assembly thread-invariant here, shares are clean" if lo < 15 \
                   else "** assembly moved even at low thread counts -- investigate **"
            print(f"  -> elem.tangent ms spread, 1-4 threads: {lo:.1f}%  {flag}")
        if hi is not None:
            print(f"  -> elem.tangent ms spread, full ladder:  {hi:.1f}%  "
                  f"(rises at 8/16 = MKL spin-wait contention with the serial "
                  f"element loop; same in both arms)")
        cross = [t for t, r in rows if r["elem_tangent"] > r["linear_solve"]]
        if cross:
            print(f"  -> crossover at {min(cross)} thread(s): the element tangent "
                  f"overtakes the solve.")
        else:
            print("  -> NO crossover on this ladder: the solve stays the largest loop.")
        # ADR-75b's threading gate was >40% for the loop under attack
        g40 = [t for t, r in rows if r["elem_tangent"] > 40.0]
        if g40:
            print(f"  -> elem.tangent clears ADR-75b's >40% threading gate at "
                  f"{min(g40)}+ thread(s).")


if __name__ == "__main__":
    main()
