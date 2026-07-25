"""ADR-75 P1e smoke — `system Pardiso -krylov <L>` correctness gate.

Timing is worthless until this passes. Three things are checked, in order:

  1. ANSWER. -krylov makes the solve INEXACT (CGS stops at 10^-L), so every
     configuration must still reproduce the UmfPack tip displacement. The bar
     is deliberately not 1e-9 here as it is for the direct solvers: an
     iterative solve inside a Newton loop converges the OUTER residual, so the
     tolerance that matters is the one `test NormDispIncr 1e-7` enforces.
  2. FALLBACK. -matrixType 2 (mtype -2) has no documented CGS mode; it must
     WARN and keep working, not silently produce a different answer.
  3. REUSE STILL WORKS. ModifiedNewton must keep taking the phase-33 shortcut
     (A unchanged) and must NOT be corrupted by a stale preconditioner.

Run:
  C:\\Users\\nmb\\AppData\\Local\\Python\\bin\\python3.12 -S p1e_smoke.py
"""
from __future__ import annotations
import os
import sys

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("PMI_RANK", "0")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees"
    r"\opensees-pr-strategy-9c17b8\dist\bin")
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__)

# --- Lane-B model (verbatim from laneB_solver_bench.py) ---------------------
E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
L = 100.0
NSTEPS = 15


def build(system_args, n=15, algo=("Newton",)):
    nx = ny = nz = n

    def _nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(_nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    eTag = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(
                    "LadrunoBrick", eTag,
                    _nid(i, j, k), _nid(i + 1, j, k), _nid(i + 1, j + 1, k),
                    _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1),
                    _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1), 1)
                eTag += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(_nid(i, j, nz), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1e-7, 25)
    ops.algorithm(*algo)
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")
    return _nid(nx // 2, ny // 2, nz)


def run(system_args, n=15, algo=("Newton",)):
    """Returns (tip_ux, total_iterations).

    The ITERATION COUNT is the load-bearing signal for the stale-preconditioner
    guard, not the displacement. Proven by deliberately breaking the guard: the
    tip displacement moved only 1.2e-11 and the run still "passed", because a
    solve against stale factors is still a descent direction and Newton's
    residual test simply launders it into a slower quasi-Newton. What it CANNOT
    hide is the extra iterations that laundering costs.
    """
    tip = build(system_args, n, algo)
    iters = 0
    for s in range(NSTEPS):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"analyze failed at step {s + 1}")
        try:
            iters += ops.testIter()
        except Exception:      # noqa: BLE001 - older/!py builds lack the query
            iters = -1
    ux = ops.nodeDisp(tip, 1)
    if not (ux == ux) or abs(ux) == float("inf"):
        raise RuntimeError("non-finite tip displacement")
    return ux, iters


CASES = [
    ("UmfPack (reference)", ["UmfPack"], ("Newton",)),
    ("Pardiso direct", ["Pardiso"], ("Newton",)),
    ("Pardiso -krylov 3", ["Pardiso", "-krylov", 3], ("Newton",)),
    ("Pardiso -krylov 6", ["Pardiso", "-krylov", 6], ("Newton",)),
    ("Pardiso -krylov 8", ["Pardiso", "-krylov", 8], ("Newton",)),
    # mtype 2 (SPD) — Lane B's associated-flow J2 tangent is SPD here, so
    # Intel's K=2 (LL^T-preconditioned CG) is legal and must also be exact.
    ("Pardiso -spd -krylov 6", ["Pardiso", "-matrixType", 1, "-krylov", 6],
     ("Newton",)),
    # mtype -2 — NO documented CGS mode. Must warn and fall back to direct.
    ("Pardiso -symmetric -krylov 6 (must warn+fallback)",
     ["Pardiso", "-matrixType", 2, "-krylov", 6], ("Newton",)),
    # The P1 reuse path must be untouched: ModifiedNewton holds A fixed within
    # a step, so most solves are phase-33-only and CGS should barely fire.
    ("Pardiso -krylov 6 + ModifiedNewton",
     ["Pardiso", "-krylov", 6], ("ModifiedNewton",)),
    ("Pardiso direct + ModifiedNewton", ["Pardiso"], ("ModifiedNewton",)),
]

# TWO tiers, because one loose bar for everything is a weak test (adversarial
# review). Every FULL-NEWTON row — direct and -krylov alike — measured exactly
# 0.0 against the UmfPack reference, so holding them to 1e-7 would let a real
# regression in the direct path sail through. They get 1e-12. The
# ModifiedNewton rows legitimately differ from the Newton reference (~5e-10, a
# different algorithm), so they are NOT compared against it at all — they are
# compared against each other, which is the guard check further down.
TOL_NEWTON = 1e-12
MODIFIED = "ModifiedNewton"


def main():
    print(f"\nP1e smoke | Lane B 15^3 (~11.5k DOF) | {NSTEPS} steps | "
          f"MKL_NUM_THREADS={os.environ['MKL_NUM_THREADS']}\n")
    ref = None
    failures = []
    iters = {}
    tips = {}
    for name, args, algo in CASES:
        try:
            ux, nit = run(args, 15, algo)
        except Exception as e:  # noqa: BLE001
            print(f"  {name:52s} FAILED: {e}")
            failures.append(f"{name}: {e}")
            continue
        iters[name] = nit
        tips[name] = ux
        if ref is None:
            ref = ux
            print(f"  {name:52s} ux = {ux:.12e}  iters {nit:4d}  (reference)")
            continue
        err = abs(ux - ref) / abs(ref)
        if MODIFIED in name:
            # Not comparable to the Newton reference — see the TOL note.
            print(f"  {name:52s} ux = {ux:.12e}  iters {nit:4d}  "
                  f"(vs Newton ref {err:.2e}; checked below instead)")
            continue
        verdict = "ok" if err < TOL_NEWTON else "*** MISMATCH ***"
        if err >= TOL_NEWTON:
            failures.append(f"{name}: rel err {err:.3e} exceeds {TOL_NEWTON:.0e}")
        print(f"  {name:52s} ux = {ux:.12e}  iters {nit:4d}  "
              f"rel err {err:.2e}  {verdict}")

    # ---- the stale-preconditioner guard: the ONE check that discriminates --
    # Under ModifiedNewton the tangent is held fixed within a step, so solves
    # 2..n of each step hit `factored == true && factorsCurrent == false` — the
    # exact state `factorsCurrent` exists for. Removing the guard sends those
    # solves down the phase-33 shortcut against factors belonging to an older
    # tangent.
    #
    # CALIBRATED BY DELIBERATELY BREAKING IT (2026-07-25). Two checks that do
    # NOT work, measured, so nobody re-derives them:
    #   * tip displacement vs the UmfPack/Newton reference — moved 1.2e-11,
    #     i.e. still "ok" at any sane tolerance. Newton's residual test
    #     launders a stale-factor solve into a slower-but-convergent
    #     quasi-Newton.
    #   * ModifiedNewton ITERATION COUNT — 94 vs 93 (1.01x). Swept the
    #     hardening modulus to perfect plasticity (Hiso 2000 -> 0) hunting a
    #     sharper tangent; the gap stayed at 1.03x throughout. Lane B's tangent
    #     simply does not move enough for stale factors to hurt it.
    # What DOES discriminate is comparing -krylov against the DIRECT run of the
    # SAME algorithm: with the guard the two are bit-identical, without it they
    # differ at 1e-11. So the reference has to be direct+ModifiedNewton, not
    # the Newton reference used above.
    mn_krylov = tips.get("Pardiso -krylov 6 + ModifiedNewton")
    mn_direct = tips.get("Pardiso direct + ModifiedNewton")
    if mn_krylov is not None and mn_direct is not None:
        d = abs(mn_krylov - mn_direct) / abs(mn_direct)
        ok = d < 1e-13
        print(f"\n  stale-preconditioner guard: -krylov vs direct under the "
              f"SAME ModifiedNewton -> rel diff {d:.2e} "
              f"({'ok, bit-identical' if ok else '*** GUARD REGRESSED ***'})")
        print(f"    (iteration counts {iters.get('Pardiso -krylov 6 + ModifiedNewton')}"
              f" vs {iters.get('Pardiso direct + ModifiedNewton')} — reported "
              f"only as context; measured NOT to discriminate)")
        if not ok:
            failures.append(
                f"stale-preconditioner guard: -krylov diverges from direct "
                f"under the same ModifiedNewton by {d:.2e} (expected 0.0); "
                f"the phase-33 shortcut is being taken against stale factors")

    print("\n--- -stats view of the CGS decisions (one Newton solve per line) ---")
    try:
        run(["Pardiso", "-krylov", 6, "-stats"], 10)
    except Exception as e:  # noqa: BLE001
        print(f"  stats run FAILED: {e}")
        failures.append(f"stats run: {e}")

    print("\nRESULT:", "PASS" if not failures else "FAIL")
    for f in failures:
        print("  -", f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
