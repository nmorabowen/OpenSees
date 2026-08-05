"""ADR-75 P1e gate — is `-krylov` (iparm[3] CGS) actually faster than direct?

Measures Lane B under FULL NEWTON, which is the regime the `factored` flag
cannot help: A changes at every iteration, so the P1 reuse gate never fires and
every solve pays a full numeric factorization. `-krylov` is the only lever that
reaches this case.

Two axes matter and pull in OPPOSITE directions, so both are swept:
  * MODEL SIZE. Factorization is ~O(N^1.5-2) while a CGS iteration is a sparse
    matvec + two triangular solves, ~O(nnz). The bigger the model, the more a
    handful of iterations should beat a refactorization.
  * THREADS. Numeric factorization parallelizes well; triangular solves do not.
    So threading makes the DIRECT path cheaper relative to CGS, and a win at
    1 thread can evaporate at 4. Measuring only one thread count would be the
    ADR-40b mistake (a dominance verdict that does not survive re-measurement).

Metric: wall seconds for the full NSTEPS analyze(), median of interleaved
rounds (A/B/A/B...), per PHASE0_RECIPE.md. Every configuration must reproduce
the direct-solver tip displacement or its timing is meaningless.

Run (per thread count):
  set MKL_NUM_THREADS=1 && python3.12 -S p1e_krylov_bench.py
  set MKL_NUM_THREADS=4 && python3.12 -S p1e_krylov_bench.py
Env: OPS_BENCH_SIZES (comma list of grid n), OPS_BENCH_REPEATS.
"""
from __future__ import annotations
import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("MKL_NUM_THREADS", "1"))
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees"
    r"\opensees-pr-strategy-9c17b8\dist\bin")
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"

HERE = Path(__file__).resolve().parent
THREADS = os.environ["MKL_NUM_THREADS"]
REPEATS = int(os.environ.get("OPS_BENCH_REPEATS", "3"))
SIZES = [int(s) for s in os.environ.get("OPS_BENCH_SIZES", "15,20,25").split(",")]
# NSTEPS matters more than it looks: CGS never refactors while it keeps winning,
# so the preconditioner AGES — it is still the step-1 factorization. A short run
# can therefore overstate the steady-state win. OPS_BENCH_NSTEPS sweeps it.
NSTEPS = int(os.environ.get("OPS_BENCH_NSTEPS", "15"))

E, nu = 200000.0, 0.3
Kbulk = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, bpar, Hiso = 250.0, 0.0, 0.0, 2000.0
LEN = 100.0

CONFIGS = [
    ("direct",   ["Pardiso"]),                  # the P1/P1d baseline
    ("krylov3",  ["Pardiso", "-krylov", 3]),    # loose  (eps 1e-3)
    ("krylov6",  ["Pardiso", "-krylov", 6]),    # Intel's own example value
    ("krylov8",  ["Pardiso", "-krylov", 8]),    # tight — expect more fallbacks
    # Do the two reuse levers COMPOSE? P1d's half-storage is worth ~1.25x on its
    # own and is a different mechanism (less data to factor, not fewer
    # factorizations), so the product is the number a symmetric production deck
    # would actually see. Legal only for -matrixType 1: Intel's K=2 is
    # documented for symmetric POSITIVE DEFINITE, and mtype -2 has no CGS mode.
    ("spd",         ["Pardiso", "-matrixType", 1]),
    ("spd+krylov6", ["Pardiso", "-matrixType", 1, "-krylov", 6]),
]
# The big rungs cost minutes per run; OPS_BENCH_CONFIGS trims to the two that
# decide the question ("direct,krylov6") without re-paying for the L sweep.
_ONLY = [s.strip() for s in os.environ.get("OPS_BENCH_CONFIGS", "").split(",") if s.strip()]
if _ONLY:
    CONFIGS = [(n, a) for (n, a) in CONFIGS if n in _ONLY]


def build(system_args, n):
    nx = ny = nz = n

    def _nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, Kbulk, G, "-iso", "voce", s0, Qinf, bpar, Hiso)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(_nid(i, j, k), i * LEN, j * LEN, k * LEN)
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
    ops.algorithm("Newton")            # FULL Newton — the target regime
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")
    return _nid(nx // 2, ny // 2, nz)


def run_once(system_args, n):
    tip = build(system_args, n)
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"analyze failed at step {s + 1}")
    wall = time.perf_counter() - t0
    ux = ops.nodeDisp(tip, 1)
    if not (ux == ux) or abs(ux) == float("inf"):
        raise RuntimeError("non-finite tip displacement")
    return wall, ux


def main():
    print(f"P1e -krylov bench | MKL_NUM_THREADS={THREADS} | REPEATS={REPEATS} | "
          f"sizes {SIZES} | {NSTEPS} steps, full Newton")
    out = {"threads": THREADS, "repeats": REPEATS, "nsteps": NSTEPS,
           "pyd": ops.__file__, "sizes": {}}

    for n in SIZES:
        ndof = 3 * (n + 1) ** 2 * n     # base plane fully fixed
        print(f"\n=== grid {n}^3  ({n**3} bricks, ~{ndof/1000:.1f}k DOF) ===")
        samples = {name: [] for name, _ in CONFIGS}
        tips = {}
        try:
            for r in range(REPEATS):
                for name, args in CONFIGS:
                    wall, ux = run_once(args, n)
                    samples[name].append(wall)
                    tips[name] = ux
        except Exception as e:  # noqa: BLE001
            print(f"  size {n} FAILED: {e}")
            out["sizes"][str(n)] = {"failed": str(e)}
            continue

        base = statistics.median(samples["direct"])
        ref = tips["direct"]
        rows = {}
        print(f"  {'config':10s} {'median_s':>10s} {'vs direct':>10s}   tip rel err")
        for name, _ in CONFIGS:
            med = statistics.median(samples[name])
            err = abs(tips[name] - ref) / abs(ref)
            rows[name] = {"median_s": med, "ratio_vs_direct": med / base,
                          "tip_ux": tips[name], "tip_rel_err": err,
                          "samples": sorted(samples[name])}
            flag = "" if err < 1e-7 else "  *** MISMATCH ***"
            print(f"  {name:10s} {med:10.3f} {med/base:9.3f}x   {err:.2e}{flag}")
        out["sizes"][str(n)] = {"ndof": ndof, "results": rows}

    # Size list in the name: a later run with different OPS_BENCH_SIZES used to
    # CLOBBER an earlier one at the same thread count (the 25^3 rung was lost to
    # the composition sweep exactly once — the numbers survived only in the
    # console log).
    dest = HERE / f"p1e_krylov_t{THREADS}_n{'-'.join(str(s) for s in SIZES)}.json"
    dest.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
