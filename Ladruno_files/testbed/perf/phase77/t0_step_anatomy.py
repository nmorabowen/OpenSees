"""ADR-77 T0 — anatomy of one IMPLICIT TRANSIENT step.

Every ADR-75 number is statics-shaped. The transient path is different: it
assembles through loops the static path does not have. Confirmed by reading the
source (2026-07-26):

  IncrementalIntegrator::formTangent   -> elem.tangent ONLY (no DOF_Group loop)
  TransientIntegrator::formTangent     -> dof.tangent + elem.tangent

so the nodal mass/damping assembly exists only here, and until this ADR's T0
instrumentation it was UNTIMED (the P3 deep hook covers the element loop alone).
This script measures the five-loop split that gap was hiding:

    dof.tangent    nodal M/C into A          (transient-only, new scope)
    elem.tangent   element c1*K+c2*C+c3*M    (deep, per-classTag)
    elem.residual  element resisting force   (deep, per-classTag)
    dof.residual   nodal inertia/damping/load into B   (new scope)
    linearSolve    factor + back-solve
    update         Domain::update (state update / material return maps)

Gate G0 forbids transplanting the L3-0 statics profile (element 35.8% -> 74.9%
under PARDISO) onto the transient step. This measures it instead, on both
solvers, so the flip is observed rather than assumed.

Bench stack per ADR-77 §5: system Pardiso + numberer LadrunoParallel + P3
profiler. Oracle runs are 1-thread (threaded PARDISO is not byte-reproducible,
ADR-75 P1f).

Run:
  OPS_PYD=<worktree>\\dist\\bin MKL_NUM_THREADS=1 python3.12 -S t0_step_anatomy.py

Env knobs:
  OPS_T0_N        cube edge elements (default 15  -> ~11.5k DOF)
  OPS_T0_STEPS    transient steps timed (default 6)
  OPS_T0_REPEATS  repeats per arm, median reported (default 3)
  OPS_T0_SOLVERS  comma list subset of "UmfPack,Pardiso"
  OPS_T0_NUMBERER "RCM" (default)

NUMBERER GOTCHA (measured 2026-07-26, ADR-77 T0): the fork's parallel numberer is
**not available in a serial build**. The registered verbs are `LadrunoParallelRCM`
/ `LadrunoParallelPlain` (ADR-30 §140's "numberer LadrunoParallel" is doc
shorthand and is rejected outright), and BOTH factories are wrapped in
`#ifdef _PARALLEL_INTERPRETERS` — `OpenSeesCommands.cpp:4774` returns 0 in serial.
A null numberer makes OpenSees fall back to RCM with only a warning, so a serial
bench that asks for the Ladruno numberer silently measures RCM. Same family as
the ADR-75 P1d "mistyped `system` costs you the SOLVER silently" trap. Serial T0
therefore runs RCM by design, and the Ladruno numberer is a T4 (parallel) knob.
"""
from __future__ import annotations

import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", os.environ["MKL_NUM_THREADS"])
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees"
    r"\.claude\worktrees\pardisio-profiling-0a03b1\dist\bin",
)
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402

HERE = Path(__file__).resolve().parent
N = int(os.environ.get("OPS_T0_N", "15"))
NSTEPS = int(os.environ.get("OPS_T0_STEPS", "6"))
REPEATS = int(os.environ.get("OPS_T0_REPEATS", "3"))
NUMBERER = os.environ.get("OPS_T0_NUMBERER", "RCM")
THREADS = os.environ["MKL_NUM_THREADS"]

# Same material/geometry constants as the ADR-75 Lane-B deck so the two studies
# are comparable; only the analysis type and the added density differ.
E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
RHO = 7.85e-9          # steel, N-mm-s units
L = 100.0
DT = 2.0e-4


def build(n, system_args, numberer):
    """Lane-B cube, transient. Returns (tipNode, nele, ndof)."""
    nx = ny = nz = n

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso,
                   "-rho", RHO)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(
                    "LadrunoBrick", e,
                    nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                    nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1),
                    1)
                e += 1
    # A step load released into a dynamic response: Newton must actually iterate.
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), 4.0e5, 0.0, -4.0e5)

    ops.constraints("Plain")
    ops.numberer(numberer)
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    ndof = (nnode - (nx + 1) * (ny + 1)) * 3
    return nid(nx // 2, ny // 2, nz), e - 1, ndof


def run_once(n, system_args, numberer, profile_tag=None):
    tip, nele, ndof = build(n, system_args, numberer)
    if profile_tag:
        ops.profiler("start", "-deep", "-perStep")
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1, DT) != 0:
            raise RuntimeError(f"analyze failed at step {s+1}")
    wall = time.perf_counter() - t0
    if profile_tag:
        ops.profiler("stop")
        out = HERE / f"t0_{profile_tag}_n{n}.h5"
        ops.profiler("report", str(out), "-run", profile_tag)
        ops.profiler("reset")
    return wall, ops.nodeDisp(tip, 3), nele, ndof


def main():
    solvers = [s.strip() for s in
               os.environ.get("OPS_T0_SOLVERS", "UmfPack,Pardiso").split(",") if s.strip()]
    print(f"ADR-77 T0 step anatomy | n={N} steps={NSTEPS} repeats={REPEATS} "
          f"dt={DT} numberer={NUMBERER} MKL_NUM_THREADS={THREADS}")
    out = {"n": N, "steps": NSTEPS, "repeats": REPEATS, "dt": DT,
           "numberer": NUMBERER, "threads": THREADS, "arms": {}}
    ref = None
    for name in solvers:
        samples, uz = [], None
        nele = ndof = None
        for r in range(REPEATS):
            w, uz, nele, ndof = run_once(N, [name], NUMBERER)
            samples.append(w)
        med = statistics.median(samples)
        # one extra profiled run per arm (profiling perturbs timing, so it is
        # never one of the timed samples)
        run_once(N, [name], NUMBERER, profile_tag=name.lower())
        if ref is None:
            ref = uz
        rel = abs(uz - ref) / abs(ref) if ref else 0.0
        out["arms"][name] = {"wall_s": med, "samples": sorted(samples),
                             "uz": uz, "rel_err_vs_first": rel,
                             "nele": nele, "ndof": ndof}
        print(f"  {name:<9} wall={med:8.3f}s  uz={uz: .9e}  relerr={rel:.1e}  "
              f"nele={nele} ndof={ndof}")
    if "UmfPack" in out["arms"] and "Pardiso" in out["arms"]:
        spd = out["arms"]["UmfPack"]["wall_s"] / out["arms"]["Pardiso"]["wall_s"]
        out["speedup_pardiso"] = spd
        print(f"  speedup (Umf/Pardiso) = {spd:.2f}x")
    (HERE / "t0_step_anatomy.json").write_text(json.dumps(out, indent=2))
    print(f"wrote {HERE / 't0_step_anatomy.json'} + per-arm .h5 profiles")


if __name__ == "__main__":
    main()
