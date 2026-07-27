"""ADR-77 T2 — the G2 headroom on a NODAL-MASS deck (the one caveat left open).

T0 killed G2 (the constant-M/C tangent cache) at 3.2-4.4% headroom against a
>=10% gate, and T0b confirmed it fails at every thread count (4.43% -> 6.85%).
But T0 §3 recorded one caveat that could reopen it, and this stage exists to
settle it rather than let it sit:

    that deck carries its mass in the MATERIAL (`LadrunoJ2 -rho`), so the
    DOF_Group tangent loop measured almost nothing (0.27% of step). A model
    with LUMPED NODAL MASS (`mass` commands -- i.e. most building models) puts
    real work in that loop.

So the experiment holds geometry, material, load, solver and step count fixed
and moves ONLY where the mass lives. Anything that changes is attributable to
that one variable; that is the whole design.

  arm A  material mass   `-rho 7.85e-9`, no `mass` commands   (= the T0 deck)
  arm B  nodal mass      `-rho 0` + `mass` on every free node, same total mass
  arm C  nodal + Rayleigh(alphaM, betaK0)   C is CONSTANT here
  arm D  nodal + Rayleigh(alphaM, betaKcurr) C is STATE-DEPENDENT here

C and D exist as a PAIR because they differ on the cache's own precondition, and
an earlier version of this file conflated them (it set betaKinit while its
comment claimed betaK -- corrected 2026-07-26). ADR-76 §5.1: with betaK on the
CURRENT stiffness, C = alphaM*M + betaK*Kt and Kt changes every iteration as the
material yields, so C is NOT cacheable. With betaK0 on the INITIAL stiffness,
C = alphaM*M + betaK0*Ki, and LadrunoBrick CACHES Ki, so C is fully constant.
Measuring only one of them would answer only half the question.

Why Rayleigh matters here at all (checked in source, Element.cpp): getDamp()
calls `addMatrix(0.0, this->getMass(), alphaM)` when alphaM != 0 -- so with any
alphaM-proportional damping the element mass matrix is formed TWICE per element
per iteration, once for the M term and once inside C. That doubling is the
entire reason the headroom moves.

G2 headroom = the work a constant-M/C cache could in principle delete:
  dof.tangent + dof.residual        (nodal M/C into A and B -- transient-only)
  + brick.inertia under elem.tangent and elem.residual   (element M re-formed)

Gate: >= 10% of step wall under PARDISO. Run at 4 threads, the T0b operating
point, since the headroom SHARE is largest when the step is shortest.

Run:
  OPS_PYD=<worktree>\\dist\\bin MKL_NUM_THREADS=4 python3.12 -S t2_nodal_mass.py
"""
from __future__ import annotations

import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "4")
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
N = int(os.environ.get("OPS_T2_N", "15"))
NSTEPS = int(os.environ.get("OPS_T2_STEPS", "4"))
REPEATS = int(os.environ.get("OPS_T2_REPEATS", "2"))
SYSTEM = os.environ.get("OPS_T2_SYSTEM", "Pardiso")
THREADS = os.environ["MKL_NUM_THREADS"]

E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
RHO = 7.85e-9
L = 100.0
DT = 2.0e-4

# Rayleigh coefficients shared by arms C and D; which SLOT they go in is the
# whole difference between the two (see the module docstring).
RAY_AM, RAY_BK = 0.5, 1.0e-6


def build(n, mass_mode, system_args):
    nx = ny = nz = n

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    rho = RHO if mass_mode == "material" else 0.0
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso,
                   "-rho", rho)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)

    if mass_mode != "material":
        # Same TOTAL mass as arm A, lumped equally on every node, so the two
        # decks carry equivalent inertia and only its ROUTE differs.
        nnode = (nx + 1) * (ny + 1) * (nz + 1)
        m = RHO * (nx * ny * nz) * (L ** 3) / nnode
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    ops.mass(nid(i, j, k), m, m, m)

    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element("LadrunoBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1),
                            nid(i, j+1, k+1), 1)
                e += 1
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), 4.0e5, 0.0, -4.0e5)

    # ops.rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm)
    if mass_mode == "nodal_rayleigh":
        ops.rayleigh(RAY_AM, 0.0, RAY_BK, 0.0)      # betaK0: C is CONSTANT
    elif mass_mode == "nodal_rayleigh_kt":
        ops.rayleigh(RAY_AM, RAY_BK, 0.0, 0.0)      # betaKcurr: C is STATE-DEPENDENT

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-7, 100)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    return nid(nx // 2, ny // 2, nz), e - 1, (nnode - (nx+1)*(ny+1)) * 3


def run_once(mass_mode, profile_tag=None):
    tip, nele, ndof = build(N, mass_mode, [SYSTEM])
    if profile_tag:
        ops.profiler("start", "-deep", "-perStep")
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1, DT) != 0:
            raise RuntimeError(f"analyze failed at step {s+1} ({mass_mode})")
    wall = time.perf_counter() - t0
    if profile_tag:
        ops.profiler("stop")
        ops.profiler("report", str(HERE / f"t2_{profile_tag}_n{N}.h5"),
                     "-run", profile_tag)
        ops.profiler("reset")
    return wall, ops.nodeDisp(tip, 3), nele, ndof


ARMS = [
    ("material", "A: material mass (-rho), = the T0 deck"),
    ("nodal", "B: nodal mass (`mass` cmd), -rho 0"),
    ("nodal_rayleigh", "C: nodal + rayleigh betaK0 (C constant)"),
    ("nodal_rayleigh_kt", "D: nodal + rayleigh betaKcurr (C varies)"),
]


def main():
    print(f"ADR-77 T2 — G2 headroom on a nodal-mass deck | n={N} steps={NSTEPS} "
          f"system={SYSTEM} threads={THREADS}")
    out = {"n": N, "steps": NSTEPS, "system": SYSTEM, "threads": THREADS,
           "arms": {}}
    for key, label in ARMS:
        samples = []
        for _ in range(REPEATS):
            w, uz, nele, ndof = run_once(key)
            samples.append(w)
        med = statistics.median(samples)
        run_once(key, profile_tag=key)
        out["arms"][key] = {"label": label, "wall_s": med, "uz": uz,
                            "nele": nele, "ndof": ndof}
        print(f"  {label:<42} wall={med:7.3f}s  uz={uz: .9e}")
    (HERE / "t2_nodal_mass.json").write_text(json.dumps(out, indent=2))
    print("wrote t2_nodal_mass.json + per-arm .h5 profiles")


if __name__ == "__main__":
    main()
