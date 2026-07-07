"""Lane-D bench, parameterized (ADR-40 sanity re-run, 2026-07-07).

Same model as Ladruno_files/testbed/perf/phase0/laneD_model.py (live-wave fix
included). Env knobs:
  OPS_DIST   - dir containing opensees.pyd (REQUIRED)
  CDL_FLAGS  - extra integrator flags, space-separated (e.g. "-commitSolveState")
  ALGO_FLAGS - extra algorithm Linear flags (e.g. "-factorOnce")
  PROF_OUT   - h5 report path; empty/unset = profiler OFF (clean wall)
Prints one machine-parseable RESULT line.
"""
import sys, os, time, math

DIST = os.environ["OPS_DIST"]
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
os.environ["PMI_RANK"] = "0"
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import opensees as ops
_got = os.path.normcase(os.path.normpath(ops.__file__))
_want = os.path.normcase(os.path.normpath(DIST))
assert _got.startswith(_want), f"WRONG PYD: {ops.__file__}"

CDL_FLAGS = os.environ.get("CDL_FLAGS", "").split()
ALGO_FLAGS = os.environ.get("ALGO_FLAGS", "").split()
PROF_OUT = os.environ.get("PROF_OUT", "")

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
NX, NY, NZ = 10, 10, 10
H = 1.0

def nid(i, j, k):
    return 1 + i + (NX + 1) * (j + (NY + 1) * k)

for k in range(NZ + 1):
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.node(nid(i, j, k), i * H, j * H, k * H)
for j in range(NY + 1):
    for i in range(NX + 1):
        ops.fix(nid(i, j, 0), 1, 1, 1)

E, nu, RHO = 200.0e9, 0.3, 7850.0
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 250.0e6, 0.0, 0.0, 2.0e9,
               "-rho", RHO)

eid = 0
for k in range(NZ):
    for j in range(NY):
        for i in range(NX):
            eid += 1
            ops.element("LadrunoBrick", eid,
                        nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                        nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1),
                        1, "-lumped")

ops.timeSeries("Path", 1, "-dt", 1.0e-3, "-values", 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
ops.pattern("Plain", 1, 1)
for j in range(NY + 1):
    for i in range(NX + 1):
        ops.load(nid(i, j, NZ), 0.0, 0.0, -1.0e7)

ops.constraints("Transformation")
ops.numberer("Plain")
ops.system("Diagonal")
ops.test("NormDispIncr", 1e-12, 1)
ops.algorithm("Linear", *ALGO_FLAGS)
ops.integrator("CentralDifferenceLadruno", "-cfl", *CDL_FLAGS)
ops.analysis("Transient")

c_dil = math.sqrt((K + 4.0 * G / 3.0) / RHO)
dt_est = 0.5 * H / c_dil
ops.analyze(1, dt_est)
dtcr = ops.criticalTimeStep()
dt = 0.8 * dtcr if (dtcr and dtcr > 0) else dt_est

NSTEPS = 2500
if PROF_OUT:
    if os.path.exists(PROF_OUT):
        os.remove(PROF_OUT)
    ops.profiler("start", "-deep", "-perStep")
t0 = time.perf_counter()
nfail = 0
for _ in range(NSTEPS):
    if ops.analyze(1, dt) != 0:
        nfail += 1
        break
t1 = time.perf_counter()
if PROF_OUT:
    ops.profiler("stop")
    ops.profiler("report", PROF_OUT)

top = nid(NX // 2, NY // 2, NZ)
uz = ops.nodeDisp(top, 3)
print(f"RESULT wall={t1-t0:.3f} nfail={nfail} dt={dt:.6e} uz={uz:.17e} "
      f"dist={DIST} cdl='{' '.join(CDL_FLAGS)}' algo='{' '.join(ALGO_FLAGS)}' "
      f"prof={'on' if PROF_OUT else 'off'}")
