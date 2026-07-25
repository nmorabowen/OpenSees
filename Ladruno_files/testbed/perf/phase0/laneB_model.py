"""
ADR-40 Phase-0 Lane B: implicit-static 3D solid on the fork's own
LadrunoBrick (classTag 33002) + LadrunoJ2 (nDMaterial).
Mesh: 12x12x12 regular hex (1728 elems, 2197 nodes, ~6591 DOF).
Profiles ONLY the analyze() block.
"""
import sys, os, time

# Ladruno ADR-75b (L3-0): DIST/OUT/system are env-overridable so this exact model
# can be re-run against a different build or solver without editing (or forking) it.
# Defaults reproduce the original ADR-40b run verbatim.
DIST = os.environ.get("OPS_PYD") or r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
os.environ["PMI_RANK"] = "0"
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("OPS_PROF_OUT") or os.path.join(HERE, "laneB_prof.h5")
if os.path.exists(OUT):
    os.remove(OUT)

# ---- units: N, mm, MPa ----
E  = 200000.0        # MPa (steel)
nu = 0.3
K  = E / (3.0 * (1.0 - 2.0 * nu))    # bulk
G  = E / (2.0 * (1.0 + nu))          # shear
s0   = 250.0         # yield stress MPa
Qinf = 0.0           # no Voce saturation
b    = 0.0           # no Voce rate
Hiso = 2000.0        # linear isotropic hardening modulus MPa

nx = ny = nz = 15
L  = 100.0           # element edge (mm) -> 1500 mm cube

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)

ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso)

# ---- nodes ----
def nid(i, j, k):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)

for k in range(nz + 1):
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(nid(i, j, k), i * L, j * L, k * L)

# ---- fixed base (k=0) ----
for j in range(ny + 1):
    for i in range(nx + 1):
        ops.fix(nid(i, j, 0), 1, 1, 1)

# ---- elements: std formulation (default) ----
eTag = 1
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            n1 = nid(i,     j,     k)
            n2 = nid(i + 1, j,     k)
            n3 = nid(i + 1, j + 1, k)
            n4 = nid(i,     j + 1, k)
            n5 = nid(i,     j,     k + 1)
            n6 = nid(i + 1, j,     k + 1)
            n7 = nid(i + 1, j + 1, k + 1)
            n8 = nid(i,     j + 1, k + 1)
            ops.element("LadrunoBrick", eTag, n1, n2, n3, n4, n5, n6, n7, n8, 1)
            eTag += 1
nElem = eTag - 1
nNode = (nx + 1) * (ny + 1) * (nz + 1)
nFixed = (nx + 1) * (ny + 1)
nDOF = (nNode - nFixed) * 3
print(f"nElem={nElem} nNode={nNode} nDOF~={nDOF}")

# ---- loads on top face (k=nz): lateral (x) + vertical down (z) ----
top_nodes = [nid(i, j, nz) for j in range(ny + 1) for i in range(nx + 1)]
Fx = 4.0e5   # N per top node (lateral, drives shear yield)
Fz = -4.0e5  # N per top node (compression)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for nd in top_nodes:
    ops.load(nd, Fx, 0.0, Fz)

# ---- ADR-40 rank-2/8 target path ----
ops.constraints("Plain")
ops.numberer("RCM")
# Ladruno ADR-75b (L3-0): OPS_SYSTEM lets the same model run on a different solver so
# the element/loop FRACTIONS can be re-measured after the Lane-1 solver win. Values are
# passed as ints where numeric -- a STRING option value fails OPS_GetIntInput, the factory
# returns 0, and OpenSees silently falls back to ProfileSPDLinSOE (banked P1d trap).
_sys = (os.environ.get("OPS_SYSTEM") or "UmfPack").split()
ops.system(*[int(a) if a.lstrip("-").isdigit() else a for a in _sys])
print("system:", _sys)
ops.test("NormDispIncr", 1e-7, 25)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0 / 15.0)
ops.analysis("Static")

# ---- timed block: profile ONLY analyze ----
# Ladruno ADR-75b (L3-0): OPS_PROF_FLAGS lets the SAME run be repeated without -deep so the
# deep-scope instrumentation tax can be MEASURED rather than argued about (ADR-40b banked
# ~0.5 us per scope instance). Default is the original "-deep -perStep".
_pf = (os.environ.get("OPS_PROF_FLAGS", "-deep -perStep")).split()
ops.profiler("start", *_pf)
print("prof flags:", _pf)
t0 = time.perf_counter()
nsteps = 15
done = 0
for s in range(nsteps):
    ok = ops.analyze(1)
    if ok != 0:
        print(f"step {s+1} FAILED (ok={ok})")
        break
    done += 1
wall = time.perf_counter() - t0
ops.profiler("stop")
ops.profiler("report", OUT)
print(f"steps done={done}/{nsteps}  wall={wall:.2f}s  report={os.path.exists(OUT)}  {OUT}")
