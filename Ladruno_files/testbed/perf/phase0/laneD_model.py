"""ADR-40 Phase-0 lane D (EXPLICIT): CentralDifferenceLadruno on a LadrunoBrick
J2 hex block, impulsive top load. Profiles ONLY the analyze block.

Run: pythoncore-3.12-64\python.exe -S laneD_model.py
"""
import sys, os, time, math

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
os.environ["PMI_RANK"] = "0"
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "laneD_prof.h5")
if os.path.exists(OUT):
    os.remove(OUT)

# ---------------- model ----------------
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)

# --- geometry: Nx x Ny x Nz hex block ---
NX, NY, NZ = 10, 10, 10          # 1000 elements, 11^3 = 1331 nodes
H = 1.0                          # element edge length (m)

def nid(i, j, k):                # node id from grid index
    return 1 + i + (NX + 1) * (j + (NY + 1) * k)

for k in range(NZ + 1):
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.node(nid(i, j, k), i * H, j * H, k * H)

# fixed base (k = 0 plane), all 3 dof
for j in range(NY + 1):
    for i in range(NX + 1):
        ops.fix(nid(i, j, 0), 1, 1, 1)

# --- material: LadrunoJ2, linear isotropic hardening (Voce with Qinf=b=0, Hiso>0) ---
# steel-like: E=200 GPa, nu=0.3 -> K = E/(3(1-2nu)), G = E/(2(1+nu)); rho=7850 kg/m^3
E   = 200.0e9
nu  = 0.3
K   = E / (3.0 * (1.0 - 2.0 * nu))
G   = E / (2.0 * (1.0 + nu))
RHO = 7850.0
s0    = 250.0e6      # yield stress
Hiso  = 2.0e9        # linear isotropic hardening modulus
ops.nDMaterial("LadrunoJ2", 1, K, G,
               "-iso", "voce", s0, 0.0, 0.0, Hiso,   # Qinf=0, b=0 -> pure linear Hiso
               "-rho", RHO)

# --- elements: LadrunoBrick std, lumped mass (rho from material) ---
# node order for a hex: 8 corners of the cell (i..i+1, j..j+1, k..k+1)
eid = 0
for k in range(NZ):
    for j in range(NY):
        for i in range(NX):
            eid += 1
            n1 = nid(i,   j,   k)
            n2 = nid(i+1, j,   k)
            n3 = nid(i+1, j+1, k)
            n4 = nid(i,   j+1, k)
            n5 = nid(i,   j,   k+1)
            n6 = nid(i+1, j,   k+1)
            n7 = nid(i+1, j+1, k+1)
            n8 = nid(i,   j+1, k+1)
            ops.element("LadrunoBrick", eid, n1, n2, n3, n4, n5, n6, n7, n8, 1,
                        "-lumped")

nElem = eid
nNode = (NX + 1) * (NY + 1) * (NZ + 1)
print(f"model built: {nElem} elements, {nNode} nodes")

# --- impulsive top load: short triangular pulse on the top face (k = NZ) ---
# push in -z on every top node
ops.timeSeries("Path", 1, "-dt", 1.0e-5, "-values",
               0.0, 1.0, 0.0, 0.0, 0.0, 0.0)   # rise then fall over ~2e-5 s, then zero
ops.pattern("Plain", 1, 1)
Pnode = -1.0e7   # N per top node, -z
for j in range(NY + 1):
    for i in range(NX + 1):
        ops.load(nid(i, j, NZ), 0.0, 0.0, Pnode)

# ---------------- analysis objects (explicit recipe from CDL tests) ----------------
ops.constraints("Transformation")
ops.numberer("Plain")
ops.system("Diagonal")
ops.test("NormDispIncr", 1e-12, 1)
ops.algorithm("Linear")
ops.integrator("CentralDifferenceLadruno", "-cfl")
ops.analysis("Transient")

# --- critical dt: estimate c = sqrt((K+4G/3)/rho) (dilatational), dt = 0.5 h/c ---
c_dil = math.sqrt((K + 4.0 * G / 3.0) / RHO)
dt_est = 0.5 * H / c_dil
print(f"c_dil = {c_dil:.1f} m/s, h/c = {H/c_dil:.3e} s, dt_est(0.5 h/c) = {dt_est:.3e} s")

# prime one step to trigger dt_cr compute, then query criticalTimeStep()
ops.analyze(1, dt_est)
try:
    dtcr = ops.criticalTimeStep()
    print(f"criticalTimeStep() = {dtcr:.3e} s")
    dt = 0.8 * dtcr if (dtcr and dtcr > 0) else dt_est
except Exception as e:
    print("criticalTimeStep() unavailable:", e)
    dt = dt_est
print(f"using dt = {dt:.3e} s")

NSTEPS = 2500

# ---------------- profile ONLY the analyze block ----------------
ops.profiler("start", "-deep", "-perStep")
t0 = time.perf_counter()
nfail = 0
for _ in range(NSTEPS):
    if ops.analyze(1, dt) != 0:
        nfail += 1
        break
t1 = time.perf_counter()
ops.profiler("stop")
ops.profiler("report", OUT)

top = nid(NX // 2, NY // 2, NZ)
print(f"analyze block: {NSTEPS} steps, {nfail} failures, wall {t1 - t0:.2f} s")
print(f"top-center disp_z = {ops.nodeDisp(top, 3):.6e} m")
print("report written:", os.path.exists(OUT), OUT)
