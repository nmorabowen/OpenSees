"""Lane A (FIBER FRAME) Phase-0 measurement.
Nonlinear fiber-section 2D moment frame, forceBeamColumn + Concrete02 + Steel02.
Profile ONLY the pushover.
Run: pythoncore-3.12-64\\python.exe -S laneA_fiber_frame.py
"""
import sys, os, time

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
os.environ["PMI_RANK"] = "0"
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "laneA_fiber_frame.h5")
if os.path.exists(OUT):
    os.remove(OUT)

# ---------------------------------------------------------------- units: kip, in
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

# ---- geometry -------------------------------------------------------------
NSTORY = 30
NBAY = 8
H = 144.0   # story height in
L = 240.0   # bay width in
nX = NBAY + 1
nY = NSTORY + 1

def nid(ix, iy):
    return iy * nX + ix + 1

for iy in range(nY):
    for ix in range(nX):
        ops.node(nid(ix, iy), ix * L, iy * H)

# base fixity
for ix in range(nX):
    ops.fix(nid(ix, 0), 1, 1, 1)

# ---- materials ------------------------------------------------------------
# Concrete02: tag  fpc  epsc0  fpcu  epsU  lambda  ft  Ets
fc = -5.0
# confined core (higher strength/ductility)
ops.uniaxialMaterial("Concrete02", 1, -6.0, -0.004, -1.2, -0.020, 0.1, 0.6, 60.0)
# unconfined cover (small residual + longer softening tail = state-det friendly)
ops.uniaxialMaterial("Concrete02", 2, -5.0, -0.002, -0.5, -0.012, 0.1, 0.5, 50.0)
# Steel02: tag  Fy  E  b  R0 cR1 cR2
ops.uniaxialMaterial("Steel02", 3, 60.0, 29000.0, 0.01, 18.0, 0.925, 0.15)

# ---- fiber sections (rectangular RC) --------------------------------------
# section builder: width b, depth d, cover c, nbar top/bot, section tag
def rc_section(tag, b, d, cover, nbar_face):
    ops.section("Fiber", tag)
    yc = d / 2.0
    zc = b / 2.0
    yci = yc - cover
    zci = zc - cover
    # confined core patch (many fibers)
    ops.patch("rect", 1, 10, 8, -yci, -zci, yci, zci)
    # cover patches (unconfined) — 4 strips
    ops.patch("rect", 2, 10, 2, yci, -zc, yc, zc)     # top
    ops.patch("rect", 2, 10, 2, -yc, -zc, -yci, zc)   # bottom
    ops.patch("rect", 2, 2, 2, -yci, zci, yci, zc)    # right
    ops.patch("rect", 2, 2, 2, -yci, -zc, yci, -zci)  # left
    # steel layers top and bottom
    As = 1.0  # in^2 per bar (#9-ish)
    ops.layer("straight", 3, nbar_face, As, yci, -zci, yci, zci)   # top
    ops.layer("straight", 3, nbar_face, As, -yci, -zci, -yci, zci) # bottom

# column section 24x24, beam section 20x30
rc_section(1, 24.0, 24.0, 2.0, 4)   # column
rc_section(2, 20.0, 30.0, 2.0, 4)   # beam

# ---- integration + transforms ---------------------------------------------
ops.beamIntegration("Lobatto", 1, 1, 5)   # column
ops.beamIntegration("Lobatto", 2, 2, 5)   # beam
ops.geomTransf("PDelta", 1)   # columns
ops.geomTransf("Linear", 2)   # beams

# ---- elements -------------------------------------------------------------
etag = 1
col_eles = []
# columns
for iy in range(NSTORY):
    for ix in range(nX):
        ops.element("forceBeamColumn", etag, nid(ix, iy), nid(ix, iy + 1), 1, 1)
        col_eles.append(etag)
        etag += 1
# beams
for iy in range(1, nY):
    for ix in range(NBAY):
        ops.element("forceBeamColumn", etag, nid(ix, iy), nid(ix + 1, iy), 2, 2)
        etag += 1

nele = etag - 1
ndof = 3 * (nY * nX - nX)  # free dofs (base fixed)
print(f"nodes={nY*nX} elements={nele} (cols={len(col_eles)}) free_dof~={ndof}")

# ---- gravity load ---------------------------------------------------------
W = 15.0  # kip per beam-node tributary point load (downward)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for iy in range(1, nY):
    for ix in range(nX):
        ops.load(nid(ix, iy), 0.0, -W, 0.0)

ops.system("UmfPack")
ops.numberer("RCM")
ops.constraints("Transformation")
ops.test("NormDispIncr", 1e-8, 25)
ops.algorithm("Newton")
ops.integrator("LoadControl", 0.1)
ops.analysis("Static")

# --- gravity (NOT profiled) ---
t0 = time.perf_counter()
ok = ops.analyze(10)
print("gravity ok:", ok, f"({time.perf_counter()-t0:.2f}s)")
ops.loadConst("-time", 0.0)

# ---- pushover: inverted-triangle lateral pattern, DisplacementControl at roof
ops.pattern("Plain", 2, 1)
for iy in range(1, nY):
    Fh = float(iy)  # triangular
    for ix in range(nX):
        ops.load(nid(ix, iy), Fh, 0.0, 0.0)

roof_ctrl = nid(0, nY - 1)   # top-left node controls
Htot = NSTORY * H
nsteps = 150
target_drift = 0.012
dU = target_drift * Htot / nsteps  # per step
ops.integrator("DisplacementControl", roof_ctrl, 1, dU)

print(f"pushover: ctrl node {roof_ctrl}, dU={dU:.4f} in/step, target roof disp={target_drift*Htot:.2f} in")

# --- pushover (PROFILED) ---
ops.profiler("start", "-deep", "-perStep")
t0 = time.perf_counter()
done = 0
for i in range(nsteps):
    ok = ops.analyze(1)
    if ok != 0:
        # retry: substep + fall back to KrylovNewton
        ops.integrator("DisplacementControl", roof_ctrl, 1, dU / 8.0)
        ops.algorithm("KrylovNewton")
        sub_ok = 0
        for _ in range(8):
            if ops.analyze(1) != 0:
                sub_ok = -1
                break
        ops.integrator("DisplacementControl", roof_ctrl, 1, dU)
        ops.algorithm("Newton")
        if sub_ok != 0:
            print(f"  FAILED at step {i}, stopping")
            break
    done += 1
wall = time.perf_counter() - t0
ops.profiler("stop")
ops.profiler("report", OUT)

roof_disp = ops.nodeDisp(roof_ctrl, 1)
print(f"pushover done: {done} steps, roof disp={roof_disp:.3f} in "
      f"(drift={roof_disp/Htot*100:.2f}%), wall={wall:.2f}s")
print("report written:", os.path.exists(OUT), OUT)
