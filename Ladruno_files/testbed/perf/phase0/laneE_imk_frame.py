"""Lane E (CONCENTRATED-PLASTICITY / IMK) Phase-0 measurement.
2D moment frame, every beam+column a LadrunoIMKBeam2d (classTag 33004) with a
modified-IMK Bilin hinge at each end (elastic interior in series). CYCLIC
DisplacementControl pushover so the hinges load/unload/reload -> that is what
forces the 2x2 internal hinge Newton to iterate. Profile ONLY the cyclic block.

Run: pythoncore-3.12-64\\python.exe -S laneE_imk_frame.py
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
OUT = os.path.join(HERE, "laneE_imk_frame.h5")
if os.path.exists(OUT):
    os.remove(OUT)

# ---------------------------------------------------------------- units: kip, in
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

# ---- geometry -------------------------------------------------------------
NSTORY = 20
NBAY = 10
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

# ---- elastic member properties (kip, in) ----------------------------------
E = 29000.0
# columns: W14x-ish
Acol, Icol = 30.0, 1500.0
# beams: W24x-ish
Abeam, Ibeam = 20.0, 2700.0

# ---- IMK (Bilin, modified-Ibarra-Medina-Krawinkler) hinge laws ------------
# Rotational moment(kip-in)-rotation(rad) backbone with cyclic deterioration.
# Params mirror the fork's own test (test_ladrunoIMKBeam2d_element.py test 6):
#   Ke, AsPos, AsNeg, My+, My-, LamdaS,LamdaC,LamdaA,LamdaK,
#   cS,cC,cA,cK, thetaP+,thetaP-, thetaPc+,thetaPc-, ResP+,ResP-, thetaU+,thetaU-, D+,D-
def bilin_hinge(tag, My, Ke):
    # Modified-IMK backbone kept in the STRAIN-HARDENING regime for the drift
    # range of interest (thetaP large -> cap onset well past 1.5% drift) so the
    # static cyclic pushover converges; cyclic deterioration (finite Lambda)
    # still active so the 2x2 hinge Newton iterates on every load/unload/reload.
    ops.uniaxialMaterial("Bilin", tag,
                         Ke, 0.05, 0.05, My, -My,
                         1.0e6, 1.0e6, 1.0e6, 1.0e6,   # Lambda huge -> cyclic deter. OFF
                         1.0, 1.0, 1.0, 1.0,   # c exponents
                         0.50, 0.50,           # thetaP (pre-cap plastic rot, large)
                         0.50, 0.50,           # thetaPc (post-cap)
                         0.5, 0.5,             # residual strength ratio
                         1.00, 1.00,           # thetaU (ultimate rot, large)
                         1.0, 1.0)             # rate of cyclic deterioration D

# hinge stiffness: n*6EI/L (stiff spring, standard IMK modeling); yield moments
# sized so hinges yield progressively over the 0.3->1.5% drift protocol while
# staying in the strain-hardening branch (stable static cyclic response).
Kcol = 10.0 * 6.0 * E * Icol / H
Kbeam = 10.0 * 6.0 * E * Ibeam / L
My_col = 14.0e3     # kip-in (strong-column / weak-beam -> spread plasticity, no soft story)
My_beam = 5.0e3     # kip-in
bilin_hinge(101, My_col, Kcol)    # column hinge law
bilin_hinge(102, My_beam, Kbeam)  # beam hinge law

# ---- transforms -----------------------------------------------------------
ops.geomTransf("Linear", 1)   # columns (linear: isolate hinge-Newton cost, no PDelta softening)
ops.geomTransf("Linear", 2)   # beams

# ---- elements: LadrunoIMKBeam, hinges BOTH ends ---------------------------
# element('LadrunoIMKBeam', tag, iNode, jNode, A, E, Iz, transfTag, '-matZ', imkTag)
etag = 1
ncol = 0
for iy in range(NSTORY):
    for ix in range(nX):
        ops.element("LadrunoIMKBeam", etag, nid(ix, iy), nid(ix, iy + 1),
                    Acol, E, Icol, 1, "-matZ", 101)
        etag += 1
        ncol += 1
nbeam = 0
for iy in range(1, nY):
    for ix in range(NBAY):
        ops.element("LadrunoIMKBeam", etag, nid(ix, iy), nid(ix + 1, iy),
                    Abeam, E, Ibeam, 2, "-matZ", 102)
        etag += 1
        nbeam += 1

nele = etag - 1
ndof = 3 * (nY * nX - nX)
print(f"nodes={nY*nX} elements={nele} (cols={ncol}, beams={nbeam}) free_dof~={ndof}")

# ---- gravity load (NOT profiled) ------------------------------------------
W = 15.0  # kip per node
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

t0 = time.perf_counter()
ok = ops.analyze(10)
print("gravity ok:", ok, f"({time.perf_counter()-t0:.2f}s)")
ops.loadConst("-time", 0.0)

# ---- cyclic lateral pushover (PROFILED) -----------------------------------
# inverted-triangle lateral force pattern; DisplacementControl at the roof,
# reversed cycles growing to ~1.5% roof drift so hinges load/unload/reload.
ops.pattern("Plain", 2, 1)
for iy in range(1, nY):
    Fh = float(iy)
    for ix in range(nX):
        ops.load(nid(ix, iy), Fh, 0.0, 0.0)

roof_ctrl = nid(0, nY - 1)
roof_ctrl_ele = 1  # first base column — a representative heavily-yielding element
Htot = NSTORY * H

# cyclic protocol: growing reversed cycles (drift ratios), 1 cycle each
drift_ratios = [0.003, 0.006, 0.010, 0.015]
targets = []
for dr in drift_ratios:
    amp = dr * Htot
    targets += [amp, -amp, amp, -amp]  # 2 reversed cycles per amplitude
targets += [0.0]

dU_base = 0.015 * Htot / 100.0  # ~0.43 in nominal step (fine incr for stable reversals)

iter_hist = {}   # global-Newton-iterations -> count of steps
peak_hinge = [0.0]

def _track():
    it = ops.testIter()
    iter_hist[it] = iter_hist.get(it, 0) + 1
    hr = ops.eleResponse(roof_ctrl_ele, "hingeRotation")  # a representative element
    if hr:
        peak_hinge[0] = max(peak_hinge[0], max(abs(x) for x in hr))

def _adaptive(step, depth=0):
    """Advance ONE nominal increment 'step' with recursive bisection on failure.
    Returns number of profiled analyze() calls done, or -1 if unrecoverable."""
    ops.integrator("DisplacementControl", roof_ctrl, 1, step)
    if ops.analyze(1) == 0:
        _track()
        return 1
    if depth >= 5:
        return -1
    # bisect: two half-steps
    n1 = _adaptive(step / 2.0, depth + 1)
    if n1 < 0:
        return -1
    n2 = _adaptive(step / 2.0, depth + 1)
    if n2 < 0:
        return -1
    return n1 + n2

def run_to(target, cur):
    """DisplacementControl from cur to target with adaptive bisection. Returns
    (new_cur, n_analyze_calls, ok)."""
    step = dU_base if target >= cur else -dU_base
    nstep = max(1, int(round(abs(target - cur) / dU_base)))
    done = 0
    for _ in range(nstep):
        n = _adaptive(step)
        if n < 0:
            return cur, done, -1
        done += n
        cur += step
    return cur, done, 0

ops.profiler("start", "-deep", "-perStep")
t0 = time.perf_counter()
cur = 0.0
total = 0
for tgt in targets:
    cur, done, ok = run_to(tgt, cur)
    total += done
    if ok != 0:
        print(f"  FAILED heading to {tgt:.2f} in (roof), stopping")
        break
wall = time.perf_counter() - t0
ops.profiler("stop")
ops.profiler("report", OUT)

# yield diagnostics: max hinge rotation across all elements (theta_y ~ My/Ke)
max_hinge = 0.0
for e in range(1, nele + 1):
    hr = ops.eleResponse(e, "hingeRotation")
    if hr:
        max_hinge = max(max_hinge, max(abs(x) for x in hr))
theta_y_col = My_col / Kcol
theta_y_beam = My_beam / Kbeam

roof_disp = ops.nodeDisp(roof_ctrl, 1)
print(f"cyclic done: {total} steps, roof disp={roof_disp:.3f} in "
      f"(drift={roof_disp/Htot*100:.3f}%), max target drift={max(drift_ratios)*100:.2f}%, "
      f"wall={wall:.2f}s")
print(f"global Newton-iters-per-step histogram (iters:steps): "
      f"{dict(sorted(iter_hist.items()))}")
print(f"peak hinge rotation (representative top-beam) = {peak_hinge[0]:.5f} rad")
print(f"max hinge rotation (all eles, final state) ={max_hinge:.5f} rad ; "
      f"theta_y col={theta_y_col:.5f} beam={theta_y_beam:.5f}"
      f"  -> peak ductility col~{peak_hinge[0]/theta_y_col:.1f}x beam~{peak_hinge[0]/theta_y_beam:.1f}x")
print("report written:", os.path.exists(OUT), OUT)
