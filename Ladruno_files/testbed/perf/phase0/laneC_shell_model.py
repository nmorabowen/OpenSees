"""Lane C (ADR-68 T2): ASDShellQ4 cantilever wall driving LadrunoJ2 through the
PLATE-FIBER (condensed sigma_22=0) path via section LayeredShell.

Usage:
  pythoncore-3.12-64\\python.exe -S laneC_shell_model.py plastic  <out.h5>
  pythoncore-3.12-64\\python.exe -S laneC_shell_model.py elastic  <out.h5>

'plastic' -> LadrunoJ2 (condensed Newton on eps_22, T2 target).
'elastic' -> ElasticIsotropic nDMaterial (same PlateFiber wrapping, no return map)
             => the analyze-time diff (plastic - elastic) IS the J2-condensation cost.
"""
import sys, os, time

MODE = sys.argv[1] if len(sys.argv) > 1 else "plastic"
OUT  = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), f"laneC_{MODE}.h5")

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
os.environ["PMI_RANK"] = "0"
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__, "| MODE:", MODE)

if os.path.exists(OUT):
    os.remove(OUT)

# ---------------------------------------------------------------- geometry
# Cantilever wall in the X-Z plane (vertical), fixed at Z=0, meshed NxN.
# ASDShellQ4 needs ndm=3, ndf=6.
NX, NZ = 40, 40            # 40x40 = 1600 elements, 41x41 = 1681 nodes
LX, LZ = 3.0, 3.0         # metres  (units: N, m, Pa)
THICK  = 0.20             # 200 mm wall
NLAYER = 8                # >=5 fibers through thickness (condensed path each)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)

nid = {}
tagn = 1
for j in range(NZ + 1):
    for i in range(NX + 1):
        x = LX * i / NX
        z = LZ * j / NZ
        ops.node(tagn, x, 0.0, z)
        nid[(i, j)] = tagn
        tagn += 1

# fix base row (Z=0) fully
for i in range(NX + 1):
    ops.fix(nid[(i, 0)], 1, 1, 1, 1, 1, 1)

# ---------------------------------------------------------------- material/section
E  = 200.0e9      # steel-like, Pa
nu = 0.30
rho = 7850.0
K  = E / (3.0 * (1.0 - 2.0 * nu))
G  = E / (2.0 * (1.0 + nu))

if MODE == "plastic":
    # LadrunoJ2: sig0=250 MPa, LINEAR isotropic hardening (Hiso), no Voce sat, no kin.
    #   nDMaterial LadrunoJ2 tag K G -iso voce s0 Qinf b Hiso  -rho rho
    fy   = 250.0e6
    Hiso = 2.0e9      # linear hardening modulus (~1% of E)
    ops.nDMaterial("LadrunoJ2", 1, K, G,
                   "-iso", "voce", fy, 0.0, 0.0, Hiso,
                   "-rho", rho)
elif MODE == "elastic":
    # matched elastic PlateFiber wrapping: same E,nu -> identical shell tangent
    # assembly, but the material return is a closed-form elastic evaluation
    # (still runs the condensed eps_22 Newton? NO -> ElasticIsotropic PlateFiber
    #  is a direct plane-stress form). This isolates constitutive cost.
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
else:
    raise SystemExit("MODE must be plastic|elastic")

# section LayeredShell tag nLayers matTag totalThickness  (2-arg distribute form)
ops.section("LayeredShell", 1, NLAYER, 1, THICK)

# ---------------------------------------------------------------- elements
etag = 1
for j in range(NZ):
    for i in range(NX):
        n1 = nid[(i, j)]
        n2 = nid[(i + 1, j)]
        n3 = nid[(i + 1, j + 1)]
        n4 = nid[(i, j + 1)]
        ops.element("ASDShellQ4", etag, n1, n2, n3, n4, 1)
        etag += 1
nelem = etag - 1
nnode = tagn - 1
print(f"elements: {nelem} ASDShellQ4 | nodes: {nnode} | layers: {NLAYER}")

# ---------------------------------------------------------------- load
# In-plane shear (X at top edge) + out-of-plane push (Y at top edge) combined,
# so a large fraction of Gauss points across the thickness yield.
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
topX_total   = 1.28e7   # N in-plane, spread on top edge (pushes broad yield)
topY_total   = 3.0e6    # N out-of-plane
ntop = NX + 1
for i in range(NX + 1):
    nd = nid[(i, NZ)]
    ops.load(nd, topX_total / ntop, topY_total / ntop, 0.0, 0.0, 0.0, 0.0)

# ---------------------------------------------------------------- analysis
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormDispIncr", 1e-7, 25)
ops.algorithm("Newton")
NSTEPS = 15
ops.integrator("LoadControl", 1.0 / NSTEPS)
ops.analysis("Static")

ndof = ops.systemSize()
print("system size (DOF):", ndof)

# ------- sanity: confirm the condensed path actually fires (plastic mode) -------
if MODE == "plastic":
    # one probe step, then read a layer stress to prove setTrialStrain ran
    pass

# ---------------------------------------------------------------- PROFILE
ops.profiler("start", "-deep", "-perStep")
t0 = time.perf_counter()
done = 0
for s in range(NSTEPS):
    ok = ops.analyze(1)
    if ok != 0:
        print(f"  step {s+1} FAILED (ok={ok}) -> stopping")
        break
    done += 1
wall = time.perf_counter() - t0
ops.profiler("stop")
ops.profiler("report", OUT)

print(f"analyze block: {done}/{NSTEPS} steps in {wall:.2f} s  | report: {os.path.exists(OUT)} {OUT}")

# tip displacement to confirm real deformation/plasticity
tipn = nid[(NX, NZ)]
ux = ops.nodeDisp(tipn, 1)
uy = ops.nodeDisp(tipn, 2)
print(f"tip disp: ux={ux:.4e} m  uy={uy:.4e} m")
