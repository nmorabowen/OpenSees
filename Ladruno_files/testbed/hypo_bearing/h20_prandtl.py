"""The quadratic-hex discriminator for the TIMs limit-load ceiling.

WHAT THIS MEASURES.  The TIMs element-technology campaign is measuring the
Prandtl-Reissner collapse load -- the one bearing-capacity problem whose limit
analysis has coincident upper and lower bounds, q_u = q0 * N_q on WEIGHTLESS
soil -- with different element technologies, all on the same problem:

    LadrunoBrick -formulation bbar   (LINEAR hex, B-bar)      ~0.99   plateau
    TenNodeTetrahedron               (quadratic Lagrange tet)  0.2828
    BezierTet10 std                  (quadratic Bernstein tet) 0.3075
    BezierTet10 -bbar                                          ~0.49-0.55

So B-bar buys ~0.31 -> ~0.50 on a quadratic tet, but the linear hex goes
0.40 -> 0.99.  Three explanations survive that data set:

    (a) SIMPLEX GEOMETRY -- one volumetric constraint per element is simply
        weaker on tets (constraint counting is adverse), so any tet tops out;
    (b) QUADRATIC ORDER  -- something about quadratic elements under
        near-incompressible perfectly-plastic flow;
    (c) the BEZIER/BERNSTEIN BASIS specifically.

(c) is already nearly excluded: at equal order and equal geometry the Bernstein
tet (0.3075) and the Lagrange tet (0.2828) agree to within 8 %, which is far
smaller than the gap being explained.

THIS SCRIPT RUNS THE MISSING CELL: `LadrunoBrick20`, the fork's 20-node
serendipity hex -- QUADRATIC and a HEX and Lagrange -- on the identical problem.

    quadratic hex reaches ~1.0     ->  the ceiling is (a) SIMPLEX GEOMETRY
    quadratic hex tops out too     ->  the ceiling is (b) QUADRATIC ORDER

NOTE ON THE VOLUMETRIC TREATMENT.  LadrunoBrick20 has no B-bar option, and it
does not need one: its two formulations are `std` (full 27-pt Gauss) and `uri`
(uniform reduced 2x2x2), and H20 + 2x2x2 IS the classical near-incompressible
workhorse -- the C3D20R of Abaqus.  Selective/reduced integration on a
quadratic hex plays the role B-bar plays on the linear hex, so `uri` is the
leg that is technology-matched to the 0.99 baseline and `std` is the locked
control that shows how much of the answer the volumetric treatment is worth.
That makes this a 2x2 design (linear/quadratic x locked/relieved) rather than
a single point.

EVERY PROBLEM PARAMETER IS THE REFERENCE DRIVER'S, UNCHANGED: phi_txc = 20 deg,
q0 = 10 kPa, gamma = 0 (weightless -- this is what makes the oracle exact),
nu = 0.45 (NOT lowered: at a low nu this cone starts at 95 % of yield and the
leg measures nothing but its own initial condition), the Chen & Han PLANE-STRAIN
cone match, non-associated (rho_bar = 0) as the gate leg with associated as the
falsification control, displacement-controlled push to s/B, the D16 adaptive
ladder and its three guards, the 2 %-drop truncation, and the plateau test.
The mesh is `dp_strip.py`'s structured graded strip grid, so it is the SAME
grid the 0.99 baseline was measured on -- only the element changes.

THE LOAD IS THE TRAP.  A uniform pressure on a quadratic serendipity face has
NEGATIVE corner weights (-A/12 at the four corners, +A/3 at the four
mid-edges).  Lumping it by tributary area -- which is right for an H8 face and
is what every mesh generator in this testbed hands you -- is WRONG here, and it
is exactly the defect that made a Bezier deck diverge at first yield.  This
script therefore (i) computes the face load by 3x3 Gauss integration of the Q8
shape functions over the real face geometry, (ii) asserts it against the
analytic rectangle weights, (iii) asserts the base reaction reproduces the
applied total to 1e-6, and (iv) offers a `patch` leg that puts an
ElasticIsotropic material on the same mesh and checks that EVERY Gauss point
sees the exact 1-D state sigma_zz = -q0, sigma_xx = sigma_yy = -K0 q0.  Check
(iii) only tests the SUM and a lumped load passes it; check (iv) is the one
that tests the distribution.

Run:  py -3.12 h20_prandtl.py [leg ...] [--h0 H] [--form uri|std] [--budget N]
      legs: patch | modes | nonassoc | assoc   (default: patch modes nonassoc)
"""
import argparse
import csv
import math
import os
import sys
import time

import numpy as np

_WT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
_BIN = os.path.join(_WT, "dist", "bin")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
if hasattr(os, "add_dll_directory") and os.path.isdir(_BIN):
    os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)

# H20_NO_ENGINE=1 imports the mesh / load-vector half of this module WITHOUT
# the solver, so the geometry and the consistent-load algebra can be unit-tested
# under any interpreter (dp_collapse.py's ADR79_NO_ENGINE idiom).
if os.environ.get("H20_NO_ENGINE") == "1":
    ops = None
else:
    import opensees as ops  # noqa: E402

    # an INSTALLED Ladruno hijacks `import opensees` through a site .pth; that
    # build is not the one under test.  Refuse rather than silently measure it.
    if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != \
            os.path.normcase(_BIN):
        raise SystemExit(f"WRONG ENGINE: imported {ops.__file__}, "
                         f"expected {_BIN}")

HERE = os.path.dirname(os.path.abspath(__file__))

# --- the problem (identical to r3_prandtl_tet10.py / dp_strip.py nq20) -------
PHI_TXC = 20.0            # the MILD cone: mesh-resolvable, and where the
                          # Prandtl/Vesic coefficients were fitted
Q0 = 10.0                 # kPa surcharge
GAMMA = 0.0               # WEIGHTLESS -- this is what makes the oracle exact
K_EL = 1.5e5              # kPa
NU = 0.45                 # free parameter; see the docstring
SY = 0.2                  # kPa apex regularizer, NOT a soil property
B_FOOT = 2.0
THICK = 0.5               # out-of-plane slab thickness, m (plane strain)
XLIM, ZBOT = 30.0, -20.0  # 14.5 B of clearance, 10 B of depth -- dp_strip.py's
R_GRADE = 1.35

# --- ADR 63 D16: the three adaptive-stepping guards -------------------------
DS_BASE, DS_MIN, DS_MAX = 2.0e-4, 2.0e-6, 1.0e-3
GROW_AFTER = 6            # recovery latch: clean steps before the step grows
STALL_WINDOW = 150        # steps
STALL_ADVANCE = 1.0e-4    # fraction of smax that must be covered in that window
DROP_TRUNCATE = 0.02      # truncate the curve where q falls >2 % below its max


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# --- cone / oracle ----------------------------------------------------------
def alpha_from_phi_txc(phi_deg):
    s = math.sin(math.radians(phi_deg))
    return 2.0 * s / (math.sqrt(3.0) * (3.0 - s))


def mc_from_cone(alpha):
    """Mohr-Coulomb angles matching a Drucker-Prager cone of slope alpha
    (Chen & Han).  `ps` is the PLANE-STRAIN match, which is the one this
    problem needs."""
    d = 1.0 - 12.0 * alpha ** 2
    s = 3.0 * math.sqrt(3.0) * alpha / (2.0 + math.sqrt(3.0) * alpha)
    return dict(txc=math.degrees(math.asin(s)),
                ps=math.degrees(math.atan(math.sqrt(9.0 * alpha ** 2 / d))))


def n_q(phi_deg):
    ph = math.radians(phi_deg)
    return math.exp(math.pi * math.tan(ph)) * math.tan(math.pi / 4 + ph / 2) ** 2


def davis_phi(phi_deg, psi_deg=0.0):
    """Davis (1968) reduced angle for NON-associated flow; at psi = 0 simply
    tan(phi*) = sin(phi).  Reported beside the oracle, never as the oracle."""
    p, s = math.radians(phi_deg), math.radians(psi_deg)
    return math.degrees(math.atan(math.cos(s) * math.sin(p)
                                  / (1.0 - math.sin(p) * math.sin(s))))


# --- the H20 strip mesh -----------------------------------------------------
def graded(lo, hi, n, h0, fine_at_lo):
    """n intervals over [lo,hi] whose element at the FINE end is h0 and which
    grow geometrically outward -- dp_strip.py's `graded`, inlined verbatim so
    the grid is bit-for-bit the one the LadrunoBrick baseline was measured on."""
    if abs(h0 * n - (hi - lo)) < 1e-12:
        return np.linspace(lo, hi, n + 1)
    a, b = 1.0, 4.0
    for _ in range(200):
        r = 0.5 * (a + b)
        tot = h0 * n if abs(r - 1) < 1e-12 else h0 * (r ** n - 1) / (r - 1)
        if tot < hi - lo:
            a = r
        else:
            b = r
    r = 0.5 * (a + b)
    e = np.concatenate([[0.0], np.cumsum(h0 * r ** np.arange(n))])
    e *= (hi - lo) / e[-1]
    return lo + e if fine_at_lo else hi - e[::-1]


def n_graded(length, h0, r=R_GRADE):
    return max(1, int(np.ceil(np.log(1.0 + length * (r - 1.0) / h0)
                              / np.log(r))))


# local node order == LadrunoHex20Shape.h ORDERING MEMO (Abaqus C3D20 /
# upstream Twenty_Node_Brick "Local Node Pattern"), as HALF-INDEX offsets
# (xi, eta, zeta) -> (2*di, 2*dj, 2*dk) with 0/1/2 for -1/0/+1.
H20_OFF = [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
           (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
           (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
           (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
           (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1)]
# the zeta = +1 face: Q8 order [c0 c1 c2 c3 m01 m12 m23 m30] (H20) and the
# plain Q4 corner face (H8) -- the H8 corner block of H20_OFF is exactly
# dp_strip.py's hex ordering, so the two meshes are the same grid.
TOPFACE_LOC = {1: [4, 5, 6, 7], 2: [4, 5, 6, 7, 12, 13, 14, 15]}
Q8_XI = [(-1., -1.), (1., -1.), (1., 1.), (-1., 1.),
         (0., -1.), (1., 0.), (0., 1.), (-1., 0.)]


def strip_mesh(h0=0.5, order=2):
    """dp_strip.py's graded tensor grid, as 8-node (order 1) or 20-node
    serendipity (order 2) hexes -- the SAME grid either way, so an
    element-technology comparison on it changes only the element.

    A half-index lattice (2*n-1 per axis) carries every corner AND every
    edge midpoint; a lattice site is a NODE iff at most one of its three
    half-indices is odd (all even = corner, one odd = edge midpoint, two odd =
    face centre and three odd = body centre, neither of which a serendipity
    element has).  Order 1 simply keeps the all-even sites."""
    nf = int(round(2.0 / h0))                          # elements across B
    no = n_graded(XLIM - 1.0, h0)
    nz = n_graded(-ZBOT - 3.0, h0)
    x = np.unique(np.concatenate([
        graded(-XLIM, -1.0, no, h0, False),            # fine at -1
        np.linspace(-1.0, 1.0, nf + 1),
        graded(1.0, XLIM, no, h0, True)]))             # fine at +1
    y = np.array([0.0, THICK])
    z = np.unique(np.concatenate([
        graded(ZBOT, -3.0, nz, h0, False),             # fine at -3
        np.linspace(-3.0, 0.0, int(round(3.0 / h0)) + 1)]))

    def halfgrid(v):
        out = np.empty(2 * len(v) - 1)
        out[0::2] = v
        out[1::2] = 0.5 * (v[:-1] + v[1:])
        return out

    X, Y, Z = halfgrid(x), halfgrid(y), halfgrid(z)
    NI, NJ, NK = len(X), len(Y), len(Z)
    nen = 8 if order == 1 else 20
    keep = 0 if order == 1 else 1        # max number of odd half-indices
    gid = -np.ones((NI, NJ, NK), dtype=np.int64)
    coords = []
    for i in range(NI):
        for j in range(NJ):
            for k in range(NK):
                if (i & 1) + (j & 1) + (k & 1) > keep:
                    continue                            # face / body centre
                gid[i, j, k] = len(coords)
                coords.append((X[i], Y[j], Z[k]))
    nodes = np.asarray(coords, dtype=float)

    cells = np.array(
        [[gid[2 * i + a, 2 * j + b, 2 * k + c] for (a, b, c) in H20_OFF[:nen]]
         for i in range(len(x) - 1) for j in range(len(y) - 1)
         for k in range(len(z) - 1)], dtype=np.int64)
    assert (cells >= 0).all(), "an element slot resolved to a missing node"

    # positive Jacobian: the corner sub-hex is the H8 test, and the map is
    # axis-aligned monotone by construction, so this must pass exactly.
    p = nodes[cells[:, :8]]
    vol = np.einsum("ij,ij->i", p[:, 6] - p[:, 0],
                    np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
    assert (vol > 0).all(), "inverted element in the strip mesh"

    tol = 1e-9
    top = np.abs(nodes[:, 2]) < tol
    sets = dict(top=np.where(top)[0],
                bottom=np.where(np.abs(nodes[:, 2] - ZBOT) < tol)[0],
                xface=np.where(np.abs(np.abs(nodes[:, 0]) - XLIM) < tol)[0],
                footing=np.where(top & (np.abs(nodes[:, 0]) <= 1.0 + tol))[0])
    # order 2: (nf+1) corner columns x 2 corner rows + nf mid-x columns x 2
    # corner rows + (nf+1) corner columns x 1 mid-y row = 5 nf + 3.
    # order 1: (nf+1) x 2.
    want = 2 * (nf + 1) if order == 1 else 5 * nf + 3
    assert len(sets["footing"]) == want, (len(sets["footing"]), want)
    return nodes, cells, sets, x, y, z


# --- CONSISTENT surcharge on a Q8 face (control 2) --------------------------
def _q8(xi, eta, order=2):
    """Q4 (order 1) / serendipity Q8 (order 2) shape functions + derivatives,
    Q8_XI ordering."""
    if order == 1:
        N, dN = np.zeros(4), np.zeros((4, 2))
        for a, (xa, ea) in enumerate(Q8_XI[:4]):
            N[a] = 0.25 * (1 + xa * xi) * (1 + ea * eta)
            dN[a, 0] = 0.25 * xa * (1 + ea * eta)
            dN[a, 1] = 0.25 * ea * (1 + xa * xi)
        return N, dN
    N, dN = np.zeros(8), np.zeros((8, 2))
    for a, (xa, ea) in enumerate(Q8_XI):
        if xa and ea:                                   # corner
            N[a] = 0.25 * (1 + xa * xi) * (1 + ea * eta) * (xa * xi + ea * eta - 1)
            dN[a, 0] = 0.25 * xa * (1 + ea * eta) * (2 * xa * xi + ea * eta)
            dN[a, 1] = 0.25 * ea * (1 + xa * xi) * (xa * xi + 2 * ea * eta)
        elif xa == 0:                                   # mid-edge on eta = +-1
            N[a] = 0.5 * (1 - xi * xi) * (1 + ea * eta)
            dN[a, 0] = -xi * (1 + ea * eta)
            dN[a, 1] = 0.5 * ea * (1 - xi * xi)
        else:                                           # mid-edge on xi = +-1
            N[a] = 0.5 * (1 + xa * xi) * (1 - eta * eta)
            dN[a, 0] = 0.5 * xa * (1 - eta * eta)
            dN[a, 1] = -eta * (1 + xa * xi)
    return N, dN


_G3 = (-math.sqrt(0.6), 0.0, math.sqrt(0.6))
_W3 = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)


def consistent_surcharge(nodes, cells, order=2):
    """int N_a dA over every z = 0 element face, by 3x3 Gauss on the real
    geometry.  Returns the per-node weights (m^2) whose sum is the top area.

    This is the load the ELEMENT's own basis implies.  For an H8 face it
    coincides with the tributary area (A/4 each); for an H20 face it does NOT
    -- the four corner weights are NEGATIVE (-A/12) and the mid-edges carry
    +A/3.  Same code path either way, so the two legs cannot diverge by a
    hand-written special case."""
    loc = TOPFACE_LOC[order]
    w = np.zeros(len(nodes))
    nfaces = 0
    for conn in cells:
        face = conn[loc]
        p = nodes[face]
        if np.any(np.abs(p[:, 2]) > 1e-9):
            continue
        nfaces += 1
        for xi, wx in zip(_G3, _W3):
            for eta, we in zip(_G3, _W3):
                N, dN = _q8(xi, eta, order)
                g1, g2 = dN[:, 0] @ p, dN[:, 1] @ p
                jac = np.linalg.norm(np.cross(g1, g2))
                w[face] += (wx * we * jac) * N
    return w, nfaces


def verify_surcharge(nodes, cells, w, nfaces, order=2):
    """Controls 2a/2b: the weights match the analytic rectangle values, and
    they sum to the physical top area."""
    loc = TOPFACE_LOC[order]
    area = 2 * XLIM * THICK
    assert abs(w.sum() - area) < 1e-8 * area, (w.sum(), area)
    nneg = int((w < 0).sum())
    # per-face analytic check on one representative face
    conn = next(c for c in cells if np.all(np.abs(nodes[c[loc], 2]) < 1e-9))
    p = nodes[conn[loc]]
    A = (p[:, 0].max() - p[:, 0].min()) * (p[:, 1].max() - p[:, 1].min())
    wl = np.zeros(len(loc))
    for xi, wx in zip(_G3, _W3):
        for eta, we in zip(_G3, _W3):
            N, dN = _q8(xi, eta, order)
            g1, g2 = dN[:, 0] @ p, dN[:, 1] @ p
            wl += (wx * we * np.linalg.norm(np.cross(g1, g2))) * N
    want = (np.array([A / 4] * 4) if order == 1
            else np.array([-A / 12] * 4 + [A / 3] * 4))
    err = float(np.abs(wl - want).max() / A)
    assert err < 1e-12, (wl, want)
    shape = "Q4 (= tributary)" if order == 1 else "Q8 serendipity"
    log(f"    [control 2] consistent {shape} surcharge over {nfaces} top "
        f"faces: sum {w.sum():.9f} m^2 vs {area:.9f} exact; {nneg} nodes "
        f"carry a NEGATIVE weight; per-face weights match the analytic "
        f"{'(A/4)' if order == 1 else '(-A/12, +A/3)'} to {err:.1e} of A")
    return nneg


# --- the model --------------------------------------------------------------
def build_model(nodes, cells, sets, w, form, assoc, elastic=False):
    alpha = alpha_from_phi_txc(PHI_TXC)
    rho = math.sqrt(2.0) * alpha
    g_el = 3.0 * K_EL * (1.0 - 2.0 * NU) / (2.0 * (1.0 + NU))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))
    if elastic:
        e_el = 9.0 * K_EL * g_el / (3.0 * K_EL + g_el)
        ops.nDMaterial("ElasticIsotropic", 1, e_el, NU, 0.0)
    else:
        ops.nDMaterial("DruckerPrager", 1, K_EL, g_el, SY, rho,
                       rho if assoc else 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for e, conn in enumerate(cells, start=1):
        if cells.shape[1] == 20:
            ops.element("LadrunoBrick20", e, *[int(c) + 1 for c in conn], 1,
                        "-formulation", form, "-b", 0.0, 0.0, -GAMMA)
        else:
            # the LINEAR-hex control, on the identical grid: `-geom linear`
            # and `-formulation bbar` are dp_strip.py's baseline settings.
            ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                        "-geom", "linear", "-b", 0.0, 0.0, -GAMMA,
                        "-formulation", form)
    # ONE fix() per node -- a second fix() on the same node adds a DUPLICATE
    # SP_Constraint rather than replacing it.  u_y = 0 on EVERY node is what
    # makes this plane strain, and it is exact for any element order.
    bt, xf = set(sets["bottom"].tolist()), set(sets["xface"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if n in xf else 0, 1, 0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if abs(w[n]) > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -Q0 * float(w[n]))
    return alpha, g_el


def surcharge_step(nodes, sets, w, tol):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso")      # tangent is UNSYMMETRIC when rho_bar != rho
        sysname = "Pardiso"
    except Exception:
        ops.system("UmfPack")
        sysname = "UmfPack"
    ops.test("NormUnbalance", tol, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    t0 = time.time()
    if ops.analyze(1) != 0:
        raise SystemExit("surcharge step failed")
    ops.reactions()
    rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    want = Q0 * float(w.sum())
    log(f"    [control 3] initial equilibrium {rz:.6f} vs {want:.6f} kN "
        f"({100 * (rz / want - 1):+.7f} %)   [{sysname}, {time.time() - t0:.1f}s]")
    assert abs(rz / want - 1) < 1e-6, "surcharge / basis convention"
    return sysname


# --- control 1 (the strong one): the 1-D patch test -------------------------
def leg_patch(nodes, cells, sets, w, form):
    """ElasticIsotropic on the same mesh.  Roller sides + fixed base + a
    consistent uniform surcharge admits the EXACT 1-D state, and that state
    lives in the H20 space, so a consistent load vector must reproduce it at
    every Gauss point.  A tributary-lumped load would not: it would show a
    boundary layer under the surface.  This is the check that tests the
    DISTRIBUTION; the base-reaction identity only tests the sum."""
    build_model(nodes, cells, sets, w, form, assoc=False, elastic=True)
    surcharge_step(nodes, sets, w, 1.0e-5 * Q0 * float(w.sum()))
    k0 = NU / (1.0 - NU)
    worst = dict(zz=0.0, xx=0.0, sh=0.0)
    for e in range(1, len(cells) + 1):
        st = np.asarray(ops.eleResponse(e, "stress"))
        if st.size == 0:
            raise SystemExit("no stress response from LadrunoBrick20")
        st = st.reshape(-1, 6)
        worst["zz"] = max(worst["zz"], np.abs(st[:, 2] / Q0 + 1.0).max())
        worst["xx"] = max(worst["xx"], np.abs(st[:, 0] / (k0 * Q0) + 1.0).max())
        worst["sh"] = max(worst["sh"], np.abs(st[:, 3:]).max() / Q0)
    ngp = st.shape[0]
    log(f"    [control 1] 1-D elastic patch test on {len(cells)} "
        f"{cells.shape[1]}-node elements x "
        f"{ngp} GP ({form}): max |sigma_zz/-q0 - 1| = {worst['zz']:.2e}, "
        f"max |sigma_xx/-K0 q0 - 1| = {worst['xx']:.2e}, "
        f"max |tau|/q0 = {worst['sh']:.2e}")
    ok = worst["zz"] < 1e-9 and worst["xx"] < 1e-9 and worst["sh"] < 1e-9
    log(f"    [control 1] {'PASS' if ok else '*** FAIL ***'} "
        f"-- the surcharge is consistent with the H20 basis"
        if ok else "    [control 1] *** FAIL: the surcharge is NOT consistent "
        "with the H20 basis ***")
    if not ok:
        raise SystemExit("patch test failed; every downstream number is void")
    return dict(leg="patch", ngp=ngp, worst=worst)


# --- control 4: the spurious-mode census on THIS mesh -----------------------
def leg_modes(nodes, cells, sets, w, form):
    """Count the zero-energy modes of the ASSEMBLED, RESTRAINED elastic
    stiffness of this exact mesh.

    Why this is not optional here.  H20 at 2x2x2 carries 6 spurious modes per
    element and ships with NO hourglass control by design (ADR 72 §3.2; Abaqus
    applies none to C3D20R).  The contract is that they are non-communicable in
    a solid assembly -- and the usage guidance that comes with it says `uri` is
    for meshes AT LEAST TWO ELEMENTS THICK, while this plane-strain strip is
    ONE element thick in y.  So the leg sits exactly on the documented caveat,
    and the honest move is to measure the rank rather than argue about it.  If
    the restrained stiffness has any zero eigenvalue beyond round-off, the
    `uri` collapse load is a mechanism, not a measurement."""
    build_model(nodes, cells, sets, w, form, assoc=False, elastic=True)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-8, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    ops.analyze(1)
    a = np.asarray(ops.printA("-ret"), dtype=float)
    n = int(round(math.sqrt(a.size)))
    assert n * n == a.size, a.size
    K = a.reshape(n, n)
    K = 0.5 * (K + K.T)
    ev = np.linalg.eigvalsh(K)
    tr = float(np.trace(K)) / n                    # mean diagonal = the scale
    rel = ev / tr
    nzero = int((rel < 1e-10).sum())
    log(f"    [control 4] spurious-mode census ({form}): {n} free DOF, "
        f"smallest 5 eigenvalues / (tr K / n) = "
        f"{', '.join(f'{v:.3e}' for v in rel[:5])}")
    log(f"    [control 4] zero-energy modes (rel < 1e-10): {nzero}  ->  "
        + ("PASS -- the assembly has full rank, so the reduced rule's 6 "
           "per-element modes do NOT communicate through this "
           "single-element-thick strip"
           if nzero == 0 else
           "*** FAIL: the mesh carries a mechanism; any collapse load "
           "measured on it is void ***"))
    if nzero:
        raise SystemExit("spurious mechanism present")
    return dict(leg="modes", form=form, ndof=n, lam_min=float(rel[0]))


# --- the ladder -------------------------------------------------------------
def attempts(tol):
    """A perfectly plastic collapse is where Newton is worst: the tangent goes
    singular in the mechanism mode exactly when the answer arrives.  Each
    increment gets this ladder before the step is shrunk, and any step that
    needed the relaxed tolerance is flagged in the CSV."""
    return [("KrylovNewton", tol, 25, 0),
            ("NewtonLineSearch", tol, 40, 0),
            ("KrylovNewton", 10.0 * tol, 60, 1)]


def set_attempt(a):
    algo, tol, it, _ = a
    ops.test("NormUnbalance", tol, it, 0)
    ops.algorithm(algo)


def leg_push(nodes, cells, sets, w, form, assoc, tag, budget, sfrac, tmax, h0):
    alpha = alpha_from_phi_txc(PHI_TXC)
    mc = mc_from_cone(alpha)
    nq = n_q(mc["ps"])
    q_exact = Q0 * nq
    k0 = NU / (1.0 - NU)
    m0 = (1.0 - k0) / (math.sqrt(3.0) * alpha * (1.0 + 2.0 * k0))

    ename = "LadrunoBrick20" if cells.shape[1] == 20 else "LadrunoBrick"
    log(f"=== leg {tag}: {ename} -formulation {form}, h0 = {h0} m "
        f"({int(round(2.0 / h0))} elements across B), {len(nodes)} nodes / "
        f"{3 * len(nodes)} DOF, {len(cells)} elements")
    log(f"    cone alpha = {alpha:.6f}  ->  phi_txc = {mc['txc']:.3f} deg, "
        f"phi_ps = {mc['ps']:.3f} deg;  "
        f"{'ASSOCIATED (control)' if assoc else 'non-associated (gate)'}")
    log(f"    ORACLE  q_u = q0 * N_q = {Q0} * {nq:.4f} = {q_exact:.3f} kPa "
        f"(exact: upper and lower bounds coincide)")
    log(f"    Davis psi=0 would give phi* = {davis_phi(mc['ps']):.3f} deg, "
        f"q_u = {Q0 * n_q(davis_phi(mc['ps'])):.1f} kPa -- a lower bound, "
        f"reported not used")
    log(f"    initial 1-D state at m0 = {m0:.4f} of yield "
        + ("(SAFE)" if m0 < 0.8 else "*** TOO CLOSE TO YIELD -- LEG IS VOID ***"))
    if m0 >= 0.8:
        raise SystemExit("initial state too close to yield; raise nu")
    log(f"    subdivision budget = {budget}, ds base/min/max = "
        f"{DS_BASE * 1e3}/{DS_MIN * 1e3}/{DS_MAX * 1e3} mm, s/B target = {sfrac}")

    build_model(nodes, cells, sets, w, form, assoc)
    tol = 1.0e-5 * max(Q0 * float(w.sum()), 1.0)
    surcharge_step(nodes, sets, w, tol)

    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = Q0 * float(w[sets["footing"]].sum())   # add back what we applied
    area = B_FOOT * THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = sfrac * B_FOOT
    # pseudo-time IS the imposed footing settlement (r3's Linear-series form,
    # which survives the trip through openseespy where the Path form does not)
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso")
    except Exception:
        ops.system("UmfPack")
    ops.analysis("Static")

    path = os.path.join(HERE, f"h20_{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "R_kN", "q_kPa", "ds_mm", "relaxed", "wall_s"])
    rows, ds, good, nfail, nrelax, nsub = [], DS_BASE, 0, 0, 0, 0
    verdict, t0 = "reached the target", time.time()
    lad = attempts(tol)
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > tmax:
            verdict = (f"WALL-CLOCK CAP of {tmax:.0f} s hit at "
                       f"s/B = {s_now / B_FOOT:.4f}")
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for a in lad:
            set_attempt(a)
            if ops.analyze(1) == 0:
                ok, relaxed = True, a[3]
                break
            nfail += 1
        if not ok:
            good = 0
            nsub += 1
            ds *= 0.5
            if nsub > budget:
                verdict = (f"subdivision budget of {budget} spent at "
                           f"s/B = {s_now / B_FOOT:.4f}")
                break
            if ds < DS_MIN:
                verdict = (f"no convergence at s/B = {s_now / B_FOOT:.4f} "
                           f"(every ladder rung failed at ds = {ds * 1e3:.5f} mm)")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= GROW_AFTER and ds < DS_MAX:          # the recovery latch
            ds, good = min(2 * ds, DS_MAX), 0
        ops.reactions()
        R = -sum(ops.nodeReaction(t, 3) for t in foot) + r_corr
        s = uz0 - ops.getTime()
        rows.append((s, s / B_FOOT, R, R / area, ds * 1e3, relaxed,
                     time.time() - t0))
        wr.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()
        if len(rows) > STALL_WINDOW:                     # the stall detector
            if rows[-1][0] - rows[-1 - STALL_WINDOW][0] < STALL_ADVANCE * smax:
                verdict = (f"stalled: < {100 * STALL_ADVANCE:.2f} % of the "
                           f"target in {STALL_WINDOW} steps, at "
                           f"s/B = {s / B_FOOT:.4f}")
                break
    fh.close()
    wall = time.time() - t0

    a = np.array([r[:4] for r in rows])
    s, q = a[:, 0], a[:, 3]
    # D16's truncation rule: a bottomed-out run writes its collapse as a final
    # converged step, which then reads as the answer.
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    if len(bad):
        ntr = len(q) - bad[0]
        s, q = s[:bad[0]], q[:bad[0]]
        verdict += f"; curve TRUNCATED at s/B = {s[-1] / B_FOOT:.4f} ({ntr} cut)"

    msk = s >= 0.9 * s[-1]
    t_last = np.polyfit(s[msk], q[msk], 1)[0] if msk.sum() > 2 else float("nan")
    n0 = max(4, len(s) // 50)
    t_init = np.polyfit(s[:n0], q[:n0], 1)[0]
    plateau = bool(abs(t_last) < 0.02 * abs(t_init))
    qmax = float(q.max())

    log(f"    {len(rows)} steps, {nfail} failed attempts, {nsub}/{budget} "
        f"subdivisions, {nrelax} relaxed, {wall:.0f} s   [{verdict}]")
    log(f"    q_max = {qmax:.2f} kPa at s/B = "
        f"{s[int(q.argmax())] / B_FOOT:.4f}; end s/B = {s[-1] / B_FOOT:.4f}")
    log(f"    tail dq/ds = {t_last:.1f} kPa/m = {100 * t_last / t_init:.2f} % "
        f"of initial  ->  {'PLATEAU' if plateau else 'STILL HARDENING'}")
    log(f"    q_num / q_exact = {qmax / q_exact:.4f}   "
        f"[{qmax:.2f} vs {q_exact:.2f} kPa]")
    return dict(leg=tag, elem=ename, form=form, h0=h0, assoc=assoc,
                nel=len(cells),
                ndof=3 * len(nodes), qmax=qmax, q_exact=q_exact,
                ratio=qmax / q_exact, plateau=plateau, m0=m0,
                t_frac=100 * t_last / t_init, s_end=float(s[-1]) / B_FOOT,
                nsteps=len(rows), nsub=nsub, budget=budget, wall=wall,
                verdict=verdict)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="*", default=["patch", "modes", "nonassoc"])
    ap.add_argument("--h0", type=float, default=0.5)
    ap.add_argument("--order", type=int, default=2, choices=[1, 2],
                    help="1 = LadrunoBrick (H8) control, 2 = LadrunoBrick20")
    ap.add_argument("--form", default="uri",
                    choices=["uri", "std", "bbar"],
                    help="H20: uri|std.  H8: bbar|std|uri.")
    ap.add_argument("--sfrac", type=float, default=0.25)
    ap.add_argument("--budget", type=int, default=40)
    ap.add_argument("--tmax", type=float, default=1.0e9)
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    if args.order == 2 and args.form == "bbar":
        raise SystemExit("LadrunoBrick20 has no B-bar; use uri (the C3D20R "
                         "analog) or std")

    # --- control 0: provenance ---------------------------------------------
    try:
        stamp = ops.ladrunoBuild()
    except Exception as exc:                             # pragma: no cover
        raise SystemExit(f"ladrunoBuild() unavailable ({exc}) -- this binary "
                         f"predates the build stamp; refusing to measure")
    log(f"[control 0] engine {ops.__file__}")
    log(f"[control 0] ladrunoBuild() = {stamp}")

    nodes, cells, sets, x, y, z = strip_mesh(args.h0, args.order)
    w, nfaces = consistent_surcharge(nodes, cells, args.order)
    verify_surcharge(nodes, cells, w, nfaces, args.order)
    log(f"[mesh] {len(nodes)} nodes / {3 * len(nodes)} DOF, {len(cells)} "
        f"{cells.shape[1]}-node hexes "
        f"({len(x) - 1} x {len(y) - 1} x {len(z) - 1}), footing set "
        f"{len(sets['footing'])} nodes")

    out = []
    for leg in args.legs:
        if leg == "patch":
            leg_patch(nodes, cells, sets, w, args.form)
        elif leg == "modes":
            leg_modes(nodes, cells, sets, w, args.form)
        elif leg in ("nonassoc", "assoc"):
            tag = (f"h{cells.shape[1]}_{leg}_{args.form}_"
                   f"h{str(args.h0).replace('0.', '')}" + args.suffix)
            out.append(leg_push(nodes, cells, sets, w, args.form,
                                leg == "assoc", tag, args.budget, args.sfrac,
                                args.tmax, args.h0))
        else:
            raise SystemExit(f"unknown leg {leg!r}")

    if out:
        print()
        hdr = (f"{'leg':>28} {'form':>5} {'h0':>6} {'q_num':>9} {'exact':>9} "
               f"{'num/exact':>10} {'tail %':>8} {'plateau':>8} {'s_end/B':>8} "
               f"{'subdiv':>8} {'wall_s':>8}")
        print(hdr)
        print("-" * len(hdr))
        for r in out:
            print(f"{r['leg']:>28} {r['form']:>5} {r['h0']:6.3f} "
                  f"{r['qmax']:9.2f} {r['q_exact']:9.2f} {r['ratio']:10.4f} "
                  f"{r['t_frac']:8.2f} {'yes' if r['plateau'] else 'NO':>8} "
                  f"{r['s_end']:8.4f} {str(r['nsub']) + '/' + str(r['budget']):>8} "
                  f"{r['wall']:8.0f}")


if __name__ == "__main__":
    main()
