"""Dynamic Rayleigh + inertia regression for BezierTri6, BezierTet10 and
LadrunoIMKBeam (2D/3D) -- WP-115.

These four elements had NO transient test with mass and Rayleigh damping, so
nothing would have noticed getResistingForceIncInertia dropping the inertia,
the -Q ground-motion term, or the Rayleigh force itself (LEDGER_quirks
"MUST snapshot the shared static `resid`": the #562 clobber recurred across
a whole family because the quirk had no per-element gate). WP-115 converted
all four to the snapshot idiom; these gates pin the result.

Gates:
  IMK  differential vs elasticBeamColumn carrying the SAME lumped masses as
       NODAL masses (rho*h per interior node, rho*h/2 at the tip) -- with no
       hinge materials IMK is the elastic beam, so every displacement history
       must agree. Legs: {alphaM, betaK, betaK0} x {IMK element -mass, IMK
       nodal mass} x {step tip load, UniformExcitation}. Fails in BOTH
       directions: dropped inertia, dropped Rayleigh, or a wrong -Q term all
       break the agreement. The oracle uses nodal masses on purpose: vanilla
       ElasticBeam2d with element -mass subtracts Q twice under
       UniformExcitation (LEDGER_quirks "ElasticBeam2d subtracts the
       ground-motion load Q TWICE"), so it is not a valid -Q reference.
  Bezier self-validating step-load rig (the test_ladrunoplane_dynamics
       pattern): the undamped peak must show the ~2x step overshoot (proves
       the rig resolves dynamics); a tiny betaK must move it <5% (catches
       dropped inertia); and 5% modal damping must pull it into a band a
       dropped Rayleigh force cannot reach (~1.0 undamped vs ~0.93 damped).

Break-on-purpose evidence (WP-115 plan doc): each leg was run against a
binary with the inertia add removed, and one with the Rayleigh add removed.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

XI = 0.05                     # target modal damping for the "Rayleigh is applied" legs
DAMPED_BAND = (0.80, 0.97)    # peak_damped / peak_undamped; SDOF estimate ~0.93


def _transient():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 50)
    ops.analysis("Transient")


def _omega1():
    lam = ops.eigen("-fullGenLapack", 1)[0]
    return math.sqrt(lam)


def _rayleigh_for(kind, w1, xi):
    """(alphaM, betaK, betaK0, betaKc) giving modal damping xi in mode 1."""
    if kind == "alphaM":
        return (2.0 * xi * w1, 0.0, 0.0, 0.0)
    if kind == "betaK":
        return (0.0, 2.0 * xi / w1, 0.0, 0.0)
    if kind == "betaK0":
        return (0.0, 0.0, 2.0 * xi / w1, 0.0)
    raise ValueError(kind)


# ==========================================================================
# LadrunoIMKBeam (2D and 3D): differential against elasticBeamColumn
# ==========================================================================
L_BEAM, A_SEC, E_MOD, I_SEC = 3.0, 0.02, 2.0e5, 1.5e-4
G_MOD, J_TOR = 8.0e4, 3.0e-4
RHO_L = 0.8                    # mass per unit length
NEL = 4


def _imk_model(dim, elem, mass_src, damping, load):
    """Cantilever along x, NEL elements; returns the tip-node y-displacement
    history over ~1.5 periods of mode 1."""
    ops.wipe()
    if dim == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.geomTransf("Linear", 1)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    h = L_BEAM / NEL
    for i in range(NEL + 1):
        if dim == 2:
            ops.node(i + 1, i * h, 0.0)
        else:
            ops.node(i + 1, i * h, 0.0, 0.0)
    ops.fix(1, *([1] * (3 if dim == 2 else 6)))

    m_el = RHO_L if mass_src == "element" else 0.0
    for e in range(1, NEL + 1):
        if dim == 2:
            sec = (A_SEC, E_MOD, I_SEC)
        else:
            sec = (A_SEC, E_MOD, G_MOD, J_TOR, I_SEC, I_SEC)
        name = "LadrunoIMKBeam" if elem == "imk" else "elasticBeamColumn"
        args = [name, e, e, e + 1, *sec, 1]
        if m_el > 0.0:
            args += ["-mass", m_el]
        ops.element(*args)
    if mass_src == "nodal":
        for i in range(2, NEL + 2):
            m = RHO_L * h * (0.5 if i == NEL + 1 else 1.0)
            ops.mass(i, *([m, m, 0.0] if dim == 2 else [m, m, m, 0.0, 0.0, 0.0]))

    ops.timeSeries("Constant", 1)
    if load == "step":
        ops.pattern("Plain", 1, 1)
        ops.load(NEL + 1, *([0.0, 1.0, 0.0] if dim == 2 else [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]))
    else:  # ground acceleration in y; elements carry -Q, nodes carry nodal mass
        ops.pattern("UniformExcitation", 1, 2, "-accel", 1)

    _transient()
    w1 = _omega1()
    ops.rayleigh(*_rayleigh_for(damping, w1, XI))
    nstep = 120
    dt = 1.5 * (2.0 * math.pi / w1) / nstep
    hist = []
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0, f"{elem} transient step failed"
        hist.append(ops.nodeDisp(NEL + 1, 2))
    return hist


@pytest.mark.parametrize("dim", [2, 3])
def test_imk_static_equals_elastic_beam(dim):
    """Premise of the differential: with no hinge materials IMK is the
    elastic beam (same tip deflection under a static tip load)."""
    tips = {}
    for elem in ("imk", "ebc"):
        ops.wipe()
        if dim == 2:
            ops.model("basic", "-ndm", 2, "-ndf", 3)
            ops.geomTransf("Linear", 1)
        else:
            ops.model("basic", "-ndm", 3, "-ndf", 6)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        h = L_BEAM / NEL
        for i in range(NEL + 1):
            ops.node(i + 1, *((i * h, 0.0) if dim == 2 else (i * h, 0.0, 0.0)))
        ops.fix(1, *([1] * (3 if dim == 2 else 6)))
        sec = (A_SEC, E_MOD, I_SEC) if dim == 2 else (A_SEC, E_MOD, G_MOD, J_TOR, I_SEC, I_SEC)
        name = "LadrunoIMKBeam" if elem == "imk" else "elasticBeamColumn"
        for e in range(1, NEL + 1):
            ops.element(name, e, e, e + 1, *sec, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(NEL + 1, *([0.0, 1.0, 0.0] if dim == 2 else [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]))
        ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
        assert ops.analyze(1) == 0
        tips[elem] = ops.nodeDisp(NEL + 1, 2)
    exact = 1.0 * L_BEAM ** 3 / (3.0 * E_MOD * I_SEC)
    assert math.isclose(tips["ebc"], exact, rel_tol=1e-9)
    assert math.isclose(tips["imk"], tips["ebc"], rel_tol=1e-9), tips


@pytest.mark.parametrize("load", ["step", "ground"])
@pytest.mark.parametrize("mass_src", ["element", "nodal"])
@pytest.mark.parametrize("damping", ["alphaM", "betaK", "betaK0"])
@pytest.mark.parametrize("dim", [2, 3])
def test_imk_dynamic_rayleigh_matches_elastic_beam(dim, damping, mass_src, load):
    imk = _imk_model(dim, "imk", mass_src, damping, load)
    ebc = _imk_model(dim, "ebc", "nodal", damping, load)   # see module docstring
    scale = max(abs(u) for u in ebc)
    assert scale > 0.0, "reference beam did not move -- rig is not driving the tip"
    worst = max(abs(a - b) for a, b in zip(imk, ebc))
    assert worst <= 1e-7 * scale, (
        f"IMK {dim}D [{damping}, {mass_src} mass, {load}] departs from "
        f"elasticBeamColumn by {worst / scale:.3e} of peak -- inertia, -Q or "
        "the Rayleigh force is being dropped or mis-signed in "
        "getResistingForceIncInertia / addInertiaLoadToUnbalance")


# ==========================================================================
# BezierTri6 / BezierTet10: self-validating step-load rig
# ==========================================================================
E_S, NU_S, RHO_S = 1000.0, 0.3, 8.0
NCELL = 2


def _tri6_column():
    """1 x NCELL column of unit cells, two BezierTri6 per cell."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E_S, NU_S)
    coord, key, mids = {}, {}, {}

    def node(xy):
        k = tuple(round(c, 9) for c in xy)
        if k not in key:
            t = len(coord) + 1
            key[k], coord[t] = t, xy
            ops.node(t, *xy)
        return key[k]

    def mid(a, b):
        k = frozenset((a, b))
        if k not in mids:
            ca, cb = coord[a], coord[b]
            mids[k] = node((0.5 * (ca[0] + cb[0]), 0.5 * (ca[1] + cb[1])))
        return mids[k]

    e = 0
    for j in range(NCELL):
        n1, n2 = node((0.0, float(j))), node((1.0, float(j)))
        n3, n4 = node((1.0, j + 1.0)), node((0.0, j + 1.0))
        for c1, c2, c3 in ((n1, n2, n3), (n1, n3, n4)):
            e += 1
            ops.element("BezierTri6", e, c1, c2, c3, mid(c1, c2), mid(c2, c3), mid(c3, c1),
                        1.0, "PlaneStrain", 1, "-rho", RHO_S)
    base = [t for t, (x, y) in coord.items() if abs(y) < 1e-12]
    top = [key[(0.0, float(NCELL))], key[(1.0, float(NCELL))]]
    for t in base:
        ops.fix(t, 1, 1)
    return top, 2


_EDGES_TET = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]   # e12 e23 e13 e14 e34 e24


def _tet10_column():
    """1 x 1 x NCELL column of unit cubes, six (Kuhn) BezierTet10 per cube."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_S, NU_S)
    coord, key, mids = {}, {}, {}

    def node(xyz):
        k = tuple(round(c, 9) for c in xyz)
        if k not in key:
            t = len(coord) + 1
            key[k], coord[t] = t, tuple(float(c) for c in xyz)
            ops.node(t, *coord[t])
        return key[k]

    def mid(a, b):
        k = frozenset((a, b))
        if k not in mids:
            ca, cb = coord[a], coord[b]
            mids[k] = node(tuple(0.5 * (ca[i] + cb[i]) for i in range(3)))
        return mids[k]

    def vol(v):
        p = [coord[t] for t in v]
        a = [p[1][i] - p[0][i] for i in range(3)]
        b = [p[2][i] - p[0][i] for i in range(3)]
        c = [p[3][i] - p[0][i] for i in range(3)]
        return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))

    e = 0
    for k in range(NCELL):
        for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
            p, path = [0, 0, 0], []
            path.append(node((p[0], p[1], p[2] + k)))
            for ax in perm:
                p[ax] = 1
                path.append(node((p[0], p[1], p[2] + k)))
            if vol(path) < 0.0:
                path[1], path[2] = path[2], path[1]
            conn = path + [mid(path[a], path[b]) for a, b in _EDGES_TET]
            e += 1
            ops.element("BezierTet10", e, *conn, 1, "-rho", RHO_S)
    base = [t for t, (x, y, z) in coord.items() if abs(z) < 1e-12]
    top = [key[(float(x), float(y), float(NCELL))] for x in (0, 1) for y in (0, 1)]
    for t in base:
        ops.fix(t, 1, 1, 1)
    return top, 3


_BUILDERS = {"BezierTri6": _tri6_column, "BezierTet10": _tet10_column}


def _bezier_run(kind, dynamic, damping=None, xi=0.0):
    """Tip x-displacement: static value, or the peak over ~1.5 periods of a
    step load applied at t = 0 (Constant series)."""
    top, ndf = _BUILDERS[kind]()
    ops.timeSeries("Constant" if dynamic else "Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 1.0, *([0.0] * (ndf - 1)))
    if not dynamic:
        ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
        assert ops.analyze(1) == 0
        return abs(ops.nodeDisp(top[0], 1))
    _transient()
    w1 = _omega1()
    if damping is not None:
        ops.rayleigh(*_rayleigh_for(damping, w1, xi))
    nstep = 150
    dt = 1.5 * (2.0 * math.pi / w1) / nstep
    peak = 0.0
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0, f"{kind} transient step failed"
        peak = max(peak, abs(ops.nodeDisp(top[0], 1)))
    return peak


@pytest.mark.parametrize("kind", ["BezierTri6", "BezierTet10"])
def test_bezier_dynamic_rayleigh_preserves_inertia(kind):
    static = _bezier_run(kind, dynamic=False)
    peak0 = _bezier_run(kind, dynamic=True)
    overshoot = peak0 / static
    assert 1.5 < overshoot < 2.2, (
        f"[{kind}] undamped peak/static = {overshoot:.3f}; expected the ~2x step "
        "overshoot -- rig is not resolving dynamics, the gate cannot discriminate")
    peak_tiny = _bezier_run(kind, dynamic=True, damping="betaK", xi=1.0e-4)
    rel = abs(peak_tiny - peak0) / peak0
    assert rel < 0.05, (
        f"[{kind}] a tiny betaK moved the peak by {rel:.1%} -- inertia is being "
        "dropped in getResistingForceIncInertia (the #562 clobber)")


@pytest.mark.parametrize("damping", ["alphaM", "betaK", "betaK0"])
@pytest.mark.parametrize("kind", ["BezierTri6", "BezierTet10"])
def test_bezier_rayleigh_force_is_applied(kind, damping):
    peak0 = _bezier_run(kind, dynamic=True)
    peak_d = _bezier_run(kind, dynamic=True, damping=damping, xi=XI)
    ratio = peak_d / peak0
    lo, hi = DAMPED_BAND
    assert lo < ratio < hi, (
        f"[{kind}, {damping}] 5% damping gave peak ratio {ratio:.3f}, expected "
        f"{lo}-{hi} -- ~1.0 means the Rayleigh force is not reaching the residual")
