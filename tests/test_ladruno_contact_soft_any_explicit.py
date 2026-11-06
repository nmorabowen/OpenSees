"""The contact SOFT penalties under every explicit integrator, not only CentralDifferenceLadruno.

`-soft` (node-to-segment and rigid plane, SOFT=1), `-mortar -soft` (SOFT=2, with its friction and
`-visc`) and `-edgeSoft` (edge pairs, normal and stick) replace the configured penalty by
k = s 4 m_eff / dt^2 under an explicit integrator, so the contact mode has omega dt = 2 sqrt(s).
They used to be gated by a dynamic_cast to CentralDifferenceLadruno: under the stock
CentralDifference, NewmarkExplicit, ExplicitDifference or ExplicitBathe the cast failed and the stiff
configured penalty ran silently -- exactly where an explicit contact run diverges.

The adapter now asks the integrator to form its tangent with the add*ToTang calls recording: an
integrator that requests no stiffness is explicit for the contact. The step size is
Domain::getDT(), or the last committed step for an integrator that forms the residual at the
committed time (stock CentralDifference); CentralDifferenceLadruno keeps getCurrentDeltaT().

Every case runs at the structural step dt = 1e-3 with a stiff penalty 1e8 on unit masses
(omega dt = 10): the stiff run diverges, the soft one stays bounded. Under statics and implicit
Newmark `-soft` / `-edgeSoft` stay bitwise inert.

Ported from the upstream port branch up/64 (tests/test_ladruno_contact_soft.py, the non-CDL
integrator cases), extended to the lanes the fork has and up/64 lacks: SOFT=2 friction, `-visc` on
SOFT=2 and the edge stick penalty.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

M = 1.0
DT = 1.0e-3
KN_STIFF = 1.0e8
RUNAWAY = 1.0

# the explicit integrators the old dynamic_cast gate missed, plus CDL as the regression anchor
EXPLICIT = ["NewmarkExplicit", "CentralDifference", "ExplicitDifference", "ExplicitBathe"]
EXPLICIT_ALL = EXPLICIT + ["CentralDifferenceLadruno"]
# stock CentralDifference and ExplicitDifference start from Ut-1 = Ut (no initial velocity), so a
# kick has no effect there; the load-driven cases start those at rest
AT_REST = {"CentralDifference", "ExplicitDifference"}


def _explicit(integrator):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    if integrator == "NewmarkExplicit":
        ops.integrator("NewmarkExplicit", 0.5)
    else:
        ops.integrator(integrator)
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _run(node, dof, nsteps, dt=DT):
    """Integrate; returns (diverged, max |u|)."""
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return True, maxabs
        u = ops.nodeDisp(node, dof)
        if not math.isfinite(u):
            return True, maxabs
        maxabs = max(maxabs, abs(u))
        if maxabs > RUNAWAY:
            return True, maxabs
    return False, maxabs


def _kick(integrator, v):
    return 0.0 if integrator in AT_REST else v


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _implicit_newmark(nsteps=20):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(nsteps):
        assert ops.analyze(1, DT) == 0


# ---------------------------------------------------------------- rigid plane (SOFT=1)
def _plane_loaded(soft, integrator, P=1.0e3, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    v = _kick(integrator, -0.5)
    if v != 0.0:
        ops.setNodeVel(1, 3, v, "-commit")
    _explicit(integrator)
    return _run(1, 3, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT_ALL)
def test_soft_plane_bounded(integrator):
    div_stiff, _ = _plane_loaded(0.0, integrator)
    div_soft, maxabs = _plane_loaded(0.1, integrator)
    assert div_stiff, "the stiff penalty should diverge at the structural step"
    assert not div_soft and maxabs < 0.05, maxabs


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_soft_uses_element_mass(integrator):
    """No nodal mass: the soft penalty is sized from the element mass (a truss with density)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 1.0)
    ops.fix(2, 1, 1, 1)
    ops.fix(1, 1, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 1.0e-3)
    ops.element("truss", 1, 1, 2, 1.0, 1, "-rho", 2.0)    # lumped mass 1 at node 1
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, "-soft", 0.1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    _explicit(integrator)
    div, maxabs = _run(1, 3, 3000)
    assert not div and maxabs < 0.05, maxabs


def _plane_impact(soft, integrator, v_mag=2.0, gap0=0.02, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, *extra)
    ops.setNodeVel(1, 3, -v_mag, "-commit")
    _explicit(integrator)
    v_out = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            return True, 0.0
        z = gap0 + ops.nodeDisp(1, 3)
        if not math.isfinite(z) or abs(z) > 1.0e6:
            return True, 0.0
        if z > 0.5 * gap0 and ops.nodeVel(1, 3) > 0.0:
            v_out = ops.nodeVel(1, 3)
    return False, abs(v_out / v_mag)


# ExplicitBathe (p = 0.54) damps the high frequencies by design; at omega dt = 2 sqrt(0.1) of the
# soft contact mode it removes a few percent of the impact energy (measured e = 0.925)
@pytest.mark.parametrize("integrator,tol", [("NewmarkExplicit", 0.05), ("ExplicitBathe", 0.10)])
def test_soft_plane_impact_restitution(integrator, tol):
    blew, e_soft = _plane_impact(0.1, integrator)
    assert not blew and abs(1.0 - e_soft) < tol, e_soft
    blew, e_stiff = _plane_impact(0.0, integrator)
    assert blew or e_stiff > 2.0, e_stiff


# ------------------------------------------------------------ node-to-segment (SOFT=1)
def _segment_loaded(soft, integrator, kt=0.0, mu=0.0, Pt=0.0, free_x=False, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(5, 0.0, 0.0, 0.0)
    ops.mass(5, M, M, M)
    ops.fix(5, 0 if free_x else 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    extra = ("-outward", 0.0, 0.0, 1.0) + (("-soft", soft) if soft > 0.0 else ())
    ops.contact(1, 1, 2, KN_STIFF, kt, mu, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, Pt, 0.0, -1.0e3)
    v = 0.0 if free_x else _kick(integrator, -0.5)
    if v != 0.0:
        ops.setNodeVel(5, 3, v, "-commit")
    _explicit(integrator)
    return _run(5, 1 if free_x else 3, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_soft_node_to_segment_bounded(integrator):
    assert _segment_loaded(0.0, integrator)[0]
    div, maxabs = _segment_loaded(0.1, integrator)
    assert not div and maxabs < 0.05, maxabs


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_soft_node_to_segment_stick_bounded(integrator):
    """Continuous stick (cap mu N = 5e3 > 50) with a stiff kt = 1e8: without softKt the soft normal
    alone leaves a stick mode with omega dt = 10."""
    assert _segment_loaded(0.0, integrator, kt=KN_STIFF, mu=5.0, Pt=50.0, free_x=True)[0]
    div, maxabs = _segment_loaded(0.1, integrator, kt=KN_STIFF, mu=5.0, Pt=50.0, free_x=True)
    assert not div and maxabs < 0.05, maxabs


def _segment_2d(soft, integrator, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x in ((1, -1.0), (2, 1.0)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    ops.node(5, 0.0, 0.0)
    ops.mass(5, M, M)
    ops.fix(5, 1, 0)
    ops.contactSurface(1, "-master", 2, 1, 2)
    ops.contactSurface(2, "-slave", 5)
    extra = ("-outward", 0.0, 1.0) + (("-soft", soft) if soft > 0.0 else ())
    ops.contact(1, 1, 2, KN_STIFF, 0.0, 0.0, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, -1.0e3)
    v = _kick(integrator, -0.5)
    if v != 0.0:
        ops.setNodeVel(5, 2, v, "-commit")
    _explicit(integrator)
    return _run(5, 2, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_soft_2d_node_to_segment_bounded(integrator):
    assert _segment_2d(0.0, integrator)[0]
    div, maxabs = _segment_2d(0.1, integrator)
    assert not div and maxabs < 0.05, maxabs


# --------------------------------------------------------------------- mortar (SOFT=2)
def _master_quad():
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def _matched_slave(z, fix_x=1):
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), z)
        ops.mass(t, M, M, M)
        ops.fix(t, fix_x, 1, 0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _mortar_loaded(soft, integrator, nsteps=1000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    _matched_slave(0.0)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contact(1, 1, 2, "-mortar", "-epsN", KN_STIFF, *extra, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    v = _kick(integrator, -0.5)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -25.0)
        if v != 0.0:
            ops.setNodeVel(t, 3, v, "-commit")
    _explicit(integrator)
    return _run(5, 3, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_soft2_mortar_bounded(integrator):
    # 1000 steps: the four slave nodes are free masses joined only through the facet; over longer
    # runs their bounces accumulate a rocking that is not the instability tested here (up/64)
    assert _mortar_loaded(0.0, integrator)[0]
    div, maxabs = _mortar_loaded(0.1, integrator)
    assert not div and maxabs < 0.05, maxabs


def _soft2_stick(integrator, epsT, mu=2.0, soft=0.1, P=200.0, nsteps=2000):
    """SOFT=2 friction in stick at a high mu with a tangential push below the cone. Returns
    (diverged, max |x|)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    _matched_slave(0.0, fix_x=0)
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, "-soft", soft, "-mu", mu, "-epsT", epsT,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 10.0, 0.0, -P)          # tangential 10 < mu P = 400: stick
    _explicit(integrator)
    return _run(5, 1, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_soft2_friction_stick_bounded(integrator):
    """A stiff -epsT 1e9 (omega_t dt ~ 30) stays bounded: the SOFT=2 stick penalty is sized too."""
    div, maxabs = _soft2_stick(integrator, epsT=1.0e9)
    assert not div and maxabs < 1.0e-3, maxabs


def _soft2_impact(integrator, muc, soft=0.1, v_mag=2.0, gap0=0.02, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    _matched_slave(gap0)
    visc = ("-visc", muc) if muc > 0.0 else ()
    ops.contact(1, 1, 2, "-mortar", "-epsN", KN_STIFF, "-soft", soft, *visc,
                "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 3, -v_mag, "-commit")
    _explicit(integrator)
    v_out = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            return None
        z = gap0 + ops.nodeDisp(5, 3)
        if not math.isfinite(z) or abs(z) > 1.0e3:
            return None
        if z > 0.5 * gap0 and ops.nodeVel(5, 3) > 0.0:
            v_out = ops.nodeVel(5, 3)
    return abs(v_out / v_mag)


@pytest.mark.parametrize("integrator", ["NewmarkExplicit", "ExplicitBathe"])
def test_soft2_visc_dissipates(integrator):
    """-visc on the SOFT=2 path: the soft impact rebounds near-conservatively at a stiff epsN, and
    the dashpot lowers the restitution monotonically."""
    e0 = _soft2_impact(integrator, 0.0)
    e1 = _soft2_impact(integrator, 200.0)
    e2 = _soft2_impact(integrator, 500.0)
    assert e0 is not None and e1 is not None and e2 is not None, (e0, e1, e2)
    assert abs(1.0 - e0) < 0.1, e0           # the stiff epsN injects energy (e >> 1) here
    assert e1 < e0 - 0.1, (e0, e1)
    assert e2 < e1 - 0.1, (e1, e2)


# ----------------------------------------------------------------- edge-edge (-edgeSoft)
EPS_EDGE = 1.0e8
ZC = 0.02
V0 = -2.0


def _bars(zc):
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, zc), (6, 1.0, 0.0, zc),
                       (7, 1.0, 0.0, zc + 0.3), (8, -1.0, 0.0, zc + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _edge_loaded(soft, integrator, nsteps=2000):
    """The slave bottom edge pressed onto the crossing master edge by a sustained load."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _bars(0.0)
    v = _kick(integrator, V0)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
        ops.mass(t, 1.0, 1.0, 1.0)
        if v != 0.0:
            ops.setNodeVel(t, 3, v, "-commit")
    args = [1, 1, 2, "-mortar", "-epsN", EPS_EDGE, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if soft:
        args += ["-edgeSoft", 0.1]
    ops.contact(*args)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, -100.0)
    _explicit(integrator)
    return _run(5, 3, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT)
def test_edge_soft_bounded(integrator):
    assert _edge_loaded(False, integrator)[0]
    div, maxabs = _edge_loaded(True, integrator)
    assert not div and maxabs < 0.05, maxabs


def _edge_stick(integrator, nsteps=2000):
    """The slave edge held on the master edge in stick (cap mu N = 1e3 > the 20 drag along y) with a
    stiff edge stick penalty -edgeKt 1e8 (omega_t dt ~ 7 on the m_eff = 2 edge mode)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _bars(0.0)
    for t in (5, 6):
        ops.fix(t, 1, 0, 0)
        ops.mass(t, 1.0, 1.0, 1.0)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS_EDGE, "-edgeedge", "-edgeBand", 0.15,
                "-edgeMu", 5.0, "-edgeKt", 1.0e8, "-edgeSoft", 0.1, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 10.0, -100.0)
    _explicit(integrator)
    return _run(5, 2, nsteps)


@pytest.mark.parametrize("integrator", EXPLICIT_ALL)
def test_edge_soft_stick_bounded(integrator):
    div, maxabs = _edge_stick(integrator)
    assert not div and maxabs < 1.0e-3, maxabs


# ------------------------------------------- statics and implicit transient: bitwise inert
def _plane_implicit(soft, transient):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.node(90, 0.0, 0.0, 0.0)
    ops.fix(90, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    ops.element("zeroLength", 90, 90, 1, "-mat", 90, "-dir", 3)
    ops.contactSurface(1, "-slave", 1)
    extra = (("-soft", soft) if soft > 0.0 else ()) + ("-visc", 5.0)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e6, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    if transient:
        _implicit_newmark()
    else:
        _static()
    return ops.nodeDisp(1, 3), ops.nodeVel(1, 3)


def _segment_implicit(soft, transient):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(5, 0.0, 0.0, 0.0)
    ops.mass(5, M, M, M)
    ops.fix(5, 0, 1, 0)
    ops.node(90, 0.0, 0.0, 0.0)
    ops.fix(90, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    ops.element("zeroLength", 90, 90, 5, "-mat", 90, 90, "-dir", 1, 3)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    extra = ("-outward", 0.0, 0.0, 1.0) + (("-soft", soft) if soft > 0.0 else ())
    ops.contact(1, 1, 2, 1.0e6, 1.0e6, 0.5, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 10.0, 0.0, -1.0e3)
    if transient:
        _implicit_newmark()
    else:
        _static()
    return ops.nodeDisp(5, 1), ops.nodeDisp(5, 3)


def _mortar_implicit(soft, transient):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    _matched_slave(0.0, fix_x=0)
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(90 + t, float(x), float(y), 0.0)
        ops.fix(90 + t, 1, 1, 1)
        ops.element("zeroLength", 900 + t, 90 + t, t, "-mat", 90, 90, "-dir", 1, 3)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, *extra, "-mu", 0.5, "-epsT", 1.0e6,
                "-visc", 5.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 2.0, 0.0, -25.0)
    if transient:
        _implicit_newmark()
    else:
        _static()
    return tuple(ops.nodeDisp(t, d) for t in (5, 6, 7, 8) for d in (1, 3))


def _edge_implicit(soft, transient):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _bars(-1.0e-4)
    for t in (5, 6):
        ops.fix(t, 1, 0, 0)
        ops.mass(t, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 500.0)
    for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
        ops.node(gnd, x, 0.0, -1.0e-4)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, 1, "-dir", 2, 3)
    args = [1, 1, 2, "-mortar", "-epsN", 1.0e6, "-edgeedge", "-edgeBand", 0.15,
            "-edgeMu", 0.5, "-edgeKt", 1.0e4, "-outward", 0.0, 0.0, 1.0]
    if soft:
        args += ["-edgeSoft", 0.1]
    ops.contact(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 5.0, -200.0)
    if transient:
        _implicit_newmark()
    else:
        _static()
    return tuple(ops.nodeDisp(t, d) for t in (5, 6) for d in (2, 3))


@pytest.mark.parametrize("transient", [False, True], ids=["static", "newmark"])
@pytest.mark.parametrize("case", [_plane_implicit, _segment_implicit, _mortar_implicit,
                                  _edge_implicit],
                         ids=["plane", "segment", "mortar", "edge"])
def test_soft_bitwise_inert_when_stiffness_is_assembled(case, transient):
    """Static and implicit Newmark integrators assemble the contact stiffness: -soft / -edgeSoft
    leave every result bitwise identical (the probe of the explicit detection leaves the tangent
    untouched)."""
    off = case(0.0, transient)
    on = case(0.1, transient)
    assert on == off, (on, off)
    assert all(math.isfinite(x) for x in off)
