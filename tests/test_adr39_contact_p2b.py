"""ADR-39 P2b-1 — faceted node-to-segment (NTS) penalty contact vs a single FIXED
master segment (frictionless). The narrow-phase kernel (projection + derived
outward normal + penalty + B-operator + kn*BᵀB tangent) validated oracle-first in
contact_prototypes/proto_p2b_nts.py (7/7); here it runs in compiled OpenSees.

P2b-1 scope (per the design gate): one adapter per (slave node, master segment),
master FIXED, `-kn $val`, `-outward` orientation. Deformable master + auto-kn +
Hertz + the FD-on-rotated ∂n/∂u tangent gate = P2b-2.

Setup: a master surface = ONE quad-4 (or tri-3) segment, all corner nodes FIXED;
slave node(s) in a SLAVE_NODES surface; `contact ... -outward ox oy oz` orients the
derived normal toward the slave's allowed half-space (robust to a just-penetrated
start, mirroring the P2a static gate).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6


def _quad_master(z=0.0, L=1.0, base=101, reverse=False):
    """4 FIXED corner nodes of a flat quad in the z=z plane (CCW). Returns the node
    tag list (reversed winding if requested)."""
    h = L / 2.0
    pts = [(-h, -h), (h, -h), (h, h), (-h, h)]
    tags = []
    for i, (x, y) in enumerate(pts):
        t = base + i
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
        tags.append(t)
    return tags[::-1] if reverse else tags


def test_p2b_static_penetration_quad():
    """Slave just-penetrated below a fixed quad; free z, loaded -z with P. NTS
    penalty static equilibrium penetration = P/kn (rel 1e-6) — the analytic gate."""
    P = 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _quad_master(z=0.0, L=1.0)
    ops.node(1, 0.0, 0.0, -1.0e-8)            # just penetrated (gap0 < 0)
    ops.fix(1, 1, 1, 0)                         # free z only (normal is +z)
    ops.contactSurface(10, "-master", 4, *m)   # one quad segment
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    pen = -(-1.0e-8 + ops.nodeDisp(1, 3))      # penetration depth (pos)
    expect = P / KN
    assert pen > 0 and abs(pen - expect) / expect < 1e-6, f"pen={pen} expect={expect}"


def test_p2b_static_penetration_tri():
    """Same analytic gate on a tri-3 master segment (linear shape, closed-ish proj)."""
    P = 5.0e2
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((101, 102, 103), ((-0.5, -0.5), (0.5, -0.5), (0.0, 0.6))):
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)                         # free z only (normal is +z)
    ops.contactSurface(10, "-master", 3, 101, 102, 103)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    pen = -(-1.0e-8 + ops.nodeDisp(1, 3))
    expect = P / KN
    assert pen > 0 and abs(pen - expect) / expect < 1e-6, f"tri pen={pen} expect={expect}"


def test_p2b_oblique_master():
    """Tilted master quad (normal n=(0,−sin30,cos30)); slave free in Z only, loaded
    −z. The derived normal must follow the tilt: the contact stiffness projects onto
    the free z-DOF as kn·nz², so z-penetration = P/(kn·nz²) (off-axis projection +
    B-operator). Free-z keeps the rank-1 contact non-singular (the inclined gate)."""
    P = 1.0e3
    th = math.radians(30.0)
    c, s = math.cos(th), math.sin(th)
    n = (0.0, -s, c)                          # outward normal of the x-tilted quad
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    h = 0.5
    mtags = []
    for i, (x, y0) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        # rotate (x, y0, 0) about the x-axis by th -> normal tilts into y-z
        t = 101 + i
        ops.node(t, x, y0 * c, y0 * s)
        ops.fix(t, 1, 1, 1)
        mtags.append(t)
    ops.node(1, 0.0, 0.0, -1.0e-8)            # just below the plane at the origin
    ops.fix(1, 1, 1, 0)                         # free z only
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", n[0], n[1], n[2])
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    pen_z = -(-1.0e-8 + ops.nodeDisp(1, 3))    # z-penetration depth (pos)
    expect = P / (KN * n[2] * n[2])
    assert pen_z > 0 and abs(pen_z - expect) / expect < 1e-5, f"oblique pen_z={pen_z} expect={expect}"


def _static_pen(master_tags, outward=(0.0, 0.0, 1.0), P=1.0e3):
    """helper: add the slave + contact onto an already-built master, run static,
    return the equilibrium penetration depth. (Master nodes created by the caller.)"""
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)                         # free z only (normal is +z)
    ops.contactSurface(10, "-master", 4, *master_tags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", *outward)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return -(-1.0e-8 + ops.nodeDisp(1, 3))


def test_p2b_winding_flip_identical():
    """GATE BLOCKER-1: reversing the master node winding must give an IDENTICAL
    contact response — the outward normal is DERIVED + oriented from -outward, never
    trusted from the node order. (A naive winding normal would flip → pass-through.)"""
    # forward winding
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _quad_master(z=0.0, L=1.0, base=101, reverse=False)
    pen_fwd = _static_pen(m)
    # reversed winding (same geometry, opposite node order)
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    mr = _quad_master(z=0.0, L=1.0, base=101, reverse=True)
    pen_rev = _static_pen(mr)
    assert abs(pen_fwd - pen_rev) / pen_fwd < 1e-9, f"winding changed pen: {pen_fwd} vs {pen_rev}"


def test_p2b_out_of_bounds_passes_through():
    """A slave whose projection lies OUTSIDE the (small) master segment must feel NO
    contact (bounded projection + sentinel — never a garbage force). Explicit: the
    laterally-offset slave is flung at the plane level and passes through (no rebound)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _quad_master(z=0.0, L=0.2, base=101)   # small quad near origin
    ops.node(1, 5.0, 5.0, 0.01)                # far outside the quad, just above
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.setNodeVel(1, 3, -2.0, "-commit")      # flung straight down
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    for _ in range(2000):
        assert ops.analyze(1, dt) == 0
    z = 0.01 + ops.nodeDisp(1, 3)
    assert z < -0.05, f"slave should pass through (no contact OOB), z={z}"


def test_p2b_explicit_impact_restitution():
    """Slave mass above the quad center, flung down, elastic rebound e≈1, bounded
    penetration (= v*sqrt(m/kn) = 0.002). Element-free model (no anchor) — exercises
    the B1 teardown fix too."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _quad_master(z=0.0, L=1.0, base=101)
    gap0 = 0.01
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.setNodeVel(1, 3, -2.0, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    max_pen, v_out, blew = 0.0, 0.0, False
    for _ in range(4000):
        if ops.analyze(1, dt) != 0:
            blew = True
            break
        z = gap0 + ops.nodeDisp(1, 3)
        vz = ops.nodeVel(1, 3)
        if z < 0:
            max_pen = max(max_pen, -z)
        if z > 0.5 * gap0 and vz > 0:
            v_out = vz
    e = v_out / 2.0
    assert (not blew) and e > 0.9 and e <= 1.05, f"blew={blew} e={e}"
    assert abs(max_pen - 0.002) / 0.002 < 0.05, f"pen={max_pen}"


def test_p2b_asdimplex_cross_check():
    """ORACLE (design gate): the NTS segment contact and a hand-placed
    ZeroLengthContactASDimplex node-pair must give the SAME static penetration for
    the same kn + normal (two independent contact implementations, rel 1e-6)."""
    P = 1.0e3
    # --- NTS segment contact ---
    pen_nts = None
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _quad_master(z=0.0, L=1.0, base=101)
    pen_nts = _static_pen(m, P=P)
    # --- ZeroLengthContactASDimplex node-pair (master node fixed, slave free z) ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(200, 0.0, 0.0, 0.0); ops.fix(200, 1, 1, 1)       # master
    ops.node(1, 0.0, 0.0, -1.0e-8); ops.fix(1, 1, 1, 0)        # slave, free z
    # element zeroLengthContactASDimplex tag mNd sNd Kn Kt mu -orient nx ny nz
    ops.element("zeroLengthContactASDimplex", 1, 200, 1, KN, 0.0, 0.0,
                "-orient", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    pen_zlc = -(-1.0e-8 + ops.nodeDisp(1, 3))
    assert pen_nts > 0 and pen_zlc > 0, (pen_nts, pen_zlc)
    assert abs(pen_nts - pen_zlc) / pen_zlc < 1e-6, f"NTS {pen_nts} vs ZLC {pen_zlc}"
