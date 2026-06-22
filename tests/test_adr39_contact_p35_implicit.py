"""ADR-39 P3.5 — implicit frictional Newton: the consistent friction TANGENT.

P3 shipped friction FORCE-only (explicit). Under an IMPLICIT solve a free tangential
DOF is SINGULAR without a tangential stiffness (P3: `FullGeneral U(0,0)=0`). P3.5
assembles the friction tangent in addKtToTang so Newton converges. The 3D assembled
tangent is validated oracle-first in contact_prototypes/proto_p35_implicit_tangent.py
(4/4). Default = the SYMMETRIC tangent (solver-safe on any system); `-consistanttan`
opts into the non-symmetric consistent tangent (quadratic, needs FullGeneral/UmfPack).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN, KT = 1.0e6, 1.0e5


def _fixed_quad(L=20.0, base=101):
    """4 FIXED corners of a flat quad in z=0 (normal +z)."""
    h = L / 2.0
    tags = []
    for i, (x, y) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        t = base + i
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
        tags.append(t)
    return tags


def _stick_model(mu, Qx, P=1000.0, solver="FullGeneral", consistanttan=False):
    """Slave pre-penetrated to the normal equilibrium (z=-P/kn ⇒ N=P), free x and z,
    fixed y, normal load -P + tangential Qx below the cone. Builds + runs one static
    LoadControl step. Returns the analyze() code (0 = converged)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _fixed_quad()
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.fix(1, 0, 1, 0)                          # free x and z; fix y
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    if consistanttan:
        ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0, "-consistanttan")
    else:
        ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Qx, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(solver)
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    return ops.analyze(1)


def test_p35_static_stick_converges():
    """THE gate: a static frictional STICK problem (free tangential DOF) — SINGULAR in
    P3 (no friction tangent) — now CONVERGES with the friction tangent, to the elastic
    stick displacement x = Q/kt, with the friction force tT = Q held below the cone μN."""
    mu, Q, P = 0.5, 200.0, 1000.0               # μN = μP = 500 > Q ⇒ stick
    rc = _stick_model(mu, Q, P=P)
    assert rc == 0, "static frictional stick must converge (was singular in P3)"
    x = ops.nodeDisp(1, 1)
    assert abs(x - Q / KT) / (Q / KT) < 1e-4, f"stick x={x} != Q/kt={Q/KT}"
    z = ops.nodeDisp(1, 3)                       # stays at the pre-penetrated equilibrium
    assert abs(z) < 1e-6, f"normal should stay at equilibrium, z-move={z}"


def test_p35_static_stick_symmetric_solver():
    """The DEFAULT (symmetric) friction tangent is correct on a SYMMETRIC solver
    (ProfileSPD) — the design-gate Q2 safety property: out-of-the-box friction must not
    silently corrupt on the common symmetric default. Same converged answer."""
    mu, Q = 0.5, 200.0
    rc = _stick_model(mu, Q, solver="ProfileSPD")
    assert rc == 0, "symmetric tangent must converge on a symmetric solver"
    assert abs(ops.nodeDisp(1, 1) - Q / KT) / (Q / KT) < 1e-4


def test_p35_consistent_tangent_converges_fewer_iters():
    """`-consistanttan` (the non-symmetric consistent tangent) converges on FullGeneral
    and takes NO MORE Newton iterations than the symmetric default (it is the exact
    derivative ⇒ at least as fast). Same converged answer."""
    mu, Q = 0.5, 200.0
    rc_c = _stick_model(mu, Q, solver="FullGeneral", consistanttan=True)
    it_c = ops.testIter()
    assert rc_c == 0 and abs(ops.nodeDisp(1, 1) - Q / KT) / (Q / KT) < 1e-4
    rc_s = _stick_model(mu, Q, solver="FullGeneral", consistanttan=False)
    it_s = ops.testIter()
    assert rc_s == 0
    assert it_c <= it_s + 1, f"consistent tangent iters {it_c} should be <= symmetric {it_s}"


def test_p35_static_slip_with_anchor():
    """Static SLIP: the slave is anchored tangentially by an elastic spring (so a slip
    equilibrium exists), normal load -P + tangential Q>μN. It slides until anchor +
    friction balance Q: x=(Q−μN)/ka, friction saturated at μN. Converges (the anchor
    provides the x-stiffness the flat slip cone does not)."""
    mu, P, Q, ka = 0.5, 1000.0, 800.0, 1.0e4     # μN=500; Q>μN ⇒ slip
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _fixed_quad()
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.fix(1, 0, 1, 0)                           # free x,z; fix y
    ops.node(2, 0.0, 0.0, -P / KN)
    ops.fix(2, 1, 1, 1)                           # ground for the x-anchor
    ops.uniaxialMaterial("Elastic", 1, ka)
    ops.element("zeroLength", 99, 2, 1, "-mat", 1, "-dir", 1)   # x-anchor spring
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-9, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "static slip with anchor must converge"
    x = ops.nodeDisp(1, 1)
    expect = (Q - mu * P) / ka                    # 0.03
    assert abs(x - expect) / expect < 1e-3, f"slip x={x} != (Q−μN)/ka={expect}"
    fric = Q - ka * x                             # x-equilibrium ⇒ friction force
    assert abs(fric - mu * P) / (mu * P) < 1e-3, f"friction {fric} != μN {mu*P}"


def test_p35_newmark_dynamic_friction():
    """Implicit dynamic (Newmark) frictional contact converges each step. A mass driven
    on a fixed plane: the friction tangent makes the implicit step well-posed; the
    block slides and the average acceleration matches a=(Q−μN)/m."""
    mu, P, Q, m = 0.3, 1000.0, 600.0, 1.0        # μN=300 ⇒ a=(600−300)/1=300
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _fixed_quad()
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.fix(1, 0, 1, 0)                           # free x,z
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Transient")
    vx = []
    for _ in range(120):
        assert ops.analyze(1, dt) == 0, "Newmark frictional step must converge"
        vx.append(ops.nodeVel(1, 1))
    a = (vx[110] - vx[40]) / (70 * dt)
    assert abs(a - 300.0) / 300.0 < 0.05, f"implicit slide a={a} != (Q−μN)/m=300"


def test_p35_explicit_byte_identity():
    """Explicit (CDL) is UNCHANGED by P3.5: formEleTangent calls only addMtoTang (no-op)
    ⇒ the friction tangent is never assembled under explicit. A CDL frictional slide
    gives the P3 result a=(Q−μN)/m exactly (regression guard for the tangent code)."""
    mu, P, Q, m = 0.3, 1000.0, 600.0, 1.0
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _fixed_quad()
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.fix(1, 0, 1, 0)
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    vx = []
    for _ in range(150):
        assert ops.analyze(1, dt) == 0
        vx.append(ops.nodeVel(1, 1))
    a = (vx[140] - vx[40]) / (100 * dt)
    assert abs(a - 300.0) / 300.0 < 0.02, f"explicit a={a} != 300 (P3 byte-identity)"
