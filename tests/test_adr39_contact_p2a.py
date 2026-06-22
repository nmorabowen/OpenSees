"""ADR-39 P2a — NTS penalty vs a RIGID analytical plane (frictionless).

P2a is the de-risked first mechanics rung (per the design gate): penalty force +
B-operator (B = nᵀ) + active set, NO projection kernel / auto-kn / master
compliance. One adapter per slave node, connectivity = {slave}.

Gates (rel tol per the fork's ~1e-7 precedent — penalty has an O(load/kn) baseline
so absolute 1e-12 is wrong):
  - axis-aligned penetration g = P/kn, IMPLICIT static (LoadControl+Newton), 1e-6
  - axis-aligned EXPLICIT impact — restitution e≈1 (elastic), bounded penetration
  - inclined plane EXPLICIT impact (off-axis n residual) — stable + rebounds

Each model carries a NEGLIGIBLE z-direction anchor truss (tiny EA, no z-motion in
any case → zero effect on the contact physics): a contacting body is connected to
"ground", and this also gives the model an element-backed FE_Element so the base
FE_Element shared scratch is allocated (a contact adapter is not element-backed).

Sign (PenaltySP_FE convention): getResidual drives the slave toward g=0 → a
penetrating node is pushed back along +n (restoring, never attracting).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6


def _anchor_truss(slave_tag, X):
    """Tiny z-direction truss from the slave to a fixed anchor (zero physical
    effect: no test has z-motion). Anchor node tag 9990+, truss tag 9990+."""
    a = 9990 + slave_tag
    ops.node(a, X[0], X[1], X[2] + 1.0)
    ops.fix(a, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 9999, 1.0)        # EA tiny (A below)
    ops.element("Truss", a, slave_tag, a, 1.0e-9, 9999)


def test_p2a_axis_aligned_penetration_static():
    """Plane x=0, normal +x; node free in x, fixed y,z; load -x with P.
    IMPLICIT static equilibrium penetration |x| = P/kn."""
    P = 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, -1.0e-8, 0.0, 0.0)            # start just-penetrated → contact active
    ops.fix(1, 0, 1, 1)                        # free x only
    _anchor_truss(1, [-1.0e-8, 0.0, 0.0])
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, KN)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, -P, 0.0, 0.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    pen = -(-1.0e-8 + ops.nodeDisp(1, 1))     # penetration depth (pos<0)
    expect = P / KN
    assert pen > 0 and abs(pen - expect) / expect < 1e-6, f"pen={pen} expect={expect}"


def _explicit_impact(normal, v_mag=2.0, m=1.0, gap0=0.01, nsteps=4000):
    """Mass node starts gap0 outside, flung toward the plane along -n. Returns
    (blew_up, restitution, max_normal_penetration)."""
    n = list(normal)
    nl = math.sqrt(sum(c * c for c in n)); n = [c / nl for c in n]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    X0 = [gap0 * c for c in n]
    ops.node(1, X0[0], X0[1], X0[2])
    ops.mass(1, m, m, m)
    _anchor_truss(1, X0)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, n[0], n[1], n[2], 0.0, 0.0, 0.0, KN)
    for d in range(3):
        ops.setNodeVel(1, d + 1, -v_mag * n[d], "-commit")
    dt = 0.5 * 2.0 * math.sqrt(m / KN)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    max_pen, v_out_n, blew = 0.0, 0.0, False
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            blew = True
            break
        u = [ops.nodeDisp(1, d + 1) for d in range(3)]
        v = [ops.nodeVel(1, d + 1) for d in range(3)]
        gap = sum(n[d] * (X0[d] + u[d]) for d in range(3))
        vn = sum(n[d] * v[d] for d in range(3))
        if gap < 0:
            max_pen = max(max_pen, -gap)
        if gap > 0.5 * gap0 and vn > 0:
            v_out_n = vn
    return blew, abs(v_out_n / v_mag), max_pen


def test_p2a_axis_aligned_impact_restitution():
    """Axis-aligned explicit impact: stable, elastic (e≈1), bounded penetration
    (= v·sqrt(m/kn) = 0.002 here)."""
    blew, e, max_pen = _explicit_impact([1, 0, 0])
    assert (not blew) and e <= 1.05 and max_pen < 0.02, f"blew={blew} e={e} pen={max_pen}"
    assert e > 0.9 and abs(max_pen - 0.002) / 0.002 < 0.05, f"e={e} pen={max_pen}"


def test_p2a_inclined_impact_offaxis():
    """Off-axis normal (3-4-5): explicit impact stable, rebounds along +n
    (active set on→off = release built in)."""
    blew, e, max_pen = _explicit_impact([0.6, 0.8, 0.0])
    assert (not blew) and e <= 1.05 and max_pen < 0.02, f"blew={blew} e={e} pen={max_pen}"
    assert e > 0.9, f"off-axis rebound e={e}"
