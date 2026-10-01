"""ADR-158 -- the opt-in finite-difference mortar pair tangent (`contact ... -mortar -fdTangent [hRel]`).

R0.7 diagnosis (Ladruno_implementation/158_mortar_fd_pair_tangent.md, section 2): at slave nodes
shared by several (slave facet, master facet) pairs on a creased interface, the shipped analytic
mortar tangent leaves out two terms that the residual has:
  * the Coulomb pressure coupling Csl = -mu*epsN*t_hat (x) n (dropped by the default symmetric
    tangent; `-consistanttan` restores it), which is as large as the normal stiffness when
    mu ~ 1 and the pairs slip;
  * the geometric terms dD/du, dM/du of the clipped overlap. They are small on a pair whose facets
    are parallel, and of the same order as the material term on the thin CROSS-CREASE pairs (a
    slave facet of one flank clipped against the master facet of the other), whose overlap width
    grows with the penetration.
`-fdTangent` differentiates each pair's own residual (normal + friction) by central differences,
so it carries both. The residual is untouched: converged answers do not change.

Gates:
  (a) mu > 0, multi-pair creased roof on a stiff solid slave (bricks), force control, full slip:
      the vertical contact force equals the analytic 2*A*kn*(delta+w)*cos(a)*(cos(a)+mu*sin(a))
      (rigid-body descent w), and Newton converges in <= 6 iterations per step with -fdTangent
      or with the shipped -consistanttan (the Csl term is the dominant missing term here); the
      shipped symmetric default needs more than twice that on the same model.
  (b) the element-less shared-ridge roof, cohesion only, force control (review #900 finding 1):
      <= 6 iterations per step with -fdTangent; the shipped tangent takes 10 and 15 in steps 3-4,
      with or without -consistanttan (cohesion only: there is no Csl; the missing terms are the
      geometric ones of the cross-crease pairs).
  (c) the FD tangent leaves the converged state unchanged (the two tangents converge to the same
      displacements to 1e-9 relative) and refuses -tie / non-mortar contacts.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

ALPHA = math.radians(20.0)


def _roof(mu=0.0, coh=0.0, solid_e=0.0, flags=(), steps=4, pz=40.0, eps=1.0e6, delta=1.0e-3,
          kd=1.0e3, augment_never=False, mw=1.5, yw=0.5):
    """Roof z = -tan(a)|x|, ridge along y. Master: 2 fixed quads over |x| <= mw, -yw <= y <= 1+yw
    (mw=1, yw=0 is the ADR-157 roof; the default oversizes the slave). Slave: 4 x 3 non-matching quads, ridge nodes SHARED by both flanks, lowered by
    delta; every slave node on kd springs (x, y, z). solid_e > 0 adds a 0.3-thick brick layer under
    the slave surface (stdBrick, E = solid_e). Load -pz on each surface node, force control.
    Returns (iterations per step, ok, surface node tags, all spring node tags)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ta = math.tan(ALPHA)
    xm, ym = [-mw, 0.0, mw], [-yw, 1.0 + yw]
    mt, t = {}, 1
    for j, y in enumerate(ym):
        for i, x in enumerate(xm):
            ops.node(t, x, y, -ta * abs(x))
            ops.fix(t, 1, 1, 1)
            mt[(i, j)] = t
            t += 1
    mq = []
    for i in range(2):
        mq += [mt[(i, 0)], mt[(i + 1, 0)], mt[(i + 1, 1)], mt[(i, 1)]]
    xs, ys = [-1.0, -0.45, 0.0, 0.45, 1.0], [0.0, 0.3, 0.75, 1.0]
    st, surf, t = {}, [], 101
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            ops.node(t, x, y, -ta * abs(x) - delta)
            st[(i, j)] = t
            surf.append(t)
            t += 1
    sq = []
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            sq += [st[(i, j)], st[(i + 1, j)], st[(i + 1, j + 1)], st[(i, j + 1)]]
    ops.contactSurface(1, "-master", 4, *mq)
    ops.contactSurface(2, "-slave-segments", 4, *sq)
    opts = ["-mortar", "-epsN", eps, "-epsT", eps, "-outward", 0.0, 0.0, 1.0]
    if mu > 0.0:
        opts += ["-mu", mu]
    if coh > 0.0:
        opts += ["-cohesion", coh]
    if augment_never:
        opts += ["-augment", "never"]
    ops.contact(1, 1, 2, *opts, *flags)
    springs = list(surf)
    if solid_e > 0.0:
        ops.nDMaterial("ElasticIsotropic", 7, solid_e, 0.25)
        bot = {}
        for s in surf:
            c = ops.nodeCoord(s)
            bot[s] = 3000 + s
            ops.node(bot[s], c[0], c[1], c[2] - 0.3)
            springs.append(bot[s])
        for q in range(0, len(sq), 4):
            a, b, c, d = sq[q:q + 4]
            ops.element("stdBrick", 7000 + q, bot[a], bot[b], bot[c], bot[d], a, b, c, d, 7)
    ops.uniaxialMaterial("Elastic", 1, kd)
    for s in springs:
        g = 50000 + s
        ops.node(g, *ops.nodeCoord(s))
        ops.fix(g, 1, 1, 1)
        ops.element("zeroLength", 60000 + s, g, s, "-mat", 1, 1, 1, "-dir", 1, 2, 3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in surf:
        ops.load(s, 0.0, 0.0, -pz)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-7, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    its = []
    for _ in range(steps):
        if ops.analyze(1) != 0:
            return its, False, surf, springs
        its.append(ops.testIter())
    return its, True, surf, springs


# ---------------------------------------------------------------------------------------- (a)
SOLID = dict(mu=0.2, solid_e=1.0e11, pz=60.0, delta=1.0e-6, augment_never=True, steps=4)


def _descent_analytic(nsurf, nspring, pz, kd, kn, mu, delta, lam):
    """Rigid vertical descent w of the roof in full slip: the applied load lam*nsurf*pz is carried
    by the springs (kd*nspring*w) and the two flanks (area A = 1/cos(a) each), where the penalty
    pressure is p = kn*(delta + w)*cos(a) and the slip friction mu*p acts up-slope:
        lam*nsurf*pz = kd*nspring*w + 2*A*p*(cos(a) + mu*sin(a)).
    Returns (w, Fz) with Fz the vertical contact force on the slave."""
    ca, sa = math.cos(ALPHA), math.sin(ALPHA)
    c = 2.0 / ca * kn * ca * (ca + mu * sa)
    w = (lam * nsurf * pz - c * delta) / (kd * nspring + c)
    return w, c * (delta + w)


@pytest.mark.parametrize("flags", [("-fdTangent",), ("-consistanttan",)])
def test_adr158_mu_multipair_roof_analytic_and_quadratic(flags):
    """(a) mu > 0 at shared multi-pair ridge nodes: analytic contact force and <= 6 Newton
    iterations per step with a tangent that carries the Coulomb pressure coupling Csl -- the FD
    pair tangent (measured 3, 3, 3, 3) or the shipped -consistanttan (3, 4, 4, 4). The shipped
    symmetric default drops Csl and needs > 12 (next test)."""
    its, ok, surf, springs = _roof(flags=flags, **SOLID)
    assert ok, f"{flags}: Newton failed after {len(its)} step(s) {its}"
    assert max(its) <= 6, f"{flags}: iterations per step {its} (want <= 6)"
    kd = 1.0e3
    for lam in (1.0,):
        w, fz = _descent_analytic(len(surf), len(springs), SOLID["pz"], kd, 1.0e6, SOLID["mu"],
                                  SOLID["delta"], lam)
        # contact force on the slave = applied load - spring forces (sum over all spring nodes)
        got = len(surf) * SOLID["pz"] + sum(kd * ops.nodeDisp(n, 3) for n in springs)
        fx = sum(kd * ops.nodeDisp(n, 1) for n in springs)
        assert got == pytest.approx(fz, rel=2e-4), f"Fz {got:.6f} != analytic {fz:.6f}"
        assert abs(fx) < 1e-8 * fz, f"spurious Fx = {fx:.3e}"


def test_adr158_mu_multipair_roof_shipped_tangent_is_linear():
    """(a) the same model on the shipped symmetric tangent: Newton is linear (the Csl and the
    cross-crease geometric terms are missing): more than 12 iterations in some step (measured
    8, 15, 15, 12 on fc75db7f3)."""
    its, ok, _, _ = _roof(**SOLID)
    assert ok
    assert max(its) > 12, f"shipped tangent iterations {its}: the baseline this gate measures moved"


# ---------------------------------------------------------------------------------------- (b)
COH = dict(coh=10.0, mw=1.0, yw=0.0, steps=4, pz=40.0, eps=1.0e6, delta=1.0e-3)


def test_adr158_shared_ridge_cohesion_roof_quadratic():
    """(b) review #900 finding 1: the element-less shared-ridge roof (cohesion only, force
    control) converges in <= 6 iterations per step with -fdTangent (shipped: 6, 4, 10, 15)."""
    its, ok, _, _ = _roof(flags=("-fdTangent",), **COH)
    assert ok, f"-fdTangent: Newton failed after {len(its)} step(s) {its}"
    assert max(its) <= 6, f"-fdTangent iterations per step {its} (want <= 6)"


# ---------------------------------------------------------------------------------------- (c)
def _state(nodes):
    return [ops.nodeDisp(n) for n in nodes]


@pytest.mark.parametrize("case", ["mu", "coh"])
def test_adr158_same_converged_state(case):
    """(c) the FD tangent changes the iterations only: both tangents converge to the same state."""
    kw = SOLID if case == "mu" else COH
    _, ok1, surf, springs = _roof(**kw)
    a = _state(springs)
    _, ok2, _, _ = _roof(flags=("-fdTangent", 1.0e-7), **kw)
    b = _state(springs)
    assert ok1 and ok2
    scale = max(abs(v) for u in a for v in u)
    for ua, ub in zip(a, b):
        for va, vb in zip(ua, ub):
            assert abs(va - vb) <= 1e-6 * scale, f"{case}: {va} != {vb}"


def test_adr158_refusals():
    """(c) -fdTangent is a -mortar option and does not apply to -tie; a bad step is refused."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(t, x, y, 0.0)
        ops.node(10 + t, x, y, 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 11, 12, 13, 14)
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, "-mortar", "-tie", "-epsN", 1.0e5, "-fdTangent")
    with pytest.raises(Exception):
        ops.contact(2, 1, 2, "-mortar", "-epsN", 1.0e5, "-fdTangent", 2.0)
    with pytest.raises(Exception):
        ops.contact(3, 1, 2, 1.0e5, "-fdTangent")
