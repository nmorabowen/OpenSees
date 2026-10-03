"""ADR-158 -- mortar friction tangent at shared multi-pair nodes: diagnosis + the -consistanttan recipe.

R0.7 diagnosis (Ladruno_implementation/158_mortar_tangent_diagnosis_consistanttan.md, section 2): at
slave nodes shared by several (slave facet, master facet) pairs on a creased interface, the shipped
analytic mortar tangent leaves out two terms that the residual has:
  * the Coulomb pressure coupling Csl = -mu*epsN*t_hat (x) n (dropped by the default symmetric
    tangent; the shipped opt-in `-consistanttan` restores it), which is as large as the normal
    stiffness when mu ~ 1 and the pairs slip;
  * the geometric terms dD/du, dM/du of the clipped overlap, of the same order as the material term
    on the thin CROSS-CREASE pairs. No shipped tangent carries them (the FD pair tangent that does is
    parked as a diagnostic oracle: contact_prototypes/adr158_fd_pair_tangent_oracle.patch).
The residual is untouched: the converged answers do not depend on the tangent.

Gates (the mu > 0 multi-pair creased roof on a stiff solid slave, force control, full slip; this is
also the binary Coulomb multi-pair test review #900 finding 2 asked for):
  (a) with -consistanttan: the vertical contact force equals the analytic
      2*A*kn*(delta+w)*cos(a)*(cos(a)+mu*sin(a)) (rigid-body descent w), Fx = 0, and Newton converges
      in <= 6 iterations per step (measured 3, 4, 4, 4). On the per-node path state (d63f49750) the
      shared ridge nodes read a traction from the other flank and this fails.
  (a') the shipped symmetric default on the same model is the linear baseline the recipe removes:
      it either needs more than 12 iterations in some step (8, 15, 15, 12 on Windows) or does not
      converge at all. A converged run in <= 12 iterations means the baseline moved.
  (c) -consistanttan changes the iterations only: when the default also converges, both reach the
      same state to 1e-6 relative (the residual and its fixed point are shared; the iterates are
      not, so the results agree to the solver tolerance, not bitwise).
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


def test_adr158_mu_multipair_roof_consistanttan_analytic():
    """(a) mu > 0 at shared multi-pair ridge nodes with the -consistanttan recipe: analytic contact
    force and <= 6 Newton iterations per step (measured 3, 4, 4, 4)."""
    its, ok, surf, springs = _roof(flags=("-consistanttan",), **SOLID)
    assert ok, f"-consistanttan: Newton failed after {len(its)} step(s) {its}"
    assert max(its) <= 6, f"-consistanttan: iterations per step {its} (want <= 6)"
    kd = 1.0e3
    w, fz = _descent_analytic(len(surf), len(springs), SOLID["pz"], kd, 1.0e6, SOLID["mu"],
                              SOLID["delta"], 1.0)
    # contact force on the slave = applied load - spring forces (sum over all spring nodes)
    got = len(surf) * SOLID["pz"] + sum(kd * ops.nodeDisp(n, 3) for n in springs)
    fx = sum(kd * ops.nodeDisp(n, 1) for n in springs)
    assert got == pytest.approx(fz, rel=2e-4), f"Fz {got:.6f} != analytic {fz:.6f}"
    assert abs(fx) < 1e-8 * fz, f"spurious Fx = {fx:.3e}"


def test_adr158_mu_multipair_roof_shipped_tangent_is_linear():
    """(a') the same model on the shipped symmetric tangent (Csl and the cross-crease geometric terms
    missing): Newton is linear, more than 12 iterations in some step (8, 15, 15, 12 on Windows), or
    it does not converge (the baseline is platform-sensitive, review #902 finding 1)."""
    its, ok, _, _ = _roof(**SOLID)
    assert (not ok) or max(its) > 12, (
        f"shipped tangent converged in {its}: the linear baseline this gate measures moved")


# ---------------------------------------------------------------------------------------- (c)
def _state(nodes):
    return [ops.nodeDisp(n) for n in nodes]


def test_adr158_consistanttan_same_converged_state():
    """(c) -consistanttan changes the iterations only: when the shipped default converges too, both
    tangents reach the same state (to the solver tolerance, not bitwise)."""
    _, ok1, surf, springs = _roof(**SOLID)
    a = _state(springs)
    _, ok2, _, _ = _roof(flags=("-consistanttan",), **SOLID)
    b = _state(springs)
    assert ok2, "-consistanttan failed"
    if not ok1:
        pytest.skip("the shipped default did not converge on this platform (see (a')): nothing to compare")
    scale = max(abs(v) for u in a for v in u)
    for ua, ub in zip(a, b):
        for va, vb in zip(ua, ub):
            assert abs(va - vb) <= 1e-6 * scale, f"{va} != {vb}"
