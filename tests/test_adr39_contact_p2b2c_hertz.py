"""ADR-39 P2b-2c / capstone B3 — the consistent ∂n/∂u normal contact tangent + Hertz.

P2b shipped the NTS penalty tangent as kn·BᵀB only (EXACT for a flat/fixed master where
the normal n is constant) and DEFERRED the geometric ∂n/∂u block. B3 adds the consistent
geometric tangent K_geom = kn·gN·∂²gN/∂u² (the curvature / large-sliding term), opted in
per contact via `-geomtan`. Validated oracle-first in
contact_prototypes/proto_b3_normal_tangent.py (analytic == FD ~1e-10 on a warped quad;
symmetric; slave block ≡ 0 for a flat facet) + the C++ standalone b3_cpp_check.cpp.

Gates (capstone B3 row):
  (1) implicit Newton converges QUADRATICALLY on a curved interface where kn·BᵀB-only is
      slower — fewer Newton iterations, SAME converged root.
  (2) Hertz: a rigid sphere (the faceted MASTER, fixed, lowered by δ) pressed into an
      elastic stdBrick block whose top nodes are the SLAVES reproduces the analytic Hertz
      contact radius a = √(Rδ) to 15 %, a centre-peaked pressure, and the a ∝ √δ,
      P ∝ δ^1.5 depth scaling.
  (3) flat / fixed-master path is BYTE-IDENTICAL to the shipped tangent (∂n/∂u ≡ 0 there).
  (5) the geometric normal tangent is SYMMETRIC ⇒ solver-safe on a symmetric solver.
  (6) `ladrunoContactForce` reports exactly 0.0 for a 3D pair that has RELEASED (static
      pull-off and explicit rebound).

Gate (2) history (r6 backport, jaabell/OpenSees upstream round 6): the previous Hertz
deck pinned the slaves ON the lowered sphere and made the block the master, and checked
the GEOMETRIC radius √(2Rδ). It passed only because of the 3D stale-force readout (a pair
that engaged during the Newton iterations and then released kept its last force, see
Ladruno_implementation/ledger/quirks/LEGACY-ladrunocontactforce-keeps-a-stale-reading-on-a.md):
every slave inside the seeded interference reported a force, so the "patch" was the
interference circle. With the readout fixed that deck's patch is the compliant one and
the geometric oracle is wrong. The redesigned deck (the indenter is the master, the
deformable body's nodes are the slaves, the compliant radius √(Rδ) is the oracle) is the
one validated for the upstream port: 37-69 slaves in contact across these cases.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 2.0e6


# ---------------------------------------------------------- Hertz FE harness
def _hertz_indent(nx=24, nz=8, W=2.0, H=1.0, R=20.0, delta=0.004, E=1000.0,
                  nu=0.3, kn=1.0e4, nm=41, geom=False, span=None):
    """A rigid sphere pressed into an elastic stdBrick block. The indenter is the MASTER:
    fixed quad facets on the sphere surface, lowered by δ. The SLAVES are the top nodes
    of the block, so `ladrunoContactForce` gives the nodal contact force of the deformable
    body (the adapter traction is not in nodeReaction). Returns
    (rc, sorted [(r, pressure, force)] of the slaves in contact)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx, dz = W / nx, H / nz
    nid = {}
    t = 1
    for k in range(nz + 1):
        for j in range(nx + 1):
            for i in range(nx + 1):
                ops.node(t, -W / 2 + i * dx, -W / 2 + j * dx, -H + k * dz)
                nid[(i, j, k)] = t
                if k == 0:
                    ops.fix(t, 1, 1, 1)                  # half-space base
                elif i in (0, nx) or j in (0, nx):
                    ops.fix(t, 1, 1, 0)                  # roller side walls
                t += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    eid = 1
    for k in range(nz):
        for j in range(nx):
            for i in range(nx):
                ops.element("stdBrick", eid,
                            nid[(i, j, k)], nid[(i + 1, j, k)],
                            nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)],
                            nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)],
                            nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)], 1)
                eid += 1
    if span is None:
        span = 2.0 * math.sqrt(R * delta)
    sdx = 2 * span / (nm - 1)
    mid = {}
    m = 400000
    for jj in range(nm):
        for ii in range(nm):
            x, y = -span + ii * sdx, -span + jj * sdx
            ops.node(m, x, y, (R - math.sqrt(R * R - x * x - y * y)) - delta)
            ops.fix(m, 1, 1, 1)
            mid[(ii, jj)] = m
            m += 1
    mt = []
    for jj in range(nm - 1):
        for ii in range(nm - 1):
            mt += [mid[(ii, jj)], mid[(ii + 1, jj)],
                   mid[(ii + 1, jj + 1)], mid[(ii, jj + 1)]]
    ops.contactSurface(10, "-master", 4, *mt)
    slaves, coord = [], {}
    for j in range(nx + 1):
        for i in range(nx + 1):
            x, y = -W / 2 + i * dx, -W / 2 + j * dx
            if abs(x) < span and abs(y) < span:
                s = nid[(i, j, nz)]
                slaves.append(s)
                coord[s] = math.hypot(x, y)
    ops.contactSurface(20, "-slave", *slaves)
    args = [1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 0.0, -1.0]
    if geom:
        args.append("-geomtan")
    ops.contact(*args)
    ops.constraints("LadrunoContact")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    rc = ops.analyze(1)
    data = []
    for s in slaves:
        f = ops.ladrunoContactForce(s)
        if f > 0.0:
            data.append((coord[s], f / (dx * dx), f))
    data.sort()
    return rc, data


# --------------------------------------------------------------------- helpers
def _warp_master(w, L=4.0, base=101):
    """A warped (non-planar) quad master, FIXED: nodes 0,2 lifted +w, 1,3 dropped −w ⇒
    a genuine saddle (κ_ξη ≠ 0). w=0 ⇒ a flat quad (κ=0)."""
    h = L / 2.0
    pts = [(-h, -h, +w), (h, -h, -w), (h, h, +w), (-h, h, -w)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        t = base + i
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 1, 1, 1)
        tags.append(t)
    return tags


def _slide_model(geom, w, Q=4.0e3, P=2.0e3, ka=2.0e4, solver="FullGeneral"):
    """A slave on a (warped) master, free in all 3 DOFs, tangentially anchored by stiff
    springs (so the tangential problem is well-posed), driven by a normal load −P plus a
    tangential load (Q, 0.7Q). Reaching equilibrium the contact point slides ACROSS the
    curvature ⇒ n changes between Newton iterates ⇒ the ∂n/∂u tangent earns quadratic
    convergence. Returns (rc, iters, ux, uy, uz)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _warp_master(w)
    ops.node(1, 0.0, 0.0, -1.0e-3); ops.fix(1, 0, 0, 0)
    ops.node(2, 0.0, 0.0, -1.0e-3); ops.fix(2, 1, 1, 1)
    ops.node(3, 0.0, 0.0, -1.0e-3); ops.fix(3, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, ka)
    ops.element("zeroLength", 91, 2, 1, "-mat", 1, "-dir", 1)   # x anchor
    ops.element("zeroLength", 92, 3, 1, "-mat", 1, "-dir", 2)   # y anchor
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if geom:
        args.append("-geomtan")
    ops.contact(*args)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.7 * Q, -P)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system(solver)
    ops.integrator("LoadControl", 1.0); ops.test("NormDispIncr", 1e-11, 60, 0)
    ops.algorithm("Newton"); ops.analysis("Static")
    rc = ops.analyze(1)
    return rc, ops.testIter(), ops.nodeDisp(1, 1), ops.nodeDisp(1, 2), ops.nodeDisp(1, 3)


# ------------------------------------------------------------------- gate (3)
def test_b3_flat_master_byte_identical():
    """Gate 3: on a FLAT master the geometric block is identically 0 ⇒ `-geomtan` gives
    the SAME converged displacement AND the SAME Newton iteration count as the shipped
    kn·BᵀB tangent (byte-identical behavior)."""
    rc0, it0, x0, y0, z0 = _slide_model(geom=False, w=0.0)
    rc1, it1, x1, y1, z1 = _slide_model(geom=True, w=0.0)
    assert rc0 == 0 and rc1 == 0
    assert it0 == it1, f"flat master must not change the tangent: iters {it0} vs {it1}"
    assert abs(x1 - x0) <= 1e-14 + 1e-12 * abs(x0)
    assert abs(y1 - y0) <= 1e-14 + 1e-12 * abs(y0)
    assert abs(z1 - z0) <= 1e-14 + 1e-12 * abs(z0)


# ------------------------------------------------------------------- gate (1)
def test_b3_curved_quadratic_convergence():
    """Gate 1: on a CURVED (warped) master the consistent ∂n/∂u tangent converges in
    STRICTLY FEWER Newton iterations than kn·BᵀB-only, to the SAME root."""
    rc_m, it_m, xm, ym, zm = _slide_model(geom=False, w=0.6)   # main term only
    rc_g, it_g, xg, yg, zg = _slide_model(geom=True, w=0.6)    # + geometric tangent
    assert rc_m == 0 and rc_g == 0, "both tangents must converge (same exact residual)"
    assert it_g < it_m, f"geometric tangent must be faster on a curve: geom {it_g} vs main {it_m}"
    # SAME equilibrium (only the tangent differs; the residual is identical)
    assert abs(xg - xm) < 1e-8 and abs(yg - ym) < 1e-8 and abs(zg - zm) < 1e-8, \
        f"roots drift: geom=({xg},{yg},{zg}) main=({xm},{ym},{zm})"


# ------------------------------------------------------------------- gate (5)
def test_b3_geomtan_symmetric_not_solver_corrupting():
    """Gate 5 — SYMMETRY-safety: the geometric NORMAL tangent is SYMMETRIC (the Hessian of
    the scalar gap), so a symmetric solver (ProfileSPD, which factors only the upper
    triangle) produces the IDENTICAL result to FullGeneral — proving the assembled block
    has no non-symmetric component that a symmetric solver would silently drop (unlike the
    Coulomb friction Csl). NOTE this proves SYMMETRY, not unconditional convergence: the
    geometric block is INDEFINITE (kn·gN·H, gN<0), so on a large soft-master patch it can
    shrink the Newton basin (see the module docstring + LEDGER_quirks); this is a small,
    well-conditioned problem where it converges. If the tangent had a non-symmetric term,
    ProfileSPD and FullGeneral would DIVERGE here — they don't."""
    rc_f, it_f, xf, yf, zf = _slide_model(geom=True, w=0.6, solver="FullGeneral")
    rc_s, it_s, xs, ys, zs = _slide_model(geom=True, w=0.6, solver="ProfileSPD")
    assert rc_f == 0 and rc_s == 0, "geometric tangent must converge on ProfileSPD"
    assert it_s == it_f, f"symmetric solver took a different path: {it_s} vs {it_f} iters"
    assert abs(xs - xf) < 1e-9 and abs(ys - yf) < 1e-9 and abs(zs - zf) < 1e-9


# --------------------------------------------------------------- gate (2) Hertz
def test_b3_contact_force_query_exact():
    """The `ladrunoContactForce` query returns a slave node's EXACT normal contact force
    (the adapter traction is NOT in nodeReaction). On a single slave loaded by P onto a
    fixed flat master, equilibrium ⇒ the query == P."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(-10, -10), (10, -10), (10, 10), (-10, 10)]):
        ops.node(101 + i, float(x), float(y), 0.0); ops.fix(101 + i, 1, 1, 1)
    ops.node(1, 0.0, 0.0, -1.0e-3); ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, 101, 102, 103, 104)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, 1.0e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -250.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0); ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert abs(ops.ladrunoContactForce(1) - 250.0) / 250.0 < 1e-6, ops.ladrunoContactForce(1)
    assert ops.ladrunoContactForce(99999) == 0.0   # unknown node ⇒ 0


def test_b3_hertz_contact_radius():
    """Hertz benchmark — a rigid sphere (R) indenting an elastic half-space by δ: the
    contact-patch radius reproduces the analytic Hertz radius a = √(Rδ) (rigid indenter,
    compliant body) to a coarse-mesh tolerance."""
    R, delta = 20.0, 0.004
    a_hertz = math.sqrt(R * delta)                   # 0.283
    rc, data = _hertz_indent(R=R, delta=delta)
    assert rc == 0, "Hertz indentation must converge"
    assert len(data) > 30, f"too few contact nodes ({len(data)}) to resolve the patch"
    a_fe = max(r for r, p, f in data)
    assert abs(a_fe - a_hertz) / a_hertz < 0.15, f"contact radius FE={a_fe:.4f} vs √(Rδ)={a_hertz:.4f}"


def test_b3_hertz_geomtan_converges_same_patch():
    """Exercise `-geomtan` (the B3 geometric tangent, `normalGeomTangent`) on the CURVED
    Hertz interface: the indentation converges with the geometric tangent ON and reproduces
    the SAME contact-patch radius as the shipped `kn·BᵀB` path — the residual (hence the
    converged contact geometry) is tangent-INVARIANT; only the Newton path differs."""
    rc0, d0 = _hertz_indent(geom=False)
    rc1, d1 = _hertz_indent(geom=True)
    assert rc0 == 0 and rc1 == 0, "Hertz must converge with -geomtan on this patch"
    a0 = max(r for r, p, f in d0)
    a1 = max(r for r, p, f in d1)
    assert abs(a1 - a0) < 1e-9, f"geomtan changed the converged patch: {a1} vs {a0}"


def test_b3_hertz_pressure_decreasing():
    """The contact pressure is HERTZIAN in shape — peaked at the center, decreasing to the
    edge: the mean pressure over the inner third of the patch is well above the mean over
    the outer third, and the peak sits in the inner half. (Per-node penalty pressures are
    noisy at this resolution, so the check is on radial-bin means, not node-by-node.)"""
    rc, data = _hertz_indent()
    assert rc == 0
    a = max(r for r, p, f in data)
    inner = [p for r, p, f in data if r < a / 3.0]
    outer = [p for r, p, f in data if r > 2.0 * a / 3.0]
    assert inner and outer
    mean_in = sum(inner) / len(inner)
    mean_out = sum(outer) / len(outer)
    assert mean_in > 1.5 * mean_out, f"pressure not center-peaked: inner {mean_in:.1f} vs outer {mean_out:.1f}"
    r_peak = max(data, key=lambda t: t[1])[0]
    assert r_peak < 0.5 * a, f"peak pressure at r={r_peak:.3f} not in the inner half (a={a:.3f})"


def test_b3_hertz_load_grows_with_depth():
    """The Hertz load–indentation relation: a deeper approach carries a larger load over a
    wider patch, with a ∝ √δ (radius ratio √2 to 0.30 for a doubled depth) and P ∝ δ^1.5
    (to 15 %). Both depths use the SAME master (span sized for the deeper δ), so the patch
    radius is a measured outcome, not a construction artifact."""
    R, d1v, d2v = 20.0, 0.003, 0.006
    span = 2.0 * math.sqrt(R * d2v)
    rc1, d1 = _hertz_indent(R=R, delta=d1v, span=span)
    rc2, d2 = _hertz_indent(R=R, delta=d2v, span=span)
    assert rc1 == 0 and rc2 == 0
    P1 = sum(f for r, p, f in d1)
    P2 = sum(f for r, p, f in d2)
    a1 = max(r for r, p, f in d1)
    a2 = max(r for r, p, f in d2)
    assert P2 > P1, f"deeper indentation must carry more load: {P2:.2f} vs {P1:.2f}"
    assert a2 > a1, f"deeper indentation must have a wider patch: {a2:.3f} vs {a1:.3f}"
    assert abs((a2 / a1) - math.sqrt(d2v / d1v)) < 0.30, f"radius ratio {a2/a1:.3f} vs √2={math.sqrt(2):.3f}"
    ratio = (d2v / d1v) ** 1.5
    assert abs(P2 / P1 - ratio) / ratio < 0.15, f"load ratio {P2/P1:.3f} vs (δ2/δ1)^1.5={ratio:.3f}"


# --------------------------------------------------------- gate (6) released pairs
def _flat_master(L, base=101):
    h = L / 2.0
    tags = []
    for i, (x, y) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        ops.node(base + i, x, y, 0.0)
        ops.fix(base + i, 1, 1, 1)
        tags.append(base + i)
    return tags


def test_b3_contact_force_is_zero_after_release():
    """A slave on a soft z spring pressed into a fixed facet, then pulled back out: the
    reported normal force goes from the penalty value to exactly 0.0 (before the r6
    backport it kept the last active value)."""
    KN1, KS, P, Q, g0 = 1.0e6, 1.0e3, 1.0e2, 50.0, 1.0e-8
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _flat_master(L=2.0)
    ops.node(1, 0.0, 0.0, -g0); ops.fix(1, 1, 1, 0)
    ops.node(2, 0.0, 0.0, -g0); ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 3)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN1, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0); ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(1) == 0
    u = -(P - KN1 * g0) / (KN1 + KS)                 # kn (g0 - u) + KS (-u) = P
    f_ref = KN1 * (g0 - u)
    f1 = ops.ladrunoContactForce(1)
    assert abs(f1 - f_ref) / f_ref < 1e-6, (f1, f_ref)

    ops.remove("loadPattern", 1)
    ops.timeSeries("Constant", 2); ops.pattern("Plain", 2, 2); ops.load(1, 0.0, 0.0, Q)
    assert ops.analyze(1) == 0
    z = -g0 + ops.nodeDisp(1, 3)
    assert abs(z - (-g0 + Q / KS)) < 1e-9, z         # separated, spring only
    assert ops.ladrunoContactForce(1) == 0.0


def test_b3_contact_force_is_zero_after_rebound():
    """Explicit drop on a fixed facet: once the slave has bounced clear the reported
    normal force is exactly 0.0."""
    KN1, gap0 = 1.0e6, 0.01
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _flat_master(L=1.0)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN1, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.setNodeVel(1, 3, -2.0, "-commit")
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("NewmarkExplicit", 0.5); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = math.sqrt(1.0 / KN1)
    f_max = 0.0
    for _ in range(4000):
        assert ops.analyze(1, dt) == 0
        f_max = max(f_max, ops.ladrunoContactForce(1))
    assert f_max > 0.0
    assert gap0 + ops.nodeDisp(1, 3) > 0.5 * gap0
    assert ops.ladrunoContactForce(1) == 0.0
