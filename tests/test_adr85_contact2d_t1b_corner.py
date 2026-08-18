"""ADR-85 T1b -- the 2D CORNER primitive in the compiled engine.

WHAT IS BEING GATED
    T1a settled the corner geometry in the header-only kernel
    (SRC/domain/contact/LadrunoContact2DKernel.h) and INVERTED the ADR's rev-2
    premise while doing it.  T1b is where the handler acts on that finding, and
    these are the two cases that separate a correct handler from a plausible
    one:

    CONCAVE (re-entrant) corners are the real coverage HOLE.  The region behind
    a re-entrant corner, inside the material, is refused by BOTH adjacent
    segments (the projection runs past the high end of one and before the low
    end of the other) while being ENTIRELY penetrating -- a finite-measure wedge
    at finite penetration with zero restoring force.  The T1a oracle measured
    that arc: length R*|turn| at penetration depth R, 19.7% of a sweep path at
    R = 0.25, scaling linearly.  Silent pass-through.  So the handler must
    instantiate a VERTEX PAIR there, and test_concave_corner_holds /
    test_concave_corner_arc_covered are the two halves of that claim (the
    symmetry axis with an analytic answer, then the whole arc).

    CONVEX corners are the opposite: their exterior Voronoi wedge lies entirely
    outside the material, so a convex vertex pair can never be penetrating and
    never carries force.  Their interior is covered TWICE, by both adjacent
    segment strips -- a DOUBLE-COUNT, not a hole.  test_convex_corner_slide is
    the ownership gate: sliding around the corner must stay single-owner and
    continuous.

HOW THE SLAVE IS PLACED, AND WHY IT IS NOT `sp`
    Two of the three tests need the slave at a PRESCRIBED position (a station on
    a sweep) rather than wherever equilibrium puts it, and reading a force back.
    Prescribing both DOFs with `fix` would leave the contact FE with no free
    equation, and imposing them with a pattern `sp` puts a non-homogeneous
    single-point constraint through a handler whose SP behavior is not part of
    any shipped contact gate.

    So the sweeps use an ELASTIC TETHER instead: a coincident fixed ground node
    plus a `zeroLength` with stiffness KS in both directions.  Three properties
    make it the right instrument rather than a workaround:

      * the station is a genuine static equilibrium with free DOFs -- no N = 0
        SOE, no all-constrained FE, no handler behavior anyone has to assume;
      * the tether IS the force gauge.  The slave's only other force is the
        contact, so |F_contact| = KS*|u| exactly, measured from displacement
        alone.  That matters: `ladrunoContactForce` is fed from the SEGMENT
        branch's setNtsForce (LadrunoContactFE.cpp:689) and whether a VERTEX
        pair reports through it is a T1b implementation choice, not a contract.
        The concave tests therefore never depend on the query;
      * the answer stays CLOSED FORM.  Contact and tether act along the same
        line (the contact normal), so they are two linear springs in series
        against the geometric depth d:

            KN*(d - s) = KS*s   =>   F = KEFF*d,   KEFF = KN*KS/(KN + KS)

        With KS = 1000*KN the tether moves the slave by d/1001, so the station
        is where it was placed to within 0.1%, and KEFF is known exactly rather
        than approximated.

TOLERANCES ARE RELATIVE (ADR-85 How/1)
    Forces are compared to KEFF*R or KEFF*d (deck quantities), penetrations to
    the penalty depth P/KN and to the master segment length.  There is no bare
    absolute anywhere.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6         # contact penalty
KS = 1.0e9         # tether stiffness; KS/KN = 1e3 (see the docstring)
KEFF = KN * KS / (KN + KS)      # the exact series stiffness of the two
SEED = 1.0e-8      # seeded penetration on the load-driven case


def _analysis(tol=1.0e-12, iters=40):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _tether(slave, x, y, gnd=900, mat=900, ele=900, ks=KS):
    """Coincident fixed ground + a 2-direction elastic link.  See the module
    docstring: this is the position clamp AND the force gauge."""
    ops.node(gnd, x, y)
    ops.fix(gnd, 1, 1)
    ops.uniaxialMaterial("Elastic", mat, ks)
    ops.element("zeroLength", ele, gnd, slave, "-mat", mat, mat, "-dir", 1, 2)


# --------------------------------------------------------------------------
# CONCAVE (re-entrant) corner: a symmetric V-notch
#
#   master nodes 101 = (-1, 1/2), 102 = (0, 0), 103 = (1, 1/2), all FIXED
#   surface      two 2-node segments (101,102) and (102,103) sharing vertex 102
#   material     everything BELOW the V (y < |x|/2); the free region is the
#                126.87 deg wedge above it, so the MATERIAL wedge is 233.13 deg
#                -- re-entrant, which is the case that has a hole.
#
# The notch is symmetric about x = 0 ON PURPOSE: the vertex normal is then the
# bisector of the two adjacent outward normals, which is exactly +y, so the
# vertex pair projects entirely onto one coordinate DOF and the equilibrium has
# a closed form.  An oblique notch would test the same code with an answer that
# has to be trusted rather than derived.
#
# `-outward 0 1` orients BOTH segments: the outward normals are
# (+0.447, +0.894) and (-0.447, +0.894), and both have a positive dot with +y.
# It is also required rather than optional here -- the reference-configuration
# centroid vote for a slave seeded 1e-8 from the vertex is degenerate.
# --------------------------------------------------------------------------
NOTCH = ((-1.0, 0.5), (0.0, 0.0), (1.0, 0.5))
L_SEG = math.hypot(1.0, 0.5)            # 1.118: the deck's element size

# The unowned wedge, from the T1a predicate (aPrev > 0 and aNext < 0): with
# segment tangents proportional to (2,-1) and (2,1), a direction d = (cos, sin)
# from the vertex is in the wedge iff 2*cos - sin > 0 and 2*cos + sin < 0, i.e.
# phi in (243.43 deg, 296.57 deg) -- an arc of exactly the 53.13 deg turn.
WEDGE_LO, WEDGE_HI = 243.43, 296.57


def _notch_master():
    for i, (x, y) in enumerate(NOTCH):
        ops.node(101 + i, x, y)
        ops.fix(101 + i, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102, 102, 103)


def test_concave_corner_holds():
    """A slave driven straight into the re-entrant corner is HELD.

    Deck: the symmetric V-notch above, slave seeded SEED below the vertex, x
    fixed and y free, pressed down with P over 5 LoadControl steps.

    WITHOUT a vertex pair this deck has NO stiffness at all on the slave's y
    DOF -- both segments refuse the projection, which is the pass-through the
    T1a oracle measured -- so the failure is loud (a singular tangent, or a
    displacement that runs away).  WITH one, the vertex normal is +y by
    symmetry, the pair is rank-1 on the free DOF, and the equilibrium is exact:

        kn * penetration = lambda * P    at every load step

    so the gate is the whole RAMP, not just the endpoint.  A vertex pair whose
    gap were measured to a segment line instead of to the vertex, or whose
    normal came from one adjacent segment (n_y = 0.894) rather than the
    bisector, lands 25% off and is caught by the same assertion.

    Bounded penetration is asserted twice over, in the two units that mean
    different things: against the analytic penalty depth (is the force right?)
    and against the master segment length (did it tunnel?).
    """
    P, nsteps = 1.0e3, 5
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _notch_master()
    ops.node(1, 0.0, -SEED)
    ops.fix(1, 1, 0)                       # x fixed, y free
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    _analysis()
    ops.integrator("LoadControl", 1.0 / nsteps)

    pens = []
    for k in range(nsteps):
        assert ops.analyze(1) == 0, (
            f"concave-corner static step {k} did not converge -- with no vertex "
            "pair the slave's y DOF has zero stiffness (both segments refuse "
            "the projection), which is exactly this signature")
        pens.append(SEED - ops.nodeDisp(1, 2))      # distance below the vertex

    print(f"\n[concave corner] penetration ramp {['%.6e' % q for q in pens]} "
          f"(analytic step {P / (nsteps * KN):.6e})")
    for k, pen in enumerate(pens):
        ref = (k + 1) * P / (nsteps * KN)
        assert pen > 0.0, f"step {k}: the slave was ejected (pen={pen:.3e})"
        assert abs(pen - ref) / ref < 1.0e-6, (
            f"step {k}: penetration {pen:.8e} != P_step/kn {ref:.8e}.  Either "
            "the vertex gap is not the distance to the vertex, or the vertex "
            "normal is not the adjacent-normal bisector (one adjacent segment "
            "normal would give a 25% error on this notch)")
    assert all(b > a for a, b in zip(pens, pens[1:])), (
        f"the penetration ramp is not monotone: {pens}")

    # no tunnelling: the slave stays in the corner region, not through it
    assert pens[-1] < 0.02 * L_SEG, (
        f"final penetration {pens[-1]:.3e} is more than 2% of the master "
        f"segment length {L_SEG:.3f} -- the corner is not holding")
    assert abs(ops.nodeDisp(1, 2)) < 0.05 * L_SEG, (
        f"slave y-displacement {ops.nodeDisp(1, 2):.3e} is a runaway against "
        f"the element size {L_SEG:.3f} (pass-through)")


def test_concave_corner_arc_covered():
    """The WHOLE unowned arc is covered, not just its symmetry axis.

    test_concave_corner_holds sits at phi = 270 deg, where a lucky
    implementation -- say, one that snaps an out-of-bounds projection back to
    the vertex on ONE of the two segments -- can still get the right answer.
    This sweeps the arc the T1a oracle identified as unowned
    (phi in (243.4, 296.6) deg for this notch, the 53.13 deg turn) and asserts
    the vertex contract at each station:

        gap = -||x_s - X_v||   =>   |F| = KEFF * R,   INDEPENDENT of phi

    A constant force over a swept arc is a strong signature: it is what a RADIAL
    (vertex) normal gives and what a segment normal cannot, since any segment's
    gap varies over the arc.  The displacement direction is checked too -- it
    must be radially outward from the vertex, which pins the normal itself and
    not only its magnitude.
    """
    R = 0.15                                # radius, well inside both segments
    phis = [250.0, 260.0, 270.0, 280.0, 290.0]
    assert WEDGE_LO < min(phis) and max(phis) < WEDGE_HI, "stations left the wedge"

    forces = []
    for i, phi in enumerate(phis):
        th = math.radians(phi)
        px, py = R * math.cos(th), R * math.sin(th)
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        _notch_master()
        ops.node(1, px, py)
        _tether(1, px, py)
        ops.contactSurface(20, "-slave", 1)
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
        ops.timeSeries("Linear", 1)         # no loads: the contact is the only
        ops.pattern("Plain", 1, 1)          # thing out of balance at u = 0
        _analysis()
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0, f"concave arc station phi={phi} did not converge"

        ux, uy = ops.nodeDisp(1, 1), ops.nodeDisp(1, 2)
        mag = math.hypot(ux, uy)
        forces.append(KS * mag)
        # the push must be radially OUTWARD from the vertex, i.e. along -d_hat
        assert mag > 0.0, f"phi={phi}: no force at all -- the arc is UNOWNED here"
        cos_align = -(ux * math.cos(th) + uy * math.sin(th)) / mag
        assert cos_align > 0.999, (
            f"phi={phi}: the restoring direction is {math.degrees(math.acos(max(-1.0, min(1.0, cos_align)))):.2f} "
            "deg off the outward radial -- the vertex normal is not "
            "(x_s - X_v)/||x_s - X_v||")

    ref = KEFF * R
    print(f"\n[concave arc] forces {['%.6e' % f for f in forces]} vs KEFF*R = {ref:.6e}")
    for phi, f in zip(phis, forces):
        assert abs(f - ref) / ref < 1.0e-4, (
            f"phi={phi}: |F| = {f:.8e} != KEFF*R = {ref:.8e}.  The vertex gap "
            "must be the distance to the vertex at every station; a segment "
            "owning this point would give a phi-dependent force")


# --------------------------------------------------------------------------
# CONVEX corner: the corner of a quarter-plane block
#
#   master nodes 101 = (-1, 0), 102 = (0, 0), 103 = (0, -1), all FIXED
#   segments     (101,102) with outward +y, (102,103) with outward +x
#   material     the quadrant x < 0, y < 0; the free region spans 270 deg
#
# `-outward 1 1` orients both (dot 0.707 with each).
# --------------------------------------------------------------------------
def _convex_master():
    for i, (x, y) in enumerate(((-1.0, 0.0), (0.0, 0.0), (0.0, -1.0))):
        ops.node(101 + i, x, y)
        ops.fix(101 + i, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102, 102, 103)


def _convex_station(px, py):
    """One tethered slave at (px, py) inside the convex corner.  Returns
    (|F| from the tether, the summed ladrunoContactForce, (ux, uy))."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _convex_master()
    ops.node(1, px, py)
    _tether(1, px, py)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 1.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    _analysis()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, f"convex station ({px:.4f},{py:.4f}) did not converge"
    ux, uy = ops.nodeDisp(1, 1), ops.nodeDisp(1, 2)
    return KS * math.hypot(ux, uy), ops.ladrunoContactForce(1), (ux, uy)


def test_convex_corner_slide():
    """Sliding a slave around a 90 deg convex corner stays SINGLE-OWNER and
    CONTINUOUS.

    The path is the circular arc of radius R about the corner vertex, from
    theta = 185 deg to 265 deg -- i.e. the slave travels from just under the
    top face round to just behind the right face, always inside the material,
    never on a face.  The two candidate depths are

        d1 = R*|sin theta|   (top face)      d2 = R*|cos theta|   (right face)

    THE BAND IS RULE-AGNOSTIC, AND THAT IS DELIBERATE.  The ADR fixes ownership
    as "segment precedence in surface order, then the unslacked vertex claim,
    else no owner"; a closest-feature reading of the same protocol picks the
    other segment on the second half of the arc.  Both are single-owner and both
    are continuous, and this file must not encode a guess about which one T1b
    shipped.  So the gate is the property, not the identity:

        any single owner       =>   KEFF*min(d1,d2) <= F <= KEFF*max(d1,d2)
        both segments active   =>   F = KEFF*(d1 + d2), which at theta = 225 deg
                                    is exactly 2x the upper bound

    That catches the double-count where it is largest (mid-arc) and, at the ends
    where a second pair is nearly separated and undetectable, the band is simply
    honest about not resolving it.  A hole gives F = 0 and is caught everywhere.

    CONTINUITY IS ADDITIVE, NOT A RATIO.  The task-level phrasing is "no spike
    beyond a factor of 2 between consecutive stations", but a ratio test is
    wrong on this path: near theta = 185 deg the top-face depth is 0.087R and
    the next station's is 0.259R, a legitimate ratio of 3.0 that has nothing to
    do with a spike.  Both candidate depths are Lipschitz in theta with constant
    R, so the honest bound is |dF| <= 1.5 * KEFF * R * dtheta, and an ownership
    JUMP -- the actual failure -- blows straight through it.

    The tether magnitude and the query are BOTH read.  They fail differently: a
    double count makes the query (a scalar sum of tractions) report d1 + d2
    while the tether reports the vector magnitude sqrt(d1^2 + d2^2), so the two
    together also catch a single pair with a doubled traction.
    """
    R, tol = 0.2, 0.02
    thetas = [185.0 + 10.0 * k for k in range(9)]        # 185 .. 265 deg

    # calibration: a station far from the corner, unambiguously owned by the top
    # face under every rule.  If this one is off, nothing below means anything.
    depth = 0.1
    f_teth, f_query, _u = _convex_station(-0.6, -depth)
    assert abs(f_teth - KEFF * depth) / (KEFF * depth) < 1.0e-4, (
        f"flat-face reference: tether force {f_teth:.8e} != KEFF*d "
        f"{KEFF * depth:.8e} -- the gauge itself is wrong")
    assert abs(f_query - KEFF * depth) / (KEFF * depth) < 1.0e-4, (
        f"flat-face reference: ladrunoContactForce {f_query:.8e} != KEFF*d "
        f"{KEFF * depth:.8e}")

    rows = []
    for th_deg in thetas:
        th = math.radians(th_deg)
        px, py = R * math.cos(th), R * math.sin(th)
        f_teth, f_query, _u = _convex_station(px, py)
        d1, d2 = R * abs(math.sin(th)), R * abs(math.cos(th))
        rows.append((th_deg, f_teth, f_query, KEFF * min(d1, d2), KEFF * max(d1, d2)))

    print("\n[convex slide] theta, F_tether, F_query, band_lo, band_hi")
    for th_deg, ft, fq, lo, hi in rows:
        print(f"  {th_deg:6.1f}  {ft:.6e}  {fq:.6e}  [{lo:.6e}, {hi:.6e}]")

    for th_deg, ft, fq, lo, hi in rows:
        assert fq > 0.0, (
            f"theta={th_deg}: NO contact force at a point strictly inside the "
            "material -- the convex corner has a coverage hole")
        for name, val in (("tether", ft), ("query", fq)):
            assert val >= (1.0 - tol) * lo, (
                f"theta={th_deg}: {name} force {val:.6e} below the "
                f"closest-feature depth {lo:.6e} -- the owner is shallower than "
                "either candidate")
            assert val <= (1.0 + tol) * hi, (
                f"theta={th_deg}: {name} force {val:.6e} above the deepest "
                f"single candidate {hi:.6e}.  Both adjacent segment strips "
                f"contain this point; a DOUBLE COUNT reports "
                f"{KEFF * R * (abs(math.sin(math.radians(th_deg))) + abs(math.cos(math.radians(th_deg)))):.6e}")

    dtheta = math.radians(thetas[1] - thetas[0])
    lip = 1.5 * KEFF * R * dtheta
    for (t0, f0, _q0, _l0, _h0), (t1, f1, _q1, _l1, _h1) in zip(rows, rows[1:]):
        assert abs(f1 - f0) <= lip, (
            f"force jumped {abs(f1 - f0):.6e} between theta={t0} and {t1}, "
            f"beyond the Lipschitz bound {lip:.6e} (both candidate depths vary "
            "at most R per radian).  The segment-to-segment handoff is not C0 "
            "-- an ownership switch between two candidates of different depth")
