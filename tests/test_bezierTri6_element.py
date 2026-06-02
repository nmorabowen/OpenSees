"""BezierTri6 (Ladruno quadratic Bézier triangle, classTag 33000) — T0 + T1.

Element class merge gate (per manifest) = T0 + T1 + T2a. This file delivers the
T0 battery and a T1 closed-form; the T2a convergence benchmark (Cook's membrane)
is the remaining pending piece.

T0:
  * rank / zero-energy: a free element has EXACTLY 3 zero modes (2D rigid body).
  * constant-strain patch: 2-triangle patch, linear field on the boundary,
    diagonal-midside node free -> every Gauss-point stress is the closed-form
    constant stress to 1e-9.
  * quadratic completeness: the distinguishing test for a *quadratic* element.
    A quadratic displacement field is in the element's span, so with the
    consistent body force b = -div(sigma) the patch reproduces it exactly and
    every Gauss-point strain equals the analytic (linear) strain at that GP's
    physical location.
T1:
  * uniaxial panel (displacement-controlled): total support reaction = E*eps*A
    and the lateral contraction = -nu*eps*y, both exact for a constant-strain-
    reproducing element.

Command: element('BezierTri6', tag, n1..n6, thick, type, matTag,
                 '-bodyForce', b1, b2). Straight-sided only (D9 guard): the 3
midside nodes sit at edge midpoints. NGAUSS = 3, ndf = 2.
"""
import math
from math import sin, cos, pi

import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy, check_constant_stress
from _testbed.equilibrium import assert_equilibrium

pytestmark = [pytest.mark.zone_a]

E = 1000.0
NU = 0.25
THK = 0.7
F_PS = E / (1.0 - NU * NU)          # plane-stress D scalar
G_PS = E / (2.0 * (1.0 + NU))
NGAUSS = 3

# ---- shared distorted 2-triangle patch (node 7 = diagonal midside = interior)
CORNERS = {1: (0.0, 0.0), 2: (2.0, 0.0), 3: (2.2, 1.8), 4: (0.1, 2.0)}


def _mid(a, b):
    return (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))


def _patch_coords():
    c = CORNERS
    return {
        **c,
        5: _mid(c[1], c[2]),   # edge 1-2
        6: _mid(c[2], c[3]),   # edge 2-3
        7: _mid(c[1], c[3]),   # diagonal 1-3  -> interior, free
        8: _mid(c[3], c[4]),   # edge 3-4
        9: _mid(c[4], c[1]),   # edge 4-1
    }


# element node order: 3 corners then 3 midsides nd4=mid(n1,n2), nd5=mid(n2,n3),
# nd6=mid(n3,n1) (standard FE hierarchical, per BezierTri6 + project memory)
TRI_A = (1, 2, 3, 5, 6, 7)
TRI_B = (1, 3, 4, 7, 8, 9)
BOUNDARY = [1, 2, 3, 4, 5, 6, 8, 9]   # everything except interior node 7


def _make_patch(coords, body=(0.0, 0.0)):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for n, (x, y) in coords.items():
        ops.node(n, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for e, nodes in ((1, TRI_A), (2, TRI_B)):
        ops.element("BezierTri6", e, *nodes, THK, "PlaneStress", 1,
                    "-bodyForce", body[0], body[1])


def _solve_with_prescribed(coords, ux_of, uy_of, body=(0.0, 0.0)):
    """Impose ux_of/uy_of(x,y) at the boundary nodes, node 7 free, then solve."""
    _make_patch(coords, body=body)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in BOUNDARY:
        x, y = coords[n]
        ops.sp(n, 1, ux_of(x, y))
        ops.sp(n, 2, uy_of(x, y))
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)


# ----------------------------------------------------------------- T0: rank
@pytest.mark.t0
def test_rank_zero_energy():
    coords = _patch_coords()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for n in TRI_A:
        ops.node(n, *coords[n])
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("BezierTri6", 1, *TRI_A, THK, "PlaneStress", 1)
    assert_zero_energy(list(TRI_A), ndf=2, ndm=2)


# ----------------------------------------------------------------- T0: patch
@pytest.mark.t0
def test_constant_strain_patch():
    coords = _patch_coords()
    g1, g2, g3, g4 = 1e-3, 3e-4, 2e-4, -4e-4
    eps_xx, eps_yy, gam_xy = g1, g4, g2 + g3
    _solve_with_prescribed(
        coords,
        ux_of=lambda x, y: g1 * x + g2 * y,
        uy_of=lambda x, y: g3 * x + g4 * y,
    )
    # interior node lands on the field
    x7, y7 = coords[7]
    assert abs(ops.nodeDisp(7, 1) - (g1 * x7 + g2 * y7)) <= 1e-9 + 1e-9
    target = [F_PS * (eps_xx + NU * eps_yy),
              F_PS * (NU * eps_xx + eps_yy),
              G_PS * gam_xy]
    sref = max(abs(eps_xx), abs(eps_yy), abs(gam_xy))
    for e in (1, 2):
        check_constant_stress(e, NGAUSS, target, E=E, strain_ref=sref)


# ------------------------------------------------- T0: quadratic completeness
#
# CRITICAL Bezier subtlety: BezierTri6 DOFs are Bernstein CONTROL-POINT
# COEFFICIENTS, not nodal field values. The basis is non-interpolatory at
# midside nodes, so to impose a field p(x,y) you must prescribe the control net,
# not point samples:
#   * vertex DOF  = p(vertex)
#   * midside DOF on edge (va,vb) = 2*p(mid) - 0.5*p(va) - 0.5*p(vb)   (deg-2
#     blossom / polar form). For a LINEAR p this reduces to p(mid) — which is why
#     the constant-strain patch above could prescribe field values directly.
# Prescribing p(mid) at midsides for a QUADRATIC field silently corrupts the
# test (the element then reproduces a different field). Verified empirically:
# with the control-net BCs + the consistent body force the max GP-strain error
# is ~1e-19 (machine zero); with naive p(mid) BCs it is ~1e-4.

_MIDSIDE_EDGES = {5: (1, 2), 6: (2, 3), 8: (3, 4), 9: (4, 1)}  # boundary midsides


def _bezier_control_coeff(fn, mid_xy, va_xy, vb_xy):
    """Degree-2 Bernstein control coefficient for the midside DOF of a polynomial
    field fn, from its values at the edge midpoint and the two endpoint vertices."""
    return (2.0 * fn(*mid_xy) - 0.5 * fn(*va_xy) - 0.5 * fn(*vb_xy))


@pytest.mark.t0
def test_quadratic_completeness_patch():
    coords = _patch_coords()
    a, c = 1.0e-4, -0.5e-4                       # u_x = a x^2 , u_y = c y^2
    def ux(x, y): return a * x * x
    def uy(x, y): return c * y * y
    # body force b = -div(sigma*); sigma* linear -> div sigma* = (2 a F, 2 c F)
    body = (-2.0 * a * F_PS, -2.0 * c * F_PS)

    _make_patch(coords, body=body)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for v in (1, 2, 3, 4):                        # vertex DOF = field value
        ops.sp(v, 1, ux(*coords[v]))
        ops.sp(v, 2, uy(*coords[v]))
    for m, (va, vb) in _MIDSIDE_EDGES.items():    # midside DOF = blossom coeff
        ops.sp(m, 1, _bezier_control_coeff(ux, coords[m], coords[va], coords[vb]))
        ops.sp(m, 2, _bezier_control_coeff(uy, coords[m], coords[va], coords[vb]))
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)

    # every GP strain == analytic linear strain (2 a x, 2 c y, 0) at the GP
    for e in (1, 2):
        gp = ops.eleResponse(e, "gpCoord")        # [x1,y1,x2,y2,x3,y3]
        st = ops.eleResponse(e, "strain")         # [exx,eyy,gxy] x 3
        assert gp and st, f"element {e} missing gpCoord/strain response"
        for k in range(NGAUSS):
            xg, yg = gp[2 * k], gp[2 * k + 1]
            exx, eyy, gxy = st[3 * k], st[3 * k + 1], st[3 * k + 2]
            assert abs(exx - 2.0 * a * xg) <= 1e-9 * abs(2.0 * a * xg) + 1e-12
            assert abs(eyy - 2.0 * c * yg) <= 1e-9 * abs(2.0 * c * yg) + 1e-12
            assert abs(gxy) <= 1e-10


# ----------------------------------------------------------------- T1: uniaxial
@pytest.mark.t1
def test_uniaxial_reaction_and_poisson():
    Lx, Ly, thk = 4.0, 2.0, 0.5
    delta = 0.01
    eps = delta / Lx                       # uniform axial strain
    coords = {
        1: (0.0, 0.0), 2: (Lx, 0.0), 3: (Lx, Ly), 4: (0.0, Ly),
        5: (Lx / 2, 0.0), 6: (Lx, Ly / 2), 7: (Lx / 2, Ly / 2),
        8: (Lx / 2, Ly), 9: (0.0, Ly / 2),
    }
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for n, (x, y) in coords.items():
        ops.node(n, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("BezierTri6", 1, 1, 2, 3, 5, 6, 7, thk, "PlaneStress", 1)
    ops.element("BezierTri6", 2, 1, 3, 4, 7, 8, 9, thk, "PlaneStress", 1)

    left = [1, 4, 9]                        # x = 0
    right = [2, 3, 6]                       # x = Lx
    for n in left:
        ops.fix(n, 1, 0)                    # roller: u_x = 0, u_y free
    ops.fix(1, 0, 1)                        # pin one point in y (no rigid body)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in right:
        ops.sp(n, 1, delta)                 # u_x = delta on the loaded edge
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)
    ops.reactions()

    # total axial support reaction = E * eps * A  (plane-stress uniaxial)
    Rx = sum(ops.nodeReaction(n, 1) for n in left)
    A = Ly * thk
    assert abs(abs(Rx) - E * eps * A) <= 1e-6 * E * eps * A

    # lateral Poisson contraction: u_y(node 4 at y=Ly) = -nu * eps * Ly
    assert abs(ops.nodeDisp(4, 2) - (-NU * eps * Ly)) <= 1e-9 + 1e-9 * NU * eps * Ly

    # force balance closes on the translational DOFs (no external x load -> the
    # right-edge reactions balance the left; net applied force is zero)
    assert_equilibrium(list(coords), n_force_dofs=2, applied_total=[0.0, 0.0], atol=1e-6)


# --------------------------------------------- T2a: convergence rate (smooth MMS)
#
# Cook's membrane is the usual continuum benchmark, but its clamped-corner stress
# singularity caps the DISPLACEMENT convergence rate at ~1.3 (measured here) and
# its quoted "23.96" is itself a fine-mesh FE value, not analytic -- so Cook is an
# ACCURACY benchmark, not a rate test. For the rate>1.8 gate we use a smooth
# manufactured solution on the unit square: no singularity, exact field on the
# whole boundary via the Bernstein control-net (T0 above), and a per-element
# constant body force b=-div(sigma) sampled at the centroid. The body-force
# sampling error is O(h^2), matching the quadratic element, so the max corner-
# displacement error converges at ~2.0-2.2 -- comfortably above 1.8. (Rectangular
# domain => each cell is affine => the diagonal midside node IS the geometric
# midpoint, so the straight-sided guard is satisfied with the simple grid.)

def _mms_u(x, y):
    s = sin(pi * x) * sin(pi * y)
    return s, s


def _mms_body(x, y):
    f = E / (1.0 - NU * NU)
    g = E / (2.0 * (1.0 + NU))
    S = sin(pi * x) * sin(pi * y)
    Cc = cos(pi * x) * cos(pi * y)
    v = pi * pi * ((f + g) * S - (f * NU + g) * Cc)
    return v, v


def _solve_mms(N):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    M = 2 * N

    def nid(i, j):
        return j * (M + 1) + i + 1

    for j in range(M + 1):
        for i in range(M + 1):
            ops.node(nid(i, j), i / M, j / M)

    eid = 1
    for a in range(N):
        for b in range(N):
            bl, br = nid(2 * a, 2 * b), nid(2 * a + 2, 2 * b)
            tr, tl = nid(2 * a + 2, 2 * b + 2), nid(2 * a, 2 * b + 2)
            bot, rgt, ctr = nid(2 * a + 1, 2 * b), nid(2 * a + 2, 2 * b + 1), nid(2 * a + 1, 2 * b + 1)
            top, lft = nid(2 * a + 1, 2 * b + 2), nid(2 * a, 2 * b + 1)
            for conn in ((bl, br, tr, bot, rgt, ctr), (bl, tr, tl, ctr, top, lft)):
                cx = sum(ops.nodeCoord(n)[0] for n in conn[:3]) / 3.0
                cy = sum(ops.nodeCoord(n)[1] for n in conn[:3]) / 3.0
                bx, by = _mms_body(cx, cy)
                ops.element("BezierTri6", eid, *conn, 1.0, "PlaneStress", 1,
                            "-bodyForce", bx, by)
                eid += 1

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(M + 1):
        for i in range(M + 1):
            if not (i in (0, M) or j in (0, M)):
                continue
            x, y = i / M, j / M
            if i % 2 == 0 and j % 2 == 0:                  # boundary vertex
                ux, uy = _mms_u(x, y)
            else:                                          # boundary midside -> blossom
                if i % 2 == 1:                             # horizontal edge
                    ua, ub = _mms_u((i - 1) / M, y), _mms_u((i + 1) / M, y)
                else:                                      # vertical edge
                    ua, ub = _mms_u(x, (j - 1) / M), _mms_u(x, (j + 1) / M)
                um = _mms_u(x, y)
                ux = 2 * um[0] - 0.5 * ua[0] - 0.5 * ub[0]
                uy = 2 * um[1] - 0.5 * ua[1] - 0.5 * ub[1]
            ops.sp(nid(i, j), 1, ux)
            ops.sp(nid(i, j), 2, uy)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)

    err = 0.0
    for j in range(0, M + 1, 2):
        for i in range(0, M + 1, 2):
            x, y = i / M, j / M
            ex, ey = _mms_u(x, y)
            err = max(err,
                      abs(ops.nodeDisp(nid(i, j), 1) - ex),
                      abs(ops.nodeDisp(nid(i, j), 2) - ey))
    return err


@pytest.mark.t2a
def test_mms_convergence_rate():
    e = {N: _solve_mms(N) for N in (4, 8, 16)}
    assert e[4] > e[8] > e[16], f"non-monotone convergence: {e}"
    rate = math.log(e[8] / e[16]) / math.log(2.0)
    assert rate > 1.8, f"convergence rate {rate:.2f} <= 1.8 (errors {e})"


# ----------------------------------------------- T0: characteristic length
#
# Crack-band materials (ASDConcrete3D) regularize the softening branch by
# ops_TheActiveElement->getCharacteristicLength(), read once on the first
# setTrialStrain. BezierTri6 OVERRIDES Element's min-inter-node-distance
# default (which on a quadratic element collapses to ~½ the edge length:
# corner-to-mid-edge) with an area-based element size lch = sqrt(2·A).
#
# For a right-isosceles triangle with legs of length L, A = L²/2, so the
# override returns lch = sqrt(2·L²/2) = L exactly. The inherited default
# would instead return the min node spacing L/2 — so this test both pins the
# formula AND fails loudly if anyone reverts the override (½× regression).

def _make_right_triangle(L):
    """Single BezierTri6, right-isosceles, legs L along +x and +y. Returns ele tag 1."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    c1, c2, c3 = (0.0, 0.0), (L, 0.0), (0.0, L)
    coords = {
        1: c1, 2: c2, 3: c3,
        4: _mid(c1, c2), 5: _mid(c2, c3), 6: _mid(c3, c1),
    }
    for n, (x, y) in coords.items():
        ops.node(n, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("BezierTri6", 1, 1, 2, 3, 4, 5, 6, THK, "PlaneStress", 1)
    return 1


@pytest.mark.t0
@pytest.mark.parametrize("L", [1.0, 3.0, 0.4])
def test_characteristic_length(L):
    ele = _make_right_triangle(L)
    lch = ops.eleResponse(ele, "charLength")
    assert lch, "BezierTri6 missing 'charLength' response"
    lch = lch[0]
    # area-based override: lch = sqrt(2·A) = L for the right-isosceles triangle
    assert abs(lch - L) <= 1e-12 * L + 1e-14, f"lch={lch}, expected {L}"
    # regression guard: must NOT be the inherited min-distance default (L/2)
    assert lch > 0.75 * L, f"lch={lch} looks like the ½-edge min-distance default"
