"""ASDShellQ4 + LayeredShell + PlateRebar: which in-plane direction the bar runs (WP-142, O8).

A PlateRebar layer at angle 0 lies along the SECTION x axis, and the section x axis is
whatever frame ASDShellQ4 hands its sections. That frame is:

  * with ``-local x1 x2 x3``: the given vector projected onto the element plane
    (ASDShellQ4.cpp setDomain, the ``m_local_x`` branch);
  * without ``-local``: the line from the midpoint of side 1-4 to the midpoint of side 2-3
    (the element xi direction), ``e1 = (P2+P3)/2 - (P1+P4)/2`` (ASDShellQ4.cpp:866).

Upstream OpenSees v3.7.1 declares a NEW ``Vector3Type e1`` inside that default branch, so the
outer ``e1`` stays (0,0,0), ``m_angle = acos(0) = +pi/2`` and the default section frame is
rotated 90 degrees from the element's side 1-2 (fixed upstream by b58e0795d, first in v3.8.0,
which this fork carries). A model without ``-local`` therefore puts its 0-degree bars along X
in this fork and along Y in v3.7.1, and in both it depends on the node ordering. ``-local``
takes the same branch in both versions: pin it whenever a section has PlateRebar layers
(LEDGER_quirks).

Oracle (T1, closed form). Concrete ElasticIsotropic with nu = 0, one elastic bar layer of
thickness t at angle 0, total thickness h. A uniaxial membrane strain eps along the in-plane
unit vector d (all other strains zero) gives

    N_d / eps = E_c (h - t) + E_s t cos^4(delta),    delta = angle(d, bar),

so along the bar K = E_c(h - t) + E_s t and across it K = E_c(h - t). K is measured
frame-free from the element's global nodal forces: K = sum(F . u) / (eps^2 A).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E_C = 30000.0
E_S = 200000.0
T_C = 0.10            # each concrete layer
T_S = 0.002           # bar layer
H = 2.0 * T_C + T_S
K_ACROSS = E_C * (H - T_S)
K_ALONG = K_ACROSS + E_S * T_S
EPS = 1.0e-4
REL = 1.0e-6

SQUARE = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]


def _rot(pts, phi_deg):
    c, s = math.cos(math.radians(phi_deg)), math.sin(math.radians(phi_deg))
    return [(c * x - s * y, s * x + c * y) for x, y in pts]


def _area(pts):
    a = 0.0
    for i in range(len(pts)):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % len(pts)]
        a += x0 * y1 - x1 * y0
    return 0.5 * abs(a)


def _unit(theta_deg):
    return math.cos(math.radians(theta_deg)), math.sin(math.radians(theta_deg))


def membrane_stiffness(pts, theta_deg, local=None, order=(1, 2, 3, 4)):
    """N/eps for a uniaxial membrane strain along theta (degrees from global X)."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for n, (x, y) in enumerate(pts, start=1):
        ops.node(n, x, y, 0.0)
    ops.nDMaterial('ElasticIsotropic', 1, E_C, 0.0)
    ops.uniaxialMaterial('Elastic', 11, E_S)
    ops.nDMaterial('PlateRebar', 21, 11, 0.0)
    ops.section('LayeredShell', 10, 3, 1, T_C, 21, T_S, 1, T_C)
    args = ['ASDShellQ4', 1, *order, 10]
    if local is not None:
        args += ['-local', *local]
    ops.element(*args)

    dx, dy = _unit(theta_deg)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y) in enumerate(pts, start=1):
        a = EPS * (dx * x + dy * y)            # u = eps (d . X) d
        ops.sp(n, 1, a * dx)
        ops.sp(n, 2, a * dy)
        for dof in (3, 4, 5, 6):
            ops.sp(n, dof, 0.0)
    ops.constraints('Penalty', 1.0e14, 1.0e14)
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-14, 10, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0

    f = list(ops.eleForce(1))
    nodes = list(ops.eleNodes(1))
    work = 0.0
    for a, n in enumerate(nodes):
        u = ops.nodeDisp(n)
        work += sum(f[6 * a + i] * u[i] for i in range(6))
    return work / (EPS * EPS * _area(pts))


def _k(delta_deg):
    return K_ACROSS + E_S * T_S * math.cos(math.radians(delta_deg)) ** 4


def test_oracle_discriminates():
    # the along/across gap is ~6 %, far above the 1e-6 gate
    assert (K_ALONG - K_ACROSS) / K_ALONG > 0.05


@pytest.mark.parametrize('theta_local', [30.0, -65.0])
def test_local_axis_sets_bar_direction(theta_local):
    lx, ly = _unit(theta_local)
    local = (lx, ly, 0.0)
    assert membrane_stiffness(SQUARE, theta_local, local) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(SQUARE, theta_local + 90.0, local) == pytest.approx(K_ACROSS, rel=REL)
    # an off-axis direction follows cos^4 of the angle to the bar
    assert membrane_stiffness(SQUARE, 0.0, local) == pytest.approx(_k(theta_local), rel=REL)


def test_local_axis_out_of_plane_component_is_projected():
    # (1, 0, 5) projects onto the element plane z = 0 as global X
    assert membrane_stiffness(SQUARE, 0.0, (1.0, 0.0, 5.0)) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(SQUARE, 90.0, (1.0, 0.0, 5.0)) == pytest.approx(K_ACROSS, rel=REL)


def test_default_frame_is_side_12_direction_on_a_square():
    # the fork's default (no -local): bar along global X = side 1->2.
    # Upstream v3.7.1 gives the opposite pair (bar along Y), see the module docstring.
    assert membrane_stiffness(SQUARE, 0.0) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(SQUARE, 90.0) == pytest.approx(K_ACROSS, rel=REL)


def test_default_frame_follows_node_ordering():
    # same square, connectivity started one corner later: the default bar turns 90 degrees
    order = (2, 3, 4, 1)
    assert membrane_stiffness(SQUARE, 90.0, order=order) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(SQUARE, 0.0, order=order) == pytest.approx(K_ACROSS, rel=REL)
    # -local makes the bar independent of the ordering
    local = (1.0, 0.0, 0.0)
    assert membrane_stiffness(SQUARE, 0.0, local, order=order) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(SQUARE, 90.0, local, order=order) == pytest.approx(K_ACROSS, rel=REL)


def test_default_frame_follows_rotated_element():
    phi = 40.0
    pts = _rot(SQUARE, phi)
    assert membrane_stiffness(pts, phi) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(pts, phi + 90.0) == pytest.approx(K_ACROSS, rel=REL)


def test_default_frame_is_midside_line_not_side_12_on_a_trapezoid():
    # side 1->2 runs at atan(0.2) = 11.3 deg; the 14-midpoint -> 23-midpoint line at
    # atan(0.1) = 5.7 deg. The default bar follows the midside line.
    pts = [(0.0, 0.0), (1.0, 0.2), (1.0, 1.0), (0.0, 1.0)]
    theta_mid = math.degrees(math.atan2(0.1, 1.0))
    theta_12 = math.degrees(math.atan2(0.2, 1.0))
    assert membrane_stiffness(pts, theta_mid) == pytest.approx(K_ALONG, rel=REL)
    assert membrane_stiffness(pts, theta_mid + 90.0) == pytest.approx(K_ACROSS, rel=REL)
    k12 = membrane_stiffness(pts, theta_12)
    assert k12 == pytest.approx(_k(theta_12 - theta_mid), rel=REL)
    assert abs(k12 - K_ALONG) / K_ALONG > 1.0e-4      # the two candidate frames are distinguishable
