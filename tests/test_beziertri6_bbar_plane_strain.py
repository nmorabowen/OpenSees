"""WP-114 (TIMs F17): BezierTri6 ``-bbar`` is the PLANE-STRAIN (1/2) mean-
dilatation split, not the 3D (1/3) one.

Before b42ca77d8, ``BezierTri6::computeBBarMatrix`` applied the 3D split
(B_avg + 2B)/3, (B_avg - B)/3 to a 3-row plane-strain B and dropped the eps_zz
row, so the material saw an in-plane trace of (theta + 2*theta_avg)/3. Under
isochoric plastic flow (psi = 0) that keeps all 3 point-wise volumetric
constraints of the plain T6: the -bbar flag relieved nothing, the pressure
checkerboarded inside the elements and GPs were driven to the Drucker-Prager
apex; a rigid punch walled early. The fix uses the 1/2 split (LadrunoQuad's).

Gates (each fails on the pre-fix binary, passes on the fixed one):

1. kinematics -- a quadratic displacement field on one element: with -bbar,
   eps_xx + eps_yy at EVERY GP equals the element-average dilatation of the
   plain element, and eps_xx - eps_yy, gamma_xy stay local. Pre-fix the trace
   is (theta + 2 theta_avg)/3.
2. constraint rank -- nearly incompressible elastic material (nu = 0.4999999):
   the number of bulk-dominated eigenvalues of the element stiffness is 1 with
   -bbar and 3 without (pre-fix -bbar: 3).
3. FD tangent -- analytic vs central-difference element tangent at a plastic,
   non-associated (psi = 0) DP state and an associated one (< 1e-6).
4. punch -- a rigid rough punch on plane-strain DP (phi = 33, psi = 0; deck in
   tests/wp114/punch_nonassoc.py, 280 T6) reaches s/B = 0.06 with no step-size
   wall and zero tension-cutoff/apex GPs, at >= 0.75 x the LadrunoQuad -bbar
   load of the same deck (quad g1). It does NOT fully reach the quad band: at
   s/B = 0.06 T6-bbar carries 372 kPa vs 427-451 kPa for the quad meshes
   (0.82-0.87 x), and it trends down with refinement as the quad does. That
   gap is measured, not asserted away -- see 04_bezier_elements.md (WP-114).
   Pre-fix: walls at s/B = 0.039 with 63 apex GPs.

Measured wall time on the fixed binary (desktop, 2026-09-18): ~15 s total,
dominated by the punch case.
"""

import math
import os
import sys

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, 'wp114'))
import punch_nonassoc as punch   # noqa: E402

# one straight-sided T6 on the unit right triangle; CPs at the Lagrange
# positions (a straight-sided Bezier T6 with mid-edge CPs is the same space)
_NODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0),
          4: (0.5, 0.0), 5: (0.5, 0.5), 6: (0.0, 0.5)}
_CONN = (1, 2, 3, 4, 5, 6)


def _one_element(bbar, mat_cmd):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial(*mat_cmd)
    for t, (x, y) in _NODES.items():
        ops.node(t, x, y)
    extra = ['-bbar'] if bbar else []
    ops.element('BezierTri6', 1, *_CONN, 1.0, 'PlaneStrain', 1, *extra)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 5, 0)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def _set_trial(u):
    for i, t in enumerate(_NODES):
        ops.ladrunoSetNodeTrial(t, u[2 * i], u[2 * i + 1], 0.0, 0.0, 0.0, 0.0)
    ops.ladrunoTrialResidualNorm()        # Element::update()


def _quadratic_field():
    # CP values with a dilatation that varies over the element. The check only
    # compares -bbar against std on the SAME CP vector, so it need not be the
    # Bezier interpolant of a particular field.
    u = []
    for t, (x, y) in _NODES.items():
        u += [0.01 * x * x + 0.004 * y + 0.003 * x * y,
              -0.006 * y * y + 0.002 * x + 0.005 * x * y]
    return u


_ELASTIC = ('ElasticIsotropic', 1, 1000.0, 0.3)


def _gp_strains(bbar, u):
    _one_element(bbar, _ELASTIC)
    _set_trial(u)
    e = np.array(ops.eleResponse(1, 'strains')).reshape(3, 3)   # xx, yy, gxy
    return e


def test_bbar_trace_is_the_element_average_dilatation():
    u = _quadratic_field()
    std = _gp_strains(False, u)
    bb = _gp_strains(True, u)
    theta = std[:, 0] + std[:, 1]
    theta_avg = theta.mean()                     # 3 equal-weight GPs, straight sides
    assert np.ptp(theta) > 1e-3, 'field must have a non-constant dilatation'
    tr = bb[:, 0] + bb[:, 1]
    # pre-fix this is (theta + 2*theta_avg)/3, off by (theta - theta_avg)/3
    np.testing.assert_allclose(tr, theta_avg, rtol=0, atol=1e-12)
    np.testing.assert_allclose(bb[:, 0] - bb[:, 1], std[:, 0] - std[:, 1], rtol=0, atol=1e-12)
    np.testing.assert_allclose(bb[:, 2], std[:, 2], rtol=0, atol=1e-12)


def _bulk_rank(bbar):
    _one_element(bbar, ('ElasticIsotropic', 1, 1.0, 0.4999999))
    _set_trial([0.0] * 12)
    K = np.array(ops.eleResponse(1, 'stiffness')).reshape(12, 12)
    ev = np.sort(np.abs(np.linalg.eigvals(K)))[::-1]
    # G ~ 1/3; bulk ~ 1.7e6: bulk-dominated modes sit > 1e4 x the shear ones
    return int(np.sum(ev > 1e3))


def test_volumetric_constraint_rank():
    assert _bulk_rank(False) == 3     # plain T6: 3 point-wise constraints
    assert _bulk_rank(True) == 1      # 1/2 mean dilatation: one per element


# ---- FD tangent (promoted from tests/wp114/test_beziertri6_bbar_numtangent.py)

_K, _G = 27777.78, 9259.26
_PHI = math.radians(33.0)
_SP = math.sin(_PHI)
_RHO = 2 * math.sqrt(2) * _SP / (math.sqrt(3) * (3 - _SP))
_SIGY = math.sqrt(3) * 6 * 5.0 * math.cos(_PHI) / (math.sqrt(3) * (3 - _SP))
_N1, _N2 = 4, 12


def _field(x, y, t):
    lam = min(t, _N1) / _N1
    s = max(t - _N1, 0.0) / _N2
    return (-1e-3 * lam * x + 6e-3 * s * y + 3e-3 * s * x * y,
            -1e-3 * lam * y - 1.5e-3 * s * x * x)


def _tangent_rel(rho_bar):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('DruckerPrager', 1, _K, _G, _SIGY, _RHO, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    for t, (x, y) in _NODES.items():
        ops.node(t, x, y)
    ops.element('BezierTri6', 1, *_CONN, 1.0, 'PlaneStrain', 1, '-bbar')
    times = [float(i) for i in range(_N1 + _N2 + 1)]
    ts = 1
    for t, (x, y) in _NODES.items():
        for d in (0, 1):
            ops.timeSeries('Path', ts, '-time', *times,
                           '-values', *[_field(x, y, tt)[d] for tt in times])
            ops.pattern('Plain', ts, ts)
            ops.sp(t, d + 1, 1.0)
            ts += 1
    ops.constraints('Penalty', 1e13, 1e13)
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    for _ in range(_N1 + _N2 - 1):
        assert ops.analyze(1) == 0
    uc = np.array([v for t in _NODES for v in ops.nodeDisp(t)[:2]])
    t1 = ops.getTime() + 1.0
    un = np.array([v for t in _NODES for v in _field(*_NODES[t], t1)])
    us = uc + 0.5 * (un - uc)
    _set_trial(us)
    K = np.array(ops.eleResponse(1, 'stiffness')).reshape(12, 12)
    br = [int(ops.eleResponse(1, 'material', g, 'ladrunoBranch')[0]) for g in (1, 2, 3)]
    Kfd = np.zeros((12, 12))
    h = 1e-8
    for j in range(12):
        up, um = us.copy(), us.copy()
        up[j] += h
        um[j] -= h
        _set_trial(up)
        fp = np.array(ops.eleResponse(1, 'forces'))
        _set_trial(um)
        fm = np.array(ops.eleResponse(1, 'forces'))
        Kfd[:, j] = (fp - fm) / (2 * h)
    return br, np.linalg.norm(K - Kfd) / np.linalg.norm(Kfd)


@pytest.mark.parametrize('flow', ['psi0', 'assoc'])
def test_bbar_plastic_tangent_matches_fd(flow):
    br, rel = _tangent_rel(0.0 if flow == 'psi0' else _RHO)
    assert any(b != 0 for b in br), br
    assert rel < 1e-6, rel


# ---- non-associated punch

def test_nonassociated_punch_no_wall_no_apex():
    tri = punch.run(ops, 'tri6_bbar', 1, '0', smax=0.06)
    assert 'error' not in tri, tri
    assert not tri['walled'], ('walled at s/B', tri['s_over_B_last'], 'apex', tri['n_apex'])
    assert tri['s_over_B_last'] >= 0.06 - 1e-9
    # zero cutoff/apex GPs at EVERY converged step, not only the last
    assert max(n for _, n in tri['census_hist']) == 0, tri['census_hist']
    quad = punch.run(ops, 'quad', 1, '0', smax=0.06)
    assert not quad['walled']
    ratio = tri['q_last'] / quad['q_last']
    # measured 0.82 (372/451); the quad g2 mesh gives 0.87 -- see docstring
    assert 0.75 <= ratio <= 1.0, (tri['q_last'], quad['q_last'])
