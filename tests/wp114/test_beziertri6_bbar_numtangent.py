"""WP-114 (TIMs F17), ask 2 element half: is BezierTri6 ``-bbar`` formed
correctly for plane strain under an UNSYMMETRIC (non-associated) tangent?

Two independent checks.

1. Algebra (no binary needed): the incompressibility constraint set that each
   plane-strain B / B-bar variant enforces at the 3 Gauss points of a T6.
   BezierTri6 applies the 3D split (B_avg + 2B)/3, (B_avg - B)/3 to a 3-row
   plane-strain B and drops the eps_zz row, so the in-plane trace the material
   sees is (theta + 2*theta_avg)/3.  Averaging that over the element gives
   theta_avg, so "trace = 0 at every GP" forces theta = 0 at every GP: the
   SAME kernel as the plain displacement element.  The 2D split (/2, as
   LadrunoQuad uses) leaves one constraint per element.

2. Consistency (needs the build): analytic element tangent vs a central finite
   difference of the resisting force w.r.t. nodal displacements, at a TRIAL
   state u* = u_c + du beyond a committed plastic state u_c.  Every trial is
   evaluated from the same committed material state (setTrialStrain never
   commits), so FD differentiates exactly the map whose derivative the element
   claims to return.  Elastic-state check separates formulation errors from
   material-tangent errors.

Run standalone for the report table:  python3.12 test_beziertri6_bbar_numtangent.py
"""

import math
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_DIST = os.path.normpath(os.path.join(_HERE, '..', '..', 'dist', 'bin'))
if os.path.isdir(_DIST) and _DIST not in sys.path:
    sys.path.insert(0, _DIST)


def _ops():
    try:
        import opensees as ops
    except ModuleNotFoundError:
        pytest.skip('no opensees build on sys.path')
    return ops


# ---------------------------------------------------------------------------
# 1. Algebra: constraint rank of each plane-strain B-bar variant on a T6
# ---------------------------------------------------------------------------

_GP3 = [(1 / 6, 1 / 6), (2 / 3, 1 / 6), (1 / 6, 2 / 3)]   # BezierTri6::GP3_xi
_W3 = [1 / 6, 1 / 6, 1 / 6]


def _t6_grad(xi, eta):
    """Lagrange T6 gradients on the unit right triangle (straight sides: same
    space as the Bernstein basis, so constraint ranks are basis-independent)."""
    l1, l2, l3 = 1 - xi - eta, xi, eta
    # N = [l1(2l1-1), l2(2l2-1), l3(2l3-1), 4l1l2, 4l2l3, 4l3l1]
    dl = np.array([[-1, -1], [1, 0], [0, 1]], float)
    L = [l1, l2, l3]
    g = np.zeros((6, 2))
    for a in range(3):
        g[a] = (4 * L[a] - 1) * dl[a]
    pairs = [(0, 1), (1, 2), (2, 0)]
    for k, (a, b) in enumerate(pairs):
        g[3 + k] = 4 * (L[a] * dl[b] + L[b] * dl[a])
    return g   # physical == reference for the unit triangle


def _variant_B(kind):
    grads = [_t6_grad(*p) for p in _GP3]
    gbar = sum(w * g for w, g in zip(_W3, grads)) / sum(_W3)
    Bs = []
    for g in grads:
        B = np.zeros((3, 12))
        for a in range(6):
            b1, b2 = g[a]
            c1, c2 = gbar[a]
            if kind == 'std':
                B[0, 2 * a], B[1, 2 * a + 1] = b1, b2
            elif kind == 'div3':      # BezierTri6::computeBBarMatrix
                B[0, 2 * a] = (c1 + 2 * b1) / 3; B[1, 2 * a] = (c1 - b1) / 3
                B[0, 2 * a + 1] = (c2 - b2) / 3; B[1, 2 * a + 1] = (c2 + 2 * b2) / 3
            elif kind == 'div2':      # LadrunoQuad::formB
                B[0, 2 * a] = b1 + (c1 - b1) / 2; B[1, 2 * a] = (c1 - b1) / 2
                B[0, 2 * a + 1] = (c2 - b2) / 2; B[1, 2 * a + 1] = b2 + (c2 - b2) / 2
            B[2, 2 * a], B[2, 2 * a + 1] = b2, b1
        Bs.append(B)
    return Bs


def constraint_ranks():
    out = {}
    for kind in ('std', 'div3', 'div2'):
        C = np.array([B[0] + B[1] for B in _variant_B(kind)])   # in-plane trace rows
        out[kind] = int(np.linalg.matrix_rank(C, tol=1e-10))
    return out


def test_div3_bbar_has_the_displacement_constraint_rank():
    r = constraint_ranks()
    assert r['std'] == 3
    assert r['div2'] == 1
    # the finding: the /3-drop-eps_zz variant relieves NOTHING in the
    # incompressible (psi=0 plastic flow) limit
    assert r['div3'] == r['std']


# ---------------------------------------------------------------------------
# 2. FD vs analytic element tangent
# ---------------------------------------------------------------------------

_K, _G = 27777.78, 9259.26               # kPa (same as the unsym-tangent test)
_PHI = math.radians(33.0)
_C = 5.0                                 # kPa cohesion
_SP = math.sin(_PHI)
_RHO = 2 * math.sqrt(2) * _SP / (math.sqrt(3) * (3 - _SP))          # compression cone
_SIGY = math.sqrt(3) * 6 * _C * math.cos(_PHI) / (math.sqrt(3) * (3 - _SP))
_E0 = 1.0e-3                             # isotropic pre-compression per direction
_GAM = 6.0e-3                            # final engineering shear
_BETA = 3.0e-3                           # non-homogeneous (x*y) component
_N1, _N2 = 4, 12
_PEN = 1.0e13

# unit square; tri: corners + mid-edges (Lagrange positions, straight sides)
_TRI_NODES = {1: (0, 0), 2: (1, 0), 3: (0, 1), 4: (0.5, 0), 5: (0.5, 0.5), 6: (0, 0.5)}
_TRI_CONN = (1, 2, 3, 4, 5, 6)
_QUAD_NODES = {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1)}
_QUAD_CONN = (1, 2, 3, 4)


def _field(x, y, t):
    """Prescribed nodal displacement at pseudo-time t (0..N1+N2)."""
    lam = min(t, _N1) / _N1
    s = max(t - _N1, 0.0) / _N2
    ux = -_E0 * lam * x + _GAM * s * y + _BETA * s * x * y
    uy = -_E0 * lam * y - 0.5 * _BETA * s * x * x
    return ux, uy


def _build(elem, rho_bar, sig_y):
    ops = _ops()
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('DruckerPrager', 1, _K, _G, sig_y, _RHO, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    nodes = _QUAD_NODES if elem == 'quad' else _TRI_NODES
    for tag, (x, y) in nodes.items():
        ops.node(tag, float(x), float(y))
    if elem == 'quad':
        ops.element('LadrunoQuad', 1, *_QUAD_CONN, 1, '-formulation', 'bbar',
                    '-type', 'PlaneStrain', '-thick', 1.0)
    else:
        extra = ['-bbar'] if elem == 'tri6_bbar' else []
        ops.element('BezierTri6', 1, *_TRI_CONN, 1.0, 'PlaneStrain', 1, *extra)
    times = [float(i) for i in range(_N1 + _N2 + 1)]
    ts = 1
    for tag, (x, y) in nodes.items():
        for d in (0, 1):
            vals = [_field(x, y, t)[d] for t in times]
            ops.timeSeries('Path', ts, '-time', *times, '-values', *vals)
            ops.pattern('Plain', ts, ts)
            ops.sp(tag, d + 1, 1.0)
            ts += 1
    ops.constraints('Penalty', _PEN, _PEN)
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    return ops, nodes


def _set_trial(ops, nodes, u):
    for i, tag in enumerate(nodes):
        ops.ladrunoSetNodeTrial(tag, u[2 * i], u[2 * i + 1], 0.0, 0.0, 0.0, 0.0)
    ops.ladrunoTrialResidualNorm()        # Element::update() on every element


def _force(ops):
    return np.array(ops.eleResponse(1, 'forces'))


def _stiff(ops, n):
    return np.array(ops.eleResponse(1, 'stiffness' if n == 12 else 'stiff')).reshape(n, n)


def _branches(ops, ngp):
    return [int(ops.eleResponse(1, 'material', g, 'ladrunoBranch')[0])
            for g in range(1, ngp + 1)]


def tangent_check(elem, rho_bar, sig_y, n_commit=_N1 + _N2 - 1, h=1.0e-8):
    ops, nodes = _build(elem, rho_bar, sig_y)
    for _ in range(n_commit):
        assert ops.analyze(1) == 0, (elem, rho_bar)
    n = 2 * len(nodes)
    ngp = 4 if elem == 'quad' else 3
    u_c = np.array([v for tag in nodes for v in ops.nodeDisp(tag)[:2]])
    # continue the loading path by half a step from the committed state
    t1 = ops.getTime() + 1.0
    u_next = np.array([v for tag in nodes for v in _field(*nodes[tag], t1)])
    u_star = u_c + 0.5 * (u_next - u_c)

    _set_trial(ops, nodes, u_star)
    K = _stiff(ops, n)
    br = _branches(ops, ngp)

    Kfd = np.zeros((n, n))
    for j in range(n):
        up = u_star.copy(); up[j] += h
        um = u_star.copy(); um[j] -= h
        _set_trial(ops, nodes, up); fp = _force(ops)
        _set_trial(ops, nodes, um); fm = _force(ops)
        Kfd[:, j] = (fp - fm) / (2 * h)
    _set_trial(ops, nodes, u_star)
    br_after = _branches(ops, ngp)

    nrm = np.linalg.norm(Kfd)
    sym = lambda A: 0.5 * (A + A.T)
    skw = lambda A: 0.5 * (A - A.T)
    return {
        'elem': elem, 'rho_bar': rho_bar, 'branches': br, 'branches_after': br_after,
        'rel': np.linalg.norm(K - Kfd) / nrm,
        'rel_sym': np.linalg.norm(sym(K) - sym(Kfd)) / nrm,
        'rel_skw': np.linalg.norm(skw(K) - skw(Kfd)) / nrm,
        'skew_frac': np.linalg.norm(skw(Kfd)) / nrm,
        'build': ops.ladrunoBuild(),
    }


_CASES = [(e, rb) for e in ('tri6_bbar', 'tri6_std', 'quad') for rb in ('psi0', 'assoc')]
_RB = {'psi0': 0.0, 'assoc': _RHO}


@pytest.mark.parametrize('elem,flow', _CASES)
def test_plastic_tangent_matches_fd(elem, flow):
    r = tangent_check(elem, _RB[flow], _SIGY)
    assert any(b != 0 for b in r['branches']), ('state not plastic', r)
    assert r['branches'] == r['branches_after']
    assert r['rel'] < 1.0e-5, r


@pytest.mark.parametrize('elem', ['tri6_bbar', 'tri6_std', 'quad'])
def test_elastic_tangent_matches_fd(elem):
    r = tangent_check(elem, 0.0, 1.0e6)
    assert all(b == 0 for b in r['branches']), r
    assert r['rel'] < 1.0e-6, r


if __name__ == '__main__':
    ops = _ops()
    print('ladrunoBuild:', ops.ladrunoBuild())
    print('constraint ranks (in-plane trace at 3 GPs):', constraint_ranks())
    print(f'rho={_RHO:.4f} sig_y={_SIGY:.3f}')
    hdr = f"{'case':<22}{'branches':<16}{'rel':>11}{'rel_sym':>11}{'rel_skw':>11}{'skewK':>9}"
    print(hdr)
    rows = [(f'{e}/{f}', e, _RB[f], _SIGY) for e, f in _CASES]
    rows += [(f'{e}/elastic', e, 0.0, 1.0e6) for e in ('tri6_bbar', 'tri6_std', 'quad')]
    for name, e, rb, sy in rows:
        r = tangent_check(e, rb, sy)
        print(f"{name:<22}{str(r['branches']):<16}{r['rel']:11.2e}{r['rel_sym']:11.2e}"
              f"{r['rel_skw']:11.2e}{r['skew_frac']:9.3f}")
