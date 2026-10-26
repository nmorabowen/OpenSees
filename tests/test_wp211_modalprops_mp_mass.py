"""WP-211: modalProperties counts mass on MP-retained DOFs under Transformation.

DomainModalProperties marked a DOF fixed whenever its DOF_Group ID was -1.
The Transformation handler also gives -1 to MP-retained (slave) DOFs, so a
rigid diaphragm whose mass sits on the slave nodes had totalFreeMass = 0,
every participation factor 0, and responseSpectrumAnalysis returned a zero
response (eigen itself, and the slave eigenvectors, were already right).
Now a DOF is fixed only if it is numbered -1 AND carries an SP constraint.

Model (ported from the RSA gap study, V1 / gap 6): one storey, 4 fixed-fixed
columns (kx != ky), rigidDiaphragm 3 with the master at the plan centre, the
diaphragm mass lumped at the master / on the 4 column tops / on a 5x5 grid,
under Transformation, Penalty and Lagrange.

Oracle: the closed-form 3-DOF (ux, uy, rz) diaphragm model: periods, total
mass, centre of mass, effective masses, participation factors and the
per-mode RSA base shear, plus an explicit Transformation == Penalty ==
Lagrange cross-check. Units: kN, m, t, s.
"""
import math
import sys

import numpy as np
import pytest
from _testbed import ops

pytestmark = [pytest.mark.zone_a]

A, B, H = 12.0, 8.0, 3.5
E = 25.0e6
BX, BY = 0.40, 0.55
I_XDIR = BY * BX ** 3 / 12.0          # bending inertia resisting X sway
I_YDIR = BX * BY ** 3 / 12.0
KX = 12.0 * E * I_XDIR / H ** 3       # fixed-fixed column lateral stiffness
KY = 12.0 * E * I_YDIR / H ** 3
M_TOT = 400.0
IP_SLAB = M_TOT * (A ** 2 + B ** 2) / 12.0
G_ACC = 9.81
XI = 0.05
SDS, SD1, TL = 1.0, 0.6, 8.0
TS = SD1 / SDS
T0 = 0.2 * TS
REL = 1e-6
# a 1e14 penalty is ill-conditioned for the dense LAPACK solve (~1e-5 drift on
# the eigenvalues, as WP-203 found); Penalty is held to penalty accuracy
REL_PENALTY = 1e-4
HANDLERS = ('Transformation', 'Penalty', 'Lagrange')
KINDS = ('master', 'slaves', 'grid')
MASTER = 100


def _sa(t):
    if t < T0:
        s = SDS * (0.4 + 0.6 * t / T0)
    elif t <= TS:
        s = SDS
    elif t <= TL:
        s = SD1 / t
    else:
        s = SD1 * TL / t ** 2
    return s * G_ACC


TN = np.concatenate([np.linspace(0.0, 2.0, 201), np.linspace(2.05, 10.0, 160)]).tolist()
SA = [_sa(t) for t in TN]


def _columns(alpha):
    x, y = alpha * A / 2, alpha * B / 2
    return [(-x, -y), (x, -y), (x, y), (-x, y)]


def _mass_points(kind, alpha, ecc):
    """Point masses (x, y, m, I_rz) on the diaphragm, about the plan centre."""
    if kind == 'master':
        return [(ecc[0], ecc[1], M_TOT, IP_SLAB)]
    if kind == 'slaves':
        return [(x, y, M_TOT / 4, 0.0) for x, y in _columns(alpha)]
    xs, ys = np.linspace(-A / 2, A / 2, 5), np.linspace(-B / 2, B / 2, 5)
    return [(x + ecc[0], y + ecc[1], M_TOT / 25, 0.0) for x in xs for y in ys]


def _reference(alpha, pts):
    K = np.zeros((3, 3))
    for x, y in _columns(alpha):
        a = np.array([[1, 0, -y], [0, 1, x]], float)
        K += a.T @ np.diag([KX, KY]) @ a
    M = np.zeros((3, 3))
    for x, y, m, ip in pts:
        a = np.array([[1, 0, -y], [0, 1, x]], float)
        M += m * a.T @ a
        M[2, 2] += ip
    Li = np.linalg.inv(np.linalg.cholesky(M))
    w2, v = np.linalg.eigh(Li @ K @ Li.T)
    phi = Li.T @ v                                   # M-orthonormal
    mtot = sum(p[2] for p in pts)
    com = [sum(p[0] * p[2] for p in pts) / mtot, sum(p[1] * p[2] for p in pts) / mtot]
    out = {'w2': w2, 'phi': phi, 'mtot': mtot, 'com': com, 'Gamma': {}, 'Meff': {}}
    for d, io in (('X', np.array([1.0, 0, 0])), ('Y', np.array([0, 1.0, 0]))):
        L = phi.T @ M @ io
        out['Gamma'][d] = L                          # generalized mass is 1
        out['Meff'][d] = L ** 2
    # per-mode base shear in X and Y from the column shears
    out['V'] = {}
    for d in ('X', 'Y'):
        rows = []
        for k in range(3):
            u = out['Gamma'][d][k] * _sa(2 * math.pi / math.sqrt(w2[k])) / w2[k] * phi[:, k]
            vx = sum(KX * (u[0] - y * u[2]) for x, y in _columns(alpha))
            vy = sum(KY * (u[1] + x * u[2]) for x, y in _columns(alpha))
            rows.append({'X': vx, 'Y': vy}[d])
        out['V'][d] = rows
    return out


def _build(kind, handler, alpha, ecc=(0.0, 0.0), sp_on_slave=False):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    cols = _columns(alpha)
    for i, (x, y) in enumerate(cols, start=1):
        ops.node(i, x, y, 0.0)
        ops.fix(i, 1, 1, 1, 1, 1, 1)
        ops.node(10 + i, x, y, H)
        # sp_on_slave: an SP on top of an MP-retained DOF (node 11, uy)
        ops.fix(10 + i, 0, 1 if (sp_on_slave and i == 1) else 0, 1, 1, 1, 0)
    pts = _mass_points(kind, alpha, ecc)
    master = (ecc[0], ecc[1]) if kind == 'master' else (0.0, 0.0)
    ops.node(MASTER, master[0], master[1], H)
    ops.fix(MASTER, 0, 0, 1, 1, 1, 0)
    slaves = [11, 12, 13, 14]
    if kind == 'grid':
        for n, (x, y, m, _) in enumerate(pts, start=200):
            ops.node(n, x, y, H)
            ops.fix(n, 0, 0, 1, 1, 1, 0)
            ops.mass(n, m, m, 0, 0, 0, 0)
            slaves.append(n)
    elif kind == 'slaves':
        for i, (_, _, m, _) in enumerate(pts, start=11):
            ops.mass(i, m, m, 0, 0, 0, 0)
    else:
        _, _, m, ip = pts[0]
        ops.mass(MASTER, m, m, 0, 0, 0, ip)
    ops.rigidDiaphragm(3, MASTER, *slaves)
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)
    for i in range(1, 5):
        # J ~ 0: column torsion stays out of the 3-DOF reference
        ops.element('elasticBeamColumn', i, i, 10 + i, BX * BY, E, E / 2.4, 1.0e-10,
                    I_XDIR, I_YDIR, 1)
    if handler == 'Penalty':
        ops.constraints('Penalty', 1.0e14, 1.0e14)
    else:
        ops.constraints(handler)
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1e-10, 10)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    return pts


def _base_shear(alpha, d):
    ops.reactions()
    io = {'X': 1, 'Y': 2}[d]
    return sum(-ops.nodeReaction(i, io) for i in range(1, 5))


def _run(kind, handler, alpha, ecc=(0.0, 0.0), eigen_args=('-fullGenLapack', 3)):
    _build(kind, handler, alpha, ecc)
    # ARPACK refuses the 3-equation Transformation model: dense LAPACK
    lam = ops.eigen(*eigen_args)
    mp = ops.modalProperties('-return')
    out = {'lam': list(lam), 'mp': mp, 'V': {}}
    for d in ('X', 'Y'):
        rows = []
        for k in range(1, 4):
            ops.responseSpectrumAnalysis({'X': 1, 'Y': 2}[d], '-Tn', *TN, '-Sa', *SA,
                                         '-mode', k)
            rows.append(_base_shear(alpha, d))
        out['V'][d] = rows
    return out


def _close(a, b, rel=REL, scale=None):
    a, b = list(a), list(b)
    assert len(a) == len(b), (a, b)
    s = scale if scale is not None else max(1.0, max(abs(x) for x in b))
    for x, y in zip(a, b):
        assert abs(x - y) <= rel * s, (a, b)


CASES = [(k, h, al) for al in (1.0, 0.6) for k in KINDS for h in HANDLERS]


@pytest.mark.parametrize('kind,handler,alpha', CASES)
def test_against_closed_form(kind, handler, alpha):
    rel = REL_PENALTY if handler == 'Penalty' else REL
    ref = _reference(alpha, _mass_points(kind, alpha, (0.0, 0.0)))
    r = _run(kind, handler, alpha)
    mp = r['mp']
    _close(r['lam'], ref['w2'], rel)
    _close(mp['totalMass'][:2], [M_TOT, M_TOT], rel)
    # the bug: totalFreeMass was 0 for the slaves / grid cases under Transformation
    _close(mp['totalFreeMass'][:2], [M_TOT, M_TOT], rel)
    _close(mp['centerOfMass'][:2], ref['com'], rel, scale=A)
    for d, key, gkey in (('X', 'partiMassMX', 'partiFactorMX'),
                         ('Y', 'partiMassMY', 'partiFactorMY')):
        _close(mp[key], ref['Meff'][d], rel, scale=M_TOT)
        assert sum(mp[key]) == pytest.approx(M_TOT, rel=rel)
        # Gamma's sign follows the eigenvector's; |Gamma| is sign-free
        _close([abs(g) for g in mp[gkey]], [abs(g) for g in ref['Gamma'][d]], rel,
               scale=math.sqrt(M_TOT))
        # per-mode RSA base shear: Gamma*phi is sign-invariant, so compare signed
        _close(r['V'][d], ref['V'][d], rel, scale=max(abs(v) for v in ref['V'][d]))
        assert max(abs(v) for v in r['V'][d]) > 0.0


@pytest.mark.parametrize('kind', KINDS)
def test_transformation_equals_penalty_and_lagrange(kind):
    rows = {h: _run(kind, h, 0.6) for h in HANDLERS}
    t = rows['Transformation']
    for h in ('Penalty', 'Lagrange'):
        o = rows[h]
        rel = REL_PENALTY if h == 'Penalty' else REL
        _close(t['lam'], o['lam'], rel)
        # only UX, UY, RZ of totalFreeMass and the plan centre of mass: Penalty
        # and Lagrange number the SP-fixed DOFs, so the z centre (a geometric
        # centre here, no UZ mass) and the RX/RY gyration masses built on it
        # legitimately differ between handlers
        picks = {'totalFreeMass': [0, 1, 5], 'centerOfMass': [0, 1]}
        for key in ('totalFreeMass', 'centerOfMass', 'partiMassMX', 'partiMassMY',
                    'partiMassRatiosCumuMX', 'partiMassRatiosCumuMY'):
            idx = picks.get(key, range(len(o['mp'][key])))
            a = [t['mp'][key][i] for i in idx]
            b = [o['mp'][key][i] for i in idx]
            _close(a, b, rel, scale=max(1.0, max(abs(x) for x in b)))
        for d in ('X', 'Y'):
            _close(t['V'][d], o['V'][d], rel, scale=max(abs(v) for v in o['V'][d]))


def test_eccentric_grid_centre_of_mass():
    # mass grid shifted off the master: centre of mass and the (now coupled)
    # torsional modes must still follow the closed form under Transformation
    ecc = (0.05 * A, -0.05 * B)
    ref = _reference(0.6, _mass_points('grid', 0.6, ecc))
    r = _run('grid', 'Transformation', 0.6, ecc)
    _close(r['mp']['centerOfMass'][:2], ref['com'], scale=A)
    assert ref['com'][0] == pytest.approx(ecc[0])
    _close(r['mp']['partiMassMX'], ref['Meff']['X'], scale=M_TOT)
    _close(r['mp']['partiMassMY'], ref['Meff']['Y'], scale=M_TOT)
    for d in ('X', 'Y'):
        _close(r['V'][d], ref['V'][d], scale=max(abs(v) for v in ref['V'][d]))


def test_sp_on_mp_retained_dof_stays_fixed():
    # an SP on a slave DOF that the diaphragm also retains (node 11, uy):
    # Transformation numbers it -1 and it carries an SP, so it is fixed and
    # its mass leaves the free mass; the other three slaves keep theirs.
    _build('slaves', 'Transformation', 1.0, sp_on_slave=True)
    ops.eigen('-fullGenLapack', 2)
    mp = ops.modalProperties('-return')
    assert mp['totalMass'][1] == pytest.approx(M_TOT, rel=1e-12)
    assert mp['totalFreeMass'][1] == pytest.approx(0.75 * M_TOT, rel=1e-12)
    assert mp['totalFreeMass'][0] == pytest.approx(M_TOT, rel=1e-12)


def test_plain_handler_sp_only_unchanged():
    # no MP constraints: the 4 column tops carry their own mass, Plain handler.
    # The base SP DOFs are the only -1 DOFs; their (zero) mass stays out.
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for i, (x, y) in enumerate(_columns(1.0), start=1):
        ops.node(i, x, y, 0.0)
        ops.fix(i, 1, 1, 1, 1, 1, 1)
        ops.mass(i, 7.0, 7.0, 7.0, 0, 0, 0)        # mass on fixed DOFs: never free
        ops.node(10 + i, x, y, H)
        ops.mass(10 + i, 25.0, 25.0, 0, 0, 0, 0)
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)
    for i in range(1, 5):
        ops.element('elasticBeamColumn', i, i, 10 + i, BX * BY, E, E / 2.4, 1.0e-3,
                    I_XDIR, I_YDIR, 1)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1e-10, 10)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    ops.eigen('-fullGenLapack', 6)
    mp = ops.modalProperties('-return')
    assert mp['totalMass'][0] == pytest.approx(4 * 32.0, rel=1e-12)
    assert mp['totalFreeMass'][0] == pytest.approx(4 * 25.0, rel=1e-12)
    assert mp['totalFreeMass'][2] == pytest.approx(0.0, abs=1e-12)


def test_rsa_combine_and_increment_use_corrected_gamma():
    alpha = 0.6
    ref = _reference(alpha, _mass_points('grid', alpha, (0.0, 0.0)))
    w = np.sqrt(ref['w2'])
    # reference per-mode master ux, then CQC
    u = []
    for k in range(3):
        g = ref['Gamma']['X'][k]
        u.append(g * _sa(2 * math.pi / w[k]) / ref['w2'][k] * ref['phi'][0, k])
    rho = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            b = w[j] / w[i]
            rho[i, j] = (8 * XI ** 2 * (1 + b) * b ** 1.5) / \
                ((1 - b ** 2) ** 2 + 4 * XI ** 2 * b * (1 + b) ** 2)
    u = np.array(u)
    u_cqc = math.sqrt(u @ rho @ u)
    assert u_cqc > 0.0

    _build('grid', 'Transformation', alpha)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    # -increment (WP-207): the trial increment of mode k is the stock u_k
    for k in range(1, 4):
        ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-mode', k, '-increment')
        assert ops.nodeResponse(MASTER, 1, 4) == pytest.approx(u[k - 1], rel=REL, abs=REL * abs(u_cqc))
    ops.responseSpectrumAnalysis('-revert')
    # -combine CQC (ADR-44): the combined master ux
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-combine', 'CQC', '-damp', XI)
    assert abs(ops.nodeDisp(MASTER, 1)) == pytest.approx(u_cqc, rel=REL)


# ci-coverage: local-only -- needs FEAST (MKL, Windows/oneAPI build) (WP-143)
@pytest.mark.skipif(sys.platform != 'win32', reason='FEAST requires MKL (Windows/oneAPI build)')
def test_feast_trimmed_modes_count_slave_mass():
    # WP-203 trimDomainModes path: FEAST over all 3 modes, mass on the slaves
    ref = _reference(0.6, _mass_points('grid', 0.6, (0.0, 0.0)))
    f_max = math.sqrt(ref['w2'][-1]) / (2 * math.pi)
    r = _run('grid', 'Transformation', 0.6, eigen_args=('-feast', 0.0, 1.5 * f_max))
    assert len(r['lam']) == 3
    _close(r['lam'], ref['w2'])
    _close(r['mp']['partiMassMX'], ref['Meff']['X'], scale=M_TOT)
    assert sum(r['mp']['partiMassMY']) == pytest.approx(M_TOT, rel=REL)
    _close(r['V']['X'], ref['V']['X'], scale=max(abs(v) for v in ref['V']['X']))
