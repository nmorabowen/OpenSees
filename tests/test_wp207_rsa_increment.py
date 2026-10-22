"""WP-207: responseSpectrumAnalysis -increment (per-mode RSA about a preloaded state).

Stock solveMode() REPLACES the trial displacement with the absolute modal one
and endMode() COMMITS it, so on a gravity + loadConst domain the gravity state
is gone after the first mode: the column's axial force drops to 0 and a later
eigen loses the P-Delta softening. `-increment` imposes trial = committed +
sign*u_k, updates and records, and never commits.

Model: 2-D 3-storey cantilever column (elasticBeamColumn), gravity 500 per
storey node (N1 = -1500), loadConst, eigen, modalProperties.

Oracles:
  * the committed state is the oracle for "gravity preserved": nodeDisp
    (committed), the element axial force after -revert, and a second eigen
    must equal their pre-RSA values;
  * for a Linear transf the element forces superpose exactly:
    F(gravity + u_k) - F(gravity) == F_stock(u_k), and -sign -1 negates it;
  * for a PDelta transf the same difference carries the extra N_g*Delta/L
    chord shear that the stock (gravity-free) mode cannot have.
"""
import os

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

TN = [0.02, 0.1, 0.3, 0.6, 1.0, 3.0]
SA = [0.40, 1.00, 1.20, 0.90, 0.60, 0.20]
H = 3.0
P = 500.0
E, A, I = 2.0e8, 0.04, 1.2e-3
NMODES = 3
TOP = 4
INCR_DISP = 4  # NodeResponseType IncrDisp = trial - committed


def _setup(transf='PDelta'):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for n in (2, 3, 4):
        ops.node(n, 0.0, (n - 1) * H)
        ops.mass(n, 2.0, 2.0, 0.0)
    ops.geomTransf(transf, 1)
    for e in (1, 2, 3):
        ops.element('elasticBeamColumn', e, e, e + 1, A, E, I, 1)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n in (2, 3, 4):
        ops.load(n, 0.0, -P, 0.0)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 20)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.loadConst('-time', 0.0)
    lam = ops.eigen('-fullGenLapack', NMODES)
    ops.modalProperties()
    return lam


def _rsa(*extra):
    return ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, *extra)


def _disp():
    return [ops.nodeDisp(n, d) for n in (2, 3, 4) for d in (1, 2, 3)]


def _forces():
    return [f for e in (1, 2, 3) for f in ops.eleResponse(e, 'globalForce')]


def _close(a, b, rel=1e-9, abs_=1e-9):
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x == pytest.approx(y, rel=rel, abs=abs_), (a, b)


def test_stock_destroys_gravity():
    # documents the bug -increment exists for (and that the default is unchanged)
    _setup()
    assert abs(_forces()[1]) == pytest.approx(1500.0, rel=1e-9)
    _rsa('-mode', 1)
    assert abs(_forces()[1]) < 1e-6 * 1500.0
    assert ops.nodeDisp(TOP, 2) == pytest.approx(0.0, abs=1e-12)


def test_gravity_state_preserved():
    lam0 = _setup()
    u_g = _disp()
    f_g = _forces()
    assert abs(f_g[1]) == pytest.approx(1500.0, rel=1e-9)  # element 1 axial (global Y)
    for k in range(1, NMODES + 1):
        _rsa('-mode', k, '-increment')
        # committed state untouched; the trial state carries gravity + mode
        _close(_disp(), u_g, abs_=1e-15)
        n1 = abs(_forces()[1])
        assert abs(n1 - 1500.0) < 0.5 * 1500.0, n1
        assert abs(ops.nodeResponse(TOP, 1, INCR_DISP)) + \
            abs(ops.nodeResponse(TOP, 3, INCR_DISP)) > 1e-8
    assert ops.responseSpectrumAnalysis('-revert') == 0
    _close(_disp(), u_g, abs_=1e-15)
    _close(_forces(), f_g, rel=1e-12, abs_=1e-9)
    for n in (2, 3, 4):
        for d in (1, 2, 3):
            assert ops.nodeResponse(n, d, INCR_DISP) == 0.0
    # the P-Delta-softened eigenvalues come back unchanged
    lam1 = ops.eigen('-fullGenLapack', NMODES)
    _close(lam1, lam0, rel=1e-10, abs_=0.0)


def test_revert_flag_and_full_loop():
    lam0 = _setup()
    f_g = _forces()
    _rsa('-mode', 2, '-increment', '-revert')
    _close(_forces(), f_g, rel=1e-12, abs_=1e-9)
    # full loop: last mode left as the trial state, nothing committed
    _rsa('-increment')
    assert max(abs(a - b) for a, b in zip(_forces(), f_g)) > 1e-6
    _rsa('-increment', '-revert')
    _close(_forces(), f_g, rel=1e-12, abs_=1e-9)
    _close(ops.eigen('-fullGenLapack', NMODES), lam0, rel=1e-10, abs_=0.0)


def test_linear_forces_superpose_and_sign():
    _setup('Linear')
    f_g = _forces()
    inc, neg, disp_k = [], [], []
    for k in range(1, NMODES + 1):
        _rsa('-mode', k, '-increment')
        inc.append([a - b for a, b in zip(_forces(), f_g)])
        disp_k.append([ops.nodeResponse(n, d, INCR_DISP) for n in (2, 3, 4) for d in (1, 2, 3)])
        _rsa('-mode', k, '-increment', '-sign', -1)
        neg.append([a - b for a, b in zip(_forces(), f_g)])
    ops.responseSpectrumAnalysis('-revert')
    # stock (destructive) per-mode runs last: F_stock(u_k) and u_k
    for k in range(1, NMODES + 1):
        _rsa('-mode', k)
        stock_f = _forces()
        scale = max(abs(x) for x in stock_f)
        assert scale > 1e-6
        _close(inc[k - 1], stock_f, rel=1e-9, abs_=1e-9 * scale)
        _close(neg[k - 1], [-x for x in stock_f], rel=1e-9, abs_=1e-9 * scale)
        _close(disp_k[k - 1], _disp(), rel=1e-12, abs_=1e-15)


def test_pdelta_carries_gravity_chord_shear():
    _setup('PDelta')
    f_g = _forces()
    n_g = [f_g[6 * (e - 1) + 1] for e in (1, 2, 3)]  # axial (global Y at end i)
    inc, drift = [], []
    for k in range(1, NMODES + 1):
        _rsa('-mode', k, '-increment')
        inc.append([a - b for a, b in zip(_forces(), f_g)])
        ux = [0.0] + [ops.nodeResponse(n, 1, INCR_DISP) for n in (2, 3, 4)]
        drift.append([ux[e] - ux[e - 1] for e in (1, 2, 3)])
    ops.responseSpectrumAnalysis('-revert')
    checked = 0
    for k in range(1, NMODES + 1):
        _rsa('-mode', k)
        stock_f = _forces()
        for e in (1, 2, 3):
            extra = inc[k - 1][6 * (e - 1)] - stock_f[6 * (e - 1)]  # shear (global X)
            expect = n_g[e - 1] * drift[k - 1][e - 1] / H
            if abs(expect) > 1e-8:
                assert abs(extra) == pytest.approx(abs(expect), rel=1e-4), (k, e, extra, expect)
                checked += 1
    assert checked >= 3


def test_recorder_sees_trial_state(tmp_path):
    _setup()
    path = os.path.join(str(tmp_path), 'f.out')
    ops.recorder('Element', '-file', path, '-ele', 1, 'globalForce')
    _rsa('-mode', 1, '-increment')
    trial = ops.eleResponse(1, 'globalForce')
    ops.remove('recorders')
    rows = [ln.split() for ln in open(path) if ln.strip()]
    assert len(rows) == 1
    _close([float(x) for x in rows[-1]], trial, rel=1e-6, abs_=1e-9)


@pytest.mark.parametrize("extra", [
    ('-mode', 1, '-sign', -1),             # -sign without -increment
    ('-mode', 1, '-revert'),               # -revert without -increment
    ('-mode', 1, '-increment', '-sign', 2),
    ('-increment', '-combine', 'SRSS'),
])
def test_guards(extra):
    _setup()
    with pytest.raises(Exception):
        _rsa(*extra)
