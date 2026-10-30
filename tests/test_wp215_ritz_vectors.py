"""WP-215: `ritzVectors` (load-dependent Ritz vectors) on the LadrunoRitzKernel,
and #1025 (the kernel's K-solve check is a normwise backward error).

ADR: Ladruno_implementation/213_ladruno_modal_basis_completeness_adr.md (D7, D8).

`ritzVectors $nvec -dir d.. -pattern tag.. -cycles n..` starts from M iota_d and
the nodal loads of each pattern (in the order given), runs the WP-213 kernel
with load starts and no basis until $nvec columns survive, and installs the
Rayleigh-Ritz pairs as the domain's mode set, tagged RITZ (`# Z`). Every check
is against a pure-numpy oracle that does not repeat the kernel's algorithm
(review of #1027): it builds the per-start Krylov subspace by a Euclidean
Arnoldi and solves the generalized (K, M) Rayleigh-Ritz problem on it; the
P-Delta K comes from the model data. Review round 1 also added: residualModes
after Ritz is one joint set (Psi^T K Psi diagonal), a ritzVectors call
mid-transient leaves CentralDifference and Linear -factorOnce untouched, loads
on SP-fixed DOFs are dropped, and the energy check that Penalty rows cannot
dilute.

V3(d): `ritzVectors 3 -dir 1` on the podium model: CQC base shear closer to
the 10-mode reference than 3 eigenmodes, static load participation 100 %.
Kernel C: shapes = numpy up to sign; nvec = N reproduces eigen -fullGenLapack;
Psi^T M Psi = I. Patterns, -cycles, the V1 rigid diaphragm under the three
handlers, a PDelta preload, the Z tags, RSA after Ritz, residualModes after
Ritz, the frequency-domain warning, and the #1025 Penalty regression.
"""
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from _testbed import ops
import test_wp213_residual_modes as t213

pytestmark = [pytest.mark.zone_a]

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
N3 = t213.N3
RSA = ['-Tn', *t213.TN.tolist(), '-Sa', *t213.SA.tolist()]


# ---------------------------------------------------------------------------
# numpy oracles, independent of the kernel's algorithm (review of #1027):
# a Ritz pair depends only on the subspace, so the oracle builds the block
# Krylov subspace with a EUCLIDEAN Arnoldi (np.linalg.qr, not the kernel's
# M-Gram-Schmidt) and solves the generalized (K, M) Rayleigh-Ritz problem on
# it through a Cholesky of the projected M (not the kernel's dsyevx on an
# M-orthonormal basis).
# ---------------------------------------------------------------------------
def _span_oracle(K, M, B):
    """Ritz values (ascending) and M-orthonormal Ritz vectors of span(B)."""
    Q, _ = np.linalg.qr(np.asarray(B, float))
    L = np.linalg.cholesky(Q.T @ M @ Q)
    Li = np.linalg.inv(L)
    lam, Z = np.linalg.eigh(Li @ (Q.T @ K @ Q) @ Li.T)
    return lam, Q @ Li.T @ Z


def _krylov_oracle(K, M, F, nvec, cycles=None):
    """Load-dependent Ritz pairs of the block Krylov subspace of the starts F.

    The subspace is the per-start Krylov prefixes {(K^-1 M)^b K^-1 f_s,
    b < cap_s}, taken block by block with the starts in order and cut at nvec
    columns (ADR 213 D7). Each start's chain is orthonormalized against its
    own earlier vectors only (Euclidean), so another start's directions never
    leak into it.
    """
    cols, kept = [], [0] * len(F)
    chains = [[] for _ in F]

    def orth(v, basis):
        for _ in range(2):
            for c in basis:
                v = v - (c @ v) * c
        return v

    block = [(s, np.linalg.solve(K, f)) for s, f in enumerate(F)]
    b = 0
    while block and len(cols) < nvec:
        nxt = []
        for s, v in block:
            g = orth(v, cols)
            if not np.linalg.norm(g) > 1e-12 * np.linalg.norm(v):
                continue
            cols.append(g / np.linalg.norm(g))
            c = orth(v, chains[s])
            chains[s].append(c / np.linalg.norm(c))
            kept[s] += 1
            if len(cols) >= nvec:
                break
            cap = cycles[s] if cycles and cycles[s] > 0 else nvec
            if b + 1 < cap:
                nxt.append((s, np.linalg.solve(K, M @ chains[s][-1])))
        block, b = nxt, b + 1
    lam, Psi = _span_oracle(K, M, np.column_stack(cols))
    return lam, Psi, kept




def _static_participation(Psi, lam, K, f):
    return sum((Psi[:, k] @ f)**2 / lam[k] for k in range(len(lam))) / (f @ np.linalg.solve(K, f))


def _dynamic_participation(Psi, M, f):
    m = np.diag(M)
    return sum((Psi[:, k] @ f)**2 for k in range(Psi.shape[1])) / (f @ (f / m))


def _v3_shapes(nmodes):
    return t213._shapes(range(1, nmodes + 1), range(1, N3 + 1))


def _same_up_to_sign(A, B, M, atol=1e-9):
    for j in range(A.shape[1]):
        assert abs(A[:, j] @ M @ B[:, j]) == pytest.approx(1.0, abs=atol), j


def _child(tmp_path, lines, quiet=False):
    script = tmp_path / 'child.py'
    script.write_text('\n'.join([
        "import sys",
        f"sys.path.insert(0, {str(HERE)!r})",
        "import test_wp213_residual_modes as t",
        "import test_wp215_ritz_vectors as t5",
        "from _testbed import ops",
    ] + lines) + '\n')
    env = {k: v for k, v in os.environ.items() if k != 'LADRUNO_OPENSEES_QUIET'}
    r = subprocess.run([sys.executable, str(script)], env=env, stdin=subprocess.DEVNULL,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-3000:]
    # the child's results are 'RESULT ...' lines: the C banner may reach the
    # pipe after Python's prints
    res = [l[len('RESULT '):] for l in r.stdout.splitlines() if l.startswith('RESULT ')]
    return res, r.stderr


def _pattern_top(tag=1, nodes=(N3,), values=(1.0,)):
    ops.timeSeries('Linear', tag)
    ops.pattern('Plain', tag, tag)
    for n, v in zip(nodes, values):
        ops.load(n, v, 0.0, 0.0)


# ---------------------------------------------------------------------------
# V3(d): Ritz vs eigen at equal count, podium model
# ---------------------------------------------------------------------------
def test_v3_d_ritz_beats_eigen_at_equal_count():
    m, k, K, M = t213._v3_mk()
    io = np.ones(N3)
    lam_np, psi_np, _ = _krylov_oracle(K, M, [M @ io], 3)
    ref = t213._v3_reference()
    v_np = t213._v3_base_shear_np(psi_np, lam_np)

    t213._v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    v_eig = t213._rsa_combined_base_shear()

    t213._v3_build()
    lam = ops.ritzVectors(3, '-dir', 1)
    np.testing.assert_allclose(lam, lam_np, rtol=1e-10)
    ops.modalProperties()
    v_ritz = t213._rsa_combined_base_shear()
    assert v_ritz == pytest.approx(v_np, rel=1e-8)
    assert abs(v_ritz / ref - 1) < abs(v_eig / ref - 1), (v_ritz, v_eig, ref)
    assert abs(v_ritz / ref - 1) < 0.05, (v_ritz, ref)

    # static load participation of M iota: 100 % (K^-1 M iota is in the span)
    psi = _v3_shapes(3)
    assert _static_participation(psi, lam, K, M @ io) == pytest.approx(1.0, abs=1e-10)


def test_v3_d_report_prints_participation(tmp_path):
    _, err = _child(tmp_path, ["t._v3_build(); ops.ritzVectors(3, '-dir', 1)"])
    line = next(l for l in err.splitlines() if l.startswith('ritzVectors: dir 1:'))
    stat, dyn = (float(x) for x in re.findall(r'static ([0-9.eE+-]+)%, dynamic ([0-9.eE+-]+)%',
                                                 line)[0])
    assert stat == pytest.approx(100.0, abs=1e-6), line
    m, k, K, M = t213._v3_mk()
    lam_np, psi_np, _ = _krylov_oracle(K, M, [M @ np.ones(N3)], 3)
    assert dyn == pytest.approx(100 * _dynamic_participation(psi_np, M, M @ np.ones(N3)),
                                rel=1e-6), line
    assert 'do not report the higher ones as natural periods' in err


# ---------------------------------------------------------------------------
# kernel C
# ---------------------------------------------------------------------------
def test_kernel_matches_numpy_ldr_and_is_m_orthonormal():
    _, _, K, M = t213._v3_mk()
    lam_np, psi_np, _ = _krylov_oracle(K, M, [M @ np.ones(N3)], 6)
    t213._v3_build()
    lam = ops.ritzVectors(6, '-dir', 1)
    np.testing.assert_allclose(lam, lam_np, rtol=1e-9)
    psi = _v3_shapes(6)
    _same_up_to_sign(psi, psi_np, M)
    np.testing.assert_allclose(psi.T @ M @ psi, np.eye(6), atol=1e-10)
    np.testing.assert_allclose(psi.T @ K @ psi, np.diag(lam), rtol=1e-9, atol=1e-8)


def test_full_count_reproduces_eigen():
    t213._v3_build()
    full = ops.eigen('-fullGenLapack', N3)
    t213._v3_build()
    lam = ops.ritzVectors(N3, '-dir', 1)
    np.testing.assert_allclose(lam, full, rtol=1e-10)
    _, _, _, M = t213._v3_mk()
    psi = _v3_shapes(N3)
    np.testing.assert_allclose(psi.T @ M @ psi, np.eye(N3), atol=1e-10)


def test_more_vectors_than_mass_dofs_installs_what_survives(tmp_path):
    out, err = _child(tmp_path, [
        "t._v3_build(); print('RESULT', len(ops.ritzVectors(t.N3 + 4, '-dir', 1)))"])
    assert out == [str(N3)]
    assert f'WARNING only {N3} of {N3 + 4} vectors survived' in err


# ---------------------------------------------------------------------------
# -pattern starts and -cycles
# ---------------------------------------------------------------------------
def _p_vec(nodes, values):
    p = np.zeros(N3)
    for n, v in zip(nodes, values):
        p[n - 1] += v
    return p


def test_pattern_start_is_load_dependent():
    _, _, K, M = t213._v3_mk()
    nodes, vals = (N3, 5), (1.0, -0.4)
    p = _p_vec(nodes, vals)
    lam_np, psi_np, _ = _krylov_oracle(K, M, [p], 4)
    t213._v3_build()
    _pattern_top(1, nodes, vals)
    lam = ops.ritzVectors(4, '-pattern', 1)
    np.testing.assert_allclose(lam, lam_np, rtol=1e-9)
    psi = _v3_shapes(4)
    _same_up_to_sign(psi, psi_np, M)
    # static participation of the pattern load 100 %; the dir-1 vectors differ
    assert _static_participation(psi, lam, K, p) == pytest.approx(1.0, abs=1e-10)
    lam_d, _, _ = _krylov_oracle(K, M, [M @ np.ones(N3)], 4)
    assert np.abs(lam / lam_d - 1).max() > 1e-3


def test_pattern_report_and_order(tmp_path):
    # starts are cycled in the order given: pattern first, then dir 1
    _, _, K, M = t213._v3_mk()
    p = _p_vec((N3,), (1.0,))
    F = [p, M @ np.ones(N3)]
    lam_np, psi_np, kept = _krylov_oracle(K, M, F, 5)
    out, err = _child(tmp_path, [
        "t._v3_build(); t5._pattern_top(1)",
        "print('RESULT', *ops.ritzVectors(5, '-pattern', 1, '-dir', 1))"])
    np.testing.assert_allclose([float(v) for v in out[0].split()], lam_np, rtol=1e-9)
    lines = [l for l in err.splitlines() if re.match(r'ritzVectors: (dir|pattern) \d+:', l)]
    assert [l.split(':')[1].strip() for l in lines] == ['pattern 1', 'dir 1']
    assert lines[0].split(': ')[2].startswith(f'{kept[0]} vector(s)')
    assert lines[1].split(': ')[2].startswith(f'{kept[1]} vector(s)')
    for l in lines:
        assert re.search(r'static 100(\.0+)?%', l) or re.search(r'static 99\.99999', l), l
    dyn = float(re.findall(r'dynamic ([0-9.eE+-]+)%', lines[0])[0])
    assert dyn == pytest.approx(100 * _dynamic_participation(psi_np, M, p), rel=1e-6)


def test_cycles_cap_blocks_per_start(tmp_path):
    _, _, K, M = t213._v3_mk()
    p = _p_vec((N3,), (1.0,))
    F = [M @ np.ones(N3), p]
    lam_np, psi_np, kept = _krylov_oracle(K, M, F, 6, cycles=[2, 0])
    assert kept == [2, 4]
    t213._v3_build()
    _pattern_top(1)
    lam = ops.ritzVectors(6, '-dir', 1, '-pattern', 1, '-cycles', 2, 0)
    np.testing.assert_allclose(lam, lam_np, rtol=1e-9)
    _same_up_to_sign(_v3_shapes(6), psi_np, M)
    # without -cycles the two starts alternate (3 + 3)
    lam_free, _, kept_free = _krylov_oracle(K, M, F, 6)
    assert kept_free == [3, 3]
    assert np.abs(lam_free / lam_np - 1).max() > 1e-6
    # every start capped below nvec: fewer vectors, warned
    out, err = _child(tmp_path, ["t._v3_build(); print('RESULT', len(ops.ritzVectors(6, '-dir', 1, '-cycles', 2)))"])
    assert out == ['2']
    assert 'WARNING only 2 of 6 vectors survived' in err


def test_errors():
    t213._v3_build()
    _pattern_top(1)
    for bad in ((), (0, '-dir', 1), (3,), (3, '-dir', 4), (3, '-dir', 1, 1),
                (3, '-dir', 2),                       # no free mass in Y
                (3, '-pattern', 99), (3, '-dir', 1, '-cycles', 2, 2),
                (3, '-dir', 1, '-bogus'), (3, '-dir', 1, '-tol', 2.0)):
        with pytest.raises(Exception):
            ops.ritzVectors(*bad)


def test_no_analysis_is_refused():
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(0, 0.0)
    ops.fix(0, 1)
    ops.node(1, 1.0)
    ops.mass(1, 1.0)
    ops.uniaxialMaterial('Elastic', 1, 100.0)
    ops.element('zeroLength', 1, 0, 1, '-mat', 1, '-dir', 1)
    with pytest.raises(Exception):
        ops.ritzVectors(1, '-dir', 1)
    ops.analysis('Static')
    assert ops.ritzVectors(1, '-dir', 1) == pytest.approx([100.0], rel=1e-12)


# ---------------------------------------------------------------------------
# mid-analysis: renumber only on a domain change; -factorOnce re-armed
# ---------------------------------------------------------------------------
def _v3_transient(integrator, algorithm, n1, n2, ritz_between):
    """Base-excited V3 transient; returns the roof displacement history."""
    t213._v3_build()
    dt = 0.002
    t = np.arange(400) * dt
    acc = np.sin(2 * math.pi * 1.3 * t) * t213.G_ACC * 0.2
    ops.timeSeries('Path', 11, '-dt', dt, '-values', *acc.tolist())
    ops.pattern('UniformExcitation', 11, 1, '-accel', 11)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.algorithm(*algorithm)
    ops.integrator(*integrator)
    ops.analysis('Transient')
    hist = []
    for _ in range(n1):
        assert ops.analyze(1, dt) == 0
        hist.append(ops.nodeDisp(N3, 1))
    if ritz_between:
        ops.ritzVectors(3, '-dir', 1)
    for _ in range(n2):
        assert ops.analyze(1, dt) == 0
        hist.append(ops.nodeDisp(N3, 1))
    return np.array(hist)


@pytest.mark.parametrize('integrator,algorithm', [
    # CentralDifference::domainChanged sets U(t - dt) = U(t): an unconditional
    # renumber would cut its history (review of #1027, finding 3)
    (('CentralDifference',), ('Linear',)),
    # Linear -factorOnce would reuse the kernel's K factor in place of Keff;
    # Linear::domainChanged now re-arms it
    (('Newmark', 0.5, 0.25), ('Linear', '-factorOnce')),
])
def test_ritz_mid_transient_leaves_the_analysis_alone(integrator, algorithm):
    ref = _v3_transient(integrator, algorithm, 60, 60, False)
    got = _v3_transient(integrator, algorithm, 60, 60, True)
    assert np.abs(ref).max() > 0
    np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-14 * np.abs(ref).max())


def _chain_with_massive_base(handler):
    # 3 storeys; the fixed base node 0 carries mass, so a load on it would
    # count in the dynamic participation's f_m' M^-1 f_m if it were kept
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1, 1)
    ops.mass(0, 5.0, 0.0, 0.0)
    for i, (mi, ki) in enumerate([(3.0, 900.0), (2.0, 700.0), (1.0, 500.0)], start=1):
        ops.node(i, 0.0, 0.0)
        ops.fix(i, 0, 1, 1)
        ops.mass(i, mi, 0.0, 0.0)
        ops.uniaxialMaterial('Elastic', i, ki)
        ops.element('zeroLength', i, i - 1, i, '-mat', i, '-dir', 1)
    if handler == 'Penalty':
        ops.constraints('Penalty', 1.0e14, 1.0e14)
    else:
        ops.constraints(handler)
    ops.system('UmfPack')
    ops.analysis('Static')
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3, 1.0, 0.0, 0.0)
    ops.load(0, 7.0, 0.0, 0.0)          # on the support: dropped


def test_pattern_load_on_fixed_dof_is_dropped(tmp_path):
    # review of #1027, finding 5: the same participation under every handler
    dyn, lam = {}, {}
    for h in ('Transformation', 'Penalty', 'Lagrange'):
        out, err = _child(tmp_path, [
            f"t5._chain_with_massive_base({h!r})",
            "print('RESULT', *ops.ritzVectors(2, '-pattern', 1))"])
        lam[h] = [float(v) for v in out[0].split()]
        line = next(l for l in err.splitlines() if l.startswith('ritzVectors: pattern 1:'))
        dyn[h] = float(re.findall(r'dynamic ([0-9.eE+-]+)%', line)[0])
    for h in ('Penalty', 'Lagrange'):
        np.testing.assert_allclose(lam[h], lam['Transformation'], rtol=PEN_TOL)
        assert dyn[h] == pytest.approx(dyn['Transformation'], rel=1e-5), dyn
    # the oracle: the free 3-DOF chain, f = the roof load only
    K = np.array([[1600.0, -700, 0], [-700, 1200, -500], [0, -500, 500]])
    M = np.diag([3.0, 2.0, 1.0])
    f = np.array([0.0, 0.0, 1.0])
    lam_np, psi_np, _ = _krylov_oracle(K, M, [f], 2)
    np.testing.assert_allclose(lam['Transformation'], lam_np, rtol=1e-9)
    assert dyn['Transformation'] == pytest.approx(100 * _dynamic_participation(psi_np, M, f),
                                                  rel=1e-6)


# ---------------------------------------------------------------------------
# constraints: V1 rigid diaphragm, mass on the constrained nodes
# ---------------------------------------------------------------------------
def _v1_ritz(handler, nvec, dirs):
    t213._v1_build(handler)
    lam = list(ops.ritzVectors(nvec, '-dir', *dirs))
    mp = ops.modalProperties('-return')
    # RSA only along the started directions (a direction no vector sees is
    # refused by the WP-212 guard)
    return {'lam': lam, 'kinds': list(mp['modeKind']),
            'rx': mp['partiMassRatiosCumuMX'][-1] / 100.0,
            'ry': mp['partiMassRatiosCumuMY'][-1] / 100.0,
            'vx': t213._v1_rsa(1, lam) if 1 in dirs else 0.0,
            'vy': t213._v1_rsa(2, lam) if 2 in dirs else 0.0}


# Penalty (alpha = 1e14) round-off is LAPACK-dependent: ~4e-7 with MKL, ~2e-6
# with the bundled reference LAPACK on Linux Zone-A (WP-213)
PEN_TOL = 2e-5


@pytest.mark.parametrize('handler', ['Transformation', 'Penalty', 'Lagrange'])
def test_rigid_diaphragm_full_set_is_exact(handler):
    exp = t213._v1_expected()
    r = _v1_ritz(handler, 3, (1, 2))
    tol = PEN_TOL if handler == 'Penalty' else 1e-8
    assert r['kinds'] == [2, 2, 2]
    np.testing.assert_allclose(r['lam'], exp['lam'], rtol=tol)
    assert r['rx'] == pytest.approx(1.0, abs=1e-8)
    assert r['ry'] == pytest.approx(1.0, abs=1e-8)
    assert r['vx'] == pytest.approx(exp['v1'], rel=tol)
    assert r['vy'] == pytest.approx(exp['v2'], rel=tol)


def test_rigid_diaphragm_handlers_agree_truncated():
    # 1 vector from Y (Y couples with RZ): not an eigenvector, so this
    # compares the kernel's Krylov start itself across the handlers
    K, M = t213._v1_reference()
    lam_np, _, _ = _krylov_oracle(K, M, [M @ np.array([0.0, 1.0, 0.0])], 1)
    assert np.abs(lam_np[0] / np.linalg.eigvals(np.linalg.solve(M, K)).real - 1).min() > 1e-4
    res = {h: _v1_ritz(h, 1, (2,)) for h in ('Transformation', 'Penalty', 'Lagrange')}
    for h, r in res.items():
        tol = PEN_TOL if h == 'Penalty' else 1e-8
        np.testing.assert_allclose(r['lam'], lam_np, rtol=tol)
    for key in ('lam', 'vy'):
        np.testing.assert_allclose(res['Transformation'][key], res['Lagrange'][key], rtol=1e-8)
        np.testing.assert_allclose(res['Penalty'][key], res['Lagrange'][key], rtol=PEN_TOL)


def test_transformation_recovers_constrained_nodes():
    t213._v1_build('Transformation')
    ops.ritzVectors(3, '-dir', 1, 2)
    for md in (1, 2, 3):
        ux, uy, rz = (ops.nodeEigenvector(100, md, d) for d in (1, 2, 6))
        for i, (x, y) in enumerate(t213.COLS, start=21):
            xs, ys = x + t213.ECC, y
            assert ops.nodeEigenvector(i, md, 1) == pytest.approx(ux - ys * rz, abs=1e-12)
            assert ops.nodeEigenvector(i, md, 2) == pytest.approx(uy + xs * rz, abs=1e-12)


def test_pattern_on_constrained_node_moves_to_retained():
    # a load on a slave of the diaphragm under Transformation = the same load
    # at the retained node plus its moment (T^T p), exactly as under Lagrange
    out = {}
    for h in ('Transformation', 'Lagrange'):
        t213._v1_build(h)
        ops.timeSeries('Linear', 7)
        ops.pattern('Plain', 7, 7)
        ops.load(23, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0)
        out[h] = ops.ritzVectors(2, '-pattern', 7)
    np.testing.assert_allclose(out['Transformation'], out['Lagrange'], rtol=1e-8)


# ---------------------------------------------------------------------------
# preloaded: gravity + loadConst with PDelta -> K - Kg
# ---------------------------------------------------------------------------
def _column_k_minus_kg(P):
    """K - Kg of t213._column from the model data (review of #1027: not from
    printA). Three elastic beams along y (E, A, I, h = 3), gravity P at nodes
    2..4, so the axial forces are -3P, -2P, -P; PDelta adds N/h on the
    transverse (ux) DOFs. DOFs (ux, uy, rz) of nodes 2..4; the sign of the
    transverse axis is immaterial for the eigenvalues."""
    E, A, I, h = 2.0e8, 0.04, 1.2e-3, 3.0
    K = np.zeros((9, 9))
    for e, N in ((0, -3 * P), (1, -2 * P), (2, -P)):
        ka = E * A / h
        c = E * I / h**3
        kb = c * np.array([[12, 6 * h, -12, 6 * h],
                           [6 * h, 4 * h * h, -6 * h, 2 * h * h],
                           [-12, -6 * h, 12, -6 * h],
                           [6 * h, 2 * h * h, -6 * h, 4 * h * h]])
        kb += (N / h) * np.array([[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]])
        dofs_b = [3 * (e - 1) + 0, 3 * (e - 1) + 2, 3 * e + 0, 3 * e + 2]   # ux, rz of i, j
        dofs_a = [3 * (e - 1) + 1, 3 * e + 1]                               # uy of i, j
        for a in range(4):
            for b in range(4):
                if dofs_b[a] >= 0 and dofs_b[b] >= 0:
                    K[dofs_b[a], dofs_b[b]] += kb[a, b]
        ka2 = ka * np.array([[1, -1], [-1, 1]])
        for a in range(2):
            for b in range(2):
                if dofs_a[a] >= 0 and dofs_a[b] >= 0:
                    K[dofs_a[a], dofs_a[b]] += ka2[a, b]
    M = np.diag([2.0, 0.0, 0.0] * 3)
    io = np.array([1.0, 0.0, 0.0] * 3)
    return K, M, io


def test_preloaded_pdelta_uses_k_minus_kg():
    lam = {}
    for P in (1500.0, 0.0):
        K, M, io = _column_k_minus_kg(P)
        t213._column('PDelta', P)
        e = ops.eigen('-fullGenLapack', 3)
        t213._column('PDelta', P)
        r = ops.ritzVectors(3, '-dir', 1)      # 3 = every mass DOF (rotations massless)
        np.testing.assert_allclose(r, e, rtol=1e-10)
        # the full mass space: the exact eigenvalues of the condensed model
        a = [0, 3, 6]
        b = [i for i in range(9) if i not in a]
        Kc = K[np.ix_(a, a)] - K[np.ix_(a, b)] @ np.linalg.solve(K[np.ix_(b, b)], K[np.ix_(b, a)])
        np.testing.assert_allclose(r, np.linalg.eigvalsh(Kc / 2.0), rtol=1e-9)
        # a truncated set: the Krylov oracle on the model-data K - Kg
        lam_np, _, _ = _krylov_oracle(K, M, [M @ io], 2)
        r2 = ops.ritzVectors(2, '-dir', 1)
        np.testing.assert_allclose(r2, lam_np, rtol=1e-8)
        lam[P] = r
    assert abs(lam[1500.0][0] / lam[0.0][0] - 1) > 1e-3


# ---------------------------------------------------------------------------
# mode-kind tags, RSA, residualModes after Ritz, frequency-domain warning
# ---------------------------------------------------------------------------
def test_modal_properties_tags_z_and_eigen_clears_them(tmp_path):
    t213._v3_build()
    ops.ritzVectors(3, '-dir', 1)
    assert ops.residualModes('-kinds') == [2, 2, 2]
    assert list(ops.modalProperties('-return')['modeKind']) == [2, 2, 2]
    ops.eigen('-fullGenLapack', 3)
    assert ops.residualModes('-kinds') == [0, 0, 0]
    assert 'modeKind' not in ops.modalProperties('-return')

    fname = tmp_path / 'mp.txt'
    _child(tmp_path, ["t._v3_build(); ops.ritzVectors(3, '-dir', 1)",
                      f"ops.modalProperties('-print', '-file', {str(fname)!r})"])
    txt = fname.read_text()
    assert txt.count('  # Z\n') == 18     # 3 rows x 6 mode sections
    assert ('# NOTE: modes 1-3 (Z) are Ritz vectors; only the lowest approximate natural '
            'periods (upper bounds), the higher ones are not modes: do not report their '
            'periods as natural periods.\n') in txt


def test_rsa_after_ritz_and_require_mass():
    _, _, K, M = t213._v3_mk()
    lam_np, psi_np, _ = _krylov_oracle(K, M, [M @ np.ones(N3)], 3)
    ratio = _dynamic_participation(psi_np, M, M @ np.ones(N3))
    t213._v3_build()
    ops.ritzVectors(3, '-dir', 1)
    ops.modalProperties()
    assert t213._ratio_x() == pytest.approx(ratio, abs=1e-10)
    # -requireMass reads the Ritz set's modal mass like any mode set
    ops.responseSpectrumAnalysis(1, *RSA, '-requireMass', ratio - 1e-3)
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *RSA, '-requireMass', min(1.0, ratio + 1e-3))
    # per-mode RSA = psi Gamma Sa / lambda of the numpy vectors
    u_np = t213._v3_modal_u1_np(psi_np, lam_np)
    for k in (1, 2, 3):
        ops.responseSpectrumAnalysis(1, *RSA, '-mode', k)
        assert abs(ops.nodeDisp(1, 1)) == pytest.approx(abs(u_np[k - 1]), rel=1e-8)


def test_residual_modes_after_ritz_is_one_joint_set():
    # review of #1027, finding 1: Ritz vectors are not eigenvectors, so the
    # residual vectors couple to them through K (Psi^T K R reached -442 on
    # V3 and CQC went to 1.137x). residualModes now re-solves span[Psi, R]
    # jointly. Oracle, independent of the kernel: span[Psi, r, K^-1 M r] with
    # r = iota - Psi Psi^T M iota is span{iota, (K^-1 M)^k iota, k = 1..4}
    # (Psi = the first three Krylov powers); generalized Rayleigh-Ritz on it.
    m, k, K, M = t213._v3_mk()
    io = np.ones(N3)
    cols, v = [io], io
    for _ in range(4):
        v = np.linalg.solve(K, M @ v)
        cols.append(v / np.linalg.norm(v))
    lam_np, psi_np = _span_oracle(K, M, np.column_stack(cols))

    t213._v3_build()
    ops.ritzVectors(3, '-dir', 1)
    lam = ops.residualModes('-dir', 1)
    assert ops.residualModes('-kinds') == [2] * 5
    np.testing.assert_allclose(lam, lam_np, rtol=1e-8)
    psi = _v3_shapes(5)
    _same_up_to_sign(psi, psi_np, M, atol=1e-8)
    np.testing.assert_allclose(psi.T @ M @ psi, np.eye(5), atol=1e-10)
    kk = psi.T @ K @ psi
    off = kk - np.diag(np.diag(kk))
    assert np.abs(off).max() < 1e-8 * np.abs(np.diag(kk)).max(), off
    np.testing.assert_allclose(np.diag(kk), lam, rtol=1e-9)

    ops.modalProperties()
    assert t213._ratio_x() == pytest.approx(1.0, abs=1e-10)
    ref = t213._v3_reference()
    v = t213._rsa_combined_base_shear()
    assert v == pytest.approx(t213._v3_base_shear_np(psi_np, lam_np), rel=1e-7)
    assert abs(v / ref - 1) < 0.01, (v, ref)
    ops.responseSpectrumAnalysis(1, *RSA, '-requireMass', 1.0)


def test_residual_modes_after_ritz_refuses_a_fixed_period():
    t213._v3_build()
    ops.ritzVectors(3, '-dir', 1)
    with pytest.raises(Exception):
        ops.residualModes('-dir', 1, '-period', 0.05)
    assert ops.residualModes('-kinds') == [2, 2, 2]


def test_frequency_domain_warns_on_ritz_basis(tmp_path):
    _, err = _child(tmp_path, [
        "t._v3_build(); ops.ritzVectors(3, '-dir', 1); ops.modalProperties()",
        f"ops.frequencyResponse('-freq', 0.1, 5.0, 20, '-baseAccel', '-dir', 1, '-node', 1, '-dof', 1, '-damp', 0.05, '-out', {str(tmp_path / 'frf.out')!r})",
    ])
    assert 'frequencyResponse: WARNING the mode set holds Ritz vectors' in err, err[-2000:]
    _, err = _child(tmp_path, [
        "t._v3_build(); ops.eigen('-fullGenLapack', 3); ops.modalProperties()",
        f"ops.frequencyResponse('-freq', 0.1, 5.0, 20, '-baseAccel', '-dir', 1, '-node', 1, '-dof', 1, '-damp', 0.05, '-out', {str(tmp_path / 'frf.out')!r})",
    ])
    assert 'Ritz vectors' not in err


# ---------------------------------------------------------------------------
# #1025: the K-solve check is a normwise backward error
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('solver', [('-fullGenLapack', 1), (1,)])
def test_issue_1025_penalty_one_mode_is_accepted(solver):
    # relative residual 1.1e-4 .. 1.6e-4 under alpha = 1e14 was refused as
    # "did not solve K x = f"; the backward error is ~eps
    for _ in range(2):          # ARPACK flipped between runs
        t213._penalty_tied_floors(10, 1.0, 500.0, solver)
        assert len(ops.residualModes('-dir', 1)) >= 1
    t213._penalty_tied_floors(10, 1.0, 500.0, solver)
    assert len(ops.ritzVectors(4, '-dir', 1)) == 4


def test_soe_that_cannot_hold_k_is_still_refused(tmp_path):
    # system Diagonal keeps only diag(K): the "solve" is b / diag(K), its
    # backward error is O(0.1) and its energy mismatch O(1); both commands
    # must refuse it (the energy check runs first)
    out, err = _child(tmp_path, [
        "t._v3_build(); ops.system('Diagonal'); ops.eigen('-fullGenLapack', 3)",
        "for f in (lambda: ops.residualModes('-dir', 1), lambda: ops.ritzVectors(3, '-dir', 1)):",
        "    try:",
        "        f(); print('RESULT ACCEPTED')",
        "    except Exception:",
        "        print('RESULT REFUSED')",
    ])
    assert out == ['REFUSED', 'REFUSED'], out
    assert err.count('did not solve K x = f (energy mismatch') == 2, err[-2000:]


def test_backward_error_alone_is_diluted_by_penalty_rows():
    # review of #1027, finding 4, pinned on the tied-floor Penalty matrices:
    # a solve with the structural stiffness doubled (forward error 50 %) has a
    # normwise backward error ~1e-13 because || |K||x| || ~ alpha ||x||, so the
    # kernel also checks the energy balance |x'Kx - x'b| / x'b against
    # max(1e-4, 10 eps |x|'|K||x| / x'Kx). No stock SOE returns such a wrong x
    # on demand, so the two measures are pinned here in numpy; the kernel's
    # use of the energy check is pinned by the Diagonal refusal above.
    nf, m, k, alpha = 10, 1.0, 500.0, 1.0e14
    n = 2 * nf
    Ks = np.zeros((n, n))
    for f in range(nf):
        for side in (0, 1):
            i = 2 * f + side
            Ks[i, i] += k
            if f > 0:
                j = 2 * (f - 1) + side
                Ks[j, j] += k
                Ks[i, j] -= k
                Ks[j, i] -= k
    Kp = np.zeros((n, n))
    for f in range(nf):
        a, b = 2 * f, 2 * f + 1
        Kp[a, a] += alpha
        Kp[b, b] += alpha
        Kp[a, b] -= alpha
        Kp[b, a] -= alpha
    K = Ks + Kp
    b = np.full(n, m)
    eps = np.finfo(float).eps

    def measures(x):
        be = np.linalg.norm(K @ x - b) / (np.linalg.norm(np.abs(K) @ np.abs(x)) + np.linalg.norm(b))
        xKx = x @ K @ x
        floor = 10 * eps * (np.abs(x) @ np.abs(K) @ np.abs(x)) / abs(xKx)
        return be, abs(xKx - x @ b) / abs(x @ b), max(1e-4, floor)

    good = np.linalg.solve(K, b)
    bad = np.linalg.solve(Ks * 2 + Kp, b)
    assert np.linalg.norm(bad - good) / np.linalg.norm(good) > 0.4
    be_g, en_g, tol_g = measures(good)
    be_b, en_b, tol_b = measures(bad)
    assert be_g < 1e-4 and en_g < tol_g          # accepted
    assert be_b < 1e-4                           # backward error alone: accepted
    assert en_b > tol_b                          # the energy check refuses it


def test_partitioned_domain_refusal_is_in_place():
    src = (ROOT / 'SRC' / 'analysis' / 'analysis' / 'LadrunoRitzVectors.cpp').read_text(encoding='utf-8')
    assert 'ladrunoThreadedUpdateAllowed()' in src and 'OpenSeesSP' in src


# ---------------------------------------------------------------------------
# classic Tcl parity
# ---------------------------------------------------------------------------
def test_classic_tcl_matches_python(tmp_path):
    exe = ROOT / 'dist' / 'bin' / 'OpenSees.exe'
    if not exe.exists():
        pytest.skip('dist/bin/OpenSees.exe not built')
    m, k, _, _ = t213._v3_mk()
    lines = ['model basic -ndm 2 -ndf 3', 'node 0 0 0', 'fix 0 1 1 1']
    for i in range(1, N3 + 1):
        lines += [f'node {i} 0 0', f'fix {i} 0 1 1', f'mass {i} {m[i-1]} 0 0',
                  f'uniaxialMaterial Elastic {i} {k[i-1]}',
                  f'element zeroLength {i} {i-1} {i} -mat {i} -dir 1']
    lines += ['timeSeries Linear 1', 'pattern Plain 1 1 {', f'  load {N3} 1 0 0', '}',
              'system UmfPack', 'analysis Static',
              'set out [open res.txt w]',
              'puts $out [ritzVectors 4 -dir 1 -pattern 1 -cycles 0 1]',
              'puts $out [residualModes -kinds]', 'close $out']
    (tmp_path / 'deck.tcl').write_text('\n'.join(lines) + '\n')
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET='1')
    subprocess.run([str(exe), 'deck.tcl'], cwd=tmp_path, env=env, check=True,
                   stdin=subprocess.DEVNULL, capture_output=True, timeout=120)
    got = (tmp_path / 'res.txt').read_text().split('\n')
    assert got[1].split() == ['2', '2', '2', '2']
    t213._v3_build()
    _pattern_top(1)
    np.testing.assert_allclose([float(v) for v in got[0].split()],
                               ops.ritzVectors(4, '-dir', 1, '-pattern', 1, '-cycles', 0, 1),
                               rtol=1e-12)
