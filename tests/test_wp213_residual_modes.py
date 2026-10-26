"""WP-213: `residualModes` (missing-mass pseudo-modes) on the LadrunoRitzKernel.

ADR: Ladruno_implementation/213_ladruno_modal_basis_completeness_adr.md.

`residualModes -dir d ...` appends, per direction, the Rayleigh-Ritz pairs of
span{r_d, K^-1 M r_d}, r_d = iota_d deflated against the current modes, to the
domain's mode set. Every check here is against a pure-numpy oracle.

V3 podium model (ADR §5): 10-DOF shear building, DOF 1 a podium with 40x the
mass and 200x the stiffness of a typical storey; ASCE-shaped spectrum
(SDS = 1.0 g, SD1 = 0.6 g, TL = 8 s), 5 % CQC. Reference: all 10 modes.
  (a) 3 eigenmodes: the podium's 82 % of the mass is missing, base shear ~6 %;
  (b) + residualModes -dir 1: modal mass 100 %, base shear within 2 %;
  (c) -period 0.05: Sa read at 0.05 s, lambda kept: pinned to the numpy oracle.
Kernel: shapes and Ritz values match a numpy implementation up to sign;
Psi^T M Psi = I and Psi^T M Phi = 0 to 1e-10; a full-size basis reproduces
`eigen -fullGenLapack`.
Constraints: V1's one-storey rigid diaphragm, mass on the constrained nodes,
under Transformation (needs WP-211, #1020), Penalty and Lagrange.
Preloaded: gravity + loadConst with PDelta: the kernel sees K - Kg (numpy
kernel on the K the LinearSOE holds).
"""
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

ROOT = Path(__file__).resolve().parent.parent
G_ACC = 9.81
XI = 0.05
SDS, SD1, TL = 1.0, 0.6, 8.0
TS = SD1 / SDS
T0 = 0.2 * TS


def _sa(t):
    if t < T0:
        s = SDS * (0.4 + 0.6 * t / T0)
    elif t <= TS:
        s = SDS
    elif t <= TL:
        s = SD1 / t
    else:
        s = SD1 * TL / t**2
    return s * G_ACC


TN = np.concatenate([np.linspace(0.0, 2.0, 201), np.linspace(2.05, 10.0, 160)])
SA = np.array([_sa(t) for t in TN])


def _S(T):
    # the fork interpolates the -Tn/-Sa table linearly and clamps at the ends
    return float(np.interp(T, TN, SA))


def _rho(w):
    w = np.asarray(w, float)
    r = np.empty((len(w), len(w)))
    for i in range(len(w)):
        for j in range(len(w)):
            b = w[j] / w[i]
            num = 8 * XI * XI * (1 + b) * b**1.5
            den = (1 - b * b)**2 + 4 * XI * XI * b * (1 + b * b) + 8 * XI * XI * b * b
            r[i, j] = 1.0 if den <= 0 else num / den
    return r


def _cqc(per_mode, w):
    v = np.asarray(per_mode, float)
    return math.sqrt(max(v @ _rho(w) @ v, 0.0))


# ---------------------------------------------------------------------------
# numpy oracle of the kernel (shape starts, two blocks, MGS twice, RR)
# ---------------------------------------------------------------------------
def _mgs_against(y, basis, M):
    for _ in range(2):
        for q in basis:
            y = y - (q @ M @ y) * q
    return y


def _residual_oracle(K, M, Phi, iotas, tol=1e-4):
    Q = []
    for k in range(Phi.shape[1]):
        y = _mgs_against(Phi[:, k].copy(), Q, M)
        Q.append(y / math.sqrt(y @ M @ y))
    X = []
    block = [io.astype(float) for io in iotas]
    for b in range(2):
        accepted = []
        for y in block:
            n0 = math.sqrt(y @ M @ y)
            y = _mgs_against(y, Q + X, M)
            n = math.sqrt(max(y @ M @ y, 0.0))
            if not n > tol * n0:
                continue
            X.append(y / n)
            accepted.append(X[-1])
        block = [np.linalg.solve(K, M @ x) for x in accepted]
    X = np.column_stack(X)
    lam, Z = np.linalg.eigh(X.T @ K @ X)
    return lam, X @ Z


# ---------------------------------------------------------------------------
# V3 podium model
# ---------------------------------------------------------------------------
N3 = 10
KBASE = 300.0


MODEL_MK = None   # (m, k) override of the V3 storey masses / stiffnesses


def _v3_mk():
    if MODEL_MK is not None:
        m, k = (np.array(a, float) for a in MODEL_MK)
    else:
        m = np.ones(N3)
        k = np.full(N3, KBASE)
        m[0] *= 40.0
        k[0] *= 200.0
    K = np.zeros((N3, N3))
    for i in range(N3):
        K[i, i] += k[i]
        if i + 1 < N3:
            K[i, i] += k[i + 1]
            K[i, i + 1] = K[i + 1, i] = -k[i + 1]
    return m, k, K, np.diag(m)


def _v3_build():
    m, k, _, _ = _v3_mk()
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1, 1)
    for i in range(1, N3 + 1):
        ops.node(i, 0.0, 0.0)
        ops.fix(i, 0, 1, 1)
        ops.mass(i, m[i - 1], 0.0, 0.0)
        ops.uniaxialMaterial('Elastic', i, k[i - 1])
        ops.element('zeroLength', i, i - 1, i, '-mat', i, '-dir', 1)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1e-12, 10)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def _v3_modes_np():
    _, _, K, M = _v3_mk()
    w2, v = np.linalg.eigh(np.diag(np.diag(M)**-0.5) @ K @ np.diag(np.diag(M)**-0.5))
    return w2, np.diag(np.diag(M)**-0.5) @ v


def _v3_modal_u1_np(Psi, lam, sa_periods=None):
    """Per-mode u at DOF 1: psi Gamma Sa / lambda, Sa at the override if given."""
    _, _, _, M = _v3_mk()
    io = np.ones(N3)
    w = np.sqrt(lam)
    out = []
    for j in range(len(lam)):
        T = 2 * math.pi / w[j]
        if sa_periods is not None and sa_periods[j] > 0:
            T = sa_periods[j]
        g = (Psi[:, j] @ M @ io) / (Psi[:, j] @ M @ Psi[:, j])
        out.append(Psi[0, j] * g * _S(T) / lam[j])
    return np.array(out)


def _v3_base_shear_np(Psi, lam, sa_periods=None):
    # an overridden mode is one rigid mode at T: Sa and the CQC omega use T,
    # lambda stays the mode's own (ADR 213 D3b)
    _, k, _, _ = _v3_mk()
    w = np.sqrt(np.asarray(lam, float))
    if sa_periods is not None:
        for j, T in enumerate(sa_periods):
            if T > 0:
                w[j] = 2 * math.pi / T
    return _cqc(k[0] * _v3_modal_u1_np(Psi, lam, sa_periods), w)


def _v3_reference():
    w2, phi = _v3_modes_np()
    return _v3_base_shear_np(phi, w2)


def _rsa_combined_base_shear():
    m, k, _, _ = _v3_mk()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN.tolist(), '-Sa', *SA.tolist(),
                                 '-combine', 'CQC', '-damp', XI)
    return k[0] * abs(ops.nodeDisp(1, 1))


def _shapes(modes, nodes, dof=1):
    return np.array([[ops.nodeEigenvector(n, md, dof) for md in modes] for n in nodes])


def _ratio_x():
    mp = ops.modalProperties('-return')
    return mp['partiMassRatiosCumuMX'][-1] / 100.0   # modalProperties reports %


def test_v3_a_truncated_basis_misses_the_podium():
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    assert _ratio_x() < 0.25
    ref = _v3_reference()
    v = _rsa_combined_base_shear()
    assert v / ref < 0.1, (v, ref)


def test_v3_b_residual_modes_complete_the_mass_and_base_shear():
    _v3_build()
    lam_e = ops.eigen('-fullGenLapack', 3)
    lam_r = ops.residualModes('-dir', 1)
    assert len(lam_r) == 2
    assert ops.residualModes('-kinds') == [0, 0, 0, 1, 1]
    ops.modalProperties()
    assert _ratio_x() == pytest.approx(1.0, abs=1e-10)

    # pseudo-modes = numpy oracle (values exactly, shapes up to sign)
    _, _, K, M = _v3_mk()
    w2, phi = _v3_modes_np()
    lam_np, psi_np = _residual_oracle(K, M, phi[:, :3], [np.ones(N3)])
    np.testing.assert_allclose(lam_r, lam_np, rtol=1e-9)
    psi = _shapes([4, 5], range(1, N3 + 1))
    for j in range(2):
        assert abs(psi[:, j] @ M @ psi_np[:, j]) == pytest.approx(1.0, abs=1e-9)

    # M-orthonormal, M-orthogonal to the eigenmodes, K-diagonal
    Phi = _shapes([1, 2, 3], range(1, N3 + 1))
    np.testing.assert_allclose(psi.T @ M @ psi, np.eye(2), atol=1e-10)
    gm = np.sqrt(np.diag(Phi.T @ M @ Phi))
    np.testing.assert_allclose((Phi / gm).T @ M @ psi, 0.0, atol=1e-10)
    np.testing.assert_allclose(psi.T @ K @ psi, np.diag(lam_r), rtol=1e-9, atol=1e-9)

    ref = _v3_reference()
    v = _rsa_combined_base_shear()
    assert v == pytest.approx(_v3_base_shear_np(np.column_stack([phi[:, :3], psi_np]),
                                                np.r_[w2[:3], lam_np]), rel=1e-8)
    assert abs(v / ref - 1.0) < 0.02, (v, ref)
    assert len(lam_e) == 3


def test_v3_c_fixed_period_vs_rayleigh():
    # -period T keeps the Rayleigh lambda and reads Sa(T) (ADR 213 D3b): the
    # pseudo-modes then add their static response times Sa(0.05). On V3 the
    # Rayleigh periods (0.31 s, 0.16 s) sit on the plateau and Sa(0.05) = 0.65 g
    # is on the ascending branch, so the fixed-period result lands BELOW
    # rayleigh; it is pinned to the oracle, not forced above the reference.
    _, _, K, M = _v3_mk()
    w2, phi = _v3_modes_np()
    lam_np, psi_np = _residual_oracle(K, M, phi[:, :3], [np.ones(N3)])
    P = np.column_stack([phi[:, :3], psi_np])
    lam_all = np.r_[w2[:3], lam_np]
    ref = _v3_reference()

    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.residualModes('-dir', 1)
    ops.modalProperties()
    v_ray = _rsa_combined_base_shear()

    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    lam_fix = ops.residualModes('-dir', 1, '-period', 0.05)
    np.testing.assert_allclose(lam_fix, lam_np, rtol=1e-9)   # lambda unchanged
    ops.modalProperties()
    assert _ratio_x() == pytest.approx(1.0, abs=1e-10)
    v_fix = _rsa_combined_base_shear()
    assert v_fix == pytest.approx(_v3_base_shear_np(P, lam_all, [0, 0, 0, 0.05, 0.05]),
                                  rel=1e-8)
    assert v_fix < v_ray
    print(f"\nV3(c): rayleigh {v_ray / ref:.4f} x ref, -period 0.05 {v_fix / ref:.4f} x ref")


def test_v3_c_rigid_residual_period_override_is_conservative(monkeypatch):
    # 100x stiffer storeys: the pseudo-modes' Rayleigh periods (0.031 s,
    # 0.016 s) are below 0.05 s, the case ASCE 7-22 12.9.1.1 means. On the
    # ascending branch Sa(0.05) >= Sa(T_R), so -period 0.05 is conservative.
    monkeypatch.setattr(sys.modules[__name__], 'KBASE', 3.0e4)
    _, _, K, M = _v3_mk()
    w2, phi = _v3_modes_np()
    lam_np, psi_np = _residual_oracle(K, M, phi[:, :3], [np.ones(N3)])
    assert (2 * math.pi / np.sqrt(lam_np) < 0.05).all()
    ref = _v3_reference()
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.residualModes('-dir', 1)
    ops.modalProperties()
    v_ray = _rsa_combined_base_shear()
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.residualModes('-dir', 1, '-period', 0.05)
    ops.modalProperties()
    v_fix = _rsa_combined_base_shear()
    assert v_fix == pytest.approx(
        _v3_base_shear_np(np.column_stack([phi[:, :3], psi_np]), np.r_[w2[:3], lam_np],
                          [0, 0, 0, 0.05, 0.05]), rel=1e-8)
    assert abs(v_ray / ref - 1) < 0.01
    assert v_fix >= v_ray >= ref * 0.99
    print(f"\nV3(c) stiff: rayleigh {v_ray / ref:.4f} x ref, -period 0.05 {v_fix / ref:.4f} x ref")


def test_period_override_correlates_the_pseudo_modes(monkeypatch):
    # Review repro (#1023, Opus round 1): a stiff two-storey podium. With
    # their own omegas the two pseudo-modes (0.026 s, 0.010 s) combine almost
    # as SRSS (rho ~ 0.02): 1.278 x ref. As one rigid mode at 0.05 s they are
    # fully correlated: 1.651 x ref.
    monkeypatch.setattr(sys.modules[__name__], 'MODEL_MK',
                        ([40, 40] + [1.0] * 8, [3e4 * 400, 3e4 * 100] + [3e4] * 8))
    _, _, K, M = _v3_mk()
    w2, phi = _v3_modes_np()
    lam_np, psi_np = _residual_oracle(K, M, phi[:, :3], [np.ones(N3)])
    P = np.column_stack([phi[:, :3], psi_np])
    lam_all = np.r_[w2[:3], lam_np]
    ref = _v3_reference()
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.residualModes('-dir', 1, '-period', 0.05)
    ops.modalProperties()
    v = _rsa_combined_base_shear()
    assert v == pytest.approx(_v3_base_shear_np(P, lam_all, [0, 0, 0, 0.05, 0.05]), rel=1e-8)
    assert v / ref == pytest.approx(1.651, abs=1e-3)


def test_period_override_on_every_rsa_path():
    _, _, K, M = _v3_mk()
    w2, phi = _v3_modes_np()
    lam_np, psi_np = _residual_oracle(K, M, phi[:, :3], [np.ones(N3)])
    P = np.column_stack([phi[:, :3], psi_np])
    lam_all = np.r_[w2[:3], lam_np]
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.residualModes('-dir', 1, '-period', 0.05)
    ops.modalProperties()
    u_np = _v3_modal_u1_np(P, lam_all, [0, 0, 0, 0.05, 0.05])
    for k in range(1, 6):
        # per mode
        ops.responseSpectrumAnalysis(1, '-Tn', *TN.tolist(), '-Sa', *SA.tolist(), '-mode', k)
        assert abs(ops.nodeDisp(1, 1)) == pytest.approx(abs(u_np[k - 1]), rel=1e-8)
    # -increment: trial - committed = u_k (the committed state is zero here)
    for k in (4, 5):
        ops.responseSpectrumAnalysis(1, '-Tn', *TN.tolist(), '-Sa', *SA.tolist(),
                                     '-mode', k, '-increment')
        assert abs(ops.nodeResponse(1, 1, 4)) == pytest.approx(abs(u_np[k - 1]), rel=1e-8)
        ops.responseSpectrumAnalysis('-revert')
    # without -period the same modes read Sa at their own periods
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.residualModes('-dir', 1)
    ops.modalProperties()
    u_ray = _v3_modal_u1_np(P, lam_all)
    ops.responseSpectrumAnalysis(1, '-Tn', *TN.tolist(), '-Sa', *SA.tolist(), '-mode', 4)
    assert abs(ops.nodeDisp(1, 1)) == pytest.approx(abs(u_ray[3]), rel=1e-8)


def test_wp212_mass_guard_sees_the_residual_modes():
    rsa = ['-Tn', *TN.tolist(), '-Sa', *SA.tolist()]
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *rsa, '-requireMass', 1.0)
    ops.residualModes('-dir', 1)
    ops.modalProperties()
    ops.responseSpectrumAnalysis(1, *rsa, '-requireMass', 1.0)
    ops.responseSpectrumAnalysis(1, *rsa, '-requireMass', 1.0, '-mode', 4)
    ops.responseSpectrumAnalysis(1, *rsa, '-requireMass', 1.0, '-combine', 'CQC', '-damp', XI)
    ops.responseSpectrumAnalysis(1, *rsa, '-requireMass', 1.0, '-mode', 5, '-increment')


def test_wp212_report_prints_100_percent(tmp_path):
    # the report line is suppressed under LADRUNO_OPENSEES_QUIET: run a child
    script = tmp_path / 'child.py'
    script.write_text('\n'.join([
        "import sys",
        f"sys.path.insert(0, {str(Path(__file__).parent)!r})",
        "import test_wp213_residual_modes as t",
        "from _testbed import ops",
        "t._v3_build(); ops.eigen('-fullGenLapack', 3); ops.residualModes('-dir', 1)",
        "ops.modalProperties()",
        "ops.responseSpectrumAnalysis(1, '-Tn', *t.TN.tolist(), '-Sa', *t.SA.tolist(), '-mode', 4)",
    ]) + '\n')
    env = {k: v for k, v in os.environ.items() if k != 'LADRUNO_OPENSEES_QUIET'}
    r = subprocess.run([sys.executable, str(script)], env=env, stdin=subprocess.DEVNULL,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    out = r.stdout + r.stderr
    assert 'cumulative effective mass in direction 1 = 100.00% of the free mass (5 modes)' in out, out[-2000:]


def test_replace_and_eigen_clears_tags():
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    lam1 = ops.residualModes('-dir', 1)
    # a second call without -replace: everything is already spanned, no change
    assert not ops.residualModes('-dir', 1)
    assert ops.residualModes('-kinds') == [0, 0, 0, 1, 1]
    lam2 = ops.residualModes('-dir', 1, '-replace')
    np.testing.assert_allclose(lam2, lam1, rtol=1e-9)
    assert ops.residualModes('-kinds') == [0, 0, 0, 1, 1]
    # a new eigen replaces the eigenvalue vector, so the record goes stale
    ops.eigen('-fullGenLapack', 4)
    assert ops.residualModes('-kinds') == [0, 0, 0, 0]


def test_period_rule_mismatch_needs_replace():
    # review round 1, item 2: a second call with another period rule used to
    # find nothing new to add and drop the -period silently
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    lam = ops.residualModes('-dir', 1)
    with pytest.raises(Exception):
        ops.residualModes('-dir', 1, '-period', 0.05)
    assert ops.residualModes('-kinds') == [0, 0, 0, 1, 1]
    lam2 = ops.residualModes('-dir', 1, '-period', 0.05, '-replace')
    np.testing.assert_allclose(lam2, lam, rtol=1e-9)
    ops.modalProperties()
    _, _, K, M = _v3_mk()
    w2, phi = _v3_modes_np()
    lam_np, psi_np = _residual_oracle(K, M, phi[:, :3], [np.ones(N3)])
    v = _rsa_combined_base_shear()
    assert v == pytest.approx(_v3_base_shear_np(np.column_stack([phi[:, :3], psi_np]),
                                                np.r_[w2[:3], lam_np], [0, 0, 0, 0.05, 0.05]),
                              rel=1e-8)
    # and back to Rayleigh periods the same way
    with pytest.raises(Exception):
        ops.residualModes('-dir', 1)
    ops.residualModes('-dir', 1, '-replace')


def test_state_changed_after_eigen_is_refused():
    # review round 1, item 3: eigen on the unloaded column, then gravity with
    # P-Delta: Phi is no longer an eigenbasis of the current K
    _column('PDelta', 0.0)
    ops.eigen('-fullGenLapack', 1)
    ops.timeSeries('Constant', 2)
    ops.pattern('Plain', 2, 2)
    for n in (2, 3, 4):
        ops.load(n, 0.0, -1500.0, 0.0)
    assert ops.analyze(1) == 0
    ops.loadConst('-time', 0.0)
    with pytest.raises(Exception):
        ops.residualModes('-dir', 1)
    ops.eigen('-fullGenLapack', 1)       # re-run on the current state: accepted
    assert len(ops.residualModes('-dir', 1)) == 2


def test_replace_rebuilds_every_replaced_direction():
    # review round 1, item 5: -dir 2, then -dir 1 -replace used to drop the
    # Y pseudo-modes and leave Y below 100 %
    _v1_build('Lagrange')
    ops.eigen('-fullGenLapack', 1)       # mode 1 = X, which it carries completely
    ops.residualModes('-dir', 2)
    assert ops.residualModes('-kinds') == [0, 1, 1]
    ops.residualModes('-dir', 1, '-replace')
    assert ops.residualModes('-kinds') == [0, 1, 1]
    mp = ops.modalProperties('-return')
    assert mp['partiMassRatiosCumuMY'][-1] / 100.0 == pytest.approx(1.0, abs=1e-8)


def test_noise_level_residual_is_warned(tmp_path):
    # review round 1, item 4: with -tol 1e-6 the X residual of V1 under
    # Penalty (finite-penalty coupling noise, 1.3e-11 of the mass; exact under
    # Lagrange / Transformation) is kept; it must be flagged
    script = tmp_path / 'child.py'
    script.write_text('\n'.join([
        "import sys",
        f"sys.path.insert(0, {str(Path(__file__).parent)!r})",
        "import test_wp213_residual_modes as t",
        "from _testbed import ops",
        "t._v1_build('Penalty'); ops.eigen('-fullGenLapack', 1)",
        "ops.residualModes('-dir', 1, 2, '-tol', 1e-6)",
    ]) + '\n')
    r = subprocess.run([sys.executable, str(script)], env=dict(os.environ),
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    out = r.stdout + r.stderr
    assert 'WARNING direction 1' in out and 'eigenvector noise' in out, out[-2000:]


def test_full_size_basis_equals_full_eigen():
    # 3 storeys: 1 eigenmode + span{r, K^-1 M r} is the whole space, so the
    # Ritz values are the exact eigenvalues 2 and 3
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1, 1)
    for i, (m, k) in enumerate([(3.0, 900.0), (2.0, 700.0), (1.0, 500.0)], start=1):
        ops.node(i, 0.0, 0.0)
        ops.fix(i, 0, 1, 1)
        ops.mass(i, m, 0.0, 0.0)
        ops.uniaxialMaterial('Elastic', i, k)
        ops.element('zeroLength', i, i - 1, i, '-mat', i, '-dir', 1)
    ops.system('UmfPack')
    ops.analysis('Static')
    full = ops.eigen('-fullGenLapack', 3)
    ops.eigen('-fullGenLapack', 1)
    lam = ops.residualModes('-dir', 1)
    np.testing.assert_allclose(lam, full[1:], rtol=1e-10)


def test_errors():
    _v3_build()
    with pytest.raises(Exception):
        ops.residualModes('-dir', 1)          # no eigen yet
    ops.eigen('-fullGenLapack', 3)
    for bad in (('-dir', 4), ('-dir',), ('-dir', 1, '-period', -1.0),
                ('-dir', 1, 1), ('-bogus',), ('-dir', 2)):  # dir 2: no free mass
        with pytest.raises(Exception):
            ops.residualModes(*bad)


def test_modal_response_history_uses_the_pseudo_modes(tmp_path):
    dt, n = 0.005, 800
    t = np.arange(n) * dt
    acc = (np.sin(2 * math.pi * 0.6 * t) + 0.8 * np.sin(2 * math.pi * 1.5 * t)) * G_ACC * 0.3
    m, k, _, _ = _v3_mk()

    def peak(nmodes, residual):
        _v3_build()
        ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
        ops.eigen('-fullGenLapack', nmodes)
        if residual:
            ops.residualModes('-dir', 1)
        ops.modalProperties()
        path = tmp_path / f'u_{nmodes}_{residual}.out'
        ops.recorder('Node', '-file', str(path), '-time', '-precision', 16,
                     '-node', 1, '-dof', 1, 'disp')
        ops.modalResponseHistory('-dt', dt, '-nsteps', n - 1, '-damp', XI,
                                 '-baseAccel', 1, '-dir', 1)
        ops.remove('recorders')
        return k[0] * np.abs(np.loadtxt(str(path))[:, 1]).max()

    ref = peak(N3, False)
    e3 = abs(peak(3, False) / ref - 1)
    e3r = abs(peak(3, True) / ref - 1)
    assert e3 > 0.5
    assert e3r < 0.1 and e3r < 0.2 * e3, (e3, e3r)


# ---------------------------------------------------------------------------
# constraints: V1 one-storey rigid diaphragm, mass on the constrained nodes
# ---------------------------------------------------------------------------
A_, B_, H_ = 12.0, 8.0, 3.5
E_ = 25.0e6
BX, BY = 0.40, 0.55
I_X = BY * BX**3 / 12.0
I_Y = BX * BY**3 / 12.0
KX = 12.0 * E_ * I_X / H_**3
KY = 12.0 * E_ * I_Y / H_**3
M_TOT = 400.0
ECC = 0.05 * A_
COLS = [(-A_ / 2, -B_ / 2), (A_ / 2, -B_ / 2), (A_ / 2, B_ / 2), (-A_ / 2, B_ / 2)]


def _v1_build(handler):
    # columns at COLS; the slab mass sits on four extra constrained nodes
    # offset by ECC in x (mass eccentricity -> Y and RZ couple, X does not)
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    slaves = []
    for i, (x, y) in enumerate(COLS, start=1):
        ops.node(i, x, y, 0.0)
        ops.fix(i, 1, 1, 1, 1, 1, 1)
        ops.node(10 + i, x, y, H_)
        ops.fix(10 + i, 0, 0, 1, 1, 1, 0)
        ops.node(20 + i, x + ECC, y, H_)
        ops.fix(20 + i, 0, 0, 1, 1, 1, 0)
        ops.mass(20 + i, M_TOT / 4, M_TOT / 4, 0, 0, 0, 0)
        slaves += [10 + i, 20 + i]
    ops.node(100, 0.0, 0.0, H_)
    ops.fix(100, 0, 0, 1, 1, 1, 0)
    ops.rigidDiaphragm(3, 100, *slaves)
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)
    for i in range(1, 5):
        # J ~ 0: column torsion is left out of the 3-DOF reference
        ops.element('elasticBeamColumn', i, i, 10 + i, BX * BY, E_, E_ / 2.4, 1.0e-10,
                    I_X, I_Y, 1)
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


def _v1_reference():
    # 3-DOF (ux, uy, rz) at the diaphragm centre
    K = np.zeros((3, 3))
    M = np.zeros((3, 3))
    for x, y in COLS:
        a = np.array([[1, 0, -y], [0, 1, x]], float)
        K += a.T @ np.diag([KX, KY]) @ a
        a = np.array([[1, 0, -y], [0, 1, x + ECC]], float)
        M += (M_TOT / 4) * a.T @ a
    return K, M


def _v1_rsa(direction, lams):
    rows = []
    for k in range(1, len(lams) + 1):
        ops.responseSpectrumAnalysis(direction, '-Tn', *TN.tolist(), '-Sa', *SA.tolist(),
                                     '-mode', k)
        ops.reactions()
        rows.append(-sum(ops.nodeReaction(i, direction) for i in range(1, 5)))
    return _cqc(rows, np.sqrt(lams))


def _v1_run(handler):
    _v1_build(handler)
    lam_e = ops.eigen('-fullGenLapack', 1)
    lam_r = ops.residualModes('-dir', 1, 2)
    mp = ops.modalProperties('-return')
    lams = list(lam_e) + list(lam_r)
    return {'lam': lams, 'kinds': ops.residualModes('-kinds'),
            'rx': mp['partiMassRatiosCumuMX'][-1] / 100.0,
            'ry': mp['partiMassRatiosCumuMY'][-1] / 100.0,
            'vx': _v1_rsa(1, lams), 'vy': _v1_rsa(2, lams)}


def _v1_expected():
    K, M = _v1_reference()
    L = np.linalg.cholesky(M)
    Li = np.linalg.inv(L)
    w2, v = np.linalg.eigh(Li @ K @ Li.T)
    phi = Li.T @ v
    out = {'lam': w2}
    for d, io in ((1, np.array([1.0, 0, 0])), (2, np.array([0, 1.0, 0]))):
        V = []
        for j in range(3):
            g = (phi[:, j] @ M @ io) / (phi[:, j] @ M @ phi[:, j])
            u = phi[:, j] * g * _S(2 * math.pi / math.sqrt(w2[j])) / w2[j]
            V.append((K @ u)[d - 1])
        out[f'v{d}'] = _cqc(V, np.sqrt(w2))
    return out


@pytest.mark.parametrize('handler', [
    'Penalty',
    'Lagrange',
    'Transformation',   # WP-211 (#1020) keeps the MP-retained mass
])
def test_rigid_diaphragm_mass_on_constrained_nodes(handler):
    exp = _v1_expected()
    r = _v1_run(handler)
    # X is uncoupled (mode 1 carries all of it); Y-RZ is completed exactly by
    # span{r_y, K^-1 M r_y}
    # Penalty (alpha = 1e14) is ~4e-7 off the exact eigenvalues already in
    # eigen's own mode 1; Lagrange is exact
    tol = 2e-6 if handler == 'Penalty' else 1e-8
    assert r['kinds'] == [0, 1, 1]
    np.testing.assert_allclose(r['lam'], exp['lam'], rtol=tol)
    assert r['rx'] == pytest.approx(1.0, abs=1e-8)
    assert r['ry'] == pytest.approx(1.0, abs=1e-8)
    assert r['vx'] == pytest.approx(exp['v1'], rel=tol)
    assert r['vy'] == pytest.approx(exp['v2'], rel=tol)


def test_penalty_and_lagrange_agree():
    p, l = _v1_run('Penalty'), _v1_run('Lagrange')
    for key in ('lam', 'vx', 'vy'):
        np.testing.assert_allclose(p[key], l[key], rtol=2e-6)   # penalty-limited


def test_transformation_recovers_constrained_nodes():
    # passes without WP-211: the pseudo-modes go in through
    # TransformationDOF_Group::setEigenvector, which rebuilds the slave nodes
    _v1_build('Transformation')
    ops.eigen('-fullGenLapack', 1)
    lam_r = ops.residualModes('-dir', 1, 2)
    assert len(lam_r) == 2
    exp = _v1_expected()
    np.testing.assert_allclose(lam_r, exp['lam'][1:], rtol=1e-8)
    m = M_TOT / 4
    for md in (2, 3):
        ux, uy, rz = (ops.nodeEigenvector(100, md, d) for d in (1, 2, 6))
        for i, (x, y) in enumerate(COLS, start=21):
            xs, ys = x + ECC, y
            assert ops.nodeEigenvector(i, md, 1) == pytest.approx(ux - ys * rz, abs=1e-12)
            assert ops.nodeEigenvector(i, md, 2) == pytest.approx(uy + xs * rz, abs=1e-12)
    # M-orthonormal on the slave masses
    P = np.array([[ops.nodeEigenvector(i, md, d) for md in (2, 3)]
                  for i in range(21, 25) for d in (1, 2)])
    np.testing.assert_allclose(m * P.T @ P, np.eye(2), atol=1e-10)


# ---------------------------------------------------------------------------
# preloaded: gravity + loadConst + PDelta -> the kernel sees K - Kg
# ---------------------------------------------------------------------------
def _column(transf, gravity):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for n in (2, 3, 4):
        ops.node(n, 0.0, (n - 1) * 3.0)
        ops.mass(n, 2.0, 0.0, 0.0)
    ops.geomTransf(transf, 1)
    for e in (1, 2, 3):
        ops.element('elasticBeamColumn', e, e, e + 1, 0.04, 2.0e8, 1.2e-3, 1)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n in (2, 3, 4):
        ops.load(n, 0.0, -gravity, 0.0)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 20)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.loadConst('-time', 0.0)


def test_preloaded_pdelta_kernel_sees_geometric_stiffness():
    # Oracle: the numpy kernel on the K the analysis' LinearSOE holds after
    # residualModes (printA), with M and the eigenmode read back in equation
    # numbering. The column's rotations are massless, so span{r, K^-1 M r} is
    # not the exact eigen-subspace here (r is not statically condensed on the
    # massless DOFs; ADR 213 §3): the check is the kernel itself, plus that the
    # preload changes what it sees.
    out = {}
    for P in (1500.0, 0.0):
        _column('PDelta', P)
        ops.eigen('-fullGenLapack', 1)
        lam = ops.residualModes('-dir', 1)
        n = ops.systemSize()
        K = np.array(ops.printA('-ret'), dtype=float).reshape(n, n)
        M = np.zeros((n, n))
        phi = np.zeros(n)
        io = np.zeros(n)
        for nd in (2, 3, 4):
            eq = ops.nodeDOFs(nd)
            M[eq[0], eq[0]] = 2.0
            io[eq[0]] = 1.0
            for d in range(3):
                phi[eq[d]] = ops.nodeEigenvector(nd, 1, d + 1)
        lam_np, psi_np = _residual_oracle(0.5 * (K + K.T), M, phi[:, None], [io])
        np.testing.assert_allclose(lam, lam_np, rtol=1e-8)
        psi = np.zeros((n, 2))
        for nd in (2, 3, 4):
            eq = ops.nodeDOFs(nd)
            for d in range(3):
                psi[eq[d]] = [ops.nodeEigenvector(nd, md, d + 1) for md in (2, 3)]
        for j in range(2):
            assert abs(psi[:, j] @ M @ psi_np[:, j]) == pytest.approx(1.0, abs=1e-8)
        # modal mass of the set is complete even with massless DOFs
        L = [(v @ M @ io)**2 / (v @ M @ v) for v in np.column_stack([phi, psi]).T]
        assert sum(L) == pytest.approx(6.0, rel=1e-10)
        out[P] = (lam, K)
    assert np.abs(out[1500.0][1] - out[0.0][1]).max() > 1.0   # K - Kg vs K
    assert abs(out[1500.0][0][0] / out[0.0][0][0] - 1) > 1e-3


# ---------------------------------------------------------------------------
# classic Tcl parity (OpenSees.exe goes through SRC/tcl/commands.cpp)
# ---------------------------------------------------------------------------
def test_classic_tcl_matches_python(tmp_path):
    exe = ROOT / 'dist' / 'bin' / 'OpenSees.exe'
    if not exe.exists():
        pytest.skip('dist/bin/OpenSees.exe not built')
    m, k, _, _ = _v3_mk()
    lines = ['model basic -ndm 2 -ndf 3', 'node 0 0 0', 'fix 0 1 1 1']
    for i in range(1, N3 + 1):
        lines += [f'node {i} 0 0', f'fix {i} 0 1 1', f'mass {i} {m[i-1]} 0 0',
                  f'uniaxialMaterial Elastic {i} {k[i-1]}',
                  f'element zeroLength {i} {i-1} {i} -mat {i} -dir 1']
    lines += ['system UmfPack', 'analysis Static', 'eigen -fullGenLapack 3',
              'set out [open res.txt w]',
              'puts $out [residualModes -dir 1]',
              'puts $out [residualModes -kinds]', 'close $out']
    (tmp_path / 'deck.tcl').write_text('\n'.join(lines) + '\n')
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET='1')
    subprocess.run([str(exe), 'deck.tcl'], cwd=tmp_path, env=env, check=True,
                   stdin=subprocess.DEVNULL, capture_output=True, timeout=120)
    got = (tmp_path / 'res.txt').read_text().split('\n')
    lam_tcl = [float(v) for v in got[0].split()]
    assert got[1].split() == ['0', '0', '0', '1', '1']
    _v3_build()
    ops.eigen('-fullGenLapack', 3)
    np.testing.assert_allclose(lam_tcl, ops.residualModes('-dir', 1), rtol=1e-12)
