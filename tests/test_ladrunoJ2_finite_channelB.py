"""ANALYTIC channel B (ADR-16 P3) — locks the analytic replacement for the numeric
R-perturbation channel-B in the C++ LadrunoJ2Finite consistent tangent, BEFORE the C++.
Oracle: ladrunoj2_finite_channelB_reference.

Channel B = the part of ∂σ/∂F flowing through the co-rotated backstress α̃=RαRᵀ
(εᵉᵗʳ, J fixed). The analytic form is the chain ∂R/∂F (polar-rotation derivative,
Sylvester) ∘ ∂α̃/∂R ∘ ∂τ/∂α̃ (return-map backstress sensitivity). Validated against the
shipped numeric channelB_dsigma_dF and against the full-minus-channelA decomposition.
"""
import numpy as np
import pytest

import ladrunoj2_reference as lr
import ladrunoj2_finite_native_reference as nat
import ladrunoj2_finite_channelB_reference as cb

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

STIFF = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                  C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])
SOFT = lr.Params(K=200.0, G=80.0, sig0=5.0, Qinf=3.0, bIso=20.0, Hiso=5.0,
                 C=[400.0, 200.0], gam=[150.0, 40.0])


def _rotz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _committed_plastic(params, nstep=6):
    m = nat.NativeFiniteJ2(params, corotate=True)
    for k in range(1, nstep + 1):
        s = 1.0 + 0.04 * k
        U = np.array([[s, 0.02 * k, 0.0], [0.0, 1.0/np.sqrt(s), 0.0],
                      [0.0, 0.0, 1.0/np.sqrt(s)]])
        m.setTrialF(_rotz(0.08 * k) @ U); m.commit()
    return m


def test_polar_derivative_matches_fd():
    """dR = R·Ω with Ω from (tr U·I − U) ω = axial(A−Aᵀ) matches a central FD of the
    polar rotation."""
    f = np.array([[1.2, 0.15, 0.03], [0.05, 0.92, 0.08], [0.02, 0.04, 0.97]])
    df = np.array([[0.01, -0.02, 0.005], [0.013, 0.004, -0.011], [-0.006, 0.009, 0.002]])
    h = 1.0e-7
    an = cb.polar_derivative(f, df)
    fd = (nat.polar_rotation(f + h * df) - nat.polar_rotation(f - h * df)) / (2.0 * h)
    assert np.abs(an - fd).max() < 1.0e-7


@pytest.mark.parametrize("params,Fprobe", [
    (STIFF, _rotz(0.6) @ np.array([[1.30, 0.06, 0.0], [0.0, 0.90, 0.0], [0.0, 0.0, 0.95]])),
    (SOFT,  _rotz(0.7) @ np.array([[1.35, 0.08, 0.0], [0.0, 0.88, 0.0], [0.0, 0.0, 0.94]])),
])
def test_channelB_analytic_matches_numeric(params, Fprobe):
    """The analytic channel B equals the shipped numeric R-perturbation channel B to
    finite-difference precision — across stiff (G≫σ_y) and soft (strong kinematic,
    low G, where channel B is largest) states."""
    m = _committed_plastic(params)
    Ban = cb.channelB_analytic(m, Fprobe)
    Bnum = m.channelB_dsigma_dF(Fprobe, h=1.0e-7)
    den = max(np.abs(Bnum).max(), 1.0e-30)
    assert np.abs(Bnum).max() > 1.0e-3, "channel B vanished — probe is not exercising it"
    assert np.abs(Ban - Bnum).max() / den < 1.0e-6


def test_channelB_plus_channelA_is_full_tangent():
    """Decomposition check: channel A (the R-frozen ∂σ/∂F) + analytic channel B equals
    the full ∂σ/∂F — i.e. the analytic B is exactly the missing co-rotation term."""
    m = _committed_plastic(STIFF)
    F = _rotz(0.55) @ np.array([[1.28, 0.05, 0.0], [0.0, 0.91, 0.0], [0.0, 0.0, 0.96]])
    full = m.tangent_dsigma_dF(F, h=1.0e-7)
    A = m.tangent_dsigma_dF(F, freeze_R=True, h=1.0e-7)
    Ban = cb.channelB_analytic(m, F)
    den = max(np.abs(full).max(), 1.0e-30)
    assert np.abs((A + Ban) - full).max() / den < 1.0e-5


def test_channelB_zero_when_elastic():
    """No plastic flow ⇒ the elastic predictor is backstress-independent ⇒ channel B = 0.
    (Virgin material, sub-yield strain — a loaded-then-tiny-increment state stays
    plastic, so it must be a fresh elastic probe.)"""
    m = nat.NativeFiniteJ2(STIFF, corotate=True)
    eps = 1.0e-3                                   # sub-yield (2G·eps ≈ 1.4 ≪ √⅔·sig0)
    F_el = np.array([[1.0 + eps, 0.0, 0.0], [0.0, 1.0 - 0.3 * eps, 0.0],
                     [0.0, 0.0, 1.0 - 0.3 * eps]])
    _, plastic = m.setTrialF(F_el)
    assert not plastic, "probe unexpectedly yielded — not an elastic test"
    B = cb.channelB_analytic(m, F_el)
    assert np.abs(B).max() == 0.0
