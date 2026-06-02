"""Consistent-tangent properties of the finite-strain-native combined-hardening J2
(co-rotating backstress), pinned in numpy BEFORE the C++ FiniteStrainNDMaterial
subclass is built.

The Cauchy stress depends on the trial F through two channels — (A) the trial
elastic strain (the shipped log-strain spatial tangent) and (B) the co-rotated
backstress α̃ = R(f_Δ) α_n R(f_Δ)ᵀ. These tests establish, with data, that:

  * channel B is SMALL (~2–4e-4 relative) and bounded even for saturated
    backstress — so it does not threaten Newton robustness, but it sits at the
    tangent-FD test tolerance (2e-4), so the exact tangent needs it;
  * channel B is FRAME-OBJECTIVE (its relative size is rotation-invariant);
  * the C++ recipe — analytic channel A + numeric channel B by perturbing only R
    — reproduces the full ∂σ/∂F to machine precision (no kernel re-derivation).

ADR: Ladruno_implementation/16_finite_native_j2_adr.md.
"""
import numpy as np
import pytest

import ladrunoj2_reference as lr
import ladrunoj2_finite_native_reference as nat

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

I3 = np.eye(3)
# three-term Chaboche AF kinematic + Voce+linear isotropic (same family as the
# objectivity oracle test), strong enough to build appreciable backstress.
PARAMS = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                   C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])


def _rot(axis, ang):
    a = np.asarray(axis, float); a /= np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return I3 + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)


def _drive_monotonic(nsteps, stretch, prerot=None):
    """Commit a monotonic stretch path (optionally rigidly pre-rotated by `prerot`);
    return the material and a continuing trial F that yields plastically."""
    m = nat.NativeFiniteJ2(PARAMS, corotate=True)
    Q = prerot if prerot is not None else I3
    a = np.log(stretch)
    for k in range(1, nsteps + 1):
        ak = a * k / nsteps
        U = np.diag([np.exp(ak), np.exp(-ak / 2), np.exp(-ak / 2)])
        m.setTrialF(Q @ U); m.commit()
    ak = a * (nsteps + 1) / nsteps
    Feval = Q @ np.diag([np.exp(ak), np.exp(-ak / 2), np.exp(-ak / 2)])
    return m, Feval


def test_tangent_recipe_is_exact():
    """Numpy-level design validation (NOT a C++ check): the additive decomposition
    full ∂σ/∂F == channel A (R frozen) + channel B (R-only perturbation) holds to
    machine precision. This is what LICENSES the C++ implementation strategy —
    analytic channel A (spatial_tangent_full) + a numeric channel-B correction. The
    actual C++ channel-B tensor (the F→L map + conventions) is gated separately in
    tests/test_ladrunoJ2_finite_native_tangentB_cpp.py against this oracle."""
    m, Feval = _drive_monotonic(12, 1.6)
    A_full = m.tangent_dsigma_dF(Feval, freeze_R=False)
    A_frz = m.tangent_dsigma_dF(Feval, freeze_R=True)   # channel A
    B = m.channelB_dsigma_dF(Feval)                     # channel B
    err = np.linalg.norm(A_full - (A_frz + B)) / np.linalg.norm(A_full)
    assert err < 1e-7, f"A + B did not reconstruct the full tangent (err={err:.2e})"


def test_channelB_is_small_but_present():
    """Channel B is real and material-dependent: omitting it perturbs the tangent
    by ~1e-4 (stiff metals: G≫σ_y) up to ~2e-3 (soft elasticity / strong kinematic
    hardening, as here) — i.e. it straddles the 2e-4 tangent-FD gate, so the exact
    tangent needs it — yet it stays bounded (≪1%), so it never threatens Newton
    robustness."""
    m, Feval = _drive_monotonic(15, 1.7)
    A_full = m.tangent_dsigma_dF(Feval, freeze_R=False)
    A_frz = m.tangent_dsigma_dF(Feval, freeze_R=True)
    rel = np.linalg.norm(A_full - A_frz) / np.linalg.norm(A_full)
    anorm = max(np.linalg.norm(a) for a in m.alpha)
    assert anorm > PARAMS.sig0, "test did not build appreciable backstress"
    assert 1e-4 < rel < 5e-3, f"channel B out of expected band: {rel:.2e}"


def test_channelB_is_frame_objective():
    """Tangent objectivity: the relative size of channel B is invariant under a
    rigid pre-rotation of the whole load path (a rotated state has the rotated
    tangent), confirming the consistent tangent is frame-indifferent."""
    Q = _rot([0.3, 1.0, 0.5], 0.9)
    m0, Fe0 = _drive_monotonic(12, 1.6)
    mQ, FeQ = _drive_monotonic(12, 1.6, prerot=Q)

    def relB(m, F):
        a = m.tangent_dsigma_dF(F, freeze_R=False)
        return np.linalg.norm(a - m.tangent_dsigma_dF(F, freeze_R=True)) \
            / np.linalg.norm(a)

    assert relB(m0, Fe0) == pytest.approx(relB(mQ, FeQ), rel=1e-6), \
        "channel-B fraction changed under a rigid rotation — tangent not objective"
