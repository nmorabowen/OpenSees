"""The finite-strain-native combined-hardening J2 (co-rotating backstress) oracle
is OBJECTIVE under large rotation — the property the v1 LogStrain-over-LadrunoJ2
path provably lacks (dSNPO §14.11). This locks the push-forward form
(R_Δ = polar(f_Δ)) BEFORE the C++ FiniteStrainNDMaterial-subclass v2 is built.

ADR: Ladruno_implementation/16_finite_native_j2_adr.md.
Design boundary it fixes is pinned (the OTHER way) by
tests/test_ladrunoJ2_finite.py::test_finite_combined_hardening_large_rotation_objectivity_is_v2
(a strict xfail on the v1 wrapper).
"""
import numpy as np
import pytest

import ladrunoj2_reference as lr
import ladrunoj2_finite_native_reference as nat

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

I3 = np.eye(3)
# three-term Chaboche AF kinematic + Voce+linear isotropic
PARAMS = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                   C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])


def _load_to_plastic_then_Q(corotate):
    """Load into the plastic range along a fixed stretch (committing each step),
    then superpose a finite rigid rotation Q. Frame indifference requires
    σ(Q F_def) == Q σ(F_def) Qᵀ."""
    m = nat.NativeFiniteJ2(PARAMS, corotate=corotate)
    F_def = np.array([[1.25, 0.10, 0.0], [0.0, 1.0 / np.sqrt(1.25), 0.0],
                      [0.0, 0.0, 1.0 / np.sqrt(1.25)]])
    sig_ref = None
    for a in np.linspace(0.2, 1.0, 8):
        sig_ref, _ = m.setTrialF(I3 + a * (F_def - I3)); m.commit()
    th = 0.6
    c, s = np.cos(th), np.sin(th)
    Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    sig_rot, _ = m.setTrialF(Q @ F_def)
    return sig_ref, Q, sig_rot


def test_native_objective_under_superposed_rotation():
    """THE gate: a finite rigid rotation superposed on a plastic combined-hardening
    state rotates the Cauchy stress rigidly. (v1 fails this by ~6.)"""
    sig_ref, Q, sig_rot = _load_to_plastic_then_Q(corotate=True)
    assert np.allclose(sig_rot, Q @ sig_ref @ Q.T, rtol=1e-9, atol=1e-9), (
        "co-rotating-α native J2 is not objective under superposed rotation")


def test_fixed_frame_alpha_is_NOT_objective_control():
    """Negative control: the SAME machinery with α left in a fixed frame (v1
    behaviour) is materially non-objective — proves the test has teeth and that
    co-rotation is what fixes it."""
    sig_ref, Q, sig_rot = _load_to_plastic_then_Q(corotate=False)
    err = np.abs(sig_rot - Q @ sig_ref @ Q.T).max()
    assert err > 1.0, f"expected gross non-objectivity for fixed-frame α, got {err:.2e}"


def test_native_objective_along_rotating_plastic_path():
    """Stronger: accumulate plastic flow WHILE rotating, step by step. The Cauchy
    stress along the rotating path equals the body-fixed path rotated by R at
    every committed step (full path objectivity, not just superposed-rotation)."""
    rot = nat.NativeFiniteJ2(PARAMS, corotate=True)
    body = nat.NativeFiniteJ2(PARAMS, corotate=True)
    saw_plastic = False
    for k in range(1, 10):
        s = 1.0 + 0.03 * k
        g = 0.015 * k
        U = np.array([[s, g, 0.0], [0.0, 1.0 / np.sqrt(s), 0.0],
                      [0.0, 0.0, 1.0 / np.sqrt(s)]])
        th = 0.07 * k
        c, sn = np.cos(th), np.sin(th)
        R = np.array([[c, -sn, 0.0], [sn, c, 0.0], [0.0, 0.0, 1.0]])
        sig_rot, pl = rot.setTrialF(R @ U)
        sig_body, _ = body.setTrialF(U)
        saw_plastic = saw_plastic or pl
        assert np.allclose(sig_rot, R @ sig_body @ R.T, rtol=1e-8, atol=1e-8), (
            f"native objectivity broke at step {k} (rotating plastic path)")
        rot.commit(); body.commit()
    assert saw_plastic, "rotating path never yielded"


def test_native_reduces_to_v1_without_rotation():
    """With NO material rotation (symmetric F ⇒ polar = I) the push-forward is a
    no-op, so the native model reproduces the v1 (fixed-frame) result exactly —
    the v2 changes nothing outside the large-rotation regime."""
    Ffin = np.array([[1.18, 0.02, 0.01], [0.02, 0.93, 0.015], [0.01, 0.015, 0.95]])

    def run(corotate):
        m = nat.NativeFiniteJ2(PARAMS, corotate=corotate)
        out = []
        for k in range(1, 7):
            sig, _ = m.setTrialF(I3 + (k / 6.0) * (Ffin - I3)); out.append(sig); m.commit()
        return np.array(out)

    assert np.allclose(run(True), run(False), rtol=1e-10, atol=1e-10), (
        "co-rotation changed the no-rotation path — it must be a no-op there")


def test_native_reduces_to_small_strain_objectivity():
    """At tiny strain the native combined-hardening model is objective too (the
    small-strain limit of the lift)."""
    eps = 1.0e-5
    m = nat.NativeFiniteJ2(PARAMS, corotate=True)
    F = np.array([[1.0 + eps, 0.5 * eps, 0.0],
                  [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
                  [0.0, 0.0, 1.0 + 0.2 * eps]])
    sig_ref, _ = m.setTrialF(F)
    th = 0.5
    c, s = np.cos(th), np.sin(th)
    Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    sig_rot, _ = m.setTrialF(Q @ F)
    assert np.allclose(sig_rot, Q @ sig_ref @ Q.T, rtol=1e-7, atol=1e-9)
