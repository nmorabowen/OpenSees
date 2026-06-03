"""IMPL-EX over the finite-strain-native combined-hardening J2 — locks the C++
`nDMaterial LadrunoJ2Finite -implex` behaviour BEFORE the C++ (oracle-first).
Oracle: `ladrunoj2_finite_implex_reference`; the finite path F(s)=exp(s·(D+W)) is
reused from the step-refinement reference (severe rotation+stretch).

Layers the classic Δγ-extrapolation IMPL-EX on the co-rotated finite return map: the
implicit co-rotated map is committed (finite plastic history identical to a fully-
implicit run); the reported stress freezes Δγ̃ and the committed log-strain flow
direction N_n (co-rotated by the SAME R_Δ=polar(f_Δ) as the backstress), giving the
constant SPD elastic log-strain tangent and a frame-indifferent explicit stress.
"""
import numpy as np
import pytest

import ladrunoj2_finite_native_reference as nat
import ladrunoj2_finite_implex_reference as fx
import ladrunoj2_finite_native_steprefine_reference as sr

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

PARAMS = sr.PARAMS
F_of_s = sr.F_of_s


def test_committed_history_identical_to_implicit_finite():
    """IMPL-EX commits the IMPLICIT co-rotated finite history verbatim — the committed
    bᵉ_n, ε̄ᵖ_n, α_n are byte-identical to a fully-implicit NativeFiniteJ2 run on the
    same F-path. (IMPL-EX changes only the reported stress/tangent.)"""
    N = 120
    imp = nat.NativeFiniteJ2(PARAMS, corotate=True)
    imx = fx.ImplexNativeFiniteJ2(PARAMS)
    for k in range(1, N + 1):
        F = F_of_s(k / N)
        imp.setTrialF(F); imp.commit()
        imx.setTrialF(F); imx.commit()
    assert np.abs(imp.Be_n - imx.Be_n).max() < 1e-12
    assert abs(imp.ebarP - imx.ebarP) < 1e-12
    for a, b in zip(imp.alpha, imx.alpha):
        assert np.abs(a - b).max() < 1e-12


def test_implex_converges_to_implicit_under_F_refinement():
    """In the sustained plastic regime the finite IMPL-EX stress converges to the
    implicit stress under F-step refinement — observed order ≈2, so the implex
    consistency error is Δt-controllable and vanishes in the limit."""
    def sustained(N):
        m = fx.ImplexNativeFiniteJ2(PARAMS)
        lo = int(0.6 * N)
        worst = 0.0
        for k in range(1, N + 1):
            sx = m.setTrialF(F_of_s(k / N))
            si = m.getStressImplicit()
            if k > lo:
                worst = max(worst, np.abs(sx - si).max())
            m.commit()
        return worst

    Ns = [120, 240, 480, 960]
    errs = [sustained(N) for N in Ns]
    for e_c, e_f in zip(errs[:-1], errs[1:]):
        assert e_f < e_c
        assert e_c / e_f > 3.0                       # order > 1.58
    order = np.log(errs[0] / errs[-1]) / np.log(Ns[-1] / Ns[0])
    assert 1.7 < order < 2.3, f"observed order {order:.2f} not ~2"


def test_implex_objective_under_superposed_rotation():
    """The finite IMPL-EX stress is frame-indifferent: the committed flow direction
    N_n co-rotates by R_Δ=polar(f_Δ) exactly like the backstress, so a superposed
    rigid rotation Q on a plastic state rotates σ̃ rigidly. (This is the property the
    whole §14.11 v2 co-rotation exists for, preserved by the implex stress.)"""
    m = fx.ImplexNativeFiniteJ2(PARAMS)
    for k in range(1, 9):
        m.setTrialF(F_of_s(0.7 * k / 8)); m.commit()
    F_def = F_of_s(0.7)
    sig_ref = m.setTrialF(F_def)                     # staged, not committed
    th = 0.6
    c, s = np.cos(th), np.sin(th)
    Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    sig_rot = m.setTrialF(Q @ F_def)
    assert np.allclose(sig_rot, Q @ sig_ref @ Q.T, rtol=1e-9, atol=1e-9)


def test_implex_tangent_is_constant_elastic_and_spd():
    """The finite IMPL-EX material tangent is the constant SPD elastic operator — no
    plastic h term and no backstress co-rotation 'channel B' (α drops out of the
    explicit stress)."""
    C = fx.ImplexNativeFiniteJ2(PARAMS).getTangent6x6()
    assert np.allclose(C, C.T)
    assert np.linalg.eigvalsh(C).min() > 0.0


def test_implex_reduces_to_small_strain_implex():
    """At tiny strain the finite IMPL-EX reduces to the small-strain IMPL-EX oracle
    (the lift is consistent in the infinitesimal limit)."""
    import ladrunoj2_reference as lr
    import ladrunoj2_implex_reference as ix
    base = np.array([[0.020, 0.012, 0.0],
                     [0.012, -0.008, 0.005],
                     [0.0, 0.005, -0.006]])
    scale = 1.0e-4
    N = 12
    fin = fx.ImplexNativeFiniteJ2(PARAMS)
    sml = ix.ImplexLadrunoJ2(PARAMS)
    for k in range(1, N + 1):
        eps = scale * (k / N) * base
        F = np.eye(3) + eps                          # symmetric ⇒ polar = I, F ≈ I+ε
        sfin = fin.setTrialF(F); fin.commit()
        ssml = sml.setTrialStrain(eps); sml.commitState()
        assert np.allclose(sfin, ssml, rtol=1e-5, atol=1e-6 * PARAMS.G * scale)
