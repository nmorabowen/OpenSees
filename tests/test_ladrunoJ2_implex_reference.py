"""IMPL-EX (Oliver–Huespe–Cante) oracle for combined-hardening J2 — locks the
Δγ-extrapolation scheme the C++ `nDMaterial LadrunoJ2[Finite] -implex` must match,
BEFORE the C++ (oracle-first). Driver/oracle: `ladrunoj2_implex_reference`.

The scheme (small-strain; the finite material layers the SAME idea on the co-rotated
return map): the implicit return map runs every step and is COMMITTED (so the plastic
history is identical to a fully-implicit run); the stress REPORTED to the solver
freezes an extrapolated multiplier Δγ̃ = Δγ_n·(Δt_{n+1}/Δt_n) and the committed flow
direction N_n, giving εᵖ̃ = εᵖ_n + Δγ̃·N_n and σ̃ = C:(ε − εᵖ̃). Since εᵖ̃ is a pure
history quantity, the reported tangent is the constant SPD elastic operator C.
"""
import numpy as np
import pytest

import ladrunoj2_reference as lr
import ladrunoj2_implex_reference as ix

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

PARAMS = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                   C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])
_BASE = np.array([[0.020, 0.012, 0.0],
                  [0.012, -0.008, 0.005],
                  [0.0, 0.005, -0.006]])   # deviatoric+volumetric, yields


def _monotone(t):    # smooth proportional path, single yield onset
    return t * _BASE


def test_committed_history_identical_to_implicit():
    """IMPL-EX commits the IMPLICIT return-map history verbatim — on a prescribed
    strain path the committed εᵖ, ε̄ᵖ, α are byte-identical to a fully-implicit run.
    (IMPL-EX changes only the REPORTED stress/tangent, never what is committed.)"""
    N = 200
    imp = lr.StatefulLadrunoJ2(PARAMS)
    imx = ix.ImplexLadrunoJ2(PARAMS)
    for k in range(1, N + 1):
        eps = _monotone(k / N)
        imp.setTrialStrain(eps); imp.commitState()
        imx.setTrialStrain(eps); imx.commitState()
    assert np.abs(imp.epsP - imx.epsP).max() < 1e-12
    assert abs(imp.ebarP - imx.ebarP) < 1e-12
    for a, b in zip(imp.alpha, imx.alpha):
        assert np.abs(a - b).max() < 1e-12


def test_tangent_is_constant_elastic_and_spd():
    """The reported IMPL-EX tangent is the elastic operator C — the same matrix at
    every step (elastic and plastic alike), symmetric and positive-definite. This is
    the whole point: a factor-once SPD global stiffness, no perfect-plasticity/
    softening tangent singularity."""
    C = ix.elastic_tangent_6x6(PARAMS)
    assert np.allclose(C, C.T)
    assert np.linalg.eigvalsh(C).min() > 0.0
    # constant across the (yielding) path
    m = ix.ImplexLadrunoJ2(PARAMS)
    saw_plastic = False
    for k in range(1, 81):
        m.setTrialStrain(_monotone(k / 80.0))
        assert np.allclose(m.getTangent(), C)         # never changes
        saw_plastic = saw_plastic or m._res["plastic"]
        m.commitState()
    assert saw_plastic, "path never yielded — tangent-constancy claim is vacuous"


def test_implex_converges_to_implicit_under_refinement():
    """In the sustained plastic regime the IMPL-EX (explicit) stress converges to the
    implicit stress under step refinement — observed order ≈2 (smooth proportional
    flow ⇒ the one-step extrapolation lag is second order). So the implex consistency
    error is controllable by Δt and vanishes in the limit."""
    def sustained_err(N):
        m = ix.ImplexLadrunoJ2(PARAMS)
        lo = int(0.6 * N)
        worst = 0.0
        for k in range(1, N + 1):
            sx = m.setTrialStrain(_monotone(k / N))
            si = m.getStressImplicit()
            if k > lo:
                worst = max(worst, np.abs(sx - si).max())
            m.commitState()
        return worst

    Ns = [100, 200, 400, 800]
    errs = [sustained_err(N) for N in Ns]
    for e_c, e_f in zip(errs[:-1], errs[1:]):
        assert e_f < e_c
        assert e_c / e_f > 3.0                         # order > 1.58
    order = np.log(errs[0] / errs[-1]) / np.log(Ns[-1] / Ns[0])
    assert 1.7 < order < 2.3, f"observed order {order:.2f} not ~2"
    assert errs[-1] < 1e-5            # tight absolute agreement at the finest resolution


def test_yield_onset_has_one_step_lag():
    """Documented IMPL-EX behaviour: at the FIRST plastic step the extrapolation has
    no increment to carry (Δγ_n = 0), so the reported stress is the ELASTIC trial —
    it lags the implicit (plastic) stress by one step. This transient is O(Δt) and
    vanishes under refinement; it is the price of the constant SPD tangent."""
    N = 200
    m = ix.ImplexLadrunoJ2(PARAMS)
    onset = None
    for k in range(1, N + 1):
        epsP_before = m.epsP.copy()
        dG_before = m.dGamma_n
        eps = _monotone(k / N)
        sx = m.setTrialStrain(eps)
        if m._res["plastic"] and onset is None:
            onset = dict(sx=sx.copy(), si=m.getStressImplicit().copy(),
                         elastic=ix._elastic_stress(PARAMS, eps, epsP_before),
                         dG_before=dG_before)
        m.commitState()
    assert onset is not None, "path never yielded"
    assert onset["dG_before"] == 0.0                   # extrapolation had no increment
    # implex == elastic trial (no plastic correction yet)
    assert np.abs(onset["sx"] - onset["elastic"]).max() < 1e-12
    # and it genuinely lags the implicit plastic stress
    assert np.abs(onset["sx"] - onset["si"]).max() > 1e-4
