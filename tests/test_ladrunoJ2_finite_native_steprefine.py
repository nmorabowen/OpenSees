"""Closes the ADR-16 P1 caveat for LadrunoJ2Finite: the *magnitude accuracy* of the
co-rotated Armstrong–Frederick backstress under SIMULTANEOUS large rotation AND large
stretch.

The shipped material transports the backstress each step by R_Δ = polar(f_Δ). ADR
§"Residual open item" flagged that this incremental polar transport differs at higher
order from the exact §14.11 exponential-map (continuous-spin) transport when a step
has both large spin and large stretch, and asked to *verify by step-refinement against
a fine-increment reference; only if a discrepancy appears, swap to exp-map transport.*

This is that verification (numpy oracle = "enough" per the ADR). Driver/oracle:
`ladrunoj2_finite_native_steprefine_reference`, driving the SAME return map
(`ladrunoj2_reference`) + Hencky kinematics (`logstrain_reference`) the C++ material
uses, along F(s)=exp(s·(D+W)) with NON-COAXIAL stretch-rate D and spin W (the regime
where the transport choice actually matters; a fixed-axis stretch would be degenerate).

Verdict (the assertions below): the polar scheme converges at order ≈1, and the
polar-vs-exact-transport difference is O(Δs²) — 3–4 orders of magnitude below the
overall step-discretization error. Swapping to exp-map transport would change nothing
measurable ⇒ NO swap. The shipped polar(f_Δ) transport is adequate.
"""
import numpy as np
import pytest

import ladrunoj2_finite_native_reference as nat
import ladrunoj2_finite_native_steprefine_reference as sr

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_NREF = 2048   # fine-increment reference resolution
_M = 64        # sub-steps for the "exact" continuous-spin transport (converged, see ref)


def test_polar_driver_matches_shipped_native_material():
    """The step-refinement driver's transport="polar" path IS the shipped logic:
    driving NativeFiniteJ2 (== C++ LadrunoJ2Finite) along the same F(s) gives the
    identical committed state. Without this, the convergence study below would be
    testing a different model than the one that ships."""
    N = 20
    m = nat.NativeFiniteJ2(sr.PARAMS, corotate=True)
    sig = None
    for n in range(1, N + 1):
        sig, _ = m.setTrialF(sr.F_of_s(n / N))
        m.commit()
    r = sr.run(N, "polar")
    assert np.abs(sig - r["sigma"]).max() < 1e-12
    for k in range(sr.PARAMS.nBack):
        assert np.abs(m.alpha[k] - r["alpha"][k]).max() < 1e-12
    assert abs(m.ebarP - r["ebarP"]) < 1e-12


def test_path_is_severe_noncoaxial_and_plastic():
    """Guard the test's teeth: the path really does rotate >90° while stretching, is
    (near-)isochoric, drives sustained plastic flow, and — crucially — the polar(f_Δ)
    and exact sub-stepped transport rotations genuinely DIFFER over a coarse step (a
    fixed-axis stretch would make every transport identical and the study vacuous)."""
    F1 = sr.F_of_s(1.0)
    assert abs(np.linalg.det(F1) - 1.0) < 1e-6                  # isochoric (W skew, D traceless)
    R1 = nat.polar_rotation(F1)
    angle = np.degrees(np.arccos((np.trace(R1) - 1.0) / 2.0))
    assert angle > 90.0                                         # large rotation
    U1 = R1.T @ F1
    stretches = np.linalg.eigvalsh(0.5 * (U1 + U1.T))
    assert stretches.max() > 1.15 and stretches.min() < 0.90    # large stretch

    r = sr.run(24, "polar")
    assert r["n_plastic"] == 24                                 # plastic every step

    Rp = sr._transport_rotation(0.0, 1.0 / 6.0, "polar", 1)
    Re = sr._transport_rotation(0.0, 1.0 / 6.0, "exact", 512)
    assert np.linalg.norm(Rp - Re) > 1e-5, (
        "transport rotations coincide — path is degenerate, study has no teeth")


def test_polar_scheme_converges_under_refinement():
    """The shipped polar-transport scheme CONVERGES under step refinement — the
    co-rotated AF magnitude tends to a well-defined limit (no O(1) discrepancy /
    inconsistency). Observed order ≈ 1 (backward-Euler return map dominates)."""
    v_ref = sr.state_vector(sr.run(_NREF, "polar"))
    Ns = [16, 32, 64, 128, 256]
    errs = [np.linalg.norm(sr.state_vector(sr.run(N, "polar")) - v_ref) for N in Ns]
    # monotone decreasing, each refinement at least ~order 0.77 (ratio > 1.7)
    for e_coarse, e_fine in zip(errs[:-1], errs[1:]):
        assert e_fine < e_coarse
        assert e_coarse / e_fine > 1.7
    # mid-range observed order is ~1
    order = np.log(errs[1] / errs[2]) / np.log(2.0)
    assert 0.85 < order < 1.35
    # the fine reference is genuinely converged
    assert errs[-1] / np.linalg.norm(v_ref) < 5e-3


def test_polar_and_exact_transport_share_the_same_limit():
    """polar(f_Δ) and the exact continuous-spin (exp-map) transport converge to the
    SAME state; their difference is second order in the step (≈O(Δs²)) — i.e. the
    'higher-order difference' the ADR flagged is real but vanishes faster than the
    first-order discretization error. So the two transports give the same physical
    answer in the limit."""
    Ns = [12, 24, 48, 96]
    gaps = [np.linalg.norm(sr.state_vector(sr.run(N, "polar"))
                           - sr.state_vector(sr.run(N, "exact", M=_M))) for N in Ns]
    for g_coarse, g_fine in zip(gaps[:-1], gaps[1:]):
        assert g_fine < g_coarse
        assert g_coarse / g_fine > 3.0          # order > 1.58 → clearly super-linear
    order = np.log(gaps[0] / gaps[-1]) / np.log(Ns[-1] / Ns[0])
    assert order > 1.6, f"transport gap order {order:.2f} not second-order-ish"


def test_transport_choice_is_negligible_vs_discretization():
    """THE verdict that closes P1: at every (coarse) resolution the polar-vs-exact
    *transport-choice* difference is far smaller than the *step-discretization* error
    of the scheme itself. So switching the shipped polar(f_Δ) transport to the exact
    §14.11 exp-map transport cannot meaningfully improve accuracy ⇒ no swap needed."""
    v_ref = sr.state_vector(sr.run(_NREF, "polar"))
    for N in [12, 24, 48]:
        v_polar = sr.state_vector(sr.run(N, "polar"))
        v_exact = sr.state_vector(sr.run(N, "exact", M=_M))
        self_err = np.linalg.norm(v_polar - v_ref)          # step-discretization error
        transport_err = np.linalg.norm(v_polar - v_exact)   # cost of the transport choice
        assert transport_err < 0.01 * self_err, (
            f"N={N}: transport choice ({transport_err:.2e}) is not negligible vs "
            f"discretization ({self_err:.2e})")
