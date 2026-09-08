"""ADR-97 P1 (wp/97b) — gate 6: every new failure path is LOUD, and every
refusal the ADR promised is actually enforced.

The class of defect ADR-94 found is not "wrong answer" but "wrong answer that
converges": a silent default to a different integrator, a dropped option, an
un-dispatched enum value, a return code no host element looks at.  So every new
code path in `Closest_Point` ends in `LADRUNO_MATERIAL_REFUSED` (never a bare
-1: `LadrunoBrick` compares against the sentinel and `stdBrick` drops
everything, ADR-94 B2), and every combination the ADR declared unsupported is
refused AT PARSE TIME, naming the tokens.

Covered here:

* a starved `n_max_iterations` on a plastic `Closest_Point` step fails the
  analysis instead of committing a non-converged state;
* `stdBrick` swallows that same refusal — pinned as a NEGATIVE control, because
  it is the reason every refusal test in this suite uses `LadrunoBrick`;
* `tangent_type Algorithmic` with any integrator other than `Closest_Point` is
  refused (ADR-97 D2 — a consistent tangent is defined only relative to a
  specific committed map);
* `integration_method Closest_Point` on a family P1 has not converted
  (a MIXED MohrCoulomb_YF x VonMises_PF pairing) is refused,
  naming the ADR phase that will deliver it -- MohrCoulomb x MohrCoulomb
  itself is SHIPPED by ADR-97 P2 and the assertion was inverted there;
* unknown tokens are still rejected (the ADR-94 contract);
* the ADR-94 B4 reproducer — Drucker-Prager driven into hydrostatic TENSION,
  past its apex — never commits a NaN under `Closest_Point`.

Zone-A, ~10 s.
"""
import numpy as np
import pytest

from _testbed import ops

import test_adr97_p1_smooth as S  # noqa: E402
import test_asdplastic_mctc as M  # noqa: E402

pytestmark = [pytest.mark.zone_a]

MC_PARAMS = ["YoungsModulus", M.E, "PoissonsRatio", M.NU,
             "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
             "MassDensity", 0.0]


def _mat_mc(tag, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters", *MC_PARAMS, "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


# Ladruno (ADR-97 wp/97c): the family gate is now about MIXED pairings, not about
# MohrCoulomb itself, and wp/97d ships HoekBrown.  `MohrCoulomb_YF` x `VonMises_PF` is one
# specializations the generator really registers and that P2 must still refuse.
def _mat_mc_mixed(tag, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "VonMises_PF", "LinearIsotropic3D_EL",
            "BackStress(TensorLinearHardeningFunction):",
            "Begin_Model_Parameters", *MC_PARAMS,
            "TensorLinearHardeningParameter", 0.0, "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


def _mat_hb(tag, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "HoekBrown_YF", "HoekBrown_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            "YoungsModulus", M.E, "PoissonsRatio", M.NU,
            # Ladruno (ADR-97 wp/97d): the parameter is `HB_sigci`, not
            # `HB_sigma_ci`.  With the wrong spelling this deck was rejected for
            # a MISSING PARAMETER, so the refusal assertion below passed without
            # ever exercising the family gate -- a vacuous test in both P1 and
            # P2.  Fixed here, and the assertion inverted (P3 ships HoekBrown).
            "HB_sigci", 50000.0, "HB_mb", 2.396510364418,
            "HB_s", 0.011743628457, "HB_a", 0.502840500848,
            "HB_mb_psi", 2.396510364418, "HB_ds", 0.0, "MassDensity", 0.0,
            "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


@pytest.fixture(scope="module")
def cp_available():
    if not S._constructible(lambda t: S.mat_vm(t)):
        pytest.skip("ASDPlasticMaterial3D Closest_Point not available")


@pytest.fixture(scope="module")
def mc_available():
    if not S._constructible(lambda t: _mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D MohrCoulomb not available")


# ===========================================================================
# 1. parser refusals
# ===========================================================================
def test_closest_point_with_algorithmic_is_the_positive_control(cp_available):
    """The combination the ADR ships must build -- otherwise every refusal
    below would pass for the wrong reason."""
    assert S._constructible(
        lambda t: S.mat_vm(t, method="Closest_Point", tangent="Algorithmic"))


@pytest.mark.parametrize("method", ["Backward_Euler", "Forward_Euler",
                                    "Modified_Euler_Error_Control",
                                    "Runge_Kutta_45_Error_Control"])
def test_algorithmic_is_refused_with_any_other_integrator(cp_available, method):
    """ADR-97 D2.  `tangent_type Algorithmic`'s enum value has existed since
    upstream with NO dispatch case anywhere -- it was dead, which is why nobody
    has been silently getting it.  Offering it on the cutting-plane
    `Backward_Euler` would ship a FOURTH almost-right tangent, which is exactly
    the class of defect ADR-94 M3 found."""
    assert not S._constructible(
        lambda t: S.mat_vm(t, method=method, tangent="Algorithmic")), (
        "tangent_type Algorithmic was accepted with integration_method %s -- "
        "the ADR-97 D2 cross-refusal is gone" % method)


def test_algorithmic_refusal_does_not_break_the_other_tangents(cp_available):
    """The refusal must be surgical: every tangent_type the chosen integrator
    does define still builds."""
    for tg in ("Secant", "Continuum", "Elastic",
               "Numerical_Algorithmic_FirstOrder",
               "Numerical_Algorithmic_SecondOrder"):
        assert S._constructible(
            lambda t, g=tg: S.mat_vm(t, method="Backward_Euler", tangent=g)), tg
        assert S._constructible(
            lambda t, g=tg: S.mat_vm(t, method="Closest_Point", tangent=g)), tg


def test_closest_point_is_refused_for_unconverted_families(mc_available):
    """ADR-97 D3: the closest-point map is added family by family.  A
    specialization whose yield function, plastic flow direction or hardening law
    has not opted in must be refused at parse time -- NOT silently run on the
    inert zero / finite-difference defaults, which is the ADR-94 M3 failure mode.

    Ladruno (ADR-97 wp/97c): this test used to assert that MohrCoulomb ITSELF was
    refused.  P2 ships it (a principal-stress-space multi-surface return; see
    tests/test_adr97_p2_principal.py), so the MohrCoulomb half of the assertion
    is INVERTED here -- deliberately, and recorded in the ADR-97 P2 report.  What
    the family gate must still refuse is a MIXED pairing: `MohrCoulomb_YF` with a
    non-Mohr-Coulomb plastic flow direction is covered by NO oracle (the
    principal-space return assumes both the surface and the potential are
    piecewise linear, and P1's smooth 6D map cannot use MohrCoulomb's Lode-angle
    gradient), so it stays refused.

    Ladruno (ADR-97 wp/97d): the HoekBrown half is now INVERTED too -- P3 ships
    `HoekBrown_YF` x `HoekBrown_PF` (a principal-space return to the CURVED
    surface; see tests/test_adr97_p3_hoekbrown.py).  The mixed HoekBrown
    pairings stay refused and are covered there."""
    assert S._constructible(lambda t: _mat_mc(t)), "control MC deck must build"
    assert S._constructible(
        lambda t: _mat_mc(t, opts=["integration_method", "Closest_Point"])), (
        "integration_method Closest_Point is refused for MohrCoulomb x "
        "MohrCoulomb -- ADR-97 P2 ships it")
    assert not S._constructible(
        lambda t: _mat_mc_mixed(t, opts=["integration_method",
                                         "Closest_Point"])), (
        "integration_method Closest_Point was accepted for the MIXED pairing "
        "MohrCoulomb_YF x VonMises_PF, which no oracle covers -- the ADR-97 D3 "
        "family gate is gone")
    assert S._constructible(lambda t: _mat_hb(t)), (
        "the control HoekBrown deck does not build at all -- fix its parameter "
        "list before reading anything into the assertion below")
    assert S._constructible(
        lambda t: _mat_hb(t, opts=["integration_method", "Closest_Point"])), (
        "integration_method Closest_Point is refused for HoekBrown_YF x "
        "HoekBrown_PF -- ADR-97 P3 ships it")


def test_unknown_tokens_are_still_rejected(cp_available):
    """The ADR-94 wp/94a contract: an unrecognised `integration_method` or
    `tangent_type` is an ERROR, not a silent default."""
    assert not S._constructible(
        lambda t: S.mat_vm(t, method="Closest_Pont", tangent="Secant"))
    assert not S._constructible(
        lambda t: S.mat_vm(t, method="Closest_Point", tangent="Algorythmic"))


# ===========================================================================
# 2. runtime refusals
# ===========================================================================
def _starved(ele):
    """A plastic Closest_Point step with `n_max_iterations` far too small."""
    return S.drive(lambda t: S.mat_vm(t, hiso=7000.0, niter=1),
                   S.VM_PATHS["triaxial"], nstep=10, ele=ele)


def test_starved_newton_is_refused_not_committed(cp_available):
    """`Closest_Point` returns LADRUNO_MATERIAL_REFUSED on exhaustion; the host
    element propagates it and the analysis fails.  Upstream's habit -- and
    `Backward_Euler`'s before ADR-84 P2a -- was to fall out of the loop and
    commit the non-converged state as success."""
    r = _starved("LadrunoBrick")
    print("starved n_max_iterations=1 codes:", r["codes"])
    assert any(c != 0 for c in r["codes"]), (
        "a starved Closest_Point Newton was committed as success")


def test_stdbrick_swallows_the_refusal_negative_control(cp_available):
    """PINNED negative control (ADR-94 B2).  `stdBrick` discards the material's
    return code by design, so the SAME starved deck 'succeeds' on it.  This is
    why every refusal test in the ADR-97 suite runs on `LadrunoBrick`."""
    r = _starved("stdBrick")
    print("starved on stdBrick codes:", r["codes"])
    assert all(c == 0 for c in r["codes"]), (
        "stdBrick now propagates material return codes -- good news, but the "
        "ADR-94 B2 pin and the ADR-97 host-choice rationale both need updating")


def test_strict_convergence_is_accepted_and_inert_on_a_converging_deck(
        cp_available):
    """`Closest_Point` ends every commit path in the existing
    `ladruno_strict_rejects` gate, on the same contract `Backward_Euler` has.
    On a deck that converges, turning the flag on must change nothing."""
    off = S.drive(lambda t: S.mat_vm(t, hiso=7000.0),
                  S.VM_PATHS["triaxial"], nstep=10)
    on = S.drive(lambda t: S.mat_vm(t, hiso=7000.0, strict=1),
                 S.VM_PATHS["triaxial"], nstep=10)
    assert off["codes"] == on["codes"] == [0] * 10
    d = S._rel(on["sigma"][-1], off["sigma"][-1])
    print("strict_convergence on/off rel difference = %.3e" % d)
    assert d == 0.0, "strict_convergence changed a converging Closest_Point deck"


def test_strict_convergence_still_refuses_a_starved_deck(cp_available):
    r = S.drive(lambda t: S.mat_vm(t, hiso=7000.0, niter=1, strict=1),
                S.VM_PATHS["triaxial"], nstep=10)
    assert any(c != 0 for c in r["codes"])


def test_dp_hydrostatic_tension_never_commits_nan_under_closest_point():
    """ADR-94 B4's reproducer, re-run on the new integrator.  Pure hydrostatic
    TENSION drives Drucker-Prager straight past its apex, where the flank return
    would have to make `sqrt(J2)` negative.  `Closest_Point` classifies the apex
    region in the ELASTIC metric and returns to the vertex; nothing NaN may
    reach the committed state, and nothing inadmissible either."""
    if not S._constructible(lambda t: S.mat_dp(t)):
        pytest.skip("DruckerPrager Closest_Point unavailable")
    legs = [np.array([2.4e-3, 2.4e-3, 2.4e-3, 0., 0., 0.])]   # 2x the apex strain
    r = S.drive(lambda t: S.mat_dp(t), legs, nstep=10, want=("stresses",
                                                             "BackStress"))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    assert np.all(np.isfinite(r["sigma"])), r["sigma"]
    p_apex = S.DP_XI_C / S.DP_ETA
    print("DP hydrostatic tension: final sigma = %s (apex p = %.6f)"
          % (np.array2string(r["sigma"][-1], precision=8), p_apex))
    assert abs(float(np.mean(r["sigma"][-1][:3])) - p_apex) < 1e-6 * p_apex
    for i in range(len(r["sigma"])):
        f = S._f_dp(r["sigma"][i], r["BackStress"][i])
        assert f <= 1e-6, (i, f)
