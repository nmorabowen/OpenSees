"""ADR-94 wp/94c -- the FIX gate for B5 (shear-slot convention), B4 (Drucker-
Prager gradients + the dead apex return) and M5 (unit-dependent yield tolerance).

Unlike the ``test_adr94_*`` sentinel files, nothing here pins a defect: every
assertion states what the fixed material must do, so a regression turns it red.

WHAT IS MEASURED
----------------
* **B5 -- one convention.**  A von Mises SIMPLE-SHEAR plastic path is compared
  against the closed-form radial return in
  ``Ladruno_implementation/adr94_oracle/vm_shear_oracle.py``.  Shear is the only
  place the two conventions differ, so this single number covers the yield
  normal, the flow direction, the plastic modulus ``n:(E:m)`` and the hardening
  rate ``H*sqrt(2/3 m:m)`` at once.  Pre-94c the same comparison is off by O(1).
* **B4 -- Drucker-Prager gradients.**  A sheared DP path must give the same
  converged stress under two independent tangent operators.  Before the fix the
  analytical ``Continuum`` tangent was built from a gradient with a 0.97
  relative error against a central difference, so it disagreed with the
  numerically differentiated one.
* **B4 -- apex.**  Hydrostatic tension driven past ``p = xi_c/eta`` must commit a
  finite stress that sits ON the yield surface, instead of running the flank
  return map past ``sqrt(J2) = 0``.
* **M5 -- ``f_relative_tol``.**  Default 0 must be byte-identical; switched on,
  the same physical problem must behave the same in kPa and in Pa.

TRAPS OBEYED (ADR-94 Sec. 8)
----------------------------
``LadrunoBrick`` everywhere a refusal could matter (``stdBrick`` swallows
material return codes); ``system("UmfPack")``; fully sp-prescribed rigs, which
this fork's Transformation handler accepts.
"""
import math
import os
import sys

import numpy as np
import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "Ladruno_implementation",
                                "adr94_oracle"))
import vm_shear_oracle as VMO  # noqa: E402

import test_adr94_hlist_numerics as N  # noqa: E402
import test_adr94_redblue_numerics as R  # noqa: E402
import test_adr94_hlist_hb as HB  # noqa: E402

pytestmark = [pytest.mark.zone_a]


# ---------------------------------------------------------------------------
# rig -- homogeneous SIMPLE SHEAR of the unit cube: u_x = gamma * z, u_y = u_z = 0.
# Every DOF is prescribed, so the strain field is exactly gamma_xz = gamma and
# the material point is driven, not solved for.
# ---------------------------------------------------------------------------
_TOP = (5, 6, 7, 8)
_BOT = (1, 2, 3, 4)


def _shear_build(mat_fn, nsteps, gamma, ele="LadrunoBrick"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in R._CUBE.items():
        ops.node(t, *map(float, c))
    for t in _BOT:
        ops.fix(t, 1, 1, 1)
    for t in _TOP:
        ops.fix(t, 0, 1, 1)
    mat_fn(1)
    ops.element(ele, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TOP:
        ops.sp(t, 1, gamma)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


def _drive_eps_sig(nsteps):
    """(codes, strains, stresses) of Gauss point 1 after each successful step."""
    codes, eps, sig = [], [], []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        ops.eleResponse(1, "forces")
        eps.append(list(ops.eleResponse(1, "strains"))[0:6])
        sig.append(list(ops.eleResponse(1, "stresses"))[0:6])
    return codes, np.array(eps), np.array(sig)


@pytest.fixture(scope="module")
def vm_available():
    try:
        _shear_build(lambda t: N.mat_vm(t), 2, 1.0e-4)
    except Exception as exc:                                  # pragma: no cover
        pytest.skip(f"ASDPlasticMaterial3D / VonMises / LadrunoBrick: {exc}")
    ops.wipe()


# ===========================================================================
# 1. B5 -- the decisive shear measurement
# ===========================================================================
# G = 26923; yield at tau = SQRT_2_over_3*sy/sqrt(2) = 17.32, i.e.
# gamma_y = 6.43e-4.  8e-3 is ~12x past yield.
GAMMA_END = 8.0e-3
SHEAR_STEPS = 10


@pytest.mark.t0m
def test_C1_vm_simple_shear_matches_the_radial_return_oracle(vm_available):
    """ADR-94 B5.  With ONE Voigt convention the Backward_Euler map on a pure
    simple-shear path is the closed-form radial return, to round-off.

    The two conventions differ ONLY in the shear slots, so this is the sharpest
    possible statement of the fix: if the yield normal, the flow direction, the
    plastic modulus or the hardening rate carried a stray factor of 2 anywhere,
    the disagreement here would be O(1), not O(1e-13).
    """
    _shear_build(lambda t: N.mat_vm(t, "Continuum"), SHEAR_STEPS, GAMMA_END)
    codes, eps, sig = _drive_eps_sig(SHEAR_STEPS)
    assert codes == [0] * SHEAR_STEPS, f"simple-shear path failed: {codes}"

    # The rig is a pure engineering shear on slot 5 (v13 == gamma_xz).
    assert abs(eps[-1][5] - GAMMA_END) < 1e-12, f"rig is not simple shear: {eps[-1]}"
    assert np.max(np.abs(np.delete(eps[-1], 5))) < 1e-14, (
        f"rig picked up other strain components: {eps[-1]}")
    # ... and it must actually be PLASTIC, or the test measures nothing.
    tau_elastic = (N.E_VM / (2.0 * (1.0 + N.NU_VM))) * GAMMA_END
    assert sig[-1][5] < 0.5 * tau_elastic, (
        f"shear path never yielded (tau={sig[-1][5]}, elastic {tau_elastic})")

    ref = VMO.radial_return(eps, N.E_VM, N.NU_VM, N.SY_VM, N.H_VM)
    scale = float(np.max(np.abs(ref)))
    err = float(np.max(np.abs(sig - ref))) / scale
    assert err <= 1e-10, (
        f"Backward_Euler disagrees with the closed-form radial return by "
        f"{err:.3e} relative on a simple-shear path -- the ADR-94 B5 shear-slot "
        f"convention has regressed.\nopensees={sig[-1]}\noracle  ={ref[-1]}")


@pytest.mark.t0m
def test_C1_vm_hardening_rate_is_shear_path_independent(vm_available):
    """ADR-94 B5, second half.  ``AllASDHardeningFunctions``' equivalent plastic
    strain rate is ``H*sqrt(2/3 m_ij m_ij)``; with an ENGINEERING-shear ``m`` the
    contraction must weight the shear slots by 1/2.  The pre-94c ``m.dot(m)``
    weighted them by 1, so the isotropic hardening modulus a shear path saw was
    not the one an axial path saw.

    Measured as: the yield stress implied by the committed shear stress after a
    known plastic strain must match the closed form to round-off.  (The oracle
    comparison above already covers this; this test names the number so a
    failure says WHICH half broke.)
    """
    _shear_build(lambda t: N.mat_vm(t, "Continuum"), SHEAR_STEPS, GAMMA_END)
    codes, eps, sig = _drive_eps_sig(SHEAR_STEPS)
    assert codes == [0] * SHEAR_STEPS

    c = VMO.SQRT_2_over_3
    G = N.E_VM / (2.0 * (1.0 + N.NU_VM))
    tau = sig[-1][5]
    # closed form: ||s|| = tau*sqrt(2) = c*sy, sy = sy0 + H*c*dLambda,
    # and tau = G*(gamma - 2*dLambda*tau/||s||*... ) -- solved by the oracle;
    # here just require the deviatoric norm to sit ON the current surface.
    s_norm = VMO.tensor_norm_stress(VMO.deviator(sig[-1]))
    dlam = (G * GAMMA_END * math.sqrt(2.0) - c * N.SY_VM) / (2.0 * G + c * c * N.H_VM)
    sy = N.SY_VM + dlam * N.H_VM * c
    assert abs(s_norm - c * sy) / (c * sy) <= 1e-10, (
        f"committed shear state is not on the hardened yield surface: "
        f"||s||={s_norm}, c*sy={c * sy}")


# ===========================================================================
# 2. B4 -- Drucker-Prager gradients agree with a numerical tangent
# ===========================================================================
@pytest.fixture(scope="module")
def dp_available():
    try:
        _shear_build(lambda t: HB.mat_dp(t), 2, 1.0e-5)
    except Exception as exc:                                  # pragma: no cover
        pytest.skip(f"DruckerPrager_YF / LadrunoBrick unavailable: {exc}")
    ops.wipe()


def _mat_dp_tangent(tag, tangent):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL", HB.IV_DP,
        "Begin_Model_Parameters",
        "YoungsModulus", HB.DP_E, "PoissonsRatio", HB.DP_NU,
        "DP_xi_c", HB.DP_XI_C, "DP_eta", HB.DP_ETA, "DP_etabar", HB.DP_ETABAR,
        "TensorLinearHardeningParameter", 0.0,
        "ScalarLinearHardeningParameter", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "DP_cohesion", 0.0,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", "Backward_Euler",
        "tangent_type", tangent,
        "End_Integration_Options")


@pytest.mark.t0m
def test_C2_dp_continuum_and_numerical_tangents_agree(dp_available):
    """ADR-94 B4.  ``Continuum`` builds the elasto-plastic operator from the
    ANALYTICAL ``df_dsigma_ij``; ``Numerical_Algorithmic_FirstOrder``
    differentiates the return map.  On a sheared (non-degenerate) plastic path
    the converged stresses must agree.

    They are not a redundant pair: before wp/94c the analytical gradient was
    2x too large on the three normal slots -- 0.971 relative error against a
    central difference -- so the analytical tangent described a different
    surface from the one the return map was iterating on.
    """
    out = {}
    for tangent in ("Continuum", "Numerical_Algorithmic_FirstOrder"):
        _shear_build(lambda t, g=tangent: _mat_dp_tangent(t, g), 10, 4.0e-3)
        codes, _, sig = _drive_eps_sig(10)
        assert codes == [0] * 10, f"DP shear path failed with {tangent}: {codes}"
        out[tangent] = sig

    a = out["Continuum"]
    b = out["Numerical_Algorithmic_FirstOrder"]
    err = float(np.max(np.abs(a - b))) / float(np.max(np.abs(b)))
    assert err <= 1e-6, (
        f"the analytical (Continuum) and numerically differentiated tangents "
        f"converge to different Drucker-Prager stresses ({err:.3e} relative) -- "
        f"the ADR-94 B4 gradient fix has regressed.\n{a[-1]}\n{b[-1]}")


# ===========================================================================
# 3. B4 -- the apex return is live
# ===========================================================================
@pytest.mark.t0m
def test_C3_dp_hydrostatic_tension_returns_to_the_apex(dp_available):
    """ADR-94 B4.  Pure hydrostatic TENSION driven to 2x the apex volumetric
    strain.  The deviator is exactly zero the whole way, so there is no flank to
    return to: the only admissible stress is the apex itself,
    ``sigma = (xi_c/eta) * I``.

    Before wp/94c the ``Backward_Euler`` apex call site was commented out, so
    ``check_apex_region`` (a ``return false`` stub in ``DruckerPrager_YF``
    anyway) was consulted and its answer discarded, and the flank return map ran
    past ``sqrt(J2) = 0``.
    """
    K = HB.DP_E / (3.0 * (1.0 - 2.0 * HB.DP_NU))
    ev = (HB.DP_P_APEX / K) * 2.0
    R._build(lambda t: HB.mat_dp(t), 10, ev, ev, ev)
    codes, hist = R._drive(10)

    assert codes == [0] * 10, (
        f"the hydrostatic-tension apex path no longer completes: {codes}")
    last = np.array(hist[-1], dtype=float)
    assert np.all(np.isfinite(last)), f"non-finite stress at the apex: {last}"

    # on the apex: hydrostatic, at p = xi_c/eta
    assert np.max(np.abs(last[3:])) <= 1e-8 * HB.DP_P_APEX, (
        f"apex stress is not hydrostatic (shear slots): {last}")
    assert np.max(np.abs(last[:3] - last[:3].mean())) <= 1e-8 * HB.DP_P_APEX, (
        f"apex stress is not hydrostatic (normal slots): {last}")
    p = float(last[:3].mean())
    assert abs(p - HB.DP_P_APEX) <= 1e-6 * HB.DP_P_APEX, (
        f"apex pressure {p} != xi_c/eta = {HB.DP_P_APEX}")

    # and f(sigma) <= tol: f = sqrt(J2) + eta*p - xi_c, sqrt(J2) == 0 here
    f = HB.DP_ETA * p - HB.DP_XI_C
    assert f <= 1e-6 * HB.DP_XI_C, f"committed apex stress is inadmissible: f={f}"


@pytest.mark.t0m
def test_C3_dp_hydrostatic_compression_is_untouched_by_the_apex(dp_available):
    """Control.  The apex is a TENSION vertex; hydrostatic compression of the
    same magnitude must stay elastic, i.e. the live apex branch must not fire on
    the compression side."""
    K = HB.DP_E / (3.0 * (1.0 - 2.0 * HB.DP_NU))
    ev = -(HB.DP_P_APEX / K) * 2.0
    R._build(lambda t: HB.mat_dp(t), 10, ev, ev, ev)
    codes, hist = R._drive(10)
    assert codes == [0] * 10, f"hydrostatic compression failed: {codes}"
    last = np.array(hist[-1], dtype=float)
    assert np.all(np.isfinite(last))
    p = float(last[:3].mean())
    # eps_v = 3*ev on the unit cube (each face is prescribed ev), so the elastic
    # answer is p = K*eps_v = 3*K*ev.  f = eta*p - xi_c < 0 there, i.e. elastic.
    p_elastic = 3.0 * K * ev
    assert abs(p - p_elastic) <= 1e-9 * abs(p_elastic), (
        f"hydrostatic compression is no longer elastic: p={p} != {p_elastic}")


# ===========================================================================
# 4. M5 -- f_relative_tol
# ===========================================================================
def _mat_mc_rel(tag, scale, ftol=1.0e-6, rtol=None, nit=100, strict=1):
    extra = [] if rtol is None else ["f_relative_tol", rtol]
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", R.IV_MC,
        "Begin_Model_Parameters",
        "YoungsModulus", R.MC_E * scale, "PoissonsRatio", R.MC_NU,
        "MC_phi", R.MC_PHI, "MC_c", R.MC_C * scale, "MC_psi", R.MC_PHI,
        "MC_ds", 0.0, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "f_absolute_tol", ftol,
        *extra,
        "n_max_iterations", nit,
        "strict_convergence", strict,
        "End_Integration_Options")


# 1e-7 * (c*cos(phi)) is 8.66e-6 in kPa and 8.66e-3 in Pa: the SAME relative
# tightness, and in both cases looser than the 1e-6 absolute the kPa run already
# passes at, so both unit systems must complete.
F_REL = 1.0e-7


@pytest.mark.t0m
def test_C4_f_relative_tol_default_is_inert():
    """Not passing the option and passing ``f_relative_tol 0`` must produce
    bit-identical histories: the default is OFF, and OFF means the absolute
    tolerance alone, exactly as before wp/94c."""
    ez = 0.01
    try:
        R._build(lambda t: _mat_mc_rel(t, 1.0, strict=0), 20, 0.0, 0.0, ez)
    except Exception as exc:                                  # pragma: no cover
        pytest.skip(f"MohrCoulomb_YF / LadrunoBrick unavailable: {exc}")
    codes_a, hist_a = R._drive(20)
    R._build(lambda t: _mat_mc_rel(t, 1.0, rtol=0.0, strict=0), 20, 0.0, 0.0, ez)
    codes_b, hist_b = R._drive(20)
    assert codes_a == codes_b, f"{codes_a} vs {codes_b}"
    assert hist_a == hist_b, (
        "`f_relative_tol 0` is not byte-identical to omitting the option")


@pytest.mark.t0m
def test_C4_f_relative_tol_makes_the_verdict_unit_independent():
    """ADR-94 M5, the fix.  The SAME physical problem in two unit systems (E and
    c scaled by ``UNIT_GAP``, strains identical) with ``strict_convergence 1``.

    With ``f_relative_tol`` OFF the fork's one fail-loud switch renders a
    verdict on the UNIT SYSTEM: the reference-unit run completes 20/20 and the
    x1e9 run is refused on step 1 (pinned in ``test_adr94_redblue_numerics``,
    which owns the ``UNIT_GAP`` constant this test reuses).  With it ON, both
    complete, and their stresses differ by exactly the unit factor.
    """
    ez = 0.01
    try:
        R._build(lambda t: _mat_mc_rel(t, 1.0, rtol=F_REL), 20, 0.0, 0.0, ez)
    except Exception as exc:                                  # pragma: no cover
        pytest.skip(f"MohrCoulomb_YF / LadrunoBrick unavailable: {exc}")
    codes_kpa, hist_kpa = R._drive(20)
    R._build(lambda t: _mat_mc_rel(t, R.UNIT_GAP, rtol=F_REL), 20, 0.0, 0.0, ez)
    codes_pa, hist_pa = R._drive(20)

    assert codes_kpa == [0] * 20, f"kPa run with f_relative_tol: {codes_kpa}"
    assert codes_pa == [0] * 20, (
        f"the x{R.UNIT_GAP:.0e} run (identical physics) is STILL refused with "
        f"f_relative_tol = {F_REL}: {codes_pa}")

    # and the two runs are the same physics: stresses differ by exactly x1000
    a = np.array(hist_kpa[-1]) * R.UNIT_GAP
    b = np.array(hist_pa[-1])
    assert float(np.max(np.abs(a - b))) <= 1e-6 * float(np.max(np.abs(b))), (
        f"kPa x1000 != Pa: {a} vs {b}")
