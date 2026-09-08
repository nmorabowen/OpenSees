"""ADR-94 `wp/94a-fail-loud` — the CORRECTED behaviour of the fix wave.

Unlike the nine `test_adr94_*.py` review files, which pin DEFECTS (so a fix
flips them red), every test here pins what `wp/94a-fail-loud` made true and must
stay green:

1. **Loud parser** (ADR-94 B1/H13) — a typo'd integration option, a typo'd model
   parameter, and a deck that never sets a required parameter each REFUSE the
   `nDMaterial` command instead of running a different material silently.
2. **Refused integrators** (ADR-94 M7/M8) — `Backward_Euler_LineSearch` and
   `Runge_Kutta_45_Error_Control_old` cannot be selected.
3. **Sentinel everywhere** (ADR-94 B2/B4) — a Drucker-Prager hydrostatic-tension
   run on `LadrunoBrick` no longer commits NaN or heap garbage.
4. **Strict mode reaches every integrator** (ADR-94 M2/B3) — a `Forward_Euler`
   deck under `strict_convergence 1` refuses rather than committing an
   inadmissible state.

Host-element rule (ADR-94 sec. 8): refusals are gated on `LadrunoBrick` or
`TenNodeTetrahedron`; `stdBrick` swallows every material return code and is
never used for a refusal assertion here.
"""
import math

import numpy as np
import pytest

from _testbed import ops

import test_asdplastic_mctc as M

pytestmark = [pytest.mark.zone_a]


# ===========================================================================
# rigs
# ===========================================================================
_TET = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
        5: (.5, 0, 0), 6: (.5, .5, 0), 7: (0, .5, 0), 8: (0, 0, .5),
        9: (.5, 0, .5), 10: (0, .5, .5)}
_TET_TOP = (4, 8, 9, 10)

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_XFACE = (2, 3, 6, 7)
_YFACE = (3, 4, 7, 8)
_ZFACE = (5, 6, 7, 8)
_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 7: (0, 0, 0), 8: (1, 0, 0)}

# Drucker-Prager: same numbers as `test_adr94_hlist_hb.py`'s apex probe.
DP_E, DP_NU = 1.0e6, 0.25
DP_XI_C, DP_ETA, DP_ETABAR = 1000.0, 0.3, 0.1
DP_P_APEX = DP_XI_C / DP_ETA
IV_DP = ("BackStress(TensorLinearHardeningFunction):"
         "DP_cohesion(ScalarLinearHardeningFunction):")

MC_PARAMS = ["YoungsModulus", M.E, "PoissonsRatio", M.NU,
             "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
             "MassDensity", 0.0]


def _mat_mc(tag, params=None, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            *(MC_PARAMS if params is None else params),
            "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


def _mat_dp(tag):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL", IV_DP,
        "Begin_Model_Parameters",
        "YoungsModulus", DP_E, "PoissonsRatio", DP_NU,
        "DP_xi_c", DP_XI_C, "DP_eta", DP_ETA, "DP_etabar", DP_ETABAR,
        "TensorLinearHardeningParameter", 0.0,
        "ScalarLinearHardeningParameter", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "DP_cohesion", 0.0,
        "End_Internal_Variables",
    )


def _fresh_model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


def _material_is_created(mat_fn):
    """True iff `mat_fn` produced a usable material.

    A REFUSED deck makes `nDMaterial` return null; depending on the front end
    that surfaces either as a raised exception or as a material tag that no
    element can bind to, so both are treated as "not created".
    """
    _fresh_model()
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    try:
        mat_fn(1)
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def mc_available():
    if not _material_is_created(lambda t: _mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb / LadrunoBrick not "
                    "available in this build")


@pytest.fixture(scope="module")
def dp_available():
    if not _material_is_created(lambda t: _mat_dp(t)):
        pytest.skip("ASDPlasticMaterial3D / DruckerPrager / LadrunoBrick not "
                    "available in this build")


def _brick_build(mat_fn, nsteps, ex, ey=None, ez=None, tol=1e-10, maxiter=50):
    _fresh_model()
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    for t, m in _FIX.items():
        ops.fix(t, *m)
    mat_fn(1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in _XFACE:
        ops.sp(n, 1, ex)
    if ey is not None:
        for n in _YFACE:
            ops.sp(n, 2, ey)
    if ez is not None:
        for n in _ZFACE:
            ops.sp(n, 3, ez)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, maxiter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


def _tet_build(mat_fn, nsteps, utop, tol=1e-8, maxiter=100):
    _fresh_model()
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
    mat_fn(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, utop)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, maxiter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


def _drive(nsteps):
    codes, hist = [], []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        ops.eleResponse(1, "forces")
        hist.append(list(ops.eleResponse(1, "stresses"))[0:6])
    return codes, (np.array(hist) if hist else np.zeros((0, 6)))


# ===========================================================================
# 1. loud parser (ADR-94 B1 / H13)
# ===========================================================================
@pytest.mark.t0m
def test_typo_in_integration_options_is_refused(mc_available):
    """FIXED by wp/94a.  ``Begin_Integration_Options`` was a run of independent
    ``if``s with no ``else``: an unrecognised token and its value were dropped
    silently, so a deck spelling ``strict_convergance`` produced byte-identical
    results to one that never asked for the flag at all (ADR-94 H13a).  The
    token now names itself on ``opserr`` and the material is NOT created."""
    assert _material_is_created(
        lambda t: _mat_mc(t, opts=["strict_convergence", 1])), (
        "control deck (correct spelling) must still build")

    assert not _material_is_created(
        lambda t: _mat_mc(t, opts=["strict_convergance", 1])), (
        "a misspelled integration option was accepted -- the `else` branch in "
        "OPS_AllASDPlasticMaterial3Ds.cpp's option chain is gone again")


@pytest.mark.t0m
def test_typo_in_model_parameter_is_refused(mc_available):
    """FIXED by wp/94a.  ``utuple_storage::setParameterByName``'s base case was
    a silent no-op, so ``MC_phii`` left ``MC_phi`` at its 0.0 default and the
    deck ran a frictionless material that still converged (ADR-94 H13b).  The
    setter now returns ``bool`` and the parser fails on ``false``."""
    bad = [p if p != "MC_phi" else "MC_phii" for p in MC_PARAMS]
    assert not _material_is_created(lambda t: _mat_mc(t, params=bad)), (
        "a misspelled model parameter was accepted -- setParameterByName's "
        "bool contract is gone again")


@pytest.mark.t0m
def test_missing_required_parameter_is_refused(mc_available):
    """FIXED by wp/94a (ADR-84 P2(e)).  An unset model parameter defaulted to
    0.0 with no warning.  After parsing, every name in ``getParameterNames()``
    except ``MassDensity`` and ``InitialP0`` must have been assigned."""
    # drop MC_c (cohesion) and its value
    i = MC_PARAMS.index("MC_c")
    missing = MC_PARAMS[:i] + MC_PARAMS[i + 2:]
    assert not _material_is_created(lambda t: _mat_mc(t, params=missing)), (
        "a deck that never set MC_c was accepted -- the required-parameter "
        "check is gone again")

    # MassDensity IS optional and must NOT trip the check
    j = MC_PARAMS.index("MassDensity")
    no_rho = MC_PARAMS[:j] + MC_PARAMS[j + 2:]
    assert _material_is_created(lambda t: _mat_mc(t, params=no_rho)), (
        "MassDensity is exempt from the required-parameter check and must "
        "still build")


@pytest.mark.t0m
@pytest.mark.parametrize("method", ["Backward_Euler_LineSearch",
                                     "Runge_Kutta_45_Error_Control_old"])
def test_broken_integrators_are_refused(mc_available, method):
    """FIXED by wp/94a (ADR-94 M7/M8).  ``Backward_Euler_LineSearch`` completes
    2/20 steps where plain ``Backward_Euler`` completes 20/20, ignores both
    ``n_max_iterations`` and ``strict_convergence``, and reports success for a
    strain the element never asked for; ``Runge_Kutta_45_Error_Control_old``'s
    drift check is an empty ``if`` body and its NaN guard calls ``exit()`` on
    the process.  Both are still compiled -- they may just not be selected."""
    assert not _material_is_created(
        lambda t: _mat_mc(t, opts=["integration_method", method])), (
        f"integration_method {method} was accepted -- its ADR-94 refusal is "
        f"gone")


@pytest.mark.t0m
def test_unknown_integration_method_name_is_refused(mc_available):
    """FIXED by wp/94a.  An unrecognised ``integration_method`` used to print a
    `cout` warning and then silently DEFAULT to Modified_Euler_Error_Control --
    a different integrator than either the one asked for or the deck's own
    default."""
    assert not _material_is_created(
        lambda t: _mat_mc(t, opts=["integration_method", "Backwards_Euler"]))
    assert not _material_is_created(
        lambda t: _mat_mc(t, opts=["tangent_type", "Concinuum"]))


# ===========================================================================
# 2. sentinel everywhere (ADR-94 B2 / B4)
# ===========================================================================
@pytest.mark.t0m
def test_dp_apex_hydrostatic_tension_never_commits_nan(dp_available):
    """FIXED by wp/94a (ADR-94 B4 + B2).  Pure hydrostatic TENSION straight
    through the Drucker-Prager apex ``p = xi_c/eta`` with an exactly zero
    deviator.  Before the fix this path built its pressure term from
    ``VoigtVector pressure_part; pressure_part *= 0.0;`` on uninitialised Eigen
    storage (``NaN*0 == NaN``) and committed NaN stress with ``analyze() == 0``
    on ``LadrunoBrick``, because the material's own NaN guard returned a bare
    ``-1`` that the element's sentinel-only compare dropped.

    Both halves are fixed: the pressure term is ``VoigtVector::Zero()``, and
    the guard returns ``LADRUNO_MATERIAL_REFUSED``.  So the invariant is now
    the strong one -- **no committed state is ever non-finite**, whatever the
    step codes are.  A refused step (non-zero code) is an acceptable outcome;
    a NaN commit is not.
    """
    K = DP_E / (3.0 * (1.0 - 2.0 * DP_NU))
    lam_end = 2.0 * (DP_P_APEX / (3.0 * K))

    _brick_build(lambda t: _mat_dp(t), 40, lam_end, lam_end, lam_end)
    codes, hist = _drive(40)

    assert np.all(np.isfinite(hist)), (
        f"a non-finite stress was COMMITTED on the DP apex path -- the "
        f"Eigen-init fix (VoigtVector::Zero) or the NaN-guard sentinel has "
        f"regressed. codes={codes}, hist={hist}")

    bad = [c for c in codes if c not in (0, -3)]
    assert not bad, (
        f"unexpected analyze() codes on the DP apex path: {codes} (0 = step "
        f"taken, -3 = global Newton gave up after the material refused)")


@pytest.mark.t0m
def test_dp_hydrostatic_compression_is_unaffected(dp_available):
    """Control for the probe above: the DP cone opens toward compression, so
    the same magnitude in compression stays elastic and must still complete
    every step with finite stress.  This is the no-regression half of the
    Eigen-init fix."""
    K = DP_E / (3.0 * (1.0 - 2.0 * DP_NU))
    lam_end = -2.0 * (DP_P_APEX / (3.0 * K))

    _brick_build(lambda t: _mat_dp(t), 40, lam_end, lam_end, lam_end)
    codes, hist = _drive(40)

    assert codes == [0] * 40 and hist.shape[0] == 40, (
        f"hydrostatic compression regressed: codes={codes}")
    assert np.all(np.isfinite(hist))


# ===========================================================================
# 3. strict mode reaches every integrator (ADR-94 M2 / B3)
# ===========================================================================
@pytest.mark.t0m
def test_strict_convergence_gates_forward_euler(mc_available):
    """FIXED by wp/94a (ADR-94 M2/H5).  ``strict_convergence`` used to be read
    ONLY inside ``Backward_Euler``; the same ``yf_val_start > yf_val_end`` =>
    "elastic, no correction" shortcut in ``Forward_Euler`` was unguarded, so
    the flag had zero effect and the run committed ``f_MC`` in the hundreds
    against a tolerance of ~0.1 (measured, `_adr94_hlist_R1B.md`).

    With the flag on, every committed state must now be admissible: either the
    step is refused (the tet propagates the material's return code, so
    ``analyze() != 0``) or ``f_MC <= tol``.  Never both zero-code and
    inadmissible.
    """
    _tet_build(lambda t: _mat_mc(t, opts=["integration_method", "Forward_Euler",
                                           "experimental_integrator", 1,  # Ladruno (ADR-97 wp/97f, D5)
                                           "strict_convergence", 1,
                                           "n_max_iterations", 100]),
               nsteps=20, utop=-0.02)
    codes, hist = [], []
    for _ in range(20):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        hist.append(_tet_stress())

    scale = 2.0 * M.C * math.cos(math.radians(M.PHI))
    if hist:
        arr = np.array(hist)
        tol = 1.0e-6 * max(scale, float(np.max(np.abs(arr))))
        fmc = np.array([M.f_mc(s) for s in arr])
        worst = float(np.max(fmc))
    else:
        tol, worst = 1.0e-6 * scale, float("-inf")

    assert worst <= tol, (
        f"Forward_Euler with strict_convergence=1 COMMITTED an inadmissible "
        f"state (max f_MC={worst:.4e} > tol {tol:.4e}) -- the strict gate on "
        f"the f-decreasing elastic exit has regressed. codes={codes}")


@pytest.mark.t0m
def test_strict_convergence_off_still_runs_forward_euler(mc_available):
    """No-regression twin: with the flag OFF the same deck behaves exactly as
    before the fix -- the strict helper is constant-false, so nothing in the
    flag-off path changed except the VALUE a failure site returns."""
    _tet_build(lambda t: _mat_mc(t, opts=["integration_method", "Forward_Euler",
                                           "experimental_integrator", 1,  # Ladruno (ADR-97 wp/97f, D5)
                                           "strict_convergence", 0,
                                           "n_max_iterations", 100]),
               nsteps=20, utop=-0.02)
    codes = [ops.analyze(1) for _ in range(20)]
    assert codes == [0] * 20, (
        f"flag-off Forward_Euler must be unchanged by wp/94a; codes={codes}")


def _tet_stress():
    ops.eleResponse(1, "forces")
    return np.array(list(ops.eleResponse(1, "stresses"))[0:6])
