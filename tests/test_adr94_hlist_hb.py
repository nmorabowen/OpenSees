"""ADR-94 H10 -- Hoek-Brown drift (jaabell/ASDP vs this tree) + the
DruckerPrager apex-stub, pinned with real material-point runs.

See ``Ladruno_implementation/_adr94_hb_drift.md`` for the function-by-function
diff analysis this file measures, and ``_adr94_hlist_R1C.md`` for the H10
verdict. Summary of what these tests pin:

1. HoekBrown_YF (our tree, ``e65e89203``) computes the yield surface with a
   two-branch ``if (arg > 0) ... else ...`` split (``arg = mb*sigma3/sigci +
   s``, sigma3 the geo-frame minor principal stress). jaabell's newer
   ``ASDP`` branch (``60d9b9b23``) replaced this with a smooth composite
   ``max(f_shear, f_tension)`` that clamps ``arg`` before ``pow`` and is
   provably continuous at the apex. Our two branches are NOT continuous at
   ``arg = 0``: the standard branch already reports strong violation well
   before ``arg`` reaches 0, but once the trial stress crosses into the
   ``arg <= 0`` region the formula jumps to ``sigma1 - sigma3 - sigci*s``,
   which is admissible up to a MUCH larger apparent tensile capacity
   (measured: ``sigci*s`` = 587.19 kPa, ``mb`` = 2.397 times the textbook
   tensile strength ``sigma_t = -s*sigci/mb`` = -245.02 kPa). Measured
   consequence: a uniaxial tension path locks onto the wrong (else-branch)
   plateau and then the analysis STALLS (``analyze() == -3``,
   "PLASTIC INCONSISTENCY - ELASTIC STEP!" printed by the H7 fallback)
   instead of cleanly yielding near the textbook tensile strength.
2. A compression path never reaches ``arg <= 0``, so both trees agree there
   (this file's compression/triaxial tests double as that parity check via
   the closed-form ``arg > 0`` formula, which is untouched by the drift).
3. DruckerPrager_YF::CHECK_APEX_REGION is a stub (``return false``) while
   ``yf_has_apex<DruckerPrager_YF>`` is declared true; and the
   Backward_Euler apex-return call site in ASDPlasticMaterial3D.h
   (``if constexpr (yf_has_apex<...>) if (check_apex_region(...)) { ... }``)
   has an entirely commented-out body, so apex handling is a no-op for
   EVERY yield function today, not just DP. Measured consequence: a pure
   hydrostatic-tension path through the DP apex (``p = xi_c/eta`` in the raw
   ``sigma.meanStress()`` convention the code actually uses -- the header
   comment calling it "positive in compression" does not match the
   arithmetic) produces **NaN stress that is COMMITTED with
   ``analyze() == 0``** -- a silent NaN corruption, not merely an unhandled
   apex.

DRIVER: single LadrunoBrick unit cube (propagates material return codes,
unlike stdBrick -- see LEDGER_quirks), 1/8-symmetry restraints, sp-prescribed
normal strains via LoadControl, default Backward_Euler + Secant tangent,
``system UmfPack`` (never FullGeneral -- a fully/near-fully prescribed
driver crashes it, see the ADR-94 plan Sec. 8 traps).
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---------------------------------------------------------------------------
# single-element cube driver (1/8-symmetry), shared by both materials
# ---------------------------------------------------------------------------
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 8: (1, 0, 0)}      # node 7 carries no fix
_XFACE = (2, 3, 6, 7)          # x = 1
_YFACE = (3, 4, 7, 8)          # y = 1
_ZFACE = (5, 6, 7, 8)          # z = 1


def _build(mat_fn, nsteps, eps_x, eps_y=None, eps_z=None):
    """eps_y/eps_z = None leaves that face group UNCONSTRAINED (free lateral
    traction, i.e. a true uniaxial-STRESS path); a number sp-prescribes it
    (a confined/triaxial path)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    for t, m in _FIX.items():
        ops.fix(t, *m)
    mat_fn(1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in _XFACE:
        ops.sp(n, 1, eps_x)
    if eps_y is not None:
        for n in _YFACE:
            ops.sp(n, 2, eps_y)
    if eps_z is not None:
        for n in _ZFACE:
            ops.sp(n, 3, eps_z)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    # UmfPack, NOT FullGeneral (ADR-94 plan Sec. 8: FullGeneral hard-crashes
    # a fully/near-fully prescribed N-small driver).
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


def _drive(nsteps):
    """Advance up to nsteps; stop at the first non-zero analyze() code.

    Returns (codes, hist): every analyze() code up to and including the
    first failure, and the committed GP-1 stress (6 comps) after each
    SUCCESSFUL step. A refused/failed step contributes a code but no row.
    """
    codes, hist = [], []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        ops.eleResponse(1, "forces")           # set the lazy strain/stress
        hist.append(list(ops.eleResponse(1, "stresses"))[0:6])
    return codes, np.array(hist) if hist else np.zeros((0, 6))


def _constructible(mat_fn):
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in _CUBE.items():
            ops.node(t, *map(float, c))
        mat_fn(1)
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


# ===========================================================================
# 1. Hoek-Brown -- realistic rock parameters (Hoek & Brown 2018 GSI formulas,
#    computed exactly as HoekBrown_Utils.h)
# ===========================================================================
HB_SIGMA_CI = 50000.0     # 50 MPa unconfined compressive strength (kPa)
HB_MI = 10.0              # intact-rock constant (e.g. a fine sandstone)
HB_GSI = 60.0             # Geological Strength Index (fair rock mass)
HB_D = 0.0                # undisturbed

HB_MB = HB_MI * math.exp((HB_GSI - 100.0) / (28.0 - 14.0 * HB_D))
HB_S = math.exp((HB_GSI - 100.0) / (9.0 - 3.0 * HB_D))
HB_A = 0.5 + (1.0 / 6.0) * (math.exp(-HB_GSI / 15.0) - math.exp(-20.0 / 3.0))
HB_SIGMA_T = -HB_S * HB_SIGMA_CI / HB_MB          # textbook tensile strength
HB_UCS = HB_SIGMA_CI * HB_S ** HB_A               # rock-mass UCS (arg>0 branch)
# The ELSE-branch plateau this tree's YF actually enforces in tension:
# yf = sigma1 - sigma3 - sigci*s = sigma_xx - sigci*s = 0  (uniaxial stress,
# sigma1=0, sigma3=-sigma_xx in the geo frame).
HB_ELSE_PLATEAU = HB_SIGMA_CI * HB_S

HB_E, HB_NU = 5.0e7, 0.25          # ~50 GPa rock elastic modulus, kPa

IV_HB = "BackStress(NullHardeningTensorFunction):"


def mat_hb(tag):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "HoekBrown_YF", "HoekBrown_PF", "LinearIsotropic3D_EL", IV_HB,
        "Begin_Model_Parameters",
        "YoungsModulus", HB_E, "PoissonsRatio", HB_NU,
        "HB_sigci", HB_SIGMA_CI, "HB_mb", HB_MB, "HB_s", HB_S, "HB_a", HB_A,
        "HB_mb_psi", HB_MB, "HB_ds", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
    )


@pytest.fixture(scope="module")
def hb_available():
    if not _constructible(lambda t: mat_hb(t)):
        pytest.skip("ASDPlasticMaterial3D / HoekBrown_YF / LadrunoBrick not "
                    "available in this build")


@pytest.mark.t0m
def test_H10_hb_tension_locks_onto_wrong_plateau_then_stalls(hb_available):
    """(i) Uniaxial-STRESS tension path (lateral faces free) driven to 3x the
    textbook tensile strain ``|sigma_t|/E``.

    EXPECTED (if the yield surface were continuous, as jaabell's composite
    is): the material yields near ``sigma_xx ~ |HB_SIGMA_T|`` = 245.0 kPa and
    the analysis proceeds cleanly (strain-controlled, perfectly plastic).

    OBSERVED on this tree: the standard (``arg > 0``) branch already reports
    strong violation well before ``arg`` reaches 0 (so a continuous surface
    would have capped growth much earlier), but the local Newton keeps
    correcting using the CURRENT branch's gradient and the committed stress
    instead climbs to the ELSE-branch's own zero, ``sigci*s`` = 587.19 kPa --
    2.40x (= HB_MB) the textbook tensile strength -- and then the analysis
    FAILS to converge (rc -3) right at that plateau, with
    "PLASTIC INCONSISTENCY - ELASTIC STEP!" printed by the unguarded H7
    fallback (``dLambda + deltaLambda < 0``). This is what the port to
    jaabell's smooth composite would change.
    """
    eps_t_end = 3.0 * (-HB_SIGMA_T) / HB_E
    _build(lambda t: mat_hb(t), 60, eps_t_end)
    codes, hist = _drive(60)

    assert len(hist) > 0, "no step committed at all -- driver regressed"
    last_sigma_xx = float(hist[-1, 0])

    # Pin the wrong (else-branch) plateau, not the textbook tensile strength.
    assert last_sigma_xx == pytest.approx(HB_ELSE_PLATEAU, rel=2.0e-3), (
        f"HB tension plateau drifted: measured {last_sigma_xx:.4f} kPa, "
        f"expected the else-branch zero sigci*s = {HB_ELSE_PLATEAU:.4f} kPa "
        f"(textbook tensile strength is only {-HB_SIGMA_T:.4f} kPa)")
    assert last_sigma_xx > 2.0 * (-HB_SIGMA_T), (
        "measured plateau is no longer well above the textbook tensile "
        "strength -- the discontinuity this test pins may have been fixed; "
        "if intentional, this test should be updated to CONFIRM the fix")

    # The analysis stalls (rc -3) at/after the plateau instead of continuing
    # to accept strain-controlled steps.
    assert codes[-1] == -3, (
        f"expected the driver to stall (-3) right after the else-branch "
        f"plateau; got codes={codes}. The discontinuity's downstream "
        f"symptom (H7's PLASTIC INCONSISTENCY fallback) may have changed.")


@pytest.mark.t0m
def test_H10_hb_uniaxial_compression_matches_closed_form_ucs(hb_available):
    """Compression never crosses ``arg <= 0``, so both trees agree here: the
    committed stress must land exactly on the closed-form rock-mass UCS,
    ``sigci * s**a`` (the ``arg > 0`` branch is untouched by the drift)."""
    eps_c_end = -1.5 * HB_UCS / HB_E
    _build(lambda t: mat_hb(t), 60, eps_c_end)
    codes, hist = _drive(60)

    assert len(hist) > 0
    last_sigma_xx = float(hist[-1, 0])
    assert last_sigma_xx == pytest.approx(-HB_UCS, rel=1.0e-6), (
        f"HB compression plateau drifted from the closed-form UCS: "
        f"measured {last_sigma_xx:.4f} vs -sigci*s**a = {-HB_UCS:.4f} kPa")


@pytest.mark.t0m
def test_H10_hb_triaxial_compression_stays_admissible(hb_available):
    """(ii) A genuinely triaxial (confined) compression path: axial strain
    dominant, lateral strain confined to 15% of it (both compressive).
    Confinement raises the HB envelope quickly (power-law growth, a ~ 0.5),
    so this path stays comfortably elastic/admissible -- it never
    approaches the arg=0 branch switch the tension test above exercises,
    which is itself the point: compression is where our tree and
    jaabell's agree bit-for-bit (both reduce to the same arg>0 formula)."""
    k = 0.15
    eps_ax = -1.3 * 0.0005237680297507399   # sized to approach (not reach)
                                              # the k=0.15 closed-form yield
    _build(lambda t: mat_hb(t), 60, eps_ax, k * eps_ax, k * eps_ax)
    codes, hist = _drive(60)

    assert codes == [0] * len(codes) and len(hist) == 60, (
        f"triaxial compression path no longer converges cleanly: "
        f"codes={codes}")

    s = hist[-1]
    # admissibility: recompute the standard (arg>0) HB yield function from
    # the committed stress (compression-positive geo frame, HoekBrown_YF.h
    # arithmetic) and require it stays safely negative (elastic).
    sigma_geo = sorted([-s[0], -s[1], -s[2]])
    sigma3, sigma2, sigma1 = sigma_geo
    arg = HB_MB * sigma3 / HB_SIGMA_CI + HB_S
    yf = sigma1 - sigma3 - HB_SIGMA_CI * (max(arg, 0.0) ** HB_A)
    assert yf < -0.5 * HB_UCS, (
        f"triaxial compression path drifted much closer to yield than "
        f"measured (yf={yf:.2f}); the sizing above assumed comfortably "
        f"elastic margin")


# ===========================================================================
# 2. DruckerPrager -- CHECK_APEX_REGION stub, H10's second half
# ===========================================================================
DP_E, DP_NU = 1.0e6, 0.25
DP_XI_C, DP_ETA, DP_ETABAR = 1000.0, 0.3, 0.1     # cohesion-like, friction
                                                    # slope, dilation slope
DP_P_APEX = DP_XI_C / DP_ETA                       # 3333.33 kPa

IV_DP = ("BackStress(TensorLinearHardeningFunction):"
         "DP_cohesion(ScalarLinearHardeningFunction):")


def mat_dp(tag):
    """Null-hardening in effect: both linear-hardening slopes are 0, so
    BackStress and DP_cohesion never evolve (perfectly plastic)."""
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


@pytest.fixture(scope="module")
def dp_available():
    if not _constructible(lambda t: mat_dp(t)):
        pytest.skip("ASDPlasticMaterial3D / DruckerPrager_YF / LadrunoBrick "
                    "not available in this build")


@pytest.mark.t0m
def test_H10_dp_apex_hydrostatic_tension_commits_nan(dp_available):
    """Pure hydrostatic-TENSION path (all three normal strains equal and
    growing) driven to 2x the apex strain ``lambda_apex = P_APEX / (3K)``,
    K the elastic bulk modulus -- so it passes straight through
    ``p = xi_c/eta`` with zero deviatoric stress the whole way (the exact
    degenerate case CHECK_APEX_REGION exists for).

    ``DruckerPrager_YF::CHECK_APEX_REGION`` is a stub returning ``false``
    (``DruckerPrager_YF.h:105-111``, "Implement!!!"), and even where a YF's
    ``check_apex_region`` DOES return true (HoekBrown), the
    ``ASDPlasticMaterial3D.h`` Backward_Euler call site's corrective body is
    entirely commented out (~2093-2161) -- apex handling is a no-op for
    every YF today.

    MEASURED CONSEQUENCE (not merely "unhandled"): the committed stress
    turns to **NaN partway through the ramp, and ``analyze()`` keeps
    returning 0 (success)** -- a silent NaN corruption, worse than a clean
    refusal.
    """
    K = DP_E / (3.0 * (1.0 - 2.0 * DP_NU))
    lam_apex = DP_P_APEX / (3.0 * K)
    lam_end = 2.0 * lam_apex

    _build(lambda t: mat_dp(t), 40, lam_end, lam_end, lam_end)
    codes, hist = _drive(40)

    assert len(hist) > 0, "no step committed -- driver regressed"

    nan_rows = np.where(np.any(np.isnan(hist), axis=1))[0]
    assert len(nan_rows) > 0, (
        f"expected the hydrostatic-tension apex path to commit NaN stress "
        f"on this build; got a clean history (max|sigma|="
        f"{np.nanmax(np.abs(hist)):.4g}). If CHECK_APEX_REGION or the "
        f"Backward_Euler apex call site were fixed, this test should be "
        f"rewritten to assert a CORRECT apex return instead of NaN.")

    first_nan = int(nan_rows[0])
    # Every code up to and including the NaN row's step must read as
    # "success" -- that IS the silent-corruption defect being pinned.
    assert all(c == 0 for c in codes[:first_nan + 1]), (
        f"NaN row {first_nan} was reached but analyze() did not report "
        f"success (0) throughout: codes={codes[:first_nan + 1]}")

    # Before the NaN, the material must have been tracking the elastic
    # hydrostatic path (p grows linearly with strain, no plasticity in an
    # associated-flow degenerate direction until the apex is reached).
    pre_nan = hist[:first_nan]
    if len(pre_nan):
        p_pre = pre_nan[:, 0]              # hydrostatic: sxx=syy=szz
        assert np.all(np.diff(p_pre) > 0), (
            "pre-NaN hydrostatic-tension history is not monotonically "
            "increasing -- the reproducer's assumption (clean elastic "
            "tracking up to the apex) no longer holds")


@pytest.mark.t0m
def test_H10_dp_apex_hydrostatic_compression_stays_admissible(dp_available):
    """Sanity twin of the tension probe: hydrostatic COMPRESSION of the same
    magnitude must stay elastic (the DP cone here opens toward compression,
    apex only on the tension side) -- confirms the NaN above is not a
    generic hydrostatic-loading bug."""
    K = DP_E / (3.0 * (1.0 - 2.0 * DP_NU))
    lam_apex = DP_P_APEX / (3.0 * K)
    lam_end = -2.0 * lam_apex

    _build(lambda t: mat_dp(t), 40, lam_end, lam_end, lam_end)
    codes, hist = _drive(40)

    assert codes == [0] * len(codes) and len(hist) == 40
    assert np.all(np.isfinite(hist)), "unexpected NaN on the compression side"
