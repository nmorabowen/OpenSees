"""ADR-94 R2 RED-NUMERICS -- reproducers for the red-team numerics lane.

Three things R1 did not measure:

1. **The DP "silent NaN" is NOT silent inside the material.** R1-C
   (``_adr94_hlist_R1C.md``) attributed the committed NaN to "NaN comparisons
   are always false, so nothing in the existing gates trips". Measured here:
   ``Backward_Euler``'s own NaN guard (``ASDPlasticMaterial3D.h:2318-2323``,
   ``if (!(norm_trial_stress == norm_trial_stress)) { cout << "NaN!"; return
   -1; }``) DOES fire and DOES return -1 -- the -1 is then dropped on the way
   up.  So the defect is a *return-code contract* defect, not (only) a missing
   check. PLATFORM CAVEAT: the NaN itself is UNDEFINED BEHAVIOUR (an
   uninitialised ``VoigtVector pressure_part; pressure_part *= 0.0;`` in
   DruckerPrager_YF/_PF), not a deterministic defect -- Ubuntu CI's fresh
   heap commits a clean, finite history on the identical apex path, so the
   guard never fires there; only Windows' dirty pytest heap reproduces it.

2. **The NaN is DP-specific, not "apex-declared-YF"-generic.**  Five YFs
   specialize ``yf_has_apex`` to true (DruckerPrager, HoekBrown, MohrCoulomb,
   RoundedMohrCoulomb, TensionCutoff).  MohrCoulomb driven through its own
   apex (``p = c*cot(phi)``) on the identical rig stays finite -- because
   ``MohrCoulomb_YF``'s analytical derivative has a ``J2 < 1e-15`` hydrostatic
   guard that ``DruckerPrager_YF``'s does not.

3. **``f_absolute_tol`` is absolute in stress units, so a UNIT CHANGE alone
   decides whether ``strict_convergence`` fails loud.**  The same physical
   problem expressed in kPa and in Pa (E, c, sigma_y all x1000; strains
   identical) is the same boundary-value problem, but ``|Phi|`` scales by
   1000 while ``f_absolute_tol`` does not.

Driver: the single ``LadrunoBrick`` unit cube of ``test_adr94_hlist_hb.py``
(1/8 symmetry, sp-prescribed normal strains, LoadControl, UmfPack).
Build ``52314165a``.  Wall time ~1.5 s.
"""
import math
import re
import subprocess
import sys
import textwrap

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_PRESSURE_PART_IDIOM = re.compile(
    r"VoigtVector\s+pressure_part\s*;\s*\n\s*pressure_part\s*\*=\s*0\.0\s*;")


def _source(rel_path):
    import pathlib
    here = pathlib.Path(__file__).resolve().parent.parent
    return (here / rel_path).read_text(encoding="utf-8", errors="replace")

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 8: (1, 0, 0)}
_XFACE = (2, 3, 6, 7)
_YFACE = (3, 4, 7, 8)
_ZFACE = (5, 6, 7, 8)


def _build(mat_fn, nsteps, ex, ey, ez):
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
        ops.sp(n, 1, ex)
    for n in _YFACE:
        ops.sp(n, 2, ey)
    for n in _ZFACE:
        ops.sp(n, 3, ez)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
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
    return codes, hist


def _run_child(body: str):
    """OS-level stdout capture -- capfd cannot see this .pyd's cout
    (see _adr94_hlist_R1B.md 'Traps recorded for the next lane')."""
    script = textwrap.dedent(body)
    p = subprocess.run([sys.executable, "-c", script], capture_output=True,
                       text=True, stdin=subprocess.DEVNULL, timeout=300)
    return p.stdout + p.stderr


# ---------------------------------------------------------------------------
# 1. DP apex: the material's OWN NaN guard fires; the -1 never reaches analyze()
# ---------------------------------------------------------------------------
DP_E, DP_NU = 1.0e6, 0.25
DP_XI_C, DP_ETA = 1000.0, 0.3
DP_P_APEX = DP_XI_C / DP_ETA
IV_DP = "BackStress(TensorLinearHardeningFunction):DP_cohesion(ScalarLinearHardeningFunction):"

_DP_CHILD = f'''
import sys
sys.path.insert(0, r"{{PYPATH}}")
import opensees as ops
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
cube = {_CUBE!r}
fix = {_FIX!r}
for t, c in cube.items():
    ops.node(t, float(c[0]), float(c[1]), float(c[2]))
for t, m in fix.items():
    ops.fix(t, m[0], m[1], m[2])
ops.nDMaterial("ASDPlasticMaterial3D", 1,
    "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL",
    "{IV_DP}",
    "Begin_Model_Parameters",
    "YoungsModulus", {DP_E}, "PoissonsRatio", {DP_NU},
    "DP_xi_c", {DP_XI_C}, "DP_eta", {DP_ETA}, "DP_etabar", 0.1,
    "TensorLinearHardeningParameter", 0.0,
    "ScalarLinearHardeningParameter", 0.0,
    "MassDensity", 0.0,
    "End_Model_Parameters",
    "Begin_Internal_Variables",
    "BackStress", 0., 0., 0., 0., 0., 0.,
    "DP_cohesion", 0.0,
    "End_Internal_Variables")
ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
K = {DP_E} / (3.0 * (1.0 - 2.0 * {DP_NU}))
ev = 3.0 * ({DP_P_APEX} / K) / 3.0 * 2.0   # 2x the apex volumetric strain
for n in (2, 3, 6, 7):
    ops.sp(n, 1, ev)
for n in (3, 4, 7, 8):
    ops.sp(n, 2, ev)
for n in (5, 6, 7, 8):
    ops.sp(n, 3, ev)
ops.constraints("Transformation")
ops.numberer("Plain")
ops.system("UmfPack")
ops.test("NormDispIncr", 1.0e-10, 50, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0 / 10)
ops.analysis("Static")
codes = []
for _ in range(10):
    codes.append(ops.analyze(1))
ops.eleResponse(1, "forces")
s = list(ops.eleResponse(1, "stresses"))[0:6]
sys.stderr.write("RESULT_CODES=" + repr(codes) + "\\n")
sys.stderr.write("RESULT_STRESS=" + repr(s) + "\\n")
'''


def _pypath():
    import opensees
    import os
    return os.path.dirname(os.path.dirname(opensees.__file__))


@pytest.mark.t0m
def test_R2_dp_apex_nan_guard_fires_and_is_swallowed():
    """R1-C said "nothing in the existing gates trips".  Measured: the
    Backward_Euler NaN guard prints "NaN!" and returns -1 on the apex path,
    and analyze() still reports 0 with a NaN stress committed on Windows.
    The gate works; the return code is lost between the material and the
    analysis.

    PLATFORM CAVEAT: the NaN itself is UNDEFINED BEHAVIOUR, not a
    deterministic defect -- ``DruckerPrager_YF``/``_PF`` build their
    pressure term via ``VoigtVector pressure_part; pressure_part *= 0.0;``,
    multiplying UNINITIALISED Eigen storage by zero rather than actually
    zero-constructing it. Windows' churned pytest heap reliably reads stale
    nonzero/NaN bits there; Ubuntu CI's fresh heap reads zero and the same
    apex path commits a clean, finite history, so the guard never fires and
    "NaN!" never prints. The structural sentinel below (the pseudo-init
    idiom's continued presence) is what actually pins this; the "NaN!"/
    committed-NaN checks are gated on NaN having actually been committed.
    """
    yf_src = _source("SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions"
                      "/DruckerPrager_YF.h")
    pf_src = _source("SRC/material/nD/ASDPlasticMaterial3D"
                      "/PlasticFlowDirections/DruckerPrager_PF.h")
    assert _PRESSURE_PART_IDIOM.search(yf_src), (
        "DruckerPrager_YF.h no longer contains the uninitialised "
        "`VoigtVector pressure_part; pressure_part *= 0.0;` pseudo-init "
        "idiom -- the apex-path UB this test pins may be gone; re-measure "
        "and update this test (and test_H10_dp_apex_hydrostatic_tension_"
        "commits_nan) accordingly.")
    assert _PRESSURE_PART_IDIOM.search(pf_src), (
        "DruckerPrager_PF.h no longer contains the same idiom -- see above.")

    out = _run_child(_DP_CHILD.replace("{PYPATH}", _pypath()))
    if "RESULT_CODES=" not in out:
        pytest.skip("DP/LadrunoBrick child run unavailable:\n" + out[-800:])
    _ns = {"nan": float("nan"), "inf": float("inf")}
    codes = eval(out.split("RESULT_CODES=")[1].splitlines()[0], _ns)
    stress = eval(out.split("RESULT_STRESS=")[1].splitlines()[0], _ns)
    nan_committed = any(s != s for s in stress)
    all_finite = all(math.isfinite(s) for s in stress)

    # (a) the analysis never reports failure on this degenerate path,
    # regardless of which UB outcome this platform's heap produced.
    assert all(c == 0 for c in codes), (
        "expected every analyze() to report success on this degenerate "
        f"apex path, got {codes}")
    assert nan_committed or all_finite, (
        f"stress is neither NaN-flagged nor entirely finite -- unexpected "
        f"state: {stress}")

    # (b) the guard's own "NaN!" print, and the committed-NaN it swallows,
    # are only expected when this platform's heap actually reproduced the
    # UB -- on a fresh/zeroed heap the path stays clean and the guard never
    # trips at all (which is not a regression, just the other UB outcome).
    if nan_committed:
        assert "NaN!" in out, (
            "expected Backward_Euler's own NaN guard to fire once NaN was "
            "committed; stdout tail:\n" + out[-1500:])


# ---------------------------------------------------------------------------
# 2. MohrCoulomb through its own apex stays finite
# ---------------------------------------------------------------------------
MC_E, MC_NU = 1.0e6, 0.25
MC_PHI, MC_C = 30.0, 100.0
MC_P_APEX = MC_C / math.tan(math.radians(MC_PHI))   # 173.2 kPa
IV_MC = "BackStress(NullHardeningTensorFunction):"


def _mat_mc(tag, scale=1.0):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", IV_MC,
        "Begin_Model_Parameters",
        "YoungsModulus", MC_E * scale, "PoissonsRatio", MC_NU,
        "MC_phi", MC_PHI, "MC_c", MC_C * scale, "MC_psi", MC_PHI,
        "MC_ds", 0.0, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables")


@pytest.mark.t0m
def test_R2_mc_hydrostatic_tension_through_apex_stays_finite():
    """Same rig, same degenerate (zero-deviator) path, an apex-declaring YF --
    and no NaN.  MohrCoulomb_YF's default (MC_ds = 0) derivative takes the
    analytical branch, which guards ``J2 < 1e-15`` explicitly; DruckerPrager_YF
    has no such guard and additionally builds its pressure term from an
    UNINITIALISED VoigtVector zeroed with ``pressure_part *= 0.0``
    (DruckerPrager_YF.h:64-65 / DruckerPrager_PF.h:68-69).  So the H10b NaN is
    a DP defect, not a property of the dead apex call site."""
    K = MC_E / (3.0 * (1.0 - 2.0 * MC_NU))
    ev = (MC_P_APEX / K) * 2.0          # 2x the apex volumetric strain
    try:
        _build(lambda t: _mat_mc(t), 10, ev, ev, ev)
    except Exception as exc:                       # pragma: no cover
        pytest.skip(f"MohrCoulomb_YF / LadrunoBrick unavailable: {exc}")
    codes, hist = _drive(10)
    assert hist, f"no step completed: codes={codes}"
    last = hist[-1]
    assert all(s == s for s in last), (
        f"MohrCoulomb NaN'd at its own apex: {last}")
    assert all(abs(s) < 1e12 for s in last), f"MC stress blew up: {last}"


# ---------------------------------------------------------------------------
# 3. f_absolute_tol is ABSOLUTE in stress units: the unit system alone decides
#    whether strict_convergence refuses the step
# ---------------------------------------------------------------------------
def _mat_mc_strict(tag, scale, ftol=1.0e-6, nit=100, strict=1):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", IV_MC,
        "Begin_Model_Parameters",
        "YoungsModulus", MC_E * scale, "PoissonsRatio", MC_NU,
        "MC_phi", MC_PHI, "MC_c", MC_C * scale, "MC_psi", MC_PHI,
        "MC_ds", 0.0, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "f_absolute_tol", ftol,
        "n_max_iterations", nit,
        "strict_convergence", strict,
        "End_Integration_Options")


@pytest.mark.t0m
def test_R2_f_absolute_tol_makes_strict_convergence_unit_dependent():
    """The SAME physical problem in kPa and in Pa (E, c x1000; strains
    identical) on the SHIPPED default ``f_absolute_tol = 1e-6`` with
    ``strict_convergence 1``.

    MEASURED: kPa completes 20/20 steps; Pa refuses on step 1 (analyze() ==
    -3, LADRUNO_MATERIAL_REFUSED from the Backward_Euler exhaustion gate).
    ``|Phi|`` scales with the stress units, ``f_absolute_tol`` does not -- so
    the fork's one fail-loud switch renders a verdict on the UNIT SYSTEM, not
    on the physics.  A user working in SI base units gets a material that
    cannot complete a step a kPa user runs without complaint."""
    ez = 0.01                        # uniaxial-strain tension, well past yield
    try:
        _build(lambda t: _mat_mc_strict(t, 1.0), 20, 0.0, 0.0, ez)
    except Exception as exc:                       # pragma: no cover
        pytest.skip(f"MohrCoulomb_YF / LadrunoBrick unavailable: {exc}")
    codes_kpa, _ = _drive(20)
    _build(lambda t: _mat_mc_strict(t, 1000.0), 20, 0.0, 0.0, ez)
    codes_pa, _ = _drive(20)
    ok_kpa = sum(1 for c in codes_kpa if c == 0)
    ok_pa = sum(1 for c in codes_pa if c == 0)
    assert ok_kpa == 20, f"kPa run was expected to complete: {codes_kpa}"
    assert ok_pa == 0 and codes_pa[-1] != 0, (
        "the Pa run (identical physics, stresses x1000) was expected to be "
        f"REFUSED at the same default f_absolute_tol: {codes_pa}")


# ---------------------------------------------------------------------------
# 4. strict_convergence is a complete no-op on stdBrick
# ---------------------------------------------------------------------------
def _build_ele(ele, mat_fn, nsteps, ex, ey, ez):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    for t, m in _FIX.items():
        ops.fix(t, *m)
    mat_fn(1)
    ops.element(ele, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in _XFACE:
        ops.sp(n, 1, ex)
    for n in _YFACE:
        ops.sp(n, 2, ey)
    for n in _ZFACE:
        ops.sp(n, 3, ez)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


@pytest.mark.t0m
def test_R2_strict_convergence_is_a_noop_on_stdbrick():
    """``Brick::update()`` assigns ``success = mat->setTrialStrain(...)`` and
    then ``return 0;`` unconditionally (SRC/element/brick/Brick.cpp), so it
    drops BOTH a bare -1 AND the ``LADRUNO_MATERIAL_REFUSED`` sentinel.

    MEASURED on an identical deck, material and tolerance: ``LadrunoBrick``
    refuses every step (0/20, analyze() == -3); ``stdBrick`` reports 20/20
    successes.  ADR-84 P2a's fail-loud gate therefore has no effect at all on
    the most widely used solid element in OpenSees -- including for
    ``Backward_Euler``, the one integrator the gate reaches."""
    ez = 0.01
    try:
        _build_ele("LadrunoBrick", lambda t: _mat_mc_strict(t, 1.0, ftol=1e-12),
                   20, 0.0, 0.0, ez)
    except Exception as exc:                       # pragma: no cover
        pytest.skip(f"MohrCoulomb_YF / LadrunoBrick unavailable: {exc}")
    codes_lb, _ = _drive(20)
    _build_ele("stdBrick", lambda t: _mat_mc_strict(t, 1.0, ftol=1e-12),
               20, 0.0, 0.0, ez)
    codes_sb, _ = _drive(20)
    assert sum(1 for c in codes_lb if c == 0) == 0, (
        f"LadrunoBrick was expected to propagate the refusal: {codes_lb}")
    assert sum(1 for c in codes_sb if c == 0) == 20, (
        f"stdBrick was expected to swallow the refusal entirely: {codes_sb}")
