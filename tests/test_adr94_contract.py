"""ADR-94 R5 -- host-element contract for ASDPlasticMaterial3D failure returns.

Scope (`Ladruno_implementation/94_asdplastic_review_plan.md` sec. 5, R5):
propagate `LADRUNO_MATERIAL_REFUSED` from EVERY integrator (today the sentinel
is only ever RETURNED from two `Backward_Euler` sites, 2198/2350), and check
revert semantics under a failed step / `ops.reset()` / `LadrunoBeginAugment`-
`LadrunoEndAugment`.  Method notes for the next lane, `_adr94_contract.md`.

Host contract, MEASURED (this file), consistent with the mission brief and
`_adr94_hlist_R1B.md`/H8:
  * ``TenNodeTetrahedron::update()`` -- ``success += materialPointers[i]->
    setTrialStrain(...)`` -- PROPAGATES every nonzero code (sentinel or bare).
  * ``Brick::update()`` (``stdBrick``) -- assigns ``success`` and never reads
    it; ``update()`` unconditionally ``return 0`` -- SWALLOWS every failure.
  * ``LadrunoBrick`` -- every ``setTrialStrain`` call site (6 of them:
    ``updateHypo``'s SSP-centroid and per-GP loops, ``formEAStrue``'s
    condensed/full loops) checks ONLY ``== LADRUNO_MATERIAL_REFUSED``.  A bare
    ``-1`` (FE singular-tangent, FE_sub, BE singular-tangent, BE NaN-guard,
    BE_LS exhaustion, RK45_old max-iter, ME NaN/max-iter, RK45 NaN/max-iter --
    everything except the two BE sentinel sites) does NOT match, so
    ``LadrunoBrick`` never prints its own refusal warning and never returns a
    failure code for that call -- SWALLOWED at the element level.  Measured
    below (``test_R5_...silent_at_element_level``): the GLOBAL Newton still
    catches the corrupted state on this rig via non-convergence, but that is
    an accident of residual magnitude, not a guarantee -- H5/H7 already show
    cases where a bad local state converges globally with rc=0.
"""
import os
import re
import subprocess
import sys

import numpy as np
import pytest

from _testbed import ops

import test_asdplastic_mctc as M
from test_adr94_hlist_mechanical import _run_child, _TESTS_DIR, _DIST_DIR  # noqa: F401

pytestmark = [pytest.mark.zone_a]

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_ELEMENT_DIR = os.path.join(_REPO, "SRC", "element")


def _src(*parts):
    with open(os.path.join(_REPO, *parts), encoding="utf-8", errors="replace") as fh:
        return fh.read()


# ===========================================================================
# structural propagation table
# ===========================================================================
def test_R5_ladrunobrick_checks_only_the_sentinel():
    """CONFIRMED (structural).  Every ``materialPointers[...]->setTrialStrain``
    call in ``LadrunoBrick.cpp`` is guarded by ``== LADRUNO_MATERIAL_REFUSED``,
    never by a blanket ``< 0`` or ``!= 0``.  6 call sites: ``updateHypo``'s SSP
    centroid (single-point) and per-GP loop, and ``formEAStrue``'s condensed
    (single-point) and full-integration loops (two loops, invoked from
    ``tang_flag`` 0 and !=0 branches).  A future 7th call site that checks
    ``< 0`` instead would (correctly) start propagating bare -1 too; this test
    would need updating, not silently pass.
    """
    src = _src("SRC", "element", "ladrunoBrick", "LadrunoBrick.cpp")
    calls = re.findall(r"setTrialStrain\([^)]*\)\s*(?:\n\s*)?==\s*LADRUNO_MATERIAL_REFUSED",
                        src)
    assert len(calls) == 6, (
        f"expected exactly 6 sentinel-only setTrialStrain checks in "
        f"LadrunoBrick.cpp, found {len(calls)} -- the call-site count "
        f"changed; re-verify R5's host contract before trusting this test.")
    # and NONE of them is spelled as a blanket failure test anywhere else
    blanket = re.findall(r"setTrialStrain\([^)]*\)\s*(?:\n\s*)?(?:<\s*0|!=\s*0)(?!\w)",
                          src)
    assert not blanket, (
        f"LadrunoBrick.cpp now also tests a setTrialStrain return with a "
        f"blanket comparison ({blanket!r}) -- the sentinel-only swallow may "
        f"have been fixed; re-verify R5.")


def test_R5_tennodetetrahedron_sums_raw_return_codes():
    """CONFIRMED (structural).  ``TenNodeTetrahedron::update()`` accumulates
    ``success += materialPointers[i]->setTrialStrain(strain)`` over every
    Gauss point and returns ``success`` -- ANY nonzero code (sentinel value
    -33086, or a bare -1/-3 from any integrator) survives to the caller.
    """
    src = _src("SRC", "element", "tetrahedron", "TenNodeTetrahedron.cpp")
    assert re.search(r"success\s*\+=\s*materialPointers\[i\]->setTrialStrain",
                      src), (
        "TenNodeTetrahedron.cpp no longer accumulates setTrialStrain's return "
        "code via 'success +='; re-verify R5's host contract.")


def test_R5_stdbrick_swallows_the_return_code():
    """CONFIRMED (structural), same defect as ADR-84 P2a's own finding.
    ``Brick::update()`` (the ``stdBrick`` element) writes
    ``success = materialPointers[i]->setTrialStrain(strain)`` inside the Gauss
    loop and the function unconditionally ``return 0;`` afterward --
    ``success`` is assigned and never read. Every non-zero return code from
    every integrator (sentinel or bare) is swallowed on this host.
    """
    src = _src("SRC", "element", "brick", "Brick.cpp")
    m = re.search(
        r"int\s+Brick::update\s*\(.*?\)\s*\{(.*?)\n\}\n", src, re.S)
    assert m, "could not locate Brick::update(void) body -- source layout changed"
    body = m.group(1)
    assert re.search(r"success\s*=\s*materialPointers\[i\]->setTrialStrain",
                      body), "Brick::update no longer assigns 'success' from setTrialStrain"
    # the last live (non-comment) statement in the body must be 'return 0;'
    tail = body.rstrip()
    live_lines = [ln for ln in tail.splitlines() if ln.strip()
                  and not ln.strip().startswith("//")]
    assert live_lines[-1].strip() == "return 0;" or "return 0 ;" in live_lines[-1], (
        f"Brick::update()'s last statement is now {live_lines[-1]!r}, not a "
        f"bare 'return 0;' -- the swallow may have been fixed; re-verify R5 "
        f"(and update ADR-84 P2a's own note, which relies on this too).")


# ===========================================================================
# runtime pin -- Backward_Euler_LineSearch exhaustion (bare -1, 2588), which
# H8 already measured fails 18/20 steps on this exact rig (a proven,
# reproducible bare-`-1` failure -- unlike the FE-NaN / ME-max-iter /
# RK45-dT_min branches, which this review could NOT force at runtime within
# budget on VonMises/MohrCoulomb decks; see the method notes in
# `_adr94_contract.md`).
# ===========================================================================
_TET = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
        5: (.5, 0, 0), 6: (.5, .5, 0), 7: (0, .5, 0), 8: (0, 0, .5),
        9: (.5, 0, .5), 10: (0, .5, .5)}
_TET_TOP = (4, 8, 9, 10)
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 7: (0, 0, 0), 8: (1, 0, 0)}


def _mat_mc_bels(tag):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
        "Begin_Model_Parameters", "YoungsModulus", M.E, "PoissonsRatio", M.NU,
        "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
        "MassDensity", 0.0, "End_Model_Parameters",
        "Begin_Internal_Variables", "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        "Begin_Integration_Options", "integration_method",
        "Backward_Euler_LineSearch", "n_max_iterations", 100,
        "End_Integration_Options",
    )


@pytest.mark.t0m
def test_R5_be_linesearch_exhaustion_propagates_on_tet():
    """CONFIRMED, runtime.  Reproduces H8's own measurement: on the ADR-84 MC
    tet leg, ``Backward_Euler_LineSearch`` fails after a couple of clean
    steps (its inner split-loop cannot chain substeps and returns a bare -1).
    ``TenNodeTetrahedron`` propagates that -1 into a global non-convergence
    (analyze() rc != 0).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
    _mat_mc_bels(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, -0.02)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    codes = [ops.analyze(1) for _ in range(20)]
    assert codes[0] == 0, f"expected step 1 to converge cleanly; codes={codes}"
    assert any(c != 0 for c in codes), (
        f"Backward_Euler_LineSearch completed all 20 steps with rc=0 on the "
        f"H8 reproducer rig -- this integrator may have been fixed; "
        f"re-verify H8/R5 before trusting this test. codes={codes}")


@pytest.mark.t0m
def test_R5_be_linesearch_exhaustion_silent_at_element_level_on_ladrunobrick():
    """CONFIRMED, runtime.  Same material/integrator on a ``LadrunoBrick``
    host: the global Newton STILL fails (rc != 0 on every step) because the
    corrupted/unintegrated stress unbalances the residual -- but
    ``LadrunoBrick``'s own refusal message
    ("...the material REFUSED the trial strain...", emitted only when the
    sentinel-only check at 1030-1035/1079-1082/1110-1113/1181-1184/1826-1829/
    3318-3321 fires) never appears, because a bare -1 does not match
    ``LADRUNO_MATERIAL_REFUSED``. The element itself never notices the
    refusal; only the accident of a bad-enough residual makes ``analyze()``
    report failure here. Run in a CHILD process (native cout/opserr on this
    build is not visible to capfd -- see ``_run_child``'s docstring).
    """
    script = (
        "import sys; sys.path.insert(0, r'" + _TESTS_DIR + "')\n"
        "from _testbed import ops\n"
        "import test_asdplastic_mctc as M\n"
        "_CUBE = {1:(0,0,0),2:(1,0,0),3:(1,1,0),4:(0,1,0),"
        "5:(0,0,1),6:(1,0,1),7:(1,1,1),8:(0,1,1)}\n"
        "_FIX = {1:(1,1,1),2:(0,1,1),3:(0,0,1),4:(1,0,1),"
        "5:(1,1,0),6:(0,1,0),7:(0,0,0),8:(1,0,0)}\n"
        "ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 3)\n"
        "[ops.node(i, *map(float, c)) for i, c in _CUBE.items()]\n"
        "[ops.fix(i, *m) for i, m in _FIX.items()]\n"
        "ops.nDMaterial('ASDPlasticMaterial3D', 1, 'MohrCoulomb_YF', "
        "'MohrCoulomb_PF', 'LinearIsotropic3D_EL', M.IV, "
        "'Begin_Model_Parameters', 'YoungsModulus', M.E, 'PoissonsRatio', "
        "M.NU, 'MC_phi', M.PHI, 'MC_c', M.C, 'MC_psi', M.PSI, 'MC_ds', 0.0, "
        "'MassDensity', 0.0, 'End_Model_Parameters', "
        "'Begin_Internal_Variables', 'BackStress', 0.,0.,0.,0.,0.,0., "
        "'End_Internal_Variables', 'Begin_Integration_Options', "
        "'integration_method', 'Backward_Euler_LineSearch', "
        "'n_max_iterations', 100, 'End_Integration_Options')\n"
        "ops.element('LadrunoBrick', 1, 1,2,3,4,5,6,7,8, 1)\n"
        "ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)\n"
        "for i in (2,3,6,7): ops.sp(i, 1, -0.02)\n"
        "ops.constraints('Transformation'); ops.numberer('Plain')\n"
        "ops.system('UmfPack'); ops.test('NormDispIncr', 1e-8, 100, 0)\n"
        "ops.algorithm('Newton'); ops.integrator('LoadControl', 0.05)\n"
        "ops.analysis('Static')\n"
        "codes = [ops.analyze(1) for _ in range(20)]\n"
        "print('CODES', codes)\n"
    )
    proc = _run_child(script)
    assert "CODES" in proc.stdout, (
        f"child process did not complete (stdout={proc.stdout!r}, "
        f"stderr={proc.stderr!r})")
    codes_line = [ln for ln in proc.stdout.splitlines() if ln.startswith("CODES")][0]
    codes = eval(codes_line.split("CODES", 1)[1].strip())
    assert any(c != 0 for c in codes), (
        f"expected the global Newton to fail at least once on this "
        f"LadrunoBrick BE_LS rig; codes={codes} -- re-verify before trusting "
        f"the 'silent at element level' half of this test.")
    combined = proc.stdout + proc.stderr
    assert "the material REFUSED the trial strain" not in combined, (
        "LadrunoBrick now prints its own refusal warning for a bare -1 "
        "(non-sentinel) failure -- the sentinel-only swallow may have been "
        "widened to a blanket check; update R5's host contract.")


# ===========================================================================
# H4 consequence -- ops.reset() leaves the material's committed state
# inconsistent with the (correctly) reset geometry
# ===========================================================================
def _mat_mc_plain(tag):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
        "Begin_Model_Parameters", "YoungsModulus", M.E, "PoissonsRatio", M.NU,
        "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
        "MassDensity", 0.0, "End_Model_Parameters",
        "Begin_Internal_Variables", "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
    )


def _tet_build_plain():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
    _mat_mc_plain(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, -0.02)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")


def _tet_stress():
    ops.eleResponse(1, "forces")
    return np.array(list(ops.eleResponse(1, "stresses"))[0:6])


@pytest.mark.t0m
def test_H4_reset_leaves_material_state_inconsistent_with_geometry():
    """CONFIRMED, runtime -- the R1-B rig could not observe this because it
    only checked the "not implemented" print, not a numeric consequence.

    ``Domain::revertToStart()`` correctly resets nodal trial/committed
    displacements to zero (``ops.reset()`` reports success), but
    ``ASDPlasticMaterial3D::revertToStart()`` is a documented no-op (prints
    "not implemented", returns -1, ignored by ``OPS_resetModel()``) -- the
    material's own Commit/TrialStress and Commit/TrialStrain survive
    untouched. Querying the element's stress IMMEDIATELY after reset (no
    ``analyze()`` call) exercises ``TenNodeTetrahedron``'s self-healing
    "stresses" response (H4's caveat): it recomputes strain from the just-
    zeroed nodal displacement and feeds that into the STALE material, whose
    ``CommitStrain``/``CommitStress`` still encode the pre-reset plastic
    state. The result is neither (a) ~zero, which a fully-reset model would
    report, nor (b) the pre-reset committed stress unchanged -- it is a THIRD,
    inconsistent value, proving geometry and material disagree about "reset"
    having happened.
    """
    _tet_build_plain()
    for _ in range(10):
        assert ops.analyze(1) == 0
    sig_before = _tet_stress()
    assert np.max(np.abs(sig_before)) > 1.0, "expected a genuinely plastic pre-reset state"

    ops.reset()
    sig_at_reset = _tet_stress()

    tol = 1.0e-6 * np.max(np.abs(sig_before))
    assert np.max(np.abs(sig_at_reset)) > tol, (
        f"post-reset queried stress is ~zero ({sig_at_reset}) -- "
        f"revertToStart() may have been implemented; re-verify H4.")
    assert np.max(np.abs(sig_at_reset - sig_before)) > tol, (
        f"post-reset queried stress ({sig_at_reset}) is unchanged from the "
        f"pre-reset committed stress ({sig_before}) -- the self-heal query "
        f"path may have changed; re-verify H4's caveat before trusting this "
        f"test.")


@pytest.mark.t0m
def test_H4_cutback_after_forced_global_failure_is_not_bitwise_reproducible():
    """MEASURED, runtime.  A step that fails to converge GLOBALLY (an
    impossible ``NormDispIncr`` budget, not a material refusal) leaves a dirty
    TRIAL stress that ``revertToLastCommit()``'s no-op body (H4) never clears
    -- ``StaticAnalysis::analyze()`` calls ``Domain::revertToLastCommit()``
    itself, and H4 already pins that this is a no-op at the material level.

    Measured consequence on this rig: retrying the SAME step (same
    LoadControl increment) after restoring a sane test tolerance, then
    running the remaining identical steps, reaches a FINAL committed stress
    that is close to -- but not bitwise/1e-12 equal to -- a reference run that
    never attempted the failing step. The gap (~6e-9 relative) sits at the
    level of the ``NormDispIncr`` 1e-8 convergence tolerance itself, so this
    probe cannot separate "H4's broken revert leaked a dirty trial state"
    from ordinary Newton-truncation noise; it is recorded as a measured,
    inconclusive data point (see ``_adr94_contract.md``), NOT as a confirmed
    corruption. The clean, confirmed H4 consequence is the reset() test above.
    """
    _tet_build_plain()
    ref_codes = [ops.analyze(1) for _ in range(20)]
    assert ref_codes == [0] * 20, f"reference run did not converge cleanly: {ref_codes}"
    sig_ref = _tet_stress()

    _tet_build_plain()
    codes_a = [ops.analyze(1) for _ in range(10)]
    assert codes_a == [0] * 10
    ops.test("NormDispIncr", 1.0e-14, 1, 0)   # impossible budget -> global fail
    rc_fail = ops.analyze(1)
    assert rc_fail != 0, "expected the impossible-tolerance step to fail"
    ops.test("NormDispIncr", 1e-8, 100, 0)    # restore; retry + remaining steps
    codes_b = [ops.analyze(1) for _ in range(10)]
    assert codes_b == [0] * 10
    sig_recovered = _tet_stress()

    diff = float(np.max(np.abs(sig_recovered - sig_ref)))
    scale = float(np.max(np.abs(sig_ref)))
    assert diff > 1.0e-9, (
        f"recovered vs reference stress is now bitwise-identical (diff="
        f"{diff:.3e}) -- either H4 was fixed or this probe's premise changed; "
        f"re-verify before trusting this test as 'inconclusive'.")
    assert diff / scale < 1.0e-6, (
        f"recovered vs reference stress differs by {diff:.3e} (relative "
        f"{diff / scale:.3e} of scale {scale:.3e}) -- this is well beyond "
        f"ordinary Newton-tolerance noise; H4's broken revert may be "
        f"corrupting the cutback recovery more than previously measured -- "
        f"escalate this finding instead of treating it as inconclusive.")


# ===========================================================================
# LadrunoBeginAugment / LadrunoEndAugment -- generic Domain flags, NOT
# contact-gated (they do not check for a contact handler at all)
# ===========================================================================
def test_R5_augment_commands_are_generic_not_contact_gated():
    """RECORDED (structural), per the R5 scope note ("if the commands are
    contact-only and refuse without a contact handler, record that and
    skip"). They are NOT contact-only: ``OPS_LadrunoBeginAugment``/
    ``OPS_LadrunoEndAugment`` (``OpenSeesOutputCommands.cpp``) simply call
    ``Domain::setContactAugmenting(true/false)`` unconditionally -- no check
    for a registered contact handler, idempotent, no args, cannot fail
    (always ``return 0``). Calling them on a plain (contact-free) ASDP model
    is legal and a no-op for the ASDP material's own revert/commit path: the
    flag only changes whether ``Domain::commit()`` fires recorders/bumps
    commitTag, which is orthogonal to how a Ladruno-material integrator
    computes or reports failure. Runtime confirmation below; deeper
    interaction with H12's per-tag cout diagnostics under repeated
    held-load commits is NOT measured here (R5's ½-session budget) and is
    flagged as a follow-up in ``_adr94_contract.md``.
    """
    src = _src("SRC", "interpreter", "OpenSeesOutputCommands.cpp")
    m = re.search(r"int OPS_LadrunoBeginAugment\(\)\s*\{(.*?)\n\}", src, re.S)
    assert m, "could not locate OPS_LadrunoBeginAugment -- source layout changed"
    body = m.group(1)
    assert "ContactHandler" not in body and "getContact" not in body, (
        "OPS_LadrunoBeginAugment now references a contact handler -- it may "
        "have been made contact-gated; re-verify this finding.")

    _tet_build_plain()
    assert ops.analyze(1) == 0
    ops.ladrunoBeginAugment()
    rc = ops.analyze(1)
    ops.ladrunoEndAugment()
    assert rc == 0, (
        f"a plain (contact-free) MC tet step failed inside an augment "
        f"bracket (rc={rc}) -- the commands may have gained a precondition; "
        f"re-verify.")
