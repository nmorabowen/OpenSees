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


@pytest.mark.t0m
def test_R5_be_linesearch_exhaustion_propagates_on_tet():
    """FIXED by wp/94a (ADR-94 M7) -- the runtime pin is now a REFUSAL pin.

    R5 measured this on `52314165a`: on the ADR-84 MC tet leg,
    ``Backward_Euler_LineSearch`` failed after a couple of clean steps (its
    inner split-loop cannot chain substeps and returned a bare -1), and
    ``TenNodeTetrahedron`` propagated that -1 into a global non-convergence.
    It was the ONE bare-`-1` mode this review could force reliably, which is
    why R5 built the whole host-contract table on it.

    wp/94a closes both halves at once: the split-loop exhaustion now returns
    ``LADRUNO_MATERIAL_REFUSED`` instead of a bare -1 (so every host, not just
    the tet, would hear it), and the parser refuses the integrator outright, so
    no deck can reach the site at all.  The test keeps its name and now pins
    the refusal: the material is not created and the command fails.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t_, c in _TET.items():
        ops.node(t_, *map(float, c))
    with pytest.raises(Exception):
        _mat_mc_bels(1)
    ops.wipe()

    # control: the SAME deck on plain Backward_Euler still builds and still
    # completes the leg -- the refusal is specific to the broken integrator,
    # not a regression of the MC material itself.
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t_, c in _TET.items():
        ops.node(t_, *map(float, c))
    for t_ in (1, 2, 3, 5, 6, 7):
        ops.fix(t_, 1, 1, 1)
    for t_ in _TET_TOP:
        ops.fix(t_, 1, 1, 0)
    _mat_mc_plain(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t_ in _TET_TOP:
        ops.sp(t_, 3, -0.02)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    codes = [ops.analyze(1) for _ in range(20)]
    assert codes == [0] * 20, (
        f"plain Backward_Euler must still complete the ADR-84 MC tet leg "
        f"20/20; codes={codes}")


@pytest.mark.t0m
def test_R5_be_linesearch_exhaustion_silent_at_element_level_on_ladrunobrick():
    """FIXED by wp/94a (ADR-94 B2 + M7).

    R5 measured, on `52314165a`, that the same BE_LS deck on a
    ``LadrunoBrick`` host made the global Newton fail on every step, but
    ``LadrunoBrick``'s own refusal warning never appeared: its six
    ``setTrialStrain`` checks compare ONLY ``== LADRUNO_MATERIAL_REFUSED``, and
    the split-loop exhaustion returned a bare -1.  The failure was visible only
    by the accident of an unbalanced residual.

    wp/94a fixes the material side (all 13 bare-`-1` sites now return the
    sentinel) and refuses the integrator at the parser.  The host side is
    deliberately UNCHANGED -- still sentinel-only, still the right contract --
    which is why the structural pin
    ``test_R5_ladrunobrick_checks_only_the_sentinel`` above must stay green.
    This test now pins that the deck cannot be built at all, in a CHILD process
    so the parser's own opserr text is captured (native output is invisible to
    capfd on this build -- see ``_run_child``'s docstring).
    """
    script = (
        "import sys; sys.path.insert(0, r'" + _TESTS_DIR + "')\n"
        "from _testbed import ops\n"
        "import test_asdplastic_mctc as M\n"
        "ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 3)\n"
        "try:\n"
        "    ops.nDMaterial('ASDPlasticMaterial3D', 1, 'MohrCoulomb_YF', "
        "'MohrCoulomb_PF', 'LinearIsotropic3D_EL', M.IV, "
        "'Begin_Model_Parameters', 'YoungsModulus', M.E, 'PoissonsRatio', "
        "M.NU, 'MC_phi', M.PHI, 'MC_c', M.C, 'MC_psi', M.PSI, 'MC_ds', 0.0, "
        "'MassDensity', 0.0, 'End_Model_Parameters', "
        "'Begin_Internal_Variables', 'BackStress', 0.,0.,0.,0.,0.,0., "
        "'End_Internal_Variables', 'Begin_Integration_Options', "
        "'integration_method', 'Backward_Euler_LineSearch', "
        "'n_max_iterations', 100, 'End_Integration_Options')\n"
        "    print('RESULT created')\n"
        "except Exception as exc:\n"
        "    print('RESULT refused')\n"
    )
    proc = _run_child(script)
    assert "RESULT" in proc.stdout, (
        f"child process did not complete (stdout={proc.stdout!r}, "
        f"stderr={proc.stderr!r})")
    assert "RESULT refused" in proc.stdout, (
        f"Backward_Euler_LineSearch was accepted by the parser -- ADR-94 M7's "
        f"refusal has regressed; stdout={proc.stdout[-1200:]!r}")

    combined = proc.stdout + proc.stderr
    assert "Backward_Euler_LineSearch" in combined and "ADR-94" in combined, (
        "the refusal must name the integrator AND cite ADR-94 so a user can "
        "find out why; got:\n" + combined[-1500:])


# ===========================================================================
# H4 consequence -- ops.reset() leaves the material's committed state
# inconsistent with the (correctly) reset geometry
# ===========================================================================


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
    """FIXED by wp/94b.  ``ops.reset()`` is a real reset now.

    Before the fix: ``Domain::revertToStart()`` correctly zeroed nodal
    trial/committed displacements, but ``ASDPlasticMaterial3D::revertToStart()``
    was a no-op (printed "not implemented", returned -1, discarded by both
    ``Domain::revertToStart()`` and ``OPS_resetModel()``), so the material's
    own Commit/Trial stress and strain survived.  Querying the element's stress
    immediately after ``reset()`` then exercised ``TenNodeTetrahedron``'s
    self-healing "stresses" response against the STALE material and produced a
    THIRD, inconsistent number -- neither zero nor the pre-reset value.
    wp/94b restores stress, strain, plastic strain, the internal variables
    (from a snapshot taken on the instance's first ``setTrialStrain``), the
    ``first_step`` / ``stress_set_externally`` flags and the tangent, so the
    queried stress after ``reset()`` is zero.
    """
    _tet_build_plain()
    for _ in range(10):
        assert ops.analyze(1) == 0
    sig_before = _tet_stress()
    assert np.max(np.abs(sig_before)) > 1.0, "expected a genuinely plastic pre-reset state"

    ops.reset()
    sig_at_reset = _tet_stress()

    tol = 1.0e-9 * np.max(np.abs(sig_before))
    assert np.max(np.abs(sig_at_reset)) < tol, (
        f"post-reset queried stress is {sig_at_reset}, not ~zero (pre-reset "
        f"was {sig_before}) -- revertToStart() is not restoring the material; "
        f"the wp/94b fix may have been reverted.")


@pytest.mark.t0m
def test_H4_cutback_after_forced_global_failure_recovers_within_newton_tolerance():
    """FIXED by wp/94b.  A step that fails to converge GLOBALLY (an impossible
    ``NormDispIncr`` budget, not a material refusal) used to leave a dirty
    TRIAL stress that ``revertToLastCommit()``'s commented-out body never
    cleared, even though ``StaticAnalysis::analyze()`` calls
    ``Domain::revertToLastCommit()`` on the way out.

    Before the fix this probe was INCONCLUSIVE: the recovered stress differed
    from a never-failed reference by ~6e-9 relative on Windows and ~5.4e-13 on
    Linux CI (where the TenNodeTetrahedron's self-healing "stresses" query
    erased the gap outright), which could not be separated from ordinary
    Newton-truncation noise.  With the revert implemented, the retry restarts
    from exactly the committed state, so the recovery is identical to
    floating-point round-off (measured ~1.9e-21 relative on Windows -- a
    single-ULP-level difference from operation reordering, ~10^11 tighter
    than the pre-fix ~6e-9 noise floor) rather than exactly bitwise; the
    tolerance below is a real gate, not the old coarse guard.
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
    assert diff / scale < 1e-9, (
        f"recovered vs reference stress differs by {diff:.3e} (relative "
        f"{diff / scale:.3e} of scale {scale:.3e}); with wp/94b's revert the "
        f"retry must restart from (numerically) the committed state -- this "
        f"is ~10^9 tighter than the pre-fix ~6e-9 noise floor, so a failure "
        f"here means the revert is incomplete again.")


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
