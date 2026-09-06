"""ADR-92 / P1 -- pytest battery for `-implex` on `LadrunoSANISAND`.

Companion to `test_ladruno_sanisand.py` (the ADR-86 constants) and
`test_ladruno_sanisand_integrator.py` (the ADR-86b integrator seams).  This
file is written BEFORE the P1 build exists (`Ladruno_implementation/
_adr92_p1_execution_plan.md` -- W1-W9 are "syntax-clean, NOT BUILT"), against
the source read directly out of `SRC/material/nD/LadrunoSANISAND.{h,cpp}`, the
ADR-92 P0 oracle results (`_adr92_p0_oracle_results.md`), and the plan's own
section 4 record of what the C++ author decided and what the plan got wrong.

REUSE, NOT RETYPING.  `_PARAMS` / `_XY` / `_P_ATM`, the confine-first
zero-free-DOF rig (`_build_confined` / `_drive_confined` / `_stress` /
`_reldiff`), the low-p `-Pmin` rig's constants and series generator
(`_PMIN_*`, `_c_series_at`), and the recorded ADR-86b regression number
(`_RECORDED_PR_SENSITIVITY`) are all imported from the two sibling files, not
copied -- a third copy of Gorini's constants was flagged as a drift risk in
review.

POST-FIRST-RUN NOTE (2026-09-06).  This file ran against the first P1 binary
(42dbf066e) and scored 11/15 before this pass. Two failures were this file's
own deck-design bugs, not the C++, and are fixed here: `-implex` vs `-implex`
off was compared at DIFFERENT effective substep caps (`_CAP_ADEQUATE` /
`_CAP_TIGHT` now make the cap a controlled variable, everywhere an ON/OFF
pair is compared); the two refusal tests asserted `analyze() != 0` on a
ZERO-free-DOF deck, where a refusal has no equation to fail and is therefore
invisible to the return code (measured: the material printed the refusal
warning 16 times and `analyze()` still returned 0) -- both now run on
`_build_free_dof_triaxial`, a genuinely free-DOF `LadrunoBrick` cube. The
floor-clamp test's ISOCHORIC deviatoric ramp also measured clamp count == 0
against that binary; `_build_floor_seeking_deck` replaces it with a
net-DILATING ramp (see that function's docstring for why isochoric shear was
probably too mild). The tangent-identity/`revertToLastCommit` failure from
that same run is NOT this file's fault and is not touched here.

WHAT A ZERO-FREE-DOF DECK CAN AND CANNOT SHOW UNDER `-implex`.  This matters
enough to say once, centrally, rather than in every test that runs into it.
`ManzariDafalias::commitState()` (reached through `ladrunoImplexCommit()`)
ALWAYS commits the IMPLICIT return (`this->integrate()`), never the
extrapolated `sigma~` -- "IMPL-EX changes only the REPORTED stress/tangent,
never what is committed" (the same rule the LadrunoJ2 IMPL-EX oracle documents
for that material).  On a zero-free-DOF deck, `analyze(1)` always converges
trivially (there is no free-DOF residual to fail) and therefore ALWAYS
commits in the same call, so `sigma~` -- the thing a REAL boundary-value
problem's Newton loop actually iterates on -- is never observable from Python
on such a deck; only the committed (implicit) answer is.  This is exactly why
the P0 oracle measured the tangent identity in a pure-Python re-implementation
and why gate 3 here (`test_tangent_identity_...`) instead uses a genuinely
free-DOF deck with a DELIBERATELY-FAILED Newton iteration (`maxIter=1` at an
unreachable tolerance) to catch `sigma~` before the domain commits over it --
see that test's docstring for the full mechanism and the source lines it
relies on.

Two consequences that shape several tests below:
  * a zero-free-DOF deck IS the right tool for anything that only needs the
    COMMITTED answer (byte-identity, stage-0 inertness, the floor clamp's
    latched fire-count, the DB round trip) -- and it is a STRONGER tool there,
    because there is no Newton tolerance to contaminate the comparison;
  * it is the WRONG tool for anything that needs to see `sigma~` itself
    (the tangent identity) or a mid-analysis refusal's element-level return
    code propagated all the way to `analyze()` (`stdBrick` discards it --
    "stdBrick swallows material return codes", `LEDGER_quirks`); those tests
    use `LadrunoBrick`, matching `test_ladruno_sanisand_integrator.py`'s own
    rule ("ELEMENT CHOICE IS LOAD-BEARING, NOT INCIDENTAL").

NUMBERS NOT INVENTED.  Every expected number below is one of: the recorded
ADR-86b regression (`sanint._RECORDED_PR_SENSITIVITY`), an algebraic identity
read off the source (`ladrunoImplexMeasureError`'s
`||d_sigma||^2 == ||d_dev||^2 + 3*(dp)^2`, `ladrunoImplexFreezeTangent`'s
`sigma~ = sigma_n + Ce*(d_eps - f*d_eps_p)` being exactly affine in `d_eps`
for FIXED `f`/`Ce`/`d_eps_p`), or a qualitative direction taken directly from
the P0 oracle table (`_adr92_p0_oracle_results.md` section 4: IMPL-EX-A's
error grows sharply with increment size at low confinement -- "breaks at
5e-4" at p0 = 5 kPa).  Nothing here is a number read off THIS build, because
this build does not exist yet.
"""
import math
import os
import tempfile

import pytest

from _testbed import ops

import test_ladruno_sanisand as sani
import test_ladruno_sanisand_integrator as sanint


pytestmark = [pytest.mark.zone_a]

_PARAMS = sani._PARAMS
_P_ATM = sani._P_ATM
_XY = sani._XY
_EQ_TOL = sani._EQ_TOL
_SENSITIVITY_FLOOR = sani._SENSITIVITY_FLOOR


def _matvec6(flat36, v6):
    """`flat36` is a row-major-flattened 6x6 (`test_ladrunoRCConcrete_material.py`'s
    own comment on this exact response: "6x6 row-major"). Fine even if the true
    storage were column-major, because `-implex` freezes ALL THREE tangent slots
    to the same SYMMETRIC isotropic Ce (`ladrunoImplexFreezeTangent`'s own
    comment: "the delivered operator stays Ce, which keeps it symmetric") --
    for a symmetric matrix the row-major and column-major flattenings coincide."""
    assert len(flat36) == 36, ('material tangent response was not a 6x6', flat36)
    out = [0.0] * 6
    for i in range(6):
        row = flat36[6 * i:6 * i + 6]
        out[i] = sum(row[j] * v6[j] for j in range(6))
    return out


def _vnorm(v):
    return math.sqrt(sum(x * x for x in v))


# ===========================================================================
#  Gate 5 -- byte-identity with `-implex` unset (THE LOAD-BEARING GATE)
# ===========================================================================
#
#  Two independent arguments, not one, because the second is provable from the
#  source without needing to have run this build at all, while the first is
#  the traditional "the published number has not moved" regression check.

# ---------------------------------------------------------------------------
#  THE SUBSTEP CAP IS A CONTROLLED VARIABLE, NOT A CONSTANT.
#
#  Measured on this build (2026-09-06), confine-first deck, deviatoric leg: the
#  SANISAND return needs between 1000 and 5000 ModifiedEuler substeps for that
#  increment.  Capped below it the update is REFUSED -- correctly, that is what
#  ADR-86b T1 built the cap for -- and the committed answer changes:
#
#      OFF, no cap      [-3.09421, -3.09421, -22.18805]
#      OFF, cap   200   [-1.26092, -1.26092,  -9.24638]   <- the CAP moved this
#      ON,  cap   200   [51.47979, 51.47979, -108.11220]
#      OFF, cap 100000  [-3.09421, -3.09421, -22.18805]
#      ON,  cap 100000  [-3.09421, -3.09421, -22.18805]   <- BIT-IDENTICAL
#
#  So an ON/OFF comparison at a cap the deck cannot meet measures the CAP, not
#  -implex.  A first draft of this file compared `-implex -maxSubsteps 200`
#  against a control with NO cap and read the difference as an IMPL-EX defect;
#  it was not.  Every ON/OFF pair below now carries the SAME cap, and any pair
#  that means to prove agreement uses one this deck can actually meet.
# ---------------------------------------------------------------------------
_CAP_ADEQUATE = 20000   # comfortably above the 1000-5000 this deck needs
_CAP_TIGHT = 200        # deliberately too tight -- only for refusal tests


def test_implex_unset_reproduces_recorded_sensitivity():
    """`LadrunoSANISAND` with no `-implex` token anywhere reproduces the SAME
    published ADR-86b regression number as before this ADR's C++ landed.

    `sanint._RECORDED_PR_SENSITIVITY` (9.325164e-02) is
    `test_ladruno_sanisand_integrator.py`'s own gate for "ADR-86b did not move
    the recorded material-point answer" -- it exists precisely to catch an
    unrelated change perturbing this deck.  ADR-92 P1 touches exactly the same
    three files that gate already watches (`LadrunoSANISAND{,3D,PlaneStrain}`),
    so re-running it here, in the file that is ADR-92's OWN warrant, is the
    direct evidence that `ladrunoTrialUpdate()`'s off-path --
    `this->integrate(); return this->ladrunoUpdateStatus();`, verbatim, in
    that order -- is exactly what ships when `-implex` is not given.
    """
    a = sani._drive_confined('LadrunoSANISAND', 8101, sani._OPTS_VANILLA)
    b = sani._drive_confined('LadrunoSANISAND', 8102, sani._OPTS_PR0)
    rel = sani._reldiff(a, b)
    assert abs(rel - sanint._RECORDED_PR_SENSITIVITY) <= sanint._RECORDED_TOL, (
        'the confine-first deck no longer reproduces its recorded p_residual '
        'sensitivity now that the ADR-92 P1 IMPL-EX code shares a translation '
        'unit with it. ladrunoTrialUpdate() is supposed to be BYTE-IDENTICAL '
        'to the pre-ADR-92 integrate()+ladrunoUpdateStatus() pair when -implex '
        'is not given -- if this moved, something in the new code path is '
        'reachable without the flag', rel, sanint._RECORDED_PR_SENSITIVITY)


def test_implex_on_matches_off_on_a_zero_free_dof_deck():
    """THE MOST IMPORTANT TEST IN THIS FILE.

    `-implex` ON (companion scheme 1, `-implexControl` OFF) must commit the
    BIT-IDENTICAL stress to `-implex` OFF, at every step, on the confine-first
    zero-free-DOF deck -- not merely close, `==` on the list.

    WHY THIS HAS TO HOLD, argued from the source rather than merely hoped for:
    on a deck with zero free DOFs, `analyze(1)` converges in exactly one trial
    call (the residual over an empty DOF set is 0 by construction) and
    therefore always reaches `commitState()` within that same call.
    `ladrunoImplexCommit()` (`LadrunoSANISAND.cpp`, "The companion return, at
    commitState only") calls `this->integrate()` -- ManzariDafalias's OWN
    return map -- using `mEpsilon`/`mEpsilon_n` exactly as they stand, which
    are set by the ELEMENT's `setTrialStrain` BEFORE `ladrunoTrialUpdate()` is
    even invoked and are therefore IDENTICAL whether `-implex` is on or off.
    `ladrunoImplexTrial()` (the extrapolation) never writes `mEpsilon`,
    `mK`/`mG` (its own tangent-freeze helper explicitly does not touch them),
    or any committed member -- only trial `mSigma`/`mEpsilonE`/`mCe`, all of
    which `integrate()` overwrites unconditionally on entry. So the OFF path's
    single `integrate()` call and the ON path's commit-time `integrate()` call
    start from the same committed state and the same trial strain and must
    produce the same answer -- this is a claim about the CODE PATH, provable
    without a binary, and it is exactly the claim gate 5 makes.

    This is what a genuine free-DOF deck CANNOT show as cleanly: there, the
    trial `sigma~` feeds Newton and can move the iteration path (though not,
    per `test_tantype_does_not_change_the_converged_answer`'s argument, a
    force-residual-converged EQUILIBRIUM). Zero free DOFs removes that
    variable entirely, which is why this is the strongest form of the gate 5
    claim available before the build exists.
    """
    opts_off = sani._OPTS_VANILLA
    opts_on = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    off = sani._drive_confined('LadrunoSANISAND', 8103, opts_off)
    on = sani._drive_confined('LadrunoSANISAND', 8104, opts_on)
    assert off == on, (
        '-implex (companion scheme 1, -implexControl off) committed a '
        'DIFFERENT stress than -implex off on a zero-free-DOF deck, where '
        'the committed answer is supposed to be the SAME integrate() call '
        'either way. Either ladrunoImplexCommit() is not calling integrate() '
        'on the identical (mEpsilon, mEpsilon_n) pair, or something the '
        'extrapolation trial touches is leaking into the committed return',
        off, on)


# ===========================================================================
#  Gate 5 -- stage-0 inertness: gravity and a LoadControl(0.0) hold
# ===========================================================================

def test_stage0_inertness_gravity_and_hold_is_bit_identical():
    """`-implex` is inert while `mElastFlag == 0` (`ladrunoImplexActive()`):
    a gravity-stage ramp AND a `LoadControl(0.0)` re-equilibration hold, BOTH
    taken entirely at stage 0, must commit BIT-IDENTICAL stress whether
    `-implex` is on or off.

    SOURCED FROM.  `ladrunoImplexActive()` is
    `mImplexOpt.enabled && (mElastFlag != 0)` -- an explicit conjunction, not
    just "the flag" -- and `ladrunoTrialUpdate()` takes the OFF branch
    whenever it is false, setting `mImplexTrialDone = false` so
    `commitState()` also takes the base path
    (`if (!mImplexTrialDone) return ManzariDafalias::commitState();`). Since
    `mElastFlag` is a process-wide STATIC (every Manzari-family constructor
    resets it -- ADR 86 risk 3), the deck below sets stage 0 EXPLICITLY for
    this material's own tag and never flips it, so this path is exercised for
    the whole test.

    THE HOLD.  A `LoadControl(0.0)` step re-forms and re-solves with a ZERO
    load increment -- the common gravity-hold idiom used across this fork's
    own test suite (`test_stagedStrain_material.py`,
    `test_ladrunoEmbeddedNode_element.py`, et al.). On this deck's zero
    free DOFs it is a no-op mechanically, but it is the literal shape ADR-92
    P1's plan (W5) names, so it is included rather than assumed equivalent to
    the ramp alone.
    """
    opts_off = sani._OPTS_VANILLA
    opts_on = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)

    def _elastic_only_plus_hold(tag, opts):
        sani._build('LadrunoSANISAND', tag, opts)
        ops.updateMaterialStage('-material', tag, '-stage', 0)
        for step in range(sani._N_EL):
            assert ops.analyze(1) == 0, f'gravity-stage step {step + 1} failed'
        ops.integrator('LoadControl', 0.0)
        assert ops.analyze(1) == 0, 'the LoadControl(0.0) hold failed to converge'
        return sani._stress()

    off = _elastic_only_plus_hold(8105, opts_off)
    on = _elastic_only_plus_hold(8106, opts_on)
    assert off == on, (
        '-implex moved the committed stress while the material was still on '
        'the ELASTIC stage (mElastFlag == 0), where ladrunoImplexActive() is '
        'supposed to be unconditionally false. gravity and a LoadControl(0.0) '
        'hold must be bit-identical with the flag on or off', off, on)


# ===========================================================================
#  Gate 4 -- the p_min floor clamp on sigma~, and its own diagnostic
# ===========================================================================

def _build_floor_seeking_deck(tag, opts, e_conf=None, n_dev=60, lat=1.5):
    """A confine-first zero-free-DOF deck, like `sani._drive_confined_at`, but
    with the deviatoric ramp's lateral/axial ratio pushed ABOVE the isochoric
    value (`sani._C_LAT = 0.5`, net volumetric strain zero) so the path net
    DILATES instead of merely shearing at constant volume.

    WHY THIS CHANGED FROM THE FIRST DRAFT.  The first draft reused
    `sani._drive_confined_at` verbatim (its own `lat = sani._C_LAT = 0.5`,
    isochoric).  Run against the first P1 binary, the clamp NEVER fired
    (`implexDetail` count == 0 for the whole 40-step leg) even though the
    BASE material's own low-p branch is independently measured to go
    NEGATIVE on that exact deck (`test_ladruno_sanisand.py`: p = -0.5647 kPa
    at deviatoric step 1 with `-Pmin` at the class default). The likely
    reason: `mImplexDEpsP` is exactly zero until the FIRST commit after the
    stage flip (`ladrunoImplexInitState`'s own comment -- "which is what
    makes the first plastic step a pure elastic prediction"), so `sigma~` on
    an isochoric shear's early steps is a much MILDER elastic-only
    prediction than the base's full elastoplastic return, and apparently
    mild enough here to stay clear of the floor. P0's own G3 path was not
    isochoric shear either -- it was described as "volumetric extension...
    while the point keeps flowing" (`_adr92_p0_oracle_results.md` section 5).
    `sani._c_series_at` already exposes the lat ratio as a parameter, so this
    supplies `lat > 0.5` (net dilation) rather than writing a new series
    generator; e_conf/e_ax/n_conf stay the SAME already-proven magnitudes
    the sibling `_PMIN_*` family uses (only the ratio changes), so this is
    not a numerically untested regime.
    """
    if e_conf is None:
        e_conf = sani._PMIN_E_CONF_LOW      # -> p ~ 0.1145 kPa at the flip
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    s_lat, s_ax = sani._c_series_at(e_conf, n_dev, n_conf=sani._C_N_CONF,
                                    e_ax=sani._C_E_AX, lat=lat)
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -e_conf)
            if y == 1.:
                ops.sp(n, 2, -e_conf)
    ops.pattern('Plain', 2, 2)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -e_conf)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-13, 25, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    return n_dev


def test_floor_clamp_fires_and_implexDetail_reports_it():
    """On a NET-DILATING path at the same low confinement the base material's
    own low-p branch is independently known to threaten, the NEW clamp on
    `sigma~` (W2) fires too, and `implexDetail` says so -- checked PER STEP,
    not only at the end, both for the clamp's own fire count/flag and for the
    `||d_sigma||^2 == ||d_dev||^2 + 3*(dp)^2` identity `ladrunoImplexMeasureError`
    computes (dev and I1 are orthogonal, so this must hold at every read,
    independent of whether the clamp fired on that particular step).

    SEE `_build_floor_seeking_deck` for why this is NOT the first draft's
    isochoric deck: that one measured clamp count == 0 against the first P1
    binary, and the fix is a path shape (net dilation, not constant-volume
    shear), not a weaker assertion.  The count is still not asserted as an
    exact number -- only that it is possible on a build to which this file
    has never been run -- see `_adr92_p0_oracle_results.md` section 5 for the
    same MECHANISM (not the same deck) measured directly by P0.
    """
    tag = 8107
    opts = ('-Presidual', 0.0, '-Pmin', sani._PMIN_LADRUNO, '-honorTolR', 0,
            '-implex', '-maxSubsteps', _CAP_ADEQUATE)
    n_dev = _build_floor_seeking_deck(tag, opts)

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(sani._C_N_CONF):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    max_count = 0.0
    ever_fired = False
    for step in range(n_dev):
        assert ops.analyze(1) == 0, f'deviatoric step {step + 1} failed'

        detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
        assert len(detail) == 6, ('implexDetail did not return the documented '
                                  '6-component vector', step, detail)
        total, dev, vol, fired_last, count, f_last = detail

        # The identity, ADR-92 P1 section 4.1: ||d_sigma||^2 == ||d_dev||^2 +
        # 3*(dp)^2, i.e. total^2 == dev^2 + vol^2 over the SAME denominator --
        # must hold at EVERY step, clamp firing or not.
        lhs = total * total
        rhs = dev * dev + vol * vol
        scale = max(lhs, rhs, 1.0e-30)
        assert abs(lhs - rhs) <= 1.0e-8 * scale, (
            'implexDetail\'s total/deviatoric/volumetric legs violate the '
            '||d_sigma||^2 == ||d_dev||^2 + 3*(dp)^2 identity the source '
            'claims (ladrunoImplexMeasureError)', step, detail)
        assert fired_last in (0.0, 1.0), ('implexDetail\'s fired flag was '
                                          'not boolean-valued', step, detail)

        max_count = max(max_count, count)
        ever_fired = ever_fired or bool(fired_last)

    assert max_count > 0 and ever_fired, (
        'the p_min clamp on sigma~ never fired on a NET-DILATING path at the '
        'same low confinement (e_conf = sani._PMIN_E_CONF_LOW) where the '
        'BASE material\'s own low-p branch is independently measured to go '
        'negative (test_ladruno_sanisand.py). If the clamp genuinely never '
        'engages on this shape either, either W2 is not wired into '
        'ladrunoImplexTrial(), or the path still is not aggressive enough -- '
        'raise `lat` further rather than weakening this assertion',
        max_count, ever_fired)


# ===========================================================================
#  Gate 3 -- tangent identity, on a NON-pseudo dt source
# ===========================================================================
#
#  MECHANISM.  `ManzariDafalias::commitState()` always commits the IMPLICIT
#  return (module docstring above), so a zero-free-DOF deck can never show
#  sigma~ to Python -- it is overwritten by the commit within the SAME
#  `analyze(1)` call. The only way to catch sigma~ is a genuinely free-DOF
#  deck with a Newton iteration that FAILS to converge (so the domain never
#  commits), forced by pairing `algorithm('Newton')` with an UNREACHABLE
#  convergence tolerance and `maxIter = 1`: exactly one
#  formTangent/formResidual/solve/update cycle runs (the standard OpenSees
#  Newton loop: form, solve, update U -> pushes the new trial strain into the
#  material -> test; with maxIter = 1 the loop stops after that ONE update
#  regardless of what test() says), then `test()` fails against the
#  unreachable tolerance and `analyze()` returns non-zero WITHOUT commit
#  ("Domain::update() returns the element's code ... IncrementalIntegrator::
#  update turns it into a failed step -- none of which involves solving
#  anything", `test_ladruno_sanisand_integrator.py`).
#
#  DO NOT call `ops.reset()` here.  MEASURED 2026-09-06: `ops.reset()` is
#  `Domain::revertToStart()` (`OpenSeesCommands.cpp:2539`), NOT
#  `revertToLastCommit()` -- it runs `ManzariDafalias::initialize()`, which
#  ZEROES `mSigma_n`/`mEpsilon_n`/`mAlpha_n`/`mFabric_n`, and then ends
#  `return this->update()` (`Domain.cpp:2385`), pushing one state determination
#  through the zeroed state.  With `sigma~ = 0 + Ce:0 = 0` the p_min clamp
#  fires HONESTLY and `getStress()` returns
#  `[-0.0101, -0.0101, -0.0101, -0.0, -0.0, -0.0]` -- the three IEEE negative
#  zeros are the fingerprint, since `-1.0 * (+0.0)` requires `dev(sigma~)` to be
#  EXACTLY zero, which no real load path produces.  Three tests read that as
#  "revertToLastCommit is not restoring".  It was the probe, not the material.
#
#  Nothing is needed in its place: a failed `StaticAnalysis` step ALREADY calls
#  `theDomain->revertToLastCommit()` and `theIntegrator->revertToLastStep()`
#  before `analyze()` returns non-zero (`StaticAnalysis.cpp:185` and four
#  sibling sites), so the correct revert has happened by the time we look.
#  There is no Python verb that calls `revertToLastCommit()` on its own; a
#  deliberately-failed step is the way to reach it.
#  commit) and the probe's OWN trial-strain norm.
#
#  WHY +h AND -h SHARE THE SAME f UNDER `-implexDt strain`.  `f` is frozen at
#  the FIRST trial call after a commit/revert from
#  `dt = GetNorm_Cov(mEpsilon - mEpsilon_n)` (ladrunoImplexArmStep). Each
#  probe below is driven by an isolated axial force pattern applied AFTER
#  `loadConst` has frozen everything else, so the probe's OWN strain
#  increment is (to within the Newton solve's own linearity, exact under a
#  frozen elastic tangent) a scalar multiple of one fixed direction; a norm is
#  insensitive to the SIGN of that scalar, so dt, and therefore f, is the same
#  for the +h and -h probes. With f AND Ce AND d_eps_p(n) equal for both,
#  sigma~ = sigma_n + Ce*(d_eps - f*d_eps_p) is EXACTLY affine in d_eps, so the
#  secant (sigma~(+h) - sigma~(-h)) / (eps~(+h) - eps~(-h)) reproduces Ce to
#  floating-point precision, not merely to FD truncation order -- this is why
#  the test can use a tight tolerance despite never seeing a "true" small h.
#
#  `-implexDt strain`, NOT `pseudo`.  Per the plan section 4.3 item 1: under
#  `pseudo`, `ops_Dt` is constant within a step regardless of what the trial
#  strain does, so the freeze is a no-op and this exact trap (an UNFROZEN,
#  strain-dependent f silently breaking d(sigma~)/d(eps) == Ce) would not be
#  exercised at all.  `strain` is the one source where the freeze fix does
#  real work, so it is the one this gate has to run on.

_PROBE_TAG = 8110
_PROBE_P0 = 100.0            # kPa -- away from the p_min floor; this gate is
                             # about the tangent, not the clamp (gate 4 owns that)
_PROBE_N_CONF = sanint._TX_N_CONF
_PROBE_DQ_NOMINAL = 6.0      # kPa (small fraction of _PROBE_P0), a few steps to
                             # establish a genuinely plastic committed history
_PROBE_N_HISTORY = 4
_PROBE_TOL_REL = sanint._TX_TOL_REL
_PROBE_MAXITER = sanint._TX_MAXITER
_PROBE_H_FRACTION = 1.0e-6   # the probe force as a fraction of _PROBE_DQ_NOMINAL/4


def _build_probe_triaxial(tag):
    """A LadrunoBrick drained-triaxial cube -- `sanint._build_triaxial`'s shape,
    with the ADR-92 P1 flags threaded through (that builder does not accept
    extra flags beyond the five positional optionals)."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,           # IntScheme, TanType, JacoType, TolF, TolR
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM,
                   '-implex', '-maxSubsteps', _CAP_ADEQUATE, '-implexDt', 'strain')
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
               '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    q = _PROBE_P0 / 4.0
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.load(n, -q, 0.0, 0.0)
            if y == 1.:
                ops.load(n, 0.0, -q, 0.0)
            if k == 1:
                ops.load(n, 0.0, 0.0, -q)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormUnbalance', _PROBE_TOL_REL * _PROBE_P0, _PROBE_MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _PROBE_N_CONF)
    ops.analysis('Static')


def _establish_plastic_history(tag):
    """Confine, flip, and take a few NORMAL (converged, committed) deviatoric
    steps so mImplexDEpsP(n) != 0 by the time the probe runs."""
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(_PROBE_N_CONF):
        assert ops.analyze(1) == 0, f'triaxial confinement step {step + 1} failed'
    ops.loadConst('-time', 0.0)
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    dq = _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 1.0 / _PROBE_N_HISTORY)
    for step in range(_PROBE_N_HISTORY):
        assert ops.analyze(1) == 0, f'history-building step {step + 1} failed'
    ops.loadConst('-time', 0.0)


def _probe_once(tag, sign, probe_pattern_tag, probe_ts_tag):
    """One forced-single-iteration trial: apply a TINY axial force (sign *
    h_fraction * the nominal per-step load), force exactly one failed Newton
    iteration, read back (stress, strain) BEFORE any commit, then revert.

    Returns (stress[6], strain[6]) -- element-convention (tension-positive,
    the -1.0 flip already applied by LadrunoSANISAND3D::getStress/getStrain).
    """
    h = sign * _PROBE_H_FRACTION * (_PROBE_DQ_NOMINAL / 4.0)
    ops.timeSeries('Linear', probe_ts_tag)
    ops.pattern('Plain', probe_pattern_tag, probe_ts_tag)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -h)

    ops.test('NormUnbalance', 1.0e-300, 1, 0)   # unreachable -> exactly 1 iteration, no commit
    ops.integrator('LoadControl', 1.0)
    rc = ops.analyze(1)
    assert rc != 0, (
        'the probe step CONVERGED and therefore COMMITTED -- it was supposed '
        'to fail an unreachable NormUnbalance tolerance after exactly one '
        'Newton iteration so sigma~ could be read before the commit '
        'overwrites it with the implicit return. If this converges, the '
        'harness itself is broken, not the material', rc)

    stress = list(ops.eleResponse(1, 'material', 1, 'stress'))
    strain = list(ops.eleResponse(1, 'material', 1, 'strain'))

    # NO ops.reset() -- see the block comment above; the failed step already
    # reverted, and reset() would zero the committed state and hand back the
    # clamp fingerprint instead.
    ops.remove('loadPattern', probe_pattern_tag)
    return stress, strain


def test_tangent_identity_frozen_ce_on_strain_dt_source():
    """Returned tangent (`getTangent()`, frozen to `Ce(p_n)` under `-implex`)
    reproduces a numerical `d(sigma~)/d(eps)` taken from TWO forced,
    uncommitted trial evaluations on a `-implexDt strain` deck -- see the
    long mechanism note above this test for why this is the honest way to
    reach `sigma~` from Python and why the two probes share the same frozen
    `f`.
    """
    _build_probe_triaxial(_PROBE_TAG)
    _establish_plastic_history(_PROBE_TAG)

    ce_flat = list(ops.eleResponse(1, 'material', 1, 'tangent'))
    stress_before = list(ops.eleResponse(1, 'material', 1, 'stress'))

    sp, ep = _probe_once(_PROBE_TAG, +1.0, probe_pattern_tag=90, probe_ts_tag=90)
    sm, em = _probe_once(_PROBE_TAG, -1.0, probe_pattern_tag=91, probe_ts_tag=91)

    stress_after = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert stress_before == stress_after, (
        'the committed stress moved across the two probe-and-reset cycles -- '
        'revertToLastCommit() is not fully restoring the committed state, so '
        'the two probes are not actually taken from the same base state',
        stress_before, stress_after)

    d_sigma = [a - b for a, b in zip(sp, sm)]
    d_eps = [a - b for a, b in zip(ep, em)]
    assert _vnorm(d_eps) > 0.0, (
        'the two probes produced the SAME trial strain -- the +h/-h force '
        'perturbation did not move the free axial DOFs at all', ep, em)

    predicted = _matvec6(ce_flat, d_eps)
    err = _vnorm([a - b for a, b in zip(d_sigma, predicted)])
    scale = max(_vnorm(d_sigma), 1.0e-30)
    rel = err / scale
    assert rel < 1.0e-6, (
        'the frozen tangent Ce(p_n) does not reproduce the numerical '
        'd(sigma~)/d(eps) from the two forced-uncommitted probes. Per the '
        'plan section 4.3 item 1, this specific gate only exercises the '
        '"-implexDt strain" f-freeze fix; a failure here on a fresh build is '
        'exactly the trap that fix exists for', rel, d_sigma, predicted)


# ===========================================================================
#  Parser / runtime refusals
# ===========================================================================

@pytest.mark.parametrize('scheme', [3, 5, 7, 8, 9])
def test_implex_refuses_unsupported_schemes(scheme):
    """ADR 92 D3 (as reversed by P0): `-implex` refuses IntScheme
    3/5/7/8/9 -- no error control on the substep, so the companion could
    never report a failed return. `setLadrunoImplexOptions` raises this via
    `opserr` + returns -1, which the parser turns into a hard construction
    failure (`delete theMaterial; return 0;`), which openseespy surfaces as
    an exception -- the same idiom `test_ladruno_sanisand.py` uses for
    `-honorTolR`'s out-of-range refusal.
    """
    ops.wipe()
    tag = 8200 + scheme
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                       scheme, 2, 1, 1.0e-7, 1.0e-7,
                       '-implex', '-maxSubsteps', _CAP_TIGHT)


def test_implex_scheme1_requires_maxsubsteps():
    """`-implex` on the default companion, IntScheme 1 (ModifiedEuler),
    REQUIRES `-maxSubsteps > 0` -- without it the companion at commitState has
    no way to fail rather than force-accept at `dT_min` (ADR-90 GATE U: single
    updates measured at 34 minutes)."""
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 8210, *_PARAMS, '-implex')


def test_implex_maxsubsteps_zero_is_also_refused():
    """`-maxSubsteps 0` means UNCAPPED (vanilla's own convention,
    `test_uncapped_is_byte_identical`), so it must be refused exactly like
    omitting the flag -- a deck that says `-maxSubsteps 0` explicitly must not
    be treated as having satisfied the requirement."""
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 8211, *_PARAMS,
                       '-implex', '-maxSubsteps', 0)


# ---------------------------------------------------------------------------
#  A refusal is only OBSERVABLE from analyze()'s return code on a deck with
#  genuine free DOFs.
#
#  The first draft of both refusal tests below used the confine-first
#  ZERO-free-DOF rig (`sanint._build_brick`, matching the SUBSTEP-CAP tests'
#  own idiom).  Run against the first P1 binary: the material DID refuse --
#  the D2 warning ("the pseudo-time increment is negative (-1)") printed 16
#  times -- but `analyze()` still returned 0, because with zero free DOFs
#  there is no equation for a refusal to fail, and the domain's pseudo time
#  moved on regardless.  Confirmed on a genuinely free-DOF deck (both
#  `stdBrick` and `LadrunoBrick`): the SAME refusal returns rc = -3.  So
#  every test in this file that needs to OBSERVE a refusal through
#  `analyze()`'s return code now uses `_build_free_dof_triaxial` below --
#  the same single-`LadrunoBrick` drained-triaxial shape gate 3's tangent-
#  identity probe uses, minus the `-implexDt strain` source (that source
#  cannot go negative at all -- it is a norm -- and D2's refusal is
#  specifically about the DEFAULT `pseudo` source going negative).
# ---------------------------------------------------------------------------

def _build_free_dof_triaxial(tag, extra_opts=(), p0=None):
    """A LadrunoBrick drained-triaxial cube with GENUINE free DOFs (the
    positive-face nodes are LOADED, not `sp`-prescribed) -- see the note
    above this function for why the refusal tests need this and the
    zero-free-DOF rig cannot substitute for it."""
    if p0 is None:
        p0 = _PROBE_P0
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,           # IntScheme, TanType, JacoType, TolF, TolR
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM,
                   *extra_opts)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
               '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    q = p0 / 4.0
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.load(n, -q, 0.0, 0.0)
            if y == 1.:
                ops.load(n, 0.0, -q, 0.0)
            if k == 1:
                ops.load(n, 0.0, 0.0, -q)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormUnbalance', _PROBE_TOL_REL * p0, _PROBE_MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _PROBE_N_CONF)
    ops.analysis('Static')


def test_implex_negative_pseudo_dt_is_refused_and_leaves_committed_state_unchanged():
    """D2's refusal by symptom: a NEGATIVE pseudo-time increment (`ops_Dt < 0`)
    is refused via `LADRUNO_MATERIAL_REFUSED`, and the committed state is
    unchanged across the refusal.

    HOW A NEGATIVE `ops_Dt` IS PRODUCED, without touching DisplacementControl
    or arc-length at all: `ops_Dt` is `currentTime - committedTime`
    (`ladrunoImplexArmStep`'s own comment), and `LoadControl`'s
    `deltaLambda` argument becomes exactly that increment. `LoadControl(-x)`
    for one step is therefore the simplest possible way to make `ops_Dt`
    negative under the DEFAULT (`pseudo`) dt source -- no limit point, no
    special integrator needed.

    A FREE-DOF deck, not the confine-first `sp` rig -- see the note above
    `_build_free_dof_triaxial` for the measured reason (a refusal is
    invisible to `analyze()`'s return code with zero free DOFs; the deck
    still visibly PRINTS the refusal warning, it just cannot be asserted on
    via the return code).
    """
    tag = 8212
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    _build_free_dof_triaxial(tag, opts)
    _establish_plastic_history(tag)

    before = list(ops.eleResponse(1, 'material', 1, 'stress'))

    ops.integrator('LoadControl', -1.0)
    rc = ops.analyze(1)
    assert rc != 0, (
        'a negative pseudo-time increment (LoadControl(-1.0)) was NOT '
        'refused. ADR-92 D2 requires the material to refuse this step '
        'because dt_{n+1}/dt_n is the extrapolation factor itself and a '
        'reversed load factor makes it a wrong answer that would pass every '
        'other gate', rc)

    # NO ops.reset(): the failed step's own revertToLastCommit already ran.
    after = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert before == after, (
        'the committed stress moved across a REFUSED step -- the refusal is '
        'supposed to leave the committed state untouched (revertToLastCommit '
        'restores the trial state on top of it)', before, after)


def test_implexcontrol_refuses_past_tolerance_and_leaves_committed_state_unchanged():
    """`-implexControl` refuses a step whose extrapolation error exceeds its
    tolerance, via `LADRUNO_MATERIAL_REFUSED`, and the committed state is
    unchanged across the refusal.

    A FREE-DOF deck (`_build_free_dof_triaxial`), for the same measured
    reason as the D2 test above: a zero-free-DOF deck's `analyze()` cannot
    fail on a refusal at all (no equations to fail), so it would report
    success and the committed state would move regardless of what the
    material returned -- which is exactly the failure mode the first draft
    of this test hit.

    THE DECK'S EXPECTATION IS QUALITATIVE, sourced from the P0 oracle's own
    measurement at the SAME calibrated parameters (`_adr92_p0_oracle_results.md`
    section 4, T2 table, p0 = 5 kPa): `implexError` grows from ~2.7e-3 at the
    nominal increment to 1.26 at 5x the nominal increment -- both values well
    on either side of the ADR-92 default tolerance (0.05). This test does not
    try to reproduce those exact numbers (a different deck, a different
    control variable -- FORCE, not strain) -- it confines to p0 = 5 kPa (the
    exact oracle confinement), takes `_PROBE_N_HISTORY` NOMINAL deviatoric
    steps (expected to pass), then ONE step TEN TIMES that per-step load
    (expected to fail the same way, on the same physics, at the same low
    confinement), and asserts only the qualitative outcome the oracle
    licenses: the big step is refused.
    """
    tag = 8213
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', 0.05, 0.01)
    _build_free_dof_triaxial(tag, opts, p0=5.0)
    _establish_plastic_history(tag)

    before = list(ops.eleResponse(1, 'material', 1, 'stress'))

    # The 10x jump: a THIRD pattern, on top of the confinement (pattern 1)
    # and the nominal deviatoric history (pattern 2), both already frozen by
    # _establish_plastic_history's own loadConst.
    big_dq = 10.0 * _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 3)
    ops.pattern('Plain', 3, 3)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -big_dq)
    ops.integrator('LoadControl', 1.0)
    rc = ops.analyze(1)
    assert rc != 0, (
        'a deviatoric load increment 10x the nominal one, at p0 = 5 kPa, was '
        'NOT refused by -implexControl. P0 measured implexError escalating '
        'sharply with increment size at this exact confinement (section 4 '
        'of the oracle results memo); if this genuinely never refuses on '
        'this deck, re-derive the deck rather than weakening this assertion',
        rc)

    # NO ops.reset(): the failed step's own revertToLastCommit already ran.
    after = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert before == after, (
        'the committed stress moved across an -implexControl refusal', before, after)


# ===========================================================================
#  getCopy / sendSelf / recvSelf -- the flags AND the committed d_eps_p
# ===========================================================================

def _roundtrip_implex(matcmd, tag, opts, n_cut, n_total):
    """`test_ladruno_sanisand.py`'s own `_roundtrip`, with the plastic-leg
    split parametrised and IMPL-EX opts threaded through. Returns (stress at
    the save, stress right after the restore, final stress)."""
    with tempfile.TemporaryDirectory(prefix='ladruno_sanisand_implex_',
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, 'sanisand_implex_rt')

        sani._build(matcmd, tag, opts)
        sani._elastic_leg(tag)
        ops.updateMaterialStage('-material', tag, '-stage', 1)
        for step in range(n_cut):
            assert ops.analyze(1) == 0, f'pre-save plastic step {step + 1} failed'
        mid = sani._stress()

        try:
            ops.database('File', dbpath)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f'database() unsupported in this build: {exc}')
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip('database save returned failure on this build')

        sani._build(matcmd, tag, opts)                  # fresh, uncommitted skeleton
        ops.database('File', dbpath)
        ops.restore(1)
        after = sani._stress()

        ops.updateMaterialStage('-material', tag, '-stage', 1)
        ops.wipeAnalysis()
        sani._analysis()
        for step in range(n_total - n_cut):
            assert ops.analyze(1) == 0, f'post-restore plastic step {step + 1} failed'
        out = sani._stress()
        ops.wipe()
        return mid, after, out


def test_implex_db_roundtrip_carries_flags_and_history():
    """A restored `-implex` material must keep running the SAME extrapolation
    it was saving -- both the option flags (the six-override discipline, ADR
    86 sec.4.4/ADR-92 P1 section header, "the wire grows") and the ONE new
    committed history variable, `mImplexDEpsP` (wire slots 16..21, sent
    because "an MP worker that receives a material mid-analysis and restarts
    its extrapolation from zero silently runs a different constitutive law
    from the rank beside it" -- the ADR 86 section 3 defect, in a new
    variable).

    THE ASSERTION IS BEHAVIOURAL, on the `test_db_roundtrip_carries_presidual`
    template: not "restore returned 0", but "the restored material finishes
    the remaining plastic steps on the SAME committed answer an unbroken run
    would reach". A material that restored WITHOUT its `-implex` flags (the
    plain base default) or WITHOUT its `d_eps_p` history (silently
    reinitialised to zero, matching a virgin Gauss point's first plastic
    step) would very likely diverge from the reference somewhere in the
    remaining plastic steps, because the extrapolation factor and the history
    it carries feed directly into every subsequent committed stress.
    """
    opts_on = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)

    ref = sani._drive('LadrunoSANISAND', 8300, opts_on)

    n_cut = 12
    mid, after, out = _roundtrip_implex('LadrunoSANISAND', 8301, opts_on,
                                        n_cut=n_cut, n_total=sani._N_PL)
    assert sani._reldiff(mid, after) <= _EQ_TOL, (
        'the restored -implex material did not come back on the state it '
        'was saved at', mid, after)
    assert sani._reldiff(ref, out) <= _EQ_TOL, (
        'after the round trip the -implex material no longer finishes on '
        'the reference (unbroken) answer -- the flags or the d_eps_p history '
        'did not survive sendSelf/recvSelf', ref, out)

    # Non-vacuity: an -implex OFF reference run on the SAME deck up to the
    # cut, if the round trip silently dropped -implex entirely (falling back
    # to the base default), would land here instead. It must not.
    ref_off = sani._drive('LadrunoSANISAND', 8302, sani._OPTS_VANILLA)
    gap = sani._reldiff(ref_off, ref)
    if gap < _SENSITIVITY_FLOOR:
        pytest.skip(
            'this deck does not distinguish -implex on from off '
            f'(gap {gap:.3e} < floor {_SENSITIVITY_FLOOR:.0e}), so the '
            'restored-vs-off comparison below would prove nothing; the '
            'restored-vs-reference comparison above already carries this '
            'test')
    assert sani._reldiff(ref_off, out) >= _SENSITIVITY_FLOOR, (
        'after the round trip the material is running as though -implex '
        'were off', ref_off, out)


def test_implex_getcopy_does_not_share_history_across_gauss_points():
    """`getCopy(const char*)` (the prototype -> per-Gauss-point path every
    `stdBrick`/`LadrunoBrick` construction uses) deliberately does NOT
    transfer `mImplexDEpsP` -- "a fresh integration point starts with
    d_eps_p = 0" (`LadrunoSANISAND.cpp`'s own comment on the getCopy(const
    char*) implementation). Two elements built from the SAME prototype
    material tag must therefore be INDEPENDENT: driving one through a full
    plastic history must not perturb what the other, never touched, reports.

    Two stdBrick cubes sharing ONE `LadrunoSANISAND` tag -- exactly what
    `getCopy("ThreeDimensional")` is for (each element calls it once per
    Gauss point on construction, independently). Element 1's nodes are driven
    through the confine-first path; element 2's nodes are never given an
    `sp` beyond the base fixities, so it sits at zero strain throughout.
    """
    tag = 8303
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    # element 1: nodes 1-8 (driven)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    # element 2: nodes 11-18 (never driven)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(10 + 4 * k + j + 1, x, y, float(k))

    opts = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    ops.element('stdBrick', 2, 11, 12, 13, 14, 15, 16, 17, 18, tag)

    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    # element 2: ALL three DOFs of ALL eight nodes fixed -- zero free DOFs and
    # zero prescribed strain, an inert bystander rather than merely "not
    # pushed as hard".
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(10 + 4 * k + j + 1, 1, 1, 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, sani._LAT * sani._E_AX)
            if y == 1.:
                ops.sp(n, 2, sani._LAT * sani._E_AX)
            if k == 1:
                ops.sp(n, 3, -sani._E_AX)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-13, 25, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / sani._NTOT)
    ops.analysis('Static')

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(sani._N_EL):
        assert ops.analyze(1) == 0, f'elastic-stage step {step + 1} failed'
    ops.updateMaterialStage('-material', tag, '-stage', 1)
    for step in range(sani._N_PL):
        assert ops.analyze(1) == 0, f'plastic-stage step {step + 1} failed'

    driven = list(ops.eleResponse(1, 'material', 1, 'stress'))
    bystander_strain = list(ops.eleResponse(2, 'material', 1, 'strain'))
    bystander_pstrain = list(ops.eleResponse(2, 'material', 1, 'plasticstrains'))

    assert _vnorm(driven) > 0.0, (
        'element 1 was never actually driven -- the deck is broken, not the '
        'claim under test', driven)
    assert all(abs(v) == 0.0 for v in bystander_strain), (
        'element 2 (never sp-driven) reports nonzero strain -- the deck '
        'setup, not getCopy, is at fault', bystander_strain)
    assert all(abs(v) == 0.0 for v in bystander_pstrain), (
        'element 2, sitting at zero strain the whole run, reports NONZERO '
        'plastic strain -- getCopy("ThreeDimensional") is leaking element '
        '1\'s d_eps_p history into element 2\'s Gauss-point copy, which '
        'shares only the PROTOTYPE material tag and must otherwise be fully '
        'independent', bystander_pstrain)
