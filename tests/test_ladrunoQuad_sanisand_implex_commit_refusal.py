"""WP-99 (F7) -- the COMMIT-TIME companion refusal must stop the run.

`Domain::commit()` (SRC/domain/domain/Domain.cpp) is `elePtr->commitState();`
-- bare, the return value dropped -- so no element, fork or vanilla, can refuse
a step at commit time.  Before this work package `LadrunoSANISAND`'s
commit-time companion failure (`ladrunoImplexCommit()` finding
`mSubstepCapHitInME`) therefore committed the PARTIALLY integrated state
ModifiedEuler had left at T < 1, returned `LADRUNO_MATERIAL_REFUSED` into a
caller that drops it, and the analysis walked on reporting every step
converged.

WP-99 closes that: a failed commit-time companion commits NOTHING (the trial is
restored from the committed state, `ManzariDafalias::commitState()` is skipped)
and sets a STICKY per-instance latch.  The latch makes `ladrunoTrialUpdate()`
-- the one entry both wrappers' `setTrialStrain()` uses -- return the sentinel
for every later update, which every element that FORWARDS `setTrialStrain`'s
return code turns into a failed step, which the analysis turns into a stop.

WHY `LadrunoQuad` IS THE LEAD DECK HERE.  The propagation audit that opens this
work package found the four `opserr` texts in `LadrunoSANISAND.cpp` /
`ManzariDafalias.cpp` were wrong about which elements act on a refusal: they
said "today LadrunoBrick" and listed `QuadUP` among the discarders.  Measured on
`9c2f964ea`, the forwarding set is LadrunoBrick (sentinel-filtered, per
ADR-33/34), LadrunoBrick20, LadrunoQuad / LadrunoCST / LadrunoLST,
BezierTet10 / BezierTri6 and vanilla FourNodeQuad / FourNodeQuadUP; the
discarders are Brick (= `stdBrick`), BbarBrick, BrickUP, SSPbrick, SSPquad and
LadrunoSolidShell.  `LadrunoQuad::update()` sums the material codes
(`ret += theMaterial[i]->setTrialStrain(eps)`) and so propagates ANY nonzero
value; `LadrunoBrick` filters for the sentinel specifically.  Both are covered
below, so the latch is proven element-agnostic.

DECK SHAPE.  Every deck here has GENUINE free DOFs (positive faces LOADED, not
`sp`-prescribed).  A zero-free-DOF deck cannot show a refusal at all: with no
free equations `analyze()` converges trivially and returns 0 whatever the
material said -- the measured trap recorded at length in
`test_ladruno_sanisand_implex.py`'s own module docstring.

`-maxSubsteps 2` is the forcing function: two ModifiedEuler substeps cannot
integrate a genuine plastic increment, so the commit-time companion caps on the
first plastic commit.  With `-implexControl` OFF the trial never integrates at
all (the companion probe is inside the control branch), so the cap can only
happen at COMMIT -- which is precisely the hole this file tests.

MEASURED WALL TIME: ~4 s for the whole file (three decks, a handful of steps
each).
"""
import math

import pytest

from _testbed import ops

import test_ladruno_sanisand as sani

_PARAMS = sani._PARAMS
_P_ATM = _PARAMS[8]
_XY = sani._XY

# --- deck constants --------------------------------------------------------
_P0 = 100.0          # kPa confinement
_N_CONF = 5
_TOL_REL = 1.0e-3    # x P0, a force tolerance
_MAXITER = 60
_CAP_STARVED = 2     # two substeps cannot integrate a plastic increment
_DQ_BIG = 60.0       # kPa deviator in ONE step -- comfortably plastic


# ---------------------------------------------------------------------------
#  Decks
# ---------------------------------------------------------------------------
def _build_quad(tag, extra_opts=()):
    """A single plane-strain `LadrunoQuad -formulation bbar` with free DOFs.

    Rollers on the two negative edges, the two positive edges LOADED, so the
    four free DOFs carry a real residual that a refusal can fail.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    for j, (x, y) in enumerate(_XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,      # IntScheme TanType JacoType TolF TolR
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM,
                   *extra_opts)
    ops.element('LadrunoQuad', 1, 1, 2, 3, 4, tag,
                '-thick', 1.0, '-type', 'PlaneStrain', '-formulation', 'bbar')
    for j, (x, y) in enumerate(_XY):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)

    q = _P0 / 2.0
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for j, (x, y) in enumerate(_XY):
        n = j + 1
        if x == 1.:
            ops.load(n, -q, 0.0)
        if y == 1.:
            ops.load(n, 0.0, -q)
    _analysis()


def _build_brick(tag, extra_opts=()):
    """The 3D twin: one `LadrunoBrick -formulation bbar` drained-triaxial cube
    with free DOFs -- the same shape `test_ladruno_sanisand_implex.py`'s
    refusal tests use, kept here so this file's three cases are read side by
    side."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM,
                   *extra_opts)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
                '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    q = _P0 / 4.0
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
    _analysis()


def _analysis():
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormUnbalance', _TOL_REL * _P0, _MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _N_CONF)
    ops.analysis('Static')


def _confine_and_flip(tag):
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(_N_CONF):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.loadConst('-time', 0.0)
    ops.updateMaterialStage('-material', tag, '-stage', 1)


def _add_deviatoric_pattern_2d(ts_tag=2, pat_tag=2, dq=_DQ_BIG):
    ops.timeSeries('Linear', ts_tag)
    ops.pattern('Plain', pat_tag, ts_tag)
    for j, (x, y) in enumerate(_XY):
        if y == 1.:
            ops.load(j + 1, 0.0, -dq / 2.0)
    ops.integrator('LoadControl', 1.0)


def _add_deviatoric_pattern_3d(ts_tag=2, pat_tag=2, dq=_DQ_BIG):
    ops.timeSeries('Linear', ts_tag)
    ops.pattern('Plain', pat_tag, ts_tag)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq / 4.0)
    ops.integrator('LoadControl', 1.0)


def _refusals():
    r = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    assert len(r) == 5, (
        'implexRefusals must be the 5-component vector WP-99 (F7) documents '
        '(total, signChange, control, companion, commitLatched)', r)
    return r


def _stress():
    return list(ops.eleResponse(1, 'material', 1, 'stress'))


# ===========================================================================
#  (1) LadrunoQuad -- the commit-time hole, closed
# ===========================================================================
def test_quad_commit_time_companion_refusal_latches_and_stops_the_run():
    """The whole point of WP-99, on the element whose propagation the audit
    corrected.

    Sequence: confine (stage 0), flip, then ONE big deviatoric step whose
    commit-time companion cannot integrate the increment in `-maxSubsteps 2`.

      * that step's commit is refused -- the companion bucket of
        `implexRefusals` increments and the per-instance `commitLatched` slot
        goes to 1;
      * NOTHING is committed: the committed stress after the refused step is
        the stress from BEFORE it, bit for bit;
      * the NEXT `analyze(1)` returns nonzero, because
        `ladrunoTrialUpdate()` refuses while latched and `LadrunoQuad::update`
        forwards the code (`ret += theMaterial[i]->setTrialStrain(eps)`).

    Before WP-99 the first bullet's counter moved (ADR-92's own
    `test_companion_refusal_at_commit_is_observable` asserted exactly that) but
    the partial state WAS committed and the next step ran on happily -- which
    is the defect, not the contract.

    Kills: a mutant that drops the latch (next step converges); one that
    latches but still calls `ManzariDafalias::commitState()` (committed stress
    moves); one that clears the latch on `revertToLastCommit()` (the failed
    step's own revert would un-latch and the run would continue).
    """
    tag = 9910
    _build_quad(tag, ('-implex', '-maxSubsteps', _CAP_STARVED))
    _confine_and_flip(tag)

    before_stress = _stress()
    before_ref = _refusals()
    assert before_ref[4] == 0.0, ('latched before the deviatoric step even '
                                  'started', before_ref)

    _add_deviatoric_pattern_2d()
    rc_big = ops.analyze(1)

    after_ref = _refusals()
    assert after_ref[3] - before_ref[3] > 0, (
        'the commit-time companion cap never incremented the companion slot '
        'of implexRefusals -- either -maxSubsteps 2 was not starved enough on '
        'this deck (raise _DQ_BIG) or the refusal is not being counted',
        before_ref, after_ref)
    assert after_ref[4] == 1.0, (
        'the per-instance commit-refusal latch (implexRefusals[4]) is not set '
        'after a commit-time companion failure -- without it the run walks '
        'past an invalid commit, which is the whole defect WP-99 closes',
        after_ref)

    after_stress = _stress()
    assert after_stress == before_stress, (
        'the committed stress MOVED across a refused commit. WP-99 requires '
        'the partial state ModifiedEuler left at T < 1 to be discarded and '
        'the trial restored from the committed state -- a refusing material '
        'must not also be an inventing one', before_stress, after_stress)

    rc_next = ops.analyze(1)
    assert rc_next != 0, (
        'the step AFTER a refused commit still converged. The latch is '
        'supposed to make every later setTrialStrain return '
        'LADRUNO_MATERIAL_REFUSED, and LadrunoQuad::update forwards it '
        '(ret += setTrialStrain), so the analysis must fail here. rc_big was '
        f'{rc_big}', rc_next)

    assert _stress() == before_stress, (
        'the committed stress moved on the REFUSED follow-up step too',
        before_stress, _stress())


# ===========================================================================
#  (2) -implexControl -- unchanged, and recoverable
# ===========================================================================
def test_implexcontrol_refuses_the_same_cap_at_trial_and_does_not_latch():
    """The same starved cap, with `-implexControl` ON, is caught one phase
    EARLIER -- at the trial, inside `ladrunoImplexTrial`'s companion probe
    (`LadrunoSANISAND.cpp`, the `if (mSubstepCapHitInME)` branch) -- so:

      * `analyze()` fails on THAT step, not the one after it;
      * the companion bucket still increments (the same failure, caught
        earlier, deliberately shares the bucket);
      * the commit-time latch is NOT set, because no commit happened.

    That difference is the whole argument for `-implexControl` being the
    RECOVERABLE path: a driver can halve its increment and retry, whereas the
    commit-time latch is sticky by design (the analysis already accepted the
    step whose commit failed, so there is nothing to revert to).

    WP-99 must not change this path at all.
    """
    tag = 9911
    _build_quad(tag, ('-implex', '-maxSubsteps', _CAP_STARVED,
                      '-implexControl', 0.02, 0.01))
    _confine_and_flip(tag)

    before_stress = _stress()
    before_ref = _refusals()

    _add_deviatoric_pattern_2d()
    rc = ops.analyze(1)

    assert rc != 0, (
        '-implexControl did not refuse a step whose companion cannot '
        'integrate the increment at all (-maxSubsteps 2). That refusal is at '
        'the TRIAL and must fail the step it belongs to', rc)

    after_ref = _refusals()
    assert after_ref[3] - before_ref[3] > 0, (
        'the trial-time companion cap refusal was not counted in the '
        'companion bucket', before_ref, after_ref)
    assert after_ref[4] == 0.0, (
        'the commit-time latch was set by a TRIAL-time refusal. -implexControl '
        'refuses before anything commits, so there is nothing to latch -- '
        'setting it here would turn the recoverable path into the sticky one',
        after_ref)

    assert _stress() == before_stress, (
        'the committed stress moved across a trial-time -implexControl '
        'refusal', before_stress, _stress())


# ===========================================================================
#  (3) LadrunoBrick -- the latch is element-agnostic
# ===========================================================================
def test_brick_commit_time_companion_refusal_latches_and_stops_the_run():
    """Case (1) on `LadrunoBrick`, whose propagation rule is DIFFERENT:
    it acts only on the sentinel value `LADRUNO_MATERIAL_REFUSED`
    (`LadrunoBrick.cpp` :1034/:1080/:1111/:1182/:1827/:3319), not on any
    nonzero code, because ADR-33/34 requires ASDConcrete3D's negative
    "best-state" codes NOT to fail a step.

    The latch therefore has to return the SENTINEL, not merely something
    nonzero, or it would work on LadrunoQuad and be silently swallowed here --
    exactly the asymmetry the audit found in the shipped warning texts. This
    test is what makes the claim "element-agnostic" a measurement.
    """
    tag = 9912
    _build_brick(tag, ('-implex', '-maxSubsteps', _CAP_STARVED))
    _confine_and_flip(tag)

    before_stress = _stress()
    before_ref = _refusals()

    _add_deviatoric_pattern_3d()
    ops.analyze(1)

    after_ref = _refusals()
    assert after_ref[3] - before_ref[3] > 0, (
        'the commit-time companion cap never incremented the companion slot '
        'on the LadrunoBrick deck', before_ref, after_ref)
    assert after_ref[4] == 1.0, (
        'the commit-refusal latch is not set on LadrunoBrick', after_ref)
    assert _stress() == before_stress, (
        'the committed stress moved across a refused commit on LadrunoBrick',
        before_stress, _stress())

    rc_next = ops.analyze(1)
    assert rc_next != 0, (
        'the step after a refused commit converged on LadrunoBrick. The latch '
        'must return LADRUNO_MATERIAL_REFUSED itself -- LadrunoBrick filters '
        'for that exact sentinel and ignores every other nonzero code',
        rc_next)


# ===========================================================================
#  (4) the latch survives what it must, and clears only on revertToStart
# ===========================================================================
def test_latch_is_cleared_only_by_reverttostart():
    """`revertToLastCommit()` must NOT clear the latch and `revertToStart()`
    must.

    The argument, recorded in the header member note: the analysis has already
    ACCEPTED the step whose commit failed, so "the last commit" is the corrupt
    datum itself -- there is nothing to go back to, and a driver that halves
    its increment would restart from a state the material has just said it
    could not produce.  `revertToStart()` puts the material back at step 0,
    where that state no longer exists, so it is the one place the latch is
    released.

    `ops.reset()` is `Domain::revertToStart()`; the failed `analyze(1)` in the
    first half already ran `Domain::revertToLastCommit()` through
    `StaticAnalysis.cpp`'s failure path, so the "survives a revertToLastCommit"
    half is measured by the fact that the run is still refusing afterwards.
    """
    tag = 9913
    _build_quad(tag, ('-implex', '-maxSubsteps', _CAP_STARVED))
    _confine_and_flip(tag)

    _add_deviatoric_pattern_2d()
    ops.analyze(1)                      # latches at commit
    assert _refusals()[4] == 1.0, 'setup: the latch did not arm'

    assert ops.analyze(1) != 0, (
        'setup: the latched material stopped refusing before the revert test '
        'even began')
    assert _refusals()[4] == 1.0, (
        'the latch was cleared by the failed step\'s own '
        'revertToLastCommit (StaticAnalysis.cpp calls Domain::'
        'revertToLastCommit on a failed solveCurrentStep). It must be sticky '
        'there: the corrupt commit IS the state being reverted to')

    ops.reset()                         # Domain::revertToStart
    assert _refusals()[4] == 0.0, (
        'revertToStart did not clear the latch. It is the one place the '
        'material is put back at step 0, i.e. the one place the commit the '
        'latch protects against no longer exists')


# ===========================================================================
#  (5) the response contract
# ===========================================================================
def test_implexrefusals_carries_the_commitlatched_slot():
    """`implexRefusals` grew 4 -> 5 in WP-99. Slots 0-3 stay the process-wide
    counters ADR-92 shipped; slot 4 is the only PER-INSTANCE entry -- 1 while
    THIS integration point is refusing. A recorder or a driver needs that
    distinction: the counters answer "how many", slot 4 answers "is the run
    dead".
    """
    tag = 9914
    _build_quad(tag, ('-implex', '-maxSubsteps', 20000))
    _confine_and_flip(tag)
    r = _refusals()
    assert all(math.isfinite(v) for v in r), r
    assert r[4] == 0.0, (
        'a deck with an ADEQUATE cap reports itself latched -- slot 4 must be '
        'the latch, not a counter', r)
    assert r[0] == pytest.approx(r[1] + r[2] + r[3]), (
        'implexRefusals[0] is no longer the sum of the three refusal buckets '
        '-- slot 4 must NOT be folded into the total; it is a flag, not a '
        'count', r)
