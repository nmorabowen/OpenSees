"""ADR-86b — the two SANISAND *integrator* seams, and the element that carries them.

Companion to ``test_ladruno_sanisand.py``, which fingerprints the ADR-86 PR-1/2/3
CONSTANTS (``-Presidual`` / ``-Pmin`` / ``-honorTolR``).  This file covers what
ADR-90's GATE U measured and ADR-90 §8.1 handed back to ADR 86:

  T1  ``-maxSubsteps N`` -- a substep-COUNT cap on ``ManzariDafalias::ModifiedEuler``.
  T2  ``LadrunoSANISAND``'s parser default for ``TanType``, 0 (ELASTIC) -> 2 (consistent).

WHY THESE TWO DEFECTS SURVIVED THREE PRs OF ADVERSARIAL REVIEW: every SANISAND
deck this fork owns is a **zero-free-DOF material-point cube**.  On such a deck
there are no equations, so a wrong tangent costs exactly nothing, and the load
steps are small enough that the substepped return never seizes.  Both defects
are properties of a BOUNDARY-VALUE problem, and the first BVP the fork put this
material into (ADR-90 WP-A2, a strip footing on softening sand) hit both at once:

  * at the parser default the h0 = 0.25 leg advanced **11 steps to s/B = 1.6e-4
    in 350 s**; with ``TanType 2`` and a reachable convergence test the same deck
    reached s/B = 2.0e-3 in 36 s -- about **7x**;
  * single ``analyze(1)`` calls cost **11 to 34 minutes** (worst 2056 s = 59 % of
    one leg's entire budget in a single step), while the stepping controller used
    **0 of its 80** pinned subdivisions and ended 6400-25000x above its own step
    floor.  The integrator never FAILED -- vanilla force-accepts at ``dT_min``
    (``ManzariDafalias.cpp:1649-1663``) -- so nothing upstream could react.

THE FOUR GATES HERE, and what each one alone would catch:

  1. ``test_uncapped_is_byte_identical``  -- the off-switch.  Default 0 must be
     UNCAPPED, and a huge cap must be indistinguishable from no cap.  Three
     spellings of "no cap" must give a BIT-IDENTICAL committed stress; and the
     ``test_ladruno_sanisand`` confine-first deck must still reproduce its own
     RECORDED sensitivity number 9.325164e-02, which is the check that neither
     T1 nor T2 moved the material-point answer at all.
  2. ``test_cap_engages_and_the_step_fails`` -- the load-bearing gate.  Measures
     the substep count the uncapped deck actually spends (via the new
     ``substeps`` response), sets the cap BELOW it, and requires
     ``analyze(1) != 0``.  Self-calibrating on purpose: a hardcoded cap would
     silently stop engaging the day the integrator gets cheaper, and a gate that
     stops engaging still passes.
  3. ``test_subdivision_recovers_after_the_cap_fires`` -- the reason the cap
     returns FAILURE instead of force-accepting.  The same cap, the same total
     strain, taken in half-size increments, must succeed: that is the D16-style
     step-cut the GATE U campaign could never trigger.
  4. ``test_tantype_default_is_the_consistent_tangent`` +
     ``test_tantype_2_costs_fewer_newton_iterations`` -- T2, said and measured.
     The second one is on a genuinely FREE-DOF deck, because that is the only
     place the tangent can be observed at all.

WHAT IS DELIBERATELY *NOT* ASSERTED
  * Nothing about vanilla ``ManzariDafalias``'s own ``TanType`` default.  It is
    still 0 (``ManzariDafalias.cpp:93``, verified at source) and is deliberately
    NOT changed -- every existing vanilla deck and golden file depends on it.
    The two parsers now disagree on that one slot on purpose; the material
    ECHOES the tangent it will run so the divergence cannot be silent.
  * Nothing about ``-Presidual``.  ADR-86b writes the decision up
    (``86_ladruno_sanisand_handoff.md`` §5b) and changes no default.

ELEMENT CHOICE IS LOAD-BEARING, NOT INCIDENTAL.  ``stdBrick`` DISCARDS
``setTrialStrain``'s return code (the LEDGER_quirks entry "stdBrick swallows
material return codes"), so the cap would fire, the material would return -1, and
the analysis would report success -- the exact silence being fixed.  ADR-86b also
had to repair four of ``LadrunoBrick::update()``'s five paths, which discarded it
too.  Every test below that reads a return code uses ``LadrunoBrick``.

COST, measured on the dev box: the whole file is a few seconds.  Nothing here is
marked ``slow`` -- ``tests/conftest.py`` makes that tier opt-in and no workflow
passes ``--runslow``, so a slow marker means the test never runs in CI at all.
"""
import os

import pytest

from _testbed import ops
from _testbed.subprocess_run import run_python_script

# The recorded confine-first deck, its parameters, and its published numbers.
# Imported rather than copied on purpose: gate 1's whole claim is that ADR-86b
# did not move THAT deck's answer, which only holds if it is literally that deck.
import test_ladruno_sanisand as sani


pytestmark = [pytest.mark.zone_a]

_HERE = os.path.dirname(os.path.abspath(__file__))

_PARAMS = sani._PARAMS
_P_ATM = sani._P_ATM
_XY = sani._XY


# ---------------------------------------------------------------------------
#  Gate 1 -- the off-switch
# ---------------------------------------------------------------------------

# `test_vanilla_equivalence_is_not_vacuous`'s published measurement on the
# confine-first zero-free-DOF deck: reldiff(p_r = 1.01, p_r = 0).  Reproduced
# here with the cap unset, so that a change to either seam that perturbed the
# material-point answer would show up as a moved number rather than as a quiet
# drift nobody compares against.
_RECORDED_PR_SENSITIVITY = 9.325164e-02
_RECORDED_TOL = 1.0e-6


def test_uncapped_is_byte_identical():
    """`-maxSubsteps` unset, explicitly 0, and absurdly large: the SAME answer.

    Three spellings of "no cap".  All three must return a bit-identical
    committed stress -- `==` on the list, not a tolerance, because there is
    nothing here that may legitimately differ:

      * unset          -> the parser's default 0
      * `-maxSubsteps 0` -> the same value said out loud
      * `-maxSubsteps 100000000` -> a cap no real update can reach

    WHAT EACH SPELLING CATCHES.  Dropping the first would let the default move
    without anyone noticing.  Dropping the third would let the cap be implemented
    as "fire when set at all" -- which passes an unset-vs-0 comparison perfectly.
    The guard in `ModifiedEuler` is `mMaxSubstepsInME > 0 && taken > cap`, and
    only the third leg tests the second half of that conjunction.
    """
    off = sani._drive_confined('LadrunoSANISAND', 1, sani._OPTS_VANILLA)
    zero = sani._drive_confined('LadrunoSANISAND', 2,
                                sani._OPTS_VANILLA + ('-maxSubsteps', 0))
    huge = sani._drive_confined('LadrunoSANISAND', 3,
                                sani._OPTS_VANILLA + ('-maxSubsteps', 100000000))

    assert off == zero, (
        '-maxSubsteps 0 is not the same as omitting the flag -- the default is '
        'not "uncapped" and every existing deck just changed', off, zero)
    assert off == huge, (
        'a cap far above any reachable substep count changed the answer, so the '
        'cap is firing on the flag being SET rather than on the count being '
        'exceeded', off, huge)


def test_adr86b_did_not_move_the_recorded_material_point_answer():
    """The published 9.325164e-02 still reads 9.325164e-02.

    This is the regression gate for BOTH seams at once, and it is the reason it
    runs on `test_ladruno_sanisand`'s own deck rather than on a new one.

    `TanType` moved from 0 to 2 in this class's parser.  On a zero-free-DOF deck
    that MUST be inert -- there are no equations, so the tangent is never used to
    solve anything and the committed stress is a pure function of the prescribed
    strain path.  "Must be" is an argument; this is the measurement.  If the
    number moves, the tangent is reaching the answer somewhere it should not, and
    every A/B in the sibling file is measuring two things.
    """
    a = sani._drive_confined('LadrunoSANISAND', 11, sani._OPTS_VANILLA)
    b = sani._drive_confined('LadrunoSANISAND', 12, sani._OPTS_PR0)
    rel = sani._reldiff(a, b)
    assert abs(rel - _RECORDED_PR_SENSITIVITY) <= _RECORDED_TOL, (
        'the confine-first deck no longer reproduces its recorded p_residual '
        'sensitivity; ADR-86b was supposed to be inert on a zero-free-DOF deck',
        rel, _RECORDED_PR_SENSITIVITY)


# ---------------------------------------------------------------------------
#  A LadrunoBrick twin of the confine-first zero-free-DOF deck
#
#  Same cube, same rollers, same `sp`-prescribed positive faces, so there are
#  still ZERO free equations and no Newton tolerance can enter the answer.  The
#  ONE difference from `test_ladruno_sanisand`'s driver is the element:
#  `LadrunoBrick`, because `stdBrick` discards `setTrialStrain`'s return code and
#  a cap that fires inside one is invisible to `analyze()`.
#
#  Zero free DOFs does NOT make the return code unobservable.  `Domain::update()`
#  returns the element's code, `AnalysisModel::updateDomain` returns that, and
#  `IncrementalIntegrator::update` turns it into a failed step -- none of which
#  involves solving anything.  So `analyze(1)` reports the material failure on a
#  deck with no equations at all, which is exactly what makes this gate cheap.
#
#  The series are built INSIDE the builder from explicit arguments rather than
#  from module constants: `test_ladruno_sanisand._c_series` documents the trap
#  where monkeypatching a module constant desynchronises the loop count from the
#  Path it walks, and this file takes the same path lengths as a parameter for
#  gate 3, so it cannot afford that coupling.
# ---------------------------------------------------------------------------
_B_E_CONF = 3.0e-6      # isotropic compressive strain per direction, stage 0
_B_N_CONF = 10          # stage-0 confinement steps
_B_E_AX = 5.0e-3        # further axial compressive strain, stage 1
_B_LAT = 0.5            # lateral extension / axial compression -- 0.5 = isochoric
_B_N_DEV = 40           # stage-1 deviatoric steps (gate 3 halves the increment
                        # by doubling this)


def _brick_series(n_conf, n_dev):
    """The two ramp SHAPES on a unit pseudo-time grid, plus a HOLD point.

    `PathSeries::getFactor` returns 0 at and beyond its last time point, so a
    series without the trailing hold silently UNLOADS the model on the final
    step (measured in the sibling file: sigma ~ 0 and a nonsense eta of 63.6).
    """
    s_lat = [i / n_conf for i in range(n_conf + 1)]
    s_ax = list(s_lat)
    r_lat = _B_LAT * _B_E_AX / _B_E_CONF
    r_ax = _B_E_AX / _B_E_CONF
    for i in range(1, n_dev + 1):
        s_lat.append(1.0 - r_lat * i / n_dev)
        s_ax.append(1.0 + r_ax * i / n_dev)
    s_lat.append(s_lat[-1])
    s_ax.append(s_ax[-1])
    return s_lat, s_ax


def _build_brick(tag, opts=(), n_conf=_B_N_CONF, n_dev=_B_N_DEV):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
                '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    s_lat, s_ax = _brick_series(n_conf, n_dev)
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -_B_E_CONF)
            if y == 1.:
                ops.sp(n, 2, -_B_E_CONF)
    ops.pattern('Plain', 2, 2)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -_B_E_CONF)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-13, 25, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def _confine(tag, n_conf=_B_N_CONF):
    # EXPLICIT every time: mElastFlag is static and the constructor of any
    # Manzari-family material just reset it (ADR 86 risk 3).
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(n_conf):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.updateMaterialStage('-material', tag, '-stage', 1)


def _substeps_at_gp1():
    """[substeps taken, cap hit] for the LAST update at Gauss point 1."""
    return list(ops.eleResponse(1, 'material', 1, 'substeps'))


def _push(n_dev):
    """Take `n_dev` deviatoric steps; return the index of the first that failed
    (1-based) or 0 if all of them converged."""
    for step in range(n_dev):
        if ops.analyze(1) != 0:
            return step + 1
    return 0


# ---------------------------------------------------------------------------
#  Gate 2 -- the cap engages, and the STEP fails
# ---------------------------------------------------------------------------

# The child script for `test_cap_engages_and_the_step_fails`.  Kept as a string
# rather than a file so that the gate and the code it runs cannot drift apart.
#
# It re-imports THIS module inside the child purely to reuse `_build_brick` /
# `_confine` / `_substeps_at_gp1`; that is why the module is import-safe (no
# work at import time beyond constants).
_CAP_CHILD = r'''
import os
import sys
sys.path.insert(0, sys.argv[1])          # the tests/ directory
sys.path.insert(0, sys.argv[2])          # the directory holding opensees.pyd
import test_ladruno_sanisand_integrator as t
from _testbed import ops

# Ladruno (ADR-86b review-fix): print which binary the CHILD actually loaded,
# BEFORE doing any work.  `sys.path.insert(0, sys.argv[2])` only makes the
# right pyd importable first -- a stray tests/opensees.pyd or an installed
# .pth alias could still shadow it (see BUILD_GOTCHAS.md sec.4's
# "sys.path[0] wins" trap), and a child that quietly ran a STALE binary would
# make this whole gate measure yesterday's integrator instead of today's,
# passing for the wrong reason. The parent checks this marker matches
# argv[2] before trusting anything else in the output.
print("ENGINE=%s" % os.path.abspath(ops.__file__))

# leg 1 -- how expensive IS this update, uncapped?
t._build_brick(1, t.sani._OPTS_VANILLA)
t._confine(1)
rc0 = ops.analyze(1)
taken, hit = t._substeps_at_gp1()
print("UNCAPPED_RC=%d" % rc0)
print("SUBSTEPS_TAKEN=%d" % int(taken))
print("UNCAPPED_CAP_HIT=%d" % int(hit))

# leg 2 -- cap BELOW the measured cost
cap = int(taken) - 1
t._build_brick(2, t.sani._OPTS_VANILLA + ("-maxSubsteps", cap))
t._confine(2)
rc = ops.analyze(1)
print("CAP=%d" % cap)
print("CAPPED_RC=%d" % rc)
'''


def test_cap_engages_and_the_step_fails():
    """Measure the substep cost, cap below it, require a failed step.

    RUNS IN A CHILD PROCESS, and that is not incidental.  The material's cap
    warning is throttled on a PROCESS-WIDE budget of 10 (it has to be: every
    Gauss point is its own material object, so a per-instance latch is not a
    throttle at all).  `test_subdivision_recovers_after_the_cap_fires` below
    fires the cap too, and on this deck it fires it hundreds of times.  So an
    in-process version of this gate asserts `n_warn >= 1` against a budget some
    OTHER test may already have spent -- it would pass or fail on file order,
    which is the one thing a gate must never do.  A fresh interpreter gets a
    fresh budget by construction.  `run_python_script` (#785) also carries the
    two Windows traps that make `subprocess` reliable under pytest's fd capture.

    SELF-CALIBRATING, and that is the design.  A hardcoded cap ("3 should be too
    few") is a gate that quietly stops engaging the day the integrator gets
    cheaper or the deck changes -- and a gate that stops engaging still passes,
    which is the failure mode `86_ladruno_sanisand_handoff` §6 calls "a skip and
    a pass look identical".  So leg 1 runs UNCAPPED and asks the material how
    many substeps it actually spent, and leg 2 caps below that.

    THREE THINGS ARE ASSERTED, and the third is the one that matters:
      * the uncapped run reports a substep count > 1 and `cap_hit == 0`;
      * the capped run FAILS a step -- `analyze(1) != 0`.  That code has to
        survive material -> element -> Domain::update -> the integrator; before
        ADR-86b four of `LadrunoBrick::update()`'s five paths dropped it on the
        floor, and `stdBrick` still does;
      * the material SAID SO exactly once per process budget, in `opserr`.  A
        silent failure would be a step that fails for no stated reason, which is
        worse for the user than the seizure it replaces.

    THE COMMITTED STATE IS UNTOUCHED, by construction rather than by assertion:
    `integrate()` writes only the trial members (`mSigma`, `mAlpha`, `mFabric`),
    never `mSigma_n` and friends, so a failed update cannot corrupt what the last
    `commitState` stored.  Nothing here can distinguish that from a lucky
    coincidence, so it is stated and not claimed as measured.
    """
    rc_child, text = run_python_script(
        _CAP_CHILD, argv=(_HERE, os.path.dirname(ops.__file__)),
        merge_stderr=True, timeout=600)

    assert rc_child == 0, ('the child process itself failed -- this is not the '
                           'gate firing, it is the harness', text[-3000:])

    # Ladruno (ADR-86b review-fix): a child that quietly loaded a DIFFERENT
    # `opensees.pyd` than the one this test session was invoked with would
    # measure yesterday's integrator and pass for the wrong reason -- the
    # classic "stale binary, green run" false pass (same check as
    # `test_soe_zero_free_equations.py::_run_child`). Fail loud before
    # trusting anything else the child printed.
    engine_line = next((l for l in text.splitlines() if l.startswith('ENGINE=')),
                       None)
    assert engine_line is not None, (
        'no ENGINE= marker from the child -- it did not even get past the '
        'import, so nothing below can be trusted', text[-3000:])
    child_engine_dir = os.path.dirname(engine_line.split('=', 1)[1])
    expected_engine_dir = os.path.dirname(os.path.abspath(ops.__file__))
    assert os.path.normcase(child_engine_dir) == os.path.normcase(expected_engine_dir), (
        f'the child loaded opensees from {child_engine_dir!r}, but this test '
        f'session is running against {expected_engine_dir!r} -- the child is '
        f'measuring a STALE or DIFFERENT binary and every other assertion in '
        f'this test is meaningless until this is fixed', text[-3000:])

    def _marker(name):
        for line in text.splitlines():
            if line.startswith(name + '='):
                return int(line.split('=', 1)[1])
        raise AssertionError('marker %s missing from the child output; the child '
                             'did not reach it' % name + '\n' + text[-3000:])

    assert _marker('UNCAPPED_RC') == 0, (
        'the uncapped deviatoric step should converge', text[-3000:])
    taken = _marker('SUBSTEPS_TAKEN')
    assert _marker('UNCAPPED_CAP_HIT') == 0, (
        'the uncapped run reports the cap as HIT; the flag is not being reset '
        'per update', taken)
    assert taken >= 2, (
        'this deck spends fewer than 2 ModifiedEuler substeps per update, so '
        'there is no room to cap BELOW it -- pick a larger deviatoric '
        'increment before trusting anything else in this file', taken)

    cap = _marker('CAP')
    rc = _marker('CAPPED_RC')

    assert rc != 0, (
        'the substep cap fired but analyze() reported SUCCESS. The return code '
        'is being dropped somewhere between ManzariDafalias::ModifiedEuler and '
        'the integrator -- which is exactly the silence ADR-90 GATE U ran into '
        '(0 of 80 subdivisions used across six legs). Check that the element is '
        'a LadrunoBrick: stdBrick discards setTrialStrain return codes.',
        cap, taken, rc)

    # The construction echo is in the same capture, so an EMPTY capture means the
    # capture itself did not take -- the positive control without which the
    # absence of a warning below would prove nothing.
    assert 'LadrunoSANISAND tag 2' in text, (
        'nothing at all was captured from the child, so the absence of a warning '
        'below proves nothing', text[-3000:])
    n_warn = text.count('substep cap')
    assert n_warn >= 1, (
        'the cap fired and the step failed, but the material said nothing. A '
        'step that fails for no stated reason is worse for the user than the '
        'seizure it replaces.', text[-2000:])
    # Throttled to a PROCESS budget of 10, not per instance: every Gauss point is
    # its own material object.  A FRESH process is what makes this bound
    # meaningful -- in-process it would be measuring whatever budget the rest of
    # the file had already spent.
    assert n_warn <= 10, (
        'more than the 10-line process budget of cap warnings was emitted; the '
        'throttle is per instance again, which on a real mesh is ~400k lines',
        n_warn)

    # The ELEMENT must have said so too. The material's line alone does not prove
    # the refusal crossed the material/element boundary -- `rc != 0` above proves
    # it reached the integrator, and this proves the element named where.
    n_ele = text.count('REFUSED the trial strain')
    assert 1 <= n_ele <= 10, (
        'the element did not report the refusal exactly once per failing update '
        'within its own process budget', n_ele, text[-3000:])


# ---------------------------------------------------------------------------
#  Gate 3 -- and the subdivision then recovers
# ---------------------------------------------------------------------------

def test_subdivision_recovers_after_the_cap_fires():
    """The same cap, the same total strain, half-size increments: it converges.

    THIS IS THE WHOLE POINT OF FAILING RATHER THAN FORCE-ACCEPTING.  Vanilla, at
    `dT == dT_min`, force-accepts the substep and advances `T`
    (`ManzariDafalias.cpp:1649-1663`), so the update always "succeeds" and an
    adaptive driver has nothing to react to: ADR-90 GATE U's six legs used 0 of
    their 80 pinned subdivisions while single steps ran for 34 minutes.  A cap
    that ALSO force-accepted would reproduce that exactly.

    So the contract is: a capped update that cannot finish returns failure, and a
    driver that responds the way a D16 subdivision controller does -- take the
    same strain in smaller pieces -- gets through.  Smaller increments need fewer
    substeps, so the cap stops binding.  Both legs travel the identical total
    strain (`_B_E_AX`); only the increment size differs.

    THE CAP IS MEASURED, NOT GUESSED, and a naive guess is wrong in an
    instructive way.  A cap of 1 -- the tightest value that is not "uncapped" --
    looks like the obvious choice and makes the test VACUOUS: the first
    `ModifiedEuler` attempt at `dT = 1` essentially always exceeds `TolE` and
    retries at a smaller `dT`, so **every** plastic update costs at least 2 loop
    iterations no matter how small the increment is, and the fine leg fails too.
    Measured on this deck (uncapped, first deviatoric step at Gauss point 1):

        n_dev  40 ->  3268 substeps      n_dev 160 ->  664
        n_dev  80 ->  1554               n_dev 320 ->  289

    i.e. the cost is very nearly PROPORTIONAL to the increment, which is what
    makes a step-cut a real remedy rather than a hope.  So the cap is set from
    the coarse leg's own measured cost, at 75 % of it: comfortably below the
    coarse step (which must fail) and comfortably above the fine one at ~50 %
    (which must pass).  Both margins are then asserted, so a future deck whose
    scaling is different fails LOUDLY here instead of going vacuous.
    """
    # --- measure the coarse leg's cost, uncapped --------------------------
    _build_brick(1, sani._OPTS_VANILLA)
    _confine(1)
    assert ops.analyze(1) == 0, 'uncapped coarse deviatoric step should converge'
    n_coarse = int(_substeps_at_gp1()[0])

    # --- and the halved-increment leg's ----------------------------------
    _build_brick(2, sani._OPTS_VANILLA, n_dev=2 * _B_N_DEV)
    _confine(2)
    n_fine = 0
    for _ in range(2):                       # two half-steps = one coarse step
        assert ops.analyze(1) == 0, 'uncapped fine deviatoric step should converge'
        n_fine = max(n_fine, int(_substeps_at_gp1()[0]))

    cap = int(0.75 * n_coarse)
    assert n_fine < cap < n_coarse, (
        'the substep cost no longer separates the two increment sizes with a '
        'usable margin, so this gate cannot discriminate: pick a new cap',
        dict(n_coarse=n_coarse, n_fine=n_fine, cap=cap))

    # --- the coarse leg must now FAIL ... ---------------------------------
    _build_brick(3, sani._OPTS_VANILLA + ('-maxSubsteps', cap))
    _confine(3)
    coarse_fail = _push(_B_N_DEV)
    assert coarse_fail == 1, (
        'the capped coarse leg did not fail on the step whose measured cost '
        'exceeds the cap', dict(cap=cap, n_coarse=n_coarse,
                                first_failed_step=coarse_fail))

    # --- ... and the SAME cap must let the halved increment through -------
    _build_brick(4, sani._OPTS_VANILLA + ('-maxSubsteps', cap),
                 n_dev=2 * _B_N_DEV)
    _confine(4)
    fine_fail = _push(2 * _B_N_DEV // 4)      # a quarter of the leg is plenty
    assert fine_fail == 0, (
        'halving the strain increment did NOT get past the substep cap, so the '
        'cap is not something a step-cut can recover from and it has traded a '
        'slow analysis for a dead one',
        dict(cap=cap, n_coarse=n_coarse, n_fine=n_fine,
             first_failed_fine_step=fine_fail))

    print('\nADR-86b T1 subdivision recovery (confine-first LadrunoBrick):')
    print('  uncapped substeps, coarse increment : %d' % n_coarse)
    print('  uncapped substeps, half increment   : %d' % n_fine)
    print('  cap used                            : %d' % cap)
    print('  capped coarse leg fails at step %d; capped fine leg takes %d steps'
          % (coarse_fail, 2 * _B_N_DEV // 4))


# ---------------------------------------------------------------------------
#  Gate 4 -- TanType
# ---------------------------------------------------------------------------

def test_tantype_default_is_the_consistent_tangent(capfd):
    """A deck that emits only the 18 positional parameters gets TanType 2.

    Asserted on the CONSTRUCTION ECHO, not on a behavioural difference, because
    the echo is the thing a user actually reads and the thing that makes the
    fork's divergence from vanilla's default non-silent.  `Print` carries it too
    and is checked as a second, independent channel -- `test_ladruno_sanisand`
    records that deleting `Print` outright once left the whole battery green.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('LadrunoSANISAND', 7, *_PARAMS)      # NO positional optionals
    out = capfd.readouterr()
    echo = out.err + out.out
    assert 'TanType = 2' in echo, (
        'the parser default for TanType is not 2. At 0 the material returns the '
        'ELASTIC tangent from getTangent(), so every deck that omits the '
        'positional block runs algorithm Newton as modified Newton -- measured '
        'at ~7x wall time on a boundary-value problem (ADR-90 GATE U).', echo)
    assert 'consistent' in echo, (
        'the echo names a TanType but does not say what it IS; a number the '
        'reader has to look up is not a disclosure', echo)
    assert 'maxSubsteps = 0' in echo and 'UNCAPPED' in echo, (
        'the echo does not state the substep cap, so a deck cannot be read back '
        'to find out whether one was in force', echo)


def test_explicit_positionals_still_win_over_the_new_default(capfd):
    """`... $Rho 1 0 1 1e-7 1e-7` must still give TanType 0, and 1 must give 1.

    The companion no gate of this kind can do without: an implementation that
    HARDCODED 2 instead of merely DEFAULTING to 2 would pass the test above
    perfectly, and would silently ignore every deck that asked for the elastic
    tangent on purpose -- which the ADR-90 GATE U driver's own A/B leg does, and
    which is the only way anyone can ever re-measure the cost of TanType 0 again.

    Both interior values are checked, not just 0: a "default unless 0" reading of
    the parser would pass a 0-only test.
    """
    for want in (0, 1):
        ops.wipe()
        ops.model('basic', '-ndm', 3, '-ndf', 3)
        ops.nDMaterial('LadrunoSANISAND', 8, *_PARAMS, 1, want, 1, 1.0e-7, 1.0e-7)
        out = capfd.readouterr()
        echo = out.err + out.out
        assert ('TanType = %d' % want) in echo, (
            'an explicit `TanType %d` did not survive the parser -- the new '
            'default is being applied unconditionally rather than as a default'
            % want, echo)


# ---------------------------------------------------------------------------
#  Gate 4b -- and the tangent is worth what it claims, on a FREE-DOF deck
# ---------------------------------------------------------------------------
#
#  Drained triaxial compression on ONE LadrunoBrick.  Negative faces rollered,
#  positive faces LOADED -- so the positive-face DOFs are genuinely free and
#  there is a Newton iteration to count.  That is the entire reason this deck
#  exists: on the zero-free-DOF cubes the rest of the SANISAND battery uses, the
#  tangent is never used to solve anything and TanType is unobservable, which is
#  precisely how this defect survived three PRs.
#
#  `system FullGeneral` because `mCep_Consistent` is UNSYMMETRIC under a
#  non-associated flow rule; a symmetric solver would be solving a different
#  matrix than the one the tangent describes and the comparison would be
#  meaningless.
# ---------------------------------------------------------------------------
#  THE SETTINGS ARE MEASURED, NOT PICKED.  Sweeping dq_total x n_dev x tolerance
#  on this deck (24 legs) gives, as (steps completed of n_dev, total iterations):
#
#    dq   n_dev  tol/P0    TanType 0        TanType 2
#    200    40    1e-2     (40,  539)       (40,  250)
#    200    40    1e-3     (29,  849) DIED  (40,  629)
#    200    40    1e-4     (17,  516) DIED  (38, 1204)
#    120    40    1e-2     (40,  225)       (40,   69)
#    120    40    1e-3     (40,  800)       (40,  283)   <-- pinned
#    120    40    1e-4     (32,  977) DIED  (40,  432)
#    120    80    1e-3     (80, 1192)       (80,  256)
#
#  The row pinned below is the one where BOTH legs complete every step, so the
#  comparison is 40 identical steps against 40 identical steps and not "40 easy
#  ones against 29 hard ones".  Note what the DIED rows say on their own: at a
#  tolerance of 1e-4*P0 the elastic tangent cannot finish the push at all while
#  the consistent one can -- a stronger statement than the iteration count, and
#  deliberately NOT the one gated, because a gate whose evidence is a failure is
#  a gate that passes for the wrong reason the day the failure moves.
_TX_P0 = 100.0          # isotropic confinement, kPa (well away from the low-p
                        # floor -- this test is about the tangent, not p_residual)
_TX_DQ = 120.0          # extra AXIAL pressure applied in stage 1, kPa
_TX_N_CONF = 5
_TX_N_DEV = 40
_TX_TOL_REL = 1.0e-3    # x P0; a force tolerance, per the ADR-86b T3 rule
_TX_MAXITER = 60


def _build_triaxial(tag, positionals):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *positionals,
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
                '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    # Unit cube: a uniform face pressure p is 4 equal nodal forces of p/4.
    q = _TX_P0 / 4.0
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
    # NormUnbalance scaled to a deck-intrinsic force -- the ADR-86b T3 rule.
    # NormDispIncr is unreachable on this material (the substepped return makes
    # the discrete map only piecewise smooth and the residual stalls), so a
    # displacement test here would measure its own tolerance, not the tangent.
    ops.test('NormUnbalance', _TX_TOL_REL * _TX_P0, _TX_MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _TX_N_CONF)
    ops.analysis('Static')


def _tx_state():
    """The converged answer: axial displacement of a free top node and the
    committed stress at Gauss point 1."""
    return (ops.nodeDisp(5, 3),
            list(ops.eleResponse(1, 'material', 1, 'stress')))


def _run_triaxial(tag, positionals):
    """Confine at stage 0, flip, then shear.  Returns (steps_converged,
    total Newton iterations over the CONVERGED deviatoric steps, uz, stress)."""
    _build_triaxial(tag, positionals)
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for _ in range(_TX_N_CONF):
        assert ops.analyze(1) == 0, 'triaxial confinement failed'
    ops.loadConst('-time', 0.0)
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    # stage 1: additional AXIAL pressure only -- drained triaxial compression.
    dq = _TX_DQ / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 1.0 / _TX_N_DEV)

    iters = 0
    done = 0
    for _ in range(_TX_N_DEV):
        if ops.analyze(1) != 0:
            break
        iters += ops.testIter()
        done += 1
    uz, sig = _tx_state()
    return done, iters, uz, sig


def test_tantype_2_costs_fewer_newton_iterations():
    """The consistent tangent is cheaper than the elastic one, MEASURED.

    T2 is a default change, and a default change with no measurement behind it is
    an opinion.  This is the measurement, on the only kind of deck where the
    tangent is observable at all.

    Both legs must complete ALL `_TX_N_DEV` steps -- asserted, not hoped for --
    so the two totals cover the same 40 steps and the comparison can never
    degrade into "40 easy steps against 29 hard ones".  The settings that make
    that true were swept, not guessed; the table is above `_TX_P0`.

    Measured on the dev box at the pinned settings: TanType 0 needs 800
    iterations over 40 steps (20.0 per step), TanType 2 needs 283 (7.1 per
    step) -- **2.8x**.  On a boundary-value problem the same difference showed
    up as ~7x of wall time (ADR-90 GATE U, a strip footing).

    NOT asserted: the ratio.  It depends on the deck, the tolerance and the
    load-step size -- the sweep table shows it running from 3.3x to 4.7x across
    neighbouring settings -- so pinning a number here would be a brittle
    re-measurement of this deck rather than a statement about the tangent.  The
    claim gated is the INEQUALITY, which is the claim the default change rests
    on; the numbers are printed for the record.
    """
    n0, it0, _, _ = _run_triaxial(21, (1, 0, 1, 1.0e-7, 1.0e-7))   # elastic
    n2, it2, _, _ = _run_triaxial(22, (1, 2, 1, 1.0e-7, 1.0e-7))   # consistent

    assert n0 == _TX_N_DEV and n2 == _TX_N_DEV, (
        'the pinned settings no longer get BOTH legs through all %d steps, so '
        'the two totals no longer cover the same steps and the comparison is '
        'between different amounts of work. Re-pin from the sweep table above '
        'rather than comparing unequal legs.' % _TX_N_DEV,
        dict(tantype0=(n0, it0), tantype2=(n2, it2)))

    per0 = it0 / n0
    per2 = it2 / n2
    assert per2 < per0, (
        'TanType 2 did not reduce the Newton iteration count per step. The '
        'parser default was moved 0 -> 2 on the claim that it does; if this '
        'fails, either the claim or the wiring is wrong.',
        dict(tantype0=(n0, it0, per0), tantype2=(n2, it2, per2)))

    # Reported, not asserted -- the number belongs in the PR description, and a
    # future reader should be able to see it without re-deriving it.
    print('\nADR-86b T2 iteration counts (drained triaxial, 1 LadrunoBrick):')
    print('  TanType 0 (elastic)    : %d steps, %d iterations, %.2f per step'
          % (n0, it0, per0))
    print('  TanType 2 (consistent) : %d steps, %d iterations, %.2f per step'
          % (n2, it2, per2))
    print('  speedup in iterations per step: %.2fx' % (per0 / per2))


# The measured agreement between the two tangents at the pinned settings, and the
# floor asserted against it.  MEASURED on the dev box, and the number is stated
# rather than assumed: see the test's docstring for what it is and is not.
#
# Ladruno (ADR-86b review fix, J-1): the original floor here was a flat 1e-5,
# on the theory that it should be "two orders TIGHTER than the per-step
# NormUnbalance tolerance".  That theory is wrong for an ACCUMULATED,
# path-dependent comparison: this deck is 40 LoadControl steps, and
# SANISAND's fabric/backstress evolution feeds each step's accepted state
# forward into every later step, so a per-step residual difference does not
# stay a per-step-sized answer difference -- it compounds.  Measured directly
# (sweep the deck's own `_TX_TOL_REL`):
#   tol = 1e-2 x P0  ->  d_uz = 3.43e-2   (~3.4x tol)
#   tol = 1e-3 x P0  ->  d_uz = 4.46e-3   (~4.5x tol)   <-- pinned tolerance
# i.e. the ANSWER discrepancy tracks the TOLERANCE roughly linearly, at
# somewhat MORE than 1x tol, not two orders below it. The floor is therefore
# set relative to the deck's own tolerance, not as an independent constant:
# 10x the per-step NormUnbalance tolerance, which still has >2x margin over
# the measured 4.5x-tol discrepancy while remaining tight enough to catch an
# actual TanType-2 defect (a wrong consistent tangent does not merely take a
# different path to the same equilibrium -- it converges Newton to a
# DIFFERENT stress state, which is a qualitatively larger effect than path
# accumulation noise).
_TX_ANSWER_FLOOR = 10.0 * _TX_TOL_REL


def test_tantype_does_not_change_the_converged_answer():
    """T2 is cost-warranted; this is the ANSWER warrant, and it is separate.

    Moving a default because it is CHEAPER is only legitimate if it is also the
    SAME.  The argument that it is -- "the tangent sets the iteration path, not
    the equilibrium point" -- is true only for a FORCE-residual convergence test,
    and `LEDGER_quirks` records the fork's own counter-example: under
    `NormDispIncr` / `EnergyIncr` a different `K_f` accepts at a different point,
    because those tests accept on `dU = K_f^-1 R`.  This deck uses
    `NormUnbalance`, which is exactly the case where the argument holds -- so
    this test is that argument's receipt, not a restatement of it.

    Both legs run the identical deck, the identical unsymmetric solver
    (`FullGeneral` -- `mCep_Consistent` is unsymmetric under a non-associated
    flow rule) and the identical `NormUnbalance` tolerance; only `TanType`
    differs.  Compared at the END of the push, on two independent quantities: a
    free node's axial displacement and the committed Gauss-point stress.

    THE FLOOR IS MEASURED, NOT GUESSED.  Convergence is declared at
    `NormUnbalance = 1e-3 x P0`, i.e. a residual of 0.1 against ~100 kN of
    applied load -- a relative equilibrium tolerance of 1e-3.  Two runs accepted
    anywhere inside that ball may legitimately differ, and on this deck they do:
    over the 40 accumulated LoadControl steps the per-step residual difference
    compounds through SANISAND's path-dependent fabric/backstress evolution, so
    the answer discrepancy tracks the deck's OWN tolerance (measured ~4.5x it at
    the pinned `1e-3`, ~3.4x it at `1e-2` -- see `_TX_ANSWER_FLOOR`'s comment),
    not some fraction of it. The floor is set at `10 x _TX_TOL_REL`, which is
    still comfortable margin over the measured discrepancy and the honest place
    to put it given the measurement. The measured value is printed.

    NOT asserted: bit-identity.  These are different iteration paths reaching the
    same point, not the same arithmetic, and demanding equality would be a gate
    that fails on a compiler flag.
    """
    n0, _, uz0, sig0 = _run_triaxial(31, (1, 0, 1, 1.0e-7, 1.0e-7))
    n2, _, uz2, sig2 = _run_triaxial(32, (1, 2, 1, 1.0e-7, 1.0e-7))

    assert n0 == n2 == _TX_N_DEV, (
        'the two legs did not both complete the push, so they are not at the '
        'same point on the load path and no comparison of the answer is '
        'possible', n0, n2)

    d_uz = abs(uz2 - uz0) / max(abs(uz0), 1.0e-30)
    d_sig = sani._reldiff(sig0, sig2)

    print('\nADR-86b T2 answer gate (same deck, same solver, same tolerance):')
    print('  axial displacement : %.6e vs %.6e   rel %.3e' % (uz0, uz2, d_uz))
    print('  GP1 stress reldiff : %.3e' % d_sig)
    print('  floor              : %.1e (the run converged at NormUnbalance '
          '%.0e x P0)' % (_TX_ANSWER_FLOOR, _TX_TOL_REL))

    assert d_uz < _TX_ANSWER_FLOOR and d_sig < _TX_ANSWER_FLOOR, (
        'TanType 0 and TanType 2 converged to DIFFERENT answers on the same '
        'deck. The default was moved 0 -> 2 on the premise that the tangent '
        'sets the iteration path and not the equilibrium point; if this fails, '
        'that premise is wrong here and the default change is not warranted.',
        dict(uz=(uz0, uz2, d_uz), stress_reldiff=d_sig,
             floor=_TX_ANSWER_FLOOR))
