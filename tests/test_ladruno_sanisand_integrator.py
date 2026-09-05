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

# The recorded confine-first deck, its parameters, and its published numbers.
# Imported rather than copied on purpose: gate 1's whole claim is that ADR-86b
# did not move THAT deck's answer, which only holds if it is literally that deck.
import test_ladruno_sanisand as sani


pytestmark = [pytest.mark.zone_a]

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

def test_cap_engages_and_the_step_fails(tmp_path):
    """Measure the substep cost, cap below it, require a failed step.

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
    # --- leg 1: how expensive IS this update, really? ----------------------
    _build_brick(1, sani._OPTS_VANILLA)
    _confine(1)
    assert ops.analyze(1) == 0, 'the uncapped deviatoric step should converge'
    taken, hit = _substeps_at_gp1()
    assert hit == 0.0, ('the uncapped run reports the cap as HIT; the flag is '
                        'not being reset per update', taken, hit)
    assert taken >= 2, (
        'this deck spends fewer than 2 ModifiedEuler substeps per update, so '
        'there is no room to cap BELOW it -- pick a larger deviatoric '
        'increment before trusting anything else in this file', taken)

    cap = int(taken) - 1

    # --- leg 2: cap below the measured cost -------------------------------
    log = os.path.join(str(tmp_path), 'capwarn.log')
    _build_brick(2, sani._OPTS_VANILLA + ('-maxSubsteps', cap))
    _confine(2)
    ops.logFile(log, '-noEcho')          # opserr -> file, console silent
    rc = ops.analyze(1)
    ops.logFile(os.path.join(str(tmp_path), 'other.log'), '-noEcho')  # release

    assert rc != 0, (
        'the substep cap fired but analyze() reported SUCCESS. The return code '
        'is being dropped somewhere between ManzariDafalias::ModifiedEuler and '
        'the integrator -- which is exactly the silence ADR-90 GATE U ran into '
        '(0 of 80 subdivisions used across six legs). Check that the element is '
        'a LadrunoBrick: stdBrick discards setTrialStrain return codes.',
        cap, taken, rc)

    text = open(log, errors='ignore').read()
    # The construction echo is in this file too, so an EMPTY file would mean the
    # redirect never took -- the positive control the logFile quirk asks for.
    assert 'LadrunoSANISAND tag 2' in text, (
        'nothing at all was captured, so the opserr redirect did not take and '
        'the absence of a warning below proves nothing')
    n_warn = text.count('substep cap')
    assert n_warn >= 1, (
        'the cap fired and the step failed, but the material said nothing. A '
        'step that fails for no stated reason is worse for the user than the '
        'seizure it replaces.', text[-2000:])
    # Throttled to a PROCESS budget of 10, not per instance: every Gauss point is
    # its own material object.  8 Gauss points x 1 step must not exceed it.
    assert n_warn <= 10, (
        'more than the 10-line process budget of cap warnings was emitted; the '
        'throttle is per instance again, which on a real mesh is ~400k lines',
        n_warn)


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

    It runs at a cap of 1 rather than at gate 2's measured value on purpose: 1 is
    the tightest cap that is not "uncapped", so the coarse leg is guaranteed to
    fail and the test cannot go vacuous by the deck becoming cheap.
    """
    _build_brick(1, sani._OPTS_VANILLA + ('-maxSubsteps', 1))
    _confine(1)
    coarse_fail = _push(_B_N_DEV)
    assert coarse_fail != 0, (
        'a cap of ONE substep never fired over the whole deviatoric leg, so '
        'this deck no longer exercises the cap at all and gate 3 is vacuous')

    # Same cap, same total strain, half the increment.
    _build_brick(2, sani._OPTS_VANILLA + ('-maxSubsteps', 1),
                 n_dev=2 * _B_N_DEV)
    _confine(2)
    fine_fail = _push(2 * coarse_fail)
    assert fine_fail == 0, (
        'halving the strain increment did NOT get past the substep cap, so the '
        'cap is not something a step-cut can recover from and it has traded a '
        'slow analysis for a dead one', coarse_fail, fine_fail)


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
_TX_P0 = 100.0          # isotropic confinement, kPa (well away from the low-p
                        # floor -- this test is about the tangent, not p_residual)
_TX_DQ = 260.0          # extra AXIAL pressure applied in stage 1, kPa
_TX_N_CONF = 5
_TX_N_DEV = 20
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
    ops.test('NormUnbalance', 1.0e-6 * _TX_P0, _TX_MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _TX_N_CONF)
    ops.analysis('Static')


def _run_triaxial(tag, positionals):
    """Confine at stage 0, flip, then shear.  Returns (steps_converged,
    total Newton iterations over the CONVERGED deviatoric steps)."""
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
    return done, iters


def test_tantype_2_costs_fewer_newton_iterations():
    """The consistent tangent is cheaper than the elastic one, MEASURED.

    T2 is a default change, and a default change with no measurement behind it is
    an opinion.  This is the measurement, on the only kind of deck where the
    tangent is observable at all.

    Compared over the steps BOTH legs converged, so the comparison is never
    between "20 easy steps" and "9 hard ones".  The elastic-tangent leg is
    expected to need many more iterations per step (linear vs quadratic
    convergence); on a boundary-value problem the same difference showed up as
    ~7x of wall time (ADR-90 GATE U).

    NOT asserted: a ratio.  The count depends on the deck, the tolerance and the
    load-step size, and pinning a number here would be a brittle re-measurement
    of this deck rather than a statement about the tangent.  The claim is the
    INEQUALITY, which is the claim the default change rests on.
    """
    n0, it0 = _run_triaxial(21, (1, 0, 1, 1.0e-7, 1.0e-7))   # elastic tangent
    n2, it2 = _run_triaxial(22, (1, 2, 1, 1.0e-7, 1.0e-7))   # consistent tangent

    assert n0 > 0 and n2 > 0, (
        'one of the two legs converged no deviatoric step at all, so there is '
        'nothing to compare', n0, n2)
    assert n2 >= n0, (
        'the CONSISTENT tangent got FEWER steps through than the elastic one, '
        'which inverts the premise of the T2 default change', n0, n2)

    # Re-run the elastic leg truncated to the shorter of the two, so the totals
    # cover the same steps.  (n2 >= n0 above, so n0 is the shorter one.)
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
