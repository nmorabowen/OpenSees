"""WP-108 -- `-maxSubsteps` IS honoured on IntScheme 2; the "NO EFFECT" warning
was wrong.

BACKGROUND (WP-105 / F12,
`Ladruno_files/testbed/hypo_bearing/adr92_f12/F12_intscheme2_verdict.md`):
`LadrunoSANISAND::schemeReachesModifiedEuler()` used to return `false` for
`mScheme == 2` (`INT_BackwardEuler`), so the constructor and `Print()` both
warned that `-maxSubsteps`/`-honorTolR` have NO EFFECT with IntScheme 2. That
was measured false: `ManzariDafalias::explicit_integrator`'s `switch(mScheme)`
(`ManzariDafalias.cpp` ~1070-1101) does not enumerate `INT_BackwardEuler`, so
whenever `BackwardEuler_CPPM`'s own recursive-halving retry ladder
(`ManzariDafalias.cpp` ~2472-2588) falls back to `explicit_integrator` --
on non-convergence or ladder exhaustion -- that call hits the switch's
`default:` case, which is `ModifiedEuler`. `-maxSubsteps` and `-honorTolR` are
both read at exactly one site, inside `ModifiedEuler()`, so the cap DOES bind
on scheme 2, conditionally (only on the fallback path, not on every step the
way it does for schemes 0/1). WP-105 measured this directly on a `p -> p_min`
floor path: a `-maxSubsteps 100` cap turned a run that completed all 40 steps
uncapped (up to 1282 substeps) into one that refuses partway through.

WP-108 fixes `schemeReachesModifiedEuler()` to return `true` for scheme 2. This
file is the reproduction of WP-105's own control, run against the fix:

  1. `test_no_false_no_effect_warning_on_intscheme2` -- the warning itself.
     Constructing an IntScheme 2 deck with `-maxSubsteps` set must NOT print
     "has NO EFFECT with IntScheme 2". Run against the PRE-WP-108 binary this
     assertion FAILS (the false warning prints) -- that is the negative
     control the PR description cites, not something this file can assert
     against a binary it isn't running.
  2. `test_maxsubsteps_is_honoured_on_intscheme2` -- the seam itself, not just
     the warning. Same floor path WP-105 used: a generous
     (effectively-uncapped) run must complete every step, and a starved cap
     (100, WP-105's own number) must fail a step before the uncapped run's own
     measured substep count is reached. If `schemeReachesModifiedEuler()`
     regressed back to `false`, the cap would stop being wired at all
     (`applyLadrunoConstants()` still sets `mMaxSubstepsInME = mMaxSubsteps`
     regardless -- the defect fixed here is only the WARNING's own truthiness
     -- but this second gate is kept anyway as a belt-and-suspenders check on
     the seam WP-105 measured, not merely on the string).

The floor path itself (`_run_floor` below) is adapted from WP-105/F12's own
`f12_matpoint.py::run_floor`: a prescribed volumetric-extension + deviatoric
shear path with EVERY positive-face DOF `sp`-prescribed, so there are ZERO
free equations and nothing but the constitutive integrator can make a step
fail -- exactly the same zero-free-DOF discipline the rest of the SANISAND
battery uses (see `test_ladruno_sanisand_integrator.py`'s header note).

COST: two 40-step zero-free-DOF matpoint pushes (plus a short consolidation
each); a few seconds on the dev box, well under the file's 60 s budget.
"""
import pytest

from _testbed import ops

import test_ladruno_sanisand as sani

pytestmark = [pytest.mark.zone_a]

_PARAMS = sani._PARAMS
_P_ATM = sani._P_ATM

# Same constants WP-105/F12 used (f12_matpoint.py: PMIN, PRESIDUAL, HONOR_TOLR).
_PMIN = 1.0e-4 * _P_ATM
_PRESIDUAL = 0.0
_HONOR_TOLR = 0

# IntScheme 2, TanType 2 (consistent -- irrelevant here, kept explicit so the
# positional block is unambiguous), JacoType 1, TolF/TolR 1e-10.
_SCHEME2_POSITIONALS = (2, 2, 1, 1.0e-10, 1.0e-10)

_NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
          (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]

_P0 = 100.0
_NSTEP = 40
_EVOL_MAX = 2.2e-4
_GAM_MAX = 6.0e-3


def _build_model(tag, opts=()):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for t, (x, y, z) in enumerate(_NODES, start=1):
        ops.node(t, float(x), float(y), float(z))
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 0, 0)
    for n in (1, 2, 5, 6):
        ops.fix(n, 0, 1, 0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 0, 0, 1)
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *_SCHEME2_POSITIONALS,
                   '-Presidual', _PRESIDUAL, '-Pmin', _PMIN,
                   '-honorTolR', _HONOR_TOLR, *opts)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag, '-formulation', 'bbar')
    q4 = -_P0 / 4.0
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, q4, 0.0, 0.0)
    for n in (3, 4, 7, 8):
        ops.load(n, 0.0, q4, 0.0)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, q4)


def _consolidate(tag):
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-8, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')
    assert ops.analyze(10) == 0, 'consolidation failed'
    ops.updateMaterialStage('-material', tag, '-stage', 1)
    ops.integrator('LoadControl', 0.0)
    assert ops.analyze(5) == 0, 'stage-1 re-equilibration failed'
    ops.loadConst('-time', 0.0)


def _floor_path_series(nstep=_NSTEP):
    """WP-105/F12's `run_floor` path: volumetric extension ramped to
    `_EVOL_MAX` by t=0.5 and held, deviatoric shear `_GAM_MAX` on throughout --
    the ADR-93 ring regime that drives p onto the p_min floor while the point
    keeps flowing. A Path series returns 0 past its last time point and
    LoadControl's accumulated pseudo-time overshoots 1.0 by an ulp on the last
    step, so the series is extended flat rather than left to unload on the
    final step."""
    ts = [i / float(nstep) for i in range(nstep + 1)]
    ramp = [min(2.0 * t, 1.0) for t in ts]
    ex = [_EVOL_MAX / 3.0 * r + _GAM_MAX / 3.0 * t for r, t in zip(ramp, ts)]
    ez = [_EVOL_MAX / 3.0 * r - 2.0 * _GAM_MAX / 3.0 * t for r, t in zip(ramp, ts)]
    ts = ts + [2.0]
    ex = ex + [ex[-1]]
    ez = ez + [ez[-1]]
    return ts, ex, ez


def _run_floor(tag, opts=(), nstep=_NSTEP):
    """Push the zero-free-DOF floor path; return the number of steps that
    converged (0..nstep)."""
    _build_model(tag, opts)
    _consolidate(tag)
    ts, ex, ez = _floor_path_series(nstep)
    ops.timeSeries('Path', 3, '-time', *ts, '-values', *ex)
    ops.timeSeries('Path', 4, '-time', *ts, '-values', *ez)
    ops.pattern('Plain', 3, 3)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, 1.0)
    for n in (3, 4, 7, 8):
        ops.sp(n, 2, 1.0)
    ops.pattern('Plain', 4, 4)
    for n in (5, 6, 7, 8):
        ops.sp(n, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-9, 200, 0)
    ops.algorithm('KrylovNewton')
    ops.integrator('LoadControl', 1.0 / nstep)
    ops.analysis('Static')
    done = 0
    for i in range(1, nstep + 1):
        if ops.analyze(1) != 0:
            return done
        done = i
    return done


# ---------------------------------------------------------------------------
#  Gate 1 -- the warning itself
# ---------------------------------------------------------------------------

def test_no_false_no_effect_warning_on_intscheme2(capfd):
    """Building an IntScheme 2 deck with `-maxSubsteps` set must not print the
    false "-maxSubsteps N has NO EFFECT with IntScheme 2" warning.

    Run against the pre-WP-108 binary this assertion FAILS -- the warning
    prints on construction (`LadrunoSANISAND.cpp` ~1206-1213) -- which is the
    negative control cited in the PR description; this file only runs against
    one binary at a time, so it cannot itself flip that switch.

    `-honorTolR` is left at its default (0) on purpose, so the sibling
    `-honorTolR` warning (gated on the identical predicate,
    `LadrunoSANISAND.cpp` ~1187-1199) cannot fire here regardless and cannot
    be mistaken for the one under test.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('LadrunoSANISAND', 1, *_PARAMS, *_SCHEME2_POSITIONALS,
                   '-maxSubsteps', 100)
    out = capfd.readouterr()
    echo = out.err + out.out
    assert 'LadrunoSANISAND tag 1' in echo, (
        'nothing was captured from construction, so the absence of the '
        'warning below proves nothing', echo)
    assert 'has NO EFFECT with IntScheme 2' not in echo, (
        'the false "-maxSubsteps has NO EFFECT with IntScheme 2" warning is '
        'still being printed -- schemeReachesModifiedEuler() regressed', echo)
    assert 'has NO EFFECT with IntScheme' not in echo, (
        'a NO EFFECT warning printed for some OTHER scheme number while '
        'building an IntScheme 2 deck -- the echoed scheme number is wrong',
        echo)


# ---------------------------------------------------------------------------
#  Gate 1b -- the SAME false claim, via Print(), not just the constructor echo
#
#  `LadrunoSANISAND::Print()` carries its OWN pair of `!schemeReachesModifiedEuler()`
#  checks (`LadrunoSANISAND.cpp` ~4345 and ~4383, one for `-maxSubsteps` and one
#  for `-honorTolR`) -- a SEPARATE call site from the constructor's, reached only
#  by `ops.printModel(...)`, not by construction alone. Both route through the
#  identical `schemeReachesModifiedEuler()` helper the constructor path uses, so
#  fixing the one function fixes both call sites -- this test is the proof that
#  it does, not a second fix.
# ---------------------------------------------------------------------------

def test_no_false_no_effect_note_in_printmodel(tmp_path):
    """`printModel` must not carry the false "does not route to ModifiedEuler(),
    so -maxSubsteps is INERT on this deck" NOTE for an IntScheme 2 material with
    `-maxSubsteps` set.

    Companion to `test_no_false_no_effect_warning_on_intscheme2` above, which
    only drives the CONSTRUCTOR's echo. `Print()` is a different function
    (deliberately left unguarded so any instance can be interrogated, per
    `test_ladruno_sanisand.py::test_print_states_what_it_ran`), with its own
    call sites into the same `schemeReachesModifiedEuler()` helper -- this is
    the check that fixing the helper closed BOTH doors, not just the one this
    file's other test happens to exercise.
    """
    _build_model(3, opts=('-maxSubsteps', 100))
    out = tmp_path / 'printmodel_intscheme2.out'
    ops.printModel('-file', str(out), '-ele')
    txt = out.read_text()

    assert 'LadrunoSANISAND Material, tag: 3' in txt, (
        'printModel produced no LadrunoSANISAND record at all -- the positive '
        'control this assertion needs is missing, so the absence of the NOTE '
        'below proves nothing', txt)
    assert 'does not route to' not in txt, (
        "printModel still carries the false \"does not route to ModifiedEuler()\" "
        "NOTE for IntScheme 2 -- schemeReachesModifiedEuler() regressed, or "
        "Print()'s call site was not covered by the fix", txt)
    assert 'is INERT on this deck' not in txt, (
        'printModel still claims -maxSubsteps/-honorTolR is INERT on this '
        'IntScheme 2 deck', txt)


# ---------------------------------------------------------------------------
#  Gate 2 -- the seam itself: the cap actually binds on scheme 2
# ---------------------------------------------------------------------------

def test_maxsubsteps_is_honoured_on_intscheme2():
    """WP-105/F12's own control, reproduced: a generous run completes every
    step of the floor path; a `-maxSubsteps 100` cap on the SAME path fails a
    step before that.

    Self-calibrating against the uncapped leg's own step count rather than a
    hardcoded "must fail before step N": WP-105 measured the cap biting at
    step 18 of 40 on its build, but the exact step is a property of the
    substep cost at each point along the path, not a promise this fix makes.
    What IS promised, and what would go back to being silently false if
    `schemeReachesModifiedEuler()` regressed, is that the cap can fail a step
    on this scheme AT ALL.
    """
    done_uncapped = _run_floor(1)
    assert done_uncapped == _NSTEP, (
        'the uncapped IntScheme 2 floor path did not complete all %d steps; '
        're-measure before trusting the capped leg below' % _NSTEP,
        done_uncapped)

    done_capped = _run_floor(2, opts=('-maxSubsteps', 100))
    assert done_capped < _NSTEP, (
        'a -maxSubsteps 100 cap on IntScheme 2 did not fail a single step of '
        'the same floor path the uncapped leg just completed -- the cap is '
        'not reaching this scheme, which is the exact defect WP-108 fixed '
        '(schemeReachesModifiedEuler() returning false for scheme 2)',
        done_capped)

    print('\nWP-108 IntScheme 2 -maxSubsteps control (floor path, p0=%.0f):'
          % _P0)
    print('  uncapped : %d/%d steps completed' % (done_uncapped, _NSTEP))
    print('  cap 100  : %d/%d steps completed (first failure at step %d)'
          % (done_capped, _NSTEP, done_capped + 1))
