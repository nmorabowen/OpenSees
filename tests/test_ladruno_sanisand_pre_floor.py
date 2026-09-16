"""WP-106 / ADR-93 II.1 -- `-pRe`, the ELASTIC-ONLY confinement floor.

ADR 93 section 1 row 1, source-verified: `m_Presidual` reaches ~30 plastic-side
mean-stress sites and `GetElasticModuli` at NONE of them, so `G, K -> 0` with `p`
whatever `p_r` is. `-pRe` is its mirror image -- it enters the three
`GetElasticModuli` overloads as `p + pRe` (before the `m_Pmin` clamp) and nothing
else.

What this module gates, in the order a reviewer should read it:

1. **OFF is byte-identical.** `-pRe 0` and the flag omitted give the same
   committed stress as each other AND as the pre-ADR-93 deck -- to exact
   equality, because the default re-executes no arithmetic (the parameter enters
   as `+ 0.0` and the one derived quantity recomputed for it sits behind an early
   return).
2. **ON is not vacuous.** The same deck at `-pRe 1.0` gives a DIFFERENT answer,
   measured, so assertion 1 has something to catch.
3. **It reaches every Gauss point.** The deck-level material is a prototype;
   `getCopy(const char*)` builds the working copy per integration point, and
   `-pRe` is not a WRAPPER constructor argument. If the transfer were missing,
   assertion 2 would read ZERO -- so 2 gates the clone too, and this test makes
   that reasoning explicit by running the plane-strain lane as well.
4. **It is a STIFFNESS floor, not a strength floor.** On a low-confinement
   drained triaxial the mobilised `eta / M^b` at the end of the path must move by
   far less than `-Presidual` moves it -- the floor may change the road, not the
   destination.
5. **It survives `reset()`** (`revertToStart` routes through the base
   `initialize()`, which zeroes the seam) and **crosses the FileDatastore wire**
   (slot 34 of the fork `Vector(35)`; ADR-86 section 3 is exactly this defect in
   `m_Presidual`).
6. **Refusals and the echo.** A negative value is refused at parse time; a
   REPEATED `-pRe` is refused rather than silently last-wins; a value above
   `0.1*P_atm` warns; a value at or below `p_min` prints a NOTE (the clamp
   already dominates as `p -> 0`, so the ring gains nothing); the construction
   echo names the floor and says which side of the model it touches.
7. **It is PINNED, not merely "different".** `_PIN_3D_PRE1` / `_PIN_PS_PRE1`
   fix the committed stress at `pRe = 1.0` to 1e-6 relative, on both
   `getCopy(const char*)` branches. A one-sided `> 1e-6` inequality passes for
   any non-zero wiring; the pin fixes the value, and the `_PERTURBED` leg beside
   it fixes the resolution. **One named mutant it does NOT kill** -- "the floor
   applied AFTER the `m_Pmin` clamp" -- was built and run and is byte-identical
   to the correct build, because that clamp never fires on a staged deck. See
   `test_pre_floor_value_is_pinned_not_merely_different`; the ordering is pinned
   in the SOURCE instead.
8. **The initial elastic operator scales with the floor** —
   `eigen` on an unstrained brick at `pRe = 3*P_atm` is exactly 2x the unfloored
   one, mode by mode, once the material is at stage 1. That is the only
   observable `refreshInitialElasticOperator()` has, and the test says why.
9. **It is INERT in the ELASTIC stage, and that is correct.** `mElastFlag == 0`
   selects the branch of `GetElasticModuli` in which the `sqrt(pn/P_atm)` factor
   is dropped, so the gravity stage is pressure-INDEPENDENT and there is nothing
   for a confinement floor to floor. Pinned with `-Pmin` as the control: the
   clamp is inert there too, so the property belongs to the vanilla stage, not
   to this flag.

Deck and helpers are imported from `test_ladruno_sanisand`, deliberately: this
file must exercise the SAME strain path the p_residual gates use, or "pRe moves
the answer by X and p_r by Y" is a comparison of two different things.

WALL TIME (conftest.py's `slow` contract): the default tier is **0.2 s** of test
time (14 passed); the one `slow` case,
`test_pre_floor_moves_stiffness_far_more_than_strength`, is **~104 s measured**
(103.7 s and 104.3 s on two `--durations` runs; two 1200-step low-`p` triaxial
legs through `_drained_triaxial`), so `--runslow` costs **104.6 s** for the
module.
"""
import math
import os
import tempfile

import pytest

from _testbed import ops

import test_ladruno_sanisand as S

pytestmark = [pytest.mark.zone_a]

_EQ_TOL = 1.0e-15          # exact: the OFF path re-executes no arithmetic
_PMIN = S._VANILLA_PMIN
_P_ATM = S._P_ATM

_OPTS_BASE = ("-Presidual", 0.0, "-Pmin", _PMIN, "-honorTolR", 0)


def _opts(pre=None):
    return _OPTS_BASE if pre is None else _OPTS_BASE + ("-pRe", float(pre))


# ---------------------------------------------------------------------------
#  PINS.  Measured on `cc4aa6db0`+ (this branch), printed at %.12e.
#
#  Why pinned and not just "different": every ON-side assertion in the original
#  draft of this file was a one-sided inequality, so ANY non-zero wiring passed
#  it.  The pin fixes the VALUE and, with the `_PERTURBED` leg beside it, the
#  RESOLUTION: a 1 % change in the request lands `2.778e-4` away, 278x outside
#  the pin, and 1e-6 is still 1000x looser than the run-to-run noise on this
#  deck (the OFF leg reproduces to 0.0 exactly).
#
#  Cross-platform: these are committed stresses after 25 Newton-converged steps
#  at `NormDispIncr 1e-13`, i.e. the solver tolerance -- not bit patterns.  1e-6
#  is the same bar ADR-94 settled on for cross-platform value pins.
# ---------------------------------------------------------------------------

_PIN_TOL = 1.0e-6
_MUTANT_BAR = 1.0e-5       # the pin must MISS by at least this on the mutant

#: `S._drive('LadrunoSANISAND', tag, _opts(1.0))`, normal components.
_PIN_3D_PRE1 = (-4.036730510440e+00, -4.036730510440e+00, -3.405550151613e+01)
#: `S._drive_ps(tag, _opts(1.0))` -- the PlaneStrain wrapper's 3-vector, in-plane
#: normal components.
_PIN_PS_PRE1 = (-1.018711858576e+01, -4.395168700731e+01)
#: A 1 % perturbation of the request (numerically, the upper bound of what the
#: "applied after the m_Pmin clamp" mutant could ever shift `pn` by). Used to
#: demonstrate the pin's RESOLUTION -- see the test docstring for why that
#: mutant itself is unobservable at run time.
_PERTURBED = 1.0 + _PMIN


def _worst_rel(measured, pinned):
    """max |x - pin| / |pin| over the pinned components."""
    return max(abs(x - p) / abs(p) for x, p in zip(measured, pinned))


# ---------------------------------------------------------------------------
#  1 + 2 -- OFF is identical, ON is not
# ---------------------------------------------------------------------------

def test_pre_floor_off_is_identical_and_on_is_not():
    """`-pRe 0` == flag omitted == the pre-ADR-93 deck; `-pRe 1.0` != any of them."""
    omitted = S._drive('LadrunoSANISAND', 101, _opts())
    zero = S._drive('LadrunoSANISAND', 102, _opts(0.0))
    floored = S._drive('LadrunoSANISAND', 103, _opts(1.0))

    assert S._reldiff(omitted, zero) <= _EQ_TOL, (
        '-pRe 0 is not the same deck as omitting the flag; the default is '
        'supposed to re-execute no arithmetic at all', omitted, zero)

    gap = S._reldiff(omitted, floored)
    assert gap > 1.0e-6, (
        'this deck cannot tell pRe = 0 from pRe = 1.0 kPa, so the equality '
        'above would pass with the parameter wired to nothing -- and, because '
        'every Gauss point is a getCopy(const char*) clone, it would also pass '
        f'if the clone dropped the request (gap {gap:.3e})')


def test_pre_floor_reaches_the_planestrain_lane():
    """The PlaneStrain wrapper is a second `getCopy(const char*)` branch."""
    base = S._drive_ps(201, _opts())
    floored = S._drive_ps(202, _opts(1.0))
    same = S._drive_ps(203, _opts(0.0))

    assert S._reldiff(base, same) <= _EQ_TOL, (
        'PlaneStrain: -pRe 0 moved the answer', base, same)
    assert S._reldiff(base, floored) > 1.0e-6, (
        'PlaneStrain: -pRe 1.0 did NOT move the answer -- the PlaneStrain2D '
        'branch of getCopy(const char*) is not carrying the request')


# ---------------------------------------------------------------------------
#  7 -- the VALUE is pinned, so a wrongly-wired floor cannot pass
# ---------------------------------------------------------------------------

def test_pre_floor_value_is_pinned_not_merely_different():
    """`-pRe 1.0` must reproduce a MEASURED stress, not merely move it.

    Mutants this pin kills:

    * **the floor applied to `K` only or to `G` only.** `K` is derived from `G`
      one line later, so either mutation rescales the two independently and the
      committed stress moves far more than 1e-6. (The eigen test below pins the
      same thing structurally: every mode scales by the *same* factor.)
    * **the request dropped in one `getCopy(const char*)` branch.** The 3D and
      PlaneStrain pins are separate values from separate wrappers.
    * **any magnitude error** -- the second half runs `pRe = 1 + p_min`, a
      perturbation of 2.778e-4, and asserts the pin rejects it. That fixes the
      pin's *resolution*, which is the property a one-sided `> 1e-6` inequality
      never had.

    A mutant this pin does NOT kill, and no runtime pin can:
    **`+ m_PreElastic` applied AFTER the `m_Pmin` clamp.** It was BUILT and RUN
    (all three overloads changed to `pn = ((pn - m_PreElastic) <= m_Pmin) ?
    (m_Pmin + m_PreElastic) : pn;`, binary rebuilt) and it is byte-identical to
    the correct build on this deck and on `-Pmin` 10.0 / 30.0 variants. Reason,
    measured: the clamp inside `GetElasticModuli` never fires on a staged deck.
    At stage 0 the branch that reads `pn` is not taken at all; at stage 1
    `Stress_Correction`'s low-`p` rescue holds the committed `p` at or above
    `m_Pmin` from the first plastic step (at `-Pmin 10` the five sub-`p_min`
    steps are all stage-0 ones, and `p` jumps 5.725 -> 11.15 at the flip). So
    the ordering is a PROVENANCE choice -- it matches the ADR-93 numpy oracle,
    which is where every II.1 number came from -- and not a behavioural one.
    `test_the_seam_is_at_every_getelasticmoduli_overload` pins it in the source,
    which is the only place it is observable.
    """
    floored = S._drive('LadrunoSANISAND', 161, _opts(1.0))
    worst = _worst_rel(floored[:3], _PIN_3D_PRE1)
    assert worst <= _PIN_TOL, (
        f'-pRe 1.0 no longer reproduces the pinned stress (worst {worst:.3e} '
        f'> {_PIN_TOL:.0e}); measured {floored[:3]}, pinned {_PIN_3D_PRE1}')

    mutant = S._drive('LadrunoSANISAND', 162, _opts(_PERTURBED))
    mutant_worst = _worst_rel(mutant[:3], _PIN_3D_PRE1)
    assert mutant_worst > _MUTANT_BAR, (
        f'a 1 % change in pRe is indistinguishable at this pin '
        f'({mutant_worst:.3e}) -- the pin has been loosened past the resolution '
        'it exists to provide')

    ps = S._drive_ps(163, _opts(1.0))
    ps_worst = _worst_rel(ps[:2], _PIN_PS_PRE1)
    assert ps_worst <= _PIN_TOL, (
        f'PlaneStrain: -pRe 1.0 no longer reproduces the pinned stress '
        f'(worst {ps_worst:.3e}); measured {ps[:2]}, pinned {_PIN_PS_PRE1}')

    ps_mutant = S._drive_ps(164, _opts(_PERTURBED))
    assert _worst_rel(ps_mutant[:2], _PIN_PS_PRE1) > _MUTANT_BAR, (
        'PlaneStrain: the after-the-clamp proxy passes the pin')


def test_the_seam_is_at_every_getelasticmoduli_overload():
    """A SOURCE pin, because the third overload has no live caller.

    `ManzariDafalias::GetElasticModuli(sigma, en, en1, nEStrain, cEStrain, K, G)`
    is dead in this fork -- both candidate call sites are commented out -- so
    "the floor was applied in only two of the three overloads" is unkillable by
    any runtime deck. It is still a real drift risk (a future caller would get
    the unfloored law silently), so it is pinned where it lives: the count of
    the seam line in the vanilla source.
    """
    src = os.path.join(S._ROOT, 'SRC', 'material', 'nD', 'UWmaterials',
                       'ManzariDafalias.cpp')
    with open(src, 'r', encoding='utf-8', errors='replace') as fh:
        text = fh.read()
    seam = 'double pn = one3 * GetTrace(sigma) + m_PreElastic;'
    n = text.count(seam)
    assert n == 3, (
        f'expected the ADR-93 II.1 seam at all three GetElasticModuli overloads, '
        f'found {n} -- the overloads have drifted apart')
    assert text.count('double pn = one3 * GetTrace(sigma);') == 0, (
        'an UNfloored `pn` line survives in ManzariDafalias.cpp: one of the '
        'three overloads is running the vanilla elastic law')


# ---------------------------------------------------------------------------
#  8 + 9 -- the initial elastic operator, and the elastic-stage inertness
# ---------------------------------------------------------------------------

def _eigen_brick(pre, n_modes, stage, do_reset):
    """One unstrained `stdBrick`, bottom face fully fixed, no prescribed strain.

    `eigen` here is a direct readout of `getInitialTangent()` (which returns
    `mCe`) assembled over one element -- the quantity
    `refreshInitialElasticOperator()` writes.
    """
    opts = _opts(pre)
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(S._XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', 171, *S._PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, 171)
    for j in range(4):
        ops.fix(j + 1, 1, 1, 1)
    if stage is not None:
        ops.updateMaterialStage('-material', 171, '-stage', stage)
    if do_reset:
        ops.reset()
    return ops.eigen('-fullGenLapack', n_modes)


def test_pre_floor_scales_the_initial_elastic_operator():
    """`refreshInitialElasticOperator()` is real, and this is its ONLY observable.

    The function re-derives `mK/mG/mCe/mCep/mCep_Consistent` at `p = P_atm` with
    the floor applied, because `ManzariDafalias::initialize()` fixes them inside
    the BASE constructor, before the fork's last write lands. It is measurable
    only where the vanilla elastic law actually depends on `p` -- i.e. at
    `mElastFlag == 1` (see the companion test below for why stage 0 cannot move).
    `updateMaterialStage ... 1` then `reset()` puts the material exactly there.

    `pRe = 3*P_atm` is chosen so the predicted factor is an EXACT binary number:
    `sqrt((P_atm + 3*P_atm)/P_atm) = 2`. The stiffness is linear in `G` (and `K`
    is a fixed multiple of `G`), so every eigenvalue must scale by 2 -- mode by
    mode. A floor applied to only one of `K`, `G` would rescale the isotropic and
    deviatoric parts of `C` differently and the ratios would spread.
    """
    off = _eigen_brick(None, 6, 1, True)
    on = _eigen_brick(3.0 * _P_ATM, 6, 1, True)

    for i, (a, b) in enumerate(zip(off, on)):
        ratio = b / a
        assert abs(ratio - 2.0) <= 1.0e-9, (
            f'mode {i + 1}: the initial elastic operator scaled by {ratio!r}, '
            f'not by sqrt((P_atm + 3*P_atm)/P_atm) = 2 -- '
            'refreshInitialElasticOperator() is not reaching mCe, or the floor '
            'is not being applied to K and G together')

    # ... and pin the unfloored deck itself, so the ratio cannot be satisfied by
    # two equally-wrong operators.
    assert abs(off[0] - 67443.766571864) / 67443.766571864 <= 1.0e-9, (
        f'the unfloored reference eigenvalue moved: {off[0]!r}')


def test_pre_floor_is_inert_in_the_elastic_stage_and_that_is_correct():
    """Stage 0 is bit-identical at `pRe = 0` and `pRe = 1e6`. On purpose.

    `mElastFlag` is a static on the base; at `0` all three `GetElasticModuli`
    overloads take the branch `G = G0*P_atm*(2.97-e)^2/(1+e)` -- WITHOUT the
    `sqrt(pn/P_atm)` factor. So the gravity / K0 stage runs a pressure-
    INDEPENDENT elastic law and `pn` is computed and unused. A *confinement*
    floor has nothing to floor there, and making it act would not remove a
    stiffness switch at the stage flip (vanilla's own switch is larger and is
    the point of the staged idiom) -- it would silently restate the G0
    calibration for every `-pRe` deck's gravity step.

    The `-Pmin` control is what makes this a property of the vanilla stage
    rather than a claim about this flag: the clamp lives one line below the
    seam, inside the same function, and is equally inert -- 0.0101 vs 10.0 kPa,
    a 31x move on `G` if the branch were live, changes nothing either.
    """
    def stage0(pre=None, pmin=_PMIN):
        opts = ("-Presidual", 0.0, "-Pmin", pmin, "-honorTolR", 0)
        if pre is not None:
            opts = opts + ("-pRe", float(pre))
        S._build('LadrunoSANISAND', 181, opts)
        S._elastic_leg(181)
        return S._stress()

    base = stage0()
    assert stage0(0.0) == base, '-pRe 0 moved the elastic stage'
    assert stage0(1.0e6) == base, (
        'the ELASTIC stage moved at pRe = 1e6. Either the mElastFlag == 0 '
        'branch of GetElasticModuli has changed and now reads `pn`, or -pRe '
        'has leaked outside the moduli. Both need a decision, not a new '
        'tolerance', base, stage0(1.0e6))
    assert stage0(None, 10.0) == base, (
        'the -Pmin control moved the elastic stage, so the inertness above is '
        'NOT the vanilla pressure-independent elastic law and the reasoning in '
        'this docstring no longer holds')


# ---------------------------------------------------------------------------
#  4 -- a stiffness floor, not a strength floor
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_pre_floor_moves_stiffness_far_more_than_strength():
    """At p0 = 0.01*P_atm: `-pRe` must not buy what `-Presidual` buys.

    `test_presidual_is_the_low_p_defect` measures that `p_r = 1.01` makes this
    path finish ~18 % ABOVE its own bounding surface (`eta/M^b - 1`). `-pRe` is
    claimed to add no strength, so on the same path its `eta/M^b` must stay near
    the `p_r = 0` leg's while the committed stress itself moves. The bar is set
    an order of magnitude below the p_r effect rather than at "zero": the floor
    changes the elastic predictor, so it changes WHICH point on the path the
    same number of steps reaches, and that is not a strength change.

    Wall time **~104 s measured** (103.7 s and 104.3 s on two `--durations`
    runs), hence `slow`: two
    `_drained_triaxial` legs of `_LOWP_STEPS = 1200` steps each at
    `NormDispIncr 1e-8`. (An earlier draft of this docstring said "~40 s (two
    800-step legs)" -- both numbers were wrong; the step count is read from
    `test_ladruno_sanisand._LOWP_STEPS` and the time is measured here. The
    module docstring and the `LEDGER_implementations.md` row carry it too, which
    is what `tests/conftest.py` asks of a `slow` case.)

    It is NOT shortened to a smaller budget: the step count is shared with
    `test_presidual_is_the_low_p_defect`, and the whole value of this case is
    that `-pRe` and `-Presidual` are measured on the SAME path. A cheaper path
    would make the `+18.1 %` control incomparable.
    """
    base, why = S._drained_triaxial('LadrunoSANISAND', 301, _opts())
    assert base is not None, f'pRe = 0 leg did not complete: {why}'
    floored, why = S._drained_triaxial('LadrunoSANISAND', 302, _opts(1.0))
    assert floored is not None, f'pRe = 1.0 leg did not complete: {why}'

    d_strength = abs(floored['err'] - base['err'])
    assert d_strength < 0.05, (
        f"-pRe moved eta/M^b by {d_strength:.3e}, which is the size of a "
        f"STRENGTH change; it is supposed to touch the elastic moduli only "
        f"(pRe=0: eta/M^b-1 = {base['err']:.4f}, "
        f"pRe=1: {floored['err']:.4f})")

    # ... and the same run must show the floor did SOMETHING, or the assertion
    # above is satisfied by a dead parameter.
    assert abs(floored['p'] - base['p']) / base['p'] > 1.0e-6, (
        'pRe = 1.0 kPa left the low-p path bit-identical -- the floor is not '
        'reaching GetElasticModuli on this deck')


# ---------------------------------------------------------------------------
#  5 -- reset() and the wire
# ---------------------------------------------------------------------------

def test_pre_floor_survives_revert_to_start():
    """`revertToStart -> initialize()` zeroes `m_PreElastic` in the BASE.

    `applyLadrunoConstants()` has to win the last write there exactly as it does
    for `m_Presidual`; if it did not, the deck would silently drop its floor
    mid-analysis. The control leg proves the deck can tell the two apart.
    """
    tag = 111
    S._build('LadrunoSANISAND', tag, _opts(1.0))
    S._elastic_leg(tag)
    S._plastic_leg(tag)
    first = S._stress()

    ops.reset()
    ops.setTime(0.0)
    ops.wipeAnalysis()
    S._analysis()
    S._elastic_leg(tag)
    S._plastic_leg(tag)
    second = S._stress()

    assert S._reldiff(first, second) <= 1.0e-12, (
        'the same history gave a different answer after reset() -- the base '
        'initialize() restored m_PreElastic = 0 and nothing took it back',
        first, second)

    unfloored = S._drive('LadrunoSANISAND', 112, _opts())
    assert S._reldiff(first, unfloored) > 1.0e-6, (
        'this deck cannot tell the floor from no floor, so the reset assertion '
        'above would pass with or without applyLadrunoConstants()')


def test_pre_floor_crosses_the_datastore_wire():
    """Slot 34 of the fork `Vector(35)`.

    ADR-86 section 3 is this exact defect in `m_Presidual`: a restored (or MP)
    material that runs a different constitutive law from the process beside it,
    with nothing warning. `_roundtrip` saves mid-path, wipes, rebuilds from the
    datastore and finishes; a dropped slot shows up as a final stress that
    matches the UNFLOORED leg instead of the floored one.
    """
    at_save, after_restore, final = S._roundtrip('LadrunoSANISAND', 121, _opts(1.0))
    assert S._reldiff(at_save, after_restore) <= 1.0e-12, (
        'restore did not reproduce the saved state', at_save, after_restore)

    straight = S._drive('LadrunoSANISAND', 122, _opts(1.0))
    assert S._reldiff(final, straight) <= 1.0e-12, (
        'the restored material finished somewhere else than the uninterrupted '
        'floored leg -- slot 34 (mPreElasticInput) did not cross the wire',
        final, straight)

    unfloored = S._drive('LadrunoSANISAND', 123, _opts())
    assert S._reldiff(straight, unfloored) > 1.0e-6, (
        'the two legs are indistinguishable, so the wire assertion above has '
        'nothing to catch')


# ---------------------------------------------------------------------------
#  6 -- refusals and the echo
# ---------------------------------------------------------------------------

def test_negative_pre_floor_is_refused(capfd):
    """A negative stiffness floor is not "a smaller floor"."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 131, *S._PARAMS, *_opts(-1.0))
    err = capfd.readouterr().err
    assert '-pRe must be >= 0' in err, err


def test_echo_names_the_floor(capfd):
    """ADR 86 section 4.4: every deck-level command echoes what it will run."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('LadrunoSANISAND', 141, *S._PARAMS, *_opts(1.0))
    err = capfd.readouterr().err
    assert 'pRe = 1' in err, err
    assert 'ELASTIC-ONLY floor' in err, err

    ops.nDMaterial('LadrunoSANISAND', 142, *S._PARAMS, *_opts())
    err = capfd.readouterr().err
    assert 'pRe = 0' in err, err
    assert 'no stiffness floor' in err, err


def test_repeated_pre_floor_is_refused(capfd):
    """A second `-pRe` is refused, not silently last-wins.

    Every other flag in this parser last-wins, which is survivable for a
    diagnostic switch. It is not for a constitutive constant: the echo prints
    exactly one value, so `-pRe 1 ... -pRe 5` would report a material the deck
    did not describe and nothing would say so.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 191, *S._PARAMS,
                       *(_OPTS_BASE + ("-pRe", 1.0, "-pRe", 5.0)))
    err = capfd.readouterr().err
    assert 'given more than once' in err, err

    # a synonym counts as the same flag
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 192, *S._PARAMS,
                       *(_OPTS_BASE + ("-pRe", 1.0, "-Pelastic", 5.0)))
    assert 'given more than once' in capfd.readouterr().err


def test_large_pre_floor_warns_but_is_accepted(capfd):
    """`-pRe` above `0.1*P_atm` is a whole-model stiffness statement.

    Unlike `-Presidual`, which the yield surface bounds, the floor multiplies
    `G, K` by `sqrt((p + pRe)/p)` at EVERY Gauss point -- it is not local to the
    ring. WP-106 measured 1 kPa (0.01*P_atm) as already a 7 % move on `G` at the
    live ring, so an order of magnitude above that is where the warning starts.
    A warning, not a refusal: there is no physical bound to appeal to.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('LadrunoSANISAND', 193, *S._PARAMS, *_opts(50.0))
    err = capfd.readouterr().err
    assert 'above 0.1*P_atm' in err, err
    assert 'every Gauss point' in err, err

    # ... and the documented working value does NOT warn
    ops.nDMaterial('LadrunoSANISAND', 194, *S._PARAMS, *_opts(1.0))
    assert 'above 0.1*P_atm' not in capfd.readouterr().err


def test_pre_floor_below_pmin_prints_a_note(capfd):
    """`pRe <= p_min` is the configuration where the floor buys nothing.

    Unfloored, the moduli argument is `max(p, p_min)`; floored it is
    `max(p + pRe, p_min)`. As `p -> 0` -- the free-surface ring, the whole
    subject of ADR 93 -- the floored argument tends to `max(pRe, p_min)`, so a
    `pRe` at or below `p_min` leaves the ring exactly where the clamp already
    had it while still perturbing `G` wherever `p ~ pRe`. Worst of both, and
    invisible from either value alone.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('LadrunoSANISAND', 195, *S._PARAMS,
                   *("-Presidual", 0.0, "-Pmin", 5.0, "-honorTolR", 0,
                     "-pRe", 1.0))
    err = capfd.readouterr().err
    assert 'NOTE: pRe <= p_min' in err, err

    # the documented pairing (pRe = 1.0 kPa, p_min = 0.0101 kPa) is silent
    ops.nDMaterial('LadrunoSANISAND', 196, *S._PARAMS, *_opts(1.0))
    assert 'NOTE: pRe <= p_min' not in capfd.readouterr().err


def test_pelastic_is_accepted_as_a_synonym():
    """`-Pelastic` is the spelling the ADR-93 decision memo (D2) wrote down."""
    a = S._drive('LadrunoSANISAND', 151, _OPTS_BASE + ("-pRe", 1.0))
    b = S._drive('LadrunoSANISAND', 152, _OPTS_BASE + ("-Pelastic", 1.0))
    assert S._reldiff(a, b) <= _EQ_TOL, (
        '-Pelastic did not resolve to the same parameter as -pRe', a, b)
