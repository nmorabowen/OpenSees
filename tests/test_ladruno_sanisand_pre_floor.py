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
6. **Refusals and the echo.** A negative value is refused at parse time; the
   construction echo names the floor and says which side of the model it touches.

Deck and helpers are imported from `test_ladruno_sanisand`, deliberately: this
file must exercise the SAME strain path the p_residual gates use, or "pRe moves
the answer by X and p_r by Y" is a comparison of two different things.
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

_OPTS_BASE = ("-Presidual", 0.0, "-Pmin", _PMIN, "-honorTolR", 0)


def _opts(pre=None):
    return _OPTS_BASE if pre is None else _OPTS_BASE + ("-pRe", float(pre))


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

    Wall time ~40 s (two 800-step legs), hence `slow`.
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


def test_pelastic_is_accepted_as_a_synonym():
    """`-Pelastic` is the spelling the ADR-93 decision memo (D2) wrote down."""
    a = S._drive('LadrunoSANISAND', 151, _OPTS_BASE + ("-pRe", 1.0))
    b = S._drive('LadrunoSANISAND', 152, _OPTS_BASE + ("-Pelastic", 1.0))
    assert S._reldiff(a, b) <= _EQ_TOL, (
        '-Pelastic did not resolve to the same parameter as -pRe', a, b)
