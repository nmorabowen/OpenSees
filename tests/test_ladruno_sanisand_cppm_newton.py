"""WP-130 (TIMs F18(c) + F18(d)): BackwardEuler_CPPM under a global Newton,
and the ModifiedEuler -> CPPM per-point fallback.

What is pinned here:

1. BYTE IDENTITY at the defaults.  Six IntScheme-2 decks
   (`wp130_sanisand_byteid.py`: vanilla ManzariDafalias, LadrunoSANISAND 3D
   incl. a huge-increment leg that drives the CPPM's halving ladder into its
   explicit fallback, the TIMs campaign set in plane strain incl. a reversal
   and a free-DOF global-Newton deck whose iteration counts are pinned too)
   reproduce, bit for bit, what the PRE-WP-130 binary produced.  WP-127's
   IntScheme-1 decks are pinned by its own test, unchanged.
2. THE PARSER refuses every flag combination that would be inert or is not
   qualified (-implex).
3. THE CENSUS: `substepStats` has 28 columns; at the defaults the vanilla
   silent explicit fallback is now COUNTED (F12 5.2/5.3 said it was invisible).
4. F18(c): with `-cppmOnFail refuse` a trial iterate the CPPM cannot return
   reaches analyze() as a failure (rc < 0) in bounded time, the refusal is
   counted, and nothing is integrated explicitly.
5. F18(d): the one-element fallback test -- a ModifiedEuler cap hit that
   refuses the step without `-meFallback` is carried by the CPPM with it, and
   agrees with an uncapped integration.

MEASURED WALL TIME: see the WP-130 PR (the byte-identity decks dominate).
"""
import json
import time

import pytest

from _testbed import ops
import test_ladruno_sanisand as sani
import wp127_sanisand_byteid as b127
import wp130_sanisand_byteid as b130

pytestmark = [pytest.mark.zone_a]

_NSTATS = 28
(CPPM_CALLS, CPPM_NFAIL, CPPM_HALV, CPPM_EXPL, CPPM_LOWP, CPPM_REF,
 ME_FB, ME_FB_OK, LAST_CPPM_REF, GUESS_TRIES, GUESS_OK) = range(17, 28)
CAP = 9
SUB = 2


def _stats(ele=1, gp=1):
    s = list(ops.eleResponse(ele, "material", gp, "substepStats"))
    assert len(s) == _NSTATS, s
    return s


# ---------------------------------------------------------------------------
#  1. byte identity against the pre-WP-130 binary
# ---------------------------------------------------------------------------

def test_scheme2_defaults_are_byte_identical():
    with open(b130.BASELINE) as fh:
        ref = json.load(fh)["decks"]
    got = b130.run_all()
    assert sorted(got) == sorted(ref)
    for name in ref:
        assert len(got[name]) == len(ref[name]), name
        for k, (a, b) in enumerate(zip(got[name], ref[name])):
            assert a == b, f"deck {name} row {k}: first differing entry " \
                f"{next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)}"


# ---------------------------------------------------------------------------
#  2. the parser
# ---------------------------------------------------------------------------

def _mat(*opts):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    try:
        ops.nDMaterial("LadrunoSANISAND", 1, *sani._PARAMS, *opts)
    except Exception:
        return False
    try:
        return 1 in ops.getNDMaterialTags() if hasattr(ops, "getNDMaterialTags") else True
    except Exception:
        return True


@pytest.mark.parametrize("opts", [
    (1, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "refuse"),               # scheme 1: inert
    (2, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "maybe"),                # bad token
    (2, 2, 1, 1e-7, 1e-7, "-cppmHalvings", 10),                   # out of range
    (2, 2, 1, 1e-7, 1e-7, "-cppmHalvings", -1),
    (1, 2, 1, 1e-7, 1e-7, "-cppmLineSearch", "on"),               # scheme 1, no fallback
    (1, 2, 1, 1e-7, 1e-7, "-meFallback", "cppm"),                 # no -maxSubsteps
    (2, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-meFallback", "cppm"),   # scheme 2
    (2, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-implex", "-cppmOnFail", "refuse"),
    (1, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-implex", "-meFallback", "cppm"),
])
def test_parser_refuses_inert_or_unqualified(opts):
    with pytest.raises(Exception):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("LadrunoSANISAND", 1, *sani._PARAMS, *opts)


@pytest.mark.parametrize("opts", [
    (2, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "refuse", "-cppmHalvings", 0),
    (2, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "explicit", "-cppmLineSearch", "on"),
    (1, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-meFallback", "cppm",
     "-cppmHalvings", 3, "-cppmLineSearch", "on"),
])
def test_parser_accepts(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoSANISAND", 1, *sani._PARAMS, *opts)


# ---------------------------------------------------------------------------
#  3. the census makes vanilla's silent explicit fallback visible
# ---------------------------------------------------------------------------

def _free_push(extra, lateral=50.0, push=20.0):
    """The byte-id free quad (100 kPa, loaded edges) under IntScheme 2 +
    `extra`, one push step of 0.1*push kPa.  Returns (rc, wall s, census)."""
    b130.build_free_quad(b130._CAMPAIGN_S2 + tuple(extra), lateral=lateral)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(sani._XY):
        if y == 1.:
            ops.load(j + 1, 0.0, -push)
    ops.integrator("LoadControl", 0.1)
    t = time.time()
    rc = ops.analyze(1)
    return rc, time.time() - t, _stats()


def test_default_counts_the_silent_explicit_fallback():
    """F12 5.3: a CPPM failure used to be invisible in every channel.  At the
    DEFAULTS (vanilla control flow, byte-identical) the census now shows it.
    Measured on build 428328adc: rc -3 after 31 Newton iterations in ~5.6 s,
    224 local-Newton failures, 440 half-increments, 4 silent explicit
    fallbacks, 0 refusals."""
    rc, wall, s = _free_push(())
    assert rc < 0
    assert s[CPPM_CALLS] > 0
    assert s[CPPM_NFAIL] > 0 and s[CPPM_HALV] > 0
    assert s[CPPM_EXPL] > 0, s          # the silent fallback, now counted
    assert s[CPPM_REF] == 0 and s[LAST_CPPM_REF] == 0


# ---------------------------------------------------------------------------
#  4. F18(c): refuse at once, and the refusal reaches analyze()
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("extra", [
    ("-cppmOnFail", "refuse", "-cppmHalvings", 0),
    ("-cppmOnFail", "refuse", "-cppmHalvings", 0, "-cppmStart", "explicit",
     "-cppmLineSearch", "on"),
])
def test_refusal_reaches_analyze_fast(extra):
    """Same deck and step as above: the refusal propagates (LADRUNO_MATERIAL_
    REFUSED -> quad -> Domain::update -> analyze rc < 0) on a trial iterate
    the CPPM cannot return, with nothing integrated explicitly and no halving.
    Measured 8-22 ms against ~5.6 s at the defaults."""
    rc_def, wall_def, _ = _free_push(())
    rc, wall, s = _free_push(extra)
    assert rc < 0
    assert s[CPPM_REF] == 1 and s[LAST_CPPM_REF] == 1, s
    assert s[CPPM_EXPL] == 0 and s[CPPM_HALV] == 0, s
    assert s[1] == 0, s                 # ModifiedEuler never ran
    assert wall < 1.0 and wall < 0.2 * wall_def, (wall, wall_def)


def test_refused_update_is_not_sticky():
    """The refusal flag is per update: after the failed step a much smaller
    load step from the same committed state integrates (nothing latched)."""
    rc, _, s0 = _free_push(("-cppmOnFail", "refuse", "-cppmHalvings", 0))
    assert rc < 0 and s0[CPPM_REF] == 1
    ops.integrator("LoadControl", 1.0e-4)
    assert ops.analyze(1) == 0
    s = _stats()
    # no new refusal; LAST_CPPM_REF keeps describing the last PLASTIC update
    # (a small step here can be elastic, which by design does not reset it)
    assert s[CPPM_REF] == s0[CPPM_REF], s


# ---------------------------------------------------------------------------
#  5. F18(d): the per-point ModifiedEuler -> CPPM fallback, one element
# ---------------------------------------------------------------------------

def _fallback_deck(extra, econf=3.0e-4, de=2.0e-3, steps=10):
    """Zero-free-DOF plane-strain quad, TIMs campaign set, IntScheme 1,
    confinement to p ~ 131 kPa (econf 3e-4), then `steps` steps of `de`
    axial compression / lateral extension. Returns (rcs, stress, census)."""
    opts = (1, 2, 1, 1.0e-7, 1.0e-7, "-Presidual", 0.0, "-Pmin", 0.0101,
            "-flipAlphaIn", "init") + tuple(extra)
    incs = b127._iso_dev(steps, de, 1.0)
    b127._build_ps(b127._CAMPAIGN, opts, 10, econf, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    rcs = []
    for _ in incs:
        rcs.append(ops.analyze(1))
        if rcs[-1] != 0:
            break
    return rcs, list(ops.eleResponse(1, "material", 1, "stress")), _stats()


def test_me_fallback_carries_a_capped_point():
    """THE one-element test (F10b(b)).  Measured on build 428328adc: uncapped
    ModifiedEuler takes up to 123 substeps in one update on this leg;
    `-maxSubsteps 20` refuses step 1; with `-meFallback cppm` every capped
    update is returned by the CPPM, all 10 steps commit, and the stress agrees
    with the uncapped integration to 1.3 % (the two integrators differ by
    their own discretisation error at this increment, F12 section 2)."""
    rc_ref, sig_ref, s_ref = _fallback_deck(())
    assert rc_ref == [0] * 10 and s_ref[CAP] == 0
    rc_cap, _, s_cap = _fallback_deck(("-maxSubsteps", 20))
    assert rc_cap[-1] < 0 and s_cap[CAP] >= 1 and s_cap[ME_FB] == 0
    rc_fb, sig_fb, s_fb = _fallback_deck(("-maxSubsteps", 20, "-meFallback", "cppm"))
    assert rc_fb == [0] * 10, rc_fb
    assert s_fb[ME_FB] >= 10 and s_fb[ME_FB_OK] == s_fb[ME_FB], s_fb
    assert s_fb[CAP] == s_fb[ME_FB]          # every cap hit was handed over
    assert s_fb[CPPM_REF] == 0 and s_fb[CPPM_EXPL] == 0 and s_fb[LAST_CPPM_REF] == 0
    for a, b in zip(sig_fb[:2], sig_ref[:2]):
        assert abs(a - b) <= 0.03 * abs(b), (sig_fb, sig_ref)


def test_me_fallback_refuses_when_the_cppm_fails_too():
    """The fallback CPPM never integrates explicitly (that would re-enter the
    ModifiedEuler that just failed): where it cannot return the increment the
    update is REFUSED.  Measured: at p ~ 44 kPa and 5e-3 steps the first
    plastic step fails both ways."""
    rcs, _, s = _fallback_deck(("-maxSubsteps", 20, "-meFallback", "cppm", "-cppmHalvings", 0),
                               econf=1.0e-4, de=5.0e-3)
    assert rcs[-1] < 0
    assert s[ME_FB] >= 1 and s[ME_FB_OK] < s[ME_FB] and s[CPPM_REF] >= 1, s
    assert s[CPPM_EXPL] == 0 and s[CPPM_LOWP] == 0, s
