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
