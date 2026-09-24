"""WP-104 -- `wipe` zeroes LadrunoSANISAND's process-wide IMPL-EX ledger.

THE DEFECT (observed by apeGmsh's live test
`tests/opensees/integration_ladruno/test_ladruno_sanisand_selfweight_implex_live.py`,
2026-09-15, on fork tip 634824e1f): a FRESH `nDMaterial LadrunoSANISAND` in a
NEW model -- after `ops.wipe()`, with a different material tag, before a single
`analyze()` -- reported

    implexRefusals = [9, 0, 0, 9, 0, 9]

i.e. the nine companion refusals an EARLIER model in the same openseespy
process had latched.  The values never grew during the new run; they were
simply inherited.  Slots 0-3 and 5 of `implexRefusals`, all seven slots of
`implexGuards` and `avgImplexError` live in ONE process-wide singleton
(`LadrunoImplexGlobals`, an anonymous-namespace static in
`SRC/material/nD/LadrunoSANISAND.cpp`); only slot 4 (`commitLatched`) is a
member of the material instance.  Nothing in `wipe` reached the singleton, so
the guide's end-of-run discipline -- "read the companion bucket
`implexRefusals[3]` ... it must be 0 on a control-off leg" -- passed or failed
on PROCESS ORDER: which pytest module, which apeGmsh case, ran first.

THE FIX.  `OPS_clearAllNDMaterial()` (the nD-material wipe hook every
interpreter and `OpenSees.exe` go through) now calls
`ladrunoSanisandResetImplexGlobals()`, which zeroes the whole singleton --
the refusal ledger, the post-latch count, the guard census, the error
accumulators AND the commit-round marker (a dangling address after a wipe,
which would otherwise degrade `avgImplexError` to since-process-start for the
rest of the process).  Same rule, same precedent as the ADR-69/72 energy-channel
reset in `Domain::clearAll()`: a wipe destroys every producer and every
consumer of these counters, so a total carried across it is a number about
objects that no longer exist.

WHAT IS DELIBERATELY UNCHANGED, and pinned here so a later "tidy-up" cannot
drift it: `ops.reset()` / `revertToStart()` do NOT zero the ledger.  They
rewind the SAME model, and a leg's running totals across reverts are exactly
what the `LEDGER_quirks` rule "read `implexRefusals` as DELTAS" is for.  The
10-per-process `opserr` throttles at the refusal sites are also untouched (a
warning budget is a log-volume contract for the process, not a census).

DECK.  The WP-99 starved-cap `LadrunoQuad` deck (`-maxSubsteps 2`), imported
from its own test module so the two files cannot drift: it latches a genuine
commit-time companion refusal within a handful of deviatoric steps, which is
the cheapest way to make every process-wide slot nonzero.  The classic-Tcl
twin of the first test is `tests/tcl/wp104_implex_refusals_wipe.tcl`
(runner: `test_wp104_implex_refusals_wipe_tcl.py`) -- `wipe` in
`SRC/tcl/commands.cpp` is a separate code path from the openseespy one and
both must reach the hook.

MEASURED WALL TIME: ~3 s for the file.
"""
import pytest

from _testbed import ops

import test_ladrunoQuad_sanisand_implex_commit_refusal as cr

pytestmark = [pytest.mark.zone_a]

_ZERO6 = [0.0] * 6
_ZERO7 = [0.0] * 7


def _refusals():
    return cr._refusals()


def _guards():
    g = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    assert len(g) == 7, ('implexGuards must be the 7-slot census (ADR-92 P2-9)', g)
    return g


def _avg_error():
    return float(ops.eleResponse(1, 'material', 1, 'avgImplexError')[0])


def _run_latching_model(tag):
    """Build the starved WP-99 quad, confine, flip, drive to the commit-time
    latch, then push ONE more step so the post-latch bucket (slot 5) moves
    too.  Returns the ledger at the end."""
    cr._build_quad(tag, ('-implex', '-maxSubsteps', cr._CAP_STARVED))
    cr._confine_and_flip(tag)
    cr._add_deviatoric_pattern_2d()
    cr._drive_until_latch()
    rc = ops.analyze(1)          # refused while latched -> slot 5 increments
    assert rc != 0, 'the step after a latched commit converged (WP-99 contract)'
    end = _refusals()
    assert end[3] > 0 and end[4] == 1.0 and end[5] > 0, (
        'the latching model did not populate every slot the test relies on '
        '(companion, commitLatched, latched) -- the DECK is broken, not the '
        'contract under test', end)
    assert end[0] == pytest.approx(end[1] + end[2] + end[3]), end
    return end


def _build_fresh_model(tag):
    """A different tag, an ADEQUATE cap: a model that never refuses."""
    cr._build_quad(tag, ('-implex', '-maxSubsteps', 20000))


# ===========================================================================
#  (1) the reported defect: a fresh model after wipe must read all-zero
# ===========================================================================
def test_fresh_model_after_wipe_starts_at_zero():
    """Model A latches (every process-wide slot nonzero).  `wipe`.  Model B
    -- new tag, adequate cap -- must report `[0,0,0,0,0,0]` BEFORE any
    analysis, and still all-zero after its own confinement leg, on every
    process-wide response: `implexRefusals`, `implexGuards`, `avgImplexError`.

    Kills: the tip-634824e1f binary (reads model A's totals -- the reported
    `[9,0,0,9,0,9]` shape); a fix that zeroes only the refusal ledger but not
    the guard census or the error accumulators; a fix hooked somewhere the
    openseespy `wipe` does not reach.
    """
    end_a = _run_latching_model(9940)
    assert end_a != _ZERO6

    ops.wipe()
    _build_fresh_model(9941)

    r0 = _refusals()
    assert r0 == _ZERO6, (
        'a FRESH LadrunoSANISAND in a NEW model, after ops.wipe(), inherited '
        'the previous model\'s refusal ledger. The counters are process-wide '
        'by design, but wipe destroys every producer and consumer, so they '
        'must read zero here -- otherwise an end-of-run '
        '`implexRefusals[3] == 0` assertion depends on which model ran first '
        f'in the process. previous model ended at {end_a}', r0)
    assert _guards() == _ZERO7, (
        'implexGuards survived ops.wipe() -- same singleton, same rule', _guards())
    assert _avg_error() == 0.0, (
        'avgImplexError survived ops.wipe()', _avg_error())

    cr._confine_and_flip(9941)
    assert _refusals() == _ZERO6, (
        'the fresh model\'s own confinement leg refused something -- an '
        'adequate-cap deck must not; check the deck before the fix',
        _refusals())


# ===========================================================================
#  (2) two consecutive models in one process: order-independent census
# ===========================================================================
def test_two_consecutive_models_both_start_at_zero_and_agree():
    """The property the guide's end-of-run check actually needs: the SAME
    model, run twice in one process with a `wipe` between, starts from
    `[0]*6` both times and ENDS with the same ledger both times.  Before
    WP-104 the second run started where the first ended and its end ledger
    was the double.
    """
    ends = []
    for tag in (9942, 9943):
        ops.wipe()
        cr._build_quad(tag, ('-implex', '-maxSubsteps', cr._CAP_STARVED))
        assert _refusals() == _ZERO6, (
            f'model with tag {tag} did not start from a zero ledger; '
            f'previous ends: {ends}', _refusals())
        assert _guards() == _ZERO7
        ends.append(_run_latching_model(tag))
    assert ends[0] == ends[1], (
        'two identical models in one process ended with DIFFERENT ledgers -- '
        'the second one inherited the first one\'s totals', ends)


# ===========================================================================
#  (3) the boundary: reset() is NOT a wipe and keeps the ledger
# ===========================================================================
def test_reset_keeps_the_ledger_within_a_model():
    """`ops.reset()` rewinds the same model; the ledger is a running total
    over that model's legs and must survive it (`LEDGER_quirks`: read it as
    deltas).  The per-instance latch (slot 4) IS cleared by revertToStart --
    that is WP-99's own contract, re-asserted so the two behaviours are read
    side by side.

    Kills: a fix that hooks the reset into `revertToStart()` / `initialize()`
    instead of the wipe hook -- which would silently zero the census between
    the arms of a multi-leg campaign.
    """
    end = _run_latching_model(9944)
    ops.reset()
    after = _refusals()
    assert after[:4] == end[:4] and after[5] == end[5], (
        'ops.reset() zeroed the process-wide ledger. Only wipe may do that; '
        'reset rewinds the SAME model and a campaign reads its legs as '
        'deltas over a total that must persist', end, after)
    assert after[4] == 0.0, (
        'the per-instance commit latch survived revertToStart() -- WP-99 '
        'says it is cleared there', after)
