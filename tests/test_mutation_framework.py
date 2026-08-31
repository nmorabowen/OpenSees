"""ADR-87 D2 — the mutation framework's own guard rails.

These tests do NOT measure a mutation score (that needs two builds; see
Ladruno_scripts/mutation_gate.py). They pin the two properties the framework
must have to be trustworthy at all:

  1. A NORMAL build is unmutated. `ladrunoMutation` returning "none" is what
     stops a sabotaged binary from being shipped, installed or benchmarked by
     mistake -- a mutant computes deliberately wrong physics while looking
     completely ordinary, so "you would notice" is not a safeguard.

  2. The verb EXISTS. The gate driver refuses to score a binary that cannot
     identify itself, because a silently-stale binary would otherwise make the
     gate report the opposite of the truth ("these tests cannot detect deleted
     physics") when in fact nothing was ever deleted.

Zone-A on purpose: this must run on every PR, since the risk it guards is
shipping, not physics.
"""

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def test_mutation_verb_exists():
    """The binary can identify itself as mutant-or-not."""
    assert hasattr(ops, "ladrunoMutation"), (
        "This binary has no `ladrunoMutation` verb. Either it predates ADR-87 D2 "
        "or the .pyd is stale -- see the installed-.pth hijack note in "
        "Ladruno_internal/BUILD_GOTCHAS.md."
    )


def test_default_build_is_not_mutated():
    """A build made the normal way must report no active mutation.

    If this fails, the binary under test was built with -DLADRUNO_MUTATE_FAMILY
    and its physics is deliberately WRONG. Nothing else in the suite means
    anything until that is fixed: every other green test would be green about a
    sabotaged kernel.
    """
    state = ops.ladrunoMutation()
    assert state == "none", (
        f"MUTANT BINARY: ladrunoMutation() == {state!r}. This build computes "
        "deliberately wrong physics (ADR-87 D2) and must never be shipped, "
        "installed or benchmarked. Rebuild without LADRUNO_MUTATE_FAMILY."
    )


def test_mutation_verb_is_quiet_and_ascii():
    """The verb is safe to call from a harness preamble.

    It is meant to be the first line of every gate run, so it must not print,
    must not need a model, and must return plain ASCII (the fork's console
    convention -- Windows consoles are cp1252 and a stray non-ASCII byte turns
    a diagnostic into a UnicodeEncodeError).
    """
    state = ops.ladrunoMutation()
    assert isinstance(state, str) and state
    assert state.isascii(), f"non-ASCII mutation string: {state!r}"
    # Callable repeatedly, with no model defined, with no side effects.
    assert ops.ladrunoMutation() == state
