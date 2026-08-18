"""Every ledger row must name the PR that landed it -- no `(this PR)` survivors (zone_a).

CLAUDE.md requires each build-control ledger row to record its PR. The natural
thing to type while writing the row is `(this PR)`, because the PR does not exist
yet -- and nothing has ever caught it at merge time. So it accumulates: a cleanup
pass in #758 cleared **62** of them, and then #759 and #761 each reintroduced
their own within the hour, by exactly the same reflex. Three occurrences in one
afternoon is a process gap, not carelessness, and the durable fix is a check that
runs in CI rather than a promise to remember.

WHAT TO DO WHEN THIS FAILS
    Replace the `(this PR)` in the named row with a markdown link to the PR that
    is landing it: `[#123](https://github.com/nmorabowen/OpenSees/pull/123)`.
    If the number is not known yet because the PR is not open, open it first --
    the row is only useful once it points somewhere. To resolve an OLD unstamped
    row, `git blame` the line, take the introducing commit's own `(#N)` suffix or
    its `Merge pull request #N`, and failing both, the OLDEST descendant merge on
    the ancestry path (that is the merge which landed it -- the newest is an
    unrelated later merge, a trap worth knowing).

WHY A TEST AND NOT A HOOK
    A hook only protects the machine it is installed on; this repo is worked by
    several agents across several worktrees. zone_a runs in CI on every PR, so a
    test is the one place the rule is enforced for everybody.

The check is deliberately literal-minded: it greps for the marker and reports the
file, line number and row subject. It does not try to validate that the PR number
is CORRECT -- that needs the merge history and belongs in review, not here.
"""
import io
import re
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
LEDGERS = [
    REPO / "Ladruno_implementation" / "LEDGER_vanilla_files.md",
    REPO / "Ladruno_implementation" / "LEDGER_implementations.md",
    REPO / "Ladruno_implementation" / "LEDGER_quirks.md",
]

MARKER = "(this PR"


@pytest.mark.parametrize("ledger", LEDGERS, ids=lambda p: p.name)
def test_no_unstamped_rows(ledger):
    """No row may still say `(this PR)` -- name the PR that landed it."""
    if not ledger.exists():                     # a ledger may be renamed one day
        pytest.skip(f"{ledger.name} not present")
    offenders = []
    with io.open(ledger, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            if MARKER not in line:
                continue
            # the row's subject = its first table cell, else the leading text
            cells = line.split("|")
            subject = (cells[1] if len(cells) > 2 else line).strip()
            offenders.append((n, subject[:88]))
    assert not offenders, (
        f"{ledger.name}: {len(offenders)} row(s) still say '(this PR)'. Replace "
        f"each with a link to the PR that lands it, e.g. "
        f"[#123](https://github.com/nmorabowen/OpenSees/pull/123):\n"
        + "\n".join(f"  line {n}: {s}" for n, s in offenders)
    )


def test_the_guard_can_actually_fail():
    """Falsifier on the guard: prove the marker would be detected if present.

    Without this, a typo in MARKER (or a ledger that silently moved) would make
    every row above pass vacuously -- which is the failure mode this whole file
    exists to prevent, reproduced one level up.
    """
    sample = "| `SRC/foo.cpp` | did a thing | (this PR) |\n"
    assert MARKER in sample
    # and that the real ledgers are actually being read, not skipped into a pass
    assert any(p.exists() for p in LEDGERS), "no ledger was found to check"
    assert any(p.stat().st_size > 1000 for p in LEDGERS if p.exists()), \
        "ledgers found but suspiciously small -- is the path still right?"
