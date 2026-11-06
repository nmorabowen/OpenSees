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
# WP-161: the ledgers are per-WP fragments; the LEDGER_*.md above are stubs then.
FRAGMENT_DIRS = {
    "LEDGER_vanilla_files.md": REPO / "Ladruno_implementation" / "ledger" / "vanilla",
    "LEDGER_implementations.md": REPO / "Ladruno_implementation" / "ledger" / "implementations",
    "LEDGER_quirks.md": REPO / "Ladruno_implementation" / "ledger" / "quirks",
}


def _sources(ledger):
    """The files that hold this ledger's rows: its fragments (WP-161), else the file."""
    d = FRAGMENT_DIRS[ledger.name]
    frags = sorted(p for p in d.glob("*.md") if not p.name.startswith("_")) if d.is_dir() else []
    return frags or ([ledger] if ledger.exists() else [])

MARKER = "(this PR"


@pytest.mark.parametrize("ledger", LEDGERS, ids=lambda p: p.name)
def test_no_unstamped_rows(ledger):
    """No row may still say `(this PR)` -- name the PR that landed it."""
    sources = _sources(ledger)
    if not sources:                             # a ledger may be renamed one day
        pytest.skip(f"{ledger.name} not present")
    offenders = []
    for src in sources:
        with io.open(src, encoding="utf-8") as fh:
            for n, line in enumerate(fh, 1):
                if MARKER not in line:
                    continue
                # the row's subject = its first table cell, else the leading text
                cells = line.split("|")
                subject = (cells[1] if len(cells) > 2 else line).strip()
                where = n if src == ledger else f"{src.name}:{n}"
                offenders.append((where, subject[:88]))
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
    sources = [s for p in LEDGERS for s in _sources(p)]
    assert sources, "no ledger was found to check"
    assert sum(s.stat().st_size for s in sources) > 1000, \
        "ledgers found but suspiciously small -- is the path still right?"
