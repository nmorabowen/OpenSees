#!/usr/bin/env python3
"""Viewer-tools ledger gate (WP-121). One rule, checked on a change (a PR diff), not a tree.

  V1 ledger     A change that ADDS a file under `Ladruno_tools/` (the profiler and
                monitor viewers) must also modify
                `Ladruno_implementation/LEDGER_implementations.md`.
                AGENTS.md has required this since before the viewers existed, and a
                2026-05-31 lesson ("prior PRs ... missed LEDGER_implementations' row;
                check that row explicitly", PR #56) restated it. It still recurred:
                #485 (profiler_monitor.py) and #487 (the whole monitor_viewer/) added
                tools with no ledger change. Evidence and acceptance runs:
                Ladruno_implementation/121_viewer_agent_surface.md.

What it cannot see (stays silent, never guesses): a change that only MODIFIES viewer
files (#55's overlay/leak badge, #484's fix) -- whether that needs a ledger edit is a
judgement, not a pattern. Renames are not additions; copies are.

The diff is three-dot (`base...head`): the change since the merge base, so a
multi-commit PR whose ledger edit landed in a later commit passes (#58's shape).

Usage:
    python ci/check_viewer_ledger.py                   # HEAD vs origin/ladruno
    python ci/check_viewer_ledger.py --base <rev> [--head <rev>]
Exit: 0 clean, 1 on a finding, 2 when git cannot resolve the range.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCOPE = "Ladruno_tools/"
LEDGER = "Ladruno_implementation/LEDGER_implementations.md"


def parse_name_status(raw: str) -> list[tuple[str, str]]:
    """Parse `git diff --name-status -z` output into (status letter, new path) pairs.

    Records are NUL-separated: `<status>\\0<path>\\0`, or for a rename/copy
    `<Rnnn|Cnnn>\\0<old>\\0<new>\\0`. The returned path is the post-change path.
    """
    fields = raw.split("\0")
    out: list[tuple[str, str]] = []
    i = 0
    while i < len(fields):
        status = fields[i]
        if not status:
            i += 1
            continue
        letter = status[0]
        if letter in ("R", "C"):
            if i + 2 >= len(fields):
                break
            out.append((letter, fields[i + 2]))
            i += 3
        else:
            if i + 1 >= len(fields):
                break
            out.append((letter, fields[i + 1]))
            i += 2
    return out


def find_violations(changes: list[tuple[str, str]]) -> list[str]:
    """V1: every added/copied path under SCOPE needs LEDGER in the same change."""
    added = sorted(p for s, p in changes if s in ("A", "C") and p.startswith(SCOPE))
    if not added:
        return []
    if any(p == LEDGER for _, p in changes):
        return []
    return [
        f"V1 {p}: new file under {SCOPE} but {LEDGER} is unchanged -- "
        f"add or update the tool's row in the same PR (AGENTS.md, 'Build-control ledgers')"
        for p in added
    ]


def _git(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    # stdin=DEVNULL: under pytest on Windows an inherited stdin makes the spawn fail
    # with OSError WinError 6/50 (LEDGER_quirks, "subprocess.run(...) without stdin=").
    return subprocess.run(
        ["git", "-C", str(root), *args],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def changes_between(root: Path, base: str, head: str) -> list[tuple[str, str]]:
    for rev in (base, head):
        r = _git(root, "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}")
        if r.returncode != 0:
            raise LookupError(f"cannot resolve {rev!r} to a commit in {root}")
    r = _git(root, "diff", "--name-status", "-z", "-M", "--no-color", f"{base}...{head}")
    if r.returncode != 0:
        raise LookupError(f"git diff {base}...{head} failed: {r.stderr.strip()}")
    return parse_name_status(r.stdout)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Viewer-tools ledger gate (WP-121, V1).")
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--base", default="origin/ladruno", help="base revision (default origin/ladruno)")
    ap.add_argument("--head", default="HEAD", help="head revision (default HEAD)")
    args = ap.parse_args(argv)

    try:
        changes = changes_between(args.root, args.base, args.head)
    except LookupError as e:
        print(f"check_viewer_ledger: ERROR {e}", file=sys.stderr)
        return 2

    findings = find_violations(changes)
    for f in findings:
        print(f)
    n_added = sum(1 for s, p in changes if s in ("A", "C") and p.startswith(SCOPE))
    print(f"check_viewer_ledger: {len(findings)} finding(s); {n_added} file(s) added under "
          f"{SCOPE} in {args.base}...{args.head}")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
