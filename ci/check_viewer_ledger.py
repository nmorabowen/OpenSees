#!/usr/bin/env python3
"""Viewer-tools ledger gate (WP-121). One rule, checked on a change (a diff), not a tree.

  V1 ledger     A change that brings a new SOURCE file into a tool directory under
                `Ladruno_tools/` (the profiler and monitor viewers) must add a line to
                `Ladruno_implementation/LEDGER_implementations.md` that names that tool
                directory (`Ladruno_tools/<tool>`) -- i.e. it edits that tool's row.
                AGENTS.md has required the row since before the viewers existed; a
                2026-05-31 lesson ("prior PRs ... missed LEDGER_implementations' row;
                check that row explicitly", PR #56) restated it, and it still recurred:
                #485 (profiler_monitor.py) and #487 (all of monitor_viewer/) shipped with
                no ledger change. Evidence: Ladruno_implementation/121_viewer_agent_surface.md.

What counts (Revision 2, after the adversarial review):
  * new file   = status A or C, or a rename whose destination is in a DIFFERENT tool
                 directory than its source (a move into Ladruno_tools/ or between tools).
                 A rename inside one tool directory is not new.
  * source     = .py .pyw .ts .tsx .js .jsx .mjs .cjs .html .css .bat .cmd .sh .ps1.
                 Not source (never flagged): .gitignore, README/docs, JSON/config,
                 fixtures, images. This scope replaces a waiver: editing the tool's row
                 is always possible, so a new source file with no row edit is always a
                 miss of the AGENTS.md rule, never a legitimate exception.
  * row edited = an ADDED ledger line that names `Ladruno_tools/<tool>` (not followed by
                 a name character: `profiler_viewer2` / `profiler_viewer.old` do not count) and is not just a re-added or whitespace-changed copy
                 of a removed line. Another row's edit (e.g. the PR's own WP row) does
                 not count; deleting the ledger does not count.

What it cannot see (stays silent): a change that only MODIFIES viewer files (#55, #484)
-- whether that needs a ledger edit is a judgement, not a pattern. And whether the
edited row actually describes the new file.

The diff is three-dot (`base...head`): the change since the merge base, so a multi-commit
PR whose ledger edit landed in a later commit passes (#58's shape).

Exit codes: 0 clean, 1 finding, 2 the range cannot be evaluated (unknown revision, no
merge base, git error), 3 git itself cannot be run. Exit 2 is deliberately LOUD, unlike the
playbook's "a rule that can't read a case stays silent": that rule is about one unreadable
case in a readable tree. Here an unresolvable range (a shallow clone, a missing base) means
the whole gate cannot run, and passing silently would switch it off unnoticed.

Usage:
    python ci/check_viewer_ledger.py                   # HEAD vs origin/ladruno
    python ci/check_viewer_ledger.py --base <rev> [--head <rev>]
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

SCOPE = "Ladruno_tools/"
LEDGER = "Ladruno_implementation/LEDGER_implementations.md"
SOURCE_SUFFIXES = (".py", ".pyw", ".ts", ".tsx", ".js", ".jsx", ".mjs", ".cjs",
                   ".html", ".css", ".bat", ".cmd", ".sh", ".ps1")

# (status letter, old path or None, new path)
Change = tuple[str, "str | None", str]


class GitUnavailable(RuntimeError):
    """git could not be executed at all (exit 3)."""


def parse_name_status(raw: str) -> list[Change]:
    """Parse `git diff --name-status -z` output into (letter, old, new) triples.

    Records are NUL-separated: `<status>\\0<path>\\0`, or for a rename/copy
    `<Rnnn|Cnnn>\\0<old>\\0<new>\\0`. `old` is None for non-rename/copy records.
    """
    fields = raw.split("\0")
    out: list[Change] = []
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
            out.append((letter, fields[i + 1], fields[i + 2]))
            i += 3
        else:
            if i + 1 >= len(fields):
                break
            out.append((letter, None, fields[i + 1]))
            i += 2
    return out


def parse_unified_diff(raw: str) -> tuple[list[str], list[str]]:
    """(added, removed) content lines of a unified diff; file headers skipped."""
    added: list[str] = []
    removed: list[str] = []
    for line in raw.splitlines():
        if line.startswith("+++") or line.startswith("---"):
            continue
        if line.startswith("+"):
            added.append(line[1:].rstrip("\r"))
        elif line.startswith("-"):
            removed.append(line[1:].rstrip("\r"))
    return added, removed


def tool_dir(path: str | None) -> str | None:
    """`Ladruno_tools/<first component>` for a path in scope, else None."""
    if not path or not path.startswith(SCOPE):
        return None
    first = path[len(SCOPE):].split("/", 1)[0]
    return SCOPE + first if first else None


def is_source(path: str) -> bool:
    return path.lower().endswith(SOURCE_SUFFIXES)


def new_sources(changes: list[Change]) -> list[str]:
    out = []
    for letter, old, new in changes:
        dest = tool_dir(new)
        if dest is None or not is_source(new):
            continue
        if letter in ("A", "C") or (letter == "R" and tool_dir(old) != dest):
            out.append(new)
    return sorted(out)


def _norm(line: str) -> str:
    return " ".join(line.split())


def row_edits(added: list[str], removed: list[str]) -> list[str]:
    """Added ledger lines that are not a re-add / whitespace variant of a removed line."""
    gone = {_norm(r) for r in removed}
    return [a for a in added if _norm(a) and _norm(a) not in gone]


def names_tool(line: str, tool: str) -> bool:
    return re.search(re.escape(tool) + r"(?![\w-]|\.\w)", line) is not None


def find_violations(changes: list[Change], ledger_added: list[str],
                    ledger_removed: list[str]) -> list[str]:
    """V1 findings for one change."""
    files = new_sources(changes)
    edits = row_edits(ledger_added, ledger_removed)
    findings = []
    for path in files:
        tool = tool_dir(path)
        assert tool is not None
        if any(names_tool(line, tool) for line in edits):
            continue
        findings.append(
            f"V1 {path}: new source file, but no added line in {LEDGER} names {tool} -- "
            f"add or update that tool's row in the same PR (AGENTS.md, 'Build-control ledgers')")
    return findings


def _git(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    # stdin=DEVNULL: under pytest on Windows an inherited stdin makes the spawn fail
    # with OSError WinError 6/50 (LEDGER_quirks, "subprocess.run(...) without stdin=").
    try:
        return subprocess.run(
            ["git", "-C", str(root), *args],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
    except OSError as e:
        raise GitUnavailable(f"cannot run git: {e}") from e


def collect(root: Path, base: str, head: str) -> tuple[list[Change], list[str], list[str]]:
    for rev in (base, head):
        r = _git(root, "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}")
        if r.returncode != 0:
            raise LookupError(f"cannot resolve {rev!r} to a commit in {root}")
    rng = f"{base}...{head}"
    changes = parse_name_status(_diff(root, rng, ["--name-status", "-z", "-M"]))
    added, removed = parse_unified_diff(_diff(root, rng, ["-U0"], [LEDGER]))
    return changes, added, removed


def _diff(root: Path, rng: str, opts: list[str], paths: list[str] | None = None) -> str:
    """`git diff <opts> <rng> [-- <paths>]`. A failure (e.g. no merge base) is an error,
    never an empty diff."""
    r = _git(root, "diff", "--no-color", "--no-ext-diff", *opts, rng, "--", *(paths or []))
    if r.returncode != 0:
        raise LookupError(f"git diff {rng} failed (no merge base?): {r.stderr.strip()}")
    return r.stdout


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Viewer-tools ledger gate (WP-121, V1).")
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--base", default="origin/ladruno", help="base revision (default origin/ladruno)")
    ap.add_argument("--head", default="HEAD", help="head revision (default HEAD)")
    args = ap.parse_args(argv)

    try:
        changes, added, removed = collect(args.root, args.base, args.head)
    except GitUnavailable as e:
        print(f"check_viewer_ledger: ERROR {e}", file=sys.stderr)
        return 3
    except LookupError as e:
        print(f"check_viewer_ledger: ERROR {e}", file=sys.stderr)
        return 2

    findings = find_violations(changes, added, removed)
    for f in findings:
        print(f)
    print(f"check_viewer_ledger: {len(findings)} finding(s); {len(new_sources(changes))} new "
          f"source file(s) under {SCOPE} in {args.base}...{args.head}")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
