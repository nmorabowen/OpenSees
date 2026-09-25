"""Self-test for ci/check_viewer_ledger.py (WP-121, rule V1).

One case per shape seen in the viewer PR history (the incident, the fix, and every
passing shape), plus a real-git case for the multi-commit PR (#58). Run:
    pytest -q ci/test_check_viewer_ledger.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_viewer_ledger as cvl  # noqa: E402

LEDGER = cvl.LEDGER


# --- find_violations: the rule itself ---------------------------------------------------

def test_flags_new_tool_without_ledger():
    # #487: the whole monitor_viewer/ added, ledger untouched.
    changes = [("A", "Ladruno_tools/monitor_viewer/monitor_reader.py"),
               ("A", "Ladruno_tools/monitor_viewer/README.md")]
    found = cvl.find_violations(changes)
    assert len(found) == 2
    assert all(f.startswith("V1 Ladruno_tools/monitor_viewer/") for f in found)


def test_flags_new_file_in_existing_tool_without_ledger():
    # #485 / #53: new files next to an existing tool, other files modified, no ledger.
    changes = [("A", "Ladruno_tools/profiler_viewer/profiler_monitor.py"),
               ("M", "Ladruno_tools/profiler_viewer/README.md"),
               ("M", "Ladruno_implementation/06_profiler.md")]
    assert len(cvl.find_violations(changes)) == 1


def test_passes_new_file_with_ledger_modified():
    # #50 / #51 shape: new files + the ledger row edited in the same change.
    changes = [("A", "Ladruno_tools/profiler_viewer/profiler_api.py"),
               ("M", LEDGER)]
    assert cvl.find_violations(changes) == []


def test_passes_ledger_added():
    changes = [("A", "Ladruno_tools/x/y.py"), ("A", LEDGER)]
    assert cvl.find_violations(changes) == []


def test_silent_on_modification_only():
    # #484 / #55: only existing viewer files modified. A judgement, not a pattern.
    changes = [("M", "Ladruno_tools/profiler_viewer/profiler_results.py"),
               ("M", "Ladruno_tools/profiler_viewer/frontend/src/App.tsx")]
    assert cvl.find_violations(changes) == []


def test_ignores_new_files_outside_scope():
    changes = [("A", "tests/test_new_thing.py"), ("A", "SRC/element/Foo.cpp")]
    assert cvl.find_violations(changes) == []


def test_scope_is_directory_bounded():
    changes = [("A", "Ladruno_tools_old/x.py"), ("A", "docs/Ladruno_tools/x.md")]
    assert cvl.find_violations(changes) == []


def test_flags_new_sibling_tool_dir():
    # A new viewer next to the existing two (monitor_viewer/ was exactly this in #487).
    assert len(cvl.find_violations([("A", "Ladruno_tools/energy_viewer/app.py")])) == 1


def test_rename_is_not_an_addition_but_copy_is():
    assert cvl.find_violations([("R", "Ladruno_tools/profiler_viewer/new_name.py")]) == []
    assert len(cvl.find_violations([("C", "Ladruno_tools/profiler_viewer/copy.py")])) == 1


def test_deletions_are_ignored():
    assert cvl.find_violations([("D", "Ladruno_tools/profiler_viewer/old.py")]) == []


# --- parse_name_status: git -z output ------------------------------------------------

def test_parse_plain_and_rename_records():
    raw = ("A\0Ladruno_tools/a.py\0"
           "M\0" + LEDGER + "\0"
           "R087\0Ladruno_tools/old.py\0Ladruno_tools/new.py\0"
           "C100\0Ladruno_tools/src.py\0Ladruno_tools/dup.py\0"
           "D\0Ladruno_tools/gone.py\0")
    assert cvl.parse_name_status(raw) == [
        ("A", "Ladruno_tools/a.py"),
        ("M", LEDGER),
        ("R", "Ladruno_tools/new.py"),
        ("C", "Ladruno_tools/dup.py"),
        ("D", "Ladruno_tools/gone.py"),
    ]


def test_parse_empty():
    assert cvl.parse_name_status("") == []


# --- end to end on a real git repository ------------------------------------------------

def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
         "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", *args],
        stdin=subprocess.DEVNULL, capture_output=True, text=True, check=True)
    return r.stdout.strip()


def _commit(repo: Path, rel: str, text: str, msg: str) -> str:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    _git(repo, "add", rel)
    _git(repo, "commit", "-q", "-m", msg)
    return _git(repo, "rev-parse", "HEAD")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    if shutil.which("git") is None:
        pytest.skip("git not on PATH: the end-to-end cases need a real repository")
    _git(tmp_path, "init", "-q")
    _commit(tmp_path, LEDGER, "| row |\n", "base")
    return tmp_path


def test_end_to_end_multi_commit_pr(repo: Path, capsys):
    # #58's shape: the tool lands in one commit, the ledger note in the next, same PR.
    base = _git(repo, "rev-parse", "HEAD")
    tool = _commit(repo, "Ladruno_tools/profiler_viewer/launch.py", "x\n", "tool")
    assert cvl.main(["--root", str(repo), "--base", base, "--head", tool]) == 1
    assert "V1 Ladruno_tools/profiler_viewer/launch.py" in capsys.readouterr().out
    fix = _commit(repo, LEDGER, "| row | launcher |\n", "ledger note")
    assert cvl.main(["--root", str(repo), "--base", base, "--head", fix]) == 0


def test_end_to_end_base_moved_on(repo: Path):
    # Three-dot: a ledger edit that landed on the BASE branch after the fork point
    # must not excuse the PR (and a base-only file must not be charged to it).
    fork = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "-b", "pr")
    pr_head = _commit(repo, "Ladruno_tools/new_tool/a.py", "x\n", "tool")
    _git(repo, "checkout", "-q", fork)
    _git(repo, "checkout", "-q", "-b", "base")
    base_head = _commit(repo, LEDGER, "| other row |\n", "someone else's row")
    assert cvl.main(["--root", str(repo), "--base", base_head, "--head", pr_head]) == 1


def test_unresolvable_base_is_an_error_not_a_pass(repo: Path, capsys):
    assert cvl.main(["--root", str(repo), "--base", "no-such-ref"]) == 2
    assert "cannot resolve" in capsys.readouterr().err
