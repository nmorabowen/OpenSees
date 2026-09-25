"""Self-test for ci/check_viewer_ledger.py (WP-121, rule V1).

One case per shape seen in the viewer PR history (the incidents, the compliant PRs, the
shapes the rule deliberately ignores), every hole the Revision-2 adversarial review found,
and real-git cases for what only git can show (multi-commit PRs, renames, merge bases).
Run:
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
PV = "Ladruno_tools/profiler_viewer"
MV = "Ladruno_tools/monitor_viewer"
PV_ROW = f"| **Stack profiler** ... `{PV}/*` ... | tooling | — | ... | shipped | #58 |"


def A(path):
    return ("A", None, path)


def M(path):
    return ("M", None, path)


def V(changes, added=(), removed=()):
    return cvl.find_violations(changes, list(added), list(removed))


# --- which files count ----------------------------------------------------------------

def test_flags_new_tool_without_ledger():
    # #487: the whole monitor_viewer/ added, ledger untouched. Only source files count.
    found = V([A(f"{MV}/monitor_reader.py"), A(f"{MV}/monitor_page.html"), A(f"{MV}/README.md")])
    assert [f.split(":")[0] for f in found] == [f"V1 {MV}/monitor_page.html", f"V1 {MV}/monitor_reader.py"]


def test_flags_new_file_in_existing_tool_without_ledger():
    # #485 / #53: new sources next to an existing tool, other files modified, no ledger.
    assert len(V([A(f"{PV}/profiler_monitor.py"), M(f"{PV}/README.md"),
                  M("Ladruno_implementation/06_profiler.md")])) == 1


def test_non_source_files_are_out_of_scope():
    # Review finding 4: a .gitignore, docs, config or fixture is never flagged.
    assert V([A(f"{PV}/.gitignore"), A(f"{PV}/README.md"), A(f"{PV}/fixture.json"),
              A(f"{PV}/frontend/public/icons.svg"), A(f"{PV}/frontend/package.json")]) == []


def test_every_source_suffix_counts():
    for suffix in cvl.SOURCE_SUFFIXES:
        assert len(V([A(f"{PV}/x{suffix}")])) == 1, suffix
    assert len(V([A(f"{PV}/Launcher.BAT")])) == 1  # case-insensitive


def test_silent_on_modification_only():
    # #484 / #55: only existing viewer files modified. A judgement, not a pattern.
    assert V([M(f"{PV}/profiler_results.py"), M(f"{PV}/frontend/src/App.tsx")]) == []


def test_ignores_new_files_outside_scope():
    assert V([A("tests/test_new_thing.py"), A("SRC/element/Foo.cpp")]) == []


def test_scope_is_directory_bounded():
    assert V([A("Ladruno_tools_old/x.py"), A("docs/Ladruno_tools/x.py")]) == []


def test_flags_new_sibling_tool_dir():
    # A new viewer next to the existing two (monitor_viewer/ was exactly this in #487).
    assert len(V([A("Ladruno_tools/energy_viewer/app.py")])) == 1


def test_renames_count_only_across_tool_directories():
    # Review finding 3: a move INTO Ladruno_tools/ or into another tool is new to that tool.
    assert V([("R", f"{PV}/old_name.py", f"{PV}/sub/new_name.py")]) == []
    assert len(V([("R", "scripts/x.py", "Ladruno_tools/new_tool/x.py")])) == 1
    assert len(V([("R", f"{PV}/x.py", f"{MV}/x.py")])) == 1
    assert V([("R", f"{PV}/x.py", "scripts/x.py")]) == []   # moved out: nothing new in scope


def test_copy_counts_deletion_does_not():
    assert len(V([("C", f"{PV}/src.py", f"{PV}/copy.py")])) == 1
    assert V([("D", None, f"{PV}/old.py")]) == []


# --- what counts as editing the tool's row ------------------------------------------------

def test_passes_when_an_added_line_names_the_tool():
    # #50 / #51 / #58 shape: the tool's row edited in the same change.
    assert V([A(f"{PV}/profiler_api.py"), M(LEDGER)], added=[PV_ROW + " P7 backend"],
             removed=[PV_ROW]) == []


def test_another_rows_edit_does_not_count():
    # Review finding 2: every WP-era PR adds its OWN row; that must not excuse the tool.
    wp_row = "| **WP-130 — something else** | tooling | — | `ci/foo.py` | draft PR open | — |"
    assert len(V([A(f"{MV}/new_panel.py"), M(LEDGER)], added=[wp_row])) == 1


def test_the_other_tools_row_does_not_count():
    assert len(V([A(f"{MV}/new_panel.py")], added=[PV_ROW + " edited"], removed=[PV_ROW])) == 1


def test_whitespace_only_row_change_does_not_count():
    assert len(V([A(f"{PV}/x.py")], added=[PV_ROW.replace(" ... ", "  ...  ")],
                 removed=[PV_ROW])) == 1


def test_moved_row_does_not_count():
    assert len(V([A(f"{PV}/x.py")], added=[PV_ROW], removed=[PV_ROW])) == 1


def test_deleted_ledger_does_not_count():
    assert len(V([A(f"{PV}/x.py"), ("D", None, LEDGER)], added=[], removed=[PV_ROW])) == 1


def test_tool_name_match_is_bounded():
    # `Ladruno_tools/profiler_viewer2` does not name `Ladruno_tools/profiler_viewer`.
    other = "| x | `Ladruno_tools/profiler_viewer2/*` |"
    assert len(V([A(f"{PV}/x.py")], added=[other])) == 1
    assert V([A(f"{PV}/x.py")], added=[f"| x | `{PV}/` and more |"]) == []
    assert V([A(f"{PV}/x.py")], added=[f"| x | {PV}. |"]) == []          # sentence end
    assert len(V([A(f"{PV}/x.py")], added=[f"| x | {PV}.old/ |"])) == 1


def test_top_level_file_is_its_own_tool():
    assert len(V([A("Ladruno_tools/tool.py")], added=[PV_ROW + " edited"])) == 1
    assert V([A("Ladruno_tools/tool.py")], added=["| x | `Ladruno_tools/tool.py` |"]) == []


# --- parsers --------------------------------------------------------------------------------

def test_parse_name_status_records():
    raw = ("A\0Ladruno_tools/a.py\0"
           "M\0" + LEDGER + "\0"
           "R087\0Ladruno_tools/old.py\0Ladruno_tools/new.py\0"
           "C100\0Ladruno_tools/src.py\0Ladruno_tools/dup.py\0"
           "D\0Ladruno_tools/gone.py\0")
    assert cvl.parse_name_status(raw) == [
        ("A", None, "Ladruno_tools/a.py"),
        ("M", None, LEDGER),
        ("R", "Ladruno_tools/old.py", "Ladruno_tools/new.py"),
        ("C", "Ladruno_tools/src.py", "Ladruno_tools/dup.py"),
        ("D", None, "Ladruno_tools/gone.py"),
    ]
    assert cvl.parse_name_status("") == []


def test_parse_unified_diff_skips_headers_and_crlf():
    raw = ("diff --git a/L b/L\r\n--- a/L\r\n+++ b/L\r\n@@ -3 +3,3 @@\r\n"
           "-old row\r\n+new row\r\n+--- dashes are content\r\n+++ plus-plus content\r\n")
    added, removed = cvl.parse_unified_diff(raw)
    assert removed == ["old row"]
    # A content line that itself starts with '++' is indistinguishable from a header and
    # is dropped -- conservative (it can only make V1 stricter).
    assert added == ["new row", "--- dashes are content"]


# --- end to end on a real git repository ---------------------------------------------------

def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
         "-c", "commit.gpgsign=false", "-c", "core.autocrlf=false", *args],
        stdin=subprocess.DEVNULL, capture_output=True, text=True, check=True)
    return r.stdout.strip()


def _write(repo: Path, rel: str, text: str) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def _commit(repo: Path, msg: str) -> str:
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", msg)
    return _git(repo, "rev-parse", "HEAD")


def _run(repo: Path, base: str, head: str = "HEAD") -> int:
    return cvl.main(["--root", str(repo), "--base", base, "--head", head])


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    if shutil.which("git") is None:
        pytest.skip("git not on PATH: the end-to-end cases need a real repository")
    _git(tmp_path, "init", "-q")
    _write(tmp_path, LEDGER, f"{PV_ROW}\n")
    _write(tmp_path, f"{PV}/profiler_results.py", "".join(f"line {i}\n" for i in range(40)))
    _write(tmp_path, "scripts/helper.py", "".join(f"helper {i}\n" for i in range(40)))
    _commit(tmp_path, "base")
    return tmp_path


def test_end_to_end_multi_commit_pr(repo: Path, capsys):
    # #58's shape: the tool lands in one commit, the row edit in the next, same PR.
    base = _git(repo, "rev-parse", "HEAD")
    _write(repo, f"{PV}/launch.py", "x\n")
    tool = _commit(repo, "tool")
    assert _run(repo, base, tool) == 1
    assert f"V1 {PV}/launch.py" in capsys.readouterr().out
    _write(repo, LEDGER, f"{PV_ROW} + launcher `launch.py`\n")
    fix = _commit(repo, "ledger note")
    assert _run(repo, base, fix) == 0


def test_end_to_end_rename_within_a_tool_is_not_new(repo: Path):
    # Kills `--no-renames`: without rename detection this move reads as A + D.
    base = _git(repo, "rev-parse", "HEAD")
    (repo / PV / "results").mkdir()
    _git(repo, "mv", f"{PV}/profiler_results.py", f"{PV}/results/profiler_results.py")
    _commit(repo, "reorganise")
    assert _run(repo, base) == 0


def test_end_to_end_move_into_scope_is_new(repo: Path):
    # Review finding 3: `git mv scripts/x.py Ladruno_tools/<new tool>/` was silent.
    base = _git(repo, "rev-parse", "HEAD")
    (repo / "Ladruno_tools" / "new_tool").mkdir(parents=True)
    _git(repo, "mv", "scripts/helper.py", "Ladruno_tools/new_tool/helper.py")
    _commit(repo, "adopt helper as a tool")
    assert _run(repo, base) == 1


def test_end_to_end_base_moved_on(repo: Path):
    # Three-dot: a row edit that landed on the BASE branch after the fork point must not
    # excuse the PR.
    fork = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "-b", "pr")
    _write(repo, f"{PV}/new_panel.py", "x\n")
    pr_head = _commit(repo, "tool")
    _git(repo, "checkout", "-q", fork)
    _git(repo, "checkout", "-q", "-b", "base")
    _write(repo, LEDGER, f"{PV_ROW} + someone else's edit\n")
    base_head = _commit(repo, "someone else's row edit")
    assert _run(repo, base_head, pr_head) == 1


def test_unresolvable_base_is_exit_2(repo: Path, capsys):
    assert _run(repo, "no-such-ref") == 2
    assert "cannot resolve" in capsys.readouterr().err


def test_no_merge_base_is_exit_2_not_a_pass(repo: Path, capsys):
    # Kills "ignore the failed diff": unrelated histories have no merge base, git diff
    # A...B fails, and an empty diff must not read as clean.
    main_head = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "--orphan", "unrelated")
    _git(repo, "rm", "-q", "-rf", ".")
    _write(repo, f"{PV}/x.py", "x\n")
    orphan = _commit(repo, "unrelated root")
    assert _run(repo, main_head, orphan) == 2
    assert "no merge base" in capsys.readouterr().err


def test_git_missing_is_exit_3(monkeypatch, capsys):
    # Review finding 6: no git on PATH was a traceback with exit 1 (= "finding").
    def boom(*a, **k):
        raise FileNotFoundError(2, "No such file or directory", "git")
    monkeypatch.setattr(cvl.subprocess, "run", boom)
    assert cvl.main(["--base", "origin/ladruno"]) == 3
    assert "cannot run git" in capsys.readouterr().err
