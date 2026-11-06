"""Self-test for quirk-lint rule `floating-runner` (WP-225) in ci/check_quirk_patterns.py.

One file per rule (WP-162). Run: pytest -q ci/test_quirk_floating_runner.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def _run(tmp_path, body, name=".github/workflows/ladruno.yml"):
    root = kit.tree(tmp_path, {name: body})
    return cq.check_floating_runner(root, kit.rel(root))


@pytest.mark.parametrize("line", [
    "    runs-on: ubuntu-latest",                    # the #1063 incident (Zone A)
    "    runs-on: [windows-latest]",
    "        os: [ubuntu-24.04, macos-latest]",       # a matrix entry
])
def test_flags_a_floating_image(tmp_path, line):
    out = _run(tmp_path, "jobs:\n  a:\n" + line + "\n")
    assert len(out) == 1 and out[0].startswith("floating-runner .github/workflows/ladruno.yml:3"), out


def test_passes_a_pinned_image_and_a_comment(tmp_path):
    body = "jobs:\n  a:\n    runs-on: ubuntu-24.04   # was ubuntu-latest until WP-225\n"
    assert _run(tmp_path, body) == []


def test_ignores_non_fork_workflows(tmp_path):
    # vanilla upstream workflows (build_cmake.yml) are out of scope (vanilla-footprint rule)
    assert _run(tmp_path, "    runs-on: ubuntu-latest\n", name=".github/workflows/build_cmake.yml") == []


def test_the_fork_workflows_are_pinned():
    root = Path(__file__).resolve().parent.parent
    assert cq.check_floating_runner(root, kit.rel(root)) == []
