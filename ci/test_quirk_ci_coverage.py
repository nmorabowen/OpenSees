"""Self-test for quirk-lint rule `ci-coverage` (alias L8, WP-143) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_ci_coverage.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def _run(tmp_path, body, name="tests/test_x.py"):
    root = kit.tree(tmp_path, {name: body})
    return cq.check_ci_coverage(root, kit.rel(root))


_ZA = "import os, sys, platform\nimport pytest\npytestmark = [pytest.mark.zone_a]\n"


@pytest.mark.parametrize("gate", [
    'pytestmark.append(pytest.mark.skipif(sys.platform != "win32", reason="mkl"))',   # the WP-136 incident
    'if os.name == "nt":\n    pass',                                                   # adr74's cleanup branch
    'if platform.system() == "Windows":\n    pass',
    'if sys.platform.startswith("win"):\n    pass',
    'TOL = 0.0\nif "win32" != sys.platform:\n    TOL = 1e-9',                            # reversed operands
])
def test_flags_an_undeclared_platform_branch(tmp_path, gate):
    out = _run(tmp_path, _ZA + gate + "\n")
    assert len(out) == 1 and out[0].startswith("L8 ") and "ci-coverage" in out[0], out


@pytest.mark.parametrize("kind", ["local-only", "partial", "portable", "nightly-windows", "pr-windows"])
def test_passes_a_declared_branch(tmp_path, kind):
    body = _ZA + f"# ci-coverage: {kind} because the MKL leg needs Pardiso\n" + \
        'pytestmark.append(pytest.mark.skipif(sys.platform != "win32", reason="mkl"))\n'
    assert _run(tmp_path, body) == []


def test_ignores_a_ternary_value_selection(tmp_path):
    body = _ZA + 'EXE = "OpenSees.exe" if os.name == "nt" else "OpenSees"\n'
    assert _run(tmp_path, body) == []


def test_ignores_non_zone_a_files_and_text_mentions(tmp_path):
    # zone_b / unmarked files are out of scope; a platform test named only in a
    # docstring or comment is not a branch (ast, not grep).
    assert _run(tmp_path, 'import sys\nif sys.platform != "win32":\n    pass\n') == []
    body = _ZA + '"""we used to skip when sys.platform != "win32"."""\n# os.name == "nt"\n'
    assert _run(tmp_path, body) == []


def test_rejects_unknown_kind_and_short_reason(tmp_path):
    gate = 'pytestmark.append(pytest.mark.skipif(sys.platform != "win32", reason="mkl"))\n'
    out = _run(tmp_path, _ZA + "# ci-coverage: sometimes whenever the moon is right\n" + gate)
    assert len(out) == 1 and "not one of" in out[0], out
    out = _run(tmp_path, _ZA + "# ci-coverage: local-only mkl\n" + gate)
    assert len(out) == 1 and "too short" in out[0], out


def test_flags_a_stale_annotation(tmp_path):
    out = _run(tmp_path, _ZA + "# ci-coverage: local-only the Pardiso leg used to be Windows-only\n")
    assert len(out) == 1 and out[0].startswith("W ") and "stale ci-coverage" in out[0], out


def test_reports_an_unparseable_test_file(tmp_path):
    out = _run(tmp_path, _ZA + "def broken(:\n")
    assert len(out) == 1 and "cannot parse" in out[0], out
