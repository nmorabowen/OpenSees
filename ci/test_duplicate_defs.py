"""Self-test for the duplicate-definition gate (WP-162): ci/check_duplicate_defs.py
(D1) and `ruff check --select F811`, both run by the static-gates step.

The incident: the #899 merge-up of #901 (WP-153 and WP-158 had both taken quirk
rule L9) merged ci/test_check_quirk_patterns.py with NO conflict and two
`def _l9`. INCIDENT below is that file's shape, reduced. F811 does not fire on
it (the first `_l9` is "used" by the test bodies in between), which is why D1
exists; the ruff cases pin that finding so a future ruff that does catch it
shows up here.
Run: pytest -q ci/test_duplicate_defs.py   (the ruff cases skip without ruff)
"""
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_duplicate_defs as dd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

INCIDENT = '''
def _l9(tmp_path, files):
    return ["the WP-153 revert helper"]


def test_l9_flags_an_integrator(tmp_path):
    assert _l9(tmp_path, {})


def _l9(tmp_path, body, stamped=False):
    return ["the WP-158 dead-decl helper"]


def test_l9_flags_the_wp158_incident(tmp_path):
    assert _l9(tmp_path, "x")
'''


def _d1(tmp_path, text, name="test_x.py"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return dd.find(p, tmp_path)


# ---------------------------------------------------------------- D1
def test_d1_flags_the_merged_helpers(tmp_path):
    out = _d1(tmp_path, INCIDENT)
    assert len(out) == 1 and out[0].startswith("D1 test_x.py:10:") and "'_l9'" in out[0], out


@pytest.mark.parametrize("text", [
    "class A:\n    def f(self):\n        pass\n\n    def f(self):\n        pass\n",       # a method
    "def test_a():\n    pass\n\n\ndef test_a():\n    pass\n",                          # a test name
    "class K:\n    pass\n\n\nclass K:\n    pass\n",                                     # a class
    "def outer():\n    def g():\n        pass\n    def g():\n        pass\n",          # a nested scope
])
def test_d1_flags_other_redefinitions(tmp_path, text):
    assert len(_d1(tmp_path, text)) == 1


@pytest.mark.parametrize("text", [
    "import sys\nif sys.platform == 'win32':\n    def f():\n        pass\nelse:\n    def f():\n        pass\n",
    "try:\n    from x import f\nexcept ImportError:\n    def f():\n        pass\n",
    "class P:\n    @property\n    def v(self):\n        return 1\n\n    @v.setter\n    def v(self, x):\n        pass\n",
    "from typing import overload\n@overload\ndef f(x: int) -> int: ...\n@overload\ndef f(x: str) -> str: ...\n"
    "def f(x):\n    return x\n",
    "def f():\n    pass\n\n\ndef g():\n    def f():\n        pass\n",               # different scopes
])
def test_d1_passes_alternatives_and_accessors(tmp_path, text):
    assert _d1(tmp_path, text) == []


def test_d1_waiver(tmp_path):
    ok = INCIDENT.replace("\n\ndef _l9(tmp_path, body", "\n\n# ladruno-lint: redef-ok intentional override for a fixture\ndef _l9(tmp_path, body")
    assert _d1(tmp_path, ok) == []
    short = INCIDENT.replace("\n\ndef _l9(tmp_path, body", "\n\n# ladruno-lint: redef-ok ok\ndef _l9(tmp_path, body")
    out = _d1(tmp_path, short, "test_y.py")
    assert len(out) == 1 and "too short" in out[0]


def test_d1_the_tree_is_clean():
    assert dd.main([]) == 0


# ---------------------------------------------------------------- ruff F811
def _ruff(*args, cwd=ROOT):
    return subprocess.run([sys.executable, "-m", "ruff", *args], cwd=cwd, stdin=subprocess.DEVNULL,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")


@pytest.fixture
def ruff():
    try:
        if _ruff("--version").returncode != 0:
            pytest.skip("ruff not installed (CI pins it; pip install ruff==0.16.10 to run locally)")
    except OSError:
        pytest.skip("cannot run python -m ruff")


def test_f811_misses_the_incident_which_is_why_d1_exists(tmp_path, ruff):
    (tmp_path / "test_merged.py").write_text(INCIDENT, encoding="utf-8")
    r = _ruff("check", "--no-cache", "--isolated", "--select", "F811", "test_merged.py", cwd=tmp_path)
    assert r.returncode == 0, "ruff F811 now catches the _l9 merge: D1 may be redundant -- " + r.stdout


def test_f811_flags_a_duplicated_test_name(tmp_path, ruff):
    (tmp_path / "test_twice.py").write_text("def test_w():\n    assert 1\n\n\ndef test_w():\n    assert 2\n",
                                            encoding="utf-8")
    r = _ruff("check", "--no-cache", "--isolated", "--select", "F811", "test_twice.py", cwd=tmp_path)
    assert r.returncode == 1 and "F811" in r.stdout, r.stdout + r.stderr


def test_f811_the_tree_is_clean(ruff):
    r = _ruff("check", "--no-cache", "--select", "F811", *dd.DEFAULT)
    assert r.returncode == 0, r.stdout + r.stderr
