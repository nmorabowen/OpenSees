"""Self-test for quirk-lint rule `gxx-static` (WP-199) in ci/check_quirk_patterns.py.

One file per rule (WP-162). Run: pytest -q ci/test_quirk_gxx_static.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402

HEAD = "import shutil\nimport subprocess\n"


def _run(tmp_path, body, name="tests/test_x_cpp.py"):
    root = kit.tree(tmp_path, {name: body})
    return cq.check_gxx_static(root, kit.rel(root))


@pytest.mark.parametrize("lookup", [
    'shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")',   # the common form
    "shutil.which('g++', path=env.get('PATH'))",                               # the membrane form
    'which("clang++")',                                                         # from shutil import which
])
def test_flags_a_dynamic_harness_build(tmp_path, lookup):
    body = HEAD + f"gpp = {lookup}\nsubprocess.run([gpp, '-O2', 'a.cpp', '-o', 'a.exe'])\n"
    out = _run(tmp_path, body)
    assert len(out) == 1 and out[0].startswith("gxx-static tests/test_x_cpp.py:3"), out


def test_flags_an_imported_but_unused_helper(tmp_path):
    # importing the helper is not enough: the compile command must use it
    body = HEAD + ("from _cxx import cxx_runtime_flags\n"
                   "gpp = shutil.which('g++')\nsubprocess.run([gpp, 'a.cpp', '-o', 'a.exe'])\n")
    assert len(_run(tmp_path, body)) == 1


@pytest.mark.parametrize("flags", ["*cxx_runtime_flags()", "*_cxx.cxx_runtime_flags()", "'-static'"])
def test_passes_a_static_build(tmp_path, flags):
    body = HEAD + f"gpp = shutil.which('g++')\nsubprocess.run([gpp, 'a.cpp', '-o', 'a.exe', {flags}])\n"
    assert _run(tmp_path, body) == []


def test_ignores_comments_docstrings_and_other_tools(tmp_path):
    body = ('"""Built with shutil.which("g++") in the docstring."""\n' + HEAD
            + "# gpp = shutil.which('g++')\ncl = shutil.which('cl')\n")
    assert _run(tmp_path, body) == []


def test_the_helper_module_itself_is_exempt(tmp_path):
    assert _run(tmp_path, HEAD + "x = shutil.which('g++')\n", name="tests/_cxx.py") == []


def test_the_fork_harnesses_link_statically():
    root = Path(__file__).resolve().parent.parent
    assert cq.check_gxx_static(root, kit.rel(root)) == []
