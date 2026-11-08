"""WP-224: the WP-99 commit-refusal counter and the probe-scope depth are safe under the WP-192 threaded commit.

Compiles tests/ladruno_commit_refusal_threads_check.cpp (header-only, includes SRC/material/LadrunoMaterialStatus.h)
with g++ -O2 -pthread and runs it: 20 rounds of 8 threads x 200000 ladrunoNoteCommitRefusal() must count exactly
1.6e6 each, and a refusal declared on thread B while thread A holds a LadrunoProbeCommitScope must be counted.
On the pre-WP-224 header (plain static ints) the same program reports lost increments and the suppression.
Skipped without a C++ compiler (CI Zone-A/Ubuntu has g++).
"""
import os
import shutil
import subprocess
import tempfile

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_MAT = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material"))


def test_commit_refusal_seam_is_thread_safe():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_commit_refusal_threads_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "crt.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-pthread", "-I", _SRC_MAT, src, "-o", exe],
                            capture_output=True, text=True)
        assert cc.returncode == 0, cc.stderr
        run = subprocess.run([exe], capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stdout + run.stderr
    assert "lost increments" in run.stdout and ": 0;" in run.stdout, run.stdout
    assert "suppression: no" in run.stdout, run.stdout
