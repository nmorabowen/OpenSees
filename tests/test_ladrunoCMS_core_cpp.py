"""Compile and execute the ADR-1000 pure C++ assembly/CSR kernels."""

import os
import shutil
import subprocess
import tempfile

import pytest


pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_CMS = os.path.normpath(
    os.path.join(_HERE, "..", "SRC", "system_of_eqn", "ladrunoCMS")
)


def test_ladruno_cms_pure_cpp_core():
    compiler = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if compiler is None:
        pytest.skip("no C++ compiler available")

    source = os.path.join(_HERE, "ladruno_cms_core_check.cpp")
    implementation = os.path.join(_CMS, "LadrunoCMSOptions.cpp")
    with tempfile.TemporaryDirectory() as temporary:
        executable = os.path.join(temporary, "ladruno_cms_core")
        compile_result = subprocess.run(
            [
                compiler,
                "-O2",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-pedantic",
                "-I",
                _CMS,
                source,
                implementation,
                "-o",
                executable,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert compile_result.returncode == 0, (
            "C++ compilation failed:\n" + compile_result.stderr
        )
        run_result = subprocess.run(
            [executable], capture_output=True, text=True, check=False
        )
        assert run_result.returncode == 0, (
            "C++ core checks failed:\n"
            + run_result.stdout
            + run_result.stderr
        )
        assert "core checks passed" in run_result.stdout
