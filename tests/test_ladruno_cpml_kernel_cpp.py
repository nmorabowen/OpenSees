"""Standalone C++17 verification test for LadrunoCPMLShapes.h and
LadrunoCPMLKernel.h (Phase P0 of LadrunoCPML) WITHOUT building OpenSees.

Compiles tests/ladruno_cpml_kernel_check.cpp with g++/clang++ or MSVC cl.exe,
runs `--self-check` (shape invariants, finite-difference tangent consistency,
SPD checks, static/interior/classical limits, input validation,
revertToLastCommit/revertToStart, and packState/unpackState round-trip across
all 9 topologies), runs the blind-oracle comparison against
tests/ladruno_cpml_cases.txt, and verifies all 5 CPML mutation flags are killed
by `--self-check`.
"""
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_SRC_CPML = os.path.normpath(os.path.join(_REPO_ROOT, "SRC", "element", "ladrunoCPML"))
_CASES = os.path.join(_HERE, "ladruno_cpml_cases.txt")
_UP_CASES = os.path.join(_HERE, "ladruno_cpml_up_cases.txt")
_CHECK_CPP = os.path.join(_HERE, "ladruno_cpml_kernel_check.cpp")

_MUTATIONS: List[str] = [
    "LADRUNO_MUTATION_CPML_NO_DAMPING",
    "LADRUNO_MUTATION_CPML_NO_KAPPA",
    "LADRUNO_MUTATION_CPML_NO_ALPHA",
    "LADRUNO_MUTATION_CPML_WRONG_SIGN_MEMORY",
    "LADRUNO_MUTATION_CPML_NO_CORNER_CROSS",
]


def _compile_all(td: str) -> Dict[str, str]:
    """Compile baseline + 5 mutant binaries in td (loading setup_env.bat once on Windows)."""
    targets: Dict[str, List[str]] = {"baseline": []}
    for m in _MUTATIONS:
        targets[m] = [m]

    exes: Dict[str, str] = {}
    for name in targets:
        exes[name] = os.path.join(td, f"{name}.exe" if os.name == "nt" else name)

    gpp = shutil.which("g++") or shutil.which("clang++") or shutil.which("c++")
    if gpp is not None:
        for name, defs in targets.items():
            dflags = [f"-D{d}=1" for d in defs]
            cc = subprocess.run(
                [gpp, "-O2", "-std=c++17", "-Wall", "-Wextra", *dflags, "-I", _SRC_CPML, _CHECK_CPP, "-o", exes[name]],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
            assert cc.returncode == 0, f"{gpp} compilation failed ({name}):\n{cc.stdout}\n{cc.stderr}"
        return exes

    cl = shutil.which("cl") or shutil.which("cl.exe")
    if cl is not None:
        for name, defs in targets.items():
            dflags = [f"/D{d}=1" for d in defs]
            cc = subprocess.run(
                [cl, "/nologo", "/EHsc", "/O2", "/std:c++17", *dflags, f"/I{_SRC_CPML}", _CHECK_CPP, f"/Fe:{exes[name]}"],
                cwd=td,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
            assert cc.returncode == 0, f"cl.exe compilation failed ({name}):\n{cc.stdout}\n{cc.stderr}"
        return exes

    # ci-coverage: portable Windows MSVC fallback when g++/clang++/cl are not on PATH
    if os.name == "nt":
        setup_env = os.path.join(_REPO_ROOT, "Ladruno_scripts", "setup_env.bat")
        if os.path.isfile(setup_env):
            bat_path = os.path.join(td, "compile_all_cpml.bat")
            with open(bat_path, "w", encoding="utf-8") as bf:
                bf.write("@echo off\n")
                bf.write(f'call "{setup_env}" >nul 2>&1\n')
                bf.write("if errorlevel 1 exit /b 1\n")
                for name, defs in targets.items():
                    dflags_str = " ".join(f"/D{d}=1" for d in defs)
                    bf.write(
                        f'cl /nologo /EHsc /O2 /std:c++17 {dflags_str} '
                        f'/I"{_SRC_CPML}" "{_CHECK_CPP}" /Fe:"{exes[name]}" /Fo:"{td}\\{name}.obj"\n'
                    )
                    bf.write("if errorlevel 1 exit /b 1\n")
            cc = subprocess.run(
                ["cmd.exe", "/d", "/c", bat_path],
                cwd=td,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
            if cc.returncode == 0 and all(os.path.isfile(p) for p in exes.values()):
                return exes

    pytest.skip("No C++17 compiler (g++, clang++, or cl.exe) available")
    return exes


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_ladruno_cpml_kernel_selfcheck_and_oracle():
    with tempfile.TemporaryDirectory() as td:
        exes = _compile_all(td)
        exe = exes["baseline"]

        # 1. Always verify standalone --self-check
        res_self = subprocess.run(
            [exe, "--self-check"],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
        )
        assert res_self.returncode == 0, (
            f"--self-check failed:\n{res_self.stdout}\n{res_self.stderr}"
        )
        assert "RESULT: PASS" in res_self.stdout, res_self.stdout

        # 2. Verify against oracle cases file(s) when present
        if os.path.isfile(_CASES):
            res_oracle = subprocess.run(
                [exe, _CASES],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
            assert res_oracle.returncode == 0, (
                f"Oracle comparison failed:\n{res_oracle.stdout}\n{res_oracle.stderr}"
            )
            assert "RESULT: PASS" in res_oracle.stdout, res_oracle.stdout
            m = re.search(r"\[cases\] verified (\d+) case", res_oracle.stdout)
            assert m and int(m.group(1)) >= 18, (
                f"Expected >= 18 cases verified from {_CASES}:\n{res_oracle.stdout}"
            )

        if os.path.isfile(_UP_CASES):
            res_up = subprocess.run(
                [exe, _UP_CASES],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
            assert res_up.returncode == 0, (
                f"Biot u-p oracle comparison failed:\n{res_up.stdout}\n{res_up.stderr}"
            )
            assert "RESULT: PASS" in res_up.stdout, res_up.stdout
            m_up = re.search(r"\[cases\] verified (\d+) case", res_up.stdout)
            assert m_up and int(m_up.group(1)) >= 10, (
                f"Expected >= 10 cases verified from {_UP_CASES}:\n{res_up.stdout}"
            )

        # 3. Verify all 5 CPML physics mutations are killed by --self-check
        for mut in _MUTATIONS:
            res_mut = subprocess.run(
                [exes[mut], "--self-check"],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
            assert res_mut.returncode != 0, (
                f"Mutant {mut} survived --self-check!\n{res_mut.stdout}"
            )
            assert "RESULT: FAIL" in res_mut.stdout, res_mut.stdout
