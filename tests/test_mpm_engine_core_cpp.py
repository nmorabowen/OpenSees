"""WP-232c: the explicit MPM engine core and the incremental hypo-kernel entry,
checked WITHOUT building OpenSees (the tests/test_hypo_kernel_cpp.py idiom).

* tests/hypo_incremental_kernel_check.cpp includes only
  SRC/element/solidTransformation/LadrunoHypoKernel.h and checks
  hypoIncrementFromF / hypoIncrementFromFPlaneStrain against
  tests/hypo_incremental_cases.txt (232b's oracle, an SVD polar) and against the
  ADR-79 nodal-displacement entry.
* tests/mpm_engine_check.cpp includes only
  SRC/domain/pattern/ladrunoMPM/LadrunoMPMCore.h and runs every case of
  tests/mpm_cases.txt with the interface's hypoelastic law: basis, critical dt,
  mass cutoff, refusals, particle_fields after 1/3/5 steps, the closed forms
  (translation, bar period, PIC vs FLIP, impact, incline, K0, small mass).

Neither file reads tests/mpm_reference.py (232b's oracle stays independent).
This file is deliberately NOT named test_ladruno_mpm_*: the mutation gate's MPM
suite must hold only tests that run the mutated binary, and these compile their
own unmutated copy of the core.
"""
import os
import re
import shutil
import subprocess
import tempfile

import pytest

pytestmark = pytest.mark.zone_a

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_KERNEL_DIR = os.path.join(_ROOT, "SRC", "element", "solidTransformation")
_CORE_DIR = os.path.join(_ROOT, "SRC", "domain", "pattern", "ladrunoMPM")

# Interface A4 (2026-11-13): the 20%-compression self-weight column misses its
# frozen 0.03 in the oracle too (0.072, explicit-GIMP cell-crossing noise); a
# strict xfail until the orchestrator rules on the gate.
_KNOWN_XFAIL = {"column_selfweight_2d_gimp selfweight_error"}


def _compile_and_run(src, cases, extra_inc):
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "check.exe")
        cmd = [gpp, "-O2", "-std=c++17", "-Wall", "-Wextra"]
        for d in extra_inc:
            cmd += ["-I", d]
        cc = subprocess.run(cmd + [src, "-o", exe], stdin=subprocess.DEVNULL,
                            capture_output=True, text=True, encoding="utf-8", errors="replace")
        assert cc.returncode == 0, f"compile failed:\n{cc.stderr}"
        rr = subprocess.run([exe, cases], stdin=subprocess.DEVNULL, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
        return rr.returncode, rr.stdout + rr.stderr


@pytest.mark.t1
def test_cpp_hypo_incremental_entry():
    rc, out = _compile_and_run(os.path.join(_HERE, "hypo_incremental_kernel_check.cpp"),
                               os.path.join(_HERE, "hypo_incremental_cases.txt"), [_KERNEL_DIR])
    assert rc == 0 and "RESULT: PASS" in out, out
    m = re.search(r"\[cases\] parsed (\d+) cases", out)
    assert m and int(m.group(1)) >= 26, out


@pytest.fixture(scope="module")
def engine_run():
    return _compile_and_run(os.path.join(_HERE, "mpm_engine_check.cpp"),
                            os.path.join(_HERE, "mpm_cases.txt"), [_CORE_DIR, _KERNEL_DIR])


def _fail_lines(out):
    keys = []
    for line in out.splitlines():
        if line.startswith("FAIL"):
            parts = line.split()
            keys.append(" ".join(parts[1:3]) if len(parts) >= 3 else line)
    return keys


@pytest.mark.t1
def test_cpp_mpm_engine_core_cases(engine_run):
    _, out = engine_run
    m = re.search(r"\[cases\] parsed (\d+) cases, (\d+) checks", out)
    assert m and int(m.group(1)) >= 30 and int(m.group(2)) >= 80, out
    unexpected = [k for k in _fail_lines(out) if k not in _KNOWN_XFAIL]
    assert not unexpected, "\n".join(unexpected) + "\n\n" + out


@pytest.mark.xfail(strict=True, reason="interface A4: uGIMP self-weight column at 20% "
                   "compression misses the frozen 0.03 (oracle 0.072)")
def test_cpp_mpm_engine_core_selfweight_column(engine_run):
    _, out = engine_run
    assert "column_selfweight_2d_gimp selfweight_error" not in _fail_lines(out), out
