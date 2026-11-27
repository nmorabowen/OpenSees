"""WP-261: the LadrunoIMKBeam safeguarded local hinge solve, WITHOUT building
OpenSees.

Compiles tests/ladrunoIMK_safe_solve_check.cpp (includes only
SRC/element/ladrunoIMKBeam/LadrunoIMKSafeSolve.h) with g++ and checks, on two
Steel01-clip hinges in series with the elastic interior of the test element:

  * the plain local Newton of ``ladrunoIMKSolveAxis`` (transcribed) FAILS on a
    hinge unloading from a plastic committed state: the arctan 2-cycle on the
    S-shaped residual a(v - th) - M(th). This pins the bug the fallback fixes;
  * the safeguarded solve converges there and on a 4000-case sweep of loading
    and unloading steps, and matches the exact solution found by enumerating
    the Steel01 regions to round-off;
  * Steel01's absolute DBL_EPSILON dead zone (no root at v ~ 50 eps) is
    accepted through the rotation-resolution floor, not reported;
  * on a law with no root (an upward moment jump) it returns false, with its
    last evaluation at the returned point: the contract that lets
    ``update()`` fail with the materials still consistent with the rotations.

Skipped if no C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import pytest

from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.normpath(os.path.join(_HERE, "..", "SRC", "element", "ladrunoIMKBeam"))


@pytest.fixture(scope="module")
def out():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladrunoIMK_safe_solve_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "imk_safe.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC, src, "-o", exe,
                             *cxx_runtime_flags()], capture_output=True, text=True)
        assert cc.returncode == 0, cc.stderr
        run = subprocess.run([exe], capture_output=True, text=True, timeout=120)
        assert run.returncode == 0, run.stderr
    vals = {}
    for line in run.stdout.splitlines():
        k, v = line.split()
        vals[k] = float(v)
    return vals


def test_plain_newton_cycles_on_unloading(out):
    # the bug, reproduced: if this ever converges the harness no longer
    # exercises the 2-cycle and the other checks prove less
    assert out["plain_unload_converged"] == 0
    assert out["sweep_plain_fail"] > 100


def test_safeguarded_solve_converges_and_is_exact(out):
    assert out["safe_unload_converged"] == 1
    assert out["safe_unload_exact_found"] == 1
    assert out["safe_unload_err"] < 1e-10          # in units of My/k0
    assert out["safe_unload_last_at_x"] == 1
    assert out["safe_unload_evals"] <= 20


def test_safeguarded_sweep(out):
    assert out["sweep_safe_fail"] == 0
    assert out["sweep_cases"] == 4000
    assert out["sweep_safe_max_err"] < 1e-10
    assert out["sweep_safe_max_evals"] <= 100


def test_steel01_dead_zone_is_not_a_failure(out):
    # Steel01 ignores |d strain| <= DBL_EPSILON (absolute): at v ~ 50 eps the
    # hinge equation has no root, and the k0 * eps residual must be accepted
    assert out["deadzone_converged"] == 1
    assert out["deadzone_evals"] <= 50


def test_no_root_reports_failure_consistently(out):
    assert out["noroot_converged"] == 0
    assert out["noroot_last_at_x"] == 1
    assert out["noroot_x_finite"] == 1
