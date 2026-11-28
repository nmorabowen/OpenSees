"""WP-264: LadrunoRebarBuckling's GA factor derivative by dual numbers.

Compiles tests/ladruno_rebarbuckling_dual_check.cpp (pure: the kernel header
SRC/material/uniaxial/LadrunoRebarBucklingKernel.h + LadrunoDual.h) with g++ and
asserts its four gates over 72 (lsr, reduction, e_cross) combinations x 400
strains spanning the buckled region:

  * bits: the templated factor in double equals a verbatim transcription of the
    pre-WP-264 LadrunoRebarBuckling::gaFactor bit for bit, and the value of the
    Dual<1> pass equals it too (the GA stress cannot move);
  * deriv: dfactor/deps by Dual<1> equals a step-plateau central FD to 1e-6 off
    the factor's kinks (it replaces the shipped h = 1e-8 FD in the tangent);
  * mut_beta / mut_sqrt (break-the-gate): the derivative with the beta term
    frozen, or with the sqrt chain rule dropped, must miss the FD.
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
_SRC = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material"))


@pytest.fixture(scope="module")
def gates():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_rebarbuckling_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lrbd.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-Wall",
                             "-I", os.path.join(_SRC, "nD"), "-I", os.path.join(_SRC, "uniaxial"),
                             src, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
    out = {}
    for line in rr.stdout.splitlines():
        p = line.split()
        if p and p[0] == "GATE":
            out[p[1]] = (int(p[2]), float(p[3]))
    assert set(out) == {"bits", "deriv", "mut_beta", "mut_sqrt"}, rr.stdout
    return out, rr.stdout


def test_factor_is_bit_identical_to_the_pre_wp_code(gates):
    g, log = gates
    assert g["bits"][0] == 1, log


def test_dual_dfactor_matches_plateau_fd(gates):
    g, log = gates
    assert g["deriv"][0] == 1, log
    assert g["deriv"][1] < 1e-6, log


@pytest.mark.parametrize("mutant", ["mut_beta", "mut_sqrt"])
def test_break_the_gate_mutants_are_flagged(gates, mutant):
    g, log = gates
    assert g[mutant][0] == 1, f"{mutant} was not flagged:\n{log}"
