"""CI wiring for the standalone g++ crack-band regularization gate (RC shell Phase 3b).

The kernel's regularize() (post-peak strain rescaling so g_reg = G_f0*(lch_ref/lch)) is proven
against the numpy oracle by tests/_testbed/rc_reg_gpp.cpp (self-asserting, non-zero exit on a
fracture-energy mismatch). This wrapper compiles+runs it from CI so the gate is enforced — the
same pattern as test_ladrunoRCConcrete_tensstiff_cpp.py / test_ladrunoJ2_*_cpp.py. Skipped if no
C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import pytest

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))
_GATE = os.path.join(_HERE, "_testbed", "rc_reg_gpp.cpp")


def test_reg_gpp_gate():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "rc_reg.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, _GATE, "-o", exe],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed to compile the regularization gate:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
    assert rr.returncode == 0, f"g++ regularization gate FAILED:\n{rr.stdout}\n{rr.stderr}"
    assert "GPP REG GATE: PASS" in rr.stdout, rr.stdout
    assert "no-op: max|delta|=0.00e+00" in rr.stdout, rr.stdout      # lch==lch_ref is a no-op
