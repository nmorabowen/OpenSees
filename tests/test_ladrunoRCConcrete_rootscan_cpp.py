"""CI wiring for the standalone g++ gate of Ladruno WP-247 (LadrunoRCKernel.h, driveBeamFiber root scan).

When the BeamFiber view's lateral Newton stalls, the kernel falls back to Model C's root rule: the nearest
increasing sign change of the tied-direction residual from the committed lateral strain, bisected, with the other
lateral unknowns balanced at every sample. tests/_testbed/rc_rootscan_gpp.cpp gates it on a lightly tied core fibre
(rho_t f_yt < f_cr): the crack-onset step that the pre-WP-247 kernel refused converges through the scan and only
there (R1), every step matches nested bracketing bisection from the same committed state (R2), gamma_xz with ties
in z is the same fibre (R3), a long cracking ramp converges throughout (R4), and the lateral box keeps a cover
fibre that is fully cracked in y on its root instead of eps_yy ~ 3e4 (R5). Self-asserting with a non-zero
exit. Skipped without a C++ compiler.
"""
import os
import shutil
import subprocess
import tempfile

import pytest
from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))
_GATE = os.path.join(_HERE, "_testbed", "rc_rootscan_gpp.cpp")


def test_rootscan_gpp_gate():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "rc_rs.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, _GATE, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
        assert cc.returncode == 0, f"g++ failed to compile the kernel gate:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert rr.returncode == 0, f"g++ root-scan gate FAILED:\n{rr.stdout}\n{rr.stderr}"
    out = rr.stdout
    assert "GPP ROOTSCAN GATE: PASS" in out, out
    for head in ("R1 fold step", "R2 reference", "R3 invariance", "R4 long ramp", "R5 lateral box"):
        assert head in out, out
