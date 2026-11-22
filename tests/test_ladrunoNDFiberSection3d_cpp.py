"""CI wiring for the standalone g++ gate of Ladruno WP-245 (LadrunoFiberShearKernel.h).

The shear map of LadrunoNDFiberSection3d (strain map a, force map b, the uniform / parabolic
profiles and the modelC / conjugate resultants) is header-only, so it is gated without OpenSees by
tests/_testbed/rc_section_gpp.cpp, with the RC kernel's BeamFiber view as the nonlinear fibre:
elastic closed forms (S1, S2: modelC K_VV = GA (1 + 1/(2n^2)), conjugate -> 6/5 GA), symmetry
(S3), uniform = the stock NDFiberSection3d formulas (S4), the chain rule by central FD (S5), the
Model C section assembly on a cracking path (S6), Vy/Vz invariance (S7) and the band check (S8).
Self-asserting with a non-zero exit. Skipped without a C++ compiler.
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
_GATE = os.path.join(_HERE, "_testbed", "rc_section_gpp.cpp")


def test_section_gpp_gate():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "rc_sec.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", os.path.join(_SRC, "section"),
                             "-I", os.path.join(_SRC, "nD"), _GATE, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
        assert cc.returncode == 0, f"g++ failed to compile the section gate:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert rr.returncode == 0, f"g++ section gate FAILED:\n{rr.stdout}\n{rr.stderr}"
    out = rr.stdout
    assert "GPP SECTION GATE: PASS" in out, out
    for head in ("S1 elastic closed forms", "S2 conjugate k", "S3 symmetry", "S4 uniform = stock",
                 "S5 chain rule", "S6 Model C assembly", "S7 invariance", "S8 band check"):
        assert head in out, out
