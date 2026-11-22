"""CI wiring for the standalone g++ gate of Ladruno WP-241 (LadrunoRCKernel.h).

A converged stress is tangent-independent, so no Zone-A test can see a wrong tangent. The new
kernel code -- the spatial add-on plane (eps1 / p1 from the 3x3 strain, multiplicity-aware TS
vectors), the damage crack gate and the interlock tangent spanning all six columns -- is gated
by tests/_testbed/rc_spatial_gpp.cpp (reduce-to-membrane, xy/xz invariance with membrane
break-once, FD of the eps1 dual, TS pinning identity + FD, interlock all-column removal + FD
with a non-vacuous raw column, crack-gate strain vs damage). It self-asserts with a non-zero
exit. Same pattern as test_ladrunoRCConcrete_tensstiff_cpp.py. Skipped without a C++ compiler.
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
_GATE = os.path.join(_HERE, "_testbed", "rc_spatial_gpp.cpp")


def test_spatial_gpp_gate():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "rc_sp.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, _GATE, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
        assert cc.returncode == 0, f"g++ failed to compile the kernel gate:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert rr.returncode == 0, f"g++ spatial gate FAILED:\n{rr.stdout}\n{rr.stderr}"
    out = rr.stdout
    assert "GPP SPATIAL GATE: PASS" in out, out
    for head in ("G1 reduce", "G2 invariance", "G3 dual", "G4 TS pinning", "G5 interlock (2a",
                 "G5 interlock (cyclic", "G6 compression", "G6 tension", "G7 TS onset",
                 "G7 uniaxial tension"):
        assert head in out, out
