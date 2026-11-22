"""CI wiring for the standalone g++ gate of Ladruno WP-243 (LadrunoRCKernel.h, driveBeamFiber).

The BeamFiber view's lateral solve and condensed tangent live in the header-only kernel, so they
are gated without OpenSees by tests/_testbed/rc_beamfiber_gpp.cpp: closed forms with and without
ties (B1, B2), a two-leg check against nested bracketing bisection from the same committed state
(B3), the condensed tangent vs central FD (B4), gamma_xy/gamma_xz invariance (B5), the crack gate
in unconfined compression with the strain gate's refusal as break-once evidence (B6), the tie
law (B7) and the Model C `sum` crack check (B8-B10). Self-asserting with a non-zero exit. Same pattern as test_ladrunoRCConcrete_spatial_cpp.py.
Skipped without a C++ compiler.
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
_GATE = os.path.join(_HERE, "_testbed", "rc_beamfiber_gpp.cpp")


def test_beamfiber_gpp_gate():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "rc_bf.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, _GATE, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
        assert cc.returncode == 0, f"g++ failed to compile the kernel gate:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert rr.returncode == 0, f"g++ BeamFiber gate FAILED:\n{rr.stdout}\n{rr.stderr}"
    out = rr.stdout
    assert "GPP BEAMFIBER GATE: PASS" in out, out
    for head in ("B1 elastic, no ties", "B2 elastic, ties", "B3 two-leg", "B4 condensed FD (elastic)",
                 "B4 condensed FD (cracked)", "B4 condensed FD (tie-yielded)", "B5 invariance",
                 "B6 damage gate", "B6 strain gate", "B7 tie law", "B7 cracked fibre",
                 "B8 crack-check cap", "B9 crack check on a shear ramp", "B10 crack-check invariance"):
        assert head in out, out
