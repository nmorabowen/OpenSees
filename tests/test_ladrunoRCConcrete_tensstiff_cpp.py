"""CI wiring for the standalone g++ tension-stiffening gate (RC shell Phase 3a).

The Phase-3 CONSISTENT TANGENT block in LadrunoRCKernel.h is the most error-prone code in the
feature, but a converged static stress is the residual root and therefore tangent-independent —
no OpenSees Zone-A test can catch a wrong TS tangent. The only tangent gate is the standalone
g++ driver tests/_testbed/rc_tensstiff_gpp.cpp (FD on the pinned direction + the degenerate
analytic cross-check + the floor closed-form, all self-asserting with a non-zero exit on
failure). This wrapper compiles and runs it from CI so the gate is enforced, not orphaned —
the same pattern as tests/test_ladrunoJ2_*_cpp.py. Skipped if no C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import pytest

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))
_GATE = os.path.join(_HERE, "_testbed", "rc_tensstiff_gpp.cpp")


def test_tensstiff_gpp_gate():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "rc_ts.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, _GATE, "-o", exe],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed to compile the kernel gate:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
    # the gate self-asserts: it returns non-zero on any floor/tangent mismatch.
    assert rr.returncode == 0, f"g++ tension-stiffening gate FAILED:\n{rr.stdout}\n{rr.stderr}"
    out = rr.stdout
    assert "GPP TENSSTIFF GATE: PASS" in out, out
    # the headline gates must actually have fired (not vacuously passed)
    assert "mismatches: 0" in out, out                                  # uniaxial closed-form
    assert "equibiaxial: pinned-both steps=" in out, out                # degen floor (both normals)
    assert "degen tangent analytic cross-check: worst rel-err=0.00e+00" in out, out
