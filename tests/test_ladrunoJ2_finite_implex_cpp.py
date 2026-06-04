"""Cross-check the C++ LadrunoJ2Finite IMPL-EX path against the numpy oracle WITHOUT
building OpenSees.

Compiles tests/ladrunoj2_finite_implex_check.cpp (pure — OpenSees-free kernels only,
mirroring SRC/material/nD/LadrunoJ2Finite.cpp::setTrialF in the useImplex branch),
runs it over the rotating+stretching combined-hardening plastic path, and asserts the
EXPLICIT (implex) Cauchy stress, committed accumulated plastic strain, and committed
total backstress match the independent oracle (tests/ladrunoj2_finite_implex_reference).

This guards the C++ convention bridging the 3×3-tensor oracle cannot: the engineering↔
tensor Voigt factors in the elastic explicit stress, the flow-direction recovery
N=Δεᵖ/Δγ, and the polar co-rotation of the frozen flow direction. Skipped if no C++
compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import ladrunoj2_finite_implex_reference as fx
import ladrunoj2_reference as lr

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))

PARAMS = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                   C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])


def _rotz(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _path():
    Fs = []
    for k in range(1, 10):
        s = 1.0 + 0.03 * k
        gg = 0.015 * k
        U = np.array([[s, gg, 0.0], [0.0, 1.0/np.sqrt(s), 0.0],
                      [0.0, 0.0, 1.0/np.sqrt(s)]])
        Fs.append(_rotz(0.07 * k) @ U)
    return Fs


def _sym_to_voigt(T):  # {00,11,22,01,12,02}
    return np.array([T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[1, 2], T[0, 2]])


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladrunoj2_finite_implex_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lj2fimx.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, src, "-o", exe],
                            capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    rows = {}
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        rows[int(p[0])] = [float(x) for x in p[1:]]   # sig6(6) ebarP(1) tot6(6)
    return rows


def test_cpp_implex_matches_oracle():
    cpp = _run_cpp()
    m = fx.ImplexNativeFiniteJ2(PARAMS)
    saw_plastic = False
    for k, F in enumerate(_path(), start=1):
        sig = m.setTrialF(F)                       # IMPL-EX (explicit) Cauchy stress
        saw_plastic = saw_plastic or m._plastic
        m.commit()

        sig6 = _sym_to_voigt(sig)
        tot6 = _sym_to_voigt(sum(m.alpha))
        ref = np.array(sig6.tolist() + [m.ebarP] + tot6.tolist())
        got = np.array(cpp[k])

        sscale = max(np.abs(ref[:6]).max(), 1.0)
        np.testing.assert_allclose(got[:6], ref[:6], rtol=1e-9, atol=1e-9 * sscale,
                                   err_msg=f"IMPL-EX Cauchy stress mismatch at step {k}")
        assert got[6] == pytest.approx(ref[6], rel=1e-9, abs=1e-12), \
            f"committed ebarP mismatch at step {k}"
        ascale = max(np.abs(ref[7:]).max(), 1.0)
        np.testing.assert_allclose(got[7:], ref[7:], rtol=1e-9, atol=1e-9 * ascale,
                                   err_msg=f"committed backstress mismatch at step {k}")
    assert saw_plastic, "path never yielded — scenario is vacuous"
