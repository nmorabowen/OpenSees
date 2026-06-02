"""Cross-check the C++ finite-strain-native J2 composition against the numpy
oracle WITHOUT building OpenSees.

Compiles tests/ladrunoj2_finite_native_check.cpp (pure — includes only the
OpenSees-free LadrunoJ2Kernel.h + LogStrainKernel.h, and mirrors line-for-line the
math of SRC/material/nD/LadrunoJ2Finite.cpp::setTrialF) with g++, runs it over a
rotating+stretching combined-hardening plastic path, and asserts that the committed
Cauchy stress, accumulated plastic strain, and total (co-rotated) backstress match
the independent 3×3-tensor oracle (tests/ladrunoj2_finite_native_reference.py).

This is the fast guard on the C++ convention bridging (engineering↔tensor Voigt,
stress/strain shear factors, polar co-rotation) that the numpy oracle — which works
in 3×3 tensors throughout — cannot catch. Skipped if no C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import ladrunoj2_reference as lr
import ladrunoj2_finite_native_reference as nat

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
    src = os.path.join(_HERE, "ladrunoj2_finite_native_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lj2fin.exe")
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


def test_cpp_native_matches_oracle():
    cpp = _run_cpp()
    m = nat.NativeFiniteJ2(PARAMS, corotate=True)
    saw_plastic = False
    for k, F in enumerate(_path(), start=1):
        sig, pl = m.setTrialF(F)
        m.commit()
        saw_plastic = saw_plastic or pl

        sig6 = _sym_to_voigt(sig)
        tot6 = _sym_to_voigt(sum(m.alpha))
        ref = np.array(sig6.tolist() + [m.ebarP] + tot6.tolist())
        got = np.array(cpp[k])

        sscale = max(np.abs(ref[:6]).max(), 1.0)
        np.testing.assert_allclose(got[:6], ref[:6], rtol=1e-9, atol=1e-9 * sscale,
                                   err_msg=f"Cauchy stress mismatch at step {k}")
        assert got[6] == pytest.approx(ref[6], rel=1e-9, abs=1e-12), \
            f"ebarP mismatch at step {k}"
        ascale = max(np.abs(ref[7:]).max(), 1.0)
        np.testing.assert_allclose(got[7:], ref[7:], rtol=1e-9, atol=1e-9 * ascale,
                                   err_msg=f"backstress mismatch at step {k}")
    assert saw_plastic, "path never yielded — scenario is vacuous"
