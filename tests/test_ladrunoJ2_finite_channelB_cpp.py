"""Cross-check the C++ ANALYTIC channel-B block of LadrunoJ2Finite (ADR-16 P3) against
the numpy oracle WITHOUT building OpenSees — validates the C++ transcription of the
polar-rotation derivative + return-map backstress sensitivity (and the F→L push).

Channel B is purely constitutive (no geometric term), so the assembled
cmatB_ijkl = Σ_m (∂σ_ij/∂F_km) F_lm compares directly, with zero OpenSees commit
semantics, to the numpy oracle's channelB_dsigma_dF mapped through the IDENTICAL
Σ_m·F_lm push-forward. The analytic C++ path should match the (convention-free) oracle
to far better than the numeric path, since both are exact (not FD) — an index/factor/
sign/Sylvester-sign bug in LadrunoJ2Finite.cpp's analytic channel B is O(1) and fails
this. Skipped if no C++ compiler is available.
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


def _diag_stretch(a):
    return np.diag([np.exp(a), np.exp(-a / 2), np.exp(-a / 2)])


def _rot(ax, ay, az, th):
    n = np.array([ax, ay, az], float); n /= np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def _base_and_eval():
    base = [_diag_stretch(0.05 * s) for s in (1, 2, 3)]          # must match the .cpp
    Feval = _rot(0.3, 1.0, 0.5, 0.35) @ _diag_stretch(0.20)
    return base, Feval


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladrunoj2_finite_channelB_analytic_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lj2cba.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, src, "-o", exe],
                            capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    out = {"PLASTIC": None, "FEVAL": None, "CB": None}
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "PLASTIC":
            out["PLASTIC"] = int(p[1])
        elif p[0] == "FEVAL":
            out["FEVAL"] = np.array([float(x) for x in p[1:]]).reshape(3, 3)
        elif p[0] == "CB":
            out["CB"] = np.array([float(x) for x in p[1:]]).reshape(3, 3, 3, 3)
    return out


def test_cpp_analytic_channelB_matches_oracle():
    cpp = _run_cpp()
    assert cpp["PLASTIC"] == 1, "eval step is elastic — channel B would be trivially zero"
    Fe = cpp["FEVAL"]
    cB_cpp = cpp["CB"]

    base, Feval = _base_and_eval()
    assert np.allclose(Fe, Feval, rtol=1e-12, atol=1e-12), "C++/py eval F mismatch"

    m = nat.NativeFiniteJ2(PARAMS, corotate=True)
    for F in base:
        m.setTrialF(F); m.commit()
    _, plastic = m.setTrialF(Feval)
    assert plastic, "oracle eval step not plastic"

    B = m.channelB_dsigma_dF(Feval)                  # ∂σ_ij/∂F_kl (R-only), numeric oracle
    cB_oracle = np.einsum("ijkm,lm->ijkl", B, Feval) # same Σ_m·F_lm push as the C++

    scale = np.linalg.norm(cB_oracle)
    assert scale > 1.0, f"channel-B tensor vanished (norm {scale:.2e}) — test is vacuous"
    relerr = np.linalg.norm(cB_cpp - cB_oracle) / scale
    # C++ analytic vs numpy NUMERIC oracle: the residual is the oracle's own FD error
    # (~1e-7, h=1e-7) plus transcription; an index/factor/Sylvester-sign bug is O(1).
    assert relerr < 1.0e-5, (
        f"C++ analytic channel-B tensor != oracle: rel Frobenius err {relerr:.3e}")
