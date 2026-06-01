"""Cross-check the C++ LogStrain kernel against the numpy oracle WITHOUT building
OpenSees. Compiles tests/logstrain_kernel_check.cpp (pure, includes only
SRC/material/nD/LogStrainKernel.h) with g++, runs it, and asserts its Cauchy
stress σ and material spatial tangent c match logstrain_reference.py.

This catches C++ transcription bugs (a typo g++ accepts but computes wrong) in
the riskiest code — the Jacobi spectral solver, the A.52/A.53 degeneracy-safe
tensor-log derivative, and the spatial-tangent assembly. Skipped if g++ is
absent (CI Zone-A/Ubuntu has it).
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import logstrain_reference as lr

K, G = 1500.0, 700.0
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))
_REP_I = [0, 1, 2, 0, 1, 2]
_REP_J = [0, 1, 2, 1, 2, 0]

# F cases — must match logstrain_kernel_check.cpp exactly.
_cz, _sz = np.cos(0.5), np.sin(0.5)
CASES = {
    "generic":  np.array([[1.20, 0.10, -0.05], [0.02, 0.90, 0.08], [0.03, -0.04, 1.10]]),
    "uniaxial": np.diag([1.25, 0.95, 0.90]),
    "rotation": np.array([[_cz, -_sz, 0.0], [_sz, _cz, 0.0], [0.0, 0.0, 1.0]]),
    "double":   np.diag([1.2, 1.2, 0.8]),
    "triple":   np.diag([1.1, 1.1, 1.1]),
}


def _oracle(F):
    """Oracle Cauchy σ (6-vec) and material spatial tangent c (6×6) for the
    elastic Hencky law — the same quantities the C++ kernel emits."""
    I3 = np.eye(3)
    res = lr.matisu_step(F, I3.copy(), I3.copy(), "elastic", K=K, G=G)
    sigma, a_full = res["sigma"], res["a"]           # a_full = (1/2J)[D:L:B] − geo
    geo = np.einsum("il,jk->ijkl", sigma, I3)
    c_t = a_full + geo                                # material part (1/2J)[D:L:B]
    sig6 = np.array([sigma[_REP_I[I], _REP_J[I]] for I in range(6)])
    c6 = np.array([[c_t[_REP_I[I], _REP_J[I], _REP_I[J], _REP_J[J]]
                    for J in range(6)] for I in range(6)])
    return sig6, c6


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "logstrain_kernel_check.cpp")
    out = {}
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lsk.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, src, "-o", exe],
                            capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    # parse "CASE n / SIGMA ... / C ..."
    name = None
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "CASE":
            name = p[1]
        elif p[0] == "SIGMA":
            out.setdefault(name, {})["sigma"] = np.array([float(x) for x in p[1:]])
        elif p[0] == "C":
            out[name]["c"] = np.array([float(x) for x in p[1:]]).reshape(6, 6)
    return out


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_kernel_matches_oracle():
    cpp = _run_cpp()
    assert set(cpp) == set(CASES), f"case mismatch: {set(cpp)} vs {set(CASES)}"
    for name, F in CASES.items():
        sig6, c6 = _oracle(F)
        np.testing.assert_allclose(cpp[name]["sigma"], sig6, rtol=1e-10, atol=1e-8,
                                   err_msg=f"{name}: Cauchy stress mismatch")
        np.testing.assert_allclose(cpp[name]["c"], c6, rtol=1e-9, atol=1e-6,
                                   err_msg=f"{name}: material spatial tangent mismatch")
