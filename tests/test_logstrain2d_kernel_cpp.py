"""Cross-check the C++ 2D log-strain algorithm against the numpy oracle WITHOUT
building OpenSees. Compiles tests/logstrain2d_kernel_check.cpp (pure — includes
only SRC/material/nD/LogStrainKernel.h), runs it, and asserts its in-plane Cauchy
σ and full in-plane modulus c2 match logstrain2d_reference.py for both plane
strain and plane stress.

This catches C++ transcription bugs in LogStrain2D's plane logic — the F₃₃-Newton
plane-stress solve and the static tangent condensation — the one genuinely new
piece of math beyond the (already kernel-tested) 3D adaptor. Skipped if no C++
compiler is present.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import logstrain2d_reference as l2

K, G = 1500.0, 700.0
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))

# F cases — must match logstrain2d_kernel_check.cpp exactly (row-major 2×2).
_CASES = {
    "ps_generic":   (False, np.array([[1.15, 0.08], [-0.04, 0.93]])),
    "ps_shear":     (False, np.array([[1.0, 0.20], [0.0, 1.0]])),
    "ps_equibi":    (False, np.array([[1.10, 0.0], [0.0, 1.10]])),
    "pss_generic":  (True,  np.array([[1.15, 0.08], [-0.04, 0.93]])),
    "pss_shear":    (True,  np.array([[1.0, 0.20], [0.0, 1.0]])),
    "pss_equibi":   (True,  np.array([[1.10, 0.0], [0.0, 1.10]])),
    "pss_compress": (True,  np.array([[0.90, 0.05], [0.03, 0.88]])),
}


def _oracle(stress, F2):
    I3 = np.eye(3)
    step = l2.plane_stress_step if stress else l2.plane_strain_step
    res = step(F2, I3.copy(), I3.copy(), "elastic", K=K, G=G)
    return res["sigma_voigt"], res["c2"].reshape(16)


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "logstrain2d_kernel_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lsk2.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, src, "-o", exe],
                            capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    out, name = {}, None
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "CASE":
            name = p[1]
        elif p[0] == "SIGMA":
            out.setdefault(name, {})["sigma"] = np.array([float(x) for x in p[1:]])
        elif p[0] == "C2":
            out[name]["c2"] = np.array([float(x) for x in p[1:]])
    return out


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_2d_kernel_matches_oracle():
    cpp = _run_cpp()
    assert set(cpp) == set(_CASES), f"case mismatch: {set(cpp)} vs {set(_CASES)}"
    for name, (stress, F2) in _CASES.items():
        sig, c2 = _oracle(stress, F2)
        np.testing.assert_allclose(cpp[name]["sigma"], sig, rtol=1e-10, atol=1e-8,
                                   err_msg=f"{name}: in-plane Cauchy σ mismatch")
        np.testing.assert_allclose(cpp[name]["c2"], c2, rtol=1e-9, atol=1e-6,
                                   err_msg=f"{name}: in-plane modulus c2 mismatch")
