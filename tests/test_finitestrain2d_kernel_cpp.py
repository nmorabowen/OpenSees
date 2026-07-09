"""Cross-check the C++ shared 2D finite-strain kernel against the numpy oracle
WITHOUT building OpenSees. Compiles tests/finitestrain2d_kernel_check.cpp (pure —
includes only the two kernel headers), runs it, and asserts its element internal
force f and stiffness K match finitestrain2d_reference.py for a Q4 (std + F-bar)
and a T3.

This catches C++ transcription bugs in the finite-strain element math — the F/F⁻¹
push-forward, the a=c−σδ tangent contraction, and the unsymmetric F-bar G₀ coupling
— the ADR-70 P0 novel-math surface. Skipped if no C++ compiler is present.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import finitestrain2d_reference as fs

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))
_SRC_PLANE = os.path.normpath(os.path.join(_HERE, "..", "SRC", "element", "ladrunoPlane"))
K, G = 1500.0, 700.0

# Fixed ref coords + displacements — MUST match finitestrain2d_kernel_check.cpp.
_XQ4 = np.array([[0.0, 0.0], [2.0, 0.1], [2.2, 1.8], [0.1, 1.9]])
_XT3 = np.array([[0.0, 0.0], [2.0, 0.2], [0.3, 1.6]])
_XT6 = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [1.0, -0.05], [1.05, 1.0], [-0.05, 1.0]])
_UQ4 = np.array([0.010, -0.006, 0.004, 0.012, -0.008, 0.003, 0.007, -0.011])
_UT3 = np.array([0.009, -0.005, 0.006, 0.008, -0.004, 0.010])
_UT6 = np.array([0.008, -0.005, 0.006, 0.009, -0.004, 0.007,
                 0.003, -0.006, 0.005, 0.004, -0.007, 0.008])

_CASES = {
    "Q4_std":  ("Q4", _XQ4, _UQ4, False),
    "Q4_fbar": ("Q4", _XQ4, _UQ4, True),
    "T3_std":  ("T3", _XT3, _UT3, False),
    "T6_std":  ("T6", _XT6, _UT6, False),
    "T6_fbar": ("T6", _XT6, _UT6, True),
}


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "finitestrain2d_kernel_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "fs2d.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, "-I", _SRC_PLANE,
                             src, "-o", exe], capture_output=True, text=True)
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
        elif p[0] == "F":
            out.setdefault(name, {})["f"] = np.array([float(x) for x in p[1:]])
        elif p[0] == "K":
            out[name]["K"] = np.array([float(x) for x in p[1:]])
    return out


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_finitestrain2d_kernel_matches_oracle():
    cpp = _run_cpp()
    assert set(cpp) == set(_CASES), f"case mismatch: {set(cpp)} vs {set(_CASES)}"
    for name, (elem, X, u, fbar) in _CASES.items():
        f_o, K_o = fs.assemble(elem, u, K, G, fbar=fbar, X=X)
        nd = f_o.size
        np.testing.assert_allclose(cpp[name]["f"], f_o, rtol=1e-9, atol=1e-9,
                                   err_msg=f"{name}: internal force mismatch")
        np.testing.assert_allclose(cpp[name]["K"].reshape(nd, nd), K_o,
                                   rtol=1e-8, atol=1e-6,
                                   err_msg=f"{name}: element stiffness mismatch")
