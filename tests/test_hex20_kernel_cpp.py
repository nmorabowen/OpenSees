"""Cross-check the pure C++ H20 serendipity kernel against the sympy/numpy oracle
WITHOUT building OpenSees. Compiles tests/hex20_kernel_check.cpp (includes only
the pure headers LadrunoHex20Shape.h + LadrunoMassLumping.h), runs it, and
asserts every printed quantity matches hex20_reference.py: the node/GP tables,
the serendipity N/dN, the distorted-hex detJ/B/K (27-pt & 8-pt), and the
reference-cube consistent mass, row-sum lump, HRZ lump, and volumes.

This catches any transcription bug in the HAND-CODED C++ serendipity algebra
(the header hand-differentiates N; the oracle derives via sympy.diff) — the
ADR-72 P0 novel-math surface — plus the row-sum negativity / HRZ positivity
facts that justify HRZ-only lumping. Skipped if no C++ compiler is present.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import hex20_reference as ref

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_BRICK = os.path.normpath(os.path.join(_HERE, "..", "SRC", "element", "ladrunoBrick"))
_SRC_INTEG = os.path.normpath(os.path.join(_HERE, "..", "SRC", "analysis", "integrator"))


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "hex20_kernel_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "hex20.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17",
                             "-I", _SRC_BRICK, "-I", _SRC_INTEG,
                             "-DLADRUNO_MASSLUMPING_STANDALONE",
                             src, "-o", exe], capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    return _parse(rr.stdout)


def _parse(text):
    nodexi = np.zeros((20, 3))
    gp27 = np.zeros((27, 4))
    gp8 = np.zeros((8, 4))
    N = np.zeros((5, 20))
    dN = np.zeros((5, 20, 3))
    detj = np.zeros(5)
    B = {0: np.zeros((6, 60)), 1: np.zeros((6, 60))}
    k27 = np.zeros((60, 60))
    k8 = np.zeros((60, 60))
    m = np.zeros((60, 60))
    rowsum = np.zeros(60)
    hrz = np.zeros(60)
    vol = {}
    for line in text.splitlines():
        p = line.split()
        if not p:
            continue
        tag = p[0]
        if tag == "NODEXI":
            nodexi[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "GP27":
            gp27[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "GP8":
            gp8[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "N":
            N[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "DN":
            dN[int(p[1]), int(p[2]), int(p[3])] = float(p[4])
        elif tag == "DETJ":
            detj[int(p[1])] = float(p[2])
        elif tag == "B":
            B[int(p[1])][int(p[2]), int(p[3])] = float(p[4])
        elif tag == "K27":
            k27[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "K8":
            k8[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "M":
            m[int(p[1]), int(p[2])] = float(p[3])
        elif tag == "ROWSUM":
            rowsum[int(p[1])] = float(p[2])
        elif tag == "HRZ":
            hrz[int(p[1])] = float(p[2])
        elif tag == "VOL27":
            vol[27] = float(p[1])
        elif tag == "VOL8":
            vol[8] = float(p[1])
    return dict(nodexi=nodexi, gp27=gp27, gp8=gp8, N=N, dN=dN, detj=detj, B=B,
                k27=k27, k8=k8, m=m, rowsum=rowsum, hrz=hrz, vol=vol)


@pytest.fixture(scope="module")
def cpp():
    return _run_cpp()


@pytest.mark.zone_a
@pytest.mark.t1
def test_node_and_gp_tables(cpp):
    # node + both GP tables match the oracle's NODE_XI / GP27_F / GP8_F exactly.
    np.testing.assert_allclose(cpp["nodexi"], np.array(ref.NODE_XI, dtype=float),
                               atol=1e-12, err_msg="NODE_XI mismatch")
    np.testing.assert_allclose(cpp["gp27"], ref.GP27_F, atol=1e-12,
                               err_msg="GP27 table mismatch")
    np.testing.assert_allclose(cpp["gp8"], ref.GP8_F, atol=1e-12,
                               err_msg="GP8 table mismatch")


@pytest.mark.zone_a
@pytest.mark.t1
def test_shape_functions(cpp):
    # N and dN at the 5 shared samples: hand-coded C++ vs sympy-derived oracle.
    for s, smp in enumerate(ref.XI_SAMPLES):
        No, dNo = ref.shape_numeric(smp)
        np.testing.assert_allclose(cpp["N"][s], No, atol=1e-12,
                                   err_msg=f"N mismatch @ sample {s}")
        np.testing.assert_allclose(cpp["dN"][s], dNo, atol=1e-12,
                                   err_msg=f"dN mismatch @ sample {s}")


@pytest.mark.zone_a
@pytest.mark.t1
def test_jacobian_and_B(cpp):
    # detJ at all 5 samples and B at samples 0 & 1 on the distorted hex.
    for s, smp in enumerate(ref.XI_SAMPLES):
        _, dNo = ref.shape_numeric(smp)
        _, detJo = ref.jacobian(ref.X_DISTORT, dNo)
        assert abs(cpp["detj"][s] - detJo) < 1e-12, f"detJ mismatch @ sample {s}"
        assert cpp["detj"][s] > 0.0, f"detJ not positive @ sample {s}"
    for s in (0, 1):
        _, dNo = ref.shape_numeric(ref.XI_SAMPLES[s])
        dNdx, _ = ref.cart_grad(ref.X_DISTORT, dNo)
        Bo = ref.fill_B(dNdx)
        np.testing.assert_allclose(cpp["B"][s], Bo, atol=1e-12,
                                   err_msg=f"B mismatch @ sample {s}")


@pytest.mark.zone_a
@pytest.mark.t1
def test_stiffness(cpp):
    # elastic K on the distorted hex, 27-pt and 8-pt, against the oracle.
    D = ref.iso_D(100.0, 0.25)
    K27o = ref.stiffness(ref.X_DISTORT, D, ref.GP27_F)
    K8o = ref.stiffness(ref.X_DISTORT, D, ref.GP8_F)
    np.testing.assert_allclose(cpp["k27"], K27o, rtol=1e-11, atol=1e-9,
                               err_msg="K27 mismatch")
    np.testing.assert_allclose(cpp["k8"], K8o, rtol=1e-11, atol=1e-9,
                               err_msg="K8 mismatch")


@pytest.mark.zone_a
@pytest.mark.t1
def test_consistent_mass(cpp):
    Mo = ref.consistent_mass(ref.X_CUBE, 1.0, ref.GP27_F)
    np.testing.assert_allclose(cpp["m"], Mo, rtol=1e-11, atol=1e-12,
                               err_msg="consistent mass mismatch")
    assert (np.diag(cpp["m"]) > 0).all(), "consistent-mass diagonal not positive"


@pytest.mark.zone_a
@pytest.mark.t1
def test_rowsum_negative_corners(cpp):
    # row-sum lump == M.sum(axis=1); corner (nodes 0..7) entries are -M/8 < 0 —
    # the FORBIDDEN lump that justifies HRZ-only (ADR 72 §3.5).
    Mo = ref.consistent_mass(ref.X_CUBE, 1.0, ref.GP27_F)
    np.testing.assert_allclose(cpp["rowsum"], Mo.sum(axis=1), atol=1e-12,
                               err_msg="row-sum lump mismatch")
    # total directional mass on the unit-density biunit cube = rho * V = 8.
    for a in range(8):                      # 8 corner nodes
        for d in range(3):
            v = cpp["rowsum"][3 * a + d]
            assert v < 0.0, f"row-sum corner dof {3*a+d} not negative: {v}"
            assert abs(v - (-8.0 / 8.0)) < 1e-12, f"row-sum corner != -M/8: {v}"


@pytest.mark.zone_a
@pytest.mark.t1
def test_hrz_positive_fractions(cpp):
    # HRZ (Ladruno::hrzLumpRaw on the 27-pt consistent mass) == oracle hrz_lump,
    # all POSITIVE, corner == 7/248 & mid-edge == 2/31 of total mass (ADR 72 §3.5).
    Mo = ref.consistent_mass(ref.X_CUBE, 1.0, ref.GP27_F)
    hrz_o = ref.hrz_lump(Mo)
    np.testing.assert_allclose(cpp["hrz"], hrz_o, rtol=1e-11, atol=1e-12,
                               err_msg="HRZ lump mismatch")
    assert (cpp["hrz"] > 0).all(), "HRZ produced a non-positive mass"
    # node 0 = corner, node 8 (dof 24) = first mid-edge; total dir mass = 8.
    for a in range(8):
        assert abs(cpp["hrz"][3 * a] - 8.0 * 7.0 / 248.0) < 1e-11, "HRZ corner != 7/248"
    for a in range(8, 20):
        assert abs(cpp["hrz"][3 * a] - 8.0 * 2.0 / 31.0) < 1e-11, "HRZ mid-edge != 2/31"
    # mass conservation: HRZ per-direction sum == total directional mass.
    assert abs(cpp["hrz"][0::3].sum() - 8.0) < 1e-10, "HRZ mass not conserved"


@pytest.mark.zone_a
@pytest.mark.t1
def test_volumes(cpp):
    # straight-edged geometry: both rules integrate the cube volume exactly.
    assert abs(cpp["vol"][27] - cpp["vol"][8]) < 1e-12, "VOL27 != VOL8"
    assert abs(cpp["vol"][27] - 8.0) < 1e-12, "cube volume != 8 @ 27-pt"
    assert abs(cpp["vol"][8] - 8.0) < 1e-12, "cube volume != 8 @ 8-pt"
