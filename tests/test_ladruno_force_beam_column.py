"""Zone-A test suite for `LadrunoForceBeamColumn` (`LadrunoForceBeamColumn2d` &
`LadrunoForceBeamColumn3d`, classTags 33024 & 33025), plus cross-element
validation of `LadrunoDispBeamColumn` (`-spanLoad` quartic bubble enhancement)
and `LadrunoIMKBeam` (`-shearHinge` and `-hingeBiaxial`).
"""
from __future__ import annotations

import hashlib
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "dist" / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import opensees as ops  # noqa: E402


# ===========================================================================
# 1. Reduce-to-forceBeamColumn Parity (2D & 3D, Elastic & Inelastic + Span Loads)
# ===========================================================================

def test_ladruno_force_beam_column_2d_parity_with_forceBeamColumn():
    """Without hinges or shear flexibility, `LadrunoForceBeamColumn2d` reproduces
    `forceBeamColumn` bit-for-bit (< 1e-13) under combined nodal + uniform +
    partial-trapezoidal + point span loads on an inelastic fiber section.
    """
    L = 5.0
    E, Fy, b = 2.0e11, 350.0e6, 0.02

    def run_2d(ele_type: str):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)  # propped cantilever

        ops.uniaxialMaterial("Steel01", 1, Fy, E, b)
        ops.section("Fiber", 1)
        for y in (-0.18, -0.12, -0.06, 0.06, 0.12, 0.18):
            ops.fiber(y, 0.0, 0.0025, 1)

        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_type, 1, 1, 2, 1, 1, "-mass", 120.0, "-cMass")

        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -40.0e3, 0.0, 15.0e3)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -18.0e3, 3.0e3)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -10.0e3, 1.0e3, 0.2, 0.8, -22.0e3, 4.0e3)
        ops.eleLoad("-ele", 1, "-type", "-beamPoint", -25.0e3, 0.5, 5.0e3)

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(4):
            assert ops.analyze(1) == 0

        disp2 = [ops.nodeDisp(2, d) for d in (1, 2, 3)]
        forces = list(ops.eleForce(1))
        sec_f = [list(ops.eleResponse(1, "section", i, "force")) for i in range(1, 6)]
        sec_e = [list(ops.eleResponse(1, "section", i, "deformation")) for i in range(1, 6)]
        mass = list(ops.eleResponse(1, "stiff"))
        return disp2, forces, sec_f, sec_e, mass

    d_ref, f_ref, sf_ref, se_ref, k_ref = run_2d("forceBeamColumn")
    d_lad, f_lad, sf_lad, se_lad, k_lad = run_2d("LadrunoForceBeamColumn")

    np.testing.assert_allclose(d_lad, d_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(f_lad, f_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(sf_lad, sf_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(se_lad, se_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(k_lad, k_ref, rtol=1e-13, atol=1e-13)


def test_ladruno_force_beam_column_3d_parity_with_forceBeamColumn():
    """Without hinges or shear flexibility, `LadrunoForceBeamColumn3d` reproduces
    `forceBeamColumn` bit-for-bit (< 1e-13) in 3D under biaxial + torsional +
    axial nodal and intra-span loads.
    """
    L = 4.5
    E, Fy, b = 2.0e11, 350.0e6, 0.02

    def run_3d(ele_type: str):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.fix(2, 0, 1, 1, 0, 0, 0)

        ops.uniaxialMaterial("Steel01", 1, Fy, E, b)
        ops.section("Fiber", 1, "-GJ", 2.5e6)
        for y in (-0.15, 0.15):
            for z in (-0.12, 0.12):
                ops.fiber(y, z, 0.003, 1)

        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_type, 1, 1, 2, 1, 1, "-mass", 150.0, "-cMass")

        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -30.0e3, 0.0, 0.0, 4.0e3, -10.0e3, 12.0e3)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -14.0e3, 8.0e3, 2.0e3)
        ops.eleLoad(
            "-ele", 1, "-type", "-beamUniform",
            -6.0e3, 4.0e3, 1.0e3, 0.25, 0.75, -12.0e3, 9.0e3, 3.0e3
        )
        ops.eleLoad("-ele", 1, "-type", "-beamPoint", -15.0e3, 10.0e3, 0.5, 4.0e3)

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(4):
            assert ops.analyze(1) == 0

        disp2 = [ops.nodeDisp(2, d) for d in range(1, 7)]
        forces = list(ops.eleForce(1))
        sec_f = [list(ops.eleResponse(1, "section", i, "force")) for i in range(1, 6)]
        sec_e = [list(ops.eleResponse(1, "section", i, "deformation")) for i in range(1, 6)]
        return disp2, forces, sec_f, sec_e

    d_ref, f_ref, sf_ref, se_ref = run_3d("forceBeamColumn")
    d_lad, f_lad, sf_lad, se_lad = run_3d("LadrunoForceBeamColumn")

    np.testing.assert_allclose(d_lad, d_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(f_lad, f_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(sf_lad, sf_ref, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(se_lad, se_ref, rtol=1e-13, atol=1e-13)


# ===========================================================================
# 2. Tier-2 Cohesive Hinges, Series Compliance, P-M Scaling & Shear Flexibility
# ===========================================================================

def test_ladruno_force_beam_column_2d_hinges_shear_pm_and_fd_tangent():
    """Verify `LadrunoForceBeamColumn2d` with:
      - `-hingeI`, `-hingeJ`, `-hingeIJ` series compliance and softening
      - `-shear GAy` flexibility addition `L/(GAy)` on `kvInit` and `f`
      - `-pm N0 alphaPM` parabolic yield-moment amplification under axial compression
      - Central finite-difference verification of `basicStiffness` (`responseID == 19`)
    """
    L = 4.0
    E, A, Iz = 2.0e11, 0.04, 1.6e-4
    GAy = 2.0e8
    Mc, Gf = 60.0e3, 400.0

    # 1. Elastic + Shear flexibility check: tip deflection of cantilever under tip load P
    # delta = P*L^3/(3*E*Iz) + P*L/GAy
    P_tip = -12.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1, "-shear", GAy)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, P_tip, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    uy_exact = P_tip * (L ** 3) / (3.0 * E * Iz) + P_tip * L / GAy
    assert ops.nodeDisp(2, 2) == pytest.approx(uy_exact, rel=1e-12)

    # 2. Dual end hinges (-hingeIJ) + -pm scaling + FD basicStiffness check
    def eval_2d_vq(ux, uy, rz, axial_p0=200.0e3, alpha_pm=0.8):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, Mc, Gf, "-linear")
        pm_args = ("-pm", axial_p0, alpha_pm) if axial_p0 > 0.0 else ()
        ops.element(
            "LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1,
            "-hingeIJ", 10, 10, *pm_args, "-shear", GAy
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, ux)
        ops.sp(2, 2, uy)
        ops.sp(2, 3, rz)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        v = np.array(ops.eleResponse(1, "basicDeformation"))
        q = np.array(ops.eleResponse(1, "basicForce"))
        kb = np.array(ops.eleResponse(1, "basicStiffness")).reshape(3, 3)
        hs = np.array(ops.eleResponse(1, "hingeSummary"))
        return v, q, kb, hs

    ux0, uy0, rz0 = -0.00005, -0.004, 0.014
    v0, q0, kb0, hs0 = eval_2d_vq(ux0, uy0, rz0, axial_p0=200.0e3, alpha_pm=0.8)
    v0_nopm, q0_nopm, kb0_nopm, hs0_nopm = eval_2d_vq(ux0, uy0, rz0, axial_p0=0.0, alpha_pm=0.0)
    assert abs(hs0[0]) > 0.0, "Hinge must be active"
    assert abs(hs0[1]) > abs(hs0_nopm[1]), "Axial compression with -pm must amplify hinge moment"
    # Verify basicStiffness matches central FD of basicForce w.r.t. basicDeformation
    eps = 1.0e-8
    DV = np.zeros((3, 3))
    DQ = np.zeros((3, 3))
    DV_nopm = np.zeros((3, 3))
    DQ_nopm = np.zeros((3, 3))
    for j in range(3):
        du = [0.0, 0.0, 0.0]
        du[j] = eps
        vp, qp, _, _ = eval_2d_vq(ux0 + du[0], uy0 + du[1], rz0 + du[2], axial_p0=200.0e3, alpha_pm=0.8)
        vm, qm, _, _ = eval_2d_vq(ux0 - du[0], uy0 - du[1], rz0 - du[2], axial_p0=200.0e3, alpha_pm=0.8)
        DV[:, j] = (vp - vm) / (2.0 * eps)
        DQ[:, j] = (qp - qm) / (2.0 * eps)
        vp0, qp0, _, _ = eval_2d_vq(ux0 + du[0], uy0 + du[1], rz0 + du[2], axial_p0=0.0, alpha_pm=0.0)
        vm0, qm0, _, _ = eval_2d_vq(ux0 - du[0], uy0 - du[1], rz0 - du[2], axial_p0=0.0, alpha_pm=0.0)
        DV_nopm[:, j] = (vp0 - vm0) / (2.0 * eps)
        DQ_nopm[:, j] = (qp0 - qm0) / (2.0 * eps)
    kb_fd = DQ @ np.linalg.inv(DV)
    kb_nopm_fd = DQ_nopm @ np.linalg.inv(DV_nopm)
    np.testing.assert_allclose(kb0, kb_fd, rtol=1e-5, atol=1e-1)
    np.testing.assert_allclose(kb0_nopm, kb_nopm_fd, rtol=1e-5, atol=1e-1)


def test_ladruno_force_beam_column_3d_biaxial_hinges_shear_and_fd():
    """Verify `LadrunoForceBeamColumn3d` with coupled `-hingeBiaxialIJ`,
    `-shear GAy GAz`, and central finite-difference `basicStiffness`.
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4
    GAy, GAz = 2.5e8, 2.0e8
    Mc_z, Gf_z = 60.0e3, 450.0
    Mc_y, Gf_y = 50.0e3, 450.0

    def eval_3d_vq(tip):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 88, Mc_z, Gf_z, Mc_y, Gf_y, "-linear")
        ops.element(
            "LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1,
            "-hingeBiaxialIJ", 88, 88, "-shear", GAy, GAz
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for dof, val in enumerate(tip, start=1):
            ops.sp(2, dof, val)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return (
            np.array(ops.eleResponse(1, "basicDeformation")),
            np.array(ops.eleResponse(1, "basicForce")),
            np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6),
            np.array(ops.eleResponse(1, "hingeSummary")),
        )

    base = [-0.00005, -0.005, 0.004, 0.0005, 0.012, 0.015]
    v0, q0, kb3d, hs3 = eval_3d_vq(base)
    assert abs(hs3[0]) > 0.0 and abs(hs3[4]) > 0.0, "Both biaxial hinge axes must be active"

    eps = 1.0e-8
    DV = np.zeros((6, 6))
    DQ = np.zeros((6, 6))
    for j in range(6):
        tp = list(base)
        tm = list(base)
        tp[j] += eps
        tm[j] -= eps
        vp, qp, _, _ = eval_3d_vq(tp)
        vm, qm, _, _ = eval_3d_vq(tm)
        DV[:, j] = (vp - vm) / (2.0 * eps)
        DQ[:, j] = (qp - qm) / (2.0 * eps)
    kb3d_fd = DQ @ np.linalg.inv(DV)
    np.testing.assert_allclose(kb3d, kb3d_fd, rtol=1e-5, atol=1e-1)


# ===========================================================================
# 3. LadrunoDispBeamColumn Exact Span-Load Internal Strain + Quartic Bubble
# ===========================================================================

def test_ladruno_disp_beam_column_span_load_bubble_2d_and_3d():
    """Verify `LadrunoDispBeamColumn2d` and `LadrunoDispBeamColumn3d` `-spanLoad`
    (default on):
      1. On a 1-element simply-supported beam (pin-roller, u_y = 0 at both ends)
         under a uniform span load w_y, standard Hermite interpolation (`-noSpanLoad`)
         produces ZERO midspan section curvature/moment when the pin-roller end
         rotations are constrained (fixed-fixed) and only half the exact midspan
         moment on a simply-supported span, whereas `-spanLoad` (default) reproduces
         the exact parabolic moment distribution M_z(x) = -w*x*(L-x)/2 (simply
         supported) and M_z(x) = w*L^2/12 - w*x*(L-x)/2 (fixed-fixed) at every
         Gauss point to machine precision (< 1e-11), matching `LadrunoForceBeamColumn`.
      2. On a 1-element inelastic fiber beam under uniform + point span loads,
         `-spanLoad` activates midspan fiber yielding and matches the force-based
         `LadrunoForceBeamColumn` response closely.
    """
    L = 6.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 2.0e-4, 1.5e-4, 8.0e10, 2.5e-4
    wy = -24.0e3
    wz = 18.0e3

    # --- 1. 2D Fixed-Fixed 1-element beam under uniform load wy ---
    # Under fixed-fixed BCs, nodal v = [0, 0, 0]!
    # With -noSpanLoad, section deformations e(x) = B(x)*v == 0 at all Gauss points!
    # With -spanLoad (default), section moments at all 5 Lobatto points equal exact
    # M(x) = w*L^2/12 * (1 - 6*xi + 6*xi^2), matching LadrunoForceBeamColumn2d to 1e-12!
    def run_2d_fixed_fixed(ele_cmd: str, extra_args=()):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_cmd, 1, 1, 2, 1, 1, *extra_args)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        sec_m = [ops.eleResponse(1, "section", i, "force")[1] for i in range(1, 6)]
        end_f = list(ops.eleForce(1))
        return sec_m, end_f

    m_no_span, _ = run_2d_fixed_fixed("LadrunoDispBeamColumn2d", ("-noSpanLoad",))
    assert np.allclose(m_no_span, 0.0, atol=1e-12), "Without -spanLoad, fixed-fixed 1-ele beam is blind to span load!"

    m_span, f_span = run_2d_fixed_fixed("LadrunoDispBeamColumn2d")
    m_fbc, f_fbc = run_2d_fixed_fixed("LadrunoForceBeamColumn2d")
    np.testing.assert_allclose(m_span, m_fbc, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(f_span, f_fbc, rtol=1e-12, atol=1e-9)
    # Check midspan (IP 3, xi = 0.5): M(L/2) = -wy * L^2 / 24 = +36.0e3 N*m
    assert m_span[2] == pytest.approx(-wy * (L ** 2) / 24.0, rel=1e-12)

    # --- 2. 2D Simply-Supported 1-element beam under uniform load wy ---
    def run_2d_simply_supported(ele_cmd: str, extra_args=()):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 0)
        ops.fix(2, 0, 1, 0)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_cmd, 1, 1, 2, 1, 1, *extra_args)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        sec_m = [ops.eleResponse(1, "section", i, "force")[1] for i in range(1, 6)]
        end_f = list(ops.eleForce(1))
        rot = [ops.nodeDisp(1, 3), ops.nodeDisp(2, 3)]
        return sec_m, end_f, rot

    m_ss_span, f_ss_span, r_ss_span = run_2d_simply_supported("LadrunoDispBeamColumn2d")
    m_ss_fbc, f_ss_fbc, r_ss_fbc = run_2d_simply_supported("LadrunoForceBeamColumn2d")
    np.testing.assert_allclose(m_ss_span, m_ss_fbc, rtol=1e-12, atol=1e-8)
    np.testing.assert_allclose(f_ss_span, f_ss_fbc, rtol=1e-12, atol=1e-8)
    np.testing.assert_allclose(r_ss_span, r_ss_fbc, rtol=1e-12, atol=1e-12)
    # Midspan moment M(L/2) = -wy * L^2 / 8 = +108.0e3 N*m, end moments = 0
    assert m_ss_span[0] == pytest.approx(0.0, abs=1e-8)
    assert m_ss_span[2] == pytest.approx(-wy * (L ** 2) / 8.0, rel=1e-12)
    assert m_ss_span[4] == pytest.approx(0.0, abs=1e-8)

    # --- 3. 2D Nonlinear Fiber Simply-Supported 1-element beam under uniform span load ---
    def run_2d_nonlinear_ss(ele_cmd: str, extra_args=()):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 0)
        ops.fix(2, 0, 1, 0)
        ops.uniaxialMaterial("Steel01", 1, 135.0e6, E, 0.02)
        ops.section("Fiber", 1)
        for y in (-0.15, -0.10, -0.05, 0.05, 0.10, 0.15):
            ops.fiber(y, 0.0, 0.002, 1)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_cmd, 1, 1, 2, 1, 1, *extra_args)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -32.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.2)
        ops.test("NormDispIncr", 1.0e-12, 30)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(5):
            assert ops.analyze(1) == 0
        mid_m = ops.eleResponse(1, "section", 3, "force")[1]
        mid_k = ops.eleResponse(1, "section", 3, "deformation")[1]
        rot = [ops.nodeDisp(1, 3), ops.nodeDisp(2, 3)]
        return mid_m, mid_k, rot

    mm_disp, mk_disp, r_disp = run_2d_nonlinear_ss("LadrunoDispBeamColumn2d")
    mm_fbc, mk_fbc, r_fbc = run_2d_nonlinear_ss("LadrunoForceBeamColumn2d")
    # Midspan moment on a single-element simply-supported beam equals w*L^2/8 = 144 kN*m
    assert mm_disp == pytest.approx(mm_fbc, rel=4e-2)
    assert abs(mk_disp) > abs(mm_disp) / (E * 1.4e-4) * 1.15, "Midspan section must have yielded inelastically"
    assert abs(mk_fbc) > abs(mm_fbc) / (E * 1.4e-4) * 1.15, "Midspan section must have yielded inelastically"

    # --- 4. 3D Simply-Supported 1-element beam under biaxial uniform span loads wy, wz ---
    def run_3d_ss(ele_cmd: str):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 0, 0)
        ops.fix(2, 0, 1, 1, 1, 0, 0)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_cmd, 1, 1, 2, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy, wz)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        sec_f = [list(ops.eleResponse(1, "section", i, "force")) for i in range(1, 6)]
        end_f = list(ops.eleForce(1))
        return sec_f, end_f

    sf3_disp, ef3_disp = run_3d_ss("LadrunoDispBeamColumn3d")
    sf3_fbc, ef3_fbc = run_3d_ss("LadrunoForceBeamColumn3d")
    np.testing.assert_allclose(sf3_disp, sf3_fbc, rtol=1e-12, atol=1e-8)
    np.testing.assert_allclose(ef3_disp, ef3_fbc, rtol=1e-12, atol=1e-8)


# ===========================================================================
# 4. LadrunoIMKBeam Series Shear Hinges & 3D Coupled Biaxial Hinges
# ===========================================================================

def test_ladruno_imk_beam_2d_and_3d_shear_hinges_and_biaxial_hinges():
    """Verify:
      1. `LadrunoIMKBeam2d` with `-shearHinge $matTag`: total shear compliance is
         `L^2/(12*EI) + 1/Ks_hinge`, and yielding of the shear hinge caps the
         element shear force V_y at Vy_yield while `shearHinge` response tokens
         report `[delta_v, V_y]`.
      2. `LadrunoIMKBeam` (3D) with `-shearHingeY $matY -shearHingeZ $matZ` and
         `-hingeBiaxial $ndMatTag` (coupled Mz-My `NDMaterial` hinge): biaxial
         hinge interaction reduces the uniaxial yield moment under simultaneous
         bending about y and z, and `basicStiffness` matches central finite
         differences.
    """
    L = 3.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # --- 1. 2D LadrunoIMKBeam2d with -shearHinge ---
    Vy_yield = 40.0e3
    Ks0 = 5.0e7
    b_sh = 0.01
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 0, 1)  # shear building column (fixed rotation at top & bottom)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 10, 200.0e3, 1.0e9, 0.02)
    ops.uniaxialMaterial("Steel01", 20, Vy_yield, Ks0, b_sh)
    ops.element(
        "LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-shearHinge", 20
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 2, 0.005)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.2)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    for _ in range(5):
        assert ops.analyze(1) == 0

    sh_resp = ops.eleResponse(1, "shearHinge")
    assert len(sh_resp) == 2
    dv_sh, vy_sh = sh_resp
    assert abs(dv_sh) > Vy_yield / Ks0, "Shear hinge must have yielded"
    # Shear force in element endForces[1] must equal vy_sh (up to sign)
    ef = ops.eleForce(1)
    assert abs(ef[1]) == pytest.approx(abs(vy_sh), rel=1e-10)
    assert abs(vy_sh) > Vy_yield and abs(vy_sh) < Vy_yield * 1.15

    # --- 2. 3D LadrunoIMKBeam with -hingeBiaxial and -shearHingeY / -shearHingeZ ---
    Mc_z, Gf_z = 80.0e3, 500.0
    Mc_y, Gf_y = 60.0e3, 500.0

    def eval_imk3d_vq(tip):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 50, Mc_z, Gf_z, Mc_y, Gf_y, "-linear")
        ops.uniaxialMaterial("Steel01", 60, 50.0e3, 8.0e7, 0.05)
        ops.uniaxialMaterial("Steel01", 61, 45.0e3, 7.0e7, 0.05)
        ops.element(
            "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
            "-hingeBiaxial", 50, "-shearHingeY", 60, "-shearHingeZ", 61
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for dof, val in enumerate(tip, start=1):
            ops.sp(2, dof, val)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return (
            np.array(ops.eleResponse(1, "basicDeformation")),
            np.array(ops.eleResponse(1, "basicForce")),
            np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6),
            np.array(ops.eleResponse(1, "hingeSummary")),
            np.array(ops.eleResponse(1, "shearHinge")),
        )

    base = [0.0001, -0.004, 0.003, 0.0005, 0.012, 0.015]
    v0, q0, kb3, hs3, sh3 = eval_imk3d_vq(base)
    assert len(sh3) == 4
    assert abs(hs3[0]) > 0.0 and abs(hs3[4]) > 0.0, "3D coupled biaxial hinges must be active"

    eps = 1.0e-8
    DV = np.zeros((6, 6))
    DQ = np.zeros((6, 6))
    for j in range(6):
        tp = list(base)
        tm = list(base)
        tp[j] += eps
        tm[j] -= eps
        vp, qp, _, _, _ = eval_imk3d_vq(tp)
        vm, qm, _, _, _ = eval_imk3d_vq(tm)
        DV[:, j] = (vp - vm) / (2.0 * eps)
        DQ[:, j] = (qp - qm) / (2.0 * eps)
    kb3_fd = DQ @ np.linalg.inv(DV)
    np.testing.assert_allclose(kb3, kb3_fd, rtol=1e-5, atol=1e-1)


# ===========================================================================
# 5. MPCO/HDF5 Recorder, Database Round-Trip, Tcl Subprocess, & OpenMP Threads
# ===========================================================================

def test_ladruno_force_beam_column_recorder_db_tcl_and_openmp(tmp_path):
    """Verify:
      - `.ladruno` HDF5 recorder writes `endForces` and `hingeState` for
        `LadrunoForceBeamColumn2d` and `LadrunoForceBeamColumn3d`
      - `database("File", ...)` `save`/`restore` round-trip preserves state
      - `OpenSees.exe` classic-Tcl subprocess parses and runs `LadrunoForceBeamColumn`
      - `ladrunoThreads(1/2/4)` produces bit-identical results
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # 1. Database save/restore + HDF5 recorder in 2D
    h5_2d = str(tmp_path / "fbc2d.ladruno")
    db_2d = str(tmp_path / "fbc2d_db")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 50.0e3, 300.0, "-linear")
    ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-shear", 2.0e8)
    ops.recorder("ladruno", h5_2d, "-E", "endForces", "hingeState", "-flush", 0, "-compress", 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -20.0e3, -15.0e3, 10.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0

    ef_before = list(ops.eleForce(1))
    hs_before = list(ops.eleResponse(1, "hingeSummary"))
    ops.database("File", db_2d)
    ops.save(1)
    ops.reset()
    ops.restore(1)
    np.testing.assert_allclose(ops.eleForce(1), ef_before, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(ops.eleResponse(1, "hingeSummary"), hs_before, rtol=1e-13, atol=1e-13)
    ops.wipe()
    assert Path(h5_2d).stat().st_size > 0

    # 2. Classic-Tcl subprocess smoke via OpenSees.exe
    opensees_exe = BIN / "OpenSees.exe"
    if opensees_exe.exists():
        tcl_script = tmp_path / "fbc_smoke.tcl"
        tcl_script.write_text(
            """
            model basic -ndm 2 -ndf 3
            node 1 0.0 0.0
            node 2 4.0 0.0
            fix 1 1 1 1
            geomTransf Linear 1
            section Elastic 1 2.0e11 0.04 1.6e-4
            beamIntegration Lobatto 1 1 4
            uniaxialMaterial LadrunoCohesiveHinge 77 50.0e3 300.0 -linear
            element LadrunoForceBeamColumn 1 1 2 1 1 -hinge 77 -shear 2.0e8
            timeSeries Linear 1
            pattern Plain 1 1 {
                load 2 -10.0e3 -15.0e3 10.0e3
            }
            system BandGeneral
            numberer Plain
            constraints Plain
            integrator LoadControl 1.0
            test NormDispIncr 1.0e-12 25
            algorithm Newton
            analysis Static
            if {[analyze 1] != 0} { exit 2 }
            puts "DISP2=[nodeDisp 2 2]"
            exit 0
            """,
            encoding="utf-8",
        )
        res = subprocess.run([str(opensees_exe), str(tcl_script)], capture_output=True, text=True, timeout=30)
        assert res.returncode == 0, f"OpenSees.exe failed: {res.stderr}\n{res.stdout}"
        assert "DISP2=" in (res.stdout + res.stderr)

    # 3. OpenMP multi-threaded bit-identity across 1, 2, 4 threads
    def run_omp(nt: int):
        ops.wipe()
        ops.ladrunoThreads(nt)
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.geomTransf("Linear", 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 50.0e3, 300.0, "-linear")
        for i in range(8):
            n1, n2 = 2 * i + 1, 2 * i + 2
            ops.node(n1, 0.0, float(i))
            ops.node(n2, L, float(i))
            ops.fix(n1, 1, 1, 1)
            ops.element("LadrunoForceBeamColumn", i + 1, n1, n2, 1, 1, "-hinge", 77, "-shear", 2.0e8)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for i in range(8):
            ops.sp(2 * i + 2, 1, -0.0001)
            ops.sp(2 * i + 2, 2, -0.006)
            ops.sp(2 * i + 2, 3, 0.004)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0
        h = hashlib.md5()
        for i in range(8):
            for v in ops.eleForce(i + 1):
                h.update(f"{v!r};".encode())
        ops.ladrunoThreads(1)
        return h.hexdigest()

    ref_digest = run_omp(1)
    assert run_omp(2) == ref_digest


# ===========================================================================
# 6. Regression Tests for Reviewer-Identified Edge Cases & Fixes
# ===========================================================================

def test_ladruno_force_beam_column_3d_pm_compression_and_tension_scaling():
    """Verify `LadrunoForceBeamColumn3d` `-pm P0 alpha` scales hinge moment capacity
    by `f(P) = clamp(1 - alpha*(P/P0) - (1-alpha)*(P/P0)^2, 0.05, 2.5)` where
    compression (`P < 0`) amplifies moment capacity below balanced load and
    tension (`P > 0`) reduces moment capacity, matching `LadrunoForceBeamColumn2d`
    and `LadrunoDispBeamColumn3d`.
    """
    L = 3.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4
    P0, alpha = 200.0e3, 0.8

    def run_3d_pm(P_axial: float):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("Steel01", 10, 60.0e3, 1.0e8, 0.01)
        ops.element("LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1, "-hinge", 10, "-pm", P0, alpha)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, P_axial, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.sp(2, 2, -0.04)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 0.2)
        ops.test("NormDispIncr", 1.0e-12, 30)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(5):
            assert ops.analyze(1) == 0
        return ops.eleResponse(1, "hingeSummary")[1]

    m_zero = abs(run_3d_pm(0.0))
    m_comp = abs(run_3d_pm(-50.0e3))
    m_tens = abs(run_3d_pm(+50.0e3))
    # Expected f(-50e3) = 1 - 0.8*(-0.25) - 0.2*(0.0625) = 1.1875
    # Expected f(+50e3) = 1 - 0.8*(+0.25) - 0.2*(0.0625) = 0.7875
    assert m_comp == pytest.approx(m_zero * 1.1875, rel=1e-2)
    assert m_tens == pytest.approx(m_zero * 0.7875, rel=1e-2)


def test_ladruno_force_beam_column_3d_implex_section_error_parameter_and_reset():
    """Verify:
      1. `LadrunoForceBeamColumn3d` with `-implex` (and no hinges) reports positive
         `implexError` / `avgImplexError` from section-level IMPLEX drift.
      2. `parameter` + `updateParameter` on `GAy`, `lch`, `section`, and `hinge`
         invalidates `Ki` and updates `basicStiffness`.
      3. `ops.reset()` (`revertToStart`) restores virgin elastic stiffness after
         inelastic loading.
    """
    L = 3.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # 1. Section-level -implex error on 3D LadrunoForceBeamColumn3d without hinges
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.uniaxialMaterial("Steel01", 1, 250.0e6, E, 0.01)
    ops.section("Fiber", 1, "-GJ", G * Jx)
    for y in (-0.1, 0.1):
        for z in (-0.1, 0.1):
            ops.fiber(y, z, 0.005, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.element("LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1, "-implex")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 2, -0.05)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    assert ops.eleResponse(1, "implexError")[0] > 0.0
    assert ops.eleResponse(1, "avgImplexError")[0] > 0.0

    # 2. Parameter updates and revertToStart on 3D LadrunoForceBeamColumn3d
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 50.0e3, 400.0, "-linear")
    ops.element("LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1, "-hinge", 10, "-shear", 1.0e8, 1.0e8, "-lch", 0.5)
    kb0 = np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6)
    ops.parameter(1, "element", 1, "GAy")
    ops.parameter(2, "element", 1, "lch")
    ops.updateParameter(1, 2.0e7)
    ops.updateParameter(2, 0.25)
    kb1 = np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6)
    assert kb1[1, 1] < kb0[1, 1], "Lowering GAy via updateParameter must reduce basic flexural stiffness"

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 2, -0.025)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    kb_yielded = np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6)
    assert kb_yielded[1, 1] < kb1[1, 1] * 0.5
    ops.reset()
    kb_reset = np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6)
    np.testing.assert_allclose(kb_reset, kb1, rtol=1e-10, atol=1e-6)


def test_ladruno_disp_beam_column_2pt_legendre_and_unsymmetric_rc_tbeam_span_load():
    """Verify:
      1. `LadrunoDispBeamColumn3d` with 2-point Gauss-Legendre integration under
         `-beamUniform` does not suffer roundoff division-by-zero in the quartic
         bubble (since 1 - 6*xi + 6*xi^2 vanishes at the 2 Gauss points) and
         preserves exact elastic stiffness `4*E*I/L`.
      2. Single-element unsymmetric RC T-beam (`Concrete01` + `Steel01` fibers)
         subjected to axial compression N = -100 kN and heavy span load
         w = -150 kN/m past yielding in both 2D and 3D `LadrunoDispBeamColumn`
         uses the coupled axial+flexural bubble `(beta_P, beta_Z)` so midspan
         axial force stays near N = -100 kN (no artificial axial locking) and
         midspan curvature / end rotation closely match `LadrunoForceBeamColumn`.
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # 1. 3D 2-point Gauss-Legendre stiffness under uniform span load
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 0, 0)
    ops.fix(2, 0, 1, 1, 1, 0, 0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Legendre", 1, 1, 2)
    ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", -24.0e3, -12.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    kb_leg2 = np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6)
    assert kb_leg2[1, 1] == pytest.approx(4.0 * E * Iz / L, rel=1e-12)
    assert kb_leg2[2, 2] == pytest.approx(4.0 * E * Iz / L, rel=1e-12)
    assert kb_leg2[3, 3] == pytest.approx(4.0 * E * Iy / L, rel=1e-12)
    assert kb_leg2[4, 4] == pytest.approx(4.0 * E * Iy / L, rel=1e-12)

    # 2. Single-element 2D & 3D unsymmetric cracked RC T-beam under axial compression + span load
    def run_rc_tbeam_2d(ele_cmd: str):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 5.0, 0.0)
        ops.fix(1, 1, 1, 0)
        ops.fix(2, 0, 1, 0)
        ops.uniaxialMaterial("Concrete01", 1, -30.0e6, -0.002, -6.0e6, -0.006)
        ops.uniaxialMaterial("Steel01", 2, 420.0e6, 200.0e9, 0.01)
        ops.section("Fiber", 1)
        ops.patch("rect", 1, 4, 1, 0.15, -0.50, 0.25, 0.50)   # flange
        ops.patch("rect", 1, 8, 1, -0.25, -0.15, 0.15, 0.15)  # web
        ops.fiber(-0.21, 0.0, 0.0018, 2)                      # bottom tension steel
        ops.fiber(0.21, 0.0, 0.0006, 2)                       # top compression steel
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_cmd, 1, 1, 2, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0e3, 0.0, 0.0)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -100.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.2)
        ops.test("NormDispIncr", 1.0e-10, 40)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(5):
            assert ops.analyze(1) == 0
        mid_f = ops.eleResponse(1, "section", 3, "force")
        mid_d = ops.eleResponse(1, "section", 3, "deformation")
        return mid_f, mid_d, ops.nodeDisp(1, 3)

    mf_disp2, md_disp2, rot_disp2 = run_rc_tbeam_2d("LadrunoDispBeamColumn2d")
    mf_fbc2, md_fbc2, rot_fbc2 = run_rc_tbeam_2d("LadrunoForceBeamColumn2d")
    # Midspan axial force must remain near the applied -100 kN (not locked at -600 kN!)
    assert mf_disp2[0] == pytest.approx(-100.0e3, rel=5e-2)
    assert mf_disp2[1] == pytest.approx(mf_fbc2[1], rel=2e-2)
    assert md_disp2[1] == pytest.approx(md_fbc2[1], rel=2e-2)
    assert rot_disp2 == pytest.approx(rot_fbc2, rel=1e-2)

    def run_rc_tbeam_3d(ele_cmd: str):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, 5.0, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 0, 0)
        ops.fix(2, 0, 1, 1, 1, 0, 0)
        ops.uniaxialMaterial("Concrete01", 1, -30.0e6, -0.002, -6.0e6, -0.006)
        ops.uniaxialMaterial("Steel01", 2, 420.0e6, 200.0e9, 0.01)
        ops.section("Fiber", 1, "-GJ", 1.0e8)
        ops.patch("rect", 1, 4, 2, 0.15, -0.50, 0.25, 0.50)   # flange
        ops.patch("rect", 1, 8, 2, -0.25, -0.15, 0.15, 0.15)  # web
        ops.fiber(-0.21, 0.0, 0.0018, 2)                      # bottom tension steel
        ops.fiber(0.21, 0.0, 0.0006, 2)                       # top compression steel
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_cmd, 1, 1, 2, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -100.0e3, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.2)
        ops.test("NormDispIncr", 1.0e-10, 40)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(5):
            assert ops.analyze(1) == 0
        mid_f = ops.eleResponse(1, "section", 3, "force")
        mid_d = ops.eleResponse(1, "section", 3, "deformation")
        return mid_f, mid_d, ops.nodeDisp(1, 6)

    mf_disp3, md_disp3, rot_disp3 = run_rc_tbeam_3d("LadrunoDispBeamColumn3d")
    mf_fbc3, md_fbc3, rot_fbc3 = run_rc_tbeam_3d("LadrunoForceBeamColumn3d")
    assert mf_disp3[0] == pytest.approx(-100.0e3, rel=5e-2)
    assert mf_disp3[1] == pytest.approx(mf_fbc3[1], rel=2e-2)
    assert md_disp3[1] == pytest.approx(md_fbc3[1], rel=2e-2)
    assert rot_disp3 == pytest.approx(rot_fbc3, rel=1e-2)


def test_ladruno_imk_beam_3d_simultaneous_biaxial_and_shear_hinge_softening():
    """Verify 3D `LadrunoIMKBeam` with softening `-hingeBiaxial` (`LadrunoCohesiveHingeBiaxial`)
    and softening `-shearHingeY` / `-shearHingeZ` (`LadrunoCohesiveHinge`) pushed
    simultaneously into the post-peak softening regime converges and its condensed
    `basicStiffness` matches central finite differences.
    """
    L = 3.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    def eval_softening_imk3d(tip):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 50, 80.0e3, 800.0, 60.0e3, 800.0, "-linear")
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 60, 45.0e3, 300.0, "-linear")
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 61, 40.0e3, 300.0, "-linear")
        ops.element(
            "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
            "-hingeBiaxial", 50, "-shearHingeY", 60, "-shearHingeZ", 61
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for dof, val in enumerate(tip, start=1):
            ops.sp(2, dof, val)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 0.25)
        ops.test("NormDispIncr", 1.0e-12, 30)
        ops.algorithm("Newton")
        ops.analysis("Static")
        for _ in range(4):
            assert ops.analyze(1) == 0
        return (
            np.array(ops.eleResponse(1, "basicDeformation")),
            np.array(ops.eleResponse(1, "basicForce")),
            np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6),
        )

    base = [0.0001, -0.015, 0.012, 0.0005, 0.002, 0.002]
    v0, q0, kb3 = eval_softening_imk3d(base)
    eps = 1.0e-8
    DV = np.zeros((6, 6))
    DQ = np.zeros((6, 6))
    for j in range(6):
        tp = list(base)
        tm = list(base)
        tp[j] += eps
        tm[j] -= eps
        vp, qp, _ = eval_softening_imk3d(tp)
        vm, qm, _ = eval_softening_imk3d(tm)
        DV[:, j] = (vp - vm) / (2.0 * eps)
        DQ[:, j] = (qp - qm) / (2.0 * eps)
    kb3_fd = DQ @ np.linalg.inv(DV)
    np.testing.assert_allclose(kb3, kb3_fd, rtol=1e-4, atol=1e0)


def test_ladruno_force_beam_column_subdivide_optional_factor_and_isfinite_validation():
    """Verify:
      1. `LadrunoForceBeamColumn2d` and `LadrunoForceBeamColumn3d` accept 1-argument
         `-subdivide <maxSubdivisions>` (defaulting `subdivideFactor` to `10.0`) both
         at the end of the argument list and immediately before another flag
         (`-subdivide 4 -mass 10.0`), as well as 2-argument `-subdivide 4 4.0`.
      2. Both 2D and 3D parsers reject non-finite (`nan`, `inf`, `-inf`) and
         out-of-range values for `-iter` (`tol <= 0` or non-finite) and
         `-subdivide` (`subdivideFactor <= 1.0` or non-finite).
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    for ndm, ndf in ((2, 3), (3, 6)):
        ops.wipe()
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.section("Elastic", 1, E, A, Iz)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 4)

        # Valid: 1-arg -subdivide followed by -mass flag (must not consume or misparse -mass)
        ops.element("LadrunoForceBeamColumn", 2, 1, 2, 1, 1, "-subdivide", 4, "-mass", 10.0)
        eig = ops.eigen("-fullGenLapack", 1)
        assert len(eig) == 1 and eig[0] > 0.0
        # Compare against explicit 2-arg -subdivide with the same -mass 10.0
        ops.remove("element", 2)
        ops.element("LadrunoForceBeamColumn", 3, 1, 2, 1, 1, "-subdivide", 4, 4.0, "-mass", 10.0)
        eig3 = ops.eigen("-fullGenLapack", 1)
        assert eig[0] == pytest.approx(eig3[0], rel=1e-12)
        # Valid: 1-arg -subdivide at end of command
        ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, "-subdivide", 4)
        assert set(ops.getEleTags()) == {1, 3}

        # Invalid -iter tol and -subdivide subdivideFactor (non-finite and out-of-range)
        bad_arg_tuples = [
            ("-iter", 10, float("nan")),
            ("-iter", 10, float("inf")),
            ("-iter", 10, -float("inf")),
            ("-iter", 10, 0.0),
            ("-iter", 10, -1.0e-12),
            ("-subdivide", 4, float("nan")),
            ("-subdivide", 4, float("inf")),
            ("-subdivide", 4, -float("inf")),
            ("-subdivide", 4, 1.0),
            ("-subdivide", 4, 0.5),
            ("-subdivide", 0),
        ]
        for idx, bad_args in enumerate(bad_arg_tuples, start=100):
            try:
                ops.element("LadrunoForceBeamColumn", idx, 1, 2, 1, 1, *bad_args)
            except Exception:
                pass
            assert idx not in ops.getEleTags(), f"ndm={ndm} accepted invalid {bad_args}"


def test_ladruno_frame_validation_fb_dr_enhancements():
    """Verify FB-1..FB-11 fixes, DR-1 (-spanLoad + -hingeIJ), DR-3 (-pmBal), Enhancement C (-axialRelease), and Enhancement D (-powell)."""
    # 1. FB-1: Interior -hingeLoc 0.5 1.0 under uniform bending (q1 = -M0, q2 = +M0)
    L = 2.0
    E, A, Iz = 3.0e10, 0.06, 4.5e-4
    EI = E * Iz
    Mc, Gf, Kpen = 60.0e3, 1.2e3, 1.0e11
    theta_target = 1.5 * (Mc * L) / (2.0 * EI)
    for elem_type in ("LadrunoForceBeamColumn", "LadrunoDispBeamColumn"):
        alphas = {}
        for mode in ("interior_xi05", "endpoint_hingeI"):
            ops.wipe()
            ops.model("basic", "-ndm", 2, "-ndf", 3)
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 0)
            ops.fix(2, 0, 1, 0)
            ops.section("Elastic", 1, E, A, Iz)
            ops.beamIntegration("Lobatto", 1, 1, 5)
            ops.geomTransf("Linear", 1)
            ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, Mc, Gf, "-linear", "-penalty", Kpen)
            ops.uniaxialMaterial("LadrunoCohesiveHinge", 11, 100.0 * Mc, Gf, "-linear", "-penalty", Kpen)
            if mode == "interior_xi05":
                ops.element(elem_type, 1, 1, 2, 1, 1, "-hingeIJ", 10, 11, "-hingeLoc", 0.5, 1.0)
            else:
                ops.element(elem_type, 1, 1, 2, 1, 1, "-hingeI", 10)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.sp(1, 3, -theta_target)
            ops.sp(2, 3, +theta_target)
            ops.system("BandGeneral")
            ops.numberer("Plain")
            ops.constraints("Transformation")
            ops.test("NormDispIncr", 1.0e-12, 30)
            ops.algorithm("Newton")
            ops.integrator("LoadControl", 0.05)
            ops.analysis("Static")
            assert ops.analyze(20) == 0
            hs = ops.eleResponse(1, "hingeState")
            alphas[mode] = abs(float(hs[0]))
        Hs = -(Mc * Mc) / (2.0 * Gf)
        k_pen_eff = (Kpen * Hs) / (Kpen + Hs)
        alpha_exact_00 = (0.5 * Mc + k_pen_eff * (Mc / Kpen)) / (4.0 * EI / L + k_pen_eff)
        f_el = (L / (6.0 * EI)) * np.array([[2.0, -1.0], [-1.0, 2.0]])
        f_el_j = f_el + np.array([[0.0, 0.0], [0.0, 1.0 / Kpen]])
        k_el_j = np.linalg.inv(f_el_j)
        b_05 = np.array([-0.5, 0.5])
        m0_05 = float(b_05 @ k_el_j @ np.array([-theta_target, +theta_target]))
        k_eff_05 = float(b_05 @ k_el_j @ b_05)
        alpha_exact_05 = (m0_05 - Mc + k_pen_eff * (Mc / Kpen)) / (k_eff_05 + k_pen_eff)
        assert alphas["endpoint_hingeI"] == pytest.approx(alpha_exact_00, rel=1.0e-9)
        assert alphas["interior_xi05"] == pytest.approx(alpha_exact_05, rel=1.0e-9)

    # 2. DR-1: LadrunoDispBeamColumn with both -spanLoad and -hingeIJ under uniform load
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 100.0 * Mc, Gf, "-linear", "-penalty", 1.0e15)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeIJ", 10, 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    wy = -24.0e3
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    mz_mid = float(ops.eleResponse(1, "section", 3, "force")[1])
    assert abs(mz_mid - abs(wy) * L * L / 24.0) / (abs(wy) * L * L / 24.0) < 1.0e-7

    # 3. DR-3: -pmBal P0c P0t Pbal muBal on LadrunoCohesiveHinge and LadrunoIMKBeam (2D & 3D)
    P0c, P0t, Pbal, muBal = 4000.0e3, 1200.0e3, 1500.0e3, 1.80
    for ndm, ndf in ((2, 3), (3, 6)):
        ops.wipe()
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.uniaxialMaterial("Steel01", 10, Mc, 1.0e9, 0.02)
            ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-matZi", 10, "-pmBal", P0c, P0t, Pbal, muBal)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, -Pbal, 0.0, 0.0)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.uniaxialMaterial("Steel01", 10, Mc, 1.0e9, 0.02)
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, E / 2.5, 2.0 * Iz, Iz, Iz, 1, "-matZi", 10, "-pmBal", P0c, P0t, Pbal, muBal)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, -Pbal, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        ops.loadConst("-time", 0.0)
        ops.pattern("Plain", 2, 1)
        if ndm == 2:
            ops.load(2, 0.0, 1.0, 0.0)
        else:
            ops.load(2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
        ops.integrator("DisplacementControl", 2, 2, 0.015)
        assert ops.analyze(1) == 0
        hs = ops.eleResponse(1, "hingeState")
        th_i, m_i = abs(float(hs[0])), abs(float(hs[1]))
        m_expected = muBal * (Mc + 0.02 * 1.0e9 * (th_i - Mc / 1.0e9))
        assert m_i == pytest.approx(m_expected, rel=1.0e-9)

    # 4. Enhancement C: -axialRelease on LadrunoForceBeamColumn & LadrunoDispBeamColumn (2D & 3D)
    for elem_type in ("LadrunoForceBeamColumn", "LadrunoDispBeamColumn"):
        for ndm, ndf in ((2, 3), (3, 6)):
            ops.wipe()
            ops.model("basic", "-ndm", ndm, "-ndf", ndf)
            if ndm == 2:
                ops.node(1, 0.0, 0.0)
                ops.node(2, L, 0.0)
                ops.fix(1, 1, 1, 1)
                ops.fix(2, 1, 0, 0)  # axially restrained at both ends!
                ops.geomTransf("Linear", 1)
                ops.uniaxialMaterial("Elastic", 1, E)
                ops.section("Fiber", 1)
                ops.fiber(0.15, 0.0, 0.04, 1)
                ops.fiber(-0.20, 0.0, 0.02, 1)
            else:
                ops.node(1, 0.0, 0.0, 0.0)
                ops.node(2, L, 0.0, 0.0)
                ops.fix(1, 1, 1, 1, 1, 1, 1)
                ops.fix(2, 1, 0, 1, 1, 1, 0)  # axially restrained at both ends!
                ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
                ops.uniaxialMaterial("Elastic", 1, E)
                ops.section("Fiber", 1, "-GJ", 1.0e8)
                ops.fiber(0.15, 0.0, 0.04, 1)
                ops.fiber(-0.20, 0.0, 0.02, 1)
                ops.fiber(0.0, 0.10, 0.01, 1)
                ops.fiber(0.0, -0.10, 0.01, 1)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.element(elem_type, 1, 1, 2, 1, 1, "-axialRelease")
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            if ndm == 2:
                ops.load(2, 0.0, -50.0e3, 0.0)
            else:
                ops.load(2, 0.0, -50.0e3, 0.0, 0.0, 0.0, 0.0)
            ops.system("BandGeneral")
            ops.numberer("Plain")
            ops.constraints("Plain")
            ops.test("NormDispIncr", 1.0e-12, 20)
            ops.algorithm("Newton")
            ops.integrator("LoadControl", 1.0)
            ops.analysis("Static")
            assert ops.analyze(1) == 0
            lf = ops.eleResponse(1, "localForce")
            assert abs(float(lf[0])) == pytest.approx(0.0, abs=1.0e-9)

    # 5. Enhancement D: -powell active-set compliance subtraction on LadrunoIMKBeam (2D & 3D)
    for ndm, ndf in ((2, 3), (3, 6)):
        ops.wipe()
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.uniaxialMaterial("Steel01", 10, Mc, 10.0 * (6.0 * EI / L), 0.02)
            ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10, "-powell")
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 0.0, -10.0e3, 0.0)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.uniaxialMaterial("Steel01", 10, Mc, 10.0 * (6.0 * EI / L), 0.02)
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, E / 2.5, 2.0 * Iz, Iz, Iz, 1, "-hinge", "both", "-matZ", 10, "-powell")
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 0.0, -10.0e3, 0.0, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        # Elastic cantilever tip deflection must match exact PL^3 / (3 EI) with zero penalty-flexibility inflation
        u_tip_exact = -10.0e3 * (L ** 3) / (3.0 * EI)
        assert float(ops.nodeDisp(2, 2)) == pytest.approx(u_tip_exact, rel=1.0e-10)



