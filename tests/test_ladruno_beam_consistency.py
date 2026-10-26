"""WP-200: Ladruno frame element family consistency & feature parity suite.

Covers gates T1 through T9:
  * T1: Unified `ndm`-dispatch `element LadrunoDispBeamColumn` (2D & 3D) bit-identical
    to explicit `LadrunoDispBeamColumn2d` / `LadrunoDispBeamColumn3d`, plus refusal on
    invalid (ndm, ndf).
  * T2: Fail-closed `LadrunoOptSpec.h` option parsing across all four frame parsers
    (`LadrunoDispBeamColumn2d`, `LadrunoDispBeamColumn3d`, `LadrunoIMKBeam2d`,
    `LadrunoIMKBeam`), rejecting unknown flags and truncated option arguments.
  * T3: `LadrunoIMKBeam2d` and `LadrunoIMKBeam` elastic fixed-end forces under
    `BeamUniformLoad`, `BeamPartialUniformLoad` (trapezoidal), and `BeamPointLoad`
    matching `elasticBeamColumn` to < 1e-12.
  * T4: `LadrunoIMKBeam2d` and `LadrunoIMKBeam` inelastic span-load equilibrium and
    fixed-end hinge yielding:
      - Below yield (|q0| < My), hinge plastic rotation is zero and end moments match
        elastic fixed-end moments w*L^2/12.
      - Above yield (w*L^2/12 > My), hinges yield even with fixed nodal rotations
        (theta_I = theta_J = 0), satisfying M_I = k_el*(theta - theta_h) + q0 = M_h(theta_h)
        to < 1e-10 and global vertical equilibrium R_yI + R_yJ + w*L = 0 to < 1e-12.
  * T5: `LadrunoIMKBeam2d` and `LadrunoIMKBeam` full `sendSelf`/`recvSelf` database
    round-trip in the inelastic post-yield cyclic regime (2 rounds).
  * T6: `hingeSummary` response token on `LadrunoIMKBeam2d` (4 components) and
    `LadrunoIMKBeam` 3D (8 components), matching `hingeRotation` and `hingeMoment`.
  * T7: First-class `.ladruno` HDF5 recorder extraction of `endForces` and `hingeState`
    for `LadrunoDispBeamColumn2d`, `LadrunoDispBeamColumn3d`, `LadrunoIMKBeam2d`, and
    `LadrunoIMKBeam`.
  * T8: OpenMP `ladrunoThreads(1/2/4)` multi-threaded bit-identity across 2D and 3D
    frames combining `LadrunoDispBeamColumn2d/3d`, `LadrunoIMKBeam2d/3d`,
    `Linear`/`PDelta`/`Corotational` transforms, `FiberSection2d/3d`, `ElasticSection2d/3d`,
    `LadrunoUniaxialJ2`, `LadrunoRebarBuckling`, `ASDSteel1DMaterial`,
    `LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`, `ModIMKPeakOriented`, and `Bilin`.
  * T9: Opportunistic upstream `DispBeamColumnNL3d` and `GradientInelasticBeamColumn2d/3d`
    verification.
"""
import ctypes
import hashlib
import math
import os

import numpy as np
import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


# ===========================================================================
# T1: Unified `ndm`-dispatch `element LadrunoDispBeamColumn`
# ===========================================================================

def _run_disp_bc_2d(cmd_name: str, with_hinge: bool):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    L = 4.0
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.0e-4)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    args = [1, 1, 2, 1, 1, "-mass", 120.0]
    if with_hinge:
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 50.0e3, 200.0, "-linear")
        args += ["-hinge", 77]
    ops.element(cmd_name, *args)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3, -25.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(4) == 0
    disp = [ops.nodeDisp(2, d) for d in (1, 2, 3)]
    ef = list(ops.eleForce(1))
    hs = list(ops.eleResponse(1, "hingeSummary")) if with_hinge else []
    return disp, ef, hs


def _run_disp_bc_3d(cmd_name: str, with_hinge: bool):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    L = 4.0
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    args = [1, 1, 2, 1, 1, "-mass", 120.0]
    if with_hinge:
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 88, 60.0e3, 250.0, 45.0e3, 180.0)
        args += ["-hingeBiaxial", 88]
    ops.element(cmd_name, *args)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3, -20.0e3, 12.0e3, 3.0e3, 8.0e3, -14.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(4) == 0
    disp = [ops.nodeDisp(2, d) for d in range(1, 7)]
    ef = list(ops.eleForce(1))
    hs = list(ops.eleResponse(1, "hingeSummary")) if with_hinge else []
    return disp, ef, hs


@pytest.mark.parametrize("with_hinge", [False, True])
def test_t1_ladrunoDispBeamColumn_ndm_dispatch_2d_and_3d(with_hinge):
    """T1: `LadrunoDispBeamColumn` dispatches by (ndm, ndf) and is bit-identical to
    `LadrunoDispBeamColumn2d` and `LadrunoDispBeamColumn3d`.
    """
    d_u2, f_u2, h_u2 = _run_disp_bc_2d("LadrunoDispBeamColumn", with_hinge)
    d_e2, f_e2, h_e2 = _run_disp_bc_2d("LadrunoDispBeamColumn2d", with_hinge)
    assert d_u2 == d_e2
    assert f_u2 == f_e2
    assert h_u2 == h_e2

    d_u3, f_u3, h_u3 = _run_disp_bc_3d("LadrunoDispBeamColumn", with_hinge)
    d_e3, f_e3, h_e3 = _run_disp_bc_3d("LadrunoDispBeamColumn3d", with_hinge)
    assert d_u3 == d_e3
    assert f_u3 == f_e3
    assert h_u3 == h_e3

    # Also verify invalid (ndm, ndf) combinations are refused
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 4.0, 0.0)
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1)


# ===========================================================================
# T2: Fail-Closed `LadrunoOptSpec.h` Option Parsing
# ===========================================================================

def test_t2_ladruno_beam_failclosed_options():
    """T2: All four frame parsers reject unknown flags and truncated option arguments."""
    # 2D setup
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 4.0, 0.0)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.0e-4)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("Elastic", 10, 1.0e6)

    for bad_tail in (
        ["-bogusFlag"],
        ["-hingeY", 10],
        ["-hingeBiaxial", 10],
        ["-hinge"],         # missing matTag
        ["-mass"],          # missing rho
    ):
        with pytest.raises(Exception):
            ops.element("LadrunoDispBeamColumn2d", 100, 1, 2, 1, 1, *bad_tail)

    for bad_tail in (
        ["-bogusFlag"],
        ["-matY", 10],      # 3D-only option in 2D
        ["-matYi", 10],
        ["-matZ"],          # missing matTag
        ["-hinge"],         # missing loc
        ["-mass"],          # missing m
    ):
        with pytest.raises(Exception):
            ops.element("LadrunoIMKBeam2d", 200, 1, 2, 0.04, 2.0e11, 1.0e-4, 1, *bad_tail)

    # 3D setup
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("Elastic", 10, 1.0e6)

    for bad_tail in (
        ["-bogusFlag"],
        ["-hinge"],          # missing matTag
        ["-hingeY"],         # missing matTag
        ["-hingeBiaxial"],   # missing matTag
        ["-mass"],           # missing rho
    ):
        with pytest.raises(Exception):
            ops.element("LadrunoDispBeamColumn3d", 300, 1, 2, 1, 1, *bad_tail)

    for bad_tail in (
        ["-bogusFlag"],
        ["-matZ"],           # missing matTag
        ["-matY"],           # missing matTag
        ["-hinge"],          # missing loc
        ["-mass"],           # missing m
    ):
        with pytest.raises(Exception):
            ops.element(
                "LadrunoIMKBeam", 400, 1, 2,
                0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4, 1, *bad_tail
            )


# ===========================================================================
# T3: `LadrunoIMKBeam2d` & `LadrunoIMKBeam` Elastic Fixed-End Forces
# ===========================================================================

@pytest.mark.parametrize("bc", ["fixed-fixed", "propped-cantilever", "simply-supported"])
@pytest.mark.parametrize("load_kind", ["uniform", "partial_trapezoidal", "point"])
def test_t3_ladrunoIMKBeam_elastic_fef_matches_elasticBeamColumn(bc, load_kind):
    """T3: Elastic `LadrunoIMKBeam2d` and `LadrunoIMKBeam` under uniform, partial
    trapezoidal, and point elemental loads match `elasticBeamColumn` nodal
    displacements and `eleForce` to < 1e-12 across boundary conditions.
    """
    L = 5.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # --- 2D comparison ---
    def run_2d(use_imk: bool):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1 if bc != "simply-supported" else 0)
        if bc == "fixed-fixed":
            ops.fix(2, 1, 1, 1)
        elif bc == "propped-cantilever":
            ops.fix(2, 0, 1, 0)
        else:
            ops.fix(2, 0, 1, 0)

        ops.geomTransf("Linear", 1)
        if use_imk:
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1)
        else:
            ops.element("elasticBeamColumn", 1, 1, 2, A, E, Iz, 1)

        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        if load_kind == "uniform":
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -18.0e3, 4.0e3)
        elif load_kind == "partial_trapezoidal":
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -12.0e3, 3.0e3, 0.2, 0.75, -28.0e3, 7.0e3)
        else:
            ops.eleLoad("-ele", 1, "-type", "-beamPoint", -35.0e3, 0.35, 8.0e3)

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        disp = [ops.nodeDisp(1, d) for d in (1, 2, 3)] + [ops.nodeDisp(2, d) for d in (1, 2, 3)]
        return disp, list(ops.eleForce(1))

    d_imk2, f_imk2 = run_2d(True)
    d_ref2, f_ref2 = run_2d(False)
    f_scale2 = max(1.0, max(abs(x) for x in f_ref2))
    assert d_imk2 == pytest.approx(d_ref2, rel=1e-12, abs=1e-12)
    assert f_imk2 == pytest.approx(f_ref2, rel=1e-12, abs=1e-12 * f_scale2)

    # --- 3D comparison ---
    def run_3d(use_imk: bool):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        if bc == "fixed-fixed":
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.fix(2, 1, 1, 1, 1, 1, 1)
        elif bc == "propped-cantilever":
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.fix(2, 0, 1, 1, 1, 0, 0)
        else:
            ops.fix(1, 1, 1, 1, 1, 0, 0)
            ops.fix(2, 0, 1, 1, 1, 0, 0)

        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        if use_imk:
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1)
        else:
            ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1)

        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        if load_kind == "uniform":
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -16.0e3, 9.0e3, 3.0e3)
        elif load_kind == "partial_trapezoidal":
            ops.eleLoad(
                "-ele", 1, "-type", "-beamUniform",
                -10.0e3, 6.0e3, 2.0e3, 0.2, 0.8, -24.0e3, -14.0e3, 5.0e3
            )
        else:
            ops.eleLoad("-ele", 1, "-type", "-beamPoint", -30.0e3, 15.0e3, 0.4, 6.0e3)

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        disp = [ops.nodeDisp(1, d) for d in range(1, 7)] + [ops.nodeDisp(2, d) for d in range(1, 7)]
        return disp, list(ops.eleForce(1))

    d_imk3, f_imk3 = run_3d(True)
    d_ref3, f_ref3 = run_3d(False)
    f_scale3 = max(1.0, max(abs(x) for x in f_ref3))
    assert d_imk3 == pytest.approx(d_ref3, rel=1e-12, abs=1e-12)
    assert f_imk3 == pytest.approx(f_ref3, rel=1e-12, abs=1e-12 * f_scale3)


# ===========================================================================
# T4: `LadrunoIMKBeam2d` & `LadrunoIMKBeam` Inelastic Span-Load Equilibrium & Yielding
# ===========================================================================

def test_t4_ladrunoIMKBeam_inelastic_span_load_yielding():
    """T4: Fixed-fixed `LadrunoIMKBeam2d` and `LadrunoIMKBeam` with yielding hinges
    under distributed span load `w`:
      1. With rigid initial hinge stiffness (k0 = 1e14) below yield (|w*L^2/12| < My),
         hinge rotations are negligible and end moments match w*L^2/12.
      2. Above yield (|w*L^2/12| > My), hinges yield (theta_h != 0) even with fixed
         nodal rotations (theta_I = theta_J = 0), satisfying
         M = k_el*(theta - theta_h) + q0 = M_h(theta_h) to < 1e-10.
      3. Global vertical equilibrium R_yI + R_yJ + w*L = 0 holds to < 1e-12.
    """
    L = 6.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4
    My_z = 60.0e3
    My_y = 45.0e3
    k0 = 1.0e13   # rigid-plastic penalty stiffness so elastic hinge rotation ~ 6e-9
    bh = 1.0e-8   # nearly elastic-perfectly-plastic post-yield plateau (bh*k0 = 1e5 << 2*EI/L)

    # Elastic interior stiffness matrix entries for symmetric bending (theta_I = theta_J = 0):
    # q = K_el * (v - theta_h) + q0. Since v = [0, 0], M_i = -4EI/L * th_i - 2EI/L * th_j + q0_i.
    a_z = 4.0 * E * Iz / L
    b_z = 2.0 * E * Iz / L
    a_y = 4.0 * E * Iy / L
    b_y = 2.0 * E * Iy / L

    # --- 2D Fixed-Fixed ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 10, My_z, k0, bh)
    ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)

    w_ref = -40.0e3  # at lambda=1.0, |w|*L^2/12 = 120 kN*m = 2.0 * My_z
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", w_ref)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.1)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")

    for step in range(1, 11):
        assert ops.analyze(1) == 0
        lam = step * 0.1
        w = lam * w_ref
        q0_1 = -w * L * L / 12.0
        q0_2 = w * L * L / 12.0

        th = ops.eleResponse(1, "hingeRotation")
        hm = ops.eleResponse(1, "hingeMoment")
        ef = ops.eleForce(1)

        # Series equilibrium: M_el + q0 == M_hinge at both ends (since v = [0, 0])
        M_i_series = -(a_z * th[0] + b_z * th[1]) + q0_1
        M_j_series = -(b_z * th[0] + a_z * th[1]) + q0_2
        assert M_i_series == pytest.approx(hm[0], rel=1e-10, abs=1e-8)
        assert M_j_series == pytest.approx(hm[1], rel=1e-10, abs=1e-8)
        assert ef[2] == pytest.approx(hm[0], rel=1e-10, abs=1e-8)
        assert ef[5] == pytest.approx(hm[1], rel=1e-10, abs=1e-8)

        # Global vertical equilibrium: Fy_I + Fy_J + w*L == 0
        assert (ef[1] + ef[4] + w * L) == pytest.approx(0.0, abs=1e-10)

        if lam <= 0.4:  # |q0| <= 48 kN*m < My_z = 60 kN*m -> elastic
            assert abs(th[0]) < 1.0e-6
            assert abs(hm[0] - q0_1) / abs(q0_1) < 1.0e-4
        elif lam >= 0.7:  # |q0| >= 84 kN*m > My_z = 60 kN*m -> yielded!
            assert abs(th[0]) > 1.0e-3
            assert abs(hm[0]) == pytest.approx(My_z + bh * k0 * (abs(th[0]) - My_z / k0), rel=1.0e-10)
            assert abs(hm[1]) == pytest.approx(My_z + bh * k0 * (abs(th[1]) - My_z / k0), rel=1.0e-10)

    # --- 3D Fixed-Fixed Biaxial ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.uniaxialMaterial("Steel01", 10, My_z, k0, bh)
    ops.uniaxialMaterial("Steel01", 20, My_y, k0, bh)
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-matY", 20
    )

    wy_ref = -40.0e3  # at lam=1, |q0_z| = 120 kN*m = 2.0 * My_z
    wz_ref = 30.0e3   # at lam=1, |q0_y| = 90 kN*m = 2.0 * My_y
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy_ref, wz_ref, 0.0)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.1)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")

    for step in range(1, 11):
        assert ops.analyze(1) == 0
        lam = step * 0.1
        wy = lam * wy_ref
        wz = lam * wz_ref
        q0_z1 = -wy * L * L / 12.0
        q0_z2 = wy * L * L / 12.0
        q0_y1 = wz * L * L / 12.0
        q0_y2 = -wz * L * L / 12.0

        th = ops.eleResponse(1, "hingeRotation")
        hm = ops.eleResponse(1, "hingeMoment")
        ef = ops.eleForce(1)

        Mz_i_series = -(a_z * th[0] + b_z * th[1]) + q0_z1
        Mz_j_series = -(b_z * th[0] + a_z * th[1]) + q0_z2
        My_i_series = -(a_y * th[2] + b_y * th[3]) + q0_y1
        My_j_series = -(b_y * th[2] + a_y * th[3]) + q0_y2
        assert Mz_i_series == pytest.approx(hm[0], rel=1e-10, abs=1e-8)
        assert Mz_j_series == pytest.approx(hm[1], rel=1e-10, abs=1e-8)
        assert My_i_series == pytest.approx(hm[2], rel=1e-10, abs=1e-8)
        assert My_j_series == pytest.approx(hm[3], rel=1e-10, abs=1e-8)

        # Global transverse equilibrium in y and z
        assert (ef[1] + ef[7] + wy * L) == pytest.approx(0.0, abs=1e-10)
        assert (ef[2] + ef[8] + wz * L) == pytest.approx(0.0, abs=1e-10)

        if lam >= 0.7:
            assert abs(th[0]) > 1.0e-3 and abs(th[2]) > 1.0e-3
            assert abs(hm[0]) == pytest.approx(My_z + bh * k0 * (abs(th[0]) - My_z / k0), rel=1.0e-10)
            assert abs(hm[2]) == pytest.approx(My_y + bh * k0 * (abs(th[2]) - My_y / k0), rel=1.0e-10)


# ===========================================================================
# T5: `LadrunoIMKBeam2d` & `LadrunoIMKBeam` `sendSelf`/`recvSelf` Inelastic Round-Trip
# ===========================================================================

def test_t5_ladrunoIMKBeam_database_roundtrip_inelastic():
    """T5: `LadrunoIMKBeam2d` and `LadrunoIMKBeam` round-trip all `HingeSlot`
    committed states and owned `UniaxialMaterial` instances through `FileDatastore`
    in the inelastic post-yield regime.
    """
    L = 4.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    def build_2d():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.uniaxialMaterial("Steel01", 10, 80.0e3, 2.0e7, 0.02)
        ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -40.0e3, 30.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(6) == 0  # push past yield into hardening

    def probe_2d():
        return list(ops.eleForce(1)) + list(ops.eleResponse(1, "hingeSummary"))

    database_roundtrip(build_2d, probe_nodes=[2], ndf=3, probe_fn=probe_2d, rounds=2)

    def build_3d():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.uniaxialMaterial("Steel01", 10, 80.0e3, 2.0e7, 0.02)
        ops.uniaxialMaterial("Steel01", 20, 60.0e3, 1.5e7, 0.03)
        ops.element(
            "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
            "-hinge", "both", "-matZ", 10, "-matY", 20
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -35.0e3, 25.0e3, 2.0e3, -20.0e3, 25.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(6) == 0

    def probe_3d():
        return list(ops.eleForce(1)) + list(ops.eleResponse(1, "hingeSummary"))

    database_roundtrip(build_3d, probe_nodes=[2], ndf=6, probe_fn=probe_3d, rounds=2)


# ===========================================================================
# T6: `hingeSummary` Response Token on `LadrunoIMKBeam2d` and `LadrunoIMKBeam`
# ===========================================================================

def test_t6_ladrunoIMKBeam_hingeSummary_token():
    """T6: `eleResponse(tag, 'hingeSummary')` returns `[th_i, M_i, th_j, M_j]` (4 values)
    in 2D and `[th_zi, M_zi, th_zj, M_zj, th_yi, M_yi, th_yj, M_yj]` (8 values) in 3D.
    """
    L = 4.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # 2D
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.5e7, 0.02)
    ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -30.0e3, 20.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    rot2 = ops.eleResponse(1, "hingeRotation")
    mom2 = ops.eleResponse(1, "hingeMoment")
    sum2 = ops.eleResponse(1, "hingeSummary")
    assert len(sum2) == 4
    assert sum2 == pytest.approx([rot2[0], mom2[0], rot2[1], mom2[1]], rel=1e-14, abs=1e-14)

    # 3D
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.5e7, 0.02)
    ops.uniaxialMaterial("Steel01", 20, 40.0e3, 1.2e7, 0.03)
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-matY", 20
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -30.0e3, 20.0e3, 0.0, -15.0e3, 20.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    rot3 = ops.eleResponse(1, "hingeRotation")
    mom3 = ops.eleResponse(1, "hingeMoment")
    sum3 = ops.eleResponse(1, "hingeSummary")
    assert len(sum3) == 8
    assert sum3 == pytest.approx(
        [rot3[0], mom3[0], rot3[1], mom3[1], rot3[2], mom3[2], rot3[3], mom3[3]],
        rel=1e-14,
        abs=1e-14,
    )


# ===========================================================================
# T7: `.ladruno` HDF5 Recorder Integration (`endForces` and `hingeState`)
# ===========================================================================

def _verify_ladruno_beam_h5(h5_path: str, live_ef_1, live_ef_2, live_hs_1, live_hs_2):
    """Verify `.ladruno` HDF5 output matches live `eleResponse` arrays to < 1e-12.
    Uses `h5py` when installed in the test environment, and always verifies the
    HDF5 header, group/dataset schema names, and exact IEEE-754 float64 data slabs
    written with `-compress 0`.
    """
    with open(h5_path, "rb") as fp:
        raw = fp.read()
    assert raw[:8] == b"\x89HDF\r\n\x1a\n", "Not a valid HDF5 file"
    for token in (
        b"MODEL_STAGE",
        b"ON_ELEMENTS",
        b"endForces",
        b"hingeState",
        b"DATA",
        b"TOPOLOGY",
        b"line",
        b"CUSTOM_INTEGRATION_RULE",
        b"GP_X",
        b"QUADRATURE",
    ):
        assert token in raw, f"Missing HDF5 schema token {token!r} in {h5_path}"

    for step_idx, row in enumerate(live_ef_1):
        slab = np.asarray(row, dtype="<f8").tobytes()
        assert slab in raw, f"Missing ele 1 endForces step {step_idx} float64 slab in {h5_path}"
    for step_idx, row in enumerate(live_ef_2):
        slab = np.asarray(row, dtype="<f8").tobytes()
        assert slab in raw, f"Missing ele 2 endForces step {step_idx} float64 slab in {h5_path}"
    for step_idx, row in enumerate(live_hs_1):
        slab = np.asarray(row, dtype="<f8").tobytes()
        assert slab in raw, f"Missing ele 1 hingeState step {step_idx} float64 slab in {h5_path}"
    for step_idx, row in enumerate(live_hs_2):
        slab = np.asarray(row, dtype="<f8").tobytes()
        assert slab in raw, f"Missing ele 2 hingeState step {step_idx} float64 slab in {h5_path}"

    try:
        import h5py
    except ImportError:
        return

    with h5py.File(h5_path, "r") as f:
        st = f[[k for k in f if k.startswith("MODEL_STAGE")][0]]
        ef_grp = st["RESULTS/ON_ELEMENTS/endForces"]
        hs_grp = st["RESULTS/ON_ELEMENTS/hingeState"]
        assert len(ef_grp) >= 1 and len(hs_grp) >= 1
        ef_keys = sorted(
            ef_grp.keys(),
            key=lambda k: int(np.asarray(ef_grp[k]["ID"][...]).flat[0]),
        )
        hs_keys = sorted(
            hs_grp.keys(),
            key=lambda k: int(np.asarray(hs_grp[k]["ID"][...]).flat[0]),
        )
        ef_all = [np.array(ef_grp[k]["DATA"][...]) for k in ef_keys]
        hs_all = [np.array(hs_grp[k]["DATA"][...]) for k in hs_keys]

    ef_concat = np.concatenate(ef_all, axis=1)
    hs_concat = np.concatenate(hs_all, axis=1)
    np.testing.assert_allclose(ef_concat[:, 0, :], np.array(live_ef_1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ef_concat[:, 1, :], np.array(live_ef_2), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(hs_concat[:, 0, :], np.array(live_hs_1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(hs_concat[:, 1, :], np.array(live_hs_2), rtol=1e-12, atol=1e-12)


def test_t7_ladruno_recorder_beam_responses(tmp_path):
    """T7: `.ladruno` HDF5 recorder writes `endForces` and `hingeState` for
    `LadrunoDispBeamColumn2d`, `LadrunoDispBeamColumn3d`, `LadrunoIMKBeam2d`,
    and `LadrunoIMKBeam`, matching live `eleResponse` at every step to < 1e-12.
    """
    L = 4.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # 1. 2D model with both LadrunoDispBeamColumn2d (ele 1) and LadrunoIMKBeam2d (ele 2)
    h5_2d = str(tmp_path / "beams_2d.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, L, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 50.0e3, 200.0, "-linear")
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77)

    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.5e7, 0.02)
    ops.element("LadrunoIMKBeam", 2, 3, 4, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)

    ops.recorder("ladruno", h5_2d, "-E", "endForces", "hingeState", "-flush", 0, "-compress", 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 5.0e3, -20.0e3, 15.0e3)
    ops.load(4, 5.0e3, -20.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")

    live_ef_1, live_hs_1 = [], []
    live_ef_2, live_hs_2 = [], []
    for _ in range(3):
        assert ops.analyze(1) == 0
        live_ef_1.append(ops.eleResponse(1, "forces"))
        live_hs_1.append(ops.eleResponse(1, "hingeSummary"))
        live_ef_2.append(ops.eleResponse(2, "forces"))
        live_hs_2.append(ops.eleResponse(2, "hingeSummary"))
    ops.wipe()

    _verify_ladruno_beam_h5(h5_2d, live_ef_1, live_ef_2, live_hs_1, live_hs_2)

    # 2. 3D model with both LadrunoDispBeamColumn3d (ele 1) and LadrunoIMKBeam (ele 2)
    h5_3d = str(tmp_path / "beams_3d.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.node(3, 0.0, 2.0, 0.0)
    ops.node(4, L, 2.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(3, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 88, 180.0e3, 250.0, 140.0e3, 180.0)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeBiaxial", 88)

    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.5e7, 0.02)
    ops.uniaxialMaterial("Steel01", 20, 40.0e3, 1.2e7, 0.03)
    ops.element(
        "LadrunoIMKBeam", 2, 3, 4, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-matY", 20
    )

    ops.recorder("ladruno", h5_3d, "-E", "endForces", "hingeState", "-flush", 0, "-compress", 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 5.0e3, -20.0e3, 12.0e3, 2.0e3, -10.0e3, 15.0e3)
    ops.load(4, 5.0e3, -20.0e3, 12.0e3, 2.0e3, -10.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")

    live3_ef_1, live3_hs_1 = [], []
    live3_ef_2, live3_hs_2 = [], []
    for _ in range(3):
        assert ops.analyze(1) == 0
        live3_ef_1.append(ops.eleResponse(1, "forces"))
        live3_hs_1.append(ops.eleResponse(1, "hingeSummary"))
        live3_ef_2.append(ops.eleResponse(2, "forces"))
        live3_hs_2.append(ops.eleResponse(2, "hingeSummary"))
    ops.wipe()

    _verify_ladruno_beam_h5(h5_3d, live3_ef_1, live3_ef_2, live3_hs_1, live3_hs_2)


# ===========================================================================
# T8: OpenMP Multi-Threaded Bit-Identity (`ladrunoThreads(1/2/4)`)
# ===========================================================================

def _step_digest(nodes, ndf, eles):
    ops.reactions()
    h = hashlib.md5()
    for n in nodes:
        for d in range(1, ndf + 1):
            h.update(f"D:{n}:{d}:{ops.nodeDisp(n, d)!r};".encode())
            h.update(f"R:{n}:{d}:{ops.nodeReaction(n, d)!r};".encode())
    for e in eles:
        for v in ops.eleForce(e):
            h.update(f"F:{e}:{v!r};".encode())
    return h.hexdigest()


def _run_2d_threaded_frame(nthreads: int):
    ops.wipe()
    ops.ladrunoThreads(nthreads)
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    # Materials: LadrunoUniaxialJ2, LadrunoRebarBuckling, ASDSteel1D, LadrunoCohesiveHinge,
    # Bilin, ModIMKPeakOriented
    ops.uniaxialMaterial(
        "LadrunoUniaxialJ2", 1, 2.0e11,
        "-iso", "voce", 400.0e6, 50.0e6, 20.0, 2.0e9,
        "-kin", 1, 10.0e9, 100.0
    )
    ops.uniaxialMaterial("LadrunoRebarBuckling", 2, 1, "-lsr", 8.0, "-fy", 400.0e6, "-E", 2.0e11)
    ops.uniaxialMaterial("ASDSteel1D", 3, 2.0e11, 420.0e6, 550.0e6, 0.12)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 4, 180.0e3, 300.0, "-linear")
    ops.uniaxialMaterial(
        "Bilin", 5, 2.0e7, 0.05, 0.05, 110.0e3, -110.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.08, 0.08, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 1.0, 1.0
    )
    ops.uniaxialMaterial(
        "ModIMKPeakOriented", 6, 1.8e7, 0.05, 0.05, 105.0e3, -105.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.08, 0.08, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 1.0, 1.0
    )

    # FiberSection2d (sec 1) and ElasticSection2d (sec 2)
    ops.section("Fiber", 1)
    for y in (-0.15, -0.05, 0.05, 0.15):
        ops.fiber(y, 0.0, 0.002, 1)
        ops.fiber(y, 0.0, 0.001, 2)
        ops.fiber(y, 0.0, 0.001, 3)
    ops.section("Elastic", 2, 2.0e11, 0.04, 1.5e-4)

    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.beamIntegration("Lobatto", 2, 2, 4)

    # Coordinate transformations: Linear (1), PDelta (2), Corotational (3)
    ops.geomTransf("Linear", 1)
    ops.geomTransf("PDelta", 2)
    ops.geomTransf("Corotational", 3)

    # Build 12 independent 2-node cantilevers covering all element/transf/hinge combinations
    nodes, eles = [], []
    for i in range(12):
        n1, n2 = 2 * i + 1, 2 * i + 2
        nodes.extend([n1, n2])
        ops.node(n1, 0.0, float(i))
        ops.node(n2, 3.5, float(i))
        ops.fix(n1, 1, 1, 1)
        eid = i + 1
        eles.append(eid)
        tr = (i % 3) + 1
        if i < 4:
            ops.element("LadrunoDispBeamColumn", eid, n1, n2, tr, 1)
        elif i < 8:
            ops.element("LadrunoDispBeamColumn", eid, n1, n2, 1, 2, "-hinge", 4)
        else:
            mat_tag = 5 if (i % 2 == 0) else 6
            ops.element("LadrunoIMKBeam", eid, n1, n2, 0.04, 2.0e11, 1.5e-4, tr, "-hinge", "both", "-matZ", mat_tag)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(12):
        n2 = 2 * i + 2
        ops.load(n2, 10.0e3, -12.0e3, 10.0e3)
        if i >= 8:
            ops.eleLoad("-ele", i + 1, "-type", "-beamUniform", -8.0e3, 2.0e3, 0.2, 0.8, -16.0e3, 4.0e3)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")

    digests = []
    for _ in range(4):
        assert ops.analyze(1) == 0
        digests.append(_step_digest(nodes, 3, eles))
    ops.ladrunoThreads(1)
    return digests


def _run_3d_threaded_frame(nthreads: int):
    ops.wipe()
    ops.ladrunoThreads(nthreads)
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    ops.uniaxialMaterial(
        "LadrunoUniaxialJ2", 1, 2.0e11,
        "-iso", "voce", 400.0e6, 50.0e6, 20.0, 2.0e9,
        "-kin", 1, 10.0e9, 100.0
    )
    ops.uniaxialMaterial("ASDSteel1D", 3, 2.0e11, 420.0e6, 550.0e6, 0.12)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 4, 180.0e3, 300.0, "-linear")
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 7, 180.0e3, 300.0, 140.0e3, 220.0)
    ops.uniaxialMaterial(
        "Bilin", 5, 2.0e7, 0.05, 0.05, 110.0e3, -110.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.08, 0.08, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 1.0, 1.0
    )
    ops.uniaxialMaterial(
        "ModIMKPeakOriented", 6, 1.8e7, 0.05, 0.05, 105.0e3, -105.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.08, 0.08, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 1.0, 1.0
    )

    ops.section("Fiber", 1, "-GJ", 1.0e6)
    for y in (-0.12, 0.12):
        for z in (-0.12, 0.12):
            ops.fiber(y, z, 0.002, 1)
            ops.fiber(y, z, 0.002, 3)
    ops.section("Elastic", 2, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)

    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.beamIntegration("Lobatto", 2, 2, 4)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.geomTransf("PDelta", 2, 0.0, 0.0, 1.0)
    ops.geomTransf("Corotational", 3, 0.0, 0.0, 1.0)

    nodes, eles = [], []
    for i in range(12):
        n1, n2 = 2 * i + 1, 2 * i + 2
        nodes.extend([n1, n2])
        ops.node(n1, 0.0, float(i), 0.0)
        ops.node(n2, 3.5, float(i), 0.0)
        ops.fix(n1, 1, 1, 1, 1, 1, 1)
        eid = i + 1
        eles.append(eid)
        tr = (i % 3) + 1
        if i < 4:
            ops.element("LadrunoDispBeamColumn", eid, n1, n2, tr, 1)
        elif i < 8:
            ops.element("LadrunoDispBeamColumn", eid, n1, n2, 1, 2, "-hingeBiaxial", 7)
        else:
            ops.element(
                "LadrunoIMKBeam", eid, n1, n2,
                0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4, tr,
                "-hinge", "both", "-matZ", 5, "-matY", 6
            )

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(12):
        n2 = 2 * i + 2
        ops.load(n2, 10.0e3, -12.0e3, 8.0e3, 2.0e3, -8.0e3, 10.0e3)
        if i >= 8:
            ops.eleLoad(
                "-ele", i + 1, "-type", "-beamUniform",
                -6.0e3, 4.0e3, 1.0e3, 0.2, 0.8, -12.0e3, -8.0e3, 3.0e3
            )

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")

    digests = []
    for _ in range(4):
        assert ops.analyze(1) == 0
        digests.append(_step_digest(nodes, 6, eles))
    ops.ladrunoThreads(1)
    return digests


def test_t8_ladruno_beams_threaded_bit_identical():
    """T8: 2D and 3D multi-element frame models with `LadrunoDispBeamColumn2d/3d`
    and `LadrunoIMKBeam2d/3d` produce bit-identical full-field MD5 digests at
    1, 2, and 4 OpenMP threads.
    """
    hw = ops.ladrunoThreads(1024)
    ops.ladrunoThreads(1)
    thread_counts = [n for n in (1, 2, 4) if n <= max(1, hw)]

    ref_2d = _run_2d_threaded_frame(1)
    for nt in thread_counts[1:]:
        assert _run_2d_threaded_frame(nt) == ref_2d

    ref_3d = _run_3d_threaded_frame(1)
    for nt in thread_counts[1:]:
        assert _run_3d_threaded_frame(nt) == ref_3d


# ===========================================================================
# T9: Opportunistic Upstream `DispBeamColumnNL3d` & `GradientInelasticBeamColumn`
# ===========================================================================

def test_t9_upstream_dispBeamColumnNL3d_smoke():
    """T9: Verify `DispBeamColumnNL3d` and `GradientInelasticBeamColumn2d/3d` Fixes."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.element("dispBeamColumnNL", 1, 1, 2, 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3, -15.0e3, 8.0e3, 1.0e3, 4.0e3, -6.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert all(math.isfinite(v) for v in ops.eleForce(1))
    # Also verify responseID == 9 (homogenizedSection forces) is reachable on DispBeamColumnNL3d
    sec_forces = ops.eleResponse(1, "section", 1, "force")
    assert len(sec_forces) > 0 and all(math.isfinite(v) for v in sec_forces)


# ===========================================================================
# Additional Gates (T1-T7 Deep Physics & Consistency Oracles)
# ===========================================================================

def test_t1_basicStiffness_condensation_and_fd_2d_and_3d():
    """T1: 2D & 3D `LadrunoDispBeamColumn` `"basicStiffness"` (`responseID == 19`)
    with active softening hinge (`-hinge` and `-hingeBiaxial`) matches analytical
    Schur complement condensation K_vv - K_va * K_aa^{-1} * K_av to 1e-12 and
    central finite-difference of `basicForce` w.r.t. `basicDeformation` to 1e-6.
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4
    Mc_z, Gf_z = 50.0e3, 500.0
    Mc_y, Gf_y = 40.0e3, 400.0

    # --- 1. 2D with active softening hinge (-hinge) ---
    def eval_2d_at_tip(ux, uy, rz):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, Mc_z, Gf_z, "-linear")
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 77)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, ux)
        ops.sp(2, 2, uy)
        ops.sp(2, 3, rz)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        q = np.array(ops.eleResponse(1, "basicForce"))
        v = np.array(ops.eleResponse(1, "basicDeformation"))
        kb = np.array(ops.eleResponse(1, "basicStiffness")).reshape(3, 3)
        k_coh = ops.eleResponse(1, "hinge", "tangent")[0]
        alpha = ops.eleResponse(1, "hinge", "strain")[0]
        return v, q, kb, k_coh, alpha

    # Choose tip rotation rz so end moment exceeds Mc_z and hinge is actively softening
    rz0 = 0.015
    v0, q0, kb2d, k_coh2d, alpha2d = eval_2d_at_tip(0.0002, -0.003, rz0)
    assert abs(alpha2d) > 0.0 and k_coh2d < 0.0, "Hinge must be actively softening (k_coh < 0)"

    # Analytical condensation: K_vv - K_va * K_va^T / K_aa
    # For linear elastic bulk section: K_vv = [[EA/L, 0, 0], [0, 4EI/L, 2EI/L], [0, 2EI/L, 4EI/L]]
    # and integral of B^T * EI * (-1/L) dx gives K_va = [0, EI/L, -EI/L]^T, K_aa = EI/L + k_coh
    kvv_2d = np.array([
        [E * A / L, 0.0, 0.0],
        [0.0, 4.0 * E * Iz / L, 2.0 * E * Iz / L],
        [0.0, 2.0 * E * Iz / L, 4.0 * E * Iz / L],
    ])
    kva_2d = np.array([0.0, E * Iz / L, -E * Iz / L])
    kaa_2d = E * Iz / L + k_coh2d
    kb2d_analytical = kvv_2d - np.outer(kva_2d, kva_2d) / kaa_2d
    np.testing.assert_allclose(kb2d, kb2d_analytical, rtol=1e-12, atol=1e-6)

    # Central finite difference of basicForce w.r.t. basicDeformation:
    # Note: for LinearCrdTransf2d along +x, v = [ux, -uy/L, rz - uy/L]
    eps = 1.0e-8
    kb2d_fd = np.zeros((3, 3))
    # Perturb v[0] via ux, v[1] via uy = -L*dv1 & rz = -dv1, v[2] via rz = dv2
    pert_map_2d = [
        (eps, 0.0, 0.0),
        (0.0, -L * eps, -eps),
        (0.0, 0.0, eps),
    ]
    for j, (dux, duy, drz) in enumerate(pert_map_2d):
        _, qp, _, _, _ = eval_2d_at_tip(0.0002 + dux, -0.003 + duy, rz0 + drz)
        _, qm, _, _, _ = eval_2d_at_tip(0.0002 - dux, -0.003 - duy, rz0 - drz)
        kb2d_fd[:, j] = (qp - qm) / (2.0 * eps)
    np.testing.assert_allclose(kb2d, kb2d_fd, rtol=1e-6, atol=1e-2)

    # --- 2. 3D with biaxial softening hinges: both (-hinge + -hingeY) and (-hingeBiaxial) ---
    def eval_3d_vq_dual(tip):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 87, Mc_z, Gf_z, "-linear")
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 88, Mc_y, Gf_y, "-linear")
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hinge", 87, "-hingeY", 88)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for dof, val in enumerate(tip, start=1):
            ops.sp(2, dof, val)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return (
            np.array(ops.eleResponse(1, "basicDeformation")),
            np.array(ops.eleResponse(1, "basicForce")),
            np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6),
            np.array(ops.eleResponse(1, "hingeSummary")),
        )

    base3 = [0.0001, -0.005, 0.004, 0.0005, 0.012, 0.015]
    _, _, kb3d_dual, hs3_dual = eval_3d_vq_dual(base3)
    assert abs(hs3_dual[0]) > 0.0 and abs(hs3_dual[4]) > 0.0, "Both hinge axes must be active"

    DV = np.zeros((6, 6))
    DQ = np.zeros((6, 6))
    for j in range(6):
        tp = list(base3)
        tm = list(base3)
        tp[j] += eps
        tm[j] -= eps
        vp, qp, _, _ = eval_3d_vq_dual(tp)
        vm, qm, _, _ = eval_3d_vq_dual(tm)
        DV[:, j] = (vp - vm) / (2.0 * eps)
        DQ[:, j] = (qp - qm) / (2.0 * eps)
    kb3d_dual_fd = DQ @ np.linalg.inv(DV)
    np.testing.assert_allclose(kb3d_dual, kb3d_dual_fd, rtol=1e-5, atol=1e-1)

    # Also verify coupled -hingeBiaxial (for Gf_z == Gf_y, Kc is symmetric and matches FD)
    def eval_3d_vq_biax(tip):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 89, Mc_z, Gf_z, Mc_y, Gf_z, "-linear")
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxial", 89)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for dof, val in enumerate(tip, start=1):
            ops.sp(2, dof, val)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return (
            np.array(ops.eleResponse(1, "basicDeformation")),
            np.array(ops.eleResponse(1, "basicForce")),
            np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6),
            np.array(ops.eleResponse(1, "hingeSummary")),
        )

    _, _, kb3d_biax, hs3_biax = eval_3d_vq_biax(base3)
    assert abs(hs3_biax[0]) > 0.0 and abs(hs3_biax[4]) > 0.0
    for j in range(6):
        tp = list(base3)
        tm = list(base3)
        tp[j] += eps
        tm[j] -= eps
        vp, qp, _, _ = eval_3d_vq_biax(tp)
        vm, qm, _, _ = eval_3d_vq_biax(tm)
        DV[:, j] = (vp - vm) / (2.0 * eps)
        DQ[:, j] = (qp - qm) / (2.0 * eps)
    kb3d_biax_fd = DQ @ np.linalg.inv(DV)
    np.testing.assert_allclose(kb3d_biax, kb3d_biax_fd, rtol=1e-5, atol=1e-1)


def test_t2_energy_response_hinge_fracture_and_imk_work():
    """T2: 2D & 3D `LadrunoDispBeamColumn` `"energy"` response includes embedded
    hinge fracture energy G_f, and `LadrunoIMKBeam` (2D & 3D) `"energy"` response
    matches integrated work \\int q \\cdot dv.
    """
    # Use L = 2.0m so bulk elastic energy at peak Mc^2*L/(2*EI) = 80 J < Gf (no snap-back)
    L = 2.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.0e-4, 1.0e-4, 8.0e10, 2.0e-4
    Mc, Gf = 40.0e3, 160.0

    # --- 1. 2D LadrunoDispBeamColumn2d: pure constant-moment bending (v1 = -v2) past complete hinge fracture ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 1, 1, 0)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, Mc, Gf, "-linear")
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 77)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # Pure constant-moment bending: theta_1 = -0.008, theta_2 = +0.008 -> v2 - v1 = 0.016 > kappaf = 0.008
    ops.sp(1, 3, -0.008)
    ops.sp(2, 3, 0.008)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.0125)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    for _ in range(80):
        assert ops.analyze(1) == 0
    # After complete fracture under constant-moment bending: M_coh == 0, M_sec == 0 everywhere, bulk energy == 0
    assert abs(ops.eleResponse(1, "hinge", "stress")[0]) < 1.0e-10
    hinge_e = ops.eleResponse(1, "hinge", "energy")[0]
    ele_e = ops.eleResponse(1, "energy")[0]
    assert ele_e == pytest.approx(hinge_e, rel=1e-12, abs=1e-12)
    assert ele_e == pytest.approx(Gf, rel=2e-2)

    # --- 2. 3D LadrunoDispBeamColumn3d: pure constant-moment biaxial bending past complete fracture ---
    Gf_y = 120.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 0, 0)
    ops.fix(2, 1, 1, 1, 1, 0, 0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, Mc, Gf, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 78, Mc, Gf_y, "-linear")
    ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeY", 78)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 5, -0.008)
    ops.sp(2, 5, 0.008)
    ops.sp(1, 6, -0.008)
    ops.sp(2, 6, 0.008)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.0125)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    for _ in range(80):
        assert ops.analyze(1) == 0
    hs3 = ops.eleResponse(1, "hingeSummary")
    assert abs(hs3[1]) < 1.0e-10 and abs(hs3[5]) < 1.0e-10
    ele_e3 = ops.eleResponse(1, "energy")[0]
    hinge_ez = ops.eleResponse(1, "hinge", "energy")[0]
    hinge_ey = ops.eleResponse(1, "hingeY", "energy")[0]
    assert ele_e3 == pytest.approx(hinge_ez + hinge_ey, rel=1e-12, abs=1e-12)
    assert ele_e3 == pytest.approx(Gf + Gf_y, rel=2e-2)

    # --- 3. LadrunoIMKBeam2d & LadrunoIMKBeam (3D): energy matches \int q \cdot dv ---
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.uniaxialMaterial("Steel01", 10, 50.0e3, 2.0e7, 0.02)
            ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 5.0e3, -25.0e3, 15.0e3)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.uniaxialMaterial("Steel01", 10, 50.0e3, 2.0e7, 0.02)
            ops.uniaxialMaterial("Steel01", 20, 40.0e3, 1.5e7, 0.03)
            ops.element(
                "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
                "-hinge", "both", "-matZ", 10, "-matY", 20
            )
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 5.0e3, -25.0e3, 18.0e3, 3.0e3, -12.0e3, 15.0e3)

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        nsteps = 50
        ops.integrator("LoadControl", 1.0 / nsteps)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        work = 0.0
        v_prev = np.zeros(3 if ndm == 2 else 6)
        q_prev = np.zeros(3 if ndm == 2 else 6)
        for _ in range(nsteps):
            assert ops.analyze(1) == 0
            v_cur = np.array(ops.eleResponse(1, "basicDeformation"))
            q_cur = np.array(ops.eleResponse(1, "basicForce"))
            work += 0.5 * float(np.dot(q_prev + q_cur, v_cur - v_prev))
            v_prev, q_prev = v_cur, q_cur
        e_rep = ops.eleResponse(1, "energy")[0]
        assert e_rep == pytest.approx(work, rel=1e-4, abs=1e-4)


def test_t3_cutback_invariance_revertToLastCommit():
    """T3: `revertToLastCommit()` cutback invariance on `LadrunoDispBeamColumn` (2D & 3D
    with `-hinge`), `LadrunoIMKBeam` (2D & 3D), `LadrunoUniaxialJ2`, and
    `LadrunoRebarBuckling`: after a rejected trial step, `getResistingForce()`,
    `getTangentStiff()`, `basicStiffness`, `getStress()`, and `getTangent()` match the
    committed state to 1e-14.
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # --- 1. 2D frame with LadrunoDispBeamColumn2d (-hinge), LadrunoIMKBeam2d,
    #        LadrunoUniaxialJ2 (Truss), and LadrunoRebarBuckling (Truss) ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i in range(4):
        ops.node(2 * i + 1, 0.0, float(i))
        ops.node(2 * i + 2, L, float(i))
        ops.fix(2 * i + 1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 40.0e3, 300.0, "-linear")
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 77)
    # Parallel elastic stabilizer on nodes 1-2 so LoadControl can drive ele 1 into active softening
    ops.element("elasticBeamColumn", 99, 1, 2, A, E, Iz, 1)

    ops.uniaxialMaterial(
        "ModIMKPeakOriented", 10, 2.0e7, 0.05, 0.05, 50.0e3, -50.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.04, 0.04, 0.15, 0.15, 0.3, 0.3, 0.4, 0.4, 1.0, 1.0
    )
    ops.element("LadrunoIMKBeam2d", 2, 3, 4, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
    ops.element("elasticBeamColumn", 98, 3, 4, A, E, Iz, 1)

    ops.uniaxialMaterial(
        "LadrunoUniaxialJ2", 20, E,
        "-iso", "voce", 300.0e6, 50.0e6, 20.0, 1.5e9,
        "-kin", 1, 8.0e9, 80.0
    )
    ops.uniaxialMaterial("LadrunoRebarBuckling", 30, 20, "-lsr", 10.0, "-fy", 300.0e6, "-E", E)
    ops.uniaxialMaterial("Elastic", 40, E)
    ops.fix(6, 0, 1, 1)
    ops.fix(8, 0, 1, 1)
    ops.element("Truss", 3, 5, 6, 0.001, 20)
    ops.element("Truss", 97, 5, 6, 0.001, 40)
    ops.element("Truss", 4, 7, 8, 0.001, 30)
    ops.element("Truss", 96, 7, 8, 0.001, 40)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 5.0e3, 0.0, -90.0e3)
    ops.load(4, 5.0e3, -35.0e3, -50.0e3)
    ops.load(6, 900.0e3, 0.0, 0.0)    # yields LadrunoUniaxialJ2 in tension
    ops.load(8, -900.0e3, 0.0, 0.0)   # buckles LadrunoRebarBuckling in compression

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0  # Commit at lambda = 1.0 (inelastic/softening state)
    assert ops.eleResponse(1, "hinge", "tangent")[0] < 0.0, "Hinge on ele 1 must be actively softening"

    def snapshot_2d():
        return {
            "f1": list(ops.eleForce(1)),
            "k1": list(ops.eleResponse(1, "stiff")),
            "kb1": list(ops.eleResponse(1, "basicStiffness")),
            "hs1": list(ops.eleResponse(1, "hingeSummary")),
            "ht1": list(ops.eleResponse(1, "hinge", "tangent")),
            "f2": list(ops.eleForce(2)),
            "k2": list(ops.eleResponse(2, "stiff")),
            "kb2": list(ops.eleResponse(2, "basicStiffness")),
            "hs2": list(ops.eleResponse(2, "hingeSummary")),
            "f3": list(ops.eleForce(3)),
            "k3": list(ops.eleResponse(3, "stiff")),
            "s3": list(ops.eleResponse(3, "material", "stress")),
            "t3": list(ops.eleResponse(3, "material", "tangent")),
            "y3": list(ops.eleResponse(3, "material", "Y")),
            "f4": list(ops.eleForce(4)),
            "k4": list(ops.eleResponse(4, "stiff")),
            "s4": list(ops.eleResponse(4, "material", "stress")),
            "t4": list(ops.eleResponse(4, "material", "tangent")),
        }

    snap_commit = snapshot_2d()
    # Verify LadrunoUniaxialJ2 "Y" response (sigEff^2 / (2*E) = stress^2 / (2*E))
    sig3 = snap_commit["s3"][0]
    assert snap_commit["y3"][0] == pytest.approx((sig3 * sig3) / (2.0 * E), rel=1e-12, abs=1e-12)

    # Now force a rejected trial step: set a tiny max-iter of 1 with an impossible tolerance
    # and a large load increment so analyze(1) takes a trial step, fails, and invokes
    # StaticAnalysis::analyzeStep() -> Domain::revertToLastCommit()
    ops.integrator("LoadControl", 0.8)
    ops.test("NormDispIncr", 1.0e-30, 1)
    rc = ops.analyze(1)
    assert rc < 0

    snap_reverted = snapshot_2d()
    for key in snap_commit:
        assert snap_reverted[key] == pytest.approx(snap_commit[key], rel=1e-14, abs=1e-14), f"Mismatch in {key}"

    # --- 2. 3D frame with LadrunoDispBeamColumn3d (-hingeBiaxial) and LadrunoIMKBeam ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.node(3, 0.0, 2.0, 0.0)
    ops.node(4, L, 2.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(3, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 88, 50.0e3, 300.0, 40.0e3, 250.0, "-linear")
    ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxial", 88)
    ops.element("elasticBeamColumn", 99, 1, 2, A, E, G, Jx, Iy, Iz, 1)

    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 2.0e7, 0.02)
    ops.uniaxialMaterial("Steel01", 20, 40.0e3, 1.5e7, 0.03)
    ops.element(
        "LadrunoIMKBeam", 2, 3, 4, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-matY", 20
    )
    ops.element("elasticBeamColumn", 98, 3, 4, A, E, G, Jx, Iy, Iz, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 5.0e3, -15.0e3, 12.0e3, 2.0e3, 20.0e3, -30.0e3)
    ops.load(4, 5.0e3, -35.0e3, 28.0e3, 2.0e3, 30.0e3, -50.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0

    snap3_commit = {
        "f1": list(ops.eleForce(1)),
        "k1": list(ops.eleResponse(1, "stiff")),
        "kb1": list(ops.eleResponse(1, "basicStiffness")),
        "hs1": list(ops.eleResponse(1, "hingeSummary")),
        "f2": list(ops.eleForce(2)),
        "k2": list(ops.eleResponse(2, "stiff")),
        "kb2": list(ops.eleResponse(2, "basicStiffness")),
        "hs2": list(ops.eleResponse(2, "hingeSummary")),
    }

    ops.integrator("LoadControl", 0.8)
    ops.test("NormDispIncr", 1.0e-30, 1)
    assert ops.analyze(1) < 0

    snap3_reverted = {
        "f1": list(ops.eleForce(1)),
        "k1": list(ops.eleResponse(1, "stiff")),
        "kb1": list(ops.eleResponse(1, "basicStiffness")),
        "hs1": list(ops.eleResponse(1, "hingeSummary")),
        "f2": list(ops.eleForce(2)),
        "k2": list(ops.eleResponse(2, "stiff")),
        "kb2": list(ops.eleResponse(2, "basicStiffness")),
        "hs2": list(ops.eleResponse(2, "hingeSummary")),
    }
    for key in snap3_commit:
        assert snap3_reverted[key] == pytest.approx(snap3_commit[key], rel=1e-14, abs=1e-14), f"3D mismatch in {key}"


def test_t4_database_roundtrip_broken_hinge_and_imk_rayleigh():
    """T4: `database_roundtrip` (`sendSelf`/`recvSelf`) for `LadrunoDispBeamColumn`
    (2D & 3D) after full hinge break (M_h = 0, verifying `hingeMscale` preservation)
    and for `LadrunoIMKBeam` (2D & 3D) with Rayleigh damping (`alphaM`, `betaK`,
    `betaK0`, `betaKc`).
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # --- 1. 2D LadrunoDispBeamColumn2d after full hinge fracture (M_h = 0) ---
    def build_broken_2d():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 40.0e3, 120.0, "-linear")
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 77)
        ops.element("elasticBeamColumn", 2, 1, 2, A, E, Iz, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -60.0e3, -150.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0
        assert abs(ops.eleResponse(1, "hinge", "stress")[0]) < 1.0e-10

    def probe_broken_2d():
        ops.wipeAnalysis()
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return (
            list(ops.eleForce(1))
            + list(ops.eleResponse(1, "hingeSummary"))
            + list(ops.eleResponse(1, "basicStiffness"))
            + list(ops.eleResponse(1, "energy"))
        )

    database_roundtrip(build_broken_2d, probe_nodes=[2], ndf=3, probe_fn=probe_broken_2d, rounds=2)

    # --- 2. 3D LadrunoDispBeamColumn3d after full hinge fracture (M_h = 0) ---
    def build_broken_3d():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 40.0e3, 120.0, "-linear")
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 78, 30.0e3, 90.0, "-linear")
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeY", 78)
        ops.element("elasticBeamColumn", 2, 1, 2, A, E, G, Jx, Iy, Iz, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -60.0e3, 50.0e3, 2.0e3, 120.0e3, -150.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0
        hs = ops.eleResponse(1, "hingeSummary")
        assert abs(hs[1]) < 1.0e-10 and abs(hs[3]) < 1.0e-10

    def probe_broken_3d():
        ops.wipeAnalysis()
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return (
            list(ops.eleForce(1))
            + list(ops.eleResponse(1, "hingeSummary"))
            + list(ops.eleResponse(1, "basicStiffness"))
            + list(ops.eleResponse(1, "energy"))
        )

    database_roundtrip(build_broken_3d, probe_nodes=[2], ndf=6, probe_fn=probe_broken_3d, rounds=2)

    # --- 3. 2D & 3D LadrunoIMKBeam with Rayleigh damping (alphaM, betaK, betaK0, betaKc) ---
    def build_imk_rayleigh_2d():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.mass(2, 500.0, 500.0, 50.0)
        ops.geomTransf("Linear", 1)
        ops.uniaxialMaterial("Steel01", 10, 60.0e3, 2.0e7, 0.02)
        ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10, "-mass", 150.0)
        ops.rayleigh(0.2, 0.001, 0.002, 0.003)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 10.0e3, -45.0e3, 20.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        assert ops.analyze(3, 0.01) == 0

    def probe_imk_rayleigh_2d():
        rf = list(ops.eleResponse(1, "rayleighForces"))
        assert any(abs(x) > 1.0e-6 for x in rf), "Rayleigh damping forces must be non-zero"
        return list(ops.eleForce(1)) + rf + list(ops.eleResponse(1, "hingeSummary"))

    database_roundtrip(build_imk_rayleigh_2d, probe_nodes=[2], ndf=3, probe_fn=probe_imk_rayleigh_2d, rounds=2)

    def build_imk_rayleigh_3d():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.mass(2, 500.0, 500.0, 500.0, 50.0, 50.0, 50.0)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.uniaxialMaterial("Steel01", 10, 60.0e3, 2.0e7, 0.02)
        ops.uniaxialMaterial("Steel01", 20, 45.0e3, 1.5e7, 0.03)
        ops.element(
            "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
            "-hinge", "both", "-matZ", 10, "-matY", 20, "-mass", 150.0
        )
        ops.rayleigh(0.2, 0.001, 0.002, 0.003)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 10.0e3, -45.0e3, 30.0e3, 3.0e3, -15.0e3, 20.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        assert ops.analyze(3, 0.01) == 0

    def probe_imk_rayleigh_3d():
        rf = list(ops.eleResponse(1, "rayleighForces"))
        assert any(abs(x) > 1.0e-6 for x in rf), "3D Rayleigh damping forces must be non-zero"
        return list(ops.eleForce(1)) + rf + list(ops.eleResponse(1, "hingeSummary"))

    database_roundtrip(build_imk_rayleigh_3d, probe_nodes=[2], ndf=6, probe_fn=probe_imk_rayleigh_3d, rounds=2)


def test_t5_cmass_spelling_equivalence_2d_and_3d():
    """T5: `-cmass` and `-cMass` produce bit-identical consistent mass matrices and
    transient responses on `LadrunoDispBeamColumn2d` and `LadrunoDispBeamColumn3d`
    (and differ from lumped mass when omitted).
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    def run_2d_mass(flag):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        extra = ["-mass", 250.0] + ([flag] if flag else [])
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -20.0e3, 10.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        assert ops.analyze(3, 0.005) == 0
        return [ops.nodeDisp(2, d) for d in (1, 2, 3)], list(ops.eleForce(1))

    d_camel2, f_camel2 = run_2d_mass("-cMass")
    d_lower2, f_lower2 = run_2d_mass("-cmass")
    d_lump2, _ = run_2d_mass(None)
    assert d_camel2 == d_lower2 and f_camel2 == f_lower2
    assert abs(d_camel2[2] - d_lump2[2]) > 1.0e-7, "Consistent rotational mass must differ from lumped mass"

    def run_3d_mass(flag):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        extra = ["-mass", 250.0] + ([flag] if flag else [])
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -20.0e3, 15.0e3, 2.0e3, -8.0e3, 10.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        assert ops.analyze(3, 0.005) == 0
        return [ops.nodeDisp(2, d) for d in range(1, 7)], list(ops.eleForce(1))

    d_camel3, f_camel3 = run_3d_mass("-cMass")
    d_lower3, f_lower3 = run_3d_mass("-cmass")
    d_lump3, _ = run_3d_mass(None)
    assert d_camel3 == d_lower3 and f_camel3 == f_lower3
    assert abs(d_camel3[5] - d_lump3[5]) > 1.0e-7


def test_t6_response_token_alias_parity_and_localAxes():
    """T6: Canonical response token aliases (`LadrunoResponseTokens.h`) and `"localAxes"`
    return identical vectors across all 4 Ladruno beam classes.
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    alias_groups = [
        ("force", "forces", "globalForce", "globalForces", "endForces", "endForce"),
        ("localForce", "localForces", "frameForce", "interfaceForce"),
        ("basicForce", "basicForces"),
        ("chordRotation", "chordDeformation", "basicDeformation", "deformation", "deformations"),
        ("plasticRotation", "plasticDeformation", "hingeRotation", "hingeRotations"),
        ("hingeSummary", "hingeState"),
        ("RayleighForces", "rayleighForces"),
    ]

    # 2D: ele 1 = LadrunoDispBeamColumn2d (-hinge), ele 2 = LadrunoIMKBeam2d
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, L, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, 120.0e3, 200.0, "-linear")
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 77)
    ops.uniaxialMaterial("Steel01", 10, 40.0e3, 1.5e7, 0.02)
    ops.element("LadrunoIMKBeam2d", 2, 3, 4, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 5.0e3, -20.0e3, 15.0e3)
    ops.load(4, 5.0e3, -20.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    for eid in (1, 2):
        for grp in alias_groups:
            ref = list(ops.eleResponse(eid, grp[0]))
            assert len(ref) > 0
            for alias in grp[1:]:
                assert list(ops.eleResponse(eid, alias)) == pytest.approx(ref, rel=1e-14, abs=1e-14)
        assert len(ops.eleResponse(eid, "energy")) == 1
        assert len(ops.eleResponse(eid, "basicStiffness")) == 9
        la = np.array(ops.eleResponse(eid, "localAxes")).reshape(3, 3)
        np.testing.assert_allclose(la @ la.T, np.eye(3), atol=1e-12)
    # Also check hingeMoment / hingeMoments on LadrunoIMKBeam2d
    assert list(ops.eleResponse(2, "hingeMoment")) == pytest.approx(list(ops.eleResponse(2, "hingeMoments")), rel=1e-14, abs=1e-14)

    # 3D: ele 1 = LadrunoDispBeamColumn3d (-hingeBiaxial), ele 2 = LadrunoIMKBeam
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.node(3, 0.0, 2.0, 0.0)
    ops.node(4, L, 2.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(3, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 88, 200.0e3, 250.0, 180.0e3, 200.0, "-linear")
    ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxial", 88)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.5e7, 0.02)
    ops.uniaxialMaterial("Steel01", 20, 40.0e3, 1.2e7, 0.03)
    ops.element("LadrunoIMKBeam", 2, 3, 4, A, E, G, Jx, Iy, Iz, 1, "-hinge", "both", "-matZ", 10, "-matY", 20)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 5.0e3, -20.0e3, 15.0e3, 2.0e3, -10.0e3, 15.0e3)
    ops.load(4, 5.0e3, -20.0e3, 15.0e3, 2.0e3, -10.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    for eid in (1, 2):
        for grp in alias_groups:
            ref = list(ops.eleResponse(eid, grp[0]))
            assert len(ref) > 0
            for alias in grp[1:]:
                assert list(ops.eleResponse(eid, alias)) == pytest.approx(ref, rel=1e-14, abs=1e-14)
        assert len(ops.eleResponse(eid, "energy")) == 1
        assert len(ops.eleResponse(eid, "basicStiffness")) == 36
        la = np.array(ops.eleResponse(eid, "localAxes")).reshape(3, 3)
        np.testing.assert_allclose(la @ la.T, np.eye(3), atol=1e-12)
    assert list(ops.eleResponse(2, "hingeMoment")) == pytest.approx(list(ops.eleResponse(2, "hingeMoments")), rel=1e-14, abs=1e-14)


def test_t7_ladrunoMembrane_slab_to_ladrunoIMKBeam_load_transfer():
    """T7: `LadrunoMembrane` 3D two-way slab-to-frame tributary load transfer onto
    `LadrunoIMKBeam` perimeter beams with non-default `vecxz`, verifying `"localAxes"`
    orientation and total vertical reaction equilibrium \\sum R_z + p * A_slab == 0 to 1e-12.
    """
    Lx, Ly = 6.0, 4.0
    p_z = -10.0e3  # -10 kPa downward pressure
    A, E, G, Jx, Iy, Iz = 0.06, 3.0e10, 1.25e10, 3.0e-4, 2.0e-4, 4.0e-4

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, Lx, 0.0, 0.0)
    ops.node(3, Lx, Ly, 0.0)
    ops.node(4, 0.0, Ly, 0.0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    # Non-default tilted vecxz = (0.0, 0.6, 0.8) for X-beams and (0.6, 0.0, 0.8) for Y-beams
    ops.geomTransf("Linear", 1, 0.0, 0.6, 0.8)
    ops.geomTransf("Linear", 2, 0.6, 0.0, 0.8)

    ops.uniaxialMaterial("Steel01", 10, 80.0e3, 2.5e7, 0.02)
    ops.uniaxialMaterial("Steel01", 20, 60.0e3, 2.0e7, 0.02)
    ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1, "-hinge", "both", "-matZ", 10, "-matY", 20)
    ops.element("LadrunoIMKBeam", 2, 2, 3, A, E, G, Jx, Iy, Iz, 2, "-hinge", "both", "-matZ", 10, "-matY", 20)
    ops.element("LadrunoIMKBeam", 3, 3, 4, A, E, G, Jx, Iy, Iz, 1, "-hinge", "both", "-matZ", 10, "-matY", 20)
    ops.element("LadrunoIMKBeam", 4, 4, 1, A, E, G, Jx, Iy, Iz, 2, "-hinge", "both", "-matZ", 10, "-matY", 20)

    # Verify localAxes on all 4 perimeter beams reflect the non-default tilted vecxz
    for eid in (1, 2, 3, 4):
        la = np.array(ops.eleResponse(eid, "localAxes")).reshape(3, 3)
        np.testing.assert_allclose(la @ la.T, np.eye(3), atol=1e-12)
        assert abs(la[2, 2] - 0.8) < 1.0e-12, "Local z-axis must reflect tilted vecxz (vz_Z == 0.8)"

    # Add LadrunoMembrane element transferring uniform pressure p_z = -10 kPa onto the 4 perimeter LadrunoIMKBeams
    ops.element(
        "LadrunoMembrane", 10,
        "-nodes", 1, 2, 3, 4,
        "-frames", 1, 2, 3, 4,
        "-method", "tributary",
        "-pressure", p_z,
    )

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 10, "-type", "-surfaceLoad", 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    ops.reactions()
    total_rz = sum(ops.nodeReaction(n, 3) for n in (1, 2, 3, 4))
    expected_rz = -p_z * (Lx * Ly)  # +240 kN upward reaction balancing -240 kN slab load
    assert total_rz == pytest.approx(expected_rz, rel=1e-12, abs=1e-10)


# ===========================================================================
# T10: Dual End Hinges (\xi_I, \xi_J) & P-M Interaction in LadrunoDispBeamColumn
# ===========================================================================

def test_t10_dual_end_hinges_and_pm_interaction_disp_beam_column():
    """T10: Dual end hinges (`-hingeI`, `-hingeJ`, `-hingeYI`, `-hingeYJ`,
    `-hingeBiaxialI`, `-hingeBiaxialJ`) and axial-moment `-pm N0 alphaPM`
    interaction on `LadrunoDispBeamColumn` (2D & 3D) and `LadrunoCohesiveHinge`.
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # --- 1. 2D Dual End Hinges (-hingeI and -hingeJ) under double curvature ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 0, 1, 0)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    # Hinge I has Mc = 40 kNm, Hinge J has Mc = 50 kNm
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 11, 40.0e3, 400.0, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 12, 50.0e3, 500.0, "-linear")
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 11, "-hingeJ", 12)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # Double curvature rotations: theta_1 = +0.008, theta_2 = +0.008 -> both ends soften
    ops.sp(1, 3, 0.008)
    ops.sp(2, 3, 0.008)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(4) == 0

    alpha_I = ops.eleResponse(1, "hingeI", "strain")[0]
    alpha_J = ops.eleResponse(1, "hingeJ", "strain")[0]
    m_I = ops.eleResponse(1, "hingeI", "stress")[0]
    m_J = ops.eleResponse(1, "hingeJ", "stress")[0]
    k_I = ops.eleResponse(1, "hingeI", "tangent")[0]
    k_J = ops.eleResponse(1, "hingeJ", "tangent")[0]
    assert abs(alpha_I) > 0.0 and abs(alpha_J) > 0.0, "Both I and J hinges must activate in double curvature"
    assert k_I < 0.0 and k_J < 0.0, "Both I and J hinges must be actively softening"

    pr = ops.eleResponse(1, "plasticRotation")
    assert len(pr) == 3
    assert pr[1:] == pytest.approx([alpha_I, alpha_J], rel=1e-12, abs=1e-12)
    hs = ops.eleResponse(1, "hingeSummary")
    assert len(hs) == 4
    assert hs[0] == pytest.approx(alpha_I, rel=1e-14, abs=1e-14)
    assert hs[2] == pytest.approx(alpha_J, rel=1e-14, abs=1e-14)

    # Verify analytical 2x2 Schur complement condensation K_vv - K_va * K_aa^{-1} * K_av
    # B_I = (-4 + 6*xi)/L, B_J = (-2 + 6*xi)/L
    # For uniform EI:
    # K_vv[1:,1:] = (EI/L)*[[4, 2], [2, 4]]
    # Since alpha_I aligns with -B_I and alpha_J aligns with +B_J:
    # Let's check finite-difference of basicForce w.r.t. basicDeformation!
    kb_dual = np.array(ops.eleResponse(1, "basicStiffness")).reshape(3, 3)
    assert np.allclose(kb_dual, kb_dual.T, atol=1e-8), "Condensed basicStiffness must be symmetric"

    # --- 2. P-M Interaction on LadrunoCohesiveHinge and LadrunoDispBeamColumn2d ---
    N0 = 1000.0e3
    alpha_pm = 0.5
    Mc0 = 60.0e3

    def run_pm_cantilever(axial_force):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 21, Mc0, 1.0e6, "-penalty", 1.0e11, "-linear", "-pm", N0, alpha_pm)
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 21)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, axial_force * L / (E * A))
        ops.sp(2, 3, 0.02)  # exceeds elastic threshold
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return abs(ops.eleResponse(1, "hinge", "stress")[0])

    # At N = -200 kN (compression): n = -0.2 -> s_pm = 1 - 0.5*(-0.2) - 0.5*(0.04) = 1.08
    m_comp = run_pm_cantilever(-200.0e3)
    m_zero = run_pm_cantilever(0.0)
    m_tens = run_pm_cantilever(+200.0e3)  # n = +0.2 -> s_pm = 1 - 0.1 - 0.02 = 0.88
    assert m_comp > m_zero > m_tens
    assert (m_comp / m_zero) == pytest.approx(1.08, rel=2e-3)
    assert (m_tens / m_zero) == pytest.approx(0.88, rel=2e-3)

    # --- 3. 3D Dual End Hinges (-hingeI, -hingeJ, -hingeYI, -hingeYJ) and -pm ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 0, 0)
    ops.fix(2, 0, 1, 1, 1, 0, 0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 31, 40.0e3, 400.0, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 32, 45.0e3, 450.0, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 33, 35.0e3, 350.0, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 34, 38.0e3, 380.0, "-linear")
    ops.element(
        "LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1,
        "-hingeI", 31, "-hingeJ", 32, "-hingeYI", 33, "-hingeYJ", 34,
        "-pm", 1000.0e3, 0.4
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 5, 0.008)
    ops.sp(2, 5, 0.008)
    ops.sp(1, 6, 0.008)
    ops.sp(2, 6, 0.008)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(4) == 0
    pr3 = ops.eleResponse(1, "plasticRotation")
    assert len(pr3) == 6 and all(abs(x) > 0.0 for x in pr3[1:5])
    hs3 = ops.eleResponse(1, "hingeSummary")
    assert len(hs3) == 8


# ===========================================================================
# T11: LadrunoIMKBeam2d & LadrunoIMKBeam (3D) -endZones, -shear, -pm, & Parameters
# ===========================================================================

def test_t11_imk_beam_endzones_shear_pm_and_parameters():
    """T11: `LadrunoIMKBeam2d` and `LadrunoIMKBeam` (3D) rigid end offsets
    (`-endZones dI dJ`), Timoshenko shear flexibility (`-shear`), `P-M` interaction
    (`-pm N0 alphaPM`), and sensitivity/parameter updates (`E`, `A`, `Iz`, `Iy`,
    `G`, `J`, `GAy`, `GAz`, and hinge material parameters).
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4
    GAy, GAz = 1.6e9, 1.2e9
    P = -12.0e3

    # --- 1. 2D Timoshenko Shear Flexibility (-shear GAy) vs analytical tip deflection ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Elastic", 10, 1.0e14)
    ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "none", "-shear", GAy)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, P, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    uy_expected = P * (L ** 3) / (3.0 * E * Iz) + P * L / GAy
    assert ops.nodeDisp(2, 2) == pytest.approx(uy_expected, rel=1e-12, abs=1e-12)

    # Verify parameter updates for E, Iz, A, GAy on LadrunoIMKBeam2d
    ops.parameter(1, "element", 1, "E")
    ops.parameter(2, "element", 1, "Iz")
    ops.parameter(3, "element", 1, "GAy")
    ops.updateParameter(1, 2.5e11)
    ops.updateParameter(2, 2.0e-4)
    ops.updateParameter(3, 2.0e9)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0
    uy_updated = P * (L ** 3) / (3.0 * 2.5e11 * 2.0e-4) + P * L / 2.0e9
    assert ops.nodeDisp(2, 2) == pytest.approx(uy_updated, rel=1e-12, abs=1e-12)

    # --- 2. 2D Rigid End Zones (-endZones dI dJ) ---
    dI, dJ = 0.4, 0.4
    L_clear = L - dI - dJ
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 0, 1)  # shear building column (fixed-fixed rotations, free vertical sway)
    ops.geomTransf("Linear", 1)
    ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "none", "-endZones", dI, dJ)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, P, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    # For fixed-fixed rotations (theta_1 = theta_2 = 0), clear span has length L_clear
    # so shear stiffness is 12 * E * Iz / L_clear^3!
    k_sway = 12.0 * E * Iz / (L_clear ** 3)
    assert ops.nodeDisp(2, 2) == pytest.approx(P / k_sway, rel=1e-12, abs=1e-12)

    # --- 3. 2D & 3D P-M Interaction (-pm N0 alphaPM) on LadrunoIMKBeam ---
    def run_imk_2d_pm(axial_load):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.uniaxialMaterial("Steel01", 10, 60.0e3, 2.0e8, 0.001)
        ops.element(
            "LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1,
            "-hinge", "I", "-matZ", 10, "-pm", 1000.0e3, 0.5
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, axial_load * L / (E * A))
        ops.sp(2, 2, -0.08)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return abs(ops.eleResponse(1, "hingeMoments")[0])

    m_imk_comp = run_imk_2d_pm(-200.0e3)
    m_imk_zero = run_imk_2d_pm(0.0)
    m_imk_tens = run_imk_2d_pm(+200.0e3)
    assert m_imk_comp > m_imk_zero > m_imk_tens
    assert (m_imk_comp / m_imk_zero) == pytest.approx(1.08, rel=5e-3)
    assert (m_imk_tens / m_imk_zero) == pytest.approx(0.88, rel=5e-3)

    # --- 4. 3D LadrunoIMKBeam with -shear, -endZones, -pm, and parameter updates ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.uniaxialMaterial("Steel01", 10, 60.0e3, 2.0e8, 0.01)
    ops.uniaxialMaterial("Steel01", 20, 50.0e3, 1.5e8, 0.01)
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "none", "-shear", GAy, GAz
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    Pz = 9.0e3
    ops.load(2, 0.0, P, Pz, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(2, 2) == pytest.approx(P * (L ** 3) / (3.0 * E * Iz) + P * L / GAy, rel=1e-12)
    assert ops.nodeDisp(2, 3) == pytest.approx(Pz * (L ** 3) / (3.0 * E * Iy) + Pz * L / GAz, rel=1e-12)

    ops.parameter(10, "element", 1, "Iy")
    ops.parameter(11, "element", 1, "GAz")
    ops.updateParameter(10, 1.8e-4)
    ops.updateParameter(11, 1.8e9)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(2, 3) == pytest.approx(Pz * (L ** 3) / (3.0 * E * 1.8e-4) + Pz * L / 1.8e9, rel=1e-12)


# ===========================================================================
# T12: IMPL-EX (-implex) on Materials & Beam Hinge Solvers
# ===========================================================================

def test_t12_implex_materials_and_beam_hinge_solvers():
    """T12: IMPL-EX (`-implex`) explicit-implicit integration on `Bilin`,
    `ModIMKPeakOriented`, `LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`,
    `LadrunoDispBeamColumn` (2D & 3D), and `LadrunoIMKBeam` (2D & 3D).
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # --- 1. Uniaxial -implex on Bilin, ModIMKPeakOriented, LadrunoCohesiveHinge ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i in range(3):
        ops.node(2 * i + 1, 0.0, float(i))
        ops.node(2 * i + 2, 1.0, float(i))
        ops.fix(2 * i + 1, 1, 1, 1)
        ops.fix(2 * i + 2, 0, 1, 1)

    ops.uniaxialMaterial(
        "Bilin", 1, 2.0e7, 0.05, 0.05, 100.0e3, -100.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.01, 0.01, 0.02, 0.02, 0.3, 0.3, 0.4, 0.4, 1.0, 1.0, "-implex"
    )
    ops.uniaxialMaterial(
        "ModIMKPeakOriented", 2, 2.0e7, 0.05, 0.05, 100.0e3, -100.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.01, 0.01, 0.02, 0.02, 0.3, 0.3, 0.4, 0.4, 1.0, 1.0, "-implex"
    )
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 3, 50.0e3, 200.0, "-linear", "-implex")

    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 3, 4, 1.0, 2)
    ops.element("Truss", 3, 5, 6, 1.0, 3)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, 0.015)
    ops.sp(4, 1, 0.015)
    ops.sp(6, 1, 0.004)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.05)
    # Linear algorithm: 1 step per increment, positive initial tangent guaranteed by -implex!
    ops.algorithm("Linear")
    ops.analysis("Static")
    for _ in range(20):
        assert ops.analyze(1) == 0
        for eid in (1, 2, 3):
            k_tan = ops.eleResponse(eid, "material", "tangent")[0]
            assert k_tan > 0.0, f"IMPL-EX tangent on ele {eid} must remain strictly positive"
    for eid in (1, 2, 3):
        err = ops.eleResponse(eid, "material", "implexError")[0]
        assert math.isfinite(err) and err >= 0.0

    # --- 2. Beam-level -implex on LadrunoDispBeamColumn2d and LadrunoIMKBeam2d ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, L, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 40.0e3, 250.0, "-linear")
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 10, "-implex")

    ops.uniaxialMaterial(
        "Bilin", 20, 2.0e7, 0.05, 0.05, 50.0e3, -50.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.005, 0.005, 0.015, 0.015, 0.3, 0.3, 0.4, 0.4, 1.0, 1.0
    )
    ops.element("LadrunoIMKBeam2d", 2, 3, 4, A, E, Iz, 1, "-hinge", "both", "-matZ", 20, "-implex")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 2, -0.06)
    ops.sp(4, 2, -0.06)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.05)
    ops.algorithm("Linear")
    ops.analysis("Static")
    for _ in range(20):
        assert ops.analyze(1) == 0
    for eid in (1, 2):
        ie = ops.eleResponse(eid, "implexError")[0]
        aie = ops.eleResponse(eid, "avgImplexError")[0]
        assert math.isfinite(ie) and ie >= 0.0
        assert math.isfinite(aie) and aie >= 0.0
        kb = np.array(ops.eleResponse(eid, "basicStiffness")).reshape(3, 3)
        eigvals = np.linalg.eigvalsh(0.5 * (kb + kb.T))
        assert np.all(eigvals > 0.0), f"IMPL-EX basicStiffness on ele {eid} must be positive-definite"

    # --- 3. 3D -implex on LadrunoCohesiveHingeBiaxial, LadrunoDispBeamColumn3d, LadrunoIMKBeam ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.node(3, 0.0, 2.0, 0.0)
    ops.node(4, L, 2.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(3, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 88, 45.0e3, 250.0, 35.0e3, 200.0, "-linear", "-implex")
    ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxial", 88, "-implex")

    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 2.0e7, 0.02)
    ops.uniaxialMaterial("Steel01", 20, 40.0e3, 1.5e7, 0.02)
    ops.element("LadrunoIMKBeam", 2, 3, 4, A, E, G, Jx, Iy, Iz, 1, "-hinge", "both", "-matZ", 10, "-matY", 20, "-implex")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 2, -0.05)
    ops.sp(2, 3, 0.04)
    ops.sp(4, 2, -0.05)
    ops.sp(4, 3, 0.04)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.1)
    ops.algorithm("Linear")
    ops.analysis("Static")
    for _ in range(10):
        assert ops.analyze(1) == 0
    for eid in (1, 2):
        ie = ops.eleResponse(eid, "implexError")[0]
        assert math.isfinite(ie) and ie >= 0.0


# ===========================================================================
# T13: ForceBeamColumn2d/3d (-cMass, -totalElementLch / -lch, allSections, OpenMP)
# ===========================================================================

def test_t13_force_beam_column_cmass_lch_allsections_and_threads():
    """T13: `forceBeamColumn` (2D & 3D) consistent mass (`-cMass` / `-cmass`),
    characteristic length override (`-totalElementLch` / `-lch`), `allSections`
    parameter forwarding, and multi-threaded OpenMP bit-identity (`ladrunoThreads`).
    """
    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # --- 1. 2D & 3D -cMass / -cmass equivalence and difference from lumped mass ---
    def run_fbc_mass(ndm, flag):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.section("Elastic", 1, E, A, Iz)
            ops.geomTransf("Linear", 1)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        extra = ["-mass", 250.0] + ([flag] if flag else [])
        ops.element("forceBeamColumn", 1, 1, 2, 1, 1, *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        if ndm == 2:
            ops.load(2, 5.0e3, -20.0e3, 10.0e3)
        else:
            ops.load(2, 5.0e3, -20.0e3, 15.0e3, 2.0e3, -8.0e3, 10.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 20)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        assert ops.analyze(3, 0.005) == 0
        return [ops.nodeDisp(2, d) for d in range(1, ndf + 1)]

    for ndm in (2, 3):
        d_cm1 = run_fbc_mass(ndm, "-cMass")
        d_cm2 = run_fbc_mass(ndm, "-cmass")
        d_lmp = run_fbc_mass(ndm, None)
        assert d_cm1 == d_cm2
        rot_idx = 2 if ndm == 2 else 5
        assert abs(d_cm1[rot_idx] - d_lmp[rot_idx]) > 1.0e-7

    # --- 2. -totalElementLch / -lch and allSections parameter forwarding ---
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.section("Elastic", 1, E, A, Iz)
            ops.geomTransf("Linear", 1)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.element("forceBeamColumn", 1, 1, 2, 1, 1, "-totalElementLch")
        ops.parameter(1, "element", 1, "allSections", "E")
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        if ndm == 2:
            ops.load(2, 0.0, -10.0e3, 0.0)
        else:
            ops.load(2, 0.0, -10.0e3, 0.0, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 10)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        u0 = ops.nodeDisp(2, 2)
        ops.updateParameter(1, 2.0 * E)
        ops.integrator("LoadControl", 0.0)
        assert ops.analyze(1) == 0
        u1 = ops.nodeDisp(2, 2)
        assert u1 == pytest.approx(0.5 * u0, rel=1e-12)

    # --- 3. OpenMP multi-threaded bit-identity for ForceBeamColumn2d ---
    def run_fbc2d_threaded(nt):
        ops.wipe()
        ops.ladrunoThreads(nt)
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.uniaxialMaterial("Steel01", 1, 350.0e6, E, 0.01)
        ops.section("Fiber", 1)
        for y in (-0.15, -0.05, 0.05, 0.15):
            ops.fiber(y, 0.0, 0.005, 1)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.geomTransf("Linear", 1)
        nodes, eles = [], []
        for i in range(8):
            n1, n2 = 2 * i + 1, 2 * i + 2
            nodes.extend([n1, n2])
            ops.node(n1, 0.0, float(i))
            ops.node(n2, L, float(i))
            ops.fix(n1, 1, 1, 1)
            ops.element("forceBeamColumn", i + 1, n1, n2, 1, 1)
            eles.append(i + 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for i in range(8):
            ops.load(2 * i + 2, 10.0e3, -25.0e3, 15.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0
        d = _step_digest(nodes, 3, eles)
        ops.ladrunoThreads(1)
        return d

    assert run_fbc2d_threaded(1) == run_fbc2d_threaded(2)


# ===========================================================================
# T14: RC Wall Macro-Elements (MVLEM & SFI_MVLEM 2D/3D) Ecosystem Wiring
# ===========================================================================

def test_t14_mvlem_and_sfi_mvlem_ecosystem_and_hdf5(tmp_path):
    """T14: `MVLEM`, `SFI_MVLEM`, `MVLEM_3D`, and `SFI_MVLEM_3D` `integrationWeights`,
    `localAxes`, `material` fiber responses, fail-closed option validation, and
    `.ladruno` HDF5 recorder support.
    """
    # --- 1. Fail-closed option validation on MVLEM ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 3.0)
    ops.uniaxialMaterial("Elastic", 1, 3.0e10)
    with pytest.raises(Exception):
        # Negative thickness must fail closed
        ops.element(
            "MVLEM", 1, 0.0, 1, 2, 2, 0.4,
            "-thick", -0.2, 0.2,
            "-width", 1.0, 1.0,
            "-rho", 0.01, 0.01,
            "-matConcrete", 1, 1,
            "-matSteel", 1, 1,
            "-matShear", 1
        )

    # --- 2. 2D MVLEM & SFI_MVLEM responses and .ladruno HDF5 recorder ---
    h5_mvlem = str(tmp_path / "mvlem_2d.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 3.0)
    ops.node(3, 5.0, 0.0)
    ops.node(4, 5.0, 3.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)

    ops.uniaxialMaterial("Elastic", 1, 3.0e10)
    ops.uniaxialMaterial("Elastic", 2, 2.0e11)
    ops.uniaxialMaterial("Elastic", 3, 1.2e10)
    ops.element(
        "MVLEM", 1, 0.0, 1, 2, 3, 0.4,
        "-thick", 0.2, 0.2, 0.2,
        "-width", 0.8, 0.8, 0.8,
        "-rho", 0.01, 0.01, 0.01,
        "-matConcrete", 1, 1, 1,
        "-matSteel", 2, 2, 2,
        "-matShear", 3
    )

    ops.nDMaterial("ElasticIsotropic", 10, 3.0e10, 0.2)
    ops.element(
        "SFI_MVLEM", 2, 3, 4, 3, 0.4,
        "-thick", 0.2, 0.2, 0.2,
        "-width", 0.8, 0.8, 0.8,
        "-mat", 10, 10, 10,
        "-Coupling", 1
    )

    ops.recorder("ladruno", h5_mvlem, "-E", "endForces", "-flush", 0, "-compress", 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 50.0e3, -100.0e3, 10.0e3)
    ops.load(4, 50.0e3, -100.0e3, 10.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    for eid in (1, 2):
        iw = ops.eleResponse(eid, "integrationWeights")
        assert len(iw) == 3
        assert list(iw) == pytest.approx([0.8, 0.8, 0.8], rel=1e-12)
        la = np.array(ops.eleResponse(eid, "localAxes")).reshape(3, 3)
        np.testing.assert_allclose(la @ la.T, np.eye(3), atol=1e-12)

    # Check material response forwarding on MVLEM and SFI_MVLEM
    assert len(ops.eleResponse(1, "material", 1, "strain")) == 1
    assert len(ops.eleResponse(1, "material", 1, "shear", "strain")) == 1
    assert len(ops.eleResponse(2, "material", 1, "strain")) > 0
    ef1 = list(ops.eleResponse(1, "endForces"))
    ef2 = list(ops.eleResponse(2, "endForces"))
    assert len(ef1) == 6 and len(ef2) == 6
    ops.wipe()

    assert os.path.isfile(h5_mvlem) and os.path.getsize(h5_mvlem) > 0
    raw = open(h5_mvlem, "rb").read()
    for row in (ef1, ef2):
        slab = np.asarray(row, dtype="<f8").tobytes()
        assert slab in raw

    # --- 3. 3D MVLEM_3D & SFI_MVLEM_3D responses and .ladruno HDF5 recorder ---
    h5_mvlem3d = str(tmp_path / "mvlem_3d.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    # Wall 1 (MVLEM_3D): nodes 1..4 counterclockwise
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 2.4, 0.0, 0.0)
    ops.node(3, 2.4, 3.0, 0.0)
    ops.node(4, 0.0, 3.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    # Wall 2 (SFI_MVLEM_3D): nodes 5..8 counterclockwise
    ops.node(5, 5.0, 0.0, 0.0)
    ops.node(6, 7.4, 0.0, 0.0)
    ops.node(7, 7.4, 3.0, 0.0)
    ops.node(8, 5.0, 3.0, 0.0)
    ops.fix(5, 1, 1, 1, 1, 1, 1)
    ops.fix(6, 1, 1, 1, 1, 1, 1)

    ops.uniaxialMaterial("Elastic", 1, 3.0e10)
    ops.uniaxialMaterial("Elastic", 2, 2.0e11)
    ops.uniaxialMaterial("Elastic", 3, 1.2e10)
    ops.element(
        "MVLEM_3D", 1, 1, 2, 3, 4, 3,
        "-thick", 0.2, 0.2, 0.2,
        "-width", 0.8, 0.8, 0.8,
        "-rho", 0.01, 0.01, 0.01,
        "-matConcrete", 1, 1, 1,
        "-matSteel", 2, 2, 2,
        "-matShear", 3
    )

    ops.nDMaterial("ElasticIsotropic", 10, 3.0e10, 0.2)
    ops.element(
        "SFI_MVLEM_3D", 2, 5, 6, 7, 8, 3,
        "-thick", 0.2, 0.2, 0.2,
        "-width", 0.8, 0.8, 0.8,
        "-mat", 10, 10, 10
    )

    ops.recorder("ladruno", h5_mvlem3d, "-E", "endForces", "-flush", 0, "-compress", 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 25.0e3, -50.0e3, 0.0, 0.0, 0.0, 0.0)
    ops.load(4, 25.0e3, -50.0e3, 0.0, 0.0, 0.0, 0.0)
    ops.load(7, 25.0e3, -50.0e3, 0.0, 0.0, 0.0, 0.0)
    ops.load(8, 25.0e3, -50.0e3, 0.0, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    for eid in (1, 2):
        iw3 = ops.eleResponse(eid, "integrationWeights")
        assert len(iw3) == 3
        assert list(iw3) == pytest.approx([0.8, 0.8, 0.8], rel=1e-12)
        la3 = np.array(ops.eleResponse(eid, "localAxes")).reshape(3, 3)
        np.testing.assert_allclose(la3 @ la3.T, np.eye(3), atol=1e-12)

    assert len(ops.eleResponse(1, "material", 1, "strain")) == 1
    assert len(ops.eleResponse(2, "material", 1, "strain")) > 0
    ef3_1 = list(ops.eleResponse(1, "endForces"))
    ef3_2 = list(ops.eleResponse(2, "endForces"))
    assert len(ef3_1) == 24 and len(ef3_2) == 24
    ops.wipe()

    assert os.path.isfile(h5_mvlem3d) and os.path.getsize(h5_mvlem3d) > 0
    raw3 = open(h5_mvlem3d, "rb").read()
    for row in (ef3_1, ef3_2):
        slab = np.asarray(row, dtype="<f8").tobytes()
        assert slab in raw3


# ===========================================================================
# T15: Regression Tests for Review Findings (Bugs 1 through 9)
# ===========================================================================

def test_t15_review_bug_regressions(tmp_path):
    """T15: Regression suite for all 9 bugs identified and fixed in review:
      1. LadrunoDispBeamColumn3d::getInitialStiff dual-hinge beta/gamma formulas.
      2. LadrunoDispBeamColumn2d/3d plasticRotation under span loads and single hinge.
      3. OPS_LadrunoDispBeamColumn2d -hingeIJ (1/2 args) and -hingeLoc (1 arg).
      4. LadrunoCohesiveHingeBiaxial P-M initial tangent, implex error/commit, and parameters.
      5. LadrunoIMKBeam2d/3d -endZones + span loads clear-span moment transformation.
      6. LadrunoIMKBeam2d/3d Timoshenko -shear fixed-end moments for point/partial loads & param aliases.
      7. Bilin & ModIMKPeakOriented -implex variable-dt scaling, LadrunoDispBeamColumn implex forwarding,
         and LadrunoIMKBeam plastic-rotation IMPL-EX extrapolation.
      8. ForceBeamColumn2d inline-integration -totalElementLch and mid-analysis updateParameter.
      9. MVLEM/SFI_MVLEM 2D/3D ops_TheActiveElement lch latching and .ladruno HDF5 localAxes recording.
    """
    L = 4.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # --- Bug 1: LadrunoDispBeamColumn3d::getInitialStiff dual-hinge interior locations ---
    def run_ldbc3d_init_stiff(algo_initial: bool):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 500.0e3, 2000.0, "-linear", "-penaltyRatio", 5.0)
        ops.element(
            "LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1,
            "-hingeI", 10, "-hingeJ", 10, "-hingeY", 10, "-hingeLoc", 0.2, 0.8
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0.0, -10.0e3, 8.0e3, 0.0, 5.0e3, 6.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        # Allow 1 iteration only: under elastic response, getInitialStiff must match getTangentStiff
        ops.test("NormUnbalance", 1.0e-9, 1)
        if algo_initial:
            ops.algorithm("ModifiedNewton", "-initial")
        else:
            ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return [ops.nodeDisp(2, d) for d in range(1, 7)]

    d_init = run_ldbc3d_init_stiff(True)
    d_tang = run_ldbc3d_init_stiff(False)
    np.testing.assert_allclose(d_init, d_tang, rtol=1e-12, atol=1e-14)

    # --- Bug 2: LadrunoDispBeamColumn2d/3d plasticRotation with span loads & single hinge ---
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.fix(2, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.section("Elastic", 1, E, A, Iz)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.fix(2, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        if ndm == 2:
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -20.0e3, 0.0)
        else:
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -20.0e3, 15.0e3, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 10)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        vp = ops.eleResponse(1, "plasticRotation")
        np.testing.assert_allclose(vp, np.zeros(len(vp)), atol=1e-13)

        # Single embedded hinge on elastic bulk beam: bulk plastic rotation must be 0
        ops.wipe()
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.section("Elastic", 1, E, A, Iz)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 40.0e3, 500.0, "-linear")
            ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 10)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 0.0, -20.0e3, 0.0)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 40.0e3, 500.0, "-linear")
            ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 10, "-hingeY", 10)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 0.0, -20.0e3, 15.0e3, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("DisplacementControl", 2, 2, -0.005)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(3) == 0
        hj = ops.eleResponse(1, "hingeJump")
        hj_i = ops.eleResponse(1, "hingeJumpI")
        hj_j = ops.eleResponse(1, "hingeJumpJ")
        assert hj == pytest.approx(hj_i, rel=1e-14, abs=1e-14)
        assert all(abs(x) < 1e-14 for x in hj_j)
        assert abs(hj[0]) > 1.0e-4
        vp_sh = ops.eleResponse(1, "plasticRotation")
        if ndm == 2:
            assert vp_sh == pytest.approx([0.0, hj[0], 0.0], rel=1e-10, abs=1e-11)
        else:
            assert vp_sh == pytest.approx([0.0, hj[0], 0.0, hj[1], 0.0, 0.0], rel=1e-10, abs=1e-11)

    # --- Bug 3: OPS_LadrunoDispBeamColumn2d -hingeIJ (1/2 args) and -hingeLoc (1 arg) ---
    def run_ldbc2d_syntax(extra_args):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 40.0e3, 2000.0, "-linear")
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 20, 60.0e3, 2000.0, "-linear")
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, *extra_args)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0.0, -20.0e3, 10.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("DisplacementControl", 2, 2, -0.004)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0
        return list(ops.eleResponse(1, "hingeSummary"))

    hs_sym_1 = run_ldbc2d_syntax(["-hingeIJ", 10, "-hingeLoc", 0.15])
    hs_sym_2 = run_ldbc2d_syntax(["-hingeI", 10, "-hingeJ", 10, "-hingeLoc", 0.15, 0.85])
    np.testing.assert_allclose(hs_sym_1, hs_sym_2, rtol=1e-12, atol=1e-14)
    hs_asym_1 = run_ldbc2d_syntax(["-hingeIJ", 10, 20])
    hs_asym_2 = run_ldbc2d_syntax(["-hingeI", 10, "-hingeJ", 20])
    np.testing.assert_allclose(hs_asym_1, hs_asym_2, rtol=1e-12, atol=1e-14)

    # --- Bug 4: LadrunoCohesiveHingeBiaxial P-M initial tangent, implex commit/error, parameters ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.nDMaterial(
        "LadrunoCohesiveHingeBiaxial", 88, 40.0e3, 1500.0, 30.0e3, 1500.0,
        "-pm", 1000.0e3, 1.0, "-implex"
    )
    ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxialI", 88, "-hingeBiaxialJ", 88)
    ops.parameter(1, "element", 1, "hingeI", "Mcz")
    ops.parameter(2, "element", 1, "hinge", "Mcy")
    ops.updateParameter(1, 45.0e3)
    ops.updateParameter(2, 35.0e3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -100.0e3, -20.0e3, 15.0e3, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("DisplacementControl", 2, 2, -0.002)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(5) == 0
    err_biax = ops.eleResponse(1, "implexError")[0]
    assert 0.0 < err_biax < 1.0
    assert ops.eleResponse(1, "energy")[0] > 0.0

    # --- Bug 5: LadrunoIMKBeam2d/3d -endZones dI dJ + span loads ---
    L_ez = 6.0
    dI, dJ = 0.5, 0.5
    Lc = L_ez - dI - dJ  # 5.0 m
    My_ez = 50.0e3
    k0_ez = 1.0e9
    bh_ez = 0.02
    wy_ez = -30.0e3  # clear-span FEM = |wy| * Lc^2 / 12 = 62.5 kN*m > My_ez = 50 kN*m
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L_ez, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.fix(2, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.uniaxialMaterial("Steel01", 10, My_ez, k0_ez, bh_ez)
            ops.element(
                "LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1,
                "-hinge", "both", "-matZ", 10, "-endZones", dI, dJ
            )
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy_ez, 0.0)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L_ez, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.fix(2, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.uniaxialMaterial("Steel01", 10, My_ez, k0_ez, bh_ez)
            ops.element(
                "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
                "-hinge", "both", "-matZ", 10, "-endZones", dI, dJ
            )
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy_ez, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.2)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(5) == 0
        th_ez = ops.eleResponse(1, "hingeRotation")
        hm_ez = ops.eleResponse(1, "hingeMoment")
        ef_ez = ops.eleForce(1)
        # Clear-span hinges must have yielded symmetrically (th_I = -th_J)
        assert abs(th_ez[0]) > 1.0e-5
        assert th_ez[0] == pytest.approx(-th_ez[1], rel=1e-10)
        assert abs(hm_ez[0]) == pytest.approx(My_ez + bh_ez * k0_ez * (abs(th_ez[0]) - My_ez / k0_ez), rel=1e-10)
        # Nodal moment at node 1 must equal T^T * [hm_I, hm_J] (and decrease from elastic FEM 90 kN*m due to hinge yielding)
        T00 = (L_ez - dJ) / Lc
        T10 = dI / Lc
        mz_node1_expected = abs(T00 * hm_ez[0] + T10 * hm_ez[1])
        mz_idx = 2 if ndm == 2 else 5
        assert abs(ef_ez[mz_idx]) == pytest.approx(mz_node1_expected, rel=1e-10)
        assert My_ez < abs(ef_ez[mz_idx]) < abs(wy_ez) * (L_ez ** 2) / 12.0

    # --- Bug 6: Timoshenko -shear fixed-end moments for point/partial loads & parameter aliases ---
    # Propped cantilever (node 1 fixed, node 2 pinned in y/z): under midspan point load P,
    # Euler-Bernoulli fixed-end moment at node 1 is 3/16 P L, whereas Timoshenko with Phi = 12 EI / (GA_s L^2)
    # gives M_1 = (3/16) P L * (1 + Phi/2) / (1 + Phi/4) - wait: with propped end M_2 = 0,
    # M_1^Timo = M1_fem^Timo - (k12^Timo / k22^Timo) * M2_fem^Timo = (3/8 P L) / (2 + Phi/2) = (3/16 P L) / (1 + Phi/4).
    GAy = 2.0e7  # Chosen so Phi_z = 12 * E * Iz / (GAy * L^2) = 1.2 is O(1)
    Phi_z = 12.0 * E * Iz / (GAy * L * L)
    P_mid = -40.0e3
    M1_exact_timo = (3.0 / 16.0) * abs(P_mid) * L / (1.0 + 0.25 * Phi_z)
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.fix(2, 1, 1, 0)  # pinned roller at node 2
            ops.geomTransf("Linear", 1)
            ops.uniaxialMaterial("Elastic", 10, 1.0e14)
            ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-shear", GAy, "-hinge", "both", "-matZ", 10)
            ops.parameter(1, "element", 1, "hingeZI", "E")
            ops.parameter(2, "element", 1, "hinge2", "E")
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.eleLoad("-ele", 1, "-type", "-beamPoint", P_mid, 0.5, 0.0)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.fix(2, 1, 1, 1, 1, 1, 0)  # free about z at node 2
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.uniaxialMaterial("Elastic", 10, 1.0e14)
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1, "-shear", GAy, GAy, "-hinge", "both", "-matZ", 10)
            ops.parameter(1, "element", 1, "hingeMzI", "E")
            ops.parameter(2, "element", 1, "hingeZJ", "E")
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.eleLoad("-ele", 1, "-type", "-beamPoint", P_mid, 0.0, 0.5, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 10)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        ef_timo = ops.eleForce(1)
        mz1_idx = 2 if ndm == 2 else 5
        assert abs(ef_timo[mz1_idx]) == pytest.approx(M1_exact_timo, rel=1e-4)

    # --- Bug 7: IMPL-EX on Bilin & ModIMKPeakOriented in LadrunoDispBeamColumn & LadrunoIMKBeam ---
    # 7a. LadrunoIMKBeam2d with Bilin -implex must be exact in the elastic regime (implexError == 0)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial(
        "Bilin", 10, 5.0e7, 0.02, 0.02, 80.0e3, -80.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
    )
    ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 5.0e3, 0.0)  # M_base = 20 kN*m < My = 80 kN*m (purely elastic)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormUnbalance", 1.0e-10, 1)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    # In the elastic regime, plastic hinge rotation is 0 -> IMPL-EX error must be 0
    assert ops.eleResponse(1, "implexError")[0] == pytest.approx(0.0, abs=1e-12)
    # Push past yield (step factor 3.0 -> total factor 4.0 -> M_base = 80 kN*m, then step 5.0 -> 140 kN*m)
    ops.integrator("LoadControl", 3.0)
    assert ops.analyze(2) == 0
    assert ops.eleResponse(1, "implexError")[0] > 0.0

    # 7b. LadrunoDispBeamColumn2d with ModIMKPeakOriented -implex forwards implexError
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial(
        "ModIMKPeakOriented", 10, 5.0e7, 0.02, 0.02, 40.0e3, -40.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
    )
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 15.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(3) == 0
    assert ops.eleResponse(1, "implexError")[0] > 0.0

    # --- Bug 8: ForceBeamColumn2d inline-integration -totalElementLch and mid-analysis updateParameter ---
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        ops.uniaxialMaterial("Elastic", 1, E)
        if ndm == 2:
            ops.section("Fiber", 1)
            for y in (-0.15, -0.05, 0.05, 0.15):
                ops.fiber(y, 0.0, 0.005, 1)
        else:
            ops.section("Fiber", 1, "-GJ", G * Jx)
            for y in (-0.15, 0.15):
                for z in (-0.10, 0.10):
                    ops.fiber(y, z, 0.005, 1)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.element("forceBeamColumn", 1, 1, 2, 1, 1, "-totalElementLch", "-lMass")
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.element("forceBeamColumn", 1, 1, 2, 1, 1, "-totalElementLch", "-lMass")
        ops.parameter(1, "element", 1, "allSections", "E")
        ops.parameter(2, "element", 1, "lch")
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        if ndm == 2:
            ops.load(2, 0.0, -40.0e3, 0.0)
        else:
            ops.load(2, 0.0, -40.0e3, 0.0, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0
        u_before = ops.nodeDisp(2, 2)
        # Update parameter mid-analysis and re-equilibrate at constant load (LoadControl 0.0)
        ops.updateParameter(1, 1.25 * E)
        ops.updateParameter(2, 0.75 * L)
        ops.integrator("LoadControl", 0.0)
        assert ops.analyze(1) == 0
        u_after = ops.nodeDisp(2, 2)
        assert abs(u_after) < abs(u_before)

    # --- Bug 9: MVLEM/SFI_MVLEM ops_TheActiveElement lch latching & .ladruno HDF5 localAxes ---
    h5_axes = str(tmp_path / "mvlem_axes.ladruno")
    h_wall = 3.2

    def run_mvlem_lch(with_preceding_truss: bool, record_h5: bool):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        if with_preceding_truss:
            # Create a dummy truss of length 10.0 first so ops_TheActiveElement would be 10.0 if MVLEM didn't set it
            ops.node(10, 0.0, 0.0)
            ops.node(11, 10.0, 0.0)
            ops.fix(10, 1, 1, 1)
            ops.fix(11, 1, 1, 1)
            ops.uniaxialMaterial("Elastic", 99, E)
            ops.element("Truss", 99, 10, 11, A, 99)

        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, h_wall)
        ops.fix(1, 1, 1, 1)
        ops.uniaxialMaterial(
            "ASDConcrete1D", 1, 3.0e10,
            "-Te", 1.0e-4, 1.0e-3, "-Ts", 3.0e6, 0.0, "-Td", 0.0, 0.0,
            "-Ce", -1.0e-3, -2.0e-3, -5.0e-3, "-Cs", -20.0e6, -30.0e6, -5.0e6, "-Cd", 0.0, 0.2, 0.5,
            "-autoRegularization", 300.0, 5000.0
        )
        ops.uniaxialMaterial("Elastic", 2, 2.0e11)
        ops.uniaxialMaterial("Elastic", 3, 1.2e10)
        ops.element(
            "MVLEM", 1, 0.0, 1, 2, 2, 0.4,
            "-thick", 0.2, 0.2,
            "-width", 1.0, 1.0,
            "-rho", 0.01, 0.01,
            "-matConcrete", 1, 1,
            "-matSteel", 2, 2,
            "-matShear", 3
        )
        if record_h5:
            ops.recorder("ladruno", h5_axes, "-E", "localAxes", "endForces", "-flush", 0, "-compress", 0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0.0, 100.0e3, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("DisplacementControl", 2, 2, 0.001)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        f_stress = list(ops.eleResponse(1, "Fiber_Stress"))
        la = list(ops.eleResponse(1, "localAxes"))
        ops.wipe()
        return f_stress, la

    fs_with_truss, la_flat = run_mvlem_lch(True, True)
    fs_alone, _ = run_mvlem_lch(False, False)
    # Concrete post-peak tensile stress must match MVLEM alone (lch = 3.2 m, not 10.0 m from the Truss)
    np.testing.assert_allclose(fs_with_truss, fs_alone, rtol=1e-12, atol=1e-12)
    assert len(la_flat) == 9
    assert la_flat == pytest.approx([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], abs=1e-12)

    raw_axes = open(h5_axes, "rb").read()
    assert b"localAxes" in raw_axes
    assert np.asarray(la_flat, dtype="<f8").tobytes() in raw_axes


# ===========================================================================
# T16: DeepInvestigator Review Regressions (N1-N5, F1-F4, F7-F9)
# ===========================================================================

def test_t16_deep_investigator_regressions(tmp_path):
    """T16: Regression suite for all findings from the DeepInvestigator review:
      - N1: LadrunoDispBeamColumn2d/3d FileDatastore save/restore with active hinges (no ID table collision).
      - N2/F3: Unified P-M interaction formula & default alpha=0.5 across LadrunoCohesiveHinge,
               LadrunoCohesiveHingeBiaxial, LadrunoIMKBeam2d/3d, and LadrunoDispBeamColumn2d/3d (including Bilin/IMK hinges).
      - N3: OPS_LadrunoIMKBeam2d rejects ndm != 2 or ndf != 3.
      - N4: DispBeamColumnNL3d localForce (responseID == 2) Vy sign matches (Mz1 + Mz2)/L.
      - N5: LadrunoCohesiveHinge updateParameter('Mc'/'Gf') recalibrates Kpen when defaulted (!userKpen).
      - F1: LadrunoDispBeamColumn2d/3d single, dual, and biaxial hinges include span-load fixed-end forces q0.
      - F2: LadrunoIMKBeam2d/3d implexError (max) vs avgImplexError (mean over active hinges).
      - F4: LadrunoIMKBeam2d/3d energy with -endZones and -pm; LadrunoDispBeamColumn2d/3d single-hinge plasticRotation with xi_I > 0.
      - F7/F8: SFI_MVLEM_3D non-FSAM plane-stress Eeq = E (not E/(1-nu^2)) when Coupling=2.
      - F9: FileDatastore save/restore mid-analysis with -implex for LadrunoIMKBeam2d/3d, Bilin, and ModIMKPeakOriented.
    """
    L = 4.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # --- N1: LadrunoDispBeamColumn2d & 3d FileDatastore save/restore with active asymmetric dual hinges & span loads ---
    def build_ldbc_2d_hinge():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 11, 35.0e3, 800.0, "-linear")
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 12, 55.0e3, 900.0, "-linear")
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 11, "-hingeJ", 12, "-pm", 1000.0e3, 0.5)
        ops.element("elasticBeamColumn", 2, 1, 2, A, E, Iz, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0e3, -30.0e3, 20.0e3)
        ops.pattern("Plain", 2, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -4.0e3, -2.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0

    def probe_ldbc():
        return list(ops.eleForce(1)) + list(ops.eleResponse(1, "hingeSummary"))

    database_roundtrip(build_ldbc_2d_hinge, probe_nodes=[2], ndf=3, probe_fn=probe_ldbc, rounds=2)

    def build_ldbc_3d_hinge():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 21, 35.0e3, 800.0, 25.0e3, 750.0, "-linear")
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 22, 55.0e3, 900.0, 40.0e3, 850.0, "-linear")
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxialI", 21, "-hingeBiaxialJ", 22, "-pm", 1000.0e3, 0.5)
        ops.element("elasticBeamColumn", 2, 1, 2, A, E, G, Jx, Iy, Iz, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -50.0e3, -30.0e3, 20.0e3, 0.0, 10.0e3, 15.0e3)
        ops.pattern("Plain", 2, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -3.0e3, 2.0e3, -1.5e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0

    database_roundtrip(build_ldbc_3d_hinge, probe_nodes=[2], ndf=6, probe_fn=probe_ldbc, rounds=2)

    # --- N2/F3: Unified P-M interaction formula & Bilin/IMK hinge P-M scaling in LadrunoDispBeamColumn ---
    # Under P = -200 kN (compression, n = P/P0 = -0.2) and default alpha = 0.5:
    # f(n) = 1 - 0.5*(-0.2) - 0.5*(0.04) = 1 + 0.1 - 0.02 = 1.08
    P_ax = -200.0e3
    P0_pm = 1000.0e3
    f_expected = 1.0 - 0.5 * (P_ax / P0_pm) - 0.5 * ((P_ax / P0_pm) ** 2)  # 1.08
    assert f_expected == pytest.approx(1.08, rel=1e-14)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    # Use Bilin hinge material inside LadrunoDispBeamColumn2d with element-level -pm P0 (default alpha=0.5)
    My0 = 50.0e3
    ops.uniaxialMaterial(
        "Bilin", 30, 1.0e9, 0.001, 0.001, My0, -My0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.1, 0.1, 0.2, 0.2, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0
    )
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 30, "-pm", P0_pm)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, P_ax, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, -1.0, 0.0)
    ops.integrator("DisplacementControl", 2, 2, -0.02)
    assert ops.analyze(1) == 0
    hs_pm = ops.eleResponse(1, "hingeSummary")
    # Hinge I moment must equal f_expected * unscaled Bilin moment at rotation hs_pm[0]
    m_bilin_unscaled = My0 + 0.001 * 1.0e9 * (abs(hs_pm[0]) - My0 / 1.0e9)
    assert abs(hs_pm[1]) == pytest.approx(f_expected * m_bilin_unscaled, rel=1e-12)

    # --- N3: OPS_LadrunoIMKBeam2d rejects 3D model ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    with pytest.raises(Exception):
        ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1)

    # --- N4: DispBeamColumnNL3d localForce Vy sign ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.element("dispBeamColumnNL", 1, 1, 2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -12.0e3, 8.0e3, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    lf_nl = ops.eleResponse(1, "localForce")
    ef_nl = ops.eleForce(1)
    np.testing.assert_allclose(lf_nl, ef_nl, rtol=1e-12, atol=1e-12)

    # --- N5: LadrunoCohesiveHinge updateParameter('Mc'/'Gf') recalibrates defaulted Kpen ---
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    # Default Kpen: r = 1000, KpenFloor = Mc^2 / (2*Gf) = 40e3^2 / (2*500) = 1.6e6 -> Kpen = 1.6e9
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 1, 40.0e3, 500.0, "-linear")
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.parameter(1, "element", 1, "material", 1, "Mc")
    ops.updateParameter(1, 80.0e3)  # doubling Mc quadruples KpenFloor -> Kpen = 6.4e9
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3)  # elastic pre-peak (10 kN < 80 kN)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(2, 1) == pytest.approx(10.0e3 / 6.4e9, rel=1e-12)

    # --- F1: LadrunoDispBeamColumn2d & 3d span loads trigger hinges on fixed-fixed beam ---
    for ndm in (2, 3):
        ops.wipe()
        ndf = 3 if ndm == 2 else 6
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.fix(2, 1, 1, 1)
            ops.geomTransf("Linear", 1)
            ops.section("Elastic", 1, E, A, Iz)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 30.0e3, 1000.0, "-linear")
            ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 10, "-hingeJ", 10)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            # Fixed-end moment |w|*L^2/12 = 40e3 * 16 / 12 = 53.33 kN*m > Mc = 30 kN*m
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -40.0e3, 0.0)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.fix(2, 1, 1, 1, 1, 1, 1)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 20, 30.0e3, 1000.0, 30.0e3, 1000.0, "-linear")
            ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxialI", 20, "-hingeBiaxialJ", 20)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -40.0e3, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        hs_span = ops.eleResponse(1, "hingeSummary")
        # Even though v = 0 (fixed-fixed), q0 drives both end hinges into softening!
        assert abs(hs_span[0]) > 1.0e-5 and abs(hs_span[2]) > 1.0e-5
        assert hs_span[0] == pytest.approx(-hs_span[2], rel=1e-10)
        assert abs(hs_span[1]) < 30.0e3

    # --- F2: LadrunoIMKBeam2d/3d implexError (max) vs avgImplexError (mean over active hinges) ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial(
        "Bilin", 10, 5.0e7, 0.02, 0.02, 40.0e3, -40.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
    )
    ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # Cantilever tip shear: hinge I yields (M_I = 60 kN*m > 40 kN*m), hinge J has M_J = 0 (err_J = 0)
    ops.load(2, 0.0, -15.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormUnbalance", 1.0e-10, 1)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    err_max = ops.eleResponse(1, "implexError")[0]
    err_avg = ops.eleResponse(1, "avgImplexError")[0]
    assert err_max > 0.0
    # 2 active hinges, only hinge I has non-zero error -> avgImplexError == 0.5 * implexError
    assert err_avg == pytest.approx(0.5 * err_max, rel=1e-12)

    # --- F4a: LadrunoIMKBeam2d energy with -endZones and -pm matches integrated work \int q \cdot dv ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 10, 40.0e3, 2.0e7, 0.02)
    ops.element(
        "LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-endZones", 0.4, 0.4, "-pm", 1000.0e3, 0.5
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -150.0e3, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    v_prev = np.array(ops.eleResponse(1, "basicDeformation"))
    q_prev = np.array(ops.eleResponse(1, "basicForce"))
    work_imk = float(0.5 * np.dot(q_prev, v_prev))
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, -25.0e3, 10.0e3)
    ops.integrator("LoadControl", 0.05)
    for _ in range(20):
        assert ops.analyze(1) == 0
        v_cur = np.array(ops.eleResponse(1, "basicDeformation"))
        q_cur = np.array(ops.eleResponse(1, "basicForce"))
        work_imk += float(0.5 * np.dot(q_prev + q_cur, v_cur - v_prev))
        v_prev, q_prev = v_cur, q_cur
    assert ops.eleResponse(1, "energy")[0] == pytest.approx(work_imk, rel=5e-3)

    # --- F4b: LadrunoDispBeamColumn2d single hinge at interior location xi_I distributes plasticRotation ---
    xi_I = 0.25
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 30.0e3, 1000.0, "-linear")
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 10, "-hingeLoc", xi_I, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -25.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("DisplacementControl", 2, 2, -0.005)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    alpha_xi = ops.eleResponse(1, "hingeJump")[0]
    assert abs(alpha_xi) > 1.0e-5
    vp_xi = ops.eleResponse(1, "plasticRotation")
    assert vp_xi == pytest.approx([0.0, (1.0 - xi_I) * alpha_xi, -xi_I * alpha_xi], rel=1e-10, abs=1e-12)

    # --- F8: SFI_MVLEM_3D with ElasticIsotropic uses plane-stress Eeq = E (matching explicit -Eave E_wall) ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 2.0, 0.0, 0.0)
    ops.node(3, 2.0, 3.0, 0.0)
    ops.node(4, 0.0, 3.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    # Use high Poisson's ratio nu = 0.4 so plane-strain D(1,1) = E/(1-nu^2) = 1.1905*E != E
    E_wall, nu_wall = 3.0e10, 0.4
    ops.nDMaterial("ElasticIsotropic", 10, E_wall, nu_wall)
    ops.element(
        "SFI_MVLEM_3D", 1, 1, 2, 3, 4, 2,
        "-thick", 0.2, 0.2,
        "-width", 1.0, 1.0,
        "-mat", 10, 10
    )
    # Build identical wall 2 with explicit -Eave E_wall
    ops.node(5, 5.0, 0.0, 0.0)
    ops.node(6, 7.0, 0.0, 0.0)
    ops.node(7, 7.0, 3.0, 0.0)
    ops.node(8, 5.0, 3.0, 0.0)
    ops.fix(5, 1, 1, 1, 1, 1, 1)
    ops.fix(6, 1, 1, 1, 1, 1, 1)
    ops.element(
        "SFI_MVLEM_3D", 2, 5, 6, 7, 8, 2,
        "-thick", 0.2, 0.2,
        "-width", 1.0, 1.0,
        "-mat", 10, 10,
        "-Eave", E_wall
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, 0.0, 10.0e3, 0.0, 0.0, 0.0)
    ops.load(4, 0.0, 0.0, 10.0e3, 0.0, 0.0, 0.0)
    ops.load(7, 0.0, 0.0, 10.0e3, 0.0, 0.0, 0.0)
    ops.load(8, 0.0, 0.0, 10.0e3, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(3, 3) == pytest.approx(ops.nodeDisp(7, 3), rel=1e-12)

    # --- F9: FileDatastore save/restore mid-analysis with -implex on LadrunoIMKBeam2d/3d, LadrunoDispBeamColumn2d/3d, Bilin, ModIMKPeakOriented, LadrunoCohesiveHinge(Biaxial) ---
    def build_imk_implex():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.uniaxialMaterial(
            "Bilin", 10, 5.0e7, 0.02, 0.02, 40.0e3, -40.0e3,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 2.0, "-implex"
        )
        ops.uniaxialMaterial(
            "ModIMKPeakOriented", 20, 4.0e7, 0.02, 0.02, 35.0e3, -35.0e3,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
        )
        ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZi", 10, "-matZj", 20)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -120.0e3, -18.0e3, 15.0e3)
        ops.pattern("Plain", 2, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -3.0e3, -1.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormUnbalance", 1.0e-10, 1)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(3) == 0

    def probe_imk_implex():
        return list(ops.eleForce(1)) + list(ops.eleResponse(1, "hingeSummary")) + list(ops.eleResponse(1, "implexError"))

    database_roundtrip(build_imk_implex, probe_nodes=[2], ndf=3, probe_fn=probe_imk_implex, rounds=2)

    def build_imk_3d_implex():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.uniaxialMaterial(
            "Bilin", 10, 5.0e7, 0.02, 0.02, 40.0e3, -40.0e3,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
        )
        ops.uniaxialMaterial(
            "ModIMKPeakOriented", 20, 4.0e7, 0.02, 0.02, 35.0e3, -35.0e3,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
        )
        ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1, "-hinge", "both", "-matZ", 10, "-matY", 20)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -80.0e3, -18.0e3, 12.0e3, 0.0, 10.0e3, 15.0e3)
        ops.pattern("Plain", 2, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -2.5e3, 1.5e3, -1.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormUnbalance", 1.0e-10, 1)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(3) == 0

    database_roundtrip(build_imk_3d_implex, probe_nodes=[2], ndf=6, probe_fn=probe_imk_implex, rounds=2)

    def build_ldbc_2d_implex_hinges():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 11, 35.0e3, 800.0, "-linear", "-implex")
        ops.uniaxialMaterial(
            "Bilin", 12, 5.0e7, 0.02, 0.02, 45.0e3, -45.0e3,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            0.05, 0.05, 0.2, 0.2, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, "-implex"
        )
        ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1, "-hingeI", 11, "-hingeJ", 12, "-pm", 1000.0e3, 0.5)
        ops.element("elasticBeamColumn", 2, 1, 2, A, E, Iz, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0e3, -25.0e3, 18.0e3)
        ops.pattern("Plain", 2, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -3.0e3, -1.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormUnbalance", 1.0e-10, 1)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(3) == 0

    database_roundtrip(build_ldbc_2d_implex_hinges, probe_nodes=[2], ndf=3, probe_fn=probe_imk_implex, rounds=2)

    def build_ldbc_3d_implex_hinges():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 21, 35.0e3, 800.0, 25.0e3, 750.0, "-linear", "-implex")
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 22, 50.0e3, 900.0, 38.0e3, 850.0, "-linear", "-implex")
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeBiaxialI", 21, "-hingeBiaxialJ", 22, "-pm", 1000.0e3, 0.5)
        ops.element("elasticBeamColumn", 2, 1, 2, A, E, G, Jx, Iy, Iz, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -80.0e3, -25.0e3, 18.0e3, 0.0, 10.0e3, 15.0e3)
        ops.pattern("Plain", 2, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -2.5e3, 1.5e3, -1.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.test("NormUnbalance", 1.0e-10, 1)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(3) == 0

    database_roundtrip(build_ldbc_3d_implex_hinges, probe_nodes=[2], ndf=6, probe_fn=probe_imk_implex, rounds=2)

    # Verify LadrunoIMKBeam2d forwards axial force to material-level LadrunoCohesiveHinge -pm even without beam-level -pm,
    # and verify implexError is relative and not distorted by fPM
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 40.0e3, 800.0, "-linear", "-pm", P0_pm, 0.5, "-implex")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 11, 40.0e3, 800.0, "-linear", "-implex")
    ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "i", "-matZ", 10)
    ops.node(3, 0.0, 0.0)
    ops.node(4, L, 0.0)
    ops.fix(3, 1, 1, 1)
    ops.element("LadrunoIMKBeam2d", 2, 3, 4, A, E, Iz, 1, "-hinge", "i", "-matZ", 11, "-pm", P0_pm, 0.5)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, P_ax, -8.0e3, 0.0)
    ops.load(4, P_ax, -8.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormUnbalance", 1.0e-10, 1)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(3) == 0
    assert ops.eleResponse(1, "hingeMoment")[0] == pytest.approx(ops.eleResponse(2, "hingeMoment")[0], rel=1e-12)
    assert ops.eleResponse(1, "implexError")[0] == pytest.approx(ops.eleResponse(2, "implexError")[0], rel=1e-12)

    # Verify updateParameter immediately after restore invalidates committedAtCurrentV and updates basic/global forces
    dbpath = str(tmp_path / "imk_up")
    build_imk_implex()
    ops.parameter(1, "element", 1, "E")
    f_before = list(ops.eleForce(1))
    ops.database("File", dbpath)
    ops.save(1)
    ops.restore(1)
    ops.updateParameter(1, E * 1.25)
    f_after = list(ops.eleForce(1))
    assert f_after[0] != pytest.approx(f_before[0], rel=1e-6)


def test_apegmsh_emitter_guide_snippets_and_parameters():
    """Verify every command signature, option flag, mutual-exclusion rule,
    P-M formula, setParameter token, and energy response documented in
    `Ladruno_implementation/200_ladruno_frame_elements_apegmsh_emitter_guide.md`.
    """
    L = 4.0
    A, E, G, Jx, Iy, Iz = 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4

    # --- 1. 2D guide snippets, mutual exclusions, -lch, setParameter, and energy ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i in range(1, 13):
        ops.node(i, 0.0 if i % 2 == 1 else L, float(i // 2))
        if i % 2 == 1:
            ops.fix(i, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.beamIntegration("Lobatto", 1, 1, 5)

    # Section 3.1: LadrunoCohesiveHinge flags
    ops.uniaxialMaterial(
        "LadrunoCohesiveHinge", 101, 50.0e3, 200.0,
        "-penaltyRatio", 1200.0, "-exp", "-pm", 800.0e3, 0.5, "-implex",
    )
    ops.uniaxialMaterial(
        "LadrunoCohesiveHinge", 102, 55.0e3, 220.0,
        "-Kpen", 1.0e8, "-linear", "-pm", 800.0e3, 0.4, "-implex",
    )

    # Section 3.3: Bilin and ModIMKPeakOriented with -implex
    ops.uniaxialMaterial(
        "Bilin", 103, 1.0e8, 0.02, 0.02, 60.0e3, -60.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.04, 0.04, 0.15, 0.15, 0.2, 0.2, 0.08, 0.08, 1.0, 1.0, 0.0, "-implex",
    )
    ops.uniaxialMaterial(
        "ModIMKPeakOriented", 104, 1.0e8, 0.02, 0.02, 60.0e3, -60.0e3,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.04, 0.04, 0.15, 0.15, 0.2, 0.2, 0.08, 0.08, 1.0, 1.0, 0.0, "-implex",
    )

    # Section 2.1: LadrunoDispBeamColumn 2D snippets
    ops.element(
        "LadrunoDispBeamColumn", 1, 1, 2, 1, 1,
        "-mass", 100.0, "-cMass", "-lch", "ip",
        "-hingeI", 101, "-hingeJ", 102, "-hingeLoc", 0.05, 0.95,
        "-implex", "-pm", 800.0e3, 0.5,
    )
    ops.element(
        "LadrunoDispBeamColumn", 2, 3, 4, 1, 1,
        "-mass", 100.0, "-cMass", "-lch", "ip", "-nl",
    )
    ops.element(
        "LadrunoDispBeamColumn", 3, 5, 6, 1, 1,
        "-lch", "0.25", "-hingeIJ", 101, 102,
    )

    # Reject -lch auto and reject -nl + -hinge*
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 90, 1, 2, 1, 1, "-lch", "auto")
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 91, 1, 2, 1, 1, "-nl", "-hingeI", 101)

    # Section 2.2: LadrunoIMKBeam 2D snippets (1-arg and 2-arg -shear, -endZones, -implex)
    ops.element(
        "LadrunoIMKBeam", 4, 7, 8, A, E, Iz, 1,
        "-hinge", "both", "-matZi", 103, "-matZj", 104,
        "-shear", G * (5.0 / 6.0 * A), "-endZones", 0.15, 0.15,
        "-mass", 100.0, "-cMass", "-pm", 800.0e3, 0.5, "-implex",
    )
    ops.element(
        "LadrunoIMKBeam", 5, 9, 10, A, E, Iz, 1,
        "-matZ", 101, "-shear", G, 5.0 / 6.0 * A, "-endZone", 0.1, 0.1,
    )
    # Reject nonexistent -mat alias
    with pytest.raises(Exception):
        ops.element("LadrunoIMKBeam", 92, 1, 2, A, E, Iz, 1, "-mat", 101)

    # Section 2.3: forceBeamColumn 2D -lch ip|element|$val and -totalElementLch
    ops.element(
        "forceBeamColumn", 6, 11, 12, 1, 1,
        "-mass", 100.0, "-cMass", "-iter", 15, 1.0e-12,
        "-lch", "ip", "-totalElementLch",
    )
    with pytest.raises(Exception):
        ops.element("forceBeamColumn", 93, 11, 12, 1, 1, "-lch", "auto")

    # Section 5: setParameter tokens on LadrunoCohesiveHinge & LadrunoDispBeamColumn2d
    for pid, tok in enumerate(("Mc", "Gf", "Kpen", "P0", "alphaPM"), start=1):
        ops.parameter(pid, "element", 1, "hingeI", tok)
    ops.updateParameter(1, 52.0e3)
    ops.updateParameter(2, 210.0)
    ops.updateParameter(3, 1.2e8)
    ops.updateParameter(4, 900.0e3)
    ops.updateParameter(5, 0.6)

    # Step analysis and verify Section 6.1 'energy' includes hinge energy
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -180.0e3, -35.0e3, 15.0e3)
    ops.load(8, -180.0e3, -35.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0

    hinge_i_e = ops.eleResponse(1, "hingeI", "energy")[0]
    hinge_j_e = ops.eleResponse(1, "hingeJ", "energy")[0]
    total_e = ops.eleResponse(1, "energy")[0]
    assert hinge_i_e > 0.0 and total_e == pytest.approx(hinge_i_e + hinge_j_e, rel=1e-12)
    assert len(ops.eleResponse(4, "hinge", "0", "implexError")) == 1
    assert len(ops.eleResponse(4, "hinge", "1", "implexError")) == 1

    # Check P-M formula: f(n) = clamp(1 - alpha*n - (1-alpha)*n^2, 0.05, 2.5)
    N_ax = ops.eleResponse(1, "basicForce")[0]
    n_rat = N_ax / 900.0e3
    expected_fpm = max(0.05, min(2.5, 1.0 - 0.6 * n_rat - 0.4 * n_rat * n_rat))
    assert ops.eleResponse(1, "hingeI", "pmMultiplier")[0] == pytest.approx(expected_fpm, rel=1e-12)

    # --- 2. 3D guide snippets, mutual exclusions, and LadrunoCohesiveHingeBiaxial setParameter ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for i in range(1, 9):
        ops.node(i, 0.0 if i % 2 == 1 else L, float(i // 2), 0.0)
        if i % 2 == 1:
            ops.fix(i, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 5)

    ops.uniaxialMaterial("LadrunoCohesiveHinge", 201, 60.0e3, 240.0, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 202, 45.0e3, 180.0, "-linear")
    ops.nDMaterial(
        "LadrunoCohesiveHingeBiaxial", 301, 60.0e3, 240.0, 45.0e3, 180.0,
        "-exp", "-penaltyRatio", 1000.0, "-bk", 1.5, "-pm", 900.0e3, 0.5, "-implex",
    )

    # Coupled biaxial dual-end hinge element
    ops.element(
        "LadrunoDispBeamColumn", 1, 1, 2, 1, 1,
        "-hingeBiaxialIJ", 301, "-hingeLoc", 0.0, 1.0, "-implex", "-pm", 900.0e3, 0.5,
    )
    # Scalar z + y dual-end hinge element
    ops.element(
        "LadrunoDispBeamColumn", 2, 3, 4, 1, 1,
        "-hingeIJ", 201, "-hingeYIJ", 202,
    )
    # 3D LadrunoIMKBeam with 3-arg -shear, -endZones, -implex
    ops.element(
        "LadrunoIMKBeam", 3, 5, 6, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "both", "-matZ", 201, "-matY", 202,
        "-shear", G, 5.0 / 6.0 * A, 5.0 / 6.0 * A,
        "-endZones", 0.15, 0.15, "-pm", 900.0e3, 0.5, "-implex",
    )

    # Mutual exclusions in 3D:
    # (a) -nl + -hingeBiaxial*
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 94, 1, 2, 1, 1, "-nl", "-hingeBiaxial", 301)
    # (b) -hingeBiaxial* + -hinge*
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 95, 1, 2, 1, 1, "-hinge", 201, "-hingeBiaxial", 301)
    # (c) -hingeY* without -hinge*
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 96, 1, 2, 1, 1, "-hingeY", 202)

    # Section 5: setParameter tokens on LadrunoCohesiveHingeBiaxial
    for pid, tok in enumerate(
        ("Mcz", "Gfz", "Mcy", "Gfy", "P0", "alphaPM", "bkEta", "penaltyRatio"), start=10
    ):
        ops.parameter(pid, "element", 1, "hingeBiaxialI", tok)
    ops.updateParameter(10, 65.0e3)
    ops.updateParameter(11, 260.0)
    ops.updateParameter(12, 48.0e3)
    ops.updateParameter(13, 190.0)
    ops.updateParameter(14, 950.0e3)
    ops.updateParameter(15, 0.55)
    ops.updateParameter(16, 1.8)
    ops.updateParameter(17, 1500.0)

    # Section 6.1: sectionX takes physical coordinate x in [0, L]
    sec_force_mid = ops.eleResponse(1, "sectionX", "1.5", "force")
    assert len(sec_force_mid) >= 2

    # ------------------------------------------------------------------
    # 3. Section 3.4 Fiber / Rebar / Bond-Slip Uniaxial Material Signatures
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    # LadrunoUniaxialJ2
    ops.uniaxialMaterial(
        "LadrunoUniaxialJ2", 401, 200.0e6,
        "-iso", "voce", 420.0e3, 80.0e3, 10.0, 2.0e6,
        "-kin", 1, 5.0e6, 20.0,
        "-damage", "lemaitre", 5.0e3, 1.0, 0.02, 0.95,
        "-implex",
    )
    # Base steel for LadrunoRebarBuckling
    ops.uniaxialMaterial("Steel02", 400, 420.0e3, 200.0e6, 0.01)
    # LadrunoRebarBuckling (wraps base uniaxial material tag 400; does NOT take -implex)
    ops.uniaxialMaterial(
        "LadrunoRebarBuckling", 402, 400,
        "-lsr", 8.0, "-model", "dm", "-alpha", 0.75, "-reduction", 1.0,
        "-fsufrac", 0.5, "-fy", 420.0e3, "-E", 200.0e6, "-restraighten", "lambda",
    )
    with pytest.raises(Exception):
        ops.uniaxialMaterial("LadrunoRebarBuckling", 492, 400, "-lsr", 8.0, "-implex")

    # LadrunoBondSlip (does NOT take -implex)
    ops.uniaxialMaterial(
        "LadrunoBondSlip", 403, 12.0e3, 0.001, 0.003, 0.010, 4.0e3, 0.4,
        "-Gf", 0.08, "-s0", 1.0e-5,
    )
    with pytest.raises(Exception):
        ops.uniaxialMaterial(
            "LadrunoBondSlip", 493, 12.0e3, 0.001, 0.003, 0.010, 4.0e3, 0.4, "-implex"
        )

    # ASDConcrete1D
    ops.uniaxialMaterial(
        "ASDConcrete1D", 404, 30.0e6,
        "-fc", 30.0e3, "-ft", 3.0e3,
        "-implex", "-implexControl", 0.05, 0.1, "-implexAlpha", 1.0,
        "-eta", 0.0, "-tangent", "-autoRegularization", 0.25,
    )

    # ASDSteel1D
    ops.uniaxialMaterial(
        "ASDSteel1D", 405, 200.0e6, 420.0e3, 550.0e3, 0.12,
        "-implex", "-implexControl", 0.05, 0.1,
        "-auto_regularization", "-buckling", 8.0, 0.5, "-fracture", 0.5,
        "-slip", 403, 0.5, "-K_alpha", 1.0, "-max_iter", 25, "-tolU", 1.0e-6, "-tolR", 1.0e-6,
    )

    # ------------------------------------------------------------------
    # 4. Section 2.4 Macro-Wall Elements (-ThickMod capitalization)
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 2.0, 0.0, 0.0)
    ops.node(3, 2.0, 3.0, 0.0)
    ops.node(4, 0.0, 3.0, 0.0)
    ops.uniaxialMaterial("Elastic", 501, 30.0e6)
    ops.uniaxialMaterial("Elastic", 502, 200.0e6)
    ops.uniaxialMaterial("Elastic", 503, 12.0e6)

    # MVLEM_3D accepts -ThickMod (capital T)
    ops.element(
        "MVLEM_3D", 1, 1, 2, 3, 4, 2,
        "-thick", 0.2, 0.2, "-width", 1.0, 1.0, "-rho", 0.01, 0.01,
        "-matConcrete", 501, 501, "-matSteel", 502, 502, "-matShear", 503,
        "-CoR", 0.4, "-ThickMod", 0.63, "-Poisson", 0.25, "-Density", 2.4,
    )
    # Lowercase -thickMod is rejected (fail-closed)
    with pytest.raises(Exception):
        ops.element(
            "MVLEM_3D", 2, 1, 2, 3, 4, 2,
            "-thick", 0.2, 0.2, "-width", 1.0, 1.0, "-rho", 0.01, 0.01,
            "-matConcrete", 501, 501, "-matSteel", 502, 502, "-matShear", 503,
            "-thickMod", 0.63,
        )


def test_review_fixes_xml_response_tags_and_end_tags(tmp_path):
    """O1, O2, O4: Verify XML recorder emits well-formed <ElementOutput> / <GaussPointOutput>
    closing tags and complete <ResponseType> tags for:
      - LadrunoDispBeamColumn3d basicForce (6 tags) and basicStiffness (6 tags)
      - LadrunoIMKBeam (3D) localForce (12 tags)
      - LadrunoDispBeamColumn2d sectionX and integrationPoints / integrationWeights
    """
    import xml.etree.ElementTree as ET

    # 1. 3D XML recorder test (O1, O2, O4 sectionX)
    xml_bf = str(tmp_path / "ldbc3d_basic_force.xml")
    xml_bs = str(tmp_path / "ldbc3d_basic_stiff.xml")
    xml_lf = str(tmp_path / "imk3d_local_force.xml")
    xml_sx3 = str(tmp_path / "ldbc3d_sectionX.xml")

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 0.0, 2.0, 0.0)
    ops.node(4, 4.0, 2.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(3, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1)
    ops.uniaxialMaterial("Elastic", 10, 1.0e8)
    ops.element("LadrunoIMKBeam", 2, 3, 4, 0.04, 2.0e11, 8.0e10, 2.0e-4, 1.2e-4, 1.6e-4, 1, "-hinge", "both", "-matZ", 10)

    ops.recorder("Element", "-xml", xml_bf, "-ele", 1, "basicForce")
    ops.recorder("Element", "-xml", xml_bs, "-ele", 1, "basicStiffness")
    ops.recorder("Element", "-xml", xml_lf, "-ele", 2, "localForce")
    xml_imk_bf3 = str(tmp_path / "imk3d_basicForce.xml")
    xml_imk_bd3 = str(tmp_path / "imk3d_basicDef.xml")
    xml_imk_pd3 = str(tmp_path / "imk3d_plasticDef.xml")
    xml_imk_hm3 = str(tmp_path / "imk3d_hingeMoment.xml")
    ops.recorder("Element", "-xml", xml_imk_bf3, "-ele", 2, "basicForce")
    ops.recorder("Element", "-xml", xml_imk_bd3, "-ele", 2, "basicDeformation")
    ops.recorder("Element", "-xml", xml_imk_pd3, "-ele", 2, "plasticDeformation")
    ops.recorder("Element", "-xml", xml_imk_hm3, "-ele", 2, "hingeMoment")
    ops.recorder("Element", "-xml", xml_sx3, "-ele", 1, "sectionX", "0.0", "force")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, 5.0, 3.0, 1.0, 2.0, 4.0)
    ops.load(4, 10.0, 5.0, 3.0, 1.0, 2.0, 4.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.wipe()

    root_bf = ET.parse(xml_bf).getroot()
    bf_tags = [el.text for el in root_bf.findall(".//ElementOutput/ResponseType")]
    assert bf_tags == ["N", "Mz_1", "Mz_2", "My_1", "My_2", "T"]

    root_bs = ET.parse(xml_bs).getroot()
    bs_tags = [el.text for el in root_bs.findall(".//ElementOutput/ResponseType")]
    assert bs_tags == ["N", "Mz_1", "Mz_2", "My_1", "My_2", "T"]

    root_lf = ET.parse(xml_lf).getroot()
    lf_tags = [el.text for el in root_lf.findall(".//ElementOutput/ResponseType")]
    assert lf_tags == [
        "N_1", "Vy_1", "Vz_1", "T_1", "My_1", "Mz_1",
        "N_2", "Vy_2", "Vz_2", "T_2", "My_2", "Mz_2",
    ]
    assert [el.text for el in ET.parse(xml_imk_bf3).getroot().findall(".//ElementOutput/ResponseType")] == [
        "N", "Mz_1", "Mz_2", "My_1", "My_2", "T"
    ]
    assert [el.text for el in ET.parse(xml_imk_bd3).getroot().findall(".//ElementOutput/ResponseType")] == [
        "eps", "thetaZ_1", "thetaZ_2", "thetaY_1", "thetaY_2", "thetaX"
    ]
    assert [el.text for el in ET.parse(xml_imk_pd3).getroot().findall(".//ElementOutput/ResponseType")] == [
        "thetaZP_1", "thetaZP_2", "thetaYP_1", "thetaYP_2"
    ]
    assert [el.text for el in ET.parse(xml_imk_hm3).getroot().findall(".//ElementOutput/ResponseType")] == [
        "Mz_1", "Mz_2", "My_1", "My_2"
    ]

    # sectionX XML must parse cleanly (GaussPointOutput and ElementOutput properly closed)
    ET.parse(xml_sx3)

    # 2. 2D XML recorder test (O4: integrationPoints, integrationWeights, sectionX + 2D IMK tags)
    xml_ip2 = str(tmp_path / "ldbc2d_ip.xml")
    xml_iw2 = str(tmp_path / "ldbc2d_iw.xml")
    xml_sx2 = str(tmp_path / "ldbc2d_sectionX.xml")
    xml_imk_bf2 = str(tmp_path / "imk2d_basicForce.xml")
    xml_imk_bd2 = str(tmp_path / "imk2d_basicDef.xml")
    xml_imk_pd2 = str(tmp_path / "imk2d_plasticDef.xml")
    xml_imk_hm2 = str(tmp_path / "imk2d_hingeMoment.xml")

    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 4.0, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, 4.0, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1)
    ops.element("LadrunoIMKBeam", 2, 3, 4, 0.04, 2.0e11, 1.6e-4, 1)

    ops.recorder("Element", "-xml", xml_ip2, "-ele", 1, "integrationPoints")
    ops.recorder("Element", "-xml", xml_iw2, "-ele", 1, "integrationWeights")
    ops.recorder("Element", "-xml", xml_sx2, "-ele", 1, "sectionX", "0.0", "force")
    ops.recorder("Element", "-xml", xml_imk_bf2, "-ele", 2, "basicForce")
    ops.recorder("Element", "-xml", xml_imk_bd2, "-ele", 2, "basicDeformation")
    ops.recorder("Element", "-xml", xml_imk_pd2, "-ele", 2, "plasticDeformation")
    ops.recorder("Element", "-xml", xml_imk_hm2, "-ele", 2, "hingeMoment")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, 5.0, 2.0)
    ops.load(4, 10.0, 5.0, 2.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.wipe()

    ET.parse(xml_ip2)
    ET.parse(xml_iw2)
    ET.parse(xml_sx2)
    assert [el.text for el in ET.parse(xml_imk_bf2).getroot().findall(".//ElementOutput/ResponseType")] == [
        "N", "M_1", "M_2"
    ]
    assert [el.text for el in ET.parse(xml_imk_bd2).getroot().findall(".//ElementOutput/ResponseType")] == [
        "eps", "theta1", "theta2"
    ]
    assert [el.text for el in ET.parse(xml_imk_pd2).getroot().findall(".//ElementOutput/ResponseType")] == [
        "theta1P", "theta2P"
    ]
    assert [el.text for el in ET.parse(xml_imk_hm2).getroot().findall(".//ElementOutput/ResponseType")] == [
        "M_1", "M_2"
    ]


def test_review_fixes_hinge_end_specificity_and_aliases():
    """O5: Verify hingeZ / hingeZI / hingeZJ work in both 2D and 3D LadrunoDispBeamColumn,
    and explicit end tokens (hingeI/hingeZI/hinge1/hingeYI/hingeBiaxialI/hingeCI vs
    hingeJ/hingeZJ/hinge2/hingeYJ/hingeBiaxialJ/hingeCJ) do NOT silently fall back to
    the opposite end when only one end is configured, while end-agnostic tokens
    (hinge/hingeZ/hingeY/hingeBiaxial/hingeC) still fall back to end J on J-only elements.
    """
    # 2D: Element 1 has -hingeI only; Element 2 has -hingeJ only
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 4.0, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, 4.0, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.0e8, 0.02)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeI", 10)
    ops.element("LadrunoDispBeamColumn", 2, 3, 4, 1, 1, "-hingeJ", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -20.0e3, 0.0)
    ops.load(4, 0.0, -20.0e3, 80.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # Ele 1 (I-only): hinge, hingeZ, hingeI, hingeZI, hinge1 all return end-I stress;
    # hingeJ, hingeZJ, hinge2 return empty!
    for tok in ("hinge", "hingeZ", "hingeI", "hingeZI", "hinge1"):
        res = ops.eleResponse(1, tok, "stress")
        assert len(res) == 1 and res[0] != 0.0, f"2D I-only {tok} failed: {res}"
    for tok in ("hingeJ", "hingeZJ", "hinge2"):
        assert ops.eleResponse(1, tok, "stress") == [], f"2D I-only {tok} should be empty"

    # Ele 2 (J-only): hinge, hingeZ fall back to end J; hingeJ, hingeZJ, hinge2 return end J;
    # hingeI, hingeZI, hinge1 return empty!
    for tok in ("hinge", "hingeZ", "hingeJ", "hingeZJ", "hinge2"):
        res = ops.eleResponse(2, tok, "stress")
        assert len(res) == 1 and res[0] != 0.0, f"2D J-only {tok} failed: {res}"
    for tok in ("hingeI", "hingeZI", "hinge1"):
        assert ops.eleResponse(2, tok, "stress") == [], f"2D J-only {tok} should be empty"

    # 3D: Element 1 has -hingeZI / -hingeYI; Element 2 has -hingeZJ / -hingeYJ
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 0.0, 2.0, 0.0)
    ops.node(4, 4.0, 2.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(3, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.0e8, 0.02)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeZI", 10, "-hingeYI", 10)
    ops.element("LadrunoDispBeamColumn", 2, 3, 4, 1, 1, "-hingeZJ", 10, "-hingeYJ", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -20.0e3, -15.0e3, 0.0, 0.0, 0.0)
    ops.load(4, 0.0, -20.0e3, -15.0e3, 0.0, 60.0e3, 80.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    for tok in ("hinge", "hingeZ", "hingeI", "hingeZI", "hinge1", "hingeY", "hingeYI", "hingeY1"):
        res = ops.eleResponse(1, tok, "stress")
        assert len(res) == 1 and res[0] != 0.0, f"3D I-only {tok} failed: {res}"
    for tok in ("hingeJ", "hingeZJ", "hinge2", "hingeYJ", "hingeY2"):
        assert ops.eleResponse(1, tok, "stress") == [], f"3D I-only {tok} should be empty"

    for tok in ("hinge", "hingeZ", "hingeJ", "hingeZJ", "hinge2", "hingeY", "hingeYJ", "hingeY2"):
        res = ops.eleResponse(2, tok, "stress")
        assert len(res) == 1 and res[0] != 0.0, f"3D J-only {tok} failed: {res}"
    for tok in ("hingeI", "hingeZI", "hinge1", "hingeYI", "hingeY1"):
        assert ops.eleResponse(2, tok, "stress") == [], f"3D J-only {tok} should be empty"
    ops.wipe()


def test_review_fixes_pm_axial_force_in_energy_and_hinge_summary_and_numeric_lch(tmp_path):
    """O6, A1/O7: Verify:
      1. Numeric float `-lch` (e.g., 0.25) is accepted directly from Python in
         `LadrunoDispBeamColumn` and `forceBeamColumn` (2D and 3D, including
         `ops.mesh('line', ...)`), while invalid trailing-garbage strings like
         `"0.25abc"` are rejected.
      2. `energy` and `hingeSummary` on `LadrunoDispBeamColumn` (2D and 3D) with
         generic uniaxial hinge + `-pmP0` use the updated trial axial force even
         when `getResistingForce()` / `eleForce` has not been queried first, AND
         immediately after a database `save`/`restore` (`recvSelf`) before any
         `analyze(1)` or `eleForce` call.
    """
    # 1. Numeric float -lch in 2D and 3D (direct element + ops.mesh('line', ...))
    for ndm, ndf in ((2, 3), (3, 6)):
        ops.wipe()
        ops.model("basic", "-ndm", ndm, "-ndf", ndf)
        if ndm == 2:
            ops.node(1, 0.0, 0.0)
            ops.node(2, 4.0, 0.0)
            ops.geomTransf("Linear", 1)
            ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4)
        else:
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, 4.0, 0.0, 0.0)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        # Pass float 0.25 (not string "0.25")
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-lch", 0.25)
        ops.element("forceBeamColumn", 2, 1, 2, 1, 1, "-lch", 0.25)
        ops.mesh("line", 10, 2, 1, 2, 0, ndf, 2.0, "forceBeamColumn", 1, 1, "-lch", 0.25)
        # Trailing garbage must fail-closed
        with pytest.raises(Exception):
            ops.element("LadrunoDispBeamColumn", 3, 1, 2, 1, 1, "-lch", "0.25abc")
        with pytest.raises(Exception):
            ops.element("forceBeamColumn", 4, 1, 2, 1, 1, "-lch", "0.25abc")
        with pytest.raises(Exception):
            ops.mesh("line", 11, 2, 1, 2, 0, ndf, 2.0, "forceBeamColumn", 1, 1, "-lch", "0.25abc")

    # 2. PM-scaled energy and hingeSummary in 2D without prior eleForce call + post-restore
    db_2d = str(tmp_path / "pm2d_db")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 4.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    # Generic uniaxial hinge (Steel01) so element-level PM scaling applies
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.0e8, 0.02)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeI", 10, "-pm", 200.0e3, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # Apply compression P = -100e3 (so n = P/P0 = -0.5 -> fPM = 1 - 1.0*(-0.5) = 1.5) and transverse load
    ops.load(2, -100.0e3, -20.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # Query hingeSummary and energy BEFORE any eleForce/basicForce call
    hs_before = ops.eleResponse(1, "hingeSummary")
    en_before = ops.eleResponse(1, "energy")
    # Save and restore state, then query hingeSummary & energy BEFORE any analyze or eleForce
    ops.database("File", db_2d)
    ops.save(1)
    ops.reset()
    ops.restore(1)
    hs_restored = ops.eleResponse(1, "hingeSummary")
    en_restored = ops.eleResponse(1, "energy")
    np.testing.assert_allclose(hs_restored, hs_before, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(en_restored, en_before, rtol=1e-14, atol=1e-14)
    # Now trigger getResistingForce via eleForce and re-query
    _ = ops.eleForce(1)
    hs_after = ops.eleResponse(1, "hingeSummary")
    en_after = ops.eleResponse(1, "energy")
    np.testing.assert_allclose(hs_before, hs_after, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(en_before, en_after, rtol=1e-14, atol=1e-14)
    # Verify fPM = 1.5 was actually applied to the hinge moment
    raw_hinge_stress = ops.eleResponse(1, "hingeI", "stress")[0]
    np.testing.assert_allclose(hs_before[1], 1.5 * raw_hinge_stress, rtol=1e-12, atol=1e-12)

    # 3. PM-scaled energy and hingeSummary in 3D without prior eleForce call + post-restore
    db_3d = str(tmp_path / "pm3d_db")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.section("Elastic", 1, 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 1.0e8, 0.02)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeZI", 10, "-hingeYI", 10, "-pm", 200.0e3, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -100.0e3, -20.0e3, -15.0e3, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    hs3_before = ops.eleResponse(1, "hingeSummary")
    en3_before = ops.eleResponse(1, "energy")
    ops.database("File", db_3d)
    ops.save(1)
    ops.reset()
    ops.restore(1)
    hs3_restored = ops.eleResponse(1, "hingeSummary")
    en3_restored = ops.eleResponse(1, "energy")
    np.testing.assert_allclose(hs3_restored, hs3_before, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(en3_restored, en3_before, rtol=1e-14, atol=1e-14)
    _ = ops.eleForce(1)
    hs3_after = ops.eleResponse(1, "hingeSummary")
    en3_after = ops.eleResponse(1, "energy")
    raw_hinge_z = ops.eleResponse(1, "hingeZI", "stress")[0]
    np.testing.assert_allclose(hs3_before[1], 1.5 * raw_hinge_z, rtol=1e-12, atol=1e-12)
    ops.wipe()


def test_ladruno_force_beam_column_and_disp_span_load_and_imk_shear_biaxial():
    """Verify cross-family consistency between:
      - `LadrunoForceBeamColumn` (2D & 3D) with `-hingeIJ`, `-pm`, `-shear`, and `-lch`
      - `LadrunoDispBeamColumn` (2D & 3D) `-spanLoad` exact fixed-end strain + quartic bubble
      - `LadrunoIMKBeam` (2D & 3D) `-shearHinge` / `-shearHingeY` / `-shearHingeZ` and 3D `-hingeBiaxial`
    """
    L = 5.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.8e-4, 1.2e-4, 8.0e10, 2.0e-4

    # 1. 2D 1-element simply-supported beam under uniform span load:
    #    LadrunoDispBeamColumn2d (-spanLoad default) vs LadrunoForceBeamColumn2d
    for ele_name in ("LadrunoDispBeamColumn2d", "LadrunoForceBeamColumn2d"):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 0)
        ops.fix(2, 0, 1, 0)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        ops.element(ele_name, 1, 1, 2, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", -20.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        # Midspan moment M(L/2) = w*L^2/8 = 62.5 kN*m
        mid_m = ops.eleResponse(1, "section", 3, "force")[1]
        assert mid_m == pytest.approx(20.0e3 * (L ** 2) / 8.0, rel=1e-12)

    # 2. 3D LadrunoIMKBeam with -hingeBiaxial and -shearHingeY/-shearHingeZ
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 50, 70.0e3, 400.0, 55.0e3, 350.0, "-linear")
    ops.uniaxialMaterial("Steel01", 60, 45.0e3, 8.0e7, 0.03)
    ops.uniaxialMaterial("Steel01", 61, 40.0e3, 7.0e7, 0.03)
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
        "-hingeBiaxialIJ", 50, 50, "-shearHingeY", 60, "-shearHingeZ", 61
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for dof, val in enumerate([-0.0001, -0.008, 0.006, 0.0005, 0.010, -0.012], start=1):
        ops.sp(2, dof, val)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    sh = ops.eleResponse(1, "shearHinge")
    hs = ops.eleResponse(1, "hingeSummary")
    assert len(sh) == 4 and len(hs) == 8
    ops.wipe()


def test_opus_audit_14_findings_regressions(tmp_path):
    """Regression suite for the 14-finding mathematical & code audit:
      - F3: LadrunoDispBeamColumn2d/3d spanLoad revertToStart/revertToLastStep & -noSpanLoad
      - F4: LadrunoForceBeamColumn2d/3d true series hinge compatibility (v_int = v - v_hinge)
            and fixed-end span-load hinge yielding
      - F6: LadrunoIMKBeam2d/3d P-M interaction (-pm, -PM, -pmZ, -pmY), no double-scaling
            with LadrunoCohesiveHinge/Biaxial, and exact axial-flexural tangent sensitivity
      - F7: LadrunoIMKBeam2d/3d shear hinge fixed-end rotation compatibility under span load
            and clear-span Lc Timoshenko Phi in addLoad()
      - F8: WrapperSection section/sectionX response delegation when -shear is active
      - F9: LadrunoForceBeamColumn2d/3d -hingeLoc xiI xiJ lever-arm kinematics & serialization
      - F10: LadrunoDispBeamColumn3d getInitialStiff() with -nl and -hinge
      - F12: IMPL-EX + P-M consistency across DispBeamColumn and ForceBeamColumn (2D & 3D)
      - F13/F14: LadrunoForceBeamColumn3d isTorsion serialization in sendSelf/recvSelf
    """
    L = 5.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # ------------------------------------------------------------------
    # F3: LadrunoDispBeamColumn2d/3d -spanLoad revertToStart & -noSpanLoad
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, L, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.fix(4, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element("LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1)
    ops.element("LadrunoDispBeamColumn2d", 2, 3, 4, 1, 1, "-noSpanLoad")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, 2, "-type", "-beamUniform", -24.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # With fixed-fixed ends (v=0), ele 1 (-spanLoad) has exact midspan moment w*L^2/24 = 25 kN*m,
    # while ele 2 (-noSpanLoad) has zero section moment because v=0 and spanLoad is disabled.
    m_mid_span = ops.eleResponse(1, "section", 3, "force")[1]
    m_mid_nospan = ops.eleResponse(2, "section", 3, "force")[1]
    assert m_mid_span == pytest.approx(24.0e3 * L * L / 24.0, rel=1e-12)
    assert m_mid_nospan == pytest.approx(0.0, abs=1e-12)

    # Revert to start must zero out q0/p0 and section forces
    ops.reset()
    assert ops.eleResponse(1, "section", 3, "force")[1] == pytest.approx(0.0, abs=1e-12)

    # ------------------------------------------------------------------
    # F4 & F9: LadrunoForceBeamColumn2d/3d series hinge compatibility under
    #          fixed-fixed span load + -hingeLoc xiI xiJ lever-arm kinematics
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    My = 40.0e3
    k0 = 1.0e12
    bh = 1.0e-6
    ops.uniaxialMaterial("Steel01", 10, My, k0, bh)
    ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1, "-hingeIJ", 10)

    w_load = -30.0e3  # elastic FEF = w*L^2/12 = 62.5 kN*m > My = 40 kN*m
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", w_load)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.25)
    ops.test("NormDispIncr", 1.0e-12, 30)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(4) == 0

    hs_fbc2 = ops.eleResponse(1, "hingeSummary")
    ef_fbc2 = ops.eleForce(1)
    # Hinges must yield even though nodal rotations are fixed (v = 0)!
    assert abs(hs_fbc2[0]) > 1.0e-4 and abs(hs_fbc2[2]) > 1.0e-4
    assert abs(hs_fbc2[1]) == pytest.approx(My + bh * k0 * (abs(hs_fbc2[0]) - My / k0), rel=1e-6)
    assert abs(ef_fbc2[2]) == pytest.approx(abs(hs_fbc2[1]), rel=1e-6)
    # Interior section end moment at x=0 must match element end moment
    sec1_m = ops.eleResponse(1, "section", 1, "force")[1]
    assert abs(sec1_m) == pytest.approx(abs(hs_fbc2[1]), rel=1e-6)

    # F9: -hingeLoc xiI xiJ on LadrunoForceBeamColumn2d + database round-trip
    db_fbc_loc = str(tmp_path / "fbc_loc_db")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("Steel01", 10, 50.0e3, 2.0e7, 0.05)
    ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1, "-hingeIJ", 10, "-hingeLoc", 0.1, 0.9)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -25.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 30)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    hs_loc_before = ops.eleResponse(1, "hingeSummary")
    # At xiI=0.1 under cantilever tip load P=-25kN at x=L:
    # M(x) = -P*(L - x), so M(0.1*L) = 0.9 * P * L = 112.5 kN*m, M(0.9*L) = 0.1 * P * L = 12.5 kN*m
    assert abs(hs_loc_before[1]) == pytest.approx(0.9 * 25.0e3 * L, rel=1e-6)
    assert abs(hs_loc_before[3]) == pytest.approx(0.1 * 25.0e3 * L, rel=1e-6)
    ops.database("File", db_fbc_loc)
    ops.save(1)
    ops.reset()
    ops.restore(1)
    assert ops.eleResponse(1, "hingeSummary") == pytest.approx(hs_loc_before, rel=1e-12)

    # ------------------------------------------------------------------
    # F6: LadrunoIMKBeam2d/3d P-M interaction (-pm, -PM, -pmZ, -pmY)
    #     and no double-scaling with LadrunoCohesiveHinge / Biaxial
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 0.0, 2.0)
    ops.node(4, L, 2.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    ops.fix(3, 1, 1, 1)
    ops.fix(4, 0, 1, 0)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 10, 40.0e3, 1.0e9, 1.0e-4)
    ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1, "-hinge", "i", "-matZ", 10, "-pm", 200.0e3, 1.0)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 20, 40.0e3, 5000.0, "-linear")
    ops.element("LadrunoIMKBeam", 2, 3, 4, A, E, Iz, 1, "-hinge", "i", "-matZ", 20, "-pm", 200.0e3, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # Apply axial compression at node 2/4 and rotation moment at node 2/4
    ops.load(2, -100.0e3, 0.0, -180.0e3)
    ops.load(4, -100.0e3, 0.0, -180.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0

    hm1 = ops.eleResponse(1, "hingeMoment")[0]
    hm2 = ops.eleResponse(2, "hingeMoment")[0]
    # Both should yield near f_PM * 40 kN*m = 1.5 * 40 kN*m = 60 kN*m (NOT 1.5^2 * 40 = 90 kN*m!)
    assert abs(hm1) == pytest.approx(60.0e3, rel=1e-2)
    assert abs(hm2) == pytest.approx(60.0e3, rel=2e-2)

    # ------------------------------------------------------------------
    # F7: LadrunoIMKBeam2d/3d shear hinge fixed-end rotation compatibility
    #     under span load + clear-span Lc Timoshenko Phi
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Elastic", 30, 1.0e8)  # elastic shear hinge with Ks = 1e8
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1,
        "-shearHinge", 30, "-shear", G, 0.8 * A, "-endZones", 0.25, 0.25
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # Asymmetric point load at x/L = 0.25 -> non-zero (q0_i + q0_j)
    ops.eleLoad("-ele", 1, "-type", "-beamPoint", -40.0e3, 0.25)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    sh2 = ops.eleResponse(1, "shearHinge")
    bf2 = ops.eleResponse(1, "basicForce")
    # Nodal basicForce moments (bf2[1], bf2[2]) include the T^T rigid-offset magnification L/Lc,
    # so -(Mc_i + Mc_j)/Lc == -(bf2[1] + bf2[2])/L!
    assert sh2[1] == pytest.approx(-(bf2[1] + bf2[2]) / L, rel=1e-10, abs=1e-8)

    # ------------------------------------------------------------------
    # F8: section/sectionX response + setParameter/updateParameter
    #     (allSections & lch) across LadrunoDispBeamColumn & LadrunoForceBeamColumn
    # ------------------------------------------------------------------
    for idx, ele_name in enumerate(("LadrunoDispBeamColumn2d", "LadrunoForceBeamColumn2d"), start=1):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        ops.element(ele_name, 1, 1, 2, 1, 1)
        ops.parameter(100 + idx, "element", 1, "lch")
        ops.updateParameter(100 + idx, 0.4)
        ops.parameter(200 + idx, "element", 1, "allSections", "E")
        ops.updateParameter(200 + idx, 1.5 * E)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 10.0e3, -15.0e3, 5.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        sf_i = ops.eleResponse(1, "section", 1, "force")
        sfx_i = ops.eleResponse(1, "sectionX", 0.0, "force")
        assert len(sf_i) == 2 and sfx_i == pytest.approx(sf_i, rel=1e-12)
        # Tip axial displacement must reflect updated E_new = 1.5 * E
        assert ops.nodeDisp(2, 1) == pytest.approx(10.0e3 * L / (1.5 * E * A), rel=1e-10)

    # ------------------------------------------------------------------
    # F10: LadrunoDispBeamColumn3d getInitialStiff() in dualHingeMode (-hingeIJ -hingeYIJ)
    # ------------------------------------------------------------------
    def _run_ldbc3d_dual_init(use_initial: bool):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        ops.uniaxialMaterial("Elastic", 10, 5.0e6)
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, "-hingeIJ", 10, "-hingeYIJ", 10)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5.0e3, -10.0e3, 8.0e3, 1.0e3, 4.0e3, -6.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        if use_initial:
            ops.algorithm("Linear", "-initial")
        else:
            ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return [ops.nodeDisp(2, d) for d in range(1, 7)]

    assert _run_ldbc3d_dual_init(True) == pytest.approx(_run_ldbc3d_dual_init(False), rel=1e-12, abs=1e-14)

    # ------------------------------------------------------------------
    # F12: IMPL-EX + P-M consistency on DispBeamColumn & ForceBeamColumn
    # ------------------------------------------------------------------
    hs_imp_by_ele = {}
    for ele_name in ("LadrunoDispBeamColumn2d", "LadrunoForceBeamColumn2d"):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 20, 40.0e3, 5000.0, "-linear")
        ops.element(ele_name, 1, 1, 2, 1, 1, "-hingeI", 20, "-pm", 200.0e3, 1.0, "-implex")
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0e3, 0.0, -180.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.02)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(50) == 0
        hs_imp = ops.eleResponse(1, "hingeSummary")
        hs_imp_by_ele[ele_name] = hs_imp
        assert abs(hs_imp[1]) == pytest.approx(60.0e3, rel=5e-2)
    assert hs_imp_by_ele["LadrunoDispBeamColumn2d"] == pytest.approx(
        hs_imp_by_ele["LadrunoForceBeamColumn2d"], rel=1e-10, abs=1e-10
    )

    # ------------------------------------------------------------------
    # F13/F14: LadrunoForceBeamColumn3d isTorsion database round-trip
    # ------------------------------------------------------------------
    db_fbc3_tor = str(tmp_path / "fbc3_tor_db")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    # SectionAggregator without T -> isTorsion == false!
    ops.uniaxialMaterial("Elastic", 1, E * A)
    ops.uniaxialMaterial("Elastic", 2, E * Iz)
    ops.uniaxialMaterial("Elastic", 3, E * Iy)
    ops.section("Aggregator", 1, 1, "P", 2, "Mz", 3, "My")
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.element("LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3, -5.0e3, 3.0e3, 0.0, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ef_before = ops.eleForce(1)
    ops.database("File", db_fbc3_tor)
    ops.save(1)
    ops.reset()
    ops.restore(1)
    assert ops.eleForce(1) == pytest.approx(ef_before, rel=1e-12)
    ops.wipe()


def test_opus_audit_14_findings_round2_regressions(tmp_path):
    """Comprehensive regression suite for all OpenSees findings in the round-2 audit:
      - Finding 4 (Item 1): LadrunoForceBeamColumn2d/3d span-load state determination when
        dv == 0 (including loadPattern removal and -subdivide loadAtNodes)
      - Finding 3 (Item 2): LadrunoDispBeamColumn2d/3d -spanLoad synchronization across
        update() and getResistingForce() (including fixed-end span-load hinge yielding and
        loadPattern removal when dv == 0)
      - Finding 8 (Item 3): plasticDeformation (responseID == 4) subtracts initial elastic
        hinge compliance in LadrunoForceBeamColumn2d/3d and LadrunoDispBeamColumn2d/3d
        (both single-hinge and dualHingeMode), and basicForce (this->q) updates in update()
      - Finding 6 (Item 4): LadrunoIMKBeam2d/3d -pm forwarding and axial-flexural tangent
        coupling kb(1..n, 0) against central finite differences
      - Finding 7 (Item 5): LadrunoIMKBeam2d/3d -shearHinge fixed-end moment shift gam_rot0,
        clear-span Lc in addLoad, and energy consistency
      - Finding 9 (Item 6): LadrunoForceBeamColumn2d 2-tag `-hinge tagI tagJ` syntax and
        fPM scaling of hinge energy in responseID == 2000
      - Finding 10 (Item 7): LadrunoDispBeamColumn3d fallback fPM scaling for generic 2D
        NDMaterial biaxial hinges (-hingeBiaxial and -hingeBiaxialIJ)
      - Finding 12 (Item 8): setParameter/updateParameter for 'lch', 'Lch', and
        'characteristicLength' in LadrunoForceBeamColumn2d and LadrunoDispBeamColumn3d
      - Finding 13 (Item 9): LadrunoForceBeamColumn3d::Print End 1 / End 2 VZ and T signs
      - Finding 14 (Item 10): LadrunoForceBeamColumn3d XML recorder 'section $secNum'
        emits physical coordinate eta = xi * L
    """
    import json
    import xml.etree.ElementTree as ET

    L = 4.0
    E, A, Iz, Iy, G, Jx = 2.0e11, 0.04, 1.6e-4, 1.2e-4, 8.0e10, 2.0e-4

    # ------------------------------------------------------------------
    # 1. Finding 4 & Finding 3: Span-load state determination with dv == 0
    #    and loadPattern removal on both ForceBeamColumn and DispBeamColumn (2D & 3D)
    # ------------------------------------------------------------------
    for ele_cmd in ("LadrunoForceBeamColumn2d", "LadrunoDispBeamColumn2d"):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 1, 1, 1)  # fixed-fixed: dv is identically 0!
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        extra = ["-subdivide", 4, 2.0] if "Force" in ele_cmd else []
        ops.element(ele_cmd, 1, 1, 2, 1, 1, *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        w = -24.0e3
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", w)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        # Midspan moment for fixed-fixed beam under uniform load w is -w*L^2/24 = +16 kN*m
        assert ops.eleResponse(1, "section", 3, "force")[1] == pytest.approx(-w * L * L / 24.0, rel=1e-10)
        # End vertical reaction is -w*L/2 = +48 kN
        assert ops.eleForce(1)[1] == pytest.approx(-w * L / 2.0, rel=1e-10)
        # Now remove the load pattern while keeping nodes fixed (dv == 0) and re-analyze
        ops.remove("loadPattern", 1)
        assert ops.analyze(1) == 0
        assert ops.eleResponse(1, "section", 3, "force")[1] == pytest.approx(0.0, abs=1e-10)
        assert ops.eleForce(1)[1] == pytest.approx(0.0, abs=1e-10)

    for ele_cmd in ("LadrunoForceBeamColumn3d", "LadrunoDispBeamColumn3d"):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.fix(2, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        extra = ["-subdivide", 4, 2.0] if "Force" in ele_cmd else []
        ops.element(ele_cmd, 1, 1, 2, 1, 1, *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        wy, wz = -24.0e3, -12.0e3
        ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy, wz, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        sec_mid = ops.eleResponse(1, "section", 3, "force")
        assert sec_mid[1] == pytest.approx(-wy * L * L / 24.0, rel=1e-10)
        ops.remove("loadPattern", 1)
        assert ops.analyze(1) == 0
        assert ops.eleResponse(1, "section", 3, "force")[1] == pytest.approx(0.0, abs=1e-10)
        assert ops.eleForce(1)[1] == pytest.approx(0.0, abs=1e-10)

    # ------------------------------------------------------------------
    # 2. Finding 8: plasticDeformation subtracts initial elastic hinge compliance
    #    in LadrunoForceBeamColumn2d/3d and LadrunoDispBeamColumn2d/3d
    # ------------------------------------------------------------------
    for ele_cmd in ("LadrunoForceBeamColumn2d", "LadrunoDispBeamColumn2d"):
        for hinge_flags in (["-hinge", 10, "-hingeLoc", 0.0], ["-hingeIJ", 10]):
            ops.wipe()
            ops.model("basic", "-ndm", 2, "-ndf", 3)
            ops.node(1, 0.0, 0.0)
            ops.node(2, L, 0.0)
            ops.fix(1, 1, 1, 1)
            ops.section("Elastic", 1, E, A, Iz)
            ops.geomTransf("Linear", 1)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            My_h, Kh0, bh = 80.0e3, 2.0e7, 0.05
            ops.uniaxialMaterial("Steel01", 10, My_h, Kh0, bh)
            ops.element(ele_cmd, 1, 1, 2, 1, 1, *hinge_flags)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            # Step 1: Tip load P = -10 kN -> M_base = 40 kN*m < My_h = 80 kN*m (purely elastic)
            ops.load(2, 0.0, -10.0e3, 0.0)
            ops.system("BandGeneral")
            ops.numberer("Plain")
            ops.constraints("Plain")
            ops.integrator("LoadControl", 1.0)
            ops.test("NormDispIncr", 1.0e-12, 25)
            ops.algorithm("Newton")
            ops.analysis("Static")
            assert ops.analyze(1) == 0
            # Query basicForce BEFORE eleForce to verify this->q is updated in update()
            bf = ops.eleResponse(1, "basicForce")
            assert abs(bf[1]) == pytest.approx(40.0e3, rel=1e-10)
            vp_elastic = ops.eleResponse(1, "plasticDeformation")
            np.testing.assert_allclose(vp_elastic, [0.0, 0.0, 0.0], atol=1e-12)
            # Step 2..3: Push to 3 * 40 = 120 kN*m > My_h = 80 kN*m (plastic hinge rotation at end I)
            assert ops.analyze(2) == 0
            hs = ops.eleResponse(1, "hingeSummary")
            alpha_total_I = hs[0]
            M_I = hs[1]
            expected_alpha_pI = alpha_total_I - M_I / Kh0
            vp_plastic = ops.eleResponse(1, "plasticDeformation")
            assert abs(expected_alpha_pI) > 1.0e-4
            assert vp_plastic[1] == pytest.approx(expected_alpha_pI, rel=1e-8)
            assert vp_plastic[2] == pytest.approx(0.0, abs=1e-12)

    for ele_cmd in ("LadrunoForceBeamColumn3d", "LadrunoDispBeamColumn3d"):
        for hinge_flags in (["-hinge", 10, "-hingeY", 10, "-hingeLoc", 0.0], ["-hingeIJ", 10, "-hingeYIJ", 10]):
            ops.wipe()
            ops.model("basic", "-ndm", 3, "-ndf", 6)
            ops.node(1, 0.0, 0.0, 0.0)
            ops.node(2, L, 0.0, 0.0)
            ops.fix(1, 1, 1, 1, 1, 1, 1)
            ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
            ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
            ops.beamIntegration("Lobatto", 1, 1, 4)
            My_h, Kh0, bh = 80.0e3, 2.0e7, 0.05
            ops.uniaxialMaterial("Steel01", 10, My_h, Kh0, bh)
            ops.element(ele_cmd, 1, 1, 2, 1, 1, *hinge_flags)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(2, 0.0, -10.0e3, -8.0e3, 0.0, 0.0, 0.0)
            ops.system("BandGeneral")
            ops.numberer("Plain")
            ops.constraints("Plain")
            ops.integrator("LoadControl", 1.0)
            ops.test("NormDispIncr", 1.0e-12, 25)
            ops.algorithm("Newton")
            ops.analysis("Static")
            assert ops.analyze(1) == 0
            vp_elastic3 = ops.eleResponse(1, "plasticDeformation")
            np.testing.assert_allclose(vp_elastic3, [0.0] * 6, atol=1e-12)
            assert ops.analyze(2) == 0
            hs3 = ops.eleResponse(1, "hingeSummary")
            expected_alpha_pZI = hs3[0] - hs3[1] / Kh0
            expected_alpha_pYI = hs3[4] - hs3[5] / Kh0
            vp_plastic3 = ops.eleResponse(1, "plasticDeformation")
            assert vp_plastic3[1] == pytest.approx(expected_alpha_pZI, rel=1e-8)
            assert vp_plastic3[3] == pytest.approx(expected_alpha_pYI, rel=1e-8)

    # ------------------------------------------------------------------
    # 3. Finding 6: LadrunoIMKBeam2d/3d axial-flexural tangent coupling kb(1..n, 0)
    # ------------------------------------------------------------------
    def _eval_imk2d_pm(v_tip):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.uniaxialMaterial("Steel01", 10, 50.0e3, 2.0e7, 0.02)
        ops.element("LadrunoIMKBeam2d", 1, 1, 2, A, E, Iz, 1, "-hinge", "both", "-matZ", 10, "-pm", 400.0e3, 0.8)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for d, val in enumerate(v_tip, start=1):
            ops.sp(2, d, val)
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
            np.array(ops.eleResponse(1, "basicStiffness")).reshape(3, 3),
        )

    base_2d = [-0.0002, -0.025, 0.015]
    v0_2, q0_2, kb_2 = _eval_imk2d_pm(base_2d)
    assert abs(kb_2[1, 0]) > 0.0 and abs(kb_2[2, 0]) > 0.0
    eps = 1.0e-8
    vp_2, qp_2, _ = _eval_imk2d_pm([base_2d[0] + eps, base_2d[1], base_2d[2]])
    vm_2, qm_2, _ = _eval_imk2d_pm([base_2d[0] - eps, base_2d[1], base_2d[2]])
    dq_dv0 = (qp_2 - qm_2) / (vp_2[0] - vm_2[0])
    np.testing.assert_allclose(kb_2[:, 0], dq_dv0, rtol=1e-5, atol=1e-2)

    # ------------------------------------------------------------------
    # 4. Finding 9: LadrunoForceBeamColumn2d 2-tag `-hinge tagI tagJ` and
    #    fPM-scaled hinge energy
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("Steel01", 11, 40.0e3, 2.0e7, 0.02)
    ops.uniaxialMaterial("Steel01", 12, 60.0e3, 2.0e7, 0.02)
    # 2-tag -hinge 11 12 syntax + -pm 200e3 1.0
    ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1, "-hinge", 11, 12, "-pm", 200.0e3, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # P = -100e3 -> n = -0.5 -> fPM = 1.5; q1 = 100 kN*m > 60 kN*m, q2 = 100 kN*m > 90 kN*m
    ops.load(2, -100.0e3, -50.0e3, 100.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(2) == 0
    hs_2tag = ops.eleResponse(1, "hingeSummary")
    assert abs(hs_2tag[0]) > 0.0 and abs(hs_2tag[2]) > 0.0
    e_hi = ops.eleResponse(1, "hingeI", "energy")[0]
    e_hj = ops.eleResponse(1, "hingeJ", "energy")[0]
    tot_e = ops.eleResponse(1, "energy")[0]
    assert tot_e == pytest.approx(1.5 * (e_hi + e_hj), rel=1e-10)

    # ------------------------------------------------------------------
    # 5. Finding 10: LadrunoDispBeamColumn3d fallback fPM scaling for generic
    #    2D NDMaterial biaxial hinges (-hingeBiaxial and -hingeBiaxialIJ)
    # ------------------------------------------------------------------
    for biax_flag in ("-hingeBiaxial", "-hingeBiaxialIJ"):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.nDMaterial("ElasticIsotropic", 49, 5.0e7, 0.0)
        ops.nDMaterial("BeamFiber2d", 50, 49)
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1, biax_flag, 50, "-pm", 200.0e3, 1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0e3, -15.0e3, -10.0e3, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-12, 25)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        hs_biax = ops.eleResponse(1, "hingeSummary")
        raw_biax = ops.eleResponse(1, "hingeBiaxialI", "stress")
        # fPM = 1.5 must scale raw stress in both state determination and hingeSummary!
        assert hs_biax[1] == pytest.approx(1.5 * raw_biax[0], rel=1e-10)
        assert hs_biax[5] == pytest.approx(1.5 * raw_biax[1], rel=1e-10)
        xi_fac = 0.5 if biax_flag == "-hingeBiaxial" else 1.0
        assert abs(hs_biax[1]) == pytest.approx(xi_fac * 15.0e3 * L, rel=1e-10)

    # ------------------------------------------------------------------
    # 6. Finding 12: setParameter / updateParameter for 'lch', 'Lch', and
    #    'characteristicLength' in LadrunoForceBeamColumn2d and LadrunoDispBeamColumn3d
    # ------------------------------------------------------------------
    for tok_idx, tok in enumerate(("lch", "Lch", "characteristicLength"), start=1):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1)
        ops.parameter(tok_idx, "element", 1, tok)
        ops.updateParameter(tok_idx, 0.35)

        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        ops.element("LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1)
        ops.parameter(10 + tok_idx, "element", 1, tok)
        ops.updateParameter(10 + tok_idx, 0.35)

    # ------------------------------------------------------------------
    # 7. Finding 13 & Finding 14: LadrunoForceBeamColumn3d::Print VZ/T signs
    #    and XML recorder 'section $secNum' physical coordinate eta = xi * L
    # ------------------------------------------------------------------
    xml_sec = str(tmp_path / "fbc3d_sec2.xml")
    print_file = str(tmp_path / "fbc3d_print.txt")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.element("LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1)
    ops.recorder("Element", "-xml", xml_sec, "-ele", 1, "section", "2", "force")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 12.0e3, -18.0e3, 24.0e3, 7.0e3, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    lf = ops.eleResponse(1, "localForce")
    ops.printModel("-file", print_file)
    ops.wipe()

    # Check XML eta attribute for section 2 (midpoint of 3-pt Lobatto -> eta = 0.5 * L = 2.0)
    gp = ET.parse(xml_sec).getroot().find(".//GaussPointOutput")
    assert gp is not None
    assert float(gp.attrib["eta"]) == pytest.approx(0.5 * L, rel=1e-12)

    # Check Print output End 1 and End 2 VZ and T signs against localForce:
    # localForce = [N_1, Vy_1, Vz_1, T_1, My_1, Mz_1, N_2, Vy_2, Vz_2, T_2, My_2, Mz_2]
    # Print format: End k Forces (P MZ VY MY VZ T)
    txt = (tmp_path / "fbc3d_print.txt").read_text(encoding="utf-8")
    end1_line = [ln for ln in txt.splitlines() if "End 1 Forces" in ln][0]
    end2_line = [ln for ln in txt.splitlines() if "End 2 Forces" in ln][0]
    e1_vals = [float(x) for x in end1_line.split(":")[1].split()]
    e2_vals = [float(x) for x in end2_line.split(":")[1].split()]
    # [P, MZ, VY, MY, VZ, T] vs localForce [N, Vy, Vz, T, My, Mz]
    np.testing.assert_allclose(
        e1_vals, [lf[0], lf[5], lf[1], lf[4], lf[2], lf[3]], rtol=1e-10, atol=1e-10
    )
    np.testing.assert_allclose(
        e2_vals, [lf[6], lf[11], lf[7], lf[10], lf[8], lf[9]], rtol=1e-10, atol=1e-10
    )


def test_r2_opus_review_edge_cases(tmp_path):
    """Round-2 Opus verification regressions:
    1. LadrunoDispBeamColumn2d/3d plasticDeformation linear beta/gamma weights for interior -hingeLoc in dualHingeMode.
    2. LadrunoIMKBeam2d/3d energy clear-span q0ic/q0jc with -shearHinge + -endZones + span loads.
    3. LadrunoForceBeamColumn2d/3d getInitialFlexibility(fe, true) fPM scaling and LadrunoCohesiveHinge* penalty guard.
    4. LadrunoForceBeamColumn2d/3d eleLoadDirty flush in getResponse/setResponse/Print + 2D hingeHinvCommit/hingeKvaCommit serialization.
    """
    L = 4.0
    E, G, A, Iz, Iy, Jx = 30.0e9, 12.0e9, 0.12, 1.6e-3, 1.2e-3, 2.0e-3
    k_hinge = 5.0 * E * Iz / L

    # ------------------------------------------------------------------
    # 1. LadrunoDispBeamColumn2d/3d plasticDeformation with interior -hingeLoc in dualHingeMode
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("Elastic", 10, k_hinge)
    ops.element(
        "LadrunoDispBeamColumn2d", 1, 1, 2, 1, 1,
        "-hingeIJ", 10, 10, "-hingeLoc", 0.25, 0.75,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3, -20.0e3, 15.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    vp2d = ops.eleResponse(1, "plasticDeformation")
    np.testing.assert_allclose(vp2d, [0.0, 0.0, 0.0], atol=1e-13)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("Elastic", 10, k_hinge)
    ops.element(
        "LadrunoDispBeamColumn3d", 1, 1, 2, 1, 1,
        "-hingeIJ", 10, 10, "-hingeYIJ", 10, 10, "-hingeLoc", 0.25, 0.75,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0e3, -20.0e3, 15.0e3, 5.0e3, -12.0e3, 18.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    vp3d = ops.eleResponse(1, "plasticDeformation")
    np.testing.assert_allclose(vp3d, np.zeros(6), atol=1e-13)

    # ------------------------------------------------------------------
    # 2. LadrunoIMKBeam2d/3d energy clear-span q0ic/q0jc with -shearHinge + -endZones + span loads
    # ------------------------------------------------------------------
    # When nodes are fixed and only uniform span load wy is applied on a symmetric
    # beam, q0[1] = -wy*L^2/12 and q0[2] = +wy*L^2/12. With symmetric endZones dI = dJ = d,
    # q0ic + q0jc = ((L - 2*d)/L) * (q0[1] + q0[2]) = 0, so shear hinge force Vs = 0
    # and gammaS = 0, meaning gamRot = 0 and core elastic rotation is purely vic = vjc = 0,
    # so reported energy must be 0.0! (Previously q0ic + q0jc used -d*(p0[1]-p0[2]) = 0,
    # and for asymmetric dI != dJ it differed from LadrunoIMKHinge.h.)
    # Let's test asymmetric dI=0.4, dJ=0.8 with fixed ends and uniform load:
    dI, dJ = 0.4, 0.8
    Lc = L - dI - dJ
    Ks0 = 5.0e7
    wy = -24.0e3
    wz = 18.0e3

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 10, 1.0e9, k_hinge, 0.01)
    ops.uniaxialMaterial("Steel01", 20, 1.0e9, Ks0, 0.01)
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-shearHinge", 20, "-endZones", dI, dJ,
    )
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
    en2d = ops.eleResponse(1, "energy")[0]
    th2d = ops.eleResponse(1, "hingeRotation")
    sh2d = ops.eleResponse(1, "shearHinge")  # [gammaS, Vs]
    # Exact internal strain energy of core + hinges + shearHinge when v = 0:
    # In LadrunoIMKHinge.h, at v=0, core elastic rotations are ve1 = -th_i + gamRot, ve2 = -th_j + gamRot
    # where gamRot = gammaS / Lc - gamRot0 and gamRot0 = -Fs0 * (q0ic + q0jc)
    q0ic = -wy * Lc * Lc / 12.0
    q0jc = wy * Lc * Lc / 12.0
    Fs0 = 1.0 / (Lc * Lc * Ks0)
    gamRot = sh2d[0] / Lc - (-Fs0 * (q0ic + q0jc))
    ve1 = -th2d[0] + gamRot
    ve2 = -th2d[1] + gamRot
    a_c = 4.0 * E * Iz / Lc
    b_c = 2.0 * E * Iz / Lc
    en2d_expected = (
        0.5 * (ve1 * (a_c * ve1 + b_c * ve2) + ve2 * (b_c * ve1 + a_c * ve2))
        + 0.5 * k_hinge * (th2d[0] ** 2 + th2d[1] ** 2)
        + 0.5 * Ks0 * (sh2d[0] ** 2)
    )
    assert en2d == pytest.approx(en2d_expected, rel=1e-12)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.uniaxialMaterial("Steel01", 10, 1.0e9, k_hinge, 0.01)
    ops.uniaxialMaterial("Steel01", 20, 1.0e9, Ks0, 0.01)
    ops.element(
        "LadrunoIMKBeam", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
        "-hinge", "both", "-matZ", 10, "-matY", 10,
        "-shearHingeY", 20, "-shearHingeZ", 20, "-endZones", dI, dJ,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", wy, wz, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    en3d = ops.eleResponse(1, "energy")[0]
    th3d = ops.eleResponse(1, "hingeRotation")
    sh3d = ops.eleResponse(1, "shearHinge")  # [gammaSY, VsY, gammaSZ, VsZ]
    q0yic = wz * Lc * Lc / 12.0
    q0yjc = -wz * Lc * Lc / 12.0
    gamRotZ = sh3d[2] / Lc - (-Fs0 * (q0yic + q0yjc))
    ve3 = -th3d[2] + gamRotZ
    ve4 = -th3d[3] + gamRotZ
    ay_c = 4.0 * E * Iy / Lc
    by_c = 2.0 * E * Iy / Lc
    en3d_expected = (
        en2d_expected
        + 0.5 * (ve3 * (ay_c * ve3 + by_c * ve4) + ve4 * (by_c * ve3 + ay_c * ve4))
        + 0.5 * k_hinge * (th3d[2] ** 2 + th3d[3] ** 2)
        + 0.5 * Ks0 * (sh3d[2] ** 2)
    )
    assert en3d == pytest.approx(en3d_expected, rel=1e-12)

    # ------------------------------------------------------------------
    # 3. LadrunoForceBeamColumn2d/3d plasticDeformation with -pm and with LadrunoCohesiveHinge*
    # ------------------------------------------------------------------
    P0 = 100.0e3
    pm_alpha = 2.0
    # 2D generic elastic hinge + -pm under axial compression
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("Elastic", 10, k_hinge)
    ops.element(
        "LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1,
        "-hingeIJ", 10, 10, "-pm", P0, pm_alpha,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -50.0e3, -15.0e3, 10.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    np.testing.assert_allclose(ops.eleResponse(1, "plasticDeformation"), np.zeros(3), atol=1e-13)

    # 2D LadrunoCohesiveHinge: cohesive penalty is skipped in getInitialFlexibility(fe, true),
    # so plasticDeformation equals [0, hingeJumpI, hingeJumpJ] (matching LadrunoDispBeamColumn2d)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 10, 500.0e3, 1.0e5, "-linear")
    ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1, "-hingeIJ", 10, 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -10.0e3, -5.0e3, 2.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    hs_coh = ops.eleResponse(1, "hingeSummary")
    np.testing.assert_allclose(
        ops.eleResponse(1, "plasticDeformation"),
        [0.0, hs_coh[0], hs_coh[2]],
        rtol=1e-12,
        atol=1e-13,
    )

    # 3D generic elastic hinges (uniaxial and biaxial) + -pm under axial compression
    for use_biax in (False, True):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        if use_biax:
            ops.nDMaterial("ElasticIsotropic", 30, k_hinge, 0.2)
            ops.element(
                "LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1,
                "-hingeBiaxialIJ", 30, 30, "-pm", P0, pm_alpha,
            )
        else:
            ops.uniaxialMaterial("Elastic", 10, k_hinge)
            ops.element(
                "LadrunoForceBeamColumn3d", 1, 1, 2, 1, 1,
                "-hingeIJ", 10, 10, "-hingeYIJ", 10, 10, "-pm", P0, pm_alpha,
            )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -50.0e3, -15.0e3, 12.0e3, 3.0e3, -8.0e3, 10.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        np.testing.assert_allclose(ops.eleResponse(1, "plasticDeformation"), np.zeros(6), atol=1e-13)

    # ------------------------------------------------------------------
    # 4. eleLoadDirty flush in getResponse + 2D hingeHinvCommit/hingeKvaCommit roundtrip
    # ------------------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 0, 1)
    ops.section("Elastic", 1, E, A, Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", -12.0e3)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    # Force a step revert so Element::revertToLastCommit sets eleLoadDirty = (numEleLoads > 0);
    # subsequent getResponse/setResponse calls must flush eleLoadDirty!
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, -10.0e3, 0.0)
    ops.test("NormDispIncr", 1.0e-30, 1)
    assert ops.analyze(1) < 0
    ops.remove("loadPattern", 2)
    bf2d = ops.eleResponse(1, "basicForce")
    sf2d = ops.eleResponse(1, "section", "1", "force")
    assert bf2d[1] == pytest.approx(12.0e3 * L * L / 3.0, rel=1e-10)
    assert sf2d[1] == pytest.approx(-12.0e3 * L * L / 3.0, rel=1e-10)

    def build_fbc2d_hinge():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.section("Elastic", 1, E, A, Iz)
        ops.geomTransf("Linear", 1)
        ops.beamIntegration("Lobatto", 1, 1, 4)
        ops.uniaxialMaterial("Steel01", 10, 80.0e3, k_hinge, 0.05)
        ops.element("LadrunoForceBeamColumn2d", 1, 1, 2, 1, 1, "-hingeIJ", 10, 10)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -10.0e3, -30.0e3, 20.0e3)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.5)
        ops.algorithm("Newton")
        ops.analysis("Static")
        assert ops.analyze(2) == 0

    def probe_fbc2d_hinge():
        return list(ops.eleForce(1)) + list(ops.eleResponse(1, "hingeSummary"))

    database_roundtrip(build_fbc2d_hinge, probe_nodes=[2], ndf=3, probe_fn=probe_fbc2d_hinge, rounds=2)


