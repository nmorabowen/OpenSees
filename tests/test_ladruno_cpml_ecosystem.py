"""Phase P2 Anchor 5 — Ecosystem, Lyapunov Stability, Stage 0->1 Transition,
LadrunoElementShell Contract, FileDatastore Checkpoint Serialization, Tcl Parity,
and Fail-Closed Parser Refusal Battery for ``LadrunoCPML`` (ELE_TAG_LadrunoCPML = 33023).
"""

import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import opensees as ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]


def _find_tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_tcl_exe()
needs_tcl_exe = pytest.mark.skipif(
    TCL_EXE is None,
    reason="classic OpenSees.exe not found in dist/bin",
)


def _ricker(t, f0, t0):
    tau = math.pi * f0 * (t - t0)
    u = tau * tau
    return (1.0 - 2.0 * u) * math.exp(-u)


# ===========================================================================
# Test Group 1 — EnergyBalanceRecorder Closure & Absorbed Energy Monotonicity
# ===========================================================================


def _run_2d_wave_energy_model(
    tmp_path,
    mode,
    alpha_max,
    use_rayleigh=False,
):
    """Build a 2D wave-propagation model with a FourNodeQuad ('quad') interior and a
    LadrunoCPML absorbing boundary in either '-box' or '-plane' mode.
    Attach EnergyBalanceRecorder, drive 60 transient steps under a Ricker load,
    and return (eb_data, int_p_ext, peak_step).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    E = 2.0e7
    nu = 0.25
    # Low inertia density (rho = 1e-8) paired with R0 = 0.9999 keeps the physical C-PML
    # damping rate sigma_max = (m+1)*cp/(2*L)*ln(1/R0) ~ 1160 s^-1 while bounding the
    # O(dt^2/8 * a^T M a) trapezoidal-vs-midpoint kinetic energy quadrature offset <= 3e-11.
    rho = 1.0e-8
    R0_eb = 0.9999
    thick = 1.0

    # 3x2 grid of 1m x 1m elements:
    #   x in [0, 1]: quad (FourNodeQuad) interior (2 elements)
    #   x in [1, 3]: LadrunoCPML absorbing layer (4 elements, L_pml = 2.0)
    nx, ny = 3, 2
    dx, dy = 1.0, 1.0

    def nid(ix, iy):
        return iy * (nx + 1) + ix + 1

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            ops.node(nid(ix, iy), ix * dx, iy * dy)

    # Fix outer PML boundary at x = 3.0
    for iy in range(ny + 1):
        ops.fix(nid(nx, iy), 1, 1)

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    ele_tag = 1
    cpml_tags = []
    for iy in range(ny):
        for ix in range(nx):
            n1 = nid(ix, iy)
            n2 = nid(ix + 1, iy)
            n3 = nid(ix + 1, iy + 1)
            n4 = nid(ix, iy + 1)
            if ix == 0:
                ops.element(
                    "quad",
                    ele_tag,
                    n1,
                    n2,
                    n3,
                    n4,
                    thick,
                    "PlaneStrain",
                    1,
                    0.0,
                    rho,
                    0.0,
                    0.0,
                )
            else:
                if mode == "box":
                    ops.element(
                        "LadrunoCPML",
                        ele_tag,
                        n1,
                        n2,
                        n3,
                        n4,
                        "-E",
                        E,
                        "-nu",
                        nu,
                        "-rho",
                        rho,
                        "-thick",
                        thick,
                        "-box",
                        0.0,
                        1.0,
                        0.0,
                        2.0,
                        "-L",
                        2.0,
                        2.0,
                        "-m",
                        2.0,
                        "-R0",
                        R0_eb,
                        "-kappaMax",
                        3.0,
                        "-alphaMax",
                        alpha_max,
                    )
                else:
                    ops.element(
                        "LadrunoCPML",
                        ele_tag,
                        n1,
                        n2,
                        n3,
                        n4,
                        "-E",
                        E,
                        "-nu",
                        nu,
                        "-rho",
                        rho,
                        "-thick",
                        thick,
                        "-plane",
                        "-x0",
                        1.0,
                        0.0,
                        "-n",
                        1.0,
                        0.0,
                        "-L",
                        2.0,
                        "-m",
                        2.0,
                        "-R0",
                        R0_eb,
                        "-kappaMax",
                        3.0,
                        "-alphaMax",
                        alpha_max,
                    )
                cpml_tags.append(ele_tag)
            ele_tag += 1

    if use_rayleigh:
        ops.rayleigh(0.05, 1.0e-4, 2.0e-4, 1.0e-4)

    n_steps = 60
    dt = 0.002
    f0 = 18.0
    t0 = 16.0 * dt
    # Smooth Ricker pulse that is zero at t=0 and t=dt so firstRecord rectangle rule matches trapz
    times = [k * dt for k in range(n_steps + 2)]
    vals = []
    for k, t in enumerate(times):
        if k <= 1 or k > 36:
            vals.append(0.0)
        else:
            w = math.sin(math.pi * (k - 1) / 35.0) ** 2
            vals.append(1.0e4 * w * _ricker(t, f0, t0))

    peak_step = int(np.argmax(np.abs(vals)))

    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    # Apply load at the interior / C-PML interface node (ix=1, iy=1)
    ops.load(nid(1, 1), 1.0, 0.5)

    eb_file = str(tmp_path / f"eb_{mode}_{alpha_max:.1f}_{int(use_rayleigh)}.txt")
    ops.recorder("EnergyBalance", "-file", eb_file, "-precision", 16)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 20, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    p_ext_hist = [0.0]
    for _ in range(n_steps):
        assert ops.analyze(1, dt) == 0
        step_p = 0.0
        for ctag in cpml_tags:
            resp = ops.eleResponse(ctag, "energy")
            assert resp and len(resp) == 1
            step_p += float(resp[0])
        p_ext_hist.append(step_p)

    ops.wipe()

    eb_data = np.atleast_2d(np.loadtxt(eb_file))
    p_arr = np.asarray(p_ext_hist, dtype=float)
    int_p_ext = float(np.sum(0.5 * (p_arr[:-1] + p_arr[1:]) * dt))
    return eb_data, int_p_ext, peak_step


@pytest.mark.parametrize("mode", ["box", "plane"])
@pytest.mark.parametrize("alpha_max", [0.0, math.pi * 10.0])
@pytest.mark.parametrize("use_rayleigh", [False, True])
def test_group1_energy_balance_closure_and_absorbed_work(
    tmp_path, mode, alpha_max, use_rayleigh
):
    """Test Group 1: EnergyBalanceRecorder closure |(KE + IE + DW) - EW| / max(EW) <= 1e-10
    at every step after the Ricker pulse peaks, and strictly positive time-integrated
    C-PML absorbed work int_0^T P_ext(t) dt > 0 via ops.eleResponse(tag, 'energy').
    """
    eb_data, int_p_ext, peak_step = _run_2d_wave_energy_model(
        tmp_path, mode=mode, alpha_max=alpha_max, use_rayleigh=use_rayleigh
    )
    assert eb_data.shape == (60, 6)

    ke = eb_data[:, 0]
    ie = eb_data[:, 1]
    dw = eb_data[:, 2]
    ew = eb_data[:, 3]

    max_ew = float(np.max(np.abs(ew)))
    assert max_ew > 0.0

    for step_idx in range(peak_step, 60):
        rel_err = abs((ke[step_idx] + ie[step_idx] + dw[step_idx]) - ew[step_idx]) / max_ew
        assert rel_err <= 1.0e-10, (
            f"mode={mode}, alpha={alpha_max}, rayleigh={use_rayleigh}, "
            f"step={step_idx + 1}: rel_err={rel_err:.3e} > 1e-10"
        )

    assert int_p_ext > 0.0, (
        f"Expected positive C-PML absorbed work integral, got {int_p_ext}"
    )


# ===========================================================================
# Test Group 2 — Late-Time Lyapunov Stability & Zero Secular Drift (N = 2500)
# ===========================================================================


def test_group2_late_time_lyapunov_stability_2d_face_and_corner():
    """Drive a 2D model containing interior, face, and corner LadrunoCPML elements
    (alpha_max = pi * 15.0, kappa_max = 4.0, R0 = 1e-4) with a Ricker pulse active
    for the first 100 steps and let it ring down freely to N = 2500 steps.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    E = 5.0e7
    nu = 0.25
    rho = 2000.0
    alpha_max = math.pi * 15.0
    kappa_max = 4.0
    R0 = 1.0e-4

    # 2x2 mesh on [0, 2] x [0, 2]:
    #   [0, 1] x [0, 1] = interior (elem 1, quad)
    #   [1, 2] x [0, 1] = +x face C-PML (elem 2, LadrunoCPML)
    #   [0, 1] x [1, 2] = +y face C-PML (elem 3, LadrunoCPML)
    #   [1, 2] x [1, 2] = +x/+y corner C-PML (elem 4, LadrunoCPML)
    for j in range(3):
        for i in range(3):
            nid = j * 3 + i + 1
            ops.node(nid, float(i), float(j))
            if i == 2 or j == 2:
                ops.fix(nid, 1, 1)
            elif i == 0 and j == 0:
                ops.fix(nid, 1, 1)
            elif i == 0:
                ops.fix(nid, 1, 0)
            elif j == 0:
                ops.fix(nid, 0, 1)

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    ops.element("quad", 1, 1, 2, 5, 4, 1.0, "PlaneStrain", 1, 0.0, rho, 0.0, 0.0)

    cpml_common = [
        "-E", E, "-nu", nu, "-rho", rho,
        "-box", 0.0, 1.0, 0.0, 1.0,
        "-L", 1.0, 1.0,
        "-m", 2.0,
        "-R0", R0,
        "-kappaMax", kappa_max,
        "-alphaMax", alpha_max,
    ]
    ops.element("LadrunoCPML", 2, 2, 3, 6, 5, *cpml_common)
    ops.element("LadrunoCPML", 3, 4, 5, 8, 7, *cpml_common)
    ops.element("LadrunoCPML", 4, 5, 6, 9, 8, *cpml_common)

    N_total = 2500
    N_pulse = 100
    dt = 0.0005
    f0 = 15.0
    t0 = 40 * dt

    times = [k * dt for k in range(N_pulse + 2)]
    vals = []
    for k, t in enumerate(times):
        if k == 0 or k >= N_pulse:
            vals.append(0.0)
        else:
            env = math.sin(math.pi * k / float(N_pulse)) ** 2
            vals.append(1.0e5 * env * _ricker(t, f0, t0))

    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 1.0, 1.0)
    ops.load(2, 0.5, 0.0)
    ops.load(4, 0.0, 0.5)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    free_nodes = [2, 4, 5]
    u_norms = np.zeros(N_total)
    v_norms = np.zeros(N_total)
    e_mech = np.zeros(N_total)

    for step in range(N_total):
        assert ops.analyze(1, dt) == 0
        u_vec = []
        v_vec = []
        for nd in free_nodes:
            u_vec.extend(ops.nodeDisp(nd))
            v_vec.extend(ops.nodeVel(nd))
        u_arr = np.asarray(u_vec, dtype=float)
        v_arr = np.asarray(v_vec, dtype=float)
        u_norms[step] = float(np.linalg.norm(u_arr))
        v_norms[step] = float(np.linalg.norm(v_arr))
        e_mech[step] = 0.5 * rho * float(np.dot(v_arr, v_arr)) + 0.5 * E * float(np.dot(u_arr, u_arr))

    u_peak = float(np.max(u_norms))
    v_peak = float(np.max(v_norms))
    assert u_peak > 0.0 and v_peak > 0.0

    # Windowed energy over full oscillation periods (200 steps) after pulse turns off
    # decays monotonically with zero late-time growth
    post_pulse = e_mech[N_pulse:]
    win = 200
    n_win = len(post_pulse) // win
    win_means = [float(np.mean(post_pulse[w * win : (w + 1) * win])) for w in range(n_win)]
    for w in range(len(win_means) - 1):
        assert win_means[w + 1] <= win_means[w] + 1.0e-15 * win_means[0], (
            f"Non-monotonic energy decay in window {w} -> {w + 1}: "
            f"{win_means[w]} -> {win_means[w + 1]}"
        )

    assert u_norms[-1] < 1.0e-4 * u_peak, (
        f"2D final disp norm {u_norms[-1]:.3e} >= 1e-4 * peak ({u_peak:.3e})"
    )
    assert v_norms[-1] < 1.0e-4 * v_peak, (
        f"2D final vel norm {v_norms[-1]:.3e} >= 1e-4 * peak ({v_peak:.3e})"
    )


def test_group2_late_time_lyapunov_stability_3d_h8_corner():
    """Drive a 3D H8 corner LadrunoCPML element (alpha_max = pi * 15.0, kappa_max = 4.0,
    R0 = 1e-4) for N = 2500 steps (Ricker pulse active for first 100 steps) and verify
    monotonic decay and final norms < 1e-4 of peak.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    E = 5.0e7
    nu = 0.25
    rho = 2000.0
    alpha_max = math.pi * 15.0
    kappa_max = 4.0
    R0 = 1.0e-4

    # Single 3D H8 corner element on [1, 2]^3 (outside box [0, 1]^3 on all 3 axes!)
    coords = [
        (1.0, 1.0, 1.0),
        (2.0, 1.0, 1.0),
        (2.0, 2.0, 1.0),
        (1.0, 2.0, 1.0),
        (1.0, 1.0, 2.0),
        (2.0, 1.0, 2.0),
        (2.0, 2.0, 2.0),
        (1.0, 2.0, 2.0),
    ]
    for idx, (x, y, z) in enumerate(coords, start=1):
        ops.node(idx, x, y, z)
        if idx != 1:
            ops.fix(idx, 1, 1, 1)

    ops.element(
        "LadrunoCPML",
        1,
        1, 2, 3, 4, 5, 6, 7, 8,
        "-E", E,
        "-nu", nu,
        "-rho", rho,
        "-box", 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        "-L", 1.0, 1.0, 1.0,
        "-m", 2.0,
        "-R0", R0,
        "-kappaMax", kappa_max,
        "-alphaMax", alpha_max,
    )

    N_total = 2500
    N_pulse = 100
    dt = 0.0005
    f0 = 15.0
    t0 = 40 * dt

    times = [k * dt for k in range(N_pulse + 2)]
    vals = []
    for k, t in enumerate(times):
        if k == 0 or k >= N_pulse:
            vals.append(0.0)
        else:
            env = math.sin(math.pi * k / float(N_pulse)) ** 2
            vals.append(1.0e5 * env * _ricker(t, f0, t0))

    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0, -0.7, 0.5)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    M11 = np.asarray(ops.eleResponse(1, "mass"), dtype=float).reshape(24, 24)[:3, :3]
    K11 = np.asarray(ops.eleResponse(1, "stiffInitial"), dtype=float).reshape(24, 24)[:3, :3]

    u_norms = np.zeros(N_total)
    v_norms = np.zeros(N_total)
    e_mech = np.zeros(N_total)

    for step in range(N_total):
        assert ops.analyze(1, dt) == 0
        u_arr = np.asarray(ops.nodeDisp(1), dtype=float)
        v_arr = np.asarray(ops.nodeVel(1), dtype=float)
        u_norms[step] = float(np.linalg.norm(u_arr))
        v_norms[step] = float(np.linalg.norm(v_arr))
        e_mech[step] = 0.5 * float(v_arr @ M11 @ v_arr) + 0.5 * float(u_arr @ K11 @ u_arr)

    u_peak = float(np.max(u_norms))
    v_peak = float(np.max(v_norms))
    assert u_peak > 0.0 and v_peak > 0.0

    post_pulse = e_mech[N_pulse:]
    win = 200
    n_win = len(post_pulse) // win
    win_means = [float(np.mean(post_pulse[w * win : (w + 1) * win])) for w in range(n_win)]
    for w in range(len(win_means) - 1):
        assert win_means[w + 1] <= win_means[w] + 1.0e-15 * win_means[0]

    assert u_norms[-1] < 1.0e-4 * u_peak, (
        f"3D H8 corner final disp norm {u_norms[-1]:.3e} >= 1e-4 * peak ({u_peak:.3e})"
    )
    assert v_norms[-1] < 1.0e-4 * v_peak, (
        f"3D H8 corner final vel norm {v_norms[-1]:.3e} >= 1e-4 * peak ({v_peak:.3e})"
    )


# ===========================================================================
# Test Group 3 — Static Gravity Stage 0 -> 1 Transition
# ===========================================================================


@pytest.mark.parametrize("transition_token", ["stage", "updateMaterialStage"])
@pytest.mark.parametrize("load_mode", ["direct_b", "selfWeight"])
def test_group3a_2d_gravity_settlement_and_zero_transient_jump(
    transition_token, load_mode
):
    """2D soil-box model:
    1. Static gravity settlement in -stage 0 matches an all-elastic FourNodeQuad ('quad')
       reference model to <= 1e-13.
    2. Transition to stage 1 (via 'stage' and 'updateMaterialStage' parameter tokens)
       latches (u_0, f_0) and holds static equilibrium under TransientAnalysis
       (20 steps, no dynamic load) to <= 1e-13 with zero jump at step 1.
    3. Subsequent dynamic Ricker pulse propagates and absorbs cleanly around u_0.
    """
    E, nu, rho = 4.0e7, 0.3, 2000.0
    g = -9.81
    bx, by = 0.0, rho * g

    def build_2d_model(use_cpml):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        # 3x2 grid on [0, 3] x [0, 2]
        for j in range(3):
            for i in range(4):
                nid = j * 4 + i + 1
                ops.node(nid, float(i), float(j))
                if j == 0:
                    ops.fix(nid, 1, 1)
                elif i in (0, 3):
                    ops.fix(nid, 1, 0)

        ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

        b_init = (bx, by) if load_mode == "direct_b" else (0.0, 1.0)
        ele_tag = 1
        cpml_tags = []
        all_tags = []
        for j in range(2):
            for i in range(3):
                n1 = j * 4 + i + 1
                n2 = n1 + 1
                n3 = n2 + 4
                n4 = n1 + 4
                all_tags.append(ele_tag)
                if not use_cpml or (i == 1 and j == 1):
                    ops.element(
                        "quad",
                        ele_tag,
                        n1, n2, n3, n4,
                        1.0,
                        "PlaneStrain",
                        1,
                        0.0,
                        rho,
                        b_init[0],
                        b_init[1],
                    )
                else:
                    ops.element(
                        "LadrunoCPML",
                        ele_tag,
                        n1, n2, n3, n4,
                        "-E", E,
                        "-nu", nu,
                        "-rho", rho,
                        "-box", 1.0, 2.0, 1.0, 2.0,
                        "-L", 1.0, 1.0,
                        "-stage", 0,
                        "-b", b_init[0], b_init[1],
                    )
                    cpml_tags.append(ele_tag)
                ele_tag += 1

        if load_mode == "selfWeight":
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.eleLoad("-ele", *all_tags, "-type", "-selfWeight", bx, by)

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static", "-noWarnings")
        assert ops.analyze(1) == 0
        ops.loadConst("-time", 0.0)
        return cpml_tags

    build_2d_model(use_cpml=False)
    u_ref = {nd: list(ops.nodeDisp(nd)) for nd in ops.getNodeTags()}

    cpml_tags = build_2d_model(use_cpml=True)
    u_0 = {nd: list(ops.nodeDisp(nd)) for nd in ops.getNodeTags()}

    for nd in u_ref:
        for d in range(2):
            assert abs(u_0[nd][d] - u_ref[nd][d]) <= 1.0e-13, (
                f"Node {nd} dof {d}: u_cpml_stage0={u_0[nd][d]} != u_ref={u_ref[nd][d]}"
            )

    # Transition all LadrunoCPML elements from stage 0 -> stage 1
    for idx, ctag in enumerate(cpml_tags, start=1):
        assert ops.eleResponse(ctag, "stage") == [0]
        if transition_token == "stage":
            ops.parameter(100 + idx, "element", ctag, "stage")
        else:
            ops.parameter(100 + idx, "element", ctag, "updateMaterialStage", str(ctag))
        ops.updateParameter(100 + idx, 1)
        assert ops.eleResponse(ctag, "stage") == [1]

    # Switch to TransientAnalysis with NO additional dynamic load
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    dt = 0.002
    for step in range(1, 21):
        assert ops.analyze(1, dt) == 0
        for nd in u_0:
            cur = ops.nodeDisp(nd)
            for d in range(2):
                diff = abs(cur[d] - u_0[nd][d])
                assert diff <= 1.0e-13, (
                    f"Transient jump at step {step}, node {nd}, dof {d}: |u - u0| = {diff:.3e}"
                )

    # Now apply a dynamic Ricker pulse on top of the pre-stressed gravity state
    n_dyn = 90
    times = [(20 + k) * dt for k in range(n_dyn + 2)]
    vals = [
        5.0e3 * (math.sin(math.pi * k / 20.0) ** 2) * _ricker(k * dt, 15.0, 8 * dt)
        if 0 < k < 20
        else 0.0
        for k in range(n_dyn + 2)
    ]
    ops.timeSeries("Path", 2, "-time", *times, "-values", *vals)
    ops.pattern("Plain", 2, 2)
    ops.load(10, 0.0, -1.0)  # Top-center node

    du_hist = []
    p_ext_hist = [0.0]
    for _ in range(n_dyn):
        assert ops.analyze(1, dt) == 0
        du_vec = [
            ops.nodeDisp(nd)[d] - u_0[nd][d]
            for nd in u_0
            for d in range(2)
        ]
        du_hist.append(float(np.linalg.norm(du_vec)))
        p_step = sum(float(ops.eleResponse(ctag, "energy")[0]) for ctag in cpml_tags)
        p_ext_hist.append(p_step)

    assert max(du_hist) > 1.0e-6
    assert du_hist[-1] < 0.25 * max(du_hist)
    int_p = sum(0.5 * (p_ext_hist[k] + p_ext_hist[k + 1]) * dt for k in range(n_dyn))
    assert int_p > 0.0


def test_group3a_3d_brick_self_weight_and_zero_transient_jump():
    """3D model with stdBrick interior and LadrunoCPML exterior using BrickSelfWeight:
    1. Static gravity settlement in -stage 0 matches all-stdBrick reference to <= 1e-13.
    2. Stage 0 -> 1 transition holds equilibrium under 20 transient steps to <= 1e-13.
    """
    E, nu, rho = 5.0e7, 0.25, 2000.0
    bz = -9.81 * rho

    def build_3d_stack(use_cpml):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for k in range(3):
            z = float(k)
            coords = [
                (0.0, 0.0, z),
                (1.0, 0.0, z),
                (1.0, 1.0, z),
                (0.0, 1.0, z),
            ]
            for i, (x, y, zc) in enumerate(coords, start=1):
                nid = k * 4 + i
                ops.node(nid, x, y, zc)
                if k == 0:
                    ops.fix(nid, 1, 1, 1)
                else:
                    ops.fix(nid, 1, 1, 0)

        ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
        if use_cpml:
            ops.element(
                "LadrunoCPML",
                1,
                1, 2, 3, 4, 5, 6, 7, 8,
                "-E", E,
                "-nu", nu,
                "-rho", rho,
                "-box", 0.0, 1.0, 0.0, 1.0, 1.0, 2.0,
                "-L", 1.0, 1.0, 1.0,
                "-stage", 0,
                "-b", 0.0, 0.0, bz,
            )
        else:
            ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 0.0, 0.0, bz)

        ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, 1, 0.0, 0.0, bz)

        ops.timeSeries("Linear", 1, "-factor", 2.5)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, 2, "-type", "-BrickSelfWeight")

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static", "-noWarnings")
        assert ops.analyze(1) == 0
        ops.loadConst("-time", 0.0)

    build_3d_stack(use_cpml=False)
    u_ref = {nd: list(ops.nodeDisp(nd)) for nd in range(1, 13)}

    build_3d_stack(use_cpml=True)
    u_0 = {nd: list(ops.nodeDisp(nd)) for nd in range(1, 13)}

    for nd in range(1, 13):
        for d in range(3):
            assert abs(u_0[nd][d] - u_ref[nd][d]) <= 1.0e-13

    ops.parameter(1, "element", 1, "stage")
    ops.updateParameter(1, 1)
    assert ops.eleResponse(1, "stage") == [1]

    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    for _ in range(20):
        assert ops.analyze(1, 0.002) == 0
        for nd in range(1, 13):
            cur = ops.nodeDisp(nd)
            for d in range(3):
                assert abs(cur[d] - u_0[nd][d]) <= 1.0e-13


def test_group3b_non_hijacking_of_soil_update_material_stage(capfd):
    """Verify that in a domain containing a LadrunoCPML element (tag 10, initialized
    at -stage 0, added FIRST in the domain) and a soil quad element with
    PressureIndependMultiYield material tag 99:
    1. Calling ops.updateMaterialStage('-material', 99, '-stage', 1) is NOT hijacked
       by LadrunoCPML tag 10 (MaterialStageParameter's domain scan skips tag 10 and
       updates soil material 99 without any 'no effect' warning, while tag 10 remains
       at stage == 0).
    2. Updating LadrunoCPML tag 10's own stage parameter transitions tag 10 to stage 1.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx, (x, y) in enumerate(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (2.0, 0.0), (2.0, 1.0)],
        start=1,
    ):
        ops.node(idx, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.fix(5, 1, 1)

    # Add LadrunoCPML (eleTag = 10) FIRST so MaterialStageParameter::setDomain() visits it first!
    ops.element(
        "LadrunoCPML",
        10,
        1, 2, 3, 4,
        "-E", 2.0e7,
        "-nu", 0.25,
        "-rho", 2000.0,
        "-plane",
        "-x0", 0.0, 0.0,
        "-n", 0.0, 1.0,
        "-L", 1.0,
        "-stage", 0,
    )

    # Add a soil element with PressureIndependMultiYield material tag 99 SECOND
    ops.nDMaterial("PressureIndependMultiYield", 99, 2, 2.0, 4.0e4, 1.2e5, 100.0, 0.1, 0.0, 100.0, 0.0, 20)
    ops.element("quad", 20, 2, 5, 6, 3, 1.0, "PlaneStrain", 99)

    assert ops.eleResponse(10, "stage") == [0]

    capfd.readouterr()
    ops.updateMaterialStage("-material", 99, "-stage", 1)
    err = capfd.readouterr().err
    # Must not warn "no effect with material tag 99" (which would happen if tag 99 weren't found),
    # and LadrunoCPML tag 10 must NOT have hijacked the scan or flipped its own stage!
    assert "no effect" not in err
    assert ops.eleResponse(10, "stage") == [0]

    # Also test element-level setParameter with mismatched tag 99 on element 10 is refused
    capfd.readouterr()
    ops.parameter(990, "element", 10, "updateMaterialStage", "99")
    err_mismatch = capfd.readouterr().err
    assert "no objects were able to identify parameter" in err_mismatch
    assert ops.eleResponse(10, "stage") == [0]

    # Updating tag 10's own stage parameter transitions LadrunoCPML tag 10 to stage 1
    ops.parameter(1010, "element", 10, "updateMaterialStage", "10")
    ops.updateParameter(1010, 1)
    assert ops.eleResponse(10, "stage") == [1]


# ===========================================================================
# Test Group 4 — LadrunoElementShell Contract Battery
# ===========================================================================


def test_group4_ki_cache_freeze_and_live_restore_drop(tmp_path):
    """Verify Ki is cached on first formation and stays frozen when parameter 'E' is
    halved (so ModifiedNewton -initial takes > 8 iterations), and that a live
    database restore drops Ki (WP-124 C6 oracle).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    E0 = 1.0e6
    for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
        ops.node(idx, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)

    ops.element(
        "LadrunoCPML",
        1,
        1, 2, 3, 4,
        "-E", E0,
        "-nu", 0.25,
        "-rho", 100.0,
        "-plane",
        "-x0", 0.0, 0.0,
        "-n", 1.0, 0.0,
        "-L", 1.0,
        "-stage", 0,
    )
    ops.parameter(1, "element", 1, "E")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 100.0, 50.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 200, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.integrator("LoadControl", 0.5)
    ops.analysis("Static", "-noWarnings")

    assert ops.analyze(1) == 0
    assert ops.testIter() <= 2

    # Halve E: Ki must remain frozen at E0 so ModifiedNewton -initial iterates > 8 times
    ops.updateParameter(1, 0.5 * E0)
    assert ops.analyze(1) == 0
    n_iter = ops.testIter()
    assert n_iter > 8, f"Expected > 8 iterations with frozen Ki, got {n_iter}"

    # Part 2: Live database restore drops Ki (WP-124 C6 oracle)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
        ops.node(idx, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.element(
        "LadrunoCPML",
        1,
        1, 2, 3, 4,
        "-E", E0,
        "-nu", 0.25,
        "-rho", 100.0,
        "-plane",
        "-x0", 0.0, 0.0,
        "-n", 1.0, 0.0,
        "-L", 1.0,
        "-stage", 0,
    )
    ops.parameter(1, "element", 1, "E")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 100.0, 50.0)
    db_prefix = str(tmp_path / "ki_live_cpml")
    ops.database("File", db_prefix)
    # Save BEFORE Ki is formed at E0
    ops.save(1)

    # Halve E and form Ki at 0.5 * E0
    ops.updateParameter(1, 0.5 * E0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 200, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.integrator("LoadControl", 0.5)
    ops.analysis("Static", "-noWarnings")
    assert ops.analyze(1) == 0

    # Restore into the LIVE domain: recvSelf restores E = E0 and MUST drop Ki (formed at 0.5*E0),
    # so the next ModifiedNewton -initial step re-forms Ki at E0 and converges in <= 2 iterations!
    ops.restore(1)
    assert ops.analyze(1) == 0
    assert ops.testIter() <= 2


@pytest.mark.parametrize("lumped", [False, True])
@pytest.mark.parametrize("stage", [0, 1])
def test_group4_ground_motion_inertia_load(lumped, stage):
    """Apply UniformExcitation ground acceleration a_g to a constrained LadrunoCPML
    element (both consistent and -lumped, stage 0 and stage 1) and verify the nodal
    reaction / unbalance equals +M * a_g to <= 1e-12, plus rigid-body relative
    acceleration == -a_g on an unconstrained element.
    """
    E, nu, rho = 1.0e7, 0.25, 2.5
    ag_x, ag_y = 1.8, -0.9

    # 1. Constrained element: reaction under transient step with UniformExcitation load
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    for idx, (x, y) in enumerate(coords, start=1):
        ops.node(idx, x, y)
        ops.fix(idx, 1, 1)

    args = [
        "LadrunoCPML", 1, 1, 2, 3, 4,
        "-E", E, "-nu", nu, "-rho", rho,
        "-box", 0.0, 0.0, 0.0, 1.0,
        "-L", 1.0, 1.0,
        "-kappaMax", 3.0,
        "-stage", stage,
    ]
    if lumped:
        args.append("-lumped")
    ops.element(*args)

    M = np.asarray(ops.eleResponse(1, "mass"), dtype=float).reshape(8, 8)
    ag_vec = np.array([ag_x, ag_y] * 4, dtype=float)
    expected_q = M @ ag_vec if not lumped else np.diag(M) * ag_vec

    ops.timeSeries("Constant", 1, "-factor", ag_x)
    ops.timeSeries("Constant", 2, "-factor", ag_y)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
    ops.pattern("UniformExcitation", 2, 2, "-accel", 2)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    assert ops.analyze(1, 0.001) == 0

    f_ele = np.asarray(ops.eleResponse(1, "force"), dtype=float)
    # With all nodes fixed (u = v = a = 0), resisting force is 0 + (-Qload) = +M * a_g
    assert np.max(np.abs(f_ele - expected_q)) <= 1.0e-12

    # 2. Unconstrained rigid-body translation under UniformExcitation -> relative accel == -a_g
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx, (x, y) in enumerate(coords, start=1):
        ops.node(idx, x, y)
    args_rb = [
        "LadrunoCPML", 1, 1, 2, 3, 4,
        "-E", E, "-nu", nu, "-rho", rho,
        "-box", 0.0, 1.0, 0.0, 1.0,
        "-L", 1.0, 1.0,
        "-stage", stage,
    ]
    if lumped:
        args_rb.append("-lumped")
    ops.element(*args_rb)
    ops.timeSeries("Constant", 1, "-factor", ag_x)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    for _ in range(3):
        assert ops.analyze(1, 0.005) == 0
    for nd in range(1, 5):
        acc = ops.nodeAccel(nd)
        assert abs(acc[0] - (-ag_x)) <= 1.0e-12
        assert abs(acc[1]) <= 1.0e-12


def test_group4_rayleigh_damping_all_four_factors(tmp_path):
    """Enable ops.rayleigh(0.1, 0.001, 0.002, 0.003) on LadrunoCPML, run transient
    steps, and verify C = C_cpml + alphaM*M + betaK*K_eff + betaK0*K0 + betaKc*Kc
    where Kc refreshes on commitState(), getResistingForceIncInertia() includes
    both C-PML and Rayleigh damping forces cleanly, and base Element response
    vocabulary (globalForce, dampingForce, dynamicForce, inertialForce) + XML output
    are well-formed.
    """
    import xml.etree.ElementTree as ET

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    E, nu, rho = 1.0e7, 0.25, 200.0
    for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
        ops.node(idx, x, y)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)

    ops.element(
        "LadrunoCPML",
        1,
        1, 2, 3, 4,
        "-E", E,
        "-nu", nu,
        "-rho", rho,
        "-plane",
        "-x0", 0.0, 0.0,
        "-n", 1.0, 0.0,
        "-L", 1.0,
        "-kappaMax", 2.0,
        "-alphaMax", 20.0,
    )

    alphaM, betaK, betaK0, betaKc = 0.1, 0.001, 0.002, 0.003
    ops.rayleigh(alphaM, betaK, betaK0, betaKc)

    xml_path = tmp_path / "cpml_global_force.xml"
    ops.recorder("Element", "-xml", str(xml_path), "-ele", 1, "globalForce")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1000.0, -400.0)
    ops.load(3, 800.0, 300.0)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    for _ in range(3):
        assert ops.analyze(1, 0.002) == 0

    K_eff = np.asarray(ops.eleResponse(1, "stiff"), dtype=float).reshape(8, 8)
    K_0 = np.asarray(ops.eleResponse(1, "stiffInitial"), dtype=float).reshape(8, 8)
    M = np.asarray(ops.eleResponse(1, "mass"), dtype=float).reshape(8, 8)
    C_total = np.asarray(ops.eleResponse(1, "damp"), dtype=float).reshape(8, 8)

    # K_eff at dt=0.002 must genuinely differ from K_0
    assert np.max(np.abs(K_eff - K_0)) > 1.0e-2 * np.max(np.abs(K_0))

    # Verify base Element response vocabulary (globalForce, dampingForce, dynamicForce, inertialForce)
    f_stat = np.asarray(ops.eleResponse(1, "force"), dtype=float)
    f_glob = np.asarray(ops.eleResponse(1, "globalForce"), dtype=float)
    f_damp = np.asarray(ops.eleResponse(1, "dampingForce"), dtype=float)
    f_dyn = np.asarray(ops.eleResponse(1, "dynamicForce"), dtype=float)
    f_iner = np.asarray(ops.eleResponse(1, "inertialForce"), dtype=float)
    assert f_glob.shape == (8,) and np.max(np.abs(f_glob - f_stat)) <= 1.0e-14
    assert f_damp.shape == (8,) and np.linalg.norm(f_damp) > 0.0
    assert f_dyn.shape == (8,) and f_iner.shape == (8,)

    # Verify parameter update on E refreshes Kc on commitState() while keeping K_0 frozen
    ops.parameter(10, "element", 1, "E")
    ops.updateParameter(10, 1.5 * E)
    assert ops.analyze(1, 0.002) == 0
    K_eff_after = np.asarray(ops.eleResponse(1, "stiff"), dtype=float).reshape(8, 8)
    K_0_after = np.asarray(ops.eleResponse(1, "stiffInitial"), dtype=float).reshape(8, 8)
    C_total_after = np.asarray(ops.eleResponse(1, "damp"), dtype=float).reshape(8, 8)
    assert np.max(np.abs(K_0_after - K_0)) <= 1.0e-14
    assert np.max(np.abs(K_eff_after - K_eff)) > 0.1 * np.max(np.abs(K_eff))

    # Also verify rho parameter update scales M cleanly
    ops.parameter(11, "element", 1, "rho")
    ops.updateParameter(11, 2.0 * rho)
    M_after = np.asarray(ops.eleResponse(1, "mass"), dtype=float).reshape(8, 8)
    assert np.max(np.abs(M_after - 2.0 * M)) <= 1.0e-12 * np.max(np.abs(M))

    # Flush XML recorder via wipe() and verify XML tag balance
    ops.wipe()
    ET.parse(xml_path)

    # Compare against an identical model without Rayleigh damping to isolate C_cpml
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
        ops.node(idx, x, y)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)
    ops.element(
        "LadrunoCPML",
        1,
        1, 2, 3, 4,
        "-E", E,
        "-nu", nu,
        "-rho", rho,
        "-plane",
        "-x0", 0.0, 0.0,
        "-n", 1.0, 0.0,
        "-L", 1.0,
        "-kappaMax", 2.0,
        "-alphaMax", 20.0,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1000.0, -400.0)
    ops.load(3, 800.0, 300.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    for _ in range(3):
        assert ops.analyze(1, 0.002) == 0
    C_cpml = np.asarray(ops.eleResponse(1, "damp"), dtype=float).reshape(8, 8)

    # After commitState(), Kc == K_eff
    C_expected = C_cpml + alphaM * M + (betaK + betaKc) * K_eff + betaK0 * K_0
    assert np.max(np.abs(C_total - C_expected)) <= 1.0e-11 * np.max(np.abs(C_expected))

    # And after E -> 1.5*E + commitState(), Kc refreshed to K_eff_after while K_0 stayed frozen
    ops.parameter(10, "element", 1, "E")
    ops.updateParameter(10, 1.5 * E)
    assert ops.analyze(1, 0.002) == 0
    C_cpml_after = np.asarray(ops.eleResponse(1, "damp"), dtype=float).reshape(8, 8)
    C_expected_after = C_cpml_after + alphaM * M + (betaK + betaKc) * K_eff_after + betaK0 * K_0
    assert np.max(np.abs(C_total_after - C_expected_after)) <= 1.0e-11 * np.max(np.abs(C_expected_after))


def test_group4_revert_to_last_commit_and_revert_to_start_protection():
    """Verify skipNextUpdate_ protects LadrunoCPML during Domain::revertToLastCommit()
    (triggered by an aborted trial step in DirectIntegrationAnalysis) and
    Domain::revertToStart() (triggered by ops.reset()):
    1. Run 3 transient steps, force an aborted trial step (ModifiedNewton -initial with
       impossible 1e-30 tolerance in 1 iteration, which advances trial kinematics, updates
       the element, fails convergence, and calls Domain::revertToLastCommit() ->
       Element::revertToLastCommit() + Domain::update()), and verify element force and
       3 subsequent steps are bit-identical (<= 1e-14) to a clean run that never took
       the aborted trial step.
    2. Call ops.reset() (Domain::revertToStart()) and verify re-running from t=0
       reproduces the trajectory to <= 1e-14.
    """
    def setup_and_run_3_steps():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
            ops.node(idx, x, y)
        ops.fix(1, 1, 1)
        ops.fix(4, 1, 1)
        ops.element(
            "LadrunoCPML",
            1,
            1, 2, 3, 4,
            "-E", 2.0e7,
            "-nu", 0.25,
            "-rho", 1500.0,
            "-box", 0.0, 0.0, 0.0, 0.0,
            "-L", 1.0, 1.0,
            "-kappaMax", 3.0,
            "-alphaMax", 25.0,
        )
        ops.timeSeries("Trig", 1, 0.0, 1.0, 0.02)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 5000.0, -2000.0)
        ops.load(3, 3000.0, 1500.0)

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient", "-noWarnings")
        for _ in range(3):
            assert ops.analyze(1, 0.001) == 0

    # Clean reference run: 3 steps + 3 more steps
    setup_and_run_3_steps()
    f_step3_ref = list(ops.eleResponse(1, "force"))
    traj_ref = []
    for _ in range(3):
        assert ops.analyze(1, 0.001) == 0
        traj_ref.append(
            list(ops.nodeDisp(2)) + list(ops.nodeVel(2)) + list(ops.eleResponse(1, "force"))
        )

    # Perturbed + aborted trial step (triggers Domain::revertToLastCommit -> Domain::update)
    setup_and_run_3_steps()
    ops.test("NormDispIncr", 1.0e-30, 1, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    # This 1-iteration step takes a non-zero trial increment, updates LadrunoCPML,
    # fails NormDispIncr(1e-30, 1), and calls Domain::revertToLastCommit() + Integrator::revertToLastStep()
    rc = ops.analyze(1, 0.001)
    assert rc < 0

    f_step3_rev = list(ops.eleResponse(1, "force"))
    assert np.max(np.abs(np.array(f_step3_rev) - np.array(f_step3_ref))) <= 1.0e-14

    # Restore normal Newton solver and continue for 3 steps
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Newton")
    traj_rev = []
    for _ in range(3):
        assert ops.analyze(1, 0.001) == 0
        traj_rev.append(
            list(ops.nodeDisp(2)) + list(ops.nodeVel(2)) + list(ops.eleResponse(1, "force"))
        )
    assert np.max(np.abs(np.array(traj_rev) - np.array(traj_ref))) <= 1.0e-14

    # Test ops.reset() (Domain::revertToStart()): reset domain to start and re-run 3 steps
    ops.reset()
    assert np.max(np.abs(np.array(ops.eleResponse(1, "force")))) <= 1.0e-15
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    for _ in range(3):
        assert ops.analyze(1, 0.001) == 0
    f_restart = list(ops.eleResponse(1, "force"))
    assert np.max(np.abs(np.array(f_restart) - np.array(f_step3_ref))) <= 1.0e-14


# ===========================================================================
# Test Group 5 — Complete sendSelf / recvSelf FileDatastore Checkpoint
#                Round-Trip Across Topologies (Q4, T6, H8, TET10, BTET10)
# ===========================================================================


def _build_2d_checkpoint_model():
    """2D model containing both Q4 (-box, with stage 0->1 latched u0/f0) and T6 (-plane)
    LadrunoCPML elements plus Rayleigh damping (including betaKc != 0).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    q4_coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    for i, (x, y) in enumerate(q4_coords, start=1):
        ops.node(i, x, y)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)

    t6_verts = [(3.0, 0.0), (2.0, 1.0), (2.0, 0.0)]
    t6_mids = [(2.5, 0.5), (2.0, 0.5), (2.5, 0.0)]
    for i, (x, y) in enumerate(t6_verts + t6_mids, start=5):
        ops.node(i, x, y)
    ops.fix(6, 1, 1)
    ops.fix(7, 1, 1)
    ops.fix(9, 1, 1)

    E, nu, rho = 2.5e7, 0.25, 1800.0
    ops.element(
        "LadrunoCPML",
        1,
        1, 2, 3, 4,
        "-E", E,
        "-nu", nu,
        "-rho", rho,
        "-box", 0.0, 0.0, 0.0, 0.0,
        "-L", 1.0, 1.0,
        "-stage", 0,
        "-b", 100.0, -200.0,
    )
    ops.element(
        "LadrunoCPML",
        2,
        5, 6, 7, 8, 9, 10,
        "-shape", "T6",
        "-E", E,
        "-nu", nu,
        "-rho", rho,
        "-plane",
        "-x0", 2.0, 0.0,
        "-n", 1.0, 0.0,
        "-L", 1.0,
    )
    ops.rayleigh(0.05, 1.0e-4, 1.0e-4, 2.0e-4)


def _build_3d_checkpoint_model(topo):
    """3D model for H8, TET10, or BTET10 (specifically testing 10-node elements
    that previously triggered the FileDatastore ID(10) collision).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    E, nu, rho = 3.0e7, 0.25, 2000.0

    if topo == "H8":
        coords = [
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
        ]
        fixed = [1, 4, 5, 8]
        free_probe = 6
        for i, xyz in enumerate(coords, start=1):
            ops.node(i, *xyz)
        for fn in fixed:
            ops.fix(fn, 1, 1, 1)
        ops.element(
            "LadrunoCPML",
            1,
            *range(1, 9),
            "-E", E, "-nu", nu, "-rho", rho,
            "-box", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            "-L", 1.0, 1.0, 1.0,
        )
    else:
        verts = [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 0.0),
        ]
        edges = [(0, 1), (1, 2), (2, 0), (0, 3), (2, 3), (1, 3)]
        mids = [
            tuple(0.5 * (verts[a][d] + verts[b][d]) for d in range(3))
            for a, b in edges
        ]
        all_pts = verts + mids
        for i, xyz in enumerate(all_pts, start=1):
            ops.node(i, *xyz)
        for fn in [2, 3, 4, 6, 9, 10]:
            ops.fix(fn, 1, 1, 1)
        free_probe = 1
        ops.element(
            "LadrunoCPML",
            1,
            *range(1, 11),
            "-shape", topo,
            "-E", E, "-nu", nu, "-rho", rho,
            "-box", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            "-L", 1.0, 1.0, 1.0,
        )
    return free_probe


@pytest.mark.parametrize("restore_mode", ["live", "fresh"])
def test_group5_checkpoint_roundtrip_2d_q4_and_t6_with_stage_latch(
    tmp_path, restore_mode
):
    """Test FileDatastore save(1) / restore(1) on 2D Q4 (-box, stage 0->1 latched)
    + T6 (-plane) with Rayleigh damping: 4 transient steps before save(1), 4 steps
    after restore(1), asserting bit-identical (<= 1e-14) kinematics and element forces.
    """
    db_prefix = str(tmp_path / f"cpml_2d_{restore_mode}")

    def configure_loads_and_stage():
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static", "-noWarnings")
        assert ops.analyze(1) == 0
        ops.parameter(1, "element", 1, "stage")
        ops.updateParameter(1, 1)

        ops.wipeAnalysis()
        ops.timeSeries("Trig", 1, 0.0, 1.0, 0.02)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 2000.0, -1000.0)
        ops.load(5, 1500.0, 800.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient", "-noWarnings")

    _build_2d_checkpoint_model()
    configure_loads_and_stage()

    dt = 0.001
    for _ in range(4):
        assert ops.analyze(1, dt) == 0

    ops.database("File", db_prefix)
    assert ops.save(1) in (0, None)

    ref_after = []
    for _ in range(4):
        assert ops.analyze(1, dt) == 0
        snap = (
            list(ops.nodeDisp(2))
            + list(ops.nodeVel(2))
            + list(ops.nodeAccel(2))
            + list(ops.nodeDisp(5))
            + list(ops.eleResponse(1, "force"))
            + list(ops.eleResponse(2, "force"))
        )
        ref_after.append(snap)

    if restore_mode == "fresh":
        _build_2d_checkpoint_model()
        configure_loads_and_stage()
        ops.database("File", db_prefix)

    assert ops.restore(1) in (0, None)
    ops.domainChange()

    got_after = []
    for _ in range(4):
        assert ops.analyze(1, dt) == 0
        snap = (
            list(ops.nodeDisp(2))
            + list(ops.nodeVel(2))
            + list(ops.nodeAccel(2))
            + list(ops.nodeDisp(5))
            + list(ops.eleResponse(1, "force"))
            + list(ops.eleResponse(2, "force"))
        )
        got_after.append(snap)

    assert np.max(np.abs(np.array(got_after) - np.array(ref_after))) <= 1.0e-14


@pytest.mark.parametrize("topo", ["H8", "TET10", "BTET10"])
@pytest.mark.parametrize("restore_mode", ["live", "fresh"])
def test_group5_checkpoint_roundtrip_3d_topologies(tmp_path, topo, restore_mode):
    """Test FileDatastore save(1) / restore(1) across 3D H8, TET10, and BTET10
    (verifying the 31-int fixed ID avoids the FileDatastore ID(10) collision).
    """
    db_prefix = str(tmp_path / f"cpml_3d_{topo}_{restore_mode}")

    def setup_transient(probe_node):
        ops.timeSeries("Trig", 1, 0.0, 1.0, 0.02)
        ops.pattern("Plain", 1, 1)
        ops.load(probe_node, 3000.0, -1500.0, 2000.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient", "-noWarnings")

    probe = _build_3d_checkpoint_model(topo)
    setup_transient(probe)

    dt = 0.001
    for _ in range(4):
        assert ops.analyze(1, dt) == 0

    ops.database("File", db_prefix)
    assert ops.save(1) in (0, None)

    ref_after = []
    for _ in range(4):
        assert ops.analyze(1, dt) == 0
        ref_after.append(
            list(ops.nodeDisp(probe))
            + list(ops.nodeVel(probe))
            + list(ops.nodeAccel(probe))
            + list(ops.eleResponse(1, "force"))
        )

    if restore_mode == "fresh":
        probe = _build_3d_checkpoint_model(topo)
        setup_transient(probe)
        ops.database("File", db_prefix)

    assert ops.restore(1) in (0, None)
    ops.domainChange()

    got_after = []
    for _ in range(4):
        assert ops.analyze(1, dt) == 0
        got_after.append(
            list(ops.nodeDisp(probe))
            + list(ops.nodeVel(probe))
            + list(ops.nodeAccel(probe))
            + list(ops.eleResponse(1, "force"))
        )

    assert np.max(np.abs(np.array(got_after) - np.array(ref_after))) <= 1.0e-14


# ===========================================================================
# Test Group 6 — Tcl Interpreter Parity (dist/bin/OpenSees.exe) &
#                Fail-Closed Parser Refusal Battery
# ===========================================================================


@needs_tcl_exe
def test_group6_tcl_interpreter_parity_2d_and_3d(tmp_path):
    """Run dist/bin/OpenSees.exe via subprocess.run(..., stdin=subprocess.DEVNULL)
    on generated Tcl scripts testing element LadrunoCPML in 2D (-box and -plane)
    and 3D (-box), and assert nodal displacements match opensees.pyd to <= 1e-13.
    """
    # --- 1. 2D (-box + -plane) Tcl vs Python ---
    tcl_2d = tmp_path / "run_2d.tcl"
    out_2d = tmp_path / "out_2d.txt"
    tcl_2d.write_text(
        f"""
wipe
model basic -ndm 2 -ndf 2
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 1.0
node 4 0.0 1.0
node 5 2.0 0.0
node 6 2.0 1.0
fix 1 1 1
fix 4 1 1
element LadrunoCPML 1 1 2 3 4 -E 2.5e7 -nu 0.25 -rho 1800.0 -box 0.0 0.0 0.0 1.0 -L 1.0 1.0 -m 2.0 -R0 1.0e-4 -kappaMax 3.0 -alphaMax 20.0
element LadrunoCPML 2 2 5 6 3 -E 2.5e7 -nu 0.25 -rho 1800.0 -plane -x0 1.0 0.0 -n 1.0 0.0 -L 1.0 -m 2.0 -R0 1.0e-4 -kappaMax 3.0 -alphaMax 20.0
timeSeries Linear 1
pattern Plain 1 1 {{
    load 5 1200.0 -450.0
    load 6 800.0 350.0
}}
constraints Plain
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-14 10 0
algorithm Linear
integrator Newmark 0.5 0.25
analysis Transient
analyze 5 0.001
set fp [open "{out_2d.as_posix()}" "w"]
puts $fp [format "%.17g %.17g %.17g %.17g" [nodeDisp 2 1] [nodeDisp 2 2] [nodeDisp 5 1] [nodeDisp 5 2]]
close $fp
wipe
""",
        encoding="utf-8",
    )
    res_2d = subprocess.run(
        [TCL_EXE, str(tcl_2d)],
        cwd=tmp_path,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert res_2d.returncode == 0, res_2d.stderr
    tcl_vals_2d = np.loadtxt(out_2d)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in [
        (1, (0.0, 0.0)), (2, (1.0, 0.0)), (3, (1.0, 1.0)),
        (4, (0.0, 1.0)), (5, (2.0, 0.0)), (6, (2.0, 1.0)),
    ]:
        ops.node(tag, x, y)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)
    ops.element(
        "LadrunoCPML", 1, 1, 2, 3, 4,
        "-E", 2.5e7, "-nu", 0.25, "-rho", 1800.0,
        "-box", 0.0, 0.0, 0.0, 1.0, "-L", 1.0, 1.0,
        "-m", 2.0, "-R0", 1.0e-4, "-kappaMax", 3.0, "-alphaMax", 20.0,
    )
    ops.element(
        "LadrunoCPML", 2, 2, 5, 6, 3,
        "-E", 2.5e7, "-nu", 0.25, "-rho", 1800.0,
        "-plane", "-x0", 1.0, 0.0, "-n", 1.0, 0.0, "-L", 1.0,
        "-m", 2.0, "-R0", 1.0e-4, "-kappaMax", 3.0, "-alphaMax", 20.0,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 1200.0, -450.0)
    ops.load(6, 800.0, 350.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    assert ops.analyze(5, 0.001) == 0
    py_vals_2d = np.array(list(ops.nodeDisp(2)) + list(ops.nodeDisp(5)), dtype=float)
    assert np.max(np.abs(tcl_vals_2d - py_vals_2d)) <= 1.0e-13

    # --- 2. 3D (-box) Tcl vs Python ---
    tcl_3d = tmp_path / "run_3d.tcl"
    out_3d = tmp_path / "out_3d.txt"
    tcl_3d.write_text(
        f"""
wipe
model basic -ndm 3 -ndf 3
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0
fix 1 1 1 1
fix 4 1 1 1
fix 5 1 1 1
fix 8 1 1 1
element LadrunoCPML 1 1 2 3 4 5 6 7 8 -E 3.0e7 -nu 0.25 -rho 2000.0 -box 0.0 0.0 0.0 1.0 0.0 1.0 -L 1.0 1.0 1.0 -m 2.0 -R0 1.0e-4 -kappaMax 3.0 -alphaMax 20.0
timeSeries Linear 1
pattern Plain 1 1 {{
    load 6 1500.0 -600.0 900.0
    load 7 1100.0 400.0 -500.0
}}
constraints Plain
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-14 10 0
algorithm Linear
integrator Newmark 0.5 0.25
analysis Transient
analyze 5 0.001
set fp [open "{out_3d.as_posix()}" "w"]
puts $fp [format "%.17g %.17g %.17g %.17g %.17g %.17g" [nodeDisp 6 1] [nodeDisp 6 2] [nodeDisp 6 3] [nodeDisp 7 1] [nodeDisp 7 2] [nodeDisp 7 3]]
close $fp
wipe
""",
        encoding="utf-8",
    )
    res_3d = subprocess.run(
        [TCL_EXE, str(tcl_3d)],
        cwd=tmp_path,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert res_3d.returncode == 0, res_3d.stderr
    tcl_vals_3d = np.loadtxt(out_3d)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    h8_coords = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    for idx, xyz in enumerate(h8_coords, start=1):
        ops.node(idx, *xyz)
    for fn in [1, 4, 5, 8]:
        ops.fix(fn, 1, 1, 1)
    ops.element(
        "LadrunoCPML", 1, *range(1, 9),
        "-E", 3.0e7, "-nu", 0.25, "-rho", 2000.0,
        "-box", 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,
        "-L", 1.0, 1.0, 1.0,
        "-m", 2.0, "-R0", 1.0e-4, "-kappaMax", 3.0, "-alphaMax", 20.0,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(6, 1500.0, -600.0, 900.0)
    ops.load(7, 1100.0, 400.0, -500.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    assert ops.analyze(5, 0.001) == 0
    py_vals_3d = np.array(list(ops.nodeDisp(6)) + list(ops.nodeDisp(7)), dtype=float)
    assert np.max(np.abs(tcl_vals_3d - py_vals_3d)) <= 1.0e-13


def _setup_2d_quad_nodes():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
        ops.node(idx, x, y)
    for idx in range(5, 10):
        ops.node(idx, float(idx), 0.0)


@pytest.mark.parametrize(
    "bad_args",
    [
        # Missing -box and -plane
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-L", 1.0],
        # Both -box and -plane
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-plane", "-x0", 0.0, 0.0, "-n", 1.0, 0.0],
        # Cross-mode: -box with -x0 or -n
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-x0", 0.0, 0.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-n", 1.0, 0.0],
        # Cross-mode: -plane with -Lx, -Ly, or multi-value -L
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-plane", "-x0", 0.0, 0.0, "-n", 1.0, 0.0, "-Lx", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-plane", "-x0", 0.0, 0.0, "-n", 1.0, 0.0, "-Ly", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-plane", "-x0", 0.0, 0.0, "-n", 1.0, 0.0, "-L", 1.0, 2.0],
        # Invalid material: E <= 0, nu <= -1.0, nu >= 0.499999, rho <= 0
        [1, 1, 2, 3, 4, "-E", 0.0, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-E", -5.0, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", -1.0, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.5, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 0.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", -10.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        # Invalid grading/geometry: L <= 0, m < 0, R0 <= 0, R0 >= 1, kappaMax < 1, alphaMax < 0
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 0.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", -1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-m", -0.5],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-R0", 0.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-R0", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-kappaMax", 0.9],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-alphaMax", -1.0],
        # Zero normal -n 0 0 and inverted box xmin > xmax
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-plane", "-x0", 0.0, 0.0, "-n", 0.0, 0.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 2.0, 1.0, 0.0, 1.0, "-L", 1.0],
        # Invalid node count in 2D (5 nodes) and -shape BT6 on 4 nodes
        [1, 1, 2, 3, 4, 5, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        [1, 1, 2, 3, 4, "-shape", "BT6", "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0],
        # Unknown flag
        [1, 1, 2, 3, 4, "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0, "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0, "-bogusFlag", 123],
    ],
)
def test_group6_fail_closed_parser_refusals_2d(bad_args):
    """Verify OPS_LadrunoCPML refuses invalid 2D flag combinations, materials,
    gradings, and topologies cleanly by raising OpenSeesError.
    """
    _setup_2d_quad_nodes()
    with pytest.raises(ops.OpenSeesError):
        ops.element("LadrunoCPML", *bad_args)


def test_group6_fail_closed_3d_and_setdomain_refusals(capfd):
    """Verify:
    1. 3D invalid node count (7 nodes), -shape BT6 in 3D, and -plane with -Lz raise OpenSeesError.
    2. Wrong node ndf (ndf=3 in 2D) and collinear/inverted element coordinates fail closed in
       LadrunoCPML::setDomain without crashing and zero out element matrices.
    """
    # 1. 3D parser refusals
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(1, 12):
        ops.node(i, float(i), float(i % 3), float(i % 2))

    # 7 nodes in 3D
    with pytest.raises(ops.OpenSeesError):
        ops.element(
            "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7,
            "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0,
            "-box", 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, "-L", 1.0,
        )

    # -shape BT6 in 3D
    with pytest.raises(ops.OpenSeesError):
        ops.element(
            "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7, 8,
            "-shape", "BT6",
            "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0,
            "-box", 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, "-L", 1.0,
        )

    # -plane with -Lz in 3D
    with pytest.raises(ops.OpenSeesError):
        ops.element(
            "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7, 8,
            "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0,
            "-plane", "-x0", 0.0, 0.0, 0.0, "-n", 1.0, 0.0, 0.0, "-Lz", 1.0,
        )

    # 2. Node with wrong ndf (ndf=3 in 2D) -> setDomain refusal
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for idx, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
        ops.node(idx, x, y)
    capfd.readouterr()
    ops.element(
        "LadrunoCPML", 1, 1, 2, 3, 4,
        "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0,
        "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0,
    )
    err_ndf = capfd.readouterr().err
    assert "DOFs, expected 2" in err_ndf
    assert np.max(np.abs(ops.eleResponse(1, "stiff"))) == 0.0

    # 3. Collinear / inverted 2D element coordinates -> setDomain refusal
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx, (x, y) in enumerate([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)], start=1):
        ops.node(idx, x, y)
    capfd.readouterr()
    ops.element(
        "LadrunoCPML", 1, 1, 2, 3, 4,
        "-E", 1.0e7, "-nu", 0.25, "-rho", 2000.0,
        "-box", 0.0, 1.0, 0.0, 1.0, "-L", 1.0,
    )
    err_inv = capfd.readouterr().err
    assert "failed C-PML kernel initialization" in err_inv
    assert np.max(np.abs(ops.eleResponse(1, "stiff"))) == 0.0
