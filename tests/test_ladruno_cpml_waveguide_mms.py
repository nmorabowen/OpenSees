"""Phase P2 Verification Suite (Worker 2.A) for LadrunoCPML (ELE_TAG_LadrunoCPML = 33023).

Covers:
1. Test Group 1: Exact Multi-Step Oracle Equivalence Through `opensees.pyd` Across All 18
   Reference Cases (`T3, Q4, T6, BT6, TET4, H8, TET10, BTET10, H20`), including inline
   (`-E/-nu/-rho`) vs `-mat` material specification, assembled `printA()` (`K_T_dyn`),
   `eleResponse("force")`, `eleResponse("basicForce")`, `eleResponse("energy")`, and
   in-engine finite-difference tangent verification.
2. Test Group 2 (Anchor 1): 1D P-Wave & S-Wave Waveguide Reflection Coefficient R <= 1%
   in both `-box` and `-plane` modes, `-box` vs `-plane` trajectory equivalence, and
   exact 37-degree oblique rotational frame invariance of `-plane` mode.
3. Test Group 3 (Anchor 2): Method of Manufactured Solutions (MMS) & Unstructured Multi-Element
   Patch Completeness across all 9 topologies (`T3, Q4, T6, BT6, TET4, H8, TET10, BTET10, H20`)
   in active C-PML corners/oblique layers, plus 2D spatial mesh refinement convergence
   (N in {4, 8, 16}) for unstructured `T3` and `Q4` C-PML meshes with rate p >= 1.85.
"""
from __future__ import annotations

import math
import os
from typing import Any

import numpy as np
import opensees as ops
import pytest

import ladruno_cpml_reference as ref

pytestmark = [pytest.mark.zone_a]

CASES_TXT_PATH = os.path.join(os.path.dirname(__file__), "ladruno_cpml_cases.txt")


# ============================================================================ #
# Helpers: Golden File Parser & OpenSees Element Builder                       #
# ============================================================================ #


def parse_cases_txt(path: str) -> list[dict[str, Any]]:
    """Parse `tests/ladruno_cpml_cases.txt` into structured dictionaries."""
    with open(path, "r", encoding="utf-8") as f:
        raw_lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]

    cases: list[dict[str, Any]] = []
    idx = 0
    while idx < len(raw_lines):
        toks = raw_lines[idx].split()
        assert toks[0] == "CASE", f"Expected CASE at line {idx}, got {raw_lines[idx]}"
        case_name = toks[1]
        idx += 1

        topo_toks = raw_lines[idx].split()
        assert topo_toks[0] == "Topology"
        topo = topo_toks[1]
        ndm = int(topo_toks[2])
        n_nodes = int(topo_toks[3])
        lumped = int(topo_toks[4])
        idx += 1

        mat_toks = raw_lines[idx].split()
        assert mat_toks[0] == "Material"
        E, nu, rho = float(mat_toks[1]), float(mat_toks[2]), float(mat_toks[3])
        idx += 1

        pml_toks = raw_lines[idx].split()
        assert pml_toks[0] == "PMLParams"
        m_deg, R0, kappa_max, alpha_max = (
            float(pml_toks[1]),
            float(pml_toks[2]),
            float(pml_toks[3]),
            float(pml_toks[4]),
        )
        idx += 1

        geom_toks = raw_lines[idx].split()
        assert geom_toks[0] == "GeomMode"
        geom_mode = geom_toks[1]
        n_gparams = int(geom_toks[2])
        geom_params = tuple(float(x) for x in geom_toks[3 : 3 + n_gparams])
        idx += 1

        coord_toks = raw_lines[idx].split()
        assert coord_toks[0] == "Coords"
        coords = np.array([float(x) for x in coord_toks[1:]], dtype=np.float64).reshape(
            n_nodes, ndm
        )
        idx += 1

        nm_toks = raw_lines[idx].split()
        assert nm_toks[0] == "NewmarkCoeffs"
        c1, c2, c3 = float(nm_toks[1]), float(nm_toks[2]), float(nm_toks[3])
        idx += 1

        ns_toks = raw_lines[idx].split()
        assert ns_toks[0] == "NumSteps"
        n_steps = int(ns_toks[1])
        idx += 1

        steps: list[dict[str, Any]] = []
        n_dof = n_nodes * ndm
        for _ in range(n_steps):
            st_toks = raw_lines[idx].split()
            assert st_toks[0] == "STEP"
            dt = float(st_toks[2])
            idx += 1

            step_data: dict[str, Any] = {"dt": dt}
            while raw_lines[idx] != "END_STEP":
                parts = raw_lines[idx].split()
                key = parts[0]
                vals = [float(x) for x in parts[1:]]
                if key in ("U", "V", "A", "R_int"):
                    step_data[key] = np.array(vals, dtype=np.float64)
                elif key in ("K_eff", "C_eff", "M_eff", "K_T_dyn"):
                    step_data[key] = np.array(vals, dtype=np.float64).reshape(n_dof, n_dof)
                elif key == "P_ext":
                    step_data[key] = vals[0]
                idx += 1
            idx += 1  # skip END_STEP
            steps.append(step_data)

        assert raw_lines[idx] == "END_CASE"
        idx += 1

        cases.append(
            {
                "name": case_name,
                "topo": topo,
                "ndm": ndm,
                "n_nodes": n_nodes,
                "lumped": lumped,
                "E": E,
                "nu": nu,
                "rho": rho,
                "pml": ref.PMLParams(m_deg=m_deg, R0=R0, kappa_max=kappa_max, alpha_max=alpha_max),
                "geom": ref.GeomSpec(mode=geom_mode, params=geom_params),
                "coords": coords,
                "c1": c1,
                "c2": c2,
                "c3": c3,
                "steps": steps,
            }
        )
    return cases


def build_cpml_element_args(
    ele_tag: int,
    node_tags: list[int],
    topo: str,
    ndm: int,
    E: float,
    nu: float,
    rho: float,
    pml: ref.PMLParams,
    geom: ref.GeomSpec,
    lumped: int,
    use_mat_tag: int | None = None,
    body_force: tuple[float, ...] | None = None,
) -> list[Any]:
    """Assemble arguments for `ops.element('LadrunoCPML', ...)`."""
    args: list[Any] = [ele_tag, *node_tags]
    if use_mat_tag is not None:
        args.extend(["-mat", use_mat_tag])
    else:
        args.extend(["-E", E, "-nu", nu, "-rho", rho])

    if geom.mode == "box":
        if ndm == 2:
            x_min, x_max, y_min, y_max, Lx, Ly = geom.params
            args.extend(["-box", x_min, x_max, y_min, y_max, "-Lx", Lx, "-Ly", Ly])
        else:
            x_min, x_max, y_min, y_max, z_min, z_max, Lx, Ly, Lz = geom.params
            args.extend(
                [
                    "-box",
                    x_min,
                    x_max,
                    y_min,
                    y_max,
                    z_min,
                    z_max,
                    "-Lx",
                    Lx,
                    "-Ly",
                    Ly,
                    "-Lz",
                    Lz,
                ]
            )
    elif geom.mode == "plane":
        if ndm == 2:
            x0, y0, nx, ny, L_pml = geom.params
            args.extend(["-plane", "-x0", x0, y0, "-n", nx, ny, "-L", L_pml])
        else:
            x0, y0, z0, nx, ny, nz, L_pml = geom.params
            args.extend(["-plane", "-x0", x0, y0, z0, "-n", nx, ny, nz, "-L", L_pml])
    else:
        raise ValueError(f"Unsupported geom.mode={geom.mode}")

    args.extend(
        [
            "-m",
            pml.m_deg,
            "-R0",
            pml.R0,
            "-kappaMax",
            pml.kappa_max,
            "-alphaMax",
            pml.alpha_max,
            "-shape",
            topo,
        ]
    )
    if lumped:
        args.append("-lumped")
    if body_force is not None:
        args.extend(["-body", *body_force])
    return args


def rel_err(a: np.ndarray | float, b: np.ndarray | float) -> float:
    """Compute relative norm error ||a - b|| / max(1.0, ||b||)."""
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    diff = float(np.linalg.norm(a_arr - b_arr))
    denom = max(1.0, float(np.linalg.norm(b_arr)))
    return diff / denom


def extract_oracle_r_disp(
    oracle: ref.ReferenceCPMLElement,
    R_int_ref: np.ndarray,
    M_eff_ref: np.ndarray,
    a_vec: np.ndarray,
) -> np.ndarray:
    """Extract `R_disp = f_K + f_disp = R_int - M_eff @ a - f_vel` from a ReferenceCPMLElement
    whose `evaluate_trial` was just called at `(u_vec, v_vec, a_vec)`.
    """
    ndm = oracle.ndm
    n_nodes = oracle.n_nodes
    f_vel_nodes = np.zeros((n_nodes, ndm), dtype=np.float64)
    for gp in oracle.gps:
        k0, k1, k2 = gp.kappa[0], gp.kappa[1], gp.kappa[2]
        a_vel_g = np.zeros(ndm, dtype=np.float64)
        for i in range(ndm):
            vg0 = gp.vel_g0[i].c0 * gp.vel_g0[i].x_trial + gp.vel_g0[i].c1 * gp.vel_g0[i].zeta_trial
            vg1 = gp.vel_g1[i].c0 * gp.vel_g1[i].x_trial + gp.vel_g1[i].c1 * gp.vel_g1[i].zeta_trial
            vg2 = gp.vel_g2[i].c0 * gp.vel_g2[i].x_trial + gp.vel_g2[i].c1 * gp.vel_g2[i].zeta_trial
            a_vel_g[i] = (k1 * k2) * vg0 + (k0 * k2) * vg1 + (k0 * k1) * vg2
        f_vel_nodes += (oracle.rho * gp.w_vol) * np.outer(gp.N, a_vel_g)
    return R_int_ref - (M_eff_ref @ a_vec) - f_vel_nodes.reshape(oracle.n_dof)


# ============================================================================ #
# Test Group 1: Exact Multi-Step Oracle Equivalence Across All 18 Cases        #
# ============================================================================ #


PARSED_CASES = parse_cases_txt(CASES_TXT_PATH)


def test_cases_txt_matches_reference_oracle() -> None:
    """Verify that the 18 parsed cases in `ladruno_cpml_cases.txt` cover all 9 topologies
    and match `ladruno_cpml_reference.py` to machine precision.
    """
    assert len(PARSED_CASES) == 18
    topos_seen = {c["topo"] for c in PARSED_CASES}
    assert topos_seen == {"T3", "Q4", "T6", "BT6", "TET4", "H8", "TET10", "BTET10", "H20"}

    ref_cases = ref.build_test_cases()
    assert len(ref_cases) == 18

    for pc, rc in zip(PARSED_CASES, ref_cases):
        assert pc["name"] == rc.name
        elem = ref.ReferenceCPMLElement(
            topo=pc["topo"],
            coords=pc["coords"],
            E=pc["E"],
            nu=pc["nu"],
            rho=pc["rho"],
            pml=pc["pml"],
            geom=pc["geom"],
            lumped=pc["lumped"],
        )
        for st in pc["steps"]:
            R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext = elem.evaluate_trial(
                st["U"], st["V"], st["A"], st["dt"], pc["c1"], pc["c2"], pc["c3"]
            )
            elem.commit()
            assert rel_err(R_int, st["R_int"]) <= 1.0e-12
            assert rel_err(K_eff, st["K_eff"]) <= 1.0e-12
            assert rel_err(C_eff, st["C_eff"]) <= 1.0e-12
            assert rel_err(M_eff, st["M_eff"]) <= 1.0e-12
            assert rel_err(K_T_dyn, st["K_T_dyn"]) <= 1.0e-12
            assert rel_err(P_ext, st["P_ext"]) <= 1.0e-12


@pytest.mark.parametrize("case_data", PARSED_CASES, ids=[c["name"] for c in PARSED_CASES])
@pytest.mark.parametrize("use_ndmat", [False, True], ids=["inline_mat", "ndmat_tag"])
def test_group1_multistep_oracle_equivalence_all_18_cases(
    case_data: dict[str, Any], use_ndmat: bool
) -> None:
    """Drive each of the 18 reference cases through multi-step Newmark(0.5, 0.25)
    TransientAnalysis in OpenSees (testing both inline E/nu/rho and `-mat $matTag`)
    and verify K_T_dyn (`printA`), `force`, `basicForce`, `stiff`, `damping`, `mass`,
    and `energy` against the Python reference oracle to <= 1e-11 relative error.
    """
    ndm = case_data["ndm"]
    n_nodes = case_data["n_nodes"]
    n_dof = n_nodes * ndm
    coords = case_data["coords"]
    E_in = float(case_data["E"])
    nu_in = float(case_data["nu"])
    rho_in = float(case_data["rho"])

    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)

    node_tags = list(range(1, n_nodes + 1))
    for idx, ntag in enumerate(node_tags):
        ops.node(ntag, *[float(x) for x in coords[idx]])

    mat_tag: int | None = None
    if use_ndmat:
        mat_tag = 10
        ops.nDMaterial("ElasticIsotropic", mat_tag, E_in, nu_in, rho_in)

    ele_args = build_cpml_element_args(
        ele_tag=1,
        node_tags=node_tags,
        topo=case_data["topo"],
        ndm=ndm,
        E=E_in,
        nu=nu_in,
        rho=rho_in,
        pml=case_data["pml"],
        geom=case_data["geom"],
        lumped=case_data["lumped"],
        use_mat_tag=mat_tag,
    )
    ops.element("LadrunoCPML", *ele_args)

    # Determine effective (E_oracle, nu_oracle) matched to the element's initial stiffness
    oracle = ref.ReferenceCPMLElement(
        topo=case_data["topo"],
        coords=coords,
        E=E_in,
        nu=nu_in,
        rho=rho_in,
        pml=case_data["pml"],
        geom=case_data["geom"],
        lumped=case_data["lumped"],
    )

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    u_c = np.zeros(n_dof, dtype=np.float64)
    v_c = np.zeros(n_dof, dtype=np.float64)
    a_c = np.zeros(n_dof, dtype=np.float64)

    for s_idx, step_info in enumerate(case_data["steps"]):
        dt = float(step_info["dt"])
        u_target = step_info["U"].copy()

        c1 = 1.0
        c2 = 2.0 / dt
        c3 = 4.0 / (dt * dt)

        # For a single unrestrained element where consistent mass is rank-deficient on
        # a rigid-body rotation mode (specifically TET4_box_face_z_cons), project u_target
        # onto the range of K_T_dyn so the free-DOF linear system is nonsingular on u_target.
        _, _, _, _, K_T_probe, _ = oracle.evaluate_trial(
            u_c, v_c, a_c, dt, c1, c2, c3
        )
        eigvals, eigvecs = np.linalg.eigh(0.5 * (K_T_probe + K_T_probe.T))
        max_eig = float(np.max(np.abs(eigvals)))
        null_mask = np.abs(eigvals) < 1.0e-10 * max_eig
        if np.any(null_mask):
            V_null = eigvecs[:, null_mask]
            u_target = u_target - V_null @ (V_null.T @ u_target)

        a_target = c3 * (u_target - u_c) - (4.0 / dt) * v_c - a_c
        v_target = v_c + 0.5 * dt * (a_c + a_target)

        R_int_ref, K_eff_ref, C_eff_ref, M_eff_ref, K_T_dyn_ref, P_ext_ref = (
            oracle.evaluate_trial(u_target, v_target, a_target, dt, c1, c2, c3)
        )
        R_disp_ref = extract_oracle_r_disp(oracle, R_int_ref, M_eff_ref, a_target)
        oracle.commit()

        # Apply constant nodal force F_ext = R_int_ref during this step
        ts_tag = 100 + s_idx
        pat_tag = 200 + s_idx
        ops.timeSeries("Constant", ts_tag)
        ops.pattern("Plain", pat_tag, ts_tag)
        for a_node, ntag in enumerate(node_tags):
            f_node = [float(R_int_ref[a_node * ndm + d]) for d in range(ndm)]
            ops.load(ntag, *f_node)

        rc = ops.analyze(1, dt)
        assert rc == 0, f"{case_data['name']} step {s_idx}: analyze failed with rc={rc}"

        ops.remove("loadPattern", pat_tag)
        ops.remove("timeSeries", ts_tag)

        # Read OpenSees nodal kinematics
        u_ops = np.array(
            [ops.nodeDisp(ntag, d + 1) for ntag in node_tags for d in range(ndm)],
            dtype=np.float64,
        )
        v_ops = np.array(
            [ops.nodeVel(ntag, d + 1) for ntag in node_tags for d in range(ndm)],
            dtype=np.float64,
        )
        a_ops = np.array(
            [ops.nodeAccel(ntag, d + 1) for ntag in node_tags for d in range(ndm)],
            dtype=np.float64,
        )
        if np.any(null_mask):
            V_null = eigvecs[:, null_mask]
            u_ops = u_ops - V_null @ (V_null.T @ u_ops)
            v_ops = v_ops - V_null @ (V_null.T @ v_ops)
            a_ops = a_ops - V_null @ (V_null.T @ a_ops)

        assert rel_err(u_ops, u_target) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: u mismatch {rel_err(u_ops, u_target):.3e}"
        )
        assert rel_err(v_ops, v_target) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: v mismatch {rel_err(v_ops, v_target):.3e}"
        )
        assert rel_err(a_ops, a_target) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: a mismatch {rel_err(a_ops, a_target):.3e}"
        )

        # Check assembled dynamic tangent K_T_dyn from printA("-ret")
        K_T_ops = np.array(ops.printA("-ret"), dtype=np.float64).reshape(n_dof, n_dof)
        assert rel_err(K_T_ops, K_T_dyn_ref) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: K_T_dyn mismatch {rel_err(K_T_ops, K_T_dyn_ref):.3e}"
        )

        # Check element matrices via eleResponse
        K_eff_ops = np.array(ops.eleResponse(1, "stiff"), dtype=np.float64).reshape(
            n_dof, n_dof
        )
        C_eff_ops = np.array(ops.eleResponse(1, "damping"), dtype=np.float64).reshape(
            n_dof, n_dof
        )
        M_eff_ops = np.array(ops.eleResponse(1, "mass"), dtype=np.float64).reshape(
            n_dof, n_dof
        )
        assert rel_err(K_eff_ops, K_eff_ref) <= 1.0e-11
        assert rel_err(C_eff_ops, C_eff_ref) <= 1.0e-11
        assert rel_err(M_eff_ops, M_eff_ref) <= 1.0e-11

        # Check element force (R_disp) and basicForce (R_int) and energy (P_ext)
        f_disp_ops = np.array(ops.eleResponse(1, "force"), dtype=np.float64)
        f_int_ops = np.array(ops.eleResponse(1, "basicForce"), dtype=np.float64)
        p_ext_ops = float(ops.eleResponse(1, "energy")[0])

        assert rel_err(f_disp_ops, R_disp_ref) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: force mismatch {rel_err(f_disp_ops, R_disp_ref):.3e}"
        )
        assert rel_err(f_int_ops, R_int_ref) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: basicForce mismatch {rel_err(f_int_ops, R_int_ref):.3e}"
        )
        assert rel_err(p_ext_ops, P_ext_ref) <= 1.0e-11, (
            f"{case_data['name']} step {s_idx}: energy mismatch {rel_err(p_ext_ops, P_ext_ref):.3e}"
        )

        u_c, v_c, a_c = u_target, v_target, a_target


@pytest.mark.parametrize("case_data", PARSED_CASES, ids=[c["name"] for c in PARSED_CASES])
def test_group1_opensees_finite_difference_tangent(case_data: dict[str, Any]) -> None:
    """Verify via finite differences in OpenSees that perturbing free nodal DOFs
    matches `ops.printA('-ret')` to <= 1e-7 relative error across all 18 cases.
    """
    ndm = case_data["ndm"]
    n_nodes = case_data["n_nodes"]
    n_dof = n_nodes * ndm
    coords = case_data["coords"]
    dt = float(case_data["steps"][0]["dt"])
    u0 = case_data["steps"][0]["U"].copy()

    def run_single_step_with_load(f_load: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ops.wipe()
        ops.model("basic", "-ndm", ndm, "-ndf", ndm)
        node_tags = list(range(1, n_nodes + 1))
        for idx, ntag in enumerate(node_tags):
            ops.node(ntag, *[float(x) for x in coords[idx]])

        ele_args = build_cpml_element_args(
            ele_tag=1,
            node_tags=node_tags,
            topo=case_data["topo"],
            ndm=ndm,
            E=case_data["E"],
            nu=case_data["nu"],
            rho=case_data["rho"],
            pml=case_data["pml"],
            geom=case_data["geom"],
            lumped=case_data["lumped"],
        )
        ops.element("LadrunoCPML", *ele_args)

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")

        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for a_node, ntag in enumerate(node_tags):
            ops.load(ntag, *[float(f_load[a_node * ndm + d]) for d in range(ndm)])

        rc = ops.analyze(1, dt)
        assert rc == 0
        u_out = np.array(
            [ops.nodeDisp(ntag, d + 1) for ntag in node_tags for d in range(ndm)],
            dtype=np.float64,
        )
        r_int_out = np.array(ops.eleResponse(1, "basicForce"), dtype=np.float64)
        k_t_out = np.array(ops.printA("-ret"), dtype=np.float64).reshape(n_dof, n_dof)
        return u_out, r_int_out, k_t_out

    oracle = ref.ReferenceCPMLElement(
        topo=case_data["topo"],
        coords=coords,
        E=case_data["E"],
        nu=case_data["nu"],
        rho=case_data["rho"],
        pml=case_data["pml"],
        geom=case_data["geom"],
        lumped=case_data["lumped"],
    )
    c1, c2, c3 = 1.0, 2.0 / dt, 4.0 / (dt * dt)
    _, _, _, _, K_T_probe, _ = oracle.evaluate_trial(
        np.zeros(n_dof), np.zeros(n_dof), np.zeros(n_dof), dt, c1, c2, c3
    )
    eigvals, eigvecs = np.linalg.eigh(0.5 * (K_T_probe + K_T_probe.T))
    max_eig = float(np.max(np.abs(eigvals)))
    null_mask = np.abs(eigvals) < 1.0e-10 * max_eig
    if np.any(null_mask):
        V_null = eigvecs[:, null_mask]
        u0 = u0 - V_null @ (V_null.T @ u0)

    f0, _, _, _, _, _ = oracle.evaluate_trial(
        u0, c2 * u0, c3 * u0, dt, c1, c2, c3
    )

    u_base, _, K_T_ops = run_single_step_with_load(f0)
    if np.any(null_mask):
        V_null = eigvecs[:, null_mask]
        u_base = u_base - V_null @ (V_null.T @ u_base)
    assert rel_err(u_base, u0) <= 1.0e-11

    eps = 1.0e-6
    K_fd_ops = np.zeros((n_dof, n_dof), dtype=np.float64)
    for j in range(n_dof):
        df = K_T_ops[:, j] * eps
        u_plus, r_plus, _ = run_single_step_with_load(f0 + df)
        u_minus, r_minus, _ = run_single_step_with_load(f0 - df)
        du_vec = (u_plus - u_minus) / (2.0 * eps)
        if np.any(null_mask):
            V_null = eigvecs[:, null_mask]
            du_vec = du_vec - V_null @ (V_null.T @ du_vec)
            ej_proj = np.zeros(n_dof, dtype=np.float64)
            ej_proj[j] = 1.0
            ej_proj = ej_proj - V_null @ (V_null.T @ ej_proj)
            assert rel_err(du_vec, ej_proj) <= 1.0e-7
        else:
            assert abs(du_vec[j] - 1.0) <= 1.0e-7
        K_fd_ops[:, j] = (r_plus - r_minus) / (2.0 * eps)

    assert rel_err(K_fd_ops, K_T_ops) <= 1.0e-7, (
        f"{case_data['name']}: OpenSees FD tangent error {rel_err(K_fd_ops, K_T_ops):.3e}"
    )


# ============================================================================ #
# Test Group 2 (Anchor 1): 1D P-Wave & S-Wave Waveguide Reflection (R <= 1%)   #
# ============================================================================ #


def ricker_wavelet(t: float, f0: float, t0: float) -> float:
    """Standard Ricker wavelet with peak frequency f0 and delay t0."""
    arg = math.pi * f0 * (t - t0)
    a2 = arg * arg
    return (1.0 - 2.0 * a2) * math.exp(-a2)


def run_1d_waveguide(
    wave_type: str,
    mode: str,
    theta_deg: float = 0.0,
    extended_ref: bool = False,
    use_roller_sp: bool = False,
) -> np.ndarray:
    """Run 2D plane-strain 1D waveguide simulation for P-wave or S-wave.

    Parameters
    ----------
    wave_type : "P" or "S"
    mode : "box" or "plane"
    theta_deg : rotation angle in degrees (0.0 for axis-aligned, 37.0 for oblique test)
    extended_ref : if True, build extended waveguide L_ext = 60.0 m with no PML damping;
                   if False, build L_x = 10.0 m interior + L_pml = 5.0 m C-PML layer.
    use_roller_sp : if True (axis-aligned only), enforce 1D motion via explicit `fix` rollers
                    on both top and bottom edges instead of periodic `equalDOF`.

    Returns
    -------
    u_rec : 1D numpy array of axial (P) or transverse (S) displacement at receiver x_r = 6.0 m.
    """
    E = 2.0e8
    nu = 0.25
    rho = 2000.0
    dx = 0.5
    Hy = 1.0
    Lx_int = 10.0
    Lpml = 5.0
    L_total = 60.0 if extended_ref else (Lx_int + Lpml)
    nx_elem = int(round(L_total / dx))

    dt = 0.0005
    T_sim = 0.16  # 320 steps: covers full P-wave and S-wave round-trip to x=15m and back to x_r=6m
    n_steps = int(round(T_sim / dt))
    f0 = 25.0
    t0 = 0.03

    theta = math.radians(theta_deg)
    ct = math.cos(theta)
    st = math.sin(theta)
    R_rot = np.array([[ct, -st], [st, ct]], dtype=np.float64)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    def ntag(i: int, j: int) -> int:
        return j * (nx_elem + 1) + i + 1

    for j in range(2):
        sy = float(j) * Hy
        for i in range(nx_elem + 1):
            sx = float(i) * dx
            xg = R_rot @ np.array([sx, sy], dtype=np.float64)
            ops.node(ntag(i, j), float(xg[0]), float(xg[1]))

    m_deg = 2.0
    R0 = 1.0e-6
    kappa_max = 1.0
    alpha_max = math.pi * f0

    for i in range(nx_elem):
        n1 = ntag(i, 0)
        n2 = ntag(i + 1, 0)
        n3 = ntag(i + 1, 1)
        n4 = ntag(i, 1)
        etag = i + 1

        if extended_ref:
            ops.element(
                "LadrunoCPML",
                etag,
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
                "-box",
                -100.0,
                200.0,
                -100.0,
                200.0,
                "-Lx",
                Lpml,
                "-Ly",
                Lpml,
                "-m",
                m_deg,
                "-R0",
                R0,
                "-kappaMax",
                kappa_max,
                "-alphaMax",
                alpha_max,
            )
        elif mode == "box":
            assert abs(theta_deg) <= 1.0e-12
            ops.element(
                "LadrunoCPML",
                etag,
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
                "-box",
                -10.0,
                Lx_int,
                -100.0,
                100.0,
                "-Lx",
                Lpml,
                "-Ly",
                Lpml,
                "-m",
                m_deg,
                "-R0",
                R0,
                "-kappaMax",
                kappa_max,
                "-alphaMax",
                alpha_max,
            )
        elif mode == "plane":
            x0_vec = R_rot @ np.array([Lx_int, 0.0], dtype=np.float64)
            n_vec = R_rot @ np.array([1.0, 0.0], dtype=np.float64)
            ops.element(
                "LadrunoCPML",
                etag,
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
                "-plane",
                "-x0",
                float(x0_vec[0]),
                float(x0_vec[1]),
                "-n",
                float(n_vec[0]),
                float(n_vec[1]),
                "-L",
                Lpml,
                "-m",
                m_deg,
                "-R0",
                R0,
                "-kappaMax",
                kappa_max,
                "-alphaMax",
                alpha_max,
            )
        else:
            raise ValueError(f"Unknown mode={mode}")

    # Fixed outer wall at i = nx_elem
    for j in (0, 1):
        ops.fix(ntag(nx_elem, j), 1, 1)

    for i in range(nx_elem):
        nb = ntag(i, 0)
        nt = ntag(i, 1)
        if use_roller_sp:
            assert abs(theta_deg) <= 1.0e-12
            if wave_type == "P":
                ops.fix(nb, 0, 1)
                ops.fix(nt, 0, 1)
            else:
                ops.fix(nb, 1, 0)
                ops.fix(nt, 1, 0)
        else:
            ops.equalDOF(nb, nt, 1, 2)

    ts_vals = [ricker_wavelet((step + 1) * dt, f0, t0) for step in range(n_steps + 5)]
    ts_times = [(step + 1) * dt for step in range(n_steps + 5)]
    ops.timeSeries("Path", 1, "-time", *ts_times, "-values", *ts_vals)
    ops.pattern("Plain", 1, 1)

    f_amp = 1.0e5 * (0.5 * Hy)
    if wave_type == "P":
        fx = f_amp * ct
        fy = f_amp * st
    else:
        fx = -f_amp * st
        fy = f_amp * ct

    ops.load(ntag(0, 0), fx, fy)
    ops.load(ntag(0, 1), fx, fy)

    ops.constraints("Plain" if use_roller_sp else "Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    rec_i = int(round(6.0 / dx))
    rec_node = ntag(rec_i, 0)
    u_rec = np.zeros(n_steps, dtype=np.float64)

    for step in range(n_steps):
        rc = ops.analyze(1, dt)
        assert rc == 0, f"Waveguide analyze failed at step {step}"
        ux = ops.nodeDisp(rec_node, 1)
        uy = ops.nodeDisp(rec_node, 2)
        if wave_type == "P":
            u_rec[step] = ux * ct + uy * st
        else:
            u_rec[step] = -ux * st + uy * ct

    return u_rec


def test_group2_anchor1_waveguide_reflection_and_invariance() -> None:
    """Anchor 1: Verify 1D P-wave and S-wave reflection ratio R <= 1% in `-box` and `-plane`
    modes, bit-identical `-box` vs `-plane` trajectories, and exact 37-degree rotational
    frame invariance in `-plane` mode.
    """
    # (a) 1D P-wave normal incidence (-box mode vs extended reference)
    u_p_ref = run_1d_waveguide(wave_type="P", mode="box", extended_ref=True, use_roller_sp=True)
    u_p_box = run_1d_waveguide(wave_type="P", mode="box", extended_ref=False, use_roller_sp=True)
    R_p_box = float(np.max(np.abs(u_p_box - u_p_ref)) / np.max(np.abs(u_p_ref)))
    assert R_p_box <= 0.01, f"1D P-wave (-box) reflection ratio {R_p_box:.4f} > 0.01"

    # (b) 1D S-wave normal incidence (-box mode vs extended reference)
    u_s_ref = run_1d_waveguide(wave_type="S", mode="box", extended_ref=True, use_roller_sp=True)
    u_s_box = run_1d_waveguide(wave_type="S", mode="box", extended_ref=False, use_roller_sp=True)
    R_s_box = float(np.max(np.abs(u_s_box - u_s_ref)) / np.max(np.abs(u_s_ref)))
    assert R_s_box <= 0.01, f"1D S-wave (-box) reflection ratio {R_s_box:.4f} > 0.01"

    # (c) 1D P-wave & S-wave with `-plane` mode (-plane -x0 10.0 0.0 -n 1.0 0.0 -L 5.0)
    u_p_plane = run_1d_waveguide(wave_type="P", mode="plane", extended_ref=False)
    u_s_plane = run_1d_waveguide(wave_type="S", mode="plane", extended_ref=False)
    u_p_box_eq = run_1d_waveguide(wave_type="P", mode="box", extended_ref=False)
    u_s_box_eq = run_1d_waveguide(wave_type="S", mode="box", extended_ref=False)

    assert np.allclose(u_p_plane, u_p_box_eq, rtol=1.0e-14, atol=1.0e-18), (
        "P-wave -plane vs -box trajectory mismatch"
    )
    assert np.allclose(u_s_plane, u_s_box_eq, rtol=1.0e-14, atol=1.0e-18), (
        "S-wave -plane vs -box trajectory mismatch"
    )
    assert rel_err(u_p_plane, u_p_box) <= 1.0e-12
    assert rel_err(u_s_plane, u_s_box) <= 1.0e-12

    R_p_plane = float(np.max(np.abs(u_p_plane - u_p_ref)) / np.max(np.abs(u_p_ref)))
    R_s_plane = float(np.max(np.abs(u_s_plane - u_s_ref)) / np.max(np.abs(u_s_ref)))
    assert R_p_plane <= 0.01
    assert R_s_plane <= 0.01

    # (d) Rotated oblique waveguide (theta = 37 deg) with `-plane` mode
    u_p_rot37 = run_1d_waveguide(
        wave_type="P", mode="plane", theta_deg=37.0, extended_ref=False
    )
    u_s_rot37 = run_1d_waveguide(
        wave_type="S", mode="plane", theta_deg=37.0, extended_ref=False
    )
    assert rel_err(u_p_rot37, u_p_plane) <= 1.0e-11, (
        f"Rotated 37-deg P-wave frame invariance error: {rel_err(u_p_rot37, u_p_plane):.3e}"
    )
    assert rel_err(u_s_rot37, u_s_plane) <= 1.0e-11, (
        f"Rotated 37-deg S-wave frame invariance error: {rel_err(u_s_rot37, u_s_plane):.3e}"
    )


# ============================================================================ #
# Test Group 3 (Anchor 2): MMS & Unstructured Patch Completeness (9 Topologies)#
# ============================================================================ #


def build_unstructured_patch(topo: str) -> tuple[np.ndarray, list[list[int]]]:
    """Build a skewed multi-element unstructured patch in [2.0, 3.6]^ndm (inside an
    active C-PML corner) for each of the 9 topologies.
    Returns (node_coords, elements_connectivity_0based).
    """
    if topo == "T3":
        coords = np.array(
            [
                [2.00, 2.00],
                [3.35, 2.08],
                [3.45, 3.42],
                [2.06, 3.28],
                [2.72, 2.68],
            ],
            dtype=np.float64,
        )
        conn = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
        return coords, conn

    if topo == "Q4":
        coords = np.array(
            [
                [2.00, 2.00],
                [2.75, 2.06],
                [2.82, 3.32],
                [2.04, 3.25],
                [3.52, 2.12],
                [3.48, 3.40],
            ],
            dtype=np.float64,
        )
        conn = [[0, 1, 2, 3], [1, 4, 5, 2]]
        return coords, conn

    if topo in ("T6", "BT6"):
        v_corners = np.array(
            [
                [2.00, 2.00],
                [3.38, 2.08],
                [3.44, 3.42],
                [2.06, 3.30],
            ],
            dtype=np.float64,
        )
        tri_verts = [(0, 1, 2), (0, 2, 3)]
        coords_list = [v_corners[i] for i in range(4)]
        edge_map: dict[tuple[int, int], int] = {}

        def get_mid(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in edge_map:
                edge_map[key] = len(coords_list)
                coords_list.append(0.5 * (v_corners[a] + v_corners[b]))
            return edge_map[key]

        conn = []
        for a, b, c in tri_verts:
            conn.append([a, b, c, get_mid(a, b), get_mid(b, c), get_mid(c, a)])
        return np.array(coords_list, dtype=np.float64), conn

    if topo == "TET4":
        # 5-tetrahedron decomposition of a skewed hexahedron (8 vertices, 5 skewed TET4s)
        origin = np.array([2.05, 2.05, 2.05], dtype=np.float64)
        e1 = np.array([1.25, 0.08, 0.04], dtype=np.float64)
        e2 = np.array([0.06, 1.22, 0.05], dtype=np.float64)
        e3 = np.array([0.05, 0.07, 1.28], dtype=np.float64)
        coords = np.array(
            [
                origin,
                origin + e1,
                origin + e1 + e2,
                origin + e2,
                origin + e3,
                origin + e1 + e3,
                origin + e1 + e2 + e3,
                origin + e2 + e3,
            ],
            dtype=np.float64,
        )
        conn = [
            [0, 1, 3, 4],
            [1, 2, 3, 6],
            [4, 6, 5, 1],
            [4, 7, 6, 3],
            [1, 3, 4, 6],
        ]
        return coords, conn

    if topo == "H8":
        # 2 skewed parallelepiped H8 bricks sharing face (1, 2, 6, 5)
        origin = np.array([2.00, 2.00, 2.00], dtype=np.float64)
        a1 = np.array([0.72, 0.05, 0.03], dtype=np.float64)
        a1_ext = np.array([1.48, 0.09, 0.05], dtype=np.float64)
        a2 = np.array([0.06, 1.22, 0.04], dtype=np.float64)
        a3 = np.array([0.04, 0.06, 1.24], dtype=np.float64)
        coords = np.array(
            [
                origin,
                origin + a1,
                origin + a1 + a2,
                origin + a2,
                origin + a3,
                origin + a1 + a3,
                origin + a1 + a2 + a3,
                origin + a2 + a3,
                origin + a1_ext,
                origin + a1_ext + a2,
                origin + a1_ext + a3,
                origin + a1_ext + a2 + a3,
            ],
            dtype=np.float64,
        )
        conn = [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 8, 9, 2, 5, 10, 11, 6],
        ]
        return coords, conn

    if topo in ("TET10", "BTET10"):
        # 2 skewed 10-node tetrahedra sharing face (0, 1, 2), oriented with det(J) > 0
        # (in TET10/BTET10, node 3 is the reference origin z4=1-z1-z2-z3)
        v_base = np.array(
            [
                [2.10, 2.05, 2.15],
                [3.35, 2.12, 2.10],
                [2.20, 3.30, 2.18],
                [2.25, 2.20, 3.35],
                [2.30, 2.25, 1.92],
            ],
            dtype=np.float64,
        )
        tet4_conn = [[0, 2, 1, 3], [0, 1, 2, 4]]
        coords_list = [v_base[i] for i in range(len(v_base))]
        edge_map_3d: dict[tuple[int, int], int] = {}

        def get_mid3(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in edge_map_3d:
                edge_map_3d[key] = len(coords_list)
                coords_list.append(0.5 * (v_base[a] + v_base[b]))
            return edge_map_3d[key]

        conn = []
        for n0, n1, n2, n3 in tet4_conn:
            conn.append(
                [
                    n0,
                    n1,
                    n2,
                    n3,
                    get_mid3(n0, n1),
                    get_mid3(n1, n2),
                    get_mid3(n2, n0),
                    get_mid3(n0, n3),
                    get_mid3(n2, n3),
                    get_mid3(n1, n3),
                ]
            )
        return np.array(coords_list, dtype=np.float64), conn

    if topo == "H20":
        v_base, h8_conn = build_unstructured_patch("H8")
        coords_list = [v_base[i] for i in range(len(v_base))]
        edge_map_h20: dict[tuple[int, int], int] = {}

        def get_mid_h20(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in edge_map_h20:
                edge_map_h20[key] = len(coords_list)
                coords_list.append(0.5 * (v_base[a] + v_base[b]))
            return edge_map_h20[key]

        edges_local = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ]
        conn = []
        for c8 in h8_conn:
            mids = [get_mid_h20(c8[a], c8[b]) for a, b in edges_local]
            conn.append(list(c8) + mids)
        return np.array(coords_list, dtype=np.float64), conn

    raise ValueError(f"Unsupported topology {topo}")


def eval_mms_spatial_field(
    x: np.ndarray, degree: int
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate exact spatial polynomial field U_spatial(x) and its gradient dU_spatial/dx
    for degree=1 (linear) or degree=2 (quadratic) in 2D or 3D.
    """
    ndm = len(x)
    if ndm == 2:
        X, Y = float(x[0]), float(x[1])
        u1 = 1.2e-4 + 2.3e-4 * X - 1.7e-4 * Y
        u2 = -0.8e-4 + 1.1e-4 * X + 1.9e-4 * Y
        grad = np.array(
            [
                [2.3e-4, -1.7e-4],
                [1.1e-4, 1.9e-4],
            ],
            dtype=np.float64,
        )
        if degree == 2:
            u1 += 0.6e-4 * X * X - 0.5e-4 * X * Y + 0.4e-4 * Y * Y
            u2 += -0.3e-4 * X * X + 0.7e-4 * X * Y - 0.5e-4 * Y * Y
            grad[0, 0] += 1.2e-4 * X - 0.5e-4 * Y
            grad[0, 1] += -0.5e-4 * X + 0.8e-4 * Y
            grad[1, 0] += -0.6e-4 * X + 0.7e-4 * Y
            grad[1, 1] += 0.7e-4 * X - 1.0e-4 * Y
        return np.array([u1, u2], dtype=np.float64), grad
    else:
        X, Y, Z = float(x[0]), float(x[1]), float(x[2])
        u1 = 1.1e-4 + 2.1e-4 * X - 1.3e-4 * Y + 0.9e-4 * Z
        u2 = -0.7e-4 + 0.8e-4 * X + 1.8e-4 * Y - 1.2e-4 * Z
        u3 = 0.5e-4 - 1.4e-4 * X + 0.6e-4 * Y + 1.5e-4 * Z
        grad = np.array(
            [
                [2.1e-4, -1.3e-4, 0.9e-4],
                [0.8e-4, 1.8e-4, -1.2e-4],
                [-1.4e-4, 0.6e-4, 1.5e-4],
            ],
            dtype=np.float64,
        )
        if degree == 2:
            u1 += 0.4e-4 * X * X - 0.3e-4 * X * Y + 0.25e-4 * Y * Z - 0.2e-4 * Z * Z
            u2 += -0.25e-4 * Y * Y + 0.35e-4 * X * Z + 0.2e-4 * X * Y
            u3 += 0.3e-4 * Z * Z - 0.4e-4 * Y * Z + 0.15e-4 * X * X
            grad[0, 0] += 0.8e-4 * X - 0.3e-4 * Y
            grad[0, 1] += -0.3e-4 * X + 0.25e-4 * Z
            grad[0, 2] += 0.25e-4 * Y - 0.4e-4 * Z
            grad[1, 0] += 0.35e-4 * Z + 0.2e-4 * Y
            grad[1, 1] += -0.5e-4 * Y + 0.2e-4 * X
            grad[1, 2] += 0.35e-4 * X
            grad[2, 0] += 0.3e-4 * X
            grad[2, 1] += -0.4e-4 * Z
            grad[2, 2] += 0.6e-4 * Z - 0.4e-4 * Y
        return np.array([u1, u2, u3], dtype=np.float64), grad


def compute_mms_nodal_dofs(
    topo: str, coords: np.ndarray, conn: list[list[int]], degree: int
) -> np.ndarray:
    """Compute nodal DOFs for the unstructured patch so that the finite element field
    reproduces `eval_mms_spatial_field(x, degree)` identically at every point in the patch.
    """
    n_nodes, ndm = coords.shape
    u_nodal = np.zeros((n_nodes, ndm), dtype=np.float64)
    for a in range(n_nodes):
        u_nodal[a], _ = eval_mms_spatial_field(coords[a], degree)

    if topo == "BT6" and degree == 2:
        for elem_nodes in conn:
            v0, v1, v2, m3, m4, m5 = elem_nodes
            for m_idx, (va, vb) in ((m3, (v0, v1)), (m4, (v1, v2)), (m5, (v2, v0))):
                u_mid, _ = eval_mms_spatial_field(coords[m_idx], degree)
                u_a, _ = eval_mms_spatial_field(coords[va], degree)
                u_b, _ = eval_mms_spatial_field(coords[vb], degree)
                u_nodal[m_idx] = 2.0 * u_mid - 0.5 * (u_a + u_b)
    elif topo == "BTET10" and degree == 2:
        for elem_nodes in conn:
            v0, v1, v2, v3, m4, m5, m6, m7, m8, m9 = elem_nodes
            edge_pairs = [
                (m4, (v0, v1)),
                (m5, (v1, v2)),
                (m6, (v2, v0)),
                (m7, (v0, v3)),
                (m8, (v2, v3)),
                (m9, (v1, v3)),
            ]
            for m_idx, (va, vb) in edge_pairs:
                u_mid, _ = eval_mms_spatial_field(coords[m_idx], degree)
                u_a, _ = eval_mms_spatial_field(coords[va], degree)
                u_b, _ = eval_mms_spatial_field(coords[vb], degree)
                u_nodal[m_idx] = 2.0 * u_mid - 0.5 * (u_a + u_b)

    return u_nodal.reshape(n_nodes * ndm)


@pytest.mark.parametrize(
    "topo,degree",
    [
        ("T3", 1),
        ("Q4", 1),
        ("T6", 2),
        ("BT6", 2),
        ("TET4", 1),
        ("H8", 1),
        ("TET10", 2),
        ("BTET10", 2),
        ("H20", 2),
    ],
)
@pytest.mark.parametrize("geom_mode", ["box_corner", "plane_oblique"])
def test_group3_anchor2_mms_patch_completeness_all_9_topologies(
    topo: str, degree: int, geom_mode: str
) -> None:
    """Anchor 2 (Part 1): Across all 9 topologies on unstructured/skewed multi-element
    patches inside active C-PML regions (`-box` corner and `-plane` oblique), prescribe a
    smooth manufactured displacement field u_MMS(x, t) (linear for T3/Q4/TET4/H8, quadratic
    for T6/BT6/TET10/BTET10/H20) and verify that:
      1. The FE interpolation reproduces exact u_MMS(x_g, t) and grad u_MMS(x_g, t) at every
         Gauss point across the skewed multi-element patch to <= 1e-13.
      2. When driven by the exact analytical/reference quadrature integral of the stretched
         CFS-PML PDE + RC filter ODEs across multiple time steps, OpenSees reproduces the
         exact manufactured nodal displacement, velocity, acceleration, and assembled
         internal force across the unstructured patch to <= 1e-11!
    """
    coords, conn = build_unstructured_patch(topo)
    n_nodes, ndm = coords.shape
    n_dof = n_nodes * ndm

    E = 2.4e8
    nu = 0.27
    rho = 2150.0
    pml = ref.PMLParams(m_deg=2.0, R0=1.0e-4, kappa_max=4.0, alpha_max=math.pi * 15.0)
    lumped_flag = 1 if topo in ("TET4", "H20") else 0

    if geom_mode == "box_corner":
        if ndm == 2:
            geom = ref.GeomSpec("box", (-1.0, 1.5, -1.0, 1.5, 2.5, 2.5))
        else:
            geom = ref.GeomSpec("box", (-1.0, 1.5, -1.0, 1.5, -1.0, 1.5, 2.5, 2.5, 2.5))
    else:
        if ndm == 2:
            c37 = math.cos(math.radians(37.0))
            s37 = math.sin(math.radians(37.0))
            geom = ref.GeomSpec("plane", (1.2, 1.1, c37, s37, 3.0))
        else:
            geom = ref.GeomSpec("plane", (1.0, 1.0, 1.0, 0.6, -0.5, 0.6245, 3.5))

    u_spatial_dofs = compute_mms_nodal_dofs(topo, coords, conn, degree)

    ref_elems: list[ref.ReferenceCPMLElement] = []
    for elem_nodes in conn:
        el_coords = coords[elem_nodes]
        el_u = u_spatial_dofs.reshape(n_nodes, ndm)[elem_nodes]
        gp_kin = ref.eval_gp_kinematics(topo, el_coords)
        for N, dN_dx, _w_vol, x_g in gp_kin:
            u_exact, grad_exact = eval_mms_spatial_field(x_g, degree)
            u_fe = N @ el_u
            grad_fe = el_u.T @ dN_dx
            assert rel_err(u_fe, u_exact) <= 1.0e-13
            assert rel_err(grad_fe, grad_exact) <= 1.0e-13

        ref_elems.append(
            ref.ReferenceCPMLElement(
                topo=topo,
                coords=el_coords,
                E=E,
                nu=nu,
                rho=rho,
                pml=pml,
                geom=geom,
                lumped=lumped_flag,
            )
        )

    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)
    for a in range(n_nodes):
        ops.node(a + 1, *[float(x) for x in coords[a]])

    for e_idx, elem_nodes in enumerate(conn):
        ntags = [n + 1 for n in elem_nodes]
        ele_args = build_cpml_element_args(
            ele_tag=e_idx + 1,
            node_tags=ntags,
            topo=topo,
            ndm=ndm,
            E=E,
            nu=nu,
            rho=rho,
            pml=pml,
            geom=geom,
            lumped=lumped_flag,
        )
        ops.element("LadrunoCPML", *ele_args)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    dt = 0.002
    n_steps = 4
    c1, c2, c3 = 1.0, 2.0 / dt, 4.0 / (dt * dt)

    u_c = np.zeros(n_dof, dtype=np.float64)
    v_c = np.zeros(n_dof, dtype=np.float64)
    a_c = np.zeros(n_dof, dtype=np.float64)

    for step in range(1, n_steps + 1):
        t_n = step * dt
        g_t = (t_n * 100.0) ** 2 * math.exp(-10.0 * t_n)
        u_mms = g_t * u_spatial_dofs
        a_mms = c3 * (u_mms - u_c) - (4.0 / dt) * v_c - a_c
        v_mms = v_c + 0.5 * dt * (a_c + a_mms)

        F_mms_global = np.zeros(n_dof, dtype=np.float64)
        K_T_mms_global = np.zeros((n_dof, n_dof), dtype=np.float64)

        for elem_nodes, ref_el in zip(conn, ref_elems):
            dof_idx = [n * ndm + d for n in elem_nodes for d in range(ndm)]
            u_el = u_mms[dof_idx]
            v_el = v_mms[dof_idx]
            a_el = a_mms[dof_idx]
            R_el, _, _, _, K_T_el, _ = ref_el.evaluate_trial(
                u_el, v_el, a_el, dt, c1, c2, c3
            )
            ref_el.commit()
            F_mms_global[dof_idx] += R_el
            K_T_mms_global[np.ix_(dof_idx, dof_idx)] += K_T_el

        ts_tag = 500 + step
        pat_tag = 600 + step
        ops.timeSeries("Constant", ts_tag)
        ops.pattern("Plain", pat_tag, ts_tag)
        for a in range(n_nodes):
            ops.load(a + 1, *[float(F_mms_global[a * ndm + d]) for d in range(ndm)])

        rc = ops.analyze(1, dt)
        assert rc == 0

        ops.remove("loadPattern", pat_tag)
        ops.remove("timeSeries", ts_tag)

        u_ops = np.array(
            [ops.nodeDisp(a + 1, d + 1) for a in range(n_nodes) for d in range(ndm)],
            dtype=np.float64,
        )
        v_ops = np.array(
            [ops.nodeVel(a + 1, d + 1) for a in range(n_nodes) for d in range(ndm)],
            dtype=np.float64,
        )
        a_ops = np.array(
            [ops.nodeAccel(a + 1, d + 1) for a in range(n_nodes) for d in range(ndm)],
            dtype=np.float64,
        )

        F_int_ops = np.zeros(n_dof, dtype=np.float64)
        for e_idx, elem_nodes in enumerate(conn):
            dof_idx = [n * ndm + d for n in elem_nodes for d in range(ndm)]
            f_el = np.array(ops.eleResponse(e_idx + 1, "basicForce"), dtype=np.float64)
            F_int_ops[dof_idx] += f_el

        K_T_ops = np.array(ops.printA("-ret"), dtype=np.float64).reshape(n_dof, n_dof)

        assert rel_err(u_ops, u_mms) <= 1.0e-11
        assert rel_err(v_ops, v_mms) <= 1.0e-11
        assert rel_err(a_ops, a_mms) <= 1.0e-11
        assert rel_err(F_int_ops, F_mms_global) <= 1.0e-11
        assert rel_err(K_T_ops, K_T_mms_global) <= 1.0e-11

        u_c, v_c, a_c = u_mms, v_mms, a_mms


def run_2d_mms_mesh_refinement(topo: str, N: int) -> float:
    """Run 2D non-polynomial trigonometric MMS spatial convergence check on a perturbed
    unstructured mesh of size N x N inside an active C-PML corner `[2.0, 3.0]^2`.
    """
    x_min_dom, x_max_dom = 2.0, 3.0
    h = (x_max_dom - x_min_dom) / float(N)

    E = 2.0e8
    nu = 0.25
    rho = 2000.0
    lam, mu, Vp, C_tensor = ref.isotropic_elasticity(E, nu, rho, 2)
    pml = ref.PMLParams(m_deg=2.0, R0=1.0e-4, kappa_max=3.0, alpha_max=math.pi * 15.0)
    geom = ref.GeomSpec("box", (0.0, 1.5, 0.0, 1.5, 2.0, 2.0))

    n_grid = N + 1
    coords = np.zeros((n_grid * n_grid, 2), dtype=np.float64)
    for j in range(n_grid):
        for i in range(n_grid):
            nid = j * n_grid + i
            x = x_min_dom + i * h
            y = x_min_dom + j * h
            if 0 < i < N and 0 < j < N:
                x += 0.08 * h * math.sin(2.0 * math.pi * i / N) * math.sin(
                    2.0 * math.pi * j / N
                )
                y -= 0.08 * h * math.sin(2.0 * math.pi * i / N) * math.sin(
                    2.0 * math.pi * j / N
                )
            coords[nid] = [x, y]

    conn: list[list[int]] = []
    for j in range(N):
        for i in range(N):
            n00 = j * n_grid + i
            n10 = j * n_grid + (i + 1)
            n11 = (j + 1) * n_grid + (i + 1)
            n01 = (j + 1) * n_grid + i
            if topo == "Q4":
                conn.append([n00, n10, n11, n01])
            elif topo == "T3":
                if (i + j) % 2 == 0:
                    conn.append([n00, n10, n11])
                    conn.append([n00, n11, n01])
                else:
                    conn.append([n00, n10, n01])
                    conn.append([n10, n11, n01])
            else:
                raise ValueError(f"Unsupported topo={topo}")

    def exact_spatial(x_pt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sx = math.sin(math.pi * (x_pt[0] - 2.0))
        cx = math.cos(math.pi * (x_pt[0] - 2.0))
        sy = math.sin(math.pi * (x_pt[1] - 2.0))
        cy = math.cos(math.pi * (x_pt[1] - 2.0))
        phi = sx * sy
        dphi_dx = math.pi * cx * sy
        dphi_dy = math.pi * sx * cy
        u_val = np.array([1.0e-3 * phi, 0.7e-3 * phi], dtype=np.float64)
        grad_val = np.array(
            [
                [1.0e-3 * dphi_dx, 1.0e-3 * dphi_dy],
                [0.7e-3 * dphi_dx, 0.7e-3 * dphi_dy],
            ],
            dtype=np.float64,
        )
        return u_val, grad_val

    if topo == "Q4":
        pts_1d, wts_1d = np.polynomial.legendre.leggauss(4)
        hq_rule = [
            ((float(xi), float(eta)), float(wx * wy))
            for eta, wy in zip(pts_1d, wts_1d)
            for xi, wx in zip(pts_1d, wts_1d)
        ]
        shape_fn = ref.shape_q4
    else:
        a1, b1, w1 = 0.445948490915965, 0.108103018168070, 0.111690794839005
        a2, b2, w2 = 0.091576213509771, 0.816847572980459, 0.054975871827661
        hq_rule = [
            ((a1, a1), w1),
            ((b1, a1), w1),
            ((a1, b1), w1),
            ((a2, a2), w2),
            ((b2, a2), w2),
            ((a2, b2), w2),
        ]
        shape_fn = ref.shape_t3

    dt = 0.002
    c1, c2, c3 = 1.0, 2.0 / dt, 4.0 / (dt * dt)
    n_nodes_total = len(coords)
    F_mms = np.zeros((n_nodes_total, 2), dtype=np.float64)

    for elem_nodes in conn:
        el_coords = coords[elem_nodes]
        for xi_pt, w_q in hq_rule:
            N_vec, dN_dxi = shape_fn(xi_pt)
            J = dN_dxi.T @ el_coords
            detJ = float(np.linalg.det(J))
            dN_dx = np.linalg.solve(J, dN_dxi.T).T
            w_vol = w_q * detJ
            x_g = N_vec @ el_coords

            u_g, grad_g = exact_spatial(x_g)
            v_g = c2 * u_g
            a_g = c3 * u_g

            d_pml, kappa_pml, alpha_pml = ref.evaluate_pml_profile_at_gp(x_g, Vp, pml, geom)

            D_F = np.zeros((2, 2), dtype=np.float64)
            for j_dir in range(2):
                for l_dir in range(2):
                    casc = ref.make_stiffness_cascade(2, j_dir, l_dir, d_pml, kappa_pml, alpha_pml)
                    _, D_F[j_dir, l_dir] = casc.evaluate_trial(1.0, dt)

            _, dg0 = ref.make_g_filter(d_pml[0], kappa_pml[0], alpha_pml[0]).evaluate_trial(1.0, dt)
            _, dg1 = ref.make_g_filter(d_pml[1], kappa_pml[1], alpha_pml[1]).evaluate_trial(1.0, dt)
            _, d10 = ref.FilterCascade(
                [
                    ref.make_g_filter(d_pml[0], kappa_pml[0], alpha_pml[0]),
                    ref.make_g_filter(d_pml[1], kappa_pml[1], alpha_pml[1]),
                ]
            ).evaluate_trial(1.0, dt)

            k0, k1 = kappa_pml[0], kappa_pml[1]
            D_acc = k0 * k1
            D_vel = k1 * dg0 + k0 * dg1
            D_disp = d10

            for i_comp in range(2):
                for j_dir in range(2):
                    sig_ij = 0.0
                    for k_comp in range(2):
                        for l_dir in range(2):
                            sig_ij += (
                                C_tensor[i_comp, j_dir, k_comp, l_dir]
                                * D_F[j_dir, l_dir]
                                * grad_g[k_comp, l_dir]
                            )
                    F_mms[elem_nodes, i_comp] += w_vol * sig_ij * dN_dx[:, j_dir]

            f_inert_g = rho * (D_acc * a_g + D_vel * v_g + D_disp * u_g)
            F_mms[elem_nodes, :] += w_vol * np.outer(N_vec, f_inert_g)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for nid in range(n_nodes_total):
        ops.node(nid + 1, float(coords[nid, 0]), float(coords[nid, 1]))

    for e_idx, elem_nodes in enumerate(conn):
        ntags = [n + 1 for n in elem_nodes]
        ele_args = build_cpml_element_args(
            ele_tag=e_idx + 1,
            node_tags=ntags,
            topo=topo,
            ndm=2,
            E=E,
            nu=nu,
            rho=rho,
            pml=pml,
            geom=geom,
            lumped=0,
        )
        ops.element("LadrunoCPML", *ele_args)

    for j in range(n_grid):
        for i in range(n_grid):
            if i == 0 or i == N or j == 0 or j == N:
                ops.fix(j * n_grid + i + 1, 1, 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(1, N):
        for i in range(1, N):
            nid = j * n_grid + i
            ops.load(nid + 1, float(F_mms[nid, 0]), float(F_mms[nid, 1]))

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    rc = ops.analyze(1, dt)
    assert rc == 0

    u_ops = np.array(
        [[ops.nodeDisp(nid + 1, 1), ops.nodeDisp(nid + 1, 2)] for nid in range(n_nodes_total)],
        dtype=np.float64,
    )

    err_sq = 0.0
    norm_sq = 0.0
    for elem_nodes in conn:
        el_coords = coords[elem_nodes]
        el_u = u_ops[elem_nodes]
        for xi_pt, w_q in hq_rule:
            N_vec, dN_dxi = shape_fn(xi_pt)
            J = dN_dxi.T @ el_coords
            w_vol = w_q * float(np.linalg.det(J))
            x_g = N_vec @ el_coords
            u_exact, _ = exact_spatial(x_g)
            u_fe = N_vec @ el_u
            err_sq += w_vol * float(np.sum((u_fe - u_exact) ** 2))
            norm_sq += w_vol * float(np.sum(u_exact ** 2))

    return math.sqrt(err_sq / norm_sq)


@pytest.mark.parametrize("topo", ["T3", "Q4"])
def test_group3_anchor2_mms_spatial_convergence_rate(topo: str) -> None:
    """Anchor 2 (Part 2): 2D spatial mesh refinement study (N in {4, 8, 16}) on
    unstructured `T3` and `Q4` meshes inside an active C-PML corner driven by a smooth
    non-polynomial trigonometric manufactured field, verifying L2 displacement error
    converges at rate p >= 1.85.
    """
    errs = [run_2d_mms_mesh_refinement(topo, N) for N in (4, 8, 16)]
    rate_4_8 = math.log2(errs[0] / errs[1])
    rate_8_16 = math.log2(errs[1] / errs[2])
    assert rate_4_8 >= 1.85, f"{topo} MMS convergence rate (4->8) = {rate_4_8:.3f} < 1.85"
    assert rate_8_16 >= 1.85, f"{topo} MMS convergence rate (8->16) = {rate_8_16:.3f} < 1.85"
