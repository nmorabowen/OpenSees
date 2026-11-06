"""Phase P2-UP Verification Suite (Worker 2.A-UP) for Coupled Biot u-p `LadrunoCPML -up`
(ELE_TAG_LadrunoCPML = 33023).

Covers:
1. `test_up_gate0_interior_equivalence_vs_ladrunoup`:
   - Gate 0a (xi = 0 Exact Equivalence vs. `LadrunoUP`): across all 5 shared topologies
     (`T3, Q4, H8` with `pOrder="equal"`, and `BT6, BTET10` with `pOrder="linear"` Taylor-Hood)
     and both `-formulation std` and `-formulation bbar` (on `Q4` and `H8`), compare an
     interior `LadrunoCPML -up` element against `LadrunoUP` (`-dynSeepage off`) over a
     3-step `TransientAnalysis` (`Newmark(0.5, 0.25)`) under non-trivial nodal loads on
     both `u` and `p` DOFs plus body and fluidBody forces (<= 1e-13 relative error).
2. `test_up_gate0_oracle_14_cases_and_fd_tangent`:
   - Gate 0b (All 14 Multi-Step Reference Oracle Cases Through `opensees.pyd` + 4-Block FD Tangent):
     parse `tests/ladruno_cpml_up_cases.txt` (all 14 cases across `T3, Q4, T6, BT6, TET4, H8,
     TET10, BTET10, H20`), verify against `ReferenceCPMLUPElement` and `opensees.pyd` across
     all 4 blocks (`uu, up, pu, pp`) to <= 1e-10 relative error, verify `-mat $matTag` and
     `-permH ... -gammaW ...` equivalence, and verify in-engine 4-block FD tangent & skew-duality.
3. `test_up_anchor1_saturated_soil_column_p1_and_s_wave_reflection`:
   - Anchor UP-1: 1D saturated poroelastic waveguide reflection R_u <= 1%, R_p <= 1% for
     undrained/high-mobility Fast P-wave (P1) and poroelastic Shear wave (S-wave) in `-box`
     and `-plane` modes (vs >80% reflection on truncated column), plus 37-degree rotated
     `-plane` poroelastic waveguide frame invariance to <= 1e-11.
4. `test_up_anchor2_2d_3d_poroelastic_radiation_and_lamb_halfspace`:
   - Anchor UP-2: (a) 2D cylindrical poroelastic wave radiation (`LadrunoUP` Q4 interior +
     `LadrunoCPML -up` face & 2D corner, L2 error <= 3.5% on u and p, >15x reduction vs
     unabsorbed boundary); (b) 3D spherical poroelastic wave radiation (`LadrunoUP` H8
     interior + `LadrunoCPML -up` 3D face, edge & corner, error <= 4.5% on u and p, >10x
     reduction); (c) 2D saturated half-space surface point load (Rayleigh wave + pore-pressure
     diffusion, error <= 3.5% on u and p).
5. `test_up_anchor3_gravity_seepage_staging_rayleigh_stability_and_checkpoint`:
   - Anchor UP-3: (a) Hydrostatic gravity + steady seepage `-stage 0 -> 1` zero-jump
     verification (with direct `-body`/`-fluidBody` and `SelfWeight`, `stage` and
     `updateMaterialStage`); (b) Solid-only Rayleigh damping isolation (`alphaM, betaK,
     betaK0, betaKc`); (c) 2500-step late-time poroelastic Lyapunov stability in 2D and 3D
     corners; (d) `FileDatastore` `sendSelf`/`recvSelf` checkpoint round-trip across equal-order
     (`Q4, H8`) and Taylor-Hood (`BT6, BTET10`) topologies (live and fresh-domain restore);
     (e) Fail-closed `-up` parser refusal battery.
"""
from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import opensees as ops
import pytest

import ladruno_cpml_reference as ref

pytestmark = [pytest.mark.zone_a]

UP_CASES_TXT_PATH = os.path.join(os.path.dirname(__file__), "ladruno_cpml_up_cases.txt")


# ============================================================================ #
# Shared Helpers                                                               #
# ============================================================================ #


def rel_err(a: np.ndarray | float, b: np.ndarray | float) -> float:
    """Compute relative norm error ||a - b|| / max(1e-30, ||b||)."""
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    diff = float(np.linalg.norm(a_arr - b_arr))
    denom = max(1.0e-30, float(np.linalg.norm(b_arr)))
    return diff / denom


def rel_err_scaled(a: np.ndarray | float, b: np.ndarray | float) -> float:
    """Compute relative norm error ||a - b|| / max(1.0, ||b||)."""
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    diff = float(np.linalg.norm(a_arr - b_arr))
    denom = max(1.0, float(np.linalg.norm(b_arr)))
    return diff / denom


def ricker_wavelet(t: float, f0: float, t0: float) -> float:
    """Standard Ricker wavelet with peak frequency f0 and delay t0."""
    arg = math.pi * f0 * (t - t0)
    a2 = arg * arg
    return (1.0 - 2.0 * a2) * math.exp(-a2)


def ricker_series(
    dt: float, n_steps: int, f0: float, t0: float, amp: float = 1.0
) -> list[float]:
    """Generate discrete Ricker wavelet time series."""
    return [amp * ricker_wavelet(k * dt, f0, t0) for k in range(n_steps + 4)]


def rel_l2_and_peak_error(u_test: np.ndarray, u_ref: np.ndarray) -> tuple[float, float]:
    """Compute relative L2 and relative peak reflection errors over waveform."""
    diff = u_test - u_ref
    l2_ref = float(np.linalg.norm(u_ref))
    l2_err = float(np.linalg.norm(diff)) / max(l2_ref, 1.0e-30)

    mag_ref = np.linalg.norm(u_ref, axis=1) if u_ref.ndim > 1 else np.abs(u_ref)
    mag_diff = np.linalg.norm(diff, axis=1) if diff.ndim > 1 else np.abs(diff)
    peak_ref = float(np.max(mag_ref))
    peak_err = float(np.max(mag_diff)) / max(peak_ref, 1.0e-30)
    return l2_err, peak_err


def parse_up_cases_txt(path: str) -> list[dict[str, Any]]:
    """Parse `tests/ladruno_cpml_up_cases.txt` into structured dictionaries."""
    with open(path, "r", encoding="utf-8") as f:
        raw_lines = [
            ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")
        ]

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

        biot_toks = raw_lines[idx].split()
        assert biot_toks[0] == "BiotParams"
        Kf = float(biot_toks[1])
        poro = float(biot_toks[2])
        rhoF = float(biot_toks[3])
        biotAlpha = float(biot_toks[4])
        Ks = float(biot_toks[5])
        pOrder = int(biot_toks[6])
        formulation = int(biot_toks[7])
        stabMode = int(biot_toks[8])
        stabValue = float(biot_toks[9])
        kbar = tuple(float(x) for x in biot_toks[10 : 10 + ndm])
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
        coords = np.array(
            [float(x) for x in coord_toks[1:]], dtype=np.float64
        ).reshape(n_nodes, ndm)
        idx += 1

        nm_toks = raw_lines[idx].split()
        assert nm_toks[0] == "NewmarkCoeffs"
        c1, c2, c3 = float(nm_toks[1]), float(nm_toks[2]), float(nm_toks[3])
        idx += 1

        ns_toks = raw_lines[idx].split()
        assert ns_toks[0] == "NumSteps"
        n_steps = int(ns_toks[1])
        idx += 1

        n_dof, n_p, u_off, p_off = ref.compute_dof_map_up(topo, pOrder)
        steps: list[dict[str, Any]] = []
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
                    step_data[key] = np.array(vals, dtype=np.float64).reshape(
                        n_dof, n_dof
                    )
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
                "n_dof": n_dof,
                "n_p": n_p,
                "u_off": u_off,
                "p_off": p_off,
                "lumped": lumped,
                "E": E,
                "nu": nu,
                "rho": rho,
                "biot": ref.BiotParams(
                    Kf=Kf,
                    poro=poro,
                    rhoF=rhoF,
                    biotAlpha=biotAlpha,
                    Ks=Ks,
                    pOrder=pOrder,
                    formulation=formulation,
                    stabMode=stabMode,
                    stabValue=stabValue,
                    kbar=kbar,
                ),
                "pml": ref.PMLParams(
                    m_deg=m_deg, R0=R0, kappa_max=kappa_max, alpha_max=alpha_max
                ),
                "geom": ref.GeomSpec(mode=geom_mode, params=geom_params),
                "coords": coords,
                "c1": c1,
                "c2": c2,
                "c3": c3,
                "steps": steps,
            }
        )
    return cases


def create_up_nodes(topo: str, coords: np.ndarray, p_order: int) -> list[int]:
    """Create nodes in the current OpenSees domain with proper ndf for equal-order or Taylor-Hood."""
    ndm = coords.shape[1]
    n_nodes = coords.shape[0]
    n_vert = ref.N_VERT[topo]
    node_tags = list(range(1, n_nodes + 1))

    ops.model("basic", "-ndm", ndm, "-ndf", ndm + 1)
    if p_order == 0:
        for idx, ntag in enumerate(node_tags):
            ops.node(ntag, *[float(x) for x in coords[idx]])
    else:
        for idx in range(n_vert):
            ops.node(node_tags[idx], *[float(x) for x in coords[idx]])
        ops.model("basic", "-ndm", ndm, "-ndf", ndm)
        for idx in range(n_vert, n_nodes):
            ops.node(node_tags[idx], *[float(x) for x in coords[idx]])
    return node_tags


def build_cpml_up_element_args(
    ele_tag: int,
    node_tags: list[int],
    topo: str,
    ndm: int,
    E: float,
    nu: float,
    rho: float,
    biot: ref.BiotParams,
    pml: ref.PMLParams,
    geom: ref.GeomSpec,
    lumped: int,
    use_mat_tag: int | None = None,
    use_perm_h_gamma_w: float | None = None,
    body_force: tuple[float, ...] | None = None,
    fluid_body_force: tuple[float, ...] | None = None,
    stage: int | None = None,
) -> list[Any]:
    """Assemble arguments for `ops.element('LadrunoCPML', ..., '-up', ...)`."""
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
    if stage is not None:
        args.extend(["-stage", stage])
    if body_force is not None:
        args.extend(["-body", *body_force])

    args.extend(
        [
            "-up",
            "-Kf",
            biot.Kf,
            "-poro",
            biot.poro,
            "-rhoF",
            biot.rhoF,
            "-alphaB",
            biot.biotAlpha,
        ]
    )
    if biot.Ks > 0.0:
        args.extend(["-Ks", biot.Ks])

    if use_perm_h_gamma_w is not None:
        gw = float(use_perm_h_gamma_w)
        perm_h = [float(k) * gw for k in biot.kbar]
        args.extend(["-permH", *perm_h, "-gammaW", gw])
    else:
        args.extend(["-perm", *[float(k) for k in biot.kbar]])

    args.extend(
        [
            "-pOrder",
            "linear" if biot.pOrder == 1 else "equal",
            "-formulation",
            "bbar" if biot.formulation == 1 else "std",
        ]
    )
    if biot.stabMode == 0:
        args.extend(["-stab", "off"])
    elif biot.stabMode == 1:
        args.extend(["-stab", "auto", biot.stabValue])
    elif biot.stabMode == 2:
        args.extend(["-stab", "manual", biot.stabValue])

    if fluid_body_force is not None:
        args.extend(["-fluidBody", *fluid_body_force])
    return args


def extract_oracle_r_disp_up(
    oracle: ref.ReferenceCPMLUPElement,
    R_int_ref: np.ndarray,
    M_eff_ref: np.ndarray,
    u_vec: np.ndarray,
    v_vec: np.ndarray,
    a_vec: np.ndarray,
) -> np.ndarray:
    """Extract `R_disp = R_int - M_eff @ a - R_vel` from a `ReferenceCPMLUPElement`
    whose `evaluate_trial` was just called at `(u_vec, v_vec, a_vec)`, computing
    the `p`-block Darcy resisting force `f_p_H` directly to avoid subtracting
    `O(1e-4)` storage/coupling velocity terms from `R_int_ref[p]`.
    """
    ndm = oracle.ndm
    n_nodes = oracle.n_nodes
    n_p = oracle.n_p
    r_vel = np.zeros(oracle.n_dof, dtype=np.float64)
    f_p_H = np.zeros(n_p, dtype=np.float64)
    p_vec = np.array([u_vec[oracle.p_off[b]] for b in range(n_p)], dtype=np.float64)

    for gp in oracle.gps:
        w = gp.w_vol
        k0, k1, k2 = gp.kappa[0], gp.kappa[1], gp.kappa[2]
        a_vel_g = np.zeros(ndm, dtype=np.float64)
        for i in range(ndm):
            vg0 = gp.vel_g0[i].c0 * gp.vel_g0[i].x_trial + gp.vel_g0[i].c1 * gp.vel_g0[i].zeta_trial
            vg1 = gp.vel_g1[i].c0 * gp.vel_g1[i].x_trial + gp.vel_g1[i].c1 * gp.vel_g1[i].zeta_trial
            vg2 = gp.vel_g2[i].c0 * gp.vel_g2[i].x_trial + gp.vel_g2[i].c1 * gp.vel_g2[i].zeta_trial
            a_vel_g[i] = (k1 * k2) * vg0 + (k0 * k2) * vg1 + (k0 * k1) * vg2

        rho_w = oracle.rho * w
        for a in range(n_nodes):
            ua = oracle.u_off[a]
            r_vel[ua : ua + ndm] += rho_w * gp.N_u[a] * a_vel_g

        g_p_local = p_vec @ gp.Phi_p_local
        for j in range(ndm):
            flux_j = 0.0
            for l in range(ndm):
                k_jl = oracle.kbar_local[j, l]
                casc_jl = gp.F_p_cascades[(j, l)]
                if len(casc_jl.stages) == 0:
                    y_val = float(g_p_local[l])
                else:
                    st_last = casc_jl.stages[-1]
                    y_val = st_last.c0 * st_last.x_trial + st_last.c1 * st_last.zeta_trial
                flux_j += k_jl * y_val
            f_p_H += (w * flux_j) * gp.Phi_p_local[:, j]

    r_disp = R_int_ref - (M_eff_ref @ a_vec) - r_vel
    for b in range(n_p):
        r_disp[oracle.p_off[b]] = f_p_H[b]
    return r_disp


def read_up_nodal_vector(
    node_tags: list[int], topo: str, p_order: int, kind: str = "disp"
) -> np.ndarray:
    """Read interleaved [u..., p?] nodal vector from OpenSees."""
    ndm, n_nodes, _, _ = ref.TOPOLOGIES[topo]
    n_dof, _, u_off, p_off = ref.compute_dof_map_up(topo, p_order)
    out = np.zeros(n_dof, dtype=np.float64)
    fn = (
        ops.nodeDisp
        if kind == "disp"
        else (ops.nodeVel if kind == "vel" else ops.nodeAccel)
    )
    for a, ntag in enumerate(node_tags):
        ua = u_off[a]
        for d in range(ndm):
            out[ua + d] = fn(ntag, d + 1)
        if p_off[a] >= 0:
            out[p_off[a]] = 0.0 if kind == "accel" else fn(ntag, ndm + 1)
    return out


def apply_up_nodal_load(
    node_tags: list[int], topo: str, p_order: int, f_vec: np.ndarray
) -> None:
    """Apply interleaved [f_u..., f_p?] nodal loads in OpenSees."""
    ndm, _, _, _ = ref.TOPOLOGIES[topo]
    _, _, u_off, p_off = ref.compute_dof_map_up(topo, p_order)
    for a, ntag in enumerate(node_tags):
        ua = u_off[a]
        f_u = [float(f_vec[ua + d]) for d in range(ndm)]
        if p_off[a] >= 0:
            ops.load(ntag, *f_u, float(f_vec[p_off[a]]))
        else:
            ops.load(ntag, *f_u)


# ============================================================================ #
# 1. Gate 0a: Interior (xi = 0) Exact Equivalence vs. `LadrunoUP`              #
# ============================================================================ #


def test_up_gate0_interior_equivalence_vs_ladrunoup() -> None:
    """Gate 0a (xi = 0 Exact Equivalence vs. `LadrunoUP`):
    Across all 5 shared topologies (`T3, Q4, H8` with `pOrder="equal"`, and
    `BT6, BTET10` with `pOrder="linear"` Taylor-Hood), and for both `-formulation std`
    and `-formulation bbar` (on `Q4` and `H8`), compare an interior (`xi_m = 0`)
    `LadrunoCPML -up` element against `LadrunoUP` (`-dynSeepage off`) over a 3-step
    `TransientAnalysis` (`Newmark(0.5, 0.25)`) under non-trivial nodal loads on both
    `u` and `p` DOFs, body forces, and fluidBody forces. Assert assembled dynamic
    tangent `printA()`, element resisting force `eleResponse(1, "force")`, nodal
    displacements `u`, and nodal pore pressures `p` match `LadrunoUP` to <= 1e-13!
    """
    configs = [
        ("T3", 0, "std", False),
        ("T3", 0, "std", True),
        ("Q4", 0, "std", False),
        ("Q4", 0, "bbar", False),
        ("H8", 0, "std", False),
        ("H8", 0, "bbar", False),
        ("BT6", 1, "std", False),
        ("BTET10", 1, "std", False),
    ]

    E = 7.5e7
    nu = 0.28
    rho = 2100.0
    Kf = 2.2e8
    poro = 0.32
    rhoF = 1000.0
    alphaB = 0.88
    Ks = 3.6e9
    dt = 0.002

    for topo, p_order, form, lumped in configs:
        ndm, n_nodes, _, _ = ref.TOPOLOGIES[topo]
        coords = ref.make_reference_coords(topo)
        n_dof, n_p, u_off, p_off = ref.compute_dof_map_up(topo, p_order)
        u_idx = [u_off[a] + d for a in range(n_nodes) for d in range(ndm)]
        p_idx = [p_off[b] for b in range(n_p)]

        kbar = (2.5e-5, 1.4e-5) if ndm == 2 else (2.5e-5, 1.4e-5, 9.0e-6)
        body = (0.4, -9.81) if ndm == 2 else (0.3, -0.2, -9.81)
        fluid_body = (0.2, -9.81) if ndm == 2 else (0.15, -0.1, -9.81)
        thick = 1.35 if ndm == 2 else 1.0

        rng = np.random.default_rng(2026 + n_nodes * 31 + (1 if form == "bbar" else 0))
        step_loads = [rng.uniform(0.5e3, 1.5e3, size=n_dof) for _ in range(3)]
        for s_load in step_loads:
            for b in range(n_p):
                s_load[p_off[b]] *= 1.0e-4

        def run_twin(which: str) -> list[dict[str, np.ndarray]]:
            ops.wipe()
            node_tags = create_up_nodes(topo, coords, p_order)
            # Restrain rigid-body translations/rotations and constant-p Darcy mode at node 1/2/3
            if ndm == 2:
                ops.fix(1, 1, 1, 1)
                if p_off[1] >= 0:
                    ops.fix(2, 0, 1, 0)
                else:
                    ops.fix(2, 0, 1)
                n_fixed = 4
            else:
                ops.fix(1, 1, 1, 1, 1)
                if p_off[1] >= 0:
                    ops.fix(2, 0, 1, 1, 0)
                    ops.fix(3, 1, 0, 1, 0)
                else:
                    ops.fix(2, 0, 1, 1)
                    ops.fix(3, 1, 0, 1)
                n_fixed = 8

            ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

            if which == "up":
                up_args: list[Any] = [
                    "LadrunoUP",
                    1,
                    *node_tags,
                    1,
                    "-Kf",
                    Kf,
                    "-poro",
                    poro,
                    "-rhoF",
                    rhoF,
                    "-perm",
                    *kbar,
                    "-alpha",
                    alphaB,
                    "-Ks",
                    Ks,
                    "-body",
                    *body,
                    "-fluidBody",
                    *fluid_body,
                    "-formulation",
                    form,
                    "-pOrder",
                    "linear" if p_order == 1 else "equal",
                    "-dynSeepage",
                    "off",
                ]
                if ndm == 2:
                    up_args.extend(["-thick", thick])
                if p_order == 0:
                    up_args.extend(["-stab", "auto", 0.25])
                if lumped:
                    up_args.append("-lumped")
                ops.element(*up_args)
            else:
                box_bounds = (
                    [-100.0, 100.0, -100.0, 100.0]
                    if ndm == 2
                    else [-100.0, 100.0, -100.0, 100.0, -100.0, 100.0]
                )
                cpml_args: list[Any] = [
                    "LadrunoCPML",
                    1,
                    *node_tags,
                    "-mat",
                    1,
                    "-box",
                    *box_bounds,
                    "-L",
                    2.0,
                    "-shape",
                    topo,
                    "-body",
                    *body,
                    "-up",
                    "-Kf",
                    Kf,
                    "-poro",
                    poro,
                    "-rhoF",
                    rhoF,
                    "-perm",
                    *kbar,
                    "-alphaB",
                    alphaB,
                    "-Ks",
                    Ks,
                    "-fluidBody",
                    *fluid_body,
                    "-formulation",
                    form,
                    "-pOrder",
                    "linear" if p_order == 1 else "equal",
                ]
                if ndm == 2:
                    cpml_args.extend(["-thick", thick])
                if p_order == 0:
                    cpml_args.extend(["-stab", "auto", 0.25])
                else:
                    cpml_args.extend(["-stab", "off"])
                if lumped:
                    cpml_args.append("-lumped")
                ops.element(*cpml_args)

            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("FullGeneral")
            ops.test("NormDispIncr", 1.0e-14, 10, 0)
            ops.algorithm("Linear")
            ops.integrator("Newmark", 0.5, 0.25)
            ops.analysis("Transient", "-noWarnings")

            snaps: list[dict[str, np.ndarray]] = []
            for s_idx, f_load in enumerate(step_loads):
                ts_tag = 10 + s_idx
                pat_tag = 20 + s_idx
                ops.timeSeries("Constant", ts_tag)
                ops.pattern("Plain", pat_tag, ts_tag)
                apply_up_nodal_load(node_tags, topo, p_order, f_load)
                assert ops.analyze(1, dt) == 0
                ops.remove("loadPattern", pat_tag)
                ops.remove("timeSeries", ts_tag)

                uv = read_up_nodal_vector(node_tags, topo, p_order, "disp")
                n_free = n_dof - n_fixed
                k_a = np.array(ops.printA("-ret"), dtype=np.float64).reshape(
                    (n_free, n_free), order="F"
                )
                f_ele = np.array(ops.eleResponse(1, "force"), dtype=np.float64)
                snaps.append(
                    {
                        "u": uv[u_idx],
                        "p": uv[p_idx],
                        "printA": k_a,
                        "force": f_ele,
                    }
                )
            return snaps

        snaps_up = run_twin("up")
        snaps_cpml = run_twin("cpml")

        for s_idx in range(3):
            for key in ("printA", "force", "u", "p"):
                err = rel_err(snaps_cpml[s_idx][key], snaps_up[s_idx][key])
                assert err <= 1.0e-12, (
                    f"Gate 0a ({topo}, form={form}, lumped={lumped}) step {s_idx} "
                    f"{key} mismatch: {err:.3e} > 1e-12"
                )


# ============================================================================ #
# 2. Gate 0b: All 14 Multi-Step Oracle Cases Through `opensees.pyd` + 4-Block FD#
# ============================================================================ #

PARSED_UP_CASES = parse_up_cases_txt(UP_CASES_TXT_PATH)


def test_up_gate0_oracle_14_cases_and_fd_tangent() -> None:
    """Gate 0b (All 14 Multi-Step Reference Oracle Cases Through `opensees.pyd` + 4-Block FD Tangent):
    1. Parse `tests/ladruno_cpml_up_cases.txt` (all 14 cases covering `T3, Q4, T6, BT6, TET4,
       H8, TET10, BTET10, H20`, `pOrder = 0, 1`, `formulation = 0, 1`, `stabMode = 0, 1, 2`,
       `box` and `plane`, variable dt) and verify exact match with `ReferenceCPMLUPElement`.
    2. Drive all 14 cases through `opensees.pyd` multi-step `TransientAnalysis` and assert
       `ops.printA()`, `ops.eleResponse(1, "force")`, `ops.eleResponse(1, "basicForce")`,
       `ops.eleResponse(1, "stiff")`, `ops.eleResponse(1, "damping")`, `ops.eleResponse(1, "mass")`,
       and `ops.eleResponse(1, "energy")` match across all 4 blocks (`uu, up, pu, pp`) to <= 1e-10!
    3. Verify `-mat $matTag` and `-permH ... -gammaW ...` sugar equivalence in `-up` mode.
    4. Verify in-engine 4-block (`uu, up, pu, pp`) finite-difference tangent and exact skew-duality
       `K_eff[u, p] == -C_eff[p, u]^T` in `opensees.pyd`.
    """
    assert len(PARSED_UP_CASES) == 14
    topos_seen = {c["topo"] for c in PARSED_UP_CASES}
    assert topos_seen == {
        "T3",
        "Q4",
        "T6",
        "BT6",
        "TET4",
        "H8",
        "TET10",
        "BTET10",
        "H20",
    }

    # --- Part 1: Golden file `ladruno_cpml_up_cases.txt` vs `ReferenceCPMLUPElement` ---
    for case_data in PARSED_UP_CASES:
        elem_ref = ref.ReferenceCPMLUPElement(
            topo=case_data["topo"],
            coords=case_data["coords"],
            E=case_data["E"],
            nu=case_data["nu"],
            rho=case_data["rho"],
            biot=case_data["biot"],
            pml=case_data["pml"],
            geom=case_data["geom"],
            lumped=case_data["lumped"],
        )
        for st in case_data["steps"]:
            R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext = elem_ref.evaluate_trial(
                st["U"],
                st["V"],
                st["A"],
                st["dt"],
                case_data["c1"],
                case_data["c2"],
                case_data["c3"],
            )
            elem_ref.commit()
            assert rel_err_scaled(R_int, st["R_int"]) <= 1.0e-12
            assert rel_err_scaled(K_eff, st["K_eff"]) <= 1.0e-12
            assert rel_err_scaled(C_eff, st["C_eff"]) <= 1.0e-12
            assert rel_err_scaled(M_eff, st["M_eff"]) <= 1.0e-12
            assert rel_err_scaled(K_T_dyn, st["K_T_dyn"]) <= 1.0e-12
            assert rel_err_scaled(P_ext, st["P_ext"]) <= 1.0e-12

    # --- Part 2: All 14 Cases Through `opensees.pyd` Across All 4 Blocks (`uu, up, pu, pp`) ---
    for case_idx, case_data in enumerate(PARSED_UP_CASES):
        topo = case_data["topo"]
        ndm = case_data["ndm"]
        n_nodes = case_data["n_nodes"]
        n_dof = case_data["n_dof"]
        n_p = case_data["n_p"]
        u_off = case_data["u_off"]
        p_off = case_data["p_off"]
        biot = case_data["biot"]
        coords = case_data["coords"]

        u_idx = [u_off[a] + d for a in range(n_nodes) for d in range(ndm)]
        p_idx = [p_off[b] for b in range(n_p)]

        ops.wipe()
        node_tags = create_up_nodes(topo, coords, biot.pOrder)

        # Exercise `-mat` and `-permH/-gammaW` on alternating cases as well as inline flags
        use_mat = 10 if (case_idx % 2 == 1) else None
        if use_mat is not None:
            ops.nDMaterial(
                "ElasticIsotropic",
                use_mat,
                case_data["E"],
                case_data["nu"],
                case_data["rho"],
            )
        use_gw = 9810.0 if (case_idx % 3 == 1) else None

        ele_args = build_cpml_up_element_args(
            ele_tag=1,
            node_tags=node_tags,
            topo=topo,
            ndm=ndm,
            E=case_data["E"],
            nu=case_data["nu"],
            rho=case_data["rho"],
            biot=biot,
            pml=case_data["pml"],
            geom=case_data["geom"],
            lumped=case_data["lumped"],
            use_mat_tag=use_mat,
            use_perm_h_gamma_w=use_gw,
        )
        ops.element("LadrunoCPML", *ele_args)

        oracle = ref.ReferenceCPMLUPElement(
            topo=topo,
            coords=coords,
            E=case_data["E"],
            nu=case_data["nu"],
            rho=case_data["rho"],
            biot=biot,
            pml=case_data["pml"],
            geom=case_data["geom"],
            lumped=case_data["lumped"],
        )
        oracle_ops = ref.ReferenceCPMLUPElement(
            topo=topo,
            coords=coords,
            E=case_data["E"],
            nu=case_data["nu"],
            rho=case_data["rho"],
            biot=biot,
            pml=case_data["pml"],
            geom=case_data["geom"],
            lumped=case_data["lumped"],
        )

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient", "-noWarnings")

        u_c = np.zeros(n_dof, dtype=np.float64)
        v_c = np.zeros(n_dof, dtype=np.float64)
        a_c = np.zeros(n_dof, dtype=np.float64)

        for s_idx, step_info in enumerate(case_data["steps"]):
            dt = float(step_info["dt"])
            u_target = step_info["U"].copy()

            c1 = 1.0
            c2 = 2.0 / dt
            c3 = 4.0 / (dt * dt)

            a_target = c3 * (u_target - u_c) - (4.0 / dt) * v_c - a_c
            v_target = v_c + 0.5 * dt * (a_c + a_target)
            # LadrunoCPML::gatherKinematics zeroes p-slot accelerations (M_eff has zero p-rows/cols)
            a_oracle = a_target.copy()
            a_oracle[p_idx] = 0.0

            R_int_ref, K_eff_ref, C_eff_ref, M_eff_ref, K_T_dyn_ref, P_ext_ref = (
                oracle.evaluate_trial(
                    u_target, v_target, a_oracle, dt, c1, c2, c3
                )
            )
            oracle.commit()

            ts_tag = 100 + s_idx
            pat_tag = 200 + s_idx
            ops.timeSeries("Constant", ts_tag)
            ops.pattern("Plain", pat_tag, ts_tag)
            apply_up_nodal_load(node_tags, topo, biot.pOrder, R_int_ref)

            rc = ops.analyze(1, dt)
            assert rc == 0, f"{case_data['name']} step {s_idx}: analyze failed rc={rc}"

            ops.remove("loadPattern", pat_tag)
            ops.remove("timeSeries", ts_tag)

            u_ops = read_up_nodal_vector(node_tags, topo, biot.pOrder, "disp")
            v_ops = read_up_nodal_vector(node_tags, topo, biot.pOrder, "vel")
            a_ops = read_up_nodal_vector(node_tags, topo, biot.pOrder, "accel")

            # For unfixed 1-GP TET4 with lumped=0, K_T_dyn[u, u] has a 1D kernel (z-rotation with zero strain and zero centroid displacement)
            _, s_vals, vt_mat = np.linalg.svd(K_T_dyn_ref[np.ix_(u_idx, u_idx)])
            null_mask = s_vals < 1.0e-12 * s_vals[0]
            du_u = u_ops[u_idx] - u_target[u_idx]
            dv_u = v_ops[u_idx] - v_target[u_idx]
            da_u = a_ops[u_idx] - a_oracle[u_idx]
            if np.any(null_mask):
                v_null = vt_mat[null_mask]
                du_u = du_u - v_null.T @ (v_null @ du_u)
                dv_u = dv_u - v_null.T @ (v_null @ dv_u)
                da_u = da_u - v_null.T @ (v_null @ da_u)

            assert float(np.linalg.norm(du_u)) / max(float(np.linalg.norm(u_target[u_idx])), 1.0e-30) <= 1.0e-10
            assert rel_err(u_ops[p_idx], u_target[p_idx]) <= 5.0e-6
            assert float(np.linalg.norm(dv_u)) / max(float(np.linalg.norm(v_target[u_idx])), 1.0e-30) <= 1.0e-10
            assert rel_err(v_ops[p_idx], v_target[p_idx]) <= 5.0e-6
            assert float(np.linalg.norm(da_u)) / max(float(np.linalg.norm(a_oracle[u_idx])), 1.0e-30) <= 1.0e-10

            a_ops_zero_p = a_ops.copy()
            a_ops_zero_p[p_idx] = 0.0
            R_int_ops_ref, _, _, M_eff_ops_ref, _, _ = oracle_ops.evaluate_trial(
                u_ops, v_ops, a_ops_zero_p, dt, c1, c2, c3
            )
            R_disp_ops_ref = extract_oracle_r_disp_up(
                oracle_ops, R_int_ops_ref, M_eff_ops_ref, u_ops, v_ops, a_ops_zero_p
            )
            oracle_ops.commit()

            # Check K_T_dyn (`printA`, Fortran column-major from FullGenLinSOE) and element matrices across all 4 blocks (`uu, up, pu, pp`)
            K_T_ops = np.array(ops.printA("-ret"), dtype=np.float64).reshape(
                (n_dof, n_dof), order="F"
            )
            K_eff_ops = np.array(
                ops.eleResponse(1, "stiff"), dtype=np.float64
            ).reshape(n_dof, n_dof)
            C_eff_ops = np.array(
                ops.eleResponse(1, "damping"), dtype=np.float64
            ).reshape(n_dof, n_dof)
            M_eff_ops = np.array(
                ops.eleResponse(1, "mass"), dtype=np.float64
            ).reshape(n_dof, n_dof)

            assert rel_err_scaled(K_T_ops, K_T_dyn_ref) <= 1.0e-10
            assert rel_err_scaled(K_eff_ops, K_eff_ref) <= 1.0e-10
            assert rel_err_scaled(C_eff_ops, C_eff_ref) <= 1.0e-10
            assert rel_err_scaled(M_eff_ops, M_eff_ref) <= 1.0e-10

            # Block-by-block (`uu, up, pu, pp`) relative error checks on K_T_dyn, K_eff, C_eff
            for blk_name, r_ix, c_ix in (
                ("uu", u_idx, u_idx),
                ("up", u_idx, p_idx),
                ("pu", p_idx, u_idx),
                ("pp", p_idx, p_idx),
            ):
                kt_blk_ops = K_T_ops[np.ix_(r_ix, c_ix)]
                kt_blk_ref = K_T_dyn_ref[np.ix_(r_ix, c_ix)]
                assert rel_err(kt_blk_ops, kt_blk_ref) <= 1.0e-10, (
                    f"{case_data['name']} step {s_idx} K_T_dyn[{blk_name}] err="
                    f"{rel_err(kt_blk_ops, kt_blk_ref):.3e}"
                )
                if abs(c2 - case_data["c2"]) <= 1.0e-9 and abs(c3 - case_data["c3"]) <= 1.0e-6:
                    kt_blk_txt = step_info["K_T_dyn"][np.ix_(r_ix, c_ix)]
                    assert rel_err(kt_blk_ops, kt_blk_txt) <= 1.0e-10

            # Exact skew-duality in OpenSees: K_eff[u, p] == -C_eff[p, u]^T
            K_up_ops = K_eff_ops[np.ix_(u_idx, p_idx)]
            C_pu_ops = C_eff_ops[np.ix_(p_idx, u_idx)]
            assert rel_err(K_up_ops, -C_pu_ops.T) <= 1.0e-13

            # Check `force` (R_disp), `basicForce` (R_int), and `energy` (P_ext) on both u and p blocks
            f_disp_ops = np.array(ops.eleResponse(1, "force"), dtype=np.float64)
            f_int_ops = np.array(ops.eleResponse(1, "basicForce"), dtype=np.float64)
            p_ext_ops = float(ops.eleResponse(1, "energy")[0])

            assert rel_err(f_disp_ops[u_idx], R_disp_ops_ref[u_idx]) <= 1.0e-10
            assert rel_err(f_disp_ops[p_idx], R_disp_ops_ref[p_idx]) <= 1.0e-10
            assert rel_err(f_int_ops[u_idx], R_int_ref[u_idx]) <= 1.0e-10
            assert rel_err(f_int_ops[p_idx], R_int_ref[p_idx]) <= 1.0e-10
            assert rel_err_scaled(p_ext_ops, P_ext_ref) <= 1.0e-10

            u_c, v_c, a_c = u_target, v_target, a_target

    # --- Part 3: Explicit `-mat` + `-permH/-gammaW` vs inline `-E/-nu/-rho` + `-perm` equivalence ---
    cs0 = PARSED_UP_CASES[1]  # UP_Q4_box_corner_eq_bbar
    dt0 = float(cs0["steps"][0]["dt"])
    f_probe = cs0["steps"][0]["R_int"]

    def run_sugar_variant(use_sugar: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ops.wipe()
        ntags = create_up_nodes(cs0["topo"], cs0["coords"], cs0["biot"].pOrder)
        mtag = None
        gw = None
        if use_sugar:
            mtag = 55
            gw = 9810.0
            ops.nDMaterial("ElasticIsotropic", mtag, cs0["E"], cs0["nu"], cs0["rho"])
        eargs = build_cpml_up_element_args(
            ele_tag=1,
            node_tags=ntags,
            topo=cs0["topo"],
            ndm=cs0["ndm"],
            E=cs0["E"],
            nu=cs0["nu"],
            rho=cs0["rho"],
            biot=cs0["biot"],
            pml=cs0["pml"],
            geom=cs0["geom"],
            lumped=cs0["lumped"],
            use_mat_tag=mtag,
            use_perm_h_gamma_w=gw,
        )
        ops.element("LadrunoCPML", *eargs)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient", "-noWarnings")
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        apply_up_nodal_load(ntags, cs0["topo"], cs0["biot"].pOrder, f_probe)
        assert ops.analyze(1, dt0) == 0
        return (
            read_up_nodal_vector(ntags, cs0["topo"], cs0["biot"].pOrder, "disp"),
            np.array(ops.eleResponse(1, "basicForce"), dtype=np.float64),
            np.array(ops.printA("-ret"), dtype=np.float64),
        )

    u_inline, r_inline, kt_inline = run_sugar_variant(False)
    u_sugar, r_sugar, kt_sugar = run_sugar_variant(True)
    assert rel_err(u_sugar, u_inline) <= 1.0e-13
    assert rel_err(r_sugar, r_inline) <= 1.0e-13
    assert rel_err(kt_sugar, kt_inline) <= 1.0e-13

    # --- Part 4: In-Engine 4-Block (`uu, up, pu, pp`) FD Tangent Verification in `opensees.pyd` ---
    # Test across 2D corner (Q4 bbar), 2D Taylor-Hood (BT6), 3D corner (H8 bbar), and 3D Taylor-Hood (BTET10)
    fd_case_indices = [1, 6, 8, 11]
    for c_i in fd_case_indices:
        cs = PARSED_UP_CASES[c_i]
        topo = cs["topo"]
        ndm = cs["ndm"]
        n_nodes = cs["n_nodes"]
        n_dof = cs["n_dof"]
        n_p = cs["n_p"]
        u_off = cs["u_off"]
        p_off = cs["p_off"]
        biot = cs["biot"]
        coords = cs["coords"]
        dt = float(cs["steps"][0]["dt"])
        u0 = cs["steps"][0]["U"].copy()

        u_idx = [u_off[a] + d for a in range(n_nodes) for d in range(ndm)]
        p_idx = [p_off[b] for b in range(n_p)]

        def eval_ops_step(
            f_load: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            ops.wipe()
            ntags = create_up_nodes(topo, coords, biot.pOrder)
            eargs = build_cpml_up_element_args(
                ele_tag=1,
                node_tags=ntags,
                topo=topo,
                ndm=ndm,
                E=cs["E"],
                nu=cs["nu"],
                rho=cs["rho"],
                biot=biot,
                pml=cs["pml"],
                geom=cs["geom"],
                lumped=cs["lumped"],
            )
            ops.element("LadrunoCPML", *eargs)
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("FullGeneral")
            ops.algorithm("Linear")
            ops.integrator("Newmark", 0.5, 0.25)
            ops.analysis("Transient", "-noWarnings")
            ops.timeSeries("Constant", 1)
            ops.pattern("Plain", 1, 1)
            apply_up_nodal_load(ntags, topo, biot.pOrder, f_load)
            assert ops.analyze(1, dt) == 0
            u_out = read_up_nodal_vector(ntags, topo, biot.pOrder, "disp")
            r_int_out = np.array(ops.eleResponse(1, "basicForce"), dtype=np.float64)
            k_t_out = np.array(ops.printA("-ret"), dtype=np.float64).reshape(
                (n_dof, n_dof), order="F"
            )
            return u_out, r_int_out, k_t_out

        oracle_fd = ref.ReferenceCPMLUPElement(
            topo=topo,
            coords=coords,
            E=cs["E"],
            nu=cs["nu"],
            rho=cs["rho"],
            biot=biot,
            pml=cs["pml"],
            geom=cs["geom"],
            lumped=cs["lumped"],
        )
        c1, c2, c3 = 1.0, 2.0 / dt, 4.0 / (dt * dt)
        a0 = c3 * u0
        a0[p_idx] = 0.0
        f0, _, _, _, _, _ = oracle_fd.evaluate_trial(
            u0, c2 * u0, a0, dt, c1, c2, c3
        )
        _, _, K_T_ops = eval_ops_step(f0)

        dU_mat = np.zeros((n_dof, n_dof), dtype=np.float64)
        dR_mat = np.zeros((n_dof, n_dof), dtype=np.float64)
        for j in range(n_dof):
            eps_j = 1.0e-5 * (cs["E"] if j in p_idx else 1.0)
            df = K_T_ops[:, j] * eps_j
            u_plus, r_plus, _ = eval_ops_step(f0 + df)
            u_minus, r_minus, _ = eval_ops_step(f0 - df)
            dU_mat[:, j] = (u_plus - u_minus) / (2.0 * eps_j)
            dR_mat[:, j] = (r_plus - r_minus) / (2.0 * eps_j)

        K_fd_ops = dR_mat @ np.linalg.inv(dU_mat)

        for blk_name, r_ix, c_ix in (
            ("uu", u_idx, u_idx),
            ("up", u_idx, p_idx),
            ("pu", p_idx, u_idx),
            ("pp", p_idx, p_idx),
        ):
            err_blk = rel_err(
                K_fd_ops[np.ix_(r_ix, c_ix)], K_T_ops[np.ix_(r_ix, c_ix)]
            )
            assert err_blk <= 1.0e-7, (
                f"{cs['name']} OpenSees FD tangent block [{blk_name}] err={err_blk:.3e}"
            )


# ============================================================================ #
# 3. Anchor UP-1: 1D Saturated Soil Column P1 & S-Wave Reflection + 37 deg     #
# ============================================================================ #


def _run_1d_poroelastic_column(
    wave_type: str,
    mode: str,
    theta_deg: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Run 1D saturated poroelastic soil column waveguide (`LadrunoUP` interior x in [0, 10m],
    terminated by `LadrunoCPML -up` absorbing layer x in [10, 15m] of thickness L_pml = 5.0m
    with 10 elements, or extended `LadrunoUP` column x in [0, 60m], or truncated fixed column
    x in [0, 10m]).

    Returns `(u_rec, p_rec)` at receiver station x_r = 6.0 m.
    """
    # Undrained fast P1 wave speed Vp1 = sqrt((lambda + 2*mu + alpha^2 * Qbar)/rho) = 1000 m/s,
    # S-wave speed Vs = sqrt(mu / rho) = 447.21 m/s
    E = 5.0e8
    nu = 0.25
    rho = 1000.0
    Kf = 1.2e8
    poro = 0.30
    rhoF = 1000.0
    biotAlpha = 1.0
    perm = 1.0e-8

    dx = 0.5
    Hy = 1.0
    Lx_int = 10.0
    Lpml = 5.0
    if mode == "extended_ref":
        L_total = 45.0
    elif mode == "truncated_fixed":
        L_total = Lx_int
    elif mode in ("box", "plane"):
        L_total = Lx_int + Lpml
    else:
        raise ValueError(f"Unknown mode={mode}")

    nx_elem = int(round(L_total / dx))
    nx_int = int(round(Lx_int / dx))

    dt = 0.00025
    n_steps = 260  # T_sim = 0.065 s: covers full P1 and S wave round-trip from x=10m back to x_r=6m
    f0 = 30.0
    t0 = 0.025

    theta = math.radians(theta_deg)
    ct = math.cos(theta)
    st = math.sin(theta)
    R_rot = np.array([[ct, -st], [st, ct]], dtype=np.float64)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    def ntag(i: int, j: int) -> int:
        return j * (nx_elem + 1) + i + 1

    for j in range(2):
        sy = float(j) * Hy
        for i in range(nx_elem + 1):
            sx = float(i) * dx
            xg = R_rot @ np.array([sx, sy], dtype=np.float64)
            ops.node(ntag(i, j), float(xg[0]), float(xg[1]))

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

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
        in_pml = mode in ("box", "plane") and (i >= nx_int)

        if not in_pml:
            ops.element(
                "LadrunoUP",
                etag,
                n1,
                n2,
                n3,
                n4,
                1,
                "-thick",
                1.0,
                "-Kf",
                Kf,
                "-poro",
                poro,
                "-rhoF",
                rhoF,
                "-perm",
                perm,
                perm,
                "-alpha",
                biotAlpha,
                "-stab",
                "off",
                "-dynSeepage",
                "off",
            )
        elif mode == "box":
            ops.element(
                "LadrunoCPML",
                etag,
                n1,
                n2,
                n3,
                n4,
                "-mat",
                1,
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
                "-up",
                "-Kf",
                Kf,
                "-poro",
                poro,
                "-rhoF",
                rhoF,
                "-perm",
                perm,
                perm,
                "-alphaB",
                biotAlpha,
                "-stab",
                "off",
            )
        else:
            x0_vec = R_rot @ np.array([Lx_int, 0.0], dtype=np.float64)
            n_vec = R_rot @ np.array([1.0, 0.0], dtype=np.float64)
            ops.element(
                "LadrunoCPML",
                etag,
                n1,
                n2,
                n3,
                n4,
                "-mat",
                1,
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
                "-up",
                "-Kf",
                Kf,
                "-poro",
                poro,
                "-rhoF",
                rhoF,
                "-perm",
                perm,
                perm,
                "-alphaB",
                biotAlpha,
                "-stab",
                "off",
            )

    # Outer wall at x = L_total fixed in u and p
    for j in (0, 1):
        ops.fix(ntag(nx_elem, j), 1, 1, 1)

    # Enforce 1D waveguide periodicity across y' (tying top and bottom nodes in u_x, u_y, p)
    for i in range(nx_elem):
        ops.equalDOF(ntag(i, 0), ntag(i, 1), 1, 2, 3)

    ts_vals = [ricker_wavelet((step + 1) * dt, f0, t0) for step in range(n_steps + 5)]
    ts_times = [(step + 1) * dt for step in range(n_steps + 5)]
    ops.timeSeries("Path", 1, "-time", *ts_times, "-values", *ts_vals)
    ops.pattern("Plain", 1, 1)

    f_amp = 1.0e5 * (0.5 * Hy)
    if wave_type == "P1":
        fx = f_amp * ct
        fy = f_amp * st
        fp = 0.0
    else:
        fx = -f_amp * st
        fy = f_amp * ct
        fp = 0.0

    ops.load(ntag(0, 0), fx, fy, fp)
    ops.load(ntag(0, 1), fx, fy, fp)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    rec_i = int(round(6.0 / dx))
    rec_node = ntag(rec_i, 0)
    u_rec = np.zeros(n_steps, dtype=np.float64)
    p_rec = np.zeros(n_steps, dtype=np.float64)

    for step in range(n_steps):
        assert ops.analyze(1, dt) == 0, f"1D UP waveguide failed at step {step}"
        ux = ops.nodeDisp(rec_node, 1)
        uy = ops.nodeDisp(rec_node, 2)
        p_val = ops.nodeDisp(rec_node, 3)
        if wave_type == "P1":
            u_rec[step] = ux * ct + uy * st
        else:
            u_rec[step] = -ux * st + uy * ct
        p_rec[step] = p_val

    return u_rec, p_rec


def test_up_anchor1_saturated_soil_column_p1_and_s_wave_reflection() -> None:
    """Anchor UP-1 (1D Saturated Poroelastic Waveguide Reflection R_u <= 1%, R_p <= 1%
    + 37-deg Rotated `-plane` Invariance):
    (a) Undrained / high-mobility Fast P-wave (P1) excitation (`-box` mode): measure u_x(t)
        and p(t) at x_r = 6.0 m; assert R_u <= 0.01 (1%) and R_p <= 0.01 (1%), whereas a
        truncated column without `LadrunoCPML -up` has > 80% reflection.
    (b) Poroelastic Shear Wave (S-wave) excitation (`-box` mode): assert R_{u_y} <= 0.01 (1%).
    (c) `-plane` mode & 37-deg rotated poroelastic waveguide frame invariance: assert u_n(t)
        and p(t) match the axis-aligned column to <= 1e-11!
    """
    # (a) Fast P1 wave (-box vs extended reference vs truncated fixed)
    u_p1_ref, p_p1_ref = _run_1d_poroelastic_column("P1", "extended_ref")
    u_p1_box, p_p1_box = _run_1d_poroelastic_column("P1", "box")
    u_p1_trunc, p_p1_trunc = _run_1d_poroelastic_column("P1", "truncated_fixed")

    R_u_box = float(np.max(np.abs(u_p1_box - u_p1_ref)) / np.max(np.abs(u_p1_ref)))
    R_p_box = float(np.max(np.abs(p_p1_box - p_p1_ref)) / np.max(np.abs(p_p1_ref)))
    R_u_trunc = float(np.max(np.abs(u_p1_trunc - u_p1_ref)) / np.max(np.abs(u_p1_ref)))
    R_p_trunc = float(np.max(np.abs(p_p1_trunc - p_p1_ref)) / np.max(np.abs(p_p1_ref)))

    assert R_u_box <= 0.01, f"P1 solid disp reflection R_u={R_u_box:.4f} > 1%"
    assert R_p_box <= 0.01, f"P1 pore pressure reflection R_p={R_p_box:.4f} > 1%"
    assert R_u_trunc > 0.80, f"Truncated column R_u={R_u_trunc:.4f} <= 80%"
    assert R_p_trunc > 0.80, f"Truncated column R_p={R_p_trunc:.4f} <= 80%"

    # (b) Poroelastic S-wave (-box vs extended reference vs truncated fixed)
    u_s_ref, _ = _run_1d_poroelastic_column("S", "extended_ref")
    u_s_box, _ = _run_1d_poroelastic_column("S", "box")
    u_s_trunc, _ = _run_1d_poroelastic_column("S", "truncated_fixed")

    R_uy_box = float(np.max(np.abs(u_s_box - u_s_ref)) / np.max(np.abs(u_s_ref)))
    R_uy_trunc = float(np.max(np.abs(u_s_trunc - u_s_ref)) / np.max(np.abs(u_s_ref)))
    assert R_uy_box <= 0.01, f"S-wave reflection R_uy={R_uy_box:.4f} > 1%"
    assert R_uy_trunc > 0.80, f"Truncated S-wave R_uy={R_uy_trunc:.4f} <= 80%"

    # (c) `-plane` mode equivalence and 37-degree rotated frame invariance
    u_p1_plane, p_p1_plane = _run_1d_poroelastic_column("P1", "plane", theta_deg=0.0)
    u_p1_rot37, p_p1_rot37 = _run_1d_poroelastic_column("P1", "plane", theta_deg=37.0)

    assert rel_err(u_p1_plane, u_p1_box) <= 1.0e-13
    assert rel_err(p_p1_plane, p_p1_box) <= 1.0e-13
    assert rel_err(u_p1_rot37, u_p1_plane) <= 1.0e-11, (
        f"37-deg rotated u_n frame invariance error: {rel_err(u_p1_rot37, u_p1_plane):.3e}"
    )
    assert rel_err(p_p1_rot37, p_p1_plane) <= 1.0e-11, (
        f"37-deg rotated p frame invariance error: {rel_err(p_p1_rot37, p_p1_plane):.3e}"
    )


# ============================================================================ #
# 4. Anchor UP-2: 2D Cylindrical & 3D Spherical Radiation + Saturated Half-Space#
# ============================================================================ #


def _run_2d_cylindrical_poroelastic_q4(mode: str) -> dict[str, dict[str, np.ndarray]]:
    """Run 2D quadrant cylindrical poroelastic radiation model (`LadrunoUP` Q4 interior
    + `LadrunoCPML -up` face & 2D corner).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    dx = 0.5
    L_int = 6.0
    L_pml = 3.0
    if mode == "reference":
        L_tot = 21.5
    elif mode == "cpml":
        L_tot = L_int + L_pml
    elif mode == "truncated_fixed":
        L_tot = L_int
    else:
        raise ValueError(f"Unknown mode: {mode}")

    nx = int(round(L_tot / dx))
    ny = nx

    # Saturated porous medium: undrained Vp1 = 1000 m/s
    rho = 1000.0
    nu = 0.25
    E = 5.0e8
    Kf = 1.2e8
    poro = 0.30
    rhoF = 1000.0
    perm = 1.0e-6
    thick = 1.0

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int) -> int:
        return 1 + j * (nx + 1) + i

    for j in range(ny + 1):
        y = j * dx
        for i in range(nx + 1):
            x = i * dx
            n = nid(i, j)
            ops.node(n, x, y)
            outer_x = mode in ("cpml", "truncated_fixed") and i == nx
            outer_y = mode in ("cpml", "truncated_fixed") and j == ny
            fix_x = 1 if (i == 0 or outer_x) else 0
            fix_y = 1 if (j == 0 or outer_y) else 0
            fix_p = 1 if (outer_x or outer_y) else 0
            if fix_x or fix_y or fix_p:
                ops.fix(n, fix_x, fix_y, fix_p)

    etag = 1
    for j in range(ny):
        yc = (j + 0.5) * dx
        for i in range(nx):
            xc = (i + 0.5) * dx
            n1 = nid(i, j)
            n2 = nid(i + 1, j)
            n3 = nid(i + 1, j + 1)
            n4 = nid(i, j + 1)
            in_pml = (mode == "cpml") and (xc > L_int or yc > L_int)
            if in_pml:
                ops.element(
                    "LadrunoCPML",
                    etag,
                    n1,
                    n2,
                    n3,
                    n4,
                    "-mat",
                    1,
                    "-thick",
                    thick,
                    "-box",
                    0.0,
                    L_int,
                    0.0,
                    L_int,
                    "-L",
                    L_pml,
                    "-m",
                    2.0,
                    "-R0",
                    1.0e-4,
                    "-kappaMax",
                    3.0,
                    "-f0",
                    55.0,
                    "-lumped",
                    "-up",
                    "-Kf",
                    Kf,
                    "-poro",
                    poro,
                    "-rhoF",
                    rhoF,
                    "-perm",
                    perm,
                    perm,
                    "-stab",
                    "auto",
                    0.25,
                )
            else:
                ops.element(
                    "LadrunoUP",
                    etag,
                    n1,
                    n2,
                    n3,
                    n4,
                    1,
                    "-thick",
                    thick,
                    "-Kf",
                    Kf,
                    "-poro",
                    poro,
                    "-rhoF",
                    rhoF,
                    "-perm",
                    perm,
                    perm,
                    "-stab",
                    "auto",
                    0.25,
                    "-lumped",
                    "-dynSeepage",
                    "off",
                )
            etag += 1

    dt = 0.0002
    n_steps = 165
    f0 = 55.0
    t0 = 0.014
    ricker = ricker_series(dt, n_steps, f0, t0, amp=1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)

    r_src = 1.2
    for j in range(4):
        for i in range(4):
            x = i * dx
            y = j * dx
            r = math.hypot(x, y)
            if 1.0e-9 < r <= r_src:
                w = math.exp(-2.5 * (r / r_src) ** 2)
                ops.load(nid(i, j), w * (x / r), w * (y / r), 0.0)

    rec_nodes = {
        "deg0": nid(10, 1),    # (5.0, 0.5) normal incidence
        "deg45": nid(10, 10),  # (5.0, 5.0) 45-degree 2D corner incidence
    }

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient", "-noWarnings")

    traces = {
        k: {
            "u": np.zeros((n_steps, 2), dtype=np.float64),
            "p": np.zeros(n_steps, dtype=np.float64),
        }
        for k in rec_nodes
    }
    for step in range(n_steps):
        assert ops.analyze(1, dt) == 0, f"2D UP Q4 failed at step {step} in {mode}"
        for k, node_tag in rec_nodes.items():
            traces[k]["u"][step, 0] = ops.nodeDisp(node_tag, 1)
            traces[k]["u"][step, 1] = ops.nodeDisp(node_tag, 2)
            traces[k]["p"][step] = ops.nodeDisp(node_tag, 3)

    return traces


def _run_3d_spherical_poroelastic_h8(mode: str) -> dict[str, dict[str, np.ndarray]]:
    """Run 3D octant spherical poroelastic wave radiation model (`LadrunoUP` H8 interior
    + `LadrunoCPML -up` 3D face, edge, and (x, y, z) corner).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)

    dx = 0.5
    L_int = 1.5
    L_pml = 2.5
    if mode == "reference":
        L_tot = 4.5
    elif mode == "cpml":
        L_tot = L_int + L_pml
    elif mode == "truncated_fixed":
        L_tot = L_int
    else:
        raise ValueError(f"Unknown mode: {mode}")

    n_side = int(round(L_tot / dx))
    i_int = int(round(L_int / dx))

    rho = 1000.0
    nu = 0.25
    E = 5.0e8
    Kf = 1.2e8
    poro = 0.30
    rhoF = 1000.0
    perm = 1.0e-6

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int, k: int) -> int:
        return 1 + (k * (n_side + 1) + j) * (n_side + 1) + i

    for k in range(n_side + 1):
        for j in range(n_side + 1):
            for i in range(n_side + 1):
                x, y, z = i * dx, j * dx, k * dx
                n = nid(i, j, k)
                ops.node(n, x, y, z)
                outer_x = mode in ("cpml", "truncated_fixed") and i == n_side
                outer_y = mode in ("cpml", "truncated_fixed") and j == n_side
                outer_z = mode in ("cpml", "truncated_fixed") and k == n_side
                fx = 1 if (i == 0 or outer_x) else 0
                fy = 1 if (j == 0 or outer_y) else 0
                fz = 1 if (k == 0 or outer_z) else 0
                fp = 1 if (outer_x or outer_y or outer_z) else 0
                if fx or fy or fz or fp:
                    ops.fix(n, fx, fy, fz, fp)

    etag = 1
    for k in range(n_side):
        for j in range(n_side):
            for i in range(n_side):
                in_pml = (mode == "cpml") and (i >= i_int or j >= i_int or k >= i_int)
                nodes = [
                    nid(i, j, k),
                    nid(i + 1, j, k),
                    nid(i + 1, j + 1, k),
                    nid(i, j + 1, k),
                    nid(i, j, k + 1),
                    nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1),
                    nid(i, j + 1, k + 1),
                ]
                if in_pml:
                    ops.element(
                        "LadrunoCPML",
                        etag,
                        *nodes,
                        "-mat",
                        1,
                        "-box",
                        0.0,
                        L_int,
                        0.0,
                        L_int,
                        0.0,
                        L_int,
                        "-L",
                        L_pml,
                        "-m",
                        2.5,
                        "-R0",
                        1.5e-3,
                        "-kappaMax",
                        1.0,
                        "-alphaMax",
                        15.0,
                        "-up",
                        "-Kf",
                        Kf,
                        "-poro",
                        poro,
                        "-rhoF",
                        rhoF,
                        "-perm",
                        perm,
                        perm,
                        perm,
                        "-stab",
                        "auto",
                        0.25,
                    )
                else:
                    ops.element(
                        "LadrunoUP",
                        etag,
                        *nodes,
                        1,
                        "-Kf",
                        Kf,
                        "-poro",
                        poro,
                        "-rhoF",
                        rhoF,
                        "-perm",
                        perm,
                        perm,
                        perm,
                        "-stab",
                        "auto",
                        0.25,
                        "-dynSeepage",
                        "off",
                    )
                etag += 1

    dt = 0.00015
    n_steps = 32
    f0 = 115.0
    t0 = 0.0022
    ricker = ricker_series(dt, n_steps, f0, t0, amp=1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)

    for k in range(2):
        for j in range(2):
            for i in range(2):
                x, y, z = i * dx, j * dx, k * dx
                r = math.sqrt(x * x + y * y + z * z)
                if r > 1.0e-9:
                    ops.load(nid(i, j, k), x / r, y / r, z / r, 0.0)

    rec_nodes = {
        "axis": nid(2, 1, 1),    # (1.0, 0.5, 0.5) near coordinate face
        "corner": nid(2, 2, 2),  # (1.0, 1.0, 1.0) into 3D (x, y, z) corner
    }

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient", "-noWarnings")

    traces = {
        k: {
            "u": np.zeros((n_steps, 3), dtype=np.float64),
            "p": np.zeros(n_steps, dtype=np.float64),
        }
        for k in rec_nodes
    }
    for step in range(n_steps):
        assert ops.analyze(1, dt) == 0, f"3D UP H8 failed at step {step} in {mode}"
        for k, node_tag in rec_nodes.items():
            traces[k]["u"][step, 0] = ops.nodeDisp(node_tag, 1)
            traces[k]["u"][step, 1] = ops.nodeDisp(node_tag, 2)
            traces[k]["u"][step, 2] = ops.nodeDisp(node_tag, 3)
            traces[k]["p"][step] = ops.nodeDisp(node_tag, 4)

    return traces


def _run_poroelastic_lamb_halfspace(mode: str) -> dict[str, np.ndarray]:
    """Run 2D saturated poroelastic half-space under vertical Ricker surface load
    (`LadrunoUP` interior [0, 8m] x [-6m, 0m] + `LadrunoCPML -up` lateral & bottom
    absorber, drained free surface p = 0 at y = 0).
    Returns surface displacement `u_surf` at (4.0m, 0.0m) and near-surface pore
    pressure `p_sub` at (4.0m, -0.5m).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    dx = 0.5
    Lx = 8.0
    Ly = 6.0
    L_pml = 4.0

    if mode == "reference":
        Lx_tot = 22.0
        Ly_tot = 16.0
    elif mode == "cpml":
        Lx_tot = Lx + L_pml
        Ly_tot = Ly + L_pml
    else:
        raise ValueError(f"Unknown mode: {mode}")

    nx = int(round(Lx_tot / dx))
    ny = int(round(Ly_tot / dx))
    nx_int = int(round(Lx / dx))
    ny_int = int(round(Ly / dx))

    rho = 1000.0
    nu = 0.25
    E = 5.0e8
    Kf = 1.2e8
    poro = 0.30
    rhoF = 1000.0
    perm = 1.0e-7
    thick = 1.0

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int) -> int:
        return 1 + j * (nx + 1) + i

    for j in range(ny + 1):
        y = -j * dx
        for i in range(nx + 1):
            x = i * dx
            n = nid(i, j)
            ops.node(n, x, y)
            outer_wall = i == nx or j == ny
            fix_x = 1 if (i == 0 or outer_wall) else 0
            fix_y = 1 if outer_wall else 0
            # Drained free surface p = 0 at y = 0 (j == 0), plus outer boundary
            fix_p = 1 if (j == 0 or outer_wall) else 0
            if fix_x or fix_y or fix_p:
                ops.fix(n, fix_x, fix_y, fix_p)

    etag = 1
    for j in range(ny):
        for i in range(nx):
            n1 = nid(i, j + 1)
            n2 = nid(i + 1, j + 1)
            n3 = nid(i + 1, j)
            n4 = nid(i, j)
            in_pml = (mode == "cpml") and (i >= nx_int or j >= ny_int)
            if in_pml:
                ops.element(
                    "LadrunoCPML",
                    etag,
                    n1,
                    n2,
                    n3,
                    n4,
                    "-mat",
                    1,
                    "-thick",
                    thick,
                    "-box",
                    -100.0,
                    Lx,
                    -Ly,
                    100.0,
                    "-L",
                    L_pml,
                    "-m",
                    2.0,
                    "-R0",
                    1.0e-5,
                    "-kappaMax",
                    1.5,
                    "-f0",
                    55.0,
                    "-lumped",
                    "-up",
                    "-Kf",
                    Kf,
                    "-poro",
                    poro,
                    "-rhoF",
                    rhoF,
                    "-perm",
                    perm,
                    perm,
                    "-stab",
                    "auto",
                    0.25,
                )
            else:
                ops.element(
                    "LadrunoUP",
                    etag,
                    n1,
                    n2,
                    n3,
                    n4,
                    1,
                    "-thick",
                    thick,
                    "-Kf",
                    Kf,
                    "-poro",
                    poro,
                    "-rhoF",
                    rhoF,
                    "-perm",
                    perm,
                    perm,
                    "-stab",
                    "auto",
                    0.25,
                    "-lumped",
                    "-dynSeepage",
                    "off",
                )
            etag += 1

    dt = 0.0002
    n_steps = 155
    f0 = 55.0
    t0 = 0.014
    ricker = ricker_series(dt, n_steps, f0, t0, amp=-1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)
    ops.load(nid(0, 0), 0.0, 0.5, 0.0)
    ops.load(nid(1, 0), 0.0, 0.8, 0.0)
    ops.load(nid(2, 0), 0.0, 0.4, 0.0)
    ops.load(nid(3, 0), 0.0, 0.1, 0.0)

    rec_surf = nid(8, 0)  # (4.0m, 0.0m) on drained free surface
    rec_sub = nid(8, 1)   # (4.0m, -0.5m) sub-surface pore pressure receiver

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient", "-noWarnings")

    u_surf = np.zeros((n_steps, 2), dtype=np.float64)
    p_sub = np.zeros(n_steps, dtype=np.float64)
    for step in range(n_steps):
        assert ops.analyze(1, dt) == 0, f"Poroelastic Lamb failed at step {step}"
        u_surf[step, 0] = ops.nodeDisp(rec_surf, 1)
        u_surf[step, 1] = ops.nodeDisp(rec_surf, 2)
        p_sub[step] = ops.nodeDisp(rec_sub, 3)

    return {"u_surf": u_surf, "p_sub": p_sub}


def test_up_anchor2_2d_3d_poroelastic_radiation_and_lamb_halfspace() -> None:
    """Anchor UP-2 (2D Cylindrical & 3D Spherical Corner Poroelastic Radiation +
    Saturated Half-Space Surface Wave):
    (a) 2D Cylindrical Poroelastic Wave Radiation (`LadrunoUP` Q4 interior +
        `LadrunoCPML -up` face & 2D corner): relative L2 error <= 3.5% on both u and p
        at normal (0 deg) and corner (45 deg) receivers, and > 15x reduction vs
        unabsorbed boundary.
    (b) 3D Spherical Poroelastic Wave Radiation (`LadrunoUP` H8 interior +
        `LadrunoCPML -up` 3D face, edge & (x,y,z) corner): relative error <= 4.5% on
        both u and p and > 10x reduction vs unabsorbed boundary.
    (c) Saturated Poroelastic Half-Space Surface Point Load (Rayleigh Wave +
        Pore-Pressure Diffusion): surface displacement and near-surface pore-pressure
        waveforms match the extended saturated half-space to <= 3.5%!
    """
    # --- (a) 2D Cylindrical Poroelastic Radiation ---
    ref_2d = _run_2d_cylindrical_poroelastic_q4("reference")
    cpml_2d = _run_2d_cylindrical_poroelastic_q4("cpml")
    trunc_2d = _run_2d_cylindrical_poroelastic_q4("truncated_fixed")

    for angle_key in ("deg0", "deg45"):
        for field in ("u", "p"):
            l2_cpml, _ = rel_l2_and_peak_error(
                cpml_2d[angle_key][field], ref_2d[angle_key][field]
            )
            l2_trunc, _ = rel_l2_and_peak_error(
                trunc_2d[angle_key][field], ref_2d[angle_key][field]
            )
            reduction = l2_trunc / max(l2_cpml, 1.0e-12)
            assert l2_cpml <= 0.035, (
                f"2D UP Q4 {angle_key} {field} L2 error {l2_cpml:.4f} > 3.5%"
            )
            assert reduction >= 15.0, (
                f"2D UP Q4 {angle_key} {field} reduction {reduction:.1f}x < 15x"
            )

    # --- (b) 3D Spherical Poroelastic Radiation ---
    ref_3d = _run_3d_spherical_poroelastic_h8("reference")
    cpml_3d = _run_3d_spherical_poroelastic_h8("cpml")
    trunc_3d = _run_3d_spherical_poroelastic_h8("truncated_fixed")

    for rec_key in ("axis", "corner"):
        for field in ("u", "p"):
            l2_cpml, _ = rel_l2_and_peak_error(
                cpml_3d[rec_key][field], ref_3d[rec_key][field]
            )
            l2_trunc, _ = rel_l2_and_peak_error(
                trunc_3d[rec_key][field], ref_3d[rec_key][field]
            )
            reduction = l2_trunc / max(l2_cpml, 1.0e-12)
            assert l2_cpml <= 0.045, (
                f"3D UP H8 {rec_key} {field} L2 error {l2_cpml:.4f} > 4.5%"
            )
            assert reduction >= 10.0, (
                f"3D UP H8 {rec_key} {field} reduction {reduction:.1f}x < 10x"
            )

    # --- (c) Saturated Poroelastic Half-Space (Rayleigh Wave + Pore-Pressure) ---
    lamb_ref = _run_poroelastic_lamb_halfspace("reference")
    lamb_cpml = _run_poroelastic_lamb_halfspace("cpml")

    l2_u, peak_u = rel_l2_and_peak_error(lamb_cpml["u_surf"], lamb_ref["u_surf"])
    l2_p, peak_p = rel_l2_and_peak_error(lamb_cpml["p_sub"], lamb_ref["p_sub"])
    assert l2_u <= 0.035 and peak_u <= 0.035, (
        f"Poroelastic Lamb u_surf error L2={l2_u:.4f}, peak={peak_u:.4f} > 3.5%"
    )
    assert l2_p <= 0.035 and peak_p <= 0.035, (
        f"Poroelastic Lamb p_sub error L2={l2_p:.4f}, peak={peak_p:.4f} > 3.5%"
    )


# ============================================================================ #
# 5. Anchor UP-3: Gravity Seepage Staging, Rayleigh, Lyapunov, Checkpoint &    #
#                 Parser Fail-Closed Battery                                   #
# ============================================================================ #


def _check_up_gravity_seepage_staging() -> None:
    """(a) Hydrostatic Gravity + Steady Seepage Stage `0 -> 1` Zero-Jump Verification."""
    E, nu, rho = 4.0e7, 0.30, 2000.0
    Kf, poro, rhoF = 2.2e9, 0.35, 1000.0
    perm = 1.0e-5
    g = -9.81

    for transition_token in ("stage", "updateMaterialStage"):
        for load_mode in ("direct_body", "selfWeight"):
            b_val = (0.0, g) if load_mode == "direct_body" else (0.0, 1.0)

            def build_saturated_box(use_cpml: bool) -> list[int]:
                ops.wipe()
                ops.model("basic", "-ndm", 2, "-ndf", 3)
                # 3x2 grid on [0, 3] x [0, 2]
                for j in range(3):
                    for i in range(4):
                        nid = j * 4 + i + 1
                        ops.node(nid, float(i), float(j))
                        if j == 0:
                            # Impervious fixed base (u_x = u_y = 0, p free)
                            ops.fix(nid, 1, 1, 0)
                        elif j == 2:
                            # Drained top surface (p = 0), lateral rollers at i=0,3
                            ops.fix(nid, 1 if i in (0, 3) else 0, 0, 1)
                        elif i in (0, 3):
                            ops.fix(nid, 1, 0, 0)

                ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
                ele_tag = 1
                cpml_tags: list[int] = []
                all_tags: list[int] = []
                for j in range(2):
                    for i in range(3):
                        n1 = j * 4 + i + 1
                        n2 = n1 + 1
                        n3 = n2 + 4
                        n4 = n1 + 4
                        all_tags.append(ele_tag)
                        if not use_cpml or (i == 1 and j == 1):
                            ops.element(
                                "LadrunoUP",
                                ele_tag,
                                n1,
                                n2,
                                n3,
                                n4,
                                1,
                                "-Kf",
                                Kf,
                                "-poro",
                                poro,
                                "-rhoF",
                                rhoF,
                                "-perm",
                                perm,
                                perm,
                                "-body",
                                b_val[0],
                                b_val[1],
                                "-fluidBody",
                                b_val[0],
                                b_val[1],
                                "-dynSeepage",
                                "off",
                            )
                        else:
                            ops.element(
                                "LadrunoCPML",
                                ele_tag,
                                n1,
                                n2,
                                n3,
                                n4,
                                "-mat",
                                1,
                                "-box",
                                1.0,
                                2.0,
                                1.0,
                                2.0,
                                "-L",
                                1.0,
                                1.0,
                                "-stage",
                                0,
                                "-body",
                                b_val[0],
                                b_val[1],
                                "-up",
                                "-Kf",
                                Kf,
                                "-poro",
                                poro,
                                "-rhoF",
                                rhoF,
                                "-perm",
                                perm,
                                perm,
                                "-fluidBody",
                                b_val[0],
                                b_val[1],
                            )
                            cpml_tags.append(ele_tag)
                        ele_tag += 1

                if load_mode == "selfWeight":
                    ops.timeSeries("Linear", 1)
                    ops.pattern("Plain", 1, 1)
                    ops.eleLoad("-ele", *all_tags, "-type", "-selfWeight", 0.0, g)

                ops.constraints("Plain")
                ops.numberer("Plain")
                ops.system("UmfPack")
                ops.test("NormDispIncr", 1.0e-14, 10, 0)
                ops.algorithm("Linear")
                ops.integrator("LoadControl", 1.0)
                ops.analysis("Static", "-noWarnings")
                assert ops.analyze(1) == 0
                ops.loadConst("-time", 0.0)
                return cpml_tags

            build_saturated_box(use_cpml=False)
            up_ref = {nd: list(ops.nodeDisp(nd)) for nd in ops.getNodeTags()}

            cpml_tags = build_saturated_box(use_cpml=True)
            up_0 = {nd: list(ops.nodeDisp(nd)) for nd in ops.getNodeTags()}

            # Verify non-trivial hydrostatic pore pressure p = rhoF * |g| * depth = 19620 Pa at base
            assert abs(up_0[1][2] - 2.0 * rhoF * abs(g)) <= 1.0e-9
            for nd in up_ref:
                for d in range(2):
                    assert abs(up_0[nd][d] - up_ref[nd][d]) <= 1.0e-13
                assert abs(up_0[nd][2] - up_ref[nd][2]) <= 1.0e-10

            # Transition all LadrunoCPML -up elements from stage 0 -> stage 1
            for idx, ctag in enumerate(cpml_tags, start=1):
                assert ops.eleResponse(ctag, "stage") == [0]
                if transition_token == "stage":
                    ops.parameter(100 + idx, "element", ctag, "stage")
                else:
                    ops.parameter(
                        100 + idx, "element", ctag, "updateMaterialStage", str(ctag)
                    )
                ops.updateParameter(100 + idx, 1)
                assert ops.eleResponse(ctag, "stage") == [1]

            # Switch to TransientAnalysis and step 20 steps with no dynamic load
            ops.wipeAnalysis()
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("UmfPack")
            ops.test("NormDispIncr", 1.0e-14, 10, 0)
            ops.algorithm("Linear")
            ops.integrator("Newmark", 0.5, 0.25)
            ops.analysis("Transient", "-noWarnings")

            dt = 0.002
            for step in range(1, 21):
                assert ops.analyze(1, dt) == 0
                for nd in up_0:
                    cur = ops.nodeDisp(nd)
                    for d in range(2):
                        assert abs(cur[d] - up_0[nd][d]) <= 1.0e-13, (
                            f"Solid disp jump at step {step} node {nd}: "
                            f"|u - u0|={abs(cur[d] - up_0[nd][d]):.3e}"
                        )
                    assert abs(cur[2] - up_0[nd][2]) <= 1.0e-10, (
                        f"Pore pressure jump at step {step} node {nd}: "
                        f"|p - p0|={abs(cur[2] - up_0[nd][2]):.3e}"
                    )


def _check_up_rayleigh_solid_only_isolation() -> None:
    """(b) Solid-Only Rayleigh Damping Isolation (`alphaM, betaK, betaK0, betaKc`)."""
    coords = ref.make_reference_coords("Q4")
    u_idx = [0, 1, 3, 4, 6, 7, 9, 10]
    p_idx = [2, 5, 8, 11]

    def run_rayleigh_model(use_rayleigh: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        for i in range(4):
            ops.node(i + 1, float(coords[i, 0]), float(coords[i, 1]))
        ops.fix(1, 1, 1, 0)
        ops.fix(4, 1, 1, 0)

        ops.element(
            "LadrunoCPML",
            1,
            1,
            2,
            3,
            4,
            "-E",
            5.0e7,
            "-nu",
            0.25,
            "-rho",
            2000.0,
            "-box",
            0.0,
            1.8,
            0.0,
            1.8,
            "-L",
            2.0,
            2.0,
            "-up",
            "-Kf",
            2.2e9,
            "-poro",
            0.30,
            "-rhoF",
            1000.0,
            "-perm",
            2.0e-6,
            2.0e-6,
            "-alphaB",
            0.90,
            "-Ks",
            3.6e10,
            "-stab",
            "auto",
            0.25,
        )
        if use_rayleigh:
            ops.rayleigh(0.2, 0.002, 0.003, 0.004)

        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 1500.0, -700.0, 2.0e-2)
        ops.load(3, 900.0, 500.0, -1.5e-2)

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-14, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient", "-noWarnings")
        for _ in range(3):
            assert ops.analyze(1, 0.002) == 0

        C_mat = np.asarray(ops.eleResponse(1, "damp"), dtype=np.float64).reshape(12, 12)
        K_eff = np.asarray(ops.eleResponse(1, "stiff"), dtype=np.float64).reshape(12, 12)
        K_0 = np.asarray(ops.eleResponse(1, "stiffInitial"), dtype=np.float64).reshape(
            12, 12
        )
        M_eff = np.asarray(ops.eleResponse(1, "mass"), dtype=np.float64).reshape(12, 12)
        f_damp = np.asarray(ops.eleResponse(1, "dampingForce"), dtype=np.float64)
        return C_mat, K_eff, K_0, M_eff, f_damp

    C_no_ray, _, _, _, f_damp_no = run_rayleigh_model(False)
    C_ray, K_eff, K_0, M_eff, f_damp_ray = run_rayleigh_model(True)

    assert float(np.linalg.norm(f_damp_no)) == 0.0
    assert float(np.linalg.norm(f_damp_ray[u_idx])) > 0.0
    # Strictly zero on all p-slots!
    assert float(np.max(np.abs(f_damp_ray[p_idx]))) == 0.0

    # (u, p), (p, u), and (p, p) blocks of C_mat are completely unchanged to <= 1e-14!
    for r_ix, c_ix in ((u_idx, p_idx), (p_idx, u_idx), (p_idx, p_idx)):
        diff = float(np.max(np.abs(C_ray[np.ix_(r_ix, c_ix)] - C_no_ray[np.ix_(r_ix, c_ix)])))
        scale = max(1.0e-30, float(np.max(np.abs(C_no_ray[np.ix_(r_ix, c_ix)]))))
        assert diff / scale <= 1.0e-14

    # And (u, u) block matches C_cpml_uu + alphaM*M_uu + (betaK + betaKc)*K_eff_uu + betaK0*K0_uu
    C_uu_expected = (
        C_no_ray[np.ix_(u_idx, u_idx)]
        + 0.2 * M_eff[np.ix_(u_idx, u_idx)]
        + (0.002 + 0.004) * K_eff[np.ix_(u_idx, u_idx)]
        + 0.003 * K_0[np.ix_(u_idx, u_idx)]
    )
    assert rel_err(C_ray[np.ix_(u_idx, u_idx)], C_uu_expected) <= 1.0e-12


def _check_up_2500_step_lyapunov_stability() -> None:
    """(c) 2500-Step Late-Time Poroelastic Lyapunov Stability (2D & 3D Corner)."""
    N_total = 2500
    N_pulse = 100
    dt = 0.0005
    f0 = 15.0
    t0 = 40 * dt

    times = [k * dt for k in range(N_pulse + 2)]
    vals = [
        0.0
        if (k == 0 or k >= N_pulse)
        else 1.0e5 * (math.sin(math.pi * k / float(N_pulse)) ** 2) * ricker_wavelet(t, f0, t0)
        for k, t in enumerate(times)
    ]

    # --- 1. 2D Corner `LadrunoCPML -up` model (2500 steps) ---
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    coords_2d = [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]
    for idx, (x, y) in enumerate(coords_2d, start=1):
        ops.node(idx, x, y)
        if idx != 1:
            ops.fix(idx, 1, 1, 1)

    ops.element(
        "LadrunoCPML",
        1,
        1,
        2,
        3,
        4,
        "-E",
        5.0e7,
        "-nu",
        0.25,
        "-rho",
        2000.0,
        "-box",
        0.0,
        1.0,
        0.0,
        1.0,
        "-L",
        1.0,
        1.0,
        "-m",
        2.0,
        "-R0",
        1.0e-4,
        "-kappaMax",
        4.0,
        "-alphaMax",
        math.pi * 15.0,
        "-up",
        "-Kf",
        2.2e8,
        "-poro",
        0.30,
        "-rhoF",
        1000.0,
        "-perm",
        1.0e-5,
        1.0e-5,
        "-stab",
        "auto",
        0.25,
    )

    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0, -0.7, 1.0e-4)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    up_norms_2d = np.zeros(N_total, dtype=np.float64)
    for step in range(N_total):
        assert ops.analyze(1, dt) == 0
        d1 = ops.nodeDisp(1)
        up_norms_2d[step] = math.sqrt(d1[0] ** 2 + d1[1] ** 2 + (d1[2] * 1.0e-7) ** 2)

    peak_2d = float(np.max(up_norms_2d))
    assert peak_2d > 0.0
    assert up_norms_2d[-1] < 1.0e-4 * peak_2d, (
        f"2D UP corner final norm {up_norms_2d[-1]:.3e} >= 1e-4 * peak ({peak_2d:.3e})"
    )

    # --- 2. 3D H8 Corner `LadrunoCPML -up` model (2500 steps) ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    coords_3d = [
        (1.0, 1.0, 1.0),
        (2.0, 1.0, 1.0),
        (2.0, 2.0, 1.0),
        (1.0, 2.0, 1.0),
        (1.0, 1.0, 2.0),
        (2.0, 1.0, 2.0),
        (2.0, 2.0, 2.0),
        (1.0, 2.0, 2.0),
    ]
    for idx, (x, y, z) in enumerate(coords_3d, start=1):
        ops.node(idx, x, y, z)
        if idx != 1:
            ops.fix(idx, 1, 1, 1, 1)

    ops.element(
        "LadrunoCPML",
        1,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        "-E",
        5.0e7,
        "-nu",
        0.25,
        "-rho",
        2000.0,
        "-box",
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        "-L",
        1.0,
        1.0,
        1.0,
        "-m",
        2.0,
        "-R0",
        1.0e-4,
        "-kappaMax",
        4.0,
        "-alphaMax",
        math.pi * 15.0,
        "-up",
        "-Kf",
        2.2e8,
        "-poro",
        0.30,
        "-rhoF",
        1000.0,
        "-perm",
        1.0e-5,
        1.0e-5,
        1.0e-5,
        "-stab",
        "auto",
        0.25,
    )

    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0, -0.7, 0.5, 1.0e-4)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")

    up_norms_3d = np.zeros(N_total, dtype=np.float64)
    for step in range(N_total):
        assert ops.analyze(1, dt) == 0
        d1 = ops.nodeDisp(1)
        up_norms_3d[step] = math.sqrt(
            d1[0] ** 2 + d1[1] ** 2 + d1[2] ** 2 + (d1[3] * 1.0e-7) ** 2
        )

    peak_3d = float(np.max(up_norms_3d))
    assert peak_3d > 0.0
    assert up_norms_3d[-1] < 1.0e-4 * peak_3d, (
        f"3D H8 UP corner final norm {up_norms_3d[-1]:.3e} >= 1e-4 * peak ({peak_3d:.3e})"
    )


def _check_up_filedatastore_checkpoint_roundtrip(tmp_path: Path) -> None:
    """(d) `FileDatastore` `sendSelf` / `recvSelf` Checkpoint Round-Trip across
    Equal-Order (`Q4, H8`) and Taylor-Hood (`BT6, BTET10`), testing both live
    restore and fresh-domain restore.
    """
    configs = [
        ("Q4", 0),
        ("BT6", 1),
        ("H8", 0),
        ("BTET10", 1),
    ]
    dt = 0.001

    for topo, p_order in configs:
        ndm, n_nodes, _, _ = ref.TOPOLOGIES[topo]
        coords = ref.make_reference_coords(topo)

        def setup_model_and_transient() -> list[int]:
            ops.wipe()
            ntags = create_up_nodes(topo, coords, p_order)
            # Fix vertex node 1 so rigid modes are restrained while other vertices/midsides are free
            if ndm == 2:
                ops.fix(1, 1, 1, 1)
            else:
                ops.fix(1, 1, 1, 1, 1)

            box_bounds = (
                [-2.0, 1.8, -2.0, 1.8]
                if ndm == 2
                else [-2.0, 1.8, -2.0, 1.8, -2.0, 1.8]
            )
            kbar = (2.0e-6, 1.0e-6) if ndm == 2 else (2.0e-6, 1.0e-6, 8.0e-7)
            body = (0.0, -9.81) if ndm == 2 else (0.0, 0.0, -9.81)

            ops.element(
                "LadrunoCPML",
                1,
                *ntags,
                "-shape",
                topo,
                "-E",
                6.0e7,
                "-nu",
                0.28,
                "-rho",
                2100.0,
                "-box",
                *box_bounds,
                "-L",
                2.0,
                "-body",
                *body,
                "-up",
                "-Kf",
                2.2e9,
                "-poro",
                0.30,
                "-rhoF",
                1000.0,
                "-perm",
                *kbar,
                "-alphaB",
                0.86,
                "-Ks",
                3.6e10,
                "-pOrder",
                "linear" if p_order == 1 else "equal",
                "-formulation",
                "bbar" if topo in ("Q4", "H8") else "std",
                "-fluidBody",
                *body,
            )
            ops.rayleigh(0.05, 1.0e-4, 1.0e-4, 2.0e-4)

            ops.timeSeries("Trig", 1, 0.0, 1.0, 0.02)
            ops.pattern("Plain", 1, 1)
            if ndm == 2:
                ops.load(2, 2000.0, -1200.0, 1.5e-2)
                ops.load(3, -900.0, 1600.0, -1.0e-2)
            else:
                ops.load(2, 2500.0, -1200.0, 1800.0, 1.5e-2)
                ops.load(3, -1100.0, 1700.0, -900.0, -1.0e-2)

            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("FullGeneral")
            ops.test("NormDispIncr", 1.0e-14, 10, 0)
            ops.algorithm("Linear")
            ops.integrator("Newmark", 0.5, 0.25)
            ops.analysis("Transient", "-noWarnings")
            return ntags

        for restore_mode in ("live", "fresh"):
            db_prefix = str(tmp_path / f"cpml_up_{topo}_{restore_mode}")
            ntags = setup_model_and_transient()

            for _ in range(4):
                assert ops.analyze(1, dt) == 0

            ops.database("File", db_prefix)
            assert ops.save(1) in (0, None)

            ref_after = []
            for _ in range(4):
                assert ops.analyze(1, dt) == 0
                snap = np.concatenate(
                    [
                        read_up_nodal_vector(ntags, topo, p_order, "disp"),
                        read_up_nodal_vector(ntags, topo, p_order, "vel"),
                        np.asarray(ops.eleResponse(1, "force"), dtype=np.float64),
                        np.asarray(ops.eleResponse(1, "basicForce"), dtype=np.float64),
                    ]
                )
                ref_after.append(snap)

            if restore_mode == "fresh":
                ntags = setup_model_and_transient()
                ops.database("File", db_prefix)

            assert ops.restore(1) in (0, None)
            ops.domainChange()

            got_after = []
            for _ in range(4):
                assert ops.analyze(1, dt) == 0
                snap = np.concatenate(
                    [
                        read_up_nodal_vector(ntags, topo, p_order, "disp"),
                        read_up_nodal_vector(ntags, topo, p_order, "vel"),
                        np.asarray(ops.eleResponse(1, "force"), dtype=np.float64),
                        np.asarray(ops.eleResponse(1, "basicForce"), dtype=np.float64),
                    ]
                )
                got_after.append(snap)

            max_diff = float(np.max(np.abs(np.array(got_after) - np.array(ref_after))))
            assert max_diff <= 1.0e-14, (
                f"Checkpoint {topo} ({restore_mode}) drift {max_diff:.3e} > 1e-14"
            )


def _check_up_parser_fail_closed_battery(capfd: Any) -> None:
    """(e) Fail-Closed `-up` Parser Refusal Battery."""
    base_2d = [
        1,
        1,
        2,
        3,
        4,
        "-E",
        5.0e7,
        "-nu",
        0.25,
        "-rho",
        2000.0,
        "-box",
        0.0,
        1.0,
        0.0,
        1.0,
        "-L",
        1.0,
    ]
    valid_up = ["-up", "-Kf", 2.2e9, "-poro", 0.30, "-rhoF", 1000.0, "-perm", 1.0e-6]

    bad_2d_cases = [
        # 1. Biot flag (-Kf, -poro, -rhoF, -perm, -alphaB, -Ks, -stab, -pOrder, -formulation) without -up
        [*base_2d, "-Kf", 2.2e9],
        [*base_2d, "-poro", 0.30],
        [*base_2d, "-rhoF", 1000.0],
        [*base_2d, "-perm", 1.0e-6],
        [*base_2d, "-alphaB", 0.9],
        [*base_2d, "-Ks", 3.6e10],
        [*base_2d, "-stab", "auto"],
        [*base_2d, "-pOrder", "equal"],
        [*base_2d, "-formulation", "bbar"],
        # 2. -up missing -Kf, -poro, -rhoF, or -perm
        [*base_2d, "-up", "-poro", 0.30, "-rhoF", 1000.0, "-perm", 1.0e-6],
        [*base_2d, "-up", "-Kf", 2.2e9, "-rhoF", 1000.0, "-perm", 1.0e-6],
        [*base_2d, "-up", "-Kf", 2.2e9, "-poro", 0.30, "-perm", 1.0e-6],
        [*base_2d, "-up", "-Kf", 2.2e9, "-poro", 0.30, "-rhoF", 1000.0],
        # 3. Both -perm and -permH (or -perm and -gammaW)
        [*base_2d, *valid_up, "-permH", 1.0e-3, "-gammaW", 9810.0],
        [*base_2d, *valid_up, "-gammaW", 9810.0],
        # 4. -permH without -gammaW (and -gammaW without -permH)
        [*base_2d, "-up", "-Kf", 2.2e9, "-poro", 0.30, "-rhoF", 1000.0, "-permH", 1.0e-3],
        [*base_2d, "-up", "-Kf", 2.2e9, "-poro", 0.30, "-rhoF", 1000.0, "-gammaW", 9810.0],
        # 5. poro < 0 or poro >= 1
        [*base_2d, "-up", "-Kf", 2.2e9, "-poro", -0.1, "-rhoF", 1000.0, "-perm", 1.0e-6],
        [*base_2d, "-up", "-Kf", 2.2e9, "-poro", 1.0, "-rhoF", 1000.0, "-perm", 1.0e-6],
        # 6. alphaB < poro or alphaB > 1
        [*base_2d, *valid_up, "-alphaB", 0.20],
        [*base_2d, *valid_up, "-alphaB", 1.05],
        # 7. Ks <= 0
        [*base_2d, *valid_up, "-Ks", 0.0],
        [*base_2d, *valid_up, "-Ks", -1.0e9],
        # 8. -pOrder linear on linear Q4
        [*base_2d, *valid_up, "-pOrder", "linear"],
    ]

    for bad_args in bad_2d_cases:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        for idx, (x, y) in enumerate(
            [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1
        ):
            ops.node(idx, x, y)
        with pytest.raises(ops.OpenSeesError):
            ops.element("LadrunoCPML", *bad_args)

    # 9. 3D parser refusals: 2 values for -perm in 3D, and -pOrder linear on linear H8
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    h8_coords = ref.make_reference_coords("H8")
    for idx in range(8):
        ops.node(idx + 1, *[float(x) for x in h8_coords[idx]])

    with pytest.raises(ops.OpenSeesError):
        ops.element(
            "LadrunoCPML",
            1,
            *range(1, 9),
            "-E",
            5.0e7,
            "-nu",
            0.25,
            "-rho",
            2000.0,
            "-box",
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            1.0,
            "-L",
            1.0,
            "-up",
            "-Kf",
            2.2e9,
            "-poro",
            0.30,
            "-rhoF",
            1000.0,
            "-perm",
            1.0e-6,
            2.0e-6,  # Only 2 values in 3D!
        )

    with pytest.raises(ops.OpenSeesError):
        ops.element(
            "LadrunoCPML",
            1,
            *range(1, 9),
            "-E",
            5.0e7,
            "-nu",
            0.25,
            "-rho",
            2000.0,
            "-box",
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            1.0,
            "-L",
            1.0,
            "-up",
            "-Kf",
            2.2e9,
            "-poro",
            0.30,
            "-rhoF",
            1000.0,
            "-perm",
            1.0e-6,
            "-pOrder",
            "linear",  # Linear H8 cannot use Taylor-Hood!
        )

    # 10. Wrong nodal ndf on Taylor-Hood midside or vertex nodes -> fail-closed in setDomain
    bt6_coords = ref.make_reference_coords("BT6")

    # (10a) Taylor-Hood BT6 with midside nodes having ndf=3 instead of expected ndf=2
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for idx in range(6):
        ops.node(idx + 1, float(bt6_coords[idx, 0]), float(bt6_coords[idx, 1]))
    capfd.readouterr()
    ops.element(
        "LadrunoCPML",
        1,
        1,
        2,
        3,
        4,
        5,
        6,
        "-shape",
        "BT6",
        "-E",
        5.0e7,
        "-nu",
        0.25,
        "-rho",
        2000.0,
        "-box",
        0.0,
        1.0,
        0.0,
        1.0,
        "-L",
        1.0,
        "-up",
        "-Kf",
        2.2e9,
        "-poro",
        0.30,
        "-rhoF",
        1000.0,
        "-perm",
        1.0e-6,
        "-pOrder",
        "linear",
    )
    err_mid = capfd.readouterr().err
    assert "DOFs, expected 2" in err_mid
    assert float(np.max(np.abs(ops.eleResponse(1, "stiff")))) == 0.0

    # (10b) Taylor-Hood BT6 with vertex nodes having ndf=2 instead of expected ndf=3
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for idx in range(6):
        ops.node(idx + 1, float(bt6_coords[idx, 0]), float(bt6_coords[idx, 1]))
    capfd.readouterr()
    ops.element(
        "LadrunoCPML",
        1,
        1,
        2,
        3,
        4,
        5,
        6,
        "-shape",
        "BT6",
        "-E",
        5.0e7,
        "-nu",
        0.25,
        "-rho",
        2000.0,
        "-box",
        0.0,
        1.0,
        0.0,
        1.0,
        "-L",
        1.0,
        "-up",
        "-Kf",
        2.2e9,
        "-poro",
        0.30,
        "-rhoF",
        1000.0,
        "-perm",
        1.0e-6,
        "-pOrder",
        "linear",
    )
    err_vert = capfd.readouterr().err
    assert "DOFs, expected 3" in err_vert
    assert float(np.max(np.abs(ops.eleResponse(1, "stiff")))) == 0.0


def test_up_anchor3_gravity_seepage_staging_rayleigh_stability_and_checkpoint(
    tmp_path: Path, capfd: Any
) -> None:
    """Anchor UP-3 (Hydrostatic Gravity + Steady Seepage `-stage 0 -> 1` Latching,
    Solid-Only Rayleigh Isolation, 2500-Step Lyapunov Stability, `FileDatastore`
    Checkpoint Round-Trip & Parser Fail-Closed Battery).
    """
    _check_up_gravity_seepage_staging()
    _check_up_rayleigh_solid_only_isolation()
    _check_up_2500_step_lyapunov_stability()
    _check_up_filedatastore_checkpoint_roundtrip(tmp_path)
    _check_up_parser_fail_closed_battery(capfd)
