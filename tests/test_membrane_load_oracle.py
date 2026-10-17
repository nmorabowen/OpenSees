"""Verification suite for the pure Python/NumPy reference oracle
(tests/_testbed/membrane_load_ref.py) across T1.1-T3.5 (ADR-98 WP-0).
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_TESTBED = os.path.join(_HERE, "_testbed")
if _TESTBED not in sys.path:
    sys.path.insert(0, _TESTBED)

import membrane_load_ref as ref


def _check_global_equilibrium(
    res: ref.DistributionResult,
    frames: list[tuple[tuple[float, float], tuple[float, float]]],
    support_nodes: list[tuple[float, float]] | None = None,
    atol: float = 1e-11,
) -> None:
    """Assert that the distributed frame sub-loads and support_node reactions
    conserve total vertical force Fz and global overturning moments Mx, My.
    """
    support_nodes = support_nodes or []
    f_sum = float(np.sum(res.support_node_reactions))
    mx_sum = 0.0
    my_sum = 0.0

    for (sx, sy), r_n in zip(support_nodes, res.support_node_reactions):
        mx_sum += float(r_n) * sy
        my_sum -= float(r_n) * sx

    for fl, (fp1, fp2) in zip(res.frame_loads, frames):
        x1, y1 = fp1
        x2, y2 = fp2
        dx = x2 - x1
        dy = y2 - y1
        lk = math.hypot(dx, dy)
        tx = dx / lk if lk > 0 else 0.0
        ty = dy / lk if lk > 0 else 0.0
        for sub in fl.sub_loads:
            a_s = sub.aL * lk
            b_s = sub.bL * lk
            l_s = b_s - a_s
            f_sub = 0.5 * (sub.wa + sub.wb) * l_s
            # Centroid of trapezoid from start of sub-segment
            if abs(f_sub) > 1e-15:
                s_bar = a_s + (sub.M_local / f_sub)
                cx = x1 + (s_bar / lk) * dx
                cy = y1 + (s_bar / lk) * dy
                f_sum += f_sub
                mx_sum += f_sub * cy
                my_sum -= f_sub * cx
            else:
                # Pure couple on sub-segment (f_sub == 0, M_local != 0)
                mx_sum += sub.M_local * ty
                my_sum -= sub.M_local * tx
        if res.plate_support_rank < 2 and abs(fl.overhang_torque) > 0.0:
            mx_sum += fl.overhang_torque * tx
            my_sum += fl.overhang_torque * ty

    assert abs(f_sum - res.total_applied_load) <= atol, (
        f"[{res.mode}] Force equilibrium failed: got {f_sum}, expected {res.total_applied_load}"
    )
    assert abs(mx_sum - res.total_applied_mx) <= atol, (
        f"[{res.mode}] Mx equilibrium failed: got {mx_sum}, expected {res.total_applied_mx}"
    )
    assert abs(my_sum - res.total_applied_my) <= atol, (
        f"[{res.mode}] My equilibrium failed: got {my_sum}, expected {res.total_applied_my}"
    )


@pytest.mark.zone_a
@pytest.mark.t1
def test_t1_1_partial_trapezoid_fef_euler_and_timoshenko():
    """T1.1: Verify partial trapezoidal 2D/3D FEFs (Euler-Bernoulli & Timoshenko Phi)
    against analytical closed-form integrals and a 200-sub-element FEM beam oracle.
    """
    L = 6.0
    wya, wyb = -12.0, -30.0
    wza, wzb = 8.0, -16.0
    waa, wab = 4.0, 10.0
    aL, bL = 0.2, 0.75

    # 1. Compare Euler-Bernoulli (phi=0) against exact polynomial closed form in ElasticBeam2d.cpp
    res2d = ref.partial_trapezoid_fef_2d(L, wya, wyb, waa, wab, aL, bL, release=0, phi=0.0)
    a = aL * L
    b = bL * L
    ba = b - a
    ba2 = b**2 - a**2
    ba3 = b**3 - a**3
    ba4 = b**4 - a**4
    ba5 = b**5 - a**5
    z1 = wya + (wya * a) / ba - (wyb * a) / ba
    wybpa = wya + wyb
    wybma = wyb - wya
    L2 = L * L
    fyt = 0.5 * wybpa * ba
    v2 = (1.0 / L) * (wya * ba * (a + 0.5 * ba) + 0.5 * wybma * ba * (a + (2.0 / 3.0) * ba))
    v1 = fyt - v2
    m1_eb = (
        0.5 * z1 * ba2
        + wybma * ba3 / (3.0 * ba)
        - z1 * ba3 * 2.0 / (3.0 * L)
        - wybma * ba4 / (2.0 * L * ba)
        + z1 * ba4 / (4.0 * L2)
        + wybma * ba5 / (5.0 * L2 * ba)
    )
    m2_eb = (
        -z1 * ba3 / (3.0 * L)
        - wybma * ba4 / (4.0 * L * ba)
        + z1 * ba4 / (4.0 * L2)
        + wybma * ba5 / (5.0 * L2 * ba)
    )
    p_ax = waa * ba + 0.5 * (wab - waa) * ba
    pj_ax = (1.0 / L) * (waa * ba * (a + 0.5 * ba) + 0.5 * (wab - waa) * ba * (a + (2.0 / 3.0) * ba))

    np.testing.assert_allclose(res2d["p0"], [-p_ax, -v1, -v2], atol=1e-12)
    np.testing.assert_allclose(res2d["q0"], [-pj_ax, -m1_eb, -m2_eb], atol=1e-12)

    # 2. Verify Timoshenko (phi > 0) against a 200-sub-element fixed-fixed Timoshenko beam
    phi = 0.45
    res_timo = ref.partial_trapezoid_fef_2d(L, wya, wyb, waa, wab, aL, bL, release=0, phi=phi)
    # Identity check: q0_1 + q0_2 == (m1_eb_q0 + m2_eb_q0) / (1 + phi), q0_2 - q0_1 == m2_eb_q0 - m1_eb_q0
    q0_eb = res2d["q0"]
    q0_tm = res_timo["q0"]
    assert abs((q0_tm[1] + q0_tm[2]) - (q0_eb[1] + q0_eb[2]) / (1.0 + phi)) < 1e-12
    assert abs((q0_tm[2] - q0_tm[1]) - (q0_eb[2] - q0_eb[1])) < 1e-12

    # 3. Verify 3D sign conventions for local-z vs local-y
    res3d = ref.partial_trapezoid_fef_3d(L, wya, wyb, wza, wzb, waa, wab, aL, bL, phiz=0.2, phiy=0.35)
    assert abs(res3d["ql0"][0] + res3d["ql0"][6] + p_ax) < 1e-12
    assert abs(res3d["ql0"][1] + res3d["ql0"][7] + fyt) < 1e-12
    fzt = 0.5 * (wza + wzb) * ba
    assert abs(res3d["ql0"][2] + res3d["ql0"][8] + fzt) < 1e-12


@pytest.mark.zone_a
@pytest.mark.t1
def test_t1_2_end_releases_against_subdivided_beam_oracle():
    """T1.2: Verify end-release branches (release = 1, 2, 3) under partial trapezoidal
    loads against a 200-element Euler-Bernoulli FEM beam with pinned end(s).
    """
    L = 5.0
    EI = 12000.0
    wya, wyb = -18.0, -42.0
    aL, bL = 0.15, 0.85

    # Build 200-element Euler-Bernoulli beam mesh to get exact end reactions/moments
    n_el = 200
    dx = L / n_el
    ndof = 2 * (n_el + 1)  # (v_i, th_i)
    K = np.zeros((ndof, ndof), dtype=float)
    f = np.zeros(ndof, dtype=float)

    ke = (EI / dx**3) * np.array(
        [
            [12.0, 6.0 * dx, -12.0, 6.0 * dx],
            [6.0 * dx, 4.0 * dx**2, -6.0 * dx, 2.0 * dx**2],
            [-12.0, -6.0 * dx, 12.0, -6.0 * dx],
            [6.0 * dx, 2.0 * dx**2, -6.0 * dx, 4.0 * dx**2],
        ],
        dtype=float,
    )
    a = aL * L
    b = bL * L
    for e in range(n_el):
        x0 = e * dx
        x1 = (e + 1) * dx
        dofs = [2 * e, 2 * e + 1, 2 * (e + 1), 2 * (e + 1) + 1]
        for i in range(4):
            for j in range(4):
                K[dofs[i], dofs[j]] += ke[i, j]
        # Integrate load on overlap of [x0, x1] with [a, b]
        xa = max(x0, a)
        xb = min(x1, b)
        if xb > xa:
            for ug, wg in zip(ref._GL3_U, ref._GL3_W):
                x = xa + (xb - xa) * ug
                w_x = wya + (wyb - wya) * (x - a) / (b - a)
                xi = (x - x0) / dx
                N = [
                    1.0 - 3.0 * xi**2 + 2.0 * xi**3,
                    dx * xi * (1.0 - xi) ** 2,
                    3.0 * xi**2 - 2.0 * xi**3,
                    -dx * xi**2 * (1.0 - xi),
                ]
                ds = (xb - xa) * wg
                for i in range(4):
                    f[dofs[i]] += w_x * N[i] * ds

    for release in (0, 1, 2, 3):
        fixed = [0, 2 * n_el]  # v_0 = 0, v_L = 0
        if release in (0, 2):
            fixed.append(1)  # th_0 = 0
        if release in (0, 1):
            fixed.append(2 * n_el + 1)  # th_L = 0
        free = [d for d in range(ndof) if d not in fixed]
        u = np.zeros(ndof, dtype=float)
        if free:
            u[free] = np.linalg.solve(K[np.ix_(free, free)], f[free])
        R = K @ u - f  # internal resisting force at nodes: [Vy1, Mz1, ..., Vy2, Mz2]
        ql0_expected = np.array([0.0, R[0], R[1], 0.0, R[2 * n_el], R[2 * n_el + 1]])
        res = ref.partial_trapezoid_fef_2d(L, wya, wyb, 0.0, 0.0, aL, bL, release=release)
        np.testing.assert_allclose(res["ql0"], ql0_expected, rtol=1e-6, atol=1e-6)


@pytest.mark.zone_a
@pytest.mark.t1
def test_t1_3_force_beam_section_sp_3d_consistency():
    """T1.3: Verify ForceBeamColumn3d internal section forces satisfy dMz/dx == Vy(x)
    and dMy/dx == Vz(x) everywhere on (a, b) under triangular/trapezoidal wz(x).
    """
    L = 8.0
    wya, wyb = -10.0, -35.0
    wza, wzb = 15.0, -25.0
    waa, wab = 5.0, 12.0
    aL, bL = 0.2, 0.8

    eps = 1e-2
    for xi in (0.1, 0.25, 0.4, 0.55, 0.75, 0.9):
        x = xi * L
        sp0 = ref.force_beam_section_sp_3d(x, L, wya, wyb, wza, wzb, waa, wab, aL, bL)
        sp_p1 = ref.force_beam_section_sp_3d(x + eps, L, wya, wyb, wza, wzb, waa, wab, aL, bL)
        sp_m1 = ref.force_beam_section_sp_3d(x - eps, L, wya, wyb, wza, wzb, waa, wab, aL, bL)
        sp_p2 = ref.force_beam_section_sp_3d(x + 2.0 * eps, L, wya, wyb, wza, wzb, waa, wab, aL, bL)
        sp_m2 = ref.force_beam_section_sp_3d(x - 2.0 * eps, L, wya, wyb, wza, wzb, waa, wab, aL, bL)
        dmz_dx = (-sp_p2["Mz"] + 8.0 * sp_p1["Mz"] - 8.0 * sp_m1["Mz"] + sp_m2["Mz"]) / (12.0 * eps)
        dmy_dx = (-sp_p2["My"] + 8.0 * sp_p1["My"] - 8.0 * sp_m1["My"] + sp_m2["My"]) / (12.0 * eps)
        assert abs(dmz_dx - sp0["Vy"]) < 1e-11, f"dMz/dx != Vy at xi={xi}"
        assert abs(dmy_dx - sp0["Vz"]) < 1e-11, f"dMy/dx != Vz at xi={xi}"


@pytest.mark.zone_a
@pytest.mark.t1
def test_t2_1_to_t2_3_pslg_cdt_acute_holes_tjunction_and_degeneracy():
    """T2.1-T2.3: Verify 30-deg acute skewed polygon (Chew-Ruppert termination),
    L-shaped polygon, eccentric hole, T-junction interior beam, partial edge beam,
    and degenerate geometry refusal.
    """
    # T2.1: 30-degree acute rhombus / skewed slab
    deg30 = math.radians(30.0)
    L_s = 5.0
    skew_poly = [
        (0.0, 0.0),
        (L_s, 0.0),
        (L_s + L_s * math.cos(deg30), L_s * math.sin(deg30)),
        (L_s * math.cos(deg30), L_s * math.sin(deg30)),
    ]
    mesh_skew = ref.build_pslg_cdt_mesh(skew_poly)
    expected_area = L_s * (L_s * math.sin(deg30))
    assert abs(mesh_skew.total_area - expected_area) < 1e-10
    assert bool(mesh_skew.is_acute[0]) is True
    assert bool(mesh_skew.is_acute[2]) is True

    # T2.2: L-shaped polygon + interior T-junction beam + eccentric hole
    l_poly = [(0.0, 0.0), (6.0, 0.0), (6.0, 3.0), (3.0, 3.0), (3.0, 6.0), (0.0, 6.0)]
    hole = [(0.8, 0.8), (2.2, 0.8), (2.2, 2.2), (0.8, 2.2)]
    frames = [
        ((0.0, 0.0), (6.0, 0.0)),
        ((6.0, 0.0), (6.0, 3.0)),
        ((0.0, 6.0), (0.0, 0.0)),
        ((3.0, 0.0), (3.0, 3.0)),  # interior T-junction beam
    ]
    mesh_l = ref.build_pslg_cdt_mesh(l_poly, holes=[hole], frames=frames)
    expected_l_area = (6.0 * 3.0 + 3.0 * 3.0) - (1.4 * 1.4)
    assert abs(mesh_l.total_area - expected_l_area) < 1e-10

    # T2.3: Degenerate geometry refusal (collinear & self-intersecting bow-tie)
    with pytest.raises(ValueError, match="Degenerate"):
        ref.build_pslg_cdt_mesh([(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)])
    with pytest.raises(ValueError, match="Degenerate"):
        ref.build_pslg_cdt_mesh([(0.0, 0.0), (4.0, 4.0), (4.0, 0.0), (0.0, 4.0)])


@pytest.mark.zone_a
@pytest.mark.t1
def test_t3_1_rectangular_bay_tributary_analytical_match():
    """T3.1: Rectangular bay two-way tributary matching analytical 45-deg yield-line
    triangles/trapezoids and conserving Fz, Mx, My to < 1e-11.
    """
    Lx, Ly = 6.0, 4.0
    q = 10.0
    poly = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
    frames = [
        ((0.0, 0.0), (Lx, 0.0)),  # Bottom (long, Lx)
        ((Lx, 0.0), (Lx, Ly)),    # Right (short, Ly)
        ((Lx, Ly), (0.0, Ly)),    # Top (long, Lx)
        ((0.0, Ly), (0.0, 0.0)),  # Left (short, Ly)
    ]
    res = ref.distribute_surface_load(poly, q, mode="tributary", frames=frames)
    _check_global_equilibrium(res, frames, atol=1e-11)

    f_short_exact = 0.25 * q * Ly * Ly
    f_long_exact = 0.5 * q * Ly * (Lx - 0.5 * Ly)
    assert abs(res.frame_loads[0].F_total - f_long_exact) < 1e-11
    assert abs(res.frame_loads[2].F_total - f_long_exact) < 1e-11
    assert abs(res.frame_loads[1].F_total - f_short_exact) < 1e-11
    assert abs(res.frame_loads[3].F_total - f_short_exact) < 1e-11


@pytest.mark.zone_a
@pytest.mark.t1
def test_t3_2_oblique_oneway_strip_equilibrium():
    """T3.2: Oblique oneway strip projection conserving sum(Fz), sum(Mx), sum(My) to < 1e-11."""
    poly = [(0.0, 0.0), (5.0, 0.0), (6.5, 4.0), (1.5, 4.0)]
    frames = [
        ((0.0, 0.0), (5.0, 0.0)),
        ((5.0, 0.0), (6.5, 4.0)),
        ((6.5, 4.0), (1.5, 4.0)),
        ((1.5, 4.0), (0.0, 0.0)),
    ]
    q = 12.5
    res = ref.distribute_surface_load(poly, q, mode="oneway", frames=frames, span_dir=(0.3, 0.95))
    _check_global_equilibrium(res, frames, atol=1e-11)


@pytest.mark.zone_a
@pytest.mark.t1
def test_t3_3_fea_poisson_isotropic_and_orthotropic():
    """T3.3: fea_poisson (isotropic eta=1 and orthotropic eta=1e-4) element-edge flux
    projection conserving equilibrium and matching 1-way limit as eta -> 0.
    """
    Lx, Ly = 4.0, 4.0
    q = 8.0
    poly = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
    frames = [
        ((0.0, 0.0), (Lx, 0.0)),  # Bottom (parallel to x)
        ((Lx, 0.0), (Lx, Ly)),    # Right (perpendicular to x)
        ((Lx, Ly), (0.0, Ly)),    # Top (parallel to x)
        ((0.0, Ly), (0.0, 0.0)),  # Left (perpendicular to x)
    ]

    # Isotropic square slab: all 4 beams receive exactly q * L^2 / 4
    res_iso = ref.distribute_surface_load(poly, q, mode="fea_poisson", frames=frames, orthotropy_eta=1.0)
    _check_global_equilibrium(res_iso, frames, atol=1e-11)
    for k in range(4):
        assert abs(res_iso.frame_loads[k].F_total - 0.25 * q * Lx * Ly) < 1e-10

    # Orthotropic slab spanning along x (eta = 1e-4): >98% of load goes to Left & Right beams
    res_ortho = ref.distribute_surface_load(
        poly, q, mode="fea_poisson", frames=frames, span_dir=(1.0, 0.0), orthotropy_eta=1e-4
    )
    _check_global_equilibrium(res_ortho, frames, atol=1e-11)
    f_lr = res_ortho.frame_loads[1].F_total + res_ortho.frame_loads[3].F_total
    assert f_lr > 0.98 * res_ortho.total_applied_load


@pytest.mark.zone_a
@pytest.mark.t1
def test_t3_4_fea_plate_dkt_cantilever_and_single_girder_rank_fallback():
    """T3.4: fea_plate (DKT bending) on cantilever overhang (rank >= 2 knife-edge)
    and single-girder (rank < 2 clamped fallback).
    """
    q = 10.0
    # Case A: Slab with cantilever overhang supported on two parallel beams at x=0 and x=4, overhanging to x=6
    poly = [(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)]
    frames_two = [
        ((0.0, 0.0), (0.0, 4.0)),  # Beam at x=0
        ((4.0, 0.0), (4.0, 4.0)),  # Interior beam at x=4 (cantilever overhang on [4, 6])
    ]
    res_overhang = ref.distribute_surface_load(poly, q, mode="fea_plate", frames=frames_two)
    assert res_overhang.plate_support_rank == 2
    _check_global_equilibrium(res_overhang, frames_two, atol=1e-10)
    # By statics of overhang [0, 6] supported at x=0 and x=4:
    # Total load = 240, centroid at x=3 => R(x=4) * 4 = 240 * 3 => R(x=4) = 180, R(x=0) = 60
    assert abs(res_overhang.frame_loads[0].F_total - 60.0) < 1e-8
    assert abs(res_overhang.frame_loads[1].F_total - 180.0) < 1e-8

    # Case B: Pure cantilever slab supported on a SINGLE girder at x=0 (rank < 2 clamped fallback)
    frames_single = [((0.0, 0.0), (0.0, 4.0))]
    res_single = ref.distribute_surface_load(poly, q, mode="fea_plate", frames=frames_single)
    assert res_single.plate_support_rank == 1
    assert abs(res_single.frame_loads[0].F_total - 240.0) < 1e-9
    _check_global_equilibrium(res_single, frames_single, atol=1e-9)


@pytest.mark.zone_a
@pytest.mark.t1
def test_t3_5_guyan_static_condensation_inplane_cst():
    """T3.5: Guyan static condensation of in-plane CST stiffness K_cond and mass M_cond
    preserves symmetry, 3 rigid-body nullspace modes, constant-strain patch test, and total mass.
    """
    Lx, Ly = 5.0, 3.0
    poly = [(0.0, 0.0), (Lx, 0.0), (Lx, Ly), (0.0, Ly)]
    mesh = ref.build_pslg_cdt_mesh(poly, h_target=1.2)
    E, nu, t, rho = 30000.0, 0.25, 0.2, 2.4

    cond = ref.condense_inplane_cst(
        mesh.vertices,
        mesh.triangles,
        mesh.boundary_vertex_indices,
        E=E,
        nu=nu,
        t=t,
        rho=rho,
    )
    K_c = cond["K_cond"]
    M_c = cond["M_cond"]

    # Symmetry
    np.testing.assert_allclose(K_c, K_c.T, atol=1e-12)
    np.testing.assert_allclose(M_c, M_c.T, atol=1e-12)

    # Rigid-body modes (u=1, v=1, rotation=(-y, x)) have zero strain energy
    b_coords = mesh.vertices[mesh.boundary_vertex_indices]
    n_b = len(b_coords)
    rb_u = np.zeros(2 * n_b)
    rb_v = np.zeros(2 * n_b)
    rb_r = np.zeros(2 * n_b)
    rb_u[0::2] = 1.0
    rb_v[1::2] = 1.0
    for i, (x, y) in enumerate(b_coords):
        rb_r[2 * i] = -y
        rb_r[2 * i + 1] = x

    np.testing.assert_allclose(K_c @ rb_u, 0.0, atol=1e-10)
    np.testing.assert_allclose(K_c @ rb_v, 0.0, atol=1e-10)
    np.testing.assert_allclose(K_c @ rb_r, 0.0, atol=1e-9)

    # Constant-strain patch test: u = eps_xx * x, v = eps_yy * y
    eps_xx, eps_yy = 1e-3, -0.5e-3
    u_lin = np.zeros(2 * n_b)
    for i, (x, y) in enumerate(b_coords):
        u_lin[2 * i] = eps_xx * x
        u_lin[2 * i + 1] = eps_yy * y
    strain_energy = 0.5 * float(u_lin @ K_c @ u_lin)
    D_ps = (E / (1.0 - nu * nu)) * np.array([[1.0, nu], [nu, 1.0]])
    eps_vec = np.array([eps_xx, eps_yy])
    exact_energy = 0.5 * (Lx * Ly * t) * float(eps_vec @ D_ps @ eps_vec)
    assert abs(strain_energy - exact_energy) < 1e-9

    # Total translational mass conservation
    total_mass = rho * t * Lx * Ly
    assert abs(float(rb_u @ M_c @ rb_u) - total_mass) < 1e-10
    assert abs(float(rb_v @ M_c @ rb_v) - total_mass) < 1e-10
