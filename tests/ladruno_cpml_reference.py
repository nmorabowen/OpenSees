"""Independent NumPy reference oracle for the LadrunoCPML kernel (Phase P0, WP 0.B).

This standalone reference implementation derives all 9 topologies (T3, Q4, T6, BT6,
TET4, H8, TET10, BTET10, H20), both 'box' and 'plane' CFS-PML coordinate stretching
profiles, exact Piecewise-Linear Recursive Convolution (RCFilter1D) cascades,
consistent and HRZ-lumped mass/damping/inertial stiffness, internal resisting force
R_int, tangent matrices (K_eff, C_eff, M_eff, K_T_dyn), and external power P_ext
from first principles according to the frozen P0 mathematical contract.

Blind-Oracle Protocol:
- Reads nothing under SRC/element/ladrunoCPML/ or sibling C++ test files.
- Includes built-in --self-check and --emit <path> modes.

Usage:
    python tests/ladruno_cpml_reference.py --self-check
    python tests/ladruno_cpml_reference.py --self-check --emit tests/ladruno_cpml_cases.txt
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass

import numpy as np

# ============================================================================ #
# 1. Topologies & Quadrature Rules (9 Topologies)                              #
# ============================================================================ #

_INV_SQRT3 = 0.57735026918962576451
_TET10_A = 0.5854101966249685
_TET10_B = 0.1381966011250105


def shape_t3(xi: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    r, s = xi
    l0 = 1.0 - r - s
    N = np.array([l0, r, s], dtype=np.float64)
    dN_dxi = np.array([
        [-1.0, -1.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ], dtype=np.float64)
    return N, dN_dxi


def gp_t3() -> list[tuple[tuple[float, float], float]]:
    return [((1.0 / 3.0, 1.0 / 3.0), 0.5)]


def shape_q4(xi: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    x, y = xi
    nodes = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
    N = np.zeros(4, dtype=np.float64)
    dN_dxi = np.zeros((4, 2), dtype=np.float64)
    for a, (xa, ya) in enumerate(nodes):
        N[a] = 0.25 * (1.0 + xa * x) * (1.0 + ya * y)
        dN_dxi[a, 0] = 0.25 * xa * (1.0 + ya * y)
        dN_dxi[a, 1] = 0.25 * (1.0 + xa * x) * ya
    return N, dN_dxi


def gp_q4() -> list[tuple[tuple[float, float], float]]:
    g = _INV_SQRT3
    pts = []
    # j (eta) outer from -g to +g, i (xi) inner from -g to +g
    for eta in (-g, +g):
        for xi in (-g, +g):
            pts.append(((xi, eta), 1.0))
    return pts


def shape_t6(xi: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    s, t = xi
    l3 = 1.0 - s - t
    N = np.array([
        s * (2.0 * s - 1.0),
        t * (2.0 * t - 1.0),
        l3 * (2.0 * l3 - 1.0),
        4.0 * s * t,
        4.0 * t * l3,
        4.0 * s * l3,
    ], dtype=np.float64)
    dN_dxi = np.array([
        [4.0 * s - 1.0, 0.0],
        [0.0, 4.0 * t - 1.0],
        [1.0 - 4.0 * l3, 1.0 - 4.0 * l3],
        [4.0 * t, 4.0 * s],
        [-4.0 * t, 4.0 * (l3 - t)],
        [4.0 * (l3 - s), -4.0 * s],
    ], dtype=np.float64)
    return N, dN_dxi


def gp_t6() -> list[tuple[tuple[float, float], float]]:
    w = 1.0 / 6.0
    return [
        ((1.0 / 6.0, 1.0 / 6.0), w),
        ((2.0 / 3.0, 1.0 / 6.0), w),
        ((1.0 / 6.0, 2.0 / 3.0), w),
    ]


def shape_bt6(xi: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    l1, l2 = xi
    l0 = 1.0 - l1 - l2
    N = np.array([
        l0 * l0,
        l1 * l1,
        l2 * l2,
        2.0 * l0 * l1,
        2.0 * l1 * l2,
        2.0 * l2 * l0,
    ], dtype=np.float64)
    dN_dxi = np.array([
        [-2.0 * l0, -2.0 * l0],
        [2.0 * l1, 0.0],
        [0.0, 2.0 * l2],
        [2.0 * (l0 - l1), -2.0 * l1],
        [2.0 * l2, 2.0 * l1],
        [-2.0 * l2, 2.0 * (l0 - l2)],
    ], dtype=np.float64)
    return N, dN_dxi


def gp_bt6() -> list[tuple[tuple[float, float], float]]:
    w = 1.0 / 6.0
    return [
        ((2.0 / 3.0, 1.0 / 6.0), w),
        ((1.0 / 6.0, 2.0 / 3.0), w),
        ((1.0 / 6.0, 1.0 / 6.0), w),
    ]


def shape_tet4(xi: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    r, s, t = xi
    l0 = 1.0 - r - s - t
    N = np.array([l0, r, s, t], dtype=np.float64)
    dN_dxi = np.array([
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=np.float64)
    return N, dN_dxi


def gp_tet4() -> list[tuple[tuple[float, float, float], float]]:
    return [((0.25, 0.25, 0.25), 1.0 / 6.0)]


def shape_h8(xi: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    x, y, z = xi
    nodes = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ]
    N = np.zeros(8, dtype=np.float64)
    dN_dxi = np.zeros((8, 3), dtype=np.float64)
    for a, (xa, ya, za) in enumerate(nodes):
        fx = 1.0 + xa * x
        fy = 1.0 + ya * y
        fz = 1.0 + za * z
        N[a] = 0.125 * fx * fy * fz
        dN_dxi[a, 0] = 0.125 * xa * fy * fz
        dN_dxi[a, 1] = 0.125 * fx * ya * fz
        dN_dxi[a, 2] = 0.125 * fx * fy * za
    return N, dN_dxi


def gp_h8() -> list[tuple[tuple[float, float, float], float]]:
    g = _INV_SQRT3
    pts = []
    # k (zeta) outer, j (eta) middle, i (xi) inner
    for zeta in (-g, +g):
        for eta in (-g, +g):
            for xi in (-g, +g):
                pts.append(((xi, eta, zeta), 1.0))
    return pts


def shape_tet10(xi: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    z1, z2, z3 = xi
    z4 = 1.0 - z1 - z2 - z3
    N = np.array([
        z1 * (2.0 * z1 - 1.0),
        z2 * (2.0 * z2 - 1.0),
        z3 * (2.0 * z3 - 1.0),
        z4 * (2.0 * z4 - 1.0),
        4.0 * z1 * z2,
        4.0 * z2 * z3,
        4.0 * z3 * z1,
        4.0 * z1 * z4,
        4.0 * z3 * z4,
        4.0 * z2 * z4,
    ], dtype=np.float64)
    dz4 = -4.0 * z4 + 1.0
    dN_dxi = np.array([
        [4.0 * z1 - 1.0, 0.0, 0.0],
        [0.0, 4.0 * z2 - 1.0, 0.0],
        [0.0, 0.0, 4.0 * z3 - 1.0],
        [dz4, dz4, dz4],
        [4.0 * z2, 4.0 * z1, 0.0],
        [0.0, 4.0 * z3, 4.0 * z2],
        [4.0 * z3, 0.0, 4.0 * z1],
        [4.0 * (z4 - z1), -4.0 * z1, -4.0 * z1],
        [-4.0 * z3, -4.0 * z3, 4.0 * (z4 - z3)],
        [-4.0 * z2, 4.0 * (z4 - z2), -4.0 * z2],
    ], dtype=np.float64)
    return N, dN_dxi


def gp_tet10() -> list[tuple[tuple[float, float, float], float]]:
    a = _TET10_A
    b = _TET10_B
    w = 1.0 / 24.0
    return [
        ((a, b, b), w),
        ((b, a, b), w),
        ((b, b, a), w),
        ((b, b, b), w),
    ]


def shape_btet10(xi: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    z1, z2, z3 = xi
    z4 = 1.0 - z1 - z2 - z3
    N = np.array([
        z1 * z1,
        z2 * z2,
        z3 * z3,
        z4 * z4,
        2.0 * z1 * z2,
        2.0 * z2 * z3,
        2.0 * z3 * z1,
        2.0 * z1 * z4,
        2.0 * z3 * z4,
        2.0 * z2 * z4,
    ], dtype=np.float64)
    dN_dxi = np.array([
        [2.0 * z1, 0.0, 0.0],
        [0.0, 2.0 * z2, 0.0],
        [0.0, 0.0, 2.0 * z3],
        [-2.0 * z4, -2.0 * z4, -2.0 * z4],
        [2.0 * z2, 2.0 * z1, 0.0],
        [0.0, 2.0 * z3, 2.0 * z2],
        [2.0 * z3, 0.0, 2.0 * z1],
        [2.0 * (z4 - z1), -2.0 * z1, -2.0 * z1],
        [-2.0 * z3, -2.0 * z3, 2.0 * (z4 - z3)],
        [-2.0 * z2, 2.0 * (z4 - z2), -2.0 * z2],
    ], dtype=np.float64)
    return N, dN_dxi


def gp_btet10() -> list[tuple[tuple[float, float, float], float]]:
    return gp_tet10()


def shape_h20(xi: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    x, y, z = xi
    corners = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ]
    mids = [
        (0.0, -1.0, -1.0),
        (1.0, 0.0, -1.0),
        (0.0, 1.0, -1.0),
        (-1.0, 0.0, -1.0),
        (0.0, -1.0, 1.0),
        (1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
        (-1.0, 0.0, 1.0),
        (-1.0, -1.0, 0.0),
        (1.0, -1.0, 0.0),
        (1.0, 1.0, 0.0),
        (-1.0, 1.0, 0.0),
    ]
    N = np.zeros(20, dtype=np.float64)
    dN_dxi = np.zeros((20, 3), dtype=np.float64)
    for a, (xa, ya, za) in enumerate(corners):
        fx = 1.0 + xa * x
        fy = 1.0 + ya * y
        fz = 1.0 + za * z
        s = xa * x + ya * y + za * z - 2.0
        N[a] = 0.125 * fx * fy * fz * s
        dN_dxi[a, 0] = 0.125 * xa * fy * fz * (s + fx)
        dN_dxi[a, 1] = 0.125 * fx * ya * fz * (s + fy)
        dN_dxi[a, 2] = 0.125 * fx * fy * za * (s + fz)

    for idx, (xa, ya, za) in enumerate(mids):
        a = 8 + idx
        if xa == 0.0:
            fx = 1.0 - x * x
            fy = 1.0 + ya * y
            fz = 1.0 + za * z
            N[a] = 0.25 * fx * fy * fz
            dN_dxi[a, 0] = 0.25 * (-2.0 * x) * fy * fz
            dN_dxi[a, 1] = 0.25 * fx * ya * fz
            dN_dxi[a, 2] = 0.25 * fx * fy * za
        elif ya == 0.0:
            fx = 1.0 + xa * x
            fy = 1.0 - y * y
            fz = 1.0 + za * z
            N[a] = 0.25 * fx * fy * fz
            dN_dxi[a, 0] = 0.25 * xa * fy * fz
            dN_dxi[a, 1] = 0.25 * fx * (-2.0 * y) * fz
            dN_dxi[a, 2] = 0.25 * fx * fy * za
        else:
            fx = 1.0 + xa * x
            fy = 1.0 + ya * y
            fz = 1.0 - z * z
            N[a] = 0.25 * fx * fy * fz
            dN_dxi[a, 0] = 0.25 * xa * fy * fz
            dN_dxi[a, 1] = 0.25 * fx * ya * fz
            dN_dxi[a, 2] = 0.25 * fx * fy * (-2.0 * z)
    return N, dN_dxi


def gp_h20() -> list[tuple[tuple[float, float, float], float]]:
    return gp_h8()


TOPOLOGIES = {
    "T3": (2, 3, shape_t3, gp_t3),
    "Q4": (2, 4, shape_q4, gp_q4),
    "T6": (2, 6, shape_t6, gp_t6),
    "BT6": (2, 6, shape_bt6, gp_bt6),
    "TET4": (3, 4, shape_tet4, gp_tet4),
    "H8": (3, 8, shape_h8, gp_h8),
    "TET10": (3, 10, shape_tet10, gp_tet10),
    "BTET10": (3, 10, shape_btet10, gp_btet10),
    "H20": (3, 20, shape_h20, gp_h20),
}


def eval_gp_kinematics(
    topo: str, coords: np.ndarray
) -> list[tuple[np.ndarray, np.ndarray, float, np.ndarray]]:
    """Return list of (N, dN_dx, w_vol, x_g) for each Gauss point of `topo`."""
    ndm, n_nodes, shape_fn, gp_fn = TOPOLOGIES[topo]
    assert coords.shape == (n_nodes, ndm)
    out = []
    for xi, w_g in gp_fn():
        N, dN_dxi = shape_fn(xi)
        # J_{mp} = d x_p / d xi_m = sum_a (dN_a / d xi_m) x_{a, p}
        J = dN_dxi.T @ coords
        detJ = float(np.linalg.det(J))
        if detJ <= 0.0:
            raise ValueError(f"Non-positive Jacobian determinant {detJ} in {topo}")
        # dN_dxi.T = J @ dN_dx.T  =>  dN_dx = solve(J, dN_dxi.T).T
        dN_dx = np.linalg.solve(J, dN_dxi.T).T
        w_vol = w_g * detJ
        x_g = N @ coords
        out.append((N, dN_dx, w_vol, x_g))
    return out


# ============================================================================ #
# 2. Isotropic Elastic Constitutive Tensor & Wave Speeds                       #
# ============================================================================ #

def isotropic_elasticity(E: float, nu: float, rho: float, ndm: int) -> tuple[float, float, float, np.ndarray]:
    mu = E / (2.0 * (1.0 + nu))
    lam = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    Vp = math.sqrt((lam + 2.0 * mu) / rho)
    C = np.zeros((ndm, ndm, ndm, ndm), dtype=np.float64)
    for i in range(ndm):
        for j in range(ndm):
            for k in range(ndm):
                for l in range(ndm):
                    val = 0.0
                    if i == j and k == l:
                        val += lam
                    if i == k and j == l:
                        val += mu
                    if i == l and j == k:
                        val += mu
                    C[i, j, k, l] = val
    return lam, mu, Vp, C


# ============================================================================ #
# 3. Geometry Modes ('box' and 'plane') & CFS-PML Stretching Profiles          #
# ============================================================================ #

@dataclass
class PMLParams:
    m_deg: float
    R0: float
    kappa_max: float
    alpha_max: float


@dataclass
class GeomSpec:
    mode: str  # "box" or "plane"
    params: tuple[float, ...]


def compute_rotation_matrix(ndm: int, geom: GeomSpec) -> np.ndarray:
    if geom.mode == "box":
        return np.eye(ndm, dtype=np.float64)
    elif geom.mode == "plane":
        if ndm == 2:
            assert len(geom.params) == 5
            nx, ny = geom.params[2], geom.params[3]
            norm = math.hypot(nx, ny)
            nx /= norm
            ny /= norm
            return np.array([
                [nx, ny],
                [-ny, nx],
            ], dtype=np.float64)
        elif ndm == 3:
            assert len(geom.params) == 7
            n = np.array(geom.params[3:6], dtype=np.float64)
            n = n / np.linalg.norm(n)
            a = np.array([1.0, 0.0, 0.0], dtype=np.float64) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0], dtype=np.float64)
            e2 = np.cross(n, a)
            e2 = e2 / np.linalg.norm(e2)
            e3 = np.cross(n, e2)
            return np.vstack([n, e2, e3])
    raise ValueError(f"Unknown geom mode {geom.mode} for ndm={ndm}")


def evaluate_pml_profile_at_gp(
    x_g: np.ndarray,
    Vp: float,
    pml: PMLParams,
    geom: GeomSpec,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (d, kappa, alpha) arrays of length 3 (with m=2 inactive in 2D)."""
    ndm = len(x_g)
    xi_depth = np.zeros(3, dtype=np.float64)
    L_arr = np.ones(3, dtype=np.float64)

    if geom.mode == "box":
        if ndm == 2:
            x_min, x_max, y_min, y_max, Lx, Ly = geom.params
            xi_depth[0] = max(0.0, x_min - x_g[0], x_g[0] - x_max)
            xi_depth[1] = max(0.0, y_min - x_g[1], x_g[1] - y_max)
            L_arr[0] = Lx
            L_arr[1] = Ly
        else:
            x_min, x_max, y_min, y_max, z_min, z_max, Lx, Ly, Lz = geom.params
            xi_depth[0] = max(0.0, x_min - x_g[0], x_g[0] - x_max)
            xi_depth[1] = max(0.0, y_min - x_g[1], x_g[1] - y_max)
            xi_depth[2] = max(0.0, z_min - x_g[2], x_g[2] - z_max)
            L_arr[0] = Lx
            L_arr[1] = Ly
            L_arr[2] = Lz
    elif geom.mode == "plane":
        if ndm == 2:
            x0, y0, nx, ny, L_pml = geom.params
            n = np.array([nx, ny], dtype=np.float64)
            n = n / np.linalg.norm(n)
            x_ref = np.array([x0, y0], dtype=np.float64)
            xi_depth[0] = max(0.0, float(np.dot(n, x_g - x_ref)))
            L_arr[:] = L_pml
        else:
            x0, y0, z0, nx, ny, nz, L_pml = geom.params
            n = np.array([nx, ny, nz], dtype=np.float64)
            n = n / np.linalg.norm(n)
            x_ref = np.array([x0, y0, z0], dtype=np.float64)
            xi_depth[0] = max(0.0, float(np.dot(n, x_g - x_ref)))
            L_arr[:] = L_pml
    else:
        raise ValueError(f"Unknown geom mode {geom.mode}")

    d = np.zeros(3, dtype=np.float64)
    kappa = np.ones(3, dtype=np.float64)
    alpha = np.zeros(3, dtype=np.float64)

    for m in range(ndm):
        xi_bar = min(1.0, max(0.0, xi_depth[m] / L_arr[m]))
        if xi_bar == 0.0:
            d[m] = 0.0
            kappa[m] = 1.0
            alpha[m] = 0.0
        else:
            d0_m = -((pml.m_deg + 1.0) * Vp * math.log(pml.R0)) / (2.0 * L_arr[m])
            xm = xi_bar ** pml.m_deg
            d[m] = d0_m * xm
            kappa[m] = 1.0 + (pml.kappa_max - 1.0) * xm
            alpha[m] = pml.alpha_max * (1.0 - xi_bar)

    return d, kappa, alpha


# ============================================================================ #
# 4. Exact Piecewise-Linear Recursive Convolution (RCFilter1D) & Cascades      #
# ============================================================================ #

def rc_coeffs(beta: float, dt: float) -> tuple[float, float, float]:
    """Compute exact piecewise-linear recursive convolution coefficients (b, phi0, phi1)."""
    theta = beta * dt
    if abs(theta) < 1.0e-8:
        b = 1.0 - theta + 0.5 * (theta ** 2) - (theta ** 3) / 6.0 + (theta ** 4) / 24.0
        phi0 = dt * (0.5 - theta / 3.0 + (theta ** 2) / 8.0 - (theta ** 3) / 30.0)
        phi1 = dt * (0.5 - theta / 6.0 + (theta ** 2) / 24.0 - (theta ** 3) / 120.0)
    else:
        b = math.exp(-theta)
        denom = (beta * beta) * dt
        phi0 = (1.0 - b * (1.0 + theta)) / denom
        phi1 = (theta - 1.0 + b) / denom
    return b, phi0, phi1


class RCFilter1D:
    """Single first-order rational filter H(w) = c0 + c1 / (beta + i*w)."""

    def __init__(self, c0: float, c1: float, beta: float) -> None:
        self.c0 = float(c0)
        self.c1 = float(c1)
        self.beta = float(beta)
        self.x_commit = 0.0
        self.zeta_commit = 0.0
        self.x_trial = 0.0
        self.zeta_trial = 0.0

    def evaluate_trial(self, x_new: float, dt: float) -> tuple[float, float]:
        """Return (y_trial, D_H) at t_{n+1} without modifying committed state."""
        if self.c1 == 0.0:
            self.x_trial = x_new
            self.zeta_trial = 0.0
            return self.c0 * x_new, self.c0
        b, phi0, phi1 = rc_coeffs(self.beta, dt)
        self.x_trial = x_new
        self.zeta_trial = b * self.zeta_commit + phi0 * self.x_commit + phi1 * x_new
        y_new = self.c0 * x_new + self.c1 * self.zeta_trial
        D_H = self.c0 + self.c1 * phi1
        return y_new, D_H

    def commit(self) -> None:
        self.x_commit = self.x_trial
        self.zeta_commit = self.zeta_trial


class FilterCascade:
    """Ordered K-stage cascade of RCFilter1D stages: stage K o ... o stage 1."""

    def __init__(self, stages: list[RCFilter1D]) -> None:
        self.stages = stages

    def evaluate_trial(self, x_in: float, dt: float) -> tuple[float, float]:
        val = x_in
        D_total = 1.0
        for st in self.stages:
            val, D_st = st.evaluate_trial(val, dt)
            D_total *= D_st
        return val, D_total

    def commit(self) -> None:
        for st in self.stages:
            st.commit()


# ============================================================================ #
# 5. Elementary Filters & Gauss-Point Filter Banks                             #
# ============================================================================ #

def make_s_filter(d_m: float, kappa_m: float, alpha_m: float) -> RCFilter1D:
    return RCFilter1D(c0=kappa_m, c1=d_m, beta=alpha_m)


def make_s_inv_filter(d_m: float, kappa_m: float, alpha_m: float) -> RCFilter1D:
    hat_alpha = alpha_m + d_m / kappa_m
    return RCFilter1D(c0=1.0 / kappa_m, c1=-d_m / (kappa_m * kappa_m), beta=hat_alpha)


def make_g_filter(d_m: float, kappa_m: float, alpha_m: float) -> RCFilter1D:
    return RCFilter1D(c0=d_m, c1=-alpha_m * d_m, beta=alpha_m)


def make_phi_filter(d_m: float, kappa_m: float, alpha_m: float) -> RCFilter1D:
    return RCFilter1D(c0=0.0, c1=d_m, beta=alpha_m)


def make_stiffness_cascade(
    ndm: int, j: int, l: int, d: np.ndarray, kappa: np.ndarray, alpha: np.ndarray
) -> FilterCascade:
    """Create the exact contractual stage cascade for F_{jl}(w) = J_s(w) / (s_j(w) s_l(w))."""
    if ndm == 2:
        if j == 0 and l == 0:
            return FilterCascade([
                make_s_inv_filter(d[0], kappa[0], alpha[0]),
                make_s_filter(d[1], kappa[1], alpha[1]),
            ])
        elif j == 1 and l == 1:
            return FilterCascade([
                make_s_inv_filter(d[1], kappa[1], alpha[1]),
                make_s_filter(d[0], kappa[0], alpha[0]),
            ])
        else:
            return FilterCascade([])
    else:
        if j == 0 and l == 0:
            return FilterCascade([
                make_s_inv_filter(d[0], kappa[0], alpha[0]),
                make_s_filter(d[1], kappa[1], alpha[1]),
                make_s_filter(d[2], kappa[2], alpha[2]),
            ])
        elif j == 1 and l == 1:
            return FilterCascade([
                make_s_inv_filter(d[1], kappa[1], alpha[1]),
                make_s_filter(d[0], kappa[0], alpha[0]),
                make_s_filter(d[2], kappa[2], alpha[2]),
            ])
        elif j == 2 and l == 2:
            return FilterCascade([
                make_s_inv_filter(d[2], kappa[2], alpha[2]),
                make_s_filter(d[0], kappa[0], alpha[0]),
                make_s_filter(d[1], kappa[1], alpha[1]),
            ])
        elif (j == 0 and l == 1) or (j == 1 and l == 0):
            return FilterCascade([make_s_filter(d[2], kappa[2], alpha[2])])
        elif (j == 0 and l == 2) or (j == 2 and l == 0):
            return FilterCascade([make_s_filter(d[1], kappa[1], alpha[1])])
        elif (j == 1 and l == 2) or (j == 2 and l == 1):
            return FilterCascade([make_s_filter(d[0], kappa[0], alpha[0])])
        else:
            raise ValueError(f"Invalid (j, l)=({j}, {l})")


class GaussPointState:
    """Holds all filter cascades for one Gauss point."""

    def __init__(
        self,
        ndm: int,
        N: np.ndarray,
        dN_dx: np.ndarray,
        w_vol: float,
        x_g: np.ndarray,
        R: np.ndarray,
        d: np.ndarray,
        kappa: np.ndarray,
        alpha: np.ndarray,
    ) -> None:
        self.ndm = ndm
        self.N = N
        self.dN_dx = dN_dx
        # Local shape function gradient: Phi'_{a, j} = sum_p R_{j, p} dN_a / dx_p
        self.Phi_local = dN_dx @ R.T
        self.w_vol = w_vol
        self.x_g = x_g
        self.d = d
        self.kappa = kappa
        self.alpha = alpha

        # For each pair (j, l) and local displacement component k and coordinate deriv l:
        # We filter G'_{kl} through F_{jl}.
        self.F_cascades = {}
        for j in range(ndm):
            for l in range(ndm):
                self.F_cascades[(j, l)] = [
                    make_stiffness_cascade(ndm, j, l, d, kappa, alpha)
                    for _ in range(ndm)  # indexed by k in 0..ndm-1
                ]

        # Inertial velocity filters: g_0, g_1, g_2 for each component i in 0..ndm-1
        self.vel_g0 = [make_g_filter(d[0], kappa[0], alpha[0]) for _ in range(ndm)]
        self.vel_g1 = [make_g_filter(d[1], kappa[1], alpha[1]) for _ in range(ndm)]
        self.vel_g2 = [make_g_filter(d[2], kappa[2], alpha[2]) for _ in range(ndm)]

        # Inertial displacement cascades for each component i in 0..ndm-1:
        # (g_1 o g_0), (g_2 o g_0), (g_2 o g_1), (phi_2 o g_1 o g_0)
        self.disp_g1_g0 = [
            FilterCascade([
                make_g_filter(d[0], kappa[0], alpha[0]),
                make_g_filter(d[1], kappa[1], alpha[1]),
            ])
            for _ in range(ndm)
        ]
        self.disp_g2_g0 = [
            FilterCascade([
                make_g_filter(d[0], kappa[0], alpha[0]),
                make_g_filter(d[2], kappa[2], alpha[2]),
            ])
            for _ in range(ndm)
        ]
        self.disp_g2_g1 = [
            FilterCascade([
                make_g_filter(d[1], kappa[1], alpha[1]),
                make_g_filter(d[2], kappa[2], alpha[2]),
            ])
            for _ in range(ndm)
        ]
        self.disp_phi2_g1_g0 = [
            FilterCascade([
                make_g_filter(d[0], kappa[0], alpha[0]),
                make_g_filter(d[1], kappa[1], alpha[1]),
                make_phi_filter(d[2], kappa[2], alpha[2]),
            ])
            for _ in range(ndm)
        ]

    def commit(self) -> None:
        for cascades in self.F_cascades.values():
            for c in cascades:
                c.commit()
        for f in self.vel_g0 + self.vel_g1 + self.vel_g2:
            f.commit()
        for c in (
            self.disp_g1_g0
            + self.disp_g2_g0
            + self.disp_g2_g1
            + self.disp_phi2_g1_g0
        ):
            c.commit()


# ============================================================================ #
# 6 & 7. Full Element Reference Oracle (Stiffness + Inertial + HRZ Lumping)    #
# ============================================================================ #

class ReferenceCPMLElement:
    """Standalone reference implementation of a single LadrunoCPML element."""

    def __init__(
        self,
        topo: str,
        coords: np.ndarray,
        E: float,
        nu: float,
        rho: float,
        pml: PMLParams,
        geom: GeomSpec,
        lumped: int,
    ) -> None:
        self.topo = topo
        self.ndm, self.n_nodes, _, _ = TOPOLOGIES[topo]
        self.n_dof = self.n_nodes * self.ndm
        self.coords = np.array(coords, dtype=np.float64).reshape(self.n_nodes, self.ndm)
        self.E = float(E)
        self.nu = float(nu)
        self.rho = float(rho)
        self.pml = pml
        self.geom = geom
        self.lumped = int(lumped)

        self.lam, self.mu, self.Vp, self.C_tensor = isotropic_elasticity(
            self.E, self.nu, self.rho, self.ndm
        )
        self.R = compute_rotation_matrix(self.ndm, self.geom)

        gp_kin = eval_gp_kinematics(self.topo, self.coords)
        self.gps: list[GaussPointState] = []
        for N, dN_dx, w_vol, x_g in gp_kin:
            d, kappa, alpha = evaluate_pml_profile_at_gp(x_g, self.Vp, self.pml, self.geom)
            self.gps.append(
                GaussPointState(
                    ndm=self.ndm,
                    N=N,
                    dN_dx=dN_dx,
                    w_vol=w_vol,
                    x_g=x_g,
                    R=self.R,
                    d=d,
                    kappa=kappa,
                    alpha=alpha,
                )
            )

        # Precompute static effective mass matrix M_eff (depends only on kappa_0*kappa_1*kappa_2)
        M_scalar_cons = np.zeros((self.n_nodes, self.n_nodes), dtype=np.float64)
        for gp in self.gps:
            D_acc = float(gp.kappa[0] * gp.kappa[1] * gp.kappa[2])
            M_scalar_cons += (self.rho * gp.w_vol * D_acc) * np.outer(gp.N, gp.N)

        self.M_eff = np.zeros((self.n_dof, self.n_dof), dtype=np.float64)
        if self.lumped == 0:
            for a in range(self.n_nodes):
                for b in range(self.n_nodes):
                    val = M_scalar_cons[a, b]
                    for i in range(self.ndm):
                        self.M_eff[a * self.ndm + i, b * self.ndm + i] = val
        else:
            M_tot = float(np.sum(M_scalar_cons))
            D_sum = float(np.trace(M_scalar_cons))
            for a in range(self.n_nodes):
                m_a = M_tot * (M_scalar_cons[a, a] / D_sum)
                for i in range(self.ndm):
                    self.M_eff[a * self.ndm + i, a * self.ndm + i] = m_a

    def evaluate_trial(
        self,
        u_vec: np.ndarray,
        v_vec: np.ndarray,
        a_vec: np.ndarray,
        dt: float,
        c1: float,
        c2: float,
        c3: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        """Evaluate R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext at trial state (u, v, a)."""
        ndm = self.ndm
        n_nodes = self.n_nodes
        n_dof = self.n_dof

        u_nodes = u_vec.reshape(n_nodes, ndm)
        v_nodes = v_vec.reshape(n_nodes, ndm)
        # Local nodal displacements u'_a = R u_a -> shape (n_nodes, ndm)
        u_local = u_nodes @ self.R.T

        f_K_local = np.zeros((n_nodes, ndm), dtype=np.float64)
        K_stiff_local = np.zeros((n_nodes, ndm, n_nodes, ndm), dtype=np.float64)

        f_vel = np.zeros((n_nodes, ndm), dtype=np.float64)
        f_disp = np.zeros((n_nodes, ndm), dtype=np.float64)

        C_scalar = np.zeros((n_nodes, n_nodes), dtype=np.float64)
        K_inert_scalar = np.zeros((n_nodes, n_nodes), dtype=np.float64)

        for gp in self.gps:
            w = gp.w_vol
            Phi = gp.Phi_local  # (n_nodes, ndm)
            # G'_{kl} = sum_a u'_{a, k} Phi'_{a, l}
            G_local = u_local.T @ Phi  # (ndm, ndm)

            # Filter G'_{kl} through F_{jl}
            G_tilde = np.zeros((ndm, ndm, ndm, ndm), dtype=np.float64)  # [j, l, k] -> store in [j, l, k, 0]
            D_F = np.zeros((ndm, ndm), dtype=np.float64)
            for j in range(ndm):
                for l in range(ndm):
                    cascades_jl = gp.F_cascades[(j, l)]
                    for k in range(ndm):
                        y_val, D_val = cascades_jl[k].evaluate_trial(float(G_local[k, l]), dt)
                        G_tilde[j, l, k, 0] = y_val
                        D_F[j, l] = D_val

            # Local stiffness force f'_{K, a, i} and tangent K'_{eff, (a,i), (b,k)}
            for j in range(ndm):
                for l in range(ndm):
                    df_jl = D_F[j, l]
                    for i in range(ndm):
                        s_il = 0.0
                        for k in range(ndm):
                            c_ijkl = self.C_tensor[i, j, k, l]
                            if c_ijkl != 0.0:
                                s_il += c_ijkl * G_tilde[j, l, k, 0]
                                # Tangent contribution: w * Phi_{a,j} * C_{ijkl} * D_{F_{jl}} * Phi_{b,l}
                                K_stiff_local[:, i, :, k] += (
                                    w * c_ijkl * df_jl
                                ) * np.outer(Phi[:, j], Phi[:, l])
                        if s_il != 0.0:
                            f_K_local[:, i] += (w * s_il) * Phi[:, j]

            # Inertial velocity and displacement filter evaluations at Gauss point g
            u_g = gp.N @ u_nodes  # (ndm,)
            v_g = gp.N @ v_nodes  # (ndm,)

            k0, k1, k2 = gp.kappa[0], gp.kappa[1], gp.kappa[2]
            a_vel_g = np.zeros(ndm, dtype=np.float64)
            a_disp_g = np.zeros(ndm, dtype=np.float64)
            D_vel_g = 0.0
            D_disp_g = 0.0

            for i in range(ndm):
                vg0, dg0 = gp.vel_g0[i].evaluate_trial(float(v_g[i]), dt)
                vg1, dg1 = gp.vel_g1[i].evaluate_trial(float(v_g[i]), dt)
                vg2, dg2 = gp.vel_g2[i].evaluate_trial(float(v_g[i]), dt)
                a_vel_g[i] = (k1 * k2) * vg0 + (k0 * k2) * vg1 + (k0 * k1) * vg2
                D_vel_g = (k1 * k2) * dg0 + (k0 * k2) * dg1 + (k0 * k1) * dg2

                ug10, d10 = gp.disp_g1_g0[i].evaluate_trial(float(u_g[i]), dt)
                ug20, d20 = gp.disp_g2_g0[i].evaluate_trial(float(u_g[i]), dt)
                ug21, d21 = gp.disp_g2_g1[i].evaluate_trial(float(u_g[i]), dt)
                uphi, dphi = gp.disp_phi2_g1_g0[i].evaluate_trial(float(u_g[i]), dt)
                a_disp_g[i] = k2 * ug10 + k1 * ug20 + k0 * ug21 + uphi
                D_disp_g = k2 * d10 + k1 * d20 + k0 * d21 + dphi

            rho_w = self.rho * w
            f_vel += rho_w * np.outer(gp.N, a_vel_g)
            f_disp += rho_w * np.outer(gp.N, a_disp_g)
            NN = np.outer(gp.N, gp.N)
            C_scalar += (rho_w * D_vel_g) * NN
            K_inert_scalar += (rho_w * D_disp_g) * NN

        # Rotate local stiffness force and tangent back to global coordinates:
        # f_{K, a} = R^T f'_{K, a}
        f_K_global = (f_K_local @ self.R).reshape(n_dof)

        K_eff = np.zeros((n_dof, n_dof), dtype=np.float64)
        C_eff = np.zeros((n_dof, n_dof), dtype=np.float64)
        for a in range(n_nodes):
            for b in range(n_nodes):
                block_local = K_stiff_local[a, :, b, :]
                block_global = self.R.T @ block_local @ self.R
                K_eff[a * ndm : (a + 1) * ndm, b * ndm : (b + 1) * ndm] = block_global
                c_ab = C_scalar[a, b]
                k_in_ab = K_inert_scalar[a, b]
                for i in range(ndm):
                    K_eff[a * ndm + i, b * ndm + i] += k_in_ab
                    C_eff[a * ndm + i, b * ndm + i] = c_ab

        f_acc = self.M_eff @ a_vec
        R_int = f_K_global + f_acc + f_vel.reshape(n_dof) + f_disp.reshape(n_dof)
        K_T_dyn = c1 * K_eff + c2 * C_eff + c3 * self.M_eff
        P_ext = float(np.dot(v_vec, R_int))
        return R_int, K_eff, C_eff, self.M_eff.copy(), K_T_dyn, P_ext

    def commit(self) -> None:
        for gp in self.gps:
            gp.commit()


# ============================================================================ #
# 8 & 9. Canonical Test Cases (18 Cases Covering All 9 Topologies)             #
# ============================================================================ #

@dataclass
class CaseSpec:
    name: str
    topo: str
    lumped: int
    E: float
    nu: float
    rho: float
    pml: PMLParams
    geom: GeomSpec
    coords: np.ndarray
    c1: float
    c2: float
    c3: float
    dts: list[float]


def make_reference_coords(topo: str) -> np.ndarray:
    """Generate non-trivial (skewed/perturbed) physical coordinates with det(J) > 0."""
    if topo == "T3":
        return np.array([
            [2.0, 0.1],
            [3.2, 0.25],
            [2.15, 1.35],
        ], dtype=np.float64)
    elif topo == "Q4":
        return np.array([
            [2.0, 2.0],
            [3.3, 2.1],
            [3.4, 3.5],
            [1.9, 3.2],
        ], dtype=np.float64)
    elif topo == "T6":
        # Node 0=(1,0), Node 1=(0,1), Node 2=(0,0) in (s,t)
        v0 = np.array([3.2, 2.1], dtype=np.float64)
        v1 = np.array([2.1, 3.3], dtype=np.float64)
        v2 = np.array([2.0, 2.0], dtype=np.float64)
        m3 = 0.5 * (v0 + v1) + np.array([0.03, 0.02])
        m4 = 0.5 * (v1 + v2) + np.array([-0.02, 0.01])
        m5 = 0.5 * (v2 + v0) + np.array([0.01, -0.02])
        return np.vstack([v0, v1, v2, m3, m4, m5])
    elif topo == "BT6":
        # Node 0=(0,0), Node 1=(1,0), Node 2=(0,1)
        v0 = np.array([2.0, 2.0], dtype=np.float64)
        v1 = np.array([3.25, 2.15], dtype=np.float64)
        v2 = np.array([2.1, 3.4], dtype=np.float64)
        m3 = 0.5 * (v0 + v1)
        m4 = 0.5 * (v1 + v2)
        m5 = 0.5 * (v2 + v0)
        return np.vstack([v0, v1, v2, m3, m4, m5])
    elif topo == "TET4":
        return np.array([
            [0.1, -0.2, 2.0],
            [1.3, -0.1, 2.1],
            [0.2, 1.1, 2.05],
            [0.15, -0.05, 3.25],
        ], dtype=np.float64)
    elif topo == "H8":
        base = np.array([
            [2.0, 2.0, 2.0],
            [3.2, 2.05, 2.0],
            [3.25, 3.3, 2.05],
            [1.95, 3.2, 2.0],
            [2.05, 2.0, 3.3],
            [3.3, 2.1, 3.25],
            [3.35, 3.35, 3.4],
            [2.0, 3.25, 3.3],
        ], dtype=np.float64)
        return base
    elif topo in ("TET10", "BTET10"):
        # Node 0=(1,0,0), Node 1=(0,1,0), Node 2=(0,0,1), Node 3=(0,0,0)
        v0 = np.array([3.3, 2.1, 2.05], dtype=np.float64)
        v1 = np.array([2.05, 3.25, 2.1], dtype=np.float64)
        v2 = np.array([2.1, 2.05, 3.35], dtype=np.float64)
        v3 = np.array([2.0, 2.0, 2.0], dtype=np.float64)
        m4 = 0.5 * (v0 + v1)
        m5 = 0.5 * (v1 + v2)
        m6 = 0.5 * (v2 + v0)
        m7 = 0.5 * (v0 + v3)
        m8 = 0.5 * (v2 + v3)
        m9 = 0.5 * (v1 + v3)
        return np.vstack([v0, v1, v2, v3, m4, m5, m6, m7, m8, m9])
    elif topo == "H20":
        corners = make_reference_coords("H8")
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]
        mids = np.array([0.5 * (corners[a] + corners[b]) for a, b in edges], dtype=np.float64)
        return np.vstack([corners, mids])
    raise ValueError(f"Unknown topo {topo}")


def build_test_cases() -> list[CaseSpec]:
    pml_std = PMLParams(m_deg=2.0, R0=1.0e-4, kappa_max=5.0, alpha_max=math.pi * 15.0)
    pml_cubic = PMLParams(m_deg=3.0, R0=1.0e-5, kappa_max=4.0, alpha_max=math.pi * 12.0)

    c37 = math.cos(math.radians(37.0))
    s37 = math.sin(math.radians(37.0))
    plane_2d = GeomSpec("plane", (1.5, 1.0, c37, s37, 2.5))
    plane_3d = GeomSpec("plane", (1.0, 1.0, 1.0, 0.6, -0.5, 0.6245, 3.0))

    box_2d_corner = GeomSpec("box", (-2.0, 1.8, -2.0, 1.8, 2.0, 2.0))
    box_2d_xface = GeomSpec("box", (-2.0, 1.8, -5.0, 5.0, 2.0, 2.0))
    box_2d_inside = GeomSpec("box", (-10.0, 10.0, -10.0, 10.0, 2.0, 2.0))

    box_3d_zface = GeomSpec("box", (-5.0, 5.0, -5.0, 5.0, -2.0, 1.8, 2.0, 2.0, 2.0))
    box_3d_xz_edge = GeomSpec("box", (-2.0, 1.8, -5.0, 5.0, -2.0, 1.8, 2.0, 2.0, 2.0))
    box_3d_corner = GeomSpec("box", (-2.0, 1.8, -2.0, 1.8, -2.0, 1.8, 2.0, 2.0, 2.0))

    cases = [
        CaseSpec(
            name="T3_box_face_x_cons",
            topo="T3",
            lumped=0,
            E=8.0e7,
            nu=0.25,
            rho=2000.0,
            pml=pml_std,
            geom=box_2d_xface,
            coords=make_reference_coords("T3"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="T3_plane_oblique_lumped",
            topo="T3",
            lumped=1,
            E=1.2e8,
            nu=0.30,
            rho=2100.0,
            pml=pml_std,
            geom=plane_2d,
            coords=make_reference_coords("T3"),
            c1=1.0,
            c2=800.0,
            c3=6.4e5,
            dts=[0.0025, 0.0020, 0.0015],  # variable dt across steps
        ),
        CaseSpec(
            name="Q4_box_corner_xy_cons",
            topo="Q4",
            lumped=0,
            E=2.0e8,
            nu=0.28,
            rho=2200.0,
            pml=pml_std,
            geom=box_2d_corner,
            coords=make_reference_coords("Q4"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="Q4_box_interior_zero_d_lumped",
            topo="Q4",
            lumped=1,
            E=1.5e8,
            nu=0.22,
            rho=1950.0,
            pml=pml_std,
            geom=box_2d_inside,
            coords=make_reference_coords("Q4"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="T6_box_corner_xy_cons",
            topo="T6",
            lumped=0,
            E=2.5e8,
            nu=0.32,
            rho=2300.0,
            pml=pml_cubic,
            geom=box_2d_corner,
            coords=make_reference_coords("T6"),
            c1=1.0,
            c2=1250.0,
            c3=1.5625e6,
            dts=[0.0016, 0.0016, 0.0016],
        ),
        CaseSpec(
            name="T6_plane_oblique_lumped",
            topo="T6",
            lumped=1,
            E=1.8e8,
            nu=0.26,
            rho=2050.0,
            pml=pml_std,
            geom=plane_2d,
            coords=make_reference_coords("T6"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="BT6_box_corner_xy_cons",
            topo="BT6",
            lumped=0,
            E=2.1e8,
            nu=0.29,
            rho=2150.0,
            pml=pml_std,
            geom=box_2d_corner,
            coords=make_reference_coords("BT6"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="BT6_plane_oblique_lumped",
            topo="BT6",
            lumped=1,
            E=1.6e8,
            nu=0.24,
            rho=1980.0,
            pml=pml_std,
            geom=plane_2d,
            coords=make_reference_coords("BT6"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="TET4_box_face_z_cons",
            topo="TET4",
            lumped=0,
            E=3.0e8,
            nu=0.25,
            rho=2250.0,
            pml=pml_std,
            geom=box_3d_zface,
            coords=make_reference_coords("TET4"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="TET4_plane_oblique3d_lumped",
            topo="TET4",
            lumped=1,
            E=2.4e8,
            nu=0.30,
            rho=2100.0,
            pml=pml_std,
            geom=plane_3d,
            coords=make_reference_coords("TET4"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="H8_box_edge_xz_cons",
            topo="H8",
            lumped=0,
            E=3.5e8,
            nu=0.27,
            rho=2350.0,
            pml=pml_std,
            geom=box_3d_xz_edge,
            coords=make_reference_coords("H8"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="H8_box_corner_xyz_lumped",
            topo="H8",
            lumped=1,
            E=4.0e8,
            nu=0.33,
            rho=2400.0,
            pml=pml_cubic,
            geom=box_3d_corner,
            coords=make_reference_coords("H8"),
            c1=1.0,
            c2=900.0,
            c3=8.1e5,
            dts=[0.0015, 0.0020, 0.0025],  # variable dt across 3D corner cascade
        ),
        CaseSpec(
            name="TET10_box_corner_xyz_cons",
            topo="TET10",
            lumped=0,
            E=2.8e8,
            nu=0.26,
            rho=2200.0,
            pml=pml_std,
            geom=box_3d_corner,
            coords=make_reference_coords("TET10"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="TET10_plane_oblique3d_lumped",
            topo="TET10",
            lumped=1,
            E=3.2e8,
            nu=0.31,
            rho=2300.0,
            pml=pml_std,
            geom=plane_3d,
            coords=make_reference_coords("TET10"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="BTET10_box_corner_xyz_cons",
            topo="BTET10",
            lumped=0,
            E=2.6e8,
            nu=0.28,
            rho=2180.0,
            pml=pml_std,
            geom=box_3d_corner,
            coords=make_reference_coords("BTET10"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="BTET10_plane_oblique3d_lumped",
            topo="BTET10",
            lumped=1,
            E=3.1e8,
            nu=0.29,
            rho=2260.0,
            pml=pml_std,
            geom=plane_3d,
            coords=make_reference_coords("BTET10"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="H20_box_corner_xyz_cons",
            topo="H20",
            lumped=0,
            E=4.5e8,
            nu=0.30,
            rho=2380.0,
            pml=pml_std,
            geom=box_3d_corner,
            coords=make_reference_coords("H20"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        CaseSpec(
            name="H20_plane_oblique3d_lumped",
            topo="H20",
            lumped=1,
            E=3.8e8,
            nu=0.25,
            rho=2320.0,
            pml=pml_std,
            geom=plane_3d,
            coords=make_reference_coords("H20"),
            c1=1.0,
            c2=1000.0,
            c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
    ]
    return cases


# ============================================================================ #
# 10. Biot u-p Coupled Poroelastic CFS-PML Reference Oracle (Phase P0-UP)      #
# ============================================================================ #

N_VERT = {
    "T3": 3,
    "Q4": 4,
    "T6": 3,
    "BT6": 3,
    "TET4": 4,
    "H8": 8,
    "TET10": 4,
    "BTET10": 4,
    "H20": 8,
}

VERTEX_EDGES = {
    "T3": [(0, 1), (1, 2), (2, 0)],
    "T6": [(0, 1), (1, 2), (2, 0)],
    "BT6": [(0, 1), (1, 2), (2, 0)],
    "Q4": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "TET4": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "TET10": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "BTET10": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "H8": [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ],
    "H20": [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ],
}


@dataclass
class BiotParams:
    Kf: float
    poro: float
    rhoF: float
    biotAlpha: float
    Ks: float
    pOrder: int       # 0 = equal, 1 = linear / Taylor-Hood
    formulation: int  # 0 = std, 1 = bbar
    stabMode: int     # 0 = off, 1 = auto, 2 = manual
    stabValue: float
    kbar: tuple[float, ...]


def eval_pressure_ref_basis(topo: str, p_order: int, xi: tuple[float, ...]) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate reference pressure basis (N_p, dN_p_dxi) at parametric point xi."""
    _, _, shape_fn, _ = TOPOLOGIES[topo]
    if p_order == 0:
        return shape_fn(xi)
    if topo in ("T3", "Q4", "TET4", "H8"):
        return shape_fn(xi)
    elif topo == "T6":
        s, t = xi
        N_p = np.array([s, t, 1.0 - s - t], dtype=np.float64)
        dN_p_dxi = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [-1.0, -1.0],
        ], dtype=np.float64)
        return N_p, dN_p_dxi
    elif topo == "BT6":
        l1, l2 = xi
        N_p = np.array([1.0 - l1 - l2, l1, l2], dtype=np.float64)
        dN_p_dxi = np.array([
            [-1.0, -1.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ], dtype=np.float64)
        return N_p, dN_p_dxi
    elif topo in ("TET10", "BTET10"):
        z1, z2, z3 = xi
        N_p = np.array([z1, z2, z3, 1.0 - z1 - z2 - z3], dtype=np.float64)
        dN_p_dxi = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ], dtype=np.float64)
        return N_p, dN_p_dxi
    elif topo == "H20":
        return shape_h8(xi)
    raise ValueError(f"Unsupported topology {topo} for p_order={p_order}")


def eval_gp_kinematics_up(
    topo: str, coords: np.ndarray, p_order: int
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]]:
    """Return list of (N_u, dN_u_dx, N_p, dN_p_dx, w_vol, x_g) for each Gauss point."""
    ndm, n_nodes, shape_fn, gp_fn = TOPOLOGIES[topo]
    assert coords.shape == (n_nodes, ndm)
    out = []
    for xi, w_g in gp_fn():
        N_u, dN_u_dxi = shape_fn(xi)
        J = dN_u_dxi.T @ coords
        detJ = float(np.linalg.det(J))
        if detJ <= 0.0:
            raise ValueError(f"Non-positive Jacobian determinant {detJ} in {topo}")
        dN_u_dx = np.linalg.solve(J, dN_u_dxi.T).T
        N_p, dN_p_dxi = eval_pressure_ref_basis(topo, p_order, xi)
        dN_p_dx = np.linalg.solve(J, dN_p_dxi.T).T
        w_vol = w_g * detJ
        x_g = N_u @ coords
        out.append((N_u, dN_u_dx, N_p, dN_p_dx, w_vol, x_g))
    return out


def compute_dof_map_up(topo: str, p_order: int) -> tuple[int, int, list[int], list[int]]:
    """Return (n_dof, n_p, u_off, p_off) for per-node interleaved [u..., p?] ordering."""
    ndm, n_nodes, _, _ = TOPOLOGIES[topo]
    n_vert = N_VERT[topo]
    n_p = n_nodes if p_order == 0 else n_vert
    p_carrier = [1] * n_nodes if p_order == 0 else ([1] * n_vert + [0] * (n_nodes - n_vert))
    u_off = [0] * n_nodes
    p_off = [-1] * n_nodes
    slot = 0
    for a in range(n_nodes):
        u_off[a] = slot
        slot += ndm
        if p_carrier[a] == 1:
            p_off[a] = slot
            slot += 1
    assert slot == n_nodes * ndm + n_p
    return slot, n_p, u_off, p_off


def compute_elem_size_h(topo: str, coords: np.ndarray) -> float:
    return max(
        float(np.linalg.norm(coords[a] - coords[b]))
        for a, b in VERTEX_EDGES[topo]
    )


def make_coupling_cascade(
    ndm: int, j: int, d: np.ndarray, kappa: np.ndarray, alpha: np.ndarray
) -> FilterCascade:
    """Create the contractual stage cascade for E_j(w) = prod_{m != j} s_m(w)."""
    if ndm == 2:
        if j == 0:
            return FilterCascade([make_s_filter(d[1], kappa[1], alpha[1])])
        elif j == 1:
            return FilterCascade([make_s_filter(d[0], kappa[0], alpha[0])])
        else:
            raise ValueError(f"Invalid j={j} for ndm=2")
    else:
        if j == 0:
            return FilterCascade([
                make_s_filter(d[1], kappa[1], alpha[1]),
                make_s_filter(d[2], kappa[2], alpha[2]),
            ])
        elif j == 1:
            return FilterCascade([
                make_s_filter(d[0], kappa[0], alpha[0]),
                make_s_filter(d[2], kappa[2], alpha[2]),
            ])
        elif j == 2:
            return FilterCascade([
                make_s_filter(d[0], kappa[0], alpha[0]),
                make_s_filter(d[1], kappa[1], alpha[1]),
            ])
        else:
            raise ValueError(f"Invalid j={j} for ndm=3")


def make_storage_cascade(
    ndm: int, d: np.ndarray, kappa: np.ndarray, alpha: np.ndarray
) -> FilterCascade:
    """Create the contractual stage cascade for J_s(w) = prod_{m=0}^{ndm-1} s_m(w)."""
    if ndm == 2:
        return FilterCascade([
            make_s_filter(d[0], kappa[0], alpha[0]),
            make_s_filter(d[1], kappa[1], alpha[1]),
        ])
    else:
        return FilterCascade([
            make_s_filter(d[0], kappa[0], alpha[0]),
            make_s_filter(d[1], kappa[1], alpha[1]),
            make_s_filter(d[2], kappa[2], alpha[2]),
        ])


class GaussPointStateUP:
    """Holds all filter cascades for one Gauss point of a Biot u-p CFS-PML element."""

    def __init__(
        self,
        ndm: int,
        N_u: np.ndarray,
        dN_u_dx: np.ndarray,
        N_p: np.ndarray,
        dN_p_dx: np.ndarray,
        w_vol: float,
        x_g: np.ndarray,
        R: np.ndarray,
        d: np.ndarray,
        kappa: np.ndarray,
        alpha: np.ndarray,
    ) -> None:
        self.ndm = ndm
        self.N_u = N_u
        self.dN_u_dx = dN_u_dx
        self.Phi_u_local = dN_u_dx @ R.T  # (n_nodes, ndm)
        self.N_p = N_p
        self.dN_p_dx = dN_p_dx
        self.Phi_p_local = dN_p_dx @ R.T  # (n_p, ndm)
        self.w_vol = w_vol
        self.x_g = x_g
        self.d = d
        self.kappa = kappa
        self.alpha = alpha

        # Solid stiffness cascades F_{jl}[G'_{kl}]
        self.F_u_cascades = {}
        # Darcy cascades F_{jl}^{(p)}[g'_{p, l}]
        self.F_p_cascades = {}
        for j in range(ndm):
            for l in range(ndm):
                self.F_u_cascades[(j, l)] = [
                    make_stiffness_cascade(ndm, j, l, d, kappa, alpha)
                    for _ in range(ndm)
                ]
                self.F_p_cascades[(j, l)] = make_stiffness_cascade(ndm, j, l, d, kappa, alpha)

        # Coupling cascades E_j^{(p)}[p_g] and E_j^{(dot_u)}[dot_eps'_{jj, g}]
        self.E_p_cascades = [
            make_coupling_cascade(ndm, j, d, kappa, alpha) for j in range(ndm)
        ]
        self.E_udot_cascades = [
            make_coupling_cascade(ndm, j, d, kappa, alpha) for j in range(ndm)
        ]

        # Storage cascade J_s^{(p)}[dot_p_g]
        self.J_s_cascade = make_storage_cascade(ndm, d, kappa, alpha)

        # Solid inertial velocity and displacement cascades
        self.vel_g0 = [make_g_filter(d[0], kappa[0], alpha[0]) for _ in range(ndm)]
        self.vel_g1 = [make_g_filter(d[1], kappa[1], alpha[1]) for _ in range(ndm)]
        self.vel_g2 = [make_g_filter(d[2], kappa[2], alpha[2]) for _ in range(ndm)]

        self.disp_g1_g0 = [
            FilterCascade([
                make_g_filter(d[0], kappa[0], alpha[0]),
                make_g_filter(d[1], kappa[1], alpha[1]),
            ])
            for _ in range(ndm)
        ]
        self.disp_g2_g0 = [
            FilterCascade([
                make_g_filter(d[0], kappa[0], alpha[0]),
                make_g_filter(d[2], kappa[2], alpha[2]),
            ])
            for _ in range(ndm)
        ]
        self.disp_g2_g1 = [
            FilterCascade([
                make_g_filter(d[1], kappa[1], alpha[1]),
                make_g_filter(d[2], kappa[2], alpha[2]),
            ])
            for _ in range(ndm)
        ]
        self.disp_phi2_g1_g0 = [
            FilterCascade([
                make_g_filter(d[0], kappa[0], alpha[0]),
                make_g_filter(d[1], kappa[1], alpha[1]),
                make_phi_filter(d[2], kappa[2], alpha[2]),
            ])
            for _ in range(ndm)
        ]

    def commit(self) -> None:
        for cascades in self.F_u_cascades.values():
            for c in cascades:
                c.commit()
        for c in self.F_p_cascades.values():
            c.commit()
        for c in self.E_p_cascades + self.E_udot_cascades:
            c.commit()
        self.J_s_cascade.commit()
        for f in self.vel_g0 + self.vel_g1 + self.vel_g2:
            f.commit()
        for c in (
            self.disp_g1_g0
            + self.disp_g2_g0
            + self.disp_g2_g1
            + self.disp_phi2_g1_g0
        ):
            c.commit()


class ReferenceCPMLUPElement:
    """Standalone reference implementation of a coupled Biot u-p LadrunoCPML element."""

    def __init__(
        self,
        topo: str,
        coords: np.ndarray,
        E: float,
        nu: float,
        rho: float,
        biot: BiotParams,
        pml: PMLParams,
        geom: GeomSpec,
        lumped: int,
    ) -> None:
        self.topo = topo
        self.ndm, self.n_nodes, _, _ = TOPOLOGIES[topo]
        self.coords = np.array(coords, dtype=np.float64).reshape(self.n_nodes, self.ndm)
        self.E = float(E)
        self.nu = float(nu)
        self.rho = float(rho)
        self.biot = biot
        self.pml = pml
        self.geom = geom
        self.lumped = int(lumped)

        self.n_dof, self.n_p, self.u_off, self.p_off = compute_dof_map_up(
            self.topo, self.biot.pOrder
        )

        # Isotropic solid skeleton moduli & 1/Qbar + undrained P1 wave speed
        self.lam, self.mu, _, self.C_tensor = isotropic_elasticity(
            self.E, self.nu, self.rho, self.ndm
        )
        self.one_over_Qbar = (self.biot.poro / self.biot.Kf) + (
            (self.biot.biotAlpha - self.biot.poro) / self.biot.Ks
            if self.biot.Ks > 0.0
            else 0.0
        )
        M_undrained = self.lam + 2.0 * self.mu + (
            (self.biot.biotAlpha ** 2) / self.one_over_Qbar
            if self.one_over_Qbar > 0.0
            else 0.0
        )
        self.Vp = math.sqrt(M_undrained / self.rho)

        self.R = compute_rotation_matrix(self.ndm, self.geom)
        # Local permeability tensor kbar'_{jl} = sum_m R_{jm} kbar_m R_{lm}
        kbar_arr = np.array(self.biot.kbar, dtype=np.float64)
        assert kbar_arr.shape == (self.ndm,)
        self.kbar_local = (self.R * kbar_arr[None, :]) @ self.R.T

        # Gauss point kinematics & filter cascades
        gp_kin = eval_gp_kinematics_up(self.topo, self.coords, self.biot.pOrder)
        self.gps: list[GaussPointStateUP] = []
        for N_u, dN_u_dx, N_p, dN_p_dx, w_vol, x_g in gp_kin:
            d, kappa, alpha = evaluate_pml_profile_at_gp(x_g, self.Vp, self.pml, self.geom)
            self.gps.append(
                GaussPointStateUP(
                    ndm=self.ndm,
                    N_u=N_u,
                    dN_u_dx=dN_u_dx,
                    N_p=N_p,
                    dN_p_dx=dN_p_dx,
                    w_vol=w_vol,
                    x_g=x_g,
                    R=self.R,
                    d=d,
                    kappa=kappa,
                    alpha=alpha,
                )
            )

        # Element volume and mean local solid gradient for B-bar
        self.V_e = sum(gp.w_vol for gp in self.gps)
        self.Phi_u_bar = np.zeros((self.n_nodes, self.ndm), dtype=np.float64)
        for gp in self.gps:
            self.Phi_u_bar += (gp.w_vol / self.V_e) * gp.Phi_u_local

        # Precompute B'_{kl, (a, m), g} and tilde_Phi'^u_{a, j, g} for each GP
        self.B_local_gps: list[np.ndarray] = []
        self.Phi_u_tilde_gps: list[np.ndarray] = []
        for gp in self.gps:
            Phi_g = gp.Phi_u_local
            B_g = np.zeros((self.ndm, self.ndm, self.n_nodes, self.ndm), dtype=np.float64)
            if self.biot.formulation == 0:
                for k in range(self.ndm):
                    for l in range(self.ndm):
                        B_g[k, l, :, k] = Phi_g[:, l]
                Phi_tilde_g = Phi_g.copy()
            else:
                diff_Phi = (self.Phi_u_bar - Phi_g) / float(self.ndm)
                for k in range(self.ndm):
                    for l in range(self.ndm):
                        B_g[k, l, :, k] += Phi_g[:, l]
                        if k == l:
                            for m in range(self.ndm):
                                B_g[k, l, :, m] += diff_Phi[:, m]
                Phi_tilde_g = self.Phi_u_bar.copy()
            self.B_local_gps.append(B_g)
            self.Phi_u_tilde_gps.append(Phi_tilde_g)

        # McGann stabilization matrix H_tilde (n_p x n_p)
        self.elem_h = compute_elem_size_h(self.topo, self.coords)
        if self.biot.pOrder == 1 or self.biot.stabMode == 0:
            self.alpha_stab = 0.0
        elif self.biot.stabMode == 1:
            self.alpha_stab = self.biot.stabValue * (self.elem_h ** 2) / (self.lam + 2.0 * self.mu)
        elif self.biot.stabMode == 2:
            self.alpha_stab = self.biot.stabValue
        else:
            raise ValueError(f"Invalid stabMode={self.biot.stabMode}")

        self.H_tilde = np.zeros((self.n_p, self.n_p), dtype=np.float64)
        if self.alpha_stab != 0.0:
            for gp in self.gps:
                self.H_tilde += (gp.w_vol * self.alpha_stab) * (gp.dN_p_dx @ gp.dN_p_dx.T)

        # Precompute static effective mass matrix M_eff (n_dof x n_dof; p-rows/cols are 0)
        M_scalar_cons = np.zeros((self.n_nodes, self.n_nodes), dtype=np.float64)
        for gp in self.gps:
            D_acc = float(gp.kappa[0] * gp.kappa[1] * gp.kappa[2])
            M_scalar_cons += (self.rho * gp.w_vol * D_acc) * np.outer(gp.N_u, gp.N_u)

        self.M_eff = np.zeros((self.n_dof, self.n_dof), dtype=np.float64)
        if self.lumped == 0:
            for a in range(self.n_nodes):
                for b in range(self.n_nodes):
                    val = M_scalar_cons[a, b]
                    for d in range(self.ndm):
                        self.M_eff[self.u_off[a] + d, self.u_off[b] + d] = val
        else:
            M_tot = float(np.sum(M_scalar_cons))
            D_sum = float(np.trace(M_scalar_cons))
            for a in range(self.n_nodes):
                m_a = M_tot * (M_scalar_cons[a, a] / D_sum)
                for d in range(self.ndm):
                    self.M_eff[self.u_off[a] + d, self.u_off[a] + d] = m_a

    def evaluate_trial(
        self,
        U_vec: np.ndarray,
        V_vec: np.ndarray,
        A_vec: np.ndarray,
        dt: float,
        c1: float,
        c2: float,
        c3: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        """Evaluate R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext at trial state (U, V, A)."""
        ndm = self.ndm
        n_nodes = self.n_nodes
        n_p = self.n_p
        n_dof = self.n_dof

        # Extract nodal solid u, v and pressure p, dot_p
        u_nodes = np.zeros((n_nodes, ndm), dtype=np.float64)
        v_nodes = np.zeros((n_nodes, ndm), dtype=np.float64)
        for a in range(n_nodes):
            u_nodes[a, :] = U_vec[self.u_off[a] : self.u_off[a] + ndm]
            v_nodes[a, :] = V_vec[self.u_off[a] : self.u_off[a] + ndm]

        p_vec = np.zeros(n_p, dtype=np.float64)
        pdot_vec = np.zeros(n_p, dtype=np.float64)
        for b in range(n_p):
            p_vec[b] = U_vec[self.p_off[b]]
            pdot_vec[b] = V_vec[self.p_off[b]]

        # Local solid displacements and velocities: u'_a = R u_a, v'_a = R v_a
        u_local = u_nodes @ self.R.T
        v_local = v_nodes @ self.R.T

        # Accumulators in local/block coordinates
        f_K_local = np.zeros((n_nodes, ndm), dtype=np.float64)
        K_stiff_local = np.zeros((n_nodes, ndm, n_nodes, ndm), dtype=np.float64)

        f_u_p_local = np.zeros((n_nodes, ndm), dtype=np.float64)
        Q_eff_local = np.zeros((n_nodes, ndm, n_p), dtype=np.float64)

        f_p_H = np.zeros(n_p, dtype=np.float64)
        H_eff = np.zeros((n_p, n_p), dtype=np.float64)

        f_p_udot = np.zeros(n_p, dtype=np.float64)
        f_p_S = np.zeros(n_p, dtype=np.float64)
        S_eff = np.zeros((n_p, n_p), dtype=np.float64)

        f_vel = np.zeros((n_nodes, ndm), dtype=np.float64)
        f_disp = np.zeros((n_nodes, ndm), dtype=np.float64)
        C_scalar = np.zeros((n_nodes, n_nodes), dtype=np.float64)
        K_inert_scalar = np.zeros((n_nodes, n_nodes), dtype=np.float64)

        alpha_B = self.biot.biotAlpha

        for g_idx, gp in enumerate(self.gps):
            w = gp.w_vol
            B_g = self.B_local_gps[g_idx]           # (ndm, ndm, n_nodes, ndm)
            Phi_u_tilde = self.Phi_u_tilde_gps[g_idx]  # (n_nodes, ndm)
            Phi_p = gp.Phi_p_local                  # (n_p, ndm)
            N_p = gp.N_p                            # (n_p,)

            # 1. Solid stiffness filter F_{jl} on G'_{kl, g}
            G_local = np.einsum("klam,am->kl", B_g, u_local)  # (ndm, ndm)
            G_tilde = np.zeros((ndm, ndm, ndm), dtype=np.float64)  # [j, l, k]
            D_F = np.zeros((ndm, ndm), dtype=np.float64)
            for j in range(ndm):
                for l in range(ndm):
                    cascades_jl = gp.F_u_cascades[(j, l)]
                    for k in range(ndm):
                        y_val, D_val = cascades_jl[k].evaluate_trial(float(G_local[k, l]), dt)
                        G_tilde[j, l, k] = y_val
                        D_F[j, l] = D_val

            for i in range(ndm):
                for j in range(ndm):
                    for k in range(ndm):
                        for l in range(ndm):
                            c_ijkl = self.C_tensor[i, j, k, l]
                            if c_ijkl != 0.0:
                                f_K_local += (w * c_ijkl * G_tilde[j, l, k]) * B_g[i, j, :, :]
                                K_stiff_local += (w * c_ijkl * D_F[j, l]) * np.einsum(
                                    "ar,bm->arbm", B_g[i, j, :, :], B_g[k, l, :, :]
                                )

            # 2. Darcy filter F_{jl}^{(p)} on g'_{p, l}
            g_p_local = p_vec @ Phi_p  # (ndm,)
            g_p_tilde = np.zeros((ndm, ndm), dtype=np.float64)
            D_F_p = np.zeros((ndm, ndm), dtype=np.float64)
            for j in range(ndm):
                for l in range(ndm):
                    y_val, D_val = gp.F_p_cascades[(j, l)].evaluate_trial(float(g_p_local[l]), dt)
                    g_p_tilde[j, l] = y_val
                    D_F_p[j, l] = D_val

            for j in range(ndm):
                flux_j = 0.0
                for l in range(ndm):
                    k_jl = self.kbar_local[j, l]
                    flux_j += k_jl * g_p_tilde[j, l]
                    H_eff += (w * k_jl * D_F_p[j, l]) * np.outer(Phi_p[:, j], Phi_p[:, l])
                f_p_H += (w * flux_j) * Phi_p[:, j]

            # 3. Coupling filter E_j on u-rows (E_j^{(p)}) and p-rows (E_j^{(dot_u)})
            p_g = float(np.dot(N_p, p_vec))
            eps_dot_sum = 0.0
            for j in range(ndm):
                p_tilde_j, D_Ej = gp.E_p_cascades[j].evaluate_trial(p_g, dt)
                f_u_p_local[:, j] -= (w * alpha_B * p_tilde_j) * Phi_u_tilde[:, j]
                Q_eff_local[:, j, :] += (w * alpha_B * D_Ej) * np.outer(Phi_u_tilde[:, j], N_p)

                eps_dot_jj = float(np.dot(v_local[:, j], Phi_u_tilde[:, j]))
                eps_dot_tilde_jj, _ = gp.E_udot_cascades[j].evaluate_trial(eps_dot_jj, dt)
                eps_dot_sum += eps_dot_tilde_jj

            f_p_udot += (w * alpha_B * eps_dot_sum) * N_p

            # 4. Storage filter J_s^{(p)} on p_dot_g
            pdot_g = float(np.dot(N_p, pdot_vec))
            pdot_tilde_g, D_Js = gp.J_s_cascade.evaluate_trial(pdot_g, dt)
            f_p_S += (w * self.one_over_Qbar * pdot_tilde_g) * N_p
            S_eff += (w * self.one_over_Qbar * D_Js) * np.outer(N_p, N_p)

            # 5. Solid inertial velocity and displacement filters
            u_g = gp.N_u @ u_nodes
            v_g = gp.N_u @ v_nodes
            k0, k1, k2 = gp.kappa[0], gp.kappa[1], gp.kappa[2]
            a_vel_g = np.zeros(ndm, dtype=np.float64)
            a_disp_g = np.zeros(ndm, dtype=np.float64)
            D_vel_g = 0.0
            D_disp_g = 0.0
            for i in range(ndm):
                vg0, dg0 = gp.vel_g0[i].evaluate_trial(float(v_g[i]), dt)
                vg1, dg1 = gp.vel_g1[i].evaluate_trial(float(v_g[i]), dt)
                vg2, dg2 = gp.vel_g2[i].evaluate_trial(float(v_g[i]), dt)
                a_vel_g[i] = (k1 * k2) * vg0 + (k0 * k2) * vg1 + (k0 * k1) * vg2
                D_vel_g = (k1 * k2) * dg0 + (k0 * k2) * dg1 + (k0 * k1) * dg2

                ug10, d10 = gp.disp_g1_g0[i].evaluate_trial(float(u_g[i]), dt)
                ug20, d20 = gp.disp_g2_g0[i].evaluate_trial(float(u_g[i]), dt)
                ug21, d21 = gp.disp_g2_g1[i].evaluate_trial(float(u_g[i]), dt)
                uphi, dphi = gp.disp_phi2_g1_g0[i].evaluate_trial(float(u_g[i]), dt)
                a_disp_g[i] = k2 * ug10 + k1 * ug20 + k0 * ug21 + uphi
                D_disp_g = k2 * d10 + k1 * d20 + k0 * d21 + dphi

            rho_w = self.rho * w
            f_vel += rho_w * np.outer(gp.N_u, a_vel_g)
            f_disp += rho_w * np.outer(gp.N_u, a_disp_g)
            NN_u = np.outer(gp.N_u, gp.N_u)
            C_scalar += (rho_w * D_vel_g) * NN_u
            K_inert_scalar += (rho_w * D_disp_g) * NN_u

        # Rotate solid local forces and blocks back to global coordinates
        f_K_global = f_K_local @ self.R
        f_u_p_global = f_u_p_local @ self.R
        # Q_eff_global[a, d, b] = sum_j R_{j, d} Q_eff_local[a, j, b]
        Q_eff_global = np.einsum("jd,ajb->adb", self.R, Q_eff_local)

        K_eff = np.zeros((n_dof, n_dof), dtype=np.float64)
        C_eff = np.zeros((n_dof, n_dof), dtype=np.float64)

        # Scatter solid-solid blocks into K_eff and C_eff
        for a in range(n_nodes):
            ua = self.u_off[a]
            for b in range(n_nodes):
                ub = self.u_off[b]
                block_local = K_stiff_local[a, :, b, :]
                block_global = self.R.T @ block_local @ self.R
                K_eff[ua : ua + ndm, ub : ub + ndm] = block_global
                c_ab = C_scalar[a, b]
                k_in_ab = K_inert_scalar[a, b]
                for d in range(ndm):
                    K_eff[ua + d, ub + d] += k_in_ab
                    C_eff[ua + d, ub + d] = c_ab

        # Scatter coupling Q_eff into K_eff (u, p) as -Q_eff and C_eff (p, u) as +Q_eff^T
        for a in range(n_nodes):
            ua = self.u_off[a]
            for d in range(ndm):
                for b in range(n_p):
                    pb = self.p_off[b]
                    q_val = Q_eff_global[a, d, b]
                    K_eff[ua + d, pb] = -q_val
                    C_eff[pb, ua + d] = +q_val

        # Scatter pressure-pressure blocks into K_eff (H_eff) and C_eff (S_eff + H_tilde)
        S_plus_HT = S_eff + self.H_tilde
        for b in range(n_p):
            pb = self.p_off[b]
            for c in range(n_p):
                pc = self.p_off[c]
                K_eff[pb, pc] = H_eff[b, c]
                C_eff[pb, pc] = S_plus_HT[b, c]

        # Assemble internal force vector R_int
        R_int = self.M_eff @ A_vec
        for a in range(n_nodes):
            ua = self.u_off[a]
            R_int[ua : ua + ndm] += (
                f_K_global[a, :]
                + f_u_p_global[a, :]
                + f_vel[a, :]
                + f_disp[a, :]
            )
        f_p_total = f_p_H + f_p_udot + f_p_S + (self.H_tilde @ pdot_vec)
        for b in range(n_p):
            R_int[self.p_off[b]] += f_p_total[b]

        K_T_dyn = c1 * K_eff + c2 * C_eff + c3 * self.M_eff
        P_ext = float(np.dot(V_vec, R_int))
        return R_int, K_eff, C_eff, self.M_eff.copy(), K_T_dyn, P_ext

    def commit(self) -> None:
        for gp in self.gps:
            gp.commit()


# ============================================================================ #
# 11. Canonical Biot u-p Test Cases (14 Cases Covering All 9 Topologies)       #
# ============================================================================ #

@dataclass
class UPCaseSpec:
    name: str
    topo: str
    lumped: int
    E: float
    nu: float
    rho: float
    biot: BiotParams
    pml: PMLParams
    geom: GeomSpec
    coords: np.ndarray
    c1: float
    c2: float
    c3: float
    dts: list[float]


def build_up_test_cases() -> list[UPCaseSpec]:
    pml_std = PMLParams(m_deg=2.0, R0=1.0e-4, kappa_max=5.0, alpha_max=math.pi * 15.0)
    pml_cubic = PMLParams(m_deg=3.0, R0=1.0e-5, kappa_max=4.0, alpha_max=math.pi * 12.0)

    c37 = math.cos(math.radians(37.0))
    s37 = math.sin(math.radians(37.0))
    plane_2d = GeomSpec("plane", (1.5, 1.0, c37, s37, 2.5))
    plane_3d = GeomSpec("plane", (1.0, 1.0, 1.0, 0.6, -0.5, 0.6245, 3.0))

    box_2d_corner = GeomSpec("box", (-2.0, 1.8, -2.0, 1.8, 2.0, 2.0))
    box_2d_xface = GeomSpec("box", (-2.0, 1.8, -5.0, 5.0, 2.0, 2.0))
    box_2d_inside = GeomSpec("box", (-10.0, 10.0, -10.0, 10.0, 2.0, 2.0))

    box_3d_zface = GeomSpec("box", (-5.0, 5.0, -5.0, 5.0, -2.0, 1.8, 2.0, 2.0, 2.0))
    box_3d_corner = GeomSpec("box", (-2.0, 1.8, -2.0, 1.8, -2.0, 1.8, 2.0, 2.0, 2.0))

    cases = [
        # 1. UP_T3_box_face_eq_std
        UPCaseSpec(
            name="UP_T3_box_face_eq_std",
            topo="T3",
            lumped=0,
            E=5.0e7,
            nu=0.25,
            rho=2000.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.30, rhoF=1000.0, biotAlpha=1.0, Ks=0.0,
                pOrder=0, formulation=0, stabMode=1, stabValue=0.25,
                kbar=(1.0e-7, 1.0e-7),
            ),
            pml=pml_std,
            geom=box_2d_xface,
            coords=make_reference_coords("T3"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 2. UP_Q4_box_corner_eq_bbar
        UPCaseSpec(
            name="UP_Q4_box_corner_eq_bbar",
            topo="Q4",
            lumped=1,
            E=8.0e7,
            nu=0.30,
            rho=2100.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.35, rhoF=1000.0, biotAlpha=0.85, Ks=3.6e10,
                pOrder=0, formulation=1, stabMode=1, stabValue=0.25,
                kbar=(2.0e-7, 2.0e-7),
            ),
            pml=pml_std,
            geom=box_2d_corner,
            coords=make_reference_coords("Q4"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 3. UP_Q4_plane_oblique_aniso_eq
        UPCaseSpec(
            name="UP_Q4_plane_oblique_aniso_eq",
            topo="Q4",
            lumped=0,
            E=6.0e7,
            nu=0.28,
            rho=2050.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.28, rhoF=1000.0, biotAlpha=0.90, Ks=3.6e10,
                pOrder=0, formulation=0, stabMode=2, stabValue=1.0e-8,
                kbar=(5.0e-6, 1.2e-6),
            ),
            pml=pml_std,
            geom=plane_2d,
            coords=make_reference_coords("Q4"),
            c1=1.0, c2=800.0, c3=6.4e5,
            dts=[0.002, 0.002, 0.002],
        ),
        # 4. UP_Q4_box_interior_zero_d
        UPCaseSpec(
            name="UP_Q4_box_interior_zero_d",
            topo="Q4",
            lumped=0,
            E=4.0e7,
            nu=0.22,
            rho=1950.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.32, rhoF=1000.0, biotAlpha=1.0, Ks=0.0,
                pOrder=0, formulation=1, stabMode=1, stabValue=0.25,
                kbar=(1.5e-7, 1.5e-7),
            ),
            pml=pml_std,
            geom=box_2d_inside,
            coords=make_reference_coords("Q4"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 5. UP_T6_box_corner_th_std
        UPCaseSpec(
            name="UP_T6_box_corner_th_std",
            topo="T6",
            lumped=0,
            E=1.0e8,
            nu=0.26,
            rho=2150.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.30, rhoF=1000.0, biotAlpha=0.82, Ks=3.6e10,
                pOrder=1, formulation=0, stabMode=0, stabValue=0.0,
                kbar=(8.0e-7, 8.0e-7),
            ),
            pml=pml_cubic,
            geom=box_2d_corner,
            coords=make_reference_coords("T6"),
            c1=1.0, c2=1250.0, c3=1.5625e6,
            dts=[0.0016, 0.0016, 0.0016],
        ),
        # 6. UP_T6_plane_oblique_eq_bbar
        UPCaseSpec(
            name="UP_T6_plane_oblique_eq_bbar",
            topo="T6",
            lumped=1,
            E=9.0e7,
            nu=0.32,
            rho=2080.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.36, rhoF=1000.0, biotAlpha=1.0, Ks=0.0,
                pOrder=0, formulation=1, stabMode=1, stabValue=0.25,
                kbar=(3.0e-7, 1.0e-7),
            ),
            pml=pml_std,
            geom=plane_2d,
            coords=make_reference_coords("T6"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 7. UP_BT6_box_corner_th_std
        UPCaseSpec(
            name="UP_BT6_box_corner_th_std",
            topo="BT6",
            lumped=0,
            E=1.1e8,
            nu=0.29,
            rho=2120.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.27, rhoF=1000.0, biotAlpha=0.80, Ks=3.6e10,
                pOrder=1, formulation=1, stabMode=0, stabValue=0.0,
                kbar=(5.0e-7, 5.0e-7),
            ),
            pml=pml_std,
            geom=box_2d_corner,
            coords=make_reference_coords("BT6"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 8. UP_TET4_box_face_eq_std
        UPCaseSpec(
            name="UP_TET4_box_face_eq_std",
            topo="TET4",
            lumped=0,
            E=1.2e8,
            nu=0.25,
            rho=2150.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.30, rhoF=1000.0, biotAlpha=1.0, Ks=0.0,
                pOrder=0, formulation=0, stabMode=0, stabValue=0.0,
                kbar=(2.0e-7, 2.0e-7, 2.0e-7),
            ),
            pml=pml_std,
            geom=box_3d_zface,
            coords=make_reference_coords("TET4"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 9. UP_H8_box_corner_xyz_eq_bbar
        UPCaseSpec(
            name="UP_H8_box_corner_xyz_eq_bbar",
            topo="H8",
            lumped=1,
            E=1.5e8,
            nu=0.33,
            rho=2200.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.34, rhoF=1000.0, biotAlpha=0.88, Ks=3.6e10,
                pOrder=0, formulation=1, stabMode=1, stabValue=0.25,
                kbar=(4.0e-7, 4.0e-7, 4.0e-7),
            ),
            pml=pml_cubic,
            geom=box_3d_corner,
            coords=make_reference_coords("H8"),
            c1=1.0, c2=900.0, c3=8.1e5,
            dts=[0.0015, 0.0020, 0.0025],
        ),
        # 10. UP_H8_plane_oblique3d_aniso
        UPCaseSpec(
            name="UP_H8_plane_oblique3d_aniso",
            topo="H8",
            lumped=0,
            E=1.3e8,
            nu=0.27,
            rho=2160.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.29, rhoF=1000.0, biotAlpha=0.84, Ks=3.6e10,
                pOrder=0, formulation=0, stabMode=1, stabValue=0.25,
                kbar=(6.0e-6, 2.0e-6, 9.0e-7),
            ),
            pml=pml_std,
            geom=plane_3d,
            coords=make_reference_coords("H8"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 11. UP_TET10_box_corner_xyz_th
        UPCaseSpec(
            name="UP_TET10_box_corner_xyz_th",
            topo="TET10",
            lumped=0,
            E=1.6e8,
            nu=0.26,
            rho=2180.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.25, rhoF=1000.0, biotAlpha=0.78, Ks=3.6e10,
                pOrder=1, formulation=0, stabMode=0, stabValue=0.0,
                kbar=(1.0e-6, 1.0e-6, 1.0e-6),
            ),
            pml=pml_std,
            geom=box_3d_corner,
            coords=make_reference_coords("TET10"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 12. UP_BTET10_box_corner_xyz_th
        UPCaseSpec(
            name="UP_BTET10_box_corner_xyz_th",
            topo="BTET10",
            lumped=1,
            E=1.4e8,
            nu=0.28,
            rho=2140.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.31, rhoF=1000.0, biotAlpha=1.0, Ks=0.0,
                pOrder=1, formulation=0, stabMode=0, stabValue=0.0,
                kbar=(7.0e-7, 7.0e-7, 7.0e-7),
            ),
            pml=pml_std,
            geom=box_3d_corner,
            coords=make_reference_coords("BTET10"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 13. UP_H20_box_corner_xyz_th
        UPCaseSpec(
            name="UP_H20_box_corner_xyz_th",
            topo="H20",
            lumped=0,
            E=2.0e8,
            nu=0.30,
            rho=2200.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.28, rhoF=1000.0, biotAlpha=0.86, Ks=3.6e10,
                pOrder=1, formulation=0, stabMode=0, stabValue=0.0,
                kbar=(3.0e-6, 1.5e-6, 8.0e-7),
            ),
            pml=pml_std,
            geom=box_3d_corner,
            coords=make_reference_coords("H20"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
        # 14. UP_H20_plane_oblique3d_eq
        UPCaseSpec(
            name="UP_H20_plane_oblique3d_eq",
            topo="H20",
            lumped=1,
            E=1.8e8,
            nu=0.34,
            rho=2190.0,
            biot=BiotParams(
                Kf=2.2e9, poro=0.38, rhoF=1000.0, biotAlpha=1.0, Ks=0.0,
                pOrder=0, formulation=1, stabMode=1, stabValue=0.25,
                kbar=(2.5e-6, 1.2e-6, 6.0e-7),
            ),
            pml=pml_std,
            geom=plane_3d,
            coords=make_reference_coords("H20"),
            c1=1.0, c2=1000.0, c3=1.0e6,
            dts=[0.002, 0.002, 0.002],
        ),
    ]
    return cases


# ============================================================================ #
# Self-Checks (Single-Phase a, b, c, d + Biot u-p Cascades & 4-Block FD)       #
# ============================================================================ #

def check_partition_of_unity_and_patch() -> None:
    """(a) Verify sum_a N_a = 1, sum_a grad_x N_a = 0, and grad_x(sum_a x_a N_a) = I."""
    for topo, (ndm, _n_nodes, _, _) in TOPOLOGIES.items():
        coords = make_reference_coords(topo)
        gp_data = eval_gp_kinematics(topo, coords)
        for N, dN_dx, w_vol, _ in gp_data:
            assert w_vol > 0.0, f"{topo}: non-positive w_vol={w_vol}"
            err_pou = abs(float(np.sum(N)) - 1.0)
            assert err_pou <= 1.0e-14, f"{topo}: partition of unity err={err_pou}"
            err_zero_grad = float(np.max(np.abs(np.sum(dN_dx, axis=0))))
            assert err_zero_grad <= 1.0e-13, f"{topo}: sum grad N err={err_zero_grad}"
            grad_x = coords.T @ dN_dx
            err_patch = float(np.max(np.abs(grad_x - np.eye(ndm))))
            assert err_patch <= 1.0e-13, f"{topo}: linear patch err={err_patch}"


def check_rc_filter_z_transform() -> None:
    """(c) Verify RCFilter1D against analytical piecewise-linear z-transform + transient."""
    dt = 0.002
    omega = 2.0 * math.pi * 11.3
    z = complex(math.cos(omega * dt), math.sin(omega * dt))

    for beta in (0.0, 1.0e-11, 25.0, 140.0):
        c0 = 1.7
        c1 = -38.5
        b, phi0, phi1 = rc_coeffs(beta, dt)
        H_disc = c0 + c1 * (phi0 + phi1 * z) / (z - b)
        zeta_ss_0 = (phi0 + phi1 * z) / (z - b)

        f_re = RCFilter1D(c0, c1, beta)
        f_im = RCFilter1D(c0, c1, beta)
        # Initialize committed x^0 at n=0: x^0 = e^{i*0} = 1 + 0i, zeta^0 = 0
        f_re.x_commit = 1.0
        f_im.x_commit = 0.0
        f_re.zeta_commit = 0.0
        f_im.zeta_commit = 0.0

        for n in range(1, 40):
            x_n = z ** n
            y_re, _ = f_re.evaluate_trial(x_n.real, dt)
            y_im, _ = f_im.evaluate_trial(x_n.imag, dt)
            f_re.commit()
            f_im.commit()

            y_num = complex(y_re, y_im)
            # Analytical response with zeta^0 = 0: H_disc * z^n - c1 * zeta_ss_0 * b^n
            y_exact = H_disc * x_n - c1 * zeta_ss_0 * (b ** n)
            err = abs(y_num - y_exact) / max(1.0, abs(y_exact))
            assert err <= 1.0e-14, f"RCFilter1D z-transform mismatch at beta={beta}, n={n}: err={err}"


def check_multistage_cascade_z_transform() -> None:
    """Verify discrete z-transform transfer function for 2-stage E_j(w) and 3-stage J_s(w) cascades."""
    dt = 0.002
    omega = 2.0 * math.pi * 9.7
    z = complex(math.cos(omega * dt), math.sin(omega * dt))

    d = np.array([42.0, 55.0, 68.0], dtype=np.float64)
    kappa = np.array([1.8, 2.3, 3.1], dtype=np.float64)
    alpha = np.array([18.0, 27.0, 39.0], dtype=np.float64)

    def h_disc_stage(st: RCFilter1D) -> complex:
        if st.c1 == 0.0:
            return complex(st.c0, 0.0)
        b, phi0, phi1 = rc_coeffs(st.beta, dt)
        return st.c0 + st.c1 * (phi0 + phi1 * z) / (z - b)

    def verify_cascade_ss(cascade_re: FilterCascade, cascade_im: FilterCascade) -> None:
        # Exact steady-state initialization at n=0 for x^n = z^n:
        # Before step 1, stage m has committed input X_{m-1}^0 = prod_{r < m} H_r(z)
        # and committed state zeta_m^0 = X_{m-1}^0 * (phi0_m + phi1_m * z) / (z - b_m).
        H_prod = complex(1.0, 0.0)
        for st_re, st_im in zip(cascade_re.stages, cascade_im.stages):
            b, phi0, phi1 = rc_coeffs(st_re.beta, dt)
            x0 = H_prod
            z0 = x0 * (phi0 + phi1 * z) / (z - b)
            st_re.x_commit = x0.real
            st_im.x_commit = x0.imag
            st_re.zeta_commit = z0.real
            st_im.zeta_commit = z0.imag
            H_prod *= h_disc_stage(st_re)

        for n in range(1, 30):
            x_n = z ** n
            y_re, _ = cascade_re.evaluate_trial(x_n.real, dt)
            y_im, _ = cascade_im.evaluate_trial(x_n.imag, dt)
            cascade_re.commit()
            cascade_im.commit()
            y_num = complex(y_re, y_im)
            y_exact = H_prod * x_n
            err = abs(y_num - y_exact) / max(1.0, abs(y_exact))
            assert err <= 1.0e-13, f"Multi-stage cascade z-transform mismatch at n={n}: err={err:.3e}"

    # Check all three 2-stage E_j cascades in 3D
    for j in range(3):
        c_re = make_coupling_cascade(3, j, d, kappa, alpha)
        c_im = make_coupling_cascade(3, j, d, kappa, alpha)
        verify_cascade_ss(c_re, c_im)

    # Check 3-stage J_s cascade in 3D and 2-stage J_s cascade in 2D
    for ndm in (2, 3):
        js_re = make_storage_cascade(ndm, d, kappa, alpha)
        js_im = make_storage_cascade(ndm, d, kappa, alpha)
        verify_cascade_ss(js_re, js_im)


def check_cases_fd_and_spd() -> None:
    """(b) & (d) Verify FD Jacobian and SPD properties across all 18 cases and 3 steps."""
    rng = np.random.default_rng(20261016)
    cases = build_test_cases()
    eps = 1.0e-7

    for cs in cases:
        elem = ReferenceCPMLElement(
            topo=cs.topo,
            coords=cs.coords,
            E=cs.E,
            nu=cs.nu,
            rho=cs.rho,
            pml=cs.pml,
            geom=cs.geom,
            lumped=cs.lumped,
        )
        n_dof = elem.n_dof
        for step_idx, dt in enumerate(cs.dts):
            u0 = rng.normal(0.0, 1.0e-3, size=n_dof)
            v0 = rng.normal(0.0, 1.0e-1, size=n_dof)
            a0 = rng.normal(0.0, 1.0e1, size=n_dof)

            R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext = elem.evaluate_trial(
                u0, v0, a0, dt, cs.c1, cs.c2, cs.c3
            )
            assert abs(P_ext - float(np.dot(v0, R_int))) <= 1.0e-12 * max(1.0, abs(P_ext))

            # (d) Symmetry & SPD checks
            m_is_full_rank = (cs.lumped == 1) or (len(elem.gps) >= elem.n_nodes)
            assert np.all(np.diag(M_eff) > 0.0), f"{cs.name}: M_eff diagonal not strictly positive"
            for name, mat, mode in (
                ("K_eff", K_eff, "sym_only"),
                ("C_eff", C_eff, "psd"),
                ("M_eff", M_eff, "pd" if m_is_full_rank else "psd"),
                ("K_T_dyn", K_T_dyn, "pd" if m_is_full_rank else "sym_only"),
            ):
                norm_m = max(1.0, float(np.linalg.norm(mat, ord="fro")))
                sym_err = float(np.max(np.abs(mat - mat.T))) / norm_m
                assert sym_err <= 1.0e-14, f"{cs.name} step {step_idx} {name} sym_err={sym_err}"
                eigvals = np.linalg.eigvalsh(0.5 * (mat + mat.T))
                min_eig = float(np.min(eigvals))
                max_eig = max(1.0, float(np.max(np.abs(eigvals))))
                if mode == "pd":
                    assert min_eig > 0.0, f"{cs.name} step {step_idx} {name} not PD: min_eig={min_eig}"
                elif mode == "psd":
                    assert min_eig >= -1.0e-12 * max_eig, (
                        f"{cs.name} step {step_idx} {name} not PSD: min_eig={min_eig}"
                    )
                else:
                    # For finite dt, K_eff has O((d*dt/kappa)^2) shift on rigid rotation modes;
                    # verify it is bounded and that K_eff(dt -> 0) is PSD to machine precision.
                    assert min_eig >= -1.0e-2 * max_eig, (
                        f"{cs.name} step {step_idx} {name} unexpected negative eig: {min_eig}"
                    )
            if step_idx == 0:
                _, K_eff_0, _, _, _, _ = elem.evaluate_trial(
                    u0, v0, a0, 0.0, cs.c1, cs.c2, cs.c3
                )
                eig0 = np.linalg.eigvalsh(0.5 * (K_eff_0 + K_eff_0.T))
                assert float(np.min(eig0)) >= -1.0e-13 * max(1.0, float(np.max(np.abs(eig0)))), (
                    f"{cs.name}: K_eff(dt=0) not PSD: min_eig={np.min(eig0)}"
                )

            # (b) Central finite-difference Jacobian check under Newmark constraint (c1 = 1)
            K_fd = np.zeros((n_dof, n_dof), dtype=np.float64)
            for j in range(n_dof):
                du = np.zeros(n_dof, dtype=np.float64)
                du[j] = eps
                u_plus = u0 + du
                v_plus = v0 + cs.c2 * du
                a_plus = a0 + cs.c3 * du
                R_plus, _, _, _, _, _ = elem.evaluate_trial(
                    u_plus, v_plus, a_plus, dt, cs.c1, cs.c2, cs.c3
                )

                u_minus = u0 - du
                v_minus = v0 - cs.c2 * du
                a_minus = a0 - cs.c3 * du
                R_minus, _, _, _, _, _ = elem.evaluate_trial(
                    u_minus, v_minus, a_minus, dt, cs.c1, cs.c2, cs.c3
                )
                K_fd[:, j] = (R_plus - R_minus) / (2.0 * eps)

            # Restore trial state at (u0, v0, a0) before commit
            elem.evaluate_trial(u0, v0, a0, dt, cs.c1, cs.c2, cs.c3)

            rel_fd_err = float(np.linalg.norm(K_fd - K_T_dyn, ord="fro")) / max(
                1.0, float(np.linalg.norm(K_T_dyn, ord="fro"))
            )
            assert rel_fd_err <= 1.0e-8, (
                f"{cs.name} step {step_idx}: FD tangent mismatch rel_err={rel_fd_err:.3e}"
            )
            elem.commit()


def _sample_up_trial_state(
    rng: np.random.Generator, elem: ReferenceCPMLUPElement
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample realistic coupled (U, V, A) trial vectors with solid u and pore pressure p."""
    U = np.zeros(elem.n_dof, dtype=np.float64)
    V = np.zeros(elem.n_dof, dtype=np.float64)
    A = np.zeros(elem.n_dof, dtype=np.float64)
    for a in range(elem.n_nodes):
        ua = elem.u_off[a]
        U[ua : ua + elem.ndm] = rng.normal(0.0, 1.0e-3, size=elem.ndm)
        V[ua : ua + elem.ndm] = rng.normal(0.0, 1.0e-1, size=elem.ndm)
        A[ua : ua + elem.ndm] = rng.normal(0.0, 1.0e1, size=elem.ndm)
    for b in range(elem.n_p):
        pb = elem.p_off[b]
        U[pb] = rng.normal(0.0, 1.0e0)
        V[pb] = rng.normal(0.0, 1.0e1)
        A[pb] = rng.normal(0.0, 1.0e2)
    return U, V, A


def check_up_cases_fd_and_blocks() -> None:
    """Verify 4-block partial FD tangents, skew-duality, and Newmark FD Jacobian on all 14 UP cases."""
    rng = np.random.default_rng(20261017)
    cases = build_up_test_cases()
    eps = 1.0e-7

    for cs in cases:
        elem = ReferenceCPMLUPElement(
            topo=cs.topo,
            coords=cs.coords,
            E=cs.E,
            nu=cs.nu,
            rho=cs.rho,
            biot=cs.biot,
            pml=cs.pml,
            geom=cs.geom,
            lumped=cs.lumped,
        )
        n_dof = elem.n_dof
        u_idx = [elem.u_off[a] + d for a in range(elem.n_nodes) for d in range(elem.ndm)]
        p_idx = [elem.p_off[b] for b in range(elem.n_p)]

        for step_idx, dt in enumerate(cs.dts):
            U0, V0, A0 = _sample_up_trial_state(rng, elem)
            R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext = elem.evaluate_trial(
                U0, V0, A0, dt, cs.c1, cs.c2, cs.c3
            )
            assert abs(P_ext - float(np.dot(V0, R_int))) <= 1.0e-12 * max(1.0, abs(P_ext))

            # Verify p-rows and p-cols of M_eff are strictly zero
            assert float(np.max(np.abs(M_eff[p_idx, :]))) == 0.0
            assert float(np.max(np.abs(M_eff[:, p_idx]))) == 0.0

            # Exact skew-duality K_eff[u, p] = -C_eff[p, u]^T to <= 1e-14
            K_up = K_eff[np.ix_(u_idx, p_idx)]
            C_pu = C_eff[np.ix_(p_idx, u_idx)]
            skew_err = float(np.max(np.abs(K_up + C_pu.T))) / max(
                1.0e-30, float(np.max(np.abs(K_up)))
            )
            assert skew_err <= 1.0e-14, f"{cs.name} step {step_idx}: skew-duality err={skew_err:.3e}"

            # Partial FD derivatives with respect to U only (holding V0, A0 fixed)
            # Note: In SI units, solid stiffness K_uu ~ E ~ 1e8 causes R_u ~ 1e5,
            # whereas coupling Q_eff ~ 1 and Darcy H_eff ~ 1e-7. When perturbing p
            # by eps = 1e-7 around a state with R_u ~ 1e5 (including committed
            # solid filter history in steps >= 1), float64 precision in R_plus - R_minus
            # is bounded by machine_eps * ||R_u|| / eps ~ 1e-16 * 1e5 / 1e-7 ~ 1e-4.
            # Using eps_p = eps * cs.E for pressure perturbations (so delta_p * Q_eff
            # is commensurate with delta_u * K_uu) while keeping eps_u = eps = 1e-7
            # tests the exact (U0, V0, A0) state and committed filter history to <= 1e-7!
            eps_u = eps
            eps_p = eps * cs.E

            # 1. Perturbing u only (holding p, V, A fixed) matches (u, u) block of K_eff
            dR_du_fd = np.zeros((n_dof, len(u_idx)), dtype=np.float64)
            for col_i, j in enumerate(u_idx):
                dU = np.zeros(n_dof, dtype=np.float64)
                dU[j] = eps_u
                R_plus, _, _, _, _, _ = elem.evaluate_trial(U0 + dU, V0, A0, dt, cs.c1, cs.c2, cs.c3)
                R_minus, _, _, _, _, _ = elem.evaluate_trial(U0 - dU, V0, A0, dt, cs.c1, cs.c2, cs.c3)
                dR_du_fd[:, col_i] = (R_plus - R_minus) / (2.0 * eps_u)

            K_uu_fd = dR_du_fd[u_idx, :]
            K_uu_an = K_eff[np.ix_(u_idx, u_idx)]
            err_k_uu = float(np.linalg.norm(K_uu_fd - K_uu_an, ord="fro")) / max(
                1.0e-30, float(np.linalg.norm(K_uu_an, ord="fro"))
            )
            assert err_k_uu <= 1.0e-7, f"{cs.name} step {step_idx}: K_uu FD err={err_k_uu:.3e}"

            # 2. Perturbing p only (holding u, V, A fixed) matches (u, p) block -Q_eff and (p, p) block H_eff
            dR_dp_fd = np.zeros((n_dof, len(p_idx)), dtype=np.float64)
            for col_b, j in enumerate(p_idx):
                dU = np.zeros(n_dof, dtype=np.float64)
                dU[j] = eps_p
                R_plus, _, _, _, _, _ = elem.evaluate_trial(U0 + dU, V0, A0, dt, cs.c1, cs.c2, cs.c3)
                R_minus, _, _, _, _, _ = elem.evaluate_trial(U0 - dU, V0, A0, dt, cs.c1, cs.c2, cs.c3)
                dR_dp_fd[:, col_b] = (R_plus - R_minus) / (2.0 * eps_p)

            K_up_fd = dR_dp_fd[u_idx, :]
            err_k_up = float(np.linalg.norm(K_up_fd - K_up, ord="fro")) / max(
                1.0e-30, float(np.linalg.norm(K_up, ord="fro"))
            )
            assert err_k_up <= 1.0e-7, f"{cs.name} step {step_idx}: K_up FD err={err_k_up:.3e}"

            H_pp_fd = dR_dp_fd[p_idx, :]
            H_pp_an = K_eff[np.ix_(p_idx, p_idx)]
            err_h_pp = float(np.linalg.norm(H_pp_fd - H_pp_an, ord="fro")) / max(
                1.0e-30, float(np.linalg.norm(H_pp_an, ord="fro"))
            )
            assert err_h_pp <= 1.0e-7, f"{cs.name} step {step_idx}: H_pp FD err={err_h_pp:.3e}"

            # 3. Perturbing u_dot only matches C_eff_uu and +Q_eff^T in C_eff_pu
            dR_dv_fd = np.zeros((n_dof, len(u_idx)), dtype=np.float64)
            for col_i, j in enumerate(u_idx):
                dV = np.zeros(n_dof, dtype=np.float64)
                dV[j] = eps_u
                R_plus, _, _, _, _, _ = elem.evaluate_trial(U0, V0 + dV, A0, dt, cs.c1, cs.c2, cs.c3)
                R_minus, _, _, _, _, _ = elem.evaluate_trial(U0, V0 - dV, A0, dt, cs.c1, cs.c2, cs.c3)
                dR_dv_fd[:, col_i] = (R_plus - R_minus) / (2.0 * eps_u)

            C_uu_fd = dR_dv_fd[u_idx, :]
            C_uu_an = C_eff[np.ix_(u_idx, u_idx)]
            denom_c_uu = max(1.0, float(np.linalg.norm(C_uu_an, ord="fro")))
            err_c_uu = float(np.linalg.norm(C_uu_fd - C_uu_an, ord="fro")) / denom_c_uu
            assert err_c_uu <= 1.0e-7, f"{cs.name} step {step_idx}: C_uu FD err={err_c_uu:.3e}"

            C_pu_fd = dR_dv_fd[p_idx, :]
            err_c_pu = float(np.linalg.norm(C_pu_fd - C_pu, ord="fro")) / max(
                1.0e-30, float(np.linalg.norm(C_pu, ord="fro"))
            )
            assert err_c_pu <= 1.0e-7, f"{cs.name} step {step_idx}: C_pu FD err={err_c_pu:.3e}"

            # 4. Perturbing p_dot only matches (p, p) block S_eff + H_tilde of C_eff
            dR_dpdot_fd = np.zeros((n_dof, len(p_idx)), dtype=np.float64)
            for col_b, j in enumerate(p_idx):
                dV = np.zeros(n_dof, dtype=np.float64)
                dV[j] = eps_p
                R_plus, _, _, _, _, _ = elem.evaluate_trial(U0, V0 + dV, A0, dt, cs.c1, cs.c2, cs.c3)
                R_minus, _, _, _, _, _ = elem.evaluate_trial(U0, V0 - dV, A0, dt, cs.c1, cs.c2, cs.c3)
                dR_dpdot_fd[:, col_b] = (R_plus - R_minus) / (2.0 * eps_p)

            C_pp_fd = dR_dpdot_fd[p_idx, :]
            C_pp_an = C_eff[np.ix_(p_idx, p_idx)]
            err_c_pp = float(np.linalg.norm(C_pp_fd - C_pp_an, ord="fro")) / max(
                1.0e-30, float(np.linalg.norm(C_pp_an, ord="fro"))
            )
            assert err_c_pp <= 1.0e-7, f"{cs.name} step {step_idx}: C_pp FD err={err_c_pp:.3e}"

            # (c) Combined Newmark FD Jacobian check across all n_dof slots
            K_dyn_fd = np.zeros((n_dof, n_dof), dtype=np.float64)
            for j in range(n_dof):
                dU = np.zeros(n_dof, dtype=np.float64)
                dU[j] = eps
                R_plus, _, _, _, _, _ = elem.evaluate_trial(
                    U0 + dU, V0 + cs.c2 * dU, A0 + cs.c3 * dU, dt, cs.c1, cs.c2, cs.c3
                )
                R_minus, _, _, _, _, _ = elem.evaluate_trial(
                    U0 - dU, V0 - cs.c2 * dU, A0 - cs.c3 * dU, dt, cs.c1, cs.c2, cs.c3
                )
                K_dyn_fd[:, j] = (R_plus - R_minus) / (2.0 * eps)

            # Restore trial state at (U0, V0, A0) before commit
            elem.evaluate_trial(U0, V0, A0, dt, cs.c1, cs.c2, cs.c3)

            err_k_dyn = float(np.linalg.norm(K_dyn_fd - K_T_dyn, ord="fro")) / max(
                1.0, float(np.linalg.norm(K_T_dyn, ord="fro"))
            )
            assert err_k_dyn <= 1.0e-7, (
                f"{cs.name} step {step_idx}: Combined Newmark FD err={err_k_dyn:.3e}"
            )
            elem.commit()


def run_self_checks() -> None:
    check_partition_of_unity_and_patch()
    print("[PASS] (a) Partition of unity & linear patch reproduction across all 9 topologies.")
    check_rc_filter_z_transform()
    print("[PASS] (c) Exact discrete z-transform transfer function + transient for RCFilter1D.")
    check_multistage_cascade_z_transform()
    print("[PASS] (UP-a) Multi-stage FilterCascade discrete z-transform (2-stage E_j & 3-stage J_s <=1e-13).")
    check_cases_fd_and_spd()
    print("[PASS] (b) & (d) Finite-difference dynamic tangent (<=1e-8) and SPD checks across all 18 cases.")
    check_up_cases_fd_and_blocks()
    print("[PASS] (UP-b) & (UP-c) 4-block partial FD tangents, skew-duality (<=1e-14), and Newmark FD (<=1e-7) across all 14 UP cases.")


# ============================================================================ #
# Emit Golden Reference Files                                                  #
# ============================================================================ #

def _fmt_floats(arr: np.ndarray) -> str:
    return " ".join(f"{float(x):.17e}" for x in np.ravel(arr))


def emit_cases_file(path: str) -> None:
    rng = np.random.default_rng(20261016)
    cases = build_test_cases()

    lines: list[str] = [
        "# Autogenerated by tests/ladruno_cpml_reference.py -- DO NOT HAND-EDIT",
        f"# Total cases: {len(cases)}",
    ]

    for cs in cases:
        elem = ReferenceCPMLElement(
            topo=cs.topo,
            coords=cs.coords,
            E=cs.E,
            nu=cs.nu,
            rho=cs.rho,
            pml=cs.pml,
            geom=cs.geom,
            lumped=cs.lumped,
        )
        n_dof = elem.n_dof
        lines.append(f"CASE {cs.name}")
        lines.append(f"Topology {cs.topo} {elem.ndm} {elem.n_nodes} {cs.lumped}")
        lines.append(f"Material {cs.E:.17e} {cs.nu:.17e} {cs.rho:.17e}")
        lines.append(
            f"PMLParams {cs.pml.m_deg:.17e} {cs.pml.R0:.17e} "
            f"{cs.pml.kappa_max:.17e} {cs.pml.alpha_max:.17e}"
        )
        g_str = " ".join(f"{float(g):.17e}" for g in cs.geom.params)
        lines.append(f"GeomMode {cs.geom.mode} {len(cs.geom.params)} {g_str}")
        lines.append(f"Coords {_fmt_floats(cs.coords)}")
        lines.append(f"NewmarkCoeffs {cs.c1:.17e} {cs.c2:.17e} {cs.c3:.17e}")
        lines.append(f"NumSteps {len(cs.dts)}")

        for s, dt in enumerate(cs.dts):
            u_vec = rng.normal(0.0, 1.0e-3, size=n_dof)
            v_vec = rng.normal(0.0, 1.0e-1, size=n_dof)
            a_vec = rng.normal(0.0, 1.0e1, size=n_dof)

            R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext = elem.evaluate_trial(
                u_vec, v_vec, a_vec, dt, cs.c1, cs.c2, cs.c3
            )
            elem.commit()

            lines.append(f"STEP {s} {dt:.17e}")
            lines.append(f"U {_fmt_floats(u_vec)}")
            lines.append(f"V {_fmt_floats(v_vec)}")
            lines.append(f"A {_fmt_floats(a_vec)}")
            lines.append(f"R_int {_fmt_floats(R_int)}")
            lines.append(f"K_eff {_fmt_floats(K_eff)}")
            lines.append(f"C_eff {_fmt_floats(C_eff)}")
            lines.append(f"M_eff {_fmt_floats(M_eff)}")
            lines.append(f"K_T_dyn {_fmt_floats(K_T_dyn)}")
            lines.append(f"P_ext {P_ext:.17e}")
            lines.append("END_STEP")

        lines.append("END_CASE")

    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[EMIT] Wrote {len(cases)} reference cases to {path}")


def emit_up_cases_file(path: str) -> None:
    rng = np.random.default_rng(20261017)
    cases = build_up_test_cases()

    lines: list[str] = [
        "# Autogenerated by tests/ladruno_cpml_reference.py (--emit-up) -- DO NOT HAND-EDIT",
        f"# Total cases: {len(cases)}",
    ]

    for cs in cases:
        elem = ReferenceCPMLUPElement(
            topo=cs.topo,
            coords=cs.coords,
            E=cs.E,
            nu=cs.nu,
            rho=cs.rho,
            biot=cs.biot,
            pml=cs.pml,
            geom=cs.geom,
            lumped=cs.lumped,
        )
        lines.append(f"CASE {cs.name}")
        lines.append(f"Topology {cs.topo} {elem.ndm} {elem.n_nodes} {cs.lumped}")
        lines.append(f"Material {cs.E:.17e} {cs.nu:.17e} {cs.rho:.17e}")
        kbar_str = " ".join(f"{float(k):.17e}" for k in cs.biot.kbar)
        lines.append(
            f"BiotParams {cs.biot.Kf:.17e} {cs.biot.poro:.17e} {cs.biot.rhoF:.17e} "
            f"{cs.biot.biotAlpha:.17e} {cs.biot.Ks:.17e} {cs.biot.pOrder} "
            f"{cs.biot.formulation} {cs.biot.stabMode} {cs.biot.stabValue:.17e} {kbar_str}"
        )
        lines.append(
            f"PMLParams {cs.pml.m_deg:.17e} {cs.pml.R0:.17e} "
            f"{cs.pml.kappa_max:.17e} {cs.pml.alpha_max:.17e}"
        )
        g_str = " ".join(f"{float(g):.17e}" for g in cs.geom.params)
        lines.append(f"GeomMode {cs.geom.mode} {len(cs.geom.params)} {g_str}")
        lines.append(f"Coords {_fmt_floats(cs.coords)}")
        lines.append(f"NewmarkCoeffs {cs.c1:.17e} {cs.c2:.17e} {cs.c3:.17e}")
        lines.append(f"NumSteps {len(cs.dts)}")

        for s, dt in enumerate(cs.dts):
            u_vec, v_vec, a_vec = _sample_up_trial_state(rng, elem)
            R_int, K_eff, C_eff, M_eff, K_T_dyn, P_ext = elem.evaluate_trial(
                u_vec, v_vec, a_vec, dt, cs.c1, cs.c2, cs.c3
            )
            elem.commit()

            lines.append(f"STEP {s} {dt:.17e}")
            lines.append(f"U {_fmt_floats(u_vec)}")
            lines.append(f"V {_fmt_floats(v_vec)}")
            lines.append(f"A {_fmt_floats(a_vec)}")
            lines.append(f"R_int {_fmt_floats(R_int)}")
            lines.append(f"K_eff {_fmt_floats(K_eff)}")
            lines.append(f"C_eff {_fmt_floats(C_eff)}")
            lines.append(f"M_eff {_fmt_floats(M_eff)}")
            lines.append(f"K_T_dyn {_fmt_floats(K_T_dyn)}")
            lines.append(f"P_ext {P_ext:.17e}")
            lines.append("END_STEP")

        lines.append("END_CASE")

    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[EMIT-UP] Wrote {len(cases)} coupled Biot u-p reference cases to {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="LadrunoCPML P0 & P0-UP Reference Oracle")
    parser.add_argument("--self-check", action="store_true", help="Run built-in verification checks")
    parser.add_argument("--emit", type=str, default="", help="Output path for single-phase golden reference cases")
    parser.add_argument("--emit-up", type=str, default="", help="Output path for coupled Biot u-p golden reference cases")
    args = parser.parse_args()

    if args.self_check or (not args.emit and not args.emit_up):
        run_self_checks()
    if args.emit:
        emit_cases_file(args.emit)
    if args.emit_up:
        emit_up_cases_file(args.emit_up)
    return 0


if __name__ == "__main__":
    sys.exit(main())

