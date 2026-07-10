"""Independent numpy oracle for the LadrunoUP Biot u-p kernel (ADR-71 P0, WP 0.B).

This file is the *independent* half of the ADR-71 P0 cross-verification: it
re-derives every kernel block (Q, H, S, H-tilde, f_seep) and the DOF map from
first principles (own shape functions, own Jacobian, own quadrature, own
assembly loops) so that agreeing to <=1e-9 with the C++ `LadrunoUPKernel.h`
*is* the verification. It reads NOTHING under SRC/ and nothing of the sibling
C++ agent's files (independence protocol, plan section "Execution model").

Block definitions implemented verbatim from the FROZEN P0 interface contract
(_adr71_implementation_plan.md, "P0 interface contract" + "0.A addendum"):

    Q[(a*ndm+i), b] += dNu[a,i] * biotAlpha * Np[b]          * dv     (coupling)
    H[b,c]          += sum_i dNp[b,i] * kbar[i] * dNp[c,i]   * dv     (Darcy)
    S[b,c]          += Np[b] * oneOverQbar * Np[c]           * dv     (storage)
    HT[b,c]         += sum_i dNp[b,i] * stabAlpha * dNp[c,i] * dv     (McGann stab)
    fseep[b]        += sum_i dNp[b,i] * kbar[i] * rhoF * drive[i] * dv (seepage)

    stabAlphaAuto = alpha0 * h*h / (Ks_skel + 4*Gs_skel/3)   (McGann 2012 eq 73)
    oneOverQbar   = n/Kf + (alpha-n)/Ks_grain                (Book eq 2.17; Ks<=0 => inf)

Emitted blocks use the CONTRACTUAL GP rule per shape (0.A addendum item 1):
    T3 = 1pt (1/3,1/3) w=1/2 ; Q4 = 2x2 Gauss ; H8 = 2x2x2 Gauss ;
    BT6 = 3pt interior {(1/6,1/6),(2/3,1/6),(1/6,2/3)} w=1/6 ;
    BTET10 = 4pt (a=0.585410196624968, b=0.138196601125011) w=1/24.
Independence lives in the *implementation*, not the GP points (which are
contract data). The self-check additionally cross-checks contractual-rule
blocks against degree-overkill (collapsed-Gauss / tensor-Gauss) integration on
the provably-exact case/block combos and asserts <=1e-12 there; known-inexact
combos are documented, not asserted.

Run:
    python tests/ladruno_up_reference.py                     # self-check (hard asserts)
    python tests/ladruno_up_reference.py --emit tests/ladruno_up_cases.txt

Pure numpy + scipy only (no opensees import).
"""
from __future__ import annotations

import argparse
import itertools
import sys

import numpy as np

try:
    from scipy.integrate import solve_ivp
    from scipy.linalg import expm
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover - environment guard
    _HAVE_SCIPY = False


# ============================================================================ #
#  Shape providers.  Each *_geom returns (N, dNdparam) for the FULL geometry    #
#  basis at a parametric point; ndm_param == ndm.  Pressure bases are derived   #
#  separately (equal-order == geometry basis; Taylor-Hood == barycentric-linear #
#  on the vertex carriers), per 0.A addendum items 3 and 4.                     #
# ============================================================================ #

_S3 = 1.0 / np.sqrt(3.0)


def q4_geom(p):
    xi, eta = p
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([
        [-(1 - eta), -(1 - xi)],
        [(1 - eta), -(1 + xi)],
        [(1 + eta), (1 + xi)],
        [-(1 + eta), (1 - xi)],
    ])
    return N, dN


def t3_geom(p):
    # param (xi1, xi2), xi3 = 1 - xi1 - xi2 ; nodes CCW: (xi3), (xi1), (xi2)
    xi1, xi2 = p
    N = np.array([1.0 - xi1 - xi2, xi1, xi2])
    dN = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    return N, dN


def h8_geom(p):
    xi, eta, ze = p
    sgn = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    N = np.empty(8)
    dN = np.empty((8, 3))
    for a in range(8):
        xa, ya, za = sgn[a]
        N[a] = 0.125 * (1 + xa * xi) * (1 + ya * eta) * (1 + za * ze)
        dN[a, 0] = 0.125 * xa * (1 + ya * eta) * (1 + za * ze)
        dN[a, 1] = 0.125 * (1 + xa * xi) * ya * (1 + za * ze)
        dN[a, 2] = 0.125 * (1 + xa * xi) * (1 + ya * eta) * za
    return N, dN


def bt6_geom(p):
    # Bernstein quadratic (0.A addendum item 3).  param (xi1, xi2), xi3=1-xi1-xi2.
    # N0=xi3^2, N1=xi1^2, N2=xi2^2, N3=2 xi1 xi3, N4=2 xi1 xi2, N5=2 xi2 xi3.
    xi1, xi2 = p
    xi3 = 1.0 - xi1 - xi2
    N = np.array([xi3 * xi3, xi1 * xi1, xi2 * xi2,
                  2 * xi1 * xi3, 2 * xi1 * xi2, 2 * xi2 * xi3])
    dN = np.array([
        [-2 * xi3, -2 * xi3],
        [2 * xi1, 0.0],
        [0.0, 2 * xi2],
        [2 * (xi3 - xi1), -2 * xi1],
        [2 * xi2, 2 * xi1],
        [-2 * xi2, 2 * (xi3 - xi2)],
    ])
    return N, dN


def btet10_geom(p):
    # Bernstein quadratic tet (0.A addendum item 3).  param (L1,L2,L3), L4=1-..
    # N0..3=L1^2..L4^2 ; N4=2L1L2,N5=2L2L3,N6=2L1L3,N7=2L1L4,N8=2L3L4,N9=2L2L4.
    L1, L2, L3 = p
    L4 = 1.0 - L1 - L2 - L3
    N = np.array([L1 * L1, L2 * L2, L3 * L3, L4 * L4,
                  2 * L1 * L2, 2 * L2 * L3, 2 * L1 * L3,
                  2 * L1 * L4, 2 * L3 * L4, 2 * L2 * L4])
    dN = np.array([
        [2 * L1, 0.0, 0.0],
        [0.0, 2 * L2, 0.0],
        [0.0, 0.0, 2 * L3],
        [-2 * L4, -2 * L4, -2 * L4],
        [2 * L2, 2 * L1, 0.0],
        [0.0, 2 * L3, 2 * L2],
        [2 * L3, 0.0, 2 * L1],
        [2 * (L4 - L1), -2 * L1, -2 * L1],
        [-2 * L3, -2 * L3, 2 * (L4 - L3)],
        [-2 * L2, 2 * (L4 - L2), -2 * L2],
    ])
    return N, dN


GEOM = {"Q4": q4_geom, "T3": t3_geom, "H8": h8_geom,
        "BT6": bt6_geom, "BTET10": btet10_geom}
NNODES = {"Q4": 4, "T3": 3, "H8": 8, "BT6": 6, "BTET10": 10}
NDM = {"Q4": 2, "T3": 2, "H8": 3, "BT6": 2, "BTET10": 3}
# vertex (p-carrier) node indices for Taylor-Hood pressure
TH_VERTICES = {"BT6": [0, 1, 2], "BTET10": [0, 1, 2, 3]}


def pressure_basis(shape, porder, p):
    """Return (Np, dNp_param) for the pressure field at parametric point p.

    equal-order  : Np == geometry N (all nodes carriers).
    th (T-H)     : barycentric-linear on the vertex carriers, COMPACTED in node
                   order restricted to carriers (0.A addendum item 4).
    """
    if porder == "equal":
        N, dN = GEOM[shape](p)
        return N, dN
    # Taylor-Hood
    if shape == "BT6":
        xi1, xi2 = p
        xi3 = 1.0 - xi1 - xi2
        Np = np.array([xi3, xi1, xi2])            # node0<->xi3, node1<->xi1, node2<->xi2
        dNp = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
        return Np, dNp
    if shape == "BTET10":
        L1, L2, L3 = p
        L4 = 1.0 - L1 - L2 - L3
        Np = np.array([L1, L2, L3, L4])           # node0<->L1 ... node3<->L4
        dNp = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0], [-1.0, -1.0, -1.0]])
        return Np, dNp
    raise ValueError(f"Taylor-Hood not defined for shape {shape}")


def n_pressure(shape, porder):
    if porder == "equal":
        return NNODES[shape]
    return len(TH_VERTICES[shape])


# ============================================================================ #
#  Quadrature rules.  Parametric-domain points + weights.  Contractual weights  #
#  sum to the parametric-domain measure (tri = 1/2, tet = 1/6, quad/hex = 2^ndm #
#  over the [-1,1] cube).  Overkill rules use the same measure so a block from   #
#  either rule is directly comparable.                                          #
# ============================================================================ #

def gp_contract(shape):
    if shape == "T3":
        return [np.array([1 / 3, 1 / 3])], [0.5]
    if shape == "Q4":
        pts, ws = [], []
        for a in (-_S3, _S3):
            for b in (-_S3, _S3):
                pts.append(np.array([a, b]))
                ws.append(1.0)
        return pts, ws
    if shape == "H8":
        pts, ws = [], []
        for a in (-_S3, _S3):
            for b in (-_S3, _S3):
                for c in (-_S3, _S3):
                    pts.append(np.array([a, b, c]))
                    ws.append(1.0)
        return pts, ws
    if shape == "BT6":
        pts = [np.array([1 / 6, 1 / 6]), np.array([2 / 3, 1 / 6]),
               np.array([1 / 6, 2 / 3])]
        return pts, [1 / 6, 1 / 6, 1 / 6]
    if shape == "BTET10":
        a = 0.585410196624968
        b = 0.138196601125011
        pts = [np.array([a, b, b]), np.array([b, a, b]),
               np.array([b, b, a]), np.array([b, b, b])]
        return pts, [1 / 24, 1 / 24, 1 / 24, 1 / 24]
    raise ValueError(shape)


def _gauss1d(n):
    x, w = np.polynomial.legendre.leggauss(n)  # on [-1,1]
    return x, w


def gp_overkill(shape):
    """High-order rule integrating the (polynomial) reference integrands exactly."""
    if shape in ("Q4",):
        x, w = _gauss1d(8)
        pts, ws = [], []
        for i in range(8):
            for j in range(8):
                pts.append(np.array([x[i], x[j]]))
                ws.append(w[i] * w[j])
        return pts, ws
    if shape == "H8":
        x, w = _gauss1d(6)
        pts, ws = [], []
        for i in range(6):
            for j in range(6):
                for k in range(6):
                    pts.append(np.array([x[i], x[j], x[k]]))
                    ws.append(w[i] * w[j] * w[k])
        return pts, ws
    if shape in ("T3", "BT6"):
        # Duffy collapse of [0,1]^2 -> reference triangle; jac (1-u); measure 1/2.
        x, w = _gauss1d(12)
        u = 0.5 * (x + 1.0)   # map [-1,1] -> [0,1]
        wu = 0.5 * w
        pts, ws = [], []
        for i in range(len(u)):
            for j in range(len(u)):
                xi1 = u[i]
                xi2 = u[j] * (1.0 - u[i])
                pts.append(np.array([xi1, xi2]))
                ws.append(wu[i] * wu[j] * (1.0 - u[i]))
        return pts, ws
    if shape in ("BTET10",):
        # Duffy collapse of [0,1]^3 -> reference tet; jac (1-u)^2 (1-v); measure 1/6.
        x, w = _gauss1d(10)
        u = 0.5 * (x + 1.0)
        wu = 0.5 * w
        pts, ws = [], []
        for i in range(len(u)):
            for j in range(len(u)):
                for k in range(len(u)):
                    L1 = u[i]
                    L2 = u[j] * (1.0 - u[i])
                    L3 = u[k] * (1.0 - u[i]) * (1.0 - u[j])
                    jac = (1.0 - u[i]) ** 2 * (1.0 - u[j])
                    pts.append(np.array([L1, L2, L3]))
                    ws.append(wu[i] * wu[j] * wu[k] * jac)
        return pts, ws
    raise ValueError(shape)


# ============================================================================ #
#  Jacobian + cartesian gradients.  Convention (documented, independently       #
#  derived): Jac[i,j] = d x_i / d xi_j = (coords^T @ dNdparam) ; detJ=det(Jac) ; #
#  dNdx[a,i] = sum_j dNdparam[a,j] * inv(Jac)[j,i]  ==  dNdparam @ inv(Jac).     #
# ============================================================================ #

def jacobian(coords, dNdparam):
    Jac = coords.T @ dNdparam          # (ndm, ndm)
    detJ = np.linalg.det(Jac)
    Jinv = np.linalg.inv(Jac)
    return Jac, detJ, Jinv


def eval_gp(shape, porder, coords, p):
    """Return N, dNdx (geom), detJ, Np, dNpdx at parametric point p."""
    N, dNp_geom = GEOM[shape](p)
    Jac, detJ, Jinv = jacobian(coords, dNp_geom)
    dNdx = dNp_geom @ Jinv
    Np, dNp_par = pressure_basis(shape, porder, p)
    dNpdx = dNp_par @ Jinv
    return N, dNdx, detJ, Np, dNpdx


# ============================================================================ #
#  Block assembly (the contract, verbatim).                                     #
# ============================================================================ #

def assemble_blocks(shape, porder, coords, params, pts, ws):
    ndm = NDM[shape]
    nNu = NNODES[shape]
    nNp = n_pressure(shape, porder)
    alpha = params["biotAlpha"]
    ooq = params["oneOverQbar"]
    rhoF = params["rhoF"]
    stab = params["stabAlpha"]
    thick = params["thick"]
    kbar = np.asarray(params["kbar"], float)
    drive = np.asarray(params["drive"], float)

    Q = np.zeros((nNu * ndm, nNp))
    H = np.zeros((nNp, nNp))
    S = np.zeros((nNp, nNp))
    HT = np.zeros((nNp, nNp))
    fseep = np.zeros(nNp)
    detJ_list = []

    for p, w in zip(pts, ws):
        N, dNdx, detJ, Np, dNpdx = eval_gp(shape, porder, coords, p)
        detJ_list.append(detJ)
        dv = w * detJ * (thick if ndm == 2 else 1.0)
        # Q
        for a in range(nNu):
            for i in range(ndm):
                r = a * ndm + i
                Q[r, :] += dNdx[a, i] * alpha * Np * dv
        # H, HT (contract: sum_i dNp[b,i]*k_i*dNp[c,i])
        for b in range(nNp):
            for c in range(nNp):
                H[b, c] += np.sum(dNpdx[b, :] * kbar * dNpdx[c, :]) * dv
                HT[b, c] += np.sum(dNpdx[b, :] * stab * dNpdx[c, :]) * dv
        # S
        S += np.outer(Np, Np) * (ooq * dv)
        # fseep
        for b in range(nNp):
            fseep[b] += np.sum(dNpdx[b, :] * kbar * rhoF * drive) * dv
    return dict(Q=Q, H=H, S=S, HT=HT, fseep=fseep, detJ=np.array(detJ_list))


def dof_map(shape, porder):
    """Per-node-interleaved [u..., p?] map (contract dofMap)."""
    ndm = NDM[shape]
    nNodes = NNODES[shape]
    if porder == "equal":
        carriers = set(range(nNodes))
    else:
        carriers = set(TH_VERTICES[shape])
    uOff = [0] * nNodes
    pOff = [-1] * nNodes
    slot = 0
    for n in range(nNodes):
        uOff[n] = slot
        slot += ndm
        if n in carriers:
            pOff[n] = slot
            slot += 1
    return slot, uOff, pOff


# ============================================================================ #
#  Scalar helpers (contract).                                                   #
# ============================================================================ #

def stab_alpha_auto(h, Ks_skel, Gs_skel, alpha0):
    return alpha0 * h * h / (Ks_skel + 4.0 * Gs_skel / 3.0)


def one_over_qbar(n, Kf, alpha, Ks_grain):
    val = n / Kf
    if Ks_grain > 0.0:
        val += (alpha - n) / Ks_grain
    return val


# Element edges (vertex-index pairs) for the h definition.  On simplices every
# vertex pair IS an edge; Q4/H8 exclude the diagonals.
_EDGES = {
    "T3": [(0, 1), (1, 2), (2, 0)],
    "Q4": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "H8": [(0, 1), (1, 2), (2, 3), (3, 0),
           (4, 5), (5, 6), (6, 7), (7, 4),
           (0, 4), (1, 5), (2, 6), (3, 7)],
    "BT6": [(0, 1), (1, 2), (2, 0)],
    "BTET10": [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
}


def elem_size_h(shape, coords):
    """Largest EDGE length (ADR section 3.3 "largest element dimension").

    Contract amendment (0.E panel, FIX 1 — adjudicated by MAIN): the plan's
    "largest vertex-pair distance" mis-cited the ADR.  The B4 calibration
    (10x10 mesh of 30 m => h = 3 m = cell SIDE) reproduces only with the edge,
    not the diagonal.  Simplices unchanged (every vertex pair is an edge)."""
    hmax = 0.0
    for i, j in _EDGES[shape]:
        hmax = max(hmax, float(np.linalg.norm(coords[i] - coords[j])))
    return hmax


# ============================================================================ #
#  Independent plain-elastic B-matrix machinery (for consolidation + inf-sup).  #
#  Built from scratch here so the ODE / Schur gates do not reuse kernel code.    #
# ============================================================================ #

def elastic_D(ndm, E, nu):
    """Isotropic elasticity, Voigt: plane strain (2D) or full 3D."""
    if ndm == 2:
        c = E / ((1 + nu) * (1 - 2 * nu))
        return c * np.array([[1 - nu, nu, 0.0],
                             [nu, 1 - nu, 0.0],
                             [0.0, 0.0, (1 - 2 * nu) / 2.0]])
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] += 2 * mu
    for i in range(3, 6):
        D[i, i] = mu
    return D


def _bmat(dNdx, ndm, nNu):
    """Voigt strain-displacement matrix (engineering shear)."""
    if ndm == 2:
        B = np.zeros((3, 2 * nNu))
        for a in range(nNu):
            B[0, 2 * a] = dNdx[a, 0]
            B[1, 2 * a + 1] = dNdx[a, 1]
            B[2, 2 * a] = dNdx[a, 1]
            B[2, 2 * a + 1] = dNdx[a, 0]
        return B
    B = np.zeros((6, 3 * nNu))
    for a in range(nNu):
        bx, by, bz = dNdx[a]
        B[0, 3 * a] = bx
        B[1, 3 * a + 1] = by
        B[2, 3 * a + 2] = bz
        B[3, 3 * a] = by          # gamma_xy
        B[3, 3 * a + 1] = bx
        B[4, 3 * a + 1] = bz      # gamma_yz
        B[4, 3 * a + 2] = by
        B[5, 3 * a] = bz          # gamma_zx
        B[5, 3 * a + 2] = bx
    return B


def elem_K_Q(shape, porder, coords, E, nu, alpha, thick, pts, ws):
    """Independent K (B^T D B) and Q (B^T alpha m Np) for any shape (2D/3D)."""
    ndm = NDM[shape]
    nNu = NNODES[shape]
    nNp = n_pressure(shape, porder)
    D = elastic_D(ndm, E, nu)
    m = np.array([1.0, 1.0, 0.0]) if ndm == 2 else np.array([1, 1, 1, 0, 0, 0], float)
    K = np.zeros((ndm * nNu, ndm * nNu))
    Q = np.zeros((ndm * nNu, nNp))
    for p, w in zip(pts, ws):
        N, dNdx, detJ, Np, dNpdx = eval_gp(shape, porder, coords, p)
        assert detJ > 0.0, f"elem_K_Q: non-positive detJ {detJ:.3e} ({shape})"
        dv = w * detJ * (thick if ndm == 2 else 1.0)
        B = _bmat(dNdx, ndm, nNu)
        K += B.T @ D @ B * dv
        Q += np.outer(B.T @ (alpha * m), Np) * dv
    return K, Q


# ============================================================================ #
#  SELF-CHECK GATES                                                             #
# ============================================================================ #

def _relerr(A, B):
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    scale = max(np.max(np.abs(B)), 1e-30)
    return np.max(np.abs(A - B)) / scale


def gate_overkill(report):
    """Gate (a): contractual-rule blocks vs degree-overkill on exact combos.

    Provably-exact combos (documented derivation in comments; HT is exact
    wherever H is — identical integrand family with scalar stab vs diag k;
    FSEEP on equal-order Bezier shapes is a degree-1 integrand vs the degree-2
    contractual rules, hence exact — 0.E FIX 5):
      - T3        : Q, H, HT, FSEEP exact (affine, const dNu) ; S quadratic -> inexact.
      - Q4 affine : all exact (const Jac) ; distorted (trapezoid) inexact.
      - H8 cube   : all exact ; distorted inexact.
      - BT6 TH    : all exact (affine, Np linear, integrand <= deg2).
      - BT6 equal : H, HT, FSEEP exact ; Q (cubic), S (quartic) -> inexact.
      - BTET10 TH : all exact.
      - BTET10 eq : H, HT, FSEEP exact ; Q, S -> inexact.
    """
    print("\n[gate a] overkill-vs-contractual quadrature on provably-exact combos")
    worst = 0.0
    checks = []
    # combos: (name, shape, porder, coords, params, exact_blocks)
    baseP = dict(biotAlpha=1.0, oneOverQbar=1.3e-3, rhoF=1.0, stabAlpha=7.0e-6,
                 thick=1.7, kbar=None, drive=None)

    def mk(shape, coords, aniso=False):
        pr = dict(baseP)
        nd = NDM[shape]
        pr["kbar"] = ([1e-5, 5e-6, 2e-5][:nd] if aniso else [1e-5] * nd)
        pr["drive"] = ([0.0, -9.81, 3.0][:nd])
        if nd == 3:
            pr["thick"] = 1.0
        return pr

    # affine geometries
    tri = np.array([[0.2, 0.1], [1.3, 0.0], [0.4, 1.1]])
    quad_par = np.array([[0.0, 0.0], [1.4, 0.2], [1.7, 1.1], [0.3, 0.9]])  # parallelogram
    quad_par[2] = quad_par[1] + quad_par[3] - quad_par[0]
    cube = _unit_cube()
    v = np.array([[0.0, 0.0], [1.4, 0.2], [0.3, 1.1]])
    bt6 = _bt6_coords(v)
    # positively oriented under the pinned node<->L map (0.E FIX 2)
    tetv = np.array([[1.2, 0.1, 0.0], [0.1, 0.0, 0.0],
                     [0.2, 1.1, 0.1], [0.0, 0.2, 1.3]])
    btet = _btet10_coords(tetv)

    combos = [
        ("T3 (Q,H,HT,FSEEP exact)", "T3", "equal", tri, mk("T3", tri),
         ["Q", "H", "HT", "fseep"], ["S"]),
        ("Q4 parallelogram (all exact)", "Q4", "equal", quad_par, mk("Q4", quad_par),
         ["Q", "H", "HT", "S", "fseep"], []),
        ("H8 cube (all exact)", "H8", "equal", cube, mk("H8", cube),
         ["Q", "H", "HT", "S", "fseep"], []),
        ("BT6 TH (all exact)", "BT6", "th", bt6, mk("BT6", bt6),
         ["Q", "H", "HT", "S", "fseep"], []),
        ("BT6 equal (H,HT,FSEEP exact)", "BT6", "equal", bt6, mk("BT6", bt6),
         ["H", "HT", "fseep"], ["Q", "S"]),
        ("BTET10 TH (all exact)", "BTET10", "th", btet, mk("BTET10", btet),
         ["Q", "H", "HT", "S", "fseep"], []),
        ("BTET10 equal (H,HT,FSEEP exact)", "BTET10", "equal", btet, mk("BTET10", btet),
         ["H", "HT", "fseep"], ["Q", "S"]),
    ]
    for name, shape, porder, coords, pr, exact, inexact in combos:
        pc, wc = gp_contract(shape)
        po, wo = gp_overkill(shape)
        bc = assemble_blocks(shape, porder, coords, pr, pc, wc)
        bo = assemble_blocks(shape, porder, coords, pr, po, wo)
        line = f"  {name:32s}"
        for blk in exact:
            e = _relerr(bc[blk], bo[blk])
            worst = max(worst, e)
            line += f"  {blk}={e:.1e}"
            assert e <= 1e-12, f"EXACT combo {name}/{blk} relerr {e:.3e} > 1e-12"
        for blk in inexact:
            e = _relerr(bc[blk], bo[blk])
            line += f"  [{blk} inexact:{e:.1e}]"
        print(line)
        checks.append((name, exact, inexact))
    print(f"  -> worst exact-combo relerr = {worst:.3e} (<=1e-12 required)  PASS")
    return worst


def gate_symmetry():
    """Gate (b): symmetry of H,S,HT ; partition-of-unity pattern of Q."""
    print("\n[gate b] symmetry (H,S,HT) + Q coupling-layout / partition-of-unity")
    worst_sym = 0.0
    worst_pou = 0.0
    for shape in ("Q4", "T3", "H8", "BT6", "BTET10"):
        for porder in (("equal", "th") if shape in TH_VERTICES else ("equal",)):
            coords = _canonical_coords(shape)
            nd = NDM[shape]
            pr = dict(biotAlpha=0.83, oneOverQbar=2.0e-3, rhoF=1.0,
                      stabAlpha=5e-6, thick=(1.0 if nd == 3 else 1.9),
                      kbar=[1e-5, 6e-6, 2e-5][:nd], drive=[0.0, -9.81, 2.0][:nd])
            pc, wc = gp_contract(shape)
            b = assemble_blocks(shape, porder, coords, pr, pc, wc)
            for blk in ("H", "S", "HT"):
                s = np.max(np.abs(b[blk] - b[blk].T))
                worst_sym = max(worst_sym, s)
                assert s <= 1e-14, f"{shape}/{porder} {blk} not symmetric ({s:.2e})"
            # Q partition of unity: sum over u-nodes of dNu[a,i] == 0 -> per (i,b)
            # column-block sums vanish.  Q[a*ndm+i,b]; sum over a of each == 0.
            nNu = NNODES[shape]
            nNp = n_pressure(shape, porder)
            Q = b["Q"]
            pou = 0.0
            for i in range(nd):
                for bcol in range(nNp):
                    ssum = sum(Q[a * nd + i, bcol] for a in range(nNu))
                    pou = max(pou, abs(ssum))
            worst_pou = max(worst_pou, pou)
            assert pou <= 1e-12, f"{shape}/{porder} Q partition-of-unity {pou:.2e}"
            # sign/layout: Q must be non-trivial and rectangular (nNu*ndm x nNp)
            assert Q.shape == (nNu * nd, nNp)
            assert np.max(np.abs(Q)) > 0.0
    print(f"  max |H-H^T|,|S-S^T|,|HT-HT^T| = {worst_sym:.2e} (<=1e-14)  PASS")
    print(f"  max |sum_a Q[a,i;b]| (partition of unity) = {worst_pou:.2e}  PASS")


def gate_consolidation():
    """Gate (c): 1-element consolidation ODE from oracle blocks (semi-discrete).

    Oedometer Q4 (unit square, plane strain): confined laterally (u_x=0 all
    nodes), base fixed (u_y=0 bottom nodes), drained top (p=0 top nodes),
    sealed elsewhere.  Quasi-static (drop M,C), with the stabilized storage
    bracket per the plan's "[S+Htilde, H, Q^T]":
        K u - Q p = f_u ;   Q^T u_dot + (S+HT) p_dot + H p = 0
    Eliminate u (u_dot = K^-1 Q p_dot):
        ((S+HT) + Q^T K^-1 Q) p_dot = -H p    ->    A p_dot = -H p
    Closed form p(t) = expm(-A^-1 H t) p0 vs scipy.integrate.

    HONESTY (0.E FIX 6): this is a STRUCTURAL sanity check — it verifies the
    assembled bracket is SPD, the coupled system is dissipative (all decay
    rates real-positive), the pressure norm decays monotonically, and the
    closed form matches the integrated ODE.  It CANNOT catch a Q sign flip
    (Q enters as Q^T K^-1 Q, sign-invariant) or an S scale factor (rescales
    the rates but stays dissipative).  Block *correctness* is carried by the
    C++ cross-check (0.D), the overkill gate (a), the hydrostatic sign gate
    (g), and the PSD/orientation asserts in the case loop.
    """
    print("\n[gate c] 1-element consolidation ODE (blocks -> decay vs scipy/expm)")
    shape, porder = "Q4", "equal"
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    E, nu, alpha, thick = 1.0e4, 0.3, 1.0, 1.0
    n, Kf = 0.4, 2.2e5
    # McGann auto stabilization from THIS gate's skeleton moduli (FIX 6: HT
    # genuinely enters the assembled bracket)
    Ks = E / (3 * (1 - 2 * nu))
    Gs = E / (2 * (1 + nu))
    stab = stab_alpha_auto(elem_size_h(shape, coords), Ks, Gs, 0.25)
    pr = dict(biotAlpha=alpha, oneOverQbar=one_over_qbar(n, Kf, alpha, -1.0),
              rhoF=1.0, stabAlpha=stab, thick=thick, kbar=[1e-4, 1e-4],
              drive=[0.0, 0.0])
    pc, wc = gp_contract(shape)
    blk = assemble_blocks(shape, porder, coords, pr, pc, wc)
    S, H, HT = blk["S"], blk["H"], blk["HT"]
    K, Q = elem_K_Q(shape, porder, coords, E, nu, alpha, thick, pc, wc)

    # DOF bookkeeping: nodes 0,1,2,3 CCW.  u dofs (2 per node), p dofs (1 per node).
    # BC: u_x=0 all ; u_y=0 nodes 0,1 (bottom) ; top nodes 2,3 u_y free.
    free_u = [2 * 2 + 1, 2 * 3 + 1]          # u_y of nodes 2,3
    # p: drained (p=0) at top nodes 2,3 ; free at bottom nodes 0,1
    free_p = [0, 1]
    Kff = K[np.ix_(free_u, free_u)]
    Qfp = Q[np.ix_(free_u, free_p)]
    Sff = (S + HT)[np.ix_(free_p, free_p)]
    Hff = H[np.ix_(free_p, free_p)]

    A = Sff + Qfp.T @ np.linalg.solve(Kff, Qfp)
    # SPD / dissipativity checks
    assert np.allclose(A, A.T, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(A)) > 0, "A not SPD"
    assert np.min(np.linalg.eigvalsh(Hff)) > 0, "H_ff not SPD (needs drained node)"
    M = np.linalg.solve(A, Hff)             # p_dot = -M p
    eigM = np.linalg.eigvals(M)
    assert np.all(eigM.real > 0), f"non-dissipative eig {eigM}"
    assert np.max(np.abs(eigM.imag)) < 1e-9, "complex decay rates"

    p0 = np.array([1.0, 1.0])
    tgrid = np.linspace(0.0, 5.0 / np.min(eigM.real), 40)
    # closed form
    p_cf = np.array([expm(-M * t) @ p0 for t in tgrid])
    if _HAVE_SCIPY:
        sol = solve_ivp(lambda t, y: -(M @ y), (tgrid[0], tgrid[-1]), p0,
                        t_eval=tgrid, rtol=1e-10, atol=1e-12)
        p_num = sol.y.T
        err = np.max(np.abs(p_num - p_cf))
    else:
        p_num, err = p_cf, 0.0
    norms = np.linalg.norm(p_cf, axis=1)
    monotone = np.all(np.diff(norms) <= 1e-12)
    print(f"  decay rates (1/time): {np.sort(eigM.real)}")
    print(f"  scipy-vs-expm max diff = {err:.3e} (<=1e-8) ; monotone decay = {monotone}")
    assert err <= 1e-8, f"consolidation ODE mismatch {err:.3e}"
    assert monotone, "pressure norm not monotonically decreasing"
    print("  -> consolidation gate PASS")


def _kuhn_tets(origin, dx):
    """Six-tet Kuhn subdivision of the cube at `origin` with side dx.  Each tet
    is re-wound (v0<->v1 swap) if needed so detJ > 0 under the pinned
    node<->L map (0.E FIX 2 orientation discipline)."""
    tets = []
    for perm in itertools.permutations((0, 1, 2)):
        v = [np.array(origin, float)]
        q = np.array(origin, float)
        for ax in perm:
            q = q.copy()
            q[ax] += dx
            v.append(q)
        verts = np.array(v)
        d = np.linalg.det(np.column_stack([verts[0] - verts[3],
                                           verts[1] - verts[3],
                                           verts[2] - verts[3]]))
        if d < 0:
            verts[[0, 1]] = verts[[1, 0]]
        tets.append(verts)
    return tets


def build_patch(shape, ncells, L=1.0):
    """Structured patch of `shape` elements on [0,L]^ndm with ncells per axis.

    Returns (nodes, elems, boundary_node_set).  Triangles = single-diagonal
    (ll->ur) split of each quad cell; tets = Kuhn 6-tet split of each cube;
    Bezier shapes get straight-side edge midpoints (shared via a registry).
    """
    reg = {}
    pts_list = []

    def get(pt):
        key = tuple(round(float(c), 12) for c in pt)
        if key not in reg:
            reg[key] = len(pts_list)
            pts_list.append([float(c) for c in pt])
        return reg[key]

    dx = L / ncells
    elems = []
    if shape in ("Q4", "T3", "BT6"):
        for j in range(ncells):
            for i in range(ncells):
                ll = np.array([i * dx, j * dx])
                lr = np.array([(i + 1) * dx, j * dx])
                ur = np.array([(i + 1) * dx, (j + 1) * dx])
                ul = np.array([i * dx, (j + 1) * dx])
                if shape == "Q4":
                    elems.append([get(ll), get(lr), get(ur), get(ul)])
                    continue
                for tri in ((ll, lr, ur), (ll, ur, ul)):   # CCW both
                    if shape == "T3":
                        elems.append([get(t) for t in tri])
                    else:
                        v0, v1, v2 = tri
                        elems.append([get(v0), get(v1), get(v2),
                                      get(0.5 * (v0 + v1)),
                                      get(0.5 * (v1 + v2)),
                                      get(0.5 * (v2 + v0))])
    elif shape == "H8":
        for k in range(ncells):
            for j in range(ncells):
                for i in range(ncells):
                    corners = [(i, j, k), (i + 1, j, k), (i + 1, j + 1, k),
                               (i, j + 1, k), (i, j, k + 1), (i + 1, j, k + 1),
                               (i + 1, j + 1, k + 1), (i, j + 1, k + 1)]
                    elems.append([get((a * dx, b * dx, c * dx))
                                  for a, b, c in corners])
    elif shape == "BTET10":
        for k in range(ncells):
            for j in range(ncells):
                for i in range(ncells):
                    for verts in _kuhn_tets((i * dx, j * dx, k * dx), dx):
                        elems.append([get(pt) for pt in _btet10_coords(verts)])
    else:
        raise ValueError(shape)
    nodes = np.array(pts_list)
    tol = 1e-9
    boundary = {a for a in range(len(nodes))
                if any(abs(c) < tol or abs(c - L) < tol for c in nodes[a])}
    return nodes, elems, boundary


def _infsup_pair(shape, porder, ncells, stabilize):
    """Assemble the pressure Schur complement S_p = Q^T K^-1 Q (+HT if
    stabilize) on a structured patch, u fixed on the whole boundary, p free.
    Returns (#eigs < 1e-10*lambda_max, total p dofs)."""
    nodes, elems, boundary = build_patch(shape, ncells)
    ndm = NDM[shape]
    nn = len(nodes)
    E, nu, alpha = 1.0e4, 0.3, 1.0
    Ks = E / (3 * (1 - 2 * nu))
    Gs = E / (2 * (1 + nu))
    if porder == "equal":
        carriers = sorted({a for e in elems for a in e})
    else:
        carriers = sorted({e[i] for e in elems for i in TH_VERTICES[shape]})
    pmap = {a: s for s, a in enumerate(carriers)}
    npdof = len(carriers)
    Kg = np.zeros((ndm * nn, ndm * nn))
    Qg = np.zeros((ndm * nn, npdof))
    HTg = np.zeros((npdof, npdof))
    pts, ws = gp_contract(shape)
    for e in elems:
        coords = nodes[list(e)]
        Ke, Qe = elem_K_Q(shape, porder, coords, E, nu, alpha, 1.0, pts, ws)
        stab = (stab_alpha_auto(elem_size_h(shape, coords), Ks, Gs, 0.25)
                if stabilize else 0.0)
        pr = dict(biotAlpha=alpha, oneOverQbar=1.0e-3, rhoF=1.0, stabAlpha=stab,
                  thick=1.0, kbar=[1e-5] * ndm, drive=[0.0] * ndm)
        HTe = assemble_blocks(shape, porder, coords, pr, pts, ws)["HT"]
        u_idx = [ndm * a + i for a in e for i in range(ndm)]
        pl = ([pmap[a] for a in e] if porder == "equal"
              else [pmap[e[i]] for i in TH_VERTICES[shape]])
        Kg[np.ix_(u_idx, u_idx)] += Ke
        Qg[np.ix_(u_idx, pl)] += Qe
        HTg[np.ix_(pl, pl)] += HTe
    fixed = {ndm * a + i for a in boundary for i in range(ndm)}
    free_u = [d for d in range(ndm * nn) if d not in fixed]
    Kff = Kg[np.ix_(free_u, free_u)]
    Qf = Qg[free_u, :]
    Sp = Qf.T @ np.linalg.solve(Kff, Qf)
    if stabilize:
        Sp = Sp + HTg
    ev = np.linalg.eigvalsh(0.5 * (Sp + Sp.T))
    lam_max = np.max(np.abs(ev))
    return int(np.sum(ev < 1e-10 * lam_max)), npdof


def gate_infsup():
    """Gate (d): inf-sup smoke per ADR section 7 P0 — ALL pairs (0.E FIX 3).

    Equal-order pairs (Q4, T3, H8, BT6, BTET10): >=2 spurious near-zero Schur
    modes without stabilization; exactly 1 (the physical constant) with HT.
    Taylor-Hood pairs (BT6, BTET10): exactly 1 with NO stabilization.
    Patch sizes: 2D = 4x4 cells (T3/BT6 single-diagonal split, 32 elems);
    H8 = 4x4x4 (ADR-literal); BTET10 = 2x2x2 Kuhn (48 tets — the quadratic
    tet patch already has 125 nodes; free-u count 84 >= nNp-1 = 26, so the
    TH Schur rank is not patch-size-limited).  K from the oracle's own
    elastic B-matrix (plane strain in 2D; the D choice does not affect the
    null count).
    """
    print("\n[gate d] inf-sup smoke (structured patches, Schur S_p eigencount)")
    eq_pairs = [("Q4", 4), ("T3", 4), ("H8", 4), ("BT6", 4), ("BTET10", 2)]
    results = []
    for shape, nc in eq_pairs:
        n0, npd = _infsup_pair(shape, "equal", nc, stabilize=False)
        n1, _ = _infsup_pair(shape, "equal", nc, stabilize=True)
        print(f"  {shape:7s} equal patch={nc}^{NDM[shape]}: "
              f"no-stab {n0:3d}/{npd}  +Htilde {n1}/{npd}")
        assert n0 >= 2, f"{shape} equal no-stab: expected >=2, got {n0}"
        assert n1 == 1, f"{shape} equal stabilized: expected ==1, got {n1}"
        results.append((shape, "equal", n0, n1, npd))
    for shape, nc in (("BT6", 4), ("BTET10", 2)):
        n0, npd = _infsup_pair(shape, "th", nc, stabilize=False)
        print(f"  {shape:7s} TH    patch={nc}^{NDM[shape]}: "
              f"no-stab {n0:3d}/{npd}  (expect ==1, LBB-stable)")
        assert n0 == 1, f"{shape} TH: expected exactly 1, got {n0}"
        results.append((shape, "th", n0, None, npd))
    print("  -> inf-sup gate PASS (all pairs)")
    return results


def gate_fd_dynseepage():
    """Gate (e): FD quantification of the dropped dynamic-seepage tangent.

    Dynamic seepage drive = body - u_ddot ; the (deliberately dropped) LHS tangent
    is  G[b, a*ndm+i] = d fseep[b] / d uddot[a,i] = -INT dNp[b,i] k_i rhoF Nu[a] dv
    (chain: u_ddot(gp)_i = sum_a Nu[a] uddot[a,i]).  Verify analytic G by central
    finite differences and print its magnitude (this is the term Chan-1988 found
    insignificant and the book drops from the LHS).
    """
    print("\n[gate e] FD-quantified dropped dynamic-seepage tangent")
    shape, porder = "Q4", "equal"
    coords = np.array([[0.0, 0.0], [1.1, 0.05], [1.2, 1.0], [0.1, 0.9]])
    ndm, nNu = NDM[shape], NNODES[shape]
    nNp = n_pressure(shape, porder)
    kbar = np.array([3e-5, 1e-5])
    rhoF = 1.0
    body = np.array([0.0, -9.81])
    thick = 1.3
    pts, ws = gp_contract(shape)

    def fseep_of_uddot(uddot_flat):
        uddot = uddot_flat.reshape(nNu, ndm)
        f = np.zeros(nNp)
        for p, w in zip(pts, ws):
            N, dNdx, detJ, Np, dNpdx = eval_gp(shape, porder, coords, p)
            dv = w * detJ * thick
            u_gp = N @ uddot                    # acceleration at gp
            drive = body - u_gp
            for b in range(nNp):
                f[b] += np.sum(dNpdx[b, :] * kbar * rhoF * drive) * dv
        return f

    # analytic G
    G = np.zeros((nNp, nNu * ndm))
    for p, w in zip(pts, ws):
        N, dNdx, detJ, Np, dNpdx = eval_gp(shape, porder, coords, p)
        dv = w * detJ * thick
        for b in range(nNp):
            for a in range(nNu):
                for i in range(ndm):
                    G[b, a * ndm + i] += -dNpdx[b, i] * kbar[i] * rhoF * N[a] * dv

    # central FD
    u0 = np.zeros(nNu * ndm)
    Gfd = np.zeros_like(G)
    eps = 1e-6
    for c in range(nNu * ndm):
        up = u0.copy(); up[c] += eps
        um = u0.copy(); um[c] -= eps
        Gfd[:, c] = (fseep_of_uddot(up) - fseep_of_uddot(um)) / (2 * eps)
    fd_err = _relerr(Gfd, G)
    # magnitude relative to H block
    blk = assemble_blocks(shape, porder, coords,
                          dict(biotAlpha=1.0, oneOverQbar=1e-3, rhoF=rhoF,
                               stabAlpha=0.0, thick=thick, kbar=kbar,
                               drive=body), pts, ws)
    ratio = np.linalg.norm(G) / max(np.linalg.norm(blk["H"]), 1e-30)
    print(f"  ||G_dynseep||_F = {np.linalg.norm(G):.6e}")
    print(f"  ||G_dynseep||/||H|| = {ratio:.6e}  (dropped LHS tangent magnitude)")
    print(f"  analytic-vs-FD relerr = {fd_err:.3e} (<=1e-7)")
    assert fd_err <= 1e-7, f"dyn-seepage tangent FD mismatch {fd_err:.3e}"
    print("  -> dyn-seepage FD gate PASS")


def gate_hydrostatic():
    """Gate (g): hydrostatic identity H*p = f_seep (0.E FIX 7).

    For nodal p_a = -rhoF*g*y_a with drive = (0,-g), the interpolated pressure
    gradient is exactly rhoF*drive at every point (linear completeness of the
    isoparametric basis), so the H*p and fseep weak-form integrands coincide
    GP-by-GP — any quadrature rule, any (even distorted) geometry, any
    (diagonal) kbar.  Pins the end-to-end sign orientation of H vs fseep.
    """
    print("\n[gate g] hydrostatic identity H*p = fseep (sign orientation)")
    g = 9.81
    rhoF = 1.3
    q_unit = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    q_trap = np.array([[0.0, 0.0], [2.0, 0.0], [1.5, 1.0], [0.3, 1.0]])
    worst = 0.0
    for name, coords in (("unit square", q_unit), ("trapezoid", q_trap)):
        pr = dict(biotAlpha=1.0, oneOverQbar=1e-3, rhoF=rhoF, stabAlpha=0.0,
                  thick=1.4, kbar=[1e-5, 4e-6], drive=[0.0, -g])
        pts, ws = gp_contract("Q4")
        blk = assemble_blocks("Q4", "equal", coords, pr, pts, ws)
        p_nodal = -rhoF * g * coords[:, 1]
        resid = blk["H"] @ p_nodal - blk["fseep"]
        rel = np.max(np.abs(resid)) / max(np.max(np.abs(blk["fseep"])), 1e-30)
        worst = max(worst, rel)
        print(f"  Q4 {name:12s}: max|H*p - fseep|/max|fseep| = {rel:.3e}")
        assert rel <= 1e-12, f"hydrostatic identity broken on {name}: {rel:.3e}"
    print(f"  -> hydrostatic gate PASS (worst {worst:.2e})")


# ============================================================================ #
#  Canonical case geometry helpers                                              #
# ============================================================================ #

def _unit_cube():
    return np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def _bt6_coords(v):
    """v = 3x2 vertex coords ; append straight-side edge midpoints
    [m01, m12, m20] per node order [v0,v1,v2,m01,m12,m20]."""
    m01 = 0.5 * (v[0] + v[1])
    m12 = 0.5 * (v[1] + v[2])
    m20 = 0.5 * (v[2] + v[0])
    return np.vstack([v, m01, m12, m20])


def _btet10_coords(v):
    """v = 4x3 vertex coords ; append edge midpoints in node order
    [v0..v3, e01,e12,e02,e03,e23,e13]."""
    e01 = 0.5 * (v[0] + v[1])
    e12 = 0.5 * (v[1] + v[2])
    e02 = 0.5 * (v[0] + v[2])
    e03 = 0.5 * (v[0] + v[3])
    e23 = 0.5 * (v[2] + v[3])
    e13 = 0.5 * (v[1] + v[3])
    return np.vstack([v, e01, e12, e02, e03, e23, e13])


def _canonical_coords(shape):
    if shape == "Q4":
        return np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    if shape == "T3":
        return np.array([[0, 0], [1, 0], [0, 1]], float)
    if shape == "H8":
        return _unit_cube()
    if shape == "BT6":
        return _bt6_coords(np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]))
    if shape == "BTET10":
        # positively oriented under the pinned node<->L map (0.E FIX 2): the
        # conventionally right-handed ordering gives detJ<0 here — v0/v1 wound.
        return _btet10_coords(np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                                        [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    raise ValueError(shape)


# ============================================================================ #
#  Canonical case definitions + emission                                        #
# ============================================================================ #

# skeleton (drained) moduli for the McGann auto-stabilization alpha
_SKEL_E, _SKEL_NU, _ALPHA0 = 25000.0, 0.3, 0.25
_KS_SKEL = _SKEL_E / (3 * (1 - 2 * _SKEL_NU))
_GS_SKEL = _SKEL_E / (2 * (1 + _SKEL_NU))


def _stab_for(shape, porder, coords):
    if porder == "th":
        return 0.0
    h = elem_size_h(shape, coords)
    return stab_alpha_auto(h, _KS_SKEL, _GS_SKEL, _ALPHA0)


def _ooq_alpha1():
    # alpha = 1, Ks_grain -> inf : 1/Qbar = n/Kf
    return one_over_qbar(0.30, 2.2e6, 1.0, -1.0)


def _ooq_alpha079():
    # B3 Boone-Ingraffea parameters: n=0.19, Kf=3000, Ks=36000, alpha=0.79
    return one_over_qbar(0.19, 3000.0, 0.79, 36000.0)


def build_cases():
    """Return the canonical case list (>=12) covering the contract's spreads."""
    cases = []

    def add(name, shape, porder, coords, alpha, ooq, kbar, thick, drive, rhoF=1.0):
        stab = _stab_for(shape, porder, coords)
        params = dict(biotAlpha=alpha, oneOverQbar=ooq, rhoF=rhoF,
                      stabAlpha=stab, thick=thick, kbar=list(kbar),
                      drive=list(drive))
        cases.append((name, shape, porder, np.asarray(coords, float), params))

    a1 = _ooq_alpha1()
    a079 = _ooq_alpha079()
    iso2 = [1e-5, 1e-5]
    ani2 = [1e-5, 5e-6]
    iso3 = [1e-5, 1e-5, 1e-5]
    ani3 = [1e-5, 5e-6, 2e-5]
    d2 = [0.0, -9.81]
    d3 = [0.0, 0.0, -9.81]

    # --- Q4 (6 cases: baseline, trapezoid, alpha=0.79, aniso, thick=2.5,
    #     distorted+aniso+inclined-drive+rhoF!=1 [0.E FIX 4 a/b/c]) ---
    q_unit = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    q_trap = np.array([[0, 0], [2.0, 0], [1.5, 1.0], [0.3, 1.0]], float)
    add("Q4_unit_a1_iso_t1", "Q4", "equal", q_unit, 1.0, a1, iso2, 1.0, d2)
    add("Q4_trap_a1_iso_t1", "Q4", "equal", q_trap, 1.0, a1, iso2, 1.0, d2)
    add("Q4_unit_a079_iso_t1", "Q4", "equal", q_unit, 0.79, a079, iso2, 1.0, d2)
    add("Q4_unit_a1_aniso_t1", "Q4", "equal", q_unit, 1.0, a1, ani2, 1.0, d2)
    add("Q4_unit_a1_iso_t25", "Q4", "equal", q_unit, 1.0, a1, iso2, 2.5, d2)
    add("Q4_trap_a1_aniso_incl_rhoF25", "Q4", "equal", q_trap, 1.0, a1, ani2,
        1.0, [3.0, -9.81], rhoF=2.5)

    # --- T3 (3 cases: unit, skewed, skewed thick=2.5) ---
    t_unit = np.array([[0, 0], [1, 0], [0, 1]], float)
    t_skew = np.array([[0.2, 0.1], [1.3, 0.0], [0.4, 1.1]], float)
    add("T3_unit_a1_iso_t1", "T3", "equal", t_unit, 1.0, a1, iso2, 1.0, d2)
    add("T3_skew_a1_iso_t1", "T3", "equal", t_skew, 1.0, a1, iso2, 1.0, d2)
    add("T3_skew_a1_iso_t25", "T3", "equal", t_skew, 1.0, a1, iso2, 2.5, d2)

    # --- H8 (4 cases: cube, distorted, cube aniso,
    #     distorted+aniso+inclined-drive [0.E FIX 4 b/c in 3D]) ---
    cube = _unit_cube()
    dist = cube.copy()
    dist[6] = [1.2, 1.15, 1.1]
    dist[2] = [1.1, 1.05, -0.05]
    add("H8_cube_a1_iso", "H8", "equal", cube, 1.0, a1, iso3, 1.0, d3)
    add("H8_dist_a1_iso", "H8", "equal", dist, 1.0, a1, iso3, 1.0, d3)
    add("H8_cube_a1_aniso", "H8", "equal", cube, 1.0, a1, ani3, 1.0, d3)
    add("H8_dist_a1_aniso_incl", "H8", "equal", dist, 1.0, a1, ani3, 1.0,
        [2.0, 1.5, -9.81])

    # --- BT6 (4 cases: equal, TH, TH alpha=0.79,
    #     TH inclined-drive + rhoF!=1 [0.E FIX 4 a/b on a simplex/TH]) ---
    bv = np.array([[0.0, 0.0], [1.4, 0.2], [0.3, 1.1]])
    bt6 = _bt6_coords(bv)
    add("BT6_straight_equal_a1_iso_t1", "BT6", "equal", bt6, 1.0, a1, iso2, 1.0, d2)
    add("BT6_straight_th_a1_iso_t1", "BT6", "th", bt6, 1.0, a1, iso2, 1.0, d2)
    add("BT6_straight_th_a079_iso_t1", "BT6", "th", bt6, 0.79, a079, iso2, 1.0, d2)
    add("BT6_straight_th_a1_iso_incl_rhoF25", "BT6", "th", bt6, 1.0, a1, iso2,
        1.0, [4.0, -9.81], rhoF=2.5)

    # --- BTET10 (2 cases: equal, TH) straight-sided ---
    # vertices wound for detJ>0 under the pinned node<->L map (0.E FIX 2)
    tv = np.array([[1.2, 0.1, 0.0], [0.1, 0.0, 0.0],
                   [0.2, 1.1, 0.1], [0.0, 0.2, 1.3]])
    btet = _btet10_coords(tv)
    add("BTET10_straight_equal_a1_iso", "BTET10", "equal", btet, 1.0, a1, iso3, 1.0, d3)
    add("BTET10_straight_th_a1_iso", "BTET10", "th", btet, 1.0, a1, iso3, 1.0, d3)

    return cases


def _fmt(x):
    return format(float(x), ".17g")


def _row(vals):
    return " ".join(_fmt(v) for v in vals)


def emit(path):
    cases = build_cases()
    lines = []
    for name, shape, porder, coords, params in cases:
        ndm = NDM[shape]
        nNodes = NNODES[shape]
        pc, wc = gp_contract(shape)
        blk = assemble_blocks(shape, porder, coords, params, pc, wc)
        totalDof, uOff, pOff = dof_map(shape, porder)
        Q, H, S, HT, fseep = blk["Q"], blk["H"], blk["S"], blk["HT"], blk["fseep"]
        lines.append(f"CASE {name} SHAPE {shape} PORDER {porder} "
                     f"NDM {ndm} NNODES {nNodes}")
        lines.append("COORDS " + _row(coords.reshape(-1)))
        lines.append("PARAMS " + _row([params["biotAlpha"], params["oneOverQbar"],
                                        params["rhoF"], params["stabAlpha"],
                                        params["thick"]] + list(params["kbar"])))
        lines.append("DRIVE " + _row(params["drive"]))
        lines.append("DOFMAP " + _row([totalDof] + uOff + pOff))
        lines.append(f"Q {Q.shape[0]} {Q.shape[1]} " + _row(Q.reshape(-1)))
        lines.append(f"H {H.shape[0]} {H.shape[1]} " + _row(H.reshape(-1)))
        lines.append(f"S {S.shape[0]} {S.shape[1]} " + _row(S.reshape(-1)))
        lines.append(f"HT {HT.shape[0]} {HT.shape[1]} " + _row(HT.reshape(-1)))
        lines.append(f"FSEEP {fseep.shape[0]} " + _row(fseep))
        lines.append("END")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n[emit] wrote {len(cases)} cases -> {path}")
    return cases


def gate_affine_detj():
    """Sanity: straight-sided Bezier + parallelogram/affine shapes have constant
    detJ (validates midside-at-midpoint placement => affine map)."""
    print("\n[gate f] affine-map detJ constancy (straight Bezier / T3)")
    worst = 0.0
    for shape in ("T3", "BT6", "BTET10"):
        coords = _canonical_coords(shape)
        pc, wc = gp_contract(shape)
        blk = assemble_blocks(shape, "equal", coords,
                              dict(biotAlpha=1.0, oneOverQbar=1e-3, rhoF=1.0,
                                   stabAlpha=0.0, thick=1.0,
                                   kbar=[1e-5] * NDM[shape],
                                   drive=[0.0] * NDM[shape]), pc, wc)
        dj = blk["detJ"]
        spread = np.max(np.abs(dj - dj[0])) / max(abs(dj[0]), 1e-30)
        worst = max(worst, spread)
        print(f"  {shape:8s} detJ spread = {spread:.2e}")
        assert spread <= 1e-13, f"{shape} detJ not constant ({spread:.2e})"
    print(f"  -> affine detJ gate PASS (worst {worst:.2e})")


def self_check():
    print("=" * 74)
    print("LadrunoUP numpy oracle - self-check (ADR-71 P0, WP 0.B)")
    print("=" * 74)
    if not _HAVE_SCIPY:
        print("WARNING: scipy unavailable - consolidation gate uses expm fallback only")
    worst = gate_overkill(None)
    gate_symmetry()
    gate_affine_detj()
    gate_consolidation()
    gate_infsup()
    gate_fd_dynseepage()
    gate_hydrostatic()
    # permanent case-level asserts (0.E FIX 2): finite blocks, per-GP detJ > 0,
    # S/H/HT positive-semidefinite — an inverted/mis-wound canonical case can
    # never ship silently again.
    cases = build_cases()
    print(f"\n[cases] per-case asserts: detJ>0 per GP, S/H/HT PSD, finite blocks")
    for name, shape, porder, coords, params in cases:
        pc, wc = gp_contract(shape)
        blk = assemble_blocks(shape, porder, coords, params, pc, wc)
        for k in ("Q", "H", "S", "HT", "fseep"):
            assert np.all(np.isfinite(blk[k])), f"{name} {k} non-finite"
        assert np.all(blk["detJ"] > 0.0), \
            f"{name}: non-positive detJ at a GP: {blk['detJ']}"
        for mkey in ("S", "H", "HT"):
            M = blk[mkey]
            ev = np.linalg.eigvalsh(0.5 * (M + M.T))
            scale = max(np.max(np.abs(ev)), 1e-30)
            assert ev.min() >= -1e-12 * scale, \
                f"{name} {mkey} not PSD (min eig {ev.min():.3e})"
        if shape == "BTET10":
            vol_quad = float(np.sum(np.asarray(wc) * blk["detJ"]))
            v = coords[:4]
            vol_true = abs(np.linalg.det(
                np.column_stack([v[1] - v[0], v[2] - v[0], v[3] - v[0]]))) / 6.0
            print(f"  {name}: sum(w*detJ) = {vol_quad:.12g} "
                  f"(= +volume {vol_true:.12g})")
            assert vol_quad > 0.0 and abs(vol_quad - vol_true) <= 1e-12 * vol_true
    print(f"[cases] {len(cases)} canonical cases build cleanly, all asserts hold")
    print("\n" + "=" * 74)
    print(f"SELF-CHECK GREEN  (worst exact-combo relerr {worst:.2e} <= 1e-12)")
    print("=" * 74)


def main(argv=None):
    ap = argparse.ArgumentParser(description="LadrunoUP numpy oracle (ADR-71 P0)")
    ap.add_argument("--emit", metavar="PATH", default=None,
                    help="run self-check then emit the cases file")
    args = ap.parse_args(argv)
    self_check()
    if args.emit:
        emit(args.emit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
