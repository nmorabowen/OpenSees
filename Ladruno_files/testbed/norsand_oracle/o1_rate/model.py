"""Constitutive functions of LadrunoNORSAND for the rate oracle (equation sheet 144a).

Everything here is written COORDINATE-FREE on 3x3 tensors (the principal-space formulas of
the sheet are isotropic tensor functions; their tensor forms are used so that no eigenvector
of a (possibly repeated) spectrum is ever needed on the integration path):

  energy          (S.2)-(S.5): sigma(eps^e) and the 4th-order elastic tangent a^e, written
                  as the tensor derivative of sigma = p 1 + (2/3)(q/eps_s) e.  Its principal
                  components are (S.3) and its spin part is (S.33) (cross-checked by
                  `elastic_tangent_spectral`, which IS (S.3)+(S.33)).
  invariants      (S.1), (S.6): p, s, R, q, n^, y, dy/dsigma (tensor form of y_a).
  zeta            (S.8)-(S.11): zeta(theta), zeta_y with the corner branch (S.9).
  yield/potential (S.12)-(S.17), cap (S.35)-(S.36), vertex rule (sheet 3.2).
  CSL, pi_i*, H   (S.22)-(S.23), (S.25), (S.41).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .params import Params

SQ23 = math.sqrt(2.0 / 3.0)
SQ32 = math.sqrt(1.5)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)
IxI = np.einsum("ij,kl->ijkl", I3, I3)
ISYM = 0.5 * (np.einsum("ik,jl->ijkl", I3, I3) + np.einsum("il,jk->ijkl", I3, I3))
IDEV = ISYM - IxI / 3.0

R_TOL_REL = 1.0e-8          # vertex rule threshold R < 1e-8 |p|  (sheet 3.2)
CORNER_TOL = 1.0e-8         # |sin 3 theta| below which the corner branch (S.9) is used

# Mandel packing of symmetric tensors: [xx, yy, zz, sqrt2 xy, sqrt2 yz, sqrt2 xz]
_MI = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]
_MW = np.array([1.0, 1.0, 1.0, math.sqrt(2.0), math.sqrt(2.0), math.sqrt(2.0)])


def t2m(T):
    return np.array([T[i, j] for (i, j) in _MI]) * _MW


def m2t(m):
    T = np.empty((3, 3))
    for k, (i, j) in enumerate(_MI):
        T[i, j] = T[j, i] = m[k] / _MW[k]
    return T


def c2m(C):
    """4th-order minor-symmetric tensor -> 6x6 Mandel matrix (sigma_m = C_m eps_m)."""
    out = np.empty((6, 6))
    for a, (i, j) in enumerate(_MI):
        for b, (k, l) in enumerate(_MI):
            out[a, b] = _MW[a] * _MW[b] * C[i, j, k, l]
    return out


def outer(A, B):
    return np.einsum("ij,kl->ijkl", A, B)


def ddot(A, B):
    return float(np.sum(A * B))


# ---------------------------------------------------------------------------------------
# BA06 energy (S.4)-(S.5), tensor form
# ---------------------------------------------------------------------------------------
@dataclass
class Elastic:
    sig: np.ndarray
    p: float
    q: float
    eps_v: float
    eps_s: float
    D11: float
    D12: float
    D22: float
    q_over_es: float
    a4: np.ndarray | None


def energy(ee, P: Params, tangent=True) -> Elastic:
    ev = float(np.trace(ee))
    e = ee - (ev / 3.0) * I3
    ne = float(np.linalg.norm(e))
    es = SQ23 * ne
    E = math.exp(-(ev - P.eps_v0) / P.kappa_hat)
    p = P.p0 * E * (1.0 + 1.5 * P.alpha0 * es * es / P.kappa_hat)
    mu = P.mu0 - P.alpha0 * P.p0 * E
    q = 3.0 * mu * es
    qes = 3.0 * mu                    # q / eps_s, exact for BA06 incl. the eps_s -> 0 limit
    D11 = -p / P.kappa_hat
    D22 = 3.0 * mu
    D12 = 3.0 * P.p0 * P.alpha0 * es * E / P.kappa_hat
    sig = p * I3 + (2.0 / 3.0) * qes * e
    a4 = None
    if tangent:
        n = e / ne if ne > 0.0 else np.zeros((3, 3))
        a4 = (D11 * IxI + SQ23 * D12 * (outer(I3, n) + outer(n, I3))
              + (2.0 / 3.0) * (D22 - qes) * outer(n, n) + (2.0 / 3.0) * qes * IDEV)
    return Elastic(sig, p, q, ev, es, D11, D12, D22, qes, a4)


def elastic_tangent_spectral(ee, P: Params, tol=1e-10):
    """(S.3) principal a^e_ab assembled by (S.33) (small strain).  Used only as a check."""
    w, V = np.linalg.eigh(ee)
    el = energy(ee, P, tangent=False)
    ev = w.sum()
    ea = w - ev / 3.0
    ne = np.linalg.norm(ea)
    nh = ea / ne if ne > 0 else np.zeros(3)
    d = np.ones(3)
    ab = (el.D11 * np.outer(d, d) + SQ23 * el.D12 * (np.outer(d, nh) + np.outer(nh, d))
          + (2.0 / 3.0) * el.D22 * np.outer(nh, nh)
          + (2.0 * el.q / (3.0 * el.eps_s) if el.eps_s > 0 else (2.0 / 3.0) * el.D22)
          * (np.eye(3) - np.outer(d, d) / 3.0 - np.outer(nh, nh)))
    sa = np.array([V[:, a] @ el.sig @ V[:, a] for a in range(3)])
    C = np.zeros((3, 3, 3, 3))
    m = [np.outer(V[:, a], V[:, a]) for a in range(3)]
    for a in range(3):
        for b in range(3):
            C += ab[a, b] * outer(m[a], m[b])
            if a != b:
                mab = np.outer(V[:, a], V[:, b])
                mba = np.outer(V[:, b], V[:, a])
                if abs(w[a] - w[b]) > tol:
                    g = (sa[a] - sa[b]) / (w[a] - w[b])
                else:
                    g = ab[a, a] - ab[a, b]
                C += 0.5 * g * (outer(mab, mab) + outer(mab, mba))
    return C


def elastic_strain_from_stress(sig, P: Params):
    """Inverse of the energy: eps^e co-axial with sig (Newton on (eps_v, eps_s))."""
    p = float(np.trace(sig)) / 3.0
    s = sig - p * I3
    R = float(np.linalg.norm(s))
    q = SQ32 * R
    n = s / R if R > 0 else np.zeros((3, 3))
    if not p < 0.0:
        raise ValueError("initial mean stress must be < 0")
    ev = P.eps_v0 - P.kappa_hat * math.log(p / P.p0)
    es = q / (3.0 * P.mu0)
    for _ in range(100):
        ee = (ev / 3.0) * I3 + SQ32 * es * n
        el = energy(ee, P, tangent=False)
        r = np.array([el.p - p, el.q - q])
        if abs(r[0]) <= 1e-14 * abs(p) and abs(r[1]) <= 1e-14 * max(abs(p), q):
            break
        J = np.array([[el.D11, el.D12], [el.D12, el.D22]])
        dx = np.linalg.solve(J, -r)
        ev += dx[0]
        es += dx[1]
    return (ev / 3.0) * I3 + SQ32 * es * n


# ---------------------------------------------------------------------------------------
# zeta(theta, rho): GA (S.10) and WW (S.11), with zeta_y (S.8) and the corner branch (S.9)
# ---------------------------------------------------------------------------------------
def _ww_z(c, rho):
    A = 4.0 * (1.0 - rho * rho)
    B = 2.0 * rho - 1.0
    s = np.sqrt(A * c * c + 5.0 * rho * rho - 4.0 * rho)
    return (A * c * c + B * B) / (2.0 * (1.0 - rho * rho) * c + B * s)


def _ww_dzdc(c, rho):
    A = 4.0 * (1.0 - rho * rho)
    B = 2.0 * rho - 1.0
    if B == 0.0:          # rho = 1/2: zeta = 2c exactly (a VERTEX at theta = pi/3, see README)
        return 2.0 + 0.0 * c
    s = np.sqrt(A * c * c + 5.0 * rho * rho - 4.0 * rho)
    num = A * c * c + B * B
    den = 2.0 * (1.0 - rho * rho) * c + B * s
    dnum = 2.0 * A * c
    dden = 2.0 * (1.0 - rho * rho) + B * A * c / s
    return (dnum * den - num * dden) / (den * den)


@lru_cache(maxsize=64)
def ww_corner_constants(rho):
    """zeta''(0), zeta''(pi/3) for WW (theta-derivatives), via the exact chain rule
    zeta'' = zeta_cc sin^2 - zeta_c cos, zeta_cc by complex step (machine precision)."""
    h = 1e-30
    zc1 = float(_ww_dzdc(1.0, rho))
    zcc_half = float(np.imag(_ww_dzdc(0.5 + 1j * h, rho)) / h)
    zc_half = float(_ww_dzdc(0.5, rho))
    zpp0 = -zc1                                        # sin 0 = 0, cos 0 = 1
    zpp60 = zcc_half * 0.75 - zc_half * 0.5            # sin^2 = 3/4, cos = 1/2
    return zpp0, zpp60


def zeta_fun(theta, y, rho, kind):
    """Return (zeta, zeta_y).  theta in [0, pi/3] computed from the same y (rule (i), 3.1)."""
    if kind == "GA":
        z = ((1.0 + rho) + SQ6 * (1.0 - rho) * y) / (2.0 * rho)
        return z, SQ6 * (1.0 - rho) / (2.0 * rho)
    c = math.cos(theta)
    z = float(_ww_z(c, rho))
    if rho == 0.5:        # zeta = 2c, zeta_y = 2 sqrt6/(12c^2-3): unbounded at theta = pi/3
        den = 12.0 * c * c - 3.0
        return z, (2.0 * SQ6 / den if den > 0.0 else float("inf"))
    s3 = math.sin(3.0 * theta)
    if abs(s3) < CORNER_TOL:
        zpp0, zpp60 = ww_corner_constants(rho)
        zy = -SQ6 * zpp0 / 9.0 if theta < math.pi / 6.0 else SQ6 * zpp60 / 9.0
    else:
        zp = -math.sin(theta) * float(_ww_dzdc(c, rho))
        zy = -(2.0 / SQ6) * zp / s3
    return z, zy


# ---------------------------------------------------------------------------------------
# yield function, potential (reading A), cap, CSL, pi_i*, H
# ---------------------------------------------------------------------------------------
@dataclass
class Plastic:
    p: float
    q: float
    R: float
    theta: float
    vertex: bool
    zeta: float
    zeta_y: float
    zeta_b: float
    zeta_b_y: float
    eta: float
    F: float
    Fp: float
    Fpi: float
    f: np.ndarray          # dF/dsigma (S.14)
    qu: np.ndarray         # uncapped flow direction (S.17)
    Omega_u: float
    w: float               # cap weight (1 = no cap)
    qflow: np.ndarray      # capped flow direction (S.35)
    Omega: float           # w Omega_u (S.36)
    psi: float
    pistar: float
    B: float               # base of (S.23) (N > 0), 1 for N = 0
    H: float               # hardening modulus (S.41)
    dilatancy: float       # D = eps_p_v/eps_p_s (nan when Omega = 0)


def eta_of(p, pi, P: Params):
    """(S.12) with F_p and F_pi of (S.13)."""
    x = p / pi
    if P.N == 0.0:
        eta = P.M * (1.0 + math.log(pi / p))
    else:
        eta = (P.M / P.N) * (1.0 - (1.0 - P.N) * x ** (P.N / (1.0 - P.N)))
    Fp = (eta - P.M) / (1.0 - P.N)
    Fpi = P.M * x ** (1.0 / (1.0 - P.N))
    return eta, Fp, Fpi


def psi_i(v, pi, P: Params):
    """State parameter at the image pressure (sheet 6)."""
    if P.csl_mode == "paper":
        return v - P.v_c0 + P.lambda_tilde * math.log(-pi)
    return (v - 1.0) - P.e0 + P.lambda_c * (-pi / P.p_a) ** P.xi


def pi_star(p, Omega, psi, P: Params):
    """(S.23): returns (pi_i*, B).  B <= 0 -> (nan, B) (guard of sheet 7)."""
    cb = P.chi_bar
    if P.N == 0.0:
        return p * math.exp(SQ23 * cb * psi * Omega / P.M), 1.0
    B = 1.0 - SQ23 * cb * psi * Omega * P.N / P.M
    if B <= 0.0:
        return float("nan"), B
    return p * B ** ((P.N - 1.0) / P.N), B


def cap_weight(eta, P: Params):
    if P.cap == "none":
        return 1.0
    if P.cap == "planar":
        return 1.0 if eta >= P.c1 * P.M else 0.0
    e1, e2 = P.c1 * P.M, P.c2 * P.M
    t = min(max((eta - e1) / (e2 - e1), 0.0), 1.0)
    if P.cap_blend == "cubic":
        return t * t * (3.0 - 2.0 * t)
    return t ** 3 * (10.0 - 15.0 * t + 6.0 * t * t)


def invariants(sig):
    p = float(np.trace(sig)) / 3.0
    s = sig - p * I3
    R = float(np.linalg.norm(s))
    return p, s, R


def plastic(sig, pi, v, P: Params, w_override=None) -> Plastic:
    p, s, R = invariants(sig)
    if not (p < 0.0 and pi < 0.0):
        raise FloatingPointError(f"p = {p}, pi_i = {pi}: need both < 0")
    q = SQ32 * R
    vertex = R < R_TOL_REL * abs(p)
    if R > 0.0:
        s2 = s @ s
        tr3 = float(np.sum(s * s2))
        y = tr3 / R ** 3
        Y = 3.0 * s2 / R ** 3 - I3 / R - 3.0 * tr3 * s / R ** 5      # tensor form of (S.6)
        nh = s / R
    else:
        y = -1.0 / SQ6          # arbitrary (zeta multiplies q = 0)
        Y = np.zeros((3, 3))
        nh = np.zeros((3, 3))
    theta = math.acos(min(1.0, max(-1.0, SQ6 * y))) / 3.0
    z, zy = zeta_fun(theta, y, P.rho, P.zeta)
    zb, zby = zeta_fun(theta, y, P.rho_bar, P.zeta)
    if vertex:                   # sheet 3.2
        nh = np.zeros((3, 3))
        Y = np.zeros((3, 3))
    eta, Fp, Fpi = eta_of(p, pi, P)
    F = z * q + p * eta
    beta = P.beta
    f = (Fp / 3.0) * I3 + SQ32 * z * nh + zy * q * Y
    qu = (beta * Fp / 3.0) * I3 + SQ32 * zb * nh + zby * q * Y
    Om_u = 0.0 if vertex else math.sqrt(1.5 * zb * zb + (zby * q) ** 2 * float(np.sum(Y * Y)))
    w = cap_weight(eta, P) if w_override is None else w_override
    qflow = -I3 / 3.0 + w * (qu + I3 / 3.0)
    Om = w * Om_u
    ps = psi_i(v, pi, P)
    pst, B = pi_star(p, Om, ps, P)
    H = -Fpi * SQ23 * P.h * (pst - pi) * Om
    dil = SQ32 * float(np.trace(qflow)) / Om if Om > 0 else float("nan")
    return Plastic(p, q, R, theta, vertex, z, zy, zb, zby, eta, F, Fp, Fpi, f, qu, Om_u, w,
                   qflow, Om, ps, pst, B, H, dil)


def continuum_tangent(ee, pi, v, P: Params, plastic_branch=True):
    """a^e, or the loading-branch a^ep of (S.42).  Returns (C4, info)."""
    el = energy(ee, P)
    Ae = el.a4
    if not plastic_branch:
        return Ae, dict(den=float("nan"))
    pq = plastic(el.sig, pi, v, P)
    Aq = np.einsum("ijkl,kl->ij", Ae, pq.qflow)
    fA = np.einsum("ij,ijkl->kl", pq.f, Ae)
    den = ddot(pq.f, Aq) + pq.H
    if not den > 0.0:
        raise FloatingPointError(f"f:a^e:q + H = {den} <= 0 (loss of uniqueness, S.42)")
    return Ae - outer(Aq, fA) / den, dict(den=den, H=pq.H)
