"""O2 algorithmic oracle kernel for LadrunoNORSAND: sheet 144a §1-§10, §13-§14.

Everything here is the closed-form algebra of the equation sheet, transcribed, in
principal space (3-vectors), plus the backward-Euler return map of sheet §9 (AB06 Box 2)
and the closed-form consistent tangent of §9.3-9.5. Compression negative throughout.

The C++ kernel must reproduce this module to ~1e-10, so every tolerance and every
branch below is part of the kernel's contract (see the CONTRACT constants and README).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .params import Params

# --------------------------------------------------------------------------------------
# NUMERICAL CONTRACT (the kernel reproduces these exactly)
# --------------------------------------------------------------------------------------
RES_TOL = 1.0e-12          # local Newton: ||(r1,r2,r3, r4/|p0|)||_2 <= RES_TOL  (S.29, scaling note §9.1)
MAX_LOCAL_ITERS = 30       # local Newton iterations before refusal "local_noconv"
MAX_LINESEARCH = 10        # halvings of the Newton step when the scaled residual does not decrease
                           # (a nested failure at the trial iterate counts as "does not decrease")
PI_TOL_REL = 1.0e-12       # nested pi_i Newton: |r(pi_i)| <= PI_TOL_REL * |pi_i,n|  (S.27)
MAX_PI_ITERS = 50          # nested safeguarded-Newton iterations (after bracketing) before failure
PI_SCAN_REL = 1.0e-3       # nested: scan step |d pi_i| = PI_SCAN_REL*|pi_i,n| from pi_i,n toward the root
                           # (well below the ~0.01|pi_i,n| width of the fold of r(pi_i) in the cap ramp)
PI_SCAN_MAX = 1000         # scan steps (travel <= PI_SCAN_MAX*PI_SCAN_REL*|pi_i,n| = |pi_i,n|) before 'pi_nobracket'
MAX_SUBSTEP_HALVINGS = 8   # a refused increment is retried as 2, 4, ..., 2^8 = 256 equal sub-increments (api.step)
R_TOL_REL = 1.0e-8         # vertex rule (§3.2): R = ||xi|| < R_TOL_REL*|p|  =>  hydrostatic vertex
CORNER_SIN3T = 1.0e-8      # corner branch of (S.9): |sin 3theta| < CORNER_SIN3T
F_TRIAL_TOL_REL = 1.0e-10  # trial state is plastic iff F(sigma_tr, pi_i,n) > F_TRIAL_TOL_REL*|p0|
EPS_S_TOL = 1.0e-14        # eps^e_s below this: n_hat^e := 0 and q/eps_s := D22 in (S.3)
REPEATED_EIG_TOL = 1.0e-10 # |eps~_a - eps~_b| below this: repeated-eigenvalue limit of g_ab in (S.33)
REPEATED_STRETCH_TOL = 1.0e-10  # |lam~_a - lam~_b| below this: limit of gamma~_ab in (S.34)
PI_MAX_NEG = 0.0           # any pi_i iterate >= PI_MAX_NEG (i.e. non-negative) is a failure

SQ23 = math.sqrt(2.0 / 3.0)
SQ32 = math.sqrt(1.5)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)
ONES = np.ones(3)


class EvalError(Exception):
    """Raised INSIDE the residual evaluation only; caught by the local Newton and
    turned into a line-search backtrack or a refusal flag. Never escapes run_path."""


# --------------------------------------------------------------------------------------
# §4 zeta(theta, rho): shape functions in the y-form of (S.8)-(S.9)
# --------------------------------------------------------------------------------------
def ww_theta(theta: float, rho: float):
    """Willam-Warnke (S.11): zeta, zeta'(theta), zeta''(theta) by the quotient rule."""
    A = 4.0 * (1.0 - rho * rho)
    B = 2.0 * rho - 1.0
    c = math.cos(theta)
    s = math.sin(theta)
    Nm = A * c * c + B * B
    Nm1 = 2.0 * A * c
    Nm2 = 2.0 * A
    two = 2.0 * (1.0 - rho * rho)
    if abs(B) < 1e-15:            # rho = 1/2 exactly: the B*S term vanishes (S -> 0 at theta = pi/3).
                                  # Unreachable through validated Params (WW needs rho > 1/2, owner
                                  # decision 2026-10-01); kept so the raw function stays defined.
        Dn, Dn1, Dn2 = two * c, two, 0.0
    else:
        S = math.sqrt(A * c * c + 5.0 * rho * rho - 4.0 * rho)
        S1 = A * c / S
        S2 = A / S - A * A * c * c / S ** 3
        Dn = two * c + B * S
        Dn1 = two + B * S1
        Dn2 = B * S2
    z = Nm / Dn
    num1 = Nm1 * Dn - Nm * Dn1
    zc = num1 / Dn ** 2
    zcc = (Nm2 * Dn - Nm * Dn2) / Dn ** 2 - 2.0 * num1 * Dn1 / Dn ** 3
    z1 = -zc * s                       # d zeta / d theta
    z2 = zcc * s * s - zc * c          # d2 zeta / d theta2
    return z, z1, z2


def ga_theta(theta: float, rho: float):
    """Gudehus-Argyris (S.10) in theta."""
    z = ((1.0 + rho) + (1.0 - rho) * math.cos(3.0 * theta)) / (2.0 * rho)
    z1 = -3.0 * (1.0 - rho) * math.sin(3.0 * theta) / (2.0 * rho)
    z2 = -9.0 * (1.0 - rho) * math.cos(3.0 * theta) / (2.0 * rho)
    return z, z1, z2


def zeta_y(theta: float, rho: float, kind: str):
    """zeta, zeta_y, zeta_yy of (S.8)-(S.10) at theta (theta computed from y, sheet §3.1).
    sin3theta, cos3theta are taken from the SAME theta as zeta', zeta'' (rule (i) of §3.1)."""
    if kind == "GA":
        z, _, _ = ga_theta(theta, rho)
        return z, SQ6 * (1.0 - rho) / (2.0 * rho), 0.0
    z, z1, z2 = ww_theta(theta, rho)
    s3 = math.sin(3.0 * theta)
    if abs(s3) < CORNER_SIN3T:                       # corner branch (S.9)
        if theta < math.pi / 6.0:
            _, _, z2c = ww_theta(0.0, rho)
            return z, -SQ6 * z2c / 9.0, 0.0
        _, _, z2c = ww_theta(math.pi / 3.0, rho)
        return z, SQ6 * z2c / 9.0, 0.0
    c3 = math.cos(3.0 * theta)
    zy = -(2.0 / SQ6) / s3 * z1
    zyy = (2.0 / 3.0) / (s3 * s3) * (z2 - 3.0 * z1 * c3 / s3)
    return z, zy, zyy


def zeta_theta(theta: float, rho: float, kind: str):
    return ga_theta(theta, rho) if kind == "GA" else ww_theta(theta, rho)


# --------------------------------------------------------------------------------------
# §2 BA06 energy in principal elastic strains
# --------------------------------------------------------------------------------------
@dataclass
class Elastic:
    eps_v: float
    eps_s: float
    nhat_e: np.ndarray
    p: float
    q: float
    D11: float
    D12: float
    D22: float
    sig: np.ndarray      # principal stresses (3,)
    ae: np.ndarray       # a^e_ab (3,3)  (S.3)


def elastic(P: Params, eps_e: np.ndarray) -> Elastic:
    ev = float(eps_e.sum())
    e = eps_e - ev / 3.0
    ne = float(np.linalg.norm(e))
    es = SQ23 * ne
    nh = e / ne if es > EPS_S_TOL else np.zeros(3)
    om = -(ev - P.eps_v0) / P.kappa_hat
    E = math.exp(om)
    p = P.p0 * E * (1.0 + 1.5 * P.alpha0 / P.kappa_hat * es * es)
    q = 3.0 * (P.mu0 - P.alpha0 * P.p0 * E) * es
    D11 = -p / P.kappa_hat
    D22 = 3.0 * P.mu0 - 3.0 * P.alpha0 * P.p0 * E
    D12 = 3.0 * P.p0 * P.alpha0 * es / P.kappa_hat * E
    sig = p * ONES + SQ23 * q * nh
    ratio = q / es if es > EPS_S_TOL else D22            # eps_s -> 0 limit (S.3 note)
    ae = (D11 * np.outer(ONES, ONES)
          + SQ23 * D12 * (np.outer(ONES, nh) + np.outer(nh, ONES))
          + (2.0 / 3.0) * D22 * np.outer(nh, nh)
          + (2.0 * ratio / 3.0) * (I3 - np.outer(ONES, ONES) / 3.0 - np.outer(nh, nh)))
    return Elastic(ev, es, nh, p, q, D11, D12, D22, sig, ae)


def energy_psi(P: Params, eps_e: np.ndarray) -> float:
    """Psi (S.4), for the closed-loop check."""
    ev = float(eps_e.sum())
    e = eps_e - ev / 3.0
    es = SQ23 * float(np.linalg.norm(e))
    om = -(ev - P.eps_v0) / P.kappa_hat
    Pt = -P.p0 * P.kappa_hat * math.exp(om)
    mu = P.mu0 + P.alpha0 / P.kappa_hat * Pt
    return Pt + 1.5 * mu * es * es


def invert_elastic(P: Params, sig: np.ndarray) -> np.ndarray:
    """Principal elastic strains from principal stresses (Newton on (eps_v, eps_s) with the
    2x2 Hessian; closed form when alpha0 = 0). Used by initial_state only."""
    p = float(sig.mean())
    xi = sig - p
    R = float(np.linalg.norm(xi))
    q = SQ32 * R
    nh = xi / R if R > 0.0 else np.zeros(3)
    if p >= 0.0:
        raise ValueError("initial stress must have p < 0")
    ev = P.eps_v0 - P.kappa_hat * math.log(p / P.p0)
    es = q / (3.0 * P.mu0)
    if P.alpha0 != 0.0:
        for _ in range(100):
            el = elastic(P, ev * ONES / 3.0 + SQ32 * es * nh if es > 0 else ev * ONES / 3.0)
            r = np.array([el.p - p, el.q - q])
            if np.linalg.norm(r) <= 1e-13 * abs(p):
                break
            H = np.array([[el.D11, el.D12], [el.D12, el.D22]])
            d = np.linalg.solve(H, r)
            ev -= d[0]
            es -= d[1]
    return ev * ONES / 3.0 + SQ32 * es * nh


# --------------------------------------------------------------------------------------
# §3 invariants and their derivatives in principal space (S.1), (S.6), (S.7), vertex rule §3.2
# --------------------------------------------------------------------------------------
@dataclass
class Invariants:
    p: float
    q: float
    R: float
    vertex: bool
    nhat: np.ndarray
    y: float
    theta: float
    y_a: np.ndarray
    nhat_ab: np.ndarray
    y_ab: np.ndarray


def invariants(sig: np.ndarray) -> Invariants:
    p = float(sig.sum()) / 3.0
    xi = sig - p
    R = float(np.linalg.norm(xi))
    if R < R_TOL_REL * abs(p):
        # vertex rule (§3.2): n_hat, y_a, n_hat_ab, y_ab := 0 ; q := 0 in F (F = p eta)
        return Invariants(p, 0.0, R, True, np.zeros(3), 0.0, math.pi / 3.0,
                          np.zeros(3), np.zeros((3, 3)), np.zeros((3, 3)))
    q = SQ32 * R
    nh = xi / R
    S3 = float(np.sum(xi ** 3))
    y = S3 / R ** 3
    arg = min(1.0, max(-1.0, SQ6 * y))
    theta = math.acos(arg) / 3.0
    y_a = 3.0 * xi * xi / R ** 3 - 3.0 * S3 * xi / R ** 5 - 1.0 / R
    P1 = I3 - np.outer(ONES, ONES) / 3.0
    nh_ab = (P1 - np.outer(nh, nh)) / R
    y_ab = (6.0 * np.diag(xi) / R ** 3
            - 3.0 * S3 / R ** 5 * (P1 - 5.0 * np.outer(xi, xi) / R ** 2)
            + (np.outer(ONES, xi) + np.outer(xi, ONES)) / R ** 3
            - 9.0 * (np.outer(xi, xi * xi) + np.outer(xi * xi, xi)) / R ** 5)
    return Invariants(p, q, R, False, nh, y, theta, y_a, nh_ab, y_ab)


# --------------------------------------------------------------------------------------
# §5 yield function F, its p / pi_i derivatives (S.12)-(S.13)
# --------------------------------------------------------------------------------------
def eta_of(P: Params, p: float, pi: float) -> float:
    if P.N == 0.0:
        return P.M * (1.0 + math.log(pi / p))
    return (P.M / P.N) * (1.0 - (1.0 - P.N) * (p / pi) ** (P.N / (1.0 - P.N)))


def pi_of_eta(P: Params, p: float, eta: float) -> float:
    """Inverse of (S.12) (BA06 2.8): pi_i on the surface through (p, eta)."""
    if P.N == 0.0:
        return p * math.exp(eta / P.M - 1.0)
    d = 1.0 - eta * P.N / P.M
    if d <= 0.0:
        raise ValueError("eta >= M/N: no yield surface through this stress")
    return p * ((1.0 - P.N) / d) ** ((1.0 - P.N) / P.N)


@dataclass
class Yield:
    eta: float
    F_p: float
    F_pp: float
    F_pi: float
    F_ppi: float
    eta_p: float
    eta_pi: float


def yield_p(P: Params, p: float, pi: float) -> Yield:
    if p >= 0.0 or pi >= PI_MAX_NEG:
        raise EvalError("p_or_pi_nonneg")
    N = P.N
    eta = eta_of(P, p, pi)
    r = p / pi
    F_p = (eta - P.M) / (1.0 - N)                       # = M ln(pi/p) for N = 0
    F_pp = -(P.M / (1.0 - N)) / p * r ** (N / (1.0 - N))
    F_pi = P.M * r ** (1.0 / (1.0 - N))
    F_ppi = (P.M / ((1.0 - N) * p)) * r ** (1.0 / (1.0 - N))
    return Yield(eta, F_p, F_pp, F_pi, F_ppi, (F_p - eta) / p, F_pi / p)


# --------------------------------------------------------------------------------------
# §6 CSL: psi_i(v, pi_i) and Lambda(pi_i) (S.22)
# --------------------------------------------------------------------------------------
def csl(P: Params, v: float, pi: float):
    if P.csl_mode == "paper":
        return v - P.v_c0 + P.lam_tilde * math.log(-pi), P.lam_tilde
    t = (-pi / P.p_a) ** P.xi
    return (v - 1.0) - P.e0 + P.lam_c * t, P.lam_c * P.xi * t


# --------------------------------------------------------------------------------------
# §7 pi_i* = Pi(p, Omega, psi_i) and its derivatives (S.23)-(S.24)
# --------------------------------------------------------------------------------------
def pistar(P: Params, p: float, Om: float, psi: float):
    """returns pi*, Pi_psi, Pi_Omega. Raises EvalError('B_nonpos') on the §7 guard."""
    cb = P.chi_bar
    a = SQ23 * cb
    if P.N == 0.0:
        ps = p * math.exp(a * psi * Om / P.M)
        return ps, ps * a * Om / P.M, ps * a * psi / P.M
    B = 1.0 - a * psi * Om * P.N / P.M
    if B <= 0.0:
        raise EvalError("B_nonpos")
    ps = p * B ** ((P.N - 1.0) / P.N)
    den = P.M - a * psi * Om * P.N
    return ps, a * Om * (1.0 - P.N) * ps / den, a * psi * (1.0 - P.N) * ps / den


# --------------------------------------------------------------------------------------
# §5.3 + §10 flow vector, Hessian, Omega (S.17)-(S.21), cap (S.35)-(S.36)
# --------------------------------------------------------------------------------------
@dataclass
class Flow:
    F: float
    f_a: np.ndarray
    q_a: np.ndarray
    q_ab: np.ndarray
    q_api: np.ndarray
    Om: float
    Om_a: np.ndarray
    Om_pi: float
    w: float
    zeta: float
    zeta_bar: float
    Y: Yield


def cap_weight(P: Params, eta: float):
    """w(eta), dw/deta of (S.35). cap='none': w = 1. planar: step at eta_1 (w = 1 for eta >= eta_1,
    i.e. the cap applies for eta < chi_cap M as in BA06 2.76; choice at equality: w = 1)."""
    if P.cap == "none":
        return 1.0, 0.0
    e1 = P.c1 * P.M
    e2 = P.c2 * P.M
    if P.cap == "planar":
        return (1.0 if eta >= e1 else 0.0), 0.0
    t = (eta - e1) / (e2 - e1)
    if t <= 0.0:
        return 0.0, 0.0
    if t >= 1.0:
        return 1.0, 0.0
    S = t ** 3 * (10.0 - 15.0 * t + 6.0 * t * t)
    Sp = 30.0 * t * t * (1.0 - t) ** 2
    return S, Sp / (e2 - e1)


def omega_only(P: Params, inv: Invariants, pi: float):
    """Omega (capped), Omega_pi, eta -- the pieces the nested pi_i loop needs (cheap path)."""
    Y = yield_p(P, inv.p, pi)
    if inv.vertex:
        return 0.0, 0.0, Y
    zb, zby, _ = zeta_y(inv.theta, P.rho_bar, P.zeta)
    Sy = float(np.dot(inv.y_a, inv.y_a))
    Omu = math.sqrt(1.5 * zb * zb + (zby * inv.q) ** 2 * Sy)
    w, w_eta = cap_weight(P, Y.eta)
    return w * Omu, w_eta * Y.eta_pi * Omu, Y


def flow(P: Params, inv: Invariants, pi: float) -> Flow:
    Y = yield_p(P, inv.p, pi)
    beta = P.beta
    if inv.vertex:
        # §3.2: purely volumetric flow; F = p eta; Omega = 0
        F = inv.p * Y.eta
        f_a = Y.F_p / 3.0 * ONES
        qu_a = beta * Y.F_p / 3.0 * ONES
        qu_ab = beta * Y.F_pp / 9.0 * np.outer(ONES, ONES)
        qu_api = beta * Y.F_ppi / 3.0 * ONES
        Omu, Omu_a = 0.0, np.zeros(3)
        z = zb = 1.0
    else:
        z, zy, _ = zeta_y(inv.theta, P.rho, P.zeta)
        zb, zby, zbyy = zeta_y(inv.theta, P.rho_bar, P.zeta)
        q, nh, y_a, nh_ab, y_ab = inv.q, inv.nhat, inv.y_a, inv.nhat_ab, inv.y_ab
        F = z * q + inv.p * Y.eta
        f_a = Y.F_p / 3.0 * ONES + SQ32 * z * nh + zy * q * y_a                       # (S.14)
        qu_a = beta * Y.F_p / 3.0 * ONES + SQ32 * zb * nh + zby * q * y_a             # (S.17)
        qu_ab = (beta * Y.F_pp / 9.0 * np.outer(ONES, ONES) + SQ32 * zb * nh_ab
                 + zby * q * y_ab + zbyy * q * np.outer(y_a, y_a)
                 + SQ32 * zby * (np.outer(nh, y_a) + np.outer(y_a, nh)))              # (S.18)
        qu_api = beta * Y.F_ppi / 3.0 * ONES                                          # (S.19)
        Sy = float(np.dot(y_a, y_a))
        Omu = math.sqrt(1.5 * zb * zb + (zby * q) ** 2 * Sy)                          # (S.20)
        Omu_a = (1.5 * zb * zby * y_a
                 + zby * q * (zbyy * q * y_a + SQ32 * zby * nh) * Sy
                 + (zby * q) ** 2 * (y_ab @ y_a)) / Omu                               # (S.21)
    w, w_eta = cap_weight(P, Y.eta)
    if w == 1.0 and w_eta == 0.0:
        return Flow(F, f_a, qu_a, qu_ab, qu_api, Omu, Omu_a, 0.0, w, z, zb, Y)
    g_a = qu_a + ONES / 3.0
    q_a = -ONES / 3.0 + w * g_a
    q_ab = w * qu_ab + w_eta * Y.eta_p / 3.0 * np.outer(g_a, ONES)                    # (S.36)
    q_api = w * qu_api + w_eta * Y.eta_pi * g_a
    Om = w * Omu
    Om_a = w * Omu_a + w_eta * Y.eta_p / 3.0 * Omu * ONES
    Om_pi = w_eta * Y.eta_pi * Omu
    return Flow(F, f_a, q_a, q_ab, q_api, Om, Om_a, Om_pi, w, z, zb, Y)


# --------------------------------------------------------------------------------------
# §8 nested scalar Newton for pi_i (S.26)-(S.27), (S.37)
# --------------------------------------------------------------------------------------
def _pi_residual(P: Params, inv: Invariants, dlam: float, v: float, pi_n: float, pi: float):
    if pi >= PI_MAX_NEG:
        raise EvalError("pi_nonneg")
    k = SQ23 * P.h
    Om, Om_pi, Y = omega_only(P, inv, pi)
    psi, Lam = csl(P, v, pi)
    ps, Ppsi, POm = pistar(P, inv.p, Om, psi)
    r = pi - pi_n - k * dlam * (ps - pi) * Om
    rp = 1.0 - k * dlam * ((POm * Om_pi + Ppsi * Lam / pi - 1.0) * Om + (ps - pi) * Om_pi)  # (S.37)
    return r, rp


def solve_pi(P: Params, inv: Invariants, dlam: float, v: float, pi_n: float):
    """Nested scalar solve of r(pi_i) = 0 (S.27)/(S.37): returns (pi_i, c = r'(pi_i), iters).

    CONTRACT: the root CONTINUOUS WITH pi_i,n (the first root of r along the direction of
    decrease of |r| from pi_i,n), never a farther one.
    1. Evaluate r at pi_i,n; |r| <= tol (PI_TOL_REL*|pi_i,n|) -> pi_i = pi_i,n, 0 iterations
       (always the case for dlam = 0, and at the vertex where Omega = 0).
    2. Scan from pi_i,n in fixed steps of PI_SCAN_REL*|pi_i,n|: AWAY from zero (more negative)
       when r(pi_i,n) > 0 (hardening: pi_i* more negative than pi_i,n; r is increasing in pi_i),
       TOWARD zero when r < 0. Stop at the first sign change of r: that is the bracket. If |r|
       GROWS between two scan points before a sign change, r has passed a local extremum: the
       root continuous with pi_i,n has annihilated in a fold (cap ramp, sheet §16.3, loop gain
       >= 1) and the only remaining root is a far one -> EvalError('pi_fold'); the caller treats
       it as a rejected iterate (line-search backtrack of the step, then substepping).
       PI_SCAN_MAX steps without a sign change -> EvalError('pi_nobracket'); a scan point
       >= PI_MAX_NEG -> EvalError('pi_nonneg').
       Without the cap r is monotone in pi_i (r' = 1 + sqrt(2/3) h dlam Omega [1 - Pi_psi Lambda/pi_i] > 0
       for every admissible state), so the scan can only end in a bracket or 'pi_nobracket'.
    3. Safeguarded Newton ("rtsafe") inside that first bracket: start from the end with the smaller
       |r|; Newton step if it lands strictly inside the bracket, bisection otherwise; the bracket is
       updated by sign every iteration; converged when |r| <= tol; MAX_PI_ITERS -> EvalError('pi_noconv').
    The fixed small scan step is what makes the selection unambiguous: a factor-2 bracket from
    pi_i,n (the earlier contract) can enclose the far root of the fold while the near one exists,
    and the inner iteration then converges to the wrong (w = 1) branch.
    An undamped Newton is what AB06 61-62 describe; with the cap on, Omega depends on pi_i
    through w(eta) and the undamped iteration overshoots where w ramps (sheet §16.3), which
    the safeguard removes without changing the root."""
    tol = PI_TOL_REL * abs(pi_n)
    a = pi_n
    ra, rpa = _pi_residual(P, inv, dlam, v, pi_n, a)
    if abs(ra) <= tol:
        return a, rpa, 0
    # 2. scan for the first bracket
    d = -1.0 if ra > 0.0 else 1.0                 # pi_i < 0: -1 = away from zero
    hstep = PI_SCAN_REL * abs(pi_n)
    it = 0
    for _ in range(PI_SCAN_MAX):
        b = a + d * hstep
        rb, rpb = _pi_residual(P, inv, dlam, v, pi_n, b)
        it += 1
        if abs(rb) <= tol:
            return b, rpb, it
        if rb * ra < 0.0:
            break
        if abs(rb) > abs(ra):
            raise EvalError("pi_fold")
        a, ra, rpa = b, rb, rpb
    else:
        raise EvalError("pi_nobracket")
    # 3. safeguarded Newton inside [lo, hi]
    lo, hi = (a, b) if a < b else (b, a)          # lo more negative
    rlo = ra if a < b else rb
    if abs(ra) < abs(rb):
        x, r, rp = a, ra, rpa
    else:
        x, r, rp = b, rb, rpb
    for _ in range(MAX_PI_ITERS):
        xn = x - r / rp if rp != 0.0 else lo
        if not (lo < xn < hi):
            xn = 0.5 * (lo + hi)
        r, rp = _pi_residual(P, inv, dlam, v, pi_n, xn)
        x = xn
        it += 1
        if abs(r) <= tol:
            return x, rp, it
        if r * rlo > 0.0:
            lo, rlo = x, r
        else:
            hi = x
    raise EvalError("pi_noconv")


# --------------------------------------------------------------------------------------
# §9 residual, Jacobian, tangent
# --------------------------------------------------------------------------------------
@dataclass
class PointEval:
    """Everything at one local iterate x = (eps_e, dlam)."""
    eps_e: np.ndarray
    dlam: float
    el: Elastic
    inv: Invariants
    fl: Flow
    pi: float
    c: float
    pi_iters: int
    psi: float
    Lam: float
    ps: float
    Ppsi: float
    POm: float
    r: np.ndarray
    # sensitivities of the converged pi_i (S.28)
    P_a: np.ndarray = field(default=None)
    Pi_b: np.ndarray = field(default=None)
    Pi_lam: float = 0.0
    Pi_v: float = 0.0


def evaluate(P: Params, eps_e: np.ndarray, dlam: float, eps_tr: np.ndarray, v: float,
             pi_n: float) -> PointEval:
    el = elastic(P, eps_e)
    inv = invariants(el.sig)
    pi, c, its = solve_pi(P, inv, dlam, v, pi_n)
    fl = flow(P, inv, pi)
    psi, Lam = csl(P, v, pi)
    ps, Ppsi, POm = pistar(P, inv.p, fl.Om, psi)
    r = np.empty(4)
    r[:3] = eps_e - eps_tr + dlam * fl.q_a
    r[3] = fl.F
    pe = PointEval(eps_e.copy(), dlam, el, inv, fl, pi, c, its, psi, Lam, ps, Ppsi, POm, r)
    k = SQ23 * P.h
    pistar_a = (ps / (3.0 * inv.p)) * ONES + POm * fl.Om_a                      # (S.24)
    pe.P_a = (k * dlam / c) * (fl.Om * pistar_a + (ps - pi) * fl.Om_a)          # (S.28)/(S.37)
    pe.Pi_b = pe.P_a @ el.ae
    pe.Pi_lam = (k * fl.Om / c) * (ps - pi)
    pe.Pi_v = (k * dlam * fl.Om / c) * Ppsi
    return pe


def scaled_norm(P: Params, r: np.ndarray) -> float:
    return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2 + (r[3] / abs(P.p0)) ** 2)


def jacobian(P: Params, pe: PointEval) -> np.ndarray:
    """(S.30)."""
    fl, el = pe.fl, pe.el
    J = np.empty((4, 4))
    J[:3, :3] = I3 + pe.dlam * (fl.q_ab @ el.ae + np.outer(fl.q_api, pe.Pi_b))
    J[:3, 3] = fl.q_a + pe.dlam * fl.q_api * pe.Pi_lam
    J[3, :3] = fl.f_a @ el.ae + fl.Y.F_pi * pe.Pi_b
    J[3, 3] = fl.Y.F_pi * pe.Pi_lam
    return J


def atilde_ep(P: Params, pe: PointEval, J: np.ndarray, vfac: float) -> np.ndarray:
    """(S.31)-(S.32): a~^ep_ab = d sigma_a / d eps~_b. vfac = d v_{n+1}/d eps~_b = v_{n+1}, the converged
    specific volume of the step, in both modes (sheet §1.2, G2 owner decision 2026-10-01; was v0 in small
    strain under the linear update). The caller must pass vfac equal to the v it passed to return_map."""
    fl = pe.fl
    s = np.empty(4)
    s[:3] = pe.dlam * fl.q_api * pe.Pi_v * vfac
    s[3] = fl.Y.F_pi * pe.Pi_v * vfac
    b = np.linalg.inv(J)
    dxde = b[:, :3] - np.outer(b @ s, ONES)
    return pe.el.ae @ dxde[:3, :]


@dataclass
class ChainData:
    """Per-sub-increment sensitivities of a converged PLASTIC return map, sheet §9.6 (S.45):
    the implicit-function-theorem columns of the map (eps~, pi_i,n, v) -> (x, pi_i,n+1).
    J = A + t Pi_x^T with t = dr/dpi_i|_x = (dlam q_api, F_pi), Pi_x = (Pi_b, Pi_lam) (S.28);
    b = J^-1, u = b t, kappa = Pi_x . u, c = r'(pi_i) (S.27), w_b = sum_c Pi_c b_cb + Pi_lam b_4b.
      dx/deps~_b = b[:, b]        dpi_i,n+1/deps~_b = w_b
      dx/dpi_i,n = -u/c           dpi_i,n+1/dpi_i,n  = (1 - kappa)/c
      dx/dv      = -u Pi_v        dpi_i,n+1/dv       = (1 - kappa) Pi_v
    Vertex branch (§3.2): Pi_x = 0, Pi_v = 0, c = 1, kappa = 0 come out of the same formulas."""
    b: np.ndarray              # (4,4) = J^-1
    u: np.ndarray              # (4,)  = b t
    w: np.ndarray              # (3,)  = Pi_x^T b[:, :3]
    kappa: float
    c: float
    Pi_v: float


@dataclass
class StepResult:
    eps_e: np.ndarray          # converged principal elastic strains
    sig: np.ndarray            # converged principal stresses
    pi: float
    dlam: float
    q_a: np.ndarray            # flow direction (capped) at the converged state
    Om: float
    F_p: float
    eta: float
    psi: float
    ps: float                  # pi_i*
    D: float                   # dissipation of the step, dlam * sum sigma_a q_a  (S.38)
    atilde: np.ndarray         # a~^ep (3,3)
    plastic: bool
    vertex: bool
    cap_active: bool
    refused: bool
    reason: str
    local_iters: int
    pi_iters: int
    res_hist: list
    ae: np.ndarray = field(default=None)       # a^e_ab (S.3) at the CONVERGED eps_e (= atilde when elastic)
    chain: ChainData = field(default=None)     # §9.6 sensitivities (plastic, accepted steps only)
    eps_tr: np.ndarray = field(default=None)   # trial principal strains eps~_a of this step


def chain_data(pe: PointEval, J: np.ndarray) -> ChainData:
    """(S.45) block data from the converged iterate and its Jacobian (S.30)."""
    fl = pe.fl
    t = np.empty(4)
    t[:3] = pe.dlam * fl.q_api
    t[3] = fl.Y.F_pi
    Pi_x = np.append(pe.Pi_b, pe.Pi_lam)
    b = np.linalg.inv(J)
    u = b @ t
    kappa = float(Pi_x @ u)
    w = Pi_x @ b[:, :3]
    return ChainData(b, u, w, kappa, pe.c, pe.Pi_v)


def return_map(P: Params, eps_tr: np.ndarray, pi_n: float, v: float, vfac: float,
               force_plastic: bool = False) -> StepResult:
    """One backward-Euler step in principal space (sheet §9.1). Never raises for
    numerical trouble: refusals come back in StepResult.refused/reason."""
    # 1-2. trial
    el = elastic(P, eps_tr)
    inv = invariants(el.sig)
    try:
        fl0 = flow(P, inv, pi_n)
        psi0, _ = csl(P, v, pi_n)
    except EvalError as e:
        return StepResult(eps_tr, el.sig, pi_n, 0.0, np.zeros(3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          el.ae, False, inv.vertex, False, True, "trial_" + str(e), 0, 0, [])
    if fl0.F <= F_TRIAL_TOL_REL * abs(P.p0) and not force_plastic:
        return StepResult(eps_tr.copy(), el.sig, pi_n, 0.0, np.zeros(3), 0.0, fl0.Y.F_p, fl0.Y.eta,
                          psi0, 0.0, 0.0, el.ae, False, inv.vertex, fl0.w < 1.0, False, "", 0, 0, [fl0.F],
                          ae=el.ae, chain=None, eps_tr=eps_tr.copy())
    # 3-4. local Newton on x = (eps_e, dlam)
    x = np.append(eps_tr, 0.0)
    hist = []
    pi_total = 0
    try:
        pe = evaluate(P, x[:3], x[3], eps_tr, v, pi_n)
    except EvalError as e:
        return StepResult(eps_tr, el.sig, pi_n, 0.0, np.zeros(3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          el.ae, True, inv.vertex, False, True, "local_" + str(e), 0, 0, [])
    rn = scaled_norm(P, pe.r)
    hist.append(rn)
    pi_total += pe.pi_iters
    it = 0
    converged = rn <= RES_TOL
    while not converged and it < MAX_LOCAL_ITERS:
        J = jacobian(P, pe)
        try:
            dx = np.linalg.solve(J, pe.r)
        except np.linalg.LinAlgError:
            return StepResult(eps_tr, el.sig, pi_n, 0.0, np.zeros(3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              el.ae, True, inv.vertex, False, True, "local_singular_J", it, pi_total, hist)
        alpha = 1.0
        accepted = False
        last_err = ""
        for _ in range(MAX_LINESEARCH + 1):
            xn = x - alpha * dx
            try:
                pen = evaluate(P, xn[:3], xn[3], eps_tr, v, pi_n)
                rnn = scaled_norm(P, pen.r)
                if rnn < rn or rnn <= RES_TOL:
                    accepted = True
                    break
            except EvalError as e:
                # a nested failure (pi_fold / pi_nobracket / pi_noconv / pi_nonneg) or an
                # inadmissible iterate (p >= 0, B <= 0): the iterate is rejected and the
                # step (dlam included) is backtracked exactly like a non-decreasing residual.
                last_err = str(e)
            alpha *= 0.5
        if not accepted:
            reason = "local_linesearch" + (":" + last_err if last_err else "")
            return StepResult(eps_tr, el.sig, pi_n, 0.0, np.zeros(3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              el.ae, True, inv.vertex, False, True, reason, it, pi_total, hist)
        x, pe, rn = xn, pen, rnn
        pi_total += pe.pi_iters
        hist.append(rn)
        it += 1
        converged = rn <= RES_TOL
    if not converged:
        return StepResult(eps_tr, el.sig, pi_n, 0.0, np.zeros(3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          el.ae, True, inv.vertex, False, True, "local_noconv", it, pi_total, hist)
    if pe.dlam < 0.0:
        return StepResult(eps_tr, el.sig, pi_n, pe.dlam, pe.fl.q_a, pe.fl.Om, 0.0, 0.0, 0.0, 0.0, 0.0,
                          el.ae, True, inv.vertex, False, True, "negative_dlambda", it, pi_total, hist)
    # 5. tangent and diagnostics
    J = jacobian(P, pe)
    at = atilde_ep(P, pe, J, vfac)
    D = pe.dlam * float(np.dot(pe.el.sig, pe.fl.q_a))
    return StepResult(pe.eps_e, pe.el.sig, pe.pi, pe.dlam, pe.fl.q_a, pe.fl.Om, pe.fl.Y.F_p, pe.fl.Y.eta,
                      pe.psi, pe.ps, D, at, True, pe.inv.vertex, pe.fl.w < 1.0, False, "", it, pi_total, hist,
                      ae=pe.el.ae, chain=chain_data(pe, J), eps_tr=eps_tr.copy())


# --------------------------------------------------------------------------------------
# §9.4 / §9.5 spectral assembly of the 4th-order tangents
# --------------------------------------------------------------------------------------
def _spectral(nvec: np.ndarray, diag_ab: np.ndarray, spin_ab: np.ndarray, half: float) -> np.ndarray:
    """C = sum_ab diag_ab m^a (x) m^b + half * sum_{a!=b} spin_ab (m^ab (x) m^ab + m^ab (x) m^ba)."""
    m = [np.outer(nvec[:, a], nvec[:, a]) for a in range(3)]
    C = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            C += diag_ab[a, b] * np.einsum("ij,kl->ijkl", m[a], m[b])
            if a != b:
                mab = np.outer(nvec[:, a], nvec[:, b])
                mba = mab.T
                C += half * spin_ab[a, b] * (np.einsum("ij,kl->ijkl", mab, mab)
                                             + np.einsum("ij,kl->ijkl", mab, mba))
    return C


def tangent_small(atilde: np.ndarray, sig: np.ndarray, eps_tr: np.ndarray, nvec: np.ndarray) -> np.ndarray:
    """(S.33): the small-strain consistent tangent from a~^ep, converged sigma_a, trial eps~_a."""
    g = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            if a == b:
                continue
            d = eps_tr[a] - eps_tr[b]
            if abs(d) < REPEATED_EIG_TOL:
                g[a, b] = atilde[a, a] - atilde[a, b]
            else:
                g[a, b] = (sig[a] - sig[b]) / d
    return _spectral(nvec, atilde, g, 0.5)


# --------------------------------------------------------------------------------------
# §9.6 chained consistent tangent across sub-increments (S.45)-(S.47)
# --------------------------------------------------------------------------------------
# Column convention (sheet §9.6 "Recursion"): the six d_eps slots J = {00, 11, 22, 01, 12, 02} with
# input tensors E_J = e_k e_k (normal) and e_k e_l + e_l e_k (shear: the independent TENSOR shear
# component, eps_kl and eps_lk moved together; tr E_J = 1 or 0).
CHAIN_IJ = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))
CHAIN_E = []
for _k, _l in CHAIN_IJ:
    _E = np.zeros((3, 3))
    _E[_k, _l] += 1.0
    _E[_l, _k] += 1.0 if _k != _l else 0.0
    CHAIN_E.append(_E)
CHAIN_E = np.array(CHAIN_E)                          # (6,3,3)
CHAIN_TRE = np.array([float(np.trace(E)) for E in CHAIN_E])


def chain_start():
    """S^eps_0 = 0 (6 columns of 3x3), S^pi_0 = 0, cumulative fraction 0 (S.46 initial state)."""
    return np.zeros((6, 3, 3)), np.zeros(6), 0.0


def chain_propagate(S_eps: np.ndarray, S_pi: np.ndarray, cum_before: float, alpha: float, v_new: float,
                    res: StepResult, nvec: np.ndarray):
    """One sub-increment of (S.46). Inputs: the sensitivities of the state ENTERING the
    sub-increment (S^eps_k, S^pi_k, sum_{j<k} alpha_j), its fraction alpha_k, v_new = v_{k+1} (the
    sub-increment's own CONVERGED specific volume, v_k exp(alpha_k tr d_eps); sheet §1.2/§9.6, G2), the
    accepted StepResult of the sub-increment (eps~_a, eps^e_a, b, u, w, kappa, c, Pi_v) and the trial
    eigenvectors. Returns (S^eps_{k+1}, S^pi_{k+1}, sum_{j<=k} alpha_j).

      T_k        = S^eps_k + alpha_k E_J                         (d eps~_k / d d_eps_J)
      S^v_{k+1}  = v_{k+1} (sum_{j<=k} alpha_j) tr E_J            (closed form: v_{k+1} = v_n exp((sum_{j<=k} alpha_j) tr d_eps)
                                                                  differentiated; was v0 (...) under the linear update)
      plastic:   S^eps_{k+1} = Phi : T_k - sum_a m^a [ (u_a/c) S^pi_k + u_a Pi_v S^v_{k+1} ]
                 S^pi_{k+1}  = sum_b w_b T^_bb + ((1 - kappa)/c) S^pi_k + (1 - kappa) Pi_v S^v_{k+1}
      elastic:   S^eps_{k+1} = T_k,  S^pi_{k+1} = S^pi_k
    Phi = d eps^e_{k+1} / d eps~ in the (S.33) form: diagonal block b_ab (a, b <= 3), spin
    (eps^e_a - eps^e_b)/(eps~_a - eps~_b), limit b_aa - b_ab (tangent_small with a~ -> b, sigma -> eps^e).
    The full 3x3 column tensors are kept (the (S.33) row convention per operator; each operator
    symmetrises its own input), sheet §9.6 "Contract"."""
    cum = cum_before + alpha
    T = S_eps + alpha * CHAIN_E                                   # (6,3,3)
    if not res.plastic:
        return T, S_pi.copy(), cum
    ch = res.chain
    Phi = tangent_small(ch.b[:3, :3], res.eps_e, res.eps_tr, nvec)
    S_v = v_new * cum * CHAIN_TRE                                 # (6,)  S^v_{k+1} = v_{k+1} cum tr E_J
    m = np.array([np.outer(nvec[:, a], nvec[:, a]) for a in range(3)])   # (3,3,3)
    # per column: Phi : T_J, then the pi_i,n and v columns (eigenvalues only)
    PhiT = np.einsum("ijkl,Jkl->Jij", Phi, T)
    coef = np.outer(S_pi, ch.u[:3] / ch.c) + np.outer(S_v, ch.u[:3] * ch.Pi_v)     # (6,3): per column, per a
    S_eps_new = PhiT - np.einsum("Ja,aij->Jij", coef, m)
    That = np.einsum("ia,Jij,jb->Jab", nvec, T, nvec)             # T in the trial basis
    Tdiag = np.einsum("Jaa->Ja", That)                            # T^_bb
    S_pi_new = Tdiag @ ch.w + ((1.0 - ch.kappa) / ch.c) * S_pi + (1.0 - ch.kappa) * ch.Pi_v * S_v
    return S_eps_new, S_pi_new, cum


def chain_assemble(P: Params, res: StepResult, nvec: np.ndarray, S_eps: np.ndarray) -> np.ndarray:
    """(S.47): C = a^e(eps^e_m) : S^eps_m with a^e in the (S.33) form on the eigen-data of the FINAL
    converged eps^e_m (spin (sigma_a - sigma_b)/(eps^e_a - eps^e_b), limit a^e_aa - a^e_ab).
    Returns the 3x3x3x3 tensor with C4[:, :, k, l] = C4[:, :, l, k] = column_J / 2 on the shear
    slots (so C4 : E = column for every symmetric E, and the parity c4_to_c6 / c6_of reduction
    C4_ijkl + C4_ijlk recovers the kernel's 6x6 column exactly) and column_J on the normal slots."""
    Ae = tangent_small(res.ae, res.sig, res.eps_e, nvec)
    cols = np.einsum("ijkl,Jkl->Jij", Ae, S_eps)                  # (6,3,3)
    C = np.zeros((3, 3, 3, 3))
    for J, (k, l) in enumerate(CHAIN_IJ):
        if k == l:
            C[:, :, k, k] = cols[J]
        else:
            C[:, :, k, l] = 0.5 * cols[J]
            C[:, :, l, k] = 0.5 * cols[J]
    return C


def tangent_finite(atilde: np.ndarray, tau: np.ndarray, eps_tr: np.ndarray, nvec: np.ndarray) -> np.ndarray:
    """(S.34): a^ep = c~ + tau (+) 1 with (tau (+) 1)_ijkl = tau_jl delta_ik. eps_tr are the trial
    principal LOG stretches; atilde from (S.32) (vfac = v in both modes since G2, nothing to substitute)."""
    lam = np.exp(eps_tr)
    c = atilde - 2.0 * np.diag(tau)
    gam = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            if a == b:
                continue
            if abs(lam[a] - lam[b]) < REPEATED_STRETCH_TOL:
                gam[a, b] = 0.5 * (atilde[b, b] - atilde[b, a]) - tau[a]
            else:
                gam[a, b] = (tau[b] * lam[a] ** 2 - tau[a] * lam[b] ** 2) / (lam[b] ** 2 - lam[a] ** 2)
    ct = _spectral(nvec, c, gam, 1.0)
    T = sum(tau[a] * np.outer(nvec[:, a], nvec[:, a]) for a in range(3))
    # (tau (+) 1)_ijkl = tau_jl delta_ik
    ct += np.einsum("jl,ik->ijkl", T, I3)
    return ct
