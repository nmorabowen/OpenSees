"""Shared plumbing of the WP-144 round-3b G1 files (Zone B; NOT a test file):
    test_g1_har.py       the HAR energy option      (sheet 144a 2.3, 2.4, 13.1h, 13.11)
    test_g1_floor.py     the p' floor Pi_f          (sheet 9.7, 13.12 - 13.14b)
    test_g1_pi0_scan.py  the unified pi_i0 rule and the (S.56) refusal   (sheet 5.4, 10.2, 13.15)

TEST RULE (plan 5, restated): every expected value in those files is written from the equation sheet (a closed form
coded HERE from the sheet, or a number printed in the sheet) BEFORE the oracle runs on the path.  Nothing below calls
an oracle for an expected value; the oracles are imported only to be driven (`make_hf`) and to be asked for their
OWN operator in the few kernel-level tests that say so.

What is here
  * the TIMs HAR constants of the sheet (2.3 gate table / 13.1h): n = 1/2, g = 807.80387674, k = 1889.48104361, p_a = 101;
  * the sheet's closed forms: (S.5h) p, q; (S.5h') D11, D12, D22, q/eps_s; (S.5h'') inverse map; (S.49) BA06 floor;
    (S.50) HAR floor (n = 1/2 closed, general n by plain bisection on the sheet's bracket); (S.12) eta; (S.53) pi_i0;
    (S.56) W_ramp; (S.22) the fork-CSL specific volume for a target psi_i;
  * `make_hf(oracle, **kw)`: builds an oracle Params from the sheet's common names, mapping p_min = 'default' to the
    oracle's own spelling (O2: None; O1: 'default'), an explicit number as is;
  * small tensor helpers (principal stress from (p, q, n_hat), Mandel 6 x 6 of a 4th-order tensor).
"""
from __future__ import annotations

import math

import numpy as np

from conftest import make_params

SQ23 = math.sqrt(2.0 / 3.0)
SQ32 = math.sqrt(1.5)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)
ONES = np.ones(3)

# ------------------------------------------------------------------------------------------------------------------
# the TIMs HAR set (sheet 2.3 "Gate table", 13.1h): n = 1/2, G0 = 264.32, nu = 0.3129, e_ref = 0.6944, p_a = 101 kPa
# ------------------------------------------------------------------------------------------------------------------
K_T, G_T, N_T, PA_T = 1889.48104361, 807.80387674, 0.5, 101.0
PMIN_T = 5.0e-3 * PA_T                       # 0.505 kPa, the parser default (sheet 1.3, 9.7)
EDGE_T = 1.0 / (K_T * (1.0 - N_T))           # 1.058491699e-3, p = 0 (sheet 13.1h)

HAR_ELASTIC = dict(energy="HAR", k=K_T, g=G_T, n_e=N_T, p_a=PA_T)

# the plastic / CSL constants of K1.14b (sheet 13.14b): TIMs M, rho = rho_bar = 0.71 (WW), fork CSL at the TIMs values,
# the K2 plastic constants N 0.4, N_bar 0.2, chi -3.5, h 280.  Condition A (S.39): N_bar <= N, rho/rho_bar = 1 >= 0.75.
PLASTIC_T = dict(M=1.3309, N=0.4, N_bar=0.2, chi=-3.5, h=280.0, rho=0.71, rho_bar=0.71, zeta="WW",
                 csl_mode="fork", e0=0.83, lambda_c=0.027, xi=0.45, cap="none")
# K2 (paper) plastic constants, sheet 14, with the BA06 energy
K2_PLASTIC = dict(M=1.2, N=0.4, N_bar=0.2, chi=-3.5, h=280.0, rho=0.7, rho_bar=0.8, zeta="WW",
                  csl_mode="paper", lambda_tilde=0.0135, v_c0=1.81, cap="none")
K2_BA06 = dict(p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0)


def har_kw(**over):
    """Common-name kwargs: HAR TIMs elastic set + the K1.14b plastic/CSL constants (+ overrides)."""
    kw = dict(HAR_ELASTIC)
    kw.update(PLASTIC_T)
    kw.update(over)
    return kw


def ba06_kw(**over):
    """Common-name kwargs: BA06 K2 set (paper CSL) (+ overrides)."""
    kw = dict(K2_BA06)
    kw.update(K2_PLASTIC)
    kw.update(over)
    return kw


def make_hf(oracle: str, **kw):
    """oracle Params from common names.  p_min: 'default' -> that oracle's spelling of the parser default
    (5e-3 p_ref); a number is passed as is; absent -> NOT passed (oracle default: O2 5e-3 p_ref, O1 0)."""
    if "p_min" in kw and kw["p_min"] == "default":
        if oracle == "O2":
            kw = dict(kw)
            del kw["p_min"]
        else:
            kw = dict(kw, p_min="default")
    return make_params(oracle, **kw)


# ------------------------------------------------------------------------------------------------------------------
# HAR closed forms, written from the sheet (S.4h)-(S.5h'') (HAR05 eq 40-46, sign-mapped)
# ------------------------------------------------------------------------------------------------------------------
def har_pq(ev, es, k=K_T, g=G_T, n=N_T, pa=PA_T):
    """(S.5h): p = -p_a k(1-n) eps* w, q = 3 g p_a eps_s w, w = [k(1-n) u]^(n/(1-n)),
    u = [eps*^2 + 3 g eps_s^2/(k(1-n))]^(1/2), eps* = 1/(k(1-n)) - eps_v.  Raises ValueError for eps* <= 0."""
    kn = k * (1.0 - n)
    est = 1.0 / kn - ev
    if est <= 0.0:
        raise ValueError("eps* <= 0: outside the HAR domain")
    u = math.sqrt(est * est + 3.0 * g * es * es / kn)
    w = (kn * u) ** (n / (1.0 - n))
    return -pa * kn * est * w, 3.0 * g * pa * es * w


def har_D(p, q, k=K_T, g=G_T, n=N_T, pa=PA_T):
    """(S.5h'): D11, D12, D22 and q/eps_s from the STRESS form, varpi^2 = p^2 + k(1-n) q^2/(3g)."""
    kn = k * (1.0 - n)
    vp2 = p * p + kn * q * q / (3.0 * g)
    fac = pa * (math.sqrt(vp2) / pa) ** n
    Z = vp2 / (p * p)
    D11 = k * fac * (1.0 - n + n / Z)
    D22 = (3.0 * g / (1.0 - n)) * fac * (1.0 - n / Z)
    D12 = n * k * p * q * fac / vp2
    return D11, D12, D22, 3.0 * g * fac


def har_inverse(p, q, k=K_T, g=G_T, n=N_T, pa=PA_T):
    """(S.5h''): eps_v^e, eps_s^e, varpi from (p < 0, q >= 0)."""
    kn = k * (1.0 - n)
    vp = math.sqrt(p * p + kn * q * q / (3.0 * g))
    ev = (1.0 / kn) * (1.0 - (abs(p) / pa) ** (1.0 - n) * (abs(p) / vp) ** n)
    es = q / (3.0 * g * pa * (vp / pa) ** n)
    return ev, es, vp


def har_Kiso(p, k=K_T, n=N_T, pa=PA_T):
    """K_iso(p) = k p_a (|p|/p_a)^n (sheet 2.3 gate table normalisation)."""
    return k * pa * (abs(p) / pa) ** n


def har_psi(ev, es, k=K_T, g=G_T, n=N_T, pa=PA_T):
    """(S.4h) Psi."""
    kn = k * (1.0 - n)
    est = 1.0 / kn - ev
    u = math.sqrt(est * est + 3.0 * g * es * es / kn)
    return pa / (k * (2.0 - n)) * (kn * u) ** ((2.0 - n) / (1.0 - n))


# ------------------------------------------------------------------------------------------------------------------
# floor closed forms: (S.49) BA06, (S.50) HAR
# ------------------------------------------------------------------------------------------------------------------
def ba06_floor(es, pmin, p0=-100.0, kh=0.01, ev0=0.0, a0=0.0):
    """(S.49): eps_v,f = eps_v0 - kh ln[p_min / (|p0| (1 + 3 a0 es^2/(2 kh)))], eps'_f = 3 a0 es/(1 + 3 a0 es^2/(2 kh))."""
    fac = 1.0 + 3.0 * a0 * es * es / (2.0 * kh)
    return dict(ev_f=ev0 - kh * math.log(pmin / (abs(p0) * fac)), epsp=3.0 * a0 * es / fac)


def har_floor(es, pmin, k=K_T, g=G_T, n=N_T, pa=PA_T):
    """(S.50): x = varpi_f/p_a solves x^2 - a x^(2n) - b = 0, a = 3k(1-n) g es^2, b = (p_min/p_a)^2;
    n = 1/2 closed, x = [a + (a^2 + 4b)^(1/2)]/2; general n: the sheet's exact bracket [x_s, x_hi] bisected.
    Returns x, eps*_f = (p_min/p_a)/(k(1-n) x^n), eps_v,f = 1/(k(1-n)) - eps*_f, q_f = 3 g p_a es x^n,
    eps'_f = n p_min q_f/((1-n) p_a^2 x^2 + n p_min^2)."""
    kn = k * (1.0 - n)
    a = 3.0 * kn * g * es * es
    b = (pmin / pa) ** 2
    if n == 0.5:
        x = 0.5 * (a + math.sqrt(a * a + 4.0 * b))
    else:
        f = lambda t: t * t - a * t ** (2.0 * n) - b           # noqa: E731
        xs = (n * a) ** (1.0 / (2.0 - 2.0 * n)) if a > 0.0 else 0.0
        xhi = max((2.0 * a) ** (1.0 / (2.0 - 2.0 * n)) if a > 0.0 else 0.0, 2.0 * math.sqrt(b))
        lo, hi = xs, xhi
        assert f(lo) < 0.0 < f(hi)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        x = 0.5 * (lo + hi)
    est_f = (pmin / pa) / (kn * x ** n)
    q_f = 3.0 * g * pa * es * x ** n
    return dict(x=x, est_f=est_f, ev_f=1.0 / kn - est_f, q_f=q_f,
                epsp=n * pmin * q_f / ((1.0 - n) * pa * pa * x * x + n * pmin * pmin))


# ------------------------------------------------------------------------------------------------------------------
# surface / CSL closed forms
# ------------------------------------------------------------------------------------------------------------------
def eta_S12(M, N, p, pi):
    """(S.12)."""
    if N == 0.0:
        return M * (1.0 + math.log(pi / p))
    return (M / N) * (1.0 - (1.0 - N) * (p / pi) ** (N / (1.0 - N)))


def pi_S53(M, N, p_init, eta_star):
    """(S.53): the inverse of (S.12) at (p_init, eta*).  Refuses (ValueError) eta* >= M/N (N > 0)."""
    if N == 0.0:
        return p_init * math.exp(eta_star / M - 1.0)
    if not eta_star < M / N:
        raise ValueError("eta* >= M/N: no surface through the state")
    return p_init * ((1.0 - N) / (1.0 - eta_star * N / M)) ** ((1.0 - N) / N)


def w_ramp(N, c1, c2):
    """(S.56): W_ramp = 1 - pi_i(eta_1)/pi_i(eta_2) = 1 - [(1 - c2 N)/(1 - c1 N)]^((1-N)/N); N = 0: 1 - exp(-(c2 - c1))."""
    if N == 0.0:
        return 1.0 - math.exp(-(c2 - c1))
    return 1.0 - ((1.0 - c2 * N) / (1.0 - c1 * N)) ** ((1.0 - N) / N)


def v_for_psi(kw, pi, psi):
    """(S.22) fork CSL: psi_i = (v - 1) - e0 + lambda_c (-pi/p_a)^xi  =>  v for a target psi_i at pi_i = pi."""
    if kw["csl_mode"] == "paper":
        return psi + kw["v_c0"] - kw["lambda_tilde"] * math.log(-pi)
    return 1.0 + kw["e0"] + psi - kw["lambda_c"] * (-pi / kw["p_a"]) ** kw["xi"]


# ------------------------------------------------------------------------------------------------------------------
# tensor helpers
# ------------------------------------------------------------------------------------------------------------------
def sigma_pq(p, q, nhat):
    """Principal stress (3,) from (S.2): sigma_a = p + sqrt(2/3) q n_a, n_a with sum n_a = 0, |n| = 1."""
    nhat = np.asarray(nhat, float)
    nhat = nhat - nhat.mean()
    nhat = nhat / np.linalg.norm(nhat)
    return p * ONES + SQ23 * q * nhat


def theta_of_dir(nhat):
    """Lode angle (S.1) of a deviatoric direction (y = sum n^3, since |n| = 1)."""
    nhat = np.asarray(nhat, float)
    nhat = nhat - nhat.mean()
    nhat = nhat / np.linalg.norm(nhat)
    y = float((nhat ** 3).sum())
    return math.acos(max(-1.0, min(1.0, SQ6 * y))) / 3.0


def dir_for_theta(theta):
    """A unit deviatoric direction with the requested Lode angle, from the one-parameter family (cos t, sin t) in
    the plane of the two basis vectors (1,-1,0)/sqrt2, (1,1,-2)/sqrt6; theta is found by bisection."""
    e1 = np.array([1.0, -1.0, 0.0]) / math.sqrt(2.0)
    e2 = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)
    # t = 0 -> e2 (compression corner, theta = pi/3); t = pi/6 ... sweep t in [0, pi/3] for theta from pi/3 to 0
    f = lambda t: theta_of_dir(math.cos(t) * e2 + math.sin(t) * e1) - theta        # noqa: E731
    lo, hi = 0.0, math.pi / 3.0
    flo, fhi = f(lo), f(hi)
    assert flo * fhi < 0.0, (flo, fhi)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm < 0.0) == (flo < 0.0):
            lo = mid
        else:
            hi = mid
    t = 0.5 * (lo + hi)
    return math.cos(t) * e2 + math.sin(t) * e1


def mandel6(C4):
    """Mandel 6 x 6 of a 4th-order tensor (orthonormal basis e_i e_i, (e_i e_j + e_j e_i)/sqrt2)."""
    B = []
    for i in range(3):
        b = np.zeros((3, 3))
        b[i, i] = 1.0
        B.append(b)
    for (i, j) in ((0, 1), (1, 2), (0, 2)):
        b = np.zeros((3, 3))
        b[i, j] = b[j, i] = 1.0 / math.sqrt(2.0)
        B.append(b)
    M = np.zeros((6, 6))
    for I, bi in enumerate(B):
        for J, bj in enumerate(B):
            M[I, J] = np.einsum("ij,ijkl,kl->", bi, C4, bj)
    return M


def ok_state(oname, st):
    f = st.flags
    return (f.get("status", "ok") in ("ok", "initial")) if oname == "O1" else (not f.get("refused", False))


def why_state(oname, st):
    f = st.flags
    return f.get("status") if oname == "O1" else f.get("reason")


# ------------------------------------------------------------------------------------------------------------------
# round 3b additions (floor / pi_i0 files): WW shape (S.11), the yield function (S.12), BA06 closed forms (S.5),
# the states of the sheet's floor examples
# ------------------------------------------------------------------------------------------------------------------
def zeta_ww(theta, rho):
    """(S.11): zeta = [A c^2 + B^2]/[2(1-rho^2) c + B (A c^2 + 5 rho^2 - 4 rho)^(1/2)], c = cos(theta), A = 4(1-rho^2), B = 2 rho - 1."""
    c = math.cos(theta)
    A, B = 4.0 * (1.0 - rho * rho), 2.0 * rho - 1.0
    return (A * c * c + B * B) / (2.0 * (1.0 - rho * rho) * c + B * math.sqrt(A * c * c + 5.0 * rho * rho - 4.0 * rho))


def invariants_pq(sig):
    """(p, q, theta) of a principal (or symmetric 3 x 3) stress, (S.1); theta = 0 on the axis (q = 0)."""
    w = np.linalg.eigvalsh(np.asarray(sig, float)) if np.ndim(sig) == 2 else np.asarray(sig, float)
    p = float(w.mean())
    xi = w - p
    R = float(np.linalg.norm(xi))
    q = SQ32 * R
    if R == 0.0:
        return p, 0.0, 0.0
    y = float((xi ** 3).sum()) / R ** 3
    return p, q, math.acos(max(-1.0, min(1.0, SQ6 * y))) / 3.0


def F_yield(M, N, rho, sig, pi):
    """(S.12) F = zeta(theta, rho) q + p eta(p, pi_i) of a principal stress (WW shape)."""
    p, q, th = invariants_pq(sig)
    return zeta_ww(th, rho) * q + p * eta_S12(M, N, p, pi)


def ba06_pq(ev, es, p0=-100.0, kh=0.01, ev0=0.0, mu0=5400.0, a0=0.0):
    """(S.5): p = p0 e^w [1 + 3 a0 es^2/(2 kh)], q = 3 (mu0 - a0 p0 e^w) es, w = -(ev - ev0)/kh."""
    E = math.exp(-(ev - ev0) / kh)
    return p0 * E * (1.0 + 1.5 * a0 * es * es / kh), 3.0 * (mu0 - a0 * p0 * E) * es


def ba06_psi(ev, es, p0=-100.0, kh=0.01, ev0=0.0, mu0=5400.0, a0=0.0):
    """(S.4): Psi = Pt + 1.5 mu_e es^2, Pt = -p0 kh e^w, mu_e = mu0 + a0 Pt/kh."""
    Pt = -p0 * kh * math.exp(-(ev - ev0) / kh)
    return Pt + 1.5 * (mu0 + a0 * Pt / kh) * es * es


def invariants_eps(eps):
    """(eps_v, eps_s, n_hat) of a principal / symmetric elastic strain."""
    w = np.linalg.eigvalsh(np.asarray(eps, float)) if np.ndim(eps) == 2 else np.asarray(eps, float)
    ev = float(w.sum())
    e = w - ev / 3.0
    ne = float(np.linalg.norm(e))
    return ev, SQ23 * ne, (e / ne if ne > 0.0 else np.zeros(3))


# the off-corner deviatoric direction of the sheet's floor examples (floor_fd.py / floor_fd_har.py, sheet 9.7, 13.14b):
# theta = 0.271 (WW corners are 0 and pi/3 = 1.047)
OFFC_DIR = np.array([-1.0, -0.35, 1.35])
OFFC_DIR = (OFFC_DIR - OFFC_DIR.mean()) / np.linalg.norm(OFFC_DIR - OFFC_DIR.mean())
WET_DIR = np.array([-1.2, -0.1, 1.3])
WET_DIR = (WET_DIR - WET_DIR.mean()) / np.linalg.norm(WET_DIR - WET_DIR.mean())


def sigma_on_surface(M, N, rho, p, eta, nhat):
    """Principal stress with F = 0 at (p, eta): q = eta |p| / zeta(theta), the deviatoric direction nhat (WW)."""
    th = theta_of_dir(nhat)
    q = eta * abs(p) / zeta_ww(th, rho)
    return sigma_pq(p, q, nhat), th, q
