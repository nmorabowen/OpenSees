"""WP-232i: the manufactured solution of the MPM convergence gates (numpy only).

After Wallstedt & Guilkey (2008), on the unit square (plane strain) or the unit cube, with rollers
on every face. For each axis i (no sum):

    u_i(X, t) = A sin(2 pi X_i) s(t),      s(t) = sin(omega t),  omega = c pi

Every X_i -> u_i is odd about X_i = 0 and X_i = 1, so the normal displacement vanishes on every
face, and u_i does not depend on the other coordinates, so F is diagonal: no rotation, no shear
stress, no shear traction. The rollers are exact boundary conditions.

Derivation (by hand; sympy is not installed. tests/test_mpm_mms_reference.py checks it by
centred finite differences through an independent route: F from differences of u, the stress
from a matrix logarithm by eigen-decomposition, the divergence by differences of P).

1. Kinematics. F_ii = 1 + du_i/dX_i = 1 + 2 pi A cos(2 pi X_i) s(t); F_ij = 0 for i != j.
   In plane strain F_33 = 1 (tests/mpm_interface.md §8.3). J = prod_k F_kk.
2. The strain the engine integrates. The hypo kernel accumulates Hughes-Winget increments in
   the unrotated frame (interface §8.2): G = (F1 - Fn) F_mid^-1, d eps_hat = R_mid^T sym(G) R_mid.
   For a diagonal F, R = I and d eps_hat_ii = 2 (F1 - Fn) / (F1 + Fn) = ln(F1 / Fn) + O(dF^3),
   so the accumulated strain feed is ln F_ii up to O(dt^2) (a midpoint rule on d ln F).
3. The stress. ElasticIsotropic receives the TOTAL feed (interface §9) and returns
   sigma_hat = lambda tr(eps) I + 2 mu eps; the engine pushes it with R (= I) and uses it as the
   CAUCHY stress, with the current volume V = J V0 (interface §3, §6.3). So
       sigma_ii = lambda sum_k ln F_kk + 2 mu ln F_ii,   sigma_ij = 0 (i != j).
   (sigma_33 in plane strain is the material's business and never enters the in-plane balance.)
4. First Piola-Kirchhoff stress P = J sigma F^-T is diagonal: P_ii = J sigma_ii / F_ii
   = (prod_{k != i} F_kk) sigma_ii.
5. Divergence. (Div_X P)_i = sum_j dP_ij/dX_j = dP_ii/dX_i. The prefactor prod_{k != i} F_kk does
   not depend on X_i; sigma_ii depends on X_i only through ln F_ii, with coefficient
   (lambda + 2 mu) (the ln F_ii term appears in the trace with lambda and alone with 2 mu). So
       (Div_X P)_i = (prod_{k != i} F_kk) (lambda + 2 mu) (dF_ii/dX_i) / F_ii,
       dF_ii/dX_i = -4 pi^2 A sin(2 pi X_i) s(t).
6. Acceleration: a_i = d2u_i/dt2 = -omega^2 A sin(2 pi X_i) s(t).
7. The body force per unit mass that makes u exact (balance of linear momentum in the
   reference configuration, rho0 a = Div_X P + rho0 b):
       b_i = a_i - (Div_X P)_i / rho0.
8. Initial conditions: u(X, 0) = 0, v_i(X, 0) = omega A sin(2 pi X_i), a(X, 0) = 0, sigma = 0.
9. Quasi-static (gate G4): u_i(X, lam) = lam A sin(2 pi X_i); replace s(t) by lam and drop the
   acceleration: b_i = -(Div_X P)_i / rho0.

At t_end = one period (2 / c) the exact u is 0. The gate error is
    e = sqrt(mean_p |u_p - u(X_p, t_end)|^2) / A,  u_p = x_p - X_p.
"""
from __future__ import annotations

import math

import numpy as np

TWO_PI = 2.0 * math.pi


def lame(E: float, nu: float) -> tuple[float, float]:
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    return lam, mu


def amplitude(fmax: float) -> float:
    """A such that max |F_ii - 1| = 2 pi A = fmax."""
    return fmax / TWO_PI


def displacement(X: np.ndarray, s: float, A: float) -> np.ndarray:
    """u(X) at the time factor s = sin(omega t) (or the load factor lam). X: [N, ndm]."""
    return A * np.sin(TWO_PI * X) * s


def velocity0(X: np.ndarray, A: float, omega: float) -> np.ndarray:
    return omega * A * np.sin(TWO_PI * X)


def stretches(X: np.ndarray, s: float, A: float) -> np.ndarray:
    """F_ii, [N, ndm]."""
    return 1.0 + TWO_PI * A * np.cos(TWO_PI * X) * s


def cauchy_diag(X: np.ndarray, s: float, A: float, lam: float, mu: float) -> np.ndarray:
    """sigma_ii (in-plane / 3D diagonal), [N, ndm] (derivation step 3)."""
    lnF = np.log(stretches(X, s, A))
    return lam * lnF.sum(axis=1, keepdims=True) + 2.0 * mu * lnF


def div_P(X: np.ndarray, s: float, A: float, lam: float, mu: float) -> np.ndarray:
    """(Div_X P)_i, [N, ndm] (derivation step 5)."""
    F = stretches(X, s, A)
    dF = -(TWO_PI ** 2) * A * np.sin(TWO_PI * X) * s
    out = np.empty_like(X, dtype=float)
    nd = X.shape[1]
    for i in range(nd):
        pre = np.ones(X.shape[0])
        for k in range(nd):
            if k != i:
                pre = pre * F[:, k]
        out[:, i] = pre * (lam + 2.0 * mu) * dF[:, i] / F[:, i]
    return out


def acceleration(X: np.ndarray, t: float, A: float, omega: float) -> np.ndarray:
    return -(omega ** 2) * A * np.sin(TWO_PI * X) * math.sin(omega * t)


def body_force(X: np.ndarray, t: float, A: float, omega: float, rho0: float, lam: float,
               mu: float) -> np.ndarray:
    """b(X, t) per unit mass, dynamic (derivation step 7). [N, ndm]."""
    s = math.sin(omega * t)
    return acceleration(X, t, A, omega) - div_P(X, s, A, lam, mu) / rho0


def body_force_static(X: np.ndarray, lam_load: float, A: float, rho0: float, lam: float,
                      mu: float) -> np.ndarray:
    """b(X, lam) per unit mass, quasi-static (derivation step 9). [N, ndm]."""
    return -div_P(X, lam_load, A, lam, mu) / rho0


def error(x: np.ndarray, X: np.ndarray, s: float, A: float) -> float:
    """The gate error e (module docstring)."""
    d = (x - X) - displacement(X, s, A)
    return float(np.sqrt(np.mean(np.sum(d * d, axis=1))) / A)


def order(hs, es, last: int = 3) -> float:
    """Least-squares slope of log e against log h over the `last` finest rungs."""
    h = np.log(np.asarray(hs, dtype=float)[-last:])
    e = np.log(np.asarray(es, dtype=float)[-last:])
    return float(np.polyfit(h, e, 1)[0])
