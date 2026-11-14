"""Acoustic tensor and localization search (sheet §14, AB06 83-97).

A_ik(n) = n_j a_ijkl n_l for any 4th-order tangent a (small-strain C of (S.33) or the
finite-strain a^ep = c~ + tau(+)1 of (S.34); the caller chooses which to pass).
Search: coarse sweep of (theta, phi) in [0, pi]^2 with n = (sin theta sin phi, cos phi,
cos theta sin phi) in the GLOBAL basis, then Nelder-Mead refinement from the best
N_REFINE grid points. det A is a smooth function of (theta, phi), so this replaces AB06's
Newton (90-97) with the same minimum.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .params import Params

GRID_N = 91          # (theta, phi) grid points per axis on [0, pi]  (2 degree spacing)
N_REFINE = 3         # best grid points refined
REFINE_XATOL = 1.0e-10
REFINE_FATOL = 1.0e-14


def acoustic_tensor(a4: np.ndarray, n: np.ndarray) -> np.ndarray:
    return np.einsum("j,ijkl,l->ik", n, a4, n)


def acoustic_principal(c_ab: np.ndarray, gam_ab: np.ndarray, tau: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    """(S.44) transcription, in the principal basis (self-check against acoustic_tensor)."""
    A = np.empty((3, 3))
    sn = float(np.sum(tau * alpha ** 2))
    for a in range(3):
        for b in range(3):
            if a == b:
                A[a, a] = alpha[a] ** 2 * c_ab[a, a] + sn + sum(alpha[c] ** 2 * gam_ab[c, a] for c in range(3) if c != a)
            else:
                A[a, b] = alpha[a] * (c_ab[a, b] + gam_ab[a, b]) * alpha[b]
    return A


def _n_of(th: float, ph: float) -> np.ndarray:
    return np.array([np.sin(th) * np.sin(ph), np.cos(ph), np.cos(th) * np.sin(ph)])


def acoustic_min_det(params: Params, state, tangent4: np.ndarray, grid_n: int = GRID_N,
                     n_refine: int = N_REFINE):
    """(min_det, n_vec): the minimum over unit n of det(n.a.n) and its direction."""
    a4 = np.asarray(tangent4, float)

    def f(x):
        return float(np.linalg.det(acoustic_tensor(a4, _n_of(x[0], x[1]))))

    ths = np.linspace(0.0, np.pi, grid_n)
    vals = np.empty((grid_n, grid_n))
    for i, th in enumerate(ths):
        for j, ph in enumerate(ths):
            vals[i, j] = f((th, ph))
    flat = np.argsort(vals, axis=None)[:n_refine]
    best_v, best_x = np.inf, None
    for k in flat:
        i, j = np.unravel_index(k, vals.shape)
        r = minimize(f, x0=np.array([ths[i], ths[j]]), method="Nelder-Mead",
                     options=dict(xatol=REFINE_XATOL, fatol=REFINE_FATOL, maxiter=2000))
        if r.fun < best_v:
            best_v, best_x = r.fun, r.x
    return best_v, _n_of(best_x[0], best_x[1])
