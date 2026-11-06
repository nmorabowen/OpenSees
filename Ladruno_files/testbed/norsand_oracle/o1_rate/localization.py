"""Acoustic-tensor localization analysis (sheet 14, (S.34), (S.44)) and the K2 driver.

acoustic_min_det(params, state, tangent4, finite=True)
    finite=True  (default, the K2 case): tangent4 is dtau/d(eps) on LOG strains (the
                 small-strain constitutive tangent run on log strains).  Its principal
                 block a~_ab = n^a n^a : tangent4 : n^b n^b (principal basis of eps^e) is
                 assembled into the spatial tangent a = c~ + tau(+)1 by (S.34):
                   c~_ab = a~_ab - 2 tau_a delta_ab,
                   g~_ab = (tau_b lam_a^2 - tau_a lam_b^2)/(lam_b^2 - lam_a^2),
                   lam_a = exp(eps^e_a)  (elastic stretches; the continuum analogue of the
                   trial stretches), repeated-stretch limit (a~_bb - a~_ba)/2 - tau_a;
                 and A(n) by (S.44).  The spin part of tangent4 is NOT used.
    finite=False (AB06 Remark 3): A_ik = n_j C_ijkl n_l with C = tangent4 as given.
    Search: coarse (theta, phi) grid on [0, pi]^2 with alpha = (sin t sin f, cos f, cos t sin f)
    in the principal basis (AB06 88), then Nelder-Mead refinement from the 4 best grid
    points (in place of AB06's Newton, eqs 90-97).  Returns (min_det, n) with n the global
    unit normal.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize

from .api import run_path, tangent
from .params import Params

STRETCH_TOL = 1e-10


def _principal(state):
    w, V = np.linalg.eigh(state.eps_e)
    tau = np.array([V[:, a] @ state.sigma @ V[:, a] for a in range(3)])
    return w, V, tau


def _alpha(th, ph):
    return np.array([np.sin(th) * np.sin(ph), np.cos(ph), np.cos(th) * np.sin(ph)])


def finite_blocks(state, tangent4):
    """(c~_ab, g~_ab, tau_a, V) of (S.34)."""
    w, V, tau = _principal(state)
    at = np.einsum("ia,ja,ijkl,kb,lb->ab", V, V, tangent4, V, V)
    lam2 = np.exp(2.0 * w)
    ct = at - 2.0 * np.diag(tau)
    g = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            if a == b:
                continue
            if abs(math.exp(w[a]) - math.exp(w[b])) < STRETCH_TOL:
                g[a, b] = 0.5 * (at[b, b] - at[b, a]) - tau[a]
            else:
                g[a, b] = (tau[b] * lam2[a] - tau[a] * lam2[b]) / (lam2[b] - lam2[a])
    return ct, g, tau, V


def finite_spatial_tangent(state, tangent4):
    """Full a = c~ + tau(+)1 as a 4th-order tensor in the global basis (for FD checks)."""
    ct, g, tau, V = finite_blocks(state, tangent4)
    m = [np.outer(V[:, a], V[:, a]) for a in range(3)]
    C = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            C += ct[a, b] * np.einsum("ij,kl->ijkl", m[a], m[b])
            if a != b:
                mab = np.outer(V[:, a], V[:, b])
                mba = np.outer(V[:, b], V[:, a])
                C += g[a, b] * (np.einsum("ij,kl->ijkl", mab, mab)
                                + np.einsum("ij,kl->ijkl", mab, mba))
    T = state.sigma
    C += np.einsum("jl,ik->ijkl", T, np.eye(3))
    return C


def _det_fun_finite(ct, g, tau):
    def detA(th, ph):
        al = _alpha(th, ph)                      # (3, ...) broadcast
        a2 = al * al
        sn = np.tensordot(tau, a2, axes=(0, 0))
        A = np.empty((3, 3) + np.shape(th))
        for a in range(3):
            A[a, a] = a2[a] * ct[a, a] + sn + sum(a2[c] * g[c, a] for c in range(3) if c != a)
            for b in range(3):
                if b != a:
                    A[a, b] = al[a] * (ct[a, b] + g[a, b]) * al[b]
        return (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
                - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
                + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))
    return detA


def _det_fun_small(Cp):
    def detA(th, ph):
        al = _alpha(th, ph)
        A = np.einsum("j...,ijkl,l...->ik...", al, Cp, al)
        return (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
                - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
                + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))
    return detA


def acoustic_min_det(params: Params, state, tangent4, finite: bool = True, n_grid: int = 181):
    if finite:
        ct, g, tau, V = finite_blocks(state, tangent4)
        detA = _det_fun_finite(ct, g, tau)
    else:
        _, V, _ = _principal(state)
        Cp = np.einsum("ia,jb,kc,ld,ijkl->abcd", V, V, V, V, tangent4)
        detA = _det_fun_small(Cp)
    th = np.linspace(0.0, math.pi, n_grid)
    TH, PH = np.meshgrid(th, th, indexing="ij")
    D = detA(TH, PH)
    order = np.argsort(D, axis=None)[:4]
    best = (float(D.flat[order[0]]), float(TH.flat[order[0]]), float(PH.flat[order[0]]))
    for k in order:
        x0 = [TH.flat[k], PH.flat[k]]
        r = minimize(lambda x: float(detA(np.array(x[0]), np.array(x[1]))), x0,
                     method="Nelder-Mead",
                     options=dict(xatol=1e-12, fatol=1e-14 * max(1.0, abs(best[0])),
                                  maxiter=4000))
        if r.fun < best[0]:
            best = (float(r.fun), float(r.x[0]), float(r.x[1]))
    n = V @ _alpha(best[1], best[2])
    return best[0], n / np.linalg.norm(n)


def k2_path(params: Params, state0, n_max: int, n1: int = 10, lam1: float = 1e-3,
            lam2: float = 4e-4, rtol: float = 1e-10, extra_after: int = 2, n_grid: int = 181):
    """AB06 6.1 protocol (S.43) on log strains (v = v0 J), min det after every step.

    Returns dict with
      n_first   : first step n with min det <= 0 (criterion 'first step'),
      n_interp  : linear-interpolated zero crossing of the normalised curve (float),
      n_interp_step : round(n_interp) (criterion 'interpolated', nearest step),
      mindet, mindet_norm (normalised by the step-n1 value), normals, states, status.
    Integration stops extra_after steps past the first crossing, or at n_max.
    """
    d1 = np.diag([math.log(1.0 + lam2), math.log(1.0 - lam1), 0.0])
    d2 = np.diag([0.0, math.log(1.0 - lam2), math.log(1.0 + lam1)])
    st = state0
    states, dets, normals = [], [], []
    n_first = None
    status = "ok"
    for n in range(1, n_max + 1):
        d = d1 if n <= n1 else d2
        st = run_path(params, st, d[None], rtol=rtol, kin="log")[-1]
        st.flags["step"] = n
        if st.flags["status"] != "ok":
            status = f"step {n}: {st.flags['status']}"
            states.append(st)
            break
        C = tangent(params, st)
        md, nv = acoustic_min_det(params, st, C, finite=True, n_grid=n_grid)
        states.append(st)
        dets.append(md)
        normals.append(nv)
        if n_first is None and md <= 0.0:
            n_first = n
        if n_first is not None and n >= n_first + extra_after:
            break
    dets = np.array(dets)
    ref = dets[n1 - 1] if len(dets) >= n1 else float("nan")
    norm = dets / ref
    n_interp = None
    for k in range(1, len(dets)):
        if dets[k] <= 0.0 < dets[k - 1]:
            n_interp = k + dets[k - 1] / (dets[k - 1] - dets[k])     # steps are k, k+1 (1-based)
            break
    return dict(n_first=n_first, n_interp=n_interp,
                n_interp_step=(int(round(n_interp)) if n_interp is not None else None),
                mindet=dets, mindet_norm=norm, normals=normals, states=states, status=status)
