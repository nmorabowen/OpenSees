"""ANALYTIC channel B for LadrunoJ2Finite — the numpy oracle for replacing the numeric
R-perturbation channel-B in the C++ consistent tangent (ADR-16 P3, "analytic channel B").

Channel B is the part of ∂σ/∂F that flows through the CO-ROTATED backstress
α̃_m = R α_m Rᵀ (R = polar(f_Δ)), with the trial elastic strain εᵉᵗʳ and J held fixed.
The shipped material computes it by finite-differencing R (≈9–18 return-map calls per GP
tangent). This module computes the SAME tensor analytically as the chain

    ∂σ/∂F = (1/J) Σ_m (∂τ/∂α̃_m) : (∂α̃_m/∂R) : (∂R/∂F)

with two analytic pieces:

  (1) ∂R/∂F — the polar-rotation derivative. For f = R U (polar), a perturbation df
      gives dR = R Ω with Ω skew solving the Sylvester equation Ω U + U Ω = A − Aᵀ,
      A = Rᵀ df. In axial-vector form (ω of Ω, w of A−Aᵀ): (tr(U) I − U) ω = w — a 3×3
      SPD solve (eigenvalues tr(U)−u_i = sum of the other two stretches > 0).

  (2) ∂τ/∂α̃_m — the return-map sensitivity to the m-th INPUT backstress at fixed strain.
      With s_tr = 2G dev(εᵉᵗʳ) (α-independent), M = s_tr − Σ_p α̃_p/D_p,
      D_p = 1+√⅔ γ_p Δγ, n = M/‖M‖, s = s_tr − 2G Δγ n:
        ∂Δγ/∂α̃_m = −n/(h D_m)                          (h = the kernel hardening modulus)
        dM(dα̃_m)  = −dα̃_m/D_m + Mp·dΔγ                 (Mp = Σ_p α̃_p √⅔ γ_p/D_p²)
        dn        = (dM − n (n:dM))/‖M‖
        ds = dτ_dev = −2G[ dΔγ·n + Δγ·dn ]             (hydrostatic part is α-independent)

Verified (tests/test_ladrunoJ2_finite_channelB.py): the polar derivative matches a
central FD of polar_rotation to ~1e-9, and channelB_analytic matches the oracle's
numeric channelB_dsigma_dF to ~1e-8 across stiff/soft × first-yield/near-saturation ×
rotated plastic states.
"""
from __future__ import annotations

import numpy as np

import ladrunoj2_reference as lr
import ladrunoj2_finite_native_reference as nat

_I3 = np.eye(3)
_ROOT23 = np.sqrt(2.0 / 3.0)


def _axial(W):
    """Axial vector w of a skew tensor W (W x = w × x): w = (W21, W02, W10)."""
    return np.array([W[2, 1], W[0, 2], W[1, 0]])


def _skew_from_axial(w):
    return np.array([[0.0, -w[2], w[1]],
                     [w[2], 0.0, -w[0]],
                     [-w[1], w[0], 0.0]])


def polar_derivative(f, df):
    """dR for f = R U (polar) given df. Solves (tr U·I − U) ω = axial(A − Aᵀ),
    A = Rᵀ df; returns dR = R Ω, Ω = skew(ω)."""
    R = nat.polar_rotation(f)
    U = R.T @ f                      # symmetric SPD (right stretch)
    U = 0.5 * (U + U.T)
    A = R.T @ df
    w = _axial(A - A.T)              # axial of the skew part (×2 folded in)
    K = np.trace(U) * _I3 - U        # SPD
    omega = np.linalg.solve(K, w)
    Omega = _skew_from_axial(omega)
    return R @ Omega


def _aux(p: lr.Params, eps_tr, ebarP_n, alpha_pf):
    """Recompute the converged return-map auxiliaries (n, ‖M‖, h, D_k, Mp, s_tr, Δγ)
    needed for the analytic backstress sensitivity. Returns None if elastic."""
    res = lr.return_map(p, eps_tr, np.zeros((3, 3)), ebarP_n, alpha_pf)
    if not res["plastic"]:
        return None
    dG = res["dGamma"]
    G = p.G
    tr = np.trace(eps_tr)
    edev = eps_tr - (tr / 3.0) * _I3
    s_tr = 2.0 * G * edev
    Dk = [1.0 + _ROOT23 * p.gam[k] * dG for k in range(p.nBack)]
    M = s_tr.copy()
    Mp = np.zeros((3, 3))
    for k in range(p.nBack):
        M -= alpha_pf[k] / Dk[k]
        Mp += alpha_pf[k] * (_ROOT23 * p.gam[k] / Dk[k] ** 2)
    normM = np.sqrt(np.tensordot(M, M))
    n = M / normM
    dtheta = 2.0 * G
    for k in range(p.nBack):
        dtheta += (2.0 / 3.0) * p.C[k] / Dk[k] ** 2
    pbar = ebarP_n + _ROOT23 * dG
    h = dtheta + (2.0 / 3.0) * p.yield_slope(pbar) - np.tensordot(n, Mp)
    return dict(dG=dG, n=n, normM=normM, h=h, Dk=Dk, Mp=Mp, G=G)


def _dtau_dev(aux, dalpha_m, m):
    """dτ_dev (= ds) for a perturbation dalpha_m of the m-th INPUT backstress α̃_m."""
    n, Mp, G, dG = aux["n"], aux["Mp"], aux["G"], aux["dG"]
    Dm, h, normM = aux["Dk"][m], aux["h"], aux["normM"]
    dGamma = -np.tensordot(n, dalpha_m) / (h * Dm)          # ∂Δγ/∂α̃_m : dα̃_m
    dM = -dalpha_m / Dm + Mp * dGamma
    dn = (dM - n * np.tensordot(n, dM)) / normM
    return -2.0 * G * (dGamma * n + dG * dn)


def channelB_analytic(model: nat.NativeFiniteJ2, F):
    """Analytic channel-B ∂σ/∂F at the model's current committed state — matches the
    oracle's numeric NativeFiniteJ2.channelB_dsigma_dF(F) (εᵉᵗʳ, J held fixed)."""
    p = model.p
    Fn = model.Fn
    fD = F @ np.linalg.inv(Fn)
    R = nat.polar_rotation(fD)
    eps_tr = __import__("logstrain_reference").hencky_strain(fD @ model.Be_n @ fD.T)
    J0 = float(np.linalg.det(F))

    alpha_pf = [R @ a @ R.T for a in model.alpha]
    aux = _aux(p, eps_tr, model.ebarP, alpha_pf)
    B = np.zeros((3, 3, 3, 3))
    if aux is None:                                          # elastic ⇒ no channel B
        return B

    Fninv = np.linalg.inv(Fn)
    for k in range(3):
        for l in range(3):
            dF = np.zeros((3, 3)); dF[k, l] = 1.0
            dR = polar_derivative(fD, dF @ Fninv)
            dsig = np.zeros((3, 3))
            for m, a in enumerate(model.alpha):
                dalpha_m = dR @ a @ R.T + R @ a @ dR.T       # ∂α̃_m for this dF
                dsig += _dtau_dev(aux, dalpha_m, m)
            B[:, :, k, l] = dsig / J0
    return B
