"""Reference (oracle) implementation of the log-strain finite-strain adaptor.

Pure-numpy mirror of the spatial multiplicative (logarithmic / Hencky) MATISU
algorithm — de Souza Neto, Perić & Owen, *Computational Methods for Plasticity*
(2008). This module is the numerical oracle the C++ `LogStrainNDMaterial`
(SRC/material/nD/, forthcoming) must match, and it de-risks the #1 numerical
trap: the eigenvalue-degeneracy branches of the tensor-log derivative.

No OpenSees / openseespy dependency — numpy only — so it runs standalone in CI.

Key references (book page → PDF page = book + 25):
  * Box 14.3 (MATISU), p.596         — bᵉᵗʳ = F_Δ bᵉ_n F_Δᵀ; εᵉ = ½ln bᵉ; σ = J⁻¹τ
  * eq. (14.30), p.582               — εᵉ = ½ ln Bᵉ  (Eulerian Hencky strain)
  * Appendix A.5.2 / Box A.6, p.743  — isotropic tensor function Y(X)=Σ y(xᵢ)Eᵢ
  * eqs (A.51), (A.52), (A.53)       — Y and dY/dX with the three eigenvalue
                                       cases (all-distinct / two-equal / triple)
  * eq. (A.46)                       — dX²/dX

Conventions: 2nd-order tensors are (3,3) ndarrays; 4th-order are (3,3,3,3).
A double contraction D:H is einsum('ijkl,kl->ij', D, H).
"""
from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------- #
#  Fourth-order building blocks                                               #
# --------------------------------------------------------------------------- #
_I = np.eye(3)


def _dyad(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """(A ⊗ B)_ijkl = A_ij B_kl."""
    return np.einsum("ij,kl->ijkl", A, B)


def sym_identity() -> np.ndarray:
    """Symmetric fourth-order identity I_S: ½(δ_ik δ_jl + δ_il δ_jk)."""
    d = _I
    return 0.5 * (np.einsum("ik,jl->ijkl", d, d) + np.einsum("il,jk->ijkl", d, d))


def dX2_dX(X: np.ndarray) -> np.ndarray:
    """Derivative of X² w.r.t. X, eq. (A.46):
    [dX²/dX]_ijkl = ½(δ_ik X_lj + δ_il X_kj + δ_jl X_ik + δ_jk X_il)."""
    d = _I
    return 0.5 * (
        np.einsum("ik,lj->ijkl", d, X)
        + np.einsum("il,kj->ijkl", d, X)
        + np.einsum("jl,ik->ijkl", d, X)
        + np.einsum("jk,il->ijkl", d, X)
    )


_IS = sym_identity()
_IxI = _dyad(_I, _I)


# --------------------------------------------------------------------------- #
#  Spectral decomposition (Box A.5)                                           #
# --------------------------------------------------------------------------- #
def spectral(X: np.ndarray):
    """Eigenvalues and eigenprojections of a symmetric 3×3 tensor.

    Returns (vals, projs) with vals a (3,) array and projs a list of three
    (3,3) eigenprojections Eᵢ = eᵢ⊗eᵢ.  numpy's symmetric eigensolver is used
    (robust, and for the *forward* function Σ y(xᵢ)Eᵢ the result is independent
    of the eigenbasis chosen within any degenerate eigenspace).
    """
    vals, vecs = np.linalg.eigh(0.5 * (X + X.T))
    projs = [np.outer(vecs[:, i], vecs[:, i]) for i in range(3)]
    return vals, projs


def _multiplicity(vals: np.ndarray, rtol: float):
    """Classify eigenvalue multiplicity (Remark A.2 relative test).

    Returns one of:
      ("distinct", None)
      ("double", a)   — index a is the singleton; the other two are equal
      ("triple", None)
    """
    scale = max(1.0, float(np.max(np.abs(vals))))
    eq = [abs(vals[i] - vals[j]) <= rtol * scale for (i, j) in ((0, 1), (1, 2), (2, 0))]
    e01, e12, e20 = eq
    if not (e01 or e12 or e20):
        return "distinct", None
    if e01 and e12 and e20:
        return "triple", None
    # exactly one pair equal → the singleton is the index not in that pair
    if e01:
        return "double", 2
    if e12:
        return "double", 0
    return "double", 1  # e20


# --------------------------------------------------------------------------- #
#  Isotropic tensor function  Y(X) = Σ y(xᵢ) Eᵢ      (A.51)                    #
# --------------------------------------------------------------------------- #
def iso_function(X: np.ndarray, y) -> np.ndarray:
    """Apply scalar function `y` to a symmetric tensor spectrally (A.51).

    Valid for any multiplicity — the grouped sum over a degenerate eigenspace is
    basis-independent.
    """
    vals, projs = spectral(X)
    return sum(y(vals[i]) * projs[i] for i in range(3))


# --------------------------------------------------------------------------- #
#  Derivative of an isotropic tensor function  dY/dX   (A.52 / A.53)          #
# --------------------------------------------------------------------------- #
def iso_function_derivative(X: np.ndarray, y, yp, rtol: float = 1e-9) -> np.ndarray:
    """Closed-form fourth-order derivative dY/dX of Y(X)=Σ y(xᵢ)Eᵢ, Box A.6.

    `y`  : scalar function,  `yp` : its first derivative.
    Branches on eigenvalue multiplicity (the D3 degeneracy handling) — analytic
    limits, NOT ε-perturbation.
    """
    vals, projs = spectral(X)
    kind, a_sing = _multiplicity(vals, rtol)
    dX2 = dX2_dX(X)

    if kind == "triple":
        # all three equal: D = y'(x) I_S            (A.52, third branch)
        return yp(vals[0]) * _IS

    if kind == "double":
        # x_a (singleton) ≠ x_b = x_c (repeated):    (A.52 second branch, A.53)
        a = a_sing
        c = (a + 1) % 3  # any of the repeated pair
        xa, xc = vals[a], vals[c]
        ya, yc = y(xa), y(xc)
        ypa, ypc = yp(xa), yp(xc)
        d = xa - xc
        s1 = (ya - yc) / d**2 - ypc / d
        s2 = 2.0 * xc * (ya - yc) / d**2 - (xa + xc) / d * ypc
        s3 = 2.0 * (ya - yc) / d**3 - (ypa + ypc) / d**2
        s4 = xc * s3
        s5 = xc * s3
        s6 = xc * xc * s3
        return (
            s1 * dX2
            - s2 * _IS
            - s3 * _dyad(X, X)
            + s4 * _dyad(X, _I)
            + s5 * _dyad(_I, X)
            - s6 * _IxI
        )

    # all distinct:                                  (A.52 first branch)
    D = np.zeros((3, 3, 3, 3))
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        xa, xb, xc = vals[a], vals[b], vals[c]
        Ea, Eb, Ec = projs[a], projs[b], projs[c]
        coef = y(xa) / ((xa - xb) * (xa - xc))
        D += coef * (
            dX2
            - (xb + xc) * _IS
            - ((xa - xb) + (xa - xc)) * _dyad(Ea, Ea)
            - (xb - xc) * (_dyad(Eb, Eb) - _dyad(Ec, Ec))
        )
        D += yp(xa) * _dyad(Ea, Ea)
    return D


# --------------------------------------------------------------------------- #
#  Scalar functions of interest                                               #
# --------------------------------------------------------------------------- #
def half_log(x: float) -> float:      # εᵉ = ½ ln bᵉ
    return 0.5 * np.log(x)


def half_log_prime(x: float) -> float:
    return 0.5 / x


def hencky_strain(Be: np.ndarray) -> np.ndarray:
    """Eulerian logarithmic (Hencky) strain εᵉ = ½ ln Bᵉ, eq. (14.30)."""
    return iso_function(Be, half_log)


def dHencky_dBe(Be: np.ndarray, rtol: float = 1e-9) -> np.ndarray:
    """∂εᵉ/∂Bᵉ — the log-map derivative (the L tensor of eq. 14.101, /2)."""
    return iso_function_derivative(Be, half_log, half_log_prime, rtol)


# --------------------------------------------------------------------------- #
#  Elastic Hencky constitutive forward map (MATISU, elastic inner material)   #
# --------------------------------------------------------------------------- #
def isotropic_elastic_kirchhoff(eps: np.ndarray, K: float, G: float) -> np.ndarray:
    """Small-strain isotropic elastic law applied to the Hencky strain → τ.

    τ = K tr(ε) I + 2G dev(ε)   (Kirchhoff stress; the UNCHANGED small-strain
    elastic law, eq. 14.33 with the Hencky hyperelastic potential).
    """
    tr = np.trace(eps)
    dev = eps - (tr / 3.0) * _I
    return K * tr * _I + 2.0 * G * dev


def cauchy_from_F_elastic(F: np.ndarray, K: float, G: float,
                          Be_n: np.ndarray | None = None,
                          F_n: np.ndarray | None = None):
    """MATISU forward map for an elastic Hencky material: F → (σ, τ, J, εᵉ).

    For a single step from the reference configuration pass Be_n=I, F_n=I, so
    Bᵉᵗʳ = F Fᵀ.  General step: Bᵉᵗʳ = F_Δ Bᵉ_n F_Δᵀ with F_Δ = F F_n⁻¹.
    """
    if Be_n is None:
        Be_n = _I.copy()
    if F_n is None:
        F_n = _I.copy()
    F_delta = F @ np.linalg.inv(F_n)
    Be_tr = F_delta @ Be_n @ F_delta.T
    eps = hencky_strain(Be_tr)                  # εᵉ = ½ ln Bᵉᵗʳ
    tau = isotropic_elastic_kirchhoff(eps, K, G)  # Kirchhoff (elastic: εᵉ=εᵉᵗʳ)
    J = float(np.linalg.det(F))
    sigma = tau / J                              # Box 14.3 (iv)
    return sigma, tau, J, eps


# --------------------------------------------------------------------------- #
#  Small-strain consistent tangents (∂τ/∂εᵉᵗʳ) — the "D" of eq. (14.100)      #
# --------------------------------------------------------------------------- #
_I_DEV = _IS - (1.0 / 3.0) * _IxI       # deviatoric projector (4th order)


def elastic_tangent_D(K: float, G: float) -> np.ndarray:
    """Isotropic elastic consistent tangent ∂τ/∂εᵉ = K I⊗I + 2G I_dev."""
    return K * _IxI + 2.0 * G * _I_DEV


def j2_return_map(eps_tr: np.ndarray, alpha_n: float,
                  K: float, G: float, sig_y0: float, H: float):
    """Small-strain von-Mises (J2) radial return with linear isotropic hardening.

    Box 14.4 form: the *trial elastic* Hencky strain goes in, the updated elastic
    strain comes out (the elastic strain IS the state — no separate εᵖ stored).

    Returns (tau, eps_e, alpha, D, plastic):
      tau     Kirchhoff stress (3×3)
      eps_e   updated elastic Hencky strain (3×3)
      alpha   updated accumulated plastic strain
      D       consistent tangent ∂τ/∂εᵉᵗʳ (3,3,3,3)
      plastic bool
    Yield: Φ = ‖s‖ − √(2/3)(σ_y0 + H·α);  flow N = s/‖s‖;  Δεᵖ = Δγ N (deviatoric).
    """
    tr = np.trace(eps_tr)
    eps_dev = eps_tr - (tr / 3.0) * _I
    s_tr = 2.0 * G * eps_dev
    q = float(np.sqrt(np.tensordot(s_tr, s_tr)))      # Frobenius ‖s_tr‖
    sqrt23 = np.sqrt(2.0 / 3.0)
    f_tr = q - sqrt23 * (sig_y0 + H * alpha_n)

    if f_tr <= 0.0 or q == 0.0:                        # elastic step
        tau = K * tr * _I + s_tr
        return tau, eps_tr.copy(), alpha_n, elastic_tangent_D(K, G), False

    dgamma = f_tr / (2.0 * G + (2.0 / 3.0) * H)
    N = s_tr / q
    s = s_tr - 2.0 * G * dgamma * N
    tau = K * tr * _I + s
    eps_e = eps_tr - dgamma * N                        # deviatoric ⇒ tr unchanged
    alpha = alpha_n + sqrt23 * dgamma

    # consistent elastoplastic tangent (Simo & Hughes; dSNPO §7.4.2)
    theta = 1.0 - 2.0 * G * dgamma / q
    theta_bar = 1.0 / (1.0 + H / (3.0 * G)) - (1.0 - theta)
    D = (K * _IxI
         + 2.0 * G * theta * _I_DEV
         - 2.0 * G * theta_bar * _dyad(N, N))
    return tau, eps_e, alpha, D, True


# --------------------------------------------------------------------------- #
#  Consistent SPATIAL tangent  a  (§14.5, eqs 14.99–14.102)                   #
# --------------------------------------------------------------------------- #
def log_derivative_L(Be: np.ndarray, rtol: float = 1e-9) -> np.ndarray:
    """L = ∂ln[Bᵉᵗʳ]/∂Bᵉᵗʳ  (eq. 14.101) = 2·∂(½ln)/∂Bᵉᵗʳ."""
    return 2.0 * dHencky_dBe(Be, rtol)


def B_tensor(Be: np.ndarray) -> np.ndarray:
    """B_ijkl = δ_ik (Bᵉᵗʳ)_jl + δ_jk (Bᵉᵗʳ)_il   (eq. 14.102)."""
    d = _I
    return np.einsum("ik,jl->ijkl", d, Be) + np.einsum("jk,il->ijkl", d, Be)


def spatial_tangent_a(D: np.ndarray, Be_tr: np.ndarray,
                      sigma: np.ndarray, J: float) -> np.ndarray:
    """Closed-form spatial consistent tangent, eq. (14.99):
        a_ijkl = (1/2J) [D : L : B]_ijkl − σ_il δ_jk
    D = ∂τ/∂εᵉᵗʳ (small-strain), L = ∂ln Bᵉᵗʳ/∂Bᵉᵗʳ, B from (14.102)."""
    L = log_derivative_L(Be_tr)
    Bt = B_tensor(Be_tr)
    DLB = np.einsum("ijpq,pqrs,rskl->ijkl", D, L, Bt)
    geo = np.einsum("il,jk->ijkl", sigma, _I)          # − σ_il δ_jk
    return DLB / (2.0 * J) - geo


# --------------------------------------------------------------------------- #
#  Full MATISU step with spatial tangent (elastic or J2 inner material)       #
# --------------------------------------------------------------------------- #
def matisu_step(F: np.ndarray, Be_n: np.ndarray, F_n: np.ndarray,
                material: str = "elastic", *, alpha_n: float = 0.0, **mat):
    """One spatial-multiplicative step. Returns a dict with σ, τ, J, spatial
    tangent a, and the updated state (Bᵉ, α).  `material` ∈ {'elastic','j2'}."""
    F_delta = F @ np.linalg.inv(F_n)
    Be_tr = F_delta @ Be_n @ F_delta.T
    eps_tr = hencky_strain(Be_tr)
    if material == "elastic":
        tau = isotropic_elastic_kirchhoff(eps_tr, mat["K"], mat["G"])
        eps_e, alpha, D = eps_tr.copy(), alpha_n, elastic_tangent_D(mat["K"], mat["G"])
        plastic = False
    elif material == "j2":
        tau, eps_e, alpha, D, plastic = j2_return_map(
            eps_tr, alpha_n, mat["K"], mat["G"], mat["sig_y0"], mat["H"])
    else:
        raise ValueError(material)
    J = float(np.linalg.det(F))
    sigma = tau / J
    a = spatial_tangent_a(D, Be_tr, sigma, J)
    Be_new = iso_function(2.0 * eps_e, np.exp)         # Bᵉ = exp[2 εᵉ]
    return dict(sigma=sigma, tau=tau, J=J, a=a, D=D, Be=Be_new, alpha=alpha,
                eps_tr=eps_tr, eps_e=eps_e, plastic=plastic)
