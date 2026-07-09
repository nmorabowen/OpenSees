"""Reference (oracle) implementation of the 2D log-strain finite-strain adaptor.

Pure-numpy mirror of `LogStrain2D` (SRC/material/nD/, forthcoming) — the plane
sibling of `LogStrainNDMaterial`. It presents an order-3 PlaneStrain / PlaneStress
face over a 3D small-strain inner material, using the SAME spatial multiplicative
(Hencky) MATISU algorithm as the 3D adaptor:

  * PLANE STRAIN — the 2×2 element deformation gradient is lifted to 3×3 with
    F₃₃ = 1 (and no out-of-plane shear) and driven through the *verified* 3D
    machinery (`logstrain_reference.matisu_step`); the in-plane block of σ and of
    the spatial tangent is read back. The full 3×3 elastic left Cauchy–Green
    state bᵉ is carried, so a plastic ε₃₃ ≠ 0 (deviatoric J2 flow) is tracked
    exactly — plane strain is genuinely the 3D problem with F₃₃ = 1.

  * PLANE STRESS — the thickness stretch λ = F₃₃ is unknown; it is solved from the
    condition σ₃₃(λ) = 0 (a local Newton), then the out-of-plane component is
    statically condensed out of the spatial tangent to give the consistent 2D
    modulus. For a linear-elastic inner the condition is linear, so λ has the
    closed form ε₃₃ = −(K − 2G/3)/(K + 4G/3)·(ε₁₁+ε₂₂) — used to anchor the
    Newton in the tests.

No OpenSees / openseespy dependency — numpy + the 3D oracle only. This module is
the numerical oracle the C++ `LogStrain2D` must match, and it de-risks the one
genuinely new piece of math: the plane-stress F₃₃ condensation.

Conventions (element face): order-3 Voigt σ = [σ₁₁, σ₂₂, σ₁₂]; ENGINEERING shear
in the tangent matrix (γ₁₂ = 2ε₁₂). Fourth-order tensors are (2,2,2,2). Index 2
(zero-based) is the out-of-plane (z) direction throughout.

References (book page → PDF page = book + 25):
  * de Souza Neto, Perić & Owen (2008) Box 14.3 (MATISU), §14.7 (plane stress via
    thickness-stretch condensation), Ch. 9 (small-strain plane-stress condensation).
"""
from __future__ import annotations

import numpy as np

import logstrain_reference as ls3d

# --------------------------------------------------------------------------- #
#  Voigt (order-3) plumbing                                                    #
# --------------------------------------------------------------------------- #
# Engineering-Voigt representative (i,j) pairs for the in-plane order-3 face.
_REP2 = ((0, 0), (1, 1), (0, 1))
_I2 = np.eye(2)
Z = 2  # zero-based out-of-plane (z) index in the lifted 3D tensors


def lift_F(F2: np.ndarray, F33: float) -> np.ndarray:
    """Embed a 2×2 deformation gradient into 3×3 with a pure thickness stretch
    F₃₃ = `F33` and NO out-of-plane shear (F₁₃ = F₃₁ = F₂₃ = F₃₂ = 0)."""
    F3 = np.eye(3)
    F3[:2, :2] = F2
    F3[2, 2] = F33
    return F3


def tensor4_to_voigt3(T: np.ndarray) -> np.ndarray:
    """Map a minor-symmetric (2,2,2,2) tensor to the 3×3 ENGINEERING-Voigt matrix
    D such that σ_I = D_IJ ε_eng_J with ε_eng = [ε₁₁, ε₂₂, γ₁₂=2ε₁₂].

    For an isotropic-consistent minor-symmetric T this needs no shear factor:
    σ₁₁ = T₁₁₁₁ε₁₁ + T₁₁₂₂ε₂₂ + 2T₁₁₁₂ε₁₂ = … + T₁₁₁₂·γ₁₂, etc.
    """
    D = np.zeros((3, 3))
    for I, (i, j) in enumerate(_REP2):
        for J, (k, l) in enumerate(_REP2):
            D[I, J] = T[i, j, k, l]
    return D


# --------------------------------------------------------------------------- #
#  Plane strain — the 2×2 F lifted to 3×3 with F₃₃ = 1                          #
# --------------------------------------------------------------------------- #
def plane_strain_step(F2: np.ndarray, Be_n3: np.ndarray, F_n3: np.ndarray,
                      material: str = "elastic", *, alpha_n: float = 0.0, **mat):
    """One plane-strain MATISU step. Returns a dict with the order-3 in-plane
    Cauchy σ (Voigt), the out-of-plane σ₃₃ (exposed, not condensed), the 2D
    material spatial modulus c2 (2,2,2,2 — the seam the element consumes; element
    forms a2 = c2 − σ_il δ_jk), the total 2D spatial tangent a2 (for checks), and
    the updated 3D state (Bᵉ, α) to be carried to the next step.
    """
    F3 = lift_F(F2, 1.0)
    res = ls3d.matisu_step(F3, Be_n3, F_n3, material, alpha_n=alpha_n, **mat)
    sigma3 = res["sigma"]
    a3 = res["a"]                                   # 3D total spatial tangent (−geo incl.)

    a2 = a3[:2, :2, :2, :2].copy()                  # in-plane block = plane-strain total
    sig2 = sigma3[:2, :2]
    geo2 = np.einsum("il,jk->ijkl", sig2, _I2)      # σ_il δ_jk (in-plane)
    c2 = a2 + geo2                                  # material seam c2 (element re-subtracts geo)

    return dict(
        sigma_voigt=np.array([sigma3[0, 0], sigma3[1, 1], sigma3[0, 1]]),
        sigma_zz=float(sigma3[2, 2]),
        c2=c2, a2=a2, F3=F3,
        Be=res["Be"], alpha=res["alpha"], sigma3=sigma3, a3=a3,
        eps_tr=res["eps_tr"], plastic=res["plastic"],
    )


# --------------------------------------------------------------------------- #
#  Plane stress — thickness stretch λ = F₃₃ solved from σ₃₃ = 0                 #
# --------------------------------------------------------------------------- #
def _sigma33_of_lambda(F2, lam, Be_n3, F_n3, material, alpha_n, mat):
    F3 = lift_F(F2, lam)
    res = ls3d.matisu_step(F3, Be_n3, F_n3, material, alpha_n=alpha_n, **mat)
    return res["sigma"][2, 2], res


def plane_stress_step(F2: np.ndarray, Be_n3: np.ndarray, F_n3: np.ndarray,
                      material: str = "elastic", *, alpha_n: float = 0.0,
                      lam0: float = 1.0, tol: float = 1e-13, maxit: int = 30, **mat):
    """One plane-stress MATISU step. Solves λ = F₃₃ so that σ₃₃ = 0 (local Newton
    with a finite-difference dσ₃₃/dλ — robust for any inner), then condenses the
    out-of-plane direction out of the spatial tangent for the consistent 2D
    modulus. Returns the same dict shape as `plane_strain_step`, plus `F33` and
    `iters`. `Be_n3`/`F_n3` carry the full 3D committed state across steps.
    """
    lam = float(lam0)
    it = 0
    for it in range(1, maxit + 1):
        r, res = _sigma33_of_lambda(F2, lam, Be_n3, F_n3, material, alpha_n, mat)
        if abs(r) <= tol * max(1.0, abs(res["sigma"][0, 0]), abs(res["sigma"][1, 1])):
            break
        h = 1e-8 * max(1.0, abs(lam))
        rp, _ = _sigma33_of_lambda(F2, lam + h, Be_n3, F_n3, material, alpha_n, mat)
        drdl = (rp - r) / h
        if drdl == 0.0:
            raise RuntimeError("plane-stress Newton: singular dσ33/dλ")
        lam -= r / drdl
    else:
        raise RuntimeError(f"plane-stress Newton did not converge (r={r:.3e})")

    F3 = lift_F(F2, lam)
    res = ls3d.matisu_step(F3, Be_n3, F_n3, material, alpha_n=alpha_n, **mat)
    sigma3, a3 = res["sigma"], res["a"]

    # static condensation of the z direction enforcing δσ₃₃ = 0:
    #   δσ_ij = [a_ijkl − a_ij33 a_33kl / a_3333] δd_kl   (i,j,k,l ∈ {0,1})
    denom = a3[Z, Z, Z, Z]
    a2 = (a3[:2, :2, :2, :2]
          - np.einsum("ij,kl->ijkl", a3[:2, :2, Z, Z], a3[Z, Z, :2, :2]) / denom)
    sig2 = sigma3[:2, :2]
    geo2 = np.einsum("il,jk->ijkl", sig2, _I2)
    c2 = a2 + geo2                                  # material seam (element re-subtracts geo)

    return dict(
        sigma_voigt=np.array([sigma3[0, 0], sigma3[1, 1], sigma3[0, 1]]),
        sigma_zz=float(sigma3[2, 2]),
        c2=c2, a2=a2, F3=F3, F33=lam, iters=it,
        Be=res["Be"], alpha=res["alpha"], sigma3=sigma3, a3=a3,
        eps_tr=res["eps_tr"], plastic=res["plastic"],
    )


# --------------------------------------------------------------------------- #
#  Analytic anchors (linear-elastic inner)                                     #
# --------------------------------------------------------------------------- #
def plane_stress_eps33_closedform(eps11: float, eps22: float, K: float, G: float) -> float:
    """Closed-form out-of-plane Hencky strain for a LINEAR-elastic inner under
    plane stress: σ₃₃ = 0 ⇒ ε₃₃ = −(K − 2G/3)/(K + 4G/3)·(ε₁₁ + ε₂₂)
    (= −ν/(1−ν)·(ε₁₁+ε₂₂), the classical thickness contraction)."""
    return -(K - 2.0 * G / 3.0) / (K + 4.0 * G / 3.0) * (eps11 + eps22)


def elastic_D_plane_strain(E: float, nu: float) -> np.ndarray:
    """Textbook small-strain plane-strain elastic matrix (engineering Voigt)."""
    f = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return f * np.array([[1.0 - nu, nu, 0.0],
                         [nu, 1.0 - nu, 0.0],
                         [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0]])


def elastic_D_plane_stress(E: float, nu: float) -> np.ndarray:
    """Textbook small-strain plane-stress elastic matrix (engineering Voigt)."""
    f = E / (1.0 - nu * nu)
    return f * np.array([[1.0, nu, 0.0],
                         [nu, 1.0, 0.0],
                         [0.0, 0.0, (1.0 - nu) / 2.0]])


def KG_to_Enu(K: float, G: float):
    """Bulk/shear → Young/Poisson (3D)."""
    E = 9.0 * K * G / (3.0 * K + G)
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return E, nu
