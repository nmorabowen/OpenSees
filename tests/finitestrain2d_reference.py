"""Reference (oracle) for the shared 2D finite-strain assembly kernel (ADR 70 P0).

Pure-numpy single-element updated-Lagrangian finite-strain assembly for the 2D
continuum family (T3 / Q4 / T6), driving the LogStrain2D material oracle
(logstrain2d_reference) for σ(F) and the in-plane spatial modulus c₂(F). This is
the oracle the C++ `LadrunoFiniteStrain2DKernel` must match, and it de-risks the
one genuinely new piece of element math: the F-bar (bar-F) consistent tangent's
unsymmetric G₀ coupling.

Per Gauss point (indices i,j,k,l ∈ {0,1}; a,b node indices; engineering shear):
  * F   = I + Σ_a u_a ⊗ ∂N_a/∂X                 (2×2 deformation gradient)
  * F-bar: F̄ = (J₀/J)^{1/2} F  (n=2, NOT ⅓) — drive the material with F̄, keep the
    spatial config (g, dv) on the actual F. (element-level F-bar; a no-op on the
    constant-strain T3 where J₀≡J.)
  * material: setTrialF(F̄) → Cauchy σ (in-plane 2×2) + full c₂[2][2][2][2]
  * spatial gradients g_j^a = Σ_k ∂N_a/∂X_k (F⁻¹)_kj      (push-forward, actual F)
  * residual   f_{a,i} = ∫ σ_ij g_j^a t dA,  dA-weight = J·detJ₀·w
  * tangent    a_ijkl = c_ijkl − σ_il δ_jk,
               K_{(a,i)(b,k)} = ∫ g_j^a a_ijkl g_l^b t dA
               + [F-bar] ∫ p_i^a (g0_m^b − g_m^b) t dA,  p_i^a = Σ_j g_j^a q_ij,
                 q_ij = (1/n) Σ_p a_ijpp − ((n-1)/n) σ_ij   (n=2; the FULL a=c−σδ,
                 NOT c alone — the dSNPO eq-15.10/15.11 unsymmetric G₀ coupling)

The geometric term −σ_il δ_jk is NOT minor-symmetric in (k,l); it must enter the
full-gradient contraction (dSNPO eq 14.99). Correctness is anchored by a
finite-difference self-check: the assembled K must equal ∂f_int/∂u (unsym-aware).

No OpenSees / openseespy dependency — numpy + the material oracle only.
"""
from __future__ import annotations

import numpy as np

import logstrain2d_reference as l2

I2 = np.eye(2)
I3 = np.eye(3)


# --------------------------------------------------------------------------- #
#  Element library: reference nodal coords, GP rule, and (N, ∂N/∂ξ) per point #
# --------------------------------------------------------------------------- #
def _t3_shapes(xi, eta):
    N = np.array([1.0 - xi - eta, xi, eta])
    dN = np.array([[-1.0, 1.0, 0.0],      # ∂N/∂ξ
                   [-1.0, 0.0, 1.0]])     # ∂N/∂η
    return N, dN


def _q4_shapes(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([
        [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],   # ∂N/∂ξ
        [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]])      # ∂N/∂η
    return N, dN


def _t6_shapes(xi, eta):
    # corner 0,1,2 (L0=1-xi-eta, L1=xi, L2=eta), midside 3=(0-1),4=(1-2),5=(2-0)
    L0, L1, L2 = 1.0 - xi - eta, xi, eta
    N = np.array([L0 * (2 * L0 - 1), L1 * (2 * L1 - 1), L2 * (2 * L2 - 1),
                  4 * L0 * L1, 4 * L1 * L2, 4 * L2 * L0])
    # ∂/∂xi (L0_xi=-1, L1_xi=1, L2_xi=0), ∂/∂eta (L0_eta=-1, L1_eta=0, L2_eta=1)
    dNdxi = np.array([1 - 4 * L0, 4 * L1 - 1, 0.0,
                      4 * (L0 - L1), 4 * L2, -4 * L2])
    dNdeta = np.array([1 - 4 * L0, 0.0, 4 * L2 - 1,
                       -4 * L1, 4 * L1, 4 * (L0 - L2)])
    return N, np.array([dNdxi, dNdeta])


# GP rules: list of (xi, eta, w). Triangles integrate over area (Σw = ½).
_GP = {
    "T3": [(1 / 3, 1 / 3, 0.5)],
    "Q4": [(g1, g2, 1.0) for g1 in (-1 / np.sqrt(3), 1 / np.sqrt(3))
           for g2 in (-1 / np.sqrt(3), 1 / np.sqrt(3))],
    "T6": [(2 / 3, 1 / 6, 1 / 6), (1 / 6, 2 / 3, 1 / 6), (1 / 6, 1 / 6, 1 / 6)],
}
_SHAPES = {"T3": _t3_shapes, "Q4": _q4_shapes, "T6": _t6_shapes}
_CENTROID = {"T3": (1 / 3, 1 / 3), "Q4": (0.0, 0.0), "T6": (1 / 3, 1 / 3)}

# Reference (undeformed) nodal coordinates — mildly distorted, positive Jacobian.
REF = {
    "T3": np.array([[0.0, 0.0], [2.0, 0.2], [0.3, 1.6]]),
    "Q4": np.array([[0.0, 0.0], [2.0, 0.1], [2.2, 1.8], [0.1, 1.9]]),
    "T6": np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0],
                    [1.0, -0.05], [1.05, 1.0], [-0.05, 1.0]]),
}


def ref_grads(elem, xi, eta, X):
    """(∂N/∂X [2×nN], detJ0) at natural point (xi,eta) for reference coords X."""
    N, dNdxi = _SHAPES[elem](xi, eta)           # dNdxi: [2 × nN]
    J = dNdxi @ X                               # J[α][i] = ∂X_i/∂ξ_α  (= Jac^T)
    detJ0 = np.linalg.det(J)
    Jinv = np.linalg.inv(J)                      # Jinv[i][α] = ∂ξ_α/∂X_i
    dNdX = Jinv @ dNdxi                          # ∂N_a/∂X_i = Σ_α Jinv[i][α] ∂N_a/∂ξ_α
    return dNdX, detJ0, N


# --------------------------------------------------------------------------- #
#  Assembly                                                                    #
# --------------------------------------------------------------------------- #
def _material(F2, K, G, plane):
    step = l2.plane_stress_step if plane == "stress" else l2.plane_strain_step
    res = step(F2, I3.copy(), I3.copy(), "elastic", K=K, G=G)
    s = res["sigma_voigt"]
    sig = np.array([[s[0], s[2]], [s[2], s[1]]])
    return sig, res["c2"]


def assemble(elem, u, K, G, fbar=False, thickness=1.0, plane="strain", X=None):
    """Return (f_int [2nN], Kmat [2nN×2nN]) for nodal displacement u [2nN]."""
    if X is None:
        X = REF[elem]
    nN = X.shape[0]
    U = u.reshape(nN, 2)
    f = np.zeros(2 * nN)
    Kmat = np.zeros((2 * nN, 2 * nN))

    # centroid deformation gradient + spatial gradients (for F-bar)
    J0 = 1.0
    if fbar:
        c_xi, c_eta = _CENTROID[elem]
        dNdX0, _, _ = ref_grads(elem, c_xi, c_eta, X)
        F0 = I2 + U.T @ dNdX0.T                 # F0_ij = δ_ij + Σ_a U_a,i dNdX_j,a
        J0 = np.linalg.det(F0)
        g0 = dNdX0.T @ np.linalg.inv(F0)        # g0[a][m] = Σ_k dNdX[k][a] F0inv[k][m]

    for (xi, eta, w) in _GP[elem]:
        dNdX, detJ0, _ = ref_grads(elem, xi, eta, X)   # dNdX[i][a]
        F = I2 + U.T @ dNdX.T                          # F[i][j] = δ + Σ_a U_a,i dNdX_j,a
        J = np.linalg.det(F)
        Fmat = F
        if fbar:
            Fmat = (J0 / J) ** 0.5 * F
        sig, c2 = _material(Fmat, K, G, plane)

        Finv = np.linalg.inv(F)
        g = dNdX.T @ Finv                              # g[a][j] = Σ_k dNdX[k][a] Finv[k][j]
        dv = J * detJ0 * w * thickness

        # residual
        for a in range(nN):
            for i in range(2):
                f[2 * a + i] += np.dot(sig[i], g[a]) * dv

        # tangent: a_ijkl = c_ijkl − σ_il δ_jk
        a4 = c2 - np.einsum("il,jk->ijkl", sig, I2)
        for a in range(nN):
            for b in range(nN):
                for i in range(2):
                    for k in range(2):
                        s = 0.0
                        for j in range(2):
                            for l in range(2):
                                s += g[a][j] * a4[i][j][k][l] * g[b][l]
                        Kmat[2 * a + i, 2 * b + k] += s * dv

        if fbar:
            # unsymmetric F-bar G₀ coupling (dSNPO eq 15.10/15.11, n=2):
            #   q_ij = (1/n) Σ_p a4_ijpp − ((n-1)/n) σ_ij   (uses the FULL a4=c−σδ)
            #   K_{ai,bk} += (Σ_j g_j^a q_ij)(g0_k^b − g_k^b) dv
            q = 0.5 * np.einsum("ijpp->ij", a4) - 0.5 * sig
            p = g @ q.T                                 # p[a][i] = Σ_j g[a][j] q_ij
            for a in range(nN):
                for b in range(nN):
                    for i in range(2):
                        for m in range(2):
                            Kmat[2 * a + i, 2 * b + m] += p[a][i] * (g0[b][m] - g[b][m]) * dv

    return f, Kmat


def internal_force(elem, u, K, G, **kw):
    return assemble(elem, u, K, G, **kw)[0]


def fd_tangent(elem, u, K, G, h=1e-7, **kw):
    """∂f_int/∂u by central differences (unsymmetric-aware)."""
    n = u.size
    Kfd = np.zeros((n, n))
    for d in range(n):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (internal_force(elem, up, K, G, **kw)
                     - internal_force(elem, um, K, G, **kw)) / (2 * h)
    return Kfd
