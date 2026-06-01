"""Spatial tangent + J2 plasticity verification for the log-strain oracle.

Pure-numpy, no openseespy — runs standalone in CI Zone-A. Extends the kinematic
oracle with:

  * the consistent SPATIAL tangent a (eq. 14.99) checked against its definition
    (eq. 14.95) by finite-differencing the Kirchhoff stress w.r.t. F — both for
    the elastic Hencky law and for J2 in the plastic regime;
  * the small-strain J2 consistent tangent ∂τ/∂εᵉᵗʳ vs finite difference;
  * exact plastic incompressibility (Remark 14.4);
  * the (c)-vs-(b) discriminator: large simple shear with the log-strain model
    shows NO spurious Jaumann stress oscillation (§14.10.3).

Reference: de Souza Neto, Perić & Owen (2008), §14.4–14.5, Box 14.4, §14.10.3.
"""
import numpy as np
import pytest

import logstrain_reference as lr

_I = np.eye(3)
MAT = dict(K=1.5e3, G=7.0e2)
J2 = dict(K=1.5e3, G=7.0e2, sig_y0=10.0, H=50.0)


def _F_generic():
    return np.array([[1.20, 0.10, -0.05],
                     [0.02, 0.90, 0.08],
                     [0.03, -0.04, 1.10]])


# --------------------------------------------------------------------------- #
#  J2 small-strain consistent tangent  ∂τ/∂εᵉᵗʳ  vs finite difference         #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t0m
@pytest.mark.parametrize("scale,regime", [(1e-3, "elastic"), (3e-2, "plastic")])
def test_j2_consistent_tangent_vs_fd(scale, regime):
    """∂τ/∂εᵉᵗʳ from j2_return_map must match a finite-difference probe in both
    the elastic and the plastic regime (algorithmic tangent, fixed α_n)."""
    # a non-trivial symmetric trial elastic strain at the requested magnitude
    eps = scale * np.array([[2.0, 0.5, -0.3],
                            [0.5, -1.0, 0.4],
                            [-0.3, 0.4, -1.0]])
    alpha_n = 0.0
    tau0, _, _, D, plastic = lr.j2_return_map(eps, alpha_n, **J2)
    assert plastic == (regime == "plastic")

    h = 1e-7
    for i in range(3):
        for j in range(i, 3):                     # symmetric perturbations
            E = np.zeros((3, 3)); E[i, j] = E[j, i] = 1.0
            tp, *_ = lr.j2_return_map(eps + h * E, alpha_n, **J2)
            tm, *_ = lr.j2_return_map(eps - h * E, alpha_n, **J2)
            dtau = (tp - tm) / (2.0 * h)           # = D : E
            analytic = np.einsum("ijkl,kl->ij", D, E)
            assert np.allclose(analytic, dtau, rtol=1e-5, atol=1e-5), (
                f"{regime}: ∂τ/∂ε mismatch at ({i},{j})")


# --------------------------------------------------------------------------- #
#  SPATIAL tangent a (14.99) vs its definition (14.95) by FD of τ(F)          #
# --------------------------------------------------------------------------- #
def _a_by_fd(material, F, Be_n, F_n, alpha_n, h=1e-7, **mat):
    """a_check = (1/J) ∂τ_ij/∂F_kq · F_lq − σ_il δ_jk  (eq. 14.95), τ by FD."""
    res = lr.matisu_step(F, Be_n, F_n, material, alpha_n=alpha_n, **mat)
    J, sigma = res["J"], res["sigma"]
    dtau_dF = np.zeros((3, 3, 3, 3))
    for k in range(3):
        for q in range(3):
            E = np.zeros((3, 3)); E[k, q] = 1.0
            tp = lr.matisu_step(F + h * E, Be_n, F_n, material, alpha_n=alpha_n, **mat)["tau"]
            tm = lr.matisu_step(F - h * E, Be_n, F_n, material, alpha_n=alpha_n, **mat)["tau"]
            dtau_dF[:, :, k, q] = (tp - tm) / (2.0 * h)
    a_mat = np.einsum("ijkq,lq->ijkl", dtau_dF, F) / J
    geo = np.einsum("il,jk->ijkl", sigma, _I)
    return res["a"], a_mat - geo


@pytest.mark.zone_a
@pytest.mark.t0m
def test_spatial_tangent_elastic_vs_fd():
    """Elastic Hencky: analytic spatial tangent (14.99) == FD definition (14.95)."""
    F = _F_generic()
    a_analytic, a_fd = _a_by_fd("elastic", F, _I.copy(), _I.copy(), 0.0, **MAT)
    assert np.allclose(a_analytic, a_fd, rtol=1e-5, atol=1e-4)


@pytest.mark.zone_a
@pytest.mark.t0m
def test_spatial_tangent_j2_plastic_vs_fd():
    """J2 in the plastic regime: analytic spatial tangent == FD definition.
    Take one committed step first so the state (Bᵉ_n, α_n) is non-trivial."""
    F1 = np.diag([1.02, 0.99, 0.99])
    s1 = lr.matisu_step(F1, _I.copy(), _I.copy(), "j2", alpha_n=0.0, **J2)
    Be_n, alpha_n, F_n = s1["Be"], s1["alpha"], F1
    F2 = np.array([[1.10, 0.06, 0.0],
                   [0.04, 0.97, 0.0],
                   [0.0, 0.0, 0.96]])
    s2 = lr.matisu_step(F2, Be_n, F_n, "j2", alpha_n=alpha_n, **J2)
    assert s2["plastic"], "expected a plastic step for this check"
    a_analytic, a_fd = _a_by_fd("j2", F2, Be_n, F_n, alpha_n, **J2)
    assert np.allclose(a_analytic, a_fd, rtol=2e-5, atol=1e-3)


# --------------------------------------------------------------------------- #
#  Exact plastic incompressibility  (Remark 14.4)                             #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t1
def test_plastic_incompressibility():
    """The deviatoric (von-Mises) plastic corrector leaves the volumetric part
    of the elastic strain untouched ⇒ tr(εᵉ) = tr(εᵉᵗʳ) ⇒ det Fᵖ = 1 exactly."""
    eps_tr = np.array([[0.05, 0.02, 0.0],
                       [0.02, -0.02, 0.01],
                       [0.0, 0.01, -0.01]])
    tau, eps_e, alpha, D, plastic = lr.j2_return_map(eps_tr, 0.0, **J2)
    assert plastic
    assert np.trace(eps_e) == pytest.approx(np.trace(eps_tr), abs=1e-14)


# --------------------------------------------------------------------------- #
#  (c) vs (b): large simple shear — NO Jaumann stress oscillation (§14.10.3)  #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t1
def test_simple_shear_no_jaumann_oscillation():
    """Drive large monotonic simple shear with the log-strain J2 model. The
    Kirchhoff shear stress τ₁₂ must rise monotonically to a hardening plateau —
    NOT oscillate (rise-then-fall) as the Jaumann-rate hypoelastic model does."""
    gammas = np.linspace(0.0, 3.0, 121)           # up to very large shear
    Be_n, alpha_n, F_n = _I.copy(), 0.0, _I.copy()
    tau12 = []
    for g in gammas[1:]:
        F = np.array([[1.0, g, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]])
        s = lr.matisu_step(F, Be_n, F_n, "j2", alpha_n=alpha_n, **J2)
        Be_n, alpha_n, F_n = s["Be"], s["alpha"], F
        tau12.append(s["tau"][0, 1])
    tau12 = np.array(tau12)

    # (1) reaches and sustains the plastic regime
    assert tau12.max() > J2["sig_y0"] / np.sqrt(3.0)
    # (2) monotonically non-decreasing (linear hardening) — the discriminator:
    #     a Jaumann model would peak then DECREASE under continued shear.
    assert np.all(np.diff(tau12) > -1e-8), "spurious stress oscillation detected"
    # (3) the maximum is at the very end (sustained, not an interior peak):
    #     a Jaumann oscillation peaks in the interior, not at the largest shear.
    assert np.argmax(tau12) == len(tau12) - 1, "interior peak ⇒ oscillation"
