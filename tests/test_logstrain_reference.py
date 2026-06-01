"""Single-point verification of the log-strain reference kernel (the oracle).

Pure-numpy, no openseespy — runs standalone in CI Zone-A. Validates the spatial
multiplicative (Hencky) MATISU math the C++ `LogStrainNDMaterial` will mirror:

  * the tensor-log derivative against finite differences across ALL eigenvalue-
    degeneracy branches (the D3 trap: uniaxial→2-fold, hydrostatic→3-fold);
  * exp/log inverse consistency;
  * frame indifference and zero stress under rigid rotation;
  * the uniaxial elastic Hencky closed form.

Reference: de Souza Neto, Perić & Owen (2008), Box 14.3 + Appendix A.5/A.6.
"""
import numpy as np
import pytest

import logstrain_reference as lr


# --------------------------------------------------------------------------- #
#  helpers                                                                    #
# --------------------------------------------------------------------------- #
def _rotation():
    """A fixed, reproducible proper rotation (det = +1)."""
    A = np.array([[0.2, 1.0, -0.3],
                  [-0.7, 0.4, 0.5],
                  [0.6, -0.2, 0.9]])
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:           # ensure proper rotation
        Q[:, 0] = -Q[:, 0]
    return Q


def _spd_from(eigs):
    """Build an SPD tensor with the given eigenvalues, rotated off-diagonal."""
    Q = _rotation()
    return Q @ np.diag(eigs) @ Q.T


def _sym_basis():
    """An orthogonal-ish basis of the 6 symmetric 3×3 directions."""
    H = []
    for i in range(3):
        E = np.zeros((3, 3)); E[i, i] = 1.0
        H.append(E)
    for i, j in ((0, 1), (1, 2), (0, 2)):
        E = np.zeros((3, 3)); E[i, j] = E[j, i] = 1.0
        H.append(E)
    return H


def _fd_directional(func, X, H, eps=1e-6):
    """Central finite-difference directional derivative func'(X):H."""
    return (func(X + eps * H) - func(X - eps * H)) / (2.0 * eps)


# --------------------------------------------------------------------------- #
#  D3 — the degeneracy net: analytic dY/dX vs finite difference               #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t0m
@pytest.mark.parametrize("eigs,label", [
    ([1.3, 2.7, 4.1], "distinct"),
    ([2.5, 2.5, 5.0], "double-low"),     # two-fold (uniaxial-like)
    ([5.0, 2.5, 2.5], "double-high"),    # two-fold, singleton first
    ([3.0, 3.0, 3.0], "triple"),         # three-fold (hydrostatic)
])
@pytest.mark.parametrize("y,yp,name", [
    (lr.half_log, lr.half_log_prime, "half_log"),
    (np.exp, np.exp, "exp"),
])
def test_iso_derivative_matches_fd(eigs, label, y, yp, name):
    """The closed-form isotropic-tensor-function derivative (A.52/A.53) must
    equal the finite-difference directional derivative in every multiplicity
    branch. This is the load-bearing D3 check."""
    X = _spd_from(eigs)
    D = lr.iso_function_derivative(X, y, yp)
    for H in _sym_basis():
        analytic = np.einsum("ijkl,kl->ij", D, H)   # D : H
        fd = _fd_directional(lambda Z: lr.iso_function(Z, y), X, H)
        assert np.allclose(analytic, fd, rtol=1e-6, atol=1e-7), (
            f"{name} / {label}: dY/dX:H mismatch\n analytic=\n{analytic}\n fd=\n{fd}"
        )


@pytest.mark.zone_a
@pytest.mark.t0m
def test_derivative_minor_symmetry():
    """dY/dX must be minor-symmetric (ijkl = jikl = ijlk) for symmetric X."""
    X = _spd_from([1.2, 3.4, 0.7])
    D = lr.dHencky_dBe(X)
    assert np.allclose(D, np.swapaxes(D, 0, 1), atol=1e-12)
    assert np.allclose(D, np.swapaxes(D, 2, 3), atol=1e-12)


# --------------------------------------------------------------------------- #
#  exp / log inverse consistency                                              #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t1
@pytest.mark.parametrize("eigs", [[1.3, 2.7, 4.1], [2.5, 2.5, 5.0], [3.0, 3.0, 3.0]])
def test_exp_log_roundtrip(eigs):
    """exp(2·½ln Bᵉ) == Bᵉ  (the predictor's Bᵉ ↔ εᵉ map, eqs 14.30/14.93)."""
    Be = _spd_from(eigs)
    eps = lr.hencky_strain(Be)
    Be_back = lr.iso_function(2.0 * eps, np.exp)
    assert np.allclose(Be_back, Be, rtol=1e-12, atol=1e-12)


# --------------------------------------------------------------------------- #
#  frame indifference / objectivity                                           #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t1
def test_zero_stress_under_rigid_rotation():
    """F = Q (pure rotation) ⇒ Bᵉ = I ⇒ εᵉ = 0 ⇒ σ = 0."""
    Q = _rotation()
    sigma, tau, J, eps = lr.cauchy_from_F_elastic(Q, K=1.5e3, G=7.0e2)
    assert np.allclose(eps, 0.0, atol=1e-12)
    assert np.allclose(sigma, 0.0, atol=1e-9)
    assert J == pytest.approx(1.0, abs=1e-12)


@pytest.mark.zone_a
@pytest.mark.t1
def test_frame_indifference():
    """σ(QF) = Q σ(F) Qᵀ for any rotation Q (material objectivity)."""
    F = np.array([[1.20, 0.10, -0.05],
                  [0.00, 0.90, 0.08],
                  [0.03, -0.04, 1.10]])
    Q = _rotation()
    K, G = 1.5e3, 7.0e2
    s, _, _, _ = lr.cauchy_from_F_elastic(F, K, G)
    sQ, _, _, _ = lr.cauchy_from_F_elastic(Q @ F, K, G)
    assert np.allclose(sQ, Q @ s @ Q.T, rtol=1e-10, atol=1e-9)


# --------------------------------------------------------------------------- #
#  uniaxial elastic Hencky closed form                                        #
# --------------------------------------------------------------------------- #
@pytest.mark.zone_a
@pytest.mark.t1
def test_uniaxial_hencky_principal_stress():
    """Diagonal F = diag(l1,l2,l3): the tensor route must reproduce the
    principal Kirchhoff stress τ_i = K·(Σ ln l) + 2G·(ln l_i − ⅓ Σ ln l),
    and σ = τ / (l1 l2 l3)."""
    l = np.array([1.25, 0.95, 0.90])
    F = np.diag(l)
    K, G = 1.5e3, 7.0e2
    sigma, tau, J, eps = lr.cauchy_from_F_elastic(F, K, G)

    lnl = np.log(l)
    tr = lnl.sum()
    tau_principal = K * tr + 2.0 * G * (lnl - tr / 3.0)
    assert np.allclose(np.diag(tau), tau_principal, rtol=1e-12, atol=1e-10)
    assert np.allclose(np.diag(eps), lnl, rtol=1e-12, atol=1e-12)
    assert J == pytest.approx(float(np.prod(l)), rel=1e-12)
    assert np.allclose(np.diag(sigma), tau_principal / np.prod(l), rtol=1e-12)
    # off-diagonals must vanish for a diagonal F
    assert np.allclose(sigma - np.diag(np.diag(sigma)), 0.0, atol=1e-10)
