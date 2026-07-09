"""Self-checks of the 2D log-strain oracle (logstrain2d_reference.py) against
KNOWN PHYSICS — no OpenSees build required, numpy only.

The plane-strain path just restricts the already-verified 3D MATISU machinery, so
these tests concentrate the scrutiny on the genuinely new pieces: the plane-stress
thickness-stretch (F₃₃) Newton solve and the static condensation of the spatial
tangent. Anchors:

  * small-strain limit reproduces the textbook plane-strain / plane-stress elastic
    matrices (validates lift + block extraction + condensation together);
  * the elastic plane-stress λ = F₃₃ matches its closed form (anchors the Newton);
  * plane strain == full 3D with F₃₃ = 1 (identity);
  * finite plane-stress genuinely has σ₃₃ = 0 and σ₁₃ = σ₂₃ = 0;
  * rigid rotation is frame-objective;
  * plastic (J2) plane strain tracks ε₃₃ ≠ 0 (deviatoric out-of-plane flow).
"""
import numpy as np
import pytest

import logstrain2d_reference as l2

K, G = 1500.0, 700.0
E, NU = l2.KG_to_Enu(K, G)
I3 = np.eye(3)


def _fresh():
    return I3.copy(), I3.copy()


# --------------------------------------------------------------------------- #
#  Small-strain limit → textbook elastic matrices                             #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mode,Dref", [
    ("strain", l2.elastic_D_plane_strain(E, NU)),
    ("stress", l2.elastic_D_plane_stress(E, NU)),
])
def test_smallstrain_elastic_matrix(mode, Dref):
    eps = 1e-6
    F2 = np.eye(2) + np.array([[eps, 0.3 * eps], [0.3 * eps, -0.4 * eps]])
    Be_n, F_n = _fresh()
    step = l2.plane_strain_step if mode == "strain" else l2.plane_stress_step
    res = step(F2, Be_n, F_n, "elastic", K=K, G=G)
    D = l2.tensor4_to_voigt3(res["a2"])
    np.testing.assert_allclose(D, Dref, rtol=1e-5, atol=1e-4 * E,
                               err_msg=f"{mode}: small-strain tangent != textbook")


# --------------------------------------------------------------------------- #
#  Plane stress: λ = F₃₃ matches the closed form; σ₃₃ ≈ 0                       #
# --------------------------------------------------------------------------- #
def test_plane_stress_lambda_closed_form():
    F2 = np.array([[1.15, 0.05], [0.03, 0.92]])
    Be_n, F_n = _fresh()
    res = l2.plane_stress_step(F2, Be_n, F_n, "elastic", K=K, G=G)
    # in-plane principal Hencky strains come from the 3D trial strain diagonal
    eps11, eps22 = res["eps_tr"][0, 0], res["eps_tr"][1, 1]
    eps33 = res["eps_tr"][2, 2]
    eps33_cf = l2.plane_stress_eps33_closedform(eps11, eps22, K, G)
    assert abs(eps33 - eps33_cf) < 1e-10, (eps33, eps33_cf)
    assert abs(res["sigma_zz"]) < 1e-9


def test_plane_stress_out_of_plane_shear_zero():
    # block-diagonal F ⇒ block-diagonal Bᵉ ⇒ σ₁₃ = σ₂₃ = 0 exactly
    F2 = np.array([[1.2, 0.15], [-0.1, 0.95]])
    Be_n, F_n = _fresh()
    res = l2.plane_stress_step(F2, Be_n, F_n, "elastic", K=K, G=G)
    s = res["sigma3"]
    assert abs(s[0, 2]) < 1e-12 and abs(s[1, 2]) < 1e-12


# --------------------------------------------------------------------------- #
#  Plane strain == 3D with F₃₃ = 1 (identity)                                  #
# --------------------------------------------------------------------------- #
def test_plane_strain_equals_3d_F33_one():
    F2 = np.array([[1.1, 0.2], [-0.05, 0.9]])
    Be_n, F_n = _fresh()
    res = l2.plane_strain_step(F2, Be_n, F_n, "elastic", K=K, G=G)
    ref = ls_matisu = __import__("logstrain_reference").matisu_step(
        l2.lift_F(F2, 1.0), I3.copy(), I3.copy(), "elastic", K=K, G=G)
    np.testing.assert_allclose(
        res["sigma_voigt"],
        [ref["sigma"][0, 0], ref["sigma"][1, 1], ref["sigma"][0, 1]], rtol=0, atol=1e-14)


# --------------------------------------------------------------------------- #
#  Frame objectivity: superposed rigid rotation                               #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mode", ["strain", "stress"])
def test_rigid_rotation_objectivity(mode):
    F2 = np.array([[1.18, 0.07], [0.02, 0.93]])
    th = 0.6
    Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    step = l2.plane_strain_step if mode == "strain" else l2.plane_stress_step
    Be_n, F_n = _fresh()
    base = step(F2, Be_n, F_n, "elastic", K=K, G=G)
    Be_n, F_n = _fresh()
    rot = step(Q @ F2, Be_n, F_n, "elastic", K=K, G=G)
    # Cauchy σ must transform as Q σ Qᵀ (in-plane block)
    s0 = np.array([[base["sigma_voigt"][0], base["sigma_voigt"][2]],
                   [base["sigma_voigt"][2], base["sigma_voigt"][1]]])
    s1 = np.array([[rot["sigma_voigt"][0], rot["sigma_voigt"][2]],
                   [rot["sigma_voigt"][2], rot["sigma_voigt"][1]]])
    np.testing.assert_allclose(s1, Q @ s0 @ Q.T, rtol=1e-9, atol=1e-8)


# --------------------------------------------------------------------------- #
#  Plastic plane strain tracks a nonzero out-of-plane elastic strain          #
# --------------------------------------------------------------------------- #
def test_plane_strain_j2_tracks_eps33():
    matp = dict(K=K, G=G, sig_y0=5.0, H=50.0)
    F2 = np.array([[1.05, 0.02], [0.0, 0.98]])
    Be_n, F_n = _fresh()
    res = l2.plane_strain_step(F2, Be_n, F_n, "j2", **matp)
    assert res["plastic"], "expected a plastic step for this stretch/yield"
    # deviatoric J2 flow gives out-of-plane plastic strain ⇒ committed Bᵉ₃₃ ≠ 1
    assert abs(res["Be"][2, 2] - 1.0) > 1e-6
    # yet the kinematic F₃₃ stayed 1 (plane strain), so σ₃₃ is generally ≠ 0
    assert abs(res["sigma_zz"]) > 1e-8


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
