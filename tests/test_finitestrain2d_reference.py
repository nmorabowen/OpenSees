"""Self-checks of the P0 2D finite-strain assembly oracle (finitestrain2d_reference)
against first principles — numpy only, no OpenSees build.

The headline gate is the FINITE-DIFFERENCE consistent-tangent: the analytically
assembled element stiffness K (including the geometric −σδ term and, for F-bar, the
unsymmetric eq-15.10 G₀ coupling) must equal ∂f_int/∂u to FD accuracy, for every
element (T3/Q4/T6) in both std and F-bar formulations. That single check validates
the whole shared-kernel math the C++ must mirror. Plus: homogeneous-F patch
(self-equilibrium), rigid-body / rank at u=0, and the F-bar volumetric-locking cure.
"""
import numpy as np
import pytest

import finitestrain2d_reference as fs

K, G = 1500.0, 700.0
ELEMS = ["T3", "Q4", "T6"]


def _rand_u(elem, scale=0.03, seed=0):
    rng = np.random.default_rng(seed)
    nN = fs.REF[elem].shape[0]
    return scale * rng.standard_normal(2 * nN)


# --------------------------------------------------------------------------- #
#  Headline: analytic K == ∂f_int/∂u  (unsym-aware), std and F-bar            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("elem", ELEMS)
@pytest.mark.parametrize("fbar", [False, True])
def test_consistent_tangent_vs_fd(elem, fbar):
    u = _rand_u(elem, seed=1)
    _, Kan = fs.assemble(elem, u, K, G, fbar=fbar)
    Kfd = fs.fd_tangent(elem, u, K, G, fbar=fbar)
    scale = max(1.0, np.max(np.abs(Kan)))
    np.testing.assert_allclose(Kan, Kfd, rtol=2e-5, atol=1e-5 * scale,
                               err_msg=f"{elem} fbar={fbar}: K != dfint/du")


# --------------------------------------------------------------------------- #
#  Homogeneous-F patch: every GP sees the same F; internal forces equilibrate  #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("elem", ELEMS)
def test_finite_rotation_objectivity(elem):
    # a finite rigid rotation of the undeformed element ⇒ F=R ⇒ σ=0 ⇒ f_int=0.
    # (this is what exposed the Jacobian-transpose gradient bug; FD-consistency
    # and the patch test both missed it.)
    th = 0.4
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    X = fs.REF[elem]
    u = ((R - fs.I2) @ X.T).T.reshape(-1)
    f, _ = fs.assemble(elem, u, K, G)
    assert np.linalg.norm(f) < 1e-8, f"{elem}: not objective (||f(rigid R)||={np.linalg.norm(f):.2e})"


@pytest.mark.parametrize("elem", ELEMS)
def test_homogeneous_patch(elem):
    Ft = np.array([[1.10, 0.05], [0.02, 0.96]])
    X = fs.REF[elem]
    u = ((Ft - fs.I2) @ X.T).T.reshape(-1)     # u_a = (Ft−I) X_a
    f, _ = fs.assemble(elem, u, K, G)
    nN = X.shape[0]
    net = f.reshape(nN, 2).sum(axis=0)         # Σ_a f_a = 0 (self-equilibrium)
    np.testing.assert_allclose(net, [0.0, 0.0], atol=1e-9)


# --------------------------------------------------------------------------- #
#  Reduce-to-linear at u=0: symmetric, 3 rigid-body modes (2D)                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("elem", ELEMS)
def test_zero_disp_stiffness_structure(elem):
    X = fs.REF[elem]
    nN = X.shape[0]
    _, K0 = fs.assemble(elem, np.zeros(2 * nN), K, G)
    Kn = np.linalg.norm(K0)
    # at u=0, σ=0 ⇒ the geometric term vanishes ⇒ K0 is symmetric (minor-sym c at B=I)
    np.testing.assert_allclose(K0, K0.T, atol=1e-7 * Kn,
                               err_msg=f"{elem}: K(u=0) not symmetric")
    # the 3 infinitesimal rigid-body modes (about the centroid) are the nullspace
    Xc = X - X.mean(axis=0)
    modes = [np.tile([1.0, 0.0], nN), np.tile([0.0, 1.0], nN),
             np.column_stack([-Xc[:, 1], Xc[:, 0]]).reshape(-1)]
    for name, r in zip(("tx", "ty", "rot"), modes):
        rel = np.linalg.norm(K0 @ r) / (Kn * np.linalg.norm(r))
        assert rel < 1e-7, f"{elem}: rigid mode {name} not in nullspace (rel={rel:.1e})"


# --------------------------------------------------------------------------- #
#  F-bar cures volumetric locking (near-incompressible); T3 F-bar is a no-op   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("elem", ["Q4", "T6"])
def test_fbar_softens_near_incompressible(elem):
    Kbulk, Gsh = 1.0e6, 700.0                  # ν ≈ 0.4997
    u = _rand_u(elem, scale=0.02, seed=3)
    _, Kstd = fs.assemble(elem, u, Kbulk, Gsh, fbar=False)
    _, Kbar = fs.assemble(elem, u, Kbulk, Gsh, fbar=True)
    e_std = u @ Kstd @ u
    e_bar = u @ Kbar @ u
    assert e_bar < e_std, f"{elem}: F-bar not softer under near-incompressibility"


def test_t3_fbar_is_noop():
    # constant-strain T3 has J₀≡J ⇒ F̄=F ⇒ F-bar identical to std
    u = _rand_u("T3", seed=5)
    f0, K0 = fs.assemble("T3", u, 1.0e6, 700.0, fbar=False)
    f1, K1 = fs.assemble("T3", u, 1.0e6, 700.0, fbar=True)
    np.testing.assert_allclose(f1, f0, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(K1, K0, rtol=1e-10, atol=1e-8)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
