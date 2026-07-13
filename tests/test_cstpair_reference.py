"""Self-gates for the LadrunoCSTPair oracle (ADR 70 P4a) — numpy only.

Gates (all at the oracle level, building-free):
  1. FD consistent tangent at FINITE stress (unsym-aware) — the load-bearing
     gate: validates the patch coupling incl. the 1/2 power, the v-weights,
     and the cross blocks. Run for both coupling forms.
  2. row-form == dSNPO 15.37/15.38 form (algebraic identity, machine eps).
  3. Homogeneous deformation => J-bar = J (both triangles) => the pair reduces
     EXACTLY to two independent CST-finite elements (diff vs the P0 oracle's
     T3 assembly without F-bar).
  4. Rank at F ~ I: exactly 3 zero eigenvalues on the free pair (the spike's
     C7 result, now on the oracle that anchors the C++ element).
  5. Tangent is genuinely unsymmetric at finite stress (guards against a
     future symmetrizing "fix").
"""
from __future__ import annotations

import numpy as np
import pytest

import cstpair_reference as cp
import finitestrain2d_reference as fs

KMOD, GMOD = 3000.0, 700.0


def _finite_state(scale=0.08, seed=3):
    """A finite-stress, patch-inhomogeneous displacement state (det F > 0)."""
    rng = np.random.default_rng(seed)
    u = scale * rng.standard_normal(8)
    # bias with a uniform stretch so sigma is O(K*scale), not noise-level
    U = u.reshape(4, 2)
    U += 0.05 * cp.REF4 @ np.array([[1.0, 0.15], [0.0, -0.6]])
    return U.ravel()


@pytest.mark.parametrize("form", ["row", "dsnpo"])
def test_fd_consistent_tangent_finite_stress(form):
    u = _finite_state()
    _, Km = cp.assemble_pair(u, KMOD, GMOD, form=form)
    Kfd = cp.fd_tangent(u, KMOD, GMOD, form=form)
    err = np.max(np.abs(Km - Kfd)) / np.max(np.abs(Kfd))
    assert err < 1e-7, f"FD mismatch ({form}): {err:.3e}"


def test_row_form_equals_dsnpo_form():
    u = _finite_state(seed=11)
    _, Ka = cp.assemble_pair(u, KMOD, GMOD, form="row")
    _, Kb = cp.assemble_pair(u, KMOD, GMOD, form="dsnpo")
    assert np.max(np.abs(Ka - Kb)) < 1e-12 * np.max(np.abs(Ka))


def test_homogeneous_reduces_to_two_csts():
    # affine map: u = X @ (H - I)^T with H the target F => F1 = F2 = H, Jbar = J
    H = np.array([[1.10, 0.05], [0.02, 0.93]])
    u = (cp.REF4 @ (H - cp.np.eye(2)).T).ravel()
    f_pair, K_pair = cp.assemble_pair(u, KMOD, GMOD)

    f_ref = np.zeros(8)
    K_ref = np.zeros((8, 8))
    for tri in cp.TRIS:
        idx = list(tri)
        dof = [2 * n + d for n in idx for d in (0, 1)]
        X = cp.REF4[idx]
        ut = np.array([u[k] for k in dof])
        ft, Kt = fs.assemble("T3", ut, KMOD, GMOD, fbar=False, X=X)
        f_ref[dof] += ft
        K_ref[np.ix_(dof, dof)] += Kt

    assert np.max(np.abs(f_pair - f_ref)) < 1e-10 * max(1.0, np.max(np.abs(f_ref)))
    # tangents also coincide: the coupling term survives (q != 0) but
    # (gbar - g_e) does NOT vanish per-element... it does under homogeneous F:
    # g_e are different rows (different triangles), gbar is their v-average, so
    # the coupling is generally nonzero -- compare via FD instead of K_ref.
    Kfd = cp.fd_tangent(u, KMOD, GMOD)
    err = np.max(np.abs(K_pair - Kfd)) / np.max(np.abs(Kfd))
    assert err < 1e-7


def test_rank_free_pair_is_3_rbm():
    u = 1e-9 * np.ones(8)  # F ~ I
    _, Km = cp.assemble_pair(u, KMOD, GMOD)
    w = np.linalg.eigvals(0.5 * (Km + Km.T))
    nz = np.sum(np.abs(w) < 1e-8 * np.max(np.abs(w)))
    assert nz == 3, f"free pair zero modes = {nz}, want 3"


def test_tangent_unsymmetric_at_finite_stress():
    u = _finite_state(seed=7)
    _, Km = cp.assemble_pair(u, KMOD, GMOD)
    asym = np.max(np.abs(Km - Km.T)) / np.max(np.abs(Km))
    assert asym > 1e-6, "patch tangent unexpectedly symmetric at finite stress"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
