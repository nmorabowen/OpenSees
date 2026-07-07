"""LadrunoComplexEigen (ADR 46, tag 33019 RESERVED) — P0 kernel gate.

P0 ships ONLY the reduced-pencil QZ kernel behind the debug entry
`complexEigen -qz p Mt Ct Kt [-tol eps]`: solve the small quadratic
eigenproblem (lambda^2 Mt + lambda Ct + Kt) z = 0 via the state-space
symmetric-block pencil + LAPACK dggev, and recover (omega0, omegaD, zeta,
lambda, kind) per physical mode. The domain-coupled projection (Rayleigh
factors / assembled C from getDamp()) is P1/P2.

Gate (ADR 46 §7 P0): eigenvalues match an INDEPENDENT companion-form
numpy.linalg.eig oracle to 1e-10 on hand-built and seeded-random pencils,
plus the analytic SDOF under/over-damped and rigid classification branches.

  test_sdof_underdamped        analytic lambda = -zeta*w +- i*w*sqrt(1-zeta^2)
  test_sdof_overdamped         two real roots, kind=overdamped, zeta -> 1
  test_rigid_plus_decay        k=0: lambda=0 (rigid) + lambda=-c/m (overdamped)
  test_classical_rayleigh      Ct = aM*I + bK*diag(w^2): zeta_k analytic
  test_2dof_nonclassical       damper on ONE dof: both pairs vs oracle
  test_random_pencils_vs_numpy seeded p=3/8 pencils vs companion-form eig
  test_bad_input_rejected      wrong entry count / unknown option -> error
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.complex_eigen_ref import (
    expand_cpp_modes,
    modes_table,
    sorted_roots,
)

pytestmark = [pytest.mark.zone_a]


def qz(M, C, K, tol=None):
    """Call the C++ kernel; returns the flat -qz output list."""
    M = np.atleast_2d(np.asarray(M, dtype=float))
    C = np.atleast_2d(np.asarray(C, dtype=float))
    K = np.atleast_2d(np.asarray(K, dtype=float))
    p = M.shape[0]
    ops.wipe()
    args = ["-qz", p]
    args += [float(x) for x in M.flatten()]
    args += [float(x) for x in C.flatten()]
    args += [float(x) for x in K.flatten()]
    if tol is not None:
        args += ["-tol", float(tol)]
    out = ops.complexEigen(*args)
    if isinstance(out, float):
        out = (out,)
    return out


def assert_roots_match(flat, M, C, K, rtol=1e-10):
    got = expand_cpp_modes(flat)
    ref = sorted_roots(M, C, K)
    assert len(got) == len(ref)
    scale = max(np.max(np.abs(ref)), 1.0)
    np.testing.assert_allclose(got, ref, rtol=0.0, atol=rtol * scale)


def assert_residuals_small(flat, M, C, K, rtol=1e-9):
    """Every reported mode must satisfy its OWN quadratic residual — this
    ties the reported lambda to the reported eigenvector (catches
    conjugate-pairing mistakes the eigenvalue-only comparison cannot see)."""
    K = np.atleast_2d(np.asarray(K, dtype=float))
    scale = max(np.linalg.norm(K, ord=np.inf), 1.0)
    for m in modes_table(flat):
        assert m["resid"] <= rtol * scale, \
            f"mode lambda={m['lambda']}: resid {m['resid']} > {rtol * scale}"


# ------------------------------------------------------------------ analytic
def test_sdof_underdamped():
    w = 2.0 * math.pi * 1.7
    zeta = 0.05
    flat = qz([[1.0]], [[2.0 * zeta * w]], [[w * w]])
    modes = modes_table(flat)
    assert len(modes) == 1
    m = modes[0]
    assert m["kind"] == "underdamped"
    assert m["omega0"] == pytest.approx(w, rel=1e-10)
    assert m["zeta"] == pytest.approx(zeta, rel=1e-10)
    assert m["omegaD"] == pytest.approx(w * math.sqrt(1 - zeta**2), rel=1e-10)
    assert m["lambda"].real == pytest.approx(-zeta * w, rel=1e-10)


def test_sdof_overdamped():
    w = 3.0
    zeta = 1.5
    flat = qz([[1.0]], [[2.0 * zeta * w]], [[w * w]])
    modes = modes_table(flat)
    assert len(modes) == 2
    lam_exact = sorted([-zeta * w + w * math.sqrt(zeta**2 - 1),
                        -zeta * w - w * math.sqrt(zeta**2 - 1)])
    lam_got = sorted(m["lambda"].real for m in modes)
    for got, exact in zip(lam_got, lam_exact):
        assert got == pytest.approx(exact, rel=1e-10)
    for m in modes:
        assert m["kind"] == "overdamped"
        assert m["omegaD"] == 0.0
        assert m["zeta"] == pytest.approx(1.0, rel=1e-12)  # -Re/|lambda|


def test_rigid_plus_decay():
    c = 0.3
    flat = qz([[1.0]], [[c]], [[0.0]])
    modes = modes_table(flat)
    assert len(modes) == 2
    kinds = sorted(m["kind"] for m in modes)
    assert kinds == ["overdamped", "rigid"]
    rigid = next(m for m in modes if m["kind"] == "rigid")
    decay = next(m for m in modes if m["kind"] == "overdamped")
    assert rigid["omega0"] == 0.0
    assert abs(rigid["lambda"]) <= 1e-10 * c
    assert decay["lambda"].real == pytest.approx(-c, rel=1e-10)


# ----------------------------------------------------- classical / oracle
def test_classical_rayleigh():
    """Rayleigh-only reduced matrices are the classical case: every mode
    stays 'real' (a plain damped oscillator per mode) with the textbook
    zeta_k = aM/(2 w_k) + bK w_k / 2 — the P1 Route-A oracle, checked here
    directly on the kernel."""
    w = np.array([2.1, 5.3, 9.8, 17.2])
    aM, bK = 0.35, 0.004
    Mt = np.eye(4)
    Kt = np.diag(w**2)
    Ct = aM * Mt + bK * Kt
    flat = qz(Mt, Ct, Kt)
    modes = modes_table(flat)
    assert len(modes) == 4
    for m, wk in zip(modes, w):
        zk = aM / (2 * wk) + bK * wk / 2
        assert m["kind"] == "underdamped"
        assert m["omega0"] == pytest.approx(wk, rel=1e-10)
        assert m["zeta"] == pytest.approx(zk, rel=1e-10)
    assert_roots_match(flat, Mt, Ct, Kt)
    assert_residuals_small(flat, Mt, Ct, Kt)


def test_2dof_nonclassical():
    """Two-mass chain with a dashpot on ONE dof — the canonical
    non-classical problem (C not simultaneously diagonalizable)."""
    M = np.eye(2)
    K = np.array([[5.0, -1.0], [-1.0, 1.0]])
    C = np.array([[0.0, 0.0], [0.0, 0.4]])
    flat = qz(M, C, K)
    modes = modes_table(flat)
    assert len(modes) == 2
    assert all(m["kind"] == "underdamped" for m in modes)
    # non-classical: the two damping ratios are NOT a Rayleigh pair's
    # signature; just require both positive and distinct
    z1, z2 = modes[0]["zeta"], modes[1]["zeta"]
    assert z1 > 0 and z2 > 0 and abs(z1 - z2) > 1e-4
    assert_roots_match(flat, M, C, K)
    assert_residuals_small(flat, M, C, K)


def test_random_pencils_vs_numpy():
    rng = np.random.default_rng(4600)
    for p in (3, 8):
        for _ in range(3):
            w = np.sort(rng.uniform(1.0, 30.0, size=p))
            Mt = np.eye(p)
            Kt = np.diag(w**2)
            # symmetric PSD damping with off-diagonal (non-classical) coupling
            R = rng.uniform(-1.0, 1.0, size=(p, p))
            Ct = 0.05 * (R @ R.T) + 0.02 * np.diag(w)
            flat = qz(Mt, Ct, Kt)
            assert_roots_match(flat, Mt, Ct, Kt)
            # the residual assertion ties each reported lambda to ITS
            # eigenvector — the check that catches beta<0 conjugate
            # mis-pairing (adversarial-review finding M1)
            assert_residuals_small(flat, Mt, Ct, Kt)
            # Ct is PSD by construction, so all modes decay
            for m in modes_table(flat):
                if m["kind"] == "underdamped":
                    assert -1e-12 < m["zeta"] < 1.0


# ------------------------------------------------------------------- guards
def test_bad_input_rejected():
    with pytest.raises(Exception):
        ops.complexEigen("-qz", 2, 1.0, 2.0, 3.0)  # far too few entries
    with pytest.raises(Exception):
        ops.complexEigen("-notAMode", 1, 1.0, 1.0, 1.0)
