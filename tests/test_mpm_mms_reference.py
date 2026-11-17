"""WP-232i: an independent check of the manufactured solution's derivation (tests/mpm_mms.py).

A centred finite-difference evaluation of the residual  rho0 a - Div_X P - rho0 b  at random
(X, t), through a route that shares nothing with mpm_mms.py's closed forms except u:
  * F from centred differences of u(X, t);
  * the Cauchy stress from the matrix logarithm of V = sqrt(F F^T) (eigen-decomposition), the
    ElasticIsotropic law on that strain (the hypo kernel's integrated feed for a coaxial stretch,
    mpm_mms.py derivation steps 2-3), P = J sigma F^-T;
  * Div_X P by centred differences of P; a by centred differences of u in t.
The residual is relative to the amplitude of the balance's terms. Never against MPM. Not a
binary test (no zone_a binary needed; it is pure numpy).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import mpm_mms as M  # noqa: E402

pytestmark = pytest.mark.zone_a

RHO0, E, NU = 1.7, 3.0, 0.27
LAM, MU = M.lame(E, NU)
OMEGA = 2.0 * math.pi
HX, HT = 1e-5, 1e-4


def _u(X3, t, A, nd):
    """u in 3D embedding; plane strain leaves u_3 = 0."""
    u = np.zeros(3)
    u[:nd] = A * np.sin(2.0 * math.pi * X3[:nd]) * math.sin(OMEGA * t)
    return u


def _F(X3, t, A, nd):
    F = np.eye(3)
    for j in range(3):
        d = np.zeros(3)
        d[j] = HX
        F[:, j] += (_u(X3 + d, t, A, nd) - _u(X3 - d, t, A, nd)) / (2 * HX)
    return F


def _P(X3, t, A, nd):
    F = _F(X3, t, A, nd)
    w, Q = np.linalg.eigh(F @ F.T)
    lnV = Q @ np.diag(0.5 * np.log(w)) @ Q.T
    sig = LAM * np.trace(lnV) * np.eye(3) + 2 * MU * lnV
    return np.linalg.det(F) * sig @ np.linalg.inv(F).T


def _residual(X, t, A, nd, static_lam=None):
    X3 = np.zeros(3)
    X3[:nd] = X
    if static_lam is None:
        Pf = lambda Y: _P(Y, t, A, nd)  # noqa: E731
        acc = (_u(X3, t + HT, A, nd) - 2 * _u(X3, t, A, nd) + _u(X3, t - HT, A, nd)) / HT ** 2
        b = M.body_force(X[None, :], t, A, OMEGA, RHO0, LAM, MU)[0]
    else:
        # quasi-static: u = lam A sin(2 pi X); s(t) = sin(omega t) chosen = lam
        tq = math.asin(static_lam) / OMEGA
        Pf = lambda Y: _P(Y, tq, A, nd)  # noqa: E731
        acc = np.zeros(3)
        b = M.body_force_static(X[None, :], static_lam, A, RHO0, LAM, MU)[0]
    div = np.zeros(3)
    for j in range(3):
        d = np.zeros(3)
        d[j] = HX
        div += (Pf(X3 + d)[:, j] - Pf(X3 - d)[:, j]) / (2 * HX)
    b3 = np.zeros(3)
    b3[:nd] = b
    r = RHO0 * acc - div - RHO0 * b3
    # relative to the amplitude of the balance's terms over the domain (a local scale would
    # divide the finite-difference round-off by a vanishing term near a node of sin)
    scale = RHO0 * OMEGA ** 2 * A + (LAM + 2 * MU) * (2 * math.pi) ** 2 * A
    return np.max(np.abs(r[:nd])) / scale, np.max(np.abs(r[nd:])) if nd < 3 else 0.0


@pytest.mark.parametrize("nd", [2, 3])
@pytest.mark.parametrize("fmax", [1e-3, 0.1])
def test_mms_body_force_balances_momentum(nd, fmax):
    rng = np.random.default_rng(232 + nd)
    A = M.amplitude(fmax)
    worst = 0.0
    for _ in range(60):
        X = rng.uniform(0.0, 1.0, nd)
        t = rng.uniform(0.0, 1.0)
        rel, out_of_plane = _residual(X, t, A, nd)
        worst = max(worst, rel)
        assert out_of_plane < 1e-9 * max(1.0, E)
    assert worst <= 1e-6, worst


@pytest.mark.parametrize("nd", [2, 3])
def test_mms_quasi_static_body_force(nd):
    rng = np.random.default_rng(7 + nd)
    A = M.amplitude(1e-3)
    worst = 0.0
    for _ in range(40):
        X = rng.uniform(0.0, 1.0, nd)
        lam = rng.uniform(0.1, 1.0)
        rel, _ = _residual(X, 0.0, A, nd, static_lam=lam)
        worst = max(worst, rel)
    assert worst <= 1e-6, worst


def test_mms_closed_forms_hold_at_exact_points():
    """u vanishes on every face, and at one period (the gate's t_end)."""
    A = M.amplitude(0.1)
    X = np.array([[0.0, 0.3], [1.0, 0.7], [0.4, 0.0], [0.9, 1.0]])
    u = M.displacement(X, 1.0, A)
    assert abs(u[0, 0]) < 1e-16 and abs(u[1, 0]) < 1e-15
    assert abs(u[2, 1]) < 1e-16 and abs(u[3, 1]) < 1e-15
    assert M.error(X, X, math.sin(OMEGA * 1.0), A) < 1e-14
    # order(): an exact power law has the exact slope
    hs = [1 / 8, 1 / 16, 1 / 32, 1 / 64]
    assert abs(M.order(hs, [3 * h ** 2 for h in hs]) - 2.0) < 1e-12
