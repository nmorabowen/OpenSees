"""ADR-39 P2b-1 ORACLE — node-to-segment (NTS) penalty contact, frictionless.

Oracle-first (the fork discipline): validate the projection + derived normal +
penalty + B-operator + self-equilibrium + the consistent tangent (vs finite
difference) in pure numpy BEFORE transcribing to SRC/domain/contact/
LadrunoContactKernel.h. No OpenSees, no build.

Mechanics (Wriggers, Computational Contact Mechanics, NTS frictionless penalty):
  master segment nodes x_i, bilinear/linear shape N_i(xi,eta), xbar = Σ N_i x_i,
  tangents g_a = ∂xbar/∂xi_a, normal n = (g1×g2)/|g1×g2| (ORIENTED outward),
  gap gN = n·(x_s − xbar); penetration ⇔ gN<0 AND projection in-bounds.
  traction tn = kn·<−gN>+ (Macaulay); slave f_s=−tn n, master f_i=+N_i tn n.
  B (1×ndof) over [u_s, u_1..u_nps] = [ nᵀ | −N_1 nᵀ ... ]; r = Bᵀ tn.
  K_c = kn BᵀB + tn·(geometric/∂n∂u block); curvature O(gN·κ) dropped.

Run: <pythoncore-3.12> proto_p2b_nts.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")   # Windows console is cp1252 by default
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
TOL = 1e-9


# ---------------------------------------------------------------- shape / geom
def shape(nps, xi, eta):
    """N_i and derivatives dN/dxi, dN/deta for tri-3 (area coords) or quad-4."""
    if nps == 4:  # bilinear quad on [-1,1]^2
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dNe = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    elif nps == 3:  # linear tri, (xi,eta) area coords, N = [1-xi-eta, xi, eta]
        N = np.array([1 - xi - eta, xi, eta])
        dNx = np.array([-1.0, 1.0, 0.0])
        dNe = np.array([-1.0, 0.0, 1.0])
    else:
        raise ValueError(nps)
    return N, dNx, dNe


def interp(N, X):
    return N @ X  # xbar = Σ N_i X_i  (X: nps×3)


def tangents(nps, xi, eta, X):
    _, dNx, dNe = shape(nps, xi, eta)
    return dNx @ X, dNe @ X  # g1, g2


def normal_raw(g1, g2):
    nn = np.cross(g1, g2)
    j = np.linalg.norm(nn)
    return nn / j, j


def orient(n, xbar, c_elem):
    """GATE BLOCKER-1: outward normal from master-element centroid, never winding."""
    return n if n @ (xbar - c_elem) > 0 else -n


# ---------------------------------------------------------------- projection
def project(nps, X, xs, tol=1e-12, maxit=10):
    """BOUNDED closest-point Newton on (x_s − xbar)·g_a = 0 (GATE BLOCKER-2).
    returns (xi, eta, status): 0 conv in-bounds, 1 conv out-of-bounds, -1 no proj."""
    if nps == 4:
        xi, eta = 0.0, 0.0
    else:
        xi, eta = 1 / 3, 1 / 3
    for _ in range(maxit):
        N, dNx, dNe = shape(nps, xi, eta)
        xbar = interp(N, X)
        g1, g2 = dNx @ X, dNe @ X
        d = xs - xbar
        R = np.array([d @ g1, d @ g2])
        if np.linalg.norm(R) < tol:
            break
        # K = [[g1·g1, g1·g2],[g2·g1, g2·g2]] (drop the O(d) curvature correction;
        # de-risk = Gauss-Newton metric, SPD for non-degenerate segments)
        K = np.array([[g1 @ g1, g1 @ g2], [g2 @ g1, g2 @ g2]])
        detK = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]
        if abs(detK) < 1e-14 * (np.linalg.norm(g1) * np.linalg.norm(g2) + 1e-300):
            return xi, eta, -1
        dxi = np.linalg.solve(K, R)
        xi += dxi[0]
        eta += dxi[1]
    else:
        return xi, eta, -1  # no convergence in maxit
    # in-bounds test
    if nps == 4:
        inb = (-1 - 1e-9 <= xi <= 1 + 1e-9) and (-1 - 1e-9 <= eta <= 1 + 1e-9)
    else:
        inb = (xi >= -1e-9 and eta >= -1e-9 and xi + eta <= 1 + 1e-9)
    return xi, eta, (0 if inb else 1)


# ---------------------------------------------------------------- force / B / r
def gap_normal_B(nps, X, xs, c_elem):
    """project, derive oriented normal, gap, shape N at proj, and B (1×ndof)."""
    xi, eta, status = project(nps, X, xs)
    N, _, _ = shape(nps, xi, eta)
    xbar = interp(N, X)
    g1, g2 = tangents(nps, xi, eta, X)
    n, _ = normal_raw(g1, g2)
    n = orient(n, xbar, c_elem)
    gN = n @ (xs - xbar)
    ndof = 3 * (1 + nps)
    B = np.zeros(ndof)
    B[0:3] = n
    for i in range(nps):
        B[3 * (1 + i):3 * (2 + i)] = -N[i] * n
    return xi, eta, status, n, gN, N, B


def residual(nps, X, xs, kn, c_elem):
    """r (ndof) = Bᵀ tn, tn = kn<−gN>+ ; zero if not penetrating / out-of-bounds."""
    xi, eta, status, n, gN, N, B = gap_normal_B(nps, X, xs, c_elem)
    ndof = 3 * (1 + nps)
    if status != 0 or gN >= 0.0:
        return np.zeros(ndof), gN, status
    tn = kn * (-gN)            # >0
    return B * tn, gN, status   # r = Bᵀ tn (B is 1×ndof → vector)


def stack(xs, X):
    """DOF vector [u_s(3), u_1..u_nps(3)] base positions, for FD."""
    return np.concatenate([xs, X.reshape(-1)])


def unstack(q, nps):
    return q[0:3].copy(), q[3:].reshape(nps, 3).copy()


def residual_q(q, nps, kn, c_elem):
    xs, X = unstack(q, nps)
    r, _, _ = residual(nps, X, xs, kn, c_elem)
    return r


def tangent_fd(nps, X, xs, kn, c_elem, h=1e-7):
    """Kernel tangent (what addKtToTang assembles) = −∂r/∂q, by central FD.
    OpenSees convention (verified vs shipped P2a: r_s=+tn·n, addKtToTang=+kn·n⊗n):
    getResidual returns the resisting-force-as-unbalance r, and the assembled tangent
    is K = −∂r/∂q. So the full consistent K_c (main kn BᵀB + ∂n/∂u + curvature) = −∂r/∂q."""
    q0 = stack(xs, X)
    ndof = q0.size
    K = np.zeros((ndof, ndof))
    for j in range(ndof):
        qp = q0.copy(); qp[j] += h
        qm = q0.copy(); qm[j] -= h
        K[:, j] = (residual_q(qp, nps, kn, c_elem) - residual_q(qm, nps, kn, c_elem)) / (2 * h)
    return -K   # kernel tangent = −∂r/∂q


def tangent_main(nps, X, xs, kn, c_elem):
    """analytic MAIN term kn BᵀB (the first-order, frozen-ξ̄ tangent)."""
    xi, eta, status, n, gN, N, B = gap_normal_B(nps, X, xs, c_elem)
    if status != 0 or gN >= 0.0:
        return np.zeros((B.size, B.size)), B
    return kn * np.outer(B, B), B


# ================================================================ TEST BATTERY
def _quad(z=0.0, L=1.0):
    """flat master quad in z=z plane, side L centered at origin, CCW (normal +z)."""
    h = L / 2
    return np.array([[-h, -h, z], [h, -h, z], [h, h, z], [-h, h, z]])


def test_interior_penetration():
    """slave pressed onto the quad interior; flat normal +z; gap = −delta exact."""
    X = _quad(z=0.0)
    c_elem = np.array([0, 0, -1.0])     # element below the face → outward = +z
    delta = 0.01
    xs = np.array([0.1, -0.05, -delta])  # penetrated by delta, projects interior
    kn = 1e6
    xi, eta, status, n, gN, N, B = gap_normal_B(4, X, xs, c_elem)
    assert status == 0, f"proj status {status}"
    assert np.allclose(n, [0, 0, 1]), n
    assert abs(gN - (-delta)) < 1e-12, gN
    r, _, _ = residual(4, X, xs, kn, c_elem)
    # getResidual convention (matches shipped P2a: resid = tn·n): slave block = +tn n,
    # i.e. the penetrating slave is pushed back along +n (restoring, never attract).
    tn = kn * delta
    assert np.allclose(r[0:3], tn * n), r[0:3]
    print(f"[interior] gN={gN:.6e}  n={n}  tn={tn:.3e}  slaveF={r[0:3]}  OK")
    return X, xs, kn, c_elem


def test_self_equilibrium():
    """Σ(slave + 4 master nodal forces) = 0 and ΣM = 0 (self-equilibrated)."""
    X = _quad()
    c_elem = np.array([0, 0, -1.0])
    xs = np.array([0.13, 0.07, -0.02])
    kn = 1e6
    r, gN, status = residual(4, X, xs, kn, c_elem)
    fs = r[0:3]
    fmaster = r[3:].reshape(4, 3)
    Fsum = fs + fmaster.sum(axis=0)
    assert np.linalg.norm(Fsum) < 1e-9 * np.linalg.norm(fs), Fsum
    # moment about origin: x_s×f_s + Σ x_i×f_i  (positions at current config)
    _, _, _, n, _, N, _ = gap_normal_B(4, X, xs, c_elem)
    M = np.cross(xs, fs)
    for i in range(4):
        M = M + np.cross(X[i], fmaster[i])
    assert np.linalg.norm(M) < 1e-9 * np.linalg.norm(fs), M
    print(f"[self-eq] |ΣF|={np.linalg.norm(Fsum):.2e}  |ΣM|={np.linalg.norm(M):.2e}  OK")


def test_winding_flip_identical():
    """GATE BLOCKER-1: reverse master node order → IDENTICAL contact force (normal
    re-derived from the element centroid, never trusted from winding)."""
    X = _quad()
    c_elem = np.array([0, 0, -1.0])
    xs = np.array([0.05, -0.03, -0.015])
    kn = 1e6
    r1, _, _ = residual(4, X, xs, kn, c_elem)
    Xr = X[::-1].copy()                 # reversed winding (normal would flip to −z)
    r2, _, _ = residual(4, Xr, xs, kn, c_elem)
    # reorder r2's master blocks back to X's node order for comparison
    r2m = r2[3:].reshape(4, 3)[::-1]
    r2cmp = np.concatenate([r2[0:3], r2m.reshape(-1)])
    assert np.allclose(r1, r2cmp, atol=1e-6 * (abs(r1).max() + 1e-30)), (r1, r2cmp)
    print(f"[winding] slaveF same={np.allclose(r1[0:3], r2[0:3])}  OK (normal re-derived)")


def test_oblique_normal():
    """tilted master quad → off-axis normal; gap along n; force along n."""
    th = np.radians(30)
    Rt = np.array([[1, 0, 0], [0, np.cos(th), -np.sin(th)], [0, np.sin(th), np.cos(th)]])
    X = (_quad() @ Rt.T)
    c_elem = (np.array([0, 0, -1.0]) @ Rt.T)
    n_expected = Rt @ np.array([0, 0, 1.0])
    delta = 0.008
    xs = -delta * n_expected + np.array([0.02, 0.0, 0.0])  # penetrate by delta
    kn = 1e6
    xi, eta, status, n, gN, N, B = gap_normal_B(4, X, xs, c_elem)
    assert status == 0, status
    assert np.allclose(n, n_expected, atol=1e-9), (n, n_expected)
    r, _, _ = residual(4, X, xs, kn, c_elem)
    assert np.allclose(r[0:3], kn * delta * n, atol=1e-6 * kn * delta), r[0:3]
    print(f"[oblique30] n={n}  gN={gN:.3e}  slaveF·n={r[0:3]@n:.3e}≈{kn*delta:.3e}  OK")


def test_out_of_bounds_zero():
    """slave beside the segment (projects out of bounds) → zero force."""
    X = _quad(L=1.0)
    c_elem = np.array([0, 0, -1.0])
    xs = np.array([5.0, 5.0, -0.01])    # far outside, below plane
    kn = 1e6
    r, gN, status = residual(4, X, xs, kn, c_elem)
    assert status != 0 and np.allclose(r, 0), (status, r)
    print(f"[oob] status={status}  |r|={np.linalg.norm(r):.1e}  OK (no garbage)")


def test_tangent_main_vs_fd_slaveblock():
    """For a FIXED master (P2b-1), the SLAVE-block tangent (∂r_s/∂u_s) must match
    the main term kn(n⊗n) since ∂n/∂u_s=0 (n depends on master nodes only).
    The master-coupled blocks differ from main-only (that's the ∂n/∂u the gate wants,
    exercised in P2b-2 with a deformable/rotated master)."""
    X = _quad()
    c_elem = np.array([0, 0, -1.0])
    xs = np.array([0.07, -0.04, -0.02])
    kn = 1e6
    Kfd = tangent_fd(4, X, xs, kn, c_elem)
    Kmain, B = tangent_main(4, X, xs, kn, c_elem)
    _, _, _, n, _, _, _ = gap_normal_B(4, X, xs, c_elem)
    # slave 3×3 block
    Kfd_ss = Kfd[0:3, 0:3]
    Kana_ss = kn * np.outer(n, n)
    err = np.linalg.norm(Kfd_ss - Kana_ss) / (np.linalg.norm(Kana_ss) + 1e-30)
    assert err < 1e-5, f"slave-block tangent err {err:.2e}\nFD=\n{Kfd_ss}\nana=\n{Kana_ss}"
    # how far is the FULL FD tangent from main-only (the ∂n/∂u magnitude)?
    dnu = np.linalg.norm(Kfd - Kmain) / (np.linalg.norm(Kfd) + 1e-30)
    print(f"[tangent] slave-block FD==kn·n⊗n err={err:.2e} | full ∂n∂u rel-weight={dnu:.3f}")
    # symmetry of the FD tangent (frictionless ⇒ symmetric)
    asym = np.linalg.norm(Kfd - Kfd.T) / (np.linalg.norm(Kfd) + 1e-30)
    print(f"[tangent] FD symmetry residual={asym:.2e} (frictionless ⇒ ~0)")
    assert asym < 1e-4, f"tangent not symmetric: {asym:.2e}"


def test_tri3_interior():
    """linear tri-3 master, interior projection, flat normal."""
    X = np.array([[0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]])
    c_elem = np.array([0.3, 0.3, -1.0])
    delta = 0.005
    xs = np.array([0.3, 0.3, -delta])
    kn = 1e6
    xi, eta, status, n, gN, N, B = gap_normal_B(3, X, xs, c_elem)
    assert status == 0, status
    assert np.allclose(n, [0, 0, 1]), n
    assert abs(gN + delta) < 1e-12, gN
    r, _, _ = residual(3, X, xs, kn, c_elem)
    assert np.allclose(r[0:3], kn * delta * n), r[0:3]
    print(f"[tri3] gN={gN:.3e} n={n} slaveF={r[0:3]} OK")


if __name__ == "__main__":
    print("=== ADR-39 P2b-1 NTS penalty oracle ===")
    test_interior_penetration()
    test_self_equilibrium()
    test_winding_flip_identical()
    test_oblique_normal()
    test_out_of_bounds_zero()
    test_tangent_main_vs_fd_slaveblock()
    test_tri3_interior()
    print("ALL P2b-1 ORACLE CHECKS PASSED")
