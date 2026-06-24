"""ADR-39 P2b-2c / capstone B3 ORACLE — consistent ∂n/∂u normal contact tangent
(the geometric/curvature block deferred since P2b) + a Hertz reference.

Oracle-first (the fork discipline): derive + FD-validate the FULL frictionless NTS
penalty tangent — the shipped main term kn·BᵀB PLUS the geometric block kn·gN·∂B/∂u
(= the ∂n/∂u + ∂N_i/∂u curvature/large-sliding terms) — in pure numpy BEFORE
transcribing the closed form to SRC/analysis/handler/LadrunoContactFE.cpp. No
OpenSees, no build.

WHAT B3 ADDS (over proto_p2b_nts.py, which shipped kn·BᵀB only):
  The residual is r = Bᵀ·tn, tn = kn·<−gN>₊, B = ∂gN/∂u = [ nᵀ | −N_i nᵀ ].
  The EXACT (consistent) tangent is K = −∂r/∂u = kn·(B⊗B + gN·H), H = ∂B/∂u = ∂²gN/∂u².
  The shipped code assembles only kn·B⊗B (EXACT for a flat/fixed master where n is
  constant ⇒ H=0). B3 adds kn·gN·H — the geometric block — so implicit Newton is
  QUADRATIC on CURVED / large-sliding interfaces where n varies with the config.

DERIVATION (Wriggers, Computational Contact Mechanics §6.3; Laursen).
  Master surface x̄(ξ) = Σ N_i(ξ) x_i, tangents a_α = x̄_{,α}, normal n = (a1×a2)/|a1×a2|,
  metric m_αβ = a_α·a_β, contravariant a^α = m^{αβ} a_β, second fundamental form
  κ_αβ = n·a_{α,β} (a_{α,β} = Σ N_{i,αβ} x_i; tri-3 ⇒ 0, quad-4 ⇒ only the ξη cross
  term, and even that is 0 for a PLANAR quad ⇒ κ=0 ⇔ flat facet). gN = n·(x_s − x̄).
  The closest point ξ̄(u) solves Φ_β = (x_s−x̄)·a_β = 0 ⇒ by the envelope theorem
  ∂gN/∂ξ = 0, so the FIRST derivative B = [n | −N_i n] has no ∂ξ̄/∂u term (shipped).
  The SECOND derivative does. Implicit differentiation of Φ gives
      ∂ξ̄^α/∂u_l = H2^{αβ} (Q_β[l] + gN·P_β[l]),     H2_αβ = m_αβ − gN·κ_αβ
  with, per DOF l (a (node,component) pair) and the unit dir e_d,
      slave  d : Q_β = a_β[d],            P_β = 0
      master j,d: Q_β = −N_j a_β[d],      P_β = N_{j,β} n[d]
  and the TOTAL derivative of n (using δn = −(n·δa_α) a^α, n·δa_α = P_α + κ_αβ ∂ξ̄^β):
      ∂n/∂u_l = −[ P_γ[l] + κ_γα ∂ξ̄^α/∂u_l ] a^γ
      ∂N_i/∂u_l = N_{i,α} ∂ξ̄^α/∂u_l
  Then ∂B/∂u_l (a full ndof vector): slave rows = ∂n/∂u_l ; master-i rows =
      −(∂N_i/∂u_l) n − N_i (∂n/∂u_l). H[:,l] = ∂B/∂u_l. K = kn(B⊗B + gN·H).
  H is SYMMETRIC (Hessian of the scalar gN) ⇒ the frictionless geometric tangent is
  SYMMETRIC and solver-safe (the friction Csl non-symmetry is a SEPARATE, P3.5 story).

BYTE-IDENTITY CONTRACT (capstone B3 gate 3): for a flat facet κ=0 ⇒ ∂n/∂u_slave=0 ⇒
  the SLAVE block of H is EXACTLY 0 ⇒ for a FIXED master (only the slave DOFs free) the
  geometric block vanishes and K = kn·n⊗n, bit-for-bit the shipped tangent.

Run: <pythoncore-3.12> proto_b3_normal_tangent.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)


# ---------------------------------------------------------------- shape / geom
def shape(nps, xi, eta):
    """N_i, first derivs dNxi/dNeta, and second derivs (Nxx, Nee, Nxe) for tri-3
    (area coords, all 2nd derivs 0) or quad-4 (only the ξη cross deriv nonzero)."""
    if nps == 4:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dNe = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
        Nxx = np.zeros(4)
        Nee = np.zeros(4)
        Nxe = 0.25 * np.array([1.0, -1.0, 1.0, -1.0])  # ∂²N/∂ξ∂η
    elif nps == 3:
        N = np.array([1 - xi - eta, xi, eta])
        dNx = np.array([-1.0, 1.0, 0.0])
        dNe = np.array([-1.0, 0.0, 1.0])
        Nxx = np.zeros(3)
        Nee = np.zeros(3)
        Nxe = np.zeros(3)
    else:
        raise ValueError(nps)
    return N, dNx, dNe, Nxx, Nee, Nxe


def project(nps, X, xs, tol=1e-13, maxit=20):
    """BOUNDED closest-point Newton on (x_s − xbar)·g_a = 0 (Gauss-Newton metric).
    returns (xi, eta, status): 0 conv in-bounds, 1 conv out-of-bounds, -1 no proj."""
    xi, eta = (0.0, 0.0) if nps == 4 else (1 / 3, 1 / 3)
    for _ in range(maxit):
        N, dNx, dNe, *_ = shape(nps, xi, eta)
        xbar = N @ X
        g1, g2 = dNx @ X, dNe @ X
        d = xs - xbar
        R = np.array([d @ g1, d @ g2])
        if np.linalg.norm(R) < tol:
            break
        K = np.array([[g1 @ g1, g1 @ g2], [g2 @ g1, g2 @ g2]])
        detK = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]
        if abs(detK) < 1e-14 * (np.linalg.norm(g1) * np.linalg.norm(g2) + 1e-300):
            return xi, eta, -1
        dxi = np.linalg.solve(K, R)
        xi += dxi[0]
        eta += dxi[1]
    else:
        return xi, eta, -1
    if nps == 4:
        inb = (-1 - 1e-9 <= xi <= 1 + 1e-9) and (-1 - 1e-9 <= eta <= 1 + 1e-9)
    else:
        inb = (xi >= -1e-9 and eta >= -1e-9 and xi + eta <= 1 + 1e-9)
    return xi, eta, (0 if inb else 1)


def orient_sign(n, xbar, refdir):
    """oriented so n·refdir > 0 (refdir toward the slave's allowed half-space)."""
    return 1.0 if n @ refdir > 0 else -1.0


def gap_normal_B(nps, X, xs, refdir):
    """project → oriented n, gap gN, shape N at proj, B (1×ndof)."""
    xi, eta, status = project(nps, X, xs)
    N, dNx, dNe, *_ = shape(nps, xi, eta)
    xbar = N @ X
    g1, g2 = dNx @ X, dNe @ X
    raw = np.cross(g1, g2)
    n = raw / np.linalg.norm(raw)
    n = orient_sign(n, xbar, refdir) * n
    gN = n @ (xs - xbar)
    ndof = 3 * (1 + nps)
    B = np.zeros(ndof)
    B[0:3] = n
    for i in range(nps):
        B[3 * (1 + i):3 * (2 + i)] = -N[i] * n
    return xi, eta, status, n, gN, N, B


# ---------------------------------------------------------------- residual / FD
def residual(nps, X, xs, kn, refdir):
    xi, eta, status, n, gN, N, B = gap_normal_B(nps, X, xs, refdir)
    ndof = 3 * (1 + nps)
    if status != 0 or gN >= 0.0:
        return np.zeros(ndof)
    tn = kn * (-gN)
    return B * tn


def stack(xs, X):
    return np.concatenate([xs, X.reshape(-1)])


def unstack(q, nps):
    return q[0:3].copy(), q[3:].reshape(nps, 3).copy()


def tangent_fd(nps, X, xs, kn, refdir, h=1e-7):
    """ground-truth consistent tangent K = −∂r/∂q by central FD (re-projects every
    perturbation ⇒ captures the full ∂ξ̄/∂u, ∂n/∂u, ∂N_i/∂u dependence)."""
    q0 = stack(xs, X)
    ndof = q0.size
    K = np.zeros((ndof, ndof))
    for j in range(ndof):
        qp = q0.copy(); qp[j] += h
        qm = q0.copy(); qm[j] -= h
        xs_p, X_p = unstack(qp, nps)
        xs_m, X_m = unstack(qm, nps)
        rp = residual(nps, X_p, xs_p, kn, refdir)
        rm = residual(nps, X_m, xs_m, kn, refdir)
        K[:, j] = (rp - rm) / (2 * h)
    return -K


# ---------------------------------------------------- ANALYTIC consistent tangent
def geom_blocks(nps, X, xs, refdir):
    """kinematics at the projection: n, gN, N, B, plus the curvature data the
    geometric block needs (a_α, metric, contravariant a^α, κ_αβ)."""
    xi, eta, status, n, gN, N, B = gap_normal_B(nps, X, xs, refdir)
    _, dNx, dNe, Nxx, Nee, Nxe = shape(nps, xi, eta)
    a1, a2 = dNx @ X, dNe @ X
    a = [a1, a2]
    m = np.array([[a1 @ a1, a1 @ a2], [a2 @ a1, a2 @ a2]])
    minv = np.linalg.inv(m)
    # contravariant tangents a^α = m^{αβ} a_β
    aup = [minv[0, 0] * a1 + minv[0, 1] * a2, minv[1, 0] * a1 + minv[1, 1] * a2]
    # second derivatives of position a_{αβ} = Σ N_{i,αβ} x_i (tri ⇒ 0; quad ⇒ only ξη)
    a_xx = Nxx @ X
    a_ee = Nee @ X
    a_xe = Nxe @ X
    add = [[a_xx, a_xe], [a_xe, a_ee]]
    kappa = np.array([[n @ add[ap][bp] for bp in range(2)] for ap in range(2)])
    dN = [dNx, dNe]  # N_{i,α}
    return dict(xi=xi, eta=eta, status=status, n=n, gN=gN, N=N, B=B,
                a=a, m=m, minv=minv, aup=aup, kappa=kappa, dN=dN)


def tangent_analytic(nps, X, xs, kn, refdir, return_H=False):
    """K = kn(B⊗B + gN·H), H = ∂B/∂u the geometric Hessian (closed form)."""
    gb = geom_blocks(nps, X, xs, refdir)
    n, gN, N, B = gb["n"], gb["gN"], gb["N"], gb["B"]
    ndof = 3 * (1 + nps)
    if gb["status"] != 0 or gN >= 0.0:
        Z = np.zeros((ndof, ndof))
        return (Z, np.zeros((ndof, ndof))) if return_H else Z
    a, aup, kappa, dN = gb["a"], gb["aup"], gb["kappa"], gb["dN"]
    # H2_αβ = m_αβ − gN κ_αβ  (the projection-sensitivity metric)
    H2 = gb["m"] - gN * kappa
    H2inv = np.linalg.inv(H2)

    # per-DOF Q_β[l], P_β[l]  (l indexes the ndof DOFs: slave xyz, then master nodes)
    Q = np.zeros((ndof, 2))
    P = np.zeros((ndof, 2))
    for d in range(3):                       # slave DOF d
        for b in range(2):
            Q[d, b] = a[b][d]                # a_β·e_d
    for j in range(nps):                     # master node j
        for d in range(3):
            l = 3 * (1 + j) + d
            for b in range(2):
                Q[l, b] = -N[j] * a[b][d]
                P[l, b] = dN[b][j] * n[d]    # N_{j,β} n[d]

    # ∂ξ̄^α/∂u_l = H2^{αβ}(Q_β + gN P_β)
    zeta = np.zeros((ndof, 2))
    for l in range(ndof):
        rhs = Q[l] + gN * P[l]
        zeta[l] = H2inv @ rhs

    # ∂n/∂u_l = −Σ_γ (P_γ[l] + Σ_α κ_γα ζ^α[l]) a^γ
    dn = np.zeros((ndof, 3))
    for l in range(ndof):
        coef = np.zeros(2)
        for g in range(2):
            coef[g] = P[l, g] + kappa[g, 0] * zeta[l, 0] + kappa[g, 1] * zeta[l, 1]
        dn[l] = -(coef[0] * aup[0] + coef[1] * aup[1])

    # ∂N_i/∂u_l = N_{i,α} ζ^α[l]
    dNi = np.zeros((ndof, nps))
    for l in range(ndof):
        for i in range(nps):
            dNi[l, i] = dN[0][i] * zeta[l, 0] + dN[1][i] * zeta[l, 1]

    # H[:,l] = ∂B/∂u_l : slave rows = ∂n/∂u_l ; master-i rows = −∂N_i/∂u_l·n − N_i·∂n/∂u_l
    H = np.zeros((ndof, ndof))
    for l in range(ndof):
        col = np.zeros(ndof)
        col[0:3] = dn[l]
        for i in range(nps):
            col[3 * (1 + i):3 * (2 + i)] = -dNi[l, i] * n - N[i] * dn[l]
        H[:, l] = col

    K = kn * (np.outer(B, B) + gN * H)
    return (K, H) if return_H else K


# ================================================================ TEST BATTERY
def _quad(z=0.0, L=1.0):
    h = L / 2
    return np.array([[-h, -h, z], [h, -h, z], [h, h, z], [-h, h, z]])


def _warp_quad(w, L=1.0):
    """warped (NON-planar) quad: nodes 0,2 lifted +w, nodes 1,3 dropped −w ⇒ a genuine
    saddle (κ_ξη ≠ 0). The canonical curved-facet test for the ∂n/∂u block."""
    h = L / 2
    return np.array([[-h, -h, +w], [h, -h, -w], [h, h, +w], [-h, h, -w]])


def test_main_recovers_shipped_flat():
    """flat fixed master ⇒ geometric block H_slave = 0 EXACTLY (κ=0); the full tangent
    SLAVE block == kn·n⊗n bit-for-bit ⇒ byte-identical to the shipped kn·BᵀB (gate 3)."""
    X = _quad()
    refdir = np.array([0, 0, 1.0])
    xs = np.array([0.07, -0.04, -0.02])
    kn = 1e6
    K, H = tangent_analytic(4, X, xs, kn, refdir, return_H=True)
    _, _, _, n, _, _, _ = gap_normal_B(4, X, xs, refdir)
    Hss = H[0:3, 0:3]
    assert np.linalg.norm(Hss) < 1e-12, f"flat-master slave Hessian must be 0, got {Hss}"
    Kss = K[0:3, 0:3]
    assert np.allclose(Kss, kn * np.outer(n, n), atol=1e-6 * kn), Kss
    print(f"[flat byte-id] |H_ss|={np.linalg.norm(Hss):.1e} (==0) ; K_ss==kn·n⊗n OK")


def test_analytic_vs_fd_flat():
    """flat master, FULL tangent (incl. master-coupled ∂n/∂u from facet tilt): analytic
    == FD even though the master nodes are free (a flat facet still tilts when a single
    master node moves ⇒ ∂n/∂u_master ≠ 0)."""
    X = _quad()
    refdir = np.array([0, 0, 1.0])
    xs = np.array([0.07, -0.04, -0.02])
    kn = 1e6
    Kfd = tangent_fd(4, X, xs, kn, refdir)
    Kan = tangent_analytic(4, X, xs, kn, refdir)
    err = np.linalg.norm(Kfd - Kan) / (np.linalg.norm(Kfd) + 1e-30)
    assert err < 1e-5, f"flat full tangent err {err:.2e}\nFD=\n{Kfd}\nAN=\n{Kan}"
    print(f"[flat full] analytic vs FD rel-err={err:.2e} OK")


def test_analytic_vs_fd_curved_quad():
    """THE gate: a WARPED (curved) quad master, slave pressed in ⇒ κ≠0 ⇒ the slave
    Hessian is NONZERO; analytic consistent tangent == FD to ~1e-6."""
    X = _warp_quad(w=0.15)
    refdir = np.array([0, 0, 1.0])
    # slave above the saddle center, pressed slightly through
    xs = np.array([0.05, -0.03, -0.02])
    kn = 1e6
    gb = geom_blocks(4, X, xs, refdir)
    assert abs(gb["kappa"][0, 1]) > 1e-3, f"warped quad must have κ_ξη≠0: {gb['kappa']}"
    K, H = tangent_analytic(4, X, xs, kn, refdir, return_H=True)
    assert np.linalg.norm(H[0:3, 0:3]) > 1e-3, "curved master ⇒ nonzero slave Hessian"
    Kfd = tangent_fd(4, X, xs, kn, refdir)
    err = np.linalg.norm(Kfd - K) / (np.linalg.norm(Kfd) + 1e-30)
    assert err < 1e-5, f"curved tangent err {err:.2e}\nFD=\n{Kfd}\nAN=\n{K}"
    sym = np.linalg.norm(K - K.T) / (np.linalg.norm(K) + 1e-30)
    assert sym < 1e-10, f"frictionless geometric tangent must be symmetric: {sym:.2e}"
    print(f"[curved quad] κ_ξη={gb['kappa'][0,1]:.3f} |H_ss|={np.linalg.norm(H[0:3,0:3]):.3e} "
          f"analytic-vs-FD={err:.2e} sym={sym:.1e} OK")


def test_analytic_vs_fd_tilted_largeslide():
    """large-sliding on a flat but ROTATED master with the slave well off-center: the
    ∂n/∂u master-coupling dominates; analytic == FD (this is the regime kn·BᵀB-only
    Newton degrades — the master blocks of the tangent are wrong without H)."""
    th = np.radians(35)
    Rt = np.array([[1, 0, 0], [0, np.cos(th), -np.sin(th)], [0, np.sin(th), np.cos(th)]])
    X = _quad(L=2.0) @ Rt.T
    refdir = Rt @ np.array([0, 0, 1.0])
    n0 = Rt @ np.array([0, 0, 1.0])
    xs = np.array([0.6, -0.5, 0.0]) @ Rt.T - 0.03 * n0   # off-center + penetrating
    kn = 1e6
    Kfd = tangent_fd(4, X, xs, kn, refdir)
    Kan = tangent_analytic(4, X, xs, kn, refdir)
    err = np.linalg.norm(Kfd - Kan) / (np.linalg.norm(Kfd) + 1e-30)
    assert err < 1e-5, f"tilted/large-slide err {err:.2e}"
    print(f"[tilted slide] analytic vs FD={err:.2e} OK")


def test_analytic_vs_fd_tri3_curved_via_slide():
    """tri-3 facet: a_{αβ}=0 ⇒ κ=0 ⇒ flat, but the master blocks still carry ∂n/∂u
    (facet tilt under master motion). analytic == FD; slave Hessian == 0 (gate 3 holds
    for tri-3 too)."""
    X = np.array([[0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]])
    refdir = np.array([0, 0, 1.0])
    xs = np.array([0.3, 0.3, -0.01])
    kn = 1e6
    K, H = tangent_analytic(3, X, xs, kn, refdir, return_H=True)
    assert np.linalg.norm(H[0:3, 0:3]) < 1e-12, "tri-3 is flat ⇒ slave Hessian 0"
    Kfd = tangent_fd(3, X, xs, kn, refdir)
    err = np.linalg.norm(Kfd - K) / (np.linalg.norm(Kfd) + 1e-30)
    assert err < 1e-5, f"tri-3 tangent err {err:.2e}"
    print(f"[tri3] slave-Hess=0, analytic vs FD={err:.2e} OK")


# ----------------------------------------------------------------- Hertz reference
def hertz_sphere(P, R, Estar):
    """rigid sphere radius R pressed by load P into an elastic half-space (modulus
    Estar = E/(1−ν²)). Returns (a, p0, delta) — contact radius, peak pressure, approach.
    p(r) = p0·sqrt(1−(r/a)²)."""
    a = (3.0 * P * R / (4.0 * Estar)) ** (1.0 / 3.0)
    p0 = 3.0 * P / (2.0 * np.pi * a * a)
    delta = a * a / R
    return a, p0, delta


def hertz_p(r, a, p0):
    r = np.asarray(r, float)
    out = np.where(r < a, p0 * np.sqrt(np.clip(1.0 - (r / a) ** 2, 0.0, None)), 0.0)
    return out


def test_hertz_reference_selfconsistent():
    """the Hertz pressure integrates back to the applied load: ∫_0^a p(r) 2πr dr = P,
    and p0 = 3P/(2πa²). (The numeric C++ Hertz benchmark lands in
    tests/test_adr39_contact_p2b2c_hertz.py; this pins the analytic reference.)"""
    P, R, Estar = 100.0, 5.0, 1.0e4
    a, p0, delta = hertz_sphere(P, R, Estar)
    rr = np.linspace(0, a, 200001)
    _trap = getattr(np, "trapezoid", getattr(np, "trapz", None))   # numpy 2.0 renamed trapz
    integral = _trap(hertz_p(rr, a, p0) * 2 * np.pi * rr, rr)
    assert abs(integral - P) / P < 1e-3, f"∫p dA={integral} != P={P}"
    assert abs(p0 - 3 * P / (2 * np.pi * a * a)) < 1e-9
    print(f"[hertz] a={a:.4f} p0={p0:.4f} δ={delta:.5f}  ∫p·dA={integral:.4f}≈P={P}  OK")


if __name__ == "__main__":
    print("=== ADR-39 P2b-2c / B3 consistent ∂n/∂u normal tangent oracle ===")
    test_main_recovers_shipped_flat()
    test_analytic_vs_fd_flat()
    test_analytic_vs_fd_curved_quad()
    test_analytic_vs_fd_tilted_largeslide()
    test_analytic_vs_fd_tri3_curved_via_slide()
    test_hertz_reference_selfconsistent()
    print("ALL B3 ORACLE CHECKS PASSED")
