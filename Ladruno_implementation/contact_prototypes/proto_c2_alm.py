"""ADR-41 C2 ORACLE — frictionless commit-cycle ALM on the mortar operators.

Oracle-first (the fork discipline): validate the ENFORCEMENT layer that turns C1's
mortar geometry (D_IJ, M_IK, weighted gap g̃_I) into a SOLVED frictionless contact —
the normal multiplier λ_N, the penalty→Uzawa augmentation, and the F_c / K_c assembly —
in pure numpy BEFORE wiring it into the OpenSees handler / FE adapter / Domain (C2 is
NOT header-only; it links into the build, so de-risk the mechanics here first).

The mechanics (standard, not dual, basis — ADR-41 §How / [[_adr41_c2_design]]):
  - normalized gap   ḡ_I = g̃_I / a_I ,  a_I = Σ_J D_IJ  (the C1 row sum = ∫N_I^s dΓ)
  - nodal pressure   t_I = min(0, λ_I + epsN·ḡ_I)        (compression ≤ 0; active iff < 0)
  - contact force    F^s = +D·t , F^m = −Mᵀ·t  (along the per-facet normal n) —
                     SELF-EQUILIBRATING via partition of unity (Σ_K D_KI = a_I = Σ_L M_IL)
  - material tangent K_c = epsN · B̃ᵀ · diag(act/a) · B̃ ,  B̃ = [D, −M]  (SPD on the active set)
  - Uzawa augment    λ_I ← min(0, λ_I + epsN·ḡ_I)  (ONCE per Domain::commit — EmbeddedRebar precedent)
NO new DOFs, NO saddle-point solver: λ is Domain-side state, the system stays the slave∪master
displacement DOFs only (the eqn count is invariant across augmentations — gate T6).

Gates (capstone C2 / ADR-41 P2):
  T1  constant pressure → t_I == −p̄ and F^s is the CONSISTENT nodal load ∫N_I p̄ (reuse C1 D,M).
  T2  K_c vs FD (material/geometry-frozen) on a ROTATED non-matched config to 1e-6.
  T3  Uzawa drives ‖ḡ‖_∞ → augTol at FINITE epsN in ≤ maxAug; converged state is epsN-INDEPENDENT.
  T4  release/reopen → exact F=0 and λ→0.
  T5  mesh-refinement consistency: recovered uniform interface pressure stays constant (the
      elliptic-Hertz half-space study is DEFERRED to the C2.2 C++ phase — it needs real elastic
      coupling the mortar-mechanics oracle does not model; noted honestly).
  T6  eqn-count / active-set: system dimension invariant across augmentations (no new DOFs).

Run: <pythoncore-3.12> proto_c2_alm.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
TOL = 1e-9
_fails = 0


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ============================================================ C1 mortar core
# Mirror of proto_c1_mortar.py (validated 30/30) — re-included so this oracle is
# standalone (the house pattern: proto_a2 re-included projection rather than import).
# Only the pieces C2 consumes: assemble() -> global D, M, g̃ on a meshed interface.
def shape(nps, xi, eta):
    if nps == 4:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dNe = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    else:
        N = np.array([1 - xi - eta, xi, eta])
        dNx = np.array([-1.0, 1.0, 0.0])
        dNe = np.array([-1.0, 0.0, 1.0])
    return N, dNx, dNe


def project(nps, X, xs, tolR=1e-12, maxit=10):
    X = np.asarray(X, float)[:nps]
    xi, eta = (0.0, 0.0) if nps == 4 else (1 / 3, 1 / 3)
    conv = False
    for _ in range(maxit):
        N, dNx, dNe = shape(nps, xi, eta)
        xbar = N @ X
        g1, g2 = dNx @ X, dNe @ X
        d = xs - xbar
        R0, R1 = d @ g1, d @ g2
        if np.hypot(R0, R1) < tolR:
            conv = True
            break
        K = np.array([[g1 @ g1, g1 @ g2], [g2 @ g1, g2 @ g2]])
        det = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]
        sc = np.linalg.norm(g1) * np.linalg.norm(g2)
        if abs(det) < 1e-14 * (sc + 1e-300):
            return xi, eta, -1
        xi += (K[1, 1] * R0 - K[0, 1] * R1) / det
        eta += (-K[1, 0] * R0 + K[0, 0] * R1) / det
    if not conv:
        return xi, eta, -1
    t = 1e-9
    inb = ((-1 - t <= xi <= 1 + t and -1 - t <= eta <= 1 + t) if nps == 4
           else (xi >= -t and eta >= -t and xi + eta <= 1 + t))
    return xi, eta, (0 if inb else 1)


def normal_oriented(nps, xi, eta, X, refDir):
    X = np.asarray(X, float)[:nps]
    _, dNx, dNe = shape(nps, xi, eta)
    raw = np.cross(dNx @ X, dNe @ X)
    j = np.linalg.norm(raw)
    if j < 1e-300:
        return None
    n = raw / j
    p = n @ refDir
    if abs(p) < 1e-12 * (np.linalg.norm(refDir) + 1e-300):
        return None
    return -n if p < 0 else n


def project_full(nps, X, xs, refDir, tolR=1e-12, maxit=10):
    X = np.asarray(X, float)[:nps]
    p = {"status": -1, "n": np.zeros(3), "phi": np.zeros(4), "gap": 0.0}
    xi, eta, st = project(nps, X, xs, tolR, maxit)
    p["status"] = st
    if st < 0:
        return p
    N, _, _ = shape(nps, xi, eta)
    n = normal_oriented(nps, xi, eta, X, refDir)
    if n is None:
        p["status"] = -1
        return p
    p["n"] = n
    p["gap"] = n @ (xs - N @ X)
    p["phi"][:len(N)] = N
    return p


def signed_area(P):
    P = np.asarray(P, float)
    n = len(P)
    return 0.5 * sum(P[i, 0] * P[(i + 1) % n, 1] - P[(i + 1) % n, 0] * P[i, 1] for i in range(n))


def ensure_ccw(P):
    return P if signed_area(P) >= 0 else P[::-1].copy()


def is_convex2(poly, tol=1e-12):
    P = np.asarray(poly, float)
    n = len(P)
    if n < 3:
        return False
    sc = max(np.linalg.norm(P[(i + 1) % n] - P[i]) for i in range(n))
    if sc < 1e-300:
        return False
    for i in range(n):
        a = P[(i + 1) % n] - P[i]
        b = P[(i + 2) % n] - P[(i + 1) % n]
        if (a[0] * b[1] - a[1] * b[0]) < -tol * sc * sc:
            return False
    return True


def _iline(S, E, A, B):
    d1, d2 = E - S, B - A
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-300:
        return S.copy()
    t = ((A[0] - S[0]) * d2[1] - (A[1] - S[1]) * d2[0]) / den
    return S + t * d1


def clip_polygon(subject, clip):
    out = [np.asarray(p, float) for p in subject]
    nC = len(clip)
    for i in range(nC):
        if not out:
            break
        A, B = np.asarray(clip[i], float), np.asarray(clip[(i + 1) % nC], float)
        edge = B - A
        side = lambda P: edge[0] * (P[1] - A[1]) - edge[1] * (P[0] - A[0])
        inp, out = out, []
        S = inp[-1]
        sS = side(S)
        for E in inp:
            sE = side(E)
            if sE >= 0.0:
                if sS < 0.0:
                    out.append(_iline(S, E, A, B))
                out.append(E)
            elif sS >= 0.0:
                out.append(_iline(S, E, A, B))
            S, sS = E, sE
    return out


def dedupe(poly, tol=1e-12):
    n = len(poly)
    return [poly[i] for i in range(n) if np.linalg.norm(poly[i] - poly[(i + 1) % n]) > tol]


def tri_quadrature(order):
    if order <= 2:
        return ([(2 / 3, 1 / 6, 1 / 6), (1 / 6, 2 / 3, 1 / 6), (1 / 6, 1 / 6, 2 / 3)],
                [1 / 3, 1 / 3, 1 / 3])
    a1, b1, w1 = 0.445948490915965, 0.108103018168070, 0.223381589678011
    a2, b2, w2 = 0.091576213509771, 0.816847572980459, 0.109951743655322
    return ([(b1, a1, a1), (a1, b1, a1), (a1, a1, b1),
             (b2, a2, a2), (a2, b2, a2), (a2, a2, b2)], [w1, w1, w1, w2, w2, w2])


def inverse_isomap_2d(nps, UV, q, tolR=1e-13, maxit=20):
    if nps == 3:
        Jm = np.column_stack([UV[1] - UV[0], UV[2] - UV[0]])
        det = Jm[0, 0] * Jm[1, 1] - Jm[0, 1] * Jm[1, 0]
        sc = np.linalg.norm(Jm[:, 0]) * np.linalg.norm(Jm[:, 1])
        if abs(det) < 1e-14 * (sc + 1e-300):
            return 0.0, 0.0, False
        xi, eta = np.linalg.solve(Jm, q - UV[0])
        return xi, eta, True
    xi = eta = 0.0
    conv = False
    for _ in range(maxit):
        N, dNx, dNe = shape(4, xi, eta)
        r = N @ UV - q
        if np.linalg.norm(r) < tolR:
            conv = True
            break
        Jm = np.column_stack([dNx @ UV, dNe @ UV])
        if abs(np.linalg.det(Jm)) < 1e-300:
            break
        d = np.linalg.solve(Jm, -r)
        xi += d[0]
        eta += d[1]
    return xi, eta, conv


def facet_normal(X, nps, refDir):
    raw = (np.cross(X[1] - X[0], X[2] - X[0]) if nps == 3
           else np.cross(X[2] - X[0], X[3] - X[1]))
    n = raw / np.linalg.norm(raw)
    return -n if (n @ refDir) < 0 else n


def aux_plane(Xm, nps_m, refDir):
    x0 = Xm[:nps_m].mean(axis=0)
    n0 = facet_normal(Xm, nps_m, refDir)
    e1 = Xm[1] - Xm[0]
    e1 = e1 - (e1 @ n0) * n0
    e1 = e1 / np.linalg.norm(e1)
    return x0, n0, e1, np.cross(n0, e1)


def to2d(X3d, x0, e1, e2):
    d = X3d - x0
    return np.column_stack([d @ e1, d @ e2])


def mortar_pair(Xs, nps_s, Xm, nps_m, refDir, order=2):
    x0, n0, e1, e2 = aux_plane(Xm, nps_m, refDir)
    n_s = facet_normal(Xs, nps_s, refDir)
    cos_t = abs(n_s @ n0)
    if cos_t < 1e-12:
        return None
    UVs = ensure_ccw(to2d(Xs[:nps_s], x0, e1, e2))
    UVm = ensure_ccw(to2d(Xm[:nps_m], x0, e1, e2))
    if not is_convex2(UVs) or not is_convex2(UVm):
        return None
    UVs_node = to2d(Xs[:nps_s], x0, e1, e2)
    poly = dedupe(clip_polygon(list(UVs), list(UVm)))
    if len(poly) < 3:
        return None
    A_min = min(abs(signed_area(UVs)), abs(signed_area(UVm)))
    if abs(signed_area(poly)) < 1e-12 * A_min:
        return None
    cen = np.mean(poly, axis=0)
    bary, wts = tri_quadrature(order)
    D = np.zeros((nps_s, nps_s))
    M = np.zeros((nps_s, nps_m))
    g = np.zeros(nps_s)
    n = len(poly)
    for i in range(n):
        tri = (cen, poly[i], poly[(i + 1) % n])
        A_aux = abs(signed_area(list(tri)))
        if A_aux < 1e-12 * A_min:
            continue
        A_phys = A_aux / cos_t
        V0, V1, V2 = tri
        for (L0, L1, L2), w in zip(bary, wts):
            q = L0 * V0 + L1 * V1 + L2 * V2
            xi_s, eta_s, ok = inverse_isomap_2d(nps_s, UVs_node, q)
            if not ok:
                continue
            Ns, _, _ = shape(nps_s, xi_s, eta_s)
            x_s = Ns @ Xs[:nps_s]
            p = project_full(nps_m, Xm, x_s, refDir)
            if p["status"] < 0:
                continue
            phi = p["phi"][:nps_m]
            gN = p["n"] @ (x_s - phi @ Xm[:nps_m])
            wJ = w * A_phys
            D += wJ * np.outer(Ns, Ns)
            M += wJ * np.outer(Ns, phi)
            g += wJ * Ns * gN
    return {"D": D, "M": M, "g": g}


def assemble(s_nodes, s_facets, m_nodes, m_facets, refDir, order=2):
    Ns, Nm = len(s_nodes), len(m_nodes)
    D = np.zeros((Ns, Ns))
    M = np.zeros((Ns, Nm))
    g = np.zeros(Ns)
    for sf in s_facets:
        nps_s = len(sf)
        Xs = np.zeros((4, 3))
        Xs[:nps_s] = s_nodes[list(sf)]
        for mf in m_facets:
            nps_m = len(mf)
            Xm = np.zeros((4, 3))
            Xm[:nps_m] = m_nodes[list(mf)]
            res = mortar_pair(Xs, nps_s, Xm, nps_m, refDir, order)
            if res is None:
                continue
            for a in range(nps_s):
                I = sf[a]
                g[I] += res["g"][a]
                for b in range(nps_s):
                    D[I, sf[b]] += res["D"][a, b]
                for k in range(nps_m):
                    M[I, mf[k]] += res["M"][a, k]
    return D, M, g


def grid_nodes(nx, ny, z=0.0, x0=0.0, y0=0.0, sx=1.0, sy=1.0):
    pts = [[x0 + sx * i / nx, y0 + sy * j / ny, z]
           for j in range(ny + 1) for i in range(nx + 1)]
    return np.array(pts, float), (lambda i, j: j * (nx + 1) + i)


def quad_mesh(nx, ny, **kw):
    nodes, idx = grid_nodes(nx, ny, **kw)
    facets = [(idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1))
              for j in range(ny) for i in range(nx)]
    return nodes, facets


# ===================================================== C2: the ALM enforcement layer
def nodal_pressure(lam, gbar, epsN):
    """t_I = min(0, λ_I + epsN·ḡ_I); active mask where the argument is < 0."""
    P = lam + epsN * gbar
    act = P < 0.0
    return np.where(act, P, 0.0), act


def contact_force(D, M, t):
    """Self-equilibrating normal force (scalar/normal component): F^s = D·t, F^m = −Mᵀ·t."""
    return D @ t, -(M.T @ t)


def Kc_material(D, M, a, act, epsN):
    """Material/penalty tangent K_c = epsN · B̃ᵀ diag(act/a) B̃, B̃ = [D, −M] (slave|master)."""
    B = np.hstack([D, -M])                      # (Ns, Ns+Nm)
    W = np.where(act, 1.0 / a, 0.0)
    return epsN * (B.T @ (W[:, None] * B))


def rot_matrix(axis, ang):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    c, s = np.cos(ang), np.sin(ang)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) * c + s * K + (1 - c) * np.outer(a, a)


# ======================================================================= tests
print("ADR-41 C2 oracle — frictionless commit-cycle ALM on the mortar operators")
RD = np.array([0.0, 0.0, 1.0])

# --- T1: constant pressure -> t_I == −p̄ and F^s == the consistent nodal load ∫N_I p̄ ---
print("\nT1  constant pressure: t_I == −p̄ and F^s is the CONSISTENT nodal load (self-equilibrated)")
sN, sF = quad_mesh(2, 2)
mN, mF = quad_mesh(3, 3)
D, M, g = assemble(sN, sF, mN, mF, RD)
a = D.sum(1)
pbar = 4.1
t = np.full(len(sN), -pbar)                     # impose uniform compressive pressure
Fs, Fm = contact_force(D, M, t)
check("F^s == −p̄·a_I (∫N_I p̄, consistent nodal load)", np.allclose(Fs, -pbar * a, atol=1e-9),
      f"max|Δ|={np.abs(Fs + pbar * a).max():.2e}")
check("self-equilibrium ΣF^s + ΣF^m == 0", abs(Fs.sum() + Fm.sum()) < 1e-9,
      f"Σ={Fs.sum() + Fm.sum():.2e}")
check("total normal force == −p̄·A_interface", abs(Fs.sum() + pbar * 1.0) < 1e-9,
      f"ΣF^s={Fs.sum():.6f}")

# --- T2: K_c vs FD (geometry frozen) on a ROTATED non-matched config to 1e-6 -------
print("\nT2  K_c == FD of F_c (geometry/active-set frozen) on a ROTATED config (≤1e-6)")
R = rot_matrix([0.3, 1.0, 0.5], 0.6)
sNr = sN @ R.T
mNr = mN @ R.T
Dr, Mr, gr = assemble(sNr, sF, mNr, mF, RD @ R.T)   # rotated interface
ar = Dr.sum(1)
Ns, Nm = len(sNr), len(mNr)
epsN = 500.0
lam = np.zeros(Ns)
# a MIXED active set: a gap field linear in x (penetration on x<0.5, separation on x>0.5)
h_node = (sN[:, 0] - 0.5) * 0.12                 # per-node target gap, −0.06 … +0.06
g0 = Dr @ h_node                                 # weighted ⇒ ~half the nodes active


def Fc_of_u(us, um, frozen_act=None, frozen=(Dr, Mr, ar, g0, lam, epsN)):
    """F_c(displacements) with geometry D,M,a,g0 FROZEN (material tangent regime). With
    frozen_act, the active SET is held fixed too (so the FD does not straddle the KKT kink —
    that is exactly the within-augmentation-sweep frozen-set regime the tangent linearizes)."""
    Dl, Ml, al, g0l, laml, eps = frozen
    gb = (g0l + Dl @ us - Ml @ um) / al          # linearized normalized weighted gap
    if frozen_act is None:
        tt, ac = nodal_pressure(laml, gb, eps)
    else:
        ac = frozen_act
        tt = np.where(ac, laml + eps * gb, 0.0)
    Fs_, Fm_ = contact_force(Dl, Ml, tt)
    return np.concatenate([Fs_, Fm_]), ac


us0, um0 = np.zeros(Ns), np.zeros(Nm)
F0, act0 = Fc_of_u(us0, um0)
check("active set non-trivial (mixed: some in contact, some open)", 0 < act0.sum() < Ns,
      f"active {act0.sum()}/{Ns}")
Kc = Kc_material(Dr, Mr, ar, act0, epsN)         # analytic material tangent
# central FD of F_c w.r.t. all slave+master normal DOFs, active set HELD at act0
ndof = Ns + Nm
Kfd = np.zeros((ndof, ndof))
h = 1e-6
for j in range(ndof):
    up = np.zeros(ndof); up[j] = h
    Fp, _ = Fc_of_u(up[:Ns], up[Ns:], frozen_act=act0)
    Fm_, _ = Fc_of_u(-up[:Ns], -up[Ns:], frozen_act=act0)
    Kfd[:, j] = (Fp - Fm_) / (2 * h)
relK = np.abs(Kc - Kfd).max() / (np.abs(Kfd).max() + 1e-300)
check("K_c == FD(F_c) (≤1e-6 rel)", relK <= 1e-6, f"max rel = {relK:.2e}")
check("K_c symmetric (frictionless)", np.allclose(Kc, Kc.T, atol=1e-9))

# --- T3: Uzawa drives ‖ḡ‖→augTol at FINITE epsN; converged state epsN-INDEPENDENT --
print("\nT3  Uzawa: ‖ḡ‖_∞ → augTol in ≤ maxAug; converged penetration epsN-INDEPENDENT")


def solve_step(D, a, z0, k, f_ext, lam, epsN, newton_it=50, tol=1e-13):
    """Equilibrium k·u + D·t = f_ext for FIXED λ (rigid foundation at z=0, master fixed).
    t = min(0, λ + epsN·ḡ), ḡ = D(z0+u)/a. Active-set Newton (piecewise linear)."""
    Ns = len(z0)
    u = np.zeros(Ns)
    for _ in range(newton_it):
        gtil = D @ (z0 + u)
        gb = gtil / a
        t, act = nodal_pressure(lam, gb, epsN)
        R = k * u + D @ t - f_ext
        if np.linalg.norm(R) < tol:
            break
        W = np.where(act, 1.0 / a, 0.0)
        J = k * np.eye(Ns) + epsN * (D @ (W[:, None] * D))
        u -= np.linalg.solve(J, R)
    gtil = D @ (z0 + u)
    return u, gtil / a, act


def uzawa(D, a, z0, k, f_ext, epsN, augTol=1e-10, maxAug=60):
    lam = np.zeros(len(z0))
    naug = 0
    for naug in range(1, maxAug + 1):
        u, gb, act = solve_step(D, a, z0, k, f_ext, lam, epsN)
        pen = np.maximum(0.0, -gb)               # penetration (ḡ < 0)
        if pen.max() < augTol:
            break
        lam = np.minimum(0.0, lam + epsN * gb)   # Uzawa update (compression clamp)
    return u, gb, lam, naug


# slave 2×2 (rigid foundation at z=0); non-uniform downward load; elastic springs k
sN2, sF2 = quad_mesh(2, 2)
D2, _, _ = assemble(sN2, sF2, sN2, sF2, RD)      # D on a self-matched slave (its Gram matrix)
a2 = D2.sum(1)
z0 = np.full(len(sN2), 0.10)                     # start 0.10 above the foundation
k = 50.0
f_ext = -np.array([1, 3, 1, 3, 8, 3, 1, 3, 1], float)   # non-uniform downward (9 nodes, 2×2 grid)
u_lo, gb_lo, lam_lo, n_lo = uzawa(D2, a2, z0, k, f_ext, epsN=1e3)
u_hi, gb_hi, lam_hi, n_hi = uzawa(D2, a2, z0, k, f_ext, epsN=1e5)
check("Uzawa converges ‖penetration‖→0 (≤1e-10) at epsN=1e3", np.maximum(0, -gb_lo).max() < 1e-10,
      f"pen={np.maximum(0, -gb_lo).max():.2e} in {n_lo} aug")
check("Uzawa converges at epsN=1e5", np.maximum(0, -gb_hi).max() < 1e-10, f"in {n_hi} aug")
check("converged displacement is epsN-INDEPENDENT (1e3 vs 1e5, ≤1e-6)",
      np.allclose(u_lo, u_hi, atol=1e-6), f"max|Δu|={np.abs(u_lo - u_hi).max():.2e}")
# at convergence the pressure t → the multiplier λ (penalty term vanishes as ḡ→0); global
# equilibrium k·u + D·t = f_ext must hold, and the contact pressure must be compressive (≤0).
t_lo = np.minimum(0.0, lam_lo + 1e3 * gb_lo)
resid = np.linalg.norm(k * u_lo + D2 @ t_lo - f_ext)
check("global equilibrium k·u + D·t == f_ext at convergence (≤1e-8)", resid < 1e-8,
      f"‖R‖={resid:.2e}")
check("converged pressure == λ (penalty term vanished) and compressive (≤0)",
      np.allclose(t_lo, lam_lo, atol=1e-7) and lam_lo.max() <= 1e-12,
      f"max|t−λ|={np.abs(t_lo - lam_lo).max():.2e}")

# --- T4: release/reopen -> exact F=0 and λ→0 --------------------------------------
print("\nT4  release: flip the load to tension → active set empties, F_contact=0, λ→0")
f_tension = +np.array([1, 3, 1, 3, 8, 3, 1, 3, 1], float)   # pull UP (separate)
u_r, gb_r, lam_r, n_r = uzawa(D2, a2, z0, k, f_tension, epsN=1e3, maxAug=80)
t_r, act_r = nodal_pressure(lam_r, gb_r, 1e3)
Fc_r = D2 @ t_r
check("no active contact nodes after release", act_r.sum() == 0, f"active={act_r.sum()}")
check("contact force F == 0", np.linalg.norm(Fc_r) < 1e-9, f"‖F‖={np.linalg.norm(Fc_r):.2e}")
check("multiplier λ → 0 (clamped open)", np.linalg.norm(lam_r) < 1e-9, f"‖λ‖={np.linalg.norm(lam_r):.2e}")

# --- T5: mesh-refinement consistency (Hertz half-space study deferred to C2.2 C++) -
print("\nT5  refinement consistency: recovered uniform interface pressure stays constant")
# (Honest scope: the elliptic-Hertz p(r)=p₀√(1−(r/a)²) needs elastic half-space coupling the
#  mortar-mechanics oracle does not model; it is a C2.2 gate against a real elastic mesh. Here
#  we check the weaker mesh-objectivity surrogate: a uniform pressure recovers through D⁻¹M on
#  progressively finer non-matched meshes with no drift.)
errs = []
for (ns, nm) in [(2, 3), (3, 4), (4, 5)]:
    sn_, sf_ = quad_mesh(ns, ns)
    mn_, mf_ = quad_mesh(nm, nm)
    Dh, Mh, _ = assemble(sn_, sf_, mn_, mf_, RD)
    um = np.full(len(mn_), 5.0)
    us = np.linalg.solve(Dh, Mh @ um)            # transfer a constant master field
    errs.append(np.abs(us - 5.0).max() / 5.0)
check("uniform pressure recovered constant under refinement (all ≤1e-6)", max(errs) <= 1e-6,
      f"errs={['%.1e' % e for e in errs]}")

# --- T6: eqn-count / no-new-DOF: system dimension invariant across augmentations ---
print("\nT6  no new DOFs: system dimension is invariant across every augmentation (Uzawa, not Lagrange)")
dims = []
lam = np.zeros(len(sN2))
for _ in range(5):
    u, gb, act = solve_step(D2, a2, z0, k, f_ext, lam, 1e3)
    J = k * np.eye(len(sN2)) + 1e3 * (D2 @ (np.where(act, 1 / a2, 0)[:, None] * D2))
    dims.append(J.shape)
    lam = np.minimum(0.0, lam + 1e3 * gb)
check("Jacobian shape constant across augmentations (no λ-DOFs added/removed)",
      len(set(dims)) == 1, f"dims={set(dims)}")
check("system size == #displacement DOFs only (Uzawa keeps λ off the equations)",
      dims[0] == (len(sN2), len(sN2)))

print(f"\n{'='*68}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(C2 mortar-ALM oracle)\n{'='*68}")
sys.exit(1 if _fails else 0)
