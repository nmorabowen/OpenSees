"""ADR-62 P2 ORACLE — integral-mortar KINEMATIC mesh-tie via constraint emission.

P1 (shipped) tied each slave NODE to one master facet by point collocation
(u_s = Σ N_i(ξ̄_s) u_{m,i}). P2 ties the slave SURFACE to the master surface in the
weak (integral) sense: assemble the mortar operators over the overlap and CONDENSE
to a per-slave-node transfer

        D u_s = M u_m      ⇒      u_s = P u_m ,   P = D⁻¹ M

where (STANDARD basis, reusing the SHIPPED LadrunoMortarKernel / proto_c1 clip):
        D_IJ = ∫_Γ N_I^s N_J^s dΓ     (slave–slave, the interface consistent mass)
        M_IK = ∫_Γ N_I^s φ_K^m dΓ     (slave–master)

THE KEY RESULT this oracle proves: P is a per-slave-node row tied to MASTER nodes
ONLY (D is pre-inverted ONCE at setup over the whole interface), so each row emits as
an ordinary EQ_Constraint  u_s = Σ_k P_sk u_{m,k}  with master-only retainers —
NO slave-slave coupling, NO MP-chains. The shipped LadrunoProjectionHandler (ADR-30)
already accepts exactly this topology (dense rows, master-only retainers), so P2 needs
NO handler change and NO kernel change — only a setup-time generator + one SPD solve.

Gates (falsifiable, numeric):
  T1  D symmetric + SPD (every covered slave node has ∫N²>0); condition number sane.
  T2  partition of unity: P·1 = 1  (Σ_k P_sk = 1)  ⇔  a constant field transfers exactly.
  T3  linear completeness / CONSTANT-STRESS patch: P·(linear field on masters) = the
      linear field at the slave nodes (to ~1e-9) ⇒ a constant-strain state crosses the
      NON-MATCHING tie exactly. (This is what makes the FE patch test pass.)
  T4  no chains: every P row's support is master nodes only (emission topology the
      projection handler accepts).
  T5  uncovered-slave detection: a slave node not covered by the master surface gives a
      ~0 D row ⇒ D singular ⇒ MUST be detected and refused (the mortar OQ-3 analog).
  T6  conforming gap: the weighted gap g̃/area ≈ 0 for coincident surfaces, and grows
      for a separated slave ⇒ the "require conforming-at-interface" check.
  T7  mortar ≠ collocation: on a genuinely non-matching mesh the two operators DIFFER
      (P_mortar ≠ P_collocation) yet BOTH are linearly complete — P2 is a distinct,
      variationally-consistent operator (optimal-convergence motivation).

Run: <pythoncore-3.12> proto_p2_mortar_tie.py   (numpy only, no build)
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


# ===================================================================== geometry
# Mirror of SRC/domain/contact/{LadrunoContactProjection,LadrunoMortarKernel}.h,
# validated by proto_a2_projection.py (7/7) and proto_c1_mortar.py (24/24). REUSED
# here verbatim (the clip + integrate is NOT the thing under test in P2 — the
# CONDENSATION P=D⁻¹M and the per-node emission topology are).
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


def project(nps, X, xs, tolR=1e-12, maxit=20):
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
        detK = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]
        scale = np.linalg.norm(g1) * np.linalg.norm(g2)
        if abs(detK) < 1e-14 * (scale + 1e-300):
            return xi, eta, -1
        dxi = (K[1, 1] * R0 - K[0, 1] * R1) / detK
        deta = (-K[1, 0] * R0 + K[0, 0] * R1) / detK
        xi += dxi
        eta += deta
    if not conv:
        return xi, eta, -1
    t = 1e-9
    if nps == 4:
        inb = (-1 - t <= xi <= 1 + t) and (-1 - t <= eta <= 1 + t)
    else:
        inb = (xi >= -t) and (eta >= -t) and (xi + eta <= 1 + t)
    return xi, eta, (0 if inb else 1)


def normal_oriented(nps, xi, eta, X, refDir):
    X = np.asarray(X, float)[:nps]
    _, dNx, dNe = shape(nps, xi, eta)
    g1, g2 = dNx @ X, dNe @ X
    raw = np.cross(g1, g2)
    j = np.linalg.norm(raw)
    if j < 1e-300:
        return None
    n = raw / j
    p = n @ refDir
    if abs(p) < 1e-12 * (np.linalg.norm(refDir) + 1e-300):
        return None
    return -n if p < 0 else n


def project_full(nps, X, xs, refDir):
    X = np.asarray(X, float)[:nps]
    p = {"status": -1, "xi": 0.0, "eta": 0.0, "gap": 0.0, "n": np.zeros(3), "phi": np.zeros(4)}
    xi, eta, st = project(nps, X, xs)
    p["status"], p["xi"], p["eta"] = st, xi, eta
    if st < 0:
        return p
    N, _, _ = shape(nps, xi, eta)
    xbar = N @ X
    n = normal_oriented(nps, xi, eta, X, refDir)
    if n is None:
        p["status"] = -1
        return p
    p["n"] = n
    p["gap"] = n @ (xs - xbar)
    p["phi"][:len(N)] = N
    return p


def facet_normal(X, nps, refDir):
    raw = np.cross(X[1] - X[0], X[2] - X[0]) if nps == 3 else np.cross(X[2] - X[0], X[3] - X[1])
    n = raw / np.linalg.norm(raw)
    return -n if (n @ refDir) < 0 else n


def aux_plane(Xm, nps_m, refDir):
    x0 = Xm[:nps_m].mean(axis=0)
    n0 = facet_normal(Xm, nps_m, refDir)
    e1 = Xm[1] - Xm[0]
    e1 = e1 - (e1 @ n0) * n0
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n0, e1)
    return x0, n0, e1, e2


def to2d(X3d, x0, e1, e2):
    d = X3d - x0
    return np.column_stack([d @ e1, d @ e2])


def signed_area(poly):
    a, n = 0.0, len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def ensure_ccw(poly):
    return poly if signed_area(poly) >= 0 else poly[::-1].copy()


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
        A = np.asarray(clip[i], float)
        B = np.asarray(clip[(i + 1) % nC], float)
        edge = B - A

        def side(P):
            return edge[0] * (P[1] - A[1]) - edge[1] * (P[0] - A[0])

        inp, out, S = out, [], out[-1]
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
    out, n = [], len(poly)
    for i in range(n):
        if np.linalg.norm(poly[i] - poly[(i + 1) % n]) > tol:
            out.append(poly[i])
    return out


def tri_quadrature(order):
    if order <= 2:
        return [(2 / 3, 1 / 6, 1 / 6), (1 / 6, 2 / 3, 1 / 6), (1 / 6, 1 / 6, 2 / 3)], [1 / 3, 1 / 3, 1 / 3]
    a1, b1, w1 = 0.445948490915965, 0.108103018168070, 0.223381589678011
    a2, b2, w2 = 0.091576213509771, 0.816847572980459, 0.109951743655322
    b = [(b1, a1, a1), (a1, b1, a1), (a1, a1, b1), (b2, a2, a2), (a2, b2, a2), (a2, a2, b2)]
    return b, [w1, w1, w1, w2, w2, w2]


def is_convex2(poly, tol=1e-12):
    P, n = np.asarray(poly, float), len(poly)
    if n < 3:
        return False
    scale = max(np.linalg.norm(P[(i + 1) % n] - P[i]) for i in range(n))
    if scale < 1e-300:
        return False
    for i in range(n):
        a = P[(i + 1) % n] - P[i]
        b = P[(i + 2) % n] - P[(i + 1) % n]
        if (a[0] * b[1] - a[1] * b[0]) < -tol * scale * scale:
            return False
    return True


def inverse_isomap_2d(nps, UV, q, tolR=1e-13, maxit=20):
    if nps == 3:
        Jm = np.column_stack([UV[1] - UV[0], UV[2] - UV[0]])
        det = Jm[0, 0] * Jm[1, 1] - Jm[0, 1] * Jm[1, 0]
        scale = np.linalg.norm(Jm[:, 0]) * np.linalg.norm(Jm[:, 1])
        if abs(det) < 1e-14 * (scale + 1e-300):
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


def clip_subtris(Xs, nps_s, Xm, nps_m, refDir, area_tol=1e-12):
    x0, n0, e1, e2 = aux_plane(Xm, nps_m, refDir)
    n_s = facet_normal(Xs, nps_s, refDir)
    UVs = ensure_ccw(to2d(Xs[:nps_s], x0, e1, e2))
    UVm = ensure_ccw(to2d(Xm[:nps_m], x0, e1, e2))
    if not is_convex2(UVs) or not is_convex2(UVm):
        return None
    UVs_nodeorder = to2d(Xs[:nps_s], x0, e1, e2)
    poly = dedupe(clip_polygon(list(UVs), list(UVm)))
    if len(poly) < 3:
        return None
    A_s, A_m = abs(signed_area(UVs)), abs(signed_area(UVm))
    if abs(signed_area(poly)) < area_tol * min(A_s, A_m):
        return None
    cen = np.mean(poly, axis=0)
    subtris, n = [], len(poly)
    for i in range(n):
        tri = (cen, poly[i], poly[(i + 1) % n])
        if abs(signed_area(list(tri))) > area_tol * min(A_s, A_m):
            subtris.append(tri)
    if not subtris:
        return None
    return subtris, UVs_nodeorder, x0, n0, e1, e2, n_s


def mortar_pair(Xs, nps_s, Xm, nps_m, refDir, order=2):
    clip = clip_subtris(Xs, nps_s, Xm, nps_m, refDir)
    if clip is None:
        return None
    subtris, UVs, x0, n0, e1, e2, n_s = clip
    cos_t = abs(n_s @ n0)
    if cos_t < 1e-12:
        return None
    bary, wts = tri_quadrature(order)
    D = np.zeros((nps_s, nps_s))
    M = np.zeros((nps_s, nps_m))
    g = np.zeros(nps_s)
    area = 0.0
    for tri in subtris:
        A_phys = abs(signed_area(list(tri))) / cos_t
        V0, V1, V2 = tri
        for (L0, L1, L2), w in zip(bary, wts):
            q = L0 * V0 + L1 * V1 + L2 * V2
            xi_s, eta_s, ok = inverse_isomap_2d(nps_s, UVs, q)
            if not ok:
                continue
            Ns, _, _ = shape(nps_s, xi_s, eta_s)
            x_s = Ns @ Xs[:nps_s]
            p = project_full(nps_m, Xm, x_s, refDir)
            if p["status"] < 0:
                continue
            phi = p["phi"][:nps_m]
            x_m = phi @ Xm[:nps_m]
            wJ = w * A_phys
            area += wJ
            D += wJ * np.outer(Ns, Ns)
            M += wJ * np.outer(Ns, phi)
            g += wJ * Ns * (p["n"] @ (x_s - x_m))
    return {"D": D, "M": M, "g": g, "area": area}


def assemble(s_nodes, s_facets, m_nodes, m_facets, refDir, order=2):
    Ns, Nm = len(s_nodes), len(m_nodes)
    D = np.zeros((Ns, Ns))
    M = np.zeros((Ns, Nm))
    g = np.zeros(Ns)
    aslave = np.zeros(Ns)                       # ∫N_I dΓ per slave node (coverage)
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
                aslave[I] += res["D"][a].sum()  # Σ_J D_IJ = ∫N_I dΓ
                for b in range(nps_s):
                    D[I, sf[b]] += res["D"][a, b]
                for k in range(nps_m):
                    M[I, mf[k]] += res["M"][a, k]
    return D, M, g, aslave


def collocation_P(s_nodes, m_nodes, m_facets, refDir):
    """P1 collocation operator for the contrast gate: each slave node → its master
    facet shape weights at the projection (the nearest in-bounds facet)."""
    Ns, Nm = len(s_nodes), len(m_nodes)
    P = np.zeros((Ns, Nm))
    for I in range(Ns):
        xs = s_nodes[I]
        best = None
        for mf in m_facets:
            nps_m = len(mf)
            Xm = np.zeros((4, 3))
            Xm[:nps_m] = m_nodes[list(mf)]
            xi, eta, st = project(nps_m, Xm, xs)
            if st != 0:
                continue
            N, _, _ = shape(nps_m, xi, eta)
            xbar = N @ Xm[:nps_m]
            d = np.linalg.norm(xs - xbar)
            if best is None or d < best[0]:
                best = (d, mf, N[:nps_m])
        if best is not None:
            _, mf, N = best
            for k in range(len(mf)):
                P[I, mf[k]] += N[k]
    return P


# ============================================================== mesh generators
def grid_nodes(nx, ny, z=0.0, x0=0.0, y0=0.0, sx=1.0, sy=1.0):
    pts = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            pts.append([x0 + sx * i / nx, y0 + sy * j / ny, z])
    return np.array(pts, float), (lambda i, j: j * (nx + 1) + i)


def quad_mesh(nx, ny, **kw):
    nodes, idx = grid_nodes(nx, ny, **kw)
    facets = [(idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1))
              for j in range(ny) for i in range(nx)]
    return nodes, facets


# ===================================================================== tests
print("ADR-62 P2 oracle — integral-mortar kinematic tie (P = D⁻¹M condensation)")
RD = np.array([0.0, 0.0, 1.0])

# A genuinely NON-MATCHING flat interface in the z=0 plane: slave 3×3, master 2×2.
# No slave node coincides with a master node (3 vs 2 divisions) ⇒ true mortar.
s_nodes, s_facets = quad_mesh(3, 3, z=0.0)
m_nodes, m_facets = quad_mesh(2, 2, z=0.0)
D, M, g, aslave = assemble(s_nodes, s_facets, m_nodes, m_facets, RD, order=4)

# --- T1: D symmetric + SPD, sane conditioning ---------------------------------
print("\nT1  D = ∫N^s N^s dΓ : symmetric, SPD, well-conditioned")
sym = np.abs(D - D.T).max()
eig = np.linalg.eigvalsh(0.5 * (D + D.T))
cond = eig.max() / eig.min() if eig.min() > 0 else np.inf
check("D symmetric", sym < 1e-14, f"max|D-Dᵀ|={sym:.2e}")
check("D SPD (all covered slaves)", eig.min() > 1e-9, f"λ_min={eig.min():.3e}")
check("D well-conditioned", cond < 1e3, f"cond(D)={cond:.1f}")
check("∫dΓ over interface == 1.0 (unit square)", abs(aslave.sum() - 1.0) < 1e-9,
      f"Σ∫N_I={aslave.sum():.9f}")

# --- T2: condense P = D⁻¹M ; partition of unity (constant field) ---------------
print("\nT2  P = D⁻¹M : partition of unity (constant field transfers exactly)")
P = np.linalg.solve(D, M)
rowsum = P.sum(axis=1)
check("Σ_k P_sk == 1 for every slave node", np.allclose(rowsum, 1.0, atol=1e-9),
      f"max|ΣP-1|={np.abs(rowsum - 1).max():.2e}")

# --- T3: linear completeness / constant-stress patch --------------------------
print("\nT3  linear completeness: P·(linear field on masters) == field at slave nodes")
# three independent linear fields a + b·x + c·y over the master nodes
worst = 0.0
for (a, b, c) in [(1, 0, 0), (0.3, 1.7, 0), (0.0, 0.0, 2.5), (-1.0, 0.5, -0.8)]:
    fm = a + b * m_nodes[:, 0] + c * m_nodes[:, 1]
    fs_exact = a + b * s_nodes[:, 0] + c * s_nodes[:, 1]
    fs_mortar = P @ fm
    worst = max(worst, np.abs(fs_mortar - fs_exact).max())
check("linear field reproduced at slave nodes (constant-stress patch)", worst < 1e-9,
      f"max nodal err={worst:.2e}")

# --- T4: no chains — every P row's support is MASTER nodes only ----------------
print("\nT4  emission topology: each P row ties a slave to MASTER nodes only (no chains)")
# P is (Ns × Nm) by construction: columns index master nodes, rows index slaves.
# A non-trivial coupling exists (P not all-zero), and no row references a slave index.
nnz_per_row = (np.abs(P) > 1e-12).sum(axis=1)
check("P shape is (n_slave × n_master)", P.shape == (len(s_nodes), len(m_nodes)),
      f"shape={P.shape}")
check("every slave ties to ≥1 master (covered)", (nnz_per_row >= 1).all(),
      f"min nnz/row={int(nnz_per_row.min())}")
check("a non-matching slave couples to MULTIPLE masters (dense row)",
      nnz_per_row.max() >= 2, f"max nnz/row={int(nnz_per_row.max())}")

# --- T5: uncovered-slave detection (D singular ⇒ refuse) ----------------------
print("\nT5  uncovered slave node ⇒ ~0 D row ⇒ singular ⇒ MUST be detected")
# master surface only covers x∈[0,1]; add a slave column at x∈[1,2] (uncovered)
s2_nodes, s2_facets = quad_mesh(3, 3, z=0.0, x0=0.0, sx=2.0)   # spans x∈[0,2]
D2, M2, g2, a2 = assemble(s2_nodes, s2_facets, m_nodes, m_facets, RD, order=4)
uncovered = np.where(a2 < 1e-9)[0]
check("some slave nodes (x>1) are uncovered by the master surface", len(uncovered) > 0,
      f"#uncovered={len(uncovered)}")
# D2 is singular (uncovered rows ≈ 0) → a naive solve must NOT be attempted
eig2 = np.linalg.eigvalsh(0.5 * (D2 + D2.T))
check("uncovered coverage ⇒ D singular (λ_min≈0) ⇒ generator refuses", eig2.min() < 1e-9,
      f"λ_min={eig2.min():.3e}")

# --- T6: conforming gap g̃/area  (require coincident surfaces) -----------------
print("\nT6  weighted gap g̃: ≈0 for coincident surfaces, grows when separated")
gbar_coincident = np.abs(g[aslave > 1e-12] / aslave[aslave > 1e-12]).max()
check("coincident surfaces ⇒ g̃/area ≈ 0", gbar_coincident < 1e-9,
      f"max|g̃/a|={gbar_coincident:.2e}")
# separate the slave by 0.05 in +z ⇒ the weighted gap must reflect it
s_off = s_nodes.copy(); s_off[:, 2] += 0.05
Do, Mo, go, ao = assemble(s_off, s_facets, m_nodes, m_facets, RD, order=4)
gbar_sep = np.abs(go[ao > 1e-12] / ao[ao > 1e-12]).max()
check("separated slave ⇒ g̃/area ≈ 0.05 (the offset)", abs(gbar_sep - 0.05) < 1e-6,
      f"g̃/a={gbar_sep:.6f}")

# --- T7: mortar ≠ collocation, but both linearly complete ---------------------
print("\nT7  mortar vs collocation: distinct operators, both linearly complete")
Pcoll = collocation_P(s_nodes, m_nodes, m_facets, RD)
diff = np.abs(P - Pcoll).max()
check("P_mortar ≠ P_collocation on the non-matching mesh", diff > 1e-3,
      f"max|P_mortar-P_coll|={diff:.3e}")
# both must reproduce a linear field (linear completeness) ...
lin = 0.4 + 1.3 * m_nodes[:, 0] - 0.7 * m_nodes[:, 1]
lin_s = 0.4 + 1.3 * s_nodes[:, 0] - 0.7 * s_nodes[:, 1]
em = np.abs(P @ lin - lin_s).max()
ec = np.abs(Pcoll @ lin - lin_s).max()
check("mortar reproduces the linear field", em < 1e-9, f"err={em:.2e}")
check("collocation also reproduces the linear field", ec < 1e-9, f"err={ec:.2e}")
check("collocation partition of unity holds too", np.allclose(Pcoll.sum(1), 1.0, atol=1e-9),
      f"max|ΣP_coll-1|={np.abs(Pcoll.sum(1) - 1).max():.2e}")
print("      → both pass the constant/linear patch; mortar is the variationally-"
      "consistent (∫-weak) operator, collocation the point-wise one. The convergence-"
      "rate advantage of mortar for non-matching interfaces is the FE-level motivation.")

# ===================================================================== summary
print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'} — "
      "P2 mortar tie: P=D⁻¹M is a per-slave master-only transfer (no chains), "
      "linearly complete, with uncovered-node + conforming-gap guards.")
sys.exit(0 if _fails == 0 else 1)
