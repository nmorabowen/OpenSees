"""ADR-62 P2.1 ORACLE — DUAL (biorthogonal) basis ⇒ SPARSE mortar mesh-tie.

P2 (shipped) condenses the mortar tie with the STANDARD slave basis:
        D_IJ = ∫ N_I^s N_J^s dΓ   (slave–slave, a FULL/coupled matrix)
        M_IK = ∫ N_I^s φ_K^m dΓ
        P = D⁻¹ M
Because D⁻¹ is DENSE, every slave row of P couples to (potentially) EVERY master node ⇒
one large interface handler group `(LᵀML)` over all master DOFs. Correct, but O(interface²).

P2.1 replaces the slave TEST functions N_I^s with a DUAL (biorthogonal) basis ψ_I
(Wohlmuth 2000) chosen so that
        ∫ ψ_I N_J^s dΓ = δ_IJ ∫ N_I^s dΓ         ⇒   D_dual is DIAGONAL.
Then P = D_dual⁻¹ M_dual is LOCAL: slave node I ties ONLY to the masters under its own
facet support ⇒ sparse rows, small local handler groups, cheap projection.

THE CLEAN CONSTRUCTION (this oracle's headline): ψ = Aᵉ · N is a per-slave-FACET LINEAR
transform of the standard basis, with
        Aᵉ = diag(cᵉ) · (Dᵉ)⁻¹ ,   cᵉ_a = Σ_b Dᵉ_ab = ∫_e N_a  (row-sum of the facet mass),
so  Dᵉ_dual = Aᵉ Dᵉ = diag(cᵉ)  EXACTLY, for ANY (even distorted) facet. Crucially this uses
only the per-facet standard Dᵉ and Mᵉ the SHIPPED `LadrunoMortarKernel::integratePair` already
returns — so the C++ path needs **NO kernel change and NO handler change**: the dual transform
is a tiny per-facet solve at the `LadrunoTie` level (exactly how P2 dodged the handler change).

Gates (falsifiable, numeric):
  T1  biorthogonality: D_dual = ∫ψ N^s is DIAGONAL (off-diagonals ≈ 0); its diagonal equals
      the standard D's row-sum (∫N_I over the slave surface).
  T2  partition of unity: Σ_k P_dual,Ik = 1 (a constant field transfers exactly).
  T3  LINEAR COMPLETENESS / constant-stress patch: P_dual·(linear on masters) = the linear
      field at the slave nodes (to ~1e-9). This is what row-sum LUMPING destroys (T6).
  T4  SPARSITY: each dual P row's master support ⊆ the masters overlapping that slave node's
      facet(s); dual P has strictly FEWER nonzeros than the (dense) standard P.
  T5  SAME TIE on linear data: standard and dual P reproduce the SAME linear field (both are
      variationally consistent) even though the operators differ — dual is the sparsified sibling.
  T6  CONTRAST: row-sum LUMPING D (the naive "make D diagonal") keeps partition of unity but
      BREAKS linear completeness ⇒ FAILS the constant-stress patch. Dual is the correct sparsifier.

Run: <pythoncore-3.12> proto_p2_1_dual_mortar.py   (numpy only, no build)
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
# Verbatim from proto_p2_mortar_tie.py (itself validated by proto_a2/proto_c1). The clip
# + integrate is NOT what P2.1 tests — the DUAL CONDENSATION on top of it is.
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


def mortar_pair(Xs, nps_s, Xm, nps_m, refDir, order=4):
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
            wJ = w * A_phys
            D += wJ * np.outer(Ns, Ns)
            M += wJ * np.outer(Ns, phi)
    return {"D": D, "M": M}


# =============================================== assembly: standard / dual / lumped
def assemble_standard(s_nodes, s_facets, m_nodes, m_facets, refDir):
    """Global standard D (Ns×Ns, dense) and M (Ns×Nm)."""
    Ns, Nm = len(s_nodes), len(m_nodes)
    D, M = np.zeros((Ns, Ns)), np.zeros((Ns, Nm))
    for sf in s_facets:
        nps_s = len(sf)
        Xs = np.zeros((4, 3)); Xs[:nps_s] = s_nodes[list(sf)]
        for mf in m_facets:
            nps_m = len(mf)
            Xm = np.zeros((4, 3)); Xm[:nps_m] = m_nodes[list(mf)]
            res = mortar_pair(Xs, nps_s, Xm, nps_m, refDir)
            if res is None:
                continue
            for a in range(nps_s):
                I = sf[a]
                for b in range(nps_s):
                    D[I, sf[b]] += res["D"][a, b]
                for k in range(nps_m):
                    M[I, mf[k]] += res["M"][a, k]
    return D, M


def assemble_dual(s_nodes, s_facets, m_nodes, m_facets, refDir):
    """DUAL condensation: per slave FACET, build Dᵉ (nps×nps) + Mᵉ (nps×Nm) over the master
    facets it overlaps, form Aᵉ = diag(rowsum Dᵉ)·(Dᵉ)⁻¹, and scatter Aᵉ·Mᵉ into M_dual and
    rowsum(Dᵉ) into the DIAGONAL D_dual. Returns (D_dual_diag, M_dual). No global dense inverse."""
    Ns, Nm = len(s_nodes), len(m_nodes)
    Dd = np.zeros(Ns)                 # diagonal only
    Md = np.zeros((Ns, Nm))
    for sf in s_facets:
        nps_s = len(sf)
        Xs = np.zeros((4, 3)); Xs[:nps_s] = s_nodes[list(sf)]
        De = np.zeros((nps_s, nps_s))
        Me = np.zeros((nps_s, Nm))
        for mf in m_facets:
            nps_m = len(mf)
            Xm = np.zeros((4, 3)); Xm[:nps_m] = m_nodes[list(mf)]
            res = mortar_pair(Xs, nps_s, Xm, nps_m, refDir)
            if res is None:
                continue
            De += res["D"]
            for a in range(nps_s):
                for k in range(nps_m):
                    Me[a, mf[k]] += res["M"][a, k]
        if np.linalg.norm(De) < 1e-300:
            continue
        c = De.sum(axis=1)            # cᵉ_a = ∫_e N_a  (row-sum of the facet mass)
        Ae = np.diag(c) @ np.linalg.inv(De)
        MeDual = Ae @ Me
        for a in range(nps_s):
            I = sf[a]
            Dd[I] += c[a]
            Md[I, :] += MeDual[a, :]
    return Dd, Md


def assemble_lumped(D, M):
    """The NAIVE sparsifier for contrast: row-sum LUMP the standard D to a diagonal. Keeps
    partition of unity but is NOT biorthogonal ⇒ breaks linear completeness (T6)."""
    return D.sum(axis=1), M


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


def master_support(s_nodes, s_facets, m_nodes, m_facets, refDir):
    """For each slave NODE, the set of master nodes on facets its support overlaps (the LOCAL
    coupling the dual basis should reproduce). Used to bound sparsity in T4."""
    Ns = len(s_nodes)
    support = [set() for _ in range(Ns)]
    for sf in s_facets:
        nps_s = len(sf)
        Xs = np.zeros((4, 3)); Xs[:nps_s] = s_nodes[list(sf)]
        for mf in m_facets:
            nps_m = len(mf)
            Xm = np.zeros((4, 3)); Xm[:nps_m] = m_nodes[list(mf)]
            if mortar_pair(Xs, nps_s, Xm, nps_m, refDir) is None:
                continue
            for a in sf:
                support[a].update(mf)
    return support


# ===================================================================== tests
print("ADR-62 P2.1 oracle — DUAL (biorthogonal) basis ⇒ sparse mortar tie")
RD = np.array([0.0, 0.0, 1.0])

# genuinely NON-MATCHING flat interface: slave 3×3, master 2×2 (no coincident interior nodes)
s_nodes, s_facets = quad_mesh(3, 3, z=0.0)
m_nodes, m_facets = quad_mesh(2, 2, z=0.0)

D_std, M_std = assemble_standard(s_nodes, s_facets, m_nodes, m_facets, RD)
Dd, Md = assemble_dual(s_nodes, s_facets, m_nodes, m_facets, RD)

# --- T1: biorthogonality — D_dual diagonal, diag == standard row-sum -----------
print("\nT1  D_dual = ∫ψ N^s is DIAGONAL; diag == ∫N_I (standard D row-sum)")
# reconstruct the full dual D once (for the off-diagonal check) via ψ = Aᵉ N applied to N^s.
# Here D_dual is stored as its diagonal by construction; verify against the biorthogonality
# identity D_dual,II = Σ_J D_std,IJ (row-sum) and that the true off-diagonals vanish.
rowsum_std = D_std.sum(axis=1)
check("D_dual diagonal == standard D row-sum (∫N_I)", np.allclose(Dd, rowsum_std, atol=1e-9),
      f"max|Δ|={np.abs(Dd - rowsum_std).max():.2e}")
# explicit off-diagonal check: build D_dual_full = A(assembled) — do it per facet and confirm
# the assembled full matrix is diagonal.
Ns = len(s_nodes)
Dfull = np.zeros((Ns, Ns))
for sf in s_facets:
    nps_s = len(sf)
    Xs = np.zeros((4, 3)); Xs[:nps_s] = s_nodes[list(sf)]
    De = np.zeros((nps_s, nps_s))
    for mf in m_facets:
        nps_m = len(mf)
        Xm = np.zeros((4, 3)); Xm[:nps_m] = m_nodes[list(mf)]
        res = mortar_pair(Xs, nps_s, Xm, nps_m, RD)
        if res is not None:
            De += res["D"]
    if np.linalg.norm(De) < 1e-300:
        continue
    c = De.sum(axis=1)
    DeDual = (np.diag(c) @ np.linalg.inv(De)) @ De     # = diag(c) exactly
    for a in range(nps_s):
        for b in range(nps_s):
            Dfull[sf[a], sf[b]] += DeDual[a, b]
offdiag = Dfull - np.diag(np.diag(Dfull))
check("D_dual off-diagonals ≈ 0 (biorthogonal)", np.abs(offdiag).max() < 1e-9,
      f"max|offdiag|={np.abs(offdiag).max():.2e}")

# --- T2: partition of unity ---------------------------------------------------
print("\nT2  P_dual = D_dual⁻¹ M_dual : partition of unity")
P_dual = Md / Dd[:, None]
P_std = np.linalg.solve(D_std, M_std)
check("Σ_k P_dual,Ik == 1", np.allclose(P_dual.sum(1), 1.0, atol=1e-9),
      f"max|ΣP-1|={np.abs(P_dual.sum(1) - 1).max():.2e}")

# --- T3: linear completeness / constant-stress patch --------------------------
print("\nT3  LINEAR COMPLETENESS: P_dual·(linear on masters) == field at slave nodes")
worst = 0.0
for (a, b, c) in [(1, 0, 0), (0.3, 1.7, 0), (0.0, 0.0, 2.5), (-1.0, 0.5, -0.8)]:
    fm = a + b * m_nodes[:, 0] + c * m_nodes[:, 1]
    fs_exact = a + b * s_nodes[:, 0] + c * s_nodes[:, 1]
    worst = max(worst, np.abs(P_dual @ fm - fs_exact).max())
check("dual reproduces the linear field (constant-stress patch)", worst < 1e-9,
      f"max nodal err={worst:.2e}")

# --- T4: SPARSITY / locality --------------------------------------------------
print("\nT4  SPARSITY: dual P rows are LOCAL (⊆ facet-support masters); fewer nnz than dense P")
support = master_support(s_nodes, s_facets, m_nodes, m_facets, RD)
nnz_dual = (np.abs(P_dual) > 1e-9).sum()
nnz_std = (np.abs(P_std) > 1e-9).sum()
# every dual nonzero must lie within that slave node's overlapping-master support
within = True
for I in range(Ns):
    cols = set(np.where(np.abs(P_dual[I]) > 1e-9)[0].tolist())
    if not cols.issubset(support[I]):
        within = False
        break
check("dual P support ⊆ per-node master facet support (local)", within)
check("dual P strictly sparser than standard (dense) P", nnz_dual < nnz_std,
      f"nnz dual={nnz_dual} < std={nnz_std}")
check("standard P is dense (some row touches most masters)",
      (np.abs(P_std) > 1e-9).sum(1).max() >= len(m_nodes) - 1,
      f"max nnz/row std={int((np.abs(P_std) > 1e-9).sum(1).max())} of {len(m_nodes)}")

# --- T5: same tie on linear data ----------------------------------------------
print("\nT5  standard & dual reproduce the SAME linear field (both variationally consistent)")
lin = 0.4 + 1.3 * m_nodes[:, 0] - 0.7 * m_nodes[:, 1]
lin_s = 0.4 + 1.3 * s_nodes[:, 0] - 0.7 * s_nodes[:, 1]
es = np.abs(P_std @ lin - lin_s).max()
ed = np.abs(P_dual @ lin - lin_s).max()
check("standard reproduces the linear field", es < 1e-9, f"err={es:.2e}")
check("dual reproduces the SAME linear field", ed < 1e-9, f"err={ed:.2e}")

# --- T6: CONTRAST — row-sum LUMPING breaks linear completeness -----------------
print("\nT6  CONTRAST: row-sum LUMPED D keeps Σ=1 but BREAKS linear completeness (fails patch)")
Dl, Ml = assemble_lumped(D_std, M_std)
P_lump = Ml / Dl[:, None]
check("lumped keeps partition of unity", np.allclose(P_lump.sum(1), 1.0, atol=1e-9),
      f"max|ΣP-1|={np.abs(P_lump.sum(1) - 1).max():.2e}")
worst_l = 0.0
for (a, b, c) in [(0.3, 1.7, 0), (0.0, 0.0, 2.5), (-1.0, 0.5, -0.8)]:
    fm = a + b * m_nodes[:, 0] + c * m_nodes[:, 1]
    fs_exact = a + b * s_nodes[:, 0] + c * s_nodes[:, 1]
    worst_l = max(worst_l, np.abs(P_lump @ fm - fs_exact).max())
check("lumped FAILS the linear patch (proves dual ≠ lumping)", worst_l > 1e-4,
      f"lumped patch err={worst_l:.2e} (dual was {worst:.1e})")

# ===================================================================== summary
print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'} — "
      "P2.1 dual basis: ψ=AᵉN per-facet ⇒ D_dual DIAGONAL ⇒ P LOCAL/sparse, still linearly "
      "complete (unlike lumping); built from the same integratePair Dᵉ/Mᵉ ⇒ no kernel/handler change.")
sys.exit(0 if _fails == 0 else 1)
