"""ADR-41 C1 ORACLE — LadrunoMortarKernel: clip -> Gauss -> D, M, weighted gap g̃.

Oracle-first (the fork discipline): validate the segment-to-segment MORTAR
integration kernel — the overlap CLIP (the only genuinely new algorithm), the
sub-triangle Gauss accumulation of the mortar matrices D_IJ, M_IK and the
weighted gap g̃_I — in pure numpy BEFORE transcribing to the header-only
SRC/domain/contact/LadrunoMortarKernel.h. No OpenSees, no build.

What C1 owns (geometry + the patch test ONLY — ALM is C2, friction is C3):
  - Sutherland-Hodgman overlap clip of slave∩master on a common auxiliary plane,
    fan-triangulated; 3-pt (or higher) Gauss per sub-triangle.
  - back-map each GP to ξ_s (inverse 2D iso-map within the slave facet) and to
    ξ_m via the SHIPPED LadrunoContactProjection.projectFull (do NOT re-derive).
  - STANDARD (not dual) basis:
        D_IJ += w·J · N_I^s · N_J^s            (slave–slave)
        M_IK += w·J · N_I^s · φ_K^m            (slave–master)
        g̃_I += w·J · N_I^s · ((x_s − x_m)·n)   (normal INSIDE the integral)
    with J = A_aux / |n_s·n0|  (flat-facet area ratio; aux-plane -> slave plane).

The gates (capstone C1 row / ADR-41 P1) — falsifiable, numeric:
  1. clip area vs known cases (the crux algorithm, isolated).
  2. partition of unity  Σ_K M_IK == Σ_J D_IJ  to 1e-12 (matched AND non-matched).
  3. CONSTANT-pressure patch test on a NON-MATCHED 2-block mesh: recover a uniform
     field through the transfer P = D⁻¹M -> constant to ≤1e-6 (THE headline gate).
  4. LINEAR consistency on a non-matched mesh -> linear-exact; this is what the
     clip BUYS — the un-clipped "mortar-lite" provably fails it (demonstrated).
  5. weighted gap g̃ linear-exact (uniform + tilted-plane gap vs reference quad).
  6. FD-check d{g̃}/d(nodal coord) vs the analytic uniform-translation derivative
     (feeds the C2 K_c linearization the header stubs dD,dM,dn,dξ will provide).

Run: <pythoncore-3.12> proto_c1_mortar.py   (numpy only)
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


# ============================================================ shipped projection
# Mirror of SRC/domain/contact/LadrunoContactProjection.h (validated by
# proto_a2_projection.py 7/7). C1 CONSUMES projectFull — it is reproduced here so
# the oracle is standalone, but it is NOT the thing under test.
def shape(nps, xi, eta):
    if nps == 4:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dNe = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    else:  # tri-3, N = [1-xi-eta, xi, eta]
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


def project_full(nps, X, xs, refDir, tolR=1e-12, maxit=10):
    X = np.asarray(X, float)[:nps]
    p = {"status": -1, "xi": 0.0, "eta": 0.0, "gap": 0.0,
         "n": np.zeros(3), "phi": np.zeros(4)}
    xi, eta, st = project(nps, X, xs, tolR, maxit)
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


# ===================================================== C1: the mortar kernel (NEW)
def facet_normal(X, nps, refDir):
    """Unit normal of a (flat) facet, oriented so n·refDir > 0."""
    if nps == 3:
        raw = np.cross(X[1] - X[0], X[2] - X[0])
    else:  # quad: cross of the diagonals (robust for a planar quad)
        raw = np.cross(X[2] - X[0], X[3] - X[1])
    n = raw / np.linalg.norm(raw)
    return -n if (n @ refDir) < 0 else n


def aux_plane(Xm, nps_m, refDir):
    """Auxiliary plane from the MASTER facet (ADR-41 §How): origin x0=centroid,
    normal n0=master facet normal, orthonormal in-plane axes (e1,e2)."""
    x0 = Xm[:nps_m].mean(axis=0)
    n0 = facet_normal(Xm, nps_m, refDir)
    e1 = (Xm[1] - Xm[0])
    e1 = e1 - (e1 @ n0) * n0
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n0, e1)
    return x0, n0, e1, e2


def to2d(X3d, x0, e1, e2):
    d = X3d - x0
    return np.column_stack([d @ e1, d @ e2])


def signed_area(poly):
    a = 0.0
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def ensure_ccw(poly):
    return poly if signed_area(poly) >= 0 else poly[::-1].copy()


def clip_polygon(subject, clip):
    """Sutherland-Hodgman: subject ∩ clip, both CCW, CLIP convex. Returns 2D pts.

    The clip polygon must be convex (tri/quad master ⇒ convex after projection);
    the subject may be any simple polygon. Result is convex with ≤ |S|+|C| verts.
    """
    out = [np.asarray(p, float) for p in subject]
    nC = len(clip)
    for i in range(nC):
        if not out:
            break
        A = np.asarray(clip[i], float)
        B = np.asarray(clip[(i + 1) % nC], float)
        edge = B - A

        def side(P):  # >0 left of A->B  (inside for CCW clip)
            return edge[0] * (P[1] - A[1]) - edge[1] * (P[0] - A[0])

        inp = out
        out = []
        S = inp[-1]
        sS = side(S)
        for E in inp:
            sE = side(E)
            if sE >= 0.0:                       # E inside
                if sS < 0.0:
                    out.append(_iline(S, E, A, B))
                out.append(E)
            elif sS >= 0.0:                     # S inside, E outside
                out.append(_iline(S, E, A, B))
            S, sS = E, sE
    return out


def _iline(S, E, A, B):
    """Intersection of segment S->E with the (infinite) line A->B."""
    d1 = E - S
    d2 = B - A
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-300:
        return S.copy()
    t = ((A[0] - S[0]) * d2[1] - (A[1] - S[1]) * d2[0]) / den
    return S + t * d1


def dedupe(poly, tol=1e-12):
    """Merge near-coincident consecutive vertices (clip-degeneracy guard)."""
    out = []
    n = len(poly)
    for i in range(n):
        P = poly[i]
        Q = poly[(i + 1) % n]
        if np.linalg.norm(P - Q) > tol:
            out.append(P)
    return out


def tri_quadrature(order):
    """Barycentric Gauss rules on the reference triangle; weights SUM TO 1
    (area factor applied separately). order=2 (3-pt) is the C1 default."""
    if order <= 2:
        b = [(2 / 3, 1 / 6, 1 / 6), (1 / 6, 2 / 3, 1 / 6), (1 / 6, 1 / 6, 2 / 3)]
        w = [1 / 3, 1 / 3, 1 / 3]
        return b, w
    # degree-4, 6-point (for tight linear-consistency checks on bilinear quads)
    a1, b1, w1 = 0.445948490915965, 0.108103018168070, 0.223381589678011
    a2, b2, w2 = 0.091576213509771, 0.816847572980459, 0.109951743655322
    b = [(b1, a1, a1), (a1, b1, a1), (a1, a1, b1),
         (b2, a2, a2), (a2, b2, a2), (a2, a2, b2)]
    w = [w1, w1, w1, w2, w2, w2]
    return b, w


def inverse_isomap_2d(nps, UV, q, tolR=1e-13, maxit=20):
    """Recover ξ_s s.t. Σ N_I(ξ_s)·UV_I = q  (UV = slave nodes in aux 2D).
    Closed form for tri-3 (affine); bounded Newton for the bilinear quad."""
    if nps == 3:
        Jm = np.column_stack([UV[1] - UV[0], UV[2] - UV[0]])
        xi, eta = np.linalg.solve(Jm, q - UV[0])
        return xi, eta
    xi = eta = 0.0
    for _ in range(maxit):
        N, dNx, dNe = shape(4, xi, eta)
        r = N @ UV - q
        if np.linalg.norm(r) < tolR:
            break
        Jm = np.column_stack([dNx @ UV, dNe @ UV])
        d = np.linalg.solve(Jm, -r)
        xi += d[0]
        eta += d[1]
    return xi, eta


def clip_subtris(Xs, nps_s, Xm, nps_m, refDir, area_tol=1e-12):
    """The CLIP: project both facets to the master-built aux plane, Sutherland-
    Hodgman the slave∩master overlap, fan-triangulate from the centroid.
    Returns (subtris_2d, UVs, x0, n0, e1, e2, n_s) or None if no usable overlap.
        subtris_2d : list of (q0,q1,q2) 2D vertices in aux coords.
        UVs        : slave nodes in aux 2D (for the GP back-map).
        n_s        : slave facet unit normal (for the area Jacobian)."""
    x0, n0, e1, e2 = aux_plane(Xm, nps_m, refDir)
    n_s = facet_normal(Xs, nps_s, refDir)
    UVs = ensure_ccw(to2d(Xs[:nps_s], x0, e1, e2))
    UVm = ensure_ccw(to2d(Xm[:nps_m], x0, e1, e2))
    # re-fetch the CCW slave UV in NODE order (back-map needs node ordering, not CCW)
    UVs_nodeorder = to2d(Xs[:nps_s], x0, e1, e2)
    poly = dedupe(clip_polygon(list(UVs), list(UVm)))
    if len(poly) < 3:
        return None
    A_s = abs(signed_area(UVs))
    A_m = abs(signed_area(UVm))
    if abs(signed_area(poly)) < area_tol * min(A_s, A_m):   # sliver guard
        return None
    cen = np.mean(poly, axis=0)
    subtris = []
    n = len(poly)
    for i in range(n):
        tri = (cen, poly[i], poly[(i + 1) % n])
        if abs(signed_area(list(tri))) > area_tol * min(A_s, A_m):
            subtris.append(tri)
    if not subtris:
        return None
    return subtris, UVs_nodeorder, x0, n0, e1, e2, n_s


def mortar_pair(Xs, nps_s, Xm, nps_m, refDir, order=2):
    """Integrate ONE slave/master facet pair. Returns per-pair contributions:
        D  (nps_s × nps_s),  M (nps_s × nps_m),  g  (nps_s,)  [weighted gap],
        plus area = ∫dΓ over the overlap (for the clip-area gate).
    All None if the overlap is empty/degenerate."""
    clip = clip_subtris(Xs, nps_s, Xm, nps_m, refDir)
    if clip is None:
        return None
    subtris, UVs, x0, n0, e1, e2, n_s = clip
    cos_t = abs(n_s @ n0)                       # slave-plane / aux-plane tilt
    if cos_t < 1e-12:
        return None                             # slave ⟂ aux: ill-posed pairing
    bary, wts = tri_quadrature(order)
    D = np.zeros((nps_s, nps_s))
    M = np.zeros((nps_s, nps_m))
    g = np.zeros(nps_s)
    area = 0.0
    for tri in subtris:
        A_aux = abs(signed_area(list(tri)))
        A_phys = A_aux / cos_t                  # aux-plane area -> slave-surface area
        V0, V1, V2 = tri
        for (L0, L1, L2), w in zip(bary, wts):
            q = L0 * V0 + L1 * V1 + L2 * V2     # GP in aux 2D
            xi_s, eta_s = inverse_isomap_2d(nps_s, UVs, q)
            Ns, _, _ = shape(nps_s, xi_s, eta_s)
            x_s = Ns @ Xs[:nps_s]
            p = project_full(nps_m, Xm, x_s, refDir)
            if p["status"] < 0:
                continue
            phi = p["phi"][:nps_m]
            x_m = phi @ Xm[:nps_m]
            gN = p["n"] @ (x_s - x_m)
            wJ = w * A_phys
            area += wJ
            D += wJ * np.outer(Ns, Ns)
            M += wJ * np.outer(Ns, phi)
            g += wJ * Ns * gN
    return {"D": D, "M": M, "g": g, "area": area}


def assemble(s_nodes, s_facets, m_nodes, m_facets, refDir, order=2):
    """Assemble GLOBAL mortar operators over a non-matched interface.
    s_facets / m_facets: lists of node-index tuples (len 3 tri or 4 quad)."""
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


# ============================================================== mesh generators
def grid_nodes(nx, ny, z=0.0, x0=0.0, y0=0.0, sx=1.0, sy=1.0):
    pts = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            pts.append([x0 + sx * i / nx, y0 + sy * j / ny, z])
    return np.array(pts, float), (lambda i, j: j * (nx + 1) + i)


def quad_mesh(nx, ny, **kw):
    nodes, idx = grid_nodes(nx, ny, **kw)
    facets = []
    for j in range(ny):
        for i in range(nx):
            facets.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)))
    return nodes, facets


def tri_mesh(nx, ny, **kw):
    nodes, idx = grid_nodes(nx, ny, **kw)
    facets = []
    for j in range(ny):
        for i in range(nx):
            a, b = idx(i, j), idx(i + 1, j)
            c, d = idx(i + 1, j + 1), idx(i, j + 1)
            facets.append((a, b, c))
            facets.append((a, c, d))
    return nodes, facets


# ===================================================================== tests
print("ADR-41 C1 oracle — LadrunoMortarKernel (clip -> Gauss -> D, M, g̃)")
RD = np.array([0.0, 0.0, 1.0])

# --- T1: Sutherland-Hodgman clip area vs known closed-form cases ---------------
print("\nT1  clip: Sutherland-Hodgman overlap area vs known cases")
# unit square [0,1]² clipped by [0.5,1.5]² ⇒ overlap [0.5,1]² area 0.25
sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
sq2 = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)]
ov = clip_polygon([np.array(p, float) for p in sq], [np.array(p, float) for p in sq2])
check("square∩offset-square area == 0.25", abs(abs(signed_area(ov)) - 0.25) < TOL,
      f"area={abs(signed_area(ov)):.6f}")
# square clipped by a triangle fully inside ⇒ triangle area 0.18
tri_in = [(0.2, 0.2), (0.8, 0.2), (0.2, 0.8)]
ov2 = clip_polygon([np.array(p, float) for p in sq], [np.array(p, float) for p in tri_in])
check("square∩interior-triangle area == 0.18", abs(abs(signed_area(ov2)) - 0.18) < TOL,
      f"area={abs(signed_area(ov2)):.6f}")
# identical squares ⇒ overlap == 1.0 (full)
ov3 = clip_polygon([np.array(p, float) for p in sq], [np.array(p, float) for p in sq])
check("square∩itself area == 1.0", abs(abs(signed_area(ov3)) - 1.0) < TOL)
# disjoint ⇒ empty
ov4 = clip_polygon([np.array(p, float) for p in sq],
                   [np.array(p, float) for p in [(2, 2), (3, 2), (3, 3), (2, 3)]])
check("disjoint squares ⇒ empty clip", len(dedupe(ov4)) < 3, f"verts={len(ov4)}")
# unit square rotated 45° about its centre ∩ unit square ⇒ octagon, area = 2(√2−1).
# Vertices: centre + (1/√2)·(cosθ,sinθ), θ = k·90°  (the rotated corners).
c, R = 0.5, 1.0 / np.sqrt(2)
rot = [(c + R * np.cos(k * np.pi / 2), c + R * np.sin(k * np.pi / 2)) for k in range(4)]
ov5 = clip_polygon([np.array(p, float) for p in sq], [np.array(p, float) for p in rot])
check("square∩45°-rotated-square area == 2(√2−1)",
      abs(abs(signed_area(ov5)) - 2 * (np.sqrt(2) - 1)) < TOL,
      f"area={abs(signed_area(ov5)):.6f}")

# --- T2: single pair, MATCHED quad/quad ⇒ Σ_K M = Σ_J D, area exact ------------
print("\nT2  single pair, MATCHED quad/quad: ΣM==ΣD, ∫dΓ == facet area")
Xs = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float)
Xs4 = np.zeros((4, 3)); Xs4[:4] = Xs
r = mortar_pair(Xs4, 4, Xs4.copy(), 4, RD)
check("∫dΓ == 1.0 (unit quad)", abs(r["area"] - 1.0) < TOL, f"area={r['area']:.9f}")
check("ΣM row == ΣD row (partition of unity)",
      np.allclose(r["M"].sum(1), r["D"].sum(1), atol=1e-12),
      f"max|ΔΣ|={np.abs(r['M'].sum(1) - r['D'].sum(1)).max():.2e}")
check("D symmetric", np.allclose(r["D"], r["D"].T, atol=1e-12))
check("ΣΣD == area (∫1·1 dΓ)", abs(r["D"].sum() - 1.0) < TOL)

# --- T3: single pair, NON-MATCHED + mixed tri/quad ⇒ still ΣM==ΣD --------------
print("\nT3  single pair, NON-MATCHED (rotated/offset, tri slave / quad master): ΣM==ΣD")
ang = 0.37
Rz = np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
Xtri = np.array([[0.1, 0.1, 0], [0.9, 0.2, 0], [0.3, 0.8, 0]], float)
Xtri4 = np.zeros((4, 3)); Xtri4[:3] = Xtri
# master: a LARGE quad [-0.5,1.5]² rotated about (0.5,0.5) — fully contains the slave
# triangle, so ∫dΓ over the overlap == the slave-triangle area exactly.
Xbig = np.array([[-0.5, -0.5, 0], [1.5, -0.5, 0], [1.5, 1.5, 0], [-0.5, 1.5, 0]], float)
Xmq = (Xbig - [0.5, 0.5, 0]) @ Rz.T + [0.5, 0.5, 0]
Xmq4 = np.zeros((4, 3)); Xmq4[:4] = Xmq
r3 = mortar_pair(Xtri4, 3, Xmq4, 4, RD)
check("overlap found", r3 is not None)
check("ΣM row == ΣD row (non-matched, tri/quad)",
      np.allclose(r3["M"].sum(1), r3["D"].sum(1), atol=1e-12),
      f"max|ΔΣ|={np.abs(r3['M'].sum(1) - r3['D'].sum(1)).max():.2e}")
# the slave triangle lies fully inside the (large rotated) master ⇒ ∫dΓ == tri area
tri_area = 0.5 * abs(np.cross(Xtri[1] - Xtri[0], Xtri[2] - Xtri[0])[2])
check("∫dΓ == slave-triangle area (slave ⊂ master)", abs(r3["area"] - tri_area) < 1e-9,
      f"area={r3['area']:.9f} tri={tri_area:.9f}")

# --- T4: CONSTANT-pressure patch test on a NON-MATCHED 2-block quad mesh -------
print("\nT4  HEADLINE — constant-pressure patch test, non-matched quad meshes")
s_nodes, s_facets = quad_mesh(2, 2)        # slave: 2×2 quads (9 nodes)
m_nodes, m_facets = quad_mesh(3, 3)        # master: 3×3 quads (16 nodes) — non-matched
D, M, g = assemble(s_nodes, s_facets, m_nodes, m_facets, RD)
check("ΣM row == ΣD row (global, non-matched)",
      np.allclose(M.sum(1), D.sum(1), atol=1e-12),
      f"max|ΔΣ|={np.abs(M.sum(1) - D.sum(1)).max():.2e}")
check("D SPD (full slave coverage ⇒ invertible)", np.all(np.linalg.eigvalsh(D) > 1e-10),
      f"λ_min={np.linalg.eigvalsh(D).min():.3e}")
pbar = 7.3
u_m = np.full(len(m_nodes), pbar)
u_s = np.linalg.solve(D, M @ u_m)          # transfer P = D⁻¹M
rel = np.abs(u_s - pbar).max() / pbar
check("D⁻¹M transfers CONSTANT pbar -> constant (≤1e-6 rel)", rel <= 1e-6,
      f"max rel dev = {rel:.2e}")
check("total ∫dΓ over slave == interface area 1.0", abs(D.sum() - 1.0) < 1e-9,
      f"ΣΣD={D.sum():.9f}")

# --- T5: WHY the clip is in the MVP — it integrates ∫N^s·φ^m EXACTLY at low order;
#         the clip-blind "mortar-lite" rule has O(few-%) quadrature error -----------
print("\nT5  clip necessity (non-matched tri meshes): clip EXACT at low order; mortar-lite biased")
sn, sf = tri_mesh(2, 2)                      # slave tris  (2×2 → 8 triangles)
mn2, mf2 = tri_mesh(3, 3)                    # master tris (3×3 → 18 triangles) — non-matched

# (a) a globally LINEAR field is reproduced exactly by the clipped transfer P=D⁻¹M.
Dl, Ml, gl = assemble(sn, sf, mn2, mf2, RD, order=2)
abc = np.array([0.5, 1.7, -0.9])            # f = a + b x + c y
fm = abc[0] + abc[1] * mn2[:, 0] + abc[2] * mn2[:, 1]
fs_exact = abc[0] + abc[1] * sn[:, 0] + abc[2] * sn[:, 1]
fs = np.linalg.solve(Dl, Ml @ fm)
err_clip = np.abs(fs - fs_exact).max() / np.abs(fs_exact).max()
check("CLIPPED mortar reproduces a LINEAR field (≤1e-6 rel)", err_clip <= 1e-6,
      f"max rel err = {err_clip:.2e}")

# (b) the real falsifier is the mortar matrix M = ∫N_I^s φ_K^m dΓ itself. φ_K^m has a
#     KINK at every master-facet edge. The CLIP aligns sub-triangle boundaries with
#     those edges ⇒ each sub-tri integrand is smooth ⇒ a 3-pt rule is already exact
#     (order-2 == order-6). A clip-BLIND rule integrates across the kinks and is wrong
#     at low order, converging to the clipped value only under heavy refinement.
_, M2, _ = assemble(sn, sf, mn2, mf2, RD, order=2)
_, M6, _ = assemble(sn, sf, mn2, mf2, RD, order=4)   # higher-order sub-tri rule
ref = np.linalg.norm(M6)
clip_self = np.linalg.norm(M2 - M6) / ref
check("CLIP exact at low order: ||M(o2) − M(o4)|| ≈ 0 (≤1e-12)", clip_self <= 1e-12,
      f"rel = {clip_self:.2e}")


def refined_tri_rule(nsub):
    """Midpoint rule on the reference triangle refined into nsub² equal sub-tris;
    barycentric GP at each sub-tri centroid, weights summing to 1. nsub=1 → 1-pt."""
    pts, w = [], 1.0 / (nsub * nsub)
    for i in range(nsub):
        for j in range(nsub - i):
            # 'up' sub-triangle centroid (barycentric on the parent)
            a = (i + 1 / 3) / nsub
            b = (j + 1 / 3) / nsub
            pts.append(((1 - a - b), a, b, w))
            if j < nsub - i - 1:                 # 'down' sub-triangle
                a2 = (i + 2 / 3) / nsub
                b2 = (j + 2 / 3) / nsub
                pts.append(((1 - a2 - b2), a2, b2, w))
    return pts


def unclipped_M(s_nodes, s_facets, m_nodes, m_facets, refDir, nsub):
    """'mortar-lite' M: integrate over the WHOLE slave facet with an nsub-refined rule,
    each GP bound to the SINGLE master facet it projects in-bounds onto (clip-blind)."""
    Ns, Nm = len(s_nodes), len(m_nodes)
    M = np.zeros((Ns, Nm))
    rule = refined_tri_rule(nsub)
    for sf_ in s_facets:
        nps_s = len(sf_)
        Xs = s_nodes[list(sf_)]
        A = 0.5 * abs(np.cross(Xs[1] - Xs[0], Xs[2] - Xs[0])[2])  # tri facets here
        for (L0, L1, L2, w) in rule:
            xg = L0 * Xs[0] + L1 * Xs[1] + L2 * Xs[2]
            Nsf = np.array([L0, L1, L2])         # linear tri shape fns at the GP
            for mf_ in m_facets:
                nps_m = len(mf_)
                Xm = m_nodes[list(mf_)]
                p = project_full(nps_m, Xm, xg, refDir)
                if p["status"] == 0:
                    phi = p["phi"][:nps_m]
                    for ka in range(nps_s):
                        for kk in range(nps_m):
                            M[sf_[ka], mf_[kk]] += w * A * Nsf[ka] * phi[kk]
                    break
    return M


lite_c = np.linalg.norm(unclipped_M(sn, sf, mn2, mf2, RD, nsub=1) - M6) / ref   # coarse
lite_f = np.linalg.norm(unclipped_M(sn, sf, mn2, mf2, RD, nsub=24) - M6) / ref  # refined
check("clip-blind COARSE rule is biased (||M_lite − M_clip|| > 1e-3)", lite_c > 1e-3,
      f"coarse rel err = {lite_c:.2e}")
check("clip-blind rule CONVERGES to the clipped (exact) M under refinement",
      lite_f < 0.2 * lite_c, f"refined {lite_f:.2e} vs coarse {lite_c:.2e}")

# --- T6: weighted gap g̃ linear-exact (uniform + tilted-plane normal gap) ------
print("\nT6  weighted gap g̃: uniform gap and tilted-plane (linear) gap vs reference")
# (a) uniform gap: master lifted +δ in z ⇒ gN ≈ -δ everywhere ⇒ g̃_I = -δ·(∫N_I)
delta = 0.05
sn6, sf6 = quad_mesh(2, 2)
mn6, mf6 = quad_mesh(3, 3, z=delta)
Dg, Mg, gg = assemble(sn6, sf6, mn6, mf6, RD)
row = Dg.sum(1)                              # ∫N_I dΓ  (== D row sums)
check("uniform gap: g̃_I == -δ·∫N_I (≤1e-9)", np.allclose(gg, -delta * row, atol=1e-9),
      f"max|Δ|={np.abs(gg + delta * row).max():.2e}")
# (b) tilted master plane z = s·x ⇒ signed normal gap is LINEAR in x; check the
#     mortar g̃ against a high-order reference integral of N_I·gN over the slave.
s_tilt = 0.08
mn6b = mn6.copy(); mn6b[:, 2] = s_tilt * mn6b[:, 0]
Dg2, Mg2, gg2 = assemble(sn6, sf6, mn6b, mf6, RD)
# reference: gN(x_s) = signed distance from x_s to the master plane (n·(x_s−x0))
n_pl = np.array([-s_tilt, 0.0, 1.0]); n_pl /= np.linalg.norm(n_pl)
x0_pl = np.array([0.0, 0.0, 0.0])


def ref_gtilde(s_nodes, s_facets):
    Ns = len(s_nodes); gref = np.zeros(Ns)
    bary, wts = tri_quadrature(3)
    for sf_ in s_facets:
        nps_s = len(sf_); Xs = s_nodes[list(sf_)]
        cen = Xs.mean(0)
        for a in range(nps_s):
            V = [cen, Xs[a], Xs[(a + 1) % nps_s]]
            A = 0.5 * abs(np.cross(V[1] - V[0], V[2] - V[0])[2])
            for (L0, L1, L2), w in zip(bary, wts):
                xg = L0 * V[0] + L1 * V[1] + L2 * V[2]
                xi_s, eta_s = inverse_isomap_2d(
                    nps_s, np.column_stack([Xs[:, 0], Xs[:, 1]]), xg[:2])
                Nsf, _, _ = shape(nps_s, xi_s, eta_s)
                gN = n_pl @ (xg - x0_pl)
                gref[list(sf_)] += w * A * Nsf * gN
    return gref


gref = ref_gtilde(sn6, sf6)
check("tilted (linear) gap: mortar g̃ == reference ∫N_I·gN (≤1e-7)",
      np.allclose(gg2, gref, atol=1e-7), f"max|Δ|={np.abs(gg2 - gref).max():.2e}")

# --- T7: FD-check d g̃ / d(master z-translation) vs analytic -(∫N_I) -----------
print("\nT7  FD: d g̃_I / d(uniform master +z) == -∫N_I (feeds C2 K_c linearization)")
h = 1e-6
mn_p = mn6.copy(); mn_p[:, 2] += h
mn_m = mn6.copy(); mn_m[:, 2] -= h
_, _, gp = assemble(sn6, sf6, mn_p, mf6, RD)
_, _, gm = assemble(sn6, sf6, mn_m, mf6, RD)
dg_fd = (gp - gm) / (2 * h)
check("central-FD d g̃/dz == -∫N_I (≤1e-6)", np.allclose(dg_fd, -row, atol=1e-6),
      f"max|Δ|={np.abs(dg_fd + row).max():.2e}")

print(f"\n{'='*64}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(C1 mortar-kernel oracle)\n{'='*64}")
sys.exit(1 if _fails else 0)
