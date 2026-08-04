"""ADR-78 ORACLE — QUADRATIC (serendipity) facets for the integral-mortar mesh-tie.

ADR-62 P2/P2.1 shipped the mortar tie for tri3/quad4 facets. ADR-78 (the apeGmsh
ADR 0086 D2 "both" follow-up) extends the kernel + generator to quadratic
serendipity facets: quad8 slave/master, tri6 MASTER only. The design (ADR-78 D1):

  ALL GEOMETRY runs on the CORNER polygon (aux plane, Sutherland-Hodgman clip,
  inverse isomap, closest-point projection) — exactly the shipped code paths —
  and only the BASIS EVALUATION changes: the full serendipity N/phi are evaluated
  at the mapped Gauss points. Exact (not approximate) in scope, because a FLAT
  facet with midsides AT the edge midpoints has a serendipity map that reduces
  identically to the corner map. That precondition is a GUARD (named refusal,
  kernel status -2), never a silent skip.

Gates (falsifiable, numeric):
  T1  serendipity shape sanity: Kronecker delta at nodes; partition of unity;
      reference integrals INT N = -1/3 (quad8 corner), 4/3 (quad8 mid),
      0 (tri6 corner — THE D2 refusal reason), 1/6 (tri6 mid).
  T2  the degree-6 12-point triangle rule integrates monomials to degree 6 exactly
      (weights sum to 1; area applied separately — the kernel convention).
  T3  D/M on a NON-MATCHING quad8-on-quad4 patch (2x2 quad8 slave vs 3x3 quad4
      master, coincident plane): partition of unity Sum_K M_IK == Sum_J D_IJ.
  T4  STANDARD basis P = D^-1 M: P.1 = 1; LINEAR PATCH EXACT (P transfers any
      linear field from master nodes to slave nodes to ~1e-9).
  T5  DUAL basis (per-facet condensation, the P2.1 construction): D_dual is
      DIAGONAL with NEGATIVE corner entries (-A/12 scale — legitimate, the sign
      cancels in P); P_dual.1 = 1; linear patch exact; strictly sparser than the
      standard P.
  T6  CONTRAST: row-sum LUMPING D breaks the linear patch for quad8 (as it did
      for quad4 in proto_p2_1) — dual is the correct sparsifier here too.
  T7  refusals: (a) tri6 slave — the per-facet corner row-sum c^e is ~0 => the
      dual scaling divides by zero (structural, no tolerance rescues it);
      (b) a CURVED quad8 (midside lifted off the corner plane / off the edge
      midpoint) is REFUSED by the validity guard (status -2), not mis-integrated.
  T8  the sign-free guards (ADR-78 D3): per-facet AREA coverage detects a slave
      protruding past the master; gapL1/area recovers a uniform interface gap.
  T9  order-4 vs order-6 quadrature: the linear patch is exact under BOTH (D and
      M share Gauss points); reported deltas justify degree-6 as the quadratic
      default (OQ-3).

SCOPE (review 2026-08-04): this oracle gates the KERNEL MATH — bases, quadrature,
D/M assembly, condensation, and the guard MEASURES (area, gapL1). The generator's
refusal LOGIC (thresholds, guard ordering, -2 hard-error vs skip, master-overlap
refusal) is C++-only and gated by tests/test_ladrunoTie_mortar_quad8.py.

Q4 TODO (gated on apeGmsh ADR 0086 S1 merging): run apeGmsh's numpy kernel
(_kernel/resolvers/_mortar.py) on THIS patch geometry and compare P row-by-row
to <=1e-12. The reference P is printed at the end for that cross-check.

Run: python proto_p2_2_quad8_mortar.py   (numpy only, no build)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True, linewidth=120)
TOL = 1e-9
_fails = 0


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ===================================================================== bases
def ncorner(nps):
    return nps if nps <= 4 else nps // 2


def shape_lin(nps, xi, eta):
    """Corner (geometry) basis + derivs: tri3 area coords / quad4 bilinear."""
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


def shape_full(nps, xi, eta):
    """Full basis: tri3/quad4 pass through; tri6/quad8 serendipity (corners first)."""
    if nps <= 4:
        return shape_lin(nps, xi, eta)[0]
    if nps == 8:
        xs = np.array([-1.0, 1.0, 1.0, -1.0])
        es = np.array([-1.0, -1.0, 1.0, 1.0])
        N = np.empty(8)
        N[:4] = 0.25 * (1 + xi * xs) * (1 + eta * es) * (xi * xs + eta * es - 1)
        N[4] = 0.5 * (1 - xi * xi) * (1 - eta)   # mid 0-1
        N[5] = 0.5 * (1 + xi) * (1 - eta * eta)  # mid 1-2
        N[6] = 0.5 * (1 - xi * xi) * (1 + eta)   # mid 2-3
        N[7] = 0.5 * (1 - xi) * (1 - eta * eta)  # mid 3-0
        return N
    # tri6
    l0, l1, l2 = 1 - xi - eta, xi, eta
    return np.array([l0 * (2 * l0 - 1), l1 * (2 * l1 - 1), l2 * (2 * l2 - 1),
                     4 * l0 * l1, 4 * l1 * l2, 4 * l2 * l0])


def parent_nodes(nps):
    if nps == 8:
        return [(-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1), (1, 0), (0, 1), (-1, 0)]
    if nps == 6:
        return [(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)]
    if nps == 4:
        return [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    return [(0, 0), (1, 0), (0, 1)]


# ===================================================================== 2D helpers (kernel mirrors)
def signed_area2(P):
    P = np.asarray(P, float)
    x, y = P[:, 0], P[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def ensure_ccw(P):
    return P[::-1].copy() if signed_area2(P) < 0 else P.copy()


def is_convex2(P, tol=1e-12):
    n = len(P)
    if n < 3:
        return False
    scale = max(np.linalg.norm(P[(i + 1) % n] - P[i]) for i in range(n))
    if scale < 1e-300:
        return False
    for i in range(n):
        a = P[(i + 1) % n] - P[i]
        b = P[(i + 2) % n] - P[(i + 1) % n]
        if a[0] * b[1] - a[1] * b[0] < -tol * scale * scale:
            return False
    return True


def clip_polygon(subj, clip):
    out = [tuple(p) for p in subj]
    for e in range(len(clip)):
        A, B = clip[e], clip[(e + 1) % len(clip)]
        ex, ey = B[0] - A[0], B[1] - A[1]
        inp, out = out, []
        if not inp:
            break
        S = inp[-1]
        sS = ex * (S[1] - A[1]) - ey * (S[0] - A[0])
        for E in inp:
            sE = ex * (E[1] - A[1]) - ey * (E[0] - A[0])
            if sE >= 0.0:
                if sS < 0.0:
                    out.append(_iline(S, E, A, B))
                out.append(E)
            elif sS >= 0.0:
                out.append(_iline(S, E, A, B))
            S, sS = E, sE
    return np.array(out) if out else np.zeros((0, 2))


def _iline(S, E, A, B):
    d1 = (E[0] - S[0], E[1] - S[1])
    d2 = (B[0] - A[0], B[1] - A[1])
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-300:
        return tuple(S)
    t = ((A[0] - S[0]) * d2[1] - (A[1] - S[1]) * d2[0]) / den
    return (S[0] + t * d1[0], S[1] + t * d1[1])


def dedupe(P, tol=1e-12):
    keep = []
    n = len(P)
    for i in range(n):
        if np.linalg.norm(P[i] - P[(i + 1) % n]) > tol:
            keep.append(P[i])
    return np.array(keep) if keep else np.zeros((0, 2))


# ===================================================================== quadrature
def tri_rule(order):
    """Barycentric points + weights SUMMING TO 1 (area applied separately).
    Mirrors the kernel's order<=4 (6-pt degree-4) and >4 (12-pt degree-6) rules;
    the kernel's order<=2 3-pt rule is NOT mirrored (unused on the tie path)."""
    if order <= 4:
        a1, b1, w1 = 0.445948490915965, 0.108103018168070, 0.223381589678011
        a2, b2, w2 = 0.091576213509771, 0.816847572980459, 0.109951743655322
        bary = np.array([[b1, a1, a1], [a1, b1, a1], [a1, a1, b1],
                         [b2, a2, a2], [a2, b2, a2], [a2, a2, b2]])
        w = np.array([w1, w1, w1, w2, w2, w2])
        return bary, w
    # Dunavant degree-6, 12 points, weights sum to 1
    g = []
    w = []
    for (l1, l2, wt, perm6) in [
            (0.501426509658179, 0.249286745170910, 0.116786275726379, False),
            (0.873821971016996, 0.063089014491502, 0.050844906370207, False),
            (0.053145049844817, 0.310352451033784, 0.082851075618374, True)]:
        l3 = 1.0 - l1 - l2
        pts = {(l1, l2, l3), (l2, l3, l1), (l3, l1, l2)} if not perm6 else \
              {(l1, l2, l3), (l1, l3, l2), (l2, l1, l3), (l2, l3, l1), (l3, l1, l2), (l3, l2, l1)}
        for p in sorted(pts):
            g.append(p)
            w.append(wt)
    return np.array(g), np.array(w)


# ===================================================================== projection (corner geometry)
def project_corner(npc, Xc, xs, tolR=1e-12, maxit=20, tolP=1e-8):
    xi, eta = (0.0, 0.0) if npc == 4 else (1 / 3, 1 / 3)
    for _ in range(maxit):
        N, dNx, dNe = shape_lin(npc, xi, eta)
        xbar = N @ Xc
        g1, g2 = dNx @ Xc, dNe @ Xc
        d = xs - xbar
        R0, R1 = d @ g1, d @ g2
        if np.hypot(R0, R1) < tolR:
            return xi, eta, 0
        K00, K01, K11 = g1 @ g1, g1 @ g2, g2 @ g2
        detK = K00 * K11 - K01 * K01
        if abs(detK) < 1e-14 * (K00 * K11 + 1e-300):
            return xi, eta, -1
        dxi = (K11 * R0 - K01 * R1) / detK
        deta = (-K01 * R0 + K00 * R1) / detK
        xi += dxi
        eta += deta
        # the kernel's SCALE-FREE step-converged escape (LadrunoContactProjection
        # project(), 2026-07 large-coordinate fix) — without it the oracle silently
        # drops GPs at facet scales >~2000 length units where the C++ integrates.
        if abs(dxi) + abs(deta) < tolP:
            return xi, eta, 0
    return xi, eta, -1


def inverse_isomap2d(npc, UV, q, tolR=1e-13, maxit=20):
    if npc == 3:
        A = np.column_stack([UV[1] - UV[0], UV[2] - UV[0]])
        det = np.linalg.det(A)
        if abs(det) < 1e-14 * (np.linalg.norm(A) ** 2 + 1e-300):
            return None
        s = np.linalg.solve(A, q - UV[0])
        return s[0], s[1]
    xi = eta = 0.0
    for _ in range(maxit):
        N, dNx, dNe = shape_lin(4, xi, eta)
        r = N @ UV - q
        if np.linalg.norm(r) < tolR:
            return xi, eta
        J = np.array([[dNx @ UV[:, 0], dNe @ UV[:, 0]], [dNx @ UV[:, 1], dNe @ UV[:, 1]]])
        det = np.linalg.det(J)
        if abs(det) < 1e-300:
            return None
        step = np.linalg.solve(J, r)
        xi -= step[0]
        eta -= step[1]
        if abs(step[0]) + abs(step[1]) < 1e-10:
            return xi, eta
    return None


# ===================================================================== the pair integral (ADR-78 kernel mirror)
def check_quadratic_facet(nps, X, tol=1e-6):
    """ADR-78 D1 guard: midsides at their (3D) edge midpoints. True = valid.
    Corner coplanarity is deliberately NOT checked — the serendipity map reduces
    to the corner map for midpoint midsides even with warped corners, so a warped
    quad8 degrades exactly like the shipped warped quad4 (~0.7% area bias)."""
    npc = ncorner(nps)
    if npc == nps:
        return True
    Xc, Xm = X[:npc], X[npc:]
    edges = [(i, (i + 1) % npc) for i in range(npc)]
    scale = max(np.linalg.norm(Xc[b] - Xc[a]) for a, b in edges)
    for k, (a, b) in enumerate(edges):
        if np.linalg.norm(Xm[k] - 0.5 * (Xc[a] + Xc[b])) > tol * scale:
            return False
    return True


def facet_normal(npc, Xc, refDir):
    if npc == 3:
        raw = np.cross(Xc[1] - Xc[0], Xc[2] - Xc[0])
    else:
        raw = np.cross(Xc[2] - Xc[0], Xc[3] - Xc[1])
    j = np.linalg.norm(raw)
    if j < 1e-300:
        return None
    n = raw / j
    return -n if n @ refDir < 0 else n


def integrate_pair(nps_s, Xs, nps_m, Xm, refDir, order=6):
    """Returns dict(D, M, g, area, gapL1, status). status: 0 ok, -1 skip, -2 invalid quadratic."""
    Xs, Xm = np.asarray(Xs, float), np.asarray(Xm, float)
    out = dict(D=np.zeros((nps_s, nps_s)), M=np.zeros((nps_s, nps_m)),
               g=np.zeros(nps_s), area=0.0, gapL1=0.0, status=-1)
    if not check_quadratic_facet(nps_s, Xs) or not check_quadratic_facet(nps_m, Xm):
        out["status"] = -2
        return out
    ncs, ncm = ncorner(nps_s), ncorner(nps_m)
    Xsc, Xmc = Xs[:ncs], Xm[:ncm]

    x0 = Xmc.mean(axis=0)
    n0 = facet_normal(ncm, Xmc, refDir)
    if n0 is None:
        return out
    t = Xmc[1] - Xmc[0]
    e1 = t - (t @ n0) * n0
    if np.linalg.norm(e1) < 1e-300:
        return out
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n0, e1)
    ns = facet_normal(ncs, Xsc, refDir)
    if ns is None:
        return out
    cos_t = abs(ns @ n0)
    if cos_t < 1e-12:
        return out

    UVs = np.array([[(X - x0) @ e1, (X - x0) @ e2] for X in Xsc])
    subj = ensure_ccw(UVs)
    clip = ensure_ccw(np.array([[(X - x0) @ e1, (X - x0) @ e2] for X in Xmc]))
    if not is_convex2(subj) or not is_convex2(clip):
        return out
    poly = dedupe(clip_polygon(subj, clip))
    if len(poly) < 3:
        return out
    A_min = min(abs(signed_area2(subj)), abs(signed_area2(clip)))
    if abs(signed_area2(poly)) < 1e-12 * A_min:
        return out

    cen = poly.mean(axis=0)
    bary, w = tri_rule(order)
    ngp = 0
    for s in range(len(poly)):
        tri = np.array([cen, poly[s], poly[(s + 1) % len(poly)]])
        A_aux = abs(signed_area2(tri))
        if A_aux < 1e-12 * A_min:
            continue
        A_phys = A_aux / cos_t
        for gp in range(len(w)):
            q = bary[gp] @ tri
            iso = inverse_isomap2d(ncs, UVs, q)
            if iso is None:
                continue
            xi_s, eta_s = iso
            Ns = shape_full(nps_s, xi_s, eta_s)
            x_s = Ns @ Xs
            xi_m, eta_m, st = project_corner(ncm, Xmc, x_s)
            if st < 0:
                continue
            Nc, dNx, dNe = shape_lin(ncm, xi_m, eta_m)
            g1, g2 = dNx @ Xmc, dNe @ Xmc
            raw = np.cross(g1, g2)
            j = np.linalg.norm(raw)
            if j < 1e-300:
                continue
            n = raw / j
            if n @ refDir < 0:
                n = -n
            phi = shape_full(nps_m, xi_m, eta_m)
            x_m = phi @ Xm
            gN = n @ (x_s - x_m)
            wJ = w[gp] * A_phys
            out["area"] += wJ
            out["gapL1"] += wJ * abs(gN)
            out["g"] += wJ * Ns * gN
            out["D"] += wJ * np.outer(Ns, Ns)
            out["M"] += wJ * np.outer(Ns, phi)
            ngp += 1
    out["status"] = 0 if ngp > 0 else -1
    return out


# ===================================================================== meshes
def quad8_mesh(nx, ny, Lx=1.0, Ly=1.0, z=0.0):
    """nx x ny serendipity quad8 mesh on [0,Lx]x[0,Ly] at height z. Returns (coords, facets)."""
    nodes, coords, facets = {}, [], []

    def nid(x, y):
        key = (round(x, 12), round(y, 12))
        if key not in nodes:
            nodes[key] = len(coords)
            coords.append([x, y, z])
        return nodes[key]

    hx, hy = Lx / nx, Ly / ny
    for ey in range(ny):
        for ex in range(nx):
            x0, y0 = ex * hx, ey * hy
            c = [nid(x0, y0), nid(x0 + hx, y0), nid(x0 + hx, y0 + hy), nid(x0, y0 + hy),
                 nid(x0 + hx / 2, y0), nid(x0 + hx, y0 + hy / 2),
                 nid(x0 + hx / 2, y0 + hy), nid(x0, y0 + hy / 2)]
            facets.append(c)
    return np.array(coords, float), facets


def quad4_mesh(nx, ny, Lx=1.0, Ly=1.0, z=0.0, x1=0.0):
    nodes, coords, facets = {}, [], []

    def nid(x, y):
        key = (round(x, 12), round(y, 12))
        if key not in nodes:
            nodes[key] = len(coords)
            coords.append([x, y, z])
        return nodes[key]

    hx, hy = (Lx - x1) / nx, Ly / ny
    for ey in range(ny):
        for ex in range(nx):
            x0, y0 = x1 + ex * hx, ey * hy
            facets.append([nid(x0, y0), nid(x0 + hx, y0), nid(x0 + hx, y0 + hy), nid(x0, y0 + hy)])
    return np.array(coords, float), facets


def assemble(sc, sf, nps_s, mc, mf, nps_m, refDir, order=6):
    Ns, Nm = len(sc), len(mc)
    D, M = np.zeros((Ns, Ns)), np.zeros((Ns, Nm))
    facet = []  # per-slave-facet dict for dual + guards
    for fs in sf:
        Xs = sc[fs]
        self_pr = integrate_pair(nps_s, Xs, nps_s, Xs, refDir, order)
        De, Me = np.zeros((nps_s, nps_s)), np.zeros((nps_s, Nm))
        area_cov, gapL1 = 0.0, 0.0
        for fm in mf:
            pr = integrate_pair(nps_s, Xs, nps_m, mc[fm], refDir, order)
            if pr["status"] == -2:
                raise ValueError("invalid quadratic facet (status -2)")
            if pr["status"] != 0:
                continue
            De += pr["D"]
            Me[:, fm] += pr["M"]
            area_cov += pr["area"]
            gapL1 += pr["gapL1"]
        for a, I in enumerate(fs):
            M[I, :] += Me[a]
            for b, J in enumerate(fs):
                D[I, J] += De[a, b]
        facet.append(dict(nodes=fs, De=De, Me=Me, area_full=self_pr["area"],
                          area_cov=area_cov, gapL1=gapL1))
    return D, M, facet


def condense_standard(D, M):
    return np.linalg.solve(D, M)


def condense_dual(facet, Ns, Nm):
    Dd = np.zeros(Ns)
    Md = np.zeros((Ns, Nm))
    for f in facet:
        De, Me, fs = f["De"], f["Me"], f["nodes"]
        if np.abs(De).sum() < 1e-300:
            continue
        c = De.sum(axis=1)                      # per-facet row-sums = INT_e N_a
        Y = np.linalg.solve(De, Me)
        for a, I in enumerate(fs):
            Dd[I] += c[a]
            Md[I] += c[a] * Y[a]
    if np.any(np.abs(Dd) <= 1e-300):
        raise ZeroDivisionError("zero dual mass (the tri6-slave degeneracy)")
    return Md / Dd[:, None]


# ===================================================================== gates
print("=" * 78)
print("ADR-78 oracle — quadratic serendipity mortar (quad8 slave on quad4 master)")
print("=" * 78)

print("\nT1  serendipity shape sanity")
for nps in (6, 8):
    pn = parent_nodes(nps)
    ok_delta = all(abs(shape_full(nps, *pn[j])[i] - (1.0 if i == j else 0.0)) < 1e-14
                   for i in range(nps) for j in range(nps))
    rng = np.random.default_rng(42)
    ok_pu = True
    for _ in range(50):
        if nps == 8:
            xi, eta = rng.uniform(-1, 1, 2)
        else:
            xi = rng.uniform(0, 1)
            eta = rng.uniform(0, 1 - xi)
        ok_pu &= abs(shape_full(nps, xi, eta).sum() - 1.0) < 1e-13
    check(f"nps={nps}: Kronecker delta at nodes", ok_delta)
    check(f"nps={nps}: partition of unity", ok_pu)
# reference integrals: quad8 via tensor Gauss; tri6 via the degree-6 rule
gx, gw = np.polynomial.legendre.leggauss(4)
I8 = np.zeros(8)
for x, wx in zip(gx, gw):
    for e, we in zip(gx, gw):
        I8 += wx * we * shape_full(8, x, e)
check("quad8: INT N_corner = -1/3 (=> -A/12 of area 4)", np.allclose(I8[:4], -1 / 3, atol=1e-12),
      f"got {I8[:4]}")
check("quad8: INT N_mid = 4/3 (=> A/3)", np.allclose(I8[4:], 4 / 3, atol=1e-12))
bary, w = tri_rule(6)
I6 = np.zeros(6)
for b, wt in zip(bary, w):
    I6 += 0.5 * wt * shape_full(6, b[1], b[2])
check("tri6: INT N_corner = 0 (THE tri6-slave refusal reason)", np.allclose(I6[:3], 0, atol=1e-14),
      f"got {I6[:3]}")
check("tri6: INT N_mid = 1/6", np.allclose(I6[3:], 1 / 6, atol=1e-14))

print("\nT2  degree-6 12-point triangle rule")
check("weights sum to 1", abs(w.sum() - 1.0) < 1e-14)
ok = True
import math
for p in range(7):
    for qd in range(7 - p):
        # exact INT x^p y^q over unit triangle = p! q! / (p+q+2)!
        exact = math.factorial(p) * math.factorial(qd) / math.factorial(p + qd + 2)
        num = sum(0.5 * wt * (b[1] ** p) * (b[2] ** qd) for b, wt in zip(bary, w))
        ok &= abs(num - exact) < 1e-14
check("all monomials x^p y^q, p+q<=6, exact", ok)

print("\nT3  quad8-on-quad4 non-matching patch: D/M assembly")
refDir = np.array([0.0, 0.0, 1.0])
sc, sf = quad8_mesh(2, 2)              # 2x2 quad8 slave on [0,1]^2, h=0.5
mc, mf = quad4_mesh(3, 3)              # 3x3 quad4 master, h=1/3 (non-matching)
D, M, facet = assemble(sc, sf, 8, mc, mf, 4, refDir)
rD, rM = D.sum(axis=1), M.sum(axis=1)
check("partition of unity Sum_K M_IK == Sum_J D_IJ", np.allclose(rD, rM, atol=1e-12),
      f"max diff {np.abs(rD - rM).max():.2e}")
check("corner row-sums NEGATIVE (~ -A/12 tributary)", all(
    rD[fsl[c]] < 0 for fsl in sf for c in range(4)))
check("total slave area = 1", abs(sum(f["area_cov"] for f in facet) - 1.0) < 1e-12)

print("\nT4  STANDARD basis: linear patch")
P = condense_standard(D, M)
check("P.1 = 1", np.allclose(P.sum(axis=1), 1.0, atol=1e-9),
      f"max |rowsum-1| {np.abs(P.sum(axis=1) - 1).max():.2e}")
lin = lambda X: 0.3 + 1.7 * X[:, 0] - 0.9 * X[:, 1]
err_std = np.abs(P @ lin(mc) - lin(sc)).max()
check("linear patch exact (std)", err_std < TOL, f"err {err_std:.2e}")

print("\nT5  DUAL basis: diagonal, negative corners, sparse, patch-exact")
Pd = condense_dual(facet, len(sc), len(mc))
Dd = np.zeros(len(sc))
for f in facet:
    c = f["De"].sum(axis=1)
    for a, I in enumerate(f["nodes"]):
        Dd[I] += c[a]
corner_ids = sorted({fsl[c] for fsl in sf for c in range(4)})
mid_ids = sorted({fsl[c] for fsl in sf for c in range(4, 8)})
check("dual masses NEGATIVE at corners, positive at mids",
      all(Dd[i] < 0 for i in corner_ids) and all(Dd[i] > 0 for i in mid_ids),
      f"corner range [{Dd[corner_ids].min():.4f},{Dd[corner_ids].max():.4f}]")
check("P_dual.1 = 1", np.allclose(Pd.sum(axis=1), 1.0, atol=1e-9))
err_dual = np.abs(Pd @ lin(mc) - lin(sc)).max()
check("linear patch exact (dual)", err_dual < TOL, f"err {err_dual:.2e}")
nnz = lambda A: int((np.abs(A) > 1e-12).sum())
check("dual strictly sparser than standard", nnz(Pd) < nnz(P), f"{nnz(Pd)} vs {nnz(P)}")

print("\nT6  CONTRAST: row-sum lumping breaks the quad8 patch")
Plump = M / D.sum(axis=1)[:, None]
err_lump = np.abs(Plump @ lin(mc) - lin(sc)).max()
check("lumped patch FAILS (err >> dual err)", err_lump > 1e3 * max(err_dual, 1e-15),
      f"lumped {err_lump:.2e} vs dual {err_dual:.2e}")

print("\nT7  refusals")
# (a) tri6 slave: per-facet corner row-sum ~ 0 => dual divides by zero
tc = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
               [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0]], float)
pr6 = integrate_pair(6, tc, 6, tc, refDir)
c6 = pr6["D"].sum(axis=1)
check("tri6 slave: corner c^e ~ 0 (dual scaling divides by zero)",
      np.allclose(c6[:3], 0, atol=1e-12) and np.all(c6[3:] > 0.1),
      f"corners {c6[:3]}, mids {c6[3:]}")
# (b) curved quad8 (midside lifted) => status -2, never integrated
Xcurved = sc[sf[0]].copy()
Xcurved[4] += np.array([0.0, 0.0, 0.05])
prc = integrate_pair(8, Xcurved, 4, mc[mf[0]], refDir)
check("curved quad8 REFUSED (status -2, not mis-integrated)", prc["status"] == -2)
# midside slid along the edge is equally refused
Xslid = sc[sf[0]].copy()
Xslid[4] += np.array([0.07, 0.0, 0.0])
check("midside-off-midpoint REFUSED (status -2)",
      integrate_pair(8, Xslid, 4, mc[mf[0]], refDir)["status"] == -2)

print("\nT8  sign-free guards (ADR-78 D3, unified per OQ-2)")
# protrusion: master covers only x<=0.9 => rightmost slave facets under-covered
mc9, mf9 = quad4_mesh(3, 3, Lx=0.9)
_, _, facet9 = assemble(sc, sf, 8, mc9, mf9, 4, refDir)
ratios = [f["area_cov"] / f["area_full"] for f in facet9]
check("per-facet AREA coverage detects protrusion", min(ratios) < 1 - 1e-3 < max(ratios) + 1,
      f"min ratio {min(ratios):.4f} (facets at x>0.9 uncovered)")
# uniform gap: slave lifted by 1e-3 => gapL1/area recovers it
scL, sfL = quad8_mesh(2, 2, z=1e-3)
_, _, facetL = assemble(scL, sfL, 8, mc, mf, 4, refDir)
gbar = sum(f["gapL1"] for f in facetL) / sum(f["area_cov"] for f in facetL)
check("gapL1/area recovers a uniform 1e-3 gap", abs(gbar - 1e-3) < 1e-9, f"got {gbar:.3e}")

print("\nT9  order-4 vs order-6 quadrature")
D4, M4, facet4 = assemble(sc, sf, 8, mc, mf, 4, refDir, order=4)
P4 = condense_standard(D4, M4)
Pd4 = condense_dual(facet4, len(sc), len(mc))
e4s = np.abs(P4 @ lin(mc) - lin(sc)).max()
e4d = np.abs(Pd4 @ lin(mc) - lin(sc)).max()
check("patch exact under order 4 too (quadrature-independent)", e4s < TOL and e4d < TOL,
      f"std {e4s:.2e} dual {e4d:.2e}")
dD = np.abs(D - D4).max()
print(f"       (D drift order4->order6: {dD:.3e} — the degree-4 rule under-integrates "
      f"the degree-6 serendipity products; order 6 is the quadratic default per OQ-3)")

print("\nT10 cross-order patches: tri6 MASTER and quad8-on-quad8")
# (a) quad4 slave on TWO tri6 masters covering the unit square — exercises the
# quadratic-master phi path (serendipity evaluated at the corner-map projection),
# which T3-T6 (quad4 masters) never touch. Linear patch must be exact.
t6c = {}
t6coords = []
def _t6(x, y):
    key = (round(x, 12), round(y, 12))
    if key not in t6c:
        t6c[key] = len(t6coords)
        t6coords.append([x, y, 0.0])
    return t6c[key]
tri_a = [_t6(0, 0), _t6(1, 0), _t6(1, 1),
         _t6(0.5, 0), _t6(1, 0.5), _t6(0.5, 0.5)]
tri_b = [_t6(0, 0), _t6(1, 1), _t6(0, 1),
         _t6(0.5, 0.5), _t6(0.5, 1), _t6(0, 0.5)]
mc6 = np.array(t6coords, float)
sc4, sf4 = quad4_mesh(2, 2)                # 2x2 quad4 slave, non-matching with the tri split
D6, M6, facet6 = assemble(sc4, sf4, 4, mc6, [tri_a, tri_b], 6, refDir)
P6 = condense_standard(D6, M6)
Pd6 = condense_dual(facet6, len(sc4), len(mc6))
e6s = np.abs(P6 @ lin(mc6) - lin(sc4)).max()
e6d = np.abs(Pd6 @ lin(mc6) - lin(sc4)).max()
check("quad4-on-TRI6-MASTER linear patch exact (std + dual)", e6s < TOL and e6d < TOL,
      f"std {e6s:.2e} dual {e6d:.2e}")
# (b) quad8 slave on a NON-MATCHING quad8 master mesh — both bases quadratic.
mc8, mf8 = quad8_mesh(3, 3)
D88, M88, facet88 = assemble(sc, sf, 8, mc8, mf8, 8, refDir)
P88 = condense_standard(D88, M88)
Pd88 = condense_dual(facet88, len(sc), len(mc8))
e88s = np.abs(P88 @ lin(mc8) - lin(sc)).max()
e88d = np.abs(Pd88 @ lin(mc8) - lin(sc)).max()
check("quad8-on-QUAD8-MASTER linear patch exact (std + dual)", e88s < TOL and e88d < TOL,
      f"std {e88s:.2e} dual {e88d:.2e}")
# (c) tri3 slave on quad8 master — the remaining mixed pairing.
sc3 = np.array([[0.2, 0.2, 0], [0.8, 0.2, 0], [0.2, 0.8, 0]], float)
D38, M38, facet38 = assemble(sc3, [[0, 1, 2]], 3, mc8, mf8, 8, refDir)
P38 = condense_standard(D38, M38)
e38 = np.abs(P38 @ lin(mc8) - lin(sc3)).max()
check("tri3-on-QUAD8-MASTER linear patch exact", e38 < TOL, f"err {e38:.2e}")

print("\nT11 the L1 gap guard is strictly stronger than the signed average")
# An accordion slave (two facets tilted oppositely, ridge on the master plane):
# the SIGNED weighted gap of the shared ridge nodes cancels across facets (the
# pre-ADR-78 per-node guard signal there is ~0), while the per-facet L1 measure
# sees the full |gap| everywhere. Quantifies the OQ-2 quirks-ledger claim.
d = 5e-3
acc = np.array([[0, 0, -d], [0.5, 0, 0], [0.5, 1, 0], [0, 1, -d],      # facet A (tilted down)
                [1, 0, d], [1, 1, d]], float)                          # facet B (tilted UP)
facA = [0, 1, 2, 3]
facB = [1, 4, 5, 2]                          # shares the ridge nodes 1, 2 (z = 0)
mcA, mfA = quad4_mesh(1, 1)                  # one flat master covering [0,1]^2
gsigned = np.zeros(len(acc))
L1_over_area = []
for fs in (facA, facB):
    pr = integrate_pair(4, acc[fs], 4, mcA[mfA[0]], refDir, order=6)
    for a, I in enumerate(fs):
        gsigned[I] += pr["g"][a]
    L1_over_area.append(pr["gapL1"] / pr["area"])
ridge_cancel = max(abs(gsigned[1]), abs(gsigned[2]))
check("ridge nodes' SIGNED gap cancels (old per-node signal ~0)",
      ridge_cancel < 1e-4 * d, f"|signed| {ridge_cancel:.2e} vs d {d:.0e}")
check("per-facet L1 gap sees the tilt (guard fires at tol < ~d/2)",
      min(L1_over_area) > 0.2 * d, f"gapL1/area {L1_over_area}")

# ---------------------------------------------------------------- cross-check reference
print("\nQ4 cross-check reference (vs apeGmsh _kernel/resolvers/_mortar.py, once merged):")
print(f"  slave nodes {len(sc)}, master nodes {len(mc)}; ||P_dual||_F = {np.linalg.norm(Pd):.12f}, "
      f"||P_std||_F = {np.linalg.norm(P):.12f}")
print(f"  P_dual row 0 (slave node 0 at {sc[0][:2]}): "
      f"{[(k, round(Pd[0, k], 9)) for k in np.nonzero(np.abs(Pd[0]) > 1e-12)[0]]}")

print("\n" + "=" * 78)
print(f"{'ALL GATES PASS' if _fails == 0 else f'{_fails} GATE(S) FAILED'}")
print("=" * 78)
sys.exit(0 if _fails == 0 else 1)
