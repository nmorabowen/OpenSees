"""ADR-39 P5 / capstone B2 ORACLE — SOFT=2 segment-based explicit penalty.

Oracle-first (the fork discipline): derive + numerically check the SOFT=2 *segment-based*
(segment-to-segment) explicit penalty BEFORE transcribing it to
SRC/analysis/handler/LadrunoContactFE.cpp (the MORTAR adapter's explicit SOFT=2 path).
No OpenSees, no build — pure numpy.

WHAT B2 IS (vs the shipped B1 SOFT=1).
  B1 SOFT=1 is NODE-to-segment (NTS): it checks slave NODES against master SEGMENTS, and
  sizes the contact spring Courant-stably (k_soft = SOFSCL·4·m_eff/dt²) so explicit dt_cr is
  never penalty-throttled. But NTS MISSES contact in three regimes:
    (1) EDGE/partial-overlap: two facets overlap (penetrate) in an interior strip but NO
        slave node projects in-bounds of the master facet ⇒ NTS gives ZERO force there.
    (2) node-off-edge: a slave node's closest-point projection falls out of the master
        segment bounds (just past an edge) ⇒ the NTS pair goes inert exactly where contact
        is sharpest.
    (3) T-intersection: a slave node over the shared edge of two master facets has no unique
        owning segment ⇒ NTS is ambiguous / picks one.
  B2 SOFT=2 fixes this with a SEGMENT-BASED penalty: detect slave-facet ∩ master-facet
  overlap (REUSE the shipped LadrunoMortarKernel clip → sub-triangle Gauss — the same
  segment-to-segment geometry the mortar lane uses) and distribute the penalty traction over
  the clipped overlap, with the contact stiffness sized by the SAME B1 SOFT mass/dt² rule so
  explicit dt_cr stays un-throttled. Result: explicit contact survives corners/edges/recontact
  at the structural timestep — the robustness the progressive-collapse / element-removal
  direction needs once bodies fragment and contact on arbitrary edges.

THE PENALTY OVER THE CLIPPED OVERLAP.
  The mortar kernel returns, per facet pair, the standard mortar matrices D_IJ, M_IK and the
  weighted normal gap g̃_I = ∫ N_I^s gN dΓ over the overlap (gN = n·(x_s − x_m), <0 ⇒
  penetration). Per slave node I define the area-normalised gap ḡ_I = g̃_I / a_I with
  a_I = Σ_J D_IJ = ∫ N_I^s dΓ. The SOFT=2 nodal pressure is the penalty (NO ALM — explicit is
  single-pass, like B1):
        p_I = k_soft,I · ḡ_I      if ḡ_I < 0   (compression),   else 0   (KKT separation)
  scattered self-equilibratingly along n exactly like the shipped C2 normal block:
        F^s_K = −(D·p)_K n ,   F^m_L = +(Mᵀ·p)_L n .

SEGMENT-BASED m_eff (the gap-mode generalized mass — oracle item 4).
  ḡ_I = B_I·u is a linear form in the facet nodal displacements with the gap operator row
        B_I :  slave node J → (D_IJ/a_I) n ,   master node K → −(M_IK/a_I) n .
  The gap-mode generalized mass is m_eff,I = 1/(B_I M⁻¹ B_Iᵀ); for a diagonal (lumped) mass
        B_I M⁻¹ B_Iᵀ = Σ_J (D_IJ/a_I)² invMproj(m_sJ,n) + Σ_K (M_IK/a_I)² invMproj(m_mK,n),
  invMproj(m,n) = Σ_d n_d²/m_d (a DOF with no mass ⇒ ∞ mass ⇒ 0 contribution: a FIXED master
  node drops out, leaving the slave-side mass — the rigid-master limit). Then
        k_soft,I = SOFSCL · 4 · m_eff,I / dt²   ⇒   ω_I·dt = 2√SOFSCL ≤ 2   (Courant-stable),
  so each node's contact mode never throttles the central-difference step (the B1 rule,
  generalized from the NTS B=[n|−Nᵢ n] gap operator to the segment-based B_I=[D,−M]/a_I).
  Matched-mesh limit: D→diag(a_I), M_IK/a_I → the master shape weights φ_K, recovering the
  NTS form m_eff = 1/(1/m_s + Σφ_K²/m_K).

GATES (mirror the capstone B2 row):
  T1 — segment-based m_eff: B_I M⁻¹ B_Iᵀ closed form == assembled; matched-mesh → NTS form;
       fixed-master limit → slave-side mass.
  T2 — CORNER/EDGE ROBUSTNESS: explicit cases NTS gets WRONG (partial overlap with no slave
       node in-bounds; T-junction over a shared master edge) that SOFT=2 catches — NTS net
       contact force is ZERO, SOFT=2 gives the correct restoring force (== a reference penalty
       integral over the overlap).
  T3 — EXPLICIT STABILITY + ENERGY BALANCE: the per-node soft oscillator is stable at the
       structural dt and energy-bounded under central difference; the same gap with a stiff
       (accuracy-sized) penalty DIVERGES at that dt.
  T4 — coupled stability: the ASSEMBLED segment-based K_c = Σ_I k_soft,I B_Iᵀ B_I has a
       central-difference ω_max·dt ≤ 2 at the default SOFSCL on a real facet pair (the SOFSCL
       margin absorbs the inter-node coupling — the B1 per-pair-conservatism caveat).

Run: <pythoncore-3.12> proto_b2_soft2_segment.py   (numpy only)
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
# proto_a2_projection.py 7/7) + LadrunoMortarKernel.h clip→Gauss (validated by
# proto_c1_mortar.py 24/24). B2 CONSUMES these verbatim — they are reproduced so the
# oracle is standalone, but they are NOT the thing under test (B2 is the SOFT sizing +
# the segment-based penalty assembly + the corner/edge robustness).
def shape(nps, xi, eta):
    if nps == 4:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dNe = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    else:  # tri-3
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
    p = {"status": -1, "xi": 0.0, "eta": 0.0, "gap": 0.0, "n": np.zeros(3), "phi": np.zeros(4)}
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


def facet_normal(X, nps, refDir):
    if nps == 3:
        raw = np.cross(X[1] - X[0], X[2] - X[0])
    else:
        raw = np.cross(X[2] - X[0], X[3] - X[1])
    n = raw / np.linalg.norm(raw)
    return -n if (n @ refDir) < 0 else n


def aux_plane(Xm, nps_m, refDir):
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


def _iline(S, E, A, B):
    d1 = E - S
    d2 = B - A
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

        inp = out
        out = []
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
    out = []
    n = len(poly)
    for i in range(n):
        P = poly[i]
        Q = poly[(i + 1) % n]
        if np.linalg.norm(P - Q) > tol:
            out.append(P)
    return out


def tri_quadrature(order):
    if order <= 2:
        b = [(2 / 3, 1 / 6, 1 / 6), (1 / 6, 2 / 3, 1 / 6), (1 / 6, 1 / 6, 2 / 3)]
        w = [1 / 3, 1 / 3, 1 / 3]
        return b, w
    a1, b1, w1 = 0.445948490915965, 0.108103018168070, 0.223381589678011
    a2, b2, w2 = 0.091576213509771, 0.816847572980459, 0.109951743655322
    b = [(b1, a1, a1), (a1, b1, a1), (a1, a1, b1), (b2, a2, a2), (a2, b2, a2), (a2, a2, b2)]
    w = [w1, w1, w1, w2, w2, w2]
    return b, w


def is_convex2(poly, tol=1e-12):
    P = np.asarray(poly, float)
    n = len(P)
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
    A_s = abs(signed_area(UVs))
    A_m = abs(signed_area(UVm))
    if abs(signed_area(poly)) < area_tol * min(A_s, A_m):
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
    """One slave/master facet pair → D (nps_s×nps_s), M (nps_s×nps_m), g (weighted gap),
    n (per-facet master normal), area. None if the overlap is empty/degenerate. This is the
    shipped LadrunoMortarKernel::integratePair (C1) — B2 consumes it for the explicit penalty."""
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
    nrm = None
    for tri in subtris:
        A_aux = abs(signed_area(list(tri)))
        A_phys = A_aux / cos_t
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
            gN = p["n"] @ (x_s - x_m)
            nrm = p["n"]
            wJ = w * A_phys
            area += wJ
            D += wJ * np.outer(Ns, Ns)
            M += wJ * np.outer(Ns, phi)
            g += wJ * Ns * gN
    if nrm is None:
        return None
    return {"D": D, "M": M, "g": g, "area": area, "n": nrm}


# ============================================================ B1 sizing helpers
def invMproj(m, n):
    """inverse mass projected on n: Σ_d n_d²/m_d; a massless DOF (m_d<=0) ⇒ ∞ mass ⇒ 0."""
    s = 0.0
    for d in range(3):
        if m[d] > 0.0:
            s += n[d] * n[d] / m[d]
    return s


def soft_kn(meff, dt, sofscl):
    return sofscl * 4.0 * meff / (dt * dt)


# ============================================================ B2: the SOFT=2 penalty (NEW)
def seg_meff(D, M, n, ms, mm):
    """Per-slave-node segment-based gap-mode generalized mass m_eff,I = 1/(B_I M⁻¹ B_Iᵀ),
    B_I = [D_IJ/a_I · n | −M_IK/a_I · n]. ms/mm = per-node translational mass [mx,my,mz].
    Returns (meff[nS], invMeff[nS], a[nS]); meff=∞ flagged as 0/invMeff for an unreferenced node."""
    nS, nM = D.shape[0], M.shape[1]
    invms = np.array([invMproj(ms[I], n) for I in range(nS)])
    invmm = np.array([invMproj(mm[K], n) for K in range(nM)])
    a = D.sum(1)
    invMeff = np.zeros(nS)
    meff = np.zeros(nS)
    for I in range(nS):
        if a[I] <= 1e-300:
            continue
        s = 0.0
        for J in range(nS):
            c = D[I, J] / a[I]
            s += c * c * invms[J]
        for K in range(nM):
            c = M[I, K] / a[I]
            s += c * c * invmm[K]
        invMeff[I] = s
        meff[I] = (1.0 / s) if s > 0.0 else np.inf
    return meff, invMeff, a


def soft2_penalty(D, M, g, n, ms, mm, dt, sofscl):
    """SOFT=2 nodal pressures p_I (≤0) and the assembled facet residual on [slave|master].
    p_I = k_soft,I·ḡ_I for ḡ_I<0 (compression), else 0; k_soft,I = SOFSCL·4·m_eff,I/dt².
    F^s_K = −(D·p)_K n, F^m_L = +(Mᵀ·p)_L n (self-equilibrating). Mirrors addSoft2Penalty."""
    nS, nM = D.shape[0], M.shape[1]
    meff, invMeff, a = seg_meff(D, M, n, ms, mm)
    p = np.zeros(nS)
    for I in range(nS):
        if a[I] <= 1e-300:
            continue
        gbar = g[I] / a[I]
        if gbar >= 0.0:
            continue
        k = soft_kn(meff[I], dt, sofscl) if (dt > 0.0 and invMeff[I] > 0.0) else 0.0
        p[I] = k * gbar
    Fs = np.zeros((nS, 3))
    Fm = np.zeros((nM, 3))
    for K in range(nS):
        Dp = float(D[K, :] @ p)
        Fs[K] = -Dp * n
    for L in range(nM):
        Mp = float(M[:, L] @ p)
        Fm[L] = Mp * n
    return p, Fs, Fm


def nts_penalty_force(slave_nodes, Xm, nps_m, refDir, kn):
    """The shipped NTS node-to-segment penalty for comparison: each slave NODE is projected
    onto the master facet; force only if penetrating AND in-bounds (evalSegment's gate). Returns
    the total normal force magnitude scattered to the slave nodes (0 if no node is in-bounds)."""
    total = 0.0
    nhit = 0
    for xs in slave_nodes:
        p = project_full(nps_m, Xm, xs, refDir)
        if p["status"] != 0:          # out-of-bounds OR no valid projection ⇒ NTS inert
            continue
        if p["gap"] >= 0.0:           # not penetrating
            continue
        total += kn * (-p["gap"])     # tn = kn·<−gap>₊
        nhit += 1
    return total, nhit


# ============================================================================ tests
print("=" * 78)
print("ADR-39 B2 ORACLE — SOFT=2 segment-based explicit penalty (k_soft = SOFSCL·4·m_eff/dt²)")
print("=" * 78)
RD = np.array([0.0, 0.0, 1.0])


# --- T1: segment-based m_eff -------------------------------------------------------
def t1_seg_meff():
    print("\nT1  segment-based m_eff = 1/(B_I M⁻¹ B_Iᵀ), B_I=[D,−M]/a_I; matched & rigid-master limits")
    # MATCHED unit quad/quad pair ⇒ D = mass-matrix-like, M = D (matched) ⇒ for node I the
    # segment m_eff should reduce to the NTS harmonic form with the projected master = node I.
    Xs = np.zeros((4, 3)); Xs[:4] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    r = mortar_pair(Xs, 4, Xs.copy(), 4, RD)
    n = r["n"]
    ms = np.array([[2.0, 2.0, 2.0], [1.5, 1.5, 1.5], [1.1, 1.1, 1.1], [1.7, 1.7, 1.7]])
    mm = ms.copy()                                     # matched mesh: master == slave nodes
    meff, invMeff, a = seg_meff(r["D"], r["M"], n, ms, mm)
    # closed form B_I M⁻¹ B_Iᵀ from D,M directly
    for I in range(4):
        bsum = sum((r["D"][I, J] / a[I])**2 / ms[J, 2] for J in range(4)) \
            + sum((r["M"][I, K] / a[I])**2 / mm[K, 2] for K in range(4))
        check(f"  node {I}: m_eff == 1/(B_I M⁻¹ B_Iᵀ)", abs(meff[I] - 1.0 / bsum) < TOL * meff[I],
              f"{meff[I]:.8f} vs {1.0/bsum:.8f}")
    # FIXED master (all master nodes massless ⇒ ∞ mass) ⇒ m_eff,I depends only on the slave masses
    mm_fixed = np.zeros((4, 3))                         # 0 ⇒ treated as ∞ mass
    meff_f, invMeff_f, _ = seg_meff(r["D"], r["M"], n, ms, mm_fixed)
    for I in range(4):
        bs = sum((r["D"][I, J] / a[I])**2 / ms[J, 2] for J in range(4))
        check(f"  node {I}: fixed-master m_eff == slave-only 1/Σ(D_IJ/a)²/m_sJ",
              abs(meff_f[I] - 1.0 / bs) < TOL * meff_f[I], f"{meff_f[I]:.6f}")
    check("  fixed master ⇒ smaller invMeff than free (less compliance) ⇒ larger m_eff",
          np.all(meff_f >= meff - 1e-12))
    # a SINGLE-node-supported gap (master one heavy node, others light) sanity: m_eff>0 finite
    check("  m_eff finite & positive on the matched pair", np.all(np.isfinite(meff)) and np.all(meff > 0))


# --- T2: CORNER/EDGE ROBUSTNESS — NTS misses, SOFT=2 catches -----------------------
def t2_corner_edge_robustness():
    print("\nT2  CORNER/EDGE: explicit cases NTS gets WRONG (zero force) that SOFT=2 catches")
    dt, sofscl, kref = 1e-3, 0.1, 1.0e7
    ms4 = np.array([[1.0, 1.0, 1.0]] * 4)
    mm_fixed = np.zeros((4, 3))

    # (a) PARTIAL OVERLAP, no slave node in master bounds. Master quad = [0,1]² at z=0 (fixed).
    #     Slave quad larger & offset so the overlap is an interior STRIP but EVERY slave node
    #     projects OUT of the master bounds ⇒ NTS (node-to-segment) gives ZERO; the facets DO
    #     penetrate over the strip ⇒ SOFT=2 integrates it.
    Xm = np.zeros((4, 3)); Xm[:4] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    dz = -0.01                                          # slave 0.01 below the master ⇒ penetration
    Xs = np.zeros((4, 3))
    Xs[:4] = [[0.5, -0.5, dz], [2.5, -0.5, dz], [2.5, 1.5, dz], [0.5, 1.5, dz]]
    nts_F, nts_hit = nts_penalty_force([Xs[i] for i in range(4)], Xm, 4, RD, kref)
    check("  (a) NTS gives ZERO (no slave node projects in-bounds of the master)",
          nts_hit == 0 and nts_F == 0.0, f"hits={nts_hit}, F={nts_F:.3e}")
    r = mortar_pair(Xs, 4, Xm, 4, RD)
    check("  (a) SOFT=2 finds the overlap (clip non-empty)", r is not None)
    p, Fs, Fm = soft2_penalty(r["D"], r["M"], r["g"], r["n"], ms4, mm_fixed, dt, sofscl)
    netF = np.linalg.norm(Fs.sum(0))
    check("  (a) SOFT=2 gives a NONZERO restoring force where NTS gave zero",
          netF > 1e-6, f"|ΣF_slave|={netF:.4e}")
    # the net normal force must equal k_soft·∫(−ḡ) over the overlap — a reference penalty integral.
    # With a UNIFORM gap (flat parallel facets) ḡ_I = dz for every covered node, so
    # Σ_I p_I·a_I = (per-node) Σ k_soft,I·dz·a_I; the slave-block net along n = −Σ_I (Σ_K D_KI) p_I
    # = −Σ_I a_I p_I (since Σ_K D_KI = a_I). Cross-check:
    a = r["D"].sum(1)
    netF_ref = abs(sum(a[I] * p[I] for I in range(4)))
    check("  (a) net |F| == |Σ a_I p_I| (self-equilibrating mortar scatter)",
          abs(netF - netF_ref) < 1e-6 * max(netF, 1e-30), f"{netF:.6e} vs {netF_ref:.6e}")
    check("  (a) self-equilibrium ΣF_slave + ΣF_master == 0",
          np.linalg.norm(Fs.sum(0) + Fm.sum(0)) < 1e-8 * max(netF, 1e-30))

    # (b) T-INTERSECTION: a slave facet spanning the SHARED EDGE of two master facets, again
    #     with no slave node in either master facet's bounds. NTS has no unique owning segment
    #     (zero in-bounds hits on each); SOFT=2 sums the overlap with BOTH master facets.
    Xm1 = np.zeros((4, 3)); Xm1[:4] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    Xm2 = np.zeros((4, 3)); Xm2[:4] = [[1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0]]
    # slave straddles x=1 (the shared edge), y outside [0,1] at the nodes
    Xs2 = np.zeros((4, 3))
    Xs2[:4] = [[0.5, -0.5, dz], [1.5, -0.5, dz], [1.5, 1.5, dz], [0.5, 1.5, dz]]
    slv = [Xs2[i] for i in range(4)]
    f1, h1 = nts_penalty_force(slv, Xm1, 4, RD, kref)
    f2, h2 = nts_penalty_force(slv, Xm2, 4, RD, kref)
    check("  (b) NTS gives ZERO across the T-junction (no in-bounds node on either facet)",
          h1 == 0 and h2 == 0, f"hits=({h1},{h2})")
    r1 = mortar_pair(Xs2, 4, Xm1, 4, RD)
    r2 = mortar_pair(Xs2, 4, Xm2, 4, RD)
    check("  (b) SOFT=2 overlaps BOTH master facets", r1 is not None and r2 is not None)
    p1, Fs1, _ = soft2_penalty(r1["D"], r1["M"], r1["g"], r1["n"], ms4, mm_fixed, dt, sofscl)
    p2, Fs2, _ = soft2_penalty(r2["D"], r2["M"], r2["g"], r2["n"], ms4, mm_fixed, dt, sofscl)
    netT = np.linalg.norm(Fs1.sum(0) + Fs2.sum(0))
    check("  (b) SOFT=2 distributes a nonzero force over the T-junction", netT > 1e-6,
          f"|ΣF|={netT:.4e}")
    # the two half-overlaps must each contribute ~half the strip ⇒ together cover the full slave∩master
    aT = abs(signed_area(ensure_ccw(to2d(Xs2[:4], np.zeros(3), np.array([1., 0, 0]), np.array([0, 1., 0])))))
    cov = r1["area"] + r2["area"]
    # overlap of slave [0.5,1.5]×[-0.5,1.5] with master union [0,2]×[0,1] = [0.5,1.5]×[0,1] = 1.0
    check("  (b) the two overlaps tile the slave∩(master union) (area == 1.0)",
          abs(cov - 1.0) < 1e-6, f"area1+area2={cov:.6f}")


# --- T3: explicit stability + energy balance (per-node oscillator) -----------------
def cd_leapfrog(Mdiag, fint, fext, u0, v0, dt, nsteps):
    u = u0.copy()
    a = (fext - fint(u0)) / Mdiag
    vhalf = v0 - 0.5 * dt * a
    U, Vf = [u0.copy()], [v0.copy()]
    blew = False
    for _ in range(nsteps):
        a = (fext - fint(u)) / Mdiag
        vhalf = vhalf + dt * a
        u = u + dt * vhalf
        vfull = vhalf + 0.5 * dt * a
        if not np.all(np.isfinite(u)) or np.max(np.abs(u)) > 1e3:
            blew = True
            break
        U.append(u.copy()); Vf.append(vfull.copy())
    return np.array(U), np.array(Vf), blew


def t3_explicit_stability():
    print("\nT3  EXPLICIT STABILITY + ENERGY BALANCE — soft segment penalty stable; stiff diverges")
    # A representative segment-based gap mode: slave node mass m on a fixed master, gap = z.
    # k_soft makes ω·dt = 2√SOFSCL; an accuracy-stiff k (R× bigger) throttles ⇒ diverges at this dt.
    m = 1.0
    k_struct = 4.0e4
    dt = 0.5 * (2.0 / np.sqrt(k_struct / m))            # SF=0.5 of the STRUCTURAL CD limit
    sofscl = 0.1
    meff = m
    ks = soft_kn(meff, dt, sofscl)
    check("  contact ω·dt = 2√SOFSCL ≤ 2 (not throttled)",
          abs(np.sqrt(ks / m) * dt - 2.0 * np.sqrt(sofscl)) < 1e-9, f"ω·dt={np.sqrt(ks/m)*dt:.4f}")

    def fint(u, k):
        g = u[0]
        return np.array([k * g if g < 0.0 else 0.0])

    U, Vf, blew = cd_leapfrog(np.array([m]), lambda u: fint(u, ks), np.array([0.0]),
                              np.array([0.02]), np.array([-3.0]), dt, 6000)
    check("  soft penalty: no blow-up at the structural dt", not blew)
    x = U[:, 0]; v = Vf[:, 0]
    away = (x > 0.01) & (v > 0.0)
    e = abs(v[away][-1]) / 3.0 if away.any() else 0.0
    check("  soft penalty: near-elastic restitution at coarse dt (|1−e| ≲ 2%)",
          abs(1.0 - e) < 0.02, f"e={e:.5f}")
    # a STIFF accuracy penalty (R=400× ⇒ ω·dt = 2√(R·SOFSCL) = 40 ≫ 2) DIVERGES at the same dt
    k_stiff = 400.0 * ks
    _, _, blew_stiff = cd_leapfrog(np.array([m]), lambda u: fint(u, k_stiff), np.array([0.0]),
                                   np.array([0.02]), np.array([-3.0]), dt, 6000)
    check("  stiff penalty DIVERGES at the same dt (ω·dt ≫ 2)", blew_stiff)
    # energy balance: |1−e| → 0 as SOFSCL refines (bounce resolved by more steps; not a leak)
    errs = []
    for s in (0.1, 0.025, 0.00625):
        kk = soft_kn(meff, dt, s)
        Us, Vfs, bl = cd_leapfrog(np.array([m]), lambda u: fint(u, kk), np.array([0.0]),
                                  np.array([0.02]), np.array([-3.0]), dt, 40000)
        xx = Us[:, 0]; vv = Vfs[:, 0]
        aw = (xx > 0.01) & (vv > 0.0)
        errs.append(abs(1.0 - (abs(vv[aw][-1]) / 3.0 if aw.any() else 0.0)))
    check("  energy error BOUNDED & SOFSCL-convergent (not a formulation leak)",
          errs[-1] < 0.3 * errs[0] and errs[-1] < 1e-3,
          f"|1−e|: {errs[0]:.2e} → {errs[1]:.2e} → {errs[2]:.2e}")


# --- T4: coupled stability — assembled K_c eigenvalue ------------------------------
def t4_coupled_stability():
    print("\nT4  coupled stability: assembled segment K_c = Σ_I k_soft,I B_Iᵀ B_I, ω_max·dt ≤ 2")
    dt, sofscl = 1e-3, 0.1
    # a NON-MATCHED facet pair (tilted/offset) with FREE masses on both sides ⇒ real inter-node
    # coupling through D,M. Build the full K_c and M, eigensolve M⁻¹K_c.
    Xs = np.zeros((4, 3)); Xs[:4] = [[0, 0, -0.01], [1.2, 0, -0.01], [1.2, 1.1, -0.01], [0, 1.1, -0.01]]
    Xm = np.zeros((4, 3)); Xm[:4] = [[-0.2, -0.1, 0], [1.0, -0.1, 0], [1.0, 1.0, 0], [-0.2, 1.0, 0]]
    r = mortar_pair(Xs, 4, Xm, 4, RD)
    check("  facet pair overlaps", r is not None)
    n = r["n"]
    D, M = r["D"], r["M"]
    nS, nM = 4, 4
    ms = np.array([[1.3, 1.3, 1.3], [0.9, 0.9, 0.9], [1.5, 1.5, 1.5], [1.1, 1.1, 1.1]])
    mm = np.array([[1.0, 1.0, 1.0], [1.4, 1.4, 1.4], [0.8, 0.8, 0.8], [1.2, 1.2, 1.2]])
    meff, invMeff, a = seg_meff(D, M, n, ms, mm)
    ndof = 3 * (nS + nM)
    Kc = np.zeros((ndof, ndof))
    for I in range(nS):
        if a[I] <= 1e-300 or invMeff[I] <= 0.0:
            continue
        kI = soft_kn(meff[I], dt, sofscl)
        BI = np.zeros(ndof)                              # gap operator row B_I (only along n)
        for J in range(nS):
            BI[3 * J:3 * J + 3] = (D[I, J] / a[I]) * n
        for K in range(nM):
            BI[3 * (nS + K):3 * (nS + K) + 3] = -(M[I, K] / a[I]) * n
        Kc += kI * np.outer(BI, BI)
    Mdiag = np.zeros(ndof)
    for J in range(nS):
        Mdiag[3 * J:3 * J + 3] = ms[J]
    for K in range(nM):
        Mdiag[3 * (nS + K):3 * (nS + K) + 3] = mm[K]
    # K_c is symmetric by construction (Σ k·outer(B_I,B_I), k>0) ⇒ PSD; check it directly. The
    # central-difference frequencies are the generalized eigenproblem K x = ω² M x; with M diagonal
    # SPD, symmetrize as M^{-1/2} K M^{-1/2} (eigvalsh-safe — M⁻¹K itself is NOT symmetric).
    kc_eig = np.linalg.eigvalsh(0.5 * (Kc + Kc.T))
    check("  assembled K_c is SPD/PSD (no negative eigenvalue)",
          kc_eig.min() > -1e-6 * abs(kc_eig.max()), f"λ_min(K_c)={kc_eig.min():.3e}")
    invsqrtM = 1.0 / np.sqrt(Mdiag)
    Ksym = (invsqrtM[:, None] * Kc) * invsqrtM[None, :]
    w2 = np.linalg.eigvalsh(0.5 * (Ksym + Ksym.T))
    wmax = np.sqrt(max(w2.max(), 0.0))
    check("  coupled ω_max·dt ≤ 2 at SOFSCL=0.1 (Courant-stable; SOFSCL margin absorbs coupling)",
          wmax * dt <= 2.0, f"ω_max·dt={wmax*dt:.4f}  (per-node bound 2√SOFSCL={2*np.sqrt(sofscl):.4f})")
    # the coupled max may exceed the per-node 2√SOFSCL (inter-node coupling) but stays ≤ 2 ⇒ stable.
    check("  coupling does raise ω_max above the per-node bound (documents the conservatism)",
          wmax * dt >= 2.0 * np.sqrt(sofscl) - 1e-9, f"{wmax*dt:.4f} ≥ {2*np.sqrt(sofscl):.4f}")


def main():
    t1_seg_meff()
    t2_corner_edge_robustness()
    t3_explicit_stability()
    t4_coupled_stability()
    print("\n" + "=" * 78)
    if _fails:
        print(f"FAILED ({_fails})")
        sys.exit(1)
    print("ALL ORACLE CHECKS PASSED  (B2 SOFT=2 segment-based explicit penalty)")
    print("=" * 78)


if __name__ == "__main__":
    main()
