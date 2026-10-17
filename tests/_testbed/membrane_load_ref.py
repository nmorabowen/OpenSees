"""Pure Python/NumPy/SciPy reference oracle for ADR-98:
1. Exact partial trapezoidal beam fixed-end forces (2D/3D, Euler-Bernoulli,
   Timoshenko shear factor Phi, end releases 0..3) and ForceBeamColumn3d
   section force response sp(x).
2. 4-Stage PSLG preprocessor + Constrained Delaunay Triangulation (CDT) +
   Chew-Ruppert quality refinement with acute-angle (<60 deg) guard.
3. Four surface-to-frame load distribution engines:
   - 'tributary': CDT-based Euclidean Voronoi / Nearest-Supported-Feature
   - 'oneway': Bidirectional ray-cast strip partitioning
   - 'fea_poisson': 1-DOF/node CST Poisson solver with element-edge diffusive
     flux apportionment (never 50/50 corner split)
   - 'fea_plate': 3-DOF/node Batoz-Bathe-Ho (1980) DKT plate bending solver
     with support rank check (knife-edge w=0 for rank>=2, clamped for rank<2)
4. Moment-conserving trapezoid projection (F_seg, M_seg) -> (w_a, w_b) per
   sub-segment and per parent frame.
5. Guyan static condensation of in-plane CST diaphragm stiffness K and mass M.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

# ============================================================================
# 1. Exact Partial Trapezoidal Fixed-End Forces & Section Response
# ============================================================================

# 3-point Gauss-Legendre on [0, 1] (exact for polynomials up to degree 5)
_GL3_U = (
    0.5 * (1.0 - math.sqrt(0.6)),
    0.5,
    0.5 * (1.0 + math.sqrt(0.6)),
)
_GL3_W = (5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0)


def _apply_end_release(m1_fixed: float, m2_fixed: float, release: int) -> tuple[float, float]:
    """Apply OpenSees moment end release (0=none, 1=node I, 2=node J, 3=both)."""
    if release == 1:
        return 0.0, m2_fixed - 0.5 * m1_fixed
    if release == 2:
        return m1_fixed - 0.5 * m2_fixed, 0.0
    if release == 3:
        return 0.0, 0.0
    return m1_fixed, m2_fixed


def partial_trapezoid_fef_2d(
    L: float,
    wya: float,
    wyb: float,
    waa: float,
    wab: float,
    aL: float,
    bL: float,
    release: int = 0,
    phi: float = 0.0,
) -> dict[str, np.ndarray]:
    """Exact 2D beam fixed-end forces for a partial trapezoidal load on [aL, bL].

    Returns dict with:
      - 'p0': basic simply-supported reactions [p0_0, p0_1, p0_2] = [-P, -Vy1_ss, -Vy2_ss]
      - 'q0': basic fixed-end forces [q0_0, q0_1, q0_2] = [-PJ, Mz1, Mz2]
      - 'ql0': local 6-DOF fixed-end forces [N1, Vy1, Mz1, N2, Vy2, Mz2]
    """
    p0 = np.zeros(3, dtype=float)
    q0 = np.zeros(3, dtype=float)
    ql0 = np.zeros(6, dtype=float)
    dxi = bL - aL
    if L <= 0.0 or dxi <= 0.0:
        return {"p0": p0, "q0": q0, "ql0": ql0}

    vy1_ss = 0.0
    vy2_ss = 0.0
    mz1_fix = 0.0
    mz2_fix = 0.0
    p_ax = 0.0
    pj_ax = 0.0
    inv_1p_phi = 1.0 / (1.0 + phi)

    for ug, wg in zip(_GL3_U, _GL3_W):
        xi = aL + dxi * ug
        wy = wya + (wyb - wya) * ug
        wa = waa + (wab - waa) * ug
        ds = L * dxi * wg

        vy1_ss += wy * (1.0 - xi) * ds
        vy2_ss += wy * xi * ds
        mz1_fix -= L * wy * xi * (1.0 - xi) * (1.0 - xi + 0.5 * phi) * inv_1p_phi * ds
        mz2_fix += L * wy * xi * (1.0 - xi) * (xi + 0.5 * phi) * inv_1p_phi * ds

        p_ax += wa * ds
        pj_ax += wa * xi * ds

    mz1, mz2 = _apply_end_release(mz1_fix, mz2_fix, release)

    p0[0] = -p_ax
    p0[1] = -vy1_ss
    p0[2] = -vy2_ss

    q0[0] = -pj_ax
    q0[1] = mz1
    q0[2] = mz2

    v_mom = (mz1 + mz2) / L
    ql0[0] = -(p_ax - pj_ax)
    ql0[1] = -vy1_ss + v_mom
    ql0[2] = mz1
    ql0[3] = -pj_ax
    ql0[4] = -vy2_ss - v_mom
    ql0[5] = mz2

    return {"p0": p0, "q0": q0, "ql0": ql0}


def partial_trapezoid_fef_3d(
    L: float,
    wya: float,
    wyb: float,
    wza: float,
    wzb: float,
    waa: float,
    wab: float,
    aL: float,
    bL: float,
    releasez: int = 0,
    releasey: int = 0,
    phiz: float = 0.0,
    phiy: float = 0.0,
) -> dict[str, np.ndarray]:
    """Exact 3D beam fixed-end forces for a partial trapezoidal load on [aL, bL].

    Returns dict with:
      - 'p0': basic reactions [p0_0..p0_4] = [-P, -Vy1_ss, -Vy2_ss, -Vz1_ss, -Vz2_ss]
      - 'q0': basic fixed-end forces [q0_0..q0_4] = [-PJ, Mz1, Mz2, My1, My2]
      - 'ql0': local 12-DOF fixed-end forces [N1, Vy1, Vz1, Tx1, My1, Mz1, N2, Vy2, Vz2, Tx2, My2, Mz2]
    """
    p0 = np.zeros(5, dtype=float)
    q0 = np.zeros(5, dtype=float)
    ql0 = np.zeros(12, dtype=float)
    dxi = bL - aL
    if L <= 0.0 or dxi <= 0.0:
        return {"p0": p0, "q0": q0, "ql0": ql0}

    f2d_y = partial_trapezoid_fef_2d(L, wya, wyb, waa, wab, aL, bL, release=releasez, phi=phiz)
    # For local z loading, +wz produces +My1 at node I and -My2 at node J
    f2d_z = partial_trapezoid_fef_2d(L, wza, wzb, 0.0, 0.0, aL, bL, release=0, phi=phiy)
    my1_fix = -f2d_z["q0"][1]
    my2_fix = -f2d_z["q0"][2]
    my1, my2 = _apply_end_release(my1_fix, my2_fix, releasey)

    p_ax = -f2d_y["p0"][0]
    pj_ax = -f2d_y["q0"][0]
    vy1_ss = -f2d_y["p0"][1]
    vy2_ss = -f2d_y["p0"][2]
    vz1_ss = -f2d_z["p0"][1]
    vz2_ss = -f2d_z["p0"][2]
    mz1 = f2d_y["q0"][1]
    mz2 = f2d_y["q0"][2]

    p0[0] = -p_ax
    p0[1] = -vy1_ss
    p0[2] = -vy2_ss
    p0[3] = -vz1_ss
    p0[4] = -vz2_ss

    q0[0] = -pj_ax
    q0[1] = mz1
    q0[2] = mz2
    q0[3] = my1
    q0[4] = my2

    vy_mom = (mz1 + mz2) / L
    vz_mom = -(my1 + my2) / L

    ql0[0] = -(p_ax - pj_ax)
    ql0[1] = -vy1_ss + vy_mom
    ql0[2] = -vz1_ss + vz_mom
    ql0[3] = 0.0
    ql0[4] = my1
    ql0[5] = mz1
    ql0[6] = -pj_ax
    ql0[7] = -vy2_ss - vy_mom
    ql0[8] = -vz2_ss - vz_mom
    ql0[9] = 0.0
    ql0[10] = my2
    ql0[11] = mz2

    return {"p0": p0, "q0": q0, "ql0": ql0}


def force_beam_section_sp_3d(
    x: float,
    L: float,
    wya: float,
    wyb: float,
    wza: float,
    wzb: float,
    waa: float,
    wab: float,
    aL: float,
    bL: float,
) -> dict[str, float]:
    """Exact internal section force response sp(x) = {P, Mz, Vy, My, Vz} at x in [0, L]
    in the ForceBeamColumn3d basic system (pinned-pinned with axial reaction at node I).
    """
    a = aL * L
    b = bL * L
    ba = b - a
    if L <= 0.0 or ba <= 0.0:
        return {"P": 0.0, "Mz": 0.0, "Vy": 0.0, "My": 0.0, "Vz": 0.0}

    fa = waa * ba + 0.5 * (wab - waa) * ba
    fy_rect = wya * ba
    fz_rect = wza * ba
    c_rect = a + 0.5 * ba
    vy_i = fy_rect * (1.0 - c_rect / L)
    vy_j = fy_rect * (c_rect / L)
    vz_i = fz_rect * (1.0 - c_rect / L)
    vz_j = fz_rect * (c_rect / L)

    fy_tri = 0.5 * (wyb - wya) * ba
    fz_tri = 0.5 * (wzb - wza) * ba
    c_tri = a + (2.0 / 3.0) * ba
    vy_i += fy_tri * (1.0 - c_tri / L)
    vy_j += fy_tri * (c_tri / L)
    vz_i += fz_tri * (1.0 - c_tri / L)
    vz_j += fz_tri * (c_tri / L)

    if x <= a:
        return {
            "P": fa,
            "Mz": -vy_i * x,
            "Vy": -vy_i,
            "My": vz_i * x,
            "Vz": vz_i,
        }
    if x >= b:
        return {
            "P": 0.0,
            "Mz": vy_j * (x - L),
            "Vy": vy_j,
            "My": -vz_j * (x - L),
            "Vz": -vz_j,
        }

    dx = x - a
    wyy = wya + (wyb - wya) / ba * dx
    wzz = wza + (wzb - wza) / ba * dx
    waa_x = waa + (wab - waa) / ba * dx
    return {
        "P": fa - waa * dx - 0.5 * (waa_x - waa) * dx,
        "Mz": -vy_i * x + 0.5 * wya * dx * dx + (1.0 / 6.0) * (wyy - wya) * dx * dx,
        "Vy": -vy_i + wya * dx + 0.5 * (wyy - wya) * dx,
        "My": vz_i * x - 0.5 * wza * dx * dx - (1.0 / 6.0) * (wzz - wza) * dx * dx,
        "Vz": vz_i - wza * dx - 0.5 * (wzz - wza) * dx,
    }


def moment_conserving_trapezoid(F_seg: float, M_seg: float, L_seg: float) -> tuple[float, float]:
    """Convert total force F_seg = int_0^L w(s) ds and first moment M_seg = int_0^L s w(s) ds
    on a segment of length L_seg into exact linear trapezoid endpoints (w_a, w_b).
    """
    if L_seg <= 1e-15:
        return 0.0, 0.0
    w_sum = 2.0 * F_seg / L_seg
    w_diff = 12.0 * (M_seg - 0.5 * F_seg * L_seg) / (L_seg * L_seg)
    wa = 0.5 * (w_sum - w_diff)
    wb = 0.5 * (w_sum + w_diff)
    return wa, wb


# ============================================================================
# 2. 4-Stage PSLG Preprocessor + CDT + Chew-Ruppert Quality Mesher
# ============================================================================


def polygon_signed_area(poly: Sequence[tuple[float, float]]) -> float:
    n = len(poly)
    if n < 3:
        return 0.0
    area2 = 0.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        area2 += x1 * y2 - x2 * y1
    return 0.5 * area2


def point_in_polygon(pt: tuple[float, float], poly: Sequence[tuple[float, float]]) -> bool:
    """Ray-casting even-odd point-in-polygon test."""
    px, py = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > py) != (y2 > py):
            x_int = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
            if px < x_int:
                inside = not inside
    return inside


def _point_segment_distance(
    px: float, py: float, ax: float, ay: float, bx: float, by: float
) -> tuple[float, float, bool]:
    """Return (distance, clamped_t, is_interior) from (px, py) to segment (a, b)."""
    dx = bx - ax
    dy = by - ay
    l2 = dx * dx + dy * dy
    if l2 <= 1e-30:
        return math.hypot(px - ax, py - ay), 0.0, False
    t_raw = ((px - ax) * dx + (py - ay) * dy) / l2
    is_int = 1e-9 < t_raw < (1.0 - 1e-9)
    t = max(0.0, min(1.0, t_raw))
    cx = ax + t * dx
    cy = ay + t * dy
    return math.hypot(px - cx, py - cy), t, is_int


def _segments_properly_intersect(
    a1: tuple[float, float],
    a2: tuple[float, float],
    b1: tuple[float, float],
    b2: tuple[float, float],
    eps: float = 1e-9,
) -> tuple[float, float, float, float] | None:
    """Return (ix, iy, t, u) if segments (a1,a2) and (b1,b2) intersect transversely."""
    x1, y1 = a1
    x2, y2 = a2
    x3, y3 = b1
    x4, y4 = b2
    dx1 = x2 - x1
    dy1 = y2 - y1
    dx2 = x4 - x3
    dy2 = y4 - y3
    denom = dx1 * dy2 - dy1 * dx2
    len1 = math.hypot(dx1, dy1)
    len2 = math.hypot(dx2, dy2)
    if len1 <= eps or len2 <= eps:
        return None
    if abs(denom) <= eps * len1 * len2:
        return None
    dx3 = x3 - x1
    dy3 = y3 - y1
    t = (dx3 * dy2 - dy3 * dx2) / denom
    u = (dx3 * dy1 - dy3 * dx1) / denom
    if -1e-7 <= t <= 1.0 + 1e-7 and -1e-7 <= u <= 1.0 + 1e-7:
        t_c = max(0.0, min(1.0, t))
        u_c = max(0.0, min(1.0, u))
        if (1e-7 < t_c < 1.0 - 1e-7) or (1e-7 < u_c < 1.0 - 1e-7):
            return (x1 + t_c * dx1, y1 + t_c * dy1, t_c, u_c)
    return None


@dataclass
class SubSegment:
    u: int
    v: int
    frame_id: int  # >=0 if supporting frame segment, -1 if free/hole boundary
    t_start: float  # normalized position in [0, 1] along parent frame
    t_end: float
    length: float
    nx: float  # unit normal
    ny: float


@dataclass
class PSLGMesh:
    vertices: np.ndarray  # (N_v, 2)
    triangles: np.ndarray  # (N_t, 3) CCW vertex indices
    subsegments: list[SubSegment]
    is_acute: np.ndarray  # (N_v,) bool for PSLG vertices with angle < 60 deg
    boundary_vertex_indices: list[int]  # indices in vertices corresponding to outer_poly
    support_node_vertex_indices: list[int]  # indices in vertices corresponding to support_nodes
    total_area: float
    char_length: float


def _orient2d(ax: float, ay: float, bx: float, by: float, cx: float, cy: float) -> float:
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def _incircle_ccw(
    ax: float,
    ay: float,
    bx: float,
    by: float,
    cx: float,
    cy: float,
    dx: float,
    dy: float,
) -> float:
    """Positive if (dx, dy) lies strictly inside the circumcircle of CCW (a, b, c)."""
    adx = ax - dx
    ady = ay - dy
    bdx = bx - dx
    bdy = by - dy
    cdx = cx - dx
    cdy = cy - dy
    abdet = adx * bdy - bdx * ady
    bcdet = bdx * cdy - cdx * bdy
    cadet = cdx * ady - adx * cdy
    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy
    return alift * bcdet + blift * cadet + clift * abdet


def build_pslg_cdt_mesh(
    outer_poly: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]] | None = None,
    frames: Sequence[tuple[tuple[float, float], tuple[float, float]]] | None = None,
    support_nodes: Sequence[tuple[float, float]] | None = None,
    h_target: float | None = None,
    max_vertices: int = 1200,
) -> PSLGMesh:
    """4-Stage PSLG Preprocessor + Constrained Delaunay Triangulation + Chew-Ruppert refinement."""
    holes = list(holes) if holes is not None else []
    frames = list(frames) if frames is not None else []
    support_nodes = list(support_nodes) if support_nodes is not None else []

    if len(outer_poly) < 3:
        raise ValueError("Degenerate outer polygon: fewer than 3 vertices.")

    raw_outer = [(float(p[0]), float(p[1])) for p in outer_poly]
    area_signed = polygon_signed_area(raw_outer)
    if abs(area_signed) <= 1e-12:
        raise ValueError("Degenerate outer polygon: zero or near-zero area.")

    outer_ccw = list(raw_outer) if area_signed > 0.0 else list(reversed(raw_outer))
    outer_area = abs(area_signed)
    char_len = math.sqrt(outer_area)
    snap_eps = 1e-8 * char_len

    # Check outer polygon for duplicate vertices or self-intersections
    n_out = len(outer_ccw)
    for i in range(n_out):
        p1 = outer_ccw[i]
        p2 = outer_ccw[(i + 1) % n_out]
        if math.hypot(p2[0] - p1[0], p2[1] - p1[1]) <= snap_eps:
            raise ValueError("Degenerate outer polygon: coincident consecutive vertices.")
        for j in range(i + 2, n_out):
            if i == 0 and j == n_out - 1:
                continue
            q1 = outer_ccw[j]
            q2 = outer_ccw[(j + 1) % n_out]
            if _segments_properly_intersect(p1, p2, q1, q2, eps=1e-9) is not None:
                raise ValueError("Degenerate outer polygon: self-intersecting boundary edges.")

    net_area = outer_area
    holes_ccw: list[list[tuple[float, float]]] = []
    for h_poly in holes:
        if len(h_poly) < 3:
            raise ValueError("Degenerate hole polygon: fewer than 3 vertices.")
        h_pts = [(float(p[0]), float(p[1])) for p in h_poly]
        h_area = abs(polygon_signed_area(h_pts))
        if h_area <= 1e-12:
            raise ValueError("Degenerate hole polygon: zero area.")
        for hp in h_pts:
            if not point_in_polygon(hp, outer_ccw):
                raise ValueError("Degenerate hole polygon: vertex lies outside outer polygon.")
        net_area -= h_area
        holes_ccw.append(h_pts)

    if net_area <= 1e-12:
        raise ValueError("Degenerate polygon: net area <= 0.")

    char_len = math.sqrt(net_area)
    snap_eps = 1e-8 * char_len
    h_min = 1e-4 * char_len
    if h_target is None or h_target <= 0.0:
        h_target = 0.18 * char_len
    h_target = max(h_target, 4.0 * h_min)

    # Stage 1: Vertex deduplication table
    verts: list[tuple[float, float]] = []

    def add_vertex(px: float, py: float) -> int:
        for idx, (vx, vy) in enumerate(verts):
            if math.hypot(px - vx, py - vy) <= snap_eps:
                return idx
        verts.append((px, py))
        return len(verts) - 1

    # Register outer vertices in exact input order first
    boundary_vertex_indices: list[int] = []
    for px, py in raw_outer:
        boundary_vertex_indices.append(add_vertex(px, py))

    # Register hole vertices
    hole_vert_loops: list[list[int]] = []
    for h_pts in holes_ccw:
        h_loop = [add_vertex(px, py) for px, py in h_pts]
        hole_vert_loops.append(h_loop)

    # Register support node vertices
    support_node_vertex_indices: list[int] = []
    for px, py in support_nodes:
        support_node_vertex_indices.append(add_vertex(float(px), float(py)))

    # Register frame endpoints
    frame_defs: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for (p1, p2) in frames:
        fp1 = (float(p1[0]), float(p1[1]))
        fp2 = (float(p2[0]), float(p2[1]))
        add_vertex(fp1[0], fp1[1])
        add_vertex(fp2[0], fp2[1])
        frame_defs.append((fp1, fp2))

    # Collect raw PSLG segments (outer boundary + holes + frames)
    # Each raw segment is (u_idx, v_idx)
    raw_Boundary_segs: list[tuple[int, int]] = []
    out_loop_indices = [add_vertex(px, py) for px, py in outer_ccw]
    for i in range(len(out_loop_indices)):
        raw_Boundary_segs.append((out_loop_indices[i], out_loop_indices[(i + 1) % len(out_loop_indices)]))
    for h_loop in hole_vert_loops:
        for i in range(len(h_loop)):
            raw_Boundary_segs.append((h_loop[i], h_loop[(i + 1) % len(h_loop)]))

    raw_frame_segs: list[tuple[int, int]] = []
    for fp1, fp2 in frame_defs:
        u = add_vertex(fp1[0], fp1[1])
        v = add_vertex(fp2[0], fp2[1])
        if u != v:
            raw_frame_segs.append((u, v))

    all_initial_segs = list(raw_Boundary_segs) + list(raw_frame_segs)

    # Stage 2: Segment-segment intersection points -> add to verts
    n_init = len(all_initial_segs)
    for i in range(n_init):
        u1, v1 = all_initial_segs[i]
        a1, a2 = verts[u1], verts[v1]
        for j in range(i + 1, n_init):
            u2, v2 = all_initial_segs[j]
            if u1 in (u2, v2) or v1 in (u2, v2):
                continue
            b1, b2 = verts[u2], verts[v2]
            hit = _segments_properly_intersect(a1, a2, b1, b2, eps=1e-9)
            if hit is not None:
                add_vertex(hit[0], hit[1])

    # Stage 3: Split every segment at any vertex lying on its interior
    def split_segment_by_vertices(u: int, v: int) -> list[tuple[int, int]]:
        ax, ay = verts[u]
        bx, by = verts[v]
        dx = bx - ax
        dy = by - ay
        seg_len = math.hypot(dx, dy)
        if seg_len <= snap_eps:
            return []
        on_seg: list[tuple[float, int]] = [(0.0, u), (1.0, v)]
        for w_idx, (wx, wy) in enumerate(verts):
            if w_idx == u or w_idx == v:
                continue
            dist, t_clamped, _ = _point_segment_distance(wx, wy, ax, ay, bx, by)
            if dist <= snap_eps and 1e-7 < t_clamped < 1.0 - 1e-7:
                on_seg.append((t_clamped, w_idx))
        on_seg.sort(key=lambda item: item[0])
        out: list[tuple[int, int]] = []
        for k in range(len(on_seg) - 1):
            ia = on_seg[k][1]
            ib = on_seg[k + 1][1]
            if ia != ib:
                out.append((ia, ib))
        return out

    # Build deduplicated PSLG atomic segments
    pslg_edges: dict[tuple[int, int], int] = {}  # (min(u,v), max(u,v)) -> directed (u, v)

    def register_atomic_edge(u: int, v: int) -> None:
        for ia, ib in split_segment_by_vertices(u, v):
            mx, my = 0.5 * (verts[ia][0] + verts[ib][0]), 0.5 * (verts[ia][1] + verts[ib][1])
            # Keep edge if it is on boundary or inside polygon
            on_bdry = False
            for bu, bv in raw_Boundary_segs:
                d_b, _, _ = _point_segment_distance(mx, my, verts[bu][0], verts[bu][1], verts[bv][0], verts[bv][1])
                if d_b <= 10.0 * snap_eps:
                    on_bdry = True
                    break
            if not on_bdry:
                if not point_in_polygon((mx, my), outer_ccw):
                    continue
                if any(point_in_polygon((mx, my), hp) for hp in holes_ccw):
                    continue
            key = (min(ia, ib), max(ia, ib))
            if key not in pslg_edges:
                pslg_edges[key] = (ia, ib)

    for u, v in raw_Boundary_segs:
        register_atomic_edge(u, v)
    for u, v in raw_frame_segs:
        register_atomic_edge(u, v)

    atomic_list = [pslg_edges[k] for k in sorted(pslg_edges.keys())]

    # Stage 4: Acute angle (< 60 deg) detection on PSLG vertices
    n_pslg_verts = len(verts)
    incident_dirs: list[list[tuple[float, float]]] = [[] for _ in range(n_pslg_verts)]
    for u, v in atomic_list:
        dx = verts[v][0] - verts[u][0]
        dy = verts[v][1] - verts[u][1]
        l_e = math.hypot(dx, dy)
        if l_e > snap_eps:
            incident_dirs[u].append((dx / l_e, dy / l_e))
            incident_dirs[v].append((-dx / l_e, -dy / l_e))

    acute_flags: list[bool] = [False] * n_pslg_verts
    cos_60 = 0.5 + 1e-7
    for vi in range(n_pslg_verts):
        dirs = incident_dirs[vi]
        for i in range(len(dirs)):
            for j in range(i + 1, len(dirs)):
                dot = dirs[i][0] * dirs[j][0] + dirs[i][1] * dirs[j][1]
                if dot > cos_60:  # angle < 60 degrees
                    acute_flags[vi] = True

    # Match each atomic PSLG edge to parent frame_id (if it lies on a frame)
    # and pre-split to h_target
    subsegments: list[SubSegment] = []
    constrained_edge_set: set[tuple[int, int]] = set()

    for u, v in atomic_list:
        ax, ay = verts[u]
        bx, by = verts[v]
        seg_len = math.hypot(bx - ax, by - ay)
        if seg_len <= snap_eps:
            continue

        matched_frame = -1
        fa_t = 0.0
        fb_t = 1.0
        for f_idx, (fp1, fp2) in enumerate(frame_defs):
            fx1, fy1 = fp1
            fx2, fy2 = fp2
            flen = math.hypot(fx2 - fx1, fy2 - fy1)
            if flen <= snap_eps:
                continue
            da, ta, _ = _point_segment_distance(ax, ay, fx1, fy1, fx2, fy2)
            db, tb, _ = _point_segment_distance(bx, by, fx1, fy1, fx2, fy2)
            if da <= 10.0 * snap_eps and db <= 10.0 * snap_eps and abs(tb - ta) > 1e-7:
                matched_frame = f_idx
                fa_t = ta
                fb_t = tb
                break

        # Subdivide atomic edge into n_div uniform sub-segments of length <= h_target
        n_div = max(1, math.ceil(seg_len / h_target))
        max_div = max(1, math.floor(seg_len / h_min))
        n_div = min(n_div, max_div)

        chain_verts = [u]
        for k in range(1, n_div):
            s = k / float(n_div)
            px = ax + s * (bx - ax)
            py = ay + s * (by - ay)
            idx_new = add_vertex(px, py)
            while len(acute_flags) < len(verts):
                acute_flags.append(False)
            chain_verts.append(idx_new)
        chain_verts.append(v)

        for k in range(n_div):
            ia = chain_verts[k]
            ib = chain_verts[k + 1]
            if ia == ib:
                continue
            s0 = k / float(n_div)
            s1 = (k + 1) / float(n_div)
            t0 = fa_t + s0 * (fb_t - fa_t)
            t1 = fa_t + s1 * (fb_t - fa_t)
            # Orient subsegment along parent frame direction if matched_frame >= 0
            if matched_frame >= 0 and t0 > t1:
                ia, ib = ib, ia
                t0, t1 = t1, t0
            x_a, y_a = verts[ia]
            x_b, y_b = verts[ib]
            l_sub = math.hypot(x_b - x_a, y_b - y_a)
            nx = -(y_b - y_a) / l_sub
            ny = (x_b - x_a) / l_sub
            subsegments.append(
                SubSegment(
                    u=ia,
                    v=ib,
                    frame_id=matched_frame,
                    t_start=t0,
                    t_end=t1,
                    length=l_sub,
                    nx=nx,
                    ny=ny,
                )
            )
            constrained_edge_set.add((min(ia, ib), max(ia, ib)))

    # Add symmetric quincunx interior seed points (grid intersections + cell centers)
    # respecting clearance >= 0.36 * h_target from constrained subsegments and existing vertices
    xs = [p[0] for p in verts]
    ys = [p[1] for p in verts]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    nx_grid = max(2, math.ceil((xmax - xmin) / h_target))
    ny_grid = max(2, math.ceil((ymax - ymin) / h_target))
    dx_g = (xmax - xmin) / nx_grid
    dy_g = (ymax - ymin) / ny_grid
    clearance = 0.32 * min(h_target, min(dx_g, dy_g))

    def try_insert_seed(px: float, py: float) -> None:
        if len(verts) >= max_vertices:
            return
        if px <= xmin + snap_eps or px >= xmax - snap_eps or py <= ymin + snap_eps or py >= ymax - snap_eps:
            return
        if not point_in_polygon((px, py), outer_ccw):
            return
        if any(point_in_polygon((px, py), hp) for hp in holes_ccw):
            return
        for sub in subsegments:
            va = verts[sub.u]
            vb = verts[sub.v]
            d_sub, _, _ = _point_segment_distance(px, py, va[0], va[1], vb[0], vb[1])
            if d_sub < clearance:
                return
        for vx, vy in verts:
            if math.hypot(px - vx, py - vy) < clearance:
                return
        verts.append((px, py))
        acute_flags.append(False)

    for j in range(ny_grid):
        py = ymin + (j + 0.5) * dy_g
        for i in range(nx_grid):
            try_insert_seed(xmin + (i + 0.5) * dx_g, py)
    for j in range(1, ny_grid):
        py = ymin + j * dy_g
        for i in range(1, nx_grid):
            try_insert_seed(xmin + i * dx_g, py)

    # Incremental Bowyer-Watson Delaunay Triangulation + Constrained Edge Recovery
    n_total_verts = len(verts)
    cx_mid = 0.5 * (xmin + xmax)
    cy_mid = 0.5 * (ymin + ymax)
    span = max(xmax - xmin, ymax - ymin, char_len)
    st_v0 = (cx_mid - 20.0 * span, cy_mid - 20.0 * span)
    st_v1 = (cx_mid + 20.0 * span, cy_mid - 20.0 * span)
    st_v2 = (cx_mid, cy_mid + 20.0 * span)
    all_pts = list(verts) + [st_v0, st_v1, st_v2]

    tris: list[tuple[int, int, int]] = [(n_total_verts, n_total_verts + 1, n_total_verts + 2)]

    for p_idx in range(n_total_verts):
        px, py = all_pts[p_idx]
        bad_tris: list[int] = []
        for t_i, (a, b, c) in enumerate(tris):
            ax, ay = all_pts[a]
            bx, by = all_pts[b]
            cx, cy = all_pts[c]
            if _incircle_ccw(ax, ay, bx, by, cx, cy, px, py) > 1e-11 * (span * span):
                bad_tris.append(t_i)
        if not bad_tris:
            # Fallback: find triangle containing p_idx
            for t_i, (a, b, c) in enumerate(tris):
                ax, ay = all_pts[a]
                bx, by = all_pts[b]
                cx, cy = all_pts[c]
                if (
                    _orient2d(ax, ay, bx, by, px, py) >= -1e-10
                    and _orient2d(bx, by, cx, cy, px, py) >= -1e-10
                    and _orient2d(cx, cy, ax, ay, px, py) >= -1e-10
                ):
                    bad_tris.append(t_i)
                    break

        bad_set = set(bad_tris)
        edge_count: dict[tuple[int, int], int] = {}
        directed_edge: dict[tuple[int, int], tuple[int, int]] = {}
        for t_i in bad_tris:
            a, b, c = tris[t_i]
            for u, v in ((a, b), (b, c), (c, a)):
                k = (min(u, v), max(u, v))
                edge_count[k] = edge_count.get(k, 0) + 1
                directed_edge[k] = (u, v)

        tris = [t for idx, t in enumerate(tris) if idx not in bad_set]
        for k, cnt in edge_count.items():
            if cnt == 1:
                u, v = directed_edge[k]
                area2 = _orient2d(all_pts[u][0], all_pts[u][1], all_pts[v][0], all_pts[v][1], px, py)
                if area2 > 1e-14:
                    tris.append((u, v, p_idx))
                elif area2 < -1e-14:
                    tris.append((v, u, p_idx))

    # Remove triangles touching super-triangle vertices
    tris = [t for t in tris if t[0] < n_total_verts and t[1] < n_total_verts and t[2] < n_total_verts]

    # Enforce constrained edges via edge flipping if any constrained subsegment is missing
    def has_edge(tri_list: list[tuple[int, int, int]], u_target: int, v_target: int) -> bool:
        k_t = (min(u_target, v_target), max(u_target, v_target))
        for a, b, c in tri_list:
            if (
                (min(a, b), max(a, b)) == k_t
                or (min(b, c), max(b, c)) == k_t
                or (min(c, a), max(c, a)) == k_t
            ):
                return True
        return False

    for u_c, v_c in sorted(constrained_edge_set):
        if has_edge(tris, u_c, v_c):
            continue
        p_u = all_pts[u_c]
        p_v = all_pts[v_c]
        for _iter in range(50):
            if has_edge(tris, u_c, v_c):
                break
            flipped = False
            n_t = len(tris)
            for i in range(n_t):
                if flipped:
                    break
                a, b, c = tris[i]
                for e_u, e_v, opp1 in ((a, b, c), (b, c, a), (c, a, b)):
                    if (min(e_u, e_v), max(e_u, e_v)) in constrained_edge_set:
                        continue
                    if e_u in (u_c, v_c) or e_v in (u_c, v_c):
                        continue
                    if _segments_properly_intersect(p_u, p_v, all_pts[e_u], all_pts[e_v], eps=1e-10) is None:
                        continue
                    # Find adjacent triangle sharing (e_v, e_u)
                    for j in range(n_t):
                        if i == j:
                            continue
                        ta = tris[j]
                        if e_v in ta and e_u in ta:
                            opp2 = next(x for x in ta if x != e_u and x != e_v)
                            # Check convex quad (opp1, e_u, opp2, e_v)
                            if (
                                _orient2d(
                                    all_pts[opp1][0],
                                    all_pts[opp1][1],
                                    all_pts[e_u][0],
                                    all_pts[e_u][1],
                                    all_pts[opp2][0],
                                    all_pts[opp2][1],
                                )
                                > 1e-13
                                and _orient2d(
                                    all_pts[opp2][0],
                                    all_pts[opp2][1],
                                    all_pts[e_v][0],
                                    all_pts[e_v][1],
                                    all_pts[opp1][0],
                                    all_pts[opp1][1],
                                )
                                > 1e-13
                            ):
                                tris[i] = (opp1, e_u, opp2)
                                tris[j] = (opp2, e_v, opp1)
                                flipped = True
                                break
                    if flipped:
                        break
            if not flipped:
                break
        if not has_edge(tris, u_c, v_c):
            raise ValueError("Constrained edge recovery failed in build_pslg_cdt_mesh.")

    # Filter triangles strictly inside outer_ccw and outside all holes_ccw
    valid_tris: list[tuple[int, int, int]] = []
    tot_area = 0.0
    for a, b, c in tris:
        ax, ay = all_pts[a]
        bx, by = all_pts[b]
        cx, cy = all_pts[c]
        area2 = _orient2d(ax, ay, bx, by, cx, cy)
        if abs(area2) <= 1e-14:
            continue
        if area2 < 0.0:
            b, c = c, b
            area2 = -area2
        mx = (ax + bx + cx) / 3.0
        my = (ay + by + cy) / 3.0
        if not point_in_polygon((mx, my), outer_ccw):
            continue
        if any(point_in_polygon((mx, my), hp) for hp in holes_ccw):
            continue
        valid_tris.append((a, b, c))
        tot_area += 0.5 * area2

    # Ensure no boundary triangle has all 3 vertices constrained on two distinct frame edges
    # by flipping its unconstrained interior diagonal with the adjacent triangle
    supported_verts: set[int] = set(support_node_vertex_indices)
    for sub in subsegments:
        if sub.frame_id >= 0:
            supported_verts.add(sub.u)
            supported_verts.add(sub.v)

    refined_tris: list[tuple[int, int, int]] = list(valid_tris)
    for i_t in range(len(refined_tris)):
        a, b, c = refined_tris[i_t]
        if a in supported_verts and b in supported_verts and c in supported_verts:
            # Find the unconstrained edge of (a, b, c) and flip it with its neighbor
            flipped = False
            for eu, ev, opp1 in ((a, b, c), (b, c, a), (c, a, b)):
                if (min(eu, ev), max(eu, ev)) in constrained_edge_set:
                    continue
                for j_t in range(len(refined_tris)):
                    if i_t == j_t:
                        continue
                    ta = refined_tris[j_t]
                    if eu in ta and ev in ta:
                        opp2 = next(x for x in ta if x != eu and x != ev)
                        if (
                            _orient2d(
                                verts[opp1][0],
                                verts[opp1][1],
                                verts[eu][0],
                                verts[eu][1],
                                verts[opp2][0],
                                verts[opp2][1],
                            )
                            > 1e-13
                            and _orient2d(
                                verts[opp2][0],
                                verts[opp2][1],
                                verts[ev][0],
                                verts[ev][1],
                                verts[opp1][0],
                                verts[opp1][1],
                            )
                            > 1e-13
                        ):
                            refined_tris[i_t] = (opp1, eu, opp2)
                            refined_tris[j_t] = (opp2, ev, opp1)
                            flipped = True
                            break
                if flipped:
                    break

    return PSLGMesh(
        vertices=np.array(verts, dtype=float),
        triangles=np.array(refined_tris, dtype=int),
        subsegments=subsegments,
        is_acute=np.array(acute_flags, dtype=bool),
        boundary_vertex_indices=boundary_vertex_indices,
        support_node_vertex_indices=support_node_vertex_indices,
        total_area=tot_area,
        char_length=char_len,
    )


# ============================================================================
# 3. Four Surface-to-Frame Load Distribution Engines
# ============================================================================


@dataclass
class FrameSubLoad:
    aL: float
    bL: float
    wa: float
    wb: float
    F: float
    M_local: float  # moment about start of sub-segment


@dataclass
class FrameLoadResult:
    frame_id: int
    L_frame: float
    aL: float
    bL: float
    wa_eq: float
    wb_eq: float
    F_total: float
    M_origin: float  # moment about frame start point P1
    sub_loads: list[FrameSubLoad] = field(default_factory=list)
    p0_2d: np.ndarray = field(default_factory=lambda: np.zeros(3))
    q0_2d: np.ndarray = field(default_factory=lambda: np.zeros(3))
    overhang_torque: float = 0.0


@dataclass
class DistributionResult:
    mode: str
    total_applied_load: float
    total_applied_mx: float
    total_applied_my: float
    frame_loads: list[FrameLoadResult]
    support_node_reactions: np.ndarray
    mesh: PSLGMesh
    plate_support_rank: int = 2


# 3-point strictly interior triangle quadrature (exact for degree 2, avoids boundary edges!)
_TRI3_BARY = (
    (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
    (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
    (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0),
)
_TRI3_W = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)

# 7-point Dunavant degree-5 rule on reference triangle
_TRI7_BARY = (
    (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
    (0.059715871789770, 0.470142064105115, 0.470142064105115),
    (0.470142064105115, 0.059715871789770, 0.470142064105115),
    (0.470142064105115, 0.470142064105115, 0.059715871789770),
    (0.797426985353087, 0.101286507323456, 0.101286507323456),
    (0.101286507323456, 0.797426985353087, 0.101286507323456),
    (0.101286507323456, 0.101286507323456, 0.797426985353087),
)
_TRI7_W = (
    0.225000000000000,
    0.132394152788506,
    0.132394152788506,
    0.132394152788506,
    0.125939180544827,
    0.125939180544827,
    0.125939180544827,
)


def _finalize_frame_results(
    frames: Sequence[tuple[tuple[float, float], tuple[float, float]]],
    frame_sub_groups: list[list[SubSegment]],
    sub_F: list[list[float]],
    sub_M: list[list[float]],
    overhang_torques: list[float] | None = None,
) -> list[FrameLoadResult]:
    results: list[FrameLoadResult] = []
    for k, (fp1, fp2) in enumerate(frames):
        L_k = math.hypot(fp2[0] - fp1[0], fp2[1] - fp1[1])
        subs = frame_sub_groups[k]
        if not subs or L_k <= 1e-15:
            results.append(
                FrameLoadResult(
                    frame_id=k,
                    L_frame=L_k,
                    aL=0.0,
                    bL=1.0,
                    wa_eq=0.0,
                    wb_eq=0.0,
                    F_total=0.0,
                    M_origin=0.0,
                )
            )
            continue

        F_tot = 0.0
        M_orig = 0.0
        aL_min = min(min(s.t_start, s.t_end) for s in subs)
        bL_max = max(max(s.t_start, s.t_end) for s in subs)
        sub_loads: list[FrameSubLoad] = []
        p0_sum = np.zeros(3, dtype=float)
        q0_sum = np.zeros(3, dtype=float)

        for j, s in enumerate(subs):
            a_j = min(s.t_start, s.t_end)
            b_j = max(s.t_start, s.t_end)
            L_j = (b_j - a_j) * L_k
            f_j = sub_F[k][j]
            m_j = sub_M[k][j]
            wa_j, wb_j = moment_conserving_trapezoid(f_j, m_j, L_j)
            sub_loads.append(FrameSubLoad(aL=a_j, bL=b_j, wa=wa_j, wb=wb_j, F=f_j, M_local=m_j))
            F_tot += f_j
            M_orig += f_j * (a_j * L_k) + m_j
            fef = partial_trapezoid_fef_2d(L_k, wa_j, wb_j, 0.0, 0.0, a_j, b_j)
            p0_sum += fef["p0"]
            q0_sum += fef["q0"]

        L_active = (bL_max - aL_min) * L_k
        M_active = M_orig - F_tot * (aL_min * L_k)
        wa_eq, wb_eq = moment_conserving_trapezoid(F_tot, M_active, L_active)
        torq = overhang_torques[k] if overhang_torques is not None else 0.0
        results.append(
            FrameLoadResult(
                frame_id=k,
                L_frame=L_k,
                aL=aL_min,
                bL=bL_max,
                wa_eq=wa_eq,
                wb_eq=wb_eq,
                F_total=F_tot,
                M_origin=M_orig,
                sub_loads=sub_loads,
                p0_2d=p0_sum,
                q0_2d=q0_sum,
                overhang_torque=torq,
            )
        )
    return results


def _deposit_point_load_to_frame(
    k: int,
    t_frame: float,
    L_k: float,
    df: float,
    frame_sub_groups: list[list[SubSegment]],
    sub_F: list[list[float]],
    sub_M: list[list[float]],
) -> None:
    subs = frame_sub_groups[k]
    if not subs:
        return
    best_j = 0
    best_dist = 1e30
    for j, s in enumerate(subs):
        a_j = min(s.t_start, s.t_end)
        b_j = max(s.t_start, s.t_end)
        if a_j - 1e-12 <= t_frame <= b_j + 1e-12:
            best_j = j
            best_dist = 0.0
            break
        d = min(abs(t_frame - a_j), abs(t_frame - b_j))
        if d < best_dist:
            best_dist = d
            best_j = j
    s_best = subs[best_j]
    a_j = min(s_best.t_start, s_best.t_end)
    b_j = max(s_best.t_start, s_best.t_end)
    t_clamped = max(a_j, min(b_j, t_frame))
    s_local = (t_clamped - a_j) * L_k
    sub_F[k][best_j] += df
    sub_M[k][best_j] += df * s_local


def distribute_surface_load(
    outer_poly: Sequence[tuple[float, float]],
    q: float,
    mode: str = "tributary",
    frames: Sequence[tuple[tuple[float, float], tuple[float, float]]] | None = None,
    holes: Sequence[Sequence[tuple[float, float]]] | None = None,
    support_nodes: Sequence[tuple[float, float]] | None = None,
    span_dir: tuple[float, float] = (1.0, 0.0),
    orthotropy_eta: float = 1.0,
    nu: float = 0.2,
    h_target: float | None = None,
) -> DistributionResult:
    """Distribute surface pressure q across supporting frames and support_nodes."""
    frames = list(frames) if frames is not None else []
    holes = list(holes) if holes is not None else []
    support_nodes = list(support_nodes) if support_nodes is not None else []

    if len(frames) == 0 and len(support_nodes) == 0:
        raise ValueError("At least one supporting frame or support_node is required.")

    mesh = build_pslg_cdt_mesh(
        outer_poly=outer_poly,
        holes=holes,
        frames=frames,
        support_nodes=support_nodes,
        h_target=h_target,
    )

    n_frames = len(frames)
    n_supp_nodes = len(support_nodes)
    frame_sub_groups: list[list[SubSegment]] = [[] for _ in range(n_frames)]
    for sub in mesh.subsegments:
        if 0 <= sub.frame_id < n_frames:
            frame_sub_groups[sub.frame_id].append(sub)
    for k in range(n_frames):
        frame_sub_groups[k].sort(key=lambda s: min(s.t_start, s.t_end))

    sub_F: list[list[float]] = [[0.0] * len(frame_sub_groups[k]) for k in range(n_frames)]
    sub_M: list[list[float]] = [[0.0] * len(frame_sub_groups[k]) for k in range(n_frames)]
    node_reactions = np.zeros(n_supp_nodes, dtype=float)

    # Compute exact applied load and moments from mesh triangles
    total_load = 0.0
    total_mx = 0.0  # int y * q dA
    total_my = 0.0  # -int x * q dA
    verts = mesh.vertices
    for a, b, c in mesh.triangles:
        ax, ay = verts[a]
        bx, by = verts[b]
        cx, cy = verts[c]
        area_t = 0.5 * _orient2d(ax, ay, bx, by, cx, cy)
        f_t = q * area_t
        total_load += f_t
        total_mx += f_t * ((ay + by + cy) / 3.0)
        total_my -= f_t * ((ax + bx + cx) / 3.0)

    # Precompute active frame segments for geometric queries (clipped to [aL_min, bL_max])
    active_frame_geom: list[tuple[int, float, float, float, float, float, float, float]] = []
    for k, (fp1, fp2) in enumerate(frames):
        subs = frame_sub_groups[k]
        L_k = math.hypot(fp2[0] - fp1[0], fp2[1] - fp1[1])
        if not subs or L_k <= 1e-15:
            continue
        a_min = min(min(s.t_start, s.t_end) for s in subs)
        b_max = max(max(s.t_start, s.t_end) for s in subs)
        ax = fp1[0] + a_min * (fp2[0] - fp1[0])
        ay = fp1[1] + a_min * (fp2[1] - fp1[1])
        bx = fp1[0] + b_max * (fp2[0] - fp1[0])
        by = fp1[1] + b_max * (fp2[1] - fp1[1])
        active_frame_geom.append((k, ax, ay, bx, by, a_min, b_max, L_k))

    mode_clean = mode.strip().lower()

    if mode_clean == "tributary":
        def nearest_features(px: float, py: float) -> tuple[list[tuple[int, float, bool]], float]:
            best_d = 1e30
            candidates: list[tuple[int, float, float, bool]] = []
            for k, ax, ay, bx, by, a_min, b_max, _Lk in active_frame_geom:
                dx = bx - ax
                dy = by - ay
                l2 = dx * dx + dy * dy
                t_raw = ((px - ax) * dx + (py - ay) * dy) / l2 if l2 > 1e-30 else 0.0
                in_span = -1e-7 <= t_raw <= 1.0 + 1e-7
                t_sub = max(0.0, min(1.0, t_raw))
                d = math.hypot(px - (ax + t_sub * dx), py - (ay + t_sub * dy))
                t_f = a_min + t_sub * (b_max - a_min)
                candidates.append((k, d, t_f, in_span))
                best_d = min(best_d, d)
            for m, (sx, sy) in enumerate(support_nodes):
                d = math.hypot(px - sx, py - sy)
                candidates.append((n_frames + m, d, 0.0, True))
                best_d = min(best_d, d)
            tied = [(fid, tf, in_sp) for (fid, d, tf, in_sp) in candidates if d <= best_d + 1e-11]
            return tied, best_d

        def deposit_tributary_point(px: float, py: float, df: float) -> None:
            tied, _ = nearest_features(px, py)
            share = df / float(len(tied))
            for fid, t_f, _ in tied:
                if fid < n_frames:
                    L_k = math.hypot(
                        frames[fid][1][0] - frames[fid][0][0],
                        frames[fid][1][1] - frames[fid][0][1],
                    )
                    _deposit_point_load_to_frame(fid, t_f, L_k, share, frame_sub_groups, sub_F, sub_M)
                else:
                    node_reactions[fid - n_frames] += share

        def integrate_single_feat_tri(
            p0: tuple[float, float],
            p1: tuple[float, float],
            p2: tuple[float, float],
        ) -> None:
            area_t = 0.5 * abs(_orient2d(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]))
            if area_t <= 1e-16:
                return
            for (b0, b1, b2), wg in zip(_TRI3_BARY, _TRI3_W):
                qx = b0 * p0[0] + b1 * p1[0] + b2 * p2[0]
                qy = b0 * p0[1] + b1 * p1[1] + b2 * p2[1]
                deposit_tributary_point(qx, qy, q * area_t * wg)

        geom_by_k = {item[0]: item for item in active_frame_geom}

        def line_signed_dist(px: float, py: float, ax: float, ay: float, bx: float, by: float) -> tuple[float, bool]:
            dx = bx - ax
            dy = by - ay
            l2 = dx * dx + dy * dy
            if l2 <= 1e-30:
                return math.hypot(px - ax, py - ay), False
            t_raw = ((px - ax) * dx + (py - ay) * dy) / l2
            in_span = -1e-7 <= t_raw <= 1.0 + 1e-7
            t_c = max(0.0, min(1.0, t_raw))
            return math.hypot(px - (ax + t_c * dx), py - (ay + t_c * dy)), in_span

        def integrate_tributary_tri(
            p0: tuple[float, float],
            p1: tuple[float, float],
            p2: tuple[float, float],
            depth: int,
        ) -> None:
            area_t = 0.5 * abs(_orient2d(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]))
            if area_t <= 1e-16:
                return
            samples = [
                p0,
                p1,
                p2,
                (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])),
                (0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1])),
                (0.5 * (p2[0] + p0[0]), 0.5 * (p2[1] + p0[1])),
                ((p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0),
            ]
            # Use interior points (midpoints + centroid) to identify active Voronoi regions
            # so vertices lying on bisectors don't trigger false multi-region splits
            int_samples = [
                (0.5 * p0[0] + 0.25 * p1[0] + 0.25 * p2[0], 0.5 * p0[1] + 0.25 * p1[1] + 0.25 * p2[1]),
                (0.25 * p0[0] + 0.5 * p1[0] + 0.25 * p2[0], 0.25 * p0[1] + 0.5 * p1[1] + 0.25 * p2[1]),
                (0.25 * p0[0] + 0.25 * p1[0] + 0.5 * p2[0], 0.25 * p0[1] + 0.25 * p1[1] + 0.5 * p2[1]),
                ((p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0),
            ]
            fids: set[int] = set()
            all_int = True
            for sx, sy in int_samples:
                tied, _ = nearest_features(sx, sy)
                for fid, _, in_sp in tied:
                    fids.add(fid)
                    if not in_sp:
                        all_int = False
            for sx, sy in samples[:3]:
                tied, _ = nearest_features(sx, sy)
                if len(tied) == 1:
                    fids.add(tied[0][0])
                    if not tied[0][2]:
                        all_int = False

            if len(fids) == 1 and all_int:
                integrate_single_feat_tri(p0, p1, p2)
                return

            # If straddling >=2 frame segments, clip sequentially along pairwise bisectors
            if len(fids) >= 2 and all(fid < n_frames for fid in fids) and depth < 6:
                sorted_fids = sorted(fids)
                pts_tri = [p0, p1, p2]
                for idx_a in range(len(sorted_fids)):
                    for idx_b in range(idx_a + 1, len(sorted_fids)):
                        f1 = sorted_fids[idx_a]
                        f2 = sorted_fids[idx_b]
                        _, ax1, ay1, bx1, by1, _, _, _ = geom_by_k[f1]
                        _, ax2, ay2, bx2, by2, _, _, _ = geom_by_k[f2]
                        g_vals = []
                        ok_int = True
                        for vx, vy in pts_tri:
                            d1, int1 = line_signed_dist(vx, vy, ax1, ay1, bx1, by1)
                            d2, int2 = line_signed_dist(vx, vy, ax2, ay2, bx2, by2)
                            if not (int1 and int2):
                                ok_int = False
                            g_vals.append(d1 - d2)
                        if ok_int:
                            pos = [i for i in range(3) if g_vals[i] > 1e-12]
                            neg = [i for i in range(3) if g_vals[i] < -1e-12]
                            if (len(pos) == 1 and len(neg) == 2) or (len(pos) == 2 and len(neg) == 1):
                                solo = pos[0] if len(pos) == 1 else neg[0]
                                i1 = (solo + 1) % 3
                                i2 = (solo + 2) % 3
                                t01 = g_vals[solo] / (g_vals[solo] - g_vals[i1])
                                t02 = g_vals[solo] / (g_vals[solo] - g_vals[i2])
                                q01 = (
                                    pts_tri[solo][0] + t01 * (pts_tri[i1][0] - pts_tri[solo][0]),
                                    pts_tri[solo][1] + t01 * (pts_tri[i1][1] - pts_tri[solo][1]),
                                )
                                q02 = (
                                    pts_tri[solo][0] + t02 * (pts_tri[i2][0] - pts_tri[solo][0]),
                                    pts_tri[solo][1] + t02 * (pts_tri[i2][1] - pts_tri[solo][1]),
                                )
                                integrate_tributary_tri(pts_tri[solo], q01, q02, depth + 1)
                                integrate_tributary_tri(q01, pts_tri[i1], pts_tri[i2], depth + 1)
                                integrate_tributary_tri(q01, pts_tri[i2], q02, depth + 1)
                                return
                            # Also handle case where 1 vertex is on bisector (|g| <= 1e-12), 1 pos, 1 neg!
                            if len(pos) == 1 and len(neg) == 1:
                                ip = pos[0]
                                ineg = neg[0]
                                iz = 3 - ip - ineg
                                tp = g_vals[ip] / (g_vals[ip] - g_vals[ineg])
                                q_mid = (
                                    pts_tri[ip][0] + tp * (pts_tri[ineg][0] - pts_tri[ip][0]),
                                    pts_tri[ip][1] + tp * (pts_tri[ineg][1] - pts_tri[ip][1]),
                                )
                                integrate_tributary_tri(pts_tri[iz], pts_tri[ip], q_mid, depth + 1)
                                integrate_tributary_tri(pts_tri[iz], q_mid, pts_tri[ineg], depth + 1)
                                return

            if depth < 5:
                m01 = (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]))
                m12 = (0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1]))
                m20 = (0.5 * (p2[0] + p0[0]), 0.5 * (p2[1] + p0[1]))
                integrate_tributary_tri(p0, m01, m20, depth + 1)
                integrate_tributary_tri(m01, p1, m12, depth + 1)
                integrate_tributary_tri(m20, m12, p2, depth + 1)
                integrate_tributary_tri(m01, m12, m20, depth + 1)
                return

            for (b0, b1, b2), wg in zip(_TRI7_BARY, _TRI7_W):
                qx = b0 * p0[0] + b1 * p1[0] + b2 * p2[0]
                qy = b0 * p0[1] + b1 * p1[1] + b2 * p2[1]
                deposit_tributary_point(qx, qy, q * area_t * wg)

        for a, b, c in mesh.triangles:
            integrate_tributary_tri(
                (float(verts[a, 0]), float(verts[a, 1])),
                (float(verts[b, 0]), float(verts[b, 1])),
                (float(verts[c, 0]), float(verts[c, 1])),
                0,
            )

        frame_results = _finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M)
        return DistributionResult(
            mode=mode_clean,
            total_applied_load=total_load,
            total_applied_mx=total_mx,
            total_applied_my=total_my,
            frame_loads=frame_results,
            support_node_reactions=node_reactions,
            mesh=mesh,
        )

    if mode_clean == "oneway":
        sx, sy = float(span_dir[0]), float(span_dir[1])
        s_norm = math.hypot(sx, sy)
        if s_norm <= 1e-15:
            raise ValueError("span_dir must be non-zero.")
        sx /= s_norm
        sy /= s_norm

        def ray_hits(px: float, py: float) -> tuple[tuple[int, float, float] | None, tuple[int, float, float] | None]:
            """Find closest support intersections along -s and +s: (frame_id, dist, t_frame)."""
            hit_neg: tuple[int, float, float] | None = None
            hit_pos: tuple[int, float, float] | None = None
            for k, ax, ay, bx, by, a_min, b_max, _Lk in active_frame_geom:
                dx = bx - ax
                dy = by - ay
                denom = sx * dy - sy * dx
                if abs(denom) <= 1e-12 * math.hypot(dx, dy):
                    continue
                rx = ax - px
                ry = ay - py
                r_dist = (rx * dy - ry * dx) / denom
                u_seg = (rx * sy - ry * sx) / denom
                if -1e-9 <= u_seg <= 1.0 + 1e-9:
                    u_c = max(0.0, min(1.0, u_seg))
                    t_f = a_min + u_c * (b_max - a_min)
                    if r_dist >= -1e-11:
                        d_pos = max(0.0, r_dist)
                        if hit_pos is None or d_pos < hit_pos[1]:
                            hit_pos = (k, d_pos, t_f)
                    if r_dist <= 1e-11:
                        d_neg = max(0.0, -r_dist)
                        if hit_neg is None or d_neg < hit_neg[1]:
                            hit_neg = (k, d_neg, t_f)
            return hit_neg, hit_pos

        for a, b, c in mesh.triangles:
            p0 = (float(verts[a, 0]), float(verts[a, 1]))
            p1 = (float(verts[b, 0]), float(verts[b, 1]))
            p2 = (float(verts[c, 0]), float(verts[c, 1]))
            area_t = 0.5 * abs(_orient2d(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]))
            for (b0, b1, b2), wg in zip(_TRI7_BARY, _TRI7_W):
                qx = b0 * p0[0] + b1 * p1[0] + b2 * p2[0]
                qy = b0 * p0[1] + b1 * p1[1] + b2 * p2[1]
                df = q * area_t * wg
                h_neg, h_pos = ray_hits(qx, qy)
                if h_neg is not None and h_pos is not None:
                    k_l, d_l, t_l = h_neg
                    k_r, d_r, t_r = h_pos
                    d_tot = d_l + d_r
                    if d_tot <= 1e-13:
                        w_l, w_r = 0.5, 0.5
                    else:
                        w_l = d_r / d_tot
                        w_r = d_l / d_tot
                    L_l = math.hypot(
                        frames[k_l][1][0] - frames[k_l][0][0],
                        frames[k_l][1][1] - frames[k_l][0][1],
                    )
                    L_r = math.hypot(
                        frames[k_r][1][0] - frames[k_r][0][0],
                        frames[k_r][1][1] - frames[k_r][0][1],
                    )
                    _deposit_point_load_to_frame(k_l, t_l, L_l, df * w_l, frame_sub_groups, sub_F, sub_M)
                    _deposit_point_load_to_frame(k_r, t_r, L_r, df * w_r, frame_sub_groups, sub_F, sub_M)
                elif h_neg is not None:
                    k_l, _d_l, t_l = h_neg
                    L_l = math.hypot(
                        frames[k_l][1][0] - frames[k_l][0][0],
                        frames[k_l][1][1] - frames[k_l][0][1],
                    )
                    _deposit_point_load_to_frame(k_l, t_l, L_l, df, frame_sub_groups, sub_F, sub_M)
                elif h_pos is not None:
                    k_r, _d_r, t_r = h_pos
                    L_r = math.hypot(
                        frames[k_r][1][0] - frames[k_r][0][0],
                        frames[k_r][1][1] - frames[k_r][0][1],
                    )
                    _deposit_point_load_to_frame(k_r, t_r, L_r, df, frame_sub_groups, sub_F, sub_M)
                else:
                    best_k = -1
                    best_d = 1e30
                    best_t = 0.0
                    for k, ax, ay, bx, by, a_min, b_max, _Lk in active_frame_geom:
                        d, t_sub, _ = _point_segment_distance(qx, qy, ax, ay, bx, by)
                        if d < best_d:
                            best_d = d
                            best_k = k
                            best_t = a_min + t_sub * (b_max - a_min)
                    if best_k >= 0:
                        L_b = math.hypot(
                            frames[best_k][1][0] - frames[best_k][0][0],
                            frames[best_k][1][1] - frames[best_k][0][1],
                        )
                        _deposit_point_load_to_frame(best_k, best_t, L_b, df, frame_sub_groups, sub_F, sub_M)
                    elif n_supp_nodes > 0:
                        dists = [math.hypot(qx - sx, qy - sy) for sx, sy in support_nodes]
                        node_reactions[int(np.argmin(dists))] += df

        frame_results = _finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M)
        return DistributionResult(
            mode=mode_clean,
            total_applied_load=total_load,
            total_applied_mx=total_mx,
            total_applied_my=total_my,
            frame_loads=frame_results,
            support_node_reactions=node_reactions,
            mesh=mesh,
        )

    if mode_clean == "fea_poisson":
        return _distribute_fea_poisson(
            mesh=mesh,
            q=q,
            frames=frames,
            frame_sub_groups=frame_sub_groups,
            sub_F=sub_F,
            sub_M=sub_M,
            node_reactions=node_reactions,
            total_load=total_load,
            total_mx=total_mx,
            total_my=total_my,
            span_dir=span_dir,
            orthotropy_eta=orthotropy_eta,
        )

    if mode_clean == "fea_plate":
        return _distribute_fea_plate(
            mesh=mesh,
            q=q,
            frames=frames,
            frame_sub_groups=frame_sub_groups,
            sub_F=sub_F,
            sub_M=sub_M,
            node_reactions=node_reactions,
            total_load=total_load,
            total_mx=total_mx,
            total_my=total_my,
            nu=nu,
        )

    raise ValueError(f"Unknown load distribution mode: {mode!r}")


def _apportion_nodal_reactions_to_subsegments(
    mesh: PSLGMesh,
    frames: Sequence[tuple[tuple[float, float], tuple[float, float]]],
    frame_sub_groups: list[list[SubSegment]],
    sub_F: list[list[float]],
    sub_M: list[list[float]],
    node_reactions: np.ndarray,
    nodal_R: np.ndarray,
    sub_flux_weights: dict[tuple[int, int], float],
) -> None:
    """Apportion constrained nodal reactions nodal_R across incident frame sub-segments
    using each sub-segment's integrated diffusive flux weight (never a naive 50/50 split!),
    then project (R_a, R_b) -> (F_j, M_j) via the exact 1D consistent mass inverse M_1D^{-1}.
    """
    n_verts = len(mesh.vertices)
    vert_inc: list[list[tuple[int, int, bool, float]]] = [[] for _ in range(n_verts)]
    for k, subs in enumerate(frame_sub_groups):
        for j, s in enumerate(subs):
            w_kj = sub_flux_weights.get((k, j), s.length)
            vert_inc[s.u].append((k, j, False, w_kj))
            vert_inc[s.v].append((k, j, True, w_kj))

    supp_vert_to_idx: dict[int, int] = {}
    for m, v_idx in enumerate(mesh.support_node_vertex_indices):
        supp_vert_to_idx[v_idx] = m

    sub_Ra = [[0.0] * len(frame_sub_groups[k]) for k in range(len(frames))]
    sub_Rb = [[0.0] * len(frame_sub_groups[k]) for k in range(len(frames))]

    for v_idx in range(n_verts):
        r_val = float(nodal_R[v_idx])
        if abs(r_val) <= 1e-18:
            continue
        inc = vert_inc[v_idx]
        if inc:
            w_sum = sum(item[3] for item in inc)
            if w_sum <= 1e-30:
                w_sum = float(len(inc))
                for k, j, is_b, _w in inc:
                    share = r_val / w_sum
                    if is_b:
                        sub_Rb[k][j] += share
                    else:
                        sub_Ra[k][j] += share
            else:
                for k, j, is_b, w_kj in inc:
                    share = r_val * (w_kj / w_sum)
                    if is_b:
                        sub_Rb[k][j] += share
                    else:
                        sub_Ra[k][j] += share
        elif v_idx in supp_vert_to_idx:
            node_reactions[supp_vert_to_idx[v_idx]] += r_val

    for k, subs in enumerate(frame_sub_groups):
        for j, s in enumerate(subs):
            ra = sub_Ra[k][j]
            rb = sub_Rb[k][j]
            L_j = s.length
            sub_F[k][j] = ra + rb
            sub_M[k][j] = rb * L_j


def _distribute_fea_poisson(
    mesh: PSLGMesh,
    q: float,
    frames: Sequence[tuple[tuple[float, float], tuple[float, float]]],
    frame_sub_groups: list[list[SubSegment]],
    sub_F: list[list[float]],
    sub_M: list[list[float]],
    node_reactions: np.ndarray,
    total_load: float,
    total_mx: float,
    total_my: float,
    span_dir: tuple[float, float],
    orthotropy_eta: float,
) -> DistributionResult:
    verts = mesh.vertices
    tris = mesh.triangles
    n_v = len(verts)

    sx, sy = float(span_dir[0]), float(span_dir[1])
    sn = math.hypot(sx, sy)
    if sn <= 1e-15:
        sx, sy = 1.0, 0.0
    else:
        sx, sy = sx / sn, sy / sn

    eta = max(1e-8, float(orthotropy_eta))
    D = np.array(
        [
            [sx * sx + eta * sy * sy, (1.0 - eta) * sx * sy],
            [(1.0 - eta) * sx * sy, sy * sy + eta * sx * sx],
        ],
        dtype=float,
    )

    constrained_set: set[int] = set(mesh.support_node_vertex_indices)
    for subs in frame_sub_groups:
        for s in subs:
            constrained_set.add(s.u)
            constrained_set.add(s.v)

    K = np.zeros((n_v, n_v), dtype=float)
    f_glob = np.zeros(n_v, dtype=float)
    tri_B: list[np.ndarray] = []

    for a, b, c in tris:
        x0, y0 = verts[a]
        x1, y1 = verts[b]
        x2, y2 = verts[c]
        area2 = _orient2d(x0, y0, x1, y1, x2, y2)
        area = 0.5 * area2
        B = (1.0 / area2) * np.array(
            [
                [y1 - y2, y2 - y0, y0 - y1],
                [x2 - x1, x0 - x2, x1 - x0],
            ],
            dtype=float,
        )
        tri_B.append(B)
        Ke = area * (B.T @ D @ B)
        idx = [a, b, c]
        free_in_tri = [v_i for v_i in idx if v_i not in constrained_set]
        if free_in_tri:
            f_share = (q * area) / float(len(free_in_tri))
            for v_i in free_in_tri:
                f_glob[v_i] += f_share
        else:
            for v_i in idx:
                f_glob[v_i] += q * area / 3.0
        for i_loc in range(3):
            for j_loc in range(3):
                K[idx[i_loc], idx[j_loc]] += Ke[i_loc, j_loc]

    free_dofs = [i for i in range(n_v) if i not in constrained_set]
    w = np.zeros(n_v, dtype=float)
    if free_dofs:
        ix = np.ix_(free_dofs, free_dofs)
        w_free = np.linalg.solve(K[ix], f_glob[free_dofs])
        for i_f, dof in enumerate(free_dofs):
            w[dof] = w_free[i_f]

    nodal_R = f_glob - K @ w

    edge_to_tris: dict[tuple[int, int], list[int]] = {}
    for t_idx, (a, b, c) in enumerate(tris):
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            edge_to_tris.setdefault(key, []).append(t_idx)

    sub_flux_weights: dict[tuple[int, int], float] = {}
    for k, subs in enumerate(frame_sub_groups):
        for j, s in enumerate(subs):
            key = (min(s.u, s.v), max(s.u, s.v))
            t_list = edge_to_tris.get(key, [])
            flux_mag = 0.0
            mx_s = 0.5 * (verts[s.u, 0] + verts[s.v, 0])
            my_s = 0.5 * (verts[s.u, 1] + verts[s.v, 1])
            for t_idx in t_list:
                a, b, c = tris[t_idx]
                grad_w = tri_B[t_idx] @ np.array([w[a], w[b], w[c]], dtype=float)
                q_flux = -(D @ grad_w)
                tx = (verts[a, 0] + verts[b, 0] + verts[c, 0]) / 3.0
                ty = (verts[a, 1] + verts[b, 1] + verts[c, 1]) / 3.0
                nx, ny = s.nx, s.ny
                if (mx_s - tx) * nx + (my_s - ty) * ny < 0.0:
                    nx, ny = -nx, -ny
                fn = q_flux[0] * nx + q_flux[1] * ny
                if abs(q) > 1e-15:
                    fn_signed = fn if q > 0.0 else -fn
                    flux_mag += max(0.0, fn_signed)
            n_vec = np.array([s.nx, s.ny], dtype=float)
            d_nn = float(n_vec @ D @ n_vec)
            sub_flux_weights[(k, j)] = s.length * (flux_mag + 1e-12 * d_nn)

    _apportion_nodal_reactions_to_subsegments(
        mesh=mesh,
        frames=frames,
        frame_sub_groups=frame_sub_groups,
        sub_F=sub_F,
        sub_M=sub_M,
        node_reactions=node_reactions,
        nodal_R=nodal_R,
        sub_flux_weights=sub_flux_weights,
    )

    frame_results = _finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M)
    return DistributionResult(
        mode="fea_poisson",
        total_applied_load=total_load,
        total_applied_mx=total_mx,
        total_applied_my=total_my,
        frame_loads=frame_results,
        support_node_reactions=node_reactions,
        mesh=mesh,
    )


def _dkt_element_stiffness(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    x3: float,
    y3: float,
    Db: np.ndarray,
) -> np.ndarray:
    """9x9 Batoz-Bathe-Ho (1980) DKT plate bending stiffness matrix for DOFs
    [w1, bx1, by1, w2, bx2, by2, w3, bx3, by3] where bx = dw/dx = -thy, by = dw/dy = +thx.
    """
    x23 = x2 - x3
    y23 = y2 - y3
    x31 = x3 - x1
    y31 = y3 - y1
    x12 = x1 - x2
    y12 = y1 - y2

    area2 = x31 * y12 - x12 * y31
    area = 0.5 * area2

    # Edges 4=(2->3), 5=(3->1), 6=(1->2)
    l4_2 = x23 * x23 + y23 * y23
    l5_2 = x31 * x31 + y31 * y31
    l6_2 = x12 * x12 + y12 * y12

    a4, a5, a6 = -x23 / l4_2, -x31 / l5_2, -x12 / l6_2
    b4, b5, b6 = 0.75 * x23 * y23 / l4_2, 0.75 * x31 * y31 / l5_2, 0.75 * x12 * y12 / l6_2
    c4 = (-0.25 * x23 * x23 + 0.5 * y23 * y23) / l4_2
    c5 = (-0.25 * x31 * x31 + 0.5 * y31 * y31) / l5_2
    c6 = (-0.25 * x12 * x12 + 0.5 * y12 * y12) / l6_2
    d4, d5, d6 = -y23 / l4_2, -y31 / l5_2, -y12 / l6_2
    e4 = (0.5 * x23 * x23 - 0.25 * y23 * y23) / l4_2
    e5 = (0.5 * x31 * x31 - 0.25 * y31 * y31) / l5_2
    e6 = (0.5 * x12 * x12 - 0.25 * y12 * y12) / l6_2

    # Nodal DOFs per node i: [w_i, bx_i, by_i] where bx_i = dw/dx, by_i = dw/dy
    Tx = np.zeros((6, 9), dtype=float)
    Ty = np.zeros((6, 9), dtype=float)

    # P1..P3 at corners 1, 2, 3
    Tx[0, 1] = 1.0
    Ty[0, 2] = 1.0
    Tx[1, 4] = 1.0
    Ty[1, 5] = 1.0
    Tx[2, 7] = 1.0
    Ty[2, 8] = 1.0

    # P4 on edge 4 (2->3):
    Tx[3, :] = [0.0, 0.0, 0.0, -1.5 * a4, c4, -b4, 1.5 * a4, c4, -b4]
    Ty[3, :] = [0.0, 0.0, 0.0, -1.5 * d4, -b4, e4, 1.5 * d4, -b4, e4]

    # P5 on edge 5 (3->1):
    Tx[4, :] = [1.5 * a5, c5, -b5, 0.0, 0.0, 0.0, -1.5 * a5, c5, -b5]
    Ty[4, :] = [1.5 * d5, -b5, e5, 0.0, 0.0, 0.0, -1.5 * d5, -b5, e5]

    # P6 on edge 6 (1->2):
    Tx[5, :] = [-1.5 * a6, c6, -b6, 1.5 * a6, c6, -b6, 0.0, 0.0, 0.0]
    Ty[5, :] = [-1.5 * d6, -b6, e6, 1.5 * d6, -b6, e6, 0.0, 0.0, 0.0]

    Ke = np.zeros((9, 9), dtype=float)
    for (L1, L2, L3), wg in zip(_TRI3_BARY, _TRI3_W):
        dP_dL1 = np.array(
            [
                4.0 * L1 - 1.0,
                0.0,
                -(4.0 * L3 - 1.0),
                -4.0 * L2,
                4.0 * (L3 - L1),
                4.0 * L2,
            ],
            dtype=float,
        )
        dP_dL2 = np.array(
            [
                0.0,
                4.0 * L2 - 1.0,
                -(4.0 * L3 - 1.0),
                4.0 * (L3 - L2),
                -4.0 * L1,
                4.0 * L1,
            ],
            dtype=float,
        )

        dHx_dL1 = dP_dL1 @ Tx
        dHx_dL2 = dP_dL2 @ Tx
        dHy_dL1 = dP_dL1 @ Ty
        dHy_dL2 = dP_dL2 @ Ty

        B = np.zeros((3, 9), dtype=float)
        B[0, :] = (y23 * dHx_dL1 + y31 * dHx_dL2) / area2
        B[1, :] = (-x23 * dHy_dL1 - x31 * dHy_dL2) / area2
        B[2, :] = (-x23 * dHx_dL1 - x31 * dHx_dL2 + y23 * dHy_dL1 + y31 * dHy_dL2) / area2

        Ke += (area * wg) * (B.T @ Db @ B)

    return Ke


def _distribute_fea_plate(
    mesh: PSLGMesh,
    q: float,
    frames: Sequence[tuple[tuple[float, float], tuple[float, float]]],
    frame_sub_groups: list[list[SubSegment]],
    sub_F: list[list[float]],
    sub_M: list[list[float]],
    node_reactions: np.ndarray,
    total_load: float,
    total_mx: float,
    total_my: float,
    nu: float = 0.2,
) -> DistributionResult:
    verts = mesh.vertices
    tris = mesh.triangles
    n_v = len(verts)
    ndof = 3 * n_v

    Db = (1.0 / (1.0 - nu * nu)) * np.array(
        [
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ],
        dtype=float,
    )

    K = np.zeros((ndof, ndof), dtype=float)
    f_glob = np.zeros(ndof, dtype=float)

    for a, b, c in tris:
        x1, y1 = verts[a]
        x2, y2 = verts[b]
        x3, y3 = verts[c]
        area = 0.5 * _orient2d(x1, y1, x2, y2, x3, y3)
        Ke = _dkt_element_stiffness(x1, y1, x2, y2, x3, y3, Db)
        dofs = [
            3 * a,
            3 * a + 1,
            3 * a + 2,
            3 * b,
            3 * b + 1,
            3 * b + 2,
            3 * c,
            3 * c + 1,
            3 * c + 2,
        ]
        for loc_n in (a, b, c):
            f_glob[3 * loc_n] += q * area / 3.0
        for i_l in range(9):
            di = dofs[i_l]
            for j_l in range(9):
                K[di, dofs[j_l]] += Ke[i_l, j_l]

    supp_vert_set: set[int] = set(mesh.support_node_vertex_indices)
    for subs in frame_sub_groups:
        for s in subs:
            supp_vert_set.add(s.u)
            supp_vert_set.add(s.v)
    supp_verts = sorted(supp_vert_set)

    supp_coords = verts[supp_verts]
    centered = supp_coords - np.mean(supp_coords, axis=0)
    cov = centered.T @ centered
    eigvals = np.linalg.eigvalsh(cov)
    sigma_max = math.sqrt(max(0.0, float(eigvals[-1])))
    sigma_min = math.sqrt(max(0.0, float(eigvals[0])))
    support_rank = 2 if (len(supp_verts) >= 3 and sigma_min >= 1e-5 * max(sigma_max, mesh.char_length)) else 1

    constrained_dofs: set[int] = set()
    for v_i in supp_verts:
        constrained_dofs.add(3 * v_i)  # w = 0
        if support_rank < 2:
            constrained_dofs.add(3 * v_i + 1)
            constrained_dofs.add(3 * v_i + 2)

    free_dofs = [d for d in range(ndof) if d not in constrained_dofs]
    u = np.zeros(ndof, dtype=float)
    if free_dofs:
        ix = np.ix_(free_dofs, free_dofs)
        u_free = np.linalg.solve(K[ix], f_glob[free_dofs])
        for idx_f, d in enumerate(free_dofs):
            u[d] = u_free[idx_f]

    R_full = f_glob - K @ u
    nodal_Rw = R_full[0::3]

    sub_flux_weights: dict[tuple[int, int], float] = {}
    for k, subs in enumerate(frame_sub_groups):
        for j, s in enumerate(subs):
            sub_flux_weights[(k, j)] = s.length

    _apportion_nodal_reactions_to_subsegments(
        mesh=mesh,
        frames=frames,
        frame_sub_groups=frame_sub_groups,
        sub_F=sub_F,
        sub_M=sub_M,
        node_reactions=node_reactions,
        nodal_R=nodal_Rw,
        sub_flux_weights=sub_flux_weights,
    )

    overhang_torques = [0.0] * len(frames)
    if support_rank < 2:
        # In _dkt_element_stiffness, DOF 1 is bx = dw/dx = -thy and DOF 2 is by = dw/dy = +thx.
        # For rigid rotations (thx, thy), w = y*thx - x*thy => bx = -thy, by = +thx.
        # Thus reaction moments about (+x, +y) are Mx = +R_by (DOF 2) and My = -R_bx (DOF 1)!
        for k, (fp1, fp2) in enumerate(frames):
            subs = frame_sub_groups[k]
            if not subs:
                continue
            dx = fp2[0] - fp1[0]
            dy = fp2[1] - fp1[1]
            lk = math.hypot(dx, dy)
            if lk <= 1e-15:
                continue
            tx, ty = dx / lk, dy / lk
            f_verts = {s.u for s in subs} | {s.v for s in subs}
            torq_k = 0.0
            for v_i in f_verts:
                mx_i = float(R_full[3 * v_i + 2])
                my_i = -float(R_full[3 * v_i + 1])
                torq_k += mx_i * tx + my_i * ty
            overhang_torques[k] = torq_k

    frame_results = _finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M, overhang_torques)
    return DistributionResult(
        mode="fea_plate",
        total_applied_load=total_load,
        total_applied_mx=total_mx,
        total_applied_my=total_my,
        frame_loads=frame_results,
        support_node_reactions=node_reactions,
        mesh=mesh,
        plate_support_rank=support_rank,
    )


# ============================================================================
# 4. Guyan Static Condensation of In-Plane CST Diaphragm Stiffness & Mass
# ============================================================================


def condense_inplane_cst(
    vertices: np.ndarray,
    triangles: np.ndarray,
    boundary_node_indices: Sequence[int],
    E: float,
    nu: float,
    t: float,
    rho: float = 0.0,
) -> dict[str, np.ndarray]:
    """Assemble 2D plane-stress CST stiffness K and consistent mass M on (vertices, triangles)
    and statically condense interior Steiner nodes onto boundary_node_indices via Guyan reduction.
    Steiner vertices lying on the outer polygon edges between consecutive boundary_node_indices
    are kinematically interpolated from the edge endpoints so the condensed N-node polygon
    satisfies the constant-strain patch test to machine precision.
    """
    n_v = len(vertices)
    ndof = 2 * n_v
    D_ps = (E / (1.0 - nu * nu)) * np.array(
        [
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ],
        dtype=float,
    )

    K = np.zeros((ndof, ndof), dtype=float)
    M = np.zeros((ndof, ndof), dtype=float)

    for a, b, c in triangles:
        x0, y0 = vertices[a]
        x1, y1 = vertices[b]
        x2, y2 = vertices[c]
        area2 = _orient2d(x0, y0, x1, y1, x2, y2)
        area = 0.5 * area2
        b0, b1, b2 = (y1 - y2) / area2, (y2 - y0) / area2, (y0 - y1) / area2
        c0, c1, c2 = (x2 - x1) / area2, (x0 - x2) / area2, (x1 - x0) / area2

        B = np.array(
            [
                [b0, 0.0, b1, 0.0, b2, 0.0],
                [0.0, c0, 0.0, c1, 0.0, c2],
                [c0, b0, c1, b1, c2, b2],
            ],
            dtype=float,
        )
        Ke = (t * area) * (B.T @ D_ps @ B)
        Me = np.zeros((6, 6), dtype=float)
        if rho > 0.0:
            m_fac = rho * t * area / 12.0
            for i_n in range(3):
                for j_n in range(3):
                    val = m_fac * (2.0 if i_n == j_n else 1.0)
                    Me[2 * i_n, 2 * j_n] = val
                    Me[2 * i_n + 1, 2 * j_n + 1] = val

        dofs = [2 * a, 2 * a + 1, 2 * b, 2 * b + 1, 2 * c, 2 * c + 1]
        for i_l in range(6):
            di = dofs[i_l]
            for j_l in range(6):
                dj = dofs[j_l]
                K[di, dj] += Ke[i_l, j_l]
                M[di, dj] += Me[i_l, j_l]

    b_nodes = list(boundary_node_indices)
    n_b = len(b_nodes)
    b_pos = {node_id: idx for idx, node_id in enumerate(b_nodes)}

    # Classify non-corner vertices into edge-slaved boundary Steiner nodes vs strictly interior nodes
    edge_slaves: dict[int, tuple[int, int, float]] = {}
    i_nodes: list[int] = []

    for v in range(n_v):
        if v in b_pos:
            continue
        vx, vy = float(vertices[v, 0]), float(vertices[v, 1])
        slaved = False
        for idx_b in range(n_b):
            u_b = b_nodes[idx_b]
            w_b = b_nodes[(idx_b + 1) % n_b]
            ax, ay = float(vertices[u_b, 0]), float(vertices[u_b, 1])
            bx, by = float(vertices[w_b, 0]), float(vertices[w_b, 1])
            seg_len = math.hypot(bx - ax, by - ay)
            dist, s_param, is_int = _point_segment_distance(vx, vy, ax, ay, bx, by)
            if is_int and dist <= 1e-7 * max(1.0, seg_len):
                edge_slaves[v] = (idx_b, (idx_b + 1) % n_b, s_param)
                slaved = True
                break
        if not slaved:
            i_nodes.append(v)

    # Build kinematic master-slave matrix P mapping [b_dofs; i_dofs] -> all 2*n_v DOFs
    n_i = len(i_nodes)
    i_pos = {node_id: idx for idx, node_id in enumerate(i_nodes)}
    n_master = 2 * (n_b + n_i)
    P = np.zeros((ndof, n_master), dtype=float)

    for v in range(n_v):
        if v in b_pos:
            col = 2 * b_pos[v]
            P[2 * v, col] = 1.0
            P[2 * v + 1, col + 1] = 1.0
        elif v in edge_slaves:
            ib0, ib1, s_val = edge_slaves[v]
            c0 = 2 * ib0
            c1 = 2 * ib1
            P[2 * v, c0] = 1.0 - s_val
            P[2 * v + 1, c0 + 1] = 1.0 - s_val
            P[2 * v, c1] = s_val
            P[2 * v + 1, c1 + 1] = s_val
        else:
            col = 2 * (n_b + i_pos[v])
            P[2 * v, col] = 1.0
            P[2 * v + 1, col + 1] = 1.0

    K_red = P.T @ K @ P
    M_red = P.T @ M @ P

    nb_dof = 2 * n_b
    K_bb = K_red[:nb_dof, :nb_dof]
    M_bb = M_red[:nb_dof, :nb_dof]
    if n_i == 0:
        K_bb = 0.5 * (K_bb + K_bb.T)
        M_bb = 0.5 * (M_bb + M_bb.T)
        return {"K_cond": K_bb, "M_cond": M_bb, "G": np.zeros((0, nb_dof))}

    K_bi = K_red[:nb_dof, nb_dof:]
    K_ib = K_red[nb_dof:, :nb_dof]
    K_ii = K_red[nb_dof:, nb_dof:]

    M_bi = M_red[:nb_dof, nb_dof:]
    M_ib = M_red[nb_dof:, :nb_dof]
    M_ii = M_red[nb_dof:, nb_dof:]

    G = -np.linalg.solve(K_ii, K_ib)
    K_cond = K_bb + K_bi @ G
    K_cond = 0.5 * (K_cond + K_cond.T)

    M_cond = M_bb + M_bi @ G + G.T @ M_ib + G.T @ M_ii @ G
    M_cond = 0.5 * (M_cond + M_cond.T)

    return {"K_cond": K_cond, "M_cond": M_cond, "G": G}

