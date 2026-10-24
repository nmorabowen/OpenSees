// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

#ifndef LadrunoMembraneKernel_h
#define LadrunoMembraneKernel_h

#ifndef NOMINMAX
#define NOMINMAX
#endif
#if defined(_WIN32)
#undef min
#undef max
#endif

// ============================================================================
// LadrunoMembraneKernel.h — Zero-dependency C++17 Computational Geometry &
// Surface-to-Frame Load Distribution + In-Plane CST Guyan Condensation Kernel
// (Normative implementation of ADR-98)
// ============================================================================

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace LadrunoMembraneKernel {

struct Point2D {
    double x = 0.0;
    double y = 0.0;
};

struct Segment2D {
    Point2D p1;
    Point2D p2;
};

// ============================================================================
// Dense Linear Algebra Helpers (Cholesky / Gaussian Elimination with Pivoting)
// ============================================================================

inline std::vector<double> solve_linear_system(
    int n,
    const std::vector<double>& A_row_major,
    const std::vector<double>& b_rhs,
    int nrhs = 1)
{
    if (n == 0) {
        return {};
    }
    std::vector<double> LU = A_row_major;
    std::vector<int> piv_order(n);
    for (int i = 0; i < n; ++i) {
        piv_order[i] = i;
    }

    for (int k = 0; k < n; ++k) {
        int piv = k;
        double max_val = std::abs(LU[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            double v = std::abs(LU[i * n + k]);
            if (v > max_val) {
                max_val = v;
                piv = i;
            }
        }
        if (max_val <= 1e-30) {
            throw std::runtime_error("Singular matrix in LadrunoMembrane::solve_linear_system.");
        }
        if (piv != k) {
            std::swap(piv_order[k], piv_order[piv]);
            for (int j = 0; j < n; ++j) {
                std::swap(LU[k * n + j], LU[piv * n + j]);
            }
        }
        double diag = LU[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            double fac = LU[i * n + k] / diag;
            LU[i * n + k] = fac;
            if (std::abs(fac) > 0.0) {
                for (int j = k + 1; j < n; ++j) {
                    LU[i * n + j] -= fac * LU[k * n + j];
                }
            }
        }
    }

    auto lu_solve = [&](const std::vector<double>& rhs) -> std::vector<double> {
        std::vector<double> Y(n * nrhs, 0.0);
        for (int i = 0; i < n; ++i) {
            int p_i = piv_order[i];
            for (int c = 0; c < nrhs; ++c) {
                double s = rhs[p_i * nrhs + c];
                for (int j = 0; j < i; ++j) {
                    s -= LU[i * n + j] * Y[j * nrhs + c];
                }
                Y[i * nrhs + c] = s;
            }
        }
        for (int i = n - 1; i >= 0; --i) {
            double diag = LU[i * n + i];
            for (int c = 0; c < nrhs; ++c) {
                double s = Y[i * nrhs + c];
                for (int j = i + 1; j < n; ++j) {
                    s -= LU[i * n + j] * Y[j * nrhs + c];
                }
                Y[i * nrhs + c] = s / diag;
            }
        }
        return Y;
    };

    std::vector<double> X = lu_solve(b_rhs);

    // One step of iterative refinement to eliminate O(kappa * eps) drift on 4th-order DKT matrices
    std::vector<double> R(n * nrhs, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < nrhs; ++c) {
            long double acc = static_cast<long double>(b_rhs[i * nrhs + c]);
            for (int j = 0; j < n; ++j) {
                acc -= static_cast<long double>(A_row_major[i * n + j]) * static_cast<long double>(X[j * nrhs + c]);
            }
            R[i * nrhs + c] = static_cast<double>(acc);
        }
    }
    std::vector<double> dX = lu_solve(R);
    for (std::size_t idx = 0; idx < X.size(); ++idx) {
        X[idx] += dX[idx];
    }
    return X;
}

inline std::pair<double, double> sym_2x2_eigenvalues(double a00, double a01, double a11)
{
    double tr = a00 + a11;
    double diff = a00 - a11;
    double disc = std::hypot(diff, 2.0 * a01);
    double l_min = 0.5 * (tr - disc);
    double l_max = 0.5 * (tr + disc);
    return {l_min, l_max};
}

// ============================================================================
// 1. Exact Partial Trapezoidal Fixed-End Forces & Section Response
// ============================================================================

struct FEF2DResult {
    std::array<double, 3> p0{};
    std::array<double, 3> q0{};
    std::array<double, 6> ql0{};
};

struct FEF3DResult {
    std::array<double, 5> p0{};
    std::array<double, 5> q0{};
    std::array<double, 12> ql0{};
};

struct SectionSP3D {
    double P = 0.0;
    double Mz = 0.0;
    double Vy = 0.0;
    double My = 0.0;
    double Vz = 0.0;
};

inline std::pair<double, double> apply_end_release(double m1_fixed, double m2_fixed, int release)
{
    if (release == 1) {
        return {0.0, m2_fixed - 0.5 * m1_fixed};
    }
    if (release == 2) {
        return {m1_fixed - 0.5 * m2_fixed, 0.0};
    }
    if (release == 3) {
        return {0.0, 0.0};
    }
    return {m1_fixed, m2_fixed};
}

inline FEF2DResult partial_trapezoid_fef_2d(
    double L,
    double wya,
    double wyb,
    double waa,
    double wab,
    double aL,
    double bL,
    int release = 0,
    double phi = 0.0)
{
    FEF2DResult res{};
    double dxi = bL - aL;
    if (L <= 0.0 || dxi <= 0.0) {
        return res;
    }

    const double s06 = std::sqrt(0.6);
    const std::array<double, 3> gl_u = {0.5 * (1.0 - s06), 0.5, 0.5 * (1.0 + s06)};
    const std::array<double, 3> gl_w = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};

    double vy1_ss = 0.0;
    double vy2_ss = 0.0;
    double mz1_fix = 0.0;
    double mz2_fix = 0.0;
    double p_ax = 0.0;
    double pj_ax = 0.0;
    double inv_1p_phi = 1.0 / (1.0 + phi);

    for (int g = 0; g < 3; ++g) {
        double ug = gl_u[g];
        double wg = gl_w[g];
        double xi = aL + dxi * ug;
        double wy = wya + (wyb - wya) * ug;
        double wa = waa + (wab - waa) * ug;
        double ds = L * dxi * wg;

        vy1_ss += wy * (1.0 - xi) * ds;
        vy2_ss += wy * xi * ds;
        mz1_fix -= L * wy * xi * (1.0 - xi) * (1.0 - xi + 0.5 * phi) * inv_1p_phi * ds;
        mz2_fix += L * wy * xi * (1.0 - xi) * (xi + 0.5 * phi) * inv_1p_phi * ds;

        p_ax += wa * ds;
        pj_ax += wa * xi * ds;
    }

    auto [mz1, mz2] = apply_end_release(mz1_fix, mz2_fix, release);

    res.p0 = {-p_ax, -vy1_ss, -vy2_ss};
    res.q0 = {-pj_ax, mz1, mz2};

    double v_mom = (mz1 + mz2) / L;
    res.ql0 = {
        -(p_ax - pj_ax),
        -vy1_ss + v_mom,
        mz1,
        -pj_ax,
        -vy2_ss - v_mom,
        mz2
    };
    return res;
}

inline FEF3DResult partial_trapezoid_fef_3d(
    double L,
    double wya,
    double wyb,
    double wza,
    double wzb,
    double waa,
    double wab,
    double aL,
    double bL,
    int releasez = 0,
    int releasey = 0,
    double phiz = 0.0,
    double phiy = 0.0)
{
    FEF3DResult res{};
    double dxi = bL - aL;
    if (L <= 0.0 || dxi <= 0.0) {
        return res;
    }

    FEF2DResult f2d_y = partial_trapezoid_fef_2d(L, wya, wyb, waa, wab, aL, bL, releasez, phiz);
    FEF2DResult f2d_z = partial_trapezoid_fef_2d(L, wza, wzb, 0.0, 0.0, aL, bL, 0, phiy);

    double my1_fix = -f2d_z.q0[1];
    double my2_fix = -f2d_z.q0[2];
    auto [my1, my2] = apply_end_release(my1_fix, my2_fix, releasey);

    double p_ax = -f2d_y.p0[0];
    double pj_ax = -f2d_y.q0[0];
    double vy1_ss = -f2d_y.p0[1];
    double vy2_ss = -f2d_y.p0[2];
    double vz1_ss = -f2d_z.p0[1];
    double vz2_ss = -f2d_z.p0[2];
    double mz1 = f2d_y.q0[1];
    double mz2 = f2d_y.q0[2];

    res.p0 = {-p_ax, -vy1_ss, -vy2_ss, -vz1_ss, -vz2_ss};
    res.q0 = {-pj_ax, mz1, mz2, my1, my2};

    double vy_mom = (mz1 + mz2) / L;
    double vz_mom = -(my1 + my2) / L;

    res.ql0 = {
        -(p_ax - pj_ax),
        -vy1_ss + vy_mom,
        -vz1_ss + vz_mom,
        0.0,
        my1,
        mz1,
        -pj_ax,
        -vy2_ss - vy_mom,
        -vz2_ss - vz_mom,
        0.0,
        my2,
        mz2
    };
    return res;
}

inline SectionSP3D force_beam_section_sp_3d(
    double x,
    double L,
    double wya,
    double wyb,
    double wza,
    double wzb,
    double waa,
    double wab,
    double aL,
    double bL)
{
    SectionSP3D sp{};
    double a = aL * L;
    double b = bL * L;
    double ba = b - a;
    if (L <= 0.0 || ba <= 0.0) {
        return sp;
    }

    double fa = waa * ba + 0.5 * (wab - waa) * ba;
    double fy_rect = wya * ba;
    double fz_rect = wza * ba;
    double c_rect = a + 0.5 * ba;
    double vy_i = fy_rect * (1.0 - c_rect / L);
    double vy_j = fy_rect * (c_rect / L);
    double vz_i = fz_rect * (1.0 - c_rect / L);
    double vz_j = fz_rect * (c_rect / L);

    double fy_tri = 0.5 * (wyb - wya) * ba;
    double fz_tri = 0.5 * (wzb - wza) * ba;
    double c_tri = a + (2.0 / 3.0) * ba;
    vy_i += fy_tri * (1.0 - c_tri / L);
    vy_j += fy_tri * (c_tri / L);
    vz_i += fz_tri * (1.0 - c_tri / L);
    vz_j += fz_tri * (c_tri / L);

    if (x <= a) {
        sp.P = fa;
        sp.Mz = -vy_i * x;
        sp.Vy = -vy_i;
        sp.My = vz_i * x;
        sp.Vz = vz_i;
        return sp;
    }
    if (x >= b) {
        sp.P = 0.0;
        sp.Mz = vy_j * (x - L);
        sp.Vy = vy_j;
        sp.My = -vz_j * (x - L);
        sp.Vz = -vz_j;
        return sp;
    }

    double dx = x - a;
    double wyy = wya + (wyb - wya) / ba * dx;
    double wzz = wza + (wzb - wza) / ba * dx;
    double waa_x = waa + (wab - waa) / ba * dx;

    sp.P = fa - waa * dx - 0.5 * (waa_x - waa) * dx;
    sp.Mz = -vy_i * x + 0.5 * wya * dx * dx + (1.0 / 6.0) * (wyy - wya) * dx * dx;
    sp.Vy = -vy_i + wya * dx + 0.5 * (wyy - wya) * dx;
    sp.My = vz_i * x - 0.5 * wza * dx * dx - (1.0 / 6.0) * (wzz - wza) * dx * dx;
    sp.Vz = vz_i - wza * dx - 0.5 * (wzz - wza) * dx;
    return sp;
}

inline std::pair<double, double> moment_conserving_trapezoid(double F_seg, double M_seg, double L_seg)
{
    if (L_seg <= 1e-15) {
        return {0.0, 0.0};
    }
    double w_sum = 2.0 * F_seg / L_seg;
    double w_diff = 12.0 * (M_seg - 0.5 * F_seg * L_seg) / (L_seg * L_seg);
    return {0.5 * (w_sum - w_diff), 0.5 * (w_sum + w_diff)};
}

// ============================================================================
// 2. 4-Stage PSLG Preprocessor + CDT + Chew-Ruppert Quality Mesher
// ============================================================================

inline double orient2d(double ax, double ay, double bx, double by, double cx, double cy)
{
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

inline double incircle_ccw(
    double ax, double ay,
    double bx, double by,
    double cx, double cy,
    double dx, double dy)
{
    double adx = ax - dx;
    double ady = ay - dy;
    double bdx = bx - dx;
    double bdy = by - dy;
    double cdx = cx - dx;
    double cdy = cy - dy;
    double abdet = adx * bdy - bdx * ady;
    double bcdet = bdx * cdy - cdx * bdy;
    double cadet = cdx * ady - adx * cdy;
    double alift = adx * adx + ady * ady;
    double blift = bdx * bdx + bdy * bdy;
    double clift = cdx * cdx + cdy * cdy;
    return alift * bcdet + blift * cadet + clift * abdet;
}

inline double polygon_signed_area(const std::vector<Point2D>& poly)
{
    std::size_t n = poly.size();
    if (n < 3) {
        return 0.0;
    }
    double area2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p1 = poly[i];
        const auto& p2 = poly[(i + 1) % n];
        area2 += p1.x * p2.y - p2.x * p1.y;
    }
    return 0.5 * area2;
}

inline bool point_in_polygon(Point2D pt, const std::vector<Point2D>& poly)
{
    bool inside = false;
    std::size_t n = poly.size();
    for (std::size_t i = 0; i < n; ++i) {
        double x1 = poly[i].x;
        double y1 = poly[i].y;
        double x2 = poly[(i + 1) % n].x;
        double y2 = poly[(i + 1) % n].y;
        if ((y1 > pt.y) != (y2 > pt.y)) {
            double x_int = x1 + (pt.y - y1) * (x2 - x1) / (y2 - y1);
            if (pt.x < x_int) {
                inside = !inside;
            }
        }
    }
    return inside;
}

struct PointSegDist {
    double dist = 0.0;
    double t = 0.0;
    bool is_interior = false;
};

inline PointSegDist point_segment_distance(
    double px, double py,
    double ax, double ay,
    double bx, double by)
{
    double dx = bx - ax;
    double dy = by - ay;
    double l2 = dx * dx + dy * dy;
    if (l2 <= 1e-30) {
        return {std::hypot(px - ax, py - ay), 0.0, false};
    }
    double t_raw = ((px - ax) * dx + (py - ay) * dy) / l2;
    bool is_int = (t_raw > 1e-9 && t_raw < 1.0 - 1e-9);
    double t = std::max(0.0, std::min(1.0, t_raw));
    double cx = ax + t * dx;
    double cy = ay + t * dy;
    return {std::hypot(px - cx, py - cy), t, is_int};
}

struct SegIntersection {
    bool hit = false;
    double ix = 0.0;
    double iy = 0.0;
    double t = 0.0;
    double u = 0.0;
};

inline SegIntersection segments_properly_intersect(
    Point2D a1, Point2D a2,
    Point2D b1, Point2D b2,
    double eps = 1e-9)
{
    double dx1 = a2.x - a1.x;
    double dy1 = a2.y - a1.y;
    double dx2 = b2.x - b1.x;
    double dy2 = b2.y - b1.y;
    double denom = dx1 * dy2 - dy1 * dx2;
    double len1 = std::hypot(dx1, dy1);
    double len2 = std::hypot(dx2, dy2);
    if (len1 <= eps || len2 <= eps) {
        return {};
    }
    if (std::abs(denom) <= eps * len1 * len2) {
        return {};
    }
    double dx3 = b1.x - a1.x;
    double dy3 = b1.y - a1.y;
    double t = (dx3 * dy2 - dy3 * dx2) / denom;
    double u = (dx3 * dy1 - dy3 * dx1) / denom;
    if (t >= -1e-7 && t <= 1.0 + 1e-7 && u >= -1e-7 && u <= 1.0 + 1e-7) {
        double t_c = std::max(0.0, std::min(1.0, t));
        double u_c = std::max(0.0, std::min(1.0, u));
        if ((t_c > 1e-7 && t_c < 1.0 - 1e-7) || (u_c > 1e-7 && u_c < 1.0 - 1e-7)) {
            return {true, a1.x + t_c * dx1, a1.y + t_c * dy1, t_c, u_c};
        }
    }
    return {};
}

struct SubSegment {
    int u = 0;
    int v = 0;
    int frame_id = -1;
    double t_start = 0.0;
    double t_end = 1.0;
    double length = 0.0;
    double nx = 0.0;
    double ny = 0.0;
};

struct PSLGMesh {
    std::vector<Point2D> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::vector<SubSegment> subsegments;
    std::vector<bool> is_acute;
    std::vector<int> boundary_vertex_indices;
    std::vector<int> support_node_vertex_indices;
    double total_area = 0.0;
    double char_length = 0.0;
};

inline PSLGMesh build_pslg_cdt_mesh(
    const std::vector<Point2D>& outer_poly,
    const std::vector<std::vector<Point2D>>& holes = {},
    const std::vector<Segment2D>& frames = {},
    const std::vector<Point2D>& support_nodes = {},
    double h_target = 0.0,
    int max_vertices = 1200)
{
    if (outer_poly.size() < 3) {
        throw std::invalid_argument("Degenerate outer polygon: fewer than 3 vertices.");
    }

    std::vector<Point2D> raw_outer = outer_poly;
    double area_signed = polygon_signed_area(raw_outer);
    if (std::abs(area_signed) <= 1e-12) {
        throw std::invalid_argument("Degenerate outer polygon: zero or near-zero area.");
    }

    std::vector<Point2D> outer_ccw = raw_outer;
    if (area_signed < 0.0) {
        std::reverse(outer_ccw.begin(), outer_ccw.end());
    }
    double outer_area = std::abs(area_signed);
    double char_len = std::sqrt(outer_area);
    double snap_eps = 1e-8 * char_len;

    std::size_t n_out = outer_ccw.size();
    for (std::size_t i = 0; i < n_out; ++i) {
        Point2D p1 = outer_ccw[i];
        Point2D p2 = outer_ccw[(i + 1) % n_out];
        if (std::hypot(p2.x - p1.x, p2.y - p1.y) <= snap_eps) {
            throw std::invalid_argument("Degenerate outer polygon: coincident consecutive vertices.");
        }
        for (std::size_t j = i + 2; j < n_out; ++j) {
            if (i == 0 && j == n_out - 1) {
                continue;
            }
            Point2D q1 = outer_ccw[j];
            Point2D q2 = outer_ccw[(j + 1) % n_out];
            if (segments_properly_intersect(p1, p2, q1, q2, 1e-9).hit) {
                throw std::invalid_argument("Degenerate outer polygon: self-intersecting boundary edges.");
            }
        }
    }

    double net_area = outer_area;
    std::vector<std::vector<Point2D>> holes_ccw;
    for (const auto& h_poly : holes) {
        if (h_poly.size() < 3) {
            throw std::invalid_argument("Degenerate hole polygon: fewer than 3 vertices.");
        }
        double h_area = std::abs(polygon_signed_area(h_poly));
        if (h_area <= 1e-12) {
            throw std::invalid_argument("Degenerate hole polygon: zero area.");
        }
        for (const auto& hp : h_poly) {
            if (!point_in_polygon(hp, outer_ccw)) {
                throw std::invalid_argument("Degenerate hole polygon: vertex lies outside outer polygon.");
            }
        }
        net_area -= h_area;
        holes_ccw.push_back(h_poly);
    }

    if (net_area <= 1e-12) {
        throw std::invalid_argument("Degenerate polygon: net area <= 0.");
    }

    char_len = std::sqrt(net_area);
    snap_eps = 1e-8 * char_len;
    double h_min = 1e-4 * char_len;
    if (h_target <= 0.0) {
        h_target = 0.18 * char_len;
    }
    h_target = std::max(h_target, 4.0 * h_min);

    std::vector<Point2D> verts;
    auto add_vertex = [&](double px, double py) -> int {
        for (std::size_t idx = 0; idx < verts.size(); ++idx) {
            if (std::hypot(px - verts[idx].x, py - verts[idx].y) <= snap_eps) {
                return static_cast<int>(idx);
            }
        }
        verts.push_back({px, py});
        return static_cast<int>(verts.size()) - 1;
    };

    std::vector<int> boundary_vertex_indices;
    for (const auto& p : raw_outer) {
        boundary_vertex_indices.push_back(add_vertex(p.x, p.y));
    }

    std::vector<std::vector<int>> hole_vert_loops;
    for (const auto& h_pts : holes_ccw) {
        std::vector<int> h_loop;
        for (const auto& p : h_pts) {
            h_loop.push_back(add_vertex(p.x, p.y));
        }
        hole_vert_loops.push_back(h_loop);
    }

    std::vector<int> support_node_vertex_indices;
    for (const auto& p : support_nodes) {
        support_node_vertex_indices.push_back(add_vertex(p.x, p.y));
    }

    for (const auto& f : frames) {
        add_vertex(f.p1.x, f.p1.y);
        add_vertex(f.p2.x, f.p2.y);
    }

    std::vector<std::pair<int, int>> raw_Boundary_segs;
    std::vector<int> out_loop_indices;
    for (const auto& p : outer_ccw) {
        out_loop_indices.push_back(add_vertex(p.x, p.y));
    }
    for (std::size_t i = 0; i < out_loop_indices.size(); ++i) {
        raw_Boundary_segs.push_back({out_loop_indices[i], out_loop_indices[(i + 1) % out_loop_indices.size()]});
    }
    for (const auto& h_loop : hole_vert_loops) {
        for (std::size_t i = 0; i < h_loop.size(); ++i) {
            raw_Boundary_segs.push_back({h_loop[i], h_loop[(i + 1) % h_loop.size()]});
        }
    }

    std::vector<std::pair<int, int>> raw_frame_segs;
    for (const auto& f : frames) {
        int u = add_vertex(f.p1.x, f.p1.y);
        int v = add_vertex(f.p2.x, f.p2.y);
        if (u != v) {
            raw_frame_segs.push_back({u, v});
        }
    }

    std::vector<std::pair<int, int>> all_initial_segs = raw_Boundary_segs;
    all_initial_segs.insert(all_initial_segs.end(), raw_frame_segs.begin(), raw_frame_segs.end());

    // Stage 2: Segment-segment intersections
    std::size_t n_init = all_initial_segs.size();
    for (std::size_t i = 0; i < n_init; ++i) {
        int u1 = all_initial_segs[i].first;
        int v1 = all_initial_segs[i].second;
        Point2D a1 = verts[u1];
        Point2D a2 = verts[v1];
        for (std::size_t j = i + 1; j < n_init; ++j) {
            int u2 = all_initial_segs[j].first;
            int v2 = all_initial_segs[j].second;
            if (u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2) {
                continue;
            }
            Point2D b1 = verts[u2];
            Point2D b2 = verts[v2];
            auto hit = segments_properly_intersect(a1, a2, b1, b2, 1e-9);
            if (hit.hit) {
                add_vertex(hit.ix, hit.iy);
            }
        }
    }

    // Stage 3: Split every segment at any vertex lying on its interior
    auto split_segment_by_vertices = [&](int u, int v) -> std::vector<std::pair<int, int>> {
        double ax = verts[u].x;
        double ay = verts[u].y;
        double bx = verts[v].x;
        double by = verts[v].y;
        double seg_len = std::hypot(bx - ax, by - ay);
        if (seg_len <= snap_eps) {
            return {};
        }
        std::vector<std::pair<double, int>> on_seg = {{0.0, u}, {1.0, v}};
        for (std::size_t w_idx = 0; w_idx < verts.size(); ++w_idx) {
            int wi = static_cast<int>(w_idx);
            if (wi == u || wi == v) {
                continue;
            }
            auto psd = point_segment_distance(verts[wi].x, verts[wi].y, ax, ay, bx, by);
            if (psd.dist <= snap_eps && psd.t > 1e-7 && psd.t < 1.0 - 1e-7) {
                on_seg.push_back({psd.t, wi});
            }
        }
        std::stable_sort(on_seg.begin(), on_seg.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
        std::vector<std::pair<int, int>> out;
        for (std::size_t k = 0; k + 1 < on_seg.size(); ++k) {
            int ia = on_seg[k].second;
            int ib = on_seg[k + 1].second;
            if (ia != ib) {
                out.push_back({ia, ib});
            }
        }
        return out;
    };

    std::map<std::pair<int, int>, std::pair<int, int>> pslg_edges;
    auto register_atomic_edge = [&](int u, int v) {
        for (auto [ia, ib] : split_segment_by_vertices(u, v)) {
            double mx = 0.5 * (verts[ia].x + verts[ib].x);
            double my = 0.5 * (verts[ia].y + verts[ib].y);
            bool on_bdry = false;
            for (auto [bu, bv] : raw_Boundary_segs) {
                auto psd = point_segment_distance(mx, my, verts[bu].x, verts[bu].y, verts[bv].x, verts[bv].y);
                if (psd.dist <= 10.0 * snap_eps) {
                    on_bdry = true;
                    break;
                }
            }
            if (!on_bdry) {
                if (!point_in_polygon({mx, my}, outer_ccw)) {
                    continue;
                }
                bool in_hole = false;
                for (const auto& hp : holes_ccw) {
                    if (point_in_polygon({mx, my}, hp)) {
                        in_hole = true;
                        break;
                    }
                }
                if (in_hole) {
                    continue;
                }
            }
            std::pair<int, int> key = {std::min(ia, ib), std::max(ia, ib)};
            if (pslg_edges.find(key) == pslg_edges.end()) {
                pslg_edges[key] = {ia, ib};
            }
        }
    };

    for (auto [u, v] : raw_Boundary_segs) {
        register_atomic_edge(u, v);
    }
    for (auto [u, v] : raw_frame_segs) {
        register_atomic_edge(u, v);
    }

    std::vector<std::pair<int, int>> atomic_list;
    for (const auto& kv : pslg_edges) {
        atomic_list.push_back(kv.second);
    }

    // Stage 4: Acute angle (< 60 deg) detection on PSLG vertices
    std::size_t n_pslg_verts = verts.size();
    std::vector<std::vector<Point2D>> incident_dirs(n_pslg_verts);
    for (auto [u, v] : atomic_list) {
        double dx = verts[v].x - verts[u].x;
        double dy = verts[v].y - verts[u].y;
        double l_e = std::hypot(dx, dy);
        if (l_e > snap_eps) {
            incident_dirs[u].push_back({dx / l_e, dy / l_e});
            incident_dirs[v].push_back({-dx / l_e, -dy / l_e});
        }
    }

    std::vector<bool> acute_flags(n_pslg_verts, false);
    const double cos_60 = 0.5 + 1e-7;
    for (std::size_t vi = 0; vi < n_pslg_verts; ++vi) {
        const auto& dirs = incident_dirs[vi];
        for (std::size_t i = 0; i < dirs.size(); ++i) {
            for (std::size_t j = i + 1; j < dirs.size(); ++j) {
                double dot = dirs[i].x * dirs[j].x + dirs[i].y * dirs[j].y;
                if (dot > cos_60) {
                    acute_flags[vi] = true;
                }
            }
        }
    }

    std::vector<SubSegment> subsegments;
    std::set<std::pair<int, int>> constrained_edge_set;

    for (auto [u, v] : atomic_list) {
        double ax = verts[u].x;
        double ay = verts[u].y;
        double bx = verts[v].x;
        double by = verts[v].y;
        double seg_len = std::hypot(bx - ax, by - ay);
        if (seg_len <= snap_eps) {
            continue;
        }

        int matched_frame = -1;
        double fa_t = 0.0;
        double fb_t = 1.0;
        for (std::size_t f_idx = 0; f_idx < frames.size(); ++f_idx) {
            double fx1 = frames[f_idx].p1.x;
            double fy1 = frames[f_idx].p1.y;
            double fx2 = frames[f_idx].p2.x;
            double fy2 = frames[f_idx].p2.y;
            double flen = std::hypot(fx2 - fx1, fy2 - fy1);
            if (flen <= snap_eps) {
                continue;
            }
            auto da = point_segment_distance(ax, ay, fx1, fy1, fx2, fy2);
            auto db = point_segment_distance(bx, by, fx1, fy1, fx2, fy2);
            if (da.dist <= 10.0 * snap_eps && db.dist <= 10.0 * snap_eps && std::abs(db.t - da.t) > 1e-7) {
                matched_frame = static_cast<int>(f_idx);
                fa_t = da.t;
                fb_t = db.t;
                break;
            }
        }

        int n_div = std::max(1, static_cast<int>(std::ceil(seg_len / h_target)));
        int max_div = std::max(1, static_cast<int>(std::floor(seg_len / h_min)));
        n_div = std::min(n_div, max_div);

        std::vector<int> chain_verts = {u};
        for (int k = 1; k < n_div; ++k) {
            double s = static_cast<double>(k) / static_cast<double>(n_div);
            double px = ax + s * (bx - ax);
            double py = ay + s * (by - ay);
            int idx_new = add_vertex(px, py);
            while (acute_flags.size() < verts.size()) {
                acute_flags.push_back(false);
            }
            chain_verts.push_back(idx_new);
        }
        chain_verts.push_back(v);

        for (int k = 0; k < n_div; ++k) {
            int ia = chain_verts[k];
            int ib = chain_verts[k + 1];
            if (ia == ib) {
                continue;
            }
            double s0 = static_cast<double>(k) / static_cast<double>(n_div);
            double s1 = static_cast<double>(k + 1) / static_cast<double>(n_div);
            double t0 = fa_t + s0 * (fb_t - fa_t);
            double t1 = fa_t + s1 * (fb_t - fa_t);
            if (matched_frame >= 0 && t0 > t1) {
                std::swap(ia, ib);
                std::swap(t0, t1);
            }
            double x_a = verts[ia].x;
            double y_a = verts[ia].y;
            double x_b = verts[ib].x;
            double y_b = verts[ib].y;
            double l_sub = std::hypot(x_b - x_a, y_b - y_a);
            double nx = -(y_b - y_a) / l_sub;
            double ny = (x_b - x_a) / l_sub;
            subsegments.push_back({ia, ib, matched_frame, t0, t1, l_sub, nx, ny});
            constrained_edge_set.insert({std::min(ia, ib), std::max(ia, ib)});
        }
    }

    double xmin = verts[0].x, xmax = verts[0].x;
    double ymin = verts[0].y, ymax = verts[0].y;
    for (const auto& p : verts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }

    int nx_grid = std::max(2, static_cast<int>(std::ceil((xmax - xmin) / h_target)));
    int ny_grid = std::max(2, static_cast<int>(std::ceil((ymax - ymin) / h_target)));
    double dx_g = (xmax - xmin) / static_cast<double>(nx_grid);
    double dy_g = (ymax - ymin) / static_cast<double>(ny_grid);
    double clearance = 0.32 * std::min(h_target, std::min(dx_g, dy_g));

    auto try_insert_seed = [&](double px, double py) {
        if (static_cast<int>(verts.size()) >= max_vertices) {
            return;
        }
        if (px <= xmin + snap_eps || px >= xmax - snap_eps || py <= ymin + snap_eps || py >= ymax - snap_eps) {
            return;
        }
        if (!point_in_polygon({px, py}, outer_ccw)) {
            return;
        }
        for (const auto& hp : holes_ccw) {
            if (point_in_polygon({px, py}, hp)) {
                return;
            }
        }
        for (const auto& sub : subsegments) {
            Point2D va = verts[sub.u];
            Point2D vb = verts[sub.v];
            if (point_segment_distance(px, py, va.x, va.y, vb.x, vb.y).dist < clearance) {
                return;
            }
        }
        for (const auto& v : verts) {
            if (std::hypot(px - v.x, py - v.y) < clearance) {
                return;
            }
        }
        verts.push_back({px, py});
        acute_flags.push_back(false);
    };

    for (int j = 0; j < ny_grid; ++j) {
        double py = ymin + (static_cast<double>(j) + 0.5) * dy_g;
        for (int i = 0; i < nx_grid; ++i) {
            try_insert_seed(xmin + (static_cast<double>(i) + 0.5) * dx_g, py);
        }
    }
    for (int j = 1; j < ny_grid; ++j) {
        double py = ymin + static_cast<double>(j) * dy_g;
        for (int i = 1; i < nx_grid; ++i) {
            try_insert_seed(xmin + static_cast<double>(i) * dx_g, py);
        }
    }

    int n_total_verts = static_cast<int>(verts.size());
    double cx_mid = 0.5 * (xmin + xmax);
    double cy_mid = 0.5 * (ymin + ymax);
    double span = std::max({xmax - xmin, ymax - ymin, char_len});
    std::vector<Point2D> all_pts = verts;
    all_pts.push_back({cx_mid - 20.0 * span, cy_mid - 20.0 * span});
    all_pts.push_back({cx_mid + 20.0 * span, cy_mid - 20.0 * span});
    all_pts.push_back({cx_mid, cy_mid + 20.0 * span});

    std::vector<std::array<int, 3>> tris = {{n_total_verts, n_total_verts + 1, n_total_verts + 2}};

    for (int p_idx = 0; p_idx < n_total_verts; ++p_idx) {
        double px = all_pts[p_idx].x;
        double py = all_pts[p_idx].y;
        std::vector<int> bad_tris;
        for (std::size_t t_i = 0; t_i < tris.size(); ++t_i) {
            auto [a, b, c] = tris[t_i];
            if (incircle_ccw(
                    all_pts[a].x, all_pts[a].y,
                    all_pts[b].x, all_pts[b].y,
                    all_pts[c].x, all_pts[c].y,
                    px, py) > 1e-11 * (span * span)) {
                bad_tris.push_back(static_cast<int>(t_i));
            }
        }
        if (bad_tris.empty()) {
            for (std::size_t t_i = 0; t_i < tris.size(); ++t_i) {
                auto [a, b, c] = tris[t_i];
                if (orient2d(all_pts[a].x, all_pts[a].y, all_pts[b].x, all_pts[b].y, px, py) >= -1e-10 &&
                    orient2d(all_pts[b].x, all_pts[b].y, all_pts[c].x, all_pts[c].y, px, py) >= -1e-10 &&
                    orient2d(all_pts[c].x, all_pts[c].y, all_pts[a].x, all_pts[a].y, px, py) >= -1e-10) {
                    bad_tris.push_back(static_cast<int>(t_i));
                    break;
                }
            }
        }

        std::vector<bool> is_bad(tris.size(), false);
        std::vector<std::pair<int, int>> edge_order;
        std::map<std::pair<int, int>, int> edge_count;
        std::map<std::pair<int, int>, std::pair<int, int>> directed_edge;

        for (int t_i : bad_tris) {
            is_bad[t_i] = true;
            auto [a, b, c] = tris[t_i];
            std::array<std::pair<int, int>, 3> edges = {{{a, b}, {b, c}, {c, a}}};
            for (auto [u, v] : edges) {
                std::pair<int, int> k = {std::min(u, v), std::max(u, v)};
                if (edge_count.find(k) == edge_count.end()) {
                    edge_order.push_back(k);
                    edge_count[k] = 1;
                } else {
                    edge_count[k] += 1;
                }
                directed_edge[k] = {u, v};
            }
        }

        std::vector<std::array<int, 3>> next_tris;
        next_tris.reserve(tris.size() + edge_order.size());
        for (std::size_t idx = 0; idx < tris.size(); ++idx) {
            if (!is_bad[idx]) {
                next_tris.push_back(tris[idx]);
            }
        }
        for (const auto& k : edge_order) {
            if (edge_count[k] == 1) {
                auto [u, v] = directed_edge[k];
                double area2 = orient2d(all_pts[u].x, all_pts[u].y, all_pts[v].x, all_pts[v].y, px, py);
                if (area2 > 1e-14) {
                    next_tris.push_back({u, v, p_idx});
                } else if (area2 < -1e-14) {
                    next_tris.push_back({v, u, p_idx});
                }
            }
        }
        tris = std::move(next_tris);
    }

    // Remove super-triangle incidents
    {
        std::vector<std::array<int, 3>> filtered;
        for (const auto& t : tris) {
            if (t[0] < n_total_verts && t[1] < n_total_verts && t[2] < n_total_verts) {
                filtered.push_back(t);
            }
        }
        tris = std::move(filtered);
    }

    auto has_edge = [&](int u_t, int v_t) -> bool {
        std::pair<int, int> kt = {std::min(u_t, v_t), std::max(u_t, v_t)};
        for (auto [a, b, c] : tris) {
            if (std::make_pair(std::min(a, b), std::max(a, b)) == kt ||
                std::make_pair(std::min(b, c), std::max(b, c)) == kt ||
                std::make_pair(std::min(c, a), std::max(c, a)) == kt) {
                return true;
            }
        }
        return false;
    };

    for (auto [u_c, v_c] : constrained_edge_set) {
        if (has_edge(u_c, v_c)) {
            continue;
        }
        Point2D p_u = all_pts[u_c];
        Point2D p_v = all_pts[v_c];
        for (int iter = 0; iter < 50; ++iter) {
            if (has_edge(u_c, v_c)) {
                break;
            }
            bool flipped = false;
            std::size_t n_t = tris.size();
            for (std::size_t i = 0; i < n_t && !flipped; ++i) {
                auto [a, b, c] = tris[i];
                std::array<std::array<int, 3>, 3> cand = {{{a, b, c}, {b, c, a}, {c, a, b}}};
                for (auto [e_u, e_v, opp1] : cand) {
                    if (constrained_edge_set.count({std::min(e_u, e_v), std::max(e_u, e_v)}) > 0) {
                        continue;
                    }
                    if (e_u == u_c || e_u == v_c || e_v == u_c || e_v == v_c) {
                        continue;
                    }
                    if (!segments_properly_intersect(p_u, p_v, all_pts[e_u], all_pts[e_v], 1e-10).hit) {
                        continue;
                    }
                    for (std::size_t j = 0; j < n_t; ++j) {
                        if (i == j) {
                            continue;
                        }
                        auto ta = tris[j];
                        bool has_eu = (ta[0] == e_u || ta[1] == e_u || ta[2] == e_u);
                        bool has_ev = (ta[0] == e_v || ta[1] == e_v || ta[2] == e_v);
                        if (has_eu && has_ev) {
                            int opp2 = ta[0];
                            for (int x : ta) {
                                if (x != e_u && x != e_v) {
                                    opp2 = x;
                                    break;
                                }
                            }
                            if (orient2d(all_pts[opp1].x, all_pts[opp1].y, all_pts[e_u].x, all_pts[e_u].y, all_pts[opp2].x, all_pts[opp2].y) > 1e-13 &&
                                orient2d(all_pts[opp2].x, all_pts[opp2].y, all_pts[e_v].x, all_pts[e_v].y, all_pts[opp1].x, all_pts[opp1].y) > 1e-13) {
                                tris[i] = {opp1, e_u, opp2};
                                tris[j] = {opp2, e_v, opp1};
                                flipped = true;
                                break;
                            }
                        }
                    }
                    if (flipped) {
                        break;
                    }
                }
            }
            if (!flipped) {
                break;
            }
        }
        if (!has_edge(u_c, v_c)) {
            throw std::runtime_error("Constrained edge recovery failed in LadrunoMembrane::build_pslg_cdt_mesh.");
        }
    }

    std::vector<std::array<int, 3>> valid_tris;
    double tot_area = 0.0;
    for (auto [a, b, c] : tris) {
        double ax = all_pts[a].x, ay = all_pts[a].y;
        double bx = all_pts[b].x, by = all_pts[b].y;
        double cx = all_pts[c].x, cy = all_pts[c].y;
        double area2 = orient2d(ax, ay, bx, by, cx, cy);
        if (std::abs(area2) <= 1e-14) {
            continue;
        }
        if (area2 < 0.0) {
            std::swap(b, c);
            area2 = -area2;
        }
        double mx = (ax + bx + cx) / 3.0;
        double my = (ay + by + cy) / 3.0;
        if (!point_in_polygon({mx, my}, outer_ccw)) {
            continue;
        }
        bool in_hole = false;
        for (const auto& hp : holes_ccw) {
            if (point_in_polygon({mx, my}, hp)) {
                in_hole = true;
                break;
            }
        }
        if (in_hole) {
            continue;
        }
        valid_tris.push_back({a, b, c});
        tot_area += 0.5 * area2;
    }

    std::set<int> supported_verts(support_node_vertex_indices.begin(), support_node_vertex_indices.end());
    for (const auto& sub : subsegments) {
        if (sub.frame_id >= 0) {
            supported_verts.insert(sub.u);
            supported_verts.insert(sub.v);
        }
    }

    std::vector<std::array<int, 3>> refined_tris = valid_tris;
    for (std::size_t i_t = 0; i_t < refined_tris.size(); ++i_t) {
        auto [a, b, c] = refined_tris[i_t];
        if (supported_verts.count(a) && supported_verts.count(b) && supported_verts.count(c)) {
            bool flipped = false;
            std::array<std::array<int, 3>, 3> cand = {{{a, b, c}, {b, c, a}, {c, a, b}}};
            for (auto [eu, ev, opp1] : cand) {
                if (constrained_edge_set.count({std::min(eu, ev), std::max(eu, ev)}) > 0) {
                    continue;
                }
                for (std::size_t j_t = 0; j_t < refined_tris.size(); ++j_t) {
                    if (i_t == j_t) {
                        continue;
                    }
                    auto ta = refined_tris[j_t];
                    bool has_eu = (ta[0] == eu || ta[1] == eu || ta[2] == eu);
                    bool has_ev = (ta[0] == ev || ta[1] == ev || ta[2] == ev);
                    if (has_eu && has_ev) {
                        int opp2 = ta[0];
                        for (int x : ta) {
                            if (x != eu && x != ev) {
                                opp2 = x;
                                break;
                            }
                        }
                        if (orient2d(verts[opp1].x, verts[opp1].y, verts[eu].x, verts[eu].y, verts[opp2].x, verts[opp2].y) > 1e-13 &&
                            orient2d(verts[opp2].x, verts[opp2].y, verts[ev].x, verts[ev].y, verts[opp1].x, verts[opp1].y) > 1e-13) {
                            refined_tris[i_t] = {opp1, eu, opp2};
                            refined_tris[j_t] = {opp2, ev, opp1};
                            flipped = true;
                            break;
                        }
                    }
                }
                if (flipped) {
                    break;
                }
            }
        }
    }

    return PSLGMesh{
        verts,
        refined_tris,
        subsegments,
        acute_flags,
        boundary_vertex_indices,
        support_node_vertex_indices,
        tot_area,
        char_len
    };
}

// ============================================================================
// 3. Four Surface-to-Frame Load Distribution Engines
// ============================================================================

struct FrameSubLoad {
    double aL = 0.0;
    double bL = 1.0;
    double wa = 0.0;
    double wb = 0.0;
    double F = 0.0;
    double M_local = 0.0;
    int u = -1;
    int v = -1;
};

struct FrameLoadResult {
    int frame_id = 0;
    double L_frame = 0.0;
    double aL = 0.0;
    double bL = 1.0;
    double wa_eq = 0.0;
    double wb_eq = 0.0;
    double F_total = 0.0;
    double M_origin = 0.0;
    std::vector<FrameSubLoad> sub_loads;
    std::array<double, 3> p0_2d{};
    std::array<double, 3> q0_2d{};
    double overhang_torque = 0.0;
    double end_F_a = 0.0; // Ladruno WP-209
    double end_F_b = 0.0; // Ladruno WP-209
};

struct YieldLineSegment2D {
    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    int fidA = -1;
    int fidB = -1;
};

struct DistributionResult {
    std::string mode;
    double total_applied_load = 0.0;
    double total_applied_mx = 0.0;
    double total_applied_my = 0.0;
    std::vector<FrameLoadResult> frame_loads;
    std::vector<double> support_node_reactions;
    PSLGMesh mesh;
    int plate_support_rank = 2;
    std::vector<int> vertex_flags;
    std::vector<double> nodal_field;
    std::vector<int> triangle_support_id;
    std::vector<Point2D> triangle_flux_2d;
    std::vector<YieldLineSegment2D> yield_lines;
};

inline constexpr std::array<std::array<double, 3>, 3> TRI3_BARY = {{
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
    {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}
}};
inline constexpr std::array<double, 3> TRI3_W = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

inline constexpr std::array<std::array<double, 3>, 7> TRI7_BARY = {{
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
    {0.059715871789770, 0.470142064105115, 0.470142064105115},
    {0.470142064105115, 0.059715871789770, 0.470142064105115},
    {0.470142064105115, 0.470142064105115, 0.059715871789770},
    {0.797426985353087, 0.101286507323456, 0.101286507323456},
    {0.101286507323456, 0.797426985353087, 0.101286507323456},
    {0.101286507323456, 0.101286507323456, 0.797426985353087}
}};
inline constexpr std::array<double, 7> TRI7_W = {
    0.225000000000000,
    0.132394152788506,
    0.132394152788506,
    0.132394152788506,
    0.125939180544827,
    0.125939180544827,
    0.125939180544827
};

inline std::vector<FrameLoadResult> finalize_frame_results(
    const std::vector<Segment2D>& frames,
    const std::vector<std::vector<SubSegment>>& frame_sub_groups,
    const std::vector<std::vector<double>>& sub_F,
    const std::vector<std::vector<double>>& sub_M,
    const std::vector<double>& overhang_torques = {},
    const std::vector<double>& end_F_a = {},
    const std::vector<double>& end_F_b = {})
{
    std::vector<FrameLoadResult> results;
    results.reserve(frames.size());
    for (std::size_t k = 0; k < frames.size(); ++k) {
        const auto& fp = frames[k];
        double L_k = std::hypot(fp.p2.x - fp.p1.x, fp.p2.y - fp.p1.y);
        const auto& subs = frame_sub_groups[k];
        if (subs.empty() || L_k <= 1e-15) {
            FrameLoadResult empty_res{};
            empty_res.frame_id = static_cast<int>(k);
            empty_res.L_frame = L_k;
            empty_res.aL = 0.0;
            empty_res.bL = 1.0;
            results.push_back(empty_res);
            continue;
        }

        double F_tot = 0.0;
        double M_orig = 0.0;
        double aL_min = 1.0;
        double bL_max = 0.0;
        for (const auto& s : subs) {
            aL_min = std::min(aL_min, std::min(s.t_start, s.t_end));
            bL_max = std::max(bL_max, std::max(s.t_start, s.t_end));
        }

        std::vector<FrameSubLoad> sub_loads;
        std::array<double, 3> p0_sum = {0.0, 0.0, 0.0};
        std::array<double, 3> q0_sum = {0.0, 0.0, 0.0};

        for (std::size_t j = 0; j < subs.size(); ++j) {
            const auto& s = subs[j];
            double a_j = std::min(s.t_start, s.t_end);
            double b_j = std::max(s.t_start, s.t_end);
            double L_j = (b_j - a_j) * L_k;
            double f_j = sub_F[k][j];
            double m_j = sub_M[k][j];
#ifdef LADRUNO_MUTATE_MEMBRANE
            f_j = 0.0;
            m_j = 0.0;
#endif
            auto [wa_j, wb_j] = moment_conserving_trapezoid(f_j, m_j, L_j);
            int u_oriented = (s.t_start <= s.t_end) ? s.u : s.v;
            int v_oriented = (s.t_start <= s.t_end) ? s.v : s.u;
            sub_loads.push_back({a_j, b_j, wa_j, wb_j, f_j, m_j, u_oriented, v_oriented});
            F_tot += f_j;
            M_orig += f_j * (a_j * L_k) + m_j;
            auto fef = partial_trapezoid_fef_2d(L_k, wa_j, wb_j, 0.0, 0.0, a_j, b_j);
            for (int c = 0; c < 3; ++c) {
                p0_sum[c] += fef.p0[c];
                q0_sum[c] += fef.q0[c];
            }
        }

        double fa_end = (k < end_F_a.size()) ? end_F_a[k] : 0.0; // Ladruno WP-209
        double fb_end = (k < end_F_b.size()) ? end_F_b[k] : 0.0; // Ladruno WP-209
#ifdef LADRUNO_MUTATE_MEMBRANE
        fa_end = 0.0;
        fb_end = 0.0;
#endif
        F_tot += fa_end + fb_end;
        M_orig += fa_end * (aL_min * L_k) + fb_end * (bL_max * L_k);
        p0_sum[0] -= (fa_end + fb_end);
        p0_sum[1] -= fa_end * (1.0 - aL_min) + fb_end * (1.0 - bL_max);
        p0_sum[2] -= fa_end * aL_min + fb_end * bL_max;

        double L_active = (bL_max - aL_min) * L_k;
        double M_active = M_orig - F_tot * (aL_min * L_k);
        auto [wa_eq, wb_eq] = moment_conserving_trapezoid(F_tot, M_active, L_active);
        double torq = (k < overhang_torques.size()) ? overhang_torques[k] : 0.0;

        results.push_back(FrameLoadResult{
            static_cast<int>(k),
            L_k,
            aL_min,
            bL_max,
            wa_eq,
            wb_eq,
            F_tot,
            M_orig,
            std::move(sub_loads),
            p0_sum,
            q0_sum,
            torq,
            fa_end,
            fb_end
        });
    }
    return results;
}

inline void deposit_point_load_to_frame(
    int k,
    double t_frame,
    double L_k,
    double df,
    const std::vector<std::vector<SubSegment>>& frame_sub_groups,
    std::vector<std::vector<double>>& sub_F,
    std::vector<std::vector<double>>& sub_M)
{
    const auto& subs = frame_sub_groups[k];
    if (subs.empty()) {
        return;
    }
    std::size_t best_j = 0;
    double best_dist = 1e30;
    for (std::size_t j = 0; j < subs.size(); ++j) {
        double a_j = std::min(subs[j].t_start, subs[j].t_end);
        double b_j = std::max(subs[j].t_start, subs[j].t_end);
        if (t_frame >= a_j - 1e-12 && t_frame <= b_j + 1e-12) {
            best_j = j;
            best_dist = 0.0;
            break;
        }
        double d = std::min(std::abs(t_frame - a_j), std::abs(t_frame - b_j));
        if (d < best_dist) {
            best_dist = d;
            best_j = j;
        }
    }
    const auto& s_best = subs[best_j];
    double a_j = std::min(s_best.t_start, s_best.t_end);
    double b_j = std::max(s_best.t_start, s_best.t_end);
    double t_clamped = std::max(a_j, std::min(b_j, t_frame));
    double s_local = (t_clamped - a_j) * L_k;
    sub_F[k][best_j] += df;
    sub_M[k][best_j] += df * s_local;
}

inline void apportion_nodal_reactions_to_subsegments(
    const PSLGMesh& mesh,
    const std::vector<Segment2D>& frames,
    const std::vector<std::vector<SubSegment>>& frame_sub_groups,
    std::vector<std::vector<double>>& sub_F,
    std::vector<std::vector<double>>& sub_M,
    std::vector<double>& node_reactions,
    const std::vector<double>& nodal_R,
    const std::map<std::pair<int, int>, double>& sub_flux_weights)
{
    std::size_t n_verts = mesh.vertices.size();
    struct IncEntry {
        int k;
        int j;
        bool is_b;
        double w;
    };
    std::vector<std::vector<IncEntry>> vert_inc(n_verts);
    for (std::size_t k = 0; k < frame_sub_groups.size(); ++k) {
        for (std::size_t j = 0; j < frame_sub_groups[k].size(); ++j) {
            const auto& s = frame_sub_groups[k][j];
            auto it = sub_flux_weights.find({static_cast<int>(k), static_cast<int>(j)});
            double w_kj = (it != sub_flux_weights.end()) ? it->second : s.length;
            vert_inc[s.u].push_back({static_cast<int>(k), static_cast<int>(j), false, w_kj});
            vert_inc[s.v].push_back({static_cast<int>(k), static_cast<int>(j), true, w_kj});
        }
    }

    std::map<int, int> supp_vert_to_idx;
    for (std::size_t m = 0; m < mesh.support_node_vertex_indices.size(); ++m) {
        supp_vert_to_idx[mesh.support_node_vertex_indices[m]] = static_cast<int>(m);
    }

    std::vector<std::vector<double>> sub_Ra(frames.size());
    std::vector<std::vector<double>> sub_Rb(frames.size());
    for (std::size_t k = 0; k < frames.size(); ++k) {
        sub_Ra[k].assign(frame_sub_groups[k].size(), 0.0);
        sub_Rb[k].assign(frame_sub_groups[k].size(), 0.0);
    }

    for (std::size_t v_idx = 0; v_idx < n_verts; ++v_idx) {
        double r_val = nodal_R[v_idx];
        if (std::abs(r_val) <= 1e-18) {
            continue;
        }
        const auto& inc = vert_inc[v_idx];
        if (!inc.empty()) {
            double w_sum = 0.0;
            for (const auto& item : inc) {
                w_sum += item.w;
            }
            if (w_sum <= 1e-30) {
                w_sum = static_cast<double>(inc.size());
                for (const auto& item : inc) {
                    double share = r_val / w_sum;
                    if (item.is_b) {
                        sub_Rb[item.k][item.j] += share;
                    } else {
                        sub_Ra[item.k][item.j] += share;
                    }
                }
            } else {
                for (const auto& item : inc) {
                    double share = r_val * (item.w / w_sum);
                    if (item.is_b) {
                        sub_Rb[item.k][item.j] += share;
                    } else {
                        sub_Ra[item.k][item.j] += share;
                    }
                }
            }
        } else {
            auto it = supp_vert_to_idx.find(static_cast<int>(v_idx));
            if (it != supp_vert_to_idx.end()) {
                node_reactions[it->second] += r_val;
            }
        }
    }

    for (std::size_t k = 0; k < frame_sub_groups.size(); ++k) {
        for (std::size_t j = 0; j < frame_sub_groups[k].size(); ++j) {
            double ra = sub_Ra[k][j];
            double rb = sub_Rb[k][j];
            double L_j = frame_sub_groups[k][j].length;
            sub_F[k][j] = ra + rb;
            sub_M[k][j] = rb * L_j;
        }
    }
}

inline std::array<double, 81> dkt_element_stiffness(
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    const std::array<double, 9>& Db)
{
    double x23 = x2 - x3, y23 = y2 - y3;
    double x31 = x3 - x1, y31 = y3 - y1;
    double x12 = x1 - x2, y12 = y1 - y2;

    double area2 = x31 * y12 - x12 * y31;
    double area = 0.5 * area2;

    double l4_2 = x23 * x23 + y23 * y23;
    double l5_2 = x31 * x31 + y31 * y31;
    double l6_2 = x12 * x12 + y12 * y12;

    double a4 = -x23 / l4_2, a5 = -x31 / l5_2, a6 = -x12 / l6_2;
    double b4 = 0.75 * x23 * y23 / l4_2, b5 = 0.75 * x31 * y31 / l5_2, b6 = 0.75 * x12 * y12 / l6_2;
    double c4 = (-0.25 * x23 * x23 + 0.5 * y23 * y23) / l4_2;
    double c5 = (-0.25 * x31 * x31 + 0.5 * y31 * y31) / l5_2;
    double c6 = (-0.25 * x12 * x12 + 0.5 * y12 * y12) / l6_2;
    double d4 = -y23 / l4_2, d5 = -y31 / l5_2, d6 = -y12 / l6_2;
    double e4 = (0.5 * x23 * x23 - 0.25 * y23 * y23) / l4_2;
    double e5 = (0.5 * x31 * x31 - 0.25 * y31 * y31) / l5_2;
    double e6 = (0.5 * x12 * x12 - 0.25 * y12 * y12) / l6_2;

    std::array<std::array<double, 9>, 6> Tx{};
    std::array<std::array<double, 9>, 6> Ty{};

    Tx[0][1] = 1.0; Ty[0][2] = 1.0;
    Tx[1][4] = 1.0; Ty[1][5] = 1.0;
    Tx[2][7] = 1.0; Ty[2][8] = 1.0;

    Tx[3] = {0.0, 0.0, 0.0, -1.5 * a4, c4, -b4, 1.5 * a4, c4, -b4};
    Ty[3] = {0.0, 0.0, 0.0, -1.5 * d4, -b4, e4, 1.5 * d4, -b4, e4};

    Tx[4] = {1.5 * a5, c5, -b5, 0.0, 0.0, 0.0, -1.5 * a5, c5, -b5};
    Ty[4] = {1.5 * d5, -b5, e5, 0.0, 0.0, 0.0, -1.5 * d5, -b5, e5};

    Tx[5] = {-1.5 * a6, c6, -b6, 1.5 * a6, c6, -b6, 0.0, 0.0, 0.0};
    Ty[5] = {-1.5 * d6, -b6, e6, 1.5 * d6, -b6, e6, 0.0, 0.0, 0.0};

    std::array<double, 81> Ke{};
    for (int g = 0; g < 3; ++g) {
        double L1 = TRI3_BARY[g][0];
        double L2 = TRI3_BARY[g][1];
        double L3 = TRI3_BARY[g][2];
        double wg = TRI3_W[g];

        std::array<double, 6> dP_dL1 = {
            4.0 * L1 - 1.0,
            0.0,
            -(4.0 * L3 - 1.0),
            -4.0 * L2,
            4.0 * (L3 - L1),
            4.0 * L2
        };
        std::array<double, 6> dP_dL2 = {
            0.0,
            4.0 * L2 - 1.0,
            -(4.0 * L3 - 1.0),
            4.0 * (L3 - L2),
            -4.0 * L1,
            4.0 * L1
        };

        std::array<double, 9> dHx_dL1{}, dHx_dL2{}, dHy_dL1{}, dHy_dL2{};
        for (int k = 0; k < 6; ++k) {
            for (int j = 0; j < 9; ++j) {
                dHx_dL1[j] += dP_dL1[k] * Tx[k][j];
                dHx_dL2[j] += dP_dL2[k] * Tx[k][j];
                dHy_dL1[j] += dP_dL1[k] * Ty[k][j];
                dHy_dL2[j] += dP_dL2[k] * Ty[k][j];
            }
        }

        std::array<std::array<double, 9>, 3> B{};
        for (int j = 0; j < 9; ++j) {
            B[0][j] = (y23 * dHx_dL1[j] + y31 * dHx_dL2[j]) / area2;
            B[1][j] = (-x23 * dHy_dL1[j] - x31 * dHy_dL2[j]) / area2;
            B[2][j] = (-x23 * dHx_dL1[j] - x31 * dHx_dL2[j] + y23 * dHy_dL1[j] + y31 * dHy_dL2[j]) / area2;
        }

        double w_fac = area * wg;
        std::array<std::array<double, 9>, 3> DbB{};
        for (int r = 0; r < 3; ++r) {
            for (int j = 0; j < 9; ++j) {
                DbB[r][j] = Db[r * 3 + 0] * B[0][j] + Db[r * 3 + 1] * B[1][j] + Db[r * 3 + 2] * B[2][j];
            }
        }
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) {
                double s = B[0][i] * DbB[0][j] + B[1][i] * DbB[1][j] + B[2][i] * DbB[2][j];
                Ke[i * 9 + j] += w_fac * s;
            }
        }
    }
    return Ke;
}

inline DistributionResult distribute_surface_load(
    const std::vector<Point2D>& outer_poly,
    double q,
    const std::string& mode = "tributary",
    const std::vector<Segment2D>& frames = {},
    const std::vector<std::vector<Point2D>>& holes = {},
    const std::vector<Point2D>& support_nodes = {},
    Point2D span_dir = {1.0, 0.0},
    double orthotropy_eta = 1.0,
    double nu = 0.2,
    double h_target = 0.0)
{
    if (frames.empty() && support_nodes.empty()) {
        throw std::invalid_argument("At least one supporting frame or support_node is required.");
    }

    PSLGMesh mesh = build_pslg_cdt_mesh(outer_poly, holes, frames, support_nodes, h_target);

    std::size_t n_frames = frames.size();
    std::size_t n_supp_nodes = support_nodes.size();
    std::vector<std::vector<SubSegment>> frame_sub_groups(n_frames);
    for (const auto& sub : mesh.subsegments) {
        if (sub.frame_id >= 0 && static_cast<std::size_t>(sub.frame_id) < n_frames) {
            frame_sub_groups[sub.frame_id].push_back(sub);
        }
    }
    for (std::size_t k = 0; k < n_frames; ++k) {
        std::stable_sort(
            frame_sub_groups[k].begin(),
            frame_sub_groups[k].end(),
            [](const SubSegment& a, const SubSegment& b) {
                return std::min(a.t_start, a.t_end) < std::min(b.t_start, b.t_end);
            });
    }

    std::vector<std::vector<double>> sub_F(n_frames);
    std::vector<std::vector<double>> sub_M(n_frames);
    for (std::size_t k = 0; k < n_frames; ++k) {
        sub_F[k].assign(frame_sub_groups[k].size(), 0.0);
        sub_M[k].assign(frame_sub_groups[k].size(), 0.0);
    }
    std::vector<double> node_reactions(n_supp_nodes, 0.0);

    double total_load = 0.0;
    double total_mx = 0.0;
    double total_my = 0.0;
    const auto& verts = mesh.vertices;
    for (auto [a, b, c] : mesh.triangles) {
        double ax = verts[a].x, ay = verts[a].y;
        double bx = verts[b].x, by = verts[b].y;
        double cx = verts[c].x, cy = verts[c].y;
        double area_t = 0.5 * orient2d(ax, ay, bx, by, cx, cy);
        double f_t = q * area_t;
        total_load += f_t;
        total_mx += f_t * ((ay + by + cy) / 3.0);
        total_my -= f_t * ((ax + bx + cx) / 3.0);
    }

    struct ActiveFrameGeom {
        int k;
        double ax, ay, bx, by, a_min, b_max, L_k;
    };
    std::vector<ActiveFrameGeom> active_frame_geom;
    std::map<int, ActiveFrameGeom> geom_by_k;
    for (std::size_t k = 0; k < n_frames; ++k) {
        const auto& subs = frame_sub_groups[k];
        double L_k = std::hypot(frames[k].p2.x - frames[k].p1.x, frames[k].p2.y - frames[k].p1.y);
        if (subs.empty() || L_k <= 1e-15) {
            continue;
        }
        double a_min = 1.0, b_max = 0.0;
        for (const auto& s : subs) {
            a_min = std::min(a_min, std::min(s.t_start, s.t_end));
            b_max = std::max(b_max, std::max(s.t_start, s.t_end));
        }
        double ax = frames[k].p1.x + a_min * (frames[k].p2.x - frames[k].p1.x);
        double ay = frames[k].p1.y + a_min * (frames[k].p2.y - frames[k].p1.y);
        double bx = frames[k].p1.x + b_max * (frames[k].p2.x - frames[k].p1.x);
        double by = frames[k].p1.y + b_max * (frames[k].p2.y - frames[k].p1.y);
        ActiveFrameGeom afg{static_cast<int>(k), ax, ay, bx, by, a_min, b_max, L_k};
        active_frame_geom.push_back(afg);
        geom_by_k[static_cast<int>(k)] = afg;
    }

    std::string mode_clean;
    for (char ch : mode) {
        if (!std::isspace(static_cast<unsigned char>(ch))) {
            mode_clean.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
        }
    }

    // Build vertex role flags:
    // 0 = interior free Steiner node, 1 = outer polygon boundary, 2 = hole boundary,
    // 3 = constrained on frame segment, 4 = constrained on point support node.
    std::vector<int> vertex_flags(verts.size(), 0);
    double bnd_tol = 1e-6 * std::max(1.0, mesh.char_length);
    for (std::size_t v = 0; v < verts.size(); ++v) {
        double vx = verts[v].x, vy = verts[v].y;
        for (std::size_t i = 0; i < outer_poly.size(); ++i) {
            const auto& pA = outer_poly[i];
            const auto& pB = outer_poly[(i + 1) % outer_poly.size()];
            if (point_segment_distance(vx, vy, pA.x, pA.y, pB.x, pB.y).dist <= bnd_tol) {
                vertex_flags[v] = std::max(vertex_flags[v], 1);
                break;
            }
        }
        for (const auto& h_loop : holes) {
            for (std::size_t i = 0; i < h_loop.size(); ++i) {
                const auto& pA = h_loop[i];
                const auto& pB = h_loop[(i + 1) % h_loop.size()];
                if (point_segment_distance(vx, vy, pA.x, pA.y, pB.x, pB.y).dist <= bnd_tol) {
                    vertex_flags[v] = std::max(vertex_flags[v], 2);
                    break;
                }
            }
        }
    }
    for (const auto& sub : mesh.subsegments) {
        if (sub.frame_id >= 0) {
            vertex_flags[sub.u] = std::max(vertex_flags[sub.u], 3);
            vertex_flags[sub.v] = std::max(vertex_flags[sub.v], 3);
        }
    }
    for (int v_s : mesh.support_node_vertex_indices) {
        if (v_s >= 0 && static_cast<std::size_t>(v_s) < vertex_flags.size()) {
            vertex_flags[v_s] = 4;
        }
    }

    struct TiedFeat {
        int fid;
        double tf;
        bool in_sp;
        double cx;
        double cy;
    };

    auto nearest_features = [&](double px, double py) -> std::pair<std::vector<TiedFeat>, double> {
        double best_d = 1e30;
        struct Cand {
            int fid;
            double d;
            double tf;
            bool in_sp;
            double cx;
            double cy;
        };
        std::vector<Cand> candidates;
        for (const auto& afg : active_frame_geom) {
            double dx = afg.bx - afg.ax;
            double dy = afg.by - afg.ay;
            double l2 = dx * dx + dy * dy;
            double t_raw = (l2 > 1e-30) ? (((px - afg.ax) * dx + (py - afg.ay) * dy) / l2) : 0.0;
            bool in_span = (t_raw >= -1e-7 && t_raw <= 1.0 + 1e-7);
            double t_sub = std::max(0.0, std::min(1.0, t_raw));
            double cx = afg.ax + t_sub * dx;
            double cy = afg.ay + t_sub * dy;
            double d = std::hypot(px - cx, py - cy);
            double t_f = afg.a_min + t_sub * (afg.b_max - afg.a_min);
            candidates.push_back({afg.k, d, t_f, in_span, cx, cy});
            if (d < best_d) {
                best_d = d;
            }
        }
        for (std::size_t m = 0; m < n_supp_nodes; ++m) {
            double cx = support_nodes[m].x;
            double cy = support_nodes[m].y;
            double d = std::hypot(px - cx, py - cy);
            candidates.push_back({static_cast<int>(n_frames + m), d, 0.0, true, cx, cy});
            if (d < best_d) {
                best_d = d;
            }
        }
        std::vector<TiedFeat> tied;
        for (const auto& c : candidates) {
            if (c.d <= best_d + 1e-11) {
                tied.push_back({c.fid, c.tf, c.in_sp, c.cx, c.cy});
            }
        }
        return {tied, best_d};
    };

    auto feature_point_dist = [&](int fid, double px, double py) -> std::pair<double, bool> {
        if (fid >= 0 && fid < static_cast<int>(n_frames)) {
            auto it = geom_by_k.find(fid);
            if (it != geom_by_k.end()) {
                const auto& g = it->second;
                double dx = g.bx - g.ax;
                double dy = g.by - g.ay;
                double l2 = dx * dx + dy * dy;
                if (l2 <= 1e-30) {
                    return {std::hypot(px - g.ax, py - g.ay), false};
                }
                double t_raw = ((px - g.ax) * dx + (py - g.ay) * dy) / l2;
                bool in_span = (t_raw >= -1e-7 && t_raw <= 1.0 + 1e-7);
                double t_c = std::max(0.0, std::min(1.0, t_raw));
                return {std::hypot(px - (g.ax + t_c * dx), py - (g.ay + t_c * dy)), in_span};
            }
        } else if (fid >= static_cast<int>(n_frames)) {
            std::size_t m = static_cast<std::size_t>(fid - static_cast<int>(n_frames));
            if (m < n_supp_nodes) {
                return {std::hypot(px - support_nodes[m].x, py - support_nodes[m].y), true};
            }
        }
        return {1e30, false};
    };

    std::vector<YieldLineSegment2D> yield_lines;
    auto record_yield_segment = [&](Point2D pA, Point2D pB, int fidA, int fidB) {
        double seg_len = std::hypot(pB.x - pA.x, pB.y - pA.y);
        if (seg_len <= 1e-7 * std::max(1.0, mesh.char_length)) {
            return;
        }
        double mx = 0.5 * (pA.x + pB.x);
        double my = 0.5 * (pA.y + pB.y);
        auto [dA, _a] = feature_point_dist(fidA, mx, my);
        auto [dB, _b] = feature_point_dist(fidB, mx, my);
        auto [_, best_d] = nearest_features(mx, my);
        double tol_d = 2e-3 * std::max(1.0, mesh.char_length);
        if (std::abs(dA - dB) <= tol_d && dA <= best_d + tol_d) {
            int fMin = std::min(fidA, fidB);
            int fMax = std::max(fidA, fidB);
            yield_lines.push_back({pA.x, pA.y, pB.x, pB.y, fMin, fMax});
        }
    };

    auto extract_yield_lines_on_tri = [&](auto& self, Point2D p0, Point2D p1, Point2D p2, int depth) -> void {
        double area_t = 0.5 * std::abs(orient2d(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y));
        if (area_t <= 1e-16) {
            return;
        }
        std::array<Point2D, 3> corners = {p0, p1, p2};
        std::array<Point2D, 4> samples = {{
            {0.5 * p0.x + 0.25 * p1.x + 0.25 * p2.x, 0.5 * p0.y + 0.25 * p1.y + 0.25 * p2.y},
            {0.25 * p0.x + 0.5 * p1.x + 0.25 * p2.x, 0.25 * p0.y + 0.5 * p1.y + 0.25 * p2.y},
            {0.25 * p0.x + 0.25 * p1.x + 0.5 * p2.x, 0.25 * p0.y + 0.25 * p1.y + 0.5 * p2.y},
            {(p0.x + p1.x + p2.x) / 3.0, (p0.y + p1.y + p2.y) / 3.0}
        }};
        std::set<int> fids;
        for (const auto& s : samples) {
            auto [tied, _] = nearest_features(s.x, s.y);
            for (const auto& tf : tied) fids.insert(tf.fid);
        }
        for (const auto& s : corners) {
            auto [tied, _] = nearest_features(s.x, s.y);
            if (tied.size() == 1) fids.insert(tied[0].fid);
        }
        if (fids.size() <= 1) {
            return;
        }
        if (fids.size() >= 2 && depth < 4) {
            std::vector<int> sorted_fids(fids.begin(), fids.end());
            for (std::size_t ia = 0; ia < sorted_fids.size(); ++ia) {
                for (std::size_t ib = ia + 1; ib < sorted_fids.size(); ++ib) {
                    int fA = sorted_fids[ia];
                    int fB = sorted_fids[ib];
                    std::array<double, 3> g_vals{};
                    for (int i = 0; i < 3; ++i) {
                        auto [d1, _1] = feature_point_dist(fA, corners[i].x, corners[i].y);
                        auto [d2, _2] = feature_point_dist(fB, corners[i].x, corners[i].y);
                        g_vals[i] = d1 - d2;
                    }
                    std::vector<int> pos, neg;
                    for (int i = 0; i < 3; ++i) {
                        if (g_vals[i] > 1e-12) pos.push_back(i);
                        if (g_vals[i] < -1e-12) neg.push_back(i);
                    }
                    if ((pos.size() == 1 && neg.size() == 2) || (pos.size() == 2 && neg.size() == 1)) {
                        int solo = (pos.size() == 1) ? pos[0] : neg[0];
                        int i1 = (solo + 1) % 3;
                        int i2 = (solo + 2) % 3;
                        double t01 = g_vals[solo] / (g_vals[solo] - g_vals[i1]);
                        double t02 = g_vals[solo] / (g_vals[solo] - g_vals[i2]);
                        Point2D q01 = {
                            corners[solo].x + t01 * (corners[i1].x - corners[solo].x),
                            corners[solo].y + t01 * (corners[i1].y - corners[solo].y)
                        };
                        Point2D q02 = {
                            corners[solo].x + t02 * (corners[i2].x - corners[solo].x),
                            corners[solo].y + t02 * (corners[i2].y - corners[solo].y)
                        };
                        record_yield_segment(q01, q02, fA, fB);
                        self(self, corners[solo], q01, q02, depth + 1);
                        self(self, q01, corners[i1], corners[i2], depth + 1);
                        self(self, q01, corners[i2], q02, depth + 1);
                        return;
                    }
                    if (pos.size() == 1 && neg.size() == 1) {
                        int ip = pos[0];
                        int ineg = neg[0];
                        int iz = 3 - ip - ineg;
                        double tp = g_vals[ip] / (g_vals[ip] - g_vals[ineg]);
                        Point2D q_mid = {
                            corners[ip].x + tp * (corners[ineg].x - corners[ip].x),
                            corners[ip].y + tp * (corners[ineg].y - corners[ip].y)
                        };
                        record_yield_segment(corners[iz], q_mid, fA, fB);
                        self(self, corners[iz], corners[ip], q_mid, depth + 1);
                        self(self, corners[iz], q_mid, corners[ineg], depth + 1);
                        return;
                    }
                }
            }
        }
    };

    if (active_frame_geom.size() + n_supp_nodes >= 2) {
        for (auto [a, b, c] : mesh.triangles) {
            extract_yield_lines_on_tri(extract_yield_lines_on_tri, verts[a], verts[b], verts[c], 0);
        }
    }

    if (mode_clean == "tributary") {
        std::vector<double> end_F_a(n_frames, 0.0); // Ladruno WP-209
        std::vector<double> end_F_b(n_frames, 0.0); // Ladruno WP-209
        auto deposit_tributary_point = [&](double px, double py, double df) {
            auto [tied, _] = nearest_features(px, py);
            double share = df / static_cast<double>(tied.size());
            for (const auto& tf : tied) {
                if (tf.fid < static_cast<int>(n_frames)) {
                    if (!tf.in_sp) { // Ladruno WP-209
                        double a_min = 0.0;
                        double b_max = 1.0;
                        for (const auto& afg : active_frame_geom) {
                            if (afg.k == tf.fid) {
                                a_min = afg.a_min;
                                b_max = afg.b_max;
                                break;
                            }
                        }
                        if (tf.tf <= 0.5 * (a_min + b_max)) {
                            end_F_a[tf.fid] += share;
                        } else {
                            end_F_b[tf.fid] += share;
                        }
                    } else {
                        double L_k = std::hypot(
                            frames[tf.fid].p2.x - frames[tf.fid].p1.x,
                            frames[tf.fid].p2.y - frames[tf.fid].p1.y);
                        deposit_point_load_to_frame(tf.fid, tf.tf, L_k, share, frame_sub_groups, sub_F, sub_M);
                    }
                } else {
                    node_reactions[tf.fid - static_cast<int>(n_frames)] += share;
                }
            }
        };

        auto integrate_single_feat_tri = [&](Point2D p0, Point2D p1, Point2D p2) {
            double area_t = 0.5 * std::abs(orient2d(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y));
            if (area_t <= 1e-16) {
                return;
            }
            for (int g = 0; g < 3; ++g) {
                double b0 = TRI3_BARY[g][0], b1 = TRI3_BARY[g][1], b2 = TRI3_BARY[g][2];
                double qx = b0 * p0.x + b1 * p1.x + b2 * p2.x;
                double qy = b0 * p0.y + b1 * p1.y + b2 * p2.y;
                deposit_tributary_point(qx, qy, q * area_t * TRI3_W[g]);
            }
        };

        auto line_signed_dist = [&](double px, double py, double ax, double ay, double bx, double by) -> std::pair<double, bool> {
            double dx = bx - ax;
            double dy = by - ay;
            double l2 = dx * dx + dy * dy;
            if (l2 <= 1e-30) {
                return {std::hypot(px - ax, py - ay), false};
            }
            double t_raw = ((px - ax) * dx + (py - ay) * dy) / l2;
            bool in_span = (t_raw >= -1e-7 && t_raw <= 1.0 + 1e-7);
            double t_c = std::max(0.0, std::min(1.0, t_raw));
            return {std::hypot(px - (ax + t_c * dx), py - (ay + t_c * dy)), in_span};
        };

        auto integrate_tributary_tri = [&](auto& self, Point2D p0, Point2D p1, Point2D p2, int depth) -> void {
            double area_t = 0.5 * std::abs(orient2d(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y));
            if (area_t <= 1e-16) {
                return;
            }
            std::array<Point2D, 3> corners = {p0, p1, p2};
            std::array<Point2D, 4> int_samples = {{
                {0.5 * p0.x + 0.25 * p1.x + 0.25 * p2.x, 0.5 * p0.y + 0.25 * p1.y + 0.25 * p2.y},
                {0.25 * p0.x + 0.5 * p1.x + 0.25 * p2.x, 0.25 * p0.y + 0.5 * p1.y + 0.25 * p2.y},
                {0.25 * p0.x + 0.25 * p1.x + 0.5 * p2.x, 0.25 * p0.y + 0.25 * p1.y + 0.5 * p2.y},
                {(p0.x + p1.x + p2.x) / 3.0, (p0.y + p1.y + p2.y) / 3.0}
            }};
            std::set<int> fids;
            bool all_int = true;
            for (const auto& s : int_samples) {
                auto [tied, _] = nearest_features(s.x, s.y);
                for (const auto& tf : tied) {
                    fids.insert(tf.fid);
                    if (!tf.in_sp) {
                        all_int = false;
                    }
                }
            }
            for (const auto& s : corners) {
                auto [tied, _] = nearest_features(s.x, s.y);
                if (tied.size() == 1) {
                    fids.insert(tied[0].fid);
                    if (!tied[0].in_sp) {
                        all_int = false;
                    }
                }
            }

            if (fids.size() == 1 && all_int) {
                integrate_single_feat_tri(p0, p1, p2);
                return;
            }

            bool all_frames = true;
            for (int fid : fids) {
                if (fid >= static_cast<int>(n_frames)) {
                    all_frames = false;
                }
            }

            if (fids.size() >= 2 && all_frames && depth < 6) {
                std::vector<int> sorted_fids(fids.begin(), fids.end());
                for (std::size_t idx_a = 0; idx_a < sorted_fids.size(); ++idx_a) {
                    for (std::size_t idx_b = idx_a + 1; idx_b < sorted_fids.size(); ++idx_b) {
                        const auto& g1 = geom_by_k[sorted_fids[idx_a]];
                        const auto& g2 = geom_by_k[sorted_fids[idx_b]];
                        std::array<double, 3> g_vals{};
                        bool ok_int = true;
                        for (int i = 0; i < 3; ++i) {
                            auto [d1, int1] = line_signed_dist(corners[i].x, corners[i].y, g1.ax, g1.ay, g1.bx, g1.by);
                            auto [d2, int2] = line_signed_dist(corners[i].x, corners[i].y, g2.ax, g2.ay, g2.bx, g2.by);
                            if (!(int1 && int2)) {
                                ok_int = false;
                            }
                            g_vals[i] = d1 - d2;
                        }
                        if (ok_int) {
                            std::vector<int> pos, neg;
                            for (int i = 0; i < 3; ++i) {
                                if (g_vals[i] > 1e-12) pos.push_back(i);
                                if (g_vals[i] < -1e-12) neg.push_back(i);
                            }
                            if ((pos.size() == 1 && neg.size() == 2) || (pos.size() == 2 && neg.size() == 1)) {
                                int solo = (pos.size() == 1) ? pos[0] : neg[0];
                                int i1 = (solo + 1) % 3;
                                int i2 = (solo + 2) % 3;
                                double t01 = g_vals[solo] / (g_vals[solo] - g_vals[i1]);
                                double t02 = g_vals[solo] / (g_vals[solo] - g_vals[i2]);
                                Point2D q01 = {
                                    corners[solo].x + t01 * (corners[i1].x - corners[solo].x),
                                    corners[solo].y + t01 * (corners[i1].y - corners[solo].y)
                                };
                                Point2D q02 = {
                                    corners[solo].x + t02 * (corners[i2].x - corners[solo].x),
                                    corners[solo].y + t02 * (corners[i2].y - corners[solo].y)
                                };
                                self(self, corners[solo], q01, q02, depth + 1);
                                self(self, q01, corners[i1], corners[i2], depth + 1);
                                self(self, q01, corners[i2], q02, depth + 1);
                                return;
                            }
                            if (pos.size() == 1 && neg.size() == 1) {
                                int ip = pos[0];
                                int ineg = neg[0];
                                int iz = 3 - ip - ineg;
                                double tp = g_vals[ip] / (g_vals[ip] - g_vals[ineg]);
                                Point2D q_mid = {
                                    corners[ip].x + tp * (corners[ineg].x - corners[ip].x),
                                    corners[ip].y + tp * (corners[ineg].y - corners[ip].y)
                                };
                                self(self, corners[iz], corners[ip], q_mid, depth + 1);
                                self(self, corners[iz], q_mid, corners[ineg], depth + 1);
                                return;
                            }
                        }
                    }
                }
            }

            if (depth < 5) {
                Point2D m01 = {0.5 * (p0.x + p1.x), 0.5 * (p0.y + p1.y)};
                Point2D m12 = {0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y)};
                Point2D m20 = {0.5 * (p2.x + p0.x), 0.5 * (p2.y + p0.y)};
                self(self, p0, m01, m20, depth + 1);
                self(self, m01, p1, m12, depth + 1);
                self(self, m20, m12, p2, depth + 1);
                self(self, m01, m12, m20, depth + 1);
                return;
            }

            for (int g = 0; g < 7; ++g) {
                double b0 = TRI7_BARY[g][0], b1 = TRI7_BARY[g][1], b2 = TRI7_BARY[g][2];
                double qx = b0 * p0.x + b1 * p1.x + b2 * p2.x;
                double qy = b0 * p0.y + b1 * p1.y + b2 * p2.y;
                deposit_tributary_point(qx, qy, q * area_t * TRI7_W[g]);
            }
        };

        for (auto [a, b, c] : mesh.triangles) {
            integrate_tributary_tri(integrate_tributary_tri, verts[a], verts[b], verts[c], 0);
        }

        std::vector<double> nodal_field(verts.size(), 0.0);
        for (std::size_t v = 0; v < verts.size(); ++v) {
            auto [_, best_d] = nearest_features(verts[v].x, verts[v].y);
            nodal_field[v] = (best_d < 1e20) ? (q * best_d) : 0.0;
        }

        std::vector<int> tri_supp(mesh.triangles.size(), 0);
        std::vector<Point2D> tri_flux(mesh.triangles.size(), {0.0, 0.0});
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
            auto [a, b, c] = mesh.triangles[t];
            double cx = (verts[a].x + verts[b].x + verts[c].x) / 3.0;
            double cy = (verts[a].y + verts[b].y + verts[c].y) / 3.0;
            auto [tied, _] = nearest_features(cx, cy);
            if (!tied.empty()) {
                tri_supp[t] = tied[0].fid;
                double vx = tied[0].cx - cx;
                double vy = tied[0].cy - cy;
                double vn = std::hypot(vx, vy);
                if (vn > 1e-15) {
                    tri_flux[t] = {std::abs(q) * (vx / vn), std::abs(q) * (vy / vn)};
                }
            }
        }

        auto frame_results = finalize_frame_results(
            frames, frame_sub_groups, sub_F, sub_M, {}, end_F_a, end_F_b); // Ladruno WP-209
        return DistributionResult{
            mode_clean,
            total_load,
            total_mx,
            total_my,
            std::move(frame_results),
            std::move(node_reactions),
            std::move(mesh),
            2,
            std::move(vertex_flags),
            std::move(nodal_field),
            std::move(tri_supp),
            std::move(tri_flux),
            std::move(yield_lines)
        };
    }

    if (mode_clean == "oneway") {
        double sx = span_dir.x, sy = span_dir.y;
        double s_norm = std::hypot(sx, sy);
        if (s_norm <= 1e-15) {
            throw std::invalid_argument("span_dir must be non-zero.");
        }
        sx /= s_norm;
        sy /= s_norm;

        struct RayHit {
            bool valid = false;
            int k = -1;
            double dist = 1e30;
            double t_f = 0.0;
        };

        auto ray_hits = [&](double px, double py) -> std::pair<RayHit, RayHit> {
            RayHit hit_neg{}, hit_pos{};
            for (const auto& afg : active_frame_geom) {
                double dx = afg.bx - afg.ax;
                double dy = afg.by - afg.ay;
                double denom = sx * dy - sy * dx;
                if (std::abs(denom) <= 1e-12 * std::hypot(dx, dy)) {
                    continue;
                }
                double rx = afg.ax - px;
                double ry = afg.ay - py;
                double r_dist = (rx * dy - ry * dx) / denom;
                double u_seg = (rx * sy - ry * sx) / denom;
                if (u_seg >= -1e-9 && u_seg <= 1.0 + 1e-9) {
                    double u_c = std::max(0.0, std::min(1.0, u_seg));
                    double t_f = afg.a_min + u_c * (afg.b_max - afg.a_min);
                    if (r_dist >= -1e-11) {
                        double d_pos = std::max(0.0, r_dist);
                        if (!hit_pos.valid || d_pos < hit_pos.dist) {
                            hit_pos = {true, afg.k, d_pos, t_f};
                        }
                    }
                    if (r_dist <= 1e-11) {
                        double d_neg = std::max(0.0, -r_dist);
                        if (!hit_neg.valid || d_neg < hit_neg.dist) {
                            hit_neg = {true, afg.k, d_neg, t_f};
                        }
                    }
                }
            }
            return {hit_neg, hit_pos};
        };

        std::vector<int> tri_supp(mesh.triangles.size(), 0);
        std::vector<Point2D> tri_flux(mesh.triangles.size(), {0.0, 0.0});

        for (std::size_t t_idx = 0; t_idx < mesh.triangles.size(); ++t_idx) {
            auto [a, b, c] = mesh.triangles[t_idx];
            Point2D p0 = verts[a], p1 = verts[b], p2 = verts[c];
            double area_t = 0.5 * std::abs(orient2d(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y));
            double cx = (p0.x + p1.x + p2.x) / 3.0;
            double cy = (p0.y + p1.y + p2.y) / 3.0;
            auto [hc_neg, hc_pos] = ray_hits(cx, cy);
            if (hc_neg.valid && hc_pos.valid) {
                if (hc_neg.dist <= hc_pos.dist) {
                    tri_supp[t_idx] = hc_neg.k;
                    tri_flux[t_idx] = {-std::abs(q) * sx, -std::abs(q) * sy};
                } else {
                    tri_supp[t_idx] = hc_pos.k;
                    tri_flux[t_idx] = {std::abs(q) * sx, std::abs(q) * sy};
                }
            } else if (hc_neg.valid) {
                tri_supp[t_idx] = hc_neg.k;
                tri_flux[t_idx] = {-std::abs(q) * sx, -std::abs(q) * sy};
            } else if (hc_pos.valid) {
                tri_supp[t_idx] = hc_pos.k;
                tri_flux[t_idx] = {std::abs(q) * sx, std::abs(q) * sy};
            } else {
                auto [tied, _] = nearest_features(cx, cy);
                if (!tied.empty()) {
                    tri_supp[t_idx] = tied[0].fid;
                }
            }

            for (int g = 0; g < 7; ++g) {
                double b0 = TRI7_BARY[g][0], b1 = TRI7_BARY[g][1], b2 = TRI7_BARY[g][2];
                double qx = b0 * p0.x + b1 * p1.x + b2 * p2.x;
                double qy = b0 * p0.y + b1 * p1.y + b2 * p2.y;
                double df = q * area_t * TRI7_W[g];
                auto [h_neg, h_pos] = ray_hits(qx, qy);
                if (h_neg.valid && h_pos.valid) {
                    double d_tot = h_neg.dist + h_pos.dist;
                    double w_l = (d_tot <= 1e-13) ? 0.5 : (h_pos.dist / d_tot);
                    double w_r = (d_tot <= 1e-13) ? 0.5 : (h_neg.dist / d_tot);
                    double L_l = std::hypot(frames[h_neg.k].p2.x - frames[h_neg.k].p1.x, frames[h_neg.k].p2.y - frames[h_neg.k].p1.y);
                    double L_r = std::hypot(frames[h_pos.k].p2.x - frames[h_pos.k].p1.x, frames[h_pos.k].p2.y - frames[h_pos.k].p1.y);
                    deposit_point_load_to_frame(h_neg.k, h_neg.t_f, L_l, df * w_l, frame_sub_groups, sub_F, sub_M);
                    deposit_point_load_to_frame(h_pos.k, h_pos.t_f, L_r, df * w_r, frame_sub_groups, sub_F, sub_M);
                } else if (h_neg.valid) {
                    double L_l = std::hypot(frames[h_neg.k].p2.x - frames[h_neg.k].p1.x, frames[h_neg.k].p2.y - frames[h_neg.k].p1.y);
                    deposit_point_load_to_frame(h_neg.k, h_neg.t_f, L_l, df, frame_sub_groups, sub_F, sub_M);
                } else if (h_pos.valid) {
                    double L_r = std::hypot(frames[h_pos.k].p2.x - frames[h_pos.k].p1.x, frames[h_pos.k].p2.y - frames[h_pos.k].p1.y);
                    deposit_point_load_to_frame(h_pos.k, h_pos.t_f, L_r, df, frame_sub_groups, sub_F, sub_M);
                } else {
                    int best_k = -1;
                    double best_d = 1e30;
                    double best_t = 0.0;
                    for (const auto& afg : active_frame_geom) {
                        auto psd = point_segment_distance(qx, qy, afg.ax, afg.ay, afg.bx, afg.by);
                        if (psd.dist < best_d) {
                            best_d = psd.dist;
                            best_k = afg.k;
                            best_t = afg.a_min + psd.t * (afg.b_max - afg.a_min);
                        }
                    }
                    if (best_k >= 0) {
                        double L_b = std::hypot(frames[best_k].p2.x - frames[best_k].p1.x, frames[best_k].p2.y - frames[best_k].p1.y);
                        deposit_point_load_to_frame(best_k, best_t, L_b, df, frame_sub_groups, sub_F, sub_M);
                    } else if (n_supp_nodes > 0) {
                        std::size_t best_m = 0;
                        double min_d = 1e30;
                        for (std::size_t m = 0; m < n_supp_nodes; ++m) {
                            double d = std::hypot(qx - support_nodes[m].x, qy - support_nodes[m].y);
                            if (d < min_d) {
                                min_d = d;
                                best_m = m;
                            }
                        }
                        node_reactions[best_m] += df;
                    }
                }
            }
        }

        std::vector<double> nodal_field(verts.size(), 0.0);
        for (std::size_t v = 0; v < verts.size(); ++v) {
            auto [h_neg, h_pos] = ray_hits(verts[v].x, verts[v].y);
            if (h_neg.valid && h_pos.valid) {
                nodal_field[v] = 0.5 * q * h_neg.dist * h_pos.dist;
            } else if (h_neg.valid) {
                nodal_field[v] = 0.5 * q * h_neg.dist * h_neg.dist;
            } else if (h_pos.valid) {
                nodal_field[v] = 0.5 * q * h_pos.dist * h_pos.dist;
            }
        }

        auto frame_results = finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M);
        return DistributionResult{
            mode_clean,
            total_load,
            total_mx,
            total_my,
            std::move(frame_results),
            std::move(node_reactions),
            std::move(mesh),
            2,
            std::move(vertex_flags),
            std::move(nodal_field),
            std::move(tri_supp),
            std::move(tri_flux),
            std::move(yield_lines)
        };
    }

    if (mode_clean == "fea_poisson") {
        std::size_t n_v = verts.size();
        double sx = span_dir.x, sy = span_dir.y;
        double sn = std::hypot(sx, sy);
        if (sn <= 1e-15) {
            sx = 1.0; sy = 0.0;
        } else {
            sx /= sn; sy /= sn;
        }
        double eta = std::max(1e-8, orthotropy_eta);
        double d00 = sx * sx + eta * sy * sy;
        double d01 = (1.0 - eta) * sx * sy;
        double d11 = sy * sy + eta * sx * sx;

        std::set<int> constrained_set(mesh.support_node_vertex_indices.begin(), mesh.support_node_vertex_indices.end());
        for (const auto& subs : frame_sub_groups) {
            for (const auto& s : subs) {
                constrained_set.insert(s.u);
                constrained_set.insert(s.v);
            }
        }

        std::vector<double> K(n_v * n_v, 0.0);
        std::vector<double> f_glob(n_v, 0.0);
        std::vector<std::array<double, 6>> tri_B(mesh.triangles.size());

        for (std::size_t t_i = 0; t_i < mesh.triangles.size(); ++t_i) {
            auto [a, b, c] = mesh.triangles[t_i];
            double x0 = verts[a].x, y0 = verts[a].y;
            double x1 = verts[b].x, y1 = verts[b].y;
            double x2 = verts[c].x, y2 = verts[c].y;
            double area2 = orient2d(x0, y0, x1, y1, x2, y2);
            double area = 0.5 * area2;
            std::array<double, 6> B = {
                (y1 - y2) / area2, (y2 - y0) / area2, (y0 - y1) / area2,
                (x2 - x1) / area2, (x0 - x2) / area2, (x1 - x0) / area2
            };
            tri_B[t_i] = B;

            std::array<int, 3> idx = {a, b, c};
            std::vector<int> free_in_tri;
            for (int v_i : idx) {
                if (constrained_set.find(v_i) == constrained_set.end()) {
                    free_in_tri.push_back(v_i);
                }
            }
            if (!free_in_tri.empty()) {
                double f_share = (q * area) / static_cast<double>(free_in_tri.size());
                for (int v_i : free_in_tri) {
                    f_glob[v_i] += f_share;
                }
            } else {
                for (int v_i : idx) {
                    f_glob[v_i] += q * area / 3.0;
                }
            }

            for (int i_loc = 0; i_loc < 3; ++i_loc) {
                double bx_i = B[i_loc];
                double by_i = B[3 + i_loc];
                double dbx_i = d00 * bx_i + d01 * by_i;
                double dby_i = d01 * bx_i + d11 * by_i;
                for (int j_loc = 0; j_loc < 3; ++j_loc) {
                    double ke_ij = area * (dbx_i * B[j_loc] + dby_i * B[3 + j_loc]);
                    K[idx[i_loc] * n_v + idx[j_loc]] += ke_ij;
                }
            }
        }

        std::vector<int> free_dofs;
        for (std::size_t i = 0; i < n_v; ++i) {
            if (constrained_set.find(static_cast<int>(i)) == constrained_set.end()) {
                free_dofs.push_back(static_cast<int>(i));
            }
        }

        std::vector<double> w(n_v, 0.0);
        if (!free_dofs.empty()) {
            int nf = static_cast<int>(free_dofs.size());
            std::vector<double> K_ff(nf * nf, 0.0);
            std::vector<double> f_f(nf, 0.0);
            for (int r = 0; r < nf; ++r) {
                f_f[r] = f_glob[free_dofs[r]];
                for (int c = 0; c < nf; ++c) {
                    K_ff[r * nf + c] = K[free_dofs[r] * n_v + free_dofs[c]];
                }
            }
            auto w_free = solve_linear_system(nf, K_ff, f_f, 1);
            for (int r = 0; r < nf; ++r) {
                w[free_dofs[r]] = w_free[r];
            }
        }

        std::vector<double> nodal_R(n_v, 0.0);
        for (std::size_t i = 0; i < n_v; ++i) {
            double kw = 0.0;
            for (std::size_t j = 0; j < n_v; ++j) {
                kw += K[i * n_v + j] * w[j];
            }
            nodal_R[i] = f_glob[i] - kw;
        }

        std::map<std::pair<int, int>, std::vector<int>> edge_to_tris;
        for (std::size_t t_idx = 0; t_idx < mesh.triangles.size(); ++t_idx) {
            auto [a, b, c] = mesh.triangles[t_idx];
            edge_to_tris[{std::min(a, b), std::max(a, b)}].push_back(static_cast<int>(t_idx));
            edge_to_tris[{std::min(b, c), std::max(b, c)}].push_back(static_cast<int>(t_idx));
            edge_to_tris[{std::min(c, a), std::max(c, a)}].push_back(static_cast<int>(t_idx));
        }

        std::map<std::pair<int, int>, double> sub_flux_weights;
        for (std::size_t k = 0; k < frame_sub_groups.size(); ++k) {
            for (std::size_t j = 0; j < frame_sub_groups[k].size(); ++j) {
                const auto& s = frame_sub_groups[k][j];
                std::pair<int, int> key = {std::min(s.u, s.v), std::max(s.u, s.v)};
                double flux_mag = 0.0;
                double mx_s = 0.5 * (verts[s.u].x + verts[s.v].x);
                double my_s = 0.5 * (verts[s.u].y + verts[s.v].y);
                auto it = edge_to_tris.find(key);
                if (it != edge_to_tris.end()) {
                    for (int t_idx : it->second) {
                        auto [a, b, c] = mesh.triangles[t_idx];
                        const auto& B = tri_B[t_idx];
                        double gwx = B[0] * w[a] + B[1] * w[b] + B[2] * w[c];
                        double gwy = B[3] * w[a] + B[4] * w[b] + B[5] * w[c];
                        double qfx = -(d00 * gwx + d01 * gwy);
                        double qfy = -(d01 * gwx + d11 * gwy);
                        double tx = (verts[a].x + verts[b].x + verts[c].x) / 3.0;
                        double ty = (verts[a].y + verts[b].y + verts[c].y) / 3.0;
                        double nx = s.nx, ny = s.ny;
                        if ((mx_s - tx) * nx + (my_s - ty) * ny < 0.0) {
                            nx = -nx;
                            ny = -ny;
                        }
                        double fn = qfx * nx + qfy * ny;
                        if (std::abs(q) > 1e-15) {
                            double fn_signed = (q > 0.0) ? fn : -fn;
                            flux_mag += std::max(0.0, fn_signed);
                        }
                    }
                }
                double d_nn = s.nx * (d00 * s.nx + d01 * s.ny) + s.ny * (d01 * s.nx + d11 * s.ny);
                sub_flux_weights[{static_cast<int>(k), static_cast<int>(j)}] = s.length * (flux_mag + 1e-12 * d_nn);
            }
        }

        apportion_nodal_reactions_to_subsegments(
            mesh, frames, frame_sub_groups, sub_F, sub_M, node_reactions, nodal_R, sub_flux_weights);

        std::vector<int> tri_supp(mesh.triangles.size(), 0);
        std::vector<Point2D> tri_flux(mesh.triangles.size(), {0.0, 0.0});
        for (std::size_t t_idx = 0; t_idx < mesh.triangles.size(); ++t_idx) {
            auto [a, b, c] = mesh.triangles[t_idx];
            const auto& B = tri_B[t_idx];
            double gwx = B[0] * w[a] + B[1] * w[b] + B[2] * w[c];
            double gwy = B[3] * w[a] + B[4] * w[b] + B[5] * w[c];
            double qfx = -(d00 * gwx + d01 * gwy);
            double qfy = -(d01 * gwx + d11 * gwy);
            tri_flux[t_idx] = {qfx, qfy};
            double cx = (verts[a].x + verts[b].x + verts[c].x) / 3.0;
            double cy = (verts[a].y + verts[b].y + verts[c].y) / 3.0;
            auto [tied, _] = nearest_features(cx, cy);
            if (!tied.empty()) {
                tri_supp[t_idx] = tied[0].fid;
            }
        }

        auto frame_results = finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M);
        return DistributionResult{
            "fea_poisson",
            total_load,
            total_mx,
            total_my,
            std::move(frame_results),
            std::move(node_reactions),
            std::move(mesh),
            2,
            std::move(vertex_flags),
            std::move(w),
            std::move(tri_supp),
            std::move(tri_flux),
            std::move(yield_lines)
        };
    }

    if (mode_clean == "fea_plate") {
        std::size_t n_v = verts.size();
        std::size_t ndof = 3 * n_v;
        double fac = 1.0 / (1.0 - nu * nu);
        std::array<double, 9> Db = {
            fac, fac * nu, 0.0,
            fac * nu, fac, 0.0,
            0.0, 0.0, fac * 0.5 * (1.0 - nu)
        };

        std::vector<double> K(ndof * ndof, 0.0);
        std::vector<double> f_glob(ndof, 0.0);

        for (auto [a, b, c] : mesh.triangles) {
            double x1 = verts[a].x, y1 = verts[a].y;
            double x2 = verts[b].x, y2 = verts[b].y;
            double x3 = verts[c].x, y3 = verts[c].y;
            double area = 0.5 * orient2d(x1, y1, x2, y2, x3, y3);
            auto Ke = dkt_element_stiffness(x1, y1, x2, y2, x3, y3, Db);
            std::array<int, 9> dofs = {
                3 * a, 3 * a + 1, 3 * a + 2,
                3 * b, 3 * b + 1, 3 * b + 2,
                3 * c, 3 * c + 1, 3 * c + 2
            };
            f_glob[3 * a] += q * area / 3.0;
            f_glob[3 * b] += q * area / 3.0;
            f_glob[3 * c] += q * area / 3.0;
            for (int i_l = 0; i_l < 9; ++i_l) {
                int di = dofs[i_l];
                for (int j_l = 0; j_l < 9; ++j_l) {
                    K[di * ndof + dofs[j_l]] += Ke[i_l * 9 + j_l];
                }
            }
        }

        std::set<int> supp_vert_set(mesh.support_node_vertex_indices.begin(), mesh.support_node_vertex_indices.end());
        for (const auto& subs : frame_sub_groups) {
            for (const auto& s : subs) {
                supp_vert_set.insert(s.u);
                supp_vert_set.insert(s.v);
            }
        }
        std::vector<int> supp_verts(supp_vert_set.begin(), supp_vert_set.end());

        double mx_c = 0.0, my_c = 0.0;
        for (int v_i : supp_verts) {
            mx_c += verts[v_i].x;
            my_c += verts[v_i].y;
        }
        mx_c /= static_cast<double>(supp_verts.size());
        my_c /= static_cast<double>(supp_verts.size());
        double c00 = 0.0, c01 = 0.0, c11 = 0.0;
        for (int v_i : supp_verts) {
            double dx = verts[v_i].x - mx_c;
            double dy = verts[v_i].y - my_c;
            c00 += dx * dx;
            c01 += dx * dy;
            c11 += dy * dy;
        }
        auto [ev_min, ev_max] = sym_2x2_eigenvalues(c00, c01, c11);
        double sigma_max = std::sqrt(std::max(0.0, ev_max));
        double sigma_min = std::sqrt(std::max(0.0, ev_min));
        int support_rank = (supp_verts.size() >= 3 && sigma_min >= 1e-5 * std::max(sigma_max, mesh.char_length)) ? 2 : 1;

        std::set<int> constrained_dofs;
        for (int v_i : supp_verts) {
            constrained_dofs.insert(3 * v_i);
            if (support_rank < 2) {
                constrained_dofs.insert(3 * v_i + 1);
                constrained_dofs.insert(3 * v_i + 2);
            }
        }

        std::vector<int> free_dofs;
        for (std::size_t d = 0; d < ndof; ++d) {
            if (constrained_dofs.find(static_cast<int>(d)) == constrained_dofs.end()) {
                free_dofs.push_back(static_cast<int>(d));
            }
        }

        std::vector<double> u(ndof, 0.0);
        if (!free_dofs.empty()) {
            int nf = static_cast<int>(free_dofs.size());
            std::vector<double> K_ff(nf * nf, 0.0);
            std::vector<double> f_f(nf, 0.0);
            for (int r = 0; r < nf; ++r) {
                f_f[r] = f_glob[free_dofs[r]];
                for (int c = 0; c < nf; ++c) {
                    K_ff[r * nf + c] = K[free_dofs[r] * ndof + free_dofs[c]];
                }
            }
            auto u_free = solve_linear_system(nf, K_ff, f_f, 1);
            for (int r = 0; r < nf; ++r) {
                u[free_dofs[r]] = u_free[r];
            }
        }

        std::vector<double> R_full(ndof, 0.0);
        for (std::size_t i = 0; i < ndof; ++i) {
            double ku = 0.0;
            for (std::size_t j = 0; j < ndof; ++j) {
                ku += K[i * ndof + j] * u[j];
            }
            R_full[i] = f_glob[i] - ku;
        }

        std::vector<double> nodal_Rw(n_v, 0.0);
        std::vector<double> nodal_w(n_v, 0.0);
        for (std::size_t v_i = 0; v_i < n_v; ++v_i) {
            nodal_Rw[v_i] = R_full[3 * v_i];
            nodal_w[v_i] = u[3 * v_i];
        }

        std::map<std::pair<int, int>, double> sub_flux_weights;
        for (std::size_t k = 0; k < frame_sub_groups.size(); ++k) {
            for (std::size_t j = 0; j < frame_sub_groups[k].size(); ++j) {
                sub_flux_weights[{static_cast<int>(k), static_cast<int>(j)}] = frame_sub_groups[k][j].length;
            }
        }

        apportion_nodal_reactions_to_subsegments(
            mesh, frames, frame_sub_groups, sub_F, sub_M, node_reactions, nodal_Rw, sub_flux_weights);

        std::vector<double> overhang_torques(n_frames, 0.0);
        if (support_rank < 2) {
            for (std::size_t k = 0; k < n_frames; ++k) {
                const auto& subs = frame_sub_groups[k];
                if (subs.empty()) {
                    continue;
                }
                double dx = frames[k].p2.x - frames[k].p1.x;
                double dy = frames[k].p2.y - frames[k].p1.y;
                double lk = std::hypot(dx, dy);
                if (lk <= 1e-15) {
                    continue;
                }
                double tx = dx / lk, ty = dy / lk;
                std::set<int> f_verts;
                for (const auto& s : subs) {
                    f_verts.insert(s.u);
                    f_verts.insert(s.v);
                }
                double torq_k = 0.0;
                for (int v_i : f_verts) {
                    double mx_i = R_full[3 * v_i + 2];
                    double my_i = -R_full[3 * v_i + 1];
                    torq_k += mx_i * tx + my_i * ty;
                }
                overhang_torques[k] = torq_k;
            }
        }

        std::vector<int> tri_supp(mesh.triangles.size(), 0);
        std::vector<Point2D> tri_flux(mesh.triangles.size(), {0.0, 0.0});
        for (std::size_t t_idx = 0; t_idx < mesh.triangles.size(); ++t_idx) {
            auto [a, b, c] = mesh.triangles[t_idx];
            double x0 = verts[a].x, y0 = verts[a].y;
            double x1 = verts[b].x, y1 = verts[b].y;
            double x2 = verts[c].x, y2 = verts[c].y;
            double area2 = orient2d(x0, y0, x1, y1, x2, y2);
            double gwx = ((y1 - y2) * nodal_w[a] + (y2 - y0) * nodal_w[b] + (y0 - y1) * nodal_w[c]) / area2;
            double gwy = ((x2 - x1) * nodal_w[a] + (x0 - x2) * nodal_w[b] + (x1 - x0) * nodal_w[c]) / area2;
            double cx = (x0 + x1 + x2) / 3.0;
            double cy = (y0 + y1 + y2) / 3.0;
            auto [tied, _] = nearest_features(cx, cy);
            if (!tied.empty()) {
                tri_supp[t_idx] = tied[0].fid;
            }
            double gn = std::hypot(gwx, gwy);
            if (gn > 1e-15) {
                double sgn = (q >= 0.0) ? -1.0 : 1.0;
                tri_flux[t_idx] = {sgn * std::abs(q) * (gwx / gn), sgn * std::abs(q) * (gwy / gn)};
            } else if (!tied.empty()) {
                double vx = tied[0].cx - cx;
                double vy = tied[0].cy - cy;
                double vn = std::hypot(vx, vy);
                if (vn > 1e-15) {
                    tri_flux[t_idx] = {std::abs(q) * (vx / vn), std::abs(q) * (vy / vn)};
                }
            }
        }

        auto frame_results = finalize_frame_results(frames, frame_sub_groups, sub_F, sub_M, overhang_torques);
        return DistributionResult{
            "fea_plate",
            total_load,
            total_mx,
            total_my,
            std::move(frame_results),
            std::move(node_reactions),
            std::move(mesh),
            support_rank,
            std::move(vertex_flags),
            std::move(nodal_w),
            std::move(tri_supp),
            std::move(tri_flux),
            std::move(yield_lines)
        };
    }

    throw std::invalid_argument("Unknown load distribution mode: " + mode);
}

// ============================================================================
// 4. Guyan Static Condensation of In-Plane CST Diaphragm Stiffness & Mass
// ============================================================================

struct CondensationResult {
    int nb_dof = 0;
    int ni_dof = 0;
    std::vector<double> K_cond; // nb_dof x nb_dof row-major
    std::vector<double> M_cond; // nb_dof x nb_dof row-major
    std::vector<double> G;      // ni_dof x nb_dof row-major
};

inline CondensationResult condense_inplane_cst(
    const std::vector<Point2D>& vertices,
    const std::vector<std::array<int, 3>>& triangles,
    const std::vector<int>& boundary_node_indices,
    double E,
    double nu,
    double t,
    double rho = 0.0)
{
    int n_v = static_cast<int>(vertices.size());
    int ndof = 2 * n_v;
    double fac = E / (1.0 - nu * nu);
    std::array<double, 9> D_ps = {
        fac, fac * nu, 0.0,
        fac * nu, fac, 0.0,
        0.0, 0.0, fac * 0.5 * (1.0 - nu)
    };

    std::vector<double> K(ndof * ndof, 0.0);
    std::vector<double> M(ndof * ndof, 0.0);

    for (auto [a, b, c] : triangles) {
        double x0 = vertices[a].x, y0 = vertices[a].y;
        double x1 = vertices[b].x, y1 = vertices[b].y;
        double x2 = vertices[c].x, y2 = vertices[c].y;
        double area2 = orient2d(x0, y0, x1, y1, x2, y2);
        double area = 0.5 * area2;
        double b0 = (y1 - y2) / area2, b1 = (y2 - y0) / area2, b2 = (y0 - y1) / area2;
        double c0 = (x2 - x1) / area2, c1 = (x0 - x2) / area2, c2 = (x1 - x0) / area2;

        std::array<std::array<double, 6>, 3> B = {{
            {b0, 0.0, b1, 0.0, b2, 0.0},
            {0.0, c0, 0.0, c1, 0.0, c2},
            {c0, b0, c1, b1, c2, b2}
        }};

        std::array<int, 6> dofs = {2 * a, 2 * a + 1, 2 * b, 2 * b + 1, 2 * c, 2 * c + 1};
        for (int i_l = 0; i_l < 6; ++i_l) {
            for (int j_l = 0; j_l < 6; ++j_l) {
                double s = 0.0;
                for (int r = 0; r < 3; ++r) {
                    for (int col = 0; col < 3; ++col) {
                        s += B[r][i_l] * D_ps[r * 3 + col] * B[col][j_l];
                    }
                }
                K[dofs[i_l] * ndof + dofs[j_l]] += (t * area) * s;
            }
        }
        if (rho > 0.0) {
            double m_fac = rho * t * area / 12.0;
            for (int i_n = 0; i_n < 3; ++i_n) {
                for (int j_n = 0; j_n < 3; ++j_n) {
                    double val = m_fac * ((i_n == j_n) ? 2.0 : 1.0);
                    M[dofs[2 * i_n] * ndof + dofs[2 * j_n]] += val;
                    M[dofs[2 * i_n + 1] * ndof + dofs[2 * j_n + 1]] += val;
                }
            }
        }
    }

    int n_b = static_cast<int>(boundary_node_indices.size());
    std::map<int, int> b_pos;
    for (int idx = 0; idx < n_b; ++idx) {
        b_pos[boundary_node_indices[idx]] = idx;
    }

    struct EdgeSlave {
        int ib0;
        int ib1;
        double s;
    };
    std::map<int, EdgeSlave> edge_slaves;
    std::vector<int> i_nodes;

    for (int v = 0; v < n_v; ++v) {
        if (b_pos.find(v) != b_pos.end()) {
            continue;
        }
        double vx = vertices[v].x;
        double vy = vertices[v].y;
        bool slaved = false;
        for (int idx_b = 0; idx_b < n_b; ++idx_b) {
            int u_b = boundary_node_indices[idx_b];
            int w_b = boundary_node_indices[(idx_b + 1) % n_b];
            double ax = vertices[u_b].x, ay = vertices[u_b].y;
            double bx = vertices[w_b].x, by = vertices[w_b].y;
            double seg_len = std::hypot(bx - ax, by - ay);
            auto psd = point_segment_distance(vx, vy, ax, ay, bx, by);
            if (psd.is_interior && psd.dist <= 1e-7 * std::max(1.0, seg_len)) {
                edge_slaves[v] = {idx_b, (idx_b + 1) % n_b, psd.t};
                slaved = true;
                break;
            }
        }
        if (!slaved) {
            i_nodes.push_back(v);
        }
    }

    int n_i = static_cast<int>(i_nodes.size());
    std::map<int, int> i_pos;
    for (int idx = 0; idx < n_i; ++idx) {
        i_pos[i_nodes[idx]] = idx;
    }

    int n_master = 2 * (n_b + n_i);
    std::vector<double> P(ndof * n_master, 0.0);
    for (int v = 0; v < n_v; ++v) {
        auto it_b = b_pos.find(v);
        if (it_b != b_pos.end()) {
            int col = 2 * it_b->second;
            P[(2 * v) * n_master + col] = 1.0;
            P[(2 * v + 1) * n_master + (col + 1)] = 1.0;
        } else {
            auto it_s = edge_slaves.find(v);
            if (it_s != edge_slaves.end()) {
                int c0 = 2 * it_s->second.ib0;
                int c1 = 2 * it_s->second.ib1;
                double s_val = it_s->second.s;
                P[(2 * v) * n_master + c0] = 1.0 - s_val;
                P[(2 * v + 1) * n_master + (c0 + 1)] = 1.0 - s_val;
                P[(2 * v) * n_master + c1] = s_val;
                P[(2 * v + 1) * n_master + (c1 + 1)] = s_val;
            } else {
                int col = 2 * (n_b + i_pos[v]);
                P[(2 * v) * n_master + col] = 1.0;
                P[(2 * v + 1) * n_master + (col + 1)] = 1.0;
            }
        }
    }

    // Compute K_red = P^T K P and M_red = P^T M P
    std::vector<double> KP(ndof * n_master, 0.0);
    std::vector<double> MP(ndof * n_master, 0.0);
    for (int i = 0; i < ndof; ++i) {
        for (int k = 0; k < ndof; ++k) {
            double kik = K[i * ndof + k];
            double mik = M[i * ndof + k];
            if (kik == 0.0 && mik == 0.0) continue;
            for (int j = 0; j < n_master; ++j) {
                double pkj = P[k * n_master + j];
                if (pkj != 0.0) {
                    KP[i * n_master + j] += kik * pkj;
                    MP[i * n_master + j] += mik * pkj;
                }
            }
        }
    }

    std::vector<double> K_red(n_master * n_master, 0.0);
    std::vector<double> M_red(n_master * n_master, 0.0);
    for (int k = 0; k < ndof; ++k) {
        for (int i = 0; i < n_master; ++i) {
            double pki = P[k * n_master + i];
            if (pki == 0.0) continue;
            for (int j = 0; j < n_master; ++j) {
                K_red[i * n_master + j] += pki * KP[k * n_master + j];
                M_red[i * n_master + j] += pki * MP[k * n_master + j];
            }
        }
    }

    int nb_dof = 2 * n_b;
    int ni_dof = 2 * n_i;
    CondensationResult out{};
    out.nb_dof = nb_dof;
    out.ni_dof = ni_dof;
    out.K_cond.assign(nb_dof * nb_dof, 0.0);
    out.M_cond.assign(nb_dof * nb_dof, 0.0);

    for (int i = 0; i < nb_dof; ++i) {
        for (int j = 0; j < nb_dof; ++j) {
            out.K_cond[i * nb_dof + j] = K_red[i * n_master + j];
            out.M_cond[i * nb_dof + j] = M_red[i * n_master + j];
        }
    }

    if (ni_dof == 0) {
        for (int i = 0; i < nb_dof; ++i) {
            for (int j = i + 1; j < nb_dof; ++j) {
                double ks = 0.5 * (out.K_cond[i * nb_dof + j] + out.K_cond[j * nb_dof + i]);
                out.K_cond[i * nb_dof + j] = ks;
                out.K_cond[j * nb_dof + i] = ks;
                double ms = 0.5 * (out.M_cond[i * nb_dof + j] + out.M_cond[j * nb_dof + i]);
                out.M_cond[i * nb_dof + j] = ms;
                out.M_cond[j * nb_dof + i] = ms;
            }
        }
        return out;
    }

    std::vector<double> K_ii(ni_dof * ni_dof, 0.0);
    std::vector<double> neg_K_ib(ni_dof * nb_dof, 0.0);
    for (int i = 0; i < ni_dof; ++i) {
        for (int j = 0; j < ni_dof; ++j) {
            K_ii[i * ni_dof + j] = K_red[(nb_dof + i) * n_master + (nb_dof + j)];
        }
        for (int j = 0; j < nb_dof; ++j) {
            neg_K_ib[i * nb_dof + j] = -K_red[(nb_dof + i) * n_master + j];
        }
    }

    out.G = solve_linear_system(ni_dof, K_ii, neg_K_ib, nb_dof);

    // K_cond = K_bb + K_bi * G
    for (int i = 0; i < nb_dof; ++i) {
        for (int k = 0; k < ni_dof; ++k) {
            double k_bik = K_red[i * n_master + (nb_dof + k)];
            for (int j = 0; j < nb_dof; ++j) {
                out.K_cond[i * nb_dof + j] += k_bik * out.G[k * nb_dof + j];
            }
        }
    }

    // M_cond = M_bb + M_bi * G + G^T * M_ib + G^T * M_ii * G
    std::vector<double> M_ii_G(ni_dof * nb_dof, 0.0);
    for (int i = 0; i < ni_dof; ++i) {
        for (int k = 0; k < ni_dof; ++k) {
            double m_iik = M_red[(nb_dof + i) * n_master + (nb_dof + k)];
            for (int j = 0; j < nb_dof; ++j) {
                M_ii_G[i * nb_dof + j] += m_iik * out.G[k * nb_dof + j];
            }
        }
    }

    for (int i = 0; i < nb_dof; ++i) {
        for (int k = 0; k < ni_dof; ++k) {
            double m_bik = M_red[i * n_master + (nb_dof + k)];
            double g_ki = out.G[k * nb_dof + i];
            for (int j = 0; j < nb_dof; ++j) {
                double m_ib_kj = M_red[(nb_dof + k) * n_master + j];
                out.M_cond[i * nb_dof + j] += m_bik * out.G[k * nb_dof + j]
                                           + g_ki * m_ib_kj
                                           + g_ki * M_ii_G[k * nb_dof + j];
            }
        }
    }

    // Symmetrize K_cond and M_cond
    for (int i = 0; i < nb_dof; ++i) {
        for (int j = i + 1; j < nb_dof; ++j) {
            double ks = 0.5 * (out.K_cond[i * nb_dof + j] + out.K_cond[j * nb_dof + i]);
            out.K_cond[i * nb_dof + j] = ks;
            out.K_cond[j * nb_dof + i] = ks;
            double ms = 0.5 * (out.M_cond[i * nb_dof + j] + out.M_cond[j * nb_dof + i]);
            out.M_cond[i * nb_dof + j] = ms;
            out.M_cond[j * nb_dof + i] = ms;
        }
    }
#ifdef LADRUNO_MUTATE_MEMBRANE
    std::fill(out.K_cond.begin(), out.K_cond.end(), 0.0);
#endif

    return out;
}

} // namespace LadrunoMembraneKernel

#endif // LadrunoMembraneKernel_h
