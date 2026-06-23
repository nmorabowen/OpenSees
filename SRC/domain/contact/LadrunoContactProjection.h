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

// ADR-39 P2b / ADR-41 A2 — LadrunoContactProjection: header-only, OpenSees-free
// closest-point projection geometry for node-to-segment (NTS) and segment-to-segment
// (mortar) contact. EXTRACTED from LadrunoContactKernel.h (PR #354) so BOTH the NTS
// penalty path AND the ADR-41 mortar Gauss-point loop consume ONE projection — no
// duplicated geometry (the staleness trap this fork keeps fighting). NO OpenSees types
// — only raw double[] + <cmath>. Mirrors the LadrunoJ2Kernel.h discipline.
//
// The NTS path needs {x̄, gap, n, φ}; the mortar GP loop ADDITIONALLY needs the
// covariant tangents t1=g1, t2=g2 and the surface metric g_ab = g_a·g_b (the slip
// return uses R = g·r, NOT I₂ — ADR-41 §FrictionalLaw / capstone ROAD-2). So:
//   evalSegment()  — NTS one-shot: project + oriented normal + gap; PENETRATING-only bool.
//   projectFull()  — mortar: geometry for ANY converged projection (gap as DATA, no
//                    penetration/in-bounds gate; the mortar clip already bounds the GP).
//
// Bit-for-bit with the shipped LadrunoContactKernel.h projection (this IS that code,
// moved): bounded closest-point Newton (cap maxIt; reject |detK| < 1e-14·|g1||g2| on a
// degenerate/collinear segment), DERIVED winding-immune oriented normal (design-gate
// BLOCKER-1/2). Oracle: contact_prototypes/proto_a2_projection.py.
//
// See Ladruno_implementation/{39_ladruno_contact_domain,41_ladruno_mortar_alm_contact,
// 48_ladruno_contact_capstone}_adr.md.

#ifndef LadrunoContactProjection_h
#define LadrunoContactProjection_h

#include <cmath>

namespace LadrunoContactProjection {

// ----------------------------------------------------------------- small vec3
inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline void cross3(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }

// ------------------------------------------------------- shape fns + derivs
// nps = 3 (linear tri, (xi,eta) area coords) or 4 (bilinear quad on [-1,1]^2).
// Fills N[nps], dNxi[nps], dNeta[nps].
inline void shape(int nps, double xi, double eta,
                  double N[4], double dNxi[4], double dNeta[4]) {
    if (nps == 4) {
        N[0] = 0.25*(1-xi)*(1-eta);  N[1] = 0.25*(1+xi)*(1-eta);
        N[2] = 0.25*(1+xi)*(1+eta);  N[3] = 0.25*(1-xi)*(1+eta);
        dNxi[0] = -0.25*(1-eta);  dNxi[1] =  0.25*(1-eta);
        dNxi[2] =  0.25*(1+eta);  dNxi[3] = -0.25*(1+eta);
        dNeta[0] = -0.25*(1-xi);  dNeta[1] = -0.25*(1+xi);
        dNeta[2] =  0.25*(1+xi);  dNeta[3] =  0.25*(1-xi);
    } else { // tri-3
        N[0] = 1 - xi - eta;  N[1] = xi;  N[2] = eta;  N[3] = 0.0;
        dNxi[0] = -1.0;  dNxi[1] = 1.0;  dNxi[2] = 0.0;  dNxi[3] = 0.0;
        dNeta[0] = -1.0; dNeta[1] = 0.0; dNeta[2] = 1.0; dNeta[3] = 0.0;
    }
}

// xbar = Σ N_i X_i ; tangents g1 = ∂x̄/∂ξ, g2 = ∂x̄/∂η.  X is nps×3.
inline void interp(int nps, const double N[4], const double X[4][3], double xbar[3]) {
    for (int d = 0; d < 3; d++) {
        xbar[d] = 0.0;
        for (int i = 0; i < nps; i++) xbar[d] += N[i] * X[i][d];
    }
}
inline void tangents(int nps, const double dNxi[4], const double dNeta[4],
                     const double X[4][3], double g1[3], double g2[3]) {
    for (int d = 0; d < 3; d++) {
        g1[d] = 0.0; g2[d] = 0.0;
        for (int i = 0; i < nps; i++) { g1[d] += dNxi[i]*X[i][d]; g2[d] += dNeta[i]*X[i][d]; }
    }
}

// in-bounds test on the parent domain (small tol slack).
inline bool inBounds(int nps, double xi, double eta) {
    const double t = 1e-9;
    if (nps == 4) return (xi >= -1 - t && xi <= 1 + t && eta >= -1 - t && eta <= 1 + t);
    return (xi >= -t && eta >= -t && xi + eta <= 1 + t);
}

// ----------------------------------------------------- bounded projection
// BOUNDED closest-point Newton on (x_s − x̄)·g_α = 0 (design-gate BLOCKER-2):
// cap maxIt; reject if |detK| < 1e-14*|g1||g2| (degenerate/collinear segment).
// Gauss-Newton metric K = [[g1·g1,g1·g2],[g2·g1,g2·g2]] (drops the O(d) curvature
// correction — SPD for non-degenerate segments; the de-risk path).
// returns: 0 converged in-bounds, 1 converged out-of-bounds, -1 no valid projection.
inline int project(int nps, const double X[4][3], const double xs[3],
                   double &xi, double &eta, double tolR, int maxIt) {
    xi = (nps == 4) ? 0.0 : (1.0/3.0);
    eta = (nps == 4) ? 0.0 : (1.0/3.0);
    double N[4], dNxi[4], dNeta[4], xbar[3], g1[3], g2[3];
    bool converged = false;
    for (int it = 0; it < maxIt; it++) {
        shape(nps, xi, eta, N, dNxi, dNeta);
        interp(nps, N, X, xbar);
        tangents(nps, dNxi, dNeta, X, g1, g2);
        double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
        double R0 = dot3(d, g1), R1 = dot3(d, g2);
        if (std::sqrt(R0*R0 + R1*R1) < tolR) { converged = true; break; }
        double K00 = dot3(g1,g1), K01 = dot3(g1,g2), K11 = dot3(g2,g2);
        double detK = K00*K11 - K01*K01;
        double scale = norm3(g1) * norm3(g2);
        if (std::fabs(detK) < 1e-14 * (scale + 1e-300)) return -1;  // degenerate
        double dxi  = ( K11*R0 - K01*R1) / detK;
        double deta = (-K01*R0 + K00*R1) / detK;
        xi += dxi; eta += deta;
    }
    if (!converged) return -1;
    return inBounds(nps, xi, eta) ? 0 : 1;
}

// derived outward normal, oriented so n·refDir > 0 (refDir points toward the slave's
// allowed half-space). Orienting by a fixed DIRECTION (not the live slave position)
// keeps the sign correct even after the slave penetrates. Winding-immune: reversing
// the node order flips the raw cross product but the refDir test flips it back
// (design-gate BLOCKER-1). refDir need not be unit.
inline bool normalOriented(int nps, double xi, double eta, const double X[4][3],
                           const double refDir[3], double n[3]) {
    double N[4], dNxi[4], dNeta[4], xbar[3], g1[3], g2[3], raw[3];
    shape(nps, xi, eta, N, dNxi, dNeta);
    interp(nps, N, X, xbar);
    tangents(nps, dNxi, dNeta, X, g1, g2);
    cross3(g1, g2, raw);
    double j = norm3(raw);
    if (j < 1e-300) return false;
    for (int d = 0; d < 3; d++) n[d] = raw[d] / j;
    // FAIL-SAFE: if refDir is (numerically) perpendicular to n the outward sense is
    // ambiguous (e.g. the slave reference point lies in the segment plane). Refuse
    // the pair rather than guess a sign — the caller should pass an explicit
    // -outward direction. (Gate H2.)
    double proj = dot3(n, refDir);
    if (std::fabs(proj) < 1e-12 * (norm3(refDir) + 1e-300)) return false;
    if (proj < 0.0) for (int d = 0; d < 3; d++) n[d] = -n[d];
    return true;
}

// One-shot segment evaluation: project xs onto the segment, derive the oriented
// normal (toward refDir), and the gap. Returns the ACTIVE flag (penetrating +
// in-bounds). On active, fills gap (<0), n[3], N[nps] (shape at the projection).
inline bool evalSegment(int nps, const double X[4][3], const double xs[3],
                        const double refDir[3], double &gap, double n[3], double N[4],
                        double tolR = 1e-12, int maxIt = 10) {
    double xi, eta;
    int st = project(nps, X, xs, xi, eta, tolR, maxIt);
    if (st != 0) return false;                       // oob or no valid projection
    if (!normalOriented(nps, xi, eta, X, refDir, n)) return false;
    double dN1[4], dN2[4], xbar[3];
    shape(nps, xi, eta, N, dN1, dN2);
    interp(nps, N, X, xbar);
    double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    gap = dot3(n, d);
    return (gap < 0.0);                              // penetrating only
}

// ------------------------------------------------ rich projection (mortar GP loop)
// Everything the segment-to-segment Gauss-point loop needs at one slave point. Unlike
// evalSegment(), projectFull() does NOT gate on penetration or in-bounds — the mortar
// overlap clip has already restricted the GP to the slave/master overlap, so gap and
// status are returned as DATA for the caller to use (weighted gap g̃, active set, etc.).
struct Projection {
    int    status;       // 0 conv in-bounds, 1 conv out-of-bounds, -1 no valid projection
    double xi, eta;      // x̄ master parametric coords at the closest point
    double gap;          // gN = n·(xs − x̄)  (signed; < 0 penetration)
    double n[3];         // oriented unit normal (toward refDir)
    double t1[3], t2[3]; // covariant tangents g1 = ∂x̄/∂ξ, g2 = ∂x̄/∂η  (NOT unit)
    double g[2][2];      // covariant surface metric g_ab = g_a · g_b  (the slip uses R=g·r)
    double phi[4];       // master shape fns φ_i at the projection (φ[3]=0 for tri-3)
};

// returns status (also in p.status): 0 conv in-bounds, 1 conv out-of-bounds, -1 no
// valid projection (degenerate segment) OR ambiguous oriented normal (refDir ⟂ n).
inline int projectFull(int nps, const double X[4][3], const double xs[3],
                       const double refDir[3], Projection &p,
                       double tolR = 1e-12, int maxIt = 10) {
    // zero the struct so a failed/garbage projection never leaks stale geometry
    p.gap = 0.0;
    for (int d = 0; d < 3; d++) { p.n[d] = 0.0; p.t1[d] = 0.0; p.t2[d] = 0.0; }
    p.g[0][0] = p.g[0][1] = p.g[1][0] = p.g[1][1] = 0.0;
    for (int i = 0; i < 4; i++) p.phi[i] = 0.0;

    int st = project(nps, X, xs, p.xi, p.eta, tolR, maxIt);
    p.status = st;
    if (st < 0) return st;                            // no valid projection (sentinel)

    double dNxi[4], dNeta[4], g1[3], g2[3], xbar[3];
    shape(nps, p.xi, p.eta, p.phi, dNxi, dNeta);
    interp(nps, p.phi, X, xbar);
    tangents(nps, dNxi, dNeta, X, g1, g2);
    if (!normalOriented(nps, p.xi, p.eta, X, refDir, p.n)) { p.status = -1; return -1; }

    for (int d = 0; d < 3; d++) { p.t1[d] = g1[d]; p.t2[d] = g2[d]; }
    p.g[0][0] = dot3(g1, g1);  p.g[0][1] = dot3(g1, g2);
    p.g[1][0] = p.g[0][1];     p.g[1][1] = dot3(g2, g2);
    double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    p.gap = dot3(p.n, d);
    return st;                                        // 0 in-bounds, 1 out-of-bounds
}

} // namespace LadrunoContactProjection

#endif
