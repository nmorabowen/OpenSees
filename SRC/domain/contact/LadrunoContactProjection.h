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
#include "LadrunoContactNormalField.h"   // ADR-63 #4a: smoothNormal() blend for evalSegmentSmooth()

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

// ADR-63 #4a — smoothed-normal one-shot segment evaluation. Identical to evalSegment() EXCEPT the
// normal is the averaged nodal-normal field interpolated at the projection (n_smooth = Σ N_i·n_i,
// normalized) instead of the per-segment faceted normalOriented(). The closest-point projection
// (x̄, ξ, η) AND the in-bounds + penetration gate are UNCHANGED (the FACET projection) ⇒ the active
// set is byte-identical to the faceted path: smoothing changes the force DIRECTION, not WHICH pairs
// are active, so R5 slide-off is untouched (ADR-63 BLOCKER-PROJECTION-CONSISTENCY). nodalNorm = the
// segment's nps nodal normals (from LadrunoContactNormalField::nodalNormals). refDir is the FALLBACK
// orientation: on a degenerate query-point blend (‖Σ N_i·n_i‖→0 — antiparallel corner normals across
// a folded element, or a node left zero by a degenerate nodal normal; BLOCKER-FALLBACK (b)) fall back
// to the faceted normalOriented(refDir) ⇒ never a garbage normal. Returns ACTIVE (penetrating + in-bounds).
inline bool evalSegmentSmooth(int nps, const double X[4][3], const double xs[3],
                              const double nodalNorm[4][3], const double refDir[3],
                              double &gap, double n[3], double N[4],
                              double tolR = 1e-12, int maxIt = 10) {
    double xi, eta;
    int st = project(nps, X, xs, xi, eta, tolR, maxIt);
    if (st != 0) return false;                       // oob or no valid projection (same gate as evalSegment)
    double dN1[4], dN2[4], xbar[3];
    shape(nps, xi, eta, N, dN1, dN2);
    // smoothed normal from the nodal field; fall back to the faceted oriented normal on a degenerate blend
    if (!LadrunoContactNormalField::smoothNormal(nps, N, nodalNorm, n)) {
        if (!normalOriented(nps, xi, eta, X, refDir, n)) return false;
    }
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

// ----------------------------------------- consistent ∂n/∂u normal tangent (B3)
// Second derivatives of the shape functions (CONSTANT in ξ,η for tri-3 / quad-4):
//   tri-3 (linear) ⇒ all 0; quad-4 (bilinear) ⇒ only the ξη cross term N_{i,ξη}≠0.
inline void shapeD2(int nps, double Nxx[4], double Nee[4], double Nxe[4]) {
    for (int i = 0; i < 4; i++) { Nxx[i] = 0.0; Nee[i] = 0.0; Nxe[i] = 0.0; }
    if (nps == 4) {
        // ∂²/∂ξ∂η of 0.25(1±ξ)(1±η): +/− 0.25 with the node-corner sign pattern.
        Nxe[0] =  0.25; Nxe[1] = -0.25; Nxe[2] =  0.25; Nxe[3] = -0.25;
    }
    // tri-3: zero (a planar facet; κ ≡ 0)
}

// The GEOMETRIC block of the consistent frictionless NTS tangent: K_geom = kn·gN·H,
// H = ∂B/∂u = ∂²gN/∂u² (the ∂n/∂u + ∂N_i/∂u curvature / large-sliding terms the
// shipped kn·BᵀB main term drops). Derived + FD-validated oracle-first in
// contact_prototypes/proto_b3_normal_tangent.py (analytic == FD ~1e-10 on a warped
// quad / tilted large-slide; SYMMETRIC to 1e-17; slave block ≡ 0 for a flat facet ⇒
// byte-identical to the shipped tangent for a flat/fixed master).
//
// Fills Kgeom[ndof][ndof] (ndof = 3·(1+nps) ≤ 15) and returns true iff ACTIVE
// (projection in-bounds AND penetrating). The caller adds fact·Kgeom to its tangent.
// X = master segment CURRENT coords (nps×3), xs = slave CURRENT position, refDir
// orients the normal (n·refDir > 0). Self-contained: re-projects (deterministic ⇒
// the SAME ξ̄/n/gap the residual used). Frictionless / normal-law only.
inline bool normalGeomTangent(int nps, const double X[4][3], const double xs[3],
                              const double refDir[3], double kn, double Kgeom[15][15]) {
    const int ndof = 3 * (1 + nps);
    for (int a = 0; a < ndof; a++)
        for (int b = 0; b < ndof; b++) Kgeom[a][b] = 0.0;

    double xi, eta;
    if (project(nps, X, xs, xi, eta, 1e-12, 10) != 0) return false;   // oob / no proj

    double N[4], dNxi[4], dNeta[4], Nxx[4], Nee[4], Nxe[4];
    shape(nps, xi, eta, N, dNxi, dNeta);
    shapeD2(nps, Nxx, Nee, Nxe);
    double xbar[3], a1[3], a2[3];
    interp(nps, N, X, xbar);
    tangents(nps, dNxi, dNeta, X, a1, a2);
    double raw[3]; cross3(a1, a2, raw);
    double jn = norm3(raw);
    if (jn < 1e-300) return false;
    double n[3]; for (int d = 0; d < 3; d++) n[d] = raw[d] / jn;
    double pr = dot3(n, refDir);
    if (std::fabs(pr) < 1e-12 * (norm3(refDir) + 1e-300)) return false;  // ambiguous
    if (pr < 0.0) for (int d = 0; d < 3; d++) n[d] = -n[d];
    double dvec[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    double gN = dot3(n, dvec);
    if (gN >= 0.0) return false;                       // not penetrating

    // covariant metric m_αβ, inverse, contravariant tangents a^α = m^{αβ} a_β
    double m00 = dot3(a1,a1), m01 = dot3(a1,a2), m11 = dot3(a2,a2);
    double detm = m00*m11 - m01*m01;
    if (std::fabs(detm) < 1e-300) return false;
    double mi00 = m11/detm, mi01 = -m01/detm, mi11 = m00/detm;
    double au1[3], au2[3];
    for (int d = 0; d < 3; d++) {
        au1[d] = mi00*a1[d] + mi01*a2[d];
        au2[d] = mi01*a1[d] + mi11*a2[d];
    }
    // second fundamental form κ_αβ = n·a_{α,β}, a_{α,β} = Σ N_{i,αβ} X_i
    double axx[3] = {0,0,0}, aee[3] = {0,0,0}, axe[3] = {0,0,0};
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) {
            axx[d] += Nxx[i]*X[i][d]; aee[d] += Nee[i]*X[i][d]; axe[d] += Nxe[i]*X[i][d];
        }
    double k00 = dot3(n,axx), k01 = dot3(n,axe), k11 = dot3(n,aee);   // κ (k10=k01)
    // projection-sensitivity metric H2_αβ = m_αβ − gN·κ_αβ, and its inverse
    double H200 = m00 - gN*k00, H201 = m01 - gN*k01, H211 = m11 - gN*k11;
    double detH2 = H200*H211 - H201*H201;
    if (std::fabs(detH2) < 1e-300) return false;
    double Hi00 = H211/detH2, Hi01 = -H201/detH2, Hi11 = H200/detH2;

    // per-DOF projected components Q_β[l], P_β[l] (l = slave xyz, then master nodes)
    double Q0[15], Q1[15], P0[15], P1[15];
    for (int l = 0; l < ndof; l++) { Q0[l]=Q1[l]=P0[l]=P1[l]=0.0; }
    for (int d = 0; d < 3; d++) { Q0[d] = a1[d]; Q1[d] = a2[d]; }      // slave: P=0
    for (int j = 0; j < nps; j++)
        for (int d = 0; d < 3; d++) {
            int l = 3*(1+j) + d;
            Q0[l] = -N[j]*a1[d]; Q1[l] = -N[j]*a2[d];
            P0[l] = dNxi[j]*n[d]; P1[l] = dNeta[j]*n[d];
        }
    // ∂ξ̄^α/∂u_l = H2^{αβ}(Q_β + gN·P_β); ∂n/∂u_l = −(c0 a^1 + c1 a^2),
    //   c_g = P_g[l] + κ_g0 ζ0 + κ_g1 ζ1. Then assemble H[:,l] and scale by kn·gN.
    for (int l = 0; l < ndof; l++) {
        double r0 = Q0[l] + gN*P0[l], r1 = Q1[l] + gN*P1[l];
        double Z0 = Hi00*r0 + Hi01*r1, Z1 = Hi01*r0 + Hi11*r1;
        double c0 = P0[l] + k00*Z0 + k01*Z1;
        double c1 = P1[l] + k01*Z0 + k11*Z1;
        double dn[3];
        for (int d = 0; d < 3; d++) dn[d] = -(c0*au1[d] + c1*au2[d]);
        double w = kn * gN;
        for (int d = 0; d < 3; d++) Kgeom[d][l] = w * dn[d];           // slave rows
        for (int i = 0; i < nps; i++) {
            double dNil = dNxi[i]*Z0 + dNeta[i]*Z1;
            for (int d = 0; d < 3; d++)
                Kgeom[3*(1+i)+d][l] = w * (-dNil*n[d] - N[i]*dn[d]);   // master-i rows
        }
    }
    return true;
}

} // namespace LadrunoContactProjection

#endif
