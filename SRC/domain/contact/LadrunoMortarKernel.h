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

// ADR-41 C1 — LadrunoMortarKernel: header-only, OpenSees-free segment-to-segment
// MORTAR integration. The clip → sub-triangle Gauss → mortar matrices D_IJ, M_IK and
// the weighted gap g̃_I, for ONE slave/master facet pair. NO OpenSees types — only raw
// double[] + <cmath> (and the SHIPPED LadrunoContactProjection); numpy-oracle-testable
// without building OpenSees. Mirrors the LadrunoJ2Kernel.h / LadrunoFrictionKernel.h
// discipline.
//
// This is the first NET-NEW mortar mechanics in the fork (A1/A2 were extractions). The
// hard part is the OVERLAP CLIP — the only genuinely new algorithm (Sutherland-Hodgman
// of slave∩master on a common auxiliary plane, fan-triangulated). The projection of each
// Gauss point onto the master facet REUSES LadrunoContactProjection::projectFull — the
// geometry is NOT re-derived here (the fork's one-projection / no-staleness rule).
//
// STANDARD (not dual) basis — ADR-47 owns the dual/biorthogonal basis:
//   D_IJ += w·J · N_I^s · N_J^s              (slave–slave mortar matrix)
//   M_IK += w·J · N_I^s · φ_K^m              (slave–master mortar matrix)
//   g̃_I += w·J · N_I^s · ( (x_s − x_m)·n )   (weighted gap, NORMAL INSIDE the integral —
//                                             ADR-41 fixed D1's flat-facet factoring)
// with J = A_aux / |n_s·n0|  (flat-facet area ratio: aux-plane cell → slave-surface area;
// for a flat slave facet this ratio is constant and exact). C1 is geometry + the patch
// test ONLY: it returns D, M, g̃ and the per-GP geometry — NO traction, NO λ, NO ALM
// (enforcement = C2; friction via LadrunoFrictionKernel = C3).
//
// SCOPE / known limitations (ADR-41 C1 adversarial review — see [[LEDGER_quirks]]):
//   - FLAT, CONVEX tri/quad facets only. Non-convex (concave/bow-tie) facets are REFUSED
//     (integratePair returns -1 via isConvex2) — they would otherwise integrate over the
//     convex hull / truncate the clip. A NON-PLANAR (warped) slave facet is NOT refused
//     but the constant `J = A_aux/|n_s·n0|` then biases the area measure (≈0.7% at a
//     0.3·edge out-of-plane lift); the partition-of-unity and patch-test gates still pass
//     because they reduce to Σφ=1 (measure-independent). The exact fix — a per-sub-triangle
//     slave Jacobian |g1×g2| — is C2/ADR-47 work; flagged, not solved, in C1.
//
// Gates (proved oracle-first in contact_prototypes/proto_c1_mortar.py, 24/24):
//   - Sutherland-Hodgman clip area vs known cases (incl. the 2(√2−1) octagon).
//   - partition of unity  Σ_K M_IK == Σ_J D_IJ  to 1e-12 (matched AND non-matched).
//   - CONSTANT-pressure patch test on a NON-MATCHED 2-block mesh: P=D⁻¹M transfers a
//     uniform field to ≤1e-6 (THE headline falsifier).
//   - clip is EXACT at low order (M(o2)==M(o4)); the clip-blind rule is O(40%) biased and
//     only converges under heavy refinement — why the clip is in the MVP, not deferred.
//   - weighted gap g̃ linear-exact (uniform + tilted-plane gap); FD d g̃/d(coord) check.
//
// See Ladruno_implementation/{41_ladruno_mortar_alm_contact,48_ladruno_contact_capstone}
// _adr.md and _adr41_c1_design.md.

#ifndef LadrunoMortarKernel_h
#define LadrunoMortarKernel_h

#include <cmath>
#include "LadrunoContactProjection.h"

namespace LadrunoMortarKernel {

namespace LP = LadrunoContactProjection;   // shipped projection (projectFull, shape, ...)

// limits: a facet has ≤4 nodes; a convex tri/quad ∩ tri/quad clip has ≤8 vertices; the
// fan triangulation has ≤8 sub-triangles; the GP rule has ≤6 points per sub-triangle.
static const int MAXV  = 8;    // max clip-polygon vertices
static const int MAXGP = 6;    // max Gauss points per sub-triangle

// ------------------------------------------------------------------ result of one pair
struct PairResult {
    int    status;     // 0 contribution accumulated, -1 empty/degenerate overlap (skip)
    int    ngp;        // number of Gauss points integrated (0 ⇒ no contribution)
    double area;       // ∫ dΓ over the overlap (slave-surface measure)
    double D[4][4];    // slave–slave mortar  D_IJ   (rows/cols = slave nodes)
    double M[4][4];    // slave–master mortar M_IK   (rows = slave, cols = master nodes)
    double g[4];       // weighted gap g̃_I           (one per slave node)
};

// =============================================================== small 2D helpers
inline double signedArea2(const double P[][2], int n) {
    double a = 0.0;
    for (int i = 0; i < n; i++) { int j = (i + 1) % n; a += P[i][0]*P[j][1] - P[j][0]*P[i][1]; }
    return 0.5 * a;
}

// reverse the winding of an n-vertex polygon in place (to enforce CCW).
inline void reverse2(double P[][2], int n) {
    for (int i = 0, j = n - 1; i < j; i++, j--) {
        double tx = P[i][0], ty = P[i][1];
        P[i][0] = P[j][0]; P[i][1] = P[j][1];
        P[j][0] = tx;      P[j][1] = ty;
    }
}
inline void ensureCCW(double P[][2], int n) { if (signedArea2(P, n) < 0.0) reverse2(P, n); }

// convexity test for a CCW polygon: every consecutive turn must be a (weak) left turn.
// The mortar path REQUIRES convex projected facets — Sutherland-Hodgman needs a convex
// CLIP polygon, and a non-convex SUBJECT would be silently integrated over its convex
// hull (wrong area) / truncate the clip buffer (diverging from the oracle). A flat,
// convex tri/quad facet always projects to a convex polygon; a concave/bow-tie or
// strongly-warped facet does not ⇒ the caller REFUSES the pair (gate, ADR-41 C1 review).
inline bool isConvex2(const double P[][2], int n, double tol = 1e-12) {
    if (n < 3) return false;
    double scale = 0.0;                          // edge-length scale for a relative tol
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double ex = P[j][0]-P[i][0], ey = P[j][1]-P[i][1];
        double e = std::sqrt(ex*ex + ey*ey);
        if (e > scale) scale = e;
    }
    if (scale < 1e-300) return false;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n, k = (i + 2) % n;
        double ax = P[j][0]-P[i][0], ay = P[j][1]-P[i][1];
        double bx = P[k][0]-P[j][0], by = P[k][1]-P[j][1];
        if (ax*by - ay*bx < -tol * scale * scale) return false;   // a right turn ⇒ concave
    }
    return true;
}

// intersection of segment S->E with the infinite line A->B (Sutherland-Hodgman edge).
inline void iline(const double S[2], const double E[2], const double A[2],
                  const double B[2], double out[2]) {
    double d1x = E[0]-S[0], d1y = E[1]-S[1];
    double d2x = B[0]-A[0], d2y = B[1]-A[1];
    double den = d1x*d2y - d1y*d2x;
    if (std::fabs(den) < 1e-300) { out[0] = S[0]; out[1] = S[1]; return; }
    double t = ((A[0]-S[0])*d2y - (A[1]-S[1])*d2x) / den;
    out[0] = S[0] + t*d1x;  out[1] = S[1] + t*d1y;
}

// Sutherland-Hodgman: subject ∩ clip, both CCW, CLIP CONVEX (tri/quad ⇒ convex). Writes
// the (convex) intersection into out[][2], returns its vertex count (≤ ns+nc, ≤ MAXV).
inline int clipPolygon(const double subj[][2], int ns, const double clip[][2], int nc,
                       double out[][2]) {
    double buf[2][MAXV + 4][2];
    int cur = 0, n = ns;
    for (int i = 0; i < ns; i++) { buf[cur][i][0] = subj[i][0]; buf[cur][i][1] = subj[i][1]; }
    for (int e = 0; e < nc; e++) {
        const double *A = clip[e];
        const double *B = clip[(e + 1) % nc];
        double ex = B[0]-A[0], ey = B[1]-A[1];
        int nxt = 1 - cur, m = 0;
        if (n == 0) { cur = nxt; break; }
        const double *S = buf[cur][n - 1];
        double sS = ex*(S[1]-A[1]) - ey*(S[0]-A[0]);   // >0 ⇒ inside (left of A->B)
        for (int i = 0; i < n && m < MAXV + 3; i++) {
            const double *E = buf[cur][i];
            double sE = ex*(E[1]-A[1]) - ey*(E[0]-A[0]);
            if (sE >= 0.0) {
                if (sS < 0.0) { iline(S, E, A, B, buf[nxt][m]); m++; }
                buf[nxt][m][0] = E[0]; buf[nxt][m][1] = E[1]; m++;
            } else if (sS >= 0.0) {
                iline(S, E, A, B, buf[nxt][m]); m++;
            }
            S = E; sS = sE;
        }
        cur = nxt; n = m;
    }
    if (n > MAXV) n = MAXV;
    for (int i = 0; i < n; i++) { out[i][0] = buf[cur][i][0]; out[i][1] = buf[cur][i][1]; }
    return n;
}

// merge near-coincident consecutive vertices; returns the surviving count.
inline int dedupe(double P[][2], int n, double tol = 1e-12) {
    double tmp[MAXV][2]; int m = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double dx = P[i][0]-P[j][0], dy = P[i][1]-P[j][1];
        if (std::sqrt(dx*dx + dy*dy) > tol) { tmp[m][0] = P[i][0]; tmp[m][1] = P[i][1]; m++; }
    }
    for (int i = 0; i < m; i++) { P[i][0] = tmp[i][0]; P[i][1] = tmp[i][1]; }
    return m;
}

// =============================================================== facet geometry
// unit normal of a (flat) facet, oriented so n·refDir > 0. quad: cross of the diagonals.
inline bool facetNormal(int nps, const double X[4][3], const double refDir[3], double n[3]) {
    double a[3], b[3], raw[3];
    if (nps == 3) {
        for (int d = 0; d < 3; d++) { a[d] = X[1][d]-X[0][d]; b[d] = X[2][d]-X[0][d]; }
    } else {
        for (int d = 0; d < 3; d++) { a[d] = X[2][d]-X[0][d]; b[d] = X[3][d]-X[1][d]; }
    }
    LP::cross3(a, b, raw);
    double j = LP::norm3(raw);
    if (j < 1e-300) return false;
    for (int d = 0; d < 3; d++) n[d] = raw[d] / j;
    if (LP::dot3(n, refDir) < 0.0) for (int d = 0; d < 3; d++) n[d] = -n[d];
    return true;
}

// auxiliary plane from the MASTER facet: origin x0=centroid, normal n0, orthonormal axes.
inline bool auxPlane(int nps, const double X[4][3], const double refDir[3],
                     double x0[3], double n0[3], double e1[3], double e2[3]) {
    for (int d = 0; d < 3; d++) {
        x0[d] = 0.0;
        for (int i = 0; i < nps; i++) x0[d] += X[i][d];
        x0[d] /= nps;
    }
    if (!facetNormal(nps, X, refDir, n0)) return false;
    double t[3];
    for (int d = 0; d < 3; d++) t[d] = X[1][d] - X[0][d];
    double tn = LP::dot3(t, n0);
    for (int d = 0; d < 3; d++) e1[d] = t[d] - tn * n0[d];
    double e1n = LP::norm3(e1);
    if (e1n < 1e-300) return false;
    for (int d = 0; d < 3; d++) e1[d] /= e1n;
    LP::cross3(n0, e1, e2);
    return true;
}

inline void to2d(const double X[3], const double x0[3], const double e1[3],
                 const double e2[3], double uv[2]) {
    double d[3] = { X[0]-x0[0], X[1]-x0[1], X[2]-x0[2] };
    uv[0] = LP::dot3(d, e1);  uv[1] = LP::dot3(d, e2);
}

// recover ξ_s s.t. Σ N_I(ξ_s)·UV_I = q  (UV = slave nodes in aux 2D). tri-3: affine;
// quad: bounded inverse-bilinear Newton (same pattern as LadrunoContactProjection).
// Returns TRUE on convergence (residual < tolR / non-degenerate affine map); FALSE lets
// the caller SKIP that GP instead of integrating a wrong root (ADR-41 C1 review). With
// the convex-facet guard in integratePair this cannot fail in C1, but the de-risk flag
// mirrors LadrunoContactProjection::project — no silent last-iterate.
inline bool inverseIsomap2D(int nps, const double UV[4][2], const double q[2],
                            double &xi, double &eta, double tolR = 1e-13, int maxIt = 20) {
    if (nps == 3) {
        double a00 = UV[1][0]-UV[0][0], a01 = UV[2][0]-UV[0][0];
        double a10 = UV[1][1]-UV[0][1], a11 = UV[2][1]-UV[0][1];
        double det = a00*a11 - a01*a10;
        double scale = std::sqrt((a00*a00+a10*a10)*(a01*a01+a11*a11));
        if (std::fabs(det) < 1e-14 * (scale + 1e-300)) { xi = eta = 0.0; return false; }
        double rx = q[0]-UV[0][0], ry = q[1]-UV[0][1];
        xi  = ( a11*rx - a01*ry) / det;
        eta = (-a10*rx + a00*ry) / det;
        return true;
    }
    xi = 0.0; eta = 0.0;
    double N[4], dNxi[4], dNeta[4];
    bool conv = false;
    for (int it = 0; it < maxIt; it++) {
        LP::shape(4, xi, eta, N, dNxi, dNeta);
        double rx = -q[0], ry = -q[1], j00 = 0, j01 = 0, j10 = 0, j11 = 0;
        for (int i = 0; i < 4; i++) {
            rx += N[i]*UV[i][0];  ry += N[i]*UV[i][1];
            j00 += dNxi[i]*UV[i][0];  j01 += dNeta[i]*UV[i][0];
            j10 += dNxi[i]*UV[i][1];  j11 += dNeta[i]*UV[i][1];
        }
        if (std::sqrt(rx*rx + ry*ry) < tolR) { conv = true; break; }
        double det = j00*j11 - j01*j10;
        if (std::fabs(det) < 1e-300) break;
        xi  -= ( j11*rx - j01*ry) / det;
        eta -= (-j10*rx + j00*ry) / det;
    }
    return conv;
}

// barycentric triangle Gauss rule, weights SUM TO 1 (area applied separately). order<=2
// ⇒ 3-pt (degree 2); else the 6-pt degree-4 rule. Returns the point count.
inline int triQuadrature(int order, double bary[MAXGP][3], double w[MAXGP]) {
    if (order <= 2) {
        const double b[3][3] = {{2./3,1./6,1./6},{1./6,2./3,1./6},{1./6,1./6,2./3}};
        for (int i = 0; i < 3; i++) { for (int k = 0; k < 3; k++) bary[i][k] = b[i][k]; w[i] = 1./3; }
        return 3;
    }
    const double a1 = 0.445948490915965, b1 = 0.108103018168070, w1 = 0.223381589678011;
    const double a2 = 0.091576213509771, b2 = 0.816847572980459, w2 = 0.109951743655322;
    const double b[6][3] = {{b1,a1,a1},{a1,b1,a1},{a1,a1,b1},{b2,a2,a2},{a2,b2,a2},{a2,a2,b2}};
    const double ww[6]   = {w1,w1,w1,w2,w2,w2};
    for (int i = 0; i < 6; i++) { for (int k = 0; k < 3; k++) bary[i][k] = b[i][k]; w[i] = ww[i]; }
    return 6;
}

// =============================================================== the mortar pair
// Integrate ONE slave/master facet pair (the ADR-39 broad phase supplies the pairing).
// Accumulates D, M, g̃ into `out` (zeroed here). Returns 0 if a contribution was added,
// -1 if the overlap is empty / a sliver / degenerate (caller simply skips it).
inline int integratePair(int nps_s, const double Xs[4][3],
                         int nps_m, const double Xm[4][3],
                         const double refDir[3], PairResult &out, int order = 2) {
    out.status = -1; out.ngp = 0; out.area = 0.0;
    for (int i = 0; i < 4; i++) { out.g[i] = 0.0;
        for (int j = 0; j < 4; j++) { out.D[i][j] = 0.0; out.M[i][j] = 0.0; } }

    double x0[3], n0[3], e1[3], e2[3], n_s[3];
    if (!auxPlane(nps_m, Xm, refDir, x0, n0, e1, e2)) return -1;
    if (!facetNormal(nps_s, Xs, refDir, n_s)) return -1;
    double cos_t = std::fabs(LP::dot3(n_s, n0));
    if (cos_t < 1e-12) return -1;                       // slave ⟂ aux plane: ill-posed

    // slave nodes in aux 2D — NODE order (for the GP back-map) and a CCW copy (clip subj)
    double UVs[4][2], subj[4][2], clip[4][2];
    for (int i = 0; i < nps_s; i++) { to2d(Xs[i], x0, e1, e2, UVs[i]);
        subj[i][0] = UVs[i][0]; subj[i][1] = UVs[i][1]; }
    for (int i = 0; i < nps_m; i++) to2d(Xm[i], x0, e1, e2, clip[i]);
    ensureCCW(subj, nps_s);
    ensureCCW(clip, nps_m);
    // C1 scope = flat, CONVEX tri/quad facets. A concave/bow-tie or strongly-warped facet
    // projects to a non-convex polygon ⇒ Sutherland-Hodgman would silently integrate the
    // slave over its convex hull and truncate the clip buffer (diverging from the oracle).
    // REFUSE the pair rather than return a wrong area (gate, ADR-41 C1 review).
    if (!isConvex2(subj, nps_s) || !isConvex2(clip, nps_m)) return -1;

    double poly[MAXV][2];
    int np = clipPolygon(subj, nps_s, clip, nps_m, poly);
    np = dedupe(poly, np);
    if (np < 3) return -1;

    double A_s = std::fabs(signedArea2(subj, nps_s));
    double A_m = std::fabs(signedArea2(clip, nps_m));
    double A_min = (A_s < A_m) ? A_s : A_m;
    double A_poly = std::fabs(signedArea2(poly, np));
    if (A_poly < 1e-12 * A_min) return -1;              // sliver guard

    double cen[2] = {0.0, 0.0};
    for (int i = 0; i < np; i++) { cen[0] += poly[i][0]; cen[1] += poly[i][1]; }
    cen[0] /= np; cen[1] /= np;

    double bary[MAXGP][3], w[MAXGP];
    int ngpr = triQuadrature(order, bary, w);

    // fan from the centroid → sub-triangles; Gauss each
    for (int s = 0; s < np; s++) {
        const double *V1 = poly[s];
        const double *V2 = poly[(s + 1) % np];
        double triP[3][2] = { {cen[0],cen[1]}, {V1[0],V1[1]}, {V2[0],V2[1]} };
        double A_aux = std::fabs(signedArea2(triP, 3));
        if (A_aux < 1e-12 * A_min) continue;
        double A_phys = A_aux / cos_t;                  // aux-plane area → slave-surface area
        for (int gp = 0; gp < ngpr; gp++) {
            double q[2] = {
                bary[gp][0]*triP[0][0] + bary[gp][1]*triP[1][0] + bary[gp][2]*triP[2][0],
                bary[gp][0]*triP[0][1] + bary[gp][1]*triP[1][1] + bary[gp][2]*triP[2][1] };
            double xi_s, eta_s;
            if (!inverseIsomap2D(nps_s, UVs, q, xi_s, eta_s)) continue;  // back-map failed
            double Ns[4], dN1[4], dN2[4];
            LP::shape(nps_s, xi_s, eta_s, Ns, dN1, dN2);
            double x_s[3] = {0,0,0};
            for (int i = 0; i < nps_s; i++) for (int d = 0; d < 3; d++) x_s[d] += Ns[i]*Xs[i][d];

            LP::Projection p;
            int st = LP::projectFull(nps_m, Xm, x_s, refDir, p);
            if (st < 0) continue;                        // no valid projection: skip GP
            double x_m[3] = {0,0,0};
            for (int k = 0; k < nps_m; k++) for (int d = 0; d < 3; d++) x_m[d] += p.phi[k]*Xm[k][d];
            double dx[3] = { x_s[0]-x_m[0], x_s[1]-x_m[1], x_s[2]-x_m[2] };
            double gN = LP::dot3(p.n, dx);

            double wJ = w[gp] * A_phys;
            out.area += wJ;
            for (int i = 0; i < nps_s; i++) {
                out.g[i] += wJ * Ns[i] * gN;
                for (int j = 0; j < nps_s; j++) out.D[i][j] += wJ * Ns[i] * Ns[j];
                for (int k = 0; k < nps_m; k++) out.M[i][k] += wJ * Ns[i] * p.phi[k];
            }
            out.ngp++;
        }
    }
    out.status = (out.ngp > 0) ? 0 : -1;
    return out.status;
}

// =============================================================== linearization (C2)
// STUB — the ALM tangent K_c needs ∂{D,M,g̃}/∂(nodal coords) and ∂{n,ξ}. ADR-41 C1 is
// geometry + the patch test only; C2 assembles K_c and consumes these derivatives (the
// oracle FD-checks them — proto_c1_mortar.py T7). Left unimplemented on purpose so C1
// stays purely additive; do NOT wire a half-derivative here (the staleness trap).
struct PairLinearization {
    // dD[i][j][a][d], dM[i][k][a][d], dn[d][a][e], dxi[a][e]: derivatives w.r.t. the
    // a-th node's d-th coordinate. Sizes settle in C2; placeholder flag for now.
    int implemented;   // 0 = not yet (C2)
};
inline int linearizePair(int /*nps_s*/, const double /*Xs*/[4][3],
                         int /*nps_m*/, const double /*Xm*/[4][3],
                         const double /*refDir*/[3], PairLinearization &lin) {
    lin.implemented = 0;
    return -1;          // ADR-41 C2
}

} // namespace LadrunoMortarKernel

#endif
