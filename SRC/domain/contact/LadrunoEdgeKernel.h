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

// ADR-57 E0 — LadrunoEdgeKernel: header-only, OpenSees-free segment-to-segment
// (edge-edge) closest-point geometry for the perpendicular edge-edge narrow phase
// (the cos_t->0 case the mortar face-overlap clip degenerates on — two facets meeting
// edge-on, no in-plane overlap to integrate, and NTS misses the interior-line contact).
// NO OpenSees types — only raw double[] + <cmath>; numpy-oracle-testable build-free.
// Mirrors the LadrunoContactProjection.h / LadrunoMortarKernel.h discipline.
//
// E0 owns the GEOMETRY ONLY (no traction, no friction, no tangent — those are E2/E3/E4):
//   closestPtSegSeg() — ERICSON'S EXACT branch ladder (RTCD 5.1.9 ClosestPtSegmentSegment):
//     the conditional clamp ladder over all 8 Voronoi sub-regions, NOT a single
//     back-substitution (the ADR-57 3-reviewer gate Lens-A A-1). Returns the closest
//     parameters (s,t), the closest points, the degeneracy gauge denom, and a MARGIN'D
//     interior flag (a near-vertex s*=0.999 is NOT interior — Lens-A A-2).
//   edgeNormal()      — common-perpendicular UNIT normal n = (d1 x d2)/||.|| (defined only
//     for non-parallel edges; ill-conditioned ~1/sin(theta_edge) — gated by tau_par).
//   edgeGap()         — signed gap gN = w.n with a BODY-FIXED (committed) sign, NOT the
//     self-referential w.n rule that masks interpenetration through contact (Lens-A A-4).
//   bOperator()       — the first variation dgN = B.[dp1,dq1,dp2,dq2] (envelope theorem
//     => B is EXACT at the closest point; the E2 penalty force + main tangent build on it).
//
// Degeneracy gauge: denom = a*e - b*b = ||d1||^2 ||d2||^2 sin^2(theta_edge) >= 0 governs
// BOTH closest-point uniqueness AND normal conditioning. The PARALLEL refusal is
// denom < tau_par*a*e <=> sin^2(theta_edge) < tau_par; tau_par is conditioning-justified
// (TAU_PAR_DEFAULT = sin^2(15deg) ~ 0.0670), NOT the draft's 5deg. A zero-length edge
// (a collapsed facet edge after element removal — the collapse motivation) is DEGENERATE.
//
// Oracle (bit-for-bit): contact_prototypes/proto_e0_closest_point.py (28/28) +
// contact_prototypes/e0_cpp_check.cpp. See Ladruno_implementation/57_ladruno_edge_edge_contact_adr.md.

#ifndef LadrunoEdgeKernel_h
#define LadrunoEdgeKernel_h

#include <cmath>

namespace LadrunoEdgeKernel {

// status codes
enum Status { EE_OK = 0, EE_PARALLEL = -1, EE_DEGENERATE = -2 };

// gauges (match the oracle)
static const double TAU_PAR_DEFAULT = 0.0669872981077807;  // sin^2(15deg): parallel/conditioning gate
static const double TAU_LEN         = 1e-9;                 // squared-length zero-edge floor
static const double EE_MARGIN       = 1e-6;                 // parametric strict-interior margin delta

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
inline double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }

// ------------------------------------------------------ closest point result
struct ClosestResult {
    int    status;       // EE_OK / EE_PARALLEL / EE_DEGENERATE
    double s, t;         // closest parameters on seg1 (slave), seg2 (master), each in [0,1]
    double c1[3], c2[3]; // closest points  c1 = p1 + s*d1 (slave), c2 = p2 + t*d2 (master)
    double denom;        // a*e - b*b = ||d1||^2 ||d2||^2 sin^2(theta_edge)
    double a, e;         // squared edge lengths
    bool   interior;     // margin'd: s,t in [MARGIN, 1-MARGIN] (a true edge-edge crossing)
};

// closest point between seg1 = p1->q1 (slave) and seg2 = p2->q2 (master).
// Ericson RTCD 5.1.9 ClosestPtSegmentSegment + the fork parallel/degenerate refusal.
// On EE_PARALLEL the Ericson s=0 fallback (s,t,c1,c2) is STILL returned so the handler
// can run the clip cross-check before dropping the pair (no silent lost contact, Lens-C).
inline ClosestResult closestPtSegSeg(const double p1[3], const double q1[3],
                                     const double p2[3], const double q2[3],
                                     double tau_par = TAU_PAR_DEFAULT) {
    double d1[3], d2[3], r[3];
    for (int k = 0; k < 3; k++) { d1[k] = q1[k]-p1[k]; d2[k] = q2[k]-p2[k]; r[k] = p1[k]-p2[k]; }
    double a = dot3(d1, d1);   // ||seg1||^2
    double e = dot3(d2, d2);   // ||seg2||^2

    ClosestResult o;
    o.status = EE_OK; o.s = 0.0; o.t = 0.0; o.denom = 0.0; o.a = a; o.e = e; o.interior = false;
    for (int k = 0; k < 3; k++) { o.c1[k] = p1[k]; o.c2[k] = p2[k]; }

    // zero-length / coincident-node guard
    if (a <= TAU_LEN || e <= TAU_LEN) { o.status = EE_DEGENERATE; return o; }

    double c = dot3(d1, r);
    double b = dot3(d1, d2);
    double f = dot3(d2, r);
    double denom = a*e - b*b;          // >= 0
    o.denom = denom;

    // Ericson's EXACT branch ladder (conditional, NOT one pass)
    double s = (denom > 1e-300) ? clamp01((b*f - c*e) / denom) : 0.0;
    double t = (b*s + f) / e;
    if (t < 0.0)      { t = 0.0; s = clamp01(-c / a); }
    else if (t > 1.0) { t = 1.0; s = clamp01((b - c) / a); }

    o.s = s; o.t = t;
    for (int k = 0; k < 3; k++) { o.c1[k] = p1[k] + s*d1[k]; o.c2[k] = p2[k] + t*d2[k]; }
    o.interior = (s >= EE_MARGIN && s <= 1.0 - EE_MARGIN &&
                  t >= EE_MARGIN && t <= 1.0 - EE_MARGIN);

    // parallel / ill-conditioned-normal refusal: sin^2(theta_edge) < tau_par
    if (denom < tau_par * a * e) o.status = EE_PARALLEL;
    return o;
}

// common-perpendicular UNIT normal n = (d1 x d2)/||.||. returns false if degenerate
// (parallel/zero edges => ||d1 x d2|| ~ 0). d1,d2 are the edge direction vectors.
inline bool edgeNormal(const double d1[3], const double d2[3], double n[3]) {
    double m[3];
    cross3(d1, d2, m);
    double nrm = norm3(m);
    if (nrm < 1e-300) return false;
    for (int k = 0; k < 3; k++) n[k] = m[k] / nrm;
    return true;
}

// signed edge gap gN = w.n, w = c1 - c2 (slave - master). The sign is BODY-FIXED:
//   committedSign in {+1,-1} => applied verbatim (the shipped, interpenetration-safe path).
//   committedSign == 0 with signRef != 0 => orient n so n.signRef >= 0; *outSign returns
//     that sign to COMMIT once at first engagement (then pass it back as committedSign).
// (The self-referential n.w rule is intentionally NOT offered here — it masks penetration.)
// Writes the oriented normal into n[3]; returns gN. *outSign (optional) gets the chosen sign.
inline double edgeGap(const double c1[3], const double c2[3], const double d1[3],
                      const double d2[3], int committedSign,
                      const double *signRef, double n[3], int *outSign = 0) {
    if (!edgeNormal(d1, d2, n)) { if (outSign) *outSign = 0; return 0.0; }
    double w[3];
    for (int k = 0; k < 3; k++) w[k] = c1[k] - c2[k];
    int sign;
    if (committedSign == 1 || committedSign == -1) {
        sign = committedSign;
    } else if (signRef != 0) {
        sign = (dot3(n, signRef) >= 0.0) ? 1 : -1;
    } else {
        sign = 1;  // no anchor supplied: leave as the raw cross-product orientation
    }
    for (int k = 0; k < 3; k++) n[k] *= sign;
    if (outSign) *outSign = sign;
    return dot3(w, n);
}

// first variation dgN = B . [dp1,dq1,dp2,dq2]. Fills B[4][3] (rows = p1,q1,p2,q2).
// Exact at the closest point (envelope theorem: w _|_ tangents, w || n). The E2 penalty
// residual f = tN * B and main tangent K = epsN * B^T B build directly on this.
inline void bOperator(const double n[3], double s, double t, double B[4][3]) {
    for (int k = 0; k < 3; k++) {
        B[0][k] = (1.0 - s) * n[k];   // p1 (slave node 0)
        B[1][k] =        s  * n[k];   // q1 (slave node 1)
        B[2][k] = -(1.0 - t) * n[k];  // p2 (master node 0)
        B[3][k] =       -t  * n[k];   // q2 (master node 1)
    }
}

} // namespace LadrunoEdgeKernel

#endif
