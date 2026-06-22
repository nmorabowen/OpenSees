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

// ADR-39 P2b — LadrunoContactKernel: header-only, OpenSees-free pure-function
// node-to-segment (NTS) penalty contact math. Mirrors LadrunoJ2Kernel.h so both
// OPS_Analysis (the FE adapter) and OPS_Domain (the surface) can include it with no
// link inversion. NO OpenSees types — only raw double[] + <cmath>.
//
// Validated oracle-first in Ladruno_implementation/contact_prototypes/proto_p2b_nts.py
// (7/7: projection, winding-flip-immune normal, penalty, B-operator, self-equilibrium,
// symmetric kn*BᵀB + ∂n/∂u tangent vs FD). 3D, master segment = tri-3 or quad-4.
//
// Conventions (match shipped P2a, verified): penetration ⇔ gap < 0; traction
// tn = kn*<−gap>₊ (Macaulay); getResidual returns +tn over the gap operator B so a
// penetrating slave is pushed back along +n (restoring, never attract). The assembled
// tangent is K = −∂r/∂q = +kn*BᵀB (main term; ∂n/∂u block is P2b-2). Outward normal is
// DERIVED + oriented from a reference point, NEVER trusted from node winding
// (design-gate BLOCKER-1).
//
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p2b_design.md.

#ifndef LadrunoContactKernel_h
#define LadrunoContactKernel_h

#include <cmath>

namespace LadrunoContactKernel {

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

// penalty traction tn = kn*<−gap>₊ . Internal Macaulay clamp (gate KMF-3): returns
// 0 for gap>=0 so a stray call can never produce a non-physical (adhesive) traction.
inline double traction(double kn, double gap) { return (gap < 0.0) ? kn * (-gap) : 0.0; }

// ------------------------------------------------ friction (ADR-39 P3, Coulomb)
// tangential component of v w.r.t. the unit normal n: v − (v·n) n. Used to map a
// relative-position vector onto the contact tangent plane.
inline void tangentPart(const double v[3], const double n[3], double out[3]) {
    double vn = dot3(v, n);
    for (int d = 0; d < 3; d++) out[d] = v[d] - vn * n[d];
}

// Clean Coulomb DIRECT return map (NOT the ASDimplex "apparent damage" form, which
// carries the quarantined shadowed-dE / ddata-OOB bugs). Structurally 1D plasticity
// on the tangential traction with a pressure-dependent (Coulomb) yield radius μN.
// Validated oracle-first in proto_p3_friction.py (6/6).
//
//   tT*  = kt·(gTeff − gpT)        trial elastic tangential traction (tangent plane)
//   cap  = μN                       Coulomb cone radius  (N = kn·<−gap>₊ ≥ 0)
//   STICK (‖tT*‖ ≤ cap): tT = tT*,  gpT unchanged
//   SLIP  (‖tT*‖ > cap): n̂ = tT*/‖tT*‖, dλ = (‖tT*‖−cap)/kt,
//                        tT = cap·n̂,  gpT += dλ·n̂
//
// SIGN (design-gate BLOCKER-1): the trial direction n̂* FOLLOWS the slave's tangential
// MOTION, so assembling +tT would ACCELERATE the slave (energy injection,
// a=g(sinθ+μcosθ)). The kernel therefore returns the FORCE the contact APPLIES to the
// slave, ALREADY NEGATED: tFric = −tT. The FE then mirrors the normal block exactly
// (resid_slave += tFric, resid_master_i += −N_i·tFric) so friction OPPOSES the motion
// (incline ⇒ a=g(sinθ−μcosθ)). The negation lives in this ONE auditable place.
//
// gpTtrial is a PURE function of committed state (gpT + dλ·n̂, set not +=), so it is
// idempotent across CDL's firstStep double getResidual (design-gate BLOCKER-2). The
// near-zero guard is PHYSICALLY scaled: slip requires ‖tT*‖ > cap > 0, so n̂ is always
// normalizable in the slip branch (no denormal 0/0); cap ≤ 0 (separated or μ=0) ⇒
// stick with the elastic traction. The IMPL-EX multiplier-extrapolation variant is
// retained only in the oracle for the P3.5 implicit leg (explicit discards the
// tangent, so IMPL-EX buys nothing here while costing a one-step onset overshoot).
// returns true if slipping (diagnostic), false if stick/inactive.
inline bool frictionReturnMap(const double gTeff[3], const double gpT[3],
                              double N, double kt, double mu,
                              double tFric[3], double gpTtrial[3]) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double cap  = mu * N;
    double norm = norm3(tTtr);
    if (cap <= 0.0 || norm <= cap) {                 // stick (or inactive)
        for (int d = 0; d < 3; d++) { tFric[d] = -tTtr[d]; gpTtrial[d] = gpT[d]; }
        return false;
    }
    double inv  = 1.0 / norm;                        // norm > cap > 0 ⇒ safe
    double dlam = (norm - cap) / kt;
    for (int d = 0; d < 3; d++) {
        double nh   = tTtr[d] * inv;
        tFric[d]    = -cap * nh;
        gpTtrial[d] = gpT[d] + dlam * nh;
    }
    return true;
}

} // namespace LadrunoContactKernel

#endif
