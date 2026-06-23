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

// ADR-39 P2b — LadrunoContactKernel: header-only, OpenSees-free pure-function NTS
// penalty NORMAL-LAW + COULOMB FRICTION math (traction, return map, consistent tangent).
// The closest-point PROJECTION geometry was EXTRACTED to LadrunoContactProjection.h
// (ADR-41 A2, capstone Track A2) — ONE projection for both the NTS path and the mortar
// Gauss-point loop. Mirrors LadrunoJ2Kernel.h so both OPS_Analysis (the FE adapter) and
// OPS_Domain (the surface) can include it with no link inversion. NO OpenSees types —
// only raw double[] + <cmath>.
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
#include "LadrunoContactProjection.h"   // ADR-41 A2: projection geometry extracted here

namespace LadrunoContactKernel {

// Geometry (shape / bounded projection / oriented normal / evalSegment / vec3 helpers)
// now lives in LadrunoContactProjection.h (ADR-41 A2) — ONE projection for the NTS path
// AND the mortar GP loop. The normal-law + friction below need only the vec3 helpers,
// imported here so the call sites read unchanged:
using LadrunoContactProjection::dot3;
using LadrunoContactProjection::norm3;

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

// ADR-39 P3.5 — friction tangent SLAVE block K_ss = D_TT·P_t + d_TN⊗n  (3×3), for the
// assembled IMPLICIT tangent K_fric = Gᵀ K_ss G (G=[I|−N_i I], scattered by the FE).
// Validated in proto_p35_implicit_tangent.py (slip==FD non-sym; stick==kt·P_t; the
// committed-N/symmetric variant drops d_TN). consistent=false ⇒ DROP d_TN⊗n (the
// SYMMETRIC default — solver-safe on any system; a symmetric SOE silently drops the
// lower triangle, so a non-sym default would corrupt the solve, design-gate Q2); true
// ⇒ include it (full consistent tangent: quadratic, NON-symmetric, needs FullGeneral/
// UmfPack/BandGen). gTeff = slip from engagement; gpT = committed plastic slip; n = unit
// normal; N = kn·<−gap>₊. Stick K_ss = kt·P_t; slip uses the idempotent D_TT·P_t = D_TT.
inline void frictionTangentBlock(const double gTeff[3], const double gpT[3],
                                 const double n[3], double N, double kn, double kt,
                                 double mu, bool consistent, double Kss[3][3]) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double cap = mu * N, nrm = norm3(tTtr);
    double Pt[3][3];                                 // tangent-plane projector I − n⊗n
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) Pt[i][j] = (i == j ? 1.0 : 0.0) - n[i]*n[j];
    if (cap <= 0.0 || nrm <= cap) {                  // STICK: K_ss = kt·P_t
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] = kt * Pt[i][j];
        return;
    }
    double inv = 1.0 / nrm, nh[3];                   // SLIP
    for (int d = 0; d < 3; d++) nh[d] = tTtr[d] * inv;
    double s = cap * kt * inv;                        // μN·kt/‖tT*‖
    for (int i = 0; i < 3; i++)                       // (μN·kt/‖tT*‖)(P_t − n̂⊗n̂)
        for (int j = 0; j < 3; j++) Kss[i][j] = s * (Pt[i][j] - nh[i]*nh[j]);
    if (consistent) {                                // + d_TN⊗n , d_TN = −μ·kn·n̂
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] += (-mu * kn * nh[i]) * n[j];
    }
}

} // namespace LadrunoContactKernel

#endif
