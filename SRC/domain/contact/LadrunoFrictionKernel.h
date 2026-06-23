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

// ADR-39 P3/P3.5 / ADR-41 A1 — LadrunoFrictionKernel: header-only, OpenSees-free
// Coulomb/Tresca contact-friction math (return map + consistent tangent). EXTRACTED
// from LadrunoContactKernel.h (PRs #360/#361) so ONE friction return map serves both
// the NTS penalty path AND the ADR-41 mortar path — no duplicated constitutive code
// (the staleness trap). NO OpenSees types — only raw double[] + <cmath>; standalone
// (its own vec helpers) so it is numpy-oracle-testable without building OpenSees.
// Mirrors the LadrunoJ2Kernel.h discipline.
//
// Placement note (capstone Track A1): ADR-41 originally sketched this at
// SRC/material/nD/ "beside LadrunoJ2Kernel.h", but the friction kernel is
// contact-specific (consumes the contact normal pressure N) and is consumed by
// OPS_Analysis (LadrunoContactFE) + OPS_Domain — so it lives with the other contact
// kernels in SRC/domain/contact/, whose dir is already on the consumers' include path
// (proven by the A2 projection extraction). Co-located: projection / friction / normal-law.
//
// UNIFIED CONE (ADR-41 / capstone, sharpened vs Abaqus TG §5.2.3):
//   cap = min( μ·N + c ,  τmax )         N = kn·<−gap>₊ ≥ 0 (contact pressure), c = cohesion
//   - pure Coulomb:        c=0, τmax≤0 (no cap)            → cap = μN
//   - Coulomb + cohesion:  c>0                              → cap = μN + c
//   - Coulomb + cap:       τmax>0 finite                    → cap = min(μN+c, τmax)
//   - Tresca (constant):   μ=0, c=τ_Tresca, τmax≤0          → cap = c
// The pressure-coupling tangent block Csl (∂t_T/∂g_N) is ACTIVE only where cap depends
// on N — i.e. μ on the Coulomb branch, 0 when the τmax cap binds (pressure-independent)
// or when separated (N≤0). This is the min()-selected non-symmetric coupling.
//
// BIT-FOR-BIT with the shipped P3/P3.5 NTS friction: the two extra params
// (cohesion, tauMax) default to (0, 0) and the N>0 gate reproduces μN exactly, so the
// existing call sites (which pass neither) are unchanged. Validated oracle-first in
// contact_prototypes/proto_a1_friction.py (Css/Csl vs FD to 1e-6; both min() branches;
// Tresca/cohesion; bit-for-bit vs the pre-extraction map).
//
// See Ladruno_implementation/{39_ladruno_contact_domain,41_ladruno_mortar_alm_contact,
// 48_ladruno_contact_capstone}_adr.md.

#ifndef LadrunoFrictionKernel_h
#define LadrunoFrictionKernel_h

#include <cmath>

namespace LadrunoFrictionKernel {

// standalone vec helpers (kept local so the kernel needs only <cmath>)
inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }

// tangential component of v w.r.t. the unit normal n: v − (v·n) n. Maps a relative-
// position vector onto the contact tangent plane (the slip measure).
inline void tangentPart(const double v[3], const double n[3], double out[3]) {
    double vn = dot3(v, n);
    for (int d = 0; d < 3; d++) out[d] = v[d] - vn * n[d];
}

// Coulomb/Tresca/cohesion critical shear cap (the unified cone), pressure N = kn·<−gap>₊.
//   capC = μN + c   (only in contact, N>0; separated ⇒ 0) ;  cap = min(capC, τmax)
inline double frictionCap(double N, double mu, double cohesion, double tauMax) {
    double capC = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    return (tauMax > 0.0 && tauMax < capC) ? tauMax : capC;
}

// Clean Coulomb/Tresca DIRECT return map (NOT the ASDimplex "apparent damage" form).
// Structurally 1D plasticity on the tangential traction with the unified yield radius
// cap = min(μN+c, τmax).
//
//   tT*  = kt·(gTeff − gpT)        trial elastic tangential traction (tangent plane)
//   STICK (‖tT*‖ ≤ cap): tT = tT*,  gpT unchanged
//   SLIP  (‖tT*‖ > cap): n̂ = tT*/‖tT*‖, dλ = (‖tT*‖−cap)/kt,  tT = cap·n̂,  gpT += dλ·n̂
//
// SIGN (design-gate BLOCKER-1): n̂* FOLLOWS the slave's tangential MOTION, so assembling
// +tT would ACCELERATE the slave (energy injection). The kernel returns the force the
// contact APPLIES to the slave, ALREADY NEGATED: tFric = −tT, so friction OPPOSES the
// motion. gpTtrial is a PURE function of committed state (set, not +=) ⇒ idempotent
// across CDL's firstStep double getResidual (design-gate BLOCKER-2). The slip-branch n̂
// is always normalizable (‖tT*‖ > cap ≥ 0). returns true if slipping (diagnostic).
inline bool frictionReturnMap(const double gTeff[3], const double gpT[3],
                              double N, double kt, double mu,
                              double tFric[3], double gpTtrial[3],
                              double cohesion = 0.0, double tauMax = 0.0) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double cap  = frictionCap(N, mu, cohesion, tauMax);
    double norm = norm3(tTtr);
    if (cap <= 0.0 || norm <= cap) {                 // stick (or inactive)
        for (int d = 0; d < 3; d++) { tFric[d] = -tTtr[d]; gpTtrial[d] = gpT[d]; }
        return false;
    }
    double inv  = 1.0 / norm;                        // norm > cap ≥ 0 ⇒ safe
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
// consistent=false ⇒ DROP d_TN⊗n (the SYMMETRIC default — solver-safe on any system;
// a symmetric SOE silently drops the lower triangle, so a non-sym default would corrupt
// the solve, design-gate Q2); true ⇒ include it (full consistent tangent: quadratic,
// NON-symmetric; needs FullGeneral/UmfPack/BandGen). The d_TN (Csl) coupling = −(∂cap/∂N)
// ·kn·n̂ with ∂cap/∂N = μ on the Coulomb branch, 0 when the τmax cap binds or separated
// (the min()-selected coupling). Stick K_ss = kt·P_t; slip uses the idempotent D_TT·P_t.
inline void frictionTangentBlock(const double gTeff[3], const double gpT[3],
                                 const double n[3], double N, double kn, double kt,
                                 double mu, bool consistent, double Kss[3][3],
                                 double cohesion = 0.0, double tauMax = 0.0) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double capC   = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    bool   capped = (tauMax > 0.0 && tauMax < capC);   // τmax cap binds ⇒ ∂cap/∂N = 0
    double cap    = capped ? tauMax : capC;
    double nrm = norm3(tTtr);
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
    double s = cap * kt * inv;                        // cap·kt/‖tT*‖
    for (int i = 0; i < 3; i++)                       // (cap·kt/‖tT*‖)(P_t − n̂⊗n̂)
        for (int j = 0; j < 3; j++) Kss[i][j] = s * (Pt[i][j] - nh[i]*nh[j]);
    if (consistent) {                                // + d_TN⊗n , d_TN = −(∂cap/∂N)·kn·n̂
        double dCap_dN = (capped || N <= 0.0) ? 0.0 : mu;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] += (-dCap_dN * kn * nh[i]) * n[j];
    }
}

} // namespace LadrunoFrictionKernel

#endif
