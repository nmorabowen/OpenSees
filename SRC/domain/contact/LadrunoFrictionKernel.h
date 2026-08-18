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
    if (cap <= 0.0) {
        // ZERO cone radius while the caller's friction guard let the return map run:
        // the unified cone min(μN+c, τmax) is 0 (reachable via a τmax-only config —
        // μ=c=0, min(0,τmax)=0 — or a separated caller) ⇒ FREE SLIP: zero tangential
        // traction, the slip absorbs the whole tangential motion. The OLD branch
        // returned the RAW ELASTIC traction here ("byte-identity is the guard's job"),
        // which silently turned `-tauMax`-without-`-mu`/`-cohesion` into an UNBOUNDED
        // elastic bond (review 2026-07 HIGH-2); the command surface now also REFUSES
        // that config, so this branch is the defense-in-depth for future callers.
        for (int d = 0; d < 3; d++) { tFric[d] = 0.0; gpTtrial[d] = gTeff[d]; }
        return true;
    }
    if (norm <= cap) {                               // stick
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
    if (cap <= 0.0) {                                // FREE SLIP (zero cone): K_ss = 0
        // must mirror frictionReturnMap's cap≤0 free-slip branch (zero traction, slip
        // absorbs the motion ⇒ zero tangential stiffness) — the OLD stick tangent here
        // was consistent with the OLD raw-elastic force, i.e. a consistently-assembled
        // unbounded bond that Newton happily converged to (review 2026-07 HIGH-2).
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] = 0.0;
        return;
    }
    if (nrm <= cap) {                                // STICK: K_ss = kt·P_t
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

// ==========================================================================
// Ladruno ADR-85 T2 -- the 2D NTS scalar friction kernel (SS How/4). Added BESIDE the
// 3D functions above, SHARING frictionCap (the unified cone) -- no duplicated cone
// logic. The existing 3D frictionReturnMap/frictionTangentBlock are UNTOUCHED above
// this marker (byte-identity-gated, shipped).
//
// DECISION (measured, not a style preference -- see Ladruno_implementation/
// _adr85_t2_design.md S1, "THE PRE-T2 EXPERIMENT"): a 2D contact's tangent space is
// ONE-dimensional (t_hat = perp(n)), so the ADR-41 covariant-metric 3D return map
// degenerates trivially -- but evaluating the SHIPPED 3D kernel degenerately (a 3rd
// component pinned to 0) was measured, not assumed, against an independent scalar
// path over 200 000 randomized configs (34% parked astride the cone radius) with
// every stick/slip disagreement adjudicated in EXACT rational arithmetic:
//   (1) threshold correctness -- the degenerate 3D path takes the WRONG stick/slip
//       branch in 149/167 (89.2%) of adjudicated disagreements, vs 18/167 (10.8%)
//       for the scalar path, whose only error is the SINGLE unavoidable rounding of
//       fl(kt*(s-sp)); the 3D path adds avoidable sqrt(tx^2+ty^2)-vs-fabs and
//       normalization rounding on top;
//   (2) the dominant finding -- storing the slip as a GLOBAL-FRAME VECTOR (the 3D
//       kernel's gpT[3]/gT0[3]) causes CATASTROPHIC CANCELLATION in the slip
//       INCREMENT: median relative error 1.5e-2 (max 100%) in the |ds|/|s| ~
//       1e-16..1e-14 regime a converging Newton iteration actually lives in, because
//       the tangential origin and the current slip are each rounded on the GLOBAL
//       axes before their componentwise difference is taken. A tangential SCALAR
//       store (this file) is exactly ZERO relative error there (Sterbenz's lemma
//       applies cleanly to a single close subtraction; it does not reconstruct
//       cleanly across a vector's separately-rounded components) -- this is the
//       ADR-41 C3 displacement-not-position crux resurfacing in a new frame: the
//       slip must be KEPT as a scalar in the tangential frame, never differenced
//       componentwise as a global-frame vector (the BINDING implementation
//       consequence T2 wiring must honor);
//   (3) in a 1-D tangent space the slip-branch tangential tangent (Pt - n_hat(x)n_hat)
//       is IDENTICALLY ZERO; the degenerate 3D path reaches exact zero in only 28.9%
//       of slip cases (residual <= 4.5e-16*kt otherwise) -- so this kernel WRITES THE
//       LITERAL 0.0 rather than compute a tangent-plane projector that is already
//       known to collapse to nothing;
//   (4) leakage -- the degenerate 3D path leaks friction force into the NORMAL
//       direction up to 3.85e-8 relative (|tF.n|/|tF|), silently perturbing the
//       normal equilibrium the frictionless T1b lane already gates; the scalar path
//       is structurally incapable of a normal-direction component (there is none).
// Exactness and simplicity point the SAME way (no trade-off to arbitrate): scalar,
// ~20 lines, wins every measured axis.
//
// SEMANTICS mirror the 3D kernel EXACTLY otherwise (ADR-85 T2 work item 1):
//   * the cap<=0.0 FREE-SLIP branch (zero traction, slip absorbs the whole motion --
//     the frictionReturnMap HIGH-2 defense-in-depth, reproduced verbatim);
//   * the ALREADY-NEGATED sign convention: tFric = -tT, so the returned force OPPOSES
//     the slave's tangential motion (design-gate BLOCKER-1);
//   * spTrial (gpTtrial) a PURE function of COMMITTED state (set, never +=), so it is
//     idempotent across CDL's firstStep double getResidual (design-gate BLOCKER-2).
//
// std::fabs / a literal sgn(), NEVER sqrt(t*t) or a normalize (t/|t|): finding (1)'s
// mechanism is exactly sqrt-vs-fabs and normalization rounding; in 1-D there is only
// a SIGN to recover (the magnitude is already |tTtr|), so a normalize is not just
// slower, it is the reintroduced source of the measured error.
inline double sgn1D(double x) { return (x < 0.0) ? -1.0 : 1.0; }

// The 2D scalar return map (ADR-85 T2 work item 1). s/sp are the tangential relative-
// displacement and the committed plastic slip, BOTH scalars in the t_hat = perp(n)
// frame (never a 2-vector -- see the file-level note above); tFric/spTrial mirror
// frictionReturnMap's tFric[3]/gpTtrial[3] one-for-one, degree-of-freedom for
// degree-of-freedom. Returns true iff slipping (diagnostic, mirrors the 3D map).
inline bool returnMap1D(double s, double sp, double N, double kt, double mu,
                        double &tFric, double &spTrial,
                        double cohesion = 0.0, double tauMax = 0.0)
{
    const double tTtr = kt * (s - sp);        // trial elastic SCALAR tangential traction
    const double cap  = frictionCap(N, mu, cohesion, tauMax);   // the SHARED unified cone
    const double atr  = std::fabs(tTtr);      // NEVER sqrt(tTtr*tTtr) -- finding (1)
    if (cap <= 0.0) {
        // mirrors frictionReturnMap's cap<=0.0 free-slip defense-in-depth (review
        // 2026-07 HIGH-2): zero cone radius => zero tangential traction, the slip
        // absorbs the whole tangential motion (never the raw elastic trial).
        tFric = 0.0; spTrial = s;
        return true;
    }
    if (atr <= cap) {                          // STICK
        tFric = -tTtr; spTrial = sp;
        return false;
    }
    const double sg   = sgn1D(tTtr);           // direction ONLY -- never a normalize
    const double dlam = (atr - cap) / kt;      // atr > cap >= 0 => safe
    tFric   = -cap * sg;
    spTrial = sp + dlam * sg;
    return true;
}

// The 2D scalar friction tangent (ADR-85 T2 work item 1/3; the tangentBlock1D name in
// the phase brief). Kss = d(tFric)/d(s) restricted to the tangential direction (the
// STICK/SLIP block); dTN = the SCALAR pressure-coupling coefficient of ADR-85 How/5,
// -(dcap/dN)*kn -- the caller (LadrunoContactFE::addFrictionTang2D) scatters it as the
// scalar-times-n_hat outer term dTN*(t_hat (x) n) the ADR calls for -- dTN ALREADY carries
// kn, so the caller must NOT reapply it -- since in 2D the
// full cross term is a 2x2 outer product of the ONE tangential direction with the
// normal direction, not a 3x3 block. consistent=false (the default) leaves dTN at 0.0
// and the caller therefore stays SYMMETRIC (design-gate Q2: a symmetric SOE silently
// drops an asymmetric contribution, so a non-symmetric default would corrupt the solve).
inline void tangentBlock1D(double s, double sp, double N, double kn, double kt, double mu,
                           bool consistent, double &Kss, double &dTN,
                           double cohesion = 0.0, double tauMax = 0.0)
{
    const double tTtr   = kt * (s - sp);
    const double capC   = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    const bool   capped = (tauMax > 0.0 && tauMax < capC);   // tauMax cap binds => dcap/dN = 0
    const double cap    = capped ? tauMax : capC;
    const double atr    = std::fabs(tTtr);
    dTN = 0.0;
    if (cap <= 0.0) {                          // FREE SLIP: mirrors frictionTangentBlock's Kss=0
        Kss = 0.0;
        return;
    }
    if (atr <= cap) {                          // STICK: exactly kt -- no projector exists in 1-D
        Kss = kt;
        return;
    }
    // SLIP: the 1-D slip-branch tangential tangent IS the (Pt - n_hat(x)n_hat) block, and in a
    // 1-D tangent space that projector is IDENTICALLY ZERO (finding (3) above) -- write the
    // literal, never compute it (the design brief's explicit requirement).
    Kss = 0.0;
    if (consistent) {
        const double dCap_dN = (capped || N <= 0.0) ? 0.0 : mu;   // the min()-selected coupling
        const double sg = sgn1D(tTtr);
        dTN = -dCap_dN * kn * sg;
    }
}

} // namespace LadrunoFrictionKernel

#endif
