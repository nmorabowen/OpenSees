/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// LadrunoConcrete3DKernel — the PURE numerical core of LadrunoConcrete3D, a CDPM2-grade
// solid-concrete plastic-damage model. Plain doubles + <cmath>, NO OpenSees dependency, so the
// SAME verified core can be:
//   * unit-tested standalone (g++) via tests/_testbed/concrete3d_kernel_check.cpp (surface +
//     hardening identities; the full oracle-numeric-dump diff is the P1 build-PR deliverable),
//   * delegated to by the small-strain  nDMaterial LadrunoConcrete3D (classTag 33017), and
//   * lifted to finite strain by  nDMaterial LogStrain  (LogStrainNDMaterial 33010) feeding
//     LadrunoBrick -geom finite  (isotropic plastic-damage is objective under large rotation).
//
// Model (ADR Ladruno_implementation/31_ladruno_concrete3d_adr.md):
//   * Plasticity: effective-stress Menetrey-Willam 3-invariant surface (friction m0, eccentricity
//     e, Willam-Warnke Lode r(theta,e)); NON-ASSOCIATED flow + dilatancy Df; hardening qh1/qh2
//     driven by a confinement-aware ductility measure x(sigma); semi-implicit return map with a
//     dedicated APEX/Lode-corner sub-algorithm.
//   * Damage: dual scalar omega_t/omega_c on the spectral tension/compression split; crack-band
//     regularized in BOTH Gf and Gc; unilateral (crack-closure) recovery.
//   * Robustness (one tangent-agnostic kernel): Tier-1 implicit-accurate return; Tier-2 IMPL-EX
//     (freeze plastic state + damage => SPD secant); Tier-3 explicit (no tangent). Duvaut-Lions
//     viscous eta at the PLASTIC level.
//
// STATUS: P0/P1 COMPLETE (surface + return map + analytic consistent tangent). The CDPM2 yield
// surface (yieldF, Eq.18 with the [1-qh1] cap), hardening laws (qh1Of/qh2Of Eq.30-31, ductilityXh
// Eq.33-36), the SEMI-IMPLICIT return map (returnMapPrincipal perfect-plastic 3-Newton w/ analytic
// 3x3 Jacobian; returnMapHardening 4-Newton w/ analytic 4x4 Jacobian; hydrostatic-tension apex;
// honest convergence), the spectral full-tensor lift (returnMapTensor), and the NON-SYMMETRIC
// analytic consistent tangent (consistentTangent) are IMPLEMENTED and verified byte-for-byte
// against the numpy oracle (tests/_testbed/concrete3d_ref.py) via the g++ oracle-numeric-dump diff
// (tests/_testbed/concrete3d_kernel_check.cpp): perfect-plastic stress+tangent to machine precision,
// hardening to the oracle's numerical-Jacobian reference floor (~1e-8). KNOWN GAP (inherited from
// the oracle, handoff §6): the low-kappa_p tiny-surface / off-meridian first-yield APEX regime needs
// sub-stepping + a dedicated apex sub-algorithm.
// P2/P3a/P3b DAMAGE COMPLETE in the kernel: the dual scalar omega_t/omega_c update (damagedUpdate,
// Eq.1 spectral split + automatic unilateral re-split, where the stress PEAK/softening comes from)
// AND the P3b ANALYTIC dual-projector DAMAGED consistent tangent (damagedTangent: D_dam:C_eff -
// sig_t(x)d_omega_t - sig_c(x)d_omega_c) — returnMap returns this damaged tangent when doTangent.
// Both verified against the oracle (NOMINAL stress ~1e-14; damaged tangent == a numerical central
// diff of the SAME C++ stress ~1e-7 AND == the oracle damaged_tangent_analytic ~1e-7). STILL TO COME
// (later phases): robustness tiers (P3 IMPL-EX), finite-strain out-contract (P4), confined-fiber
// condensation (§4.6), the cyclic beta_c temper (P2f). NO OpenSees nDMaterial wrapper yet (classTag
// 33017 reserved).
//
// Tensor convention: symmetric tensors stored as 6 TENSOR components in the canonical OpenSees
// order {00,11,22,01,12,02} (off-diagonal slots hold TRUE tensor components, not engineering —
// the engineering<->tensor conversion is the caller's responsibility, as in LadrunoJ2.cpp).
// Sign convention: COMPRESSION NEGATIVE; fc, ft are positive magnitudes.

#ifndef LadrunoConcrete3DKernel_h
#define LadrunoConcrete3DKernel_h

#include <cmath>

namespace Ladruno {
namespace Concrete3D {

// ---------------------------------------------------------------------------
// Material parameters (POD; the nDMaterial wrapper fills + serializes these)
// ---------------------------------------------------------------------------
struct Params {
    // elasticity
    double E = 0.0;
    double nu = 0.0;
    double rho = 0.0;
    // strengths (POSITIVE magnitudes)
    double fc = 0.0;
    double ft = 0.0;
    // surface shape
    double e = 0.0;     // eccentricity (deviatoric out-of-roundness); see eccentricityFromKupfer
    double m0 = 0.0;    // friction; = m0Of(fc,ft,e)
    // fracture energies (crack-band regularized)
    double Gf = 0.0;    // tensile
    double Gc = 0.0;    // compressive (WEAKEST-calibrated knob — ADR §6)
    // flow / ductility (P1)
    double Df = 0.0;    // dilatancy (non-associated flow)
    // damage (P2) softening ductility (Eq.56: x_s = 1 + (As-1) Rs)
    double As = 2.0;    // confinement-ductility amplitude (=x_s in uniaxial compression)
    // hardening (P1, CDPM2 Eqs. 30-36)
    double qh0 = 0.3;   // initial yield fraction (qh1 at kp=0)
    double Hp = 0.5;    // hardening modulus (qh1 slope at kp=1^-, qh2 slope for kp>1)
    double Ah = 0.08;   // ductility-measure params (Eq.33; literature defaults, recalibrate)
    double Bh = 0.003;
    double Ch = 2.0;
    double Dh = 1.0e-6;
    // rate / robustness
    double eta = 0.0;            // Duvaut-Lions viscosity (0 => inviscid, byte-identical)
    bool   implex = false;       // Tier-2 (IMPL-EX)
    double implexRmax = 2.0;     // IMPL-EX extrapolation time-ratio cap (review ALG-2; r=dt/dt_n clamped [0,rmax])
    // regularization
    double lch = 1.0;            // parent-element characteristic length
    double lch_ref = 1.0;        // reference length of the input Gf/Gc
};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
static const double SQRT3 = 1.7320508075688772;
static const double SQRT6 = 2.449489742783178;
static const double SQRT1_5 = 1.224744871391589;

// ---------------------------------------------------------------------------
// Stress invariants. sig = {s00,s11,s22,s01,s12,s02} TENSOR components.
// Returns xi (=I1/sqrt3), rho (=sqrt(2 J2)), theta (Lode angle in [0,pi/3]).
// ---------------------------------------------------------------------------
inline void invariants(const double sig[6], double& xi, double& rho, double& theta)
{
    const double I1 = sig[0] + sig[1] + sig[2];
    const double p = I1 / 3.0;
    const double d0 = sig[0] - p, d1 = sig[1] - p, d2 = sig[2] - p;
    const double s01 = sig[3], s12 = sig[4], s02 = sig[5];
    const double J2 = 0.5 * (d0 * d0 + d1 * d1 + d2 * d2) + s01 * s01 + s12 * s12 + s02 * s02;
    const double J3 = d0 * (d1 * d2 - s12 * s12)
                    - s01 * (s01 * d2 - s12 * s02)
                    + s02 * (s01 * s12 - d1 * s02);
    xi = I1 / SQRT3;
    rho = (J2 > 0.0) ? std::sqrt(2.0 * J2) : 0.0;
    if (J2 <= 1.0e-300) {
        theta = 0.0;
    } else {
        double c3 = (3.0 * SQRT3 / 2.0) * J3 / std::pow(J2, 1.5);
        if (c3 > 1.0) c3 = 1.0; else if (c3 < -1.0) c3 = -1.0;
        theta = std::acos(c3) / 3.0;
    }
}

// ---------------------------------------------------------------------------
// Willam-Warnke elliptic Lode function r(theta,e), e in (0.5,1].
// r(0)=1/e (tensile meridian), r(pi/3)=1 (compressive meridian) => r_c/r_t = e.
// Convex for e in [0.5,1] (Willam-Warnke 1975). NOTE: a naive per-sextant 1/r second-derivative
// ("g+g''") polar test reports false violations near theta~52deg for larger e and is NOT the
// correct convexity criterion for the elliptic interpolation — do not "fix" a non-bug.
// ---------------------------------------------------------------------------
inline double lodeR(double theta, double e)
{
    const double ct = std::cos(theta);
    const double oneMe2 = 1.0 - e * e;
    const double num = 4.0 * oneMe2 * ct * ct + (2.0 * e - 1.0) * (2.0 * e - 1.0);
    double rad = 4.0 * oneMe2 * ct * ct + 5.0 * e * e - 4.0 * e;
    if (rad < 0.0) rad = 0.0;
    const double den = 2.0 * oneMe2 * ct + (2.0 * e - 1.0) * std::sqrt(rad);
    return num / den;
}

inline double m0Of(double fc, double ft, double e)
{
    return 3.0 * (fc * fc - ft * ft) / (fc * ft) * e / (e + 1.0);
}

// elastic moduli from (E, nu) — the oracle's make_material convention.
inline double bulkK(const Params& mp)  { return mp.E / (3.0 * (1.0 - 2.0 * mp.nu)); }
inline double shearG(const Params& mp) { return mp.E / (2.0 * (1.0 + mp.nu)); }

// ---------------------------------------------------------------------------
// CDPM2 yield function f_p — Grassl et al. 2013 IJSS Eq. (18). qh1=qh2=1 => FAILURE surface
// Eq. (21) = Menetrey-Willam 1995. Note xi/(sqrt3 fc) = (I1/3)/fc = sigma_V/fc (Eq.12).
// ---------------------------------------------------------------------------
inline double yieldF(const double sig[6], const Params& mp, double qh1 = 1.0, double qh2 = 1.0)
{
    double xi, rho, theta;
    invariants(sig, xi, rho, theta);
    const double r = lodeR(theta, mp.e);
    const double sigV_fc = xi / (SQRT3 * mp.fc);                 // sigma_V/fc
    const double AV = rho / (SQRT6 * mp.fc) + sigV_fc;           // hardening-cap base
    const double RR = rho * r / (SQRT6 * mp.fc) + sigV_fc;       // m0-friction bracket (Lode r)
    const double quad = SQRT1_5 * rho / mp.fc;                   // sqrt(3/2) rho/fc
    const double cap = (1.0 - qh1) * AV * AV + quad;
    return cap * cap + mp.m0 * qh1 * qh1 * qh2 * RR - (qh1 * qh1) * (qh2 * qh2);  // Eq.(18)
}

// CDPM2 hardening building blocks (Grassl et al. 2013 Eqs. 30-36). qh1: qh0->1 over kp in [0,1]
// (Eq.30); qh2: 1 then 1+Hp(kp-1) (Eq.31); ductility xh (Eq.33) with Rh=-sigV/fc-1/3 (Eq.34) =>
// more ductile under compression. (Used by the P1 hardening return map; mirrors the numpy oracle.)
inline double qh1Of(double kp, double qh0, double Hp)
{
    if (kp < 1.0)
        return qh0 + (1.0 - qh0) * (kp * kp * kp - 3.0 * kp * kp + 3.0 * kp)
            - Hp * (kp * kp * kp - 3.0 * kp * kp + 2.0 * kp);
    return 1.0;
}
inline double qh2Of(double kp, double Hp) { return kp < 1.0 ? 1.0 : 1.0 + Hp * (kp - 1.0); }
inline double ductilityXh(double sigV, double fc, double Ah, double Bh, double Ch, double Dh)
{
    const double Rh = -sigV / fc - 1.0 / 3.0;                    // Eq.(34)
    if (Rh >= 0.0)
        return Ah - (Ah - Bh) * std::exp(-Rh / Ch);             // Eq.(33) upper
    const double Eh = Bh - Dh;                                   // Eq.(35)
    const double Fh = (Bh - Dh) * Ch / (Ah - Bh);               // Eq.(36)
    return Eh * std::exp(Rh / Fh) + Dh;                          // Eq.(33) lower
}

// ---------------------------------------------------------------------------
// MW surface value directly from invariants (theta enters via the FROZEN Lode r).
//   yfInv     : failure surface Eq.21 (qh1=qh2=1) — used by the perfect-plastic map.
//   yfInvHard : Eq.18 with hardening qh1(kp), qh2(kp) — used by the hardening map.
// (Mirror the oracle's _yf_inv / _yf_inv_hard byte-for-byte.)
// ---------------------------------------------------------------------------
inline double yfInv(double xi, double rho, double r, const Params& mp)
{
    return 1.5 * rho * rho / (mp.fc * mp.fc)
         + mp.m0 * (rho * r / (SQRT6 * mp.fc) + xi / (SQRT3 * mp.fc)) - 1.0;
}

inline double yfInvHard(double xi, double rho, double r, double kp, const Params& mp)
{
    const double q1 = qh1Of(kp, mp.qh0, mp.Hp);
    const double q2 = qh2Of(kp, mp.Hp);
    const double sigV_fc = xi / (SQRT3 * mp.fc);
    const double AV = rho / (SQRT6 * mp.fc) + sigV_fc;
    const double RR = rho * r / (SQRT6 * mp.fc) + sigV_fc;
    const double quad = SQRT1_5 * rho / mp.fc;
    const double cap = (1.0 - q1) * AV * AV + quad;
    return cap * cap + mp.m0 * q1 * q1 * q2 * RR - (q1 * q1) * (q2 * q2);
}

// kp-derivatives of the hardening functions (kp<1 branch; 0 for kp>=1). Eq.30-31.
inline double dqh1OfdKp(double kp, double qh0, double Hp)
{
    if (kp >= 1.0) return 0.0;
    return (1.0 - qh0) * (3.0 * kp * kp - 6.0 * kp + 3.0)
         - Hp * (3.0 * kp * kp - 6.0 * kp + 2.0);
}
inline double dqh2OfdKp(double kp, double Hp) { return kp < 1.0 ? 0.0 : Hp; }

// d(ductility xh)/d(sigV) — analytic, piecewise across Rh=0 (Eq.33-36). Rh=-sigV/fc-1/3.
inline double dDuctilityXhdSigV(double sigV, double fc, double Ah, double Bh, double Ch, double Dh)
{
    const double Rh = -sigV / fc - 1.0 / 3.0;
    const double dRh = -1.0 / fc;
    if (Rh >= 0.0)
        return (Ah - Bh) / Ch * std::exp(-Rh / Ch) * dRh;            // d/dRh of upper branch * dRh
    const double Eh = Bh - Dh;
    const double Fh = (Bh - Dh) * Ch / (Ah - Bh);
    return Eh / Fh * std::exp(Rh / Fh) * dRh;                        // lower branch
}

// ---------------------------------------------------------------------------
// Convenience: eccentricity that makes the equibiaxial strength hit fcc/fc=target (ADR 4.1b).
// Inline bisection mirror of the oracle's eccentricity_from_kupfer (+ equibiaxial_strength).
// e is CLAMPED to the open convexity band (0.5, 1] — a user-supplied e<=0.5 is rejected/clamped
// by the wrapper; this routine searches strictly inside (0.5, 1).
// ---------------------------------------------------------------------------
inline double equibiaxialStrength(double fc, double ft, double e)
{
    // |sigma| at biaxial compression sigma=(-b,-b,0) crossing the failure surface (f increases w/ b).
    double lo = 0.5 * fc, hi = 3.0 * fc;
    Params mp; mp.fc = fc; mp.ft = ft; mp.e = e; mp.m0 = m0Of(fc, ft, e);
    auto fOfB = [&](double b) { double s[6] = {-b, -b, 0, 0, 0, 0}; return yieldF(s, mp); };
    double flo = fOfB(lo);
    for (int it = 0; it < 200; ++it) {
        double mid = 0.5 * (lo + hi), fm = fOfB(mid);
        if (flo * fm <= 0.0) hi = mid; else { lo = mid; flo = fm; }
        if (hi - lo < 1.0e-12 * fc) break;
    }
    return 0.5 * (lo + hi);
}

inline double eccentricityFromKupfer(double fc, double ft, double targetFccRatio = 1.16)
{
    auto g = [&](double e) { return equibiaxialStrength(fc, ft, e) / fc - targetFccRatio; };
    double lo = 0.5 + 1.0e-6, hi = 1.0 - 1.0e-9;
    double glo = g(lo);
    for (int it = 0; it < 200; ++it) {
        double mid = 0.5 * (lo + hi), gm = g(mid);
        if (glo * gm <= 0.0) hi = mid; else { lo = mid; glo = gm; }
        if (hi - lo < 1.0e-12) break;
    }
    return 0.5 * (lo + hi);
}

// ===========================================================================
// P2 DAMAGE kinematics (CDPM2 §2.3, Grassl et al. 2013) — ports the oracle verbatim.
//   equivStrainGeneral  Eq.37 : ε̃ (== σ̄/E uniaxial tension, == ε0 on the failure surface)
//   alphaCompression    Eq.46 : 0 (tension) .. 1 (compression)
//   damageDrivers              : ε̃, α_c, and the softening ductility x_s (Eq.56-57)
//   solveOmegaBracketed        : the implicit (1-ω)D = f·exp(-(kd1+ω·kd2)/eps_f) root, bisection-
//                                safeguarded so a non-monotone F never clamp-stalls to 0 (PR #261)
// sig_pr = 3 EFFECTIVE principal stresses (the damage drivers are frame-invariant).
// ===========================================================================
inline double equivStrainGeneral(const double sig_pr[3], const Params& mp)
{
    const double eps0 = mp.ft / mp.E;
    const double sv[6] = { sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0 };
    double xi, rho, theta; invariants(sv, xi, rho, theta);
    const double sigV = xi / SQRT3;
    const double A = rho * lodeR(theta, mp.e) / (SQRT6 * mp.fc) + sigV / mp.fc;
    double rad = (eps0 * eps0 * mp.m0 * mp.m0 / 4.0) * A * A
               + 3.0 * eps0 * eps0 * rho * rho / (2.0 * mp.fc * mp.fc);
    if (rad < 0.0) rad = 0.0;
    return (eps0 * mp.m0 / 2.0) * A + std::sqrt(rad);
}

inline double alphaCompression(const double sig_pr[3])
{
    double nrm2 = 0.0, num = 0.0;
    for (int i = 0; i < 3; ++i) {
        nrm2 += sig_pr[i] * sig_pr[i];
        const double spc = sig_pr[i] < 0.0 ? sig_pr[i] : 0.0;   // negative part
        num += spc * sig_pr[i];
    }
    if (nrm2 <= 1.0e-300) return 0.0;
    return num / nrm2;
}

inline void damageDrivers(const double sig_pr[3], const Params& mp, double& et, double& ac, double& xs)
{
    et = equivStrainGeneral(sig_pr, mp);
    ac = alphaCompression(sig_pr);
    const double sv[6] = { sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0 };
    double xi, rho, theta; invariants(sv, xi, rho, theta);
    const double sigV = xi / SQRT3;
    const double Rs = (sigV <= 0.0 && rho > 1.0e-12) ? (-SQRT6 * sigV / rho) : 0.0;   // Eq.57
    xs = 1.0 + (mp.As - 1.0) * Rs;                                                     // Eq.56
}

inline double solveOmegaBracketed(double kd1, double kd2, double sig_eff, double f, double eps_f)
{
    auto Fof = [&](double w) { return (1.0 - w) * sig_eff - f * std::exp(-(kd1 + w * kd2) / eps_f); };
    double lo = 0.0, hi = 1.0;
    if (Fof(lo) <= 0.0) return 0.0;   // not damage-loading
    if (Fof(hi) >= 0.0) return 1.0;   // fully damaged
    double w = 0.5;
    for (int it = 0; it < 100; ++it) {
        const double Fw = Fof(w);
        if (std::fabs(Fw) < 1.0e-13 * (f + 1.0)) break;
        if (Fw > 0.0) lo = w; else hi = w;
        const double dF = -sig_eff + f * std::exp(-(kd1 + w * kd2) / eps_f) * kd2 / eps_f;
        const double wn = (dF != 0.0) ? (w - Fw / dF) : w;       // safeguarded Newton...
        w = (lo < wn && wn < hi) ? wn : 0.5 * (lo + hi);         // ...bisection fallback
    }
    return w < 0.0 ? 0.0 : (w > 1.0 ? 1.0 : w);
}

// ===========================================================================
// Return-map result bundles (principal space). The "internals" (xi,rho,dlam,r,
// rho_tr, s_tr) are exposed so the consistent-tangent assembly can form the
// principal Jacobian d(sig_principal)/d(sig_trial_principal) analytically by the
// implicit-function theorem on the SAME residual the Newton solved.
// ===========================================================================
struct PrincipalResult {
    double sp[3];        // returned principal stresses (paired index-wise with the trial)
    double kp;           // updated kappa_p (== input for the perfect-plastic map)
    bool   plastic;
    bool   converged;
    bool   apex;
    double f_after;      // INDEPENDENT f at the returned stress (honest convergence)
    // internals at the converged invariant state (radial, non-apex):
    double xi, rho, dlam, r, rho_tr;
    double s_tr[3];      // trial deviatoric principal stresses
};

// ---------------------------------------------------------------------------
// Perfect-plastic principal return (Grassl Eq.21 + Lode-independent flow Eq.22-25).
// 3-unknown (xi,rho,dlam) Newton with ANALYTIC 3x3 Jacobian; theta frozen at trial
// (radial deviatoric return); hydrostatic-tension apex fallback. Mirrors the oracle
// return_map_principal.
// ---------------------------------------------------------------------------
inline PrincipalResult returnMapPrincipal(const double sigTr[3], const Params& mp,
                                          double tol = 1.0e-12)
{
    PrincipalResult R;
    const double fc = mp.fc, m0 = mp.m0, Df = mp.Df, K = bulkK(mp), G = shearG(mp);
    double sv[6] = {sigTr[0], sigTr[1], sigTr[2], 0, 0, 0};
    double xi_tr, rho_tr, th_tr;
    invariants(sv, xi_tr, rho_tr, th_tr);
    R.kp = 0.0; R.apex = false;
    const double r = lodeR(th_tr, mp.e);
    R.r = r; R.rho_tr = rho_tr;
    const double p_tr = xi_tr / SQRT3;
    for (int a = 0; a < 3; ++a) R.s_tr[a] = sigTr[a] - p_tr;

    const double f_tr = yieldF(sv, mp);
    if (f_tr <= tol * fc) {
        for (int a = 0; a < 3; ++a) R.sp[a] = sigTr[a];
        R.plastic = false; R.converged = true; R.f_after = f_tr;
        R.xi = xi_tr; R.rho = rho_tr; R.dlam = 0.0;
        return R;
    }

    const double m_v = Df * m0 / (SQRT3 * fc);     // dg/dxi (Lode-independent flow)
    double xi = xi_tr, rho = rho_tr, dlam = 0.0;
    bool apex = false, converged = false;
    for (int it = 0; it < 100; ++it) {
        const double m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc);   // dg/drho (NO Lode r)
        const double R1 = xi - xi_tr + 3.0 * K * dlam * m_v;
        const double R2 = rho - rho_tr + 2.0 * G * dlam * m_s;
        const double R3 = yfInv(xi, rho, r, mp);
        if (std::fabs(R1) + std::fabs(R2) + std::fabs(R3) < tol * (fc + 1.0)) { converged = true; break; }
        // analytic 3x3 Jacobian (matches the oracle)
        const double J00 = 1.0,                       J01 = 0.0,                                   J02 = 3.0 * K * m_v;
        const double J10 = 0.0,                       J11 = 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), J12 = 2.0 * G * m_s;
        const double J20 = m0 / (SQRT3 * fc),         J21 = 3.0 * rho / (fc * fc) + m0 * r / (SQRT6 * fc), J22 = 0.0;
        // solve J dx = -R (3x3 by cofactors)
        const double det = J00 * (J11 * J22 - J12 * J21) - J01 * (J10 * J22 - J12 * J20) + J02 * (J10 * J21 - J11 * J20);
        const double b0 = -R1, b1 = -R2, b2 = -R3;
        const double dxi  = ( b0 * (J11 * J22 - J12 * J21) - J01 * (b1 * J22 - J12 * b2) + J02 * (b1 * J21 - J11 * b2) ) / det;
        const double drho = ( J00 * (b1 * J22 - J12 * b2) - b0 * (J10 * J22 - J12 * J20) + J02 * (J10 * b2 - b1 * J20) ) / det;
        const double dl   = ( J00 * (J11 * b2 - b1 * J21) - J01 * (J10 * b2 - b1 * J20) + b0 * (J10 * J21 - J11 * J20) ) / det;
        xi += dxi; rho += drho; dlam += dl;
        if (rho < 0.0) { apex = true; break; }
    }

    R.xi = xi; R.rho = rho; R.dlam = dlam;
    if (apex) {
        xi = SQRT3 * fc / m0;                 // hydrostatic-tension vertex (f(xi,0)=0)
        const double p_new = xi / SQRT3;
        for (int a = 0; a < 3; ++a) R.sp[a] = p_new;
        R.xi = xi; R.rho = 0.0;
        R.dlam = (m_v > 0.0) ? (xi_tr - xi) / (3.0 * K * m_v) : 0.0;
        R.apex = true; R.plastic = true;
        const double f_apex = yfInv(xi, 0.0, r, mp);
        R.f_after = f_apex;
        R.converged = std::fabs(f_apex) < tol * (fc + 1.0);
        return R;
    }

    const double p_new = xi / SQRT3;
    const double dev_scale = (rho_tr > 0.0) ? rho / rho_tr : 0.0;
    for (int a = 0; a < 3; ++a) R.sp[a] = R.s_tr[a] * dev_scale + p_new;
    double svn[6] = {R.sp[0], R.sp[1], R.sp[2], 0, 0, 0};
    R.f_after = yieldF(svn, mp);
    R.plastic = true; R.converged = converged;
    return R;
}

// ---------------------------------------------------------------------------
// Hardening principal return (Grassl Eq.18 + Eq.30-36). 4-unknown (xi,rho,dlam,kp)
// Newton with ANALYTIC 4x4 Jacobian (the oracle uses a NUMERICAL Jacobian — this is
// the build-PR deliverable); theta frozen; hydrostatic apex fallback; HONEST
// convergence (independent f at the returned stress w/ its own Lode angle).
// ---------------------------------------------------------------------------
inline PrincipalResult returnMapHardening(const double sigTr[3], const Params& mp, double kp_n,
                                          double tol = 1.0e-11)
{
    PrincipalResult R;
    const double fc = mp.fc, m0 = mp.m0, Df = mp.Df, K = bulkK(mp), G = shearG(mp);
    double sv[6] = {sigTr[0], sigTr[1], sigTr[2], 0, 0, 0};
    double xi_tr, rho_tr, th_tr;
    invariants(sv, xi_tr, rho_tr, th_tr);
    const double r = lodeR(th_tr, mp.e);
    R.r = r; R.rho_tr = rho_tr; R.apex = false;
    const double p_tr = xi_tr / SQRT3;
    for (int a = 0; a < 3; ++a) R.s_tr[a] = sigTr[a] - p_tr;

    const double f_tr = yfInvHard(xi_tr, rho_tr, r, kp_n, mp);
    if (f_tr <= tol * fc) {
        for (int a = 0; a < 3; ++a) R.sp[a] = sigTr[a];
        R.kp = kp_n; R.plastic = false; R.converged = true; R.f_after = f_tr;
        R.xi = xi_tr; R.rho = rho_tr; R.dlam = 0.0;
        return R;
    }

    const double cos2 = (2.0 * std::cos(th_tr)) * (2.0 * std::cos(th_tr));
    const double m_v = Df * m0 / (SQRT3 * fc);
    double xi = xi_tr, rho = rho_tr, dlam = 0.0, kp = kp_n;
    bool apex = false, converged = false;
    for (int it = 0; it < 100; ++it) {
        const double m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc);
        const double mnorm = std::sqrt(m_v * m_v + m_s * m_s);
        const double sigV = xi / SQRT3;
        const double xh = ductilityXh(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
        const double R1 = xi - xi_tr + 3.0 * K * dlam * m_v;
        const double R2 = rho - rho_tr + 2.0 * G * dlam * m_s;
        const double R3 = yfInvHard(xi, rho, r, kp, mp);
        const double R4 = kp - kp_n - dlam * mnorm / xh * cos2;
        if (std::fabs(R1) < tol * fc && std::fabs(R2) < tol * fc
            && std::fabs(R3) < tol && std::fabs(R4) < tol) { converged = true; break; }

        // analytic 4x4 Jacobian
        const double q1 = qh1Of(kp, mp.qh0, mp.Hp), q2 = qh2Of(kp, mp.Hp);
        const double dq1 = dqh1OfdKp(kp, mp.qh0, mp.Hp), dq2 = dqh2OfdKp(kp, mp.Hp);
        const double sigV_fc = xi / (SQRT3 * fc);
        const double AV = rho / (SQRT6 * fc) + sigV_fc;
        const double RR = rho * r / (SQRT6 * fc) + sigV_fc;
        const double quad = SQRT1_5 * rho / fc;
        const double cap = (1.0 - q1) * AV * AV + quad;
        // dR3
        const double dAV_dxi = 1.0 / (SQRT3 * fc), dAV_drho = 1.0 / (SQRT6 * fc);
        const double dcap_dxi = (1.0 - q1) * 2.0 * AV * dAV_dxi;
        const double dcap_drho = (1.0 - q1) * 2.0 * AV * dAV_drho + SQRT1_5 / fc;
        const double dcap_dkp = -dq1 * AV * AV;
        const double dRR_dxi = 1.0 / (SQRT3 * fc), dRR_drho = r / (SQRT6 * fc);
        const double dR3_dxi  = 2.0 * cap * dcap_dxi + m0 * q1 * q1 * q2 * dRR_dxi;
        const double dR3_drho = 2.0 * cap * dcap_drho + m0 * q1 * q1 * q2 * dRR_drho;
        const double dq1sq_q2 = 2.0 * q1 * q2 * dq1 + q1 * q1 * dq2;
        const double dq1sq_q2sq = 2.0 * q1 * q2 * q2 * dq1 + 2.0 * q1 * q1 * q2 * dq2;
        const double dR3_dkp = 2.0 * cap * dcap_dkp + m0 * RR * dq1sq_q2 - dq1sq_q2sq;
        // dR4
        const double dms_drho = 3.0 / (fc * fc);
        const double dmnorm_drho = (m_s * dms_drho) / mnorm;
        const double dxh_dsigV = dDuctilityXhdSigV(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
        const double dxh_dxi = dxh_dsigV / SQRT3;
        const double g_val = mnorm / xh;
        const double dg_dxi  = -mnorm / (xh * xh) * dxh_dxi;
        const double dg_drho = dmnorm_drho / xh;
        const double J[4][4] = {
            { 1.0, 0.0, 3.0 * K * m_v, 0.0 },
            { 0.0, 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), 2.0 * G * m_s, 0.0 },
            { dR3_dxi, dR3_drho, 0.0, dR3_dkp },
            { -dlam * cos2 * dg_dxi, -dlam * cos2 * dg_drho, -g_val * cos2, 1.0 }
        };
        double b[4] = { -R1, -R2, -R3, -R4 };
        // 4x4 solve via Gauss elimination w/ partial pivot
        double M[4][5];
        for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = J[i][j]; M[i][4] = b[i]; }
        for (int c = 0; c < 4; ++c) {
            int piv = c; for (int rr = c + 1; rr < 4; ++rr) if (std::fabs(M[rr][c]) > std::fabs(M[piv][c])) piv = rr;
            for (int j = 0; j < 5; ++j) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
            for (int rr = 0; rr < 4; ++rr) if (rr != c) { double f = M[rr][c] / M[c][c]; for (int j = c; j < 5; ++j) M[rr][j] -= f * M[c][j]; }
        }
        xi  += M[0][4] / M[0][0];
        rho += M[1][4] / M[1][1];
        dlam+= M[2][4] / M[2][2];
        kp  += M[3][4] / M[3][3];
        if (rho < 0.0) { apex = true; break; }
    }

    R.xi = xi; R.rho = rho; R.dlam = dlam; R.kp = kp;
    if (apex) {
        // 1-D apex: f(xi,0;kp)=0 by bisection (mirror the oracle)
        double lo = 0.0, hi = SQRT3 * fc / m0 * 4.0;
        for (int it = 0; it < 80; ++it) {
            double mid = 0.5 * (lo + hi);
            if (yfInvHard(lo, 0.0, r, kp, mp) * yfInvHard(mid, 0.0, r, kp, mp) <= 0.0) hi = mid; else lo = mid;
        }
        xi = 0.5 * (lo + hi);
        const double p_new = xi / SQRT3;
        for (int a = 0; a < 3; ++a) R.sp[a] = p_new;
        R.xi = xi; R.rho = 0.0; R.apex = true;
    } else {
        const double p_new = xi / SQRT3;
        const double dev_scale = (rho_tr > 0.0) ? rho / rho_tr : 0.0;
        for (int a = 0; a < 3; ++a) R.sp[a] = R.s_tr[a] * dev_scale + p_new;
    }
    // HONEST convergence: recompute f at the returned stress with ITS OWN Lode angle.
    double svn[6] = {R.sp[0], R.sp[1], R.sp[2], 0, 0, 0};
    R.f_after = yieldF(svn, mp, qh1Of(kp, mp.qh0, mp.Hp), qh2Of(kp, mp.Hp));
    // ADMISSIBILITY + SAFETY (PR #249 adversarial-review fix). A valid plastic return needs
    // dlam>=0 and a NON-DECREASING hardening variable (kp>=kp_n). The semi-implicit hardening
    // Newton can overshoot to rho<0; the apex bisection then returns a SIGN-FLIPPED hydrostatic-
    // tension vertex with kp<kp_n (or kp<0) for a deep-COMPRESSION trial — the documented KNOWN
    // GAP (handoff §6: low-kp / off-meridian first-yield + the apex needs sub-stepping + a
    // dedicated apex sub-algorithm). Be HONEST (never report converged for an inadmissible/off-
    // surface state — the old `(converged||apex)&&|f|<tol` lied here) and SAFE (never hand back
    // the garbage iterate — kp<kp_n would poison the committed history): fall back to the elastic
    // predictor so the caller's status!=0 cuts the step. This deliberately diverges from the
    // numpy oracle's (equally-arbitrary) apex teleport — the kernel is the safe reference here.
    const bool admissible = std::isfinite(R.f_after) && dlam >= -1.0e-12 && kp >= kp_n - 1.0e-12;
    R.converged = (converged || apex) && std::fabs(R.f_after) < 1.0e-7 * (fc + 1.0) && admissible;
    if (!R.converged) {
        for (int a = 0; a < 3; ++a) R.sp[a] = sigTr[a];   // safe fallback = elastic predictor
        R.kp = kp_n; R.xi = xi_tr; R.rho = rho_tr; R.dlam = 0.0; R.apex = false;
    }
    R.plastic = true;
    return R;
}

// §4.6 — lateral mixed-control residual for the confined-fiber / triaxial driver.
//   mode 0 free (sigma_lat=0), 1 active (sigma_lat=-p), 2 passive (sigma_lat=-sigma_hoop(eps_lat)).
// Returns the residual whose root (over the lateral strains) the condensation Newton drives to 0.
inline double lateralResidual(double sigmaLat, double epsLat, int mode, double p,
                              double (*hoopLaw)(double))
{
    if (mode == 0) return sigmaLat;
    if (mode == 1) return sigmaLat + p;
    if (mode == 2) return sigmaLat + (hoopLaw ? hoopLaw(epsLat) : 0.0);
    return sigmaLat;
}

// ===========================================================================
// 3x3 symmetric eigensolver (cyclic Jacobi). w[a] paired with eigenvector V[*][a]
// (columns). Robust + orthonormal V — the pairing (sp_a <-> w_a <-> V[:,a]) is what
// the radial deviatoric return preserves, so NO sorting is needed (the MW invariants
// are symmetric in the three principal values).
// ===========================================================================
inline void eig3sym(const double A[3][3], double w[3], double V[3][3])
{
    double a[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) { a[i][j] = A[i][j]; V[i][j] = (i == j) ? 1.0 : 0.0; }
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        if (off < 1.0e-300) break;
        for (int p = 0; p < 2; ++p) for (int q = p + 1; q < 3; ++q) {
            if (std::fabs(a[p][q]) < 1.0e-300) continue;
            const double app = a[p][p], aqq = a[q][q], apq = a[p][q];
            const double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
            const double c = std::cos(phi), s = std::sin(phi);
            for (int k = 0; k < 3; ++k) {
                const double akp = a[k][p], akq = a[k][q];
                a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
            }
            for (int k = 0; k < 3; ++k) {
                const double apk = a[p][k], aqk = a[q][k];
                a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
            }
            for (int k = 0; k < 3; ++k) {
                const double vkp = V[k][p], vkq = V[k][q];
                V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
            }
        }
    }
    for (int i = 0; i < 3; ++i) w[i] = a[i][i];
}

inline void voigtToMat(const double v[6], double M[3][3])
{
    M[0][0] = v[0]; M[1][1] = v[1]; M[2][2] = v[2];
    M[0][1] = M[1][0] = v[3]; M[1][2] = M[2][1] = v[4]; M[0][2] = M[2][0] = v[5];
}
inline void matToVoigt(const double M[3][3], double v[6])
{
    v[0] = M[0][0]; v[1] = M[1][1]; v[2] = M[2][2];
    v[3] = M[0][1]; v[4] = M[1][2]; v[5] = M[0][2];
}

// Elastic tensor in the ORACLE convention (off-diagonal slots = TRUE tensor components;
// C[i][i]=2G for the shear rows, i.e. dsig_ij = 2G deps_ij). This matches
// tests/_testbed/concrete3d_ref.py::elastic_C so the C++ tangent numbers diff 1:1.
inline void elasticC(const Params& mp, double C[6][6])
{
    const double K = bulkK(mp), G = shearG(mp); const double lam = K - 2.0 * G / 3.0;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) C[i][j] = 0.0;
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) C[i][j] = lam; C[i][i] += 2.0 * G; }
    for (int i = 3; i < 6; ++i) C[i][i] = 2.0 * G;
}

// Trial stress = sig_n + C:deps (oracle elastic_pred_tensor; tensor-component deps).
inline void elasticPredTensor(const double sig_n[6], const double deps[6], const Params& mp, double sig_tr[6])
{
    const double K = bulkK(mp), G = shearG(mp); const double lam = K - 2.0 * G / 3.0;
    const double tr = deps[0] + deps[1] + deps[2];
    for (int i = 0; i < 6; ++i) sig_tr[i] = sig_n[i] + 2.0 * G * deps[i];
    for (int i = 0; i < 3; ++i) sig_tr[i] += lam * tr;
}

// ---------------------------------------------------------------------------
// Principal Jacobian D[a][b] = d(sp_a)/d(w_b) by the implicit-function theorem on
// the SAME residual the inner Newton solved (perfect-plastic OR hardening). This is
// the analytic backbone of the consistent tangent. The single Lode directional
// gradient dr/dw (corner-singular in closed form: dr/dtheta->0 cancels 1/sin3theta
// to a finite limit) is taken by a robust scalar central difference — an isolated,
// well-conditioned micro-derivative; the rest is closed-form. The whole 6x6 is then
// FD-verified against the numpy oracle (gate floors 1e-6).
// ---------------------------------------------------------------------------
inline double lodeOfPrincipal(const double w[3], double e)
{
    double sv[6] = {w[0], w[1], w[2], 0, 0, 0}, xi, rho, th;
    invariants(sv, xi, rho, th);
    return lodeR(th, e);
}

inline void principalJacobian(const double w[3], const PrincipalResult& pr, const Params& mp,
                              bool hardening, double D[3][3])
{
    const double fc = mp.fc, m0 = mp.m0, K = bulkK(mp), G = shearG(mp);
    const double rho_tr = pr.rho_tr, rho = pr.rho, dlam = pr.dlam, r = pr.r, xi = pr.xi, kp = pr.kp;
    const double m_v = mp.Df * m0 / (SQRT3 * fc);
    const double m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc);

    // dr/dw_b : scalar central difference of the (smooth-in-w) Lode r.
    double dr_dw[3];
    {
        const double h = 1.0e-6 * (fc + 1.0);
        for (int b = 0; b < 3; ++b) {
            double wp[3] = {w[0], w[1], w[2]}, wm[3] = {w[0], w[1], w[2]};
            wp[b] += h; wm[b] -= h;
            dr_dw[b] = (lodeOfPrincipal(wp, mp.e) - lodeOfPrincipal(wm, mp.e)) / (2.0 * h);
        }
    }

    // Inner Newton Jacobian J_u (rows R1,R2,R3[,R4]) at the converged state, and the
    // trial-derivative columns dR/dw_b. Then d(xi,rho,...)/dw_b = -J_u^{-1} dR/dw_b.
    // We only need dxi/dw_b and drho/dw_b for D[a][b].
    if (!hardening) {
        // 3x3 J_u
        const double J[3][3] = {
            { 1.0, 0.0, 3.0 * K * m_v },
            { 0.0, 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), 2.0 * G * m_s },
            { m0 / (SQRT3 * fc), 3.0 * rho / (fc * fc) + m0 * r / (SQRT6 * fc), 0.0 }
        };
        const double det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0]) + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        for (int b = 0; b < 3; ++b) {
            const double s_tr_b = pr.s_tr[b];
            // dR/dw_b = [ -dxi_tr/dw_b ; -drho_tr/dw_b ; dR3/dr * dr/dw_b ]
            const double rhs0 = -(1.0 / SQRT3);
            const double rhs1 = -((rho_tr > 0.0) ? s_tr_b / rho_tr : 0.0);
            const double rhs2 = (m0 * rho / (SQRT6 * fc)) * dr_dw[b];
            // solve J u = -rhs  (so u = d(xi,rho,dlam)/dw_b)
            const double c0 = -rhs0, c1 = -rhs1, c2 = -rhs2;
            const double dxi  = ( c0*(J[1][1]*J[2][2]-J[1][2]*J[2][1]) - J[0][1]*(c1*J[2][2]-J[1][2]*c2) + J[0][2]*(c1*J[2][1]-J[1][1]*c2) ) / det;
            const double drho = ( J[0][0]*(c1*J[2][2]-J[1][2]*c2) - c0*(J[1][0]*J[2][2]-J[1][2]*J[2][0]) + J[0][2]*(J[1][0]*c2-c1*J[2][0]) ) / det;
            const double dscale = (rho_tr > 0.0) ? (drho * rho_tr - rho * (s_tr_b / rho_tr)) / (rho_tr * rho_tr) : 0.0;
            const double scale = (rho_tr > 0.0) ? rho / rho_tr : 0.0;
            for (int a = 0; a < 3; ++a) {
                const double ds_tr_a = (a == b ? 1.0 : 0.0) - 1.0 / 3.0;
                D[a][b] = ds_tr_a * scale + pr.s_tr[a] * dscale + (1.0 / SQRT3) * dxi;
            }
        }
    } else {
        // 4x4 J_u (xi,rho,dlam,kp); reuse the analytic entries from returnMapHardening.
        double sv[6] = {w[0], w[1], w[2], 0, 0, 0}, xitr, rhotr, thtr;
        invariants(sv, xitr, rhotr, thtr);
        const double cos2 = (2.0 * std::cos(thtr)) * (2.0 * std::cos(thtr));
        const double mnorm = std::sqrt(m_v * m_v + m_s * m_s);
        const double sigV = xi / SQRT3;
        const double xh = ductilityXh(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
        const double q1 = qh1Of(kp, mp.qh0, mp.Hp), q2 = qh2Of(kp, mp.Hp);
        const double dq1 = dqh1OfdKp(kp, mp.qh0, mp.Hp), dq2 = dqh2OfdKp(kp, mp.Hp);
        const double sigV_fc = xi / (SQRT3 * fc);
        const double AV = rho / (SQRT6 * fc) + sigV_fc;
        const double RR = rho * r / (SQRT6 * fc) + sigV_fc;
        const double quad = SQRT1_5 * rho / fc;
        const double cap = (1.0 - q1) * AV * AV + quad;
        const double dcap_dxi = (1.0 - q1) * 2.0 * AV * (1.0 / (SQRT3 * fc));
        const double dcap_drho = (1.0 - q1) * 2.0 * AV * (1.0 / (SQRT6 * fc)) + SQRT1_5 / fc;
        const double dcap_dkp = -dq1 * AV * AV;
        const double dR3_dxi  = 2.0 * cap * dcap_dxi + m0 * q1 * q1 * q2 * (1.0 / (SQRT3 * fc));
        const double dR3_drho = 2.0 * cap * dcap_drho + m0 * q1 * q1 * q2 * (r / (SQRT6 * fc));
        const double dq1sq_q2 = 2.0 * q1 * q2 * dq1 + q1 * q1 * dq2;
        const double dq1sq_q2sq = 2.0 * q1 * q2 * q2 * dq1 + 2.0 * q1 * q1 * q2 * dq2;
        const double dR3_dkp = 2.0 * cap * dcap_dkp + m0 * RR * dq1sq_q2 - dq1sq_q2sq;
        const double dmnorm_drho = (m_s * (3.0 / (fc * fc))) / mnorm;
        const double dxh_dxi = dDuctilityXhdSigV(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh) / SQRT3;
        const double g_val = mnorm / xh;
        const double dg_dxi  = -mnorm / (xh * xh) * dxh_dxi;
        const double dg_drho = dmnorm_drho / xh;
        const double Ju[4][4] = {
            { 1.0, 0.0, 3.0 * K * m_v, 0.0 },
            { 0.0, 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), 2.0 * G * m_s, 0.0 },
            { dR3_dxi, dR3_drho, 0.0, dR3_dkp },
            { -dlam * cos2 * dg_dxi, -dlam * cos2 * dg_drho, -g_val * cos2, 1.0 }
        };
        // dr depends on theta (frozen) which depends on w; the R4 cos2 term ALSO depends on
        // theta(w). For the principal-block Jacobian we include the dominant trial couplings
        // (xi_tr,rho_tr) analytically and the r/theta coupling via dr/dw (R3) + dcos2/dw (R4).
        // dcos2/dw_b via the same scalar-FD route as dr/dw (cheap, robust).
        double dcos2_dw[3];
        {
            const double h = 1.0e-6 * (fc + 1.0);
            for (int b = 0; b < 3; ++b) {
                double wp[3] = {w[0], w[1], w[2]}, wm[3] = {w[0], w[1], w[2]};
                wp[b] += h; wm[b] -= h;
                double svp[6] = {wp[0],wp[1],wp[2],0,0,0}, svm[6] = {wm[0],wm[1],wm[2],0,0,0};
                double x1,r1,t1,x2,r2,t2; invariants(svp,x1,r1,t1); invariants(svm,x2,r2,t2);
                const double cp = (2.0*std::cos(t1))*(2.0*std::cos(t1));
                const double cm = (2.0*std::cos(t2))*(2.0*std::cos(t2));
                dcos2_dw[b] = (cp - cm) / (2.0 * h);
            }
        }
        for (int b = 0; b < 3; ++b) {
            const double s_tr_b = pr.s_tr[b];
            const double rhs0 = -(1.0 / SQRT3);
            const double rhs1 = -((rhotr > 0.0) ? s_tr_b / rhotr : 0.0);
            // dR3/dr = d(yfInvHard)/dr = m0 q1^2 q2 * rho/(sqrt6 fc)  (q1^2 q2 == 1 only when
            // perfect-plastic; omitting it was a hardening-only tangent bug).
            const double rhs2 = (m0 * q1 * q1 * q2 * rho / (SQRT6 * fc)) * dr_dw[b];
            const double rhs3 = -dlam * g_val * dcos2_dw[b];            // dR4/dcos2 * dcos2/dw_b
            double rhs[4] = { rhs0, rhs1, rhs2, rhs3 };
            // solve Ju u = -rhs  (Gauss w/ pivot)
            double M[4][5];
            for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = Ju[i][j]; M[i][4] = -rhs[i]; }
            for (int c = 0; c < 4; ++c) {
                int piv = c; for (int rr = c + 1; rr < 4; ++rr) if (std::fabs(M[rr][c]) > std::fabs(M[piv][c])) piv = rr;
                for (int j = 0; j < 5; ++j) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
                for (int rr = 0; rr < 4; ++rr) if (rr != c) { double f = M[rr][c] / M[c][c]; for (int j = c; j < 5; ++j) M[rr][j] -= f * M[c][j]; }
            }
            const double dxi  = M[0][4] / M[0][0];
            const double drho = M[1][4] / M[1][1];
            const double dscale = (rhotr > 0.0) ? (drho * rhotr - rho * (s_tr_b / rhotr)) / (rhotr * rhotr) : 0.0;
            const double scale = (rhotr > 0.0) ? rho / rhotr : 0.0;
            for (int a = 0; a < 3; ++a) {
                const double ds_tr_a = (a == b ? 1.0 : 0.0) - 1.0 / 3.0;
                D[a][b] = ds_tr_a * scale + pr.s_tr[a] * dscale + (1.0 / SQRT3) * dxi;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Consistent (algorithmic) tangent dsigma/depsilon (6x6, oracle tensor convention).
// Spectral lift of the principal Jacobian D[a][b]: dsigma/dsig_tr (isotropic
// tensor-function derivative, de Souza Neto-Peric-Owen 2008 Box A.6) then : C.
// NON-SYMMETRIC for non-associated flow (and ~2% even associated, from the spectral
// recompose) => Tier-1 needs an unsymmetric solver UNCONDITIONALLY.
// ---------------------------------------------------------------------------
inline void consistentTangent(const double sig_tr[6], const double w[3], const double V[3][3],
                              const PrincipalResult& pr, const Params& mp, bool hardening,
                              double Dtan6[6][6])
{
    // Elastic step, the non-converged safe fallback (returnMapHardening reset to the elastic
    // predictor), or a converged apex => the elastic operator. NOTE (PR #249 review): at a true
    // apex the physical tangent collapses toward zero (the stress is pinned at the vertex — an
    // ATTRACTING PLATEAU, NOT measure-zero), so this elastic surrogate is ~K too STIFF there and
    // can slow the global Newton. That is acceptable for P1 because the apex is the deferred
    // KNOWN-GAP regime (handoff §6) and the return map flags it (converged=false ⇒ step-cut);
    // the dedicated rank-deficient apex tangent lands with the apex sub-algorithm (P2+).
    if (!pr.plastic || !pr.converged || pr.apex) { elasticC(mp, Dtan6); return; }
    // eigenprojections E_a = e_a (x) e_a (rank-1)
    double E[3][3][3];
    for (int a = 0; a < 3; ++a) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) E[a][i][j] = V[i][a] * V[j][a];
    double D[3][3];
    principalJacobian(w, pr, mp, hardening, D);

    // dsigma/dsig_tr as a 4th-order tensor 𝔻_ijkl
    double Dt[3][3][3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) {
        double val = 0.0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) val += D[a][b] * E[a][i][j] * E[b][k][l];
        Dt[i][j][k][l] = val;
    }
    // spin part: sum_{a!=b} gamma_ab * 0.5*(E_a_ik E_b_jl + E_a_il E_b_jk)
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
        if (a == b) continue;
        double dw = w[a] - w[b];
        double gamma;
        if (std::fabs(dw) > 1.0e-9 * (mp.fc + 1.0)) gamma = (pr.sp[a] - pr.sp[b]) / dw;
        else gamma = D[a][a] - D[a][b];          // l'Hopital limit (repeated eigenvalues)
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l)
            Dt[i][j][k][l] += gamma * 0.5 * (E[a][i][k] * E[b][j][l] + E[a][i][l] * E[b][j][k]);
    }

    // contract with the elastic operator: dsigma/deps_ij = 𝔻_ijkl C_klmn ... but C is
    // isotropic so 𝔻:C in tensor form. Build C as a 4th-order tensor (lam dij dkl + 2G I_sym).
    const double K = bulkK(mp), G = shearG(mp); const double lam = K - 2.0 * G / 3.0;
    double Ct[3][3][3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) {
        const double dij = (i == j) ? 1.0 : 0.0, dkl = (k == l) ? 1.0 : 0.0;
        const double dik = (i == k) ? 1.0 : 0.0, djl = (j == l) ? 1.0 : 0.0, dil = (i == l) ? 1.0 : 0.0, djk = (j == k) ? 1.0 : 0.0;
        Ct[i][j][k][l] = lam * dij * dkl + G * (dik * djl + dil * djk);
    }
    double T[3][3][3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int m = 0; m < 3; ++m) for (int n = 0; n < 3; ++n) {
        double v = 0.0;
        for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) v += Dt[i][j][k][l] * Ct[k][l][m][n];
        T[i][j][m][n] = v;
    }
    // pack to the 6x6 ORACLE Voigt convention: rows/cols {00,11,22,01,12,02}; the column
    // (mn) sum over the full tensor double-counts off-diagonal pairs, so the shear COLUMNS
    // carry a factor 2 (deps_mn=deps_nm). Rows take the symmetric component.
    const int I[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};
    for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) {
        const int i = I[A][0], j = I[A][1], m = I[B][0], n = I[B][1];
        double v = 0.5 * (T[i][j][m][n] + T[j][i][m][n]);
        if (B >= 3) v += 0.5 * (T[i][j][n][m] + T[j][i][n][m]);   // off-diag column => +(nm)
        Dtan6[A][B] = v;
    }
}

// ===========================================================================
// Committed history. Small-strain plastic-damage state. (P2 will add the damage
// kappas; the projector frame is recomputed from sig_tr each step, not stored.)
// ===========================================================================
struct State {
    double eps[6] = {0, 0, 0, 0, 0, 0};   // committed total strain (tensor comps)
    double sig[6] = {0, 0, 0, 0, 0, 0};   // committed NOMINAL stress (tensor comps)
    double sigEff[6] = {0, 0, 0, 0, 0, 0}; // committed EFFECTIVE stress (drives the next return + the damage plastic strain)
    double kp = 0.0;                       // committed kappa_p
    // P2 damage history (CDPM2 §2.3): et_max = max equiv-strain (Eq.43); kd*1/kd*2 the plastic /
    // damage-scaled inelastic-strain parts (Eq.44/45/48/49) feeding eps_i = kd1 + omega*kd2 (Eq.52).
    double et_max = 0.0;
    double kdt1 = 0.0, kdt2 = 0.0;         // tensile
    double kdc = 0.0, kdc1 = 0.0, kdc2 = 0.0;  // compressive (kdc = the alpha_c-weighted history)
    // P3 Tier-2 IMPL-EX bookkeeping (committed IMPLICIT damage + the per-variable increments and the
    // committed dt, for the next step's extrapolation x~ = x_n + (dt/dt_n)*dx_n). Unused when !implex.
    double wt = 0.0, wc = 0.0;             // committed IMPLICIT dual damage
    double dwt = 0.0, dwc = 0.0;           // committed implicit damage increments
    double depl[6] = {0, 0, 0, 0, 0, 0};   // committed implicit plastic-strain increment (tensor)
    double dt_n = 0.0;                      // committed time step
};

// ---------------------------------------------------------------------------
// Elastic (compliance) plastic strain eps_p = eps - C^-1 : sig_eff (tensor Voigt), isotropic
// closed form (no linear solve) — mirrors the oracle _plastic_strain6.
// ---------------------------------------------------------------------------
inline void plasticStrain6(const double sig_eff[6], const double eps[6], const Params& mp, double epl[6])
{
    const double E = mp.E, nu = mp.nu;
    epl[0] = eps[0] - (sig_eff[0] - nu * (sig_eff[1] + sig_eff[2])) / E;
    epl[1] = eps[1] - (sig_eff[1] - nu * (sig_eff[0] + sig_eff[2])) / E;
    epl[2] = eps[2] - (sig_eff[2] - nu * (sig_eff[0] + sig_eff[1])) / E;
    for (int k = 3; k < 6; ++k) epl[k] = eps[k] - sig_eff[k] * (1.0 + nu) / E;  // eps_ij = sig_ij/(2G)
}

// ---------------------------------------------------------------------------
// P2 dual-damage NOMINAL stress update (mirrors the oracle damaged_step_tensor EXACTLY): from the
// committed history `in` + the new EFFECTIVE stress/strain, accumulate the CDPM2 damage drivers,
// solve omega_t/omega_c (bracketed, physical-floored), and recompose the nominal stress via the
// spectral split  sigma = (1-omega_t)<sig_bar>+ + (1-omega_c)<sig_bar>-  (Eq.1). Writes the new
// damage history into `out`. The split is recomputed from the CONVERGED effective stress every step
// => automatic, tier-independent unilateral crack-closure (ADR §4.3 BLOCKING).
// ---------------------------------------------------------------------------
inline void damagedUpdate(const Params& mp, const State& in, const double sig_eff[6],
                          const double eps_new[6], State& out,
                          double* wtOut = nullptr, double* wcOut = nullptr)
{
    const double eps0 = mp.ft / mp.E;
    const double eps_f  = mp.Gf / (mp.ft * mp.lch);
    const double eps_fc = mp.Gc / (mp.fc * mp.lch);

    double A[3][3], w[3], V[3][3];
    voigtToMat(sig_eff, A);
    eig3sym(A, w, V);

    double et, ac, xs; damageDrivers(w, mp, et, ac, xs);

    // plastic-strain increment (tensor Frobenius norm ||M||_F = sqrt(sum M_ij^2))
    double epl[6], epl_n[6];
    plasticStrain6(sig_eff,   eps_new,   mp, epl);
    plasticStrain6(in.sigEff, in.eps,    mp, epl_n);
    double Md[3][3], Mn[3][3];
    voigtToMat(epl, Md); voigtToMat(epl_n, Mn);
    double dnorm2 = 0.0;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        const double d = Md[i][j] - Mn[i][j]; dnorm2 += d * d;
    }
    const double dnorm = std::sqrt(dnorm2);

    const double et_max_n = in.et_max;
    const double det_raw = (et - et_max_n > 0.0) ? (et - et_max_n) : 0.0;
    const double lo = et_max_n > eps0 ? et_max_n : eps0;
    const double above = (et - lo > 0.0) ? (et - lo) : 0.0;       // increment ABOVE the onset eps0
    const bool loading = det_raw > 0.0 && et > eps0;

    double kdt1 = in.kdt1, kdt2 = in.kdt2, kdc = in.kdc, kdc1 = in.kdc1, kdc2 = in.kdc2;
    if (loading) {
        kdt2 += above / xs;          kdt1 += dnorm / xs;                 // Eq.45 / Eq.44 (no alpha_c)
        kdc  += ac * above;          kdc2 += ac * above / xs;            // Eq.47 / Eq.49
        kdc1 += ac * dnorm / xs;                                         // Eq.48 (beta_c = 1, monotonic)
    }
    const double et_max = et_max_n > et ? et_max_n : et;

    // softening drives = the extreme effective principal of each sign; physical FLOOR (review-fix):
    // never solve omega on a numerical-residual stress (~1e-10 MPa) -> flips wt 0<->1 in compression.
    double Dt = w[0]; for (int i = 1; i < 3; ++i) if (w[i] > Dt) Dt = w[i]; if (Dt < 0.0) Dt = 0.0;
    double mn = w[0]; for (int i = 1; i < 3; ++i) if (w[i] < mn) mn = w[i];
    double Dc = -mn; if (Dc < 0.0) Dc = 0.0;
    const double wt = (et_max > eps0 && Dt > 1.0e-6 * mp.ft)
                    ? solveOmegaBracketed(kdt1, kdt2, Dt, mp.ft, eps_f) : 0.0;
    const double wc = (kdc > 0.0 && Dc > 1.0e-6 * mp.fc)
                    ? solveOmegaBracketed(kdc1, kdc2, Dc, mp.fc, eps_fc) : 0.0;

    // nominal principal stresses (Eq.1) then recompose on the SAME eigenvectors
    double sp[3];
    for (int i = 0; i < 3; ++i) {
        const double st = w[i] > 0.0 ? w[i] : 0.0;
        const double sc = w[i] < 0.0 ? w[i] : 0.0;
        sp[i] = (1.0 - wt) * st + (1.0 - wc) * sc;
    }
    double S[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = 0.0; for (int a = 0; a < 3; ++a) v += V[i][a] * sp[a] * V[j][a];
        S[i][j] = v;
    }
    matToVoigt(S, out.sig);
    out.et_max = et_max; out.kdt1 = kdt1; out.kdt2 = kdt2;
    out.kdc = kdc; out.kdc1 = kdc1; out.kdc2 = kdc2;
    if (wtOut) *wtOut = wt;   // expose the damage variables for the wrapper's recorders (read-only)
    if (wcOut) *wcOut = wc;
}

// ===========================================================================
// P3b — ANALYTIC dual-projector DAMAGED consistent tangent. Supporting machinery
// (ports the oracle's isotropic_tangent / _dscalar_dsig / damaged_tangent_analytic
// VERBATIM), then the tangent itself. The C++ public entry point returnMap upgrades
// its P1 EFFECTIVE tangent to this damaged tangent when doTangent is requested.
// ===========================================================================

// de Souza Neto-Peric-Owen 2008 Box A.6 derivative dY/dX of an isotropic symmetric-tensor
// function Y = sum_a y(lam_a) E_a, given eigenvalues `lam`, eigenvectors `V` (columns), and the
// per-eigenvalue values yv=y(lam_a), ypv=y'(lam_a):
//     (dY/dX : S) = sum_a ypv_a (E_a:S) E_a + sum_{a!=b} G_ab E_a S E_b,
//     G_ab = (yv_a - yv_b)/(lam_a - lam_b)  (-> ypv_a as lam_b -> lam_a, l'Hopital).
// Mirrors the oracle isotropic_tangent byte-for-byte: operate on real 3x3 matrices, build each
// 6x6 column by applying to the basis Voigt tensor e_j, pack with matToVoigt {00,11,22,01,12,02}.
inline void isotropicTangent(const double lam[3], const double V[3][3],
                             const double yv[3], const double ypv[3], double D6[6][6],
                             double tol = 1.0e-8)
{
    double E[3][3][3];
    for (int a = 0; a < 3; ++a) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
        E[a][i][j] = V[i][a] * V[j][a];
    double G[3][3];
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
        if (a == b) { G[a][b] = 0.0; continue; }
        const double dl = lam[a] - lam[b];
        G[a][b] = (std::fabs(dl) < tol) ? ypv[a] : (yv[a] - yv[b]) / dl;
    }
    for (int j = 0; j < 6; ++j) {
        double Sj[6] = {0, 0, 0, 0, 0, 0}; Sj[j] = 1.0;
        double S[3][3]; voigtToMat(Sj, S);
        double out[3][3];
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) out[i][k] = 0.0;
        for (int a = 0; a < 3; ++a) {                                   // sum_a ypv_a (E_a:S) E_a
            double EaS = 0.0;
            for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) EaS += E[a][i][k] * S[i][k];
            const double c = ypv[a] * EaS;
            for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) out[i][k] += c * E[a][i][k];
        }
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {       // sum_{a!=b} G_ab E_a S E_b
            if (a == b) continue;
            double ES[3][3];
            for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) {
                double v = 0.0; for (int m = 0; m < 3; ++m) v += E[a][i][m] * S[m][k]; ES[i][k] = v;
            }
            double M[3][3];
            for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) {
                double v = 0.0; for (int m = 0; m < 3; ++m) v += ES[i][m] * E[b][m][k]; M[i][k] = v;
            }
            for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) out[i][k] += G[a][b] * M[i][k];
        }
        double out6[6]; matToVoigt(out, out6);
        for (int i = 0; i < 6; ++i) D6[i][j] = out6[i];
    }
}

// 6x6 inverse by Gauss-Jordan w/ partial pivot — the elastic compliance C^-1 for the
// plastic-strain-rate chain term d(eps_p)/d(eps) = I - C^-1 C_eff. Returns false if singular.
inline bool invert6(const double A[6][6], double Inv[6][6])
{
    double M[6][12];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) { M[i][j] = A[i][j]; M[i][j + 6] = (i == j) ? 1.0 : 0.0; }
    for (int c = 0; c < 6; ++c) {
        int piv = c; for (int r = c + 1; r < 6; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
        if (std::fabs(M[piv][c]) < 1.0e-300) return false;
        for (int j = 0; j < 12; ++j) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
        const double d = M[c][c];
        for (int j = 0; j < 12; ++j) M[c][j] /= d;
        for (int r = 0; r < 6; ++r) if (r != c) {
            const double f = M[r][c];
            for (int j = 0; j < 12; ++j) M[r][j] -= f * M[c][j];
        }
    }
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Inv[i][j] = M[i][j + 6];
    return true;
}

// Scalar damage-driver as a function of the full effective stress (eigendecompose -> principals).
//   which 0 = equiv strain (Eq.37), 1 = softening ductility xs (Eq.56-57), 2 = alpha_c (Eq.46).
inline double scalarDriver(int which, const double sig6[6], const Params& mp)
{
    double A[3][3], w[3], V[3][3]; voigtToMat(sig6, A); eig3sym(A, w, V);
    if (which == 0) return equivStrainGeneral(w, mp);
    if (which == 1) { double et, ac, xs; damageDrivers(w, mp, et, ac, xs); return xs; }
    return alphaCompression(w);
}

// Isolated micro-FD gradient d(scalarDriver)/d(sig) per Voigt component (oracle _dscalar_dsig,
// h=1e-6). These are PER-COMPONENT grads of a scalar-of-stress => NO [1,1,1,2,2,2] tensor weight
// (unlike the eigenprojection / ||Deps_p|| gradients below, which DO carry it — LEDGER quirk).
inline void dscalarDsig(int which, const double sig6[6], const Params& mp, double g[6], double h = 1.0e-6)
{
    for (int k = 0; k < 6; ++k) {
        double dp[6], dm[6];
        for (int i = 0; i < 6; ++i) { dp[i] = sig6[i]; dm[i] = sig6[i]; }
        dp[k] += h; dm[k] -= h;
        g[k] = (scalarDriver(which, dp, mp) - scalarDriver(which, dm, mp)) / (2.0 * h);
    }
}

// ---------------------------------------------------------------------------
// P3b ANALYTIC dual-projector DAMAGED consistent tangent (mirrors the oracle
// damaged_tangent_analytic VERBATIM):
//     C = D_dam : C_eff  -  sig_t (x) dω_t/dε  -  sig_c (x) dω_c/dε        (ADR §4.3 MAJOR)
// D_dam = isotropicTangent spectral derivative of the per-principal damaged stress with ω FROZEN
// (Box A.6); the two rank-1 updates carry the ω-sensitivity via the IFT on the bracketed ω-solve
// F(ω)=(1-ω)D - f exp(-(kd1+ω kd2)/eps_f)=0 (H=dF/dω=D[(1-ω)kd2/eps_f - 1]), chained through the
// CDPM2 damage histories. `Ceff` = the P1 EFFECTIVE consistent tangent (from returnMapTensor);
// `sig_eff` = the converged effective stress; `in` = the committed damage history. RECOMPUTES the
// same drivers as damagedUpdate (self-contained, exactly as the oracle does), so the ω used here is
// byte-identical to the update's. FD-verified == the P2d numerical reference ~1e-10.
// Voigt-WEIGHT quirk (LEDGER): the TENSOR gradients (eigenprojection Emax/Emin, ||Deps_p||) carry
// the [1,1,1,2,2,2] double-contraction weight W6 BEFORE Ceff^T; the per-component micro-FD scalar
// grads (det/dxs/dac) do NOT (already per-component).
// ---------------------------------------------------------------------------
inline void damagedTangent(const Params& mp, const State& in, const double sig_eff[6],
                           const double eps_new[6], const double Ceff[6][6], double D6[6][6])
{
    static const double W6[6] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    const double eps0   = mp.ft / mp.E;
    const double eps_f  = mp.Gf / (mp.ft * mp.lch);
    const double eps_fc = mp.Gc / (mp.fc * mp.lch);

    double A[3][3], w[3], V[3][3];
    voigtToMat(sig_eff, A);
    eig3sym(A, w, V);
    double et, ac, xs; damageDrivers(w, mp, et, ac, xs);

    // plastic-strain INCREMENT (tensor Voigt) + its Frobenius norm (== ddot6 with W6)
    double epl[6], epl_n[6], depl[6];
    plasticStrain6(sig_eff,   eps_new, mp, epl);
    plasticStrain6(in.sigEff, in.eps,  mp, epl_n);
    for (int i = 0; i < 6; ++i) depl[i] = epl[i] - epl_n[i];
    double dnorm2 = 0.0;
    for (int i = 0; i < 6; ++i) dnorm2 += W6[i] * depl[i] * depl[i];
    const double dnorm = std::sqrt(dnorm2);

    const double et_max_n = in.et_max;
    const double det_raw = (et - et_max_n > 0.0) ? (et - et_max_n) : 0.0;
    const double loref = et_max_n > eps0 ? et_max_n : eps0;
    const double above = (et - loref > 0.0) ? (et - loref) : 0.0;
    const bool loading = det_raw > 0.0 && et > eps0;

    double kdt1 = in.kdt1, kdt2 = in.kdt2, kdc = in.kdc, kdc1 = in.kdc1, kdc2 = in.kdc2;
    if (loading) {
        kdt2 += above / xs;          kdt1 += dnorm / xs;
        kdc  += ac * above;          kdc2 += ac * above / xs;
        kdc1 += ac * dnorm / xs;
    }
    const double et_max2 = et_max_n > et ? et_max_n : et;

    // extreme effective principals + their eigenprojections (argmax/argmin; eig3sym is unsorted)
    int imax = 0, imin = 0;
    for (int i = 1; i < 3; ++i) { if (w[i] > w[imax]) imax = i; if (w[i] < w[imin]) imin = i; }
    const double Dt = w[imax] > 0.0 ? w[imax] : 0.0;
    const double Dc = (-w[imin]) > 0.0 ? -w[imin] : 0.0;
    const double wt = (et_max2 > eps0 && Dt > 1.0e-6 * mp.ft)
                    ? solveOmegaBracketed(kdt1, kdt2, Dt, mp.ft, eps_f) : 0.0;
    const double wc = (kdc > 0.0 && Dc > 1.0e-6 * mp.fc)
                    ? solveOmegaBracketed(kdc1, kdc2, Dc, mp.fc, eps_fc) : 0.0;

    // D_dam = spectral derivative of the per-principal damaged stress with ω FROZEN
    double yv[3], ypv[3];
    for (int a = 0; a < 3; ++a) {
        const double st = w[a] > 0.0 ? w[a] : 0.0;
        const double sc = w[a] < 0.0 ? w[a] : 0.0;
        yv[a]  = (1.0 - wt) * st + (1.0 - wc) * sc;
        ypv[a] = (w[a] > 0.0) ? (1.0 - wt) : (1.0 - wc);
    }
    double Ddam[6][6];
    isotropicTangent(w, V, yv, ypv, Ddam);

    // sig_t = <sig_bar>+ , sig_c = <sig_bar>- recomposed on the SAME eigenvectors
    double sig_t[6], sig_c[6];
    {
        double St[3][3], Sc[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
            double vt = 0.0, vc = 0.0;
            for (int a = 0; a < 3; ++a) {
                const double sp = w[a];
                vt += V[i][a] * (sp > 0.0 ? sp : 0.0) * V[j][a];
                vc += V[i][a] * (sp < 0.0 ? sp : 0.0) * V[j][a];
            }
            St[i][j] = vt; Sc[i][j] = vc;
        }
        matToVoigt(St, sig_t); matToVoigt(Sc, sig_c);
    }

    // C = D_dam @ C_eff
    double C[6][6];
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
        double v = 0.0; for (int k = 0; k < 6; ++k) v += Ddam[i][k] * Ceff[k][j];
        C[i][j] = v;
    }

    // --- chain-rule gradient pieces (each d(.)/dε, 6-vector) ---
    // Ceff^T @ g  (d(scalar of sig_eff)/dε = (d sig_eff/dε)^T @ (d scalar/d sig_eff))
    auto CeffT = [&](const double g[6], double out[6]) {
        for (int k = 0; k < 6; ++k) { double v = 0.0; for (int i = 0; i < 6; ++i) v += Ceff[i][k] * g[i]; out[k] = v; }
    };
    double g_et[6], g_xs[6], g_ac[6], det_deps[6], dxs_deps[6], dac_deps[6];
    dscalarDsig(0, sig_eff, mp, g_et);   CeffT(g_et, det_deps);
    dscalarDsig(1, sig_eff, mp, g_xs);   CeffT(g_xs, dxs_deps);
    dscalarDsig(2, sig_eff, mp, g_ac);   CeffT(g_ac, dac_deps);

    // eigenprojection gradients (WITH W6 tensor weight): dDt/dε = Ceff^T (W6 . Emax), dDc/dε = -Ceff^T (W6 . Emin)
    double dDt_deps[6], dDc_deps[6];
    {
        double Em[3][3], En[3][3], Em6[6], En6[6], tmp[6];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
            Em[i][j] = V[i][imax] * V[j][imax];
            En[i][j] = V[i][imin] * V[j][imin];
        }
        matToVoigt(Em, Em6); matToVoigt(En, En6);
        for (int i = 0; i < 6; ++i) tmp[i] = W6[i] * Em6[i];
        if (Dt > 0.0) CeffT(tmp, dDt_deps); else for (int i = 0; i < 6; ++i) dDt_deps[i] = 0.0;
        for (int i = 0; i < 6; ++i) tmp[i] = W6[i] * En6[i];
        if (Dc > 0.0) { CeffT(tmp, dDc_deps); for (int i = 0; i < 6; ++i) dDc_deps[i] = -dDc_deps[i]; }
        else for (int i = 0; i < 6; ++i) dDc_deps[i] = 0.0;
    }

    // ||Deps_p|| gradient: depl_deps = I - C^-1 Ceff ; dnorm/dε = depl_deps^T (W6 . depl) / dnorm
    double dnorm_deps[6];
    if (dnorm > 1.0e-14) {
        double Cel[6][6], Cinv[6][6], CinvCeff[6][6];
        elasticC(mp, Cel); invert6(Cel, Cinv);
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            double v = 0.0; for (int k = 0; k < 6; ++k) v += Cinv[i][k] * Ceff[k][j];
            CinvCeff[i][j] = v;
        }
        double x[6]; for (int i = 0; i < 6; ++i) x[i] = W6[i] * depl[i];
        for (int k = 0; k < 6; ++k) {
            double v = 0.0;
            for (int i = 0; i < 6; ++i) {
                const double depl_deps_ik = ((i == k) ? 1.0 : 0.0) - CinvCeff[i][k];
                v += depl_deps_ik * x[i];
            }
            dnorm_deps[k] = v / dnorm;
        }
    } else for (int i = 0; i < 6; ++i) dnorm_deps[i] = 0.0;

    // dkd*/dε under loading (else zero)
    double dkdt1[6], dkdt2[6], dkdc1[6], dkdc2[6];
    if (loading) {
        const double ixs = 1.0 / xs, ixs2 = ixs * ixs;
        for (int i = 0; i < 6; ++i) {
            dkdt2[i] = det_deps[i] * ixs - above * dxs_deps[i] * ixs2;
            dkdt1[i] = dnorm_deps[i] * ixs - dnorm * dxs_deps[i] * ixs2;
            dkdc2[i] = (dac_deps[i] * above + ac * det_deps[i]) * ixs - ac * above * dxs_deps[i] * ixs2;
            dkdc1[i] = (dac_deps[i] * dnorm + ac * dnorm_deps[i]) * ixs - ac * dnorm * dxs_deps[i] * ixs2;
        }
    } else for (int i = 0; i < 6; ++i) { dkdt1[i] = dkdt2[i] = dkdc1[i] = dkdc2[i] = 0.0; }

    // ω via IFT (only when interior 0<ω<1; clamped/inactive ω is insensitive => dω=0)
    double dwt[6], dwc[6];
    if (wt > 0.0 && wt < 1.0) {
        const double Ht = Dt * ((1.0 - wt) * kdt2 / eps_f - 1.0);
        const double a0 = -(1.0 - wt) / Ht;
        const double a1 = -(1.0 - wt) * Dt / (eps_f * Ht);
        const double a2 = -(1.0 - wt) * Dt * wt / (eps_f * Ht);
        for (int i = 0; i < 6; ++i) dwt[i] = a0 * dDt_deps[i] + a1 * dkdt1[i] + a2 * dkdt2[i];
    } else for (int i = 0; i < 6; ++i) dwt[i] = 0.0;
    if (wc > 0.0 && wc < 1.0) {
        const double Hc = Dc * ((1.0 - wc) * kdc2 / eps_fc - 1.0);
        const double a0 = -(1.0 - wc) / Hc;
        const double a1 = -(1.0 - wc) * Dc / (eps_fc * Hc);
        const double a2 = -(1.0 - wc) * Dc * wc / (eps_fc * Hc);
        for (int i = 0; i < 6; ++i) dwc[i] = a0 * dDc_deps[i] + a1 * dkdc1[i] + a2 * dkdc2[i];
    } else for (int i = 0; i < 6; ++i) dwc[i] = 0.0;

    // assemble: C - sig_t (x) dω_t - sig_c (x) dω_c
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
        D6[i][j] = C[i][j] - sig_t[i] * dwt[j] - sig_c[i] * dwc[j];
}

// ---------------------------------------------------------------------------
// Atomic incremental return (mirrors the oracle return_map_tensor exactly):
// elastic-predict sig_tr = sig_n + C:deps -> eigendecompose -> principal return
// (radial, preserves V) -> recompose. Consistent tangent on request.
//   hardening=true  : full CDPM2 (qh1/qh2/kp).  false : perfect-plastic failure surface.
// ---------------------------------------------------------------------------
inline int returnMapTensor(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                           bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                           bool doTangent)
{
    double sig_tr[6];
    elasticPredTensor(sig_n, deps, mp, sig_tr);
    double A[3][3], w[3], V[3][3];
    voigtToMat(sig_tr, A);
    eig3sym(A, w, V);

    PrincipalResult pr = hardening ? returnMapHardening(w, mp, kp_n) : returnMapPrincipal(w, mp);
    kp_new = hardening ? pr.kp : kp_n;

    // recompose sig = V diag(sp) V^T (radial return preserves V)
    double S[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = 0.0;
        for (int a = 0; a < 3; ++a) v += V[i][a] * pr.sp[a] * V[j][a];
        S[i][j] = v;
    }
    matToVoigt(S, sig_new);

    // SAFETY NET (PR #249 review): never propagate a non-finite stress into the global solve.
    // Catches the NaN/Inf that degenerate params produce (e.g. ft=0/fc=ft => m0<=0 => the apex
    // xi=sqrt3 fc/m0 blows up). The kernel ASSUMES valid params (fc>0, ft>0, 0<=nu<0.5, m0>0) —
    // the nDMaterial wrapper enforces that (handoff §5); this is the last-line defensive guard.
    bool finite = true;
    for (int i = 0; i < 6; ++i) if (!std::isfinite(sig_new[i])) finite = false;
    if (!finite) {
        for (int i = 0; i < 6; ++i) sig_new[i] = sig_tr[i];   // elastic predictor (best available)
        kp_new = kp_n;
        if (doTangent) elasticC(mp, Dtan6);
        return 2;
    }

    if (doTangent) consistentTangent(sig_tr, w, V, pr, mp, hardening, Dtan6);
    return pr.converged ? 0 : 2;   // 0 OK, 2 = no-converge (honest flag)
}

// ===========================================================================
// Public small-strain entry point. Computes the TRIAL response (stress, tangent)
// from the COMMITTED history `in`; writes the new (uncommitted) state to `out`.
// The caller commits by copying out->in on commitState.
//   sigEffImplicit = the UNDAMAGED effective stress (== sigma at P1; the LogStrain
//     b^e fix per ADR R3 — kept separate so a future damage/IMPL-EX choice never
//     corrupts the finite-strain recovery).
//   dt : current time increment (ops_Dt). Only used for the Tier-2 IMPL-EX extrapolation ratio
//        r = dt/dt_n; ignored for Tier-1 (mp.implex == false).
// Returns: 0 converged, 2 no-converge (honest off-surface flag).
//
// TIER-1 (mp.implex == false): reports the IMPLICIT nominal stress + the analytic dual-projector
//   damaged tangent (P3b), exactly as before.
// TIER-2 (mp.implex == true, ADR §4.4): ALSO runs the implicit solve (for the committed internal
//   variables + the next-step extrapolation source), but REPORTS the EXPLICIT stress assembled from
//   EXTRAPOLATED internals (the plastic-strain increment + the dual damage frozen at the committed
//   rate, ratio clamped to [0, implexRmax]) and the degraded-elastic SECANT D_dam(w~):C0. The secant
//   is symmetric-part SPD only in SINGLE-SIGN principal regimes (review NUM-1) — it is still the exact
//   d(sigma)/d(strain) of the reported explicit stress. sigEffImplicit stays the IMPLICIT effective
//   stress regardless of tier (the LogStrain b^e contract, ADR R3).
// ===========================================================================
inline int returnMap(const Params& mp, const double strain[6], const State& in, State& out,
                     double sigma[6], double sigEffImplicit[6], double Dtan6[6][6],
                     bool doTangent, double dt = 0.0, bool hardening = true,
                     double* wtOut = nullptr, double* wcOut = nullptr)
{
    double deps[6];
    for (int i = 0; i < 6; ++i) deps[i] = strain[i] - in.eps[i];
    // (1) IMPLICIT EFFECTIVE-stress return from the committed EFFECTIVE state (NOT the nominal sig).
    double sig_eff[6], kp_new;
    int status = returnMapTensor(mp, in.sigEff, deps, in.kp, hardening, sig_eff, kp_new, Dtan6, doTangent);
    for (int i = 0; i < 6; ++i) { out.eps[i] = strain[i]; out.sigEff[i] = sig_eff[i]; sigEffImplicit[i] = sig_eff[i]; }
    out.kp = kp_new;
    // (2) IMPLICIT P2 dual-damage NOMINAL stress (writes out.sig + the damage history). Unilateral by re-split.
    double wt_impl = 0.0, wc_impl = 0.0;
    damagedUpdate(mp, in, sig_eff, strain, out, &wt_impl, &wc_impl);
    // carry the committed IMPLICIT damage + the per-variable increments for the next IMPL-EX extrapolation
    out.wt = wt_impl; out.wc = wc_impl;
    out.dwt = wt_impl - in.wt; out.dwc = wc_impl - in.wc;
    { double epl[6], epl_n[6];
      plasticStrain6(out.sigEff, out.eps, mp, epl);
      plasticStrain6(in.sigEff,  in.eps,  mp, epl_n);
      for (int i = 0; i < 6; ++i) out.depl[i] = epl[i] - epl_n[i]; }
    out.dt_n = (dt > 0.0) ? dt : in.dt_n;
    if (wtOut) *wtOut = wt_impl;   // recorders show the IMPLICIT (accurate) damage in both tiers
    if (wcOut) *wcOut = wc_impl;

    if (!mp.implex) {
        // ---- TIER-1: report the implicit nominal + the analytic damaged tangent (P3b) ----
        for (int i = 0; i < 6; ++i) sigma[i] = out.sig[i];
        if (doTangent && status == 0) {
            double Ceff[6][6];
            for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Ceff[i][j] = Dtan6[i][j];
            damagedTangent(mp, in, sig_eff, strain, Ceff, Dtan6);
        }
        return status;
    }

    // ---- TIER-2 IMPL-EX: report the EXPLICIT extrapolated stress + the degraded-elastic SECANT ----
    double r = (in.dt_n > 0.0 && dt > 0.0) ? (dt / in.dt_n) : 0.0;   // forward-only; floored at 0
    if (r > mp.implexRmax) r = mp.implexRmax;                        // clamp (review ALG-2/NUM-2/NUM-3)
    double wt_x = in.wt + r * in.dwt; if (wt_x < 0.0) wt_x = 0.0; if (wt_x > 1.0 - 1.0e-12) wt_x = 1.0 - 1.0e-12;
    double wc_x = in.wc + r * in.dwc; if (wc_x < 0.0) wc_x = 0.0; if (wc_x > 1.0 - 1.0e-12) wc_x = 1.0 - 1.0e-12;
    double deps_eff[6];
    for (int i = 0; i < 6; ++i) deps_eff[i] = deps[i] - r * in.depl[i];   // frozen plastic-strain increment
    double sig_bar_x[6];
    elasticPredTensor(in.sigEff, deps_eff, mp, sig_bar_x);               // LINEAR in deps => elastic tangent
    double A[3][3], w[3], V[3][3]; voigtToMat(sig_bar_x, A); eig3sym(A, w, V);
    double sp[3];
    for (int i = 0; i < 3; ++i) {
        const double st = w[i] > 0.0 ? w[i] : 0.0, sc = w[i] < 0.0 ? w[i] : 0.0;
        sp[i] = (1.0 - wt_x) * st + (1.0 - wc_x) * sc;                   // Eq.1 with FROZEN damage
    }
    double S[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = 0.0; for (int a = 0; a < 3; ++a) v += V[i][a] * sp[a] * V[j][a];
        S[i][j] = v;
    }
    matToVoigt(S, sigma);
    if (doTangent) {
        double yv[3], ypv[3];
        for (int i = 0; i < 3; ++i) {
            yv[i] = (1.0 - wt_x) * (w[i] > 0.0 ? w[i] : 0.0) + (1.0 - wc_x) * (w[i] < 0.0 ? w[i] : 0.0);
            ypv[i] = (w[i] > 0.0) ? (1.0 - wt_x) : (1.0 - wc_x);
        }
        double Ddam[6][6], C0[6][6]; isotropicTangent(w, V, yv, ypv, Ddam); elasticC(mp, C0);
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            double v = 0.0; for (int k = 0; k < 6; ++k) v += Ddam[i][k] * C0[k][j];
            Dtan6[i][j] = v;                                             // D_dam(w~) : C0
        }
    }
    return status;
}

} // namespace Concrete3D
} // namespace Ladruno

#endif // LadrunoConcrete3DKernel_h
