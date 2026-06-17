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
// STATUS: P0/P1 SURFACE + HARDENING LAWS. The CDPM2 yield surface (yieldF, Eq.18 with the [1-qh1]
// cap) and the hardening building blocks (qh1Of/qh2Of Eq.30-31, ductilityXh Eq.33-36) are
// IMPLEMENTED and mirror the numpy oracle 1:1 (gated: surface normalization, Kupfer fcc/fc, the
// hardening unit identities, and the oracle's 4-unknown hardening return map). The C++ return map
// (the actual stress update, P1-build), dual damage (P2), robustness tiers (P3), finite-strain
// out-contract (P4), and confined-fiber condensation (§4.6) are STUBS with ADR-phase markers.
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
    // hardening (P1, CDPM2 Eqs. 30-36)
    double qh0 = 0.3;   // initial yield fraction (qh1 at kp=0)
    double Hp = 0.5;    // hardening modulus (qh1 slope at kp=1^-, qh2 slope for kp>1)
    double Ah = 0.08;   // ductility-measure params (Eq.33; literature defaults, recalibrate)
    double Bh = 0.003;
    double Ch = 2.0;
    double Dh = 1.0e-6;
    // rate / robustness
    double eta = 0.0;            // Duvaut-Lions viscosity (0 => inviscid, byte-identical)
    bool   implex = false;       // Tier-2
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

// Convenience: eccentricity that makes the equibiaxial strength hit fcc/fc=target (ADR 4.1b).
// Mirror of the oracle's bisection; the wrapper calls this once at setup.
// Ladruno-NOTE: DECLARATION-ONLY at P0/P1 (undefined symbol). A TU that merely uses the surface
// (yieldF/m0Of/invariants/lodeR) compiles+links fine; any TU that CALLS this fails to link until
// the build PR provides the definition (.cpp mirror of the oracle bisection) — that is intended,
// not a regression.
double eccentricityFromKupfer(double fc, double ft, double targetFccRatio /* =1.16 */);

// ===========================================================================
// STUBS — implemented in later phases. Signatures fix the kernel contract now.
// ===========================================================================

// P1 — semi-implicit Menetrey-Willam return map (non-associated flow + ductility + apex handler).
// out: sigma (nominal), sigEffImplicit (UNDAMAGED, true backward-Euler — consumed by LogStrain
//      so IMPL-EX/Tier choice never corrupts the finite-strain b^e recovery — ADR §4.4 BLOCKING),
//      Dtan6 (6x6, NON-SYMMETRIC for non-associated flow => Tier-1 needs an unsymmetric solver).
// dScaleOverride: <0 => Tier-1 implicit; >=0 => Tier-2 IMPL-EX (freeze plastic state + damage).
struct State; // committed history (plastic strain, kappa_p, kappa_dt/dc, projector frame, dtime)
int returnMap(const Params& mp, const double strain[6], State& st,
              double sigma[6], double sigEffImplicit[6], double Dtan6[6][6],
              bool doTangent, double dScaleOverride /* <0 = implicit */);

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

} // namespace Concrete3D
} // namespace Ladruno

#endif // LadrunoConcrete3DKernel_h
