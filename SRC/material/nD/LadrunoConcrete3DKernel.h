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
    // P2h compression->tension damage-coupling temper (the -ctTemper modes): scales the tensile
    // damage-history (kdt1/kdt2) accumulation by a weight w_t. 0=none (w_t=1, literal CDPM2,
    // byte-identical); 1=alphat (w_t=1-alpha_c); 2=proj (fraction of the plastic-strain increment along
    // POSITIVE effective-stress directions). Both temper modes => w_t~0 in compression (no tension
    // pre-damage). See tensileDamageWeight.
    int    ctTemper = 0;
    // PV20 tension->compression damage-coupling temper (the -tcTemper modes; mirror of ctTemper). 0 = none (literal
    // CDPM2 Eq.47/48, byte-identical; the kernel/oracle default). 2 = proj — the nDMaterial DEFAULT: the kdc1 plastic
    // measure is the part of ||d eps_p|| along NON-tensile effective-stress directions, and the CDPM2 compressive drive
    // eqc is fed by eps_tilde of <sig_bar>- (tensile principals zeroed). Identical to none whenever no effective
    // principal is tensile. See compressiveDamageWeight / compressiveEquivStrain.
    int    tcTemper = 0;
    // Tensile softening law (WP concrete3d-oracle-diagnosis). 0 = legacy exponential sigma = ft exp(-eps_i/eps_f),
    // eps_f = Gf/(ft lch), with the legacy kdt2 = (kdt - eps0)/x_s history (byte-identical to the pre-2026-09
    // kernel; the kernel/oracle default so every fixture stays pinned). 1 = CDPM2 BILINEAR (Grassl 2013
    // Eq.51/58/59: s1 = 0.3 ft at wf1 = 0.15 wf, wf = 4.444 Gf/ft, w = lch*eps_i) with the literal Eq.44/45
    // histories — the nDMaterial wrapper DEFAULT (-tensionLaw bilinear|exp).
    int    tensionLaw = 0;
    // Compressive softening strain eps_fc used DIRECTLY when > 0 (the wrapper's -epsFc, or its Gc-energy
    // calibration calibrateEpsFc). 0 => legacy eps_fc = Gc/(fc lch).
    double epsFc = 0.0;
    // Plastic potential (WP concrete3d-flow-potential, B1). 0 = legacy v1 flow (m_v = Df m0/(sqrt3 fc), qh1=1-shaped
    // m_s — always dilatant; the kernel/oracle default so every fixture stays pinned). 1 = the FULL CDPM2 potential
    // (Grassl 2013 Eq.22-29: [1-qh1] cap + m_g(sigV,kp), Df = CDPM2 dilation constant) — the nDMaterial DEFAULT
    // (-flowPotential cdpm2|legacy).
    int    flowPotential = 0;
    // Return-map sub-incrementation depth (B1): if the direct return fails, halve the strain increment down to
    // 2^-maxSubIncr (OOFEM performPlasticityReturn). 0 => the direct return only (byte-identical legacy).
    int    maxSubIncr = 0;
    // TOTAL attempt budget for the sub-increment loop in returnMapTensor (successes + failures), WP
    // concrete3d-hang-diagnosis: a GP that only ever succeeds near the 2^-maxSubIncr floor alternates
    // fail/succeed and can legally spend up to ~2 * 2^maxSubIncr attempts (of up to 100 Newton iterations
    // each) reaching done=1.0 -- at the CDPM2 wrapper default maxSubIncr=10 that is ~2048 attempts per
    // material point per global Newton iteration, which is the multi-hour analyze(1) hang (nothing ever
    // fails, so nothing is logged and the step is never cut). Independent of maxSubIncr; only matters when
    // maxSubIncr > 0. Default 64 bounds the worst case to 64 * 100 Newton iterations per call.
    int    maxSubAttempts = 64;
    // Sub-incrementation MODE (only when maxSubIncr > 0, hardening map). 0 = DETERMINISTIC (DEFAULT, WP
    // concrete3d-hang-diagnosis #877 follow-up): the piece count n = clamp(ceil(f_trial / subIncrC), 1,
    // subIncrMaxPieces) is a function of the trial overshoot only and is ALWAYS applied, with a deterministic
    // ladder n -> 2n -> 4n on a piece failure; refuse only after the ladder. 1 = ADAPTIVE: the direct return first,
    // then halving/doubling on failure (bounded by maxSubAttempts) -- the pre-follow-up path, opt-in so old numbers
    // stay reproducible (nDMaterial -subIncr adaptive).
    int    subIncrMode = 0;
    double subIncrC = 0.3;
    int    subIncrMaxPieces = 64;
    int    subIncrForceN = 0;   // > 0 pins the deterministic piece count n (finite-difference legs of the algorithmic tangent use the CENTRAL n)
    int    subIncrRescue = 1;   // deterministic mode only: after the ladder fails, try the ADAPTIVE path (bounded by maxSubAttempts) before refusing
    // DEAD POINTS (WP concrete3d-hang-diagnosis, owner decision 2026-09-28): the residual-strength fraction (1 - omega) at or
    // below which a point's crack (omega_t) or crush (omega_c) is treated as fully open. Committed damage >= omegaDead:
    //   * omega_t: TENSION CUTOFF ON THE PLASTIC FLOW -- the tensile spectral part of the trial effective stress is carried
    //     ELASTICALLY (no flow, no kappa_p growth from tension; returns to zero at eps_p on unloading) and the return map
    //     runs on the compressive remainder only (a cracked point keeps its compressive strut);
    //   * omega_c: STRICT FREEZE -- a crushed point carries nothing: kappa_p, the plastic strain and the damage histories
    //     are frozen, the effective stress is elastic on the fixed plastic strain, both damages go to the floor OMEGA_MAX,
    //     nominal = (1-OMEGA_MAX)*sig_eff, tangent = (1-OMEGA_MAX)*C (the jump to the floor costs at most
    //     (1-omegaDead)*|sig_eff| of nominal stress; freezing at the committed omega instead was measured to break the
    //     Gc-calibration gate, 62 % off, because the residual then grows with the elastic sig_eff).
    // Without it kappa_p and sig_eff run away on a point that has already lost its strength (K&R coarse: kappa_p 3.4e4,
    // sig_eff 813 MPa at fc = 24; G5 band: kappa_p 3.2e4, sig_eff 779 MPa at omega_t = 0.9993) until the return map cannot
    // integrate them; the nominal residual (1-omega)*sig_eff is then a spurious fraction of ft. omegaDead >= 1 disables
    // (parser: -deadThreshold takes [0.99, 1); -noDead sets 2.0, the A/B knob that reproduces the pre-treatment behaviour).
    // The default is the smallest threshold that clears the measured runaway states (ExplicitBathe tension softening on
    // the legacy exp law only reaches omega_t = 0.99856 when its return map refuses). Parser flag -deadThreshold.
    double omegaDead = 0.998;
    // Compressive damage drive (B2, WP concrete3d-damage-drive). 0 = legacy fork drive ((1-wc)(-sig_min) = fc exp(..),
    // histories from the onset only; the kernel/oracle default). 1 = CDPM2 Eq.47-49/53/55 (OOFEM computeDamage /
    // computeDamageParamCompression): eqc += alpha_c d(eps_tilde), kappa_dc = max eqc, kdc2 from the start, kdc1 with
    // the post-onset fraction, (1-wc) E kappa_dc = ft exp(..) — the nDMaterial DEFAULT (-compressionDrive cdpm2|legacy).
    int    compDrive = 0;
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
// B5: units-free acceptance tolerance on the (dimensionless) yield function after a return. Was 1e-7*(fc+1):
// 3.1e-6 in MPa (fc = 30) but 0.3 in Pa (fc = 3e6), so an SI model accepted returns 30 % off the surface. The
// value equals the MPa one at fc = 30 (every fixture unchanged).
static const double F_TOL_HONEST = 3.1e-6;

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

inline void voigtToMat(const double v[6], double M[3][3]);   // fwd decl (defined below; used by 'proj')

// P2h compression->tension damage-coupling TEMPER weight w_t (mirror the oracle tensile_damage_weight):
// scales the tensile damage-history accumulation. 0=none (w_t=1, byte-identical); 1=alphat (w_t=1-alpha_c);
// 2=proj (fraction of the plastic-strain increment depl[6] along POSITIVE effective-stress principal
// directions: project depl into the stress eigenframe V, keep the diagonal entries where w_stress>0).
inline double tensileDamageWeight(const Params& mp, double ac, const double depl6[6],
                                  const double w_stress[3], const double V[3][3])
{
    if (mp.ctTemper == 1) {                                  // alphat
        const double w = 1.0 - ac;
        return w > 0.0 ? w : 0.0;
    }
    if (mp.ctTemper == 2) {                                  // proj
        double M[3][3]; voigtToMat(depl6, M);
        double nrm2 = 0.0;
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) nrm2 += M[i][j] * M[i][j];
        if (nrm2 <= 1.0e-300) return 1.0;
        const double floor = 1.0e-6 * mp.ft;
        double tens = 0.0;
        for (int a = 0; a < 3; ++a) {
            if (w_stress[a] <= floor) continue;             // tensile-stress directions only
            double d = 0.0;                                  // (V^T M V)_aa = plastic strain along axis a
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) d += V[i][a] * M[i][j] * V[j][a];
            tens += d * d;
        }
        return std::sqrt(tens) / std::sqrt(nrm2);
    }
    return 1.0;                                              // none (literal CDPM2)
}

// PV20 -tcTemper proj (mirror the oracle compressive_damage_weight): kdc1 (Eq.48) weight w_c = ||Phi depl Phi|| /
// ||depl|| in the effective-stress eigenframe, phi_a = tcPhi (1 on non-tensile directions incl. a 1e-3 ft dead zone,
// 0 on a crack direction; continuous). EXACTLY 1.0 when no effective principal is > TC_DEAD ft (every compressive backbone, the Gc
// table). Literal CDPM2 counts the crack-opening plastic strain of a cracked-then-compressed state (RC panel in shear)
// as crushing history: the strut softened at ~0.2 fc in PV20 (tau 1.76 -> 0.03 MPa vs 4.26 in the test).
// phi_a = 1 (sigma_bar_a <= TC_DEAD ft) .. 0 (sigma_bar_a >= TC_BAND ft), linear between (oracle _tc_phi)
static const double TC_DEAD = 1.0e-3, TC_BAND = 0.05;
inline double tcPhi(const Params& mp, double s)
{
    const double v = 1.0 - (s - TC_DEAD * mp.ft) / ((TC_BAND - TC_DEAD) * mp.ft);
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}
inline double compressiveDamageWeight(const Params& mp, const double depl6[6],
                                      const double w_stress[3], const double V[3][3])
{
    if (mp.tcTemper != 2) return 1.0;
    double mx = w_stress[0]; for (int a = 1; a < 3; ++a) if (w_stress[a] > mx) mx = w_stress[a];
    if (mx <= TC_DEAD * mp.ft) return 1.0;
    double M[3][3]; voigtToMat(depl6, M);
    double nrm2 = 0.0;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) nrm2 += M[i][j] * M[i][j];
    if (nrm2 <= 1.0e-300) return 1.0;
    double phi[3];
    for (int a = 0; a < 3; ++a) phi[a] = tcPhi(mp, w_stress[a]);
    double p2 = 0.0;                                         // ||Phi (V^T M V) Phi||_F^2
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
        double d = 0.0;
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) d += V[i][a] * M[i][j] * V[j][b];
        const double q = phi[a] * d * phi[b];
        p2 += q * q;
    }
    return std::sqrt(p2) / std::sqrt(nrm2);
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

// PV20 -tcTemper proj (mirror the oracle compressive_equiv_strain): the equivalent strain feeding the CDPM2
// compressive history (Eq.47) is eps_tilde of the COMPRESSIVE part (tensile principals scaled by tcPhi -> zeroed on
// a crack); `et` (the full eps_tilde) otherwise and whenever no effective principal is > TC_DEAD ft. Stops the hardened crack stress (qh2 ft, qh2 ~ 3 across a crack)
// from pushing kappa_dc past eps0 while the strut is at ~0.2 fc.
inline double compressiveEquivStrain(const Params& mp, const double w_stress[3], double et)
{
    if (mp.tcTemper != 2) return et;
    double mx = w_stress[0]; for (int a = 1; a < 3; ++a) if (w_stress[a] > mx) mx = w_stress[a];
    if (mx <= TC_DEAD * mp.ft) return et;
    double wc[3];
    for (int a = 0; a < 3; ++a) wc[a] = w_stress[a] < 0.0 ? w_stress[a] : tcPhi(mp, w_stress[a]) * w_stress[a];
    return equivStrainGeneral(wc, mp);
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

// ---------------------------------------------------------------------------
// CDPM2 BILINEAR tension law + literal Eq.45 history (mirror of the oracle _bilinear_sigma /
// _solve_omega_bilinear / _omega_t / _tension_hist_update / _eps_fc byte-for-byte). See Params::tensionLaw.
// ---------------------------------------------------------------------------
static const double BILIN_S1 = 0.3, BILIN_W1 = 0.15;
static const double OMEGA_TAN_FLOOR = 1.0e-6;                // residual (1-omega) in the damaged TANGENT only
// B2: omega is capped at OMEGA_MAX for the STRESS too. With omega = 1 exactly (bilinear law past wf) every fully
// cracked state carries IDENTICALLY zero stress, i.e. a whole plateau of far states is an exact root of any
// equilibrium (OOFEM con2dpm2 at one sub-step: the element Newton jumped at step 5 from the physical lateral root
// +2.04e-3 to a fully cracked +1.21e-2 with sigma == 0). The 1e-6 residual of the effective stress removes the
// plateau root. Only omega > 1 - 1e-6 is affected (every fixture regenerates byte-identically).
static const double OMEGA_MAX = 1.0 - 1.0e-6;
static const double BILIN_GF = 0.5 * (0.15 + 0.3);           // Gf/(ft wf) = 0.225  =>  wf = 4.444 Gf/ft

inline double epsFcOf(const Params& mp) { return mp.epsFc > 0.0 ? mp.epsFc : mp.Gc / (mp.fc * mp.lch); }

inline double bilinearSigma(double wcr, double ft, double wf, double* slope = nullptr)
{
    const double wf1 = BILIN_W1 * wf, s1 = BILIN_S1 * ft;
    if (wcr <= wf1) { if (slope) *slope = -(ft - s1) / wf1; return ft - (ft - s1) * wcr / wf1; }
    if (wcr <= wf)  { if (slope) *slope = -s1 / (wf - wf1);  return s1 * (wf - wcr) / (wf - wf1); }
    if (slope) *slope = 0.0;
    return 0.0;
}

inline double solveOmegaBilinear(double kd1, double kd2, double D, double ft, double Gf, double lch)
{
    const double wf = Gf / (BILIN_GF * ft);
    const double wf1 = BILIN_W1 * wf, s1 = BILIN_S1 * ft, h = lch;
    if (D <= ft) return 0.0;
    double den = D * wf1 + (s1 - ft) * kd2 * h;
    if (den > 0.0) {
        const double w = ((D - ft) * wf1 - (s1 - ft) * kd1 * h) / den;
        if (0.0 <= w && w <= 1.0 && h * (kd1 + w * kd2) <= wf1) return w;
    }
    den = D * (wf - wf1) - s1 * kd2 * h;
    if (den > 0.0) {
        const double w = (D * (wf - wf1) + s1 * (kd1 * h - wf)) / den;
        const double wi = h * (kd1 + w * kd2);
        if (0.0 <= w && w <= 1.0 && wf1 < wi && wi <= wf) return w;
    }
    if (h * kd1 >= wf) return 1.0;
    auto F = [&](double w) { return (1.0 - w) * D - bilinearSigma(h * (kd1 + w * kd2), ft, wf); };
    double lo = 0.0, hi = 1.0;
    if (F(lo) <= 0.0) return 0.0;
    if (F(hi) >= 0.0) return 1.0;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) break;
        if (F(mid) > 0.0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

inline double omegaT(const Params& mp, double kdt1, double kdt2, double D)
{
    if (mp.tensionLaw != 1) return solveOmegaBracketed(kdt1, kdt2, D, mp.ft, mp.Gf / (mp.ft * mp.lch));   // legacy: uncapped
    const double w = solveOmegaBilinear(kdt1, kdt2, D, mp.ft, mp.Gf, mp.lch);
    return w < OMEGA_MAX ? w : OMEGA_MAX;
}

// Compressive damage for the active drive (mirror of the oracle _omega_c); <= OMEGA_MAX.
inline double omegaC(const Params& mp, double kdc, double kdc1, double kdc2, double sigcMax, double eps_fc)
{
    if (mp.compDrive != 1)                                   // legacy drive: uncapped (byte-identical)
        return (kdc > 0.0 && sigcMax > 1.0e-6 * mp.fc) ? solveOmegaBracketed(kdc1, kdc2, sigcMax, mp.fc, eps_fc) : 0.0;
    const double w = (kdc > mp.ft / mp.E) ? solveOmegaBracketed(kdc1, kdc2, mp.E * kdc, mp.ft, eps_fc) : 0.0;
    return w < OMEGA_MAX ? w : OMEGA_MAX;
}

// CDPM2 compressive histories after one step (B2; mirror of the oracle _comp_hist_update_cdpm2).
inline void compHistUpdateCdpm2(const Params& mp, double& kdc, double& kdc1, double& kdc2, double& eqc,
                                double etp, double et, double ac, double bc, double dnorm, double xs)
{
    const double eps0 = mp.ft / mp.E;
    const double eqcNew = eqc + ac * (et - etp);
    if (eqcNew > kdc) {
        const double d = eqcNew - kdc;
        if (eqcNew > eps0) {
            const double frac = (kdc >= eps0) ? 1.0 : (eqcNew - eps0) / d;
            kdc1 += ac * bc * frac * dnorm / xs;
        }
        kdc2 += d / xs;
        kdc = eqcNew;
    }
    eqc = eqcNew;
}

inline void tensionHistUpdate(const Params& mp, double& kdt1, double& kdt2, double et, double et_max_n,
                              double dnorm, double xs, double wtw)
{
    const double eps0 = mp.ft / mp.E;
    const double det_raw = (et - et_max_n > 0.0) ? (et - et_max_n) : 0.0;
    if (mp.tensionLaw == 1) {
        if (det_raw > 0.0) {
            kdt2 += wtw * det_raw / xs;
            if (et > eps0) {
                const double frac = (et_max_n >= eps0) ? 1.0 : (et - eps0) / det_raw;
                kdt1 += wtw * frac * dnorm / xs;
            }
        }
        return;
    }
    if (det_raw > 0.0 && et > eps0) {
        const double lo = et_max_n > eps0 ? et_max_n : eps0;
        const double above = (et - lo > 0.0) ? (et - lo) : 0.0;
        kdt2 += wtw * above / xs;
        kdt1 += wtw * dnorm / xs;
    }
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
// VERTEX (hydrostatic-axis) return — mirrors the oracle _vertex_kappa / cdpm2_vertex_potential_grad /
// return_map_vertex byte-for-byte (same bisection, same acceptance test).
//   vertexKappa: kp consistent with the vertex projection (Grassl Eq.32): kp_n + ||d eps_p||/xh(sigV),
//     ||d eps_p||_F = sqrt((sigV_tr-sigV)^2/(3K^2) + (rho_tr/2G)^2) (volumetric + deviatoric, the same
//     Frobenius norm as the regular map's R4), theta=pi/3 on the axis => (2cos)^2 = 1 (OOFEM's choice).
//     NB OOFEM computeTempKappa uses 1/9 instead of 1/3 on the volumetric term (1/sqrt3 of the norm its
//     own regular return uses) — the one deliberate difference to OOFEM on the axis.
//   cdpm2VertexPotentialGrad: (dg/dsigV, dg/drho) of the FULL CDPM2 potential (Eq.22-29) at rho=0; used
//     ONLY for the compressive-vertex cone test (the v1 regular flow m_v=Df m0/(sqrt3 fc) has no cap
//     term => its volumetric flow is always dilatant => it has NO compressive cone). Df clamped >0.5.
//   returnMapVertex: F(sigV)=f_p(sigV,0;kp(sigV))=0 by bisection in (sigV_tr,0); the yield function is
//     regular on the axis (Lode r multiplies rho=0). Accept iff dlam=(sigV_tr-sigV)/(K dg/dsigV)>=0 and
//     rho_tr <= 2G dlam dg/drho (trial inside the cone of normals).
// ---------------------------------------------------------------------------
inline double vertexKappa(double kp_n, double sigV_tr, double rho_tr, double sigV, const Params& mp)
{
    const double K = bulkK(mp), G = shearG(mp);
    const double d = sigV_tr - sigV, g = rho_tr / (2.0 * G);
    const double eq = std::sqrt(d * d / (3.0 * K * K) + g * g);
    return kp_n + eq / ductilityXh(sigV, mp.fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
}

inline void cdpm2VertexPotentialGrad(double sigV, double kp, const Params& mp, double& dgs, double& dgr)
{
    const double fc = mp.fc, ft = mp.ft, m0 = mp.m0;
    const double Df = mp.Df > 0.5 + 1.0e-6 ? mp.Df : 0.5 + 1.0e-6;
    const double q1 = qh1Of(kp, mp.qh0, mp.Hp), q2 = qh2Of(kp, mp.Hp);
    const double AG = 3.0 * ft * q2 / fc + m0 / 2.0;
    const double BG = q2 / 3.0 * (1.0 + ft / fc) / (std::log(AG) + std::log(Df + 1.0) - std::log(2.0 * Df - 1.0)
                                                  - std::log(3.0 * q2 + m0 / 2.0));
    const double mQ = AG * std::exp((sigV - ft * q2 / 3.0) / (fc * BG));
    const double Bl = sigV / fc;
    const double Al = (1.0 - q1) * Bl * Bl;
    dgs = 4.0 * (1.0 - q1) / fc * Al * Bl + q1 * q1 * mQ / fc;
    dgr = Al / (SQRT6 * fc) * (4.0 * (1.0 - q1) * Bl + 6.0) + m0 * q1 * q1 / (SQRT6 * fc);
}

// ---------------------------------------------------------------------------
// FULL CDPM2 PLASTIC POTENTIAL (B1; mirror of the oracle cdpm2_potential_derivs / flow_grad_jac). Grassl 2013
// Eq.22-29 in the OOFEM ConcreteDPM2 form (sigV = I1/3):
//   g = Al^2 + qh1^2 (m0 rho/(sqrt6 fc) + m_g/fc), Al = (1-qh1) Bl^2 + sqrt(3/2) rho/fc, Bl = sigV/fc + rho/(sqrt6 fc),
//   m_g = A_g B_g fc e^R, R = (sigV - qh2 ft/3)/(B_g fc), A_g = 3 ft qh2/fc + m0/2,
//   B_g = qh2/3 (1+ft/fc) / (ln A_g + ln(Df+1) - ln(2Df-1) - ln(3 qh2 + m0/2)).
// Returns dg/dsigV (gs), dg/drho (gr) and their (sigV, rho, kp) derivatives (the Hessian rows the analytic
// return-map Jacobian and consistent tangent need). Df clamped > 0.5; R capped at 700 (overflow guard for
// far-off Newton iterates only).
// ---------------------------------------------------------------------------
inline void cdpm2PotentialDerivs(double sigV, double rho, double kp, const Params& mp,
                                 double& gs, double& gr, double dgs[3], double dgr[3])
{
    const double fc = mp.fc, ft = mp.ft, m0 = mp.m0;
    const double Df = mp.Df > 0.5 + 1.0e-6 ? mp.Df : 0.5 + 1.0e-6;
    const double q1 = qh1Of(kp, mp.qh0, mp.Hp), q2 = qh2Of(kp, mp.Hp);
    const double dq1 = dqh1OfdKp(kp, mp.qh0, mp.Hp), dq2 = dqh2OfdKp(kp, mp.Hp);
    const double a = 1.0 - q1, c = 1.0 + ft / fc;
    const double AG = 3.0 * ft * q2 / fc + m0 / 2.0, AG_k = 3.0 * ft * dq2 / fc;
    const double L = std::log(AG) + std::log(Df + 1.0) - std::log(2.0 * Df - 1.0) - std::log(3.0 * q2 + m0 / 2.0);
    const double L_k = AG_k / AG - 3.0 * dq2 / (3.0 * q2 + m0 / 2.0);
    const double BG = q2 / 3.0 * c / L;
    const double BG_k = (dq2 / 3.0 * c * L - q2 / 3.0 * c * L_k) / (L * L);
    const double X = sigV - ft * q2 / 3.0;
    double Rg = X / (fc * BG); if (Rg > 700.0) Rg = 700.0;
    const double eR = std::exp(Rg), mQ = AG * eR;
    const double R_s = 1.0 / (fc * BG);
    const double R_k = -(ft * dq2 / 3.0) / (fc * BG) - X * BG_k / (fc * BG * BG);
    const double mQ_s = mQ * R_s, mQ_k = AG_k * eR + mQ * R_k;
    const double Bl = sigV / fc + rho / (SQRT6 * fc), Bl_s = 1.0 / fc, Bl_r = 1.0 / (SQRT6 * fc);
    const double Al = a * Bl * Bl + SQRT1_5 * rho / fc;
    const double Al_s = 2.0 * a * Bl * Bl_s, Al_r = 2.0 * a * Bl * Bl_r + SQRT1_5 / fc, Al_k = -dq1 * Bl * Bl;
    gs = 4.0 * a * Al * Bl / fc + q1 * q1 * mQ / fc;
    gr = Al / (SQRT6 * fc) * (4.0 * a * Bl + 6.0) + m0 * q1 * q1 / (SQRT6 * fc);
    dgs[0] = 4.0 * a * (Al_s * Bl + Al * Bl_s) / fc + q1 * q1 * mQ_s / fc;
    dgs[1] = 4.0 * a * (Al_r * Bl + Al * Bl_r) / fc;
    dgs[2] = (-4.0 * dq1 * Al * Bl + 4.0 * a * Al_k * Bl) / fc + (2.0 * q1 * dq1 * mQ + q1 * q1 * mQ_k) / fc;
    dgr[0] = (Al_s * (4.0 * a * Bl + 6.0) + Al * 4.0 * a * Bl_s) / (SQRT6 * fc);
    dgr[1] = (Al_r * (4.0 * a * Bl + 6.0) + Al * 4.0 * a * Bl_r) / (SQRT6 * fc);
    dgr[2] = (Al_k * (4.0 * a * Bl + 6.0) - Al * 4.0 * dq1 * Bl) / (SQRT6 * fc) + 2.0 * m0 * q1 * dq1 / (SQRT6 * fc);
}

// (m_v, m_s) = (dg/dxi, dg/drho) of the CDPM2 potential in the (xi, rho) frame + d/d(xi, rho, kp).
inline void cdpm2FlowGradJac(double xi, double rho, double kp, const Params& mp,
                             double& m_v, double& m_s, double dmv[3], double dms[3])
{
    double gs, gr, dgs[3], dgr[3];
    cdpm2PotentialDerivs(xi / SQRT3, rho, kp, mp, gs, gr, dgs, dgr);
    m_v = gs / SQRT3; m_s = gr;
    dmv[0] = dgs[0] / 3.0; dmv[1] = dgs[1] / SQRT3; dmv[2] = dgs[2] / SQRT3;
    dms[0] = dgr[0] / SQRT3; dms[1] = dgr[1]; dms[2] = dgr[2];
}

inline bool returnMapVertex(double sigV_tr, double rho_tr, const Params& mp, double kp_n,
                            double& sigV, double& kp, double& dlam)
{
    const double fc = mp.fc, K = bulkK(mp), G = shearG(mp);
    bool tension;
    if (sigV_tr > 0.0) tension = true;
    else if (sigV_tr < 0.0 && kp_n < 1.0) tension = false;
    else return false;
    auto F = [&](double s) { return yfInvHard(SQRT3 * s, 0.0, 1.0, vertexKappa(kp_n, sigV_tr, rho_tr, s, mp), mp); };
    if (!(F(sigV_tr) > 0.0 && F(0.0) < 0.0)) return false;   // (also rejects NaN)
    double a = tension ? 0.0 : sigV_tr, b = tension ? sigV_tr : 0.0;   // a < b
    double Fa = F(a);
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (a + b);
        if (mid <= a || mid >= b) break;
        const double Fm = F(mid);
        if (Fm == 0.0) { a = b = mid; break; }
        if ((Fm > 0.0) == (Fa > 0.0)) { a = mid; Fa = Fm; } else b = mid;
    }
    const double s = 0.5 * (a + b);
    const double k = vertexKappa(kp_n, sigV_tr, rho_tr, s, mp);
    double dgs, dgr;
    if (tension && mp.flowPotential == 1) {          // B1: the SAME potential as the regular map
        double d1[3], d2[3]; cdpm2PotentialDerivs(s, 0.0, k, mp, dgs, dgr, d1, d2);
    }
    else if (tension) { dgs = mp.Df * mp.m0 / fc; dgr = mp.m0 / (SQRT6 * fc); }
    else cdpm2VertexPotentialGrad(s, k, mp, dgs, dgr);
    if (dgs == 0.0) return false;
    const double dl = (sigV_tr - s) / (K * dgs);
    if (!(dl >= 0.0) || rho_tr > 2.0 * G * dl * dgr * (1.0 + 1.0e-10) + 1.0e-14 * fc) return false;
    sigV = s; kp = k; dlam = dl;
    return true;
}

// ---------------------------------------------------------------------------
// Hardening principal return (Grassl Eq.18 + Eq.30-36). 4-unknown (xi,rho,dlam,kp)
// Newton with ANALYTIC 4x4 Jacobian (the oracle uses a NUMERICAL Jacobian — this is
// the build-PR deliverable); theta frozen; hydrostatic apex fallback; HONEST
// convergence (independent f at the returned stress w/ its own Lode angle).
// ---------------------------------------------------------------------------
// Residual R[4] (and, if J != nullptr, the analytic 4x4 Jacobian) of the hardening return system with the FULL
// CDPM2 plastic potential (B1): R1 = xi - xi_tr + 3K dlam m_v, R2 = rho - rho_tr + 2G dlam m_s, R3 = f_p (frozen
// Lode r), R4 = kp - kp_n - dlam ||m||/xh cos2. Used by the globalized Newton (line search needs R alone).
//   legacyFlow = true (WP concrete3d-hang-diagnosis #877 follow-up): the SAME system with the legacy v1 flow
//   (m_v = Df m0/(sqrt3 fc) constant, m_s = 3 rho/fc^2 + m0/(sqrt6 fc)) so the globalized Newton can be the
//   rescue for flowPotential = legacy too; the Hessian terms are dmv = 0 and dms = (0, 3/fc^2, 0), which
//   reproduces the legacy-branch Jacobian rows of the plain newton() exactly.
inline void cdpm2HardeningResidual(double xi, double rho, double dlam, double kp, double xi_tr, double rho_tr,
                                   double kp_n, double r, double cos2, const Params& mp, double R[4],
                                   double J[4][4] = nullptr, bool legacyFlow = false)
{
    const double fc = mp.fc, m0 = mp.m0, K = bulkK(mp), G = shearG(mp);
    double m_v, m_s, dmv[3], dms[3];
    if (legacyFlow) {
        m_v = mp.Df * m0 / (SQRT3 * fc);
        m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc);
        dmv[0] = dmv[1] = dmv[2] = 0.0;
        dms[0] = 0.0; dms[1] = 3.0 / (fc * fc); dms[2] = 0.0;
    } else {
        cdpm2FlowGradJac(xi, rho, kp, mp, m_v, m_s, dmv, dms);
    }
    const double mnorm = std::sqrt(m_v * m_v + m_s * m_s);
    const double sigV = xi / SQRT3;
    const double xh = ductilityXh(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
    R[0] = xi - xi_tr + 3.0 * K * dlam * m_v;
    R[1] = rho - rho_tr + 2.0 * G * dlam * m_s;
    R[2] = yfInvHard(xi, rho, r, kp, mp);
    R[3] = kp - kp_n - dlam * mnorm / xh * cos2;
    if (!J) return;
    const double q1 = qh1Of(kp, mp.qh0, mp.Hp), q2 = qh2Of(kp, mp.Hp);
    const double dq1 = dqh1OfdKp(kp, mp.qh0, mp.Hp), dq2 = dqh2OfdKp(kp, mp.Hp);
    const double sigV_fc = xi / (SQRT3 * fc);
    const double AV = rho / (SQRT6 * fc) + sigV_fc;
    const double RR = rho * r / (SQRT6 * fc) + sigV_fc;
    const double cap = (1.0 - q1) * AV * AV + SQRT1_5 * rho / fc;
    const double dcap_dxi = (1.0 - q1) * 2.0 * AV / (SQRT3 * fc);
    const double dcap_drho = (1.0 - q1) * 2.0 * AV / (SQRT6 * fc) + SQRT1_5 / fc;
    const double dcap_dkp = -dq1 * AV * AV;
    const double dq1sq_q2 = 2.0 * q1 * q2 * dq1 + q1 * q1 * dq2;
    const double dq1sq_q2sq = 2.0 * q1 * q2 * q2 * dq1 + 2.0 * q1 * q1 * q2 * dq2;
    const double dxh_dxi = dDuctilityXhdSigV(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh) / SQRT3;
    double dmn[3];
    for (int j = 0; j < 3; ++j) dmn[j] = (m_v * dmv[j] + m_s * dms[j]) / mnorm;
    J[0][0] = 1.0 + 3.0 * K * dlam * dmv[0]; J[0][1] = 3.0 * K * dlam * dmv[1];
    J[0][2] = 3.0 * K * m_v;                 J[0][3] = 3.0 * K * dlam * dmv[2];
    J[1][0] = 2.0 * G * dlam * dms[0];       J[1][1] = 1.0 + 2.0 * G * dlam * dms[1];
    J[1][2] = 2.0 * G * m_s;                 J[1][3] = 2.0 * G * dlam * dms[2];
    J[2][0] = 2.0 * cap * dcap_dxi + m0 * q1 * q1 * q2 / (SQRT3 * fc);
    J[2][1] = 2.0 * cap * dcap_drho + m0 * q1 * q1 * q2 * r / (SQRT6 * fc);
    J[2][2] = 0.0;
    J[2][3] = 2.0 * cap * dcap_dkp + m0 * RR * dq1sq_q2 - dq1sq_q2sq;
    J[3][0] = -dlam * cos2 * (dmn[0] / xh - mnorm / (xh * xh) * dxh_dxi);
    J[3][1] = -dlam * cos2 * dmn[1] / xh;
    J[3][2] = -mnorm / xh * cos2;
    J[3][3] = 1.0 - dlam * cos2 * dmn[2] / xh;
}

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
    // The regular semi-implicit Newton. clampRho=false: the ORIGINAL scheme (an iterate crossing the
    // hydrostatic axis, rho<0, aborts with apex=true -> vertex candidate). clampRho=true: the RETRY used
    // only after a rejected vertex (the trial is OUTSIDE the cone of normals, so a regular rho>0 solution
    // exists and the abort was a Newton overshoot): rho is clamped to >=0 and iteration continues (OOFEM
    // performRegularReturn's max(rho,0)).
    const bool cdpm2Flow = (mp.flowPotential == 1);
    auto newton = [&](bool clampRho) {
    xi = xi_tr; rho = rho_tr; dlam = 0.0; kp = kp_n; apex = false; converged = false;
    for (int it = 0; it < 100; ++it) {
        // plastic potential gradient (B1): the full CDPM2 potential, or the legacy v1 flow (byte-identical)
        double m_vv = m_v, m_s, dmv[3] = {0.0, 0.0, 0.0}, dms[3] = {0.0, 0.0, 0.0};
        if (cdpm2Flow) cdpm2FlowGradJac(xi, rho, kp, mp, m_vv, m_s, dmv, dms);
        else m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc);
        const double mnorm = std::sqrt(m_vv * m_vv + m_s * m_s);
        const double sigV = xi / SQRT3;
        const double xh = ductilityXh(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
        const double R1 = xi - xi_tr + 3.0 * K * dlam * m_vv;
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
        double J[4][4] = {
            { 1.0, 0.0, 3.0 * K * m_v, 0.0 },
            { 0.0, 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), 2.0 * G * m_s, 0.0 },
            { dR3_dxi, dR3_drho, 0.0, dR3_dkp },
            { -dlam * cos2 * dg_dxi, -dlam * cos2 * dg_drho, -g_val * cos2, 1.0 }
        };
        if (cdpm2Flow) {   // B1: rows 1, 2, 4 with the potential Hessian (row 3, the yield function, unchanged)
            double dmn[3];
            for (int j = 0; j < 3; ++j) dmn[j] = (m_vv * dmv[j] + m_s * dms[j]) / mnorm;
            J[0][0] = 1.0 + 3.0 * K * dlam * dmv[0]; J[0][1] = 3.0 * K * dlam * dmv[1];
            J[0][2] = 3.0 * K * m_vv;                J[0][3] = 3.0 * K * dlam * dmv[2];
            J[1][0] = 2.0 * G * dlam * dms[0];       J[1][1] = 1.0 + 2.0 * G * dlam * dms[1];
            J[1][2] = 2.0 * G * m_s;                 J[1][3] = 2.0 * G * dlam * dms[2];
            J[3][0] = -dlam * cos2 * (dmn[0] / xh - mnorm / (xh * xh) * dxh_dxi);
            J[3][1] = -dlam * cos2 * dmn[1] / xh;
            J[3][2] = -g_val * cos2;
            J[3][3] = 1.0 - dlam * cos2 * dmn[2] / xh;
        }
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
        if (rho < 0.0) { if (!clampRho) { apex = true; break; } rho = 0.0; }
    }
    };
    // B1 GLOBALIZED Newton for the CDPM2 potential (mirror of the oracle _newton_glob): OOFEM's projections
    // (rho >= 0, dlam >= 0, kp >= kp_n) + a backtracking line search on the scaled residual. The plain Newton fails
    // sporadically on far trials (iterates cross the kp = 1 kink / go kp < 0) and the sub-incremented fallback
    // then made the stress a DISCONTINUOUS function of the strain (spurious element-Newton roots: OOFEM con2dpm2
    // at one sub-step gave -3.22 MPa + return-map warnings in the C++ build). A vertex solution cannot satisfy R2
    // with rho pinned at 0 => stop early (rho stuck at 0 for 3 iterations) and flag apex => vertex candidate.
    //
    // WP concrete3d-hang-diagnosis review #877, defect 2 (MAJOR), option A (mirror of the oracle _newton_glob):
    // the ORIGINAL scheme projected (rho, dlam, kp) onto their admissible ranges on EVERY line-search trial
    // iterate, not just the accepted one. In a TENSION-dominated trial at kappa_p < 1 (m0*RR > 1, the hardening
    // system is locally INDEFINITE there -- df/dkappa_p > 0) that per-iterate projection repeatedly pins the
    // iterate back onto the same clamped face: the line search bottoms out at a=1/64 nearly every step and the
    // loop burns its 100-iteration budget (~800 residual evaluations) before falling through to the plain Newton
    // anyway, which is also what made the direct and sub-incremented returns land on different states near
    // first cracking (stress discontinuous in strain). Fix: the line search evaluates the UNPROJECTED iterate;
    // the physically-required projection (rho>=0, dlam>=0, kp>=kp_n) is applied exactly once, to the FINAL
    // returned iterate (on convergence and on the rho-stuck apex exit). The caller's admissibility gate
    // (dlam>=-1e-12, kp>=kp_n-1e-12, on-surface f_after) remains the honesty check on whatever root is found.
    auto newtonGlob = [&]() {
        xi = xi_tr; rho = rho_tr; dlam = 0.0; kp = kp_n; apex = false; converged = false;
        double Rr[4];
        cdpm2HardeningResidual(xi, rho, dlam, kp, xi_tr, rho_tr, kp_n, r, cos2, mp, Rr, nullptr, !cdpm2Flow);
        int stuck = 0;
        auto project = [&]() {
            if (rho < 0.0) rho = 0.0;
            if (dlam < 0.0) dlam = 0.0;
            if (kp < kp_n) kp = kp_n;
        };
        for (int it = 0; it < 100; ++it) {
            if (std::fabs(Rr[0]) < tol * fc && std::fabs(Rr[1]) < tol * fc
                && std::fabs(Rr[2]) < tol && std::fabs(Rr[3]) < tol) {
                // ADMISSIBILITY on the UNPROJECTED root (review #877 minor 1; mirror of the oracle _newton_glob): the caller's
                // gate sees the projected dlam / kp (admissible by construction), so a root with dlam < 0 or kp < kp_n used to
                // be accepted as its clamp. Outside the cone it is a non-convergence and falls through to the plain scheme /
                // vertex return.
                const bool admissibleRoot = (dlam >= -1.0e-12) && (kp >= kp_n - 1.0e-12);
                project(); converged = admissibleRoot; return;
            }
            double Rj[4], J[4][4];
            cdpm2HardeningResidual(xi, rho, dlam, kp, xi_tr, rho_tr, kp_n, r, cos2, mp, Rj, J, !cdpm2Flow);
            double M[4][5];
            for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = J[i][j]; M[i][4] = -Rr[i]; }
            for (int c = 0; c < 4; ++c) {
                int piv = c; for (int rr = c + 1; rr < 4; ++rr) if (std::fabs(M[rr][c]) > std::fabs(M[piv][c])) piv = rr;
                for (int j = 0; j < 5; ++j) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
                if (M[c][c] == 0.0) return;
                for (int rr = 0; rr < 4; ++rr) if (rr != c) { double f = M[rr][c] / M[c][c]; for (int j = c; j < 5; ++j) M[rr][j] -= f * M[c][j]; }
            }
            double step[4];
            for (int i = 0; i < 4; ++i) { step[i] = M[i][4] / M[i][i]; if (!std::isfinite(step[i])) return; }
            const double sc[4] = {1.0 / fc, 1.0 / fc, 1.0, 1.0};
            double n0 = 0.0; for (int i = 0; i < 4; ++i) n0 += (Rr[i] * sc[i]) * (Rr[i] * sc[i]);
            n0 = std::sqrt(n0);
            double a = 1.0, un[4], Rn[4];
            for (;;) {
                un[0] = xi + a * step[0]; un[1] = rho + a * step[1]; un[2] = dlam + a * step[2]; un[3] = kp + a * step[3];
                // UNPROJECTED (option A): no per-iterate clamp; see the block comment above newtonGlob.
                cdpm2HardeningResidual(un[0], un[1], un[2], un[3], xi_tr, rho_tr, kp_n, r, cos2, mp, Rn, nullptr, !cdpm2Flow);
                double nn = 0.0; bool fin = true;
                for (int i = 0; i < 4; ++i) { fin = fin && std::isfinite(Rn[i]); nn += (Rn[i] * sc[i]) * (Rn[i] * sc[i]); }
                if ((fin && std::sqrt(nn) < (1.0 - 1.0e-4 * a) * n0) || a < 1.0 / 64.0) break;
                a *= 0.5;
            }
            xi = un[0]; rho = un[1]; dlam = un[2]; kp = un[3];
            for (int i = 0; i < 4; ++i) Rr[i] = Rn[i];
            stuck = (rho <= 0.0) ? stuck + 1 : 0;
            if (stuck >= 3) { project(); apex = true; return; }
        }
    };
    // cdpm2: the globalized Newton first; if it fails (not an axis overshoot), fall back to the plain scheme
    // (whose rho<0 abort feeds the vertex test + the clamped retry below) before giving up.
    // legacy (WP concrete3d-hang-diagnosis #877 follow-up): the direct plain Newton runs first EXACTLY as
    // before (every converging case, hence every pinned fixture, is byte-identical); only a non-convergent,
    // non-apex outcome hands over to the globalized Newton (legacy-flow residual, unprojected line search) as
    // a rescue, and only then to the honest failure. Measured on the plain legacy Newton alone: 10-13 % of
    // ordinary tension-dominated first-cracking increments (expansive lateral strain, what an FE Newton
    // iterate produces) failed the return map and silently fell to the elastic trial.
    if (cdpm2Flow && xi_tr > 0.0) {
        // TENSION-dominated trial (sigma_V_trial > 0): plain Newton FIRST (mirror of the oracle). Measured at the
        // first-crack step (virgin sigma_xx = 2 MPa, kappa_p < 1): newtonGlob burns 660-900 residual evaluations
        // failing (line search bottoms out where the hardening system is locally indefinite) and the plain Newton
        // then converges in ~25 iterations to the SAME state; newtonGlob (option A) is the rescue. If the plain
        // scheme aborted on an axis overshoot (apex) and the rescue does not converge either, keep the apex
        // verdict so the vertex test below still runs.
        newton(false);
        if (!converged) { const bool plainApex = apex; newtonGlob(); if (!converged && !apex) apex = plainApex; }
    }
    else if (cdpm2Flow) { newtonGlob(); if (!converged && !apex) newton(false); }
    else { newton(false); if (!converged && !apex) newtonGlob(); }
    const bool overshot = apex;

    R.xi = xi; R.rho = rho; R.dlam = dlam; R.kp = kp; R.plastic = true;
    if (!apex) {
        const double p_new = xi / SQRT3;
        const double dev_scale = (rho_tr > 0.0) ? rho / rho_tr : 0.0;
        for (int a = 0; a < 3; ++a) R.sp[a] = R.s_tr[a] * dev_scale + p_new;
        // HONEST convergence: recompute f at the returned stress with ITS OWN Lode angle.
        double svn[6] = {R.sp[0], R.sp[1], R.sp[2], 0, 0, 0};
        R.f_after = yieldF(svn, mp, qh1Of(kp, mp.qh0, mp.Hp), qh2Of(kp, mp.Hp));
        // ADMISSIBILITY (PR #249 adversarial-review fix): a valid plastic return needs dlam>=0 and a
        // NON-DECREASING hardening variable (kp>=kp_n). Never report converged for an inadmissible/
        // off-surface state.
        const bool admissible = std::isfinite(R.f_after) && dlam >= -1.0e-12 && kp >= kp_n - 1.0e-12;
        if (converged && std::fabs(R.f_after) < F_TOL_HONEST && admissible) {
            R.converged = true;
            return R;
        }
    }
    // VERTEX RETURN (WP concrete3d-oracle-diagnosis; mirrors the oracle return_map_vertex). Reached when
    // the regular (radial) return overshot the hydrostatic axis (rho<0), did not converge, or landed
    // inadmissible. The OLD branch projected EVERY such trial onto the hydrostatic-TENSION vertex with the
    // kp of the ABORTED Newton iterate: a deep-compression trial was sign-flipped to tension (then rejected
    // by the PR #249 gate => ELASTIC fallback — hydrostatic compression never yielded, OOFEM con2dpm3), and
    // the tension-vertex kp was an arbitrary iterate value (hydrostatic tension step-size dependent, no
    // convergence under refinement, con2dpm4). returnMapVertex solves f(sigV, rho=0; kp(sigV)) = 0 with kp
    // CONSISTENT with the vertex plastic strain, on the TENSION vertex (sigV_tr>0) or — while the [1-qh1]
    // cap closes the surface (kp_n<1) — the COMPRESSION vertex, accepted only inside the cone of plastic-
    // potential normals. Otherwise: SAFE honest failure (elastic predictor, converged=false => the
    // caller's status!=0 cuts the step) — the PR #249 contract is unchanged.
    {
        double sV = 0.0, kpv = kp_n, dlv = 0.0;
        if (returnMapVertex(xi_tr / SQRT3, rho_tr, mp, kp_n, sV, kpv, dlv)) {
            double svv[6] = {sV, sV, sV, 0, 0, 0};
            const double fv = yieldF(svv, mp, qh1Of(kpv, mp.qh0, mp.Hp), qh2Of(kpv, mp.Hp));
            if (std::isfinite(fv) && std::fabs(fv) < F_TOL_HONEST && kpv >= kp_n - 1.0e-12) {
                for (int a = 0; a < 3; ++a) R.sp[a] = sV;
                R.xi = SQRT3 * sV; R.rho = 0.0; R.dlam = dlv; R.kp = kpv;
                R.apex = true; R.converged = true; R.f_after = fv;
                return R;
            }
        }
    }
    // Regular RETRY with rho clamped at 0 (only after an axis overshoot whose vertex was rejected).
    if (overshot) {
        newton(true);
        if (converged && rho > 0.0) {
            const double p_new = xi / SQRT3;
            const double dev_scale = (rho_tr > 0.0) ? rho / rho_tr : 0.0;
            for (int a = 0; a < 3; ++a) R.sp[a] = R.s_tr[a] * dev_scale + p_new;
            double svn[6] = {R.sp[0], R.sp[1], R.sp[2], 0, 0, 0};
            R.f_after = yieldF(svn, mp, qh1Of(kp, mp.qh0, mp.Hp), qh2Of(kp, mp.Hp));
            const bool admissible = std::isfinite(R.f_after) && dlam >= -1.0e-12 && kp >= kp_n - 1.0e-12;
            if (std::fabs(R.f_after) < F_TOL_HONEST && admissible) {
                R.xi = xi; R.rho = rho; R.dlam = dlam; R.kp = kp; R.apex = false; R.converged = true;
                return R;
            }
        }
    }
    for (int a = 0; a < 3; ++a) R.sp[a] = sigTr[a];   // safe fallback = elastic predictor
    R.kp = kp_n; R.xi = xi_tr; R.rho = rho_tr; R.dlam = 0.0; R.apex = false;
    R.f_after = yfInvHard(xi_tr, rho_tr, r, kp_n, mp);
    R.converged = false;
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
// Hardening-history sensitivities of ONE direct return (needed to chain sub-incremented pieces exactly): with
// sigma_new = R(sigma_tr, kappa_n) and kappa_new = K(sigma_tr, kappa_n),
//   Sk = d sigma_new / d kappa_n   (tensor-Voigt vector),
//   R  = d kappa_new / d sigma_tr  (row over the tensor-Voigt sigma_tr; shear entries carry the factor 2),
//   Kk = d kappa_new / d kappa_n.
// Defaults are the ELASTIC piece (Sk = 0, R = 0, Kk = 1). Fed by the same implicit-function solve as the principal
// Jacobian (one extra right-hand side, dR4/dkappa_n = -1), so it costs one more back-substitution.
// ---------------------------------------------------------------------------
struct PieceSens {
    double Sk[6] = {0, 0, 0, 0, 0, 0};
    double R[6]  = {0, 0, 0, 0, 0, 0};
    double Kk = 1.0;
};

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
                              bool hardening, double D[3][3], double* hs = nullptr)
{
    // hs (optional, 7 doubles, hardening branch): [0..2] d sp_a/d kappa_n, [3..5] d kappa_new/d w_b, [6] d kappa_new/d kappa_n
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
        double Ju[4][4] = {
            { 1.0, 0.0, 3.0 * K * m_v, 0.0 },
            { 0.0, 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), 2.0 * G * m_s, 0.0 },
            { dR3_dxi, dR3_drho, 0.0, dR3_dkp },
            { -dlam * cos2 * dg_dxi, -dlam * cos2 * dg_drho, -g_val * cos2, 1.0 }
        };
        double gvAct = g_val;          // ||m||/xh of the ACTIVE potential (the R4 trial-cos2 derivative below)
        if (mp.flowPotential == 1) {   // B1: the CDPM2-potential rows (mirror of the returnMapHardening Jacobian)
            double mv, ms, dmv[3], dms[3], dmn[3];
            cdpm2FlowGradJac(xi, rho, kp, mp, mv, ms, dmv, dms);
            const double mn = std::sqrt(mv * mv + ms * ms), gv = mn / xh;
            for (int j = 0; j < 3; ++j) dmn[j] = (mv * dmv[j] + ms * dms[j]) / mn;
            Ju[0][0] = 1.0 + 3.0 * K * dlam * dmv[0]; Ju[0][1] = 3.0 * K * dlam * dmv[1];
            Ju[0][2] = 3.0 * K * mv;                  Ju[0][3] = 3.0 * K * dlam * dmv[2];
            Ju[1][0] = 2.0 * G * dlam * dms[0];       Ju[1][1] = 1.0 + 2.0 * G * dlam * dms[1];
            Ju[1][2] = 2.0 * G * ms;                  Ju[1][3] = 2.0 * G * dlam * dms[2];
            Ju[3][0] = -dlam * cos2 * (dmn[0] / xh - mn / (xh * xh) * dxh_dxi);
            Ju[3][1] = -dlam * cos2 * dmn[1] / xh;
            Ju[3][2] = -gv * cos2;
            Ju[3][3] = 1.0 - dlam * cos2 * dmn[2] / xh;
            gvAct = gv;
        }
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
            const double rhs3 = -dlam * gvAct * dcos2_dw[b];            // dR4/dcos2 * dcos2/dw_b
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
            if (hs) hs[3 + b] = M[3][4] / M[3][3];                     // d kappa_new / d w_b
            const double dscale = (rhotr > 0.0) ? (drho * rhotr - rho * (s_tr_b / rhotr)) / (rhotr * rhotr) : 0.0;
            const double scale = (rhotr > 0.0) ? rho / rhotr : 0.0;
            for (int a = 0; a < 3; ++a) {
                const double ds_tr_a = (a == b ? 1.0 : 0.0) - 1.0 / 3.0;
                D[a][b] = ds_tr_a * scale + pr.s_tr[a] * dscale + (1.0 / SQRT3) * dxi;
            }
        }
        if (hs) {
            // kappa_n column: R4 = kappa - kappa_n - ..., so J u = +e4 (u = d(xi,rho,dlam,kappa)/d kappa_n)
            double M[4][5];
            for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = Ju[i][j]; M[i][4] = (i == 3) ? 1.0 : 0.0; }
            for (int c = 0; c < 4; ++c) {
                int piv = c; for (int rr = c + 1; rr < 4; ++rr) if (std::fabs(M[rr][c]) > std::fabs(M[piv][c])) piv = rr;
                for (int j = 0; j < 5; ++j) { double tt = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = tt; }
                for (int rr = 0; rr < 4; ++rr) if (rr != c) { double f = M[rr][c] / M[c][c]; for (int j = c; j < 5; ++j) M[rr][j] -= f * M[c][j]; }
            }
            const double dxi_k = M[0][4] / M[0][0], drho_k = M[1][4] / M[1][1];
            for (int a = 0; a < 3; ++a) hs[a] = pr.s_tr[a] * ((rhotr > 0.0) ? drho_k / rhotr : 0.0) + dxi_k / SQRT3;
            hs[6] = M[3][4] / M[3][3];
        }
    }
}

// ---------------------------------------------------------------------------
// Principal Jacobian of the hardening VERTEX return (returnMapVertex). sp_a = sigV for all a, where
// F(sigV; sigV_tr, rho_tr) = f(sigV, 0; kp(sigV, sigV_tr, rho_tr)) = 0 (implicit-function theorem):
//   dsigV/dw_b = -(F_sigVtr * 1/3 + F_rhotr * s_tr_b/rho_tr) / F_sigV ,
//   F_x = f_x + f_kp kp_x ;  kp = kp_n + eq/xh(sigV), eq = sqrt((sigV_tr-sigV)^2/(3K^2) + (rho_tr/2G)^2).
// The rho_tr/rho_tr factor cancels analytically (kp_rhotr * s_b/rho_tr = s_b/(4G^2 eq xh)), so an exactly
// hydrostatic trial is regular. Every row equal => zero deviatoric stiffness (the stress is pinned to the
// axis) — the rank-deficient apex tangent handoff §6 listed as owed. Returns false on a degenerate state.
// ---------------------------------------------------------------------------
inline bool vertexPrincipalJacobian(const double w[3], const PrincipalResult& pr, const Params& mp, double D[3][3],
                                    double* hs = nullptr)
{
    const double fc = mp.fc, m0 = mp.m0, K = bulkK(mp), G = shearG(mp);
    const double sigV = pr.xi / SQRT3, kp = pr.kp;
    const double sigV_tr = (w[0] + w[1] + w[2]) / 3.0;
    const double d = sigV_tr - sigV, g = pr.rho_tr / (2.0 * G);
    const double eq = std::sqrt(d * d / (3.0 * K * K) + g * g);
    const double xh = ductilityXh(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
    const double dxh = dDuctilityXhdSigV(sigV, fc, mp.Ah, mp.Bh, mp.Ch, mp.Dh);
    if (!(eq > 0.0) || !(xh > 0.0)) return false;
    const double kp_s  = (-d / (3.0 * K * K)) / (eq * xh) - eq * dxh / (xh * xh);   // dkp/dsigV
    const double kp_st = ( d / (3.0 * K * K)) / (eq * xh);                           // dkp/dsigV_tr
    // f on the axis: f = cap^2 + m0 q1^2 q2 s - q1^2 q2^2, cap = (1-q1) s^2, s = sigV/fc
    const double q1 = qh1Of(kp, mp.qh0, mp.Hp), q2 = qh2Of(kp, mp.Hp);
    const double dq1 = dqh1OfdKp(kp, mp.qh0, mp.Hp), dq2 = dqh2OfdKp(kp, mp.Hp);
    const double s = sigV / fc, cap = (1.0 - q1) * s * s;
    const double f_s = (2.0 * cap * (1.0 - q1) * 2.0 * s + m0 * q1 * q1 * q2) / fc;
    const double f_k = 2.0 * cap * (-dq1 * s * s) + m0 * s * (2.0 * q1 * q2 * dq1 + q1 * q1 * dq2)
                     - (2.0 * q1 * q2 * q2 * dq1 + 2.0 * q1 * q1 * q2 * dq2);
    const double F_s = f_s + f_k * kp_s;
    if (!(std::fabs(F_s) > 0.0)) return false;
    for (int b = 0; b < 3; ++b) {
        const double dev_b = f_k * pr.s_tr[b] / (4.0 * G * G * eq * xh);         // F_rhotr * s_b/rho_tr
        const double dsdw = -(f_k * kp_st / 3.0 + dev_b) / F_s;
        for (int a = 0; a < 3; ++a) D[a][b] = dsdw;
        // d kappa_new / d w_b = kp_s dsigV/dw_b + kp_st/3 + s_b/(4 G^2 eq xh)
        if (hs) hs[3 + b] = kp_s * dsdw + kp_st / 3.0 + pr.s_tr[b] / (4.0 * G * G * eq * xh);
    }
    if (hs) {                                   // kappa_n enters kp = kp_n + eq/xh with unit weight: F_kn = f_k
        const double dsdk = -f_k / F_s;
        for (int a = 0; a < 3; ++a) hs[a] = dsdk;
        hs[6] = 1.0 + kp_s * dsdk;
    }
    return true;
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
                              double Dtan6[6][6], PieceSens* ps = nullptr)
{
    if (ps) *ps = PieceSens();   // elastic / fallback default
    // Elastic step, the non-converged safe fallback (returnMapHardening reset to the elastic
    // predictor), or a converged apex => the elastic operator. NOTE (PR #249 review): at a true
    // apex the physical tangent collapses toward zero (the stress is pinned at the vertex — an
    // ATTRACTING PLATEAU, NOT measure-zero), so this elastic surrogate is ~K too STIFF there and
    // can slow the global Newton. That is acceptable for P1 because the apex is the deferred
    // KNOWN-GAP regime (handoff §6) and the return map flags it (converged=false ⇒ step-cut);
    // the dedicated rank-deficient apex tangent lands with the apex sub-algorithm (P2+).
    // (WP concrete3d-oracle-diagnosis) the HARDENING vertex return now has its own analytic principal
    // Jacobian (vertexPrincipalJacobian: rank-1, purely volumetric, zero deviatoric stiffness); the
    // perfect-plastic apex (returnMapPrincipal) keeps the elastic surrogate.
    if (!pr.plastic || !pr.converged || (pr.apex && !hardening)) { elasticC(mp, Dtan6); return; }
    // eigenprojections E_a = e_a (x) e_a (rank-1)
    double E[3][3][3];
    for (int a = 0; a < 3; ++a) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) E[a][i][j] = V[i][a] * V[j][a];
    double D[3][3];
    double hs[7] = {0, 0, 0, 0, 0, 0, 1.0};
    double* hsp = (ps && hardening) ? hs : nullptr;
    if (pr.apex) {
        if (!vertexPrincipalJacobian(w, pr, mp, D, hsp)) { elasticC(mp, Dtan6); return; }
    } else {
        principalJacobian(w, pr, mp, hardening, D, hsp);
    }
    if (hsp) {
        double Sm[3][3], Gm[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
            double s1 = 0.0, s2 = 0.0;
            for (int a = 0; a < 3; ++a) { s1 += hs[a] * E[a][i][j]; s2 += hs[3 + a] * E[a][i][j]; }
            Sm[i][j] = s1; Gm[i][j] = s2;
        }
        matToVoigt(Sm, ps->Sk);
        ps->R[0] = Gm[0][0]; ps->R[1] = Gm[1][1]; ps->R[2] = Gm[2][2];
        ps->R[3] = 2.0 * Gm[0][1]; ps->R[4] = 2.0 * Gm[1][2]; ps->R[5] = 2.0 * Gm[0][2];
        ps->Kk = hs[6];
    }

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
    // P2g — MONOTONE (no-heal) cyclic damage: the running MAX of each channel's effective drive stress.
    // omega is solved against these (not the live drive), so on an elastic unload (live drive drops,
    // histories frozen) omega stays FIXED (no healing) and the nominal stress unloads along the damage
    // secant (1-omega)*sig_bar. On any monotonic path max == live => byte-identical to the pre-P2g kernel.
    double sigtMax = 0.0, sigcMax = 0.0;
    // B2 CDPM2 compressive drive: the compressive equivalent strain eqc = sum alpha_c d(eps_tilde) (Eq.47) and the
    // previous step's eps_tilde (unused by the legacy drive; zero-initialized => legacy byte-identical).
    double eqc = 0.0, etPrev = 0.0;
    // diagnostic of the LAST return (not history): 0 direct, n >= 2 sub-incremented pieces, -1 final failure
    int subInfo = 0;
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

// forward decl (damagedTangent's d beta_c/dε does a composite micro-FD through the return map, and
// returnMapTensor is defined further down)
inline int returnMapTensor(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                           bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                           bool doTangent, int* subInfo = nullptr);

// CDPM2 beta_c (Eq.50, P2f): the factor scaling the PLASTIC-strain part of the compressive-damage
// driver kappa_dc1 (Eq.48). beta_c = ft*qh2(kp)*sqrt(2/3) / (rho_bar*sqrt(1+2*Df^2)), rho_bar = sqrt(2 J2)
// of the EFFECTIVE stress. Mirror of the oracle beta_c(). In monotonic compression ~ft/(fc*sqrt(1+2Df^2))
// << 1, so it makes compression markedly MORE DUCTILE than the beta_c=1 simplification (faithful CDPM2).
// rho_bar->0 guard (hydrostatic) + clamp [0,1] (a plastic-contribution fraction <=1; inactive in the
// damaging regime where rho_bar is large, so the clamp never binds and the analytic tangent stays smooth).
inline double betaC(const double sig_pr[3], double kp, const Params& mp)
{
    double sv[6] = { sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0 };
    double xi, rho, theta; invariants(sv, xi, rho, theta);
    if (rho <= 1.0e-12) return 1.0;
    double bc = mp.ft * qh2Of(kp, mp.Hp) * std::sqrt(2.0 / 3.0) / (rho * std::sqrt(1.0 + 2.0 * mp.Df * mp.Df));
    if (bc < 0.0) bc = 0.0; if (bc > 1.0) bc = 1.0;
    return bc;
}

// ---------------------------------------------------------------------------
// P2 dual-damage NOMINAL stress update (mirrors the oracle damaged_step_tensor EXACTLY): from the
// committed history `in` + the new EFFECTIVE stress/strain, accumulate the CDPM2 damage drivers,
// solve omega_t/omega_c (bracketed, physical-floored), and recompose the nominal stress via the
// spectral split  sigma = (1-omega_t)<sig_bar>+ + (1-omega_c)<sig_bar>-  (Eq.1). Writes the new
// damage history into `out`. The split is recomputed from the CONVERGED effective stress every step
// => automatic, tier-independent unilateral crack-closure (ADR §4.3 BLOCKING). `kp` = the current
// (post-return) hardening variable, for beta_c (Eq.50).
// ---------------------------------------------------------------------------
inline void damagedUpdate(const Params& mp, const State& in, const double sig_eff[6], double kp,
                          const double eps_new[6], State& out,
                          double* wtOut = nullptr, double* wcOut = nullptr)
{
    const double eps0 = mp.ft / mp.E;
    const double eps_fc = epsFcOf(mp);

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
    double dnc = dnorm;
    {
        double depl[6]; for (int i = 0; i < 6; ++i) depl[i] = epl[i] - epl_n[i];
        const double wtw = tensileDamageWeight(mp, ac, depl, w, V);     // P2h ctTemper weight (1 if none)
        tensionHistUpdate(mp, kdt1, kdt2, et, et_max_n, dnorm, xs, wtw); // Eq.44/45 (law-dependent)
        dnc = dnorm * compressiveDamageWeight(mp, depl, w, V);          // PV20 tcTemper (== dnorm for none)
    }
    const double etc = compressiveEquivStrain(mp, w, et);                // PV20 tcTemper (== et for none)
    double eqc = in.eqc;
    if (mp.compDrive == 1) {                                             // B2: CDPM2 Eq.47-49
        compHistUpdateCdpm2(mp, kdc, kdc1, kdc2, eqc, in.etPrev, etc, ac, betaC(w, kp, mp), dnc, xs);
    } else if (loading) {
        kdc  += ac * above;          kdc2 += ac * above / xs;            // Eq.47 / Eq.49
        kdc1 += ac * betaC(w, kp, mp) * dnc / xs;                        // Eq.48 with the full CDPM2 beta_c (Eq.50, P2f)
    }
    const double et_max = et_max_n > et ? et_max_n : et;

    // P2i — multiaxial-consistent TENSILE drive E*et (Eq.37 equivalent strain) instead of the extreme
    // tensile principal, gated by the presence of a real tensile principal (reduces to the extreme
    // principal in uniaxial tension, E*et == sig_bar_t; in biaxial/triaxial tension E*et > the extreme
    // principal => damage onsets at a lower per-principal stress, the CDPM2-consistent envelope). The
    // COMPRESSIVE drive stays the extreme principal: et is ft-scaled (== eps0 on ANY failure surface), so
    // E*et could never reach fc and would never onset wc. Physical FLOOR (review-fix): never solve omega
    // on a numerical-residual stress (~1e-10 MPa) -> flips wt 0<->1 in compression.
    double maxw = w[0]; for (int i = 1; i < 3; ++i) if (w[i] > maxw) maxw = w[i];
    double Dt = (maxw > 1.0e-6 * mp.ft) ? mp.E * et : 0.0;
    double mn = w[0]; for (int i = 1; i < 3; ++i) if (w[i] < mn) mn = w[i];
    double Dc = -mn; if (Dc < 0.0) Dc = 0.0;
    // P2g — drive omega with the MONOTONE running max (no heal on unload); max == live on monotonic paths.
    const double sigtMax = in.sigtMax > Dt ? in.sigtMax : Dt;
    const double sigcMax = in.sigcMax > Dc ? in.sigcMax : Dc;
    const double wt = (et_max > eps0 && sigtMax > 1.0e-6 * mp.ft)
                    ? omegaT(mp, kdt1, kdt2, sigtMax) : 0.0;
    const double wc = omegaC(mp, kdc, kdc1, kdc2, sigcMax, eps_fc);

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
    out.sigtMax = sigtMax; out.sigcMax = sigcMax;   // P2g monotone drive history
    out.eqc = eqc; out.etPrev = etc;                // B2 CDPM2 compressive drive history (PV20: eps_tilde_c)
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
    if (which == 3) return compressiveEquivStrain(mp, w, equivStrainGeneral(w, mp));   // PV20 tcTemper
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
// Return-map evaluation at the committed state for the composite micro-FDs of the damaged tangent. For a live point this is
// exactly returnMapTensor(mp, in.sigEff, d, in.kp, true, ...). For a TENSION-DEAD point (returnMap re-bases `in` on the
// trial: eps = new strain, sigEff = sig_tr, so the FD increment d is a perturbation of the trial) the perturbed trial is
// split spectrally, the tensile part is carried elastically and the return map runs on the compressive remainder -- the
// same map returnMap applies, so the FD tracks the tangent actually being assembled.
//
// PIECE COUNT PINNED: the deterministic map is discontinuous where n = ceil(f_tr/c) changes, so a +/- leg that lands on
// the other side of an n boundary would differentiate the jump. Both legs therefore use the CENTRAL evaluation's n
// (rmPiecesFD at the central increment, passed as nForce); if the central point itself sits exactly on a boundary the
// legs simply follow the central side, which is the branch returnMap took for the reported stress.
inline int detPieces(const Params& mp, const double sig_n[6], const double deps[6], double kp_n);

inline int rmPiecesFD(const Params& mp, const State& in, const double d[6])
{
    if (!(in.wt >= mp.omegaDead)) return detPieces(mp, in.sigEff, d, in.kp);
    double sigTr[6]; elasticPredTensor(in.sigEff, d, mp, sigTr);
    double A[3][3], w[3], V[3][3]; voigtToMat(sigTr, A); eig3sym(A, w, V);
    double sm[3]; for (int a = 0; a < 3; ++a) sm[a] = w[a] < 0.0 ? w[a] : 0.0;
    double S[3][3], minus[6];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = 0.0; for (int a = 0; a < 3; ++a) v += V[i][a] * sm[a] * V[j][a];
        S[i][j] = v;
    }
    matToVoigt(S, minus);
    const double zero[6] = {0, 0, 0, 0, 0, 0};
    return detPieces(mp, minus, zero, in.kp);
}

inline void rmForFD(const Params& mp, const State& in, const double d[6], int nForce, double sb[6], double& kp)
{
    double dum[6][6];
    Params q = mp; q.subIncrForceN = nForce;
    if (!(in.wt >= mp.omegaDead)) { returnMapTensor(q, in.sigEff, d, in.kp, true, sb, kp, dum, false); return; }
    double sigTr[6]; elasticPredTensor(in.sigEff, d, mp, sigTr);
    double A[3][3], w[3], V[3][3]; voigtToMat(sigTr, A); eig3sym(A, w, V);
    double sm[3], sq[3];
    for (int a = 0; a < 3; ++a) { sm[a] = w[a] < 0.0 ? w[a] : 0.0; sq[a] = w[a] > 0.0 ? w[a] : 0.0; }
    double S[3][3], Q[3][3], minus[6], plus[6];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = 0.0, u = 0.0;
        for (int a = 0; a < 3; ++a) { v += V[i][a] * sm[a] * V[j][a]; u += V[i][a] * sq[a] * V[j][a]; }
        S[i][j] = v; Q[i][j] = u;
    }
    matToVoigt(S, minus); matToVoigt(Q, plus);
    const double zero[6] = {0, 0, 0, 0, 0, 0};
    returnMapTensor(q, minus, zero, in.kp, true, sb, kp, dum, false);
    for (int i = 0; i < 6; ++i) sb[i] += plus[i];
}

inline void damagedTangent(const Params& mp, const State& in, const double sig_eff[6],
                           const double eps_new[6], double kp_new, const double Ceff[6][6], double D6[6][6])
{
    static const double W6[6] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    const double eps0   = mp.ft / mp.E;
    const double eps_f  = mp.Gf / (mp.ft * mp.lch);
    const double eps_fc = epsFcOf(mp);
    // piece count of the CENTRAL evaluation, pinned on every micro-FD leg below (see rmForFD)
    int nFD = 1;
    { double dC[6]; for (int i = 0; i < 6; ++i) dC[i] = eps_new[i] - in.eps[i]; nFD = rmPiecesFD(mp, in, dC); }
    const bool bilin = (mp.tensionLaw == 1);

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

    const double bc = betaC(w, kp_new, mp);              // Eq.50 (P2f): scales the kdc1 plastic part
    const double wtw = tensileDamageWeight(mp, ac, depl, w, V);   // P2h ctTemper weight (1 if none)
    const double wcw = compressiveDamageWeight(mp, depl, w, V);   // PV20 tcTemper weight (1 if none)
    const double dnc = dnorm * wcw;                               // kdc1 plastic measure (== dnorm for none)
    const double etc = compressiveEquivStrain(mp, w, et);         // kdc drive (== et for none)
    double kdt1 = in.kdt1, kdt2 = in.kdt2, kdc = in.kdc, kdc1 = in.kdc1, kdc2 = in.kdc2;
    tensionHistUpdate(mp, kdt1, kdt2, et, et_max_n, dnorm, xs, wtw);
    const bool cdc = (mp.compDrive == 1);
    const double kdc_n = in.kdc, eqc_n = in.eqc, etp_n = in.etPrev;
    double eqcNew = eqc_n;
    if (cdc) {                                                           // B2: CDPM2 Eq.47-49
        compHistUpdateCdpm2(mp, kdc, kdc1, kdc2, eqcNew, etp_n, etc, ac, bc, dnc, xs);
    } else if (loading) {
        kdc  += ac * above;          kdc2 += ac * above / xs;
        kdc1 += ac * bc * dnc / xs;
    }
    const bool cAdv = cdc && (eqcNew > kdc_n);
    const double et_max2 = et_max_n > et ? et_max_n : et;

    // extreme effective principals + their eigenprojections (argmax/argmin; eig3sym is unsorted)
    int imax = 0, imin = 0;
    for (int i = 1; i < 3; ++i) { if (w[i] > w[imax]) imax = i; if (w[i] < w[imin]) imin = i; }
    const double Dt = (w[imax] > 1.0e-6 * mp.ft) ? mp.E * et : 0.0;   // P2i: E*et tensile drive (Eq.37)
    const double Dc = (-w[imin]) > 0.0 ? -w[imin] : 0.0;
    // P2g — MONOTONE drive (mirror damagedUpdate). Solve omega against the running max; tLoading/cLoading
    // mark whether each channel is ADVANCING its max (== loading). On UNLOAD (live drive < committed max)
    // the drive is frozen and the histories are frozen, so d(omega)/d(eps)=0 and the tangent collapses to
    // the SPD damage secant D_dam:C_eff. On loading max == live => byte-identical to the pre-P2g tangent.
    const double sigtMax = in.sigtMax > Dt ? in.sigtMax : Dt;
    const double sigcMax = in.sigcMax > Dc ? in.sigcMax : Dc;
    const bool tLoading = Dt >= in.sigtMax;
    const bool cLoading = Dc >= in.sigcMax;
    const double wt = (et_max2 > eps0 && sigtMax > 1.0e-6 * mp.ft)
                    ? omegaT(mp, kdt1, kdt2, sigtMax) : 0.0;
    const double wc = omegaC(mp, kdc, kdc1, kdc2, sigcMax, eps_fc);

    // D_dam = spectral derivative of the per-principal damaged stress with ω FROZEN
    // RESIDUAL TANGENT STIFFNESS (WP concrete3d-oracle-diagnosis): the bilinear tension law reaches omega_t = 1
    // EXACTLY at w = wf, so on a fully open crack every TENSILE principal direction (incl. the uniaxial-tension
    // laterals, whose effective stress sits on the Macaulay kink) has ZERO tangent stiffness => a singular
    // global system (single-element tension: NaN / runaway lateral strains, thousands of return-map
    // warnings). The TANGENT keeps (1-omega) >= OMEGA_TAN_FLOOR; the STRESS is untouched (exact zero), so
    // the dissipated energy and every stress fixture are unchanged — only omega within 1e-6 of 1 is affected.
    const double kT = (1.0 - wt) > OMEGA_TAN_FLOOR ? (1.0 - wt) : OMEGA_TAN_FLOOR;
    const double kC = (1.0 - wc) > OMEGA_TAN_FLOOR ? (1.0 - wc) : OMEGA_TAN_FLOOR;
    double yv[3], ypv[3];
    for (int a = 0; a < 3; ++a) {
        const double st = w[a] > 0.0 ? w[a] : 0.0;
        const double sc = w[a] < 0.0 ? w[a] : 0.0;
        yv[a]  = kT * st + kC * sc;
        ypv[a] = (w[a] > 0.0) ? kT : kC;
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

    // (C3c, cost) dω_t/dε and dω_c/dε are NON-zero only for an interior damage 0<ω<1 (see the IFT
    // branches below: a clamped or inactive ω is insensitive). Everything from here to the assembly
    // feeds ONLY dwt/dwc -- three micro-FD scalar gradients (36 eigendecompositions), the dnorm
    // gradient, and under loading two more composite FDs through the return map. Skipping it when
    // neither ω is interior (every ELASTIC point) leaves dwt = dwc = 0, i.e. the SAME assembly
    // arithmetic and a bit-identical tangent. Re-applied on the damage-drive kernel: the interior test
    // is the union of the IFT branch conditions (wt < OMEGA_MAX && bilin; wc < OMEGA_MAX && cdc; wc < 1
    // legacy), so `< 1.0` is a safe superset.
    double dwt[6] = {0, 0, 0, 0, 0, 0}, dwc[6] = {0, 0, 0, 0, 0, 0};
    const bool needDw = (wt > 0.0 && wt < 1.0) || (wc > 0.0 && wc < 1.0);
    if (needDw) {
    // --- chain-rule gradient pieces (each d(.)/dε, 6-vector) ---
    // Ceff^T @ g  (d(scalar of sig_eff)/dε = (d sig_eff/dε)^T @ (d scalar/d sig_eff))
    auto CeffT = [&](const double g[6], double out[6]) {
        for (int k = 0; k < 6; ++k) { double v = 0.0; for (int i = 0; i < 6; ++i) v += Ceff[i][k] * g[i]; out[k] = v; }
    };
    double g_et[6], g_xs[6], g_ac[6], det_deps[6], dxs_deps[6], dac_deps[6];
    dscalarDsig(0, sig_eff, mp, g_et);   CeffT(g_et, det_deps);
    dscalarDsig(1, sig_eff, mp, g_xs);   CeffT(g_xs, dxs_deps);
    dscalarDsig(2, sig_eff, mp, g_ac);   CeffT(g_ac, dac_deps);

    // drive gradients. P2i: the TENSILE drive is E*et, so d(Dt)/dε = E * det_deps (the equiv-strain
    // gradient already assembled above), NOT the extreme-principal eigenprojection. det_deps is a
    // per-Voigt-component micro-FD grad ⇒ it carries NO W6 weight (the W6 quirk applies to the tensor
    // eigenprojection only). The COMPRESSIVE drive stays the extreme principal: dDc/dε = -Ceff^T (W6 . Emin).
    // (P2g) each is frozen (zero) on unload so the -sig(x)d(omega) rank-update vanishes ⇒ SPD secant.
    double dDt_deps[6], dDc_deps[6];
    {
        if (Dt > 0.0 && tLoading) for (int i = 0; i < 6; ++i) dDt_deps[i] = mp.E * det_deps[i];
        else for (int i = 0; i < 6; ++i) dDt_deps[i] = 0.0;
        double En[3][3], En6[6], tmp[6];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) En[i][j] = V[i][imin] * V[j][imin];
        matToVoigt(En, En6);
        for (int i = 0; i < 6; ++i) tmp[i] = W6[i] * En6[i];
        if (Dc > 0.0 && cLoading) { CeffT(tmp, dDc_deps); for (int i = 0; i < 6; ++i) dDc_deps[i] = -dDc_deps[i]; }
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

    // d(beta_c)/dε (P2f): beta_c depends on BOTH rho_bar(sig_bar) AND qh2(kp), so a single composite
    // micro-FD THROUGH the return map captures the full gradient (mirror of the oracle; same step as the
    // numerical reference so the FD truncation correlates). Only needed under compressive loading.
    double dbc_deps[6] = {0,0,0,0,0,0};
    const bool needBc = cdc ? (cAdv && eqcNew > eps0) : loading;
    if (needBc && bc > 0.0) {
        double deps[6]; for (int i = 0; i < 6; ++i) deps[i] = eps_new[i] - in.eps[i];
        const double base = mp.fc / mp.E;
        for (int j = 0; j < 6; ++j) {
            const double hh = 1.0e-6 * (std::fabs(deps[j]) + base);
            double dp[6], dm[6]; for (int i = 0; i < 6; ++i) { dp[i] = deps[i]; dm[i] = deps[i]; }
            dp[j] += hh; dm[j] -= hh;
            double sbp[6], sbm[6], kpp, kpm, dum[6][6];
            rmForFD(mp, in, dp, nFD, sbp, kpp);
            rmForFD(mp, in, dm, nFD, sbm, kpm);
            double Ap[3][3], wp[3], Vp[3][3], Am[3][3], wm[3], Vm[3][3];
            voigtToMat(sbp, Ap); eig3sym(Ap, wp, Vp);
            voigtToMat(sbm, Am); eig3sym(Am, wm, Vm);
            dbc_deps[j] = (betaC(wp, kpp, mp) - betaC(wm, kpm, mp)) / (2.0 * hh);
        }
    }

    // d(w_t)/dε for the ctTemper modes (P2h): none -> 0 (unchanged tangent); alphat -> -d(alpha_c)/dε
    // (analytic, reuses dac_deps); proj -> composite micro-FD through the return map (w_t = the tensile-
    // stress-projected plastic-strain fraction). Only under loading.
    double dwtw_deps[6] = {0,0,0,0,0,0};
    const bool tAdv = bilin ? (det_raw > 0.0) : loading;         // the tensile history advances this step
    if (tAdv && mp.ctTemper == 1) {
        for (int i = 0; i < 6; ++i) dwtw_deps[i] = -dac_deps[i];
    } else if (tAdv && mp.ctTemper == 2) {
        double deps[6]; for (int i = 0; i < 6; ++i) deps[i] = eps_new[i] - in.eps[i];
        const double base = mp.fc / mp.E;
        double epln[6]; plasticStrain6(in.sigEff, in.eps, mp, epln);
        for (int j = 0; j < 6; ++j) {
            const double hh = 1.0e-6 * (std::fabs(deps[j]) + base);
            double dp[6], dm[6]; for (int i = 0; i < 6; ++i) { dp[i] = deps[i]; dm[i] = deps[i]; }
            dp[j] += hh; dm[j] -= hh;
            double sbp[6], sbm[6], kpp, kpm, dum[6][6];
            rmForFD(mp, in, dp, nFD, sbp, kpp);
            rmForFD(mp, in, dm, nFD, sbm, kpm);
            double Ap[3][3], wp[3], Vp[3][3], Am[3][3], wm[3], Vm[3][3];
            voigtToMat(sbp, Ap); eig3sym(Ap, wp, Vp);
            voigtToMat(sbm, Am); eig3sym(Am, wm, Vm);
            double eplp[6], eplm[6], dplp[6], dplm[6], epsp[6], epsm[6];
            for (int i = 0; i < 6; ++i) { epsp[i] = in.eps[i] + dp[i]; epsm[i] = in.eps[i] + dm[i]; }
            plasticStrain6(sbp, epsp, mp, eplp);
            plasticStrain6(sbm, epsm, mp, eplm);
            for (int i = 0; i < 6; ++i) { dplp[i] = eplp[i] - epln[i]; dplm[i] = eplm[i] - epln[i]; }
            const double wpv = tensileDamageWeight(mp, alphaCompression(wp), dplp, wp, Vp);
            const double wmv = tensileDamageWeight(mp, alphaCompression(wm), dplm, wm, Vm);
            dwtw_deps[j] = (wpv - wmv) / (2.0 * hh);
        }
    }

    // PV20 tcTemper proj: d(dnorm w_c)/dε (w_c by composite micro-FD through the return map, as the ctTemper proj
    // weight) and d(eps_tilde_c)/dε. none (or no tensile effective principal) -> dnc_deps == dnorm_deps,
    // detc_deps == det_deps (byte-identical).
    double dnc_deps[6], detc_deps[6];
    for (int i = 0; i < 6; ++i) { dnc_deps[i] = dnorm_deps[i]; detc_deps[i] = det_deps[i]; }
    {
        double mxw = w[0]; for (int a = 1; a < 3; ++a) if (w[a] > mxw) mxw = w[a];
        if (mp.tcTemper == 2 && mxw > TC_DEAD * mp.ft) {
            double deps[6]; for (int i = 0; i < 6; ++i) deps[i] = eps_new[i] - in.eps[i];
            const double base = mp.fc / mp.E;
            double epln[6]; plasticStrain6(in.sigEff, in.eps, mp, epln);
            double dwcw[6];
            for (int j = 0; j < 6; ++j) {
                const double hh = 1.0e-6 * (std::fabs(deps[j]) + base);
                double dp[6], dm[6]; for (int i = 0; i < 6; ++i) { dp[i] = deps[i]; dm[i] = deps[i]; }
                dp[j] += hh; dm[j] -= hh;
                double sbp[6], sbm[6], kpp, kpm, dum[6][6];
                rmForFD(mp, in, dp, nFD, sbp, kpp);
                rmForFD(mp, in, dm, nFD, sbm, kpm);
                double Ap[3][3], wp[3], Vp[3][3], Am[3][3], wm[3], Vm[3][3];
                voigtToMat(sbp, Ap); eig3sym(Ap, wp, Vp);
                voigtToMat(sbm, Am); eig3sym(Am, wm, Vm);
                double eplp[6], eplm[6], dplp[6], dplm[6], epsp[6], epsm[6];
                for (int i = 0; i < 6; ++i) { epsp[i] = in.eps[i] + dp[i]; epsm[i] = in.eps[i] + dm[i]; }
                plasticStrain6(sbp, epsp, mp, eplp);
                plasticStrain6(sbm, epsm, mp, eplm);
                for (int i = 0; i < 6; ++i) { dplp[i] = eplp[i] - epln[i]; dplm[i] = eplm[i] - epln[i]; }
                dwcw[j] = (compressiveDamageWeight(mp, dplp, wp, Vp) - compressiveDamageWeight(mp, dplm, wm, Vm))
                        / (2.0 * hh);
            }
            for (int i = 0; i < 6; ++i) dnc_deps[i] = wcw * dnorm_deps[i] + dnorm * dwcw[i];
            double g_etc[6]; dscalarDsig(3, sig_eff, mp, g_etc); CeffT(g_etc, detc_deps);
        }
    }

    // dkd*/dε under loading (else zero).  kdc1 = ac * bc * dnorm / xs (Eq.48 with beta_c) => product rule.
    // kdt2 = w_t * above / xs, kdt1 = w_t * dnorm / xs (Eq.45/44 with the ctTemper weight) => product rule.
    double dkdt1[6], dkdt2[6], dkdc1[6], dkdc2[6];
    for (int i = 0; i < 6; ++i) { dkdt1[i] = dkdt2[i] = 0.0; }
    if (bilin && det_raw > 0.0) {
        // literal Eq.45/44: kdt2 += w det_raw/xs (history advancing), kdt1 += w frac dnorm/xs (past onset)
        const double ixs = 1.0 / xs, ixs2 = ixs * ixs;
        const bool cross = et_max_n < eps0;
        const double frac = cross ? (et - eps0) / det_raw : 1.0;
        const double dfr = cross ? (eps0 - et_max_n) / (det_raw * det_raw) : 0.0;   // d frac / d et
        for (int i = 0; i < 6; ++i) {
            dkdt2[i] = wtw * (det_deps[i] * ixs - det_raw * dxs_deps[i] * ixs2) + (det_raw * ixs) * dwtw_deps[i];
            if (et > eps0)
                dkdt1[i] = wtw * (frac * dnorm_deps[i] * ixs + dnorm * dfr * det_deps[i] * ixs
                                  - frac * dnorm * dxs_deps[i] * ixs2) + (frac * dnorm * ixs) * dwtw_deps[i];
        }
    }
    if (loading) {
        const double ixs = 1.0 / xs, ixs2 = ixs * ixs;
        for (int i = 0; i < 6; ++i) {
            if (!bilin) {
                dkdt2[i] = wtw * (det_deps[i] * ixs - above * dxs_deps[i] * ixs2) + (above * ixs) * dwtw_deps[i];
                dkdt1[i] = wtw * (dnorm_deps[i] * ixs - dnorm * dxs_deps[i] * ixs2) + (dnorm * ixs) * dwtw_deps[i];
            }
            dkdc2[i] = (dac_deps[i] * above + ac * det_deps[i]) * ixs - ac * above * dxs_deps[i] * ixs2;
            dkdc1[i] = (dac_deps[i] * bc * dnc + ac * dbc_deps[i] * dnc + ac * bc * dnc_deps[i]) * ixs
                     - ac * bc * dnc * dxs_deps[i] * ixs2;
        }
    } else for (int i = 0; i < 6; ++i) { if (!bilin) { dkdt1[i] = dkdt2[i] = 0.0; } dkdc1[i] = dkdc2[i] = 0.0; }
    double dkdc[6] = {0, 0, 0, 0, 0, 0};
    if (cdc) {   // B2: kappa_dc = eqc_n + ac (et - etp_n) while advancing; kdc2 += d/xs; kdc1 += ac bc frac dnorm/xs
        for (int i = 0; i < 6; ++i) { dkdc1[i] = dkdc2[i] = 0.0; }
        if (cAdv) {
            const double dd = eqcNew - kdc_n, ixs = 1.0 / xs;
            for (int i = 0; i < 6; ++i) {
                dkdc[i] = ac * detc_deps[i] + (etc - etp_n) * dac_deps[i];
                dkdc2[i] = dkdc[i] * ixs - dd * dxs_deps[i] * ixs * ixs;
            }
            if (eqcNew > eps0) {
                const bool cross = kdc_n < eps0;
                const double frac = cross ? (eqcNew - eps0) / dd : 1.0;
                const double dfr = cross ? (eps0 - kdc_n) / (dd * dd) : 0.0;
                const double inc = ac * bc * frac * dnc * ixs;
                for (int i = 0; i < 6; ++i)
                    dkdc1[i] = (dac_deps[i] * bc * frac * dnc + ac * dbc_deps[i] * frac * dnc
                                + ac * bc * dfr * dkdc[i] * dnc + ac * bc * frac * dnc_deps[i]) * ixs
                             - inc * dxs_deps[i] * ixs;
            }
        }
    }

    // ω via IFT (only when interior 0<ω<1; clamped/inactive ω is insensitive => dω=0)
    // (P2g) D = the MONOTONE drive sigtMax/sigcMax; dDt_deps/dDc_deps are already zeroed on unload, so an
    // unloading channel contributes d(omega)=0 (secant). On loading D == live drive => unchanged.
    if (wt > 0.0 && wt < OMEGA_MAX && bilin) {
        // IFT on F(w) = (1-w)D - sigma(h(kd1+w kd2)): F_w = -D - sigma' h kd2 (sigma' = active branch slope)
        const double wf_b = mp.Gf / (BILIN_GF * mp.ft);
        double slope = 0.0; bilinearSigma(mp.lch * (kdt1 + wt * kdt2), mp.ft, wf_b, &slope);
        const double spd = slope * mp.lch;
        const double Fw = -sigtMax - spd * kdt2;
        for (int i = 0; i < 6; ++i)
            dwt[i] = (-(1.0 - wt) / Fw) * dDt_deps[i] + (spd / Fw) * dkdt1[i] + (spd * wt / Fw) * dkdt2[i];
    } else if (wt > 0.0 && wt < 1.0) {
        const double Ht = sigtMax * ((1.0 - wt) * kdt2 / eps_f - 1.0);
        const double a0 = -(1.0 - wt) / Ht;
        const double a1 = -(1.0 - wt) * sigtMax / (eps_f * Ht);
        const double a2 = -(1.0 - wt) * sigtMax * wt / (eps_f * Ht);
        for (int i = 0; i < 6; ++i) dwt[i] = a0 * dDt_deps[i] + a1 * dkdt1[i] + a2 * dkdt2[i];
    } else for (int i = 0; i < 6; ++i) dwt[i] = 0.0;
    if (wc > 0.0 && wc < OMEGA_MAX && cdc) {     // B2: D = E kappa_dc (the IFT does not depend on f = ft)
        const double Dcc = mp.E * kdc;
        const double Hc = Dcc * ((1.0 - wc) * kdc2 / eps_fc - 1.0);
        const double a0 = -(1.0 - wc) / Hc;
        const double a1 = -(1.0 - wc) * Dcc / (eps_fc * Hc);
        const double a2 = -(1.0 - wc) * Dcc * wc / (eps_fc * Hc);
        for (int i = 0; i < 6; ++i) dwc[i] = a0 * mp.E * dkdc[i] + a1 * dkdc1[i] + a2 * dkdc2[i];
    } else if (wc > 0.0 && wc < 1.0) {
        const double Hc = sigcMax * ((1.0 - wc) * kdc2 / eps_fc - 1.0);
        const double a0 = -(1.0 - wc) / Hc;
        const double a1 = -(1.0 - wc) * sigcMax / (eps_fc * Hc);
        const double a2 = -(1.0 - wc) * sigcMax * wc / (eps_fc * Hc);
        for (int i = 0; i < 6; ++i) dwc[i] = a0 * dDc_deps[i] + a1 * dkdc1[i] + a2 * dkdc2[i];
    } else for (int i = 0; i < 6; ++i) dwc[i] = 0.0;

    }   // needDw (C3c)

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
inline int returnMapTensor1(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                            bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                            bool doTangent, PieceSens* ps = nullptr);

inline int returnMapTensorAdaptive(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                                   bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                                   bool doTangent, int* subInfo);

// Dimensionless yield-function value f_tr of the ELASTIC TRIAL (the same f_tr returnMapHardening tests).
inline double trialOvershoot(const Params& mp, const double sig_n[6], const double deps[6], double kp_n)
{
    double sig_tr[6]; elasticPredTensor(sig_n, deps, mp, sig_tr);
    double A[3][3], w[3], V[3][3]; voigtToMat(sig_tr, A); eig3sym(A, w, V);
    double sv[6] = {w[0], w[1], w[2], 0.0, 0.0, 0.0}, xi, rho, th;
    invariants(sv, xi, rho, th);
    return yfInvHard(xi, rho, lodeR(th, mp.e), kp_n, mp);
}

// The deterministic piece count n = clamp(ceil(f_tr/c), 1, nmax) of a (state, increment); Params::subIncrForceN > 0 pins it.
inline int detPieces(const Params& mp, const double sig_n[6], const double deps[6], double kp_n)
{
    if (mp.subIncrForceN > 0) return mp.subIncrForceN;
    const double c = (mp.subIncrC > 0.0) ? mp.subIncrC : 0.3;
    const int nmaxP = (mp.subIncrMaxPieces > 0) ? mp.subIncrMaxPieces : 64;
    const double fTr = trialOvershoot(mp, sig_n, deps, kp_n);
    int n = 1;
    if (fTr > 0.0) { const double q = std::ceil(fTr / c); n = (q >= (double)nmaxP) ? nmaxP : (q < 1.0 ? 1 : (int)q); }
    return n;
}

// DETERMINISTIC sub-incrementation (mirror of the oracle _return_map_tensor_det): n = clamp(ceil(f_tr/c), 1, nmax)
// equal pieces, each a direct return, chained from the committed state, ALWAYS applied -- one (state, increment)
// always takes the same path. A piece failure redoes the whole chain with 2n and then 4n pieces; the honest failure
// (the direct-return fallback, status != 0) only after that. Structurally bounded: n + 2n + 4n <= 7 nmax direct
// returns. n = 1 (f_tr <= c, incl. elastic trials) is one direct return, byte-identical to the direct map.
// subInfo: 0 = direct, L >= 2 = the ladder level (piece count) that succeeded, -1 = final failure.
//
// DISCONTINUITY (honest statement; the earlier text here said the map was discontinuous only at ladder failures, which is
// false): a chain of n pieces and a chain of n + 1 pieces are two different, each consistent, integrations of the same
// increment, so sigma_eff and kappa_p JUMP where ceil(f_tr/c) changes, and again wherever the ladder switches level. The
// jump is the discretization difference between the two integrations -- measured at 1e-19-apart strains on the
// reviewer's probes: |dsigma_eff| 0.06-0.40 MPa (0.4-1.5 % of |sigma_eff|), kappa_p up to ~2.5 near first cracking.
// Chosen over the failure-driven adaptive path because that one is discontinuous at every attempt boundary and depends on
// the Newton iterate's noise (1.8 % vs 13.3 % measured), and because the deterministic map is a function of (state,
// increment) only. Consumers that pin equality of nominally identical Gauss points (EAS alpha, hosting parity) must
// allow the jump. The reported TANGENT is the CHAIN's (accumulated forward through the pieces, see the loop), not the
// last piece's; Params::subIncrForceN pins n for finite-difference references of it.
inline int returnMapTensorDet(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                              double sig_new[6], double& kp_new, double Dtan6[6][6], bool doTangent, int* subInfo)
{
    const int n = detPieces(mp, sig_n, deps, kp_n);
    if (n == 1) {
        const int st = returnMapTensor1(mp, sig_n, deps, kp_n, true, sig_new, kp_new, Dtan6, doTangent);
        if (st == 0) { if (subInfo) *subInfo = 0; return 0; }
    }
    const int levels[3] = { n, 2 * n, 4 * n };
    for (int li = (n == 1 ? 1 : 0); li < 3; ++li) {
        const int L = levels[li];
        double s[6], k = kp_n, sn[6], kn, Dt[6][6], G[6][6], Gn[6][6], Cinv[6][6], C0m[6][6], g[6], gn[6];
        PieceSens ps;
        for (int i = 0; i < 6; ++i) s[i] = sig_n[i];
        if (doTangent) {
            for (int i = 0; i < 6; ++i) { g[i] = 0.0; for (int j = 0; j < 6; ++j) G[i][j] = 0.0; }
            elasticC(mp, C0m); invert6(C0m, Cinv);
        }
        bool ok = true;
        for (int p = 0; p < L && ok; ++p) {
            double d[6]; for (int i = 0; i < 6; ++i) d[i] = deps[i] / L;
            if (returnMapTensor1(mp, s, d, k, true, sn, kn, Dt, doTangent, doTangent ? &ps : nullptr) != 0) ok = false;
            else {
                for (int i = 0; i < 6; ++i) s[i] = sn[i];
                k = kn;
                if (doTangent) {
                    // CHAIN tangent, accumulated forward through the pieces. Piece p is sig_p = R(sig_(p-1) + C d/L, kappa_(p-1)),
                    // kappa_p = K(same arguments), with the direct-return sensitivities Dt (= d sig_p/d(d/L)), Sk, R, Kk:
                    //   d sig_p/d sig_(p-1) = Dt C^-1,   d sig_p/d kappa_(p-1) = Sk,
                    //   d kappa_p/d sig_(p-1) = R,       d kappa_p/d kappa_(p-1) = Kk,     d kappa_p/d(d/L) = R C
                    // and, with G = d sig/d(deps) (6x6) and g = d kappa/d(deps) (1x6), G_0 = 0, g_0 = 0:
                    //   G_p = (Dt C^-1) G_(p-1) + Sk (x) g_(p-1) + Dt/L ,   g_p = R G_(p-1) + Kk g_(p-1) + (R C)/L .
                    // This is the EXACT derivative of the fixed-n chain (the kappa history coupling included); the only
                    // approximations left are the Lode-angle scalar central differences inside the principal Jacobian
                    // and, in a ladder/rescue state, the chain that was actually integrated.
                    double A[6][6], RC[6];
                    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
                        double v = 0.0; for (int q = 0; q < 6; ++q) v += Dt[i][q] * Cinv[q][j];
                        A[i][j] = v;
                    }
                    for (int j = 0; j < 6; ++j) { double v = 0.0; for (int q = 0; q < 6; ++q) v += ps.R[q] * C0m[q][j]; RC[j] = v; }
                    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
                        double v = Dt[i][j] / L + ps.Sk[i] * g[j]; for (int q = 0; q < 6; ++q) v += A[i][q] * G[q][j];
                        Gn[i][j] = v;
                    }
                    for (int j = 0; j < 6; ++j) {
                        double v = ps.Kk * g[j] + RC[j] / L; for (int q = 0; q < 6; ++q) v += ps.R[q] * G[q][j];
                        gn[j] = v;
                    }
                    for (int i = 0; i < 6; ++i) { g[i] = gn[i]; for (int j = 0; j < 6; ++j) G[i][j] = Gn[i][j]; }
                }
            }
        }
        if (ok) {
            for (int i = 0; i < 6; ++i) sig_new[i] = s[i];
            kp_new = k;
            if (doTangent) for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Dtan6[i][j] = G[i][j];
            if (subInfo) *subInfo = L;
            return 0;
        }
    }
    // ladder exhausted: last resort = the ADAPTIVE path (direct, then halving; bounded by maxSubAttempts), then the
    // honest failure. Reached only in the small share of states the fixed-n chains cannot integrate at all, where the
    // alternative is a step cut. (The deterministic map is NOT continuous overall: see the DISCONTINUITY note above.)
    if (mp.subIncrRescue)
        return returnMapTensorAdaptive(mp, sig_n, deps, kp_n, true, sig_new, kp_new, Dtan6, doTangent, subInfo);
    const int st0 = returnMapTensor1(mp, sig_n, deps, kp_n, true, sig_new, kp_new, Dtan6, doTangent);
    if (subInfo) *subInfo = (st0 == 0) ? 0 : -1;
    return st0;
}

// Tensor return with OOFEM-style SUB-INCREMENTATION (B1; mirror of the oracle return_map_tensor). The direct
// return first; if it fails and mp.maxSubIncr > 0 (hardening map), the strain increment is halved and the
// sub-increments integrated in sequence from the committed state (doubling back, to at most 2x the size that
// just succeeded, after each success) down to 2^-maxSubIncr. The reported tangent is the LAST sub-increment's
// consistent tangent (an approximation of the sub-stepped algorithmic tangent). maxSubIncr = 0 => byte-identical
// to the direct return.
//
// ATTEMPT BUDGET (WP concrete3d-hang-diagnosis, 2026-09-27): mp.maxSubAttempts caps the TOTAL number of
// returnMapTensor1 calls in this loop (successes + failures). Without it, a GP that only ever succeeds at the
// 2^-maxSubIncr floor alternates fail/succeed and can legally run ~2 * 2^maxSubIncr attempts (~2048 at the
// CDPM2 wrapper's maxSubIncr=10) of up to 100 Newton iterations each, PER material point PER global Newton
// iteration -- with nothing ever failing, so nothing is logged and the step is never cut. That is the observed
// analyze(1) hang (C3/B1 vecchio_shim, sheikh_uzumeri confined column): 100% CPU, no opserr output, no
// recorder output, for 30-40+ minutes on a single step. Once the budget is exhausted, return the honest
// failure (st0, the direct-return fallback) exactly like any other non-convergence, so the caller's status
// != 0 cuts the step. maxSubIncr = 0 stays byte-identical (the budget only applies inside this loop).
inline int returnMapTensor(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                           bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                           bool doTangent, int* subInfo)
{
    // subInfo (diagnostic, optional): 0 = direct return, n >= 2 = sub-incremented in n pieces, -1 = FINAL failure
    if (mp.maxSubIncr > 0 && hardening && mp.subIncrMode == 0)
        return returnMapTensorDet(mp, sig_n, deps, kp_n, sig_new, kp_new, Dtan6, doTangent, subInfo);
    return returnMapTensorAdaptive(mp, sig_n, deps, kp_n, hardening, sig_new, kp_new, Dtan6, doTangent, subInfo);
}

// ADAPTIVE sub-incrementation (subIncrMode = 1, and the deterministic mode's last resort): direct return first, on
// failure halving/doubling down to 2^-maxSubIncr, bounded by maxSubAttempts. See the comments above.
inline int returnMapTensorAdaptive(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                                   bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                                   bool doTangent, int* subInfo)
{
    const int st0 = returnMapTensor1(mp, sig_n, deps, kp_n, hardening, sig_new, kp_new, Dtan6, doTangent);
    if (subInfo) *subInfo = (st0 == 0) ? 0 : -1;
    if (st0 == 0 || mp.maxSubIncr <= 0 || !hardening) return st0;
    const int maxAttempts = (mp.maxSubAttempts > 0) ? mp.maxSubAttempts : 64;
    int pieces = 0, attempts = 0;
    double s[6], k = kp_n, done = 0.0, frac = 0.5;
    const double floorFrac = std::ldexp(1.0, -mp.maxSubIncr);
    for (int i = 0; i < 6; ++i) s[i] = sig_n[i];
    double sn[6], kn, Dt[6][6];
    while (done < 1.0) {
        if (attempts >= maxAttempts) return st0;        // honest failure: budget exhausted, keep the fallback
        ++attempts;
        const double f = (frac < 1.0 - done) ? frac : 1.0 - done;
        double d[6]; for (int i = 0; i < 6; ++i) d[i] = deps[i] * f;
        const int st = returnMapTensor1(mp, s, d, k, true, sn, kn, Dt, doTangent && (done + f >= 1.0));
        if (st == 0) {
            for (int i = 0; i < 6; ++i) s[i] = sn[i];
            k = kn; done += f; ++pieces;
            frac = (2.0 * f < 1.0) ? 2.0 * f : 1.0;     // at most 2x the size that just succeeded (not the stale frac)
        } else {
            frac *= 0.5;
            if (frac < floorFrac) return st0;          // honest failure: keep the direct-return fallback
        }
    }
    for (int i = 0; i < 6; ++i) sig_new[i] = s[i];
    kp_new = k;
    if (doTangent) for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Dtan6[i][j] = Dt[i][j];
    if (subInfo) *subInfo = pieces;
    return 0;
}

inline int returnMapTensor1(const Params& mp, const double sig_n[6], const double deps[6], double kp_n,
                            bool hardening, double sig_new[6], double& kp_new, double Dtan6[6][6],
                            bool doTangent, PieceSens* ps)
{
    if (ps) *ps = PieceSens();
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

    if (doTangent) consistentTangent(sig_tr, w, V, pr, mp, hardening, Dtan6, ps);
    return pr.converged ? 0 : 2;   // 0 OK, 2 = no-converge (honest flag)
}

// ===========================================================================
// DEAD-AWARE EFFECTIVE-stress return, shared by returnMap and driveConfinedFiber (review M2: the BeamFiber view called
// returnMapTensor/damagedUpdate directly and so never saw the dead-point treatment). Precondition: the point is not
// CRUSHED (inRaw.wc < mp.omegaDead; a crushed point is frozen by the caller). Decided on the COMMITTED omega_t:
// CRACKED (omega_t dead): tension cutoff on the plastic flow. sig_tr = sig_n + C:deps is split spectrally; the tensile part
// is carried elastically, the return map runs on the compressive remainder with a zero increment. The committed state is
// re-based on the trial (eps = new strain, sigEff = sig_tr) so the plastic-strain increment the damage update sees is the
// compressive return's alone. (Absorbing the tension into the plastic strain instead was measured and rejected: the
// permanent strain locks a full-stiffness compression in on unloading.)
// Outputs: inEff -> the committed state to hand to damagedUpdate/damagedTangent (inRaw itself, or the re-based copy in
// cutBuf); cutT; sig_eff/kp_new/Dtan6 = the EFFECTIVE stress, kappa_p and effective tangent; return = the return-map status.
// Includes the Duvaut-Lions relaxation (only when !implex && eta > 0 && dt > 0).
// ===========================================================================
inline int effectiveReturn(const Params& mp, const double strain[6], const State& inRaw, bool hardening, double dt,
                           bool doTangent, State& cutBuf, const State*& inEff, bool& cutT,
                           double sig_eff[6], double& kp_new, double Dtan6[6][6], int* subInfo)
{
    double deps[6];
    for (int i = 0; i < 6; ++i) deps[i] = strain[i] - inRaw.eps[i];
    cutT = (inRaw.wt >= mp.omegaDead);
    double cutW[3] = {0, 0, 0}, cutV[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double sigPlus[6] = {0, 0, 0, 0, 0, 0}, sigMinus[6] = {0, 0, 0, 0, 0, 0};
    if (cutT) {
        double sigTr[6]; elasticPredTensor(inRaw.sigEff, deps, mp, sigTr);
        double A[3][3]; voigtToMat(sigTr, A); eig3sym(A, cutW, cutV);
        double sm[3], sq[3];
        for (int a = 0; a < 3; ++a) { sm[a] = cutW[a] < 0.0 ? cutW[a] : 0.0; sq[a] = cutW[a] > 0.0 ? cutW[a] : 0.0; }
        double S[3][3], Q[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
            double v = 0.0, u = 0.0;
            for (int a = 0; a < 3; ++a) { v += cutV[i][a] * sm[a] * cutV[j][a]; u += cutV[i][a] * sq[a] * cutV[j][a]; }
            S[i][j] = v; Q[i][j] = u;
        }
        cutBuf = inRaw;
        matToVoigt(S, sigMinus); matToVoigt(Q, sigPlus);
        for (int i = 0; i < 6; ++i) { cutBuf.sigEff[i] = sigTr[i]; cutBuf.eps[i] = strain[i]; deps[i] = 0.0; }
    }
    inEff = cutT ? &cutBuf : &inRaw;
    const State& in = *inEff;
    // (1) IMPLICIT EFFECTIVE-stress return from the committed EFFECTIVE state (NOT the nominal sig).
    int status = returnMapTensor(mp, cutT ? sigMinus : in.sigEff, deps, in.kp, hardening, sig_eff, kp_new, Dtan6, doTangent,
                                 subInfo);
    if (cutT) for (int i = 0; i < 6; ++i) sig_eff[i] += sigPlus[i];   // tensile part carried elastically
    // (1b) Duvaut-Lions viscoplastic relaxation at the PLASTIC level (ADR §4.4; oracle PR #316). Relax the
    //   inviscid effective return + kp toward the elastic trial by beta = dt/(eta+dt) (Simo-Hughes closed
    //   form). beta < 1 only with a positive viscosity AND a positive dt; eta==0 OR dt<=0 => beta=1 =>
    //   BYTE-identical to the inviscid Tier-1 path (a missing time increment falls back to inviscid, NOT
    //   to the elastic beta->0 limit). Damage then follows from the RELAXED effective stress (downstream
    //   uses sig_eff/kp_new), and the EFFECTIVE consistent tangent blends C_eff <- (1-beta)C0 + beta C_eff
    //   (damagedTangent chains its damage linearization through this blended C_eff). v1: Tier-1 only —
    //   gated on !implex so the IMPL-EX implicit solve stays inviscid (matches the oracle scope; the
    //   -eta + -implex composition is deferred).
    if (!mp.implex && mp.eta > 0.0 && dt > 0.0) {
        const double beta = dt / (mp.eta + dt);
        double sig_tr[6];
        elasticPredTensor(in.sigEff, deps, mp, sig_tr);
        for (int i = 0; i < 6; ++i) sig_eff[i] = (1.0 - beta) * sig_tr[i] + beta * sig_eff[i];
        kp_new = (1.0 - beta) * in.kp + beta * kp_new;
        if (doTangent && status == 0) {
            double C0[6][6]; elasticC(mp, C0);
            for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
                Dtan6[i][j] = (1.0 - beta) * C0[i][j] + beta * Dtan6[i][j];
        }
    }
    if (cutT && doTangent && status == 0) {
        // sig_eff = Plus(sig_tr) + RM(Minus(sig_tr)):  d sig_eff/d eps = C0 + (Dret C0^-1 - I) Ddam- C0, with Dret = A C0 the
        // return map's own tangent at the compressive remainder and Ddam- the spectral derivative of the negative-part map.
        double yv[3], ypv[3], Ddam[6][6], C0[6][6], C0i[6][6], T1[6][6], T2[6][6];
        for (int a = 0; a < 3; ++a) { yv[a] = cutW[a] < 0.0 ? cutW[a] : 0.0; ypv[a] = cutW[a] < 0.0 ? 1.0 : 0.0; }
        isotropicTangent(cutW, cutV, yv, ypv, Ddam);
        elasticC(mp, C0); invert6(C0, C0i);
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            double v = 0.0; for (int k = 0; k < 6; ++k) v += Dtan6[i][k] * C0i[k][j];
            T1[i][j] = v - (i == j ? 1.0 : 0.0);                           // Dret C0^-1 - I
        }
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            double v = 0.0; for (int k = 0; k < 6; ++k) v += Ddam[i][k] * C0[k][j];
            T2[i][j] = v;                                                  // Ddam- C0
        }
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            double v = C0[i][j]; for (int k = 0; k < 6; ++k) v += T1[i][k] * T2[k][j];
            Dtan6[i][j] = v;
        }
    }
    return status;
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
inline int returnMap(const Params& mp, const double strain[6], const State& inRaw, State& out,
                     double sigma[6], double sigEffImplicit[6], double Dtan6[6][6],
                     bool doTangent, double dt = 0.0, bool hardening = true,
                     double* wtOut = nullptr, double* wcOut = nullptr)
{
    // (0) DEAD POINTS (WP concrete3d-hang-diagnosis #877 follow-up, owner decision 2026-09-28; see Params::omegaDead).
    // Mirror of the oracle damaged_step_tensor. Decided on the COMMITTED damage.
    double deps[6];
    for (int i = 0; i < 6; ++i) deps[i] = strain[i] - inRaw.eps[i];
    if (inRaw.wc >= mp.omegaDead) {
        // CRUSHED: strict freeze. The point is dead in every direction: both damages go to the floor (OMEGA_MAX) and stay
        // there, nominal = (1-OMEGA_MAX)*sig_eff (scalar), tangent (1-OMEGA_MAX)*C, plastic state and histories frozen.
        // Review #877 minor 2 asked to freeze at the COMMITTED omega instead (no jump). MEASURED and rejected: with the
        // committed omega (>= omegaDead) the residual (1-omega)*sig_eff GROWS with the elastic sig_eff of the frozen point and
        // the Gc-calibration gate (uniaxial compression dissipates Gc within 5 %) goes 62 % off, because the post-peak tail
        // never decays below 1 % of the peak. The jump to the floor is small in ABSOLUTE terms -- the nominal stress drops
        // by at most (1-omegaDead)*|sig_eff,dead| = 2e-3 |sig_eff| at the freeze (~0.06 MPa at fc = 30), 1e-2 at the loosest
        // admissible threshold -- which is what a dead point should do (the law's tail is exhausted to that fraction).
        double sigEffD[6]; elasticPredTensor(inRaw.sigEff, deps, mp, sigEffD);
        const double k = 1.0 - OMEGA_MAX;
        out = inRaw;
        out.wt = OMEGA_MAX; out.wc = OMEGA_MAX;
        for (int i = 0; i < 6; ++i) {
            out.eps[i] = strain[i]; out.sigEff[i] = sigEffD[i]; sigEffImplicit[i] = sigEffD[i];
            out.sig[i] = k * sigEffD[i]; sigma[i] = out.sig[i]; out.depl[i] = 0.0;
        }
        out.subInfo = 0; out.dwt = 0.0; out.dwc = 0.0;
        out.dt_n = (dt > 0.0) ? dt : inRaw.dt_n;
        if (wtOut) *wtOut = OMEGA_MAX;
        if (wcOut) *wcOut = OMEGA_MAX;
        if (doTangent) {
            const double kT = k > OMEGA_TAN_FLOOR ? k : OMEGA_TAN_FLOOR;
            double C0[6][6]; elasticC(mp, C0);
            for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Dtan6[i][j] = kT * C0[i][j];
        }
        return 0;
    }
    // (1)+(1b) dead-aware EFFECTIVE-stress return (CRACKED points: tension cutoff on the plastic flow; see effectiveReturn)
    State cutBuf; const State* inP = nullptr; bool cutT = false;
    double sig_eff[6], kp_new;
    int status = effectiveReturn(mp, strain, inRaw, hardening, dt, doTangent, cutBuf, inP, cutT, sig_eff, kp_new, Dtan6,
                                 &out.subInfo);
    const State& in = *inP;
    for (int i = 0; i < 6; ++i) { out.eps[i] = strain[i]; out.sigEff[i] = sig_eff[i]; sigEffImplicit[i] = sig_eff[i]; }
    out.kp = kp_new;
    // (2) IMPLICIT P2 dual-damage NOMINAL stress (writes out.sig + the damage history). Unilateral by re-split.
    double wt_impl = 0.0, wc_impl = 0.0;
    damagedUpdate(mp, in, sig_eff, kp_new, strain, out, &wt_impl, &wc_impl);
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
            damagedTangent(mp, in, sig_eff, strain, kp_new, Ceff, Dtan6);
        }
        return status;
    }

    // ---- TIER-2 IMPL-EX: report the EXPLICIT extrapolated stress + the degraded-elastic SECANT ----
    double r = (in.dt_n > 0.0 && dt > 0.0) ? (dt / in.dt_n) : 0.0;   // forward-only; floored at 0
    if (r > mp.implexRmax) r = mp.implexRmax;                        // clamp (review ALG-2/NUM-2/NUM-3)
    double wt_x = in.wt + r * in.dwt; if (wt_x < 0.0) wt_x = 0.0; if (wt_x > 1.0 - 1.0e-12) wt_x = 1.0 - 1.0e-12;
    double wc_x = in.wc + r * in.dwc; if (wc_x < 0.0) wc_x = 0.0; if (wc_x > 1.0 - 1.0e-12) wc_x = 1.0 - 1.0e-12;
    double deps_eff[6];
    for (int i = 0; i < 6; ++i) deps_eff[i] = (strain[i] - inRaw.eps[i]) - r * inRaw.depl[i];   // frozen plastic-strain increment
    double sig_bar_x[6];
    elasticPredTensor(inRaw.sigEff, deps_eff, mp, sig_bar_x);               // LINEAR in deps => elastic tangent
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

// The deterministic piece count that returnMap(mp, strain, inRaw, ...) uses for its (first) return map -- for tests and
// finite-difference references that must pin n on both legs (Params::subIncrForceN). 1 for a crushed (frozen) point.
inline int returnMapPieces(const Params& mp, const State& inRaw, const double strain[6])
{
    double deps[6];
    for (int i = 0; i < 6; ++i) deps[i] = strain[i] - inRaw.eps[i];
    if (inRaw.wc >= mp.omegaDead) return 1;
    if (inRaw.wt >= mp.omegaDead) {
        State rebased = inRaw;                       // rmPiecesFD reads the re-based trial state like damagedTangent does
        double sigTr[6]; elasticPredTensor(inRaw.sigEff, deps, mp, sigTr);
        for (int i = 0; i < 6; ++i) { rebased.sigEff[i] = sigTr[i]; rebased.eps[i] = strain[i]; }
        const double zero[6] = {0, 0, 0, 0, 0, 0};
        return rmPiecesFD(mp, rebased, zero);
    }
    return detPieces(mp, inRaw.sigEff, deps, inRaw.kp);
}

// ===========================================================================
// P5 — CONFINED-FIBER VIEW (§4.6, "Mander by mechanism"). The same triaxial kernel condensed into a
// 1-D beam-column fiber against the PASSIVE transverse-steel confinement. The axial strain (the retained
// BeamFiber comps {0,3,5} = axial 00 + the two transverse shears 01,02) is imposed; the lateral block
// {1,2,4} (the two lateral normals 11,22 + their in-plane shear 12) is condensed by a nested Newton so
// the EFFECTIVE lateral normal stresses balance the passive hoop and the in-plane shear is traction-free:
//      sigEff_11 + sig_hoop(eps_11) = 0,   sigEff_22 + sig_hoop(eps_22) = 0,   sigEff_12 = 0.
// The non-associated dilatancy mobilizes the hoop tension self-consistently => confined strength +
// ductility EMERGE (NO pre-baked Mander backbone). Mirrors the oracle confined_step EXACTLY (symmetric
// eps_11=eps_22, shears 0); hoopK=0 reduces to the free uniaxial-stress driver byte-for-byte.
// Circular/spiral hoops only (the symmetric two-normal condensation; rectangular ties are anisotropic).
// NB the lateral balance is on the EFFECTIVE stress (matches the shipped P5a oracle); at/near the peak
// the damage is ~0 so nominal==effective and the Mander match is unaffected (record in LEDGER_quirks).
// ---------------------------------------------------------------------------

// Passive circular-hoop spring: tension-only (slack in lateral contraction, eps_lat<=0), elastic-
// perfectly-plastic (stiffness K, yield fy). Returns the confining pressure sig_hoop >= 0 (the concrete's
// effective lateral stress balances -sig_hoop). Mirrors the oracle hoop_stress.
inline double hoopStress(double epsLat, double K, double fy)
{
    if (epsLat <= 0.0) return 0.0;
    const double s = K * epsLat;
    return (s < fy) ? s : fy;
}
// d sig_hoop / d eps_lat: K below yield (and in tension), 0 once yielded or slack. Oracle hoop_stiffness.
inline double hoopStiffness(double epsLat, double K, double fy)
{
    return (epsLat > 0.0 && K * epsLat < fy) ? K : 0.0;
}

// Solve the 3x3 system A x = b by Gauss elimination with partial pivoting. false on a singular pivot.
inline bool solve3(const double A[3][3], const double b[3], double x[3])
{
    double M[3][4];
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) M[i][j] = A[i][j]; M[i][3] = b[i]; }
    for (int c = 0; c < 3; ++c) {
        int piv = c; double best = std::fabs(M[c][c]);
        for (int r = c + 1; r < 3; ++r) if (std::fabs(M[r][c]) > best) { best = std::fabs(M[r][c]); piv = r; }
        if (best < 1.0e-300) return false;
        if (piv != c) for (int j = 0; j < 4; ++j) { double t = M[piv][j]; M[piv][j] = M[c][j]; M[c][j] = t; }
        for (int r = 0; r < 3; ++r) {
            if (r == c) continue;
            const double f = M[r][c] / M[c][c];
            for (int j = c; j < 4; ++j) M[r][j] -= f * M[c][j];
        }
    }
    for (int i = 0; i < 3; ++i) x[i] = M[i][3] / M[i][i];
    return true;
}

// Confined-fiber update. `strain` holds the imposed retained comps {0,3,5} on entry; the lateral block
// {1,2,4} starts from the committed guess and is OVERWRITTEN with the converged condensed strains. Writes
// the NOMINAL stress (all 6 comps; the wrapper reads {0,3,5}), the implicit effective stress, the new
// State, and (doTangent) the condensed consistent tangent on the retained comps:
//   dsig_R/deps_R = Cdam_RR - Cdam_RL * (Ceff_LL + Hoop)^-1 * Ceff_LR
// (Cdam = the P2 damaged tangent for the output rows; the constraint block uses the EFFECTIVE tangent
// Ceff_LL + the hoop stiffness on the lateral-normal diagonal). Reduces to a plain static condensation
// where omega->0 (Cdam->Ceff). Returns 0 converged / 2 the inviscid effective return did not converge.
inline int driveConfinedFiber(const Params& mp, double strain[6], const State& inRaw, State& out,
                              double sigma[6], double sigEffImpl[6], double Dtan6[6][6],
                              bool doTangent, double hoopK, double hoopFy, double dt = 0.0)
{
    const int L[3] = {1, 2, 4};                       // condensed lateral block (eps_11, eps_22, gamma_12)
    const double tol = 1.0e-10 * (mp.fc + 1.0);
    double sig_eff[6] = {0,0,0,0,0,0}, kp_new = inRaw.kp, Ceff[6][6];
    int status = 0;
    // DEAD POINTS (review M2; the same treatment as returnMap, via the shared effectiveReturn): a CRUSHED fibre point
    // (omega_c >= omegaDead) is frozen -- elastic on the fixed plastic strain, both damages at the floor; a CRACKED one
    // (omega_t >= omegaDead) carries its tensile effective stress elastically and returns on the compressive remainder.
    const bool crushed = (inRaw.wc >= mp.omegaDead);
    State cutBuf; const State* inP = &inRaw; bool cutT = false;
    auto effective = [&](const double eps[6], int* subInfo) {
        if (crushed) {
            double deps[6]; for (int i = 0; i < 6; ++i) deps[i] = eps[i] - inRaw.eps[i];
            elasticPredTensor(inRaw.sigEff, deps, mp, sig_eff); kp_new = inRaw.kp; elasticC(mp, Ceff);
            if (subInfo) *subInfo = 0;
            return 0;
        }
        return effectiveReturn(mp, eps, inRaw, true, 0.0, true, cutBuf, inP, cutT, sig_eff, kp_new, Ceff, subInfo);
    };
    for (int it = 0; it < 80; ++it) {                 // nested lateral Newton vs the hoop residual
        status = effective(strain, nullptr);
        double r[3] = { sig_eff[1] + hoopStress(strain[1], hoopK, hoopFy),
                        sig_eff[2] + hoopStress(strain[2], hoopK, hoopFy),
                        sig_eff[4] };
        if (std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) < tol) break;
        double J[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) J[i][j] = Ceff[L[i]][L[j]];
        J[0][0] += hoopStiffness(strain[1], hoopK, hoopFy);
        J[1][1] += hoopStiffness(strain[2], hoopK, hoopFy);
        double dx[3];
        if (!solve3(J, r, dx)) break;
        for (int j = 0; j < 3; ++j) strain[L[j]] -= dx[j];          // J dx = r  =>  step -dx
    }
    // final sync (mirror the oracle's post-loop re-evaluation): recompute the effective return at the
    // converged lateral strain so sig_eff/Ceff/kp_new always correspond to the committed strain (guards
    // the rare 80-iter-exhausted case where the last lateral update post-dates the last effective return),
    // and flag a genuinely unmet lateral balance as non-converged (status 2 => the caller cuts the step).
    {
        status = effective(strain, &out.subInfo);
        const double rf[3] = { sig_eff[1] + hoopStress(strain[1], hoopK, hoopFy),
                               sig_eff[2] + hoopStress(strain[2], hoopK, hoopFy),
                               sig_eff[4] };
        if (std::sqrt(rf[0]*rf[0] + rf[1]*rf[1] + rf[2]*rf[2]) >= 1.0e-6 * (mp.fc + 1.0)) status = 2;
    }
    const State& in = *inP;
    double Cdam[6][6];
    if (crushed) {
        // frozen point (same rule as returnMap, see the note there): histories and kappa_p as committed, both damages at the
        // floor, nominal = (1-OMEGA_MAX) sig_eff
        out = inRaw;
        const double k = 1.0 - OMEGA_MAX;
        for (int i = 0; i < 6; ++i) { out.eps[i] = strain[i]; out.sigEff[i] = sig_eff[i]; sigEffImpl[i] = sig_eff[i];
                                       out.sig[i] = k * sig_eff[i]; sigma[i] = out.sig[i]; out.depl[i] = 0.0; }
        out.wt = OMEGA_MAX; out.wc = OMEGA_MAX; out.dwt = 0.0; out.dwc = 0.0;
        out.dt_n = (dt > 0.0) ? dt : inRaw.dt_n;
        if (doTangent) { double C0[6][6]; elasticC(mp, C0); const double kT = k > OMEGA_TAN_FLOOR ? k : OMEGA_TAN_FLOOR;
                          for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Cdam[i][j] = kT * C0[i][j]; }
    } else {
        // commit the converged effective state + the IMPLICIT P2 dual-damage nominal stress (mirror returnMap)
        for (int i = 0; i < 6; ++i) { out.eps[i] = strain[i]; out.sigEff[i] = sig_eff[i]; sigEffImpl[i] = sig_eff[i]; }
        out.kp = kp_new;
        double wt_impl = 0.0, wc_impl = 0.0;
        damagedUpdate(mp, in, sig_eff, kp_new, strain, out, &wt_impl, &wc_impl);
        for (int i = 0; i < 6; ++i) sigma[i] = out.sig[i];
        out.wt = wt_impl; out.wc = wc_impl; out.dwt = wt_impl - in.wt; out.dwc = wc_impl - in.wc;
        { double epl[6], epl_n[6];
          plasticStrain6(out.sigEff, out.eps, mp, epl);
          plasticStrain6(in.sigEff,  in.eps,  mp, epl_n);
          for (int i = 0; i < 6; ++i) out.depl[i] = epl[i] - epl_n[i]; }
        out.dt_n = (dt > 0.0) ? dt : in.dt_n;
        if (doTangent && status == 0) damagedTangent(mp, in, sig_eff, strain, kp_new, Ceff, Cdam);
    }

    if (doTangent) {
        if (status != 0) { elasticC(mp, Dtan6); return status; }   // safe fallback; caller cuts the step
        double Kll[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Kll[i][j] = Ceff[L[i]][L[j]];
        Kll[0][0] += hoopStiffness(strain[1], hoopK, hoopFy);
        Kll[1][1] += hoopStiffness(strain[2], hoopK, hoopFy);
        double Tmat[3][6];                                          // Tmat = Kll^-1 * Ceff[L][:]
        for (int k = 0; k < 6; ++k) {
            double rhs[3] = { Ceff[L[0]][k], Ceff[L[1]][k], Ceff[L[2]][k] }, tcol[3];
            if (!solve3(Kll, rhs, tcol)) { tcol[0] = tcol[1] = tcol[2] = 0.0; }
            for (int i = 0; i < 3; ++i) Tmat[i][k] = tcol[i];
        }
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {   // Dtan = Cdam - Cdam[:,L] * Tmat
            double v = Cdam[i][j];
            for (int a = 0; a < 3; ++a) v -= Cdam[i][L[a]] * Tmat[a][j];
            Dtan6[i][j] = v;                                        // only retained rows/cols {0,3,5} are read
        }
    }
    return status;
}

// ===========================================================================
// Gc = PHYSICAL compressive fracture energy (WP concrete3d-oracle-diagnosis; mirrors the oracle
// compression_energy_density / calibrate_eps_fc_table / eps_fc_from_gc). The compressive law is in STRAIN
// form (eps_fc) and its plastic driver is scaled by beta_c/x_s (Eq.48/50), so the legacy eps_fc = Gc/(fc lch)
// dissipates ~an order of magnitude more than Gc per unit area. The uniaxial stress-strain response depends
// on eps_fc but NOT on lch, so g(eps_fc) = int_peak sigma d(eps - sigma/E) (post-peak energy per unit volume)
// is a material function and Gc = lch g(eps_fc): tabulate g once (at material construction, free uniaxial
// stress via driveConfinedFiber hoopK=0, the kernel's own single-step contract), invert per lch.
// ===========================================================================
static const int    GC_TABLE_N  = 8;
static const double GC_TABLE_LO = 0.01, GC_TABLE_HI = 3.0;   // eps_fc range in units of fc/E

inline double compressionEnergyDensity(const Params& mp0, double epsFc, int maxSteps = 20000, double* peakOut = nullptr)
{
    Params mp = mp0; mp.epsFc = epsFc; mp.implex = false; mp.eta = 0.0;
    const double E = mp.E, fc = mp.fc, de = fc / (20.0 * E);        // FIXED step (see the oracle docstring)
    double e11 = 0.0, sPrev = 0.0, eiPrev = 0.0, peak = 0.0, g = 0.0;
    bool post = false;
    State in;
    for (int n = 1; n <= maxSteps; ++n) {
        e11 -= de;
        double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i];
        strain[0] = e11; strain[3] = strain[5] = 0.0;
        State out; double sig[6], sigEff[6], Dt[6][6];
        driveConfinedFiber(mp, strain, in, out, sig, sigEff, Dt, false, 0.0, 1.0e30, 0.0);
        in = out;
        const double s = -sig[0], e = -e11, ei = e - s / E;
        if (s > peak && !post) peak = s;
        else if (s < peak) post = true;
        if (post) {
            g += 0.5 * (s + sPrev) * (ei - eiPrev);
            if (s < 0.01 * peak || n == maxSteps) {
                const double ds = sPrev - s;
                if (ds > 0.0) g += s * s * (ei - eiPrev) / ds;       // exponential tail beyond the stop
                break;
            }
        }
        sPrev = s; eiPrev = ei;
    }
    if (peakOut) *peakOut = peak;
    return g;
}

inline void calibrateEpsFcTable(const Params& mp, double efc[GC_TABLE_N], double g[GC_TABLE_N])
{
    const double base = mp.fc / mp.E, lo = std::log(GC_TABLE_LO), hi = std::log(GC_TABLE_HI);
    for (int k = 0; k < GC_TABLE_N; ++k) {
        efc[k] = base * std::exp(lo + (hi - lo) * k / (GC_TABLE_N - 1));
        g[k] = compressionEnergyDensity(mp, efc[k]);
    }
}

// Invert the monotone table for g = Gc/lch (log-log linear). status 0 inside, -1 below (Gc too small for this
// lch — the brittle/snap-back limit; clamped to the smallest eps_fc), +1 above (clamped to the largest).
inline double epsFcFromGc(const double efc[GC_TABLE_N], const double g[GC_TABLE_N], double Gc, double lch,
                          int* status = nullptr)
{
    const double tgt = Gc / lch;
    if (status) *status = 0;
    if (!(tgt > g[0])) { if (status) *status = -1; return efc[0]; }
    if (!(tgt < g[GC_TABLE_N - 1])) { if (status) *status = 1; return efc[GC_TABLE_N - 1]; }
    for (int k = 1; k < GC_TABLE_N; ++k) {
        if (tgt <= g[k]) {
            const double t = (std::log(tgt) - std::log(g[k - 1])) / (std::log(g[k]) - std::log(g[k - 1]));
            return std::exp(std::log(efc[k - 1]) + t * (std::log(efc[k]) - std::log(efc[k - 1])));
        }
    }
    if (status) *status = 1;
    return efc[GC_TABLE_N - 1];
}

} // namespace Concrete3D
} // namespace Ladruno

#endif // LadrunoConcrete3DKernel_h
