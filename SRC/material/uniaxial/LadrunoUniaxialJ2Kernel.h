/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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
// Created: 11/2026  (WP-262)
//
// LadrunoUniaxialJ2Kernel — the PURE numerical core of LadrunoUniaxialJ2's
// return map (combined Voce+linear isotropic + Chaboche AF kinematic, optional
// Lemaitre damage, IMPL-EX), extracted from LadrunoUniaxialJ2::setTrialStrain
// so it can be written ONCE for any scalar:
//
//   * T = P = double            the production path (setTrialStrain), which
//                               performs exactly the operations of the former
//                               in-class code (bit-identity gate:
//                               tests/ladruno_uniaxialj2_dual_check.cpp);
//   * T = P = Dual<1>           a material constant seeded => its parameter
//                               sensitivity (DDM getStressSensitivity /
//                               commitSensitivity);
//   * T = Dual<1>, P = double   the strain seeded => the consistent tangent.
//
// Plain C++17 + <math.h>, no OpenSees dependency. The kernel never prints: it
// returns a status and the class prints the warnings it always printed.
//
// Model: sig_tr = E (eps - eps^p_n), eta = sig_tr - sum_k X_k,n, scalar Newton on
// the plastic multiplier dp of
//   R(dp) = |eta| - E dp - sum_k dp (C_k - g_k s X_k,n)/(1 + g_k dp) - sig_y(p_n + dp)
// consistent tangent E h/(E + h); Lemaitre damage as the 1D reduction of the 3D
// LadrunoJ2 layer (R_v == 1 at uniaxial stress). See LadrunoUniaxialJ2.cpp.

#ifndef LadrunoUniaxialJ2Kernel_h
#define LadrunoUniaxialJ2Kernel_h

#include <math.h>
#include <type_traits>
#include <LadrunoDual.h>
#include <LadrunoHardening.h>
#include <LadrunoDamage.h>

namespace ladruno_uniaxial_j2 {

using ladruno_dual::value;

const int MAXBACK = 8;   // matches LadrunoUniaxialJ2::MAXBACK

enum { STATUS_OK = 0, STATUS_SINGULAR = 1, STATUS_NO_CONVERGE = 2 };

// material constants (P = double in production; a Dual to seed one of them)
template <class P>
struct ParamsT {
  P E, sig0, Qinf, bIso, Hiso;
  int nBack;
  P C[MAXBACK];      // AF kinematic moduli C_k
  P gam[MAXBACK];    // AF recall constants gamma_k
  // Lemaitre damage (not differentiated: plain doubles)
  bool   dmgOn;
  double dmgR, dmgS, dmgPD, dmgDc;
  bool   dmgImplex;
  double lchScale;   // = lch / lch_ref (1 = local); read from the element by the class
};

// history (committed C or trial T)
template <class T>
struct StateT {
  T plasticStrain, ebarp, dGamma, wp, D;
  T X[MAXBACK];
};

// strain: trial total strain. C: committed state; CdD: committed implicit damage
// increment (IMPL-EX extrapolation). Tr: trial state (out). stress / tangent (out).
// outIter (optional): local Newton iteration count; outResidual: final |R|.
template <class T, class P>
inline int returnMap(const ParamsT<P>& p, const T& strain,
                     const StateT<T>& C, const T& CdD,
                     StateT<T>& Tr, T& stress, T& tangent,
                     double* outResidual = 0)
{
  const int nBack = p.nBack;
  const P& E = p.E;
  auto yieldStress = [&](const T& pb) { return T(Ladruno::yieldStressVoceLinear(pb, p.sig0, p.Qinf, p.bIso, p.Hiso)); };
  auto yieldSlope  = [&](const T& pb) { return T(Ladruno::yieldSlopeVoceLinear(pb, p.sig0, p.Qinf, p.bIso, p.Hiso)); };

  // trial (elastic predictor): freeze plastic state
  T sigTr = E*(strain - C.plasticStrain);
  T Xtot = 0.0;
  for (int k = 0; k < nBack; k++) Xtot += C.X[k];
  T eta = sigTr - Xtot;                      // relative stress
  T absEta = fabs(eta);
  T syn = yieldStress(C.ebarp);              // sig_y at committed pbar
  T f = absEta - syn;

  // unit-aware tolerance (never collapses to a fixed absolute when sig0 == 0,
  // e.g. kinematic-only); scaled to the operative stress magnitude. Control
  // flow only: it reads the VALUE of a Dual.
  double stressScale = fmax(value(absEta), value(syn));
  double tolY = 1.0e-12 * fmax(stressScale, 1.0e-300);

  if (value(f) <= tolY) {
    // ---- elastic step: state frozen ----
    Tr.plasticStrain = C.plasticStrain;
    Tr.ebarp = C.ebarp;
    Tr.dGamma = 0.0;
    Tr.wp = C.wp;
    for (int k = 0; k < nBack; k++) Tr.X[k] = C.X[k];
    stress = sigTr;
    tangent = E;
    Tr.D = C.D;                              // no new damage on an elastic step
    if (p.dmgOn) {                           // prior damage still degrades sigma/tangent
      T dScale = p.dmgImplex ? T(C.D + CdD) : C.D;   // IMPL-EX: frozen extrapolated D~
      if (value(dScale) > p.dmgDc) dScale = p.dmgDc;
      if (value(dScale) < 0.0)     dScale = 0.0;
      if (value(dScale) > 0.0) { T omd = 1.0 - dScale; stress *= omd; tangent *= omd; }
    }
    return STATUS_OK;
  }

  // ---- plastic corrector ----
  int s = (value(eta) >= 0.0) ? 1 : -1;

  // initial plastic modulus (for the Newton seed): iso slope + sum C_k
  T h0 = yieldSlope(C.ebarp);
  for (int k = 0; k < nBack; k++) h0 += p.C[k];
  T dp = f / (E + h0);                       // perfectly-plastic-style seed
  if (value(dp) < 0.0) dp = 0.0;

  // scalar Newton:  R(dp) = |eta| - E dp
  //   - sum_k dp (C_k - gamma_k s X_k,n)/(1 + gamma_k dp) - sig_y(pbar_n + dp) = 0
  //   a_k = C_k - gamma_k s X_k,n   (constant in dp)
  int status = STATUS_OK;
  const int maxIter = 50;
  int iter = 0;
  T R = 1.0, dR = 0.0;
  for (; iter < maxIter; iter++) {
    T R0 = absEta - E*dp - yieldStress(C.ebarp + dp);
    T dRl = -E - yieldSlope(C.ebarp + dp);
    for (int k = 0; k < nBack; k++) {
      T ak = p.C[k] - p.gam[k]*s*C.X[k];
      T den = 1.0 + p.gam[k]*dp;
      R0 -= dp*ak/den;
      dRl -= ak/(den*den);
    }
    R = R0;
    dR = dRl;
    if (fabs(value(R)) <= tolY) break;
    if (value(dR) == 0.0) { status = STATUS_SINGULAR; break; }
    dp -= R/dR;
    if (value(dp) < 0.0) dp = 0.0;
  }
  if (iter >= maxIter) status = STATUS_NO_CONVERGE;
  if (outResidual) *outResidual = fabs(value(R));

  // Dual T only: one value-zeroed Newton step at the converged dp so its
  // derivative is the implicit-function-theorem one (see LadrunoJ2Kernel.h and
  // quirk WP-262-dual-derivative-lags-local-newton). Compiled out for double.
  if constexpr (ladruno_dual::is_dual<T>::value) {
    if (status == STATUS_OK && iter < maxIter && value(dR) != 0.0) {
      T corr = R/dR;
      ladruno_dual::zero_value(corr);
      dp -= corr;
    }
  }

  // ---- state update ----
  T pbar = C.ebarp + dp;
  Tr.ebarp = pbar;
  Tr.dGamma = dp;
  Tr.plasticStrain = C.plasticStrain + dp*s;
  for (int k = 0; k < nBack; k++)
    Tr.X[k] = (C.X[k] + p.C[k]*dp*s) / (1.0 + p.gam[k]*dp);

  // ---- stress ----
  stress = sigTr - E*dp*s;

  // ---- plastic work increment (fracture/LCF hook): dWp = sigma * deps^p ----
  Tr.wp = C.wp + stress*(dp*s);

  // ---- consistent tangent: E_alg = E h/(E+h),  h = sig_y' + sum a_k/(1+g_k dp)^2
  T h = yieldSlope(pbar);
  for (int k = 0; k < nBack; k++) {
    T ak = p.C[k] - p.gam[k]*s*C.X[k];
    T den = 1.0 + p.gam[k]*dp;
    h += ak/(den*den);
  }
  tangent = E*h/(E + h);

  // ---- Lemaitre coupled ductile damage (uniaxial: triaxiality 1/3 => R_v == 1,
  //      so Y = sigma~^2/(2E)). The effective return above is D-independent; D is
  //      recovered from the converged effective stress and the (dp,D) coupling
  //      enters the CONSISTENT TANGENT: E_nom = (1-D)E_alg - sigma~ dD/deps.
  Tr.D = C.D;
  if (p.dmgOn) {
    const T sigEff  = stress;                // effective stress (pre-degradation)
    const T EalgEff = tangent;               // effective consistent tangent
    Ladruno::DamageParams dpar;
    dpar.on = p.dmgOn; dpar.r = p.dmgR; dpar.s = p.dmgS; dpar.pD = p.dmgPD; dpar.Dc = p.dmgDc;
    const T Y = (value(E) > 0.0) ? T(sigEff*sigEff/(2.0*E)) : T(0.0);
    const T dDraw = p.lchScale * Ladruno::damageIncrement(Y, C.ebarp, pbar, dpar);
    T Dnew = C.D + dDraw;
    bool capped = false;
    if (value(Dnew) >= p.dmgDc) { Dnew = p.dmgDc; capped = true; }
    if (value(Dnew) < 0.0) Dnew = 0.0;
    Tr.D = Dnew;                             // implicit damage (committed for extrapolation)

    // IMPL-EX: degrade by the FROZEN extrapolated D~ = CD + CdD (constant tangent
    // omd*E_alg, no dD/deps softening term); else fully implicit.
    T dScale = p.dmgImplex ? T(C.D + CdD) : Dnew;
    if (value(dScale) > p.dmgDc) dScale = p.dmgDc;
    if (value(dScale) < 0.0)     dScale = 0.0;
    const T omd = 1.0 - dScale;

    T dDdeps = 0.0;
    if (!p.dmgImplex && !capped && value(dDraw) > 0.0) {
      // d(pbar)/deps = s E/(E+h) ; dY/deps = (sigma~/E) E_alg ; partials of dD
      T ddD_dY, ddD_dp;
      Ladruno::damageIncrementPartials(Y, C.ebarp, pbar, dpar, ddD_dY, ddD_dp);
      const T dpbar_deps = s*E/(E + h);
      const T dY_deps    = (sigEff/E)*EalgEff;
      dDdeps = p.lchScale * (ddD_dp*dpbar_deps + ddD_dY*dY_deps);
    }
    stress  = omd*sigEff;
    tangent = omd*EalgEff - sigEff*dDdeps;
  }

  return status;
}

} // namespace ladruno_uniaxial_j2

#endif
