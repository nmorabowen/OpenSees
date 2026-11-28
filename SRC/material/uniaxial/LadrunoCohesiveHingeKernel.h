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
// Created: 11/2026  (WP-267)
//
// LadrunoCohesiveHingeKernel — the PURE numerical core of LadrunoCohesiveHinge
// (calibration, P-M multiplier, envelope, trial state incl. IMPL-EX), extracted
// from the class so it is written ONCE for any scalar:
//
//   * T = P = double            the production path, which performs exactly the
//                               operations of the former in-class code
//                               (bit-identity gate: tests/ladruno_cohesivehinge_dual_check.cpp);
//   * T = P = Dual<1>           a material constant seeded => its parameter
//                               sensitivity (DDM getStressSensitivity /
//                               commitSensitivity);
//   * T = Dual<1>, P = double   the strain (rotation jump) seeded => the tangent.
//
// Plain C++17 + <math.h>, no OpenSees dependency. The kernel never prints: the
// calibration returns a "clamped" flag and the class prints the warning.
// The law itself is documented in LadrunoCohesiveHinge.h.

#ifndef LadrunoCohesiveHingeKernel_h
#define LadrunoCohesiveHingeKernel_h

#include <math.h>
#include <type_traits>
#include <LadrunoDual.h>

namespace ladruno_cohesive_hinge {

using ladruno_dual::value;

enum { EXP = 0, LINEAR = 1 };   // matches LadrunoCohesiveHinge::SofteningShape

// material constants (P = double in production; a Dual to seed one of them)
template <class P>
struct ParamsT {
  P Mc, Gf, Kpen;          // Kpen as stored by the class (computeDerived output)
  bool userKpen;
  int shape;
  double penaltyRatio;
  P P0, alphaPM, P0t, Pbal, muBal;
  double P_axial;          // supplied by the parent beam; not differentiated
  P kappa0, aExp, kappaf;  // derived by computeDerived()
};

// committed history (the trial kappaMax is an output of trialState)
template <class T>
struct StateT {
  T kappaMax, kappaMaxOld;
};

// Calibration (verbatim LadrunoCohesiveHinge::computeDerived, opserr removed).
// Returns 1 if a penalty at or below the floor was clamped; KpenIn / KpenFloorOut /
// KpenSafeOut (optional) receive the values the class prints in its warning.
template <class P>
inline int computeDerived(ParamsT<P>& p, double* KpenIn = 0, double* KpenFloorOut = 0,
                          double* KpenSafeOut = 0)
{
  if (p.Mc <= 0.0 || p.Gf <= 0.0) { p.kappa0 = 0.0; p.aExp = 0.0; p.kappaf = 0.0; return 0; }

  const P KpenFloor = p.Mc*p.Mc / (2.0*p.Gf);   // pre-peak energy == Gf at the floor

  // default penalty (!userKpen or Kpen<=0): penaltyRatio multiples above the floor
  if (!p.userKpen || p.Kpen <= 0.0) {
    double r = (p.penaltyRatio > 1.0) ? p.penaltyRatio : 1000.0;
    p.Kpen = r * KpenFloor;
  }

  // guarded floor: Kpen must exceed the floor or the softening branch is non-physical
  int clamped = 0;
  if (p.Kpen <= KpenFloor) {
    double rSafe = (p.penaltyRatio > 1.0) ? p.penaltyRatio : 1000.0;
    P KpenSafe = rSafe * KpenFloor;
    if (KpenIn) *KpenIn = value(p.Kpen);
    if (KpenFloorOut) *KpenFloorOut = value(KpenFloor);
    if (KpenSafeOut) *KpenSafeOut = value(KpenSafe);
    p.Kpen = KpenSafe;
    clamped = 1;
  }

  p.kappa0 = p.Mc / p.Kpen;                          // jump at activation
  P Esoft = p.Gf - p.Mc*p.Mc/(2.0*p.Kpen);           // energy left for the softening branch (>0 by floor)
  p.aExp   = p.Mc / Esoft;                           // exp: area Mc/aExp == Esoft
  p.kappaf = 2.0*Esoft / p.Mc;                       // lin: area Mc*kappaf/2 == Esoft
  return clamped;
}

// P-M multiplier f(P_axial) (verbatim LadrunoCohesiveHinge::getPMMultiplier).
// The branch and clamp tests read values; the clamped plateaus have zero slope.
template <class P>
inline P pmMultiplier(const ParamsT<P>& p)
{
  const double P_axial = p.P_axial;
  if (p.P0 <= 0.0) return P(1.0);
  // 4-parameter C^1 piecewise-parabolic RC balance-point law (-pmBal)
  if (p.Pbal > 0.0 && p.P0 > p.Pbal && p.P0t > 0.0) {
    P f = 0.0;
    if (P_axial >= 0.0 && P_axial <= p.P0t) {
      f = 1.0 - P_axial / p.P0t;
    } else if (P_axial < 0.0 && P_axial >= -p.Pbal) {
      P u = (P_axial + p.Pbal) / p.Pbal;
      f = 1.0 + (p.muBal - 1.0) * (1.0 - u * u);
    } else if (P_axial < -p.Pbal && P_axial >= -p.P0) {
      P denom = p.P0 - p.Pbal;
      P u = (P_axial + p.Pbal) / denom;
      f = p.muBal * (1.0 - u * u);
    }
    if (f < 0.05) f = 0.05;
    else if (f > 3.0) f = 3.0;
    return f;
  }
  // OpenSees beam basic axial force q(0): positive in tension, negative in compression.
  // Parabolic P-M scaling: f(n) = clamp(1 - alphaPM * n - (1 - alphaPM) * n^2, 0.05, 2.5)
  P n = P_axial / p.P0;
  P f = 1.0 - p.alphaPM * n - (1.0 - p.alphaPM) * n * n;
  if (f < 0.05) f = 0.05;
  else if (f > 2.5) f = 2.5;
  return f;
}

// initial tangent (LadrunoCohesiveHinge::getInitialTangent)
template <class P>
inline P initialTangent(const ParamsT<P>& p)
{
  if (p.P0 > 0.0) return P(p.Kpen * pmMultiplier(p));
  return p.Kpen;
}

// envelope moment Menv(km) and consistent slope dMenv/dkm at km >= 0
// (verbatim LadrunoCohesiveHinge::envelope)
template <class T, class P>
inline void envelope(const ParamsT<P>& p, const T& km, T& Menv, T& dMenv)
{
  if (km <= p.kappa0) {                      // pre-peak penalty (hinge closed)
    Menv  = p.Kpen * km;
    dMenv = p.Kpen;
  }
  else if (p.shape == LINEAR) {
    T d = km - p.kappa0;
    if (d >= p.kappaf) {                     // fully broken: zero residual moment
      Menv = 0.0; dMenv = 0.0;
    } else {
      Menv  = p.Mc * (1.0 - d/p.kappaf);
      dMenv = -p.Mc/p.kappaf;
    }
  } else {                                   // EXP
    T e = exp(-p.aExp * (km - p.kappa0));
    Menv  =  p.Mc * e;
    dMenv = -p.aExp * p.Mc * e;
  }
  if (p.P0 > 0.0) {
    P f = pmMultiplier(p);
    Menv  *= f;
    dMenv *= f;
  }
}

// Trial state (verbatim LadrunoCohesiveHinge::setTrialStrain minus the work
// integral, which the class keeps). implex / timeFactor: the IMPL-EX switch and
// the frozen step ratio the class computes from LadrunoImplexClock.
template <class T, class P>
inline void trialState(const ParamsT<P>& p, bool implex, double timeFactor,
                       const T& strain, const StateT<T>& C,
                       T& TkappaMax, T& stress, T& tangent)
{
  T akappa = fabs(strain);

  // irreversible damage driver: max |jump| ever seen
  bool loading = (akappa >= C.kappaMax);
  TkappaMax = loading ? akappa : C.kappaMax;

  if (implex) {
    // Oliver et al. (2008) IMPL-EX step extrapolation of the damage driver (see
    // LadrunoCohesiveHinge.cpp): the secant of the extrapolated committed history.
    double time_factor = timeFactor;
    T dKappaMax_n = C.kappaMax - C.kappaMaxOld;
    T kappaMax_tilde = C.kappaMax + time_factor * dKappaMax_n;
    if (kappaMax_tilde < C.kappaMax) kappaMax_tilde = C.kappaMax;

    T Menv_tilde = 0.0, dMenv_tilde = 0.0;
    envelope(p, kappaMax_tilde, Menv_tilde, dMenv_tilde);
    P fpm = (p.P0 > 0.0) ? pmMultiplier(p) : P(1.0);
    P KpenEff = p.Kpen * fpm;
    T Ksec_tilde;
    if (kappaMax_tilde > 0.0) Ksec_tilde = Menv_tilde / kappaMax_tilde;
    else                      Ksec_tilde = KpenEff;
    P Kfloor = 1.0e-6 * KpenEff;
    if (Ksec_tilde < Kfloor) Ksec_tilde = Kfloor;

    stress  = Ksec_tilde * strain;
    tangent = Ksec_tilde;
    return;
  }

  T Menv, dMenv;
  envelope(p, TkappaMax, Menv, dMenv);

  // secant-to-origin (isotropic damage); Ksec = Menv(kappaMax)/kappaMax
  P KpenEff = (p.P0 > 0.0) ? P(p.Kpen * pmMultiplier(p)) : p.Kpen;
  T Ksec;
  if (TkappaMax > 0.0) Ksec = Menv / TkappaMax;
  else                 Ksec = KpenEff;

  stress = Ksec * strain;                    // sign carried by `strain`

  if (loading && akappa > 0.0)
    tangent = dMenv;                         // consistent envelope slope (may be < 0)
  else
    tangent = Ksec;                          // unload/reload secant (>= 0)
}

} // namespace ladruno_cohesive_hinge

#endif
