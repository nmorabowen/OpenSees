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
// LadrunoBondSlipKernel — the PURE numerical core of LadrunoBondSlip (derived
// constants, MC2010 backbone, friction-slider trial state), extracted from the
// class so it is written ONCE for any scalar:
//
//   * T = P = double            the production path, which performs exactly the
//                               operations of the former in-class code
//                               (bit-identity gate: tests/ladruno_bondslip_dual_check.cpp);
//   * T = P = Dual<1>           a backbone constant seeded => its parameter
//                               sensitivity (DDM getStressSensitivity /
//                               commitSensitivity);
//   * T = Dual<1>, P = double   the slip seeded => the tangent.
//
// Plain C++17 + <math.h>, no OpenSees dependency. The law is documented in
// LadrunoBondSlip.h.

#ifndef LadrunoBondSlipKernel_h
#define LadrunoBondSlipKernel_h

#include <math.h>
#include <type_traits>
#include <LadrunoDual.h>

namespace ladruno_bond_slip {

using ladruno_dual::value;

// x^y for the backbone power law. LadrunoDual.h has pow(Dual, double) only; the
// exponent alpha is a parameter, so a seeded alpha needs pow(Dual, Dual). The
// value is ::pow(x, y) (bit-identical), d = y x^(y-1) dx + x^y ln(x) dy.
// Only called at x > 0 (the power-law branch is s0 < a <= s1, and s0/s1 > 0).
inline double powT(double x, double y) { return ::pow(x, y); }
template <int N>
inline ladruno_dual::Dual<N> powT(const ladruno_dual::Dual<N>& x, double y)
{
  return ladruno_dual::pow(x, y);
}
template <int N>
inline ladruno_dual::Dual<N> powT(const ladruno_dual::Dual<N>& x, const ladruno_dual::Dual<N>& y)
{
  ladruno_dual::Dual<N> r(::pow(x.v, y.v));
  const double kx = y.v * ::pow(x.v, y.v - 1.0);
  const double ky = r.v * ::log(x.v);
  for (int i = 0; i < N; i++) r.d[i] = kx*x.d[i] + ky*y.d[i];
  return r;
}

// backbone constants (P = double in production; a Dual to seed one of them)
template <class P>
struct ParamsT {
  P tau_max, s1, s2, s3, tau_f, alpha, Gf, s0;
  P k0;                       // derived by setDerived()
};

// committed state
template <class T>
struct StateT {
  T slip, stress, slipMaxAbs;
};

// branch of the last trial state (for the tests; the class ignores it)
enum { BR_ENVELOPE = 0, BR_ENVELOPE_FLOORED = 1, BR_UNLOAD = 2, BR_SLIDE = 3 };

// s0 default, k0, and the Gf override of s3 (verbatim LadrunoBondSlip::setDerived)
template <class P>
inline void setDerived(ParamsT<P>& p)
{
  if (p.s1 <= 0.0) p.s1 = 1.0e-12;
  if (p.s0 <= 0.0) p.s0 = 0.1 * p.s1;
  if (p.s0 > p.s1) p.s0 = p.s1;
  // initial stiffness = secant to the power-law value at s0 (D4.1 linear segment)
  P tau_s0 = p.tau_max * powT(p.s0 / p.s1, p.alpha);
  if (p.s0 > 0.0) p.k0 = tau_s0 / p.s0;
  else            p.k0 = p.tau_max / p.s1;
  // D4.3: Gf overrides s3 so the softening triangle dissipates Gf per unit area
  if (p.Gf > 0.0 && p.tau_max > p.tau_f) {
    p.s3 = p.s2 + 2.0 * p.Gf / (p.tau_max - p.tau_f);
  }
  if (p.s2 < p.s1) p.s2 = p.s1;
  if (p.s3 < p.s2) p.s3 = p.s2;
}

// backbone (a = |slip| >= 0) (verbatim LadrunoBondSlip::envelope)
template <class T, class P>
inline T envelope(const ParamsT<P>& p, const T& a)
{
  if (a <= p.s0)  return p.k0 * a;
  if (a <= p.s1)  return p.tau_max * powT(a / p.s1, p.alpha);
  if (a <= p.s2)  return T(p.tau_max);
  if (a <= p.s3)  return p.tau_max - (p.tau_max - p.tau_f) * (a - p.s2) / (p.s3 - p.s2);
  return T(p.tau_f);
}

// (verbatim LadrunoBondSlip::envelopeTangent)
template <class T, class P>
inline T envelopeTangent(const ParamsT<P>& p, const T& a)
{
  if (a <= p.s0)  return T(p.k0);
  if (a <= p.s1)  return p.tau_max * p.alpha / p.s1 * powT(a / p.s1, p.alpha - 1.0);
  if (a <= p.s2)  return T(0.0);
  if (a <= p.s3)  return T(-(p.tau_max - p.tau_f) / (p.s3 - p.s2));
  return T(0.0);
}

// Trial state (verbatim LadrunoBondSlip::setTrialStrain minus the work
// integral, which the class keeps). Returns the branch (BR_*).
template <class T, class P>
inline int trialState(const ParamsT<P>& p, const T& slip, const StateT<T>& C,
                      T& TslipMaxAbs, T& stress, T& tangent)
{
  T a = fabs(slip);
  double sgn = (slip >= 0.0) ? 1.0 : -1.0;

  int branch;
  if (a >= C.slipMaxAbs) {             // virgin loading on the envelope
    TslipMaxAbs = a;
    stress  = sgn * envelope(p, a);
    tangent = envelopeTangent(p, a);   // (>=0 ascending, <0 softening, 0 plateau)
    branch = BR_ENVELOPE;
    if (tangent == 0.0) {              // avoid a hard 0 pivot
      if (a <= p.s2) tangent = p.k0;
      else           tangent = -1.0e-12 * p.k0;
      branch = BR_ENVELOPE_FLOORED;
    }
  }
  else {                               // inside the peak: elastic, then friction
    TslipMaxAbs = C.slipMaxAbs;
    // Friction slider: k0 from the committed point, capped at the envelope stress
    // of the peak slip in either direction.
    const T tauCap = envelope(p, C.slipMaxAbs);
    stress = C.stress + p.k0 * (slip - C.slip);
    branch = BR_UNLOAD;
    if (stress > tauCap)       { stress = tauCap;  branch = BR_SLIDE; }
    else if (stress < -tauCap) { stress = -tauCap; branch = BR_SLIDE; }
    tangent = p.k0;                    // also while sliding: keeps the pivot finite
  }
  return branch;
}

} // namespace ladruno_bond_slip

#endif
