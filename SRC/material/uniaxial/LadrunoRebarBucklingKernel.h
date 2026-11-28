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
// Created: 11/2026  (WP-264)
//
// LadrunoRebarBucklingKernel — the closed-form Gomes-Appleton buckled-average factor of LadrunoRebarBuckling
// (-model ga), extracted from LadrunoRebarBuckling::gaFactor so it is written ONCE for any scalar:
//
//   * T = double                 the production stress (LadrunoRebarBuckling::gaFactor), which performs exactly the
//                                operations of the former in-class code (bit-identity gate:
//                                tests/ladruno_rebarbuckling_dual_check.cpp);
//   * T = ladruno_dual::Dual<1>  the strain seeded => dfactor/deps for the consistent tangent, exact. It replaces the
//                                central FD (h = 1e-8) the tangent used before WP-264.
//
// Ported from ReinforcingSteel::Buckled_stress_Gomes(ess, fss): the upstream stress function hardcodes its own
// beta = 1, gama = 0.1, Dft = 0.25, so only lsr, the reduction r and the anchor drive it.
//   fs_buck = sqrt(32 / (e_cross - eps)) / (3 pi lsr)
//   beta    = 1 - gama (Dft - |fs_buck - 1|) / Dft   if |fs_buck - 1| <= Dft, else 1
//   factor  = min(fs_buck, 1) beta (1 - r) + r
// Kinks (fs_buck = 1, |fs_buck - 1| = Dft, the 1e-300 denominator floor): a Dual takes the active branch (one-sided),
// where the FD it replaces straddled them.
//
// Plain C++17 + <math.h>, no OpenSees dependency.

#ifndef LadrunoRebarBucklingKernel_h
#define LadrunoRebarBucklingKernel_h

#include <math.h>
#include <LadrunoDual.h>

namespace ladruno_rebar_buckling {

template <class T>
inline T gaFactorT(const T& eps, double e_cross, double lsr, double gaReduction)
{
  const double C = 9.42477796076938;               // 3*pi (upstream constant)
  T denom = e_cross - eps;                          // > 0 in the buckled region
  if (denom < 1.0e-300) denom = T(1.0e-300);
  const T fs_buck = sqrt(32.0 / denom) / (C * lsr);
  T betaLoc = T(1.0);
  const double gamaLoc = 0.1, Dft = 0.25;
  const T stress_diff = fabs(fs_buck - 1.0);
  if (stress_diff <= Dft) betaLoc = 1.0 - gamaLoc * (Dft - stress_diff) / Dft;
  const T m = (fs_buck < 1.0) ? fs_buck : T(1.0);
  return m * betaLoc * (1.0 - gaReduction) + gaReduction;
}

// dfactor/deps at eps by one Dual<1> pass (value: the double factor, bit for bit)
inline double gaFactorDeriv(double eps, double e_cross, double lsr, double gaReduction, double* factor = nullptr)
{
  const ladruno_dual::Dual<1> f = gaFactorT(ladruno_dual::Dual<1>::variable(eps, 0), e_cross, lsr, gaReduction);
  if (factor) *factor = f.v;
  return f.d[0];
}

} // namespace ladruno_rebar_buckling

#endif
