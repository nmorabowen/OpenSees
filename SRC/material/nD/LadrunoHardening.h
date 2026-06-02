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

// Ladruno: shared isotropic-hardening law for the J2 material family.
//
// The same sig_y(pbar) is consumed by BOTH the 3D nDMaterial LadrunoJ2 and the
// uniaxial LadrunoUniaxialJ2. Keeping it in ONE header is the "oracle contract":
// the uniaxial material must reduce to the 3D model under a uniaxial stress
// state to ~1e-12, which only holds if the isotropic backbone is byte-identical
// in both. (The Armstrong-Frederick kinematic update is provably equivalent
// between the two but is written in each one's natural variable -- the tensor
// plastic multiplier dGamma in 3D, the accumulated plastic strain increment
// dpbar in 1D, related by dpbar = sqrt(2/3) dGamma -- so it is not shared here.)
//
// Header-only. Written: N. Mora-Bowen (Ladruno), 2026.
// See Ladruno_implementation/13_ladruno_uniaxial_j2_adr.md.

#ifndef LadrunoHardening_h
#define LadrunoHardening_h

#include <math.h>

namespace Ladruno {

// Voce saturation + linear isotropic hardening (Abaqus *PLASTIC form):
//   sig_y(p) = sig0 + Qinf (1 - e^{-b p}) + Hiso p
//   Qinf = Hiso = 0  => perfectly plastic ;  Qinf = 0 => linear
inline double yieldStressVoceLinear(double p, double sig0, double Qinf,
                                    double b, double Hiso)
{
  return sig0 + Qinf*(1.0 - exp(-b*p)) + Hiso*p;
}

// d sig_y / d p  (the isotropic plastic modulus contribution H_iso)
inline double yieldSlopeVoceLinear(double p, double sig0, double Qinf,
                                   double b, double Hiso)
{
  return Qinf*b*exp(-b*p) + Hiso;
}

} // namespace Ladruno

#endif
