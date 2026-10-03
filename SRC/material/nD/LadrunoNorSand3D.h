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

// LadrunoNorSand3D: the 3D (order 6, {11,22,33,12,23,13}) wrapper of LadrunoNorSand (the LadrunoSANISAND3D / PlaneStrain
// pattern). A thin class: it fixes the classTag (ND_TAG_LadrunoNorSand3D) and the dimensional view; every
// method (state cycle, tangent, sendSelf/recvSelf, responses, refusal/latch) is LadrunoNorSand's.
// LadrunoNorSand::getCopy("ThreeDimensional") hands these out, so the user
// declares ONE `nDMaterial LadrunoNorSand ...` and the element picks the view.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoNorSand3D_h
#define LadrunoNorSand3D_h

#include <classTags.h>
#include "LadrunoNorSand.h"

class LadrunoNorSand3D : public LadrunoNorSand {
 public:
  LadrunoNorSand3D();                         // null constructor (broker / recvSelf / getCopy)
  LadrunoNorSand3D(int tag, const ladruno_norsand::Params& p, const double sig0[6],
        double v0, double pi0, double dens = 0.0);
  ~LadrunoNorSand3D();

  const char* getClassType(void) const { return "LadrunoNorSand3D"; }

  NDMaterial* getCopy(void);
  using LadrunoNorSand::getCopy;   // keep getCopy(const char*) visible
};

#endif
