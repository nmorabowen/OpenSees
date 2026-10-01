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

// LadrunoNorSandPlaneStrain: the PlaneStrain (order 3, {11,22,12}, eps_33 = 0) wrapper of LadrunoNorSand (the LadrunoSANISAND3D / PlaneStrain
// pattern). A thin class: it fixes the classTag (ND_TAG_LadrunoNorSandPlaneStrain) and the dimensional view; every
// method (state cycle, tangent, sendSelf/recvSelf, responses, refusal/latch) is LadrunoNorSand's.
// LadrunoNorSand::getCopy("PlaneStrain") hands these out, so the user
// declares ONE `nDMaterial LadrunoNorSand ...` and the element picks the view.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoNorSandPlaneStrain_h
#define LadrunoNorSandPlaneStrain_h

#include <classTags.h>
#include "LadrunoNorSand.h"

class LadrunoNorSandPlaneStrain : public LadrunoNorSand {
 public:
  LadrunoNorSandPlaneStrain();                         // null constructor (broker / recvSelf / getCopy)
  LadrunoNorSandPlaneStrain(int tag, const ladruno_norsand::Params& p, const double sigma0[6],
        double v0, double pi0, double density = 0.0);
  ~LadrunoNorSandPlaneStrain();

  const char* getClassType(void) const { return "LadrunoNorSandPlaneStrain"; }

  NDMaterial* getCopy(void);
  using LadrunoNorSand::getCopy;   // keep getCopy(const char*) visible
};

#endif
