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

// Author: N. Mora-Bowen (Ladruno), 06/2026
//
// Factory for the 2D LadrunoIMKBeam2d element. Not registered as its own
// command: it is reached through OPS_LadrunoIMKBeam(), which dispatches on the
// model's ndm so a single 'LadrunoIMKBeam' command builds the 2D or 3D class.
//
// Usage (ndm 2, ndf 3):
//   element('LadrunoIMKBeam', tag, iNode, jNode, A, E, Iz, transfTag
//           [, '-hinge', <both|i|j>]      # ends for the SYMMETRIC -matZ (default both)
//           [, '-matZ', imkTag]           # strong-axis (Mz) law, both -hinge ends
//           [, '-matZi', tag] [, '-matZj', tag]   # per-end strong-axis (asymmetric)
//           [, '-mass', rho])             # mass per unit length (lumped)
//
// Mirrors the 3D factory but with the planar input signature: no G/Jx/Iy and no
// weak-axis (-matY) options. See OPS_LadrunoIMKBeam.cpp for the 3D form and
// Ladruno_implementation/14_ladruno_imk_beam.md.

#include "LadrunoIMKBeam2d.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <CrdTransf.h>
#include <UniaxialMaterial.h>

#include <string.h>

// Parse a "<flag> <matTag>" option and return the looked-up uniaxialMaterial.
// Returns 0 on any error (missing arg / bad int / tag not found).
static UniaxialMaterial *getMatArg2d(const char *flag)
{
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING LadrunoIMKBeam (2D) -- " << flag << " needs a material tag\n";
    return 0;
  }
  int mtag, num = 1;
  if (OPS_GetIntInput(&num, &mtag) < 0) {
    opserr << "WARNING LadrunoIMKBeam (2D) -- invalid " << flag << " tag\n";
    return 0;
  }
  UniaxialMaterial *m = OPS_getUniaxialMaterial(mtag);
  if (m == 0)
    opserr << "WARNING LadrunoIMKBeam (2D) -- uniaxialMaterial " << mtag
           << " (" << flag << ") not found\n";
  return m;
}

// Reached from OPS_LadrunoIMKBeam() after it has confirmed ndm==2 && ndf==3.
void *OPS_LadrunoIMKBeam2d()
{
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoIMKBeam tag iNode jNode A E Iz transfTag "
              "<-hinge both|i|j> <-matZ tag> <-mass rho>\n";
    return 0;
  }

  // tag, iNode, jNode
  int iData[3];
  int num = 3;
  if (OPS_GetIntInput(&num, iData) < 0) {
    opserr << "WARNING LadrunoIMKBeam (2D) -- invalid tag/iNode/jNode\n";
    return 0;
  }

  // A, E, Iz
  double dData[3];
  num = 3;
  if (OPS_GetDoubleInput(&num, dData) < 0) {
    opserr << "WARNING LadrunoIMKBeam (2D) -- invalid A E Iz\n";
    return 0;
  }

  // transfTag
  int transfTag;
  num = 1;
  if (OPS_GetIntInput(&num, &transfTag) < 0) {
    opserr << "WARNING LadrunoIMKBeam (2D) -- invalid transfTag\n";
    return 0;
  }
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING LadrunoIMKBeam (2D) -- CrdTransf " << transfTag << " not found\n";
    return 0;
  }

  // options
  int hingeMode = 0;  // 0 both, 1 i, 2 j  (gates the SYMMETRIC -matZ)
  UniaxialMaterial *matZ = 0;   // symmetric strong-axis law (applied per -hinge)
  UniaxialMaterial *matZi = 0, *matZj = 0;
  bool setZi = false, setZj = false;
  double rho = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-hinge") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam (2D) -- -hinge needs both|i|j\n";
        return 0;
      }
      const char *h = OPS_GetString();
      if (strcmp(h, "both") == 0)      hingeMode = 0;
      else if (strcmp(h, "i") == 0)    hingeMode = 1;
      else if (strcmp(h, "j") == 0)    hingeMode = 2;
      else {
        opserr << "WARNING LadrunoIMKBeam (2D) -- -hinge must be both|i|j\n";
        return 0;
      }

    } else if (strcmp(opt, "-matZ") == 0) {
      if ((matZ = getMatArg2d("-matZ")) == 0) return 0;
    } else if (strcmp(opt, "-matZi") == 0) {
      if ((matZi = getMatArg2d("-matZi")) == 0) return 0;  setZi = true;
    } else if (strcmp(opt, "-matZj") == 0) {
      if ((matZj = getMatArg2d("-matZj")) == 0) return 0;  setZj = true;

    } else if (strcmp(opt, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam (2D) -- -mass needs rho\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) {
        opserr << "WARNING LadrunoIMKBeam (2D) -- invalid -mass rho\n";
        return 0;
      }

    } else {
      opserr << "WARNING LadrunoIMKBeam (2D) -- unknown option '" << opt << "'\n";
      return 0;
    }
  }

  // map onto the two slots [Zi, Zj]: a per-end override wins; otherwise the
  // symmetric -matZ law applies to the ends selected by -hinge.
  bool hI = (hingeMode == 0 || hingeMode == 1);
  bool hJ = (hingeMode == 0 || hingeMode == 2);

  UniaxialMaterial *mZi = setZi ? matZi : ((matZ && hI) ? matZ : 0);
  UniaxialMaterial *mZj = setZj ? matZj : ((matZ && hJ) ? matZ : 0);

  Element *theEle = new LadrunoIMKBeam2d(iData[0], iData[1], iData[2],
                                         dData[0], dData[1], dData[2],
                                         *theTransf, mZi, mZj, rho);
  return theEle;
}
