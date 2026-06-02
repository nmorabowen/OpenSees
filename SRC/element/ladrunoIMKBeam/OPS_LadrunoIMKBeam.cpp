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
// Factory for the LadrunoIMKBeam element (Tcl + Python).
//
// Usage:
//   element('LadrunoIMKBeam', tag, iNode, jNode, A, E, G, Jx, Iy, Iz, transfTag
//           [, '-hinge', <both|i|j>]      # which ends hinge (default both)
//           [, '-matZ', imkTag]           # strong-axis (Mz) moment-rotation law
//           [, '-matY', imkTag]           # weak-axis  (My) moment-rotation law
//           [, '-mass', rho])             # mass per unit length (lumped)
//
// matZ/matY are pre-built uniaxialMaterial tags (the element is material-
// agnostic: steel -> IMKBilin/Bilin, concrete -> ModIMKPeakOriented/pinching).
// Omitting an axis material leaves that axis elastic. The column-face offset is
// supplied through the geomTransf (-jntOffset), not here. Asymmetric end laws
// (-matZ differing at i and j) are a v1.x extension.

#include "LadrunoIMKBeam.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <CrdTransf.h>
#include <UniaxialMaterial.h>

#include <string.h>

void *OPS_LadrunoIMKBeam()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 6) {
    opserr << "WARNING LadrunoIMKBeam -- ndm must be 3 and ndf must be 6\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoIMKBeam tag iNode jNode A E G Jx Iy Iz "
              "transfTag <-hinge both|i|j> <-matZ tag> <-matY tag> <-mass rho>\n";
    return 0;
  }

  // tag, iNode, jNode
  int iData[3];
  int num = 3;
  if (OPS_GetIntInput(&num, iData) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid tag/iNode/jNode\n";
    return 0;
  }

  // A, E, G, Jx, Iy, Iz
  double dData[6];
  num = 6;
  if (OPS_GetDoubleInput(&num, dData) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid A E G Jx Iy Iz\n";
    return 0;
  }

  // transfTag
  int transfTag;
  num = 1;
  if (OPS_GetIntInput(&num, &transfTag) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid transfTag\n";
    return 0;
  }
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING LadrunoIMKBeam -- CrdTransf " << transfTag << " not found\n";
    return 0;
  }

  // options
  int hingeMode = 0;  // 0 both, 1 i, 2 j
  UniaxialMaterial *matZ = 0;
  UniaxialMaterial *matY = 0;
  double rho = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-hinge") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -hinge needs both|i|j\n";
        return 0;
      }
      const char *h = OPS_GetString();
      if (strcmp(h, "both") == 0)      hingeMode = 0;
      else if (strcmp(h, "i") == 0)    hingeMode = 1;
      else if (strcmp(h, "j") == 0)    hingeMode = 2;
      else {
        opserr << "WARNING LadrunoIMKBeam -- -hinge must be both|i|j\n";
        return 0;
      }

    } else if (strcmp(opt, "-matZ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -matZ needs a material tag\n";
        return 0;
      }
      int mtag;
      num = 1;
      if (OPS_GetIntInput(&num, &mtag) < 0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -matZ tag\n";
        return 0;
      }
      matZ = OPS_getUniaxialMaterial(mtag);
      if (matZ == 0) {
        opserr << "WARNING LadrunoIMKBeam -- uniaxialMaterial " << mtag
               << " (matZ) not found\n";
        return 0;
      }

    } else if (strcmp(opt, "-matY") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -matY needs a material tag\n";
        return 0;
      }
      int mtag;
      num = 1;
      if (OPS_GetIntInput(&num, &mtag) < 0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -matY tag\n";
        return 0;
      }
      matY = OPS_getUniaxialMaterial(mtag);
      if (matY == 0) {
        opserr << "WARNING LadrunoIMKBeam -- uniaxialMaterial " << mtag
               << " (matY) not found\n";
        return 0;
      }

    } else if (strcmp(opt, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -mass needs rho\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -mass rho\n";
        return 0;
      }

    } else {
      opserr << "WARNING LadrunoIMKBeam -- unknown option '" << opt << "'\n";
      return 0;
    }
  }

  // map hinge ends x axis materials onto the four slots [Zi, Zj, Yi, Yj]
  bool hI = (hingeMode == 0 || hingeMode == 1);
  bool hJ = (hingeMode == 0 || hingeMode == 2);

  UniaxialMaterial *mZi = (matZ && hI) ? matZ : 0;
  UniaxialMaterial *mZj = (matZ && hJ) ? matZ : 0;
  UniaxialMaterial *mYi = (matY && hI) ? matY : 0;
  UniaxialMaterial *mYj = (matY && hJ) ? matY : 0;

  Element *theEle = new LadrunoIMKBeam(iData[0], iData[1], iData[2],
                                       dData[0], dData[1], dData[2], dData[3],
                                       dData[4], dData[5], *theTransf,
                                       mZi, mZj, mYi, mYj, rho);
  return theEle;
}
