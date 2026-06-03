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

// Ladruno: OPS parser for the embedded-reinforcement coupling element.
//   element LadrunoEmbeddedRebar tag rebarNode nHost h1 ... hNhost
//           -shape N1 ... NNhost  -dir dx dy [dz]
//           ( -bond matTag [-bondScale bs] | -perfect kAxial )
//           [-kt kt]
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoEmbeddedRebar.h>
#include <ID.h>
#include <Vector.h>
#include <UniaxialMaterial.h>
#include <elementAPI.h>
#include <string.h>

void* OPS_LadrunoEmbeddedRebar(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING LadrunoEmbeddedRebar: model ndm must be 2 or 3\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING insufficient args\n"
           << "Want: element LadrunoEmbeddedRebar tag rebarNode nHost h1..hNhost "
           << "-shape N1..NNhost -dir dx dy [dz] "
           << "(-bond matTag [-bondScale bs] | -perfect kAxial) [-kt kt]\n";
    return 0;
  }

  int idata[3];                    // tag, rebarNode, nHost
  int n = 3;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING LadrunoEmbeddedRebar: invalid tag/rebarNode/nHost\n";
    return 0;
  }
  int tag = idata[0], rebarNode = idata[1], nHost = idata[2];
  if (nHost < 1) {
    opserr << "WARNING LadrunoEmbeddedRebar: nHost must be >= 1\n";
    return 0;
  }

  ID hostNodes(nHost);
  for (int i = 0; i < nHost; i++) {
    int h; n = 1;
    if (OPS_GetIntInput(&n, &h) < 0) {
      opserr << "WARNING LadrunoEmbeddedRebar: invalid host node " << i << "\n";
      return 0;
    }
    hostNodes(i) = h;
  }

  Vector Nshape(nHost), dir(ndm);
  bool haveShape = false, haveDir = false;
  double kt = 1.0e12;              // default transverse penalty
  double bondScale = 1.0;
  double kAxialPerfect = 0.0;
  int bondTag = -1;
  bool perfect = false;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if (strcmp(opt, "-shape") == 0) {
      n = nHost;
      if (OPS_GetDoubleInput(&n, &Nshape(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -shape wants " << nHost << " values\n";
        return 0;
      }
      haveShape = true;
    }
    else if (strcmp(opt, "-dir") == 0) {
      n = ndm;
      if (OPS_GetDoubleInput(&n, &dir(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -dir wants " << ndm << " values\n";
        return 0;
      }
      haveDir = true;
    }
    else if (strcmp(opt, "-kt") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &kt) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -kt wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-bondScale") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &bondScale) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -bondScale wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-bond") == 0) {
      n = 1;
      if (OPS_GetIntInput(&n, &bondTag) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -bond wants a matTag\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-perfect") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &kAxialPerfect) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -perfect wants kAxial\n";
        return 0;
      }
      perfect = true;
    }
    else {
      opserr << "WARNING LadrunoEmbeddedRebar: unknown option '" << opt << "'\n";
      return 0;
    }
  }

  if (!haveShape || !haveDir) {
    opserr << "WARNING LadrunoEmbeddedRebar: -shape and -dir are required\n";
    return 0;
  }
  if (bondTag < 0 && !perfect) {
    opserr << "WARNING LadrunoEmbeddedRebar: supply -bond matTag or -perfect kAxial\n";
    return 0;
  }

  UniaxialMaterial* bondMat = 0;
  if (bondTag >= 0) {
    bondMat = OPS_getUniaxialMaterial(bondTag);
    if (bondMat == 0) {
      opserr << "WARNING LadrunoEmbeddedRebar: -bond matTag " << bondTag
             << " not found\n";
      return 0;
    }
  }

  Element* e = new LadrunoEmbeddedRebar(tag, ndm, rebarNode, hostNodes, Nshape,
                                        dir, kt, bondScale, bondMat, kAxialPerfect);
  if (e == 0) {
    opserr << "WARNING LadrunoEmbeddedRebar: could not create element\n";
    return 0;
  }
  return e;
}
