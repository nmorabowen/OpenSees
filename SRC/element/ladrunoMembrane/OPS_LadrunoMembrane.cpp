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

#include "LadrunoMembrane.h"

#include <Domain.h>
#include <NDMaterial.h>
#include <OPS_Globals.h>
#include <elementAPI.h>

#include <array>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static const char *
memGetTok(char *buf, int len)
{
  const char *s = OPS_GetStringFromAll(buf, len);
  if (s == 0) {
    return 0;
  }
  if (s != buf) {
    std::strncpy(buf, s, len - 1);
    buf[len - 1] = '\0';
  }
  return buf;
}

static int
memReadIntList(std::vector<int> &out)
{
  out.clear();
  char tok[128];
  while (OPS_GetNumRemainingInputArgs() > 0) {
    if (memGetTok(tok, sizeof(tok)) == 0) {
      break;
    }
    char *endp = 0;
    long v = std::strtol(tok, &endp, 10);
    if (endp == tok || *endp != '\0') {
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    out.push_back(static_cast<int>(v));
  }
  return static_cast<int>(out.size());
}

static int
memReadDoubleList(std::vector<double> &out, int maxCount)
{
  out.clear();
  char tok[128];
  while (OPS_GetNumRemainingInputArgs() > 0 && static_cast<int>(out.size()) < maxCount) {
    if (memGetTok(tok, sizeof(tok)) == 0) {
      break;
    }
    char *endp = 0;
    double v = std::strtod(tok, &endp);
    if (endp == tok || *endp != '\0') {
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    out.push_back(v);
  }
  return static_cast<int>(out.size());
}

void *
OPS_LadrunoMembrane(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING element LadrunoMembrane: insufficient arguments\n";
    return 0;
  }

  int numData = 1;
  int eleTag = 0;
  if (OPS_GetIntInput(&numData, &eleTag) < 0) {
    opserr << "WARNING element LadrunoMembrane: invalid eleTag\n";
    return 0;
  }

  std::vector<int> polyNodes;
  std::vector<int> frameTags;
  std::vector<int> supportNodes;
  std::vector<std::vector<int>> holeLoops;
  std::string method = "tributary";
  std::array<double, 3> spanAxis = {1.0, 0.0, 0.0};
  bool hasSpanAxis = false;
  double orthoRatio = 1.0;
  int segments = 0;
  double meshSize = 0.0;
  double thick = 0.0;
  bool hasMatTag = false;
  int matTag = 0;
  double E_direct = 0.0;
  double nu_direct = 0.2;
  double rho = 0.0;
  bool lumped = false;
  double pressure = 0.0;
  double b1 = 0.0, b2 = 0.0, b3 = 0.0;

  char tok[128];
  while (OPS_GetNumRemainingInputArgs() > 0) {
    if (memGetTok(tok, sizeof(tok)) == 0) {
      break;
    }
    if (std::strcmp(tok, "-nodes") == 0) {
      if (memReadIntList(polyNodes) < 3) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -nodes requires at least 3 node tags\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-frames") == 0) {
      if (memReadIntList(frameTags) < 1) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -frames requires at least 1 frame element tag\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-supportNodes") == 0) {
      if (memReadIntList(supportNodes) < 1) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -supportNodes requires at least 1 node tag\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-hole") == 0) {
      std::vector<int> hole;
      if (memReadIntList(hole) < 3) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -hole requires at least 3 node tags\n";
        return 0;
      }
      holeLoops.push_back(hole);
    } else if (std::strcmp(tok, "-method") == 0) {
      char mTok[128];
      if (OPS_GetNumRemainingInputArgs() < 1 || memGetTok(mTok, sizeof(mTok)) == 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -method requires a method name\n";
        return 0;
      }
      std::string mStr;
      for (int i = 0; mTok[i] != '\0'; ++i) {
        mStr.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(mTok[i]))));
      }
      if (mStr != "tributary" && mStr != "oneway" &&
          mStr != "fea_poisson" && mStr != "fea_plate") {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": unknown -method '" << mTok
               << "' (expected tributary, oneway, fea_poisson, or fea_plate)\n";
        return 0;
      }
      method = mStr;
    } else if (std::strcmp(tok, "-spanAxis") == 0) {
      std::vector<double> axVals;
      int nAx = memReadDoubleList(axVals, 3);
      if (nAx < 2) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -spanAxis requires 2 or 3 components\n";
        return 0;
      }
      spanAxis[0] = axVals[0];
      spanAxis[1] = axVals[1];
      spanAxis[2] = (nAx >= 3) ? axVals[2] : 0.0;
      hasSpanAxis = true;
    } else if (std::strcmp(tok, "-orthoRatio") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &orthoRatio) < 0 ||
          orthoRatio <= 0.0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -orthoRatio requires a positive value\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-segments") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetIntInput(&numData, &segments) < 0 ||
          segments <= 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -segments requires a positive integer\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-meshSize") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &meshSize) < 0 ||
          meshSize < 0.0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -meshSize requires a non-negative value\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-thick") == 0 || std::strcmp(tok, "-thickness") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &thick) < 0 ||
          thick <= 0.0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -thick requires a positive thickness\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-mat") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetIntInput(&numData, &matTag) < 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -mat requires an NDMaterial tag\n";
        return 0;
      }
      hasMatTag = true;
    } else if (std::strcmp(tok, "-stiff") == 0) {
      double sVals[3] = {0.0, 0.0, 0.0};
      numData = 3;
      if (OPS_GetNumRemainingInputArgs() < 3 || OPS_GetDoubleInput(&numData, sVals) < 0 ||
          sVals[0] <= 0.0 || sVals[1] <= -1.0 || sVals[1] >= 0.5 || sVals[2] <= 0.0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -stiff requires E > 0, -1 < nu < 0.5, t > 0\n";
        return 0;
      }
      E_direct = sVals[0];
      nu_direct = sVals[1];
      thick = sVals[2];
    } else if (std::strcmp(tok, "-E") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &E_direct) < 0 ||
          E_direct <= 0.0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -E requires a positive modulus\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-nu") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &nu_direct) < 0 ||
          nu_direct <= -1.0 || nu_direct >= 0.5) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -nu requires -1 < nu < 0.5\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-rho") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &rho) < 0 ||
          rho < 0.0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -rho requires a non-negative density\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-lumped") == 0) {
      lumped = true;
    } else if (std::strcmp(tok, "-pressure") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &pressure) < 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -pressure requires a numeric value\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-b1") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &b1) < 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -b1 requires a numeric value\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-b2") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &b2) < 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -b2 requires a numeric value\n";
        return 0;
      }
    } else if (std::strcmp(tok, "-b3") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &b3) < 0) {
        opserr << "WARNING element LadrunoMembrane " << eleTag
               << ": -b3 requires a numeric value\n";
        return 0;
      }
    } else {
      opserr << "WARNING element LadrunoMembrane " << eleTag
             << ": unknown option '" << tok << "'\n";
      return 0;
    }
  }

  if (polyNodes.size() < 3) {
    opserr << "WARNING element LadrunoMembrane " << eleTag
           << ": -nodes with at least 3 polygon vertices is required\n";
    return 0;
  }

  NDMaterial *theMat = nullptr;
  if (hasMatTag) {
    theMat = OPS_getNDMaterial(matTag);
    if (theMat == nullptr) {
      opserr << "WARNING element LadrunoMembrane " << eleTag
             << ": NDMaterial " << matTag << " not found\n";
      return 0;
    }
    if (thick <= 0.0) {
      opserr << "WARNING element LadrunoMembrane " << eleTag
             << ": -mat requires positive -thick\n";
      return 0;
    }
  }

  LadrunoMembrane *theEle = new LadrunoMembrane(
      eleTag,
      polyNodes,
      frameTags,
      supportNodes,
      holeLoops,
      method,
      spanAxis,
      hasSpanAxis,
      orthoRatio,
      segments,
      meshSize,
      thick,
      theMat,
      E_direct,
      nu_direct,
      rho,
      lumped,
      pressure,
      b1,
      b2,
      b3);

  Domain *theDomain = OPS_GetDomain();
  if (theDomain != nullptr) {
    if (!theEle->validateNodesAndBasis(theDomain)) {
      delete theEle;
      return 0;
    }
  }

  return theEle;
}
