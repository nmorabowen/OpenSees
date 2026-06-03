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
//
//   Host nodes — two forms (ADR 20 §9):
//     explicit : ... tag rebarNode  nHost h1 ... hNhost  ...
//     by host  : ... tag rebarNode  -host hostEleTag     ...   (host must already
//                exist; its external-node list IS the host-node list)
//
//   Shape weights N_i — two forms:
//     -shape N1 ... NNhost        (user-supplied; works for any host)
//     -xi   x1 ... x_ndm          (queried from the host element — requires -host;
//                                  the host must implement getInterpolationWeights,
//                                  e.g. LadrunoBrick (ξ,η,ζ) / BezierTet10 (L1,L2,L3))
//
//   element LadrunoEmbeddedRebar tag rebarNode {nHost h1..hN | -host eleTag}
//           {-shape N1..NN | -xi x1..x_ndm}  -dir dx dy [dz]
//           ( -bond matTag [-bondScale bs] | -perfect kAxial )
//           [-kt {kt | auto}] [-ktAlpha a]
//
//   -kt auto (ADR 20 §10.2a): resolve the transverse penalty from the host
//   element's own initial stiffness, kt = ktAlpha * max|K_host(i,i)| (default
//   ktAlpha = 1e3) — mesh/material-independent conditioning. Requires -host.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoEmbeddedRebar.h>
#include <ID.h>
#include <Vector.h>
#include <UniaxialMaterial.h>
#include <Element.h>
#include <Domain.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>   // atoi

void* OPS_LadrunoEmbeddedRebar(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING LadrunoEmbeddedRebar: model ndm must be 2 or 3\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient args\n"
           << "Want: element LadrunoEmbeddedRebar tag rebarNode "
           << "{nHost h1..hN | -host eleTag} {-shape N1..NN | -xi x1..x_ndm} "
           << "-dir dx dy [dz] (-bond matTag [-bondScale bs] | -perfect kAxial) [-kt kt]\n";
    return 0;
  }

  int idata[2];                    // tag, rebarNode
  int n = 2;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING LadrunoEmbeddedRebar: invalid tag/rebarNode\n";
    return 0;
  }
  int tag = idata[0], rebarNode = idata[1];

  // --- host nodes: either an explicit count+list, or `-host eleTag` ----------
  int nHost = 0;
  ID hostNodes;
  Element* hostEle = 0;            // set only in the -host form; needed for -xi
  int hostEleTag = -1;             // -host tag, threaded to the element for -kt auto

  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING LadrunoEmbeddedRebar: missing host spec (nHost.. or -host)\n";
    return 0;
  }
  const char* hostTok = OPS_GetString();
  if (strcmp(hostTok, "-host") == 0) {
    int eleTag; n = 1;
    if (OPS_GetIntInput(&n, &eleTag) < 0) {
      opserr << "WARNING LadrunoEmbeddedRebar: -host wants a host element tag\n";
      return 0;
    }
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) {
      opserr << "WARNING LadrunoEmbeddedRebar: no active domain for -host\n";
      return 0;
    }
    hostEle = theDomain->getElement(eleTag);
    if (hostEle == 0) {
      opserr << "WARNING LadrunoEmbeddedRebar: -host element " << eleTag
             << " not found (define the host solid before the rebar element)\n";
      return 0;
    }
    hostEleTag = eleTag;                        // remember for -kt auto
    hostNodes = hostEle->getExternalNodes();   // ID copy
    nHost = hostNodes.Size();
    if (nHost < 1) {
      opserr << "WARNING LadrunoEmbeddedRebar: -host element " << eleTag
             << " has no external nodes\n";
      return 0;
    }
  }
  else {
    // explicit form: the peeked token is nHost
    nHost = atoi(hostTok);
    if (nHost < 1) {
      opserr << "WARNING LadrunoEmbeddedRebar: nHost must be >= 1 "
             << "(or use -host eleTag); got '" << hostTok << "'\n";
      return 0;
    }
    hostNodes = ID(nHost);
    for (int i = 0; i < nHost; i++) {
      int h; n = 1;
      if (OPS_GetIntInput(&n, &h) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: invalid host node " << i << "\n";
        return 0;
      }
      hostNodes(i) = h;
    }
  }

  Vector Nshape(nHost), dir(ndm);
  bool haveShape = false, haveDir = false;
  double kt = 1.0e12;              // default transverse penalty (numeric form)
  bool ktAuto = false;            // -kt auto: resolve kt from host stiffness
  double ktAlpha = 1.0e3;         // -ktAlpha: multiplier for the auto kt
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
    else if (strcmp(opt, "-xi") == 0) {
      // query the host element for the shape-function weights at this natural
      // coordinate (ADR 20 §9 single source of truth). Requires the -host form.
      if (hostEle == 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -xi requires the -host form "
               << "(no host element to query); use -shape instead\n";
        return 0;
      }
      Vector xi(ndm);
      n = ndm;
      if (OPS_GetDoubleInput(&n, &xi(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -xi wants " << ndm
               << " natural coords\n";
        return 0;
      }
      if (hostEle->getInterpolationWeights(xi, Nshape) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: host element "
               << hostEle->getTag() << " (" << hostEle->getClassType()
               << ") does not implement getInterpolationWeights; supply -shape\n";
        return 0;
      }
      if (Nshape.Size() != nHost) {
        opserr << "WARNING LadrunoEmbeddedRebar: host returned " << Nshape.Size()
               << " weights but has " << nHost << " nodes\n";
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
      // numeric value, OR the sentinel 'auto' (resolve from host stiffness at
      // first assembly: kt = ktAlpha * max|K_host(i,i)|, ADR 20 §10.2a).
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoEmbeddedRebar: -kt wants a value or 'auto'\n";
        return 0;
      }
      const char* ktTok = OPS_GetString();
      if (strcmp(ktTok, "auto") == 0) {
        ktAuto = true;
      } else {
        kt = atof(ktTok);
      }
    }
    else if (strcmp(opt, "-ktAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &ktAlpha) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -ktAlpha wants a value\n";
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

  if (!haveShape) {
    opserr << "WARNING LadrunoEmbeddedRebar: shape weights required "
           << "(-shape N1..NN, or -xi x1..x_ndm with the -host form)\n";
    return 0;
  }
  if (!haveDir) {
    opserr << "WARNING LadrunoEmbeddedRebar: -dir is required\n";
    return 0;
  }
  if (bondTag < 0 && !perfect) {
    opserr << "WARNING LadrunoEmbeddedRebar: supply -bond matTag or -perfect kAxial\n";
    return 0;
  }
  if (ktAuto && hostEleTag < 0) {
    opserr << "WARNING LadrunoEmbeddedRebar: -kt auto requires the -host form "
           << "(no host element to read the stiffness scale from)\n";
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
                                        dir, kt, bondScale, bondMat, kAxialPerfect,
                                        hostEleTag, ktAuto, ktAlpha);
  if (e == 0) {
    opserr << "WARNING LadrunoEmbeddedRebar: could not create element\n";
    return 0;
  }
  return e;
}
