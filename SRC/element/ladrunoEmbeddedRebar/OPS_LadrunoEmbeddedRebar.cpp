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
//           [-corot {-xiB b1..b_ndm | -shapeB N1..NN}]
//           [-enforce {penalty | al}]
//
//   -enforce (ADR 20 §10.1): constraint treatment. penalty (default, Mode P) or
//   al (augmented Lagrangian — per-step Uzawa drives the perfect-bond gap -> 0 at
//   MODERATE kt; bond-slip axial untouched). nitsche is not built; transformation
//   is deferred indefinitely (both rejected at parse time).
//
//   -corot (ADR 20 §10.5): co-rotate the bar axis each step from current host
//   geometry (secant embed-point -> point B), so the axial/transverse split stays
//   frame-objective under large host rotation. Needs a point B along the bar via
//   -xiB (host query) or -shapeB (explicit weights). Default OFF (frozen -dir).
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
  // peek the host-spec token. Use OPS_GetStringFromAll (NOT OPS_GetString): in
  // openseespy the explicit-form nHost arrives as a typed int and OPS_GetString
  // rejects it ("Invalid String Input!"); GetStringFromAll stringifies any arg.
  char hostTok[128];
  OPS_GetStringFromAll(hostTok, sizeof(hostTok));
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
  bool corot = false;             // -corot: co-rotate the bar axis (ADR §10.5)
  Vector NshapeB(nHost);          // weights at point B along the bar (corot)
  bool haveB = false;
  int enforce = 0;                // -enforce: 0=penalty (default), 1=al (ADR §10.1)
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
    else if (strcmp(opt, "-corot") == 0) {
      corot = true;                 // co-rotate the bar axis each step (ADR §10.5)
    }
    else if (strcmp(opt, "-enforce") == 0) {
      // constraint-enforcement strategy (ADR 20 §10.1). penalty (default) | al;
      // nitsche not built; transformation deferred indefinitely.
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoEmbeddedRebar: -enforce wants penalty|al\n";
        return 0;
      }
      const char* mode = OPS_GetString();
      if (strcmp(mode, "penalty") == 0)      enforce = 0;
      else if (strcmp(mode, "al") == 0)      enforce = 1;
      else if (strcmp(mode, "nitsche") == 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -enforce nitsche not implemented "
               << "(ADR 20 §10.7, research-grade); use penalty|al\n";
        return 0;
      }
      else if (strcmp(mode, "transformation") == 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -enforce transformation is DEFERRED "
               << "INDEFINITELY (ADR 20 §10.1 / D2 — needs a multi-retained MPC); "
               << "use penalty|al\n";
        return 0;
      }
      else {
        opserr << "WARNING LadrunoEmbeddedRebar: unknown -enforce '" << mode
               << "' (want penalty|al)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-shapeB") == 0) {
      // explicit weights at point B along the bar (the corot secant endpoint).
      n = nHost;
      if (OPS_GetDoubleInput(&n, &NshapeB(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -shapeB wants " << nHost << " values\n";
        return 0;
      }
      haveB = true;
    }
    else if (strcmp(opt, "-xiB") == 0) {
      // query the host for the point-B weights (corot secant). Requires -host.
      if (hostEle == 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -xiB requires the -host form; "
               << "use -shapeB instead\n";
        return 0;
      }
      Vector xiB(ndm);
      n = ndm;
      if (OPS_GetDoubleInput(&n, &xiB(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: -xiB wants " << ndm
               << " natural coords\n";
        return 0;
      }
      if (hostEle->getInterpolationWeights(xiB, NshapeB) < 0) {
        opserr << "WARNING LadrunoEmbeddedRebar: host element " << hostEle->getTag()
               << " does not implement getInterpolationWeights; supply -shapeB\n";
        return 0;
      }
      if (NshapeB.Size() != nHost) {
        opserr << "WARNING LadrunoEmbeddedRebar: host returned " << NshapeB.Size()
               << " B-weights but has " << nHost << " nodes\n";
        return 0;
      }
      haveB = true;
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
      // GetStringFromAll (not OPS_GetString): a numeric -kt arrives typed in
      // openseespy and OPS_GetString would reject it.
      char ktTok[64];
      OPS_GetStringFromAll(ktTok, sizeof(ktTok));
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
  if (corot && !haveB) {
    opserr << "WARNING LadrunoEmbeddedRebar: -corot requires a point-B along the "
           << "bar (-xiB x1..x_ndm with -host, or -shapeB N1..NN)\n";
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
                                        hostEleTag, ktAuto, ktAlpha,
                                        corot, corot ? &NshapeB : 0, enforce);
  if (e == 0) {
    opserr << "WARNING LadrunoEmbeddedRebar: could not create element\n";
    return 0;
  }
  return e;
}
