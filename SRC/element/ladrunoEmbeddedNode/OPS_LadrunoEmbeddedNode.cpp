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

// Ladruno: OPS parser for the general node-to-host coupling element (ADR 23, Phase 1).
//
//   Host nodes — two forms (mirror LadrunoEmbeddedRebar / ADR 20 §9):
//     explicit : ... tag cNode  nHost h1 ... hNhost  ...
//     by host  : ... tag cNode  -host hostEleTag     ...
//
//   Shape weights N_i — two forms:
//     -shape N1 ... NNhost        (user-supplied; any host)
//     -xi   x1 ... x_ndm          (queried via host getInterpolationWeights; needs -host)
//
//   element LadrunoEmbeddedNode tag cNode {nHost h1..hN | -host eleTag}
//           {-shape N1..NN | -xi x1..x_ndm}
//           [-k {Ku | auto}] [-kAlpha a]
//           [-pressure [-kp Kp]]      # Phase 1b: also tie the pressure DOF (u-p nodes)
//           [-enforce {penalty | al}]
//           [-bipenalty {-dtcr dt | -wcap beta}]
//
//   -pressure (ADR 23 Phase 1b): also couple the constrained node's pressure DOF
//   (index ndm) to the host's interpolated pressure, g_p = p_c - sum N_i p_host,i,
//   t_p = K_p*g_p. Requires all coupled nodes to be u-p (ndf >= ndm+1); else U-only.
//   K_p via -kp (default 1e12; host pressure-block auto-scale deferred).
//
//   -k auto (ADR 23 D3): resolve the isotropic translational penalty from the host's
//   own initial stiffness, K_u = kAlpha * max|K_host(i,i)| (default kAlpha=1e3) —
//   mesh/material-independent conditioning. Requires -host. NB do NOT pass ASD's 1e18
//   as kAlpha (it is not E-scaled — condition-number blow-up).
//
//   -bipenalty (ADR 23 D5): explicit critical-step control via a mass penalty m_p
//   lumped on the constrained node. Budget: -dtcr <dt> or -wcap <beta> (needs -host).
//   Gated on -enforce penalty (ignored with al).
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoEmbeddedNode.h>
#include <ID.h>
#include <Vector.h>
#include <Element.h>
#include <Domain.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>   // atoi/atof

void* OPS_LadrunoEmbeddedNode(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING LadrunoEmbeddedNode: model ndm must be 2 or 3\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient args\n"
           << "Want: element LadrunoEmbeddedNode tag cNode "
           << "{nHost h1..hN | -host eleTag} {-shape N1..NN | -xi x1..x_ndm} "
           << "[-k {Ku|auto}] [-kAlpha a] [-enforce {penalty|al}] "
           << "[-bipenalty {-dtcr dt | -wcap beta}]\n";
    return 0;
  }

  int idata[2];                    // tag, cNode
  int n = 2;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING LadrunoEmbeddedNode: invalid tag/cNode\n";
    return 0;
  }
  int tag = idata[0], cNode = idata[1];

  // --- host nodes: explicit count+list, or `-host eleTag` --------------------
  int nHost = 0;
  ID hostNodes;
  Element* hostEle = 0;
  int hostEleTag = -1;

  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING LadrunoEmbeddedNode: missing host spec (nHost.. or -host)\n";
    return 0;
  }
  char hostTok[128];
  OPS_GetStringFromAll(hostTok, sizeof(hostTok));
  if (strcmp(hostTok, "-host") == 0) {
    int eleTag; n = 1;
    if (OPS_GetIntInput(&n, &eleTag) < 0) {
      opserr << "WARNING LadrunoEmbeddedNode: -host wants a host element tag\n";
      return 0;
    }
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) {
      opserr << "WARNING LadrunoEmbeddedNode: no active domain for -host\n";
      return 0;
    }
    hostEle = theDomain->getElement(eleTag);
    if (hostEle == 0) {
      opserr << "WARNING LadrunoEmbeddedNode: -host element " << eleTag
             << " not found (define the host element before this one)\n";
      return 0;
    }
    hostEleTag = eleTag;
    hostNodes = hostEle->getExternalNodes();
    nHost = hostNodes.Size();
    if (nHost < 1) {
      opserr << "WARNING LadrunoEmbeddedNode: -host element " << eleTag
             << " has no external nodes\n";
      return 0;
    }
  }
  else {
    nHost = atoi(hostTok);
    if (nHost < 1) {
      opserr << "WARNING LadrunoEmbeddedNode: nHost must be >= 1 "
             << "(or use -host eleTag); got '" << hostTok << "'\n";
      return 0;
    }
    hostNodes = ID(nHost);
    for (int i = 0; i < nHost; i++) {
      int h; n = 1;
      if (OPS_GetIntInput(&n, &h) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: invalid host node " << i << "\n";
        return 0;
      }
      hostNodes(i) = h;
    }
  }

  Vector Nshape(nHost);
  bool haveShape = false;
  double Ku = 1.0e12;             // default isotropic penalty (numeric form)
  bool ktAuto = false;           // -k auto
  double ktAlpha = 1.0e3;        // -kAlpha multiplier for the auto K_u
  int enforce = 0;               // 0 = penalty, 1 = al
  bool bipenalty = false;
  int bpMode = 0;                // 0 = -dtcr, 1 = -wcap
  double bpDt = 0.0, bpBeta = 0.0;
  bool bpBudgetSet = false;
  bool pressure = false;         // -pressure: opt-in UP (pressure) tie (Phase 1b)
  double Kp = 1.0e12;            // -kp: pressure penalty (auto-scale deferred)

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if (strcmp(opt, "-shape") == 0) {
      n = nHost;
      if (OPS_GetDoubleInput(&n, &Nshape(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -shape wants " << nHost << " values\n";
        return 0;
      }
      haveShape = true;
    }
    else if (strcmp(opt, "-xi") == 0) {
      if (hostEle == 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -xi requires the -host form "
               << "(no host element to query); use -shape\n";
        return 0;
      }
      Vector xi(ndm);
      n = ndm;
      if (OPS_GetDoubleInput(&n, &xi(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -xi wants " << ndm
               << " natural coords\n";
        return 0;
      }
      if (hostEle->getInterpolationWeights(xi, Nshape) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: host element " << hostEle->getTag()
               << " (" << hostEle->getClassType()
               << ") does not implement getInterpolationWeights; supply -shape\n";
        return 0;
      }
      if (Nshape.Size() != nHost) {
        opserr << "WARNING LadrunoEmbeddedNode: host returned " << Nshape.Size()
               << " weights but has " << nHost << " nodes\n";
        return 0;
      }
      haveShape = true;
    }
    else if (strcmp(opt, "-k") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoEmbeddedNode: -k wants a value or 'auto'\n";
        return 0;
      }
      char kTok[64];
      OPS_GetStringFromAll(kTok, sizeof(kTok));
      if (strcmp(kTok, "auto") == 0) ktAuto = true;
      else                           Ku = atof(kTok);
    }
    else if (strcmp(opt, "-kAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &ktAlpha) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -kAlpha wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-enforce") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoEmbeddedNode: -enforce wants penalty|al\n";
        return 0;
      }
      const char* mode = OPS_GetString();
      if (strcmp(mode, "penalty") == 0)      enforce = 0;
      else if (strcmp(mode, "al") == 0)      enforce = 1;
      else if (strcmp(mode, "nitsche") == 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -enforce nitsche not implemented "
               << "(ADR 23 / ADR 20 §10.7); use penalty|al\n";
        return 0;
      }
      else if (strcmp(mode, "transformation") == 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -enforce transformation is DEFERRED "
               << "INDEFINITELY (ADR 23 / ADR 20 §10.1); use penalty|al\n";
        return 0;
      }
      else {
        opserr << "WARNING LadrunoEmbeddedNode: unknown -enforce '" << mode
               << "' (want penalty|al)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-pressure") == 0) {
      // opt-in UP (pressure) tie (ADR 23 Phase 1b). Activated at setDomain only if all
      // coupled nodes are u-p (ndf >= ndm+1); else degrades to U-only. NB the M5
      // -rot/-pressure mutual-exclusion guard lands with -rot in Phase 2.
      pressure = true;
    }
    else if (strcmp(opt, "-kp") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &Kp) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -kp wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-bipenalty") == 0) {
      bipenalty = true;
    }
    else if (strcmp(opt, "-dtcr") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &bpDt) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -dtcr wants a target step\n";
        return 0;
      }
      bipenalty = true; bpMode = 0; bpBudgetSet = true;
    }
    else if (strcmp(opt, "-wcap") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &bpBeta) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -wcap wants a frequency ratio beta\n";
        return 0;
      }
      bipenalty = true; bpMode = 1; bpBudgetSet = true;
    }
    else {
      opserr << "WARNING LadrunoEmbeddedNode: unknown option '" << opt << "'\n";
      return 0;
    }
  }

  if (!haveShape) {
    opserr << "WARNING LadrunoEmbeddedNode: shape weights required "
           << "(-shape N1..NN, or -xi x1..x_ndm with the -host form)\n";
    return 0;
  }
  if (ktAuto && hostEleTag < 0) {
    opserr << "WARNING LadrunoEmbeddedNode: -k auto requires the -host form "
           << "(no host element to read the stiffness scale from)\n";
    return 0;
  }
  if (bipenalty && !bpBudgetSet) {
    opserr << "WARNING LadrunoEmbeddedNode: -bipenalty needs a budget — "
           << "-dtcr <dt> (explicit) or -wcap <beta> (auto, with -host)\n";
    return 0;
  }
  if (bipenalty && bpMode == 1 && hostEleTag < 0) {
    opserr << "WARNING LadrunoEmbeddedNode: -wcap requires the -host form "
           << "(no host element to read omega_host from); use -dtcr instead\n";
    return 0;
  }
  if (bipenalty && enforce == 1) {
    opserr << "WARNING LadrunoEmbeddedNode: -bipenalty ignored with -enforce al "
           << "(augmented Lagrangian needs no mass penalty)\n";
    bipenalty = false; bpBudgetSet = false;
  }

  Element* e = new LadrunoEmbeddedNode(tag, ndm, cNode, hostNodes, Nshape, Ku,
                                       hostEleTag, ktAuto, ktAlpha, enforce,
                                       bipenalty, bpMode, bpDt, bpBeta,
                                       pressure, Kp);
  if (e == 0) {
    opserr << "WARNING LadrunoEmbeddedNode: could not create element\n";
    return 0;
  }
  return e;
}
