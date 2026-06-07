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
//           [-rot [-kr {Kr|auto}] [-krAlpha a]            # Phase 2: also tie rotations
//                 {-dNdx N1x N1y [N1z] .. | -gradXi x1..x_ndm}]
//           [-normal nx ny [nz] [-orient ox oy oz]        # Phase 2b: material interface
//            -matN tag | -matT1 tag | -matT2 tag ..  [-corot]]
//           [-enforce {penalty | al}]
//           [-bipenalty {-dtcr dt | -wcap beta}]
//           [-absolute]                                   # opt out of initial-gap capture
//
//   -absolute (alias -noInitGap, ADR 23): by DEFAULT the element captures the gap at
//   activation (setDomain) and drives traction from the RELATIVE gap (g - g0), so an element
//   added to an ALREADY-DEFORMED host (staged construction) is born stress-free instead of
//   yanking the slave to chase the host's accumulated displacement. -absolute disables that
//   capture and keeps the ABSOLUTE tie u_c = sum N_i u_host (legacy v1 / a deliberate
//   snap-to-host). Capture is a no-op when the element is added at the undeformed state.
//
//   -matN/-matT1/-matT2 (ADR 23 Phase 2b, D9): turn the translational tie into a
//   material-driven INTERFACE — each local direction (normal e_0, tangents e_1/e_2)
//   carries a UNIAXIAL material instead of the bare penalty K_u (any direction without
//   a -mat* keeps the K_u penalty). Needs -normal nx ny [nz] for the frame (+ optional
//   -orient for the first tangent). Models: cohesive (softening uniaxial), unilateral
//   gap (ENT/ElasticPPGap on the normal), elastic bedding, bond. Coulomb friction only
//   APPROXIMATE (uncoupled per-direction; rigorous -> LadrunoContact).
//
//   -corot (ADR 23 Phase 2b v2, D9.1): co-rotate the material frame with the host
//   CONTINUUM rotation theta = skew(grad u) at the embedded point, so a contact normal
//   follows the deformed host (frameCur = R(theta)*frame). Material mode only; needs the
//   host gradients dN/dx (supply -dNdx/-gradXi, or -xi on a gradient-capable host — same
//   surface as -rot). The dropped dframe/du tangent term keeps the residual exact but the
//   tangent inexact (frame-objective for explicit, step-halving for implicit, may slow
//   Newton on stiff-normal large-rotation contact). Default = the v1 REFERENCE frame.
//
//   -pressure (ADR 23 Phase 1b): also couple the constrained node's pressure DOF
//   (index ndm) to the host's interpolated pressure, g_p = p_c - sum N_i p_host,i,
//   t_p = K_p*g_p. Requires all coupled nodes to be u-p (ndf >= ndm+1); else U-only.
//   K_p via -kp (default 1e12; host pressure-block auto-scale deferred).
//
//   -rot (ADR 23 Phase 2): also tie the constrained node's ROTATION DOFs to the host
//   CONTINUUM rotation theta = 1/2 curl(u) = skew(grad u) at the embedded point (3D:
//   3 rotations; 2D: the single drilling R_z). Needs the host cartesian gradients
//   dN_i/dx: supply them explicitly with -dNdx (nHost*ndm values, row order), query
//   the host with -gradXi xi.. (needs -host), or let -xi auto-query a gradient-capable
//   host (LadrunoBrick / BezierTet10). K_r via -kr (default 1e12) or -kr auto =
//   krAlpha*K_u*lch^2 (needs -host for lch; default krAlpha=1e3). MUTUALLY EXCLUSIVE
//   with -pressure (the extra DOF is rotation OR pressure — ambiguous in 2D ndf=3).
//   UR is APPROXIMATE / mesh-limited (CST/TET4 host -> single rigid-spin tie; moment-
//   critical embedments need a higher-order host like BezierTet10).
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
#include <Matrix.h>
#include <Element.h>
#include <Domain.h>
#include <UniaxialMaterial.h>   // ADR 23 Phase 2b (D9) — interface materials
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
  bool rot = false;              // -rot: opt-in UR (rotation) tie (Phase 2)
  double Kr = 1.0e12;            // -kr: rotational penalty (numeric form)
  bool krAuto = false;           // -kr auto: K_r = krAlpha·K_u·lch² (needs -host)
  double krAlpha = 1.0e3;        // -krAlpha multiplier for the auto K_r
  Matrix gradN;                  // host cartesian gradients ∂N_i/∂x_j (nHost × ndm)
  bool haveGrad = false;         // gradients supplied (-dNdx / -gradXi / -xi auto)
  Vector xiStored(ndm);          // ξ captured from -xi (reused for -rot auto gradients)
  bool haveXi = false;
  // ADR 23 Phase 2b (D9) — material-driven interface
  Vector normalDir(ndm);         // -normal: frame normal e_0
  bool haveNormal = false;
  Vector orientDir(ndm);         // -orient: first-tangent hint
  bool haveOrient = false;
  int matTag[3] = { -1, -1, -1 };  // -matN / -matT1 / -matT2 uniaxial tags
  bool corot = false;            // -corot (v2): co-rotate the material frame with the host
  bool initGapCapture = true;    // ON by default; -absolute keeps the absolute tie (no offset)

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
      xiStored = xi; haveXi = true;   // reuse ξ for -rot auto gradients (Phase 2)
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
    else if (strcmp(opt, "-rot") == 0) {
      // opt-in UR (rotation) tie (ADR 23 Phase 2). Activated at setDomain only if the
      // constrained node carries the rotation DOFs AND host gradients are present.
      rot = true;
    }
    else if (strcmp(opt, "-kr") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoEmbeddedNode: -kr wants a value or 'auto'\n";
        return 0;
      }
      char krTok[64];
      OPS_GetStringFromAll(krTok, sizeof(krTok));
      if (strcmp(krTok, "auto") == 0) krAuto = true;
      else                            Kr = atof(krTok);
    }
    else if (strcmp(opt, "-krAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &krAlpha) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -krAlpha wants a value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-dNdx") == 0) {
      // explicit host gradients ∂N_i/∂x_j (the rotation analog of -shape): nHost·ndm
      // values in row order [node0: ∂/∂x ∂/∂y (∂/∂z); node1: …]. Works for ANY host.
      Vector tmp(nHost * ndm);
      n = nHost * ndm;
      if (OPS_GetDoubleInput(&n, &tmp(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -dNdx wants " << nHost * ndm
               << " values (nHost x ndm gradients)\n";
        return 0;
      }
      gradN.resize(nHost, ndm);
      for (int i = 0; i < nHost; i++)
        for (int j = 0; j < ndm; j++)
          gradN(i, j) = tmp(i * ndm + j);
      haveGrad = true;
    }
    else if (strcmp(opt, "-gradXi") == 0) {
      // query the host for ∂N/∂x at a natural coordinate (needs -host). The gradient
      // analog of -xi; uses Element::getInterpolationGradients.
      if (hostEle == 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -gradXi requires the -host form; "
               << "use -dNdx\n";
        return 0;
      }
      Vector xg(ndm);
      n = ndm;
      if (OPS_GetDoubleInput(&n, &xg(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -gradXi wants " << ndm
               << " natural coords\n";
        return 0;
      }
      gradN.resize(nHost, ndm);
      if (hostEle->getInterpolationGradients(xg, gradN) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: host element " << hostEle->getTag()
               << " (" << hostEle->getClassType()
               << ") does not implement getInterpolationGradients; supply -dNdx\n";
        return 0;
      }
      if (gradN.noRows() != nHost) {
        opserr << "WARNING LadrunoEmbeddedNode: host returned " << gradN.noRows()
               << " gradient rows but has " << nHost << " nodes\n";
        return 0;
      }
      haveGrad = true;
    }
    else if (strcmp(opt, "-normal") == 0) {
      // D9 local-frame normal e_0 (ndm comps). Tangents are auto-built (or oriented
      // by -orient). Required when any direction carries a material on the normal.
      n = ndm;
      if (OPS_GetDoubleInput(&n, &normalDir(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -normal wants " << ndm << " comps\n";
        return 0;
      }
      haveNormal = true;
    }
    else if (strcmp(opt, "-orient") == 0) {
      n = ndm;
      if (OPS_GetDoubleInput(&n, &orientDir(0)) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: -orient wants " << ndm << " comps\n";
        return 0;
      }
      haveOrient = true;
    }
    else if (strcmp(opt, "-corot") == 0) {
      // D9 v2: co-rotate the material frame with the host continuum rotation at the
      // embedded point. Material mode only (the isotropic/penalty tie is frame-free);
      // needs host gradients ∂N/∂x (resolved below, like -rot).
      corot = true;
    }
    else if (strcmp(opt, "-absolute") == 0 || strcmp(opt, "-noInitGap") == 0) {
      // ADR 23 — opt out of the initial-gap capture: keep the ABSOLUTE tie u_c = Σ N_i u_host
      // (the legacy v1 behavior, or a deliberate snap-to-host). Default = capture ON (the
      // element is born stress-free at the current host configuration ⇒ no staged jolt).
      initGapCapture = false;
    }
    else if (strcmp(opt, "-matN") == 0 || strcmp(opt, "-matT1") == 0 ||
             strcmp(opt, "-matT2") == 0) {
      int slot = (strcmp(opt, "-matN") == 0) ? 0 : (strcmp(opt, "-matT1") == 0) ? 1 : 2;
      if (slot >= ndm) {
        opserr << "WARNING LadrunoEmbeddedNode: " << opt << " needs ndm=3 (no t2 in 2D)\n";
        return 0;
      }
      int mt; n = 1;
      if (OPS_GetIntInput(&n, &mt) < 0) {
        opserr << "WARNING LadrunoEmbeddedNode: " << opt << " wants a uniaxial tag\n";
        return 0;
      }
      matTag[slot] = mt;
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

  // ADR 23 M5/D7-1 — parse-time -rot / -pressure mutual-exclusion guard. The 2D
  // ndf=3 case is ambiguous (u_x,u_y,R_z) vs (u_x,u_y,p); ASD rejects both flags
  // together (ASDEmbeddedNodeElement.cpp:277) and so do we, as the PRIMARY mechanism
  // (the setDomain ndf precedence is only a defensive backstop).
  if (rot && pressure) {
    opserr << "WARNING LadrunoEmbeddedNode: cannot use both -rot and -pressure "
           << "(the constrained node's extra DOF is either rotation OR pressure, "
           << "ambiguous in 2D ndf=3); pick one\n";
    return 0;
  }

  // ADR 23 Phase 2 / 2b v2 — resolve the host gradients ∂N/∂x needed by BOTH the UR tie
  // (-rot) and the D9 frame co-rotation (-corot). If neither -dNdx nor -gradXi was given
  // but the host was given via -xi, query the host at that ξ (the convenience path);
  // otherwise an explicit gradient is required.
  bool needGrad = (rot || corot);
  if (needGrad && !haveGrad && haveXi && hostEle != 0) {
    gradN.resize(nHost, ndm);
    if (hostEle->getInterpolationGradients(xiStored, gradN) < 0) {
      opserr << "WARNING LadrunoEmbeddedNode: -rot/-corot with -xi but host element "
             << hostEle->getTag() << " (" << hostEle->getClassType()
             << ") does not implement getInterpolationGradients; supply -dNdx\n";
      return 0;
    }
    haveGrad = true;
  }
  if (needGrad && !haveGrad) {
    opserr << "WARNING LadrunoEmbeddedNode: -rot/-corot needs host gradients ∂N/∂x — "
           << "supply -dNdx N1x N1y [N1z] …, -gradXi ξ.. (with -host), or -xi with a "
           << "gradient-capable -host\n";
    return 0;
  }

  // ADR 23 Phase 2b (D9) — material-driven interface. matMode is on iff any -mat* was
  // given. -matN (a material on the NORMAL) requires -normal so the frame is defined;
  // a tangent material also needs the frame, so -normal is required for ANY -mat*.
  bool matMode = (matTag[0] >= 0 || matTag[1] >= 0 || matTag[2] >= 0);
  // ADR 23 Phase 2b v2 — -corot only applies to a material interface (the isotropic /
  // penalty tie is already frame-objective). Reject -corot without any -mat* at parse time.
  if (corot && !matMode) {
    opserr << "WARNING LadrunoEmbeddedNode: -corot only applies to a material interface "
           << "(supply -matN/-matT*); the isotropic/penalty tie is already frame-objective\n";
    return 0;
  }
  UniaxialMaterial* matN = 0;
  UniaxialMaterial* matT1 = 0;
  UniaxialMaterial* matT2 = 0;
  if (matMode) {
    if (!haveNormal) {
      opserr << "WARNING LadrunoEmbeddedNode: -matN/-matT* (material interface) requires "
             << "-normal nx ny [nz] to define the local frame\n";
      return 0;
    }
    UniaxialMaterial** slot[3] = { &matN, &matT1, &matT2 };
    for (int d = 0; d < ndm; d++) {
      if (matTag[d] < 0) continue;
      UniaxialMaterial* m = OPS_getUniaxialMaterial(matTag[d]);
      if (m == 0) {
        opserr << "WARNING LadrunoEmbeddedNode: interface material tag " << matTag[d]
               << " (dir " << d << ") not found\n";
        return 0;
      }
      *slot[d] = m;
    }
    if (bipenalty)
      opserr << "INFO LadrunoEmbeddedNode " << tag << ": -bipenalty with a material "
             << "interface uses k_eff = max(initial tangents); a STIFFENING material "
             << "(e.g. a re-contacting gap) breaks the closed-form dt_cr — set -dtcr "
             << "from the closed-contact stiffness or use -enforce al (ADR 23 M6).\n";
  }

  Element* e = new LadrunoEmbeddedNode(tag, ndm, cNode, hostNodes, Nshape, Ku,
                                       hostEleTag, ktAuto, ktAlpha, enforce,
                                       bipenalty, bpMode, bpDt, bpBeta,
                                       pressure, Kp,
                                       rot, Kr, krAuto, krAlpha,
                                       (rot || corot) ? &gradN : 0,
                                       matMode, corot, haveNormal ? &normalDir : 0,
                                       haveOrient ? &orientDir : 0,
                                       matN, matT1, matT2, initGapCapture);
  if (e == 0) {
    opserr << "WARNING LadrunoEmbeddedNode: could not create element\n";
    return 0;
  }
  return e;
}
