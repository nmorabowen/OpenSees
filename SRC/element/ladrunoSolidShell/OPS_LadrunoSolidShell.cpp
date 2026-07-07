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

// Author: N. Mora-Bowen (Ladruno), 07/2026
//
// Factory for the LadrunoSolidShell element (Tcl + Python) — ADR 66 P5.1.
//
// Usage:
//   element('LadrunoSolidShell', tag, n1..n8, matTag
//           [, '-nz', n]                     # through-thickness points (default 2)
//           [, '-quadz', 'gauss'|'lobatto']  # default: gauss (nz=2) / lobatto (nz>=3)
//           [, '-formulation', 'ans'|'std']  # default ans; std = locking A/B control
//           [, '-lumped']                    # row-sum/HRZ diagonal mass (the
//                                            # explicit-path mass, ADR 66 G9);
//                                            # default = consistent
//           [, '-geom', 'linear'])           # corot/finite reserved (ADR 66 P5.4)
//
// Nodes 1-4 = bottom face (counterclockwise), 5-8 = top face: the node
// ordering DEFINES the thickness direction. Consumes any 3D NDMaterial
// (LadrunoConcrete3D is the headline host).

#include "LadrunoSolidShell.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
#include <FiniteStrainNDMaterial.h>   // reject setTrialF-driven materials (P5.4)

#include <string.h>

void *OPS_LadrunoSolidShell()
{
  if (OPS_GetNumRemainingInputArgs() < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoSolidShell eleTag? n1? ... n8? matTag? "
              "<-nz n> <-quadz gauss|lobatto> <-formulation ans|std> "
              "<-lumped> <-geom linear>\n";
    return 0;
  }

  // tag + 8 node tags + matTag
  int idata[10];
  int num = 10;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data for LadrunoSolidShell\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(idata[9]);
  if (mat == 0) {
    opserr << "WARNING material not found (tag " << idata[9]
           << ") for LadrunoSolidShell element " << idata[0] << endln;
    return 0;
  }

  int nz = 2;
  bool quadzGiven = false;
  int massType = 0;                 // 0 consistent (default), 1 row-sum/HRZ lumped
  LadrunoSolidShell::QuadZ quadz = LadrunoSolidShell::QuadZ::GAUSS;
  LadrunoSolidShell::Formulation formulation =
      LadrunoSolidShell::Formulation::ANS;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-nz") == 0) {
      int n1 = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetIntInput(&n1, &nz) < 0) {
        opserr << "WARNING -nz needs an integer for LadrunoSolidShell "
               << idata[0] << endln;
        return 0;
      }
    }
    else if (strcmp(opt, "-quadz") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -quadz needs a value for LadrunoSolidShell "
               << idata[0] << endln;
        return 0;
      }
      const char *q = OPS_GetString();
      if (strcmp(q, "gauss") == 0 || strcmp(q, "Gauss") == 0)
        quadz = LadrunoSolidShell::QuadZ::GAUSS;
      else if (strcmp(q, "lobatto") == 0 || strcmp(q, "Lobatto") == 0)
        quadz = LadrunoSolidShell::QuadZ::LOBATTO;
      else {
        opserr << "WARNING unknown -quadz '" << q << "' for LadrunoSolidShell "
               << idata[0] << " (use gauss|lobatto)\n";
        return 0;
      }
      quadzGiven = true;
    }
    else if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -formulation needs a value for LadrunoSolidShell "
               << idata[0] << endln;
        return 0;
      }
      const char *f = OPS_GetString();
      if (strcmp(f, "ans") == 0 || strcmp(f, "ANS") == 0)
        formulation = LadrunoSolidShell::Formulation::ANS;
      else if (strcmp(f, "std") == 0 || strcmp(f, "standard") == 0)
        formulation = LadrunoSolidShell::Formulation::STD;
      else {
        opserr << "WARNING unknown -formulation '" << f
               << "' for LadrunoSolidShell " << idata[0] << " (use ans|std)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-lumped") == 0) {
      massType = 1;                 // row-sum == HRZ on the trilinear shape (ADR 66 G9)
    }
    else if (strcmp(opt, "-geom") == 0 || strcmp(opt, "-geometry") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -geom needs a value for LadrunoSolidShell "
               << idata[0] << endln;
        return 0;
      }
      const char *g = OPS_GetString();
      if (strcmp(g, "linear") != 0) {
        // ADR 66 D5: a single nodal-cloud corotation is ill-conditioned for
        // thin shells (cond(S) ~ (L/t)^2); corot ships GUARDED at P5.4, finite
        // with it. Refuse now rather than run silently wrong.  // Ladruno
        opserr << "WARNING LadrunoSolidShell " << idata[0] << ": -geom '" << g
               << "' is reserved (ADR 66 P5.4); v1 is -geom linear only\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING unknown option '" << opt << "' for LadrunoSolidShell "
             << idata[0] << endln;
    }
  }

  // default z-rule: gauss for nz=2, lobatto (face points) for nz>=3 (ADR 66 D3)
  if (!quadzGiven && nz >= 3)
    quadz = LadrunoSolidShell::QuadZ::LOBATTO;

  // validate the (rule, nz) combination up front — the same tables the element
  // uses (gauss 2-5, lobatto 3-7; lobatto nz=2 is the trapezoid rule and
  // cannot integrate the zeta^2 bending stiffness)
  {
    double p[LadrunoSolidShell::NZ_MAX], w[LadrunoSolidShell::NZ_MAX];
    if (nz < LadrunoSolidShell::NZ_MIN || nz > LadrunoSolidShell::NZ_MAX ||
        LadrunoSolidShell::zRuleStatic(quadz, nz, p, w) != 0) {
      opserr << "WARNING LadrunoSolidShell " << idata[0] << ": unsupported -nz "
             << nz << " for -quadz "
             << (quadz == LadrunoSolidShell::QuadZ::GAUSS ? "gauss" : "lobatto")
             << " (gauss: 2-5, lobatto: 3-7); stack elements through the "
                "thickness for finer resolution\n";
      return 0;
    }
  }

  // a finite-strain NDMaterial (setTrialF-driven) cannot be driven by the
  // small-strain setTrialStrain path (the LadrunoBrick GEOM-1 guard)
  if (dynamic_cast<FiniteStrainNDMaterial *>(mat) != 0) {
    opserr << "WARNING LadrunoSolidShell " << idata[0]
           << ": finite-strain NDMaterials (e.g. nDMaterial LogStrain) need "
              "-geom finite, which is reserved (ADR 66 P5.4)\n";
    return 0;
  }

  // Pre-validate the ThreeDimensional copy at PARSE time. The constructor
  // makes 4*nz copies via getCopy("ThreeDimensional"); a material that exists
  // but does not implement a 3D copy returns 0 there and would leave null
  // materialPointers that crash on the first form/commit. Probe once, delete
  // the probe, and fail cleanly here instead.  // Ladruno
  NDMaterial *probe = mat->getCopy("ThreeDimensional");
  if (probe == 0) {
    opserr << "WARNING LadrunoSolidShell " << idata[0] << ": material "
           << idata[9] << " does not provide a 'ThreeDimensional' copy "
              "(LadrunoSolidShell needs a full 3D NDMaterial)\n";
    return 0;
  }
  delete probe;

  return new LadrunoSolidShell(idata[0],
                               idata[1], idata[2], idata[3], idata[4],
                               idata[5], idata[6], idata[7], idata[8],
                               *mat, nz, quadz, formulation, massType);
}
