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
// Factory for the LadrunoBrick20 element (Tcl + Python). ADR 72 P1/P2.
//
// Usage:
//   element('LadrunoBrick20', tag, n1..n20, matTag
//           [, '-formulation', <std|uri>]   # default std
//           [, '-lumped']                   # HRZ lumped mass (P3; default consistent)
//           [, '-b', bx, by, bz]
//           [, '-damp', dampTag])
//
// -formulation std = full 27-pt Gauss; uri (alias: reduced) = uniform 2x2x2
// reduced integration, the C3D20R analog (P2) — prefer std for eigen /
// single-element-thick / point-loaded / soft-support cases (ADR 72 §3.2).
// '-hourglass' is a HARD error by design: the H20@2x2x2 spurious modes are
// non-communicable in solid meshes and no hourglass control exists on this
// element (ADR 72 §2.2 — deliberate asymmetry with LadrunoBrick -formulation
// uri, whose H8@1-pt modes ARE communicable and MUST be stabilized).  // Ladruno

#include "LadrunoBrick20.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
#include <Damping.h>

#include <string.h>

void *OPS_LadrunoBrick20()
{
  if (OPS_GetNumRemainingInputArgs() < 22) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoBrick20 eleTag? n1? ... n20? matTag? "
              "<-formulation std|uri> <-lumped> <-b bx by bz> <-damp dampTag>\n";
    return 0;
  }

  // tag + 20 node tags + matTag
  int idata[22];
  int num = 22;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data for LadrunoBrick20\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(idata[21]);
  if (mat == 0) {
    opserr << "WARNING material not found (tag " << idata[21]
           << ") for LadrunoBrick20 element " << idata[0] << endln;
    return 0;
  }

  // F2: pre-validate the 3D capability HERE, so a plane-stress-only (or other
  // non-3D) material is refused with a clear message and NO element is created
  // (the ctor's clone loop keeps a defensive backstop; no exit(-1) — fork
  // policy).  // Ladruno
  NDMaterial *probe3d = mat->getCopy("ThreeDimensional");
  if (probe3d == 0) {
    opserr << "WARNING LadrunoBrick20 " << idata[0] << ": material " << idata[21]
           << " (" << mat->getClassType() << ") cannot produce a "
              "'ThreeDimensional' copy - LadrunoBrick20 is a 3D solid element "
              "and needs a full 3D NDMaterial (e.g. ElasticIsotropic, LadrunoJ2, "
              "ASDConcrete3D). No element created.\n";
    return 0;
  }
  delete probe3d;

  LadrunoBrick20::Formulation formulation = LadrunoBrick20::Formulation::STD;
  double bf[3] = { 0.0, 0.0, 0.0 };
  int massType = 0;
  Damping *theDamping = 0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -formulation needs a value for LadrunoBrick20 " << idata[0] << endln;
        return 0;
      }
      const char *f = OPS_GetString();
      if (strcmp(f, "std") == 0 || strcmp(f, "standard") == 0)
        formulation = LadrunoBrick20::Formulation::STD;
      else if (strcmp(f, "uri") == 0 || strcmp(f, "reduced") == 0)
        formulation = LadrunoBrick20::Formulation::URI;
      else {
        opserr << "WARNING unknown -formulation '" << f << "' for LadrunoBrick20 "
               << idata[0] << " (use std or uri)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-hourglass") == 0 || strcmp(opt, "-hg") == 0) {
      // HARD error by design (never accepted): the H20 2x2x2 spurious modes are
      // non-communicable in solid meshes (Abaqus applies no control to C3D20R)
      // and NO hourglass control exists on this element — ADR 72 §2.2.
      opserr << "ERROR LadrunoBrick20 " << idata[0]
             << ": -hourglass is not available on this element BY DESIGN "
                "(ADR 72 §2.2): the H20@2x2x2 modes are non-communicable in "
                "solid meshes and no control exists here (unlike LadrunoBrick "
                "-formulation uri, whose H8@1-pt modes must be stabilized). "
                "Use -formulation std for eigen / single-stack / point-loaded "
                "cases.\n";
      return 0;
    }
    else if (strcmp(opt, "-lumped") == 0 || strcmp(opt, "-lump") == 0) {
      // HRZ lumped mass (ADR 72 §3.5, P3): positive-by-construction diagonal
      // lump of the 27-pt consistent mass. Row-sum lumping of an H20 gives
      // NEGATIVE corner masses, so HRZ is the ONLY lumping this element offers.  // Ladruno
      massType = 1;
    }
    else if (strcmp(opt, "-b") == 0 || strcmp(opt, "-bodyForce") == 0) {
      int n3 = 3;
      if (OPS_GetNumRemainingInputArgs() < 3 || OPS_GetDoubleInput(&n3, bf) < 0) {
        opserr << "WARNING invalid -b body force for LadrunoBrick20 " << idata[0] << endln;
        return 0;
      }
    }
    else if (strcmp(opt, "-damp") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        int dampingTag = 0;
        int n1 = 1;
        if (OPS_GetIntInput(&n1, &dampingTag) < 0) return 0;
        theDamping = OPS_getDamping(dampingTag);
        if (theDamping == 0) {
          opserr << "WARNING damping not found (tag " << dampingTag
                 << ") for LadrunoBrick20 " << idata[0] << endln;
          return 0;
        }
      } else {
        // F6: a trailing -damp with no tag is an input error, not a no-op
        opserr << "WARNING LadrunoBrick20 " << idata[0]
               << ": -damp requires a damping tag\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING unknown option '" << opt << "' for LadrunoBrick20 " << idata[0] << endln;
    }
  }

  return new LadrunoBrick20(idata[0],
                            idata[1],  idata[2],  idata[3],  idata[4],
                            idata[5],  idata[6],  idata[7],  idata[8],
                            idata[9],  idata[10], idata[11], idata[12],
                            idata[13], idata[14], idata[15], idata[16],
                            idata[17], idata[18], idata[19], idata[20],
                            *mat, formulation, bf[0], bf[1], bf[2],
                            massType, theDamping);
}
