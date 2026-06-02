/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Author: N. Mora-Bowen (Ladruno), 05/2026
//
// Factory for the LadrunoBrick element (Tcl + Python).
//
// Usage:
//   element('LadrunoBrick', tag, n1..n8, matTag
//           [, '-formulation', <std|bbar|uri|eas>]   # default std
//           [, '-hourglass', <viscous|stiffness|physical>, coeff]  # uri only
//           [, '-lumped']
//           [, '-b', bx, by, bz]
//           [, '-damp', dampTag])
//
// std + bbar + uri(stiffness|physical|viscous) + eas are all implemented and
// accepted at construction. NOTE: uri -hourglass viscous is rate-form damping
// and EXPLICIT-ONLY — it adds no hourglass stiffness, so the element tangent is
// rank-deficient under an implicit/eigen solver (use stiffness or physical
// there). -damp is only honoured by the std/bbar kernel (see the guard below).

#include "LadrunoBrick.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
#include <FiniteStrainNDMaterial.h>   // -geom finite requires a finite-strain material
#include <SolidTransformation.h>      // geometry-method ids (linear/corot/finite)
#include <Damping.h>

#include <string.h>

void *OPS_LadrunoBrick()
{
  if (OPS_GetNumRemainingInputArgs() < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoBrick eleTag? n1? ... n8? matTag? "
              "<-formulation std|bbar|uri|eas> <-hourglass type coeff> "
              "<-lumped> <-b bx by bz> <-damp dampTag>\n";
    return 0;
  }

  // tag + 8 node tags + matTag
  int idata[10];
  int num = 10;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING: invalid integer data for LadrunoBrick\n";
    return 0;
  }

  NDMaterial *mat = OPS_getNDMaterial(idata[9]);
  if (mat == 0) {
    opserr << "WARNING material not found (tag " << idata[9]
           << ") for LadrunoBrick element " << idata[0] << endln;
    return 0;
  }

  LadrunoBrick::Formulation formulation = LadrunoBrick::Formulation::STD;
  LadrunoBrick::Hourglass   hgType = LadrunoBrick::Hourglass::STIFFNESS;  // default for uri
  double hgCoeff = 0.0;
  double bf[3] = { 0.0, 0.0, 0.0 };
  int massType = 0;
  Damping *theDamping = 0;
  int geomMethodID = SolidTransformation::METHOD_LINEAR;   // -geom (default linear)

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -formulation needs a value for LadrunoBrick " << idata[0] << endln;
        return 0;
      }
      const char *f = OPS_GetString();
      if (strcmp(f, "std") == 0 || strcmp(f, "standard") == 0)
        formulation = LadrunoBrick::Formulation::STD;
      else if (strcmp(f, "bbar") == 0 || strcmp(f, "Bbar") == 0 || strcmp(f, "bBar") == 0)
        formulation = LadrunoBrick::Formulation::BBAR;
      else if (strcmp(f, "uri") == 0 || strcmp(f, "reduced") == 0)
        formulation = LadrunoBrick::Formulation::URI;
      else if (strcmp(f, "eas") == 0)
        formulation = LadrunoBrick::Formulation::EAS;
      else {
        opserr << "WARNING unknown -formulation '" << f << "' for LadrunoBrick "
               << idata[0] << " (use std|bbar|uri|eas)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-hourglass") == 0 || strcmp(opt, "-hg") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -hourglass needs a type for LadrunoBrick " << idata[0] << endln;
        return 0;
      }
      const char *h = OPS_GetString();
      if (strcmp(h, "viscous") == 0)        hgType = LadrunoBrick::Hourglass::VISCOUS;
      else if (strcmp(h, "stiffness") == 0) hgType = LadrunoBrick::Hourglass::STIFFNESS;
      else if (strcmp(h, "physical") == 0)  hgType = LadrunoBrick::Hourglass::PHYSICAL;
      else {
        opserr << "WARNING unknown -hourglass '" << h << "' for LadrunoBrick "
               << idata[0] << " (use viscous|stiffness|physical)\n";
        return 0;
      }
      // optional numeric coefficient. Read it as a NUMBER (OPS_GetDoubleInput
      // accepts Python numeric args AND Tcl numeric strings — OPS_GetString
      // would return "Invalid String Input!" for a Python float and silently
      // drop the coeff). If the next token is not a number it is the next
      // option (e.g. -lumped): GetDoubleInput advances the cursor even on
      // failure, so ResetCurrentInputArg(-1) un-gets exactly that one token
      // and the option loop re-reads it.
      if (OPS_GetNumRemainingInputArgs() > 0) {
        int n1 = 1;
        double tmp = 0.0;
        if (OPS_GetDoubleInput(&n1, &tmp) == 0)
          hgCoeff = tmp;
        else
          OPS_ResetCurrentInputArg(-1);
      }
    }
    else if (strcmp(opt, "-lumped") == 0 || strcmp(opt, "-lump") == 0) {
      massType = 1;
    }
    else if (strcmp(opt, "-geom") == 0 || strcmp(opt, "-geometry") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -geom needs a value for LadrunoBrick " << idata[0] << endln;
        return 0;
      }
      const char *g = OPS_GetString();
      if (strcmp(g, "linear") == 0)
        geomMethodID = SolidTransformation::METHOD_LINEAR;
      else if (strcmp(g, "finite") == 0)
        geomMethodID = SolidTransformation::METHOD_FINITE;
      else if (strcmp(g, "corot") == 0 || strcmp(g, "corotational") == 0)
        geomMethodID = SolidTransformation::METHOD_COROT;
      else {
        opserr << "WARNING unknown -geom '" << g << "' for LadrunoBrick "
               << idata[0] << " (use linear|corot|finite)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-b") == 0 || strcmp(opt, "-bodyForce") == 0) {
      int n3 = 3;
      if (OPS_GetNumRemainingInputArgs() < 3 || OPS_GetDoubleInput(&n3, bf) < 0) {
        opserr << "WARNING invalid -b body force for LadrunoBrick " << idata[0] << endln;
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
                 << ") for LadrunoBrick " << idata[0] << endln;
          return 0;
        }
      }
    }
    else {
      opserr << "WARNING unknown option '" << opt << "' for LadrunoBrick " << idata[0] << endln;
    }
  }

  // -damp is only wired through the std/bbar kernel; the uri/physical/eas
  // condensed single-point kernels do not apply element-level Damping. Drop it
  // with a clear diagnostic for those formulations rather than silently
  // allocating a no-op damping object that is committed every step.
  if (theDamping != 0 &&
      formulation != LadrunoBrick::Formulation::STD &&
      formulation != LadrunoBrick::Formulation::BBAR) {
    opserr << "WARNING LadrunoBrick " << idata[0]
           << ": -damp is only supported with -formulation std|bbar; ignoring "
              "the damping object for this formulation\n";
    theDamping = 0;
  }

  // -geom finite (v3): updated-Lagrangian. v3 supports the std formulation only
  // (bbar+finite = F-bar and uri/eas+finite are reserved), and requires a
  // finite-strain material (driven by setTrialF(F), e.g. nDMaterial LogStrain).
  if (geomMethodID == SolidTransformation::METHOD_FINITE) {
    if (formulation != LadrunoBrick::Formulation::STD) {
      opserr << "WARNING LadrunoBrick " << idata[0]
             << ": -geom finite currently supports only -formulation std "
                "(bbar+finite = F-bar, uri/eas+finite are reserved)\n";
      return 0;
    }
    if (dynamic_cast<FiniteStrainNDMaterial *>(mat) == 0) {
      opserr << "WARNING LadrunoBrick " << idata[0]
             << ": -geom finite requires a finite-strain NDMaterial driven by "
                "setTrialF(F) (e.g. nDMaterial LogStrain); material " << idata[9]
             << " is not one\n";
      return 0;
    }
    if (theDamping != 0) {
      opserr << "WARNING LadrunoBrick " << idata[0]
             << ": -damp is not yet supported with -geom finite (the finite "
                "assembly does not apply element damping); reserved follow-up\n";
      return 0;
    }
  }

  // -geom corot (v2): EICR small-strain corotational. v2 ships std + bbar only.
  // uri/eas under corot are a deferred follow-up (ADR 10 §6/§7): EAS condensation
  // in the corotated frame is unvalidated, and uri's PHYSICAL-hourglass path does
  // not route through the globalize seams at all (frame-inconsistent). Reject the
  // unsupported combos at parse time, mirroring the -geom finite guard.  // Ladruno (sweep #1)
  if (geomMethodID == SolidTransformation::METHOD_COROT &&
      formulation != LadrunoBrick::Formulation::STD &&
      formulation != LadrunoBrick::Formulation::BBAR) {
    opserr << "WARNING LadrunoBrick " << idata[0]
           << ": -geom corot currently supports only -formulation std|bbar "
              "(uri/eas under corot are a deferred follow-up)\n";
    return 0;
  }

  return new LadrunoBrick(idata[0],
                          idata[1], idata[2], idata[3], idata[4],
                          idata[5], idata[6], idata[7], idata[8],
                          *mat, formulation, bf[0], bf[1], bf[2],
                          massType, hgType, hgCoeff, theDamping, geomMethodID);
}
