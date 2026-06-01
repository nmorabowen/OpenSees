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
// std + bbar + uri(stiffness|physical) + eas are all implemented. Only
// uri -hourglass viscous remains refused at construction (rate-form, explicit).

#include "LadrunoBrick.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
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
      // optional numeric coefficient. Only consume the next token if it is a
      // number — otherwise it is the next option (e.g. -lumped), so push it
      // back for the option loop. (Avoids swallowing flags after -hourglass.)
      if (OPS_GetNumRemainingInputArgs() > 0) {
        const char *peek = OPS_GetString();
        char *endp = 0;
        double val = strtod(peek, &endp);
        if (endp != peek && *endp == '\0')
          hgCoeff = val;
        else
          OPS_ResetCurrentInputArg(-1);   // not numeric -> leave it for the loop
      }
    }
    else if (strcmp(opt, "-lumped") == 0 || strcmp(opt, "-lump") == 0) {
      massType = 1;
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

  return new LadrunoBrick(idata[0],
                          idata[1], idata[2], idata[3], idata[4],
                          idata[5], idata[6], idata[7], idata[8],
                          *mat, formulation, bf[0], bf[1], bf[2],
                          massType, hgType, hgCoeff, theDamping);
}
