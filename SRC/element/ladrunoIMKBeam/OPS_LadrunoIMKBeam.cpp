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

// Author: N. Mora-Bowen (Ladruno), 06/2026
//
// Factory for the LadrunoIMKBeam element (Tcl + Python).
//
// Usage:
//   element('LadrunoIMKBeam', tag, iNode, jNode, A, E, G, Jx, Iy, Iz, transfTag
//           [, '-hinge', <both|i|j>]      # ends for the SYMMETRIC -matZ/-matY (default both)
//           [, '-matZ', imkTag]           # strong-axis (Mz) law, both -hinge ends
//           [, '-matY', imkTag]           # weak-axis  (My) law, both -hinge ends
//           [, '-matZi', tag] [, '-matZj', tag]   # per-end strong-axis (asymmetric)
//           [, '-matYi', tag] [, '-matYj', tag]   # per-end weak-axis  (asymmetric)
//           [, '-mass', rho])             # mass per unit length (lumped)
//
// matZ/matY are pre-built uniaxialMaterial tags (the element is material-
// agnostic: steel -> IMKBilin/Bilin, concrete -> ModIMKPeakOriented/pinching).
// Omitting an axis material leaves that axis elastic. The column-face offset is
// supplied through the geomTransf (-jntOffset), not here.
//
// Asymmetric end laws: -matZi/-matZj/-matYi/-matYj set one end/axis explicitly
// and take precedence over the symmetric -matZ/-matY (and are independent of
// -hinge -- naming an end IS the request to hinge it there). Typical use: a
// different deteriorated backbone at each beam end, or a single-axis hinge at
// only one end. e.g. '-matZi', 5, '-matZj', 6 -> strong-axis hinges with laws
// 5 at i and 6 at j; '-matZ', 5, '-matYj', 7 -> symmetric strong-axis law 5 at
// both ends plus an extra weak-axis hinge (law 7) at j only.
//
// "Elastic end" vs "release": omitting an end's material leaves it ELASTIC (it
// still carries moment, full 4EI/L stiffness) -- this is NOT a moment release.
// There is no -release flag yet (deferred). To approximate a true pin (M==0) at
// an end today, give it a near-zero-stiffness Elastic hinge, e.g.
//   uniaxialMaterial Elastic 7 [expr 1e-5*4*$E*$Iz/$L];  ... -matZj 7
// (residual M/M_other ~ 1e-5; keep k above ~1e-8*4EI/L). See
// Ladruno_implementation/{14_ladruno_imk_beam.md sec.8, LEDGER_quirks.md}.

#include "LadrunoIMKBeam.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <CrdTransf.h>
#include <UniaxialMaterial.h>

#include <string.h>

// Parse a "<flag> <matTag>" option and return the looked-up uniaxialMaterial.
// Returns 0 on any error (missing arg / bad int / tag not found); since the
// flags are only consumed when present, a 0 return always means the command is
// malformed and the factory should abort.
static UniaxialMaterial *getMatArg(const char *flag)
{
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING LadrunoIMKBeam -- " << flag << " needs a material tag\n";
    return 0;
  }
  int mtag, num = 1;
  if (OPS_GetIntInput(&num, &mtag) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid " << flag << " tag\n";
    return 0;
  }
  UniaxialMaterial *m = OPS_getUniaxialMaterial(mtag);
  if (m == 0)
    opserr << "WARNING LadrunoIMKBeam -- uniaxialMaterial " << mtag
           << " (" << flag << ") not found\n";
  return m;
}

// Planar builder (LadrunoIMKBeam2d), defined in OPS_LadrunoIMKBeam2d.cpp. The
// single 'LadrunoIMKBeam' command dispatches here to the 2D class when the model
// is planar (ndm==2, ndf==3).
extern void *OPS_LadrunoIMKBeam2d();

void *OPS_LadrunoIMKBeam()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  // ndm-dispatch: one user-facing command, planar or spatial class chosen from
  // the model dimension (same UX as elasticBeamColumn).
  if (ndm == 2 && ndf == 3)
    return OPS_LadrunoIMKBeam2d();

  if (ndm != 3 || ndf != 6) {
    opserr << "WARNING LadrunoIMKBeam -- model must be ndm 3/ndf 6 (3D) "
              "or ndm 2/ndf 3 (2D)\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 10) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoIMKBeam tag iNode jNode A E G Jx Iy Iz "
              "transfTag <-hinge both|i|j> <-matZ tag> <-matY tag> <-mass rho>\n";
    return 0;
  }

  // tag, iNode, jNode
  int iData[3];
  int num = 3;
  if (OPS_GetIntInput(&num, iData) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid tag/iNode/jNode\n";
    return 0;
  }

  // A, E, G, Jx, Iy, Iz
  double dData[6];
  num = 6;
  if (OPS_GetDoubleInput(&num, dData) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid A E G Jx Iy Iz\n";
    return 0;
  }

  // transfTag
  int transfTag;
  num = 1;
  if (OPS_GetIntInput(&num, &transfTag) < 0) {
    opserr << "WARNING LadrunoIMKBeam -- invalid transfTag\n";
    return 0;
  }
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING LadrunoIMKBeam -- CrdTransf " << transfTag << " not found\n";
    return 0;
  }

  // options
  int hingeMode = 0;  // 0 both, 1 i, 2 j  (gates the SYMMETRIC -matZ/-matY)
  UniaxialMaterial *matZ = 0;   // symmetric strong-axis law (applied per -hinge)
  UniaxialMaterial *matY = 0;   // symmetric weak-axis law
  // per-end overrides (asymmetric end laws): take precedence over -matZ/-matY,
  // independent of -hinge (naming an end IS the request to hinge it there).
  UniaxialMaterial *matZi = 0, *matZj = 0, *matYi = 0, *matYj = 0;
  bool setZi = false, setZj = false, setYi = false, setYj = false;
  double rho = 0.0;
  int cMass = 0;
  double GAvy = 0.0, GAvz = 0.0;
  double P0 = 0.0, pmAlpha = 0.5;
  double dI = 0.0, dJ = 0.0;
  bool useImplex = false;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-hinge") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -hinge needs both|i|j|none\n";
        return 0;
      }
      const char *h = OPS_GetString();
      if (strcmp(h, "both") == 0 || strcmp(h, "BOTH") == 0)         hingeMode = 0;
      else if (strcmp(h, "i") == 0 || strcmp(h, "I") == 0)          hingeMode = 1;
      else if (strcmp(h, "j") == 0 || strcmp(h, "J") == 0)          hingeMode = 2;
      else if (strcmp(h, "none") == 0 || strcmp(h, "NONE") == 0)    hingeMode = 3;
      else {
        opserr << "WARNING LadrunoIMKBeam -- -hinge must be both|i|j|none\n";
        return 0;
      }

    } else if (strcmp(opt, "-matZ") == 0) {
      if ((matZ = getMatArg("-matZ")) == 0) return 0;
    } else if (strcmp(opt, "-matY") == 0) {
      if ((matY = getMatArg("-matY")) == 0) return 0;
    } else if (strcmp(opt, "-matZi") == 0) {
      if ((matZi = getMatArg("-matZi")) == 0) return 0;  setZi = true;
    } else if (strcmp(opt, "-matZj") == 0) {
      if ((matZj = getMatArg("-matZj")) == 0) return 0;  setZj = true;
    } else if (strcmp(opt, "-matYi") == 0) {
      if ((matYi = getMatArg("-matYi")) == 0) return 0;  setYi = true;
    } else if (strcmp(opt, "-matYj") == 0) {
      if ((matYj = getMatArg("-matYj")) == 0) return 0;  setYj = true;

    } else if (strcmp(opt, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -mass needs rho\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -mass rho\n";
        return 0;
      }

    } else if (strcmp(opt, "-cMass") == 0 || strcmp(opt, "-cmass") == 0) {
      cMass = 1;

    } else if (strcmp(opt, "-implex") == 0) {
      useImplex = true;

    } else if (strcmp(opt, "-endZones") == 0 || strcmp(opt, "-endZone") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING LadrunoIMKBeam -- -endZones needs dI dJ\n";
        return 0;
      }
      double dz[2] = {0.0, 0.0};
      num = 2;
      if (OPS_GetDoubleInput(&num, dz) < 0 || dz[0] < 0.0 || dz[1] < 0.0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -endZones dI dJ\n";
        return 0;
      }
      dI = dz[0];
      dJ = dz[1];

    } else if (strcmp(opt, "-shear") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING LadrunoIMKBeam -- -shear needs <GAy GAz> or <G Avy Avz>\n";
        return 0;
      }
      double sVals[3] = {0.0, 0.0, 0.0};
      num = 2;
      if (OPS_GetDoubleInput(&num, sVals) < 0 || sVals[0] <= 0.0 || sVals[1] <= 0.0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -shear parameters\n";
        return 0;
      }
      int nShearArgs = 2;
      if (OPS_GetNumRemainingInputArgs() > 0) {
        const char *nxt = OPS_GetString();
        bool nextIsFlag = (nxt != 0 && nxt[0] == '-' &&
                           ((nxt[1] >= 'a' && nxt[1] <= 'z') ||
                            (nxt[1] >= 'A' && nxt[1] <= 'Z')));
        OPS_ResetCurrentInputArg(-1);
        if (!nextIsFlag) {
          num = 1;
          if (OPS_GetDoubleInput(&num, &sVals[2]) < 0 || sVals[2] <= 0.0) {
            opserr << "WARNING LadrunoIMKBeam -- invalid -shear Avz\n";
            return 0;
          }
          nShearArgs = 3;
        }
      }
      if (nShearArgs == 3) {
        GAvy = sVals[0] * sVals[1];
        GAvz = sVals[0] * sVals[2];
      } else {
        GAvy = sVals[0];
        GAvz = sVals[1];
      }

    } else if (strcmp(opt, "-pm") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoIMKBeam -- -pm needs P0 [alpha]\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &P0) < 0 || P0 <= 0.0) {
        opserr << "WARNING LadrunoIMKBeam -- invalid -pm P0 (must be > 0)\n";
        return 0;
      }
      if (OPS_GetNumRemainingInputArgs() > 0) {
        const char *nxt = OPS_GetString();
        bool nextIsFlag = (nxt != 0 && nxt[0] == '-' &&
                           ((nxt[1] >= 'a' && nxt[1] <= 'z') ||
                            (nxt[1] >= 'A' && nxt[1] <= 'Z')));
        OPS_ResetCurrentInputArg(-1);
        if (!nextIsFlag) {
          num = 1;
          if (OPS_GetDoubleInput(&num, &pmAlpha) < 0) {
            opserr << "WARNING LadrunoIMKBeam -- invalid -pm alpha\n";
            return 0;
          }
        }
      }

    } else {
      opserr << "WARNING LadrunoIMKBeam -- unknown option '" << opt << "'\n";
      return 0;
    }
  }

  // map onto the four slots [Zi, Zj, Yi, Yj]: a per-end override wins; otherwise
  // the symmetric -matZ/-matY law applies to the ends selected by -hinge.
  bool hI = (hingeMode == 0 || hingeMode == 1);
  bool hJ = (hingeMode == 0 || hingeMode == 2);

  UniaxialMaterial *mZi = setZi ? matZi : ((matZ && hI) ? matZ : 0);
  UniaxialMaterial *mZj = setZj ? matZj : ((matZ && hJ) ? matZ : 0);
  UniaxialMaterial *mYi = setYi ? matYi : ((matY && hI) ? matY : 0);
  UniaxialMaterial *mYj = setYj ? matYj : ((matY && hJ) ? matY : 0);

  Element *theEle = new LadrunoIMKBeam(iData[0], iData[1], iData[2],
                                       dData[0], dData[1], dData[2], dData[3],
                                       dData[4], dData[5], *theTransf,
                                       mZi, mZj, mYi, mYj, rho, cMass,
                                       GAvy, GAvz, P0, pmAlpha, dI, dJ, useImplex);
  return theEle;
}
