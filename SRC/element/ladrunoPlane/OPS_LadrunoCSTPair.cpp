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

// Parser/factory for LadrunoCSTPair (Ladruno fork, ADR 70 P4a). Flag-based:
//   element LadrunoCSTPair tag n1 n2 n3 n4 matTag
//       <-thick t> <-type PlaneStrain> <-rho r> <-body b1 b2>
//
// The macro-element is FINITE-STRAIN ONLY (it exists solely as the T3
// F-bar-Patch volumetric cure, dSNPO §15.1.9) and PlaneStrain only, so the
// material must be a FiniteStrainND2DMaterial (e.g. LogStrain2D) and
// '-type PlaneStress' is refused. n1..n4 must be a CCW quad; the pair is
// split along the n1-n3 diagonal into (n1,n2,n3) and (n1,n3,n4).

#include <LadrunoCSTPair.h>
#include <NDMaterial.h>
#include <FiniteStrainND2DMaterial.h>
#include <elementAPI.h>
#include <string.h>

void *OPS_LadrunoCSTPair()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 2 || ndf != 2) {
    opserr << "WARNING LadrunoCSTPair -- model ndm/ndf not 2/2\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoCSTPair tag n1 n2 n3 n4 matTag "
              "<-thick t> <-type PlaneStrain> <-rho r> <-body b1 b2>\n";
    return 0;
  }

  int idata[5];
  int num = 5;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING LadrunoCSTPair -- invalid tag/node integers\n";
    return 0;
  }

  int matTag;
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING LadrunoCSTPair -- invalid matTag\n";
    return 0;
  }
  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "WARNING LadrunoCSTPair -- material " << matTag << " not found\n";
    return 0;
  }

  double thk = 1.0, rho = 0.0, b1 = 0.0, b2 = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &thk) < 0) { opserr << "WARNING LadrunoCSTPair -- bad -thick\n"; return 0; }

    } else if (strcmp(opt, "-type") == 0) {
      // PlaneStrain only: the finite volume weight dv = J·A·t holds the
      // thickness fixed — it omits the out-of-plane stretch λ = F_33 a finite
      // plane-stress state develops. Refuse anything else (family rule, ADR 70).
      const char *t = OPS_GetString();
      if (strcmp(t, "PlaneStrain") != 0 && strcmp(t, "PlaneStrain2D") != 0) {
        opserr << "WARNING LadrunoCSTPair -- PlaneStrain only (got -type " << t
               << "); the finite path omits the thickness stretch (ADR 70)\n";
        return 0;
      }

    } else if (strcmp(opt, "-rho") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) { opserr << "WARNING LadrunoCSTPair -- bad -rho\n"; return 0; }

    } else if (strcmp(opt, "-body") == 0) {
      double bb[2];
      num = 2;
      if (OPS_GetDoubleInput(&num, bb) < 0) { opserr << "WARNING LadrunoCSTPair -- bad -body\n"; return 0; }
      b1 = bb[0]; b2 = bb[1];

    } else if (strcmp(opt, "-geom") == 0) {
      // the pair has no geometry axis — it IS the finite-strain patch cure.
      // Accept "-geom finite" as a harmless no-op; refuse "linear" loudly.
      const char *gopt = OPS_GetString();
      if (strcmp(gopt, "finite") != 0) {
        opserr << "WARNING LadrunoCSTPair -- the pair is finite-strain only "
                  "(F-bar-Patch); '-geom " << gopt << "' is not available\n";
        return 0;
      }

    } else {
      opserr << "WARNING LadrunoCSTPair -- ignoring unknown option " << opt << "\n";
    }
  }

  // the patch drives both materials by F̄ (setTrialF): a FiniteStrainND2DMaterial
  // is mandatory. Reject before the element exists.
  if (dynamic_cast<FiniteStrainND2DMaterial *>(mat) == 0) {
    opserr << "WARNING LadrunoCSTPair -- needs a FiniteStrainND2DMaterial "
              "(e.g. LogStrain2D); material " << matTag << " is not one\n";
    return 0;
  }

  return new LadrunoCSTPair(idata[0], idata[1], idata[2], idata[3], idata[4],
                            *mat, thk, rho, b1, b2);
}
