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

// Parser/factory for LadrunoCST (Ladruno fork). Flag-based:
//   element LadrunoCST tag n1 n2 n3 matTag
//       <-type PlaneStrain|PlaneStress> <-thick t> <-rho r>
//       <-body b1 b2> <-pressure p>

#include <LadrunoCST.h>
#include <NDMaterial.h>
#include <elementAPI.h>
#include <string.h>

void *OPS_LadrunoCST()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 2 || ndf != 2) {
    opserr << "WARNING LadrunoCST -- model ndm/ndf not 2/2\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoCST tag n1 n2 n3 matTag "
              "<-type PlaneStrain|PlaneStress> <-thick t> <-rho r> "
              "<-body b1 b2> <-pressure p>\n";
    return 0;
  }

  int idata[4];
  int num = 4;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING LadrunoCST -- invalid tag/node integers\n";
    return 0;
  }

  int matTag;
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING LadrunoCST -- invalid matTag\n";
    return 0;
  }
  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "WARNING LadrunoCST -- material " << matTag << " not found\n";
    return 0;
  }

  double thk = 1.0, rho = 0.0, b1 = 0.0, b2 = 0.0, pressure = 0.0;
  char typeBuf[24];
  strcpy(typeBuf, "PlaneStrain");

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &thk) < 0) { opserr << "WARNING LadrunoCST -- bad -thick\n"; return 0; }

    } else if (strcmp(opt, "-type") == 0) {
      const char *t = OPS_GetString();
      if (strcmp(t, "PlaneStrain") == 0 || strcmp(t, "PlaneStrain2D") == 0)
        strcpy(typeBuf, "PlaneStrain");
      else if (strcmp(t, "PlaneStress") == 0 || strcmp(t, "PlaneStress2D") == 0)
        strcpy(typeBuf, "PlaneStress");
      else { opserr << "WARNING LadrunoCST -- unknown -type " << t << "\n"; return 0; }

    } else if (strcmp(opt, "-rho") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) { opserr << "WARNING LadrunoCST -- bad -rho\n"; return 0; }

    } else if (strcmp(opt, "-body") == 0) {
      double bb[2];
      num = 2;
      if (OPS_GetDoubleInput(&num, bb) < 0) { opserr << "WARNING LadrunoCST -- bad -body\n"; return 0; }
      b1 = bb[0]; b2 = bb[1];

    } else if (strcmp(opt, "-pressure") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &pressure) < 0) { opserr << "WARNING LadrunoCST -- bad -pressure\n"; return 0; }

    } else {
      opserr << "WARNING LadrunoCST -- ignoring unknown option " << opt << "\n";
    }
  }

  return new LadrunoCST(idata[0], idata[1], idata[2], idata[3],
                        *mat, typeBuf, thk, rho, b1, b2, pressure);
}
