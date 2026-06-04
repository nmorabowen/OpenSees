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

// Parser/factory for LadrunoQuad (Ladruno fork). Flag-based:
//   element LadrunoQuad tag n1 n2 n3 n4 matTag
//       <-formulation std|bbar|ssp|eas> <-type PlaneStrain|PlaneStress>
//       <-thick t> <-rho r> <-body b1 b2> <-pressure p>

#include <LadrunoQuad.h>
#include <NDMaterial.h>
#include <elementAPI.h>
#include <string.h>

void *OPS_LadrunoQuad()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 2 || ndf != 2) {
    opserr << "WARNING LadrunoQuad -- model ndm/ndf not 2/2\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoQuad tag n1 n2 n3 n4 matTag "
              "<-formulation std|bbar|ssp|eas> <-type PlaneStrain|PlaneStress> "
              "<-thick t> <-rho r> <-body b1 b2> <-pressure p>\n";
    return 0;
  }

  int idata[5];
  int num = 5;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING LadrunoQuad -- invalid tag/node integers\n";
    return 0;
  }

  int matTag;
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING LadrunoQuad -- invalid matTag\n";
    return 0;
  }
  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "WARNING LadrunoQuad -- material " << matTag << " not found\n";
    return 0;
  }

  // defaults
  double thk = 1.0, rho = 0.0, b1 = 0.0, b2 = 0.0, pressure = 0.0;
  char typeBuf[24];
  strcpy(typeBuf, "PlaneStrain");
  LadrunoQuad::Formulation form = LadrunoQuad::Formulation::STD;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &thk) < 0) { opserr << "WARNING LadrunoQuad -- bad -thick\n"; return 0; }

    } else if (strcmp(opt, "-type") == 0) {
      const char *t = OPS_GetString();
      if (strcmp(t, "PlaneStrain") == 0 || strcmp(t, "PlaneStrain2D") == 0)
        strcpy(typeBuf, "PlaneStrain");
      else if (strcmp(t, "PlaneStress") == 0 || strcmp(t, "PlaneStress2D") == 0)
        strcpy(typeBuf, "PlaneStress");
      else { opserr << "WARNING LadrunoQuad -- unknown -type " << t << "\n"; return 0; }

    } else if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      const char *f = OPS_GetString();
      if (strcmp(f, "std") == 0)       form = LadrunoQuad::Formulation::STD;
      else if (strcmp(f, "bbar") == 0) form = LadrunoQuad::Formulation::BBAR;
      else if (strcmp(f, "ssp") == 0)  form = LadrunoQuad::Formulation::SSP;
      else if (strcmp(f, "eas") == 0)  form = LadrunoQuad::Formulation::EAS;
      else { opserr << "WARNING LadrunoQuad -- unknown -formulation " << f << "\n"; return 0; }

    } else if (strcmp(opt, "-rho") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) { opserr << "WARNING LadrunoQuad -- bad -rho\n"; return 0; }

    } else if (strcmp(opt, "-body") == 0) {
      double bb[2];
      num = 2;
      if (OPS_GetDoubleInput(&num, bb) < 0) { opserr << "WARNING LadrunoQuad -- bad -body\n"; return 0; }
      b1 = bb[0]; b2 = bb[1];

    } else if (strcmp(opt, "-pressure") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &pressure) < 0) { opserr << "WARNING LadrunoQuad -- bad -pressure\n"; return 0; }

    } else {
      opserr << "WARNING LadrunoQuad -- ignoring unknown option " << opt << "\n";
    }
  }

  if (form == LadrunoQuad::Formulation::EAS) {
    opserr << "WARNING LadrunoQuad -- formulation 'eas' is reserved but not yet "
              "implemented (ADR 25 Phase 3); use std, bbar, or ssp\n";
    return 0;
  }

  if (form == LadrunoQuad::Formulation::BBAR && strcmp(typeBuf, "PlaneStress") == 0) {
    opserr << "WARNING LadrunoQuad -- '-formulation bbar' is for PlaneStrain only "
              "(volumetric locking is a plane-strain/incompressible issue)\n";
    return 0;
  }

  return new LadrunoQuad(idata[0], idata[1], idata[2], idata[3], idata[4],
                         *mat, typeBuf, thk, form, rho, b1, b2, pressure);
}
