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
//       <-type PlaneStrain|PlaneStress> <-geom linear|finite>
//       <-thick t> <-rho r> <-body b1 b2> <-pressure p>

#include <LadrunoCST.h>
#include <NDMaterial.h>
#include <FiniteStrainND2DMaterial.h>   // Ladruno (ADR 70): -geom finite material check
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
              "<-type PlaneStrain|PlaneStress> <-geom linear|finite> "
              "<-thick t> <-rho r> <-body b1 b2> <-pressure p>\n";
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
  double bvB1 = 0.0, bvB2 = 0.0;   // Ladruno (W2-E1): explicit bulk-viscosity coeffs (off by default)
  char typeBuf[24];
  strcpy(typeBuf, "PlaneStrain");
  LadrunoCST::Geom geom = LadrunoCST::Geom::LINEAR;   // Ladruno (ADR 70)

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

    } else if (strcmp(opt, "-geom") == 0) {   // Ladruno (ADR 70)
      const char *gopt = OPS_GetString();
      if (strcmp(gopt, "linear") == 0)      geom = LadrunoCST::Geom::LINEAR;
      else if (strcmp(gopt, "finite") == 0) geom = LadrunoCST::Geom::FINITE;
      else { opserr << "WARNING LadrunoCST -- unknown -geom " << gopt
                    << " (want linear|finite)\n"; return 0; }

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

    } else if (strcmp(opt, "-bulkViscosity") == 0 || strcmp(opt, "-bv") == 0) {
      // Ladruno (W2-E1): explicit bulk viscosity reads TWO doubles (linear b1,
      // quadratic b2). Both must be >= 0; a negative coeff is warned and ignored.
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING LadrunoCST -- -bulkViscosity needs two values (b1 b2)\n";
        return 0;
      }
      double bv[2] = { 0.0, 0.0 };
      num = 2;
      if (OPS_GetDoubleInput(&num, bv) < 0) { opserr << "WARNING LadrunoCST -- bad -bulkViscosity\n"; return 0; }
      if (bv[0] < 0.0) { opserr << "WARNING LadrunoCST -- -bulkViscosity b1 < 0; ignoring (b1=0)\n"; bv[0] = 0.0; }
      if (bv[1] < 0.0) { opserr << "WARNING LadrunoCST -- -bulkViscosity b2 < 0; ignoring (b2=0)\n"; bv[1] = 0.0; }
      bvB1 = bv[0];
      bvB2 = bv[1];

    } else {
      opserr << "WARNING LadrunoCST -- ignoring unknown option " << opt << "\n";
    }
  }

  // Ladruno (ADR 70): -geom finite is plane-strain only. formFinite's current
  // volume weight dv = J·detJ0·thickness·w holds the thickness fixed — it omits
  // the out-of-plane stretch λ = F_33 that a finite plane-stress state develops.
  // Refuse here — before the element exists — rather than assemble a silently
  // wrong finite plane-stress residual/tangent. Mirrors LadrunoQuad.
  if (geom == LadrunoCST::Geom::FINITE && strcmp(typeBuf, "PlaneStress") == 0) {
    opserr << "WARNING LadrunoCST -- '-geom finite' is PlaneStrain only "
              "(finite plane-stress omits the thickness stretch lambda in the "
              "volume weight; use -type PlaneStrain, ADR 70)\n";
    return 0;
  }

  // Ladruno (ADR 70): -geom finite drives the material by the deformation
  // gradient (setTrialF), so the material MUST be a FiniteStrainND2DMaterial
  // (e.g. LogStrain2D). Reject here — before the element exists — so a misuse
  // can never reach the assembly with a non-finite material.
  if (geom == LadrunoCST::Geom::FINITE &&
      dynamic_cast<FiniteStrainND2DMaterial *>(mat) == 0) {
    opserr << "WARNING LadrunoCST -- '-geom finite' needs a FiniteStrainND2DMaterial "
              "(e.g. LogStrain2D); material " << matTag << " is not one\n";
    return 0;
  }

  // Ladruno (ADR 70): the explicit bulk-viscosity force block lives only in the
  // small-strain resisting-force path; strip it (with a diagnostic) under -geom
  // finite so Print stays honest. Mirrors LadrunoQuad.
  if ((bvB1 > 0.0 || bvB2 > 0.0) && geom == LadrunoCST::Geom::FINITE) {
    opserr << "WARNING LadrunoCST " << idata[0]
           << " -- -bulkViscosity is not applied with -geom finite; ignoring it\n";
    bvB1 = bvB2 = 0.0;
  }

  return new LadrunoCST(idata[0], idata[1], idata[2], idata[3],
                        *mat, typeBuf, thk, geom, rho, b1, b2, pressure,
                        bvB1, bvB2);   // Ladruno (ADR 70) geom; (W2-E1) bulk-viscosity
}
