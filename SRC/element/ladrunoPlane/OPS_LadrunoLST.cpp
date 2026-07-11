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

// Parser/factory for LadrunoLST (Ladruno fork, ADR 70 P3). Flag-based:
//   element LadrunoLST tag n1 n2 n3 n4 n5 n6 matTag
//       <-type PlaneStrain|PlaneStress> <-geom linear|finite>
//       <-thick t> <-rho r> <-body b1 b2> <-pressure p>
// Node order matches upstream SixNodeTri: corners n1,n2,n3 (CCW), midsides
// n4=(1-2), n5=(2-3), n6=(3-1). std ONLY: '-formulation bbar' is refused —
// constant element-mean dilatation is rank-deficient on the T6 (two conformal
// zero-energy modes; refuted at P3, see ADR 70 §9).

#include <LadrunoLST.h>
#include <NDMaterial.h>
#include <FiniteStrainND2DMaterial.h>   // Ladruno (ADR 70): -geom finite material check
#include <elementAPI.h>
#include <string.h>

void *OPS_LadrunoLST()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 2 || ndf != 2) {
    opserr << "WARNING LadrunoLST -- model ndm/ndf not 2/2\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LadrunoLST tag n1 n2 n3 n4 n5 n6 matTag "
              "<-type PlaneStrain|PlaneStress> "
              "<-geom linear|finite> <-thick t> <-rho r> <-body b1 b2> <-pressure p>\n";
    return 0;
  }

  int idata[7];
  int num = 7;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr << "WARNING LadrunoLST -- invalid tag/node integers\n";
    return 0;
  }

  int matTag;
  num = 1;
  if (OPS_GetIntInput(&num, &matTag) < 0) {
    opserr << "WARNING LadrunoLST -- invalid matTag\n";
    return 0;
  }
  NDMaterial *mat = OPS_getNDMaterial(matTag);
  if (mat == 0) {
    opserr << "WARNING LadrunoLST -- material " << matTag << " not found\n";
    return 0;
  }

  double thk = 1.0, rho = 0.0, b1 = 0.0, b2 = 0.0, pressure = 0.0;
  double bvB1 = 0.0, bvB2 = 0.0;   // Ladruno (W2-E1): explicit bulk-viscosity coeffs (off by default)
  char typeBuf[24];
  strcpy(typeBuf, "PlaneStrain");
  LadrunoLST::Formulation form = LadrunoLST::Formulation::STD;
  LadrunoLST::Geom geom = LadrunoLST::Geom::LINEAR;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &thk) < 0) { opserr << "WARNING LadrunoLST -- bad -thick\n"; return 0; }

    } else if (strcmp(opt, "-type") == 0) {
      const char *t = OPS_GetString();
      if (strcmp(t, "PlaneStrain") == 0 || strcmp(t, "PlaneStrain2D") == 0)
        strcpy(typeBuf, "PlaneStrain");
      else if (strcmp(t, "PlaneStress") == 0 || strcmp(t, "PlaneStress2D") == 0)
        strcpy(typeBuf, "PlaneStress");
      else { opserr << "WARNING LadrunoLST -- unknown -type " << t << "\n"; return 0; }

    } else if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      const char *f = OPS_GetString();
      if (strcmp(f, "std") == 0)       form = LadrunoLST::Formulation::STD;
      else if (strcmp(f, "bbar") == 0) {
        // Ladruno (ADR 70 P3): REFUTED. Constant element-mean dilatation
        // (B-bar / centroid F-bar) is rank-deficient on the T6: the two
        // quadratic conformal modes (Re/Im z^2) have zero deviatoric strain
        // and zero MEAN dilatation, so a free element gets 5 zero-energy
        // modes (3 RBM + 2 spurious). No uncontrolled mechanism ships; the
        // T6 volumetric cure is ADR 70 P4 (F-bar-Patch / projected dilatation).
        opserr << "WARNING LadrunoLST -- '-formulation bbar' is NOT available: "
                  "constant mean-dilatation B-bar/F-bar is rank-deficient on "
                  "the T6 (2 conformal zero-energy modes, ADR 70 P3); use "
                  "LadrunoQuad bbar or wait for the nodal F-bar-Patch (P4)\n";
        return 0;
      }
      else { opserr << "WARNING LadrunoLST -- unknown -formulation " << f
                    << " (want std)\n"; return 0; }

    } else if (strcmp(opt, "-geom") == 0) {
      const char *gopt = OPS_GetString();
      if (strcmp(gopt, "linear") == 0)      geom = LadrunoLST::Geom::LINEAR;
      else if (strcmp(gopt, "finite") == 0) geom = LadrunoLST::Geom::FINITE;
      else { opserr << "WARNING LadrunoLST -- unknown -geom " << gopt
                    << " (want linear|finite)\n"; return 0; }

    } else if (strcmp(opt, "-rho") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0) { opserr << "WARNING LadrunoLST -- bad -rho\n"; return 0; }

    } else if (strcmp(opt, "-body") == 0) {
      double bb[2];
      num = 2;
      if (OPS_GetDoubleInput(&num, bb) < 0) { opserr << "WARNING LadrunoLST -- bad -body\n"; return 0; }
      b1 = bb[0]; b2 = bb[1];

    } else if (strcmp(opt, "-pressure") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &pressure) < 0) { opserr << "WARNING LadrunoLST -- bad -pressure\n"; return 0; }

    } else if (strcmp(opt, "-bulkViscosity") == 0 || strcmp(opt, "-bv") == 0) {
      // Ladruno (W2-E1): explicit bulk viscosity reads TWO doubles (linear b1,
      // quadratic b2). Both must be >= 0; a negative coeff is warned and ignored.
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING LadrunoLST -- -bulkViscosity needs two values (b1 b2)\n";
        return 0;
      }
      double bv[2] = { 0.0, 0.0 };
      num = 2;
      if (OPS_GetDoubleInput(&num, bv) < 0) { opserr << "WARNING LadrunoLST -- bad -bulkViscosity\n"; return 0; }
      if (bv[0] < 0.0) { opserr << "WARNING LadrunoLST -- -bulkViscosity b1 < 0; ignoring (b1=0)\n"; bv[0] = 0.0; }
      if (bv[1] < 0.0) { opserr << "WARNING LadrunoLST -- -bulkViscosity b2 < 0; ignoring (b2=0)\n"; bv[1] = 0.0; }
      bvB1 = bv[0];
      bvB2 = bv[1];

    } else {
      opserr << "WARNING LadrunoLST -- ignoring unknown option " << opt << "\n";
    }
  }

  // Ladruno (ADR 70): -geom finite is plane-strain only. formFinite's current
  // volume weight dv = J·detJ0·thickness·w holds the thickness fixed — it omits
  // the out-of-plane stretch λ = F_33 that a finite plane-stress state develops.
  // Refuse here — before the element exists — rather than assemble a silently
  // wrong finite plane-stress residual/tangent. Mirrors LadrunoQuad/CST.
  if (geom == LadrunoLST::Geom::FINITE && strcmp(typeBuf, "PlaneStress") == 0) {
    opserr << "WARNING LadrunoLST -- '-geom finite' is PlaneStrain only "
              "(finite plane-stress omits the thickness stretch lambda in the "
              "volume weight; use -type PlaneStrain, ADR 70)\n";
    return 0;
  }

  // Ladruno (ADR 70): -geom finite drives the material by the deformation
  // gradient (setTrialF), so the material MUST be a FiniteStrainND2DMaterial
  // (e.g. LogStrain2D). Reject here — before the element exists — so a misuse
  // can never reach the assembly with a non-finite material.
  if (geom == LadrunoLST::Geom::FINITE &&
      dynamic_cast<FiniteStrainND2DMaterial *>(mat) == 0) {
    opserr << "WARNING LadrunoLST -- '-geom finite' needs a FiniteStrainND2DMaterial "
              "(e.g. LogStrain2D); material " << matTag << " is not one\n";
    return 0;
  }

  // Ladruno (ADR 70): the explicit bulk-viscosity force block lives only in the
  // small-strain resisting-force path; strip it (with a diagnostic) under -geom
  // finite so Print stays honest. Mirrors LadrunoQuad/CST.
  if ((bvB1 > 0.0 || bvB2 > 0.0) && geom == LadrunoLST::Geom::FINITE) {
    opserr << "WARNING LadrunoLST " << idata[0]
           << " -- -bulkViscosity is not applied with -geom finite; ignoring it\n";
    bvB1 = bvB2 = 0.0;
  }

  return new LadrunoLST(idata[0], idata[1], idata[2], idata[3],
                        idata[4], idata[5], idata[6],
                        *mat, typeBuf, thk, form, geom, rho, b1, b2, pressure,
                        bvB1, bvB2);   // Ladruno (ADR 70) form/geom; (W2-E1) bulk-viscosity
}
