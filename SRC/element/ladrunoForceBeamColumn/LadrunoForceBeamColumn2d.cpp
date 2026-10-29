/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
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

// Standalone Ladruno force-based (flexibility) 2D fiber frame element
// (Ladruno WP-200) combining exact force-based equilibrium state
// determination with Tier-1 per-IP/element/user characteristic length (lch),
// Tier-2 single/dual cohesive fracture hinges with 2x2 Schur-complement
// static condensation, P-M interaction, and IMPLEX support.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <cmath>
#include <vector>

#include <Information.h>
#include <Parameter.h>
#include <LadrunoForceBeamColumn2d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <OPS_Globals.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <LadrunoResponseTokens.h>
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <LadrunoCohesiveHinge.h>
#include <Bilin.h>
#include <ModIMKPeakOriented.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <CompositeSimpsonBeamIntegration.h>

extern thread_local Element *ops_TheActiveElement;

thread_local Matrix LadrunoForceBeamColumn2d::theMatrix(6, 6);
thread_local Vector LadrunoForceBeamColumn2d::theVector(6);
thread_local double LadrunoForceBeamColumn2d::workArea[200];

thread_local Vector LadrunoForceBeamColumn2d::vsSubdivide[maxNumSections];
thread_local Vector LadrunoForceBeamColumn2d::SsrSubdivide[maxNumSections];
thread_local Matrix LadrunoForceBeamColumn2d::fsSubdivide[maxNumSections];

static bool lfbc2dTryOptionalDouble(double &out)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return false;
  char tok[64];
  const char *s = OPS_GetStringFromAll(tok, sizeof(tok));
  OPS_ResetCurrentInputArg(-1);
  if (s == 0 || s[0] == '\0') return false;
  char *endp = 0;
  (void)strtod(s, &endp);
  if (endp == s || *endp != '\0') return false;
  int one = 1;
  return (OPS_GetDoubleInput(&one, &out) == 0);
}

static bool lfbc2dTryOptionalInt(int &out)
{
  if (OPS_GetNumRemainingInputArgs() < 1) return false;
  char tok[64];
  const char *s = OPS_GetStringFromAll(tok, sizeof(tok));
  OPS_ResetCurrentInputArg(-1);
  if (s == 0 || s[0] == '\0') return false;
  char *endp = 0;
  long v = strtol(s, &endp, 10);
  if (endp == s || *endp != '\0') return false;
  int one = 1;
  if (OPS_GetIntInput(&one, &out) < 0) return false;
  (void)v;
  return true;
}

static bool lfbc2dGetHingeImplexError(UniaxialMaterial *mat, double &errOut)
{
  if (mat == 0) return false;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(mat)) {
    errOut = ch->getImplexError();
    return true;
  }
  if (Bilin *b = dynamic_cast<Bilin*>(mat)) {
    errOut = b->getImplexError();
    return true;
  }
  if (ModIMKPeakOriented *m = dynamic_cast<ModIMKPeakOriented*>(mat)) {
    errOut = m->getImplexError();
    return true;
  }
  return false;
}

static bool lfbc2dHingeIsImplex(UniaxialMaterial *mat)
{
  if (mat == 0) return false;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(mat))
    return ch->isImplex();
  if (Bilin *b = dynamic_cast<Bilin*>(mat))
    return b->isImplex();
  if (ModIMKPeakOriented *m = dynamic_cast<ModIMKPeakOriented*>(mat))
    return m->isImplex();
  return false;
}

static double lfbc2dScalePM(double N, double P0, double alpha)
{
  if (P0 <= 0.0) return 1.0;
  double n = N / P0;
  double f = 1.0 - alpha * n - (1.0 - alpha) * n * n;
  if (f < 0.05) f = 0.05;
  else if (f > 2.5) f = 2.5;
  return f;
}

static double lfbc2dScalePMDeriv(double N, double P0, double alpha)
{
  if (P0 <= 0.0) return 0.0;
  double n = N / P0;
  double f = 1.0 - alpha * n - (1.0 - alpha) * n * n;
  if (f <= 0.05 || f >= 2.5) return 0.0;
  return (-alpha - 2.0 * (1.0 - alpha) * n) / P0;
}

static double lfbc2dEvalPM(double N, double P0, double alpha,
                           double P0t, double Pbal, double muBal,
                           double &dfdN)
{
  dfdN = 0.0;
  if (P0 <= 0.0) return 1.0;
  if (Pbal > 0.0 && P0 > Pbal && P0t > 0.0 && muBal >= 1.0) {
    double P = -N; // compression positive
    double f = 1.0;
    double dfdP = 0.0;
    if (P <= 0.0) {
      double u = 1.0 + P / P0t;
      if (u <= 0.05) {
        return 0.05;
      }
      f = u;
      dfdP = 1.0 / P0t;
    } else if (Pbal > 0.0 && P <= Pbal) {
      double r = P / Pbal;
      f = 1.0 + (muBal - 1.0) * r * (2.0 - r);
      dfdP = (muBal - 1.0) * (2.0 - 2.0 * r) / Pbal;
    } else {
      double denom = P0 - Pbal;
      if (denom <= 1.0e-12 * P0) {
        return 0.05;
      }
      double r = (P - Pbal) / denom;
      f = muBal * (1.0 - r * r);
      dfdP = -2.0 * muBal * r / denom;
    }
    if (f <= 0.05) {
      dfdN = 0.0;
      return 0.05;
    }
    if (f >= 3.0) {
      dfdN = 0.0;
      return 3.0;
    }
    dfdN = -dfdP;
    return f;
  }
  dfdN = lfbc2dScalePMDeriv(N, P0, alpha);
  return lfbc2dScalePM(N, P0, alpha);
}

static void lfbc2dInvertSym2x2Guarded(double A, double B, double C, double (&Hinv)[2][2])
{
  double avg = 0.5 * (A + B);
  double hdiff = 0.5 * (A - B);
  double disc = sqrt(hdiff * hdiff + C * C);
  double lam1 = avg + disc;
  double lam2 = avg - disc;
  double v1z, v1y, v2z, v2y;
  if (fabs(C) > 1.0e-30 * (fabs(hdiff) + 1.0e-300)) {
    double dx = lam1 - B, dy = C;
    double n = sqrt(dx * dx + dy * dy);
    v1z = dx / n; v1y = dy / n;
  } else if (A >= B) {
    v1z = 1.0; v1y = 0.0;
  } else {
    v1z = 0.0; v1y = 1.0;
  }
  v2z = -v1y; v2y = v1z;

  double scale = fabs(A) + fabs(B) + 2.0 * fabs(C);
  double floor = 1.0e-8 * scale + 1.0e-300;
  double l1s = (fabs(lam1) < floor) ? ((lam1 < 0.0) ? -floor : floor) : lam1;
  double l2s = (fabs(lam2) < floor) ? ((lam2 < 0.0) ? -floor : floor) : lam2;
  double inv1 = 1.0 / l1s;
  double inv2 = 1.0 / l2s;
  Hinv[0][0] = v1z * v1z * inv1 + v2z * v2z * inv2;
  Hinv[1][1] = v1y * v1y * inv1 + v2y * v2y * inv2;
  Hinv[0][1] = v1z * v1y * inv1 + v2z * v2y * inv2;
  Hinv[1][0] = Hinv[0][1];
}

static void lfbc2dInvertGen2x2Guarded(double H00, double H01, double H10, double H11,
                                      double (&Hinv)[2][2])
{
  double scale = fabs(H00) + fabs(H11) + fabs(H01) + fabs(H10);
  double det = H00 * H11 - H01 * H10;
  double floor = (1.0e-8 * scale) * (1.0e-8 * scale) + 1.0e-300;
  if (fabs(det) >= floor) {
    double invDet = 1.0 / det;
    Hinv[0][0] =  H11 * invDet;
    Hinv[0][1] = -H01 * invDet;
    Hinv[1][0] = -H10 * invDet;
    Hinv[1][1] =  H00 * invDet;
  } else {
    lfbc2dInvertSym2x2Guarded(H00, H11, 0.5 * (H01 + H10), Hinv);
  }
}

void* OPS_LadrunoForceBeamColumn2d()
{
  if (OPS_GetNDM() != 2 || OPS_GetNDF() != 3) {
    opserr << "WARNING LadrunoForceBeamColumn2d -- model must be ndm 2/ndf 3\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING LadrunoForceBeamColumn2d: insufficient arguments; want: "
           << "eleTag iNode jNode transfTag beamIntTag ...\n";
    return 0;
  }

  int iData[5];
  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) < 0) {
    opserr << "WARNING LadrunoForceBeamColumn2d: invalid int inputs (eleTag iNode jNode)\n";
    return 0;
  }

  // Distinguish modern form (transfTag beamIntTag) from legacy inline form
  // (numSecs secTag1 ... secTagN transfTag):
  // In modern form, there are exactly 2 integer args before the end of args or the first '-flag'.
  bool legacyForm = false;
  if (OPS_GetNumRemainingInputArgs() >= 3) {
    int saveArgs = OPS_GetNumRemainingInputArgs();
    char tok1[64], tok2[64], tok3[64];
    (void)OPS_GetStringFromAll(tok1, sizeof(tok1));
    (void)OPS_GetStringFromAll(tok2, sizeof(tok2));
    const char *s3 = OPS_GetStringFromAll(tok3, sizeof(tok3));
    int consumed = saveArgs - OPS_GetNumRemainingInputArgs();
    OPS_ResetCurrentInputArg(-consumed);
    if (s3 != 0 && s3[0] != '-') {
      char *endp = 0;
      (void)strtol(s3, &endp, 10);
      if (endp != s3 && *endp == '\0') {
        legacyForm = true;
      }
    }
  }

  int transfTag = 0;
  int beamIntTag = 0;
  std::vector<int> legacySecTags;
  BeamIntegration *legacyBi = 0;

  if (!legacyForm) {
    numData = 2;
    if (OPS_GetIntInput(&numData, &iData[3]) < 0) {
      opserr << "WARNING LadrunoForceBeamColumn2d: invalid transfTag or beamIntTag\n";
      return 0;
    }
    transfTag = iData[3];
    beamIntTag = iData[4];
  } else {
    int numSec = 0;
    numData = 1;
    if (OPS_GetIntInput(&numData, &numSec) < 0 || numSec <= 0 || numSec > 20) {
      opserr << "WARNING LadrunoForceBeamColumn2d: invalid numSecs (" << numSec << ")\n";
      return 0;
    }
    if (OPS_GetNumRemainingInputArgs() < numSec + 1) {
      opserr << "WARNING LadrunoForceBeamColumn2d: insufficient section tags or transfTag\n";
      return 0;
    }
    legacySecTags.resize(numSec);
    numData = numSec;
    if (OPS_GetIntInput(&numData, legacySecTags.data()) < 0) {
      opserr << "WARNING LadrunoForceBeamColumn2d: invalid section tags\n";
      return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, &transfTag) < 0) {
      opserr << "WARNING LadrunoForceBeamColumn2d: invalid transfTag\n";
      return 0;
    }
  }

  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];

  double mass = 0.0;
  int cmass = 0;
  int maxIters = 10;
  double tol = 1.0e-12;
  int maxSubdivisions = 4;
  double subdivideFactor = 10.0;
  int lchMode = 0;
  double userLch = 0.0;
  Damping *theDamping = 0;
  UniaxialMaterial *hingeMat = 0;
  UniaxialMaterial *hingeMatJ = 0;
  int dualHingeMode = 0;
  double hingeXiI = 0.0;
  double hingeXiJ = 1.0;
  bool eleImplex = false;
  double elePmP0 = 0.0;
  double elePmAlpha = 0.5;
  double elePmP0t = 0.0;
  double elePmPbal = 0.0;
  double elePmMuBal = 1.0;
  double gavy = 0.0;
  bool axialRelease = false;
  int spanLoadMode = 1;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (type == 0) break;

    if (strcmp(type, "-cMass") == 0 || strcmp(type, "-cmass") == 0) {
      cmass = 1;
    }
    else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -mass needs a value\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &mass) < 0 || !std::isfinite(mass) || mass < 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid -mass value\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-damp") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -damp needs a damping tag\n";
        delete legacyBi;
        return 0;
      }
      int dTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &dTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid damping tag\n";
        delete legacyBi;
        return 0;
      }
      theDamping = OPS_getDamping(dTag);
      if (theDamping == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: damping " << dTag << " not found\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -iter needs maxIters and tol\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &maxIters) < 0 || maxIters <= 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid maxIters\n";
        delete legacyBi;
        return 0;
      }
      if (OPS_GetDoubleInput(&numData, &tol) < 0 || !std::isfinite(tol) || tol <= 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid tol\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-subdivide") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -subdivide needs maxSubdivisions [subdivideFactor]\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &maxSubdivisions) < 0 || maxSubdivisions <= 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid maxSubdivisions\n";
        delete legacyBi;
        return 0;
      }
      subdivideFactor = 10.0;
      double optSubFac = subdivideFactor;
      if (lfbc2dTryOptionalDouble(optSubFac)) {
        if (!std::isfinite(optSubFac) || optSubFac <= 1.0) {
          opserr << "WARNING LadrunoForceBeamColumn2d: invalid subdivideFactor (must be > 1.0)\n";
          delete legacyBi;
          return 0;
        }
        subdivideFactor = optSubFac;
      }
    }
    else if (strcmp(type, "-totalElementLch") == 0) {
      lchMode = 1;
    }
    else if (strcmp(type, "-lch") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -lch needs {ip|element|<value>}\n";
        delete legacyBi;
        return 0;
      }
      char lchBuf[64];
      const char *arg = OPS_GetStringFromAll(lchBuf, sizeof(lchBuf));
      if (strcmp(arg, "ip") == 0) {
        lchMode = 0;
      } else if (strcmp(arg, "element") == 0) {
        lchMode = 1;
      } else {
        char *endp = 0;
        double val = strtod(arg, &endp);
        if (endp == arg || *endp != '\0' || !std::isfinite(val) || val <= 0.0) {
          opserr << "WARNING LadrunoForceBeamColumn2d: invalid -lch '" << arg
                 << "' (expected 'ip', 'element', or a positive number)\n";
          delete legacyBi;
          return 0;
        }
        lchMode = 2;
        userLch = val;
      }
    }
    else if (strcmp(type, "-Integration") == 0 || strcmp(type, "-integration") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -Integration needs an integration type\n";
        delete legacyBi;
        return 0;
      }
      const char *integType = OPS_GetString();
      delete legacyBi;
      legacyBi = 0;
      if (strcmp(integType, "Lobatto") == 0)
        legacyBi = new LobattoBeamIntegration();
      else if (strcmp(integType, "Legendre") == 0)
        legacyBi = new LegendreBeamIntegration();
      else if (strcmp(integType, "Radau") == 0)
        legacyBi = new RadauBeamIntegration();
      else if (strcmp(integType, "NewtonCotes") == 0)
        legacyBi = new NewtonCotesBeamIntegration();
      else if (strcmp(integType, "Trapezoidal") == 0)
        legacyBi = new TrapezoidalBeamIntegration();
      else if (strcmp(integType, "CompositeSimpson") == 0)
        legacyBi = new CompositeSimpsonBeamIntegration();
      else {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid -Integration type '" << integType << "'\n";
        return 0;
      }
    }
    else if (strcmp(type, "-hinge") == 0 || strcmp(type, "-hingeZ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hinge needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid -hinge material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMat = OPS_getUniaxialMaterial(hTag);
      if (hingeMat == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hinge uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      int optTagJ = 0;
      if (lfbc2dTryOptionalInt(optTagJ)) {
        hingeMatJ = OPS_getUniaxialMaterial(optTagJ);
        if (hingeMatJ == 0) {
          opserr << "WARNING LadrunoForceBeamColumn2d: -hinge uniaxialMaterial " << optTagJ << " not found\n";
          delete legacyBi;
          return 0;
        }
        dualHingeMode = 1;
      }
    }
    else if (strcmp(type, "-hingeI") == 0 || strcmp(type, "-hinge1") == 0 || strcmp(type, "-hingeZI") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeI needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid -hingeI material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMat = OPS_getUniaxialMaterial(hTag);
      if (hingeMat == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeI uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeJ") == 0 || strcmp(type, "-hinge2") == 0 || strcmp(type, "-hingeZJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeJ needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid -hingeJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatJ = OPS_getUniaxialMaterial(hTag);
      if (hingeMatJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeJ uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeIJ") == 0 || strcmp(type, "-hingeZIJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeIJ needs at least one uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0, hTagJ = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: invalid -hingeIJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hTagJ = hTagI;
      int optTagJ = 0;
      if (lfbc2dTryOptionalInt(optTagJ))
        hTagJ = optTagJ;
      hingeMat = OPS_getUniaxialMaterial(hTagI);
      hingeMatJ = OPS_getUniaxialMaterial(hTagJ);
      if (hingeMat == 0 || hingeMatJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeIJ uniaxialMaterial not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeLoc") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeLoc needs xiI [xiJ]\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &hingeXiI) < 0 || !std::isfinite(hingeXiI) ||
          hingeXiI < 0.0 || hingeXiI > 1.0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -hingeLoc xiI must be in [0, 1]\n";
        delete legacyBi;
        return 0;
      }
      double optXiJ = 0.0;
      if (lfbc2dTryOptionalDouble(optXiJ)) {
        if (!std::isfinite(optXiJ) || optXiJ < 0.0 || optXiJ > 1.0) {
          opserr << "WARNING LadrunoForceBeamColumn2d: -hingeLoc xiJ must be in [0, 1]\n";
          delete legacyBi;
          return 0;
        }
        hingeXiJ = optXiJ;
      } else {
        hingeXiJ = 1.0 - hingeXiI;
      }
    }
    else if (strcmp(type, "-implex") == 0 || strcmp(type, "-IMPLEX") == 0) {
      eleImplex = true;
    }
    else if (strcmp(type, "-pm") == 0 || strcmp(type, "-PM") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          OPS_GetDoubleInput(&numData, &elePmP0) < 0 || !std::isfinite(elePmP0) || elePmP0 <= 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -pm needs P0 > 0\n";
        delete legacyBi;
        return 0;
      }
      double optAlpha = elePmAlpha;
      if (lfbc2dTryOptionalDouble(optAlpha)) {
        if (!std::isfinite(optAlpha)) {
          opserr << "WARNING LadrunoForceBeamColumn2d: invalid -pm alpha\n";
          delete legacyBi;
          return 0;
        }
        elePmAlpha = optAlpha;
      }
      elePmP0t = 0.0;
      elePmPbal = 0.0;
      elePmMuBal = 1.0;
    }
    else if (strcmp(type, "-pmBal") == 0 || strcmp(type, "-PMBal") == 0) {
      double pb[4] = {0.0, 0.0, 0.0, 1.0};
      numData = 4;
      if (OPS_GetNumRemainingInputArgs() < 4 ||
          OPS_GetDoubleInput(&numData, pb) < 0 ||
          !std::isfinite(pb[0]) || !std::isfinite(pb[1]) ||
          !std::isfinite(pb[2]) || !std::isfinite(pb[3]) ||
          pb[0] <= 0.0 || pb[1] <= 0.0 || pb[2] < 0.0 || pb[2] >= pb[0] || pb[3] < 1.0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -pmBal needs P0c>0 P0t>0 0<=Pbal<P0c muBal>=1\n";
        delete legacyBi;
        return 0;
      }
      elePmP0 = pb[0];
      elePmP0t = pb[1];
      elePmPbal = pb[2];
      elePmMuBal = pb[3];
    }
    else if (strcmp(type, "-shear") == 0 || strcmp(type, "-shearY") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          OPS_GetDoubleInput(&numData, &gavy) < 0 || !std::isfinite(gavy) || gavy <= 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: -shear needs GAvy > 0\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-axialRelease") == 0 || strcmp(type, "-decoupleAxial") == 0) {
      axialRelease = true;
    }
    else if (strcmp(type, "-spanLoad") == 0 || strcmp(type, "-spanLoadEnhancement") == 0) {
      int optFlag = 1;
      if (lfbc2dTryOptionalInt(optFlag))
        spanLoadMode = (optFlag != 0) ? 1 : 0;
      else
        spanLoadMode = 1;
    }
    else if (strcmp(type, "-noSpanLoad") == 0) {
      spanLoadMode = 0;
    }
    else {
      opserr << "WARNING LadrunoForceBeamColumn2d: unknown option '" << type << "'\n";
      delete legacyBi;
      return 0;
    }
  }

  if (hingeMatJ != 0)
    dualHingeMode = 1;

  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING LadrunoForceBeamColumn2d: coord transformation " << transfTag << " not found\n";
    delete legacyBi;
    return 0;
  }

  BeamIntegration *bi = 0;
  SectionForceDeformation **sections = 0;
  int nSecs = 0;

  if (!legacyForm) {
    BeamIntegrationRule *theRule = OPS_getBeamIntegrationRule(beamIntTag);
    if (theRule == 0) {
      opserr << "WARNING LadrunoForceBeamColumn2d: beam integration " << beamIntTag << " not found\n";
      return 0;
    }
    bi = theRule->getBeamIntegration();
    if (bi == 0) {
      opserr << "WARNING LadrunoForceBeamColumn2d: beam integration is null\n";
      return 0;
    }
    const ID &secTags = theRule->getSectionTags();
    nSecs = secTags.Size();
    if (nSecs <= 0 || nSecs > 20) {
      opserr << "WARNING LadrunoForceBeamColumn2d: number of sections (" << nSecs << ") must be in 1..20\n";
      return 0;
    }
    sections = new SectionForceDeformation *[nSecs];
    for (int i = 0; i < nSecs; i++) {
      sections[i] = OPS_getSectionForceDeformation(secTags(i));
      if (sections[i] == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: section " << secTags(i) << " not found\n";
        delete [] sections;
        return 0;
      }
    }
  } else {
    if (legacyBi == 0)
      legacyBi = new LobattoBeamIntegration();
    bi = legacyBi;
    nSecs = (int)legacySecTags.size();
    sections = new SectionForceDeformation *[nSecs];
    for (int i = 0; i < nSecs; i++) {
      sections[i] = OPS_getSectionForceDeformation(legacySecTags[i]);
      if (sections[i] == 0) {
        opserr << "WARNING LadrunoForceBeamColumn2d: section " << legacySecTags[i] << " not found\n";
        delete [] sections;
        delete legacyBi;
        return 0;
      }
    }
  }

  UniaxialMaterial *hI_tmp = 0, *hJ_tmp = 0;
  if (eleImplex || elePmP0 > 0.0) {
    auto cfgHinge = [&](UniaxialMaterial *&hPtr, UniaxialMaterial *&tmpPtr) {
      if (hPtr != 0) {
        tmpPtr = hPtr->getCopy();
        if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(tmpPtr)) {
          if (eleImplex) ch->setImplex(true);
          if (elePmP0 > 0.0) {
            if (elePmP0t > 0.0 && elePmMuBal >= 1.0)
              ch->setPMBal(elePmP0, elePmP0t, elePmPbal, elePmMuBal);
            else
              ch->setPM(elePmP0, elePmAlpha);
          }
        } else if (Bilin *b = dynamic_cast<Bilin*>(tmpPtr)) {
          if (eleImplex) b->setImplex(true);
        } else if (ModIMKPeakOriented *m = dynamic_cast<ModIMKPeakOriented*>(tmpPtr)) {
          if (eleImplex) m->setImplex(true);
        }
        hPtr = tmpPtr;
      }
    };
    cfgHinge(hingeMat, hI_tmp);
    cfgHinge(hingeMatJ, hJ_tmp);
  }

  Element *theEle = new LadrunoForceBeamColumn2d(
      eleTag, nodeI, nodeJ, nSecs, sections, *bi, *theTransf,
      mass, maxIters, tol, cmass, theDamping, lchMode, userLch,
      hingeMat, hingeMatJ, dualHingeMode, hingeXiI, hingeXiJ,
      elePmP0, elePmAlpha, maxSubdivisions, subdivideFactor, eleImplex,
      gavy, elePmP0t, elePmPbal, elePmMuBal, axialRelease, spanLoadMode);

  delete hI_tmp;
  delete hJ_tmp;
  delete [] sections;
  delete legacyBi;
  return theEle;
}

LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d(
    int tag, int nodeI, int nodeJ,
    int numSec, SectionForceDeformation **sec,
    BeamIntegration &bi,
    CrdTransf &coordTransf, double massDensPerUnitLength,
    int maxNumIters, double tolerance,
    int cm, Damping *damping,
    int lchM, double userL,
    UniaxialMaterial *hingeMat,
    UniaxialMaterial *hingeMatJ,
    int dualMode,
    double xiI, double xiJ,
    double pmP0_, double pmAlpha_,
    int maxSub, double subFac,
    bool eleImplex_,
    double gavy_,
    double pmP0t_, double pmPbal_, double pmMuBal_,
    bool axialRelease_, int spanLoadMode_)
  : Element(tag, ELE_TAG_LadrunoForceBeamColumn2d),
    numSections(0), sections(0), crdTransf(0), beamIntegr(0),
    connectedExternalNodes(2), load(NEGD),
    rho(massDensPerUnitLength), cMass(cm),
    maxIters(maxNumIters), tol(tolerance),
    initialFlag(0),
    kv(NEBD, NEBD), Se(NEBD),
    kvcommit(NEBD, NEBD), Secommit(NEBD),
    kvbulk(NEBD, NEBD), kvbulkCommit(NEBD, NEBD),
    fs(0), vs(0), Ssr(0), vscommit(0), fscommit(0), Ssrcommit(0),
    numEleLoads(0), eleLoadDirty(false), Ki(0), parameterID(0), theDamping(0),
    maxSubdivisions(maxSub), subdivideFactor(subFac),
    current_section_lch(0.0), lchMode(lchM), userLch(userL), GAvy(gavy_),
    hingeOn(0), dualHingeMode(dualMode), hingeXiI(xiI), hingeXiJ(xiJ),
    theHinge(0), theHingeJ(0),
    hingeJump(0.0), hingeJumpCommit(0.0),
    hingeJumpJ(0.0), hingeJumpJCommit(0.0),
    hingeMscale(0.0), hingeMscaleCommit(0.0),
    pmP0(pmP0_), pmAlpha(pmAlpha_),
    pmP0t(pmP0t_), pmPbal(pmPbal_), pmMuBal(pmMuBal_),
    axialRelease(axialRelease_), spanLoadMode(spanLoadMode_),
    useImplex(eleImplex_), eleImplexErr(0.0), eleImplexErrCommit(0.0),
    postRestoreSync(false), vRestoreValid(false),
    vCommitBulk(NEBD), vTrialBulk(NEBD), vTrialBulkCommit(NEBD),
    vSolved(NEBD), solvedValid(false),
    vSolvedCommit(NEBD), solvedCommitValid(false)
{
  load.Zero();
  vRestore[0] = vRestore[1] = vRestore[2] = 0.0;
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) {
      hingeHinv[m][n] = 0.0;
      hingeHinvCommit[m][n] = 0.0;
    }
    for (int a = 0; a < 3; a++) {
      hingeKva[a][m] = 0.0;
      hingeKvaCommit[a][m] = 0.0;
      hingeKav[m][a] = 0.0;
      hingeKavCommit[m][a] = 0.0;
    }
  }

  theNodes[0] = 0;
  theNodes[1] = 0;
  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;

  if (hingeMat != 0) {
    theHinge = hingeMat->getCopy();
    if (theHinge == 0) {
      opserr << "LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d - failed to copy hinge material\n";
      exit(-1);
    }
    hingeOn = 1;
  }
  if (hingeMatJ != 0) {
    theHingeJ = hingeMatJ->getCopy();
    if (theHingeJ == 0) {
      opserr << "LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d - failed to copy hingeJ material\n";
      exit(-1);
    }
    hingeOn = 1;
    dualHingeMode = 1;
  }

  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d - could not create copy of beam integration object\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy2d();
  if (crdTransf == 0) {
    opserr << "Error: LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d - could not create copy of coordinate transformation object\n";
    exit(-1);
  }

  if (damping) {
    theDamping = (*damping).getCopy();
    if (!theDamping) {
      opserr << "Error: LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d - could not create copy of damping object\n";
      exit(-1);
    }
  }

  this->setSectionPointers(numSec, sec);

  for (int i = 0; i < maxNumEleLoads; i++) {
    eleLoads[i] = 0;
    eleLoadFactors[i] = 0.0;
  }
}

LadrunoForceBeamColumn2d::LadrunoForceBeamColumn2d()
  : Element(0, ELE_TAG_LadrunoForceBeamColumn2d),
    numSections(0), sections(0), crdTransf(0), beamIntegr(0),
    connectedExternalNodes(2), load(NEGD),
    rho(0.0), cMass(0),
    maxIters(10), tol(1.0e-12),
    initialFlag(0),
    kv(NEBD, NEBD), Se(NEBD),
    kvcommit(NEBD, NEBD), Secommit(NEBD),
    kvbulk(NEBD, NEBD), kvbulkCommit(NEBD, NEBD),
    fs(0), vs(0), Ssr(0), vscommit(0), fscommit(0), Ssrcommit(0),
    numEleLoads(0), eleLoadDirty(false), Ki(0), parameterID(0), theDamping(0),
    maxSubdivisions(4), subdivideFactor(10.0),
    current_section_lch(0.0), lchMode(0), userLch(0.0), GAvy(0.0),
    hingeOn(0), dualHingeMode(0), hingeXiI(0.0), hingeXiJ(1.0),
    theHinge(0), theHingeJ(0),
    hingeJump(0.0), hingeJumpCommit(0.0),
    hingeJumpJ(0.0), hingeJumpJCommit(0.0),
    hingeMscale(0.0), hingeMscaleCommit(0.0),
    pmP0(0.0), pmAlpha(0.5),
    pmP0t(0.0), pmPbal(0.0), pmMuBal(1.0),
    axialRelease(false), spanLoadMode(1),
    useImplex(false), eleImplexErr(0.0), eleImplexErrCommit(0.0),
    postRestoreSync(false), vRestoreValid(false),
    vCommitBulk(NEBD), vTrialBulk(NEBD), vTrialBulkCommit(NEBD),
    vSolved(NEBD), solvedValid(false),
    vSolvedCommit(NEBD), solvedCommitValid(false)
{
  load.Zero();
  vRestore[0] = vRestore[1] = vRestore[2] = 0.0;
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) {
      hingeHinv[m][n] = 0.0;
      hingeHinvCommit[m][n] = 0.0;
    }
    for (int a = 0; a < 3; a++) {
      hingeKva[a][m] = 0.0;
      hingeKvaCommit[a][m] = 0.0;
      hingeKav[m][a] = 0.0;
      hingeKavCommit[m][a] = 0.0;
    }
  }
  theNodes[0] = 0;
  theNodes[1] = 0;
  for (int i = 0; i < maxNumEleLoads; i++) {
    eleLoads[i] = 0;
    eleLoadFactors[i] = 0.0;
  }
}

LadrunoForceBeamColumn2d::~LadrunoForceBeamColumn2d()
{
  if (sections) {
    for (int i = 0; i < numSections; i++)
      if (sections[i])
        delete sections[i];
    delete [] sections;
  }
  if (fs) delete [] fs;
  if (vs) delete [] vs;
  if (Ssr) delete [] Ssr;
  if (vscommit) delete [] vscommit;
  if (fscommit) delete [] fscommit;
  if (Ssrcommit) delete [] Ssrcommit;
  if (crdTransf) delete crdTransf;
  if (beamIntegr) delete beamIntegr;
  if (Ki) delete Ki;
  if (theDamping) delete theDamping;
  if (theHinge) delete theHinge;
  if (theHingeJ) delete theHingeJ;
}

int
LadrunoForceBeamColumn2d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
LadrunoForceBeamColumn2d::getExternalNodes(void)
{
  return connectedExternalNodes;
}

Node **
LadrunoForceBeamColumn2d::getNodePtrs(void)
{
  return theNodes;
}

int
LadrunoForceBeamColumn2d::getNumDOF(void)
{
  return NEGD;
}

void
LadrunoForceBeamColumn2d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }

  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  if (theNodes[0] == 0 || theNodes[1] == 0) {
    opserr << "LadrunoForceBeamColumn2d::setDomain: Node not found in domain\n";
    return;
  }

  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  if (dofNd1 != NND || dofNd2 != NND) {
    opserr << "LadrunoForceBeamColumn2d::setDomain: node DOF mismatch\n";
    return;
  }

  if (crdTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "LadrunoForceBeamColumn2d::setDomain: Error initializing coordinate transformation\n";
    return;
  }

  if (theDamping && theDamping->setDomain(theDomain, NEBD)) {
    opserr << "LadrunoForceBeamColumn2d::setDomain: Error initializing damping\n";
    return;
  }

  double L = crdTransf->getInitialLength();
  if (L <= 0.0) {
    opserr << "LadrunoForceBeamColumn2d::setDomain: Zero element length\n";
    return;
  }

  if (initialFlag == 0)
    this->initializeSectionHistoryVariables();

  this->DomainComponent::setDomain(theDomain);

  ops_TheActiveElement = this;
  this->update();
  ops_TheActiveElement = 0;
}

int
LadrunoForceBeamColumn2d::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping) {
    if (theDamping) delete theDamping;
    theDamping = (*damping).getCopy();
    if (!theDamping) {
      opserr << "LadrunoForceBeamColumn2d::setDamping -- failed to get copy of damping\n";
      return -1;
    }
    if (theDamping->setDomain(theDomain, NEBD)) {
      opserr << "LadrunoForceBeamColumn2d::setDamping -- Error initializing damping\n";
      return -2;
    }
  }
  return 0;
}

int
LadrunoForceBeamColumn2d::commitState(void)
{
  int err = 0;
  int i = 0;

  if ((err = this->Element::commitState()) != 0) {
    opserr << "LadrunoForceBeamColumn2d::commitState() - Element::commitState failed\n";
    return err;
  }

  do {
    vscommit[i] = vs[i];
    fscommit[i] = fs[i];
    Ssrcommit[i] = Ssr[i];
    err = sections[i++]->commitState();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  if ((err = crdTransf->commitState()) != 0)
    return err;

  if (theDamping && (err = theDamping->commitState()) != 0)
    return err;

  if (hingeOn) {
    if (theHinge != 0 && (err = theHinge->commitState()) != 0)
      return err;
    if (theHingeJ != 0 && (err = theHingeJ->commitState()) != 0)
      return err;
    hingeJumpCommit = hingeJump;
    hingeJumpJCommit = hingeJumpJ;
    hingeMscaleCommit = hingeMscale;
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++)
        hingeHinvCommit[m][n] = hingeHinv[m][n];
      for (int a = 0; a < 3; a++) {
        hingeKvaCommit[a][m] = hingeKva[a][m];
        hingeKavCommit[m][a] = hingeKav[m][a];
      }
    }
  }

  const Vector &v = crdTransf->getBasicTrialDisp();
  vRestore[0] = v(0);
  vRestore[1] = v(1);
  vRestore[2] = v(2);
  vRestoreValid = true;
  postRestoreSync = false;

  double betaI = 1.0 - hingeXiI, gammaI = -hingeXiI;
  double betaJ = hingeXiJ - 1.0, gammaJ = hingeXiJ;
  vCommitBulk(0) = v(0);
  vCommitBulk(1) = v(1) - (hingeOn && theHinge ? betaI * hingeJumpCommit : 0.0)
                        - (hingeOn && theHingeJ ? betaJ * hingeJumpJCommit : 0.0);
  vCommitBulk(2) = v(2) - (hingeOn && theHinge ? gammaI * hingeJumpCommit : 0.0)
                        - (hingeOn && theHingeJ ? gammaJ * hingeJumpJCommit : 0.0);

  kvcommit = kv;
  kvbulkCommit = kvbulk;
  Secommit = Se;
  vTrialBulkCommit = vTrialBulk;
  vSolvedCommit = vSolved;
  solvedCommitValid = solvedValid;
  solvedLoadsCommit = solvedLoads;
  solvedLoadDataCommit = solvedLoadData;
  eleImplexErrCommit = eleImplexErr;

  return err;
}

int
LadrunoForceBeamColumn2d::revertToLastCommit(void)
{
  int err = 0;
  int i = 0;

  double L = (crdTransf != 0) ? crdTransf->getInitialLength() : 0.0;
  double wt[maxNumSections];
  if (beamIntegr != 0 && L > 0.0)
    beamIntegr->getSectionWeights(numSections, L, wt);
  else
    for (int k = 0; k < numSections; k++) wt[k] = 1.0;

  ops_TheActiveElement = this;
  do {
    if (lchMode == 1)
      current_section_lch = L;
    else if (lchMode == 2)
      current_section_lch = userLch;
    else
      current_section_lch = wt[i] * L;

    vs[i] = vscommit[i];
    err = sections[i]->revertToLastCommit();
    sections[i]->setTrialSectionDeformation(vs[i]);
    Ssr[i] = sections[i]->getStressResultant();
    fs[i] = sections[i]->getSectionFlexibility();
    Ssrcommit[i] = Ssr[i];
    fscommit[i] = fs[i];
    i++;
  } while (err == 0 && i < numSections);
  current_section_lch = (lchMode == 2) ? userLch : L;
  ops_TheActiveElement = 0;

  if (err)
    return err;

  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;

  if (theDamping && (err = theDamping->revertToLastCommit()) != 0)
    return err;

  if (hingeOn) {
    hingeJump = hingeJumpCommit;
    hingeJumpJ = hingeJumpJCommit;
    hingeMscale = hingeMscaleCommit;
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++)
        hingeHinv[m][n] = hingeHinvCommit[m][n];
      for (int a = 0; a < 3; a++) {
        hingeKva[a][m] = hingeKvaCommit[a][m];
        hingeKav[m][a] = hingeKavCommit[m][a];
      }
    }
    if (theHinge) err += theHinge->revertToLastCommit();
    if (theHingeJ) err += theHingeJ->revertToLastCommit();
  }

  double betaI = 1.0 - hingeXiI, gammaI = -hingeXiI;
  double betaJ = hingeXiJ - 1.0, gammaJ = hingeXiJ;
  vCommitBulk(0) = vRestore[0];
  vCommitBulk(1) = vRestore[1] - (hingeOn && theHinge ? betaI * hingeJumpCommit : 0.0)
                               - (hingeOn && theHingeJ ? betaJ * hingeJumpJCommit : 0.0);
  vCommitBulk(2) = vRestore[2] - (hingeOn && theHinge ? gammaI * hingeJumpCommit : 0.0)
                               - (hingeOn && theHingeJ ? gammaJ * hingeJumpJCommit : 0.0);
  vTrialBulk = vTrialBulkCommit;
  // The trial state is the committed one again, solved at the committed basic
  // deformation and loads: restore that record, so the load re-application of the
  // next analyze does not re-solve at zero increment.
  vSolved = vSolvedCommit;
  solvedValid = solvedCommitValid;
  solvedLoads = solvedLoadsCommit;
  solvedLoadData = solvedLoadDataCommit;

  Se = Secommit;
  kv = kvcommit;
  kvbulk = kvbulkCommit;
  eleImplexErr = eleImplexErrCommit;
  eleLoadDirty = (numEleLoads > 0);
  initialFlag = 1;

  return err;
}

int
LadrunoForceBeamColumn2d::revertToStart(void)
{
  int err = 0;
  int i = 0;

  Secommit.Zero();
  kvcommit.Zero();
  kvbulkCommit.Zero();
  Se.Zero();
  kv.Zero();
  kvbulk.Zero();
  vCommitBulk.Zero();
  vTrialBulk.Zero();
  vTrialBulkCommit.Zero();
  solvedValid = false;
  solvedCommitValid = false;

  if (Ki != 0) {
    delete Ki;
    Ki = 0;
  }

  do {
    vscommit[i].Zero();
    vs[i].Zero();
    Ssr[i].Zero();
    Ssrcommit[i].Zero();
    fs[i].Zero();
    fscommit[i].Zero();
    err = sections[i++]->revertToStart();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  this->initializeSectionHistoryVariables();

  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  if (theDamping && (err = theDamping->revertToStart()) != 0)
    return err;

  if (hingeOn) {
    hingeJump = 0.0;
    hingeJumpCommit = 0.0;
    hingeJumpJ = 0.0;
    hingeJumpJCommit = 0.0;
    hingeMscale = 0.0;
    hingeMscaleCommit = 0.0;
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++) {
        hingeHinv[m][n] = 0.0;
        hingeHinvCommit[m][n] = 0.0;
      }
      for (int a = 0; a < 3; a++) {
        hingeKva[a][m] = 0.0;
        hingeKvaCommit[a][m] = 0.0;
        hingeKav[m][a] = 0.0;
        hingeKavCommit[m][a] = 0.0;
      }
    }
    if (theHinge) err += theHinge->revertToStart();
    if (theHingeJ) err += theHingeJ->revertToStart();
  }

  eleImplexErr = 0.0;
  eleImplexErrCommit = 0.0;
  postRestoreSync = false;
  vRestoreValid = false;
  eleLoadDirty = false;
  initialFlag = 0;

  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
    ops_TheActiveElement = this;
    err += this->update();
    ops_TheActiveElement = 0;
  }

  return err;
}

const Matrix &
LadrunoForceBeamColumn2d::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  static thread_local Matrix f(NEBD, NEBD);
  static thread_local Matrix kvInit(NEBD, NEBD);

  this->getInitialFlexibility(f, false);
  f.Invert(kvInit);

  if (hingeOn) {
    double betaI = 1.0 - hingeXiI, gammaI = -hingeXiI;
    double betaJ = hingeXiJ - 1.0, gammaJ = hingeXiJ;
    bool hasI = (theHinge != 0);
    bool hasJ = (theHingeJ != 0);
    double KcohI = hasI ? theHinge->getInitialTangent() : 0.0;
    double KcohJ = hasJ ? theHingeJ->getInitialTangent() : 0.0;
    double Kva0[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    for (int a = 0; a < 3; a++) {
      Kva0[a][0] = hasI ? -(kvInit(a, 1) * betaI + kvInit(a, 2) * gammaI) : 0.0;
      Kva0[a][1] = hasJ ? -(kvInit(a, 1) * betaJ + kvInit(a, 2) * gammaJ) : 0.0;
    }
    double H00 = hasI ? (-(betaI * Kva0[1][0] + gammaI * Kva0[2][0]) + KcohI) : 1.0;
    double H11 = hasJ ? (-(betaJ * Kva0[1][1] + gammaJ * Kva0[2][1]) + KcohJ) : 1.0;
    double H01 = (hasI && hasJ) ? -0.5 * ((betaI * Kva0[1][1] + gammaI * Kva0[2][1]) +
                                          (betaJ * Kva0[1][0] + gammaJ * Kva0[2][0])) : 0.0;
    double Hinv0[2][2];
    lfbc2dInvertSym2x2Guarded(H00, H11, H01, Hinv0);
    if (!hasI) { Hinv0[0][0] = 0.0; Hinv0[0][1] = 0.0; Hinv0[1][0] = 0.0; }
    if (!hasJ) { Hinv0[1][1] = 0.0; Hinv0[0][1] = 0.0; Hinv0[1][0] = 0.0; }
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        double c = 0.0;
        for (int m = 0; m < 2; m++)
          for (int n = 0; n < 2; n++)
            c += Kva0[a][m] * Hinv0[m][n] * Kva0[b][n];
        kvInit(a, b) -= c;
      }
    }
  }

  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));
  return *Ki;
}

const Matrix &
LadrunoForceBeamColumn2d::getTangentStiff(void)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  crdTransf->update();
  return crdTransf->getGlobalStiffMatrix(kv, Se);
}

void
LadrunoForceBeamColumn2d::computeReactions(double *p0)
{
  int type;
  double L = crdTransf->getInitialLength();

  for (int i = 0; i < numEleLoads; i++) {
    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wa = data(1) * loadFactor;
      double wy = data(0) * loadFactor;

      p0[0] -= wa * L;
      double V = 0.5 * wy * L;
      p0[1] -= V;
      p0[2] -= V;
    }
    else if (type == LOAD_TAG_BeamUniformMoment) {
      double mz = data(2) * loadFactor;
      p0[1] += mz;
      p0[2] -= mz;
    }
    else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
      double wya = data(0) * loadFactor;
      double wyb = data(1) * loadFactor;
      double waa = data(2) * loadFactor;
      double wab = data(3) * loadFactor;
      double a = data(4) * L;
      double b = data(5) * L;
      if (b <= a || a < 0.0 || b > L)
        continue;

      p0[0] -= waa * (b - a) + 0.5 * (wab - waa) * (b - a);
      double Fy = wya * (b - a);
      double c = a + 0.5 * (b - a);
      p0[1] -= Fy * (1.0 - c / L);
      p0[2] -= Fy * c / L;
      Fy = 0.5 * (wyb - wya) * (b - a);
      c = a + 2.0 / 3.0 * (b - a);
      p0[1] -= Fy * (1.0 - c / L);
      p0[2] -= Fy * c / L;
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0) * loadFactor;
      double N = data(1) * loadFactor;
      double aOverL = data(2);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double V1 = P * (1.0 - aOverL);
      double V2 = P * aOverL;

      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
    }
  }
}

void
LadrunoForceBeamColumn2d::computeReactionSensitivity(double *dp0dh, int gradNumber)
{
  int type;
  double L = crdTransf->getInitialLength();
  double dLdh = crdTransf->getdLdh();

  for (int i = 0; i < numEleLoads; i++) {
    const Vector &data = eleLoads[i]->getData(type, 1.0);

    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wy = data(0) * 1.0;
      double wa = data(1) * 1.0;

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwadh = sens(1);

      dp0dh[0] -= wa * dLdh + dwadh * L;
      double dVdh = 0.5 * (wy * dLdh + dwydh * L);
      dp0dh[1] -= dVdh;
      dp0dh[2] -= dVdh;
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0) * 1.0;
      double N = data(1) * 1.0;
      double aOverL = data(2);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPdh = sens(0);
      double dNdh = sens(1);
      double daLdh = sens(2);

      double dV1dh = P * (0.0 - daLdh) + dPdh * (1.0 - aOverL);
      double dV2dh = P * daLdh + dPdh * aOverL;

      dp0dh[0] -= dNdh;
      dp0dh[1] -= dV1dh;
      dp0dh[2] -= dV2dh;
    }
  }
}

const Vector &
LadrunoForceBeamColumn2d::getResistingForce(void)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  crdTransf->update();

  double p0[3];
  Vector p0Vec(p0, 3);
  p0Vec.Zero();

  if (numEleLoads > 0)
    this->computeReactions(p0);

  theVector = crdTransf->getGlobalResistingForce(Se, p0Vec);

  if (rho != 0)
    theVector.addVector(1.0, load, -1.0);

  return theVector;
}

const Vector &
LadrunoForceBeamColumn2d::getDampingForce(void)
{
  crdTransf->update();
  return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(3));
}

void
LadrunoForceBeamColumn2d::initializeSectionHistoryVariables(void)
{
  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    fs[i] = sections[i]->getInitialFlexibility();
    vs[i] = Vector(order);
    Ssr[i] = Vector(order);
    vscommit[i] = Vector(order);
    fscommit[i] = fs[i];
    Ssrcommit[i] = Ssr[i];
  }
}

int
LadrunoForceBeamColumn2d::evalBulkState(const Vector &v_tilde,
                                        Vector &Se_out, Matrix &Kbulk_out, bool commitSubdivide)
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  static thread_local Vector vr(NEBD);
  static thread_local Matrix f(NEBD, NEBD);
  static thread_local Matrix I(NEBD, NEBD);
  static thread_local Vector dv(NEBD);
  static thread_local Vector vin(NEBD);
  static thread_local Vector dSe(NEBD);
  static thread_local Vector dvToDo(NEBD);
  static thread_local Vector dvTrial(NEBD);
  static thread_local Vector SeTrial(NEBD);
  static thread_local Matrix kvTrial(NEBD, NEBD);
  static thread_local Vector SeCur(NEBD);
  static thread_local Matrix kvCur(NEBD, NEBD);
  static thread_local Vector vsCur[maxNumSections];
  static thread_local Matrix fsCur[maxNumSections];
  static thread_local Vector SsrCur[maxNumSections];

  I.Zero();
  for (int i = 0; i < NEBD; i++)
    I(i, i) = 1.0;

  if (initialFlag == 0) {
    bool kvZero = true;
    for (int r = 0; r < NEBD && kvZero; r++)
      for (int c = 0; c < NEBD; c++)
        if (kvbulk(r, c) != 0.0) { kvZero = false; break; }
    if (kvZero) {
      this->getInitialFlexibility(f, false);
      f.Solve(I, kvbulk);
      kvbulkCommit = kvbulk;
    }
  }

  // Iterate from the last TRIAL bulk state, as ForceBeamColumn2d does. Restarting
  // every call from the committed state re-solves the whole step to the element
  // tolerance each time, so the force carries a tolerance-sized error that no
  // global Newton iteration removes: the residual stalls at that floor instead of
  // converging quadratically. revertToLastCommit resets the trial state and
  // vTrialBulk to the committed ones.
  SeCur = Se;
  kvCur = kvbulk;
  for (int i = 0; i < numSections; i++) {
    vsCur[i] = vs[i];
    fsCur[i] = fs[i];
    SsrCur[i] = Ssr[i];
  }

  vin = vTrialBulk;
  dvToDo = v_tilde;
  dvToDo -= vTrialBulk;
  dvTrial = dvToDo;

  double factor = subdivideFactor;
  int numSubdivide = 1;
  bool converged = false;

  while (converged == false && numSubdivide <= maxSubdivisions) {
    for (int l = 0; l < 3; l++) {
      SeTrial = SeCur;
      kvTrial = kvCur;
      for (int i = 0; i < numSections; i++) {
        vsSubdivide[i] = vsCur[i];
        fsSubdivide[i] = fsCur[i];
        SsrSubdivide[i] = SsrCur[i];
      }

      dSe.addMatrixVector(0.0, kvTrial, dvTrial, 1.0);
      SeTrial += dSe;

      if (initialFlag != 2) {
        int numIters = maxIters;
        if (l == 1)
          numIters = 10 * maxIters;

        for (int j = 0; j < numIters; j++) {
          f.Zero();
          vr.Zero();

          for (int i = 0; i < numSections; i++) {
            int order = sections[i]->getOrder();
            const ID &code = sections[i]->getType();

            static thread_local Vector Ss;
            static thread_local Vector dSs;
            static thread_local Vector dvs;
            static thread_local Matrix fb;

            Ss.setData(workArea, order);
            dSs.setData(&workArea[order], order);
            dvs.setData(&workArea[2 * order], order);
            fb.setData(&workArea[3 * order], order, NEBD);

            double xL = xi[i];
            double xL1 = xL - 1.0;
            double wtL = wt[i] * L;

            if (lchMode == 1)
              current_section_lch = L;
            else if (lchMode == 2)
              current_section_lch = userLch;
            else
              current_section_lch = wtL;

            for (int ii = 0; ii < order; ii++) {
              switch (code(ii)) {
              case SECTION_RESPONSE_P:
                Ss(ii) = axialRelease ? 0.0 : SeTrial(0);
                break;
              case SECTION_RESPONSE_MZ:
                Ss(ii) = xL1 * SeTrial(1) + xL * SeTrial(2);
                break;
              case SECTION_RESPONSE_VY:
                Ss(ii) = oneOverL * (SeTrial(1) + SeTrial(2));
                break;
              default:
                Ss(ii) = 0.0;
                break;
              }
            }

            if (numEleLoads > 0) {
              this->computeSectionForces(Ss, i);
              if (axialRelease) {
                for (int ii = 0; ii < order; ii++)
                  if (code(ii) == SECTION_RESPONSE_P)
                    Ss(ii) = 0.0;
              }
            }

            dSs = Ss;
            dSs.addVector(1.0, SsrSubdivide[i], -1.0);

            if (l == 0) {
              dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
            } else if (l == 2) {
              if (j == 0) {
                const Matrix &fs0 = sections[i]->getInitialFlexibility();
                dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
              } else {
                dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
              }
            } else {
              const Matrix &fs0 = sections[i]->getInitialFlexibility();
              dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
            }

            if (initialFlag != 0 || j > 0)
              vsSubdivide[i] += dvs;

            if (sections[i]->setTrialSectionDeformation(vsSubdivide[i]) < 0) {
              opserr << "LadrunoForceBeamColumn2d::evalBulkState() - section failed in setTrial\n";
              return -1;
            }

            SsrSubdivide[i] = sections[i]->getStressResultant();
            fsSubdivide[i] = sections[i]->getSectionFlexibility();

            dSs = Ss;
            dSs.addVector(1.0, SsrSubdivide[i], -1.0);
            dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);

            const Matrix &fSec = fsSubdivide[i];
            fb.Zero();
            double tmp;
            for (int ii = 0; ii < order; ii++) {
              switch (code(ii)) {
              case SECTION_RESPONSE_P:
                if (!axialRelease) {
                  for (int jj = 0; jj < order; jj++)
                    fb(jj, 0) += fSec(jj, ii) * wtL;
                }
                break;
              case SECTION_RESPONSE_MZ:
                for (int jj = 0; jj < order; jj++) {
                  tmp = fSec(jj, ii) * wtL;
                  fb(jj, 1) += xL1 * tmp;
                  fb(jj, 2) += xL * tmp;
                }
                break;
              case SECTION_RESPONSE_VY:
                for (int jj = 0; jj < order; jj++) {
                  tmp = oneOverL * fSec(jj, ii) * wtL;
                  fb(jj, 1) += tmp;
                  fb(jj, 2) += tmp;
                }
                break;
              default:
                break;
              }
            }

            for (int ii = 0; ii < order; ii++) {
              switch (code(ii)) {
              case SECTION_RESPONSE_P:
                if (axialRelease) {
                  const Matrix &fs0 = sections[i]->getInitialFlexibility();
                  f(0, 0) += fs0(ii, ii) * wtL;
                } else {
                  for (int jj = 0; jj < NEBD; jj++)
                    f(0, jj) += fb(ii, jj);
                }
                break;
              case SECTION_RESPONSE_MZ:
                for (int jj = 0; jj < NEBD; jj++) {
                  tmp = fb(ii, jj);
                  f(1, jj) += xL1 * tmp;
                  f(2, jj) += xL * tmp;
                }
                break;
              case SECTION_RESPONSE_VY:
                for (int jj = 0; jj < NEBD; jj++) {
                  tmp = oneOverL * fb(ii, jj);
                  f(1, jj) += tmp;
                  f(2, jj) += tmp;
                }
                break;
              default:
                break;
              }
            }

            dvs.addVector(1.0, vsSubdivide[i], 1.0);
            for (int ii = 0; ii < order; ii++) {
              double dei = dvs(ii) * wtL;
              switch (code(ii)) {
              case SECTION_RESPONSE_P:
                if (axialRelease) {
                  const Matrix &fs0 = sections[i]->getInitialFlexibility();
                  vr(0) += fs0(ii, ii) * SeTrial(0) * wtL;
                } else {
                  vr(0) += dei;
                }
                break;
              case SECTION_RESPONSE_MZ:
                vr(1) += xL1 * dei;
                vr(2) += xL * dei;
                break;
              case SECTION_RESPONSE_VY:
                tmp = oneOverL * dei;
                vr(1) += tmp;
                vr(2) += tmp;
                break;
              default:
                break;
              }
            }
          }

          current_section_lch = (lchMode == 2) ? userLch : L;

          if (GAvy > 0.0) {
            double fsY = oneOverL / GAvy;
            f(1, 1) += fsY;
            f(1, 2) += fsY;
            f(2, 1) += fsY;
            f(2, 2) += fsY;
            double dvShear = fsY * (SeTrial(1) + SeTrial(2));
            vr(1) += dvShear;
            vr(2) += dvShear;
          }

          if (f.Solve(I, kvTrial) < 0)
            return -1;

          dv = vin;
          dv += dvTrial;
          dv -= vr;

          dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);
          double dW = dv ^ dSe;
          SeTrial += dSe;

          if (fabs(dW) < tol) {
            dvToDo -= dvTrial;
            vin += dvTrial;

            SeCur = SeTrial;
            kvCur = kvTrial;
            for (int k = 0; k < numSections; k++) {
              vsCur[k] = vsSubdivide[k];
              fsCur[k] = fsSubdivide[k];
              SsrCur[k] = SsrSubdivide[k];
            }

            if (dvToDo.Norm() <= DBL_EPSILON) {
              converged = true;
            } else {
              dvTrial = dvToDo;
              numSubdivide = 1;
            }

            j = numIters + 1;
            l = 4;
          } else {
            if (j == (numIters - 1) && (l == 2)) {
              dvTrial /= factor;
              numSubdivide++;
            }
          }
        }
      } else {
        SeCur = SeTrial;
        kvCur = kvTrial;
        converged = true;
        break;
      }
    }
  }

  if (!converged)
    return -1;

  current_section_lch = (lchMode == 2) ? userLch : L;
  Se_out = SeCur;
  Kbulk_out = kvCur;
  if (commitSubdivide) {
    for (int k = 0; k < numSections; k++) {
      vs[k] = vsCur[k];
      fs[k] = fsCur[k];
      Ssr[k] = SsrCur[k];
    }
    vTrialBulk = v_tilde;
  }
  return 0;
}

int
LadrunoForceBeamColumn2d::solveHingeState(const Vector &v)
{
  const bool hasI = (theHinge != 0);
  const bool hasJ = (theHingeJ != 0);
  const bool freezeI = postRestoreSync && hasI && lfbc2dHingeIsImplex(theHinge);
  const bool freezeJ = postRestoreSync && hasJ && lfbc2dHingeIsImplex(theHingeJ);
  const double betaI = 1.0 - hingeXiI;
  const double gammaI = -hingeXiI;
  const double betaJ = hingeXiJ - 1.0;
  const double gammaJ = hingeXiJ;

  double alphaI = hasI ? hingeJumpCommit : 0.0;
  double alphaJ = hasJ ? hingeJumpJCommit : 0.0;

  LadrunoCohesiveHinge *cohI = hasI ? dynamic_cast<LadrunoCohesiveHinge*>(theHinge) : 0;
  LadrunoCohesiveHinge *cohJ = hasJ ? dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) : 0;

  static thread_local Vector v_tilde(NEBD);
  static thread_local Vector q_bulk(NEBD);
  static thread_local Matrix Kbulk(NEBD, NEBD);

  double Hinv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double Kva[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  double Kav[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double gI = 0.0, gJ = 0.0;

  const int maxOuter = 20;
  for (int outer = 0; outer < maxOuter; outer++) {
    v_tilde(0) = v(0);
    v_tilde(1) = v(1) - (hasI ? betaI * alphaI : 0.0) - (hasJ ? betaJ * alphaJ : 0.0);
    v_tilde(2) = v(2) - (hasI ? gammaI * alphaI : 0.0) - (hasJ ? gammaJ * alphaJ : 0.0);

    int err = this->evalBulkState(v_tilde, q_bulk, Kbulk, false);
    if (err != 0)
      return err;

    for (int a = 0; a < 3; a++) {
      Kva[a][0] = hasI ? -(Kbulk(a, 1) * betaI + Kbulk(a, 2) * gammaI) : 0.0;
      Kva[a][1] = hasJ ? -(Kbulk(a, 1) * betaJ + Kbulk(a, 2) * gammaJ) : 0.0;
    }

    if (cohI != 0) cohI->setAxialForce(q_bulk(0));
    if (cohJ != 0) cohJ->setAxialForce(q_bulk(0));
    double dfPMIDN = 0.0, dfPMJDN = 0.0;
    double fPMI = (cohI == 0 && pmP0 > 0.0)
                      ? lfbc2dEvalPM(q_bulk(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dfPMIDN)
                      : 1.0;
    double fPMJ = (cohJ == 0 && pmP0 > 0.0)
                      ? lfbc2dEvalPM(q_bulk(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dfPMJDN)
                      : 1.0;

    double McohI = 0.0, KcohI = 0.0, hI = 0.0;
    gI = 0.0;
    if (hasI) {
      if (!freezeI)
        theHinge->setTrialStrain(alphaI);
      double rawM = theHinge->getStress();
      McohI = fPMI * rawM;
      KcohI = fPMI * theHinge->getTangent();
      if (!freezeI) {
        if (cohI != 0 && cohI->hasPM() && !cohI->isImplex()) {
          double fM = cohI->getPMMultiplier();
          if (fM > 0.0) gI = (cohI->getPMDerivative() / fM) * McohI;
        } else if (cohI == 0 && pmP0 > 0.0) {
          gI = dfPMIDN * rawM;
        }
      }
      double qBulkI = betaI * q_bulk(1) + gammaI * q_bulk(2);
      hI = (freezeI || initialFlag == 0) ? 0.0 : (-qBulkI + McohI);
    }

    double McohJ = 0.0, KcohJ = 0.0, hJ = 0.0;
    gJ = 0.0;
    if (hasJ) {
      if (!freezeJ)
        theHingeJ->setTrialStrain(alphaJ);
      double rawM = theHingeJ->getStress();
      McohJ = fPMJ * rawM;
      KcohJ = fPMJ * theHingeJ->getTangent();
      if (!freezeJ) {
        if (cohJ != 0 && cohJ->hasPM() && !cohJ->isImplex()) {
          double fM = cohJ->getPMMultiplier();
          if (fM > 0.0) gJ = (cohJ->getPMDerivative() / fM) * McohJ;
        } else if (cohJ == 0 && pmP0 > 0.0) {
          gJ = dfPMJDN * rawM;
        }
      }
      double qBulkJ = betaJ * q_bulk(1) + gammaJ * q_bulk(2);
      hJ = (freezeJ || initialFlag == 0) ? 0.0 : (-qBulkJ + McohJ);
    }

    double H00 = hasI ? (-(betaI * Kva[1][0] + gammaI * Kva[2][0]) + KcohI + gI * Kva[0][0]) : 1.0;
    double H11 = hasJ ? (-(betaJ * Kva[1][1] + gammaJ * Kva[2][1]) + KcohJ + gJ * Kva[0][1]) : 1.0;
    double H01 = (hasI && hasJ) ? (-(betaI * Kva[1][1] + gammaI * Kva[2][1]) + gI * Kva[0][1]) : 0.0;
    double H10 = (hasI && hasJ) ? (-(betaJ * Kva[1][0] + gammaJ * Kva[2][0]) + gJ * Kva[0][0]) : 0.0;

    lfbc2dInvertGen2x2Guarded(H00, H01, H10, H11, Hinv);
    if (!hasI) { Hinv[0][0] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }
    if (!hasJ) { Hinv[1][1] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }

    double Mhere = fabs(McohI) + fabs(McohJ) + fabs(q_bulk(1)) + fabs(q_bulk(2));
    if (Mhere > hingeMscale) hingeMscale = Mhere;

    double hMag = fabs(hI) + fabs(hJ);
    const double outerTol = 1.0e-12 * hingeMscale + 1.0e-14;
    if (hMag == 0.0 || (outer > 0 && hMag <= outerTol)) {
      break;
    }

    // Solve the inner 2x2 hinge equilibrium against the linearized bulk response
    double aI0 = alphaI;
    double aJ0 = alphaJ;
    double damp = 1.0;
    double hnormPrev = 1.0e300;
    for (int inner = 0; inner < 35; inner++) {
      double daI = alphaI - aI0;
      double daJ = alphaJ - aJ0;
      double q0_lin = q_bulk(0) + Kva[0][0] * daI + Kva[0][1] * daJ;
      double q1_lin = q_bulk(1) + Kva[1][0] * daI + Kva[1][1] * daJ;
      double q2_lin = q_bulk(2) + Kva[2][0] * daI + Kva[2][1] * daJ;

      if (cohI != 0) cohI->setAxialForce(q0_lin);
      if (cohJ != 0) cohJ->setAxialForce(q0_lin);
      dfPMIDN = 0.0;
      dfPMJDN = 0.0;
      fPMI = (cohI == 0 && pmP0 > 0.0)
                 ? lfbc2dEvalPM(q0_lin, pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dfPMIDN)
                 : 1.0;
      fPMJ = (cohJ == 0 && pmP0 > 0.0)
                 ? lfbc2dEvalPM(q0_lin, pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dfPMJDN)
                 : 1.0;

      gI = 0.0;
      if (hasI) {
        if (!freezeI) theHinge->setTrialStrain(alphaI);
        double rawM = theHinge->getStress();
        McohI = fPMI * rawM;
        KcohI = fPMI * theHinge->getTangent();
        if (!freezeI) {
          if (cohI != 0 && cohI->hasPM() && !cohI->isImplex()) {
            double fM = cohI->getPMMultiplier();
            if (fM > 0.0) gI = (cohI->getPMDerivative() / fM) * McohI;
          } else if (cohI == 0 && pmP0 > 0.0) {
            gI = dfPMIDN * rawM;
          }
        }
        hI = freezeI ? 0.0 : (-(betaI * q1_lin + gammaI * q2_lin) + McohI);
      }
      gJ = 0.0;
      if (hasJ) {
        if (!freezeJ) theHingeJ->setTrialStrain(alphaJ);
        double rawM = theHingeJ->getStress();
        McohJ = fPMJ * rawM;
        KcohJ = fPMJ * theHingeJ->getTangent();
        if (!freezeJ) {
          if (cohJ != 0 && cohJ->hasPM() && !cohJ->isImplex()) {
            double fM = cohJ->getPMMultiplier();
            if (fM > 0.0) gJ = (cohJ->getPMDerivative() / fM) * McohJ;
          } else if (cohJ == 0 && pmP0 > 0.0) {
            gJ = dfPMJDN * rawM;
          }
        }
        hJ = freezeJ ? 0.0 : (-(betaJ * q1_lin + gammaJ * q2_lin) + McohJ);
      }

      H00 = hasI ? (-(betaI * Kva[1][0] + gammaI * Kva[2][0]) + KcohI + gI * Kva[0][0]) : 1.0;
      H11 = hasJ ? (-(betaJ * Kva[1][1] + gammaJ * Kva[2][1]) + KcohJ + gJ * Kva[0][1]) : 1.0;
      H01 = (hasI && hasJ) ? (-(betaI * Kva[1][1] + gammaI * Kva[2][1]) + gI * Kva[0][1]) : 0.0;
      H10 = (hasI && hasJ) ? (-(betaJ * Kva[1][0] + gammaJ * Kva[2][0]) + gJ * Kva[0][0]) : 0.0;
      lfbc2dInvertGen2x2Guarded(H00, H01, H10, H11, Hinv);
      if (!hasI) { Hinv[0][0] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }
      if (!hasJ) { Hinv[1][1] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }

      double hInMag = fabs(hI) + fabs(hJ);
      const double innerTol = 1.0e-13 * hingeMscale + 1.0e-15;
      if (hInMag == 0.0 || (inner > 0 && hInMag <= innerTol))
        break;

      double hnorm = hInMag / (hingeMscale + 1.0e-12);
      if (hnorm > hnormPrev) damp *= 0.5;
      else damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
      hnormPrev = hnorm;

      if (hasI && !freezeI) alphaI -= damp * (Hinv[0][0] * hI + Hinv[0][1] * hJ);
      if (hasJ && !freezeJ) alphaJ -= damp * (Hinv[1][0] * hI + Hinv[1][1] * hJ);
    }
  }

  v_tilde(0) = v(0);
  v_tilde(1) = v(1) - (hasI ? betaI * alphaI : 0.0) - (hasJ ? betaJ * alphaJ : 0.0);
  v_tilde(2) = v(2) - (hasI ? gammaI * alphaI : 0.0) - (hasJ ? gammaJ * alphaJ : 0.0);
  int err = this->evalBulkState(v_tilde, q_bulk, Kbulk, true);
  if (err != 0)
    return err;

  if (cohI != 0) cohI->setAxialForce(q_bulk(0));
  if (cohJ != 0) cohJ->setAxialForce(q_bulk(0));
  if (hasI && !freezeI) theHinge->setTrialStrain(alphaI);
  if (hasJ && !freezeJ) theHingeJ->setTrialStrain(alphaJ);

  double dfPMIDN = 0.0, dfPMJDN = 0.0;
  double fPMI = (cohI == 0 && pmP0 > 0.0)
                    ? lfbc2dEvalPM(q_bulk(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dfPMIDN)
                    : 1.0;
  double fPMJ = (cohJ == 0 && pmP0 > 0.0)
                    ? lfbc2dEvalPM(q_bulk(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dfPMJDN)
                    : 1.0;
  double KcohI = hasI ? (fPMI * theHinge->getTangent()) : 0.0;
  double KcohJ = hasJ ? (fPMJ * theHingeJ->getTangent()) : 0.0;
  gI = 0.0;
  if (hasI && !freezeI) {
    if (cohI != 0) {
      double fM = cohI->getPMMultiplier();
      if (fM > 0.0) gI = (cohI->getPMDerivative() / fM) * (fPMI * theHinge->getStress());
    }
    else if (pmP0 > 0.0) gI = dfPMIDN * theHinge->getStress();
  }
  gJ = 0.0;
  if (hasJ && !freezeJ) {
    if (cohJ != 0) {
      double fM = cohJ->getPMMultiplier();
      if (fM > 0.0) gJ = (cohJ->getPMDerivative() / fM) * (fPMJ * theHingeJ->getStress());
    }
    else if (pmP0 > 0.0) gJ = dfPMJDN * theHingeJ->getStress();
  }
  for (int a = 0; a < 3; a++) {
    Kva[a][0] = hasI ? -(Kbulk(a, 1) * betaI + Kbulk(a, 2) * gammaI) : 0.0;
    Kva[a][1] = hasJ ? -(Kbulk(a, 1) * betaJ + Kbulk(a, 2) * gammaJ) : 0.0;
    Kav[0][a] = hasI ? (-(betaI * Kbulk(1, a) + gammaI * Kbulk(2, a)) + gI * Kbulk(0, a)) : 0.0;
    Kav[1][a] = hasJ ? (-(betaJ * Kbulk(1, a) + gammaJ * Kbulk(2, a)) + gJ * Kbulk(0, a)) : 0.0;
  }
  double H00 = hasI ? (-(betaI * Kva[1][0] + gammaI * Kva[2][0]) + KcohI + gI * Kva[0][0]) : 1.0;
  double H11 = hasJ ? (-(betaJ * Kva[1][1] + gammaJ * Kva[2][1]) + KcohJ + gJ * Kva[0][1]) : 1.0;
  double H01 = (hasI && hasJ) ? (-(betaI * Kva[1][1] + gammaI * Kva[2][1]) + gI * Kva[0][1]) : 0.0;
  double H10 = (hasI && hasJ) ? (-(betaJ * Kva[1][0] + gammaJ * Kva[2][0]) + gJ * Kva[0][0]) : 0.0;
  lfbc2dInvertGen2x2Guarded(H00, H01, H10, H11, Hinv);
  if (!hasI) { Hinv[0][0] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }
  if (!hasJ) { Hinv[1][1] = 0.0; Hinv[0][1] = 0.0; Hinv[1][0] = 0.0; }

  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) hingeHinv[m][n] = Hinv[m][n];
    for (int a = 0; a < 3; a++) {
      hingeKva[a][m] = Kva[a][m];
      hingeKav[m][a] = Kav[m][a];
    }
  }
  hingeJump = alphaI;
  hingeJumpJ = alphaJ;
  if (postRestoreSync) {
    for (int m = 0; m < 2; m++) {
      for (int n = 0; n < 2; n++) hingeHinvCommit[m][n] = hingeHinv[m][n];
      for (int a = 0; a < 3; a++) {
        hingeKvaCommit[a][m] = hingeKva[a][m];
        hingeKavCommit[m][a] = hingeKav[m][a];
      }
    }
  }

  Se = q_bulk;
  kvbulk = Kbulk;
  kv = Kbulk;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      double c = 0.0;
      for (int m = 0; m < 2; m++)
        for (int n = 0; n < 2; n++)
          c += hingeKva[a][m] * hingeHinv[m][n] * hingeKav[n][b];
      kv(a, b) -= c;
    }
  }

  return 0;
}

bool
LadrunoForceBeamColumn2d::solvedStateMatches(const Vector &v)
{
  if (!solvedValid || (size_t)numEleLoads != solvedLoads.size())
    return false;
  for (int i = 0; i < NEBD; i++)
    if (v(i) != vSolved(i))
      return false;
  size_t k = 0;
  for (int i = 0; i < numEleLoads; i++) {
    if (eleLoads[i] != solvedLoads[i])
      return false;
    int type;
    const Vector &data = eleLoads[i]->getData(type, eleLoadFactors[i]);
    if (k + 2 + data.Size() > solvedLoadData.size())
      return false;
    if (solvedLoadData[k++] != eleLoadFactors[i] || solvedLoadData[k++] != (double)type)
      return false;
    for (int j = 0; j < data.Size(); j++)
      if (solvedLoadData[k++] != data(j))
        return false;
  }
  return k == solvedLoadData.size();
}

void
LadrunoForceBeamColumn2d::recordSolvedState(const Vector &v)
{
  vSolved = v;
  solvedLoads.clear();
  solvedLoadData.clear();
  for (int i = 0; i < numEleLoads; i++) {
    solvedLoads.push_back(eleLoads[i]);
    int type;
    const Vector &data = eleLoads[i]->getData(type, eleLoadFactors[i]);
    solvedLoadData.push_back(eleLoadFactors[i]);
    solvedLoadData.push_back((double)type);
    for (int j = 0; j < data.Size(); j++)
      solvedLoadData.push_back(data(j));
  }
  solvedValid = true;
}

int
LadrunoForceBeamColumn2d::update(void)
{
  ops_TheActiveElement = this;

  if (initialFlag == 2)
    this->revertToLastCommit();

  crdTransf->update();
  const Vector &v = crdTransf->getBasicTrialDisp();

  static thread_local Vector dv(NEBD);
  dv = crdTransf->getBasicIncrDeltaDisp();

  if (initialFlag == 1 && dv.Norm() <= DBL_EPSILON && numEleLoads == 0 && !eleLoadDirty) {
    ops_TheActiveElement = 0;
    return 0;
  }

  // Nothing changed since the last solve: a load re-application (eleLoadDirty,
  // e.g. DisplacementControl::domainChanged's applyLoad + formUnbalance) must not
  // re-solve, or the sub-tolerance section residual left by the last solve is pushed
  // into the fibres and can flip them onto their unloading tangent.
  if (initialFlag == 1 && this->solvedStateMatches(v)) {
    eleLoadDirty = false;
    ops_TheActiveElement = 0;
    return 0;
  }

  if (postRestoreSync) {
    if (!vRestoreValid) {
      vRestore[0] = v(0); vRestore[1] = v(1); vRestore[2] = v(2);
      vRestoreValid = true;
    } else if (v(0) != vRestore[0] || v(1) != vRestore[1] || v(2) != vRestore[2]) {
      postRestoreSync = false;
    }
  }

  static thread_local Vector dvPred(NEBD);
  static thread_local Vector SePred(NEBD);
  dvPred(0) = v(0) - vRestore[0];
  dvPred(1) = v(1) - vRestore[1];
  dvPred(2) = v(2) - vRestore[2];
  SePred = Secommit;
  SePred.addMatrixVector(1.0, kvcommit, dvPred, 1.0);

  int err = 0;
  if (hingeOn) {
    err = this->solveHingeState(v);
  } else {
    err = this->evalBulkState(v, Se, kvbulk, true);
    kv = kvbulk;
  }

  if (err != 0) {
    ops_TheActiveElement = 0;
    return err;
  }

  if (initialFlag == 0 || (initialFlag == 3 && !vRestoreValid)) {
    vTrialBulkCommit = vTrialBulk;
    kvcommit = kv;
    kvbulkCommit = kvbulk;
    Secommit = Se;
    for (int k = 0; k < numSections; k++) {
      vscommit[k] = vs[k];
      fscommit[k] = fs[k];
      Ssrcommit[k] = Ssr[k];
    }
    double betaI = 1.0 - hingeXiI, gammaI = -hingeXiI;
    double betaJ = hingeXiJ - 1.0, gammaJ = hingeXiJ;
    vCommitBulk(0) = v(0);
    vCommitBulk(1) = v(1) - (hingeOn && theHinge ? betaI * hingeJump : 0.0)
                          - (hingeOn && theHingeJ ? betaJ * hingeJumpJ : 0.0);
    vCommitBulk(2) = v(2) - (hingeOn && theHinge ? gammaI * hingeJump : 0.0)
                          - (hingeOn && theHingeJ ? gammaJ * hingeJumpJ : 0.0);
  }

  eleImplexErr = 0.0;
  if (useImplex && initialFlag != 0 && initialFlag != 3) {
    static thread_local Vector dSeErr(NEBD);
    dSeErr = Se;
    dSeErr -= SePred;
    eleImplexErr = dSeErr.Norm() / (Se.Norm() + 1.0e-12);
  }
  if (hingeOn) {
    double eH = 0.0;
    if (lfbc2dGetHingeImplexError(theHinge, eH) && eH > eleImplexErr)
      eleImplexErr = eH;
    if (lfbc2dGetHingeImplexError(theHingeJ, eH) && eH > eleImplexErr)
      eleImplexErr = eH;
  }

  initialFlag = 1;
  eleLoadDirty = false;
  this->recordSolvedState(v);

  if (theDamping && theDamping->update(Se)) {
    opserr << "LadrunoForceBeamColumn2d::update: failed to update damping\n";
    ops_TheActiveElement = 0;
    return -1;
  }

  ops_TheActiveElement = 0;
  return err;
}

void
LadrunoForceBeamColumn2d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
{
  b.Zero();
  double L = crdTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:
      b(i, 1) = xi - 1.0;
      b(i, 2) = xi;
      break;
    case SECTION_RESPONSE_P:
      b(i, 0) = 1.0;
      break;
    case SECTION_RESPONSE_VY:
      b(i, 1) = b(i, 2) = 1.0 / L;
      break;
    default:
      break;
    }
  }
}

void
LadrunoForceBeamColumn2d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
{
  bp.Zero();
  double L = crdTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:
      bp(i, 1) = xi * (xi - 1.0) * L * L / 2.0;
      break;
    case SECTION_RESPONSE_P:
      bp(i, 0) = (1.0 - xi) * L;
      break;
    case SECTION_RESPONSE_VY:
      bp(i, 1) = (xi - 0.5) * L;
      break;
    default:
      break;
    }
  }
}

const Matrix &
LadrunoForceBeamColumn2d::getMass(void)
{
  theMatrix.Zero();
  double L = crdTransf->getInitialLength();
  if (rho != 0.0) {
    if (cMass == 0) {
      theMatrix(0, 0) = theMatrix(1, 1) = theMatrix(3, 3) = theMatrix(4, 4) = 0.5 * rho * L;
    } else {
      static thread_local Matrix ml(6, 6);
      ml.Zero();
      double m = rho * L / 420.0;
      ml(0, 0) = ml(3, 3) = m * 140.0;
      ml(0, 3) = ml(3, 0) = m * 70.0;

      ml(1, 1) = ml(4, 4) = m * 156.0;
      ml(1, 4) = ml(4, 1) = m * 54.0;
      ml(2, 2) = ml(5, 5) = m * 4.0 * L * L;
      ml(2, 5) = ml(5, 2) = -m * 3.0 * L * L;
      ml(1, 2) = ml(2, 1) = m * 22.0 * L;
      ml(4, 5) = ml(5, 4) = -ml(1, 2);
      ml(1, 5) = ml(5, 1) = -m * 13.0 * L;
      ml(2, 4) = ml(4, 2) = -ml(1, 5);

      theMatrix = crdTransf->getGlobalMatrixFromLocal(ml);
    }
  }
  return theMatrix;
}

void
LadrunoForceBeamColumn2d::zeroLoad(void)
{
  if (numEleLoads > 0)
    eleLoadDirty = true;
  load.Zero();
  numEleLoads = 0;
}

int
LadrunoForceBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  if (numEleLoads == maxNumEleLoads) {
    opserr << "LadrunoForceBeamColumn2d::addLoad() - maxNumEleLoads exceeded\n";
    return -1;
  }
  eleLoads[numEleLoads] = theLoad;
  eleLoadFactors[numEleLoads] = loadFactor;
  numEleLoads++;
  eleLoadDirty = true;
  postRestoreSync = false;
  return 0;
}

void
LadrunoForceBeamColumn2d::computeSectionForces(Vector &sp, int isec)
{
  if (!spanLoadMode)
    return;
  int type;
  double L = crdTransf->getInitialLength();
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double x = xi[isec] * L;

  int order = sections[isec]->getOrder();
  const ID &code = sections[isec]->getType();

  for (int i = 0; i < numEleLoads; i++) {
    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wa = data(1) * loadFactor;
      double wy = data(0) * loadFactor;
      for (int ii = 0; ii < order; ii++) {
        switch (code(ii)) {
        case SECTION_RESPONSE_P:
          sp(ii) += wa * (L - x);
          break;
        case SECTION_RESPONSE_MZ:
          sp(ii) += wy * 0.5 * x * (x - L);
          break;
        case SECTION_RESPONSE_VY:
          sp(ii) += wy * (x - 0.5 * L);
          break;
        default:
          break;
        }
      }
    }
    else if (type == LOAD_TAG_BeamUniformMoment) {
      double mz = data(2) * loadFactor;
      for (int ii = 0; ii < order; ii++) {
        if (code(ii) == SECTION_RESPONSE_VY)
          sp(ii) += mz;
      }
    }
    else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
      double wya = data(0) * loadFactor;
      double wyb = data(1) * loadFactor;
      double waa = data(2) * loadFactor;
      double wab = data(3) * loadFactor;
      double a = data(4) * L;
      double b = data(5) * L;
      if (b <= a || a < 0.0 || b > L)
        continue;

      double Fa = waa * (b - a) + 0.5 * (wab - waa) * (b - a);
      double Fy = wya * (b - a);
      double c = a + 0.5 * (b - a);
      double VI = Fy * (1.0 - c / L);
      double VJ = Fy * c / L;
      Fy = 0.5 * (wyb - wya) * (b - a);
      c = a + 2.0 / 3.0 * (b - a);
      VI += Fy * (1.0 - c / L);
      VJ += Fy * c / L;

      for (int ii = 0; ii < order; ii++) {
        if (x <= a) {
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            sp(ii) += Fa;
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) -= VI * x;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) -= VI;
            break;
          default:
            break;
          }
        } else if (x >= b) {
          switch (code(ii)) {
          case SECTION_RESPONSE_MZ:
            sp(ii) += VJ * (x - L);
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) += VJ;
            break;
          default:
            break;
          }
        } else {
          double wx = wya + (wyb - wya) / (b - a) * (x - a);
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            sp(ii) += Fa - waa * (x - a) - 0.5 * (wab - waa) / (b - a) * (x - a) * (x - a);
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) += -VI * x + wya * (x - a) * 0.5 * (x - a) + 0.5 * (wx - wya) * (x - a) * (x - a) / 3.0;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) += -VI + wya * (x - a) + 0.5 * (wx - wya) * (x - a);
            break;
          default:
            break;
          }
        }
      }
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0) * loadFactor;
      double N = data(1) * loadFactor;
      double aOverL = data(2);
      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double a = aOverL * L;
      double V1 = P * (1.0 - aOverL);
      double V2 = P * aOverL;

      for (int ii = 0; ii < order; ii++) {
        if (x <= a) {
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            sp(ii) += N;
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) -= x * V1;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) -= V1;
            break;
          default:
            break;
          }
        } else {
          switch (code(ii)) {
          case SECTION_RESPONSE_MZ:
            sp(ii) -= (L - x) * V2;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) += V2;
            break;
          default:
            break;
          }
        }
      }
    }
    else {
      opserr << "LadrunoForceBeamColumn2d::addLoad -- load type unknown for element with tag: "
             << this->getTag() << endln;
    }
  }
}

void
LadrunoForceBeamColumn2d::computeSectionForceSensitivity(Vector &dspdh, int isec, int gradNumber)
{
  int type;
  double L = crdTransf->getInitialLength();
  double dLdh = crdTransf->getdLdh();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double dxidh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dxidh);

  double x = xi[isec] * L;
  double dxdh = xi[isec] * dLdh + dxidh[isec] * L;

  int order = sections[isec]->getOrder();
  const ID &code = sections[isec]->getType();

  for (int i = 0; i < numEleLoads; i++) {
    const Vector &data = eleLoads[i]->getData(type, 1.0);
    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wy = data(0) * 1.0;
      double wa = data(1) * 1.0;
      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwadh = sens(1);
      for (int ii = 0; ii < order; ii++) {
        switch (code(ii)) {
        case SECTION_RESPONSE_P:
          dspdh(ii) += dwadh * (L - x) + wa * (dLdh - dxdh);
          break;
        case SECTION_RESPONSE_MZ:
          dspdh(ii) += 0.5 * (dwydh * x * (x - L) + wy * (dxdh * (2.0 * x - L) - x * dLdh));
          break;
        case SECTION_RESPONSE_VY:
          dspdh(ii) += dwydh * (x - 0.5 * L) + wy * (dxdh - 0.5 * dLdh);
          break;
        default:
          break;
        }
      }
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0) * 1.0;
      double N = data(1) * 1.0;
      double aOverL = data(2);
      if (aOverL < 0.0 || aOverL > 1.0)
        continue;
      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPdh = sens(0);
      double dNdh = sens(1);
      double daLdh = sens(2);
      double a = aOverL * L;
      double V1 = P * (1.0 - aOverL);
      double V2 = P * aOverL;
      double dV1dh = P * (0.0 - daLdh) + dPdh * (1.0 - aOverL);
      double dV2dh = P * daLdh + dPdh * aOverL;
      for (int ii = 0; ii < order; ii++) {
        if (x <= a) {
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            dspdh(ii) += dNdh;
            break;
          case SECTION_RESPONSE_MZ:
            dspdh(ii) -= (dxdh * V1 + x * dV1dh);
            break;
          case SECTION_RESPONSE_VY:
            dspdh(ii) -= dV1dh;
            break;
          default:
            break;
          }
        } else {
          switch (code(ii)) {
          case SECTION_RESPONSE_MZ:
            dspdh(ii) -= (dLdh - dxdh) * V2 + (L - x) * dV2dh;
            break;
          case SECTION_RESPONSE_VY:
            dspdh(ii) += dV2dh;
            break;
          default:
            break;
          }
        }
      }
    }
  }
}

int
LadrunoForceBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (cMass == 0) {
    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L;
    load(0) -= m * Raccel1(0);
    load(1) -= m * Raccel1(1);
    load(3) -= m * Raccel2(0);
    load(4) -= m * Raccel2(1);
  } else {
    static thread_local Vector Raccel(6);
    for (int i = 0; i < 3; i++) {
      Raccel(i) = Raccel1(i);
      Raccel(i + 3) = Raccel2(i);
    }
    load.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  return 0;
}

const Vector &
LadrunoForceBeamColumn2d::getResistingForceIncInertia(void)
{
  static thread_local Vector res(6);

  res = this->getResistingForce();

  if (theDamping)
    res += this->getDampingForce();

  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    if (cMass == 0) {
      double L = crdTransf->getInitialLength();
      double m = 0.5 * rho * L;
      res(0) += m * accel1(0);
      res(1) += m * accel1(1);
      res(3) += m * accel2(0);
      res(4) += m * accel2(1);
    } else {
      static thread_local Vector accel(6);
      for (int i = 0; i < 3; i++) {
        accel(i) = accel1(i);
        accel(i + 3) = accel2(i);
      }
      res.addMatrixVector(1.0, this->getMass(), accel, 1.0);
    }

    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      res += this->getRayleighDampingForces();
  } else {
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      res += this->getRayleighDampingForces();
  }

  theVector = res;
  return theVector;
}

int
LadrunoForceBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  int i, j, k;
  int loc = 0;

  static thread_local ID idData(20);
  idData.Zero();
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;
  idData(6) = crdTransf->getClassTag();
  int crdTransfDbTag = crdTransf->getDbTag();
  if (crdTransfDbTag == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag != 0)
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(7) = crdTransfDbTag;

  idData(8) = beamIntegr->getClassTag();
  int beamIntegrDbTag = beamIntegr->getDbTag();
  if (beamIntegrDbTag == 0) {
    beamIntegrDbTag = theChannel.getDbTag();
    if (beamIntegrDbTag != 0)
      beamIntegr->setDbTag(beamIntegrDbTag);
  }
  idData(9) = beamIntegrDbTag;
  idData(10) = maxSubdivisions;

  idData(11) = 0;
  idData(12) = 0;
  if (theDamping) {
    idData(11) = theDamping->getClassTag();
    int dmpDbTag = theDamping->getDbTag();
    if (dmpDbTag == 0) {
      dmpDbTag = theChannel.getDbTag();
      if (dmpDbTag != 0)
        theDamping->setDbTag(dmpDbTag);
    }
    idData(12) = dmpDbTag;
  }

  idData(13) = cMass + 2 * lchMode + (useImplex ? 8 : 0) + (axialRelease ? 16 : 0) + (spanLoadMode == 0 ? 32 : 0);
  idData(14) = hingeOn + 2 * dualHingeMode;
  idData(15) = (theHinge != 0) ? theHinge->getClassTag() : 0;
  idData(16) = 0;
  if (theHinge != 0) {
    int hDb = theHinge->getDbTag();
    if (hDb == 0) {
      hDb = theChannel.getDbTag();
      if (hDb != 0) theHinge->setDbTag(hDb);
    }
    idData(16) = hDb;
  }
  idData(17) = (theHingeJ != 0) ? theHingeJ->getClassTag() : 0;
  idData(18) = 0;
  if (theHingeJ != 0) {
    int hDbJ = theHingeJ->getDbTag();
    if (hDbJ == 0) {
      hDbJ = theChannel.getDbTag();
      if (hDbJ != 0) theHingeJ->setDbTag(hDbJ);
    }
    idData(18) = hDbJ;
  }
  idData(19) = vRestoreValid ? 1 : 0;

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send crdTrans\n";
    return -1;
  }

  if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send beamIntegr\n";
    return -1;
  }

  ID idSections(2 * numSections);
  loc = 0;
  for (i = 0; i < numSections; i++) {
    int sectClassTag = sections[i]->getClassTag();
    int sectDbTag = sections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      sections[i]->setDbTag(sectDbTag);
    }
    idSections(loc) = sectClassTag;
    idSections(loc + 1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send ID sections\n";
    return -1;
  }

  for (j = 0; j < numSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoForceBeamColumn2d::sendSelf() - section " << j << " failed to send itself\n";
      return -1;
    }
  }

  if (theHinge != 0 && theHinge->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send theHinge\n";
    return -1;
  }
  if (theHingeJ != 0 && theHingeJ->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send theHingeJ\n";
    return -1;
  }

  int secDefSize = 0;
  for (i = 0; i < numSections; i++)
    secDefSize += sections[i]->getOrder();

  Vector dData(3 + NEBD + 2 * NEBD * NEBD + secDefSize + 5 + 8 + 2 * NEBD + 1 + 19);
  loc = 0;
  dData(loc++) = rho;
  dData(loc++) = tol;
  dData(loc++) = subdivideFactor;

  for (i = 0; i < NEBD; i++)
    dData(loc++) = Secommit(i);

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      dData(loc++) = kvcommit(i, j);

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      dData(loc++) = kvbulkCommit(i, j);

  for (k = 0; k < numSections; k++)
    for (i = 0; i < sections[k]->getOrder(); i++)
      dData(loc++) = (vscommit[k])(i);

  dData(loc++) = alphaM;
  dData(loc++) = betaK;
  dData(loc++) = betaK0;
  dData(loc++) = betaKc;
  dData(loc++) = userLch;

  dData(loc++) = hingeJumpCommit;
  dData(loc++) = hingeJumpJCommit;
  dData(loc++) = hingeXiI;
  dData(loc++) = hingeXiJ;
  dData(loc++) = pmP0;
  dData(loc++) = pmAlpha;
  dData(loc++) = hingeMscaleCommit;
  dData(loc++) = eleImplexErrCommit;

  for (i = 0; i < NEBD; i++)
    dData(loc++) = vCommitBulk(i);
  for (i = 0; i < NEBD; i++)
    dData(loc++) = vRestore[i];
  dData(loc++) = GAvy;
  for (int m = 0; m < 2; m++)
    for (int n = 0; n < 2; n++)
      dData(loc++) = hingeHinvCommit[m][n];
  for (int a = 0; a < 3; a++)
    for (int m = 0; m < 2; m++)
      dData(loc++) = hingeKvaCommit[a][m];
  for (int m = 0; m < 2; m++)
    for (int a = 0; a < 3; a++)
      dData(loc++) = hingeKavCommit[m][a];
  dData(loc++) = pmP0t;
  dData(loc++) = pmPbal;
  dData(loc++) = pmMuBal;

  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf() - failed to send Vector data\n";
    return -1;
  }

  if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn2d::sendSelf -- could not send Damping\n";
    return -1;
  }

  return 0;
}

int
LadrunoForceBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  int i, j, k;

  static thread_local ID idData(20);
  if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
    opserr << "LadrunoForceBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters = idData(4);
  initialFlag = idData(5);
  maxSubdivisions = idData(10);
  cMass = idData(13) & 1;
  lchMode = (idData(13) >> 1) & 3;
  useImplex = ((idData(13) & 8) != 0);
  axialRelease = ((idData(13) & 16) != 0);
  spanLoadMode = ((idData(13) & 32) != 0) ? 0 : 1;
  hingeOn = idData(14) & 1;
  dualHingeMode = (idData(14) >> 1) & 1;

  int crdTransfClassTag = idData(6);
  int crdTransfDbTag = idData(7);
  int beamIntegrClassTag = idData(8);
  int beamIntegrDbTag = idData(9);

  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
    if (crdTransf != 0) delete crdTransf;
    crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);
    if (crdTransf == 0) {
      opserr << "LadrunoForceBeamColumn2d::recvSelf() - failed to obtain CrdTransf\n";
      return -1;
    }
  }
  crdTransf->setDbTag(crdTransfDbTag);
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)
    return -3;

  if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntegrClassTag) {
    if (beamIntegr != 0) delete beamIntegr;
    beamIntegr = theBroker.getNewBeamIntegration(beamIntegrClassTag);
    if (beamIntegr == 0) {
      opserr << "LadrunoForceBeamColumn2d::recvSelf() - failed to obtain BeamIntegration\n";
      return -1;
    }
  }
  beamIntegr->setDbTag(beamIntegrDbTag);
  if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0)
    return -3;

  ID idSections(2 * idData(3));
  int loc = 0;
  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)
    return -1;

  if (numSections != idData(3)) {
    if (numSections != 0) {
      for (i = 0; i < numSections; i++)
        delete sections[i];
      delete [] sections;
    }
    numSections = idData(3);
    if (vscommit != 0) delete [] vscommit;
    vscommit = new Vector[numSections];
    if (fscommit != 0) delete [] fscommit;
    fscommit = new Matrix[numSections];
    if (Ssrcommit != 0) delete [] Ssrcommit;
    Ssrcommit = new Vector[numSections];
    if (fs != 0) delete [] fs;
    fs = new Matrix[numSections];
    if (vs != 0) delete [] vs;
    vs = new Vector[numSections];
    if (Ssr != 0) delete [] Ssr;
    Ssr = new Vector[numSections];
    sections = new SectionForceDeformation *[numSections];
    loc = 0;
    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc + 1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0)
        return -1;
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0)
        return -1;
    }
    this->initializeSectionHistoryVariables();
  } else {
    loc = 0;
    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc + 1);
      loc += 2;
      if (sections[i]->getClassTag() != sectClassTag) {
        delete sections[i];
        sections[i] = theBroker.getNewSection(sectClassTag);
        if (sections[i] == 0)
          return -1;
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0)
        return -1;
    }
  }

  int hTagI = idData(15), hDbI = idData(16);
  if (hTagI != 0) {
    if (theHinge == 0 || theHinge->getClassTag() != hTagI) {
      if (theHinge) delete theHinge;
      theHinge = theBroker.getNewUniaxialMaterial(hTagI);
      if (theHinge == 0) return -1;
    }
    theHinge->setDbTag(hDbI);
    if (theHinge->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
  } else if (theHinge != 0) {
    delete theHinge;
    theHinge = 0;
  }

  int hTagJ = idData(17), hDbJ = idData(18);
  if (hTagJ != 0) {
    if (theHingeJ == 0 || theHingeJ->getClassTag() != hTagJ) {
      if (theHingeJ) delete theHingeJ;
      theHingeJ = theBroker.getNewUniaxialMaterial(hTagJ);
      if (theHingeJ == 0) return -1;
    }
    theHingeJ->setDbTag(hDbJ);
    if (theHingeJ->recvSelf(commitTag, theChannel, theBroker) < 0) return -1;
  } else if (theHingeJ != 0) {
    delete theHingeJ;
    theHingeJ = 0;
  }

  int secDefSize = 0;
  for (int ii = 0; ii < numSections; ii++)
    secDefSize += sections[ii]->getOrder();

  Vector dData(3 + NEBD + 2 * NEBD * NEBD + secDefSize + 5 + 8 + 2 * NEBD + 1 + 19);
  if (theChannel.recvVector(dbTag, commitTag, dData) < 0)
    return -1;

  loc = 0;
  rho = dData(loc++);
  tol = dData(loc++);
  subdivideFactor = dData(loc++);

  for (i = 0; i < NEBD; i++)
    Secommit(i) = dData(loc++);

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      kvcommit(i, j) = dData(loc++);

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      kvbulkCommit(i, j) = dData(loc++);

  kv = kvcommit;
  kvbulk = kvbulkCommit;
  Se = Secommit;

  for (k = 0; k < numSections; k++) {
    int order = sections[k]->getOrder();
    vscommit[k] = Vector(order);
    for (i = 0; i < order; i++)
      (vscommit[k])(i) = dData(loc++);
    vs[k] = vscommit[k];
    fs[k] = sections[k]->getSectionFlexibility();
    Ssr[k] = sections[k]->getStressResultant();
    fscommit[k] = fs[k];
    Ssrcommit[k] = Ssr[k];
  }

  alphaM = dData(loc++);
  betaK = dData(loc++);
  betaK0 = dData(loc++);
  betaKc = dData(loc++);
  userLch = dData(loc++);

  hingeJumpCommit = dData(loc++);
  hingeJumpJCommit = dData(loc++);
  hingeXiI = dData(loc++);
  hingeXiJ = dData(loc++);
  pmP0 = dData(loc++);
  pmAlpha = dData(loc++);
  hingeMscaleCommit = dData(loc++);
  eleImplexErrCommit = dData(loc++);

  hingeJump = hingeJumpCommit;
  hingeJumpJ = hingeJumpJCommit;
  hingeMscale = hingeMscaleCommit;
  eleImplexErr = eleImplexErrCommit;

  for (i = 0; i < NEBD; i++)
    vCommitBulk(i) = dData(loc++);
  vTrialBulk = vCommitBulk;
  vTrialBulkCommit = vCommitBulk;
  solvedValid = false;
  solvedCommitValid = false;
  for (i = 0; i < NEBD; i++)
    vRestore[i] = dData(loc++);
  GAvy = dData(loc++);
  for (int m = 0; m < 2; m++)
    for (int n = 0; n < 2; n++)
      hingeHinv[m][n] = hingeHinvCommit[m][n] = dData(loc++);
  for (int a = 0; a < 3; a++)
    for (int m = 0; m < 2; m++)
      hingeKva[a][m] = hingeKvaCommit[a][m] = dData(loc++);
  for (int m = 0; m < 2; m++)
    for (int a = 0; a < 3; a++)
      hingeKav[m][a] = hingeKavCommit[m][a] = dData(loc++);
  pmP0t = dData(loc++);
  pmPbal = dData(loc++);
  pmMuBal = dData(loc++);

  if (Ki != 0) {
    delete Ki;
    Ki = 0;
  }

  initialFlag = 2;
  vRestoreValid = (idData(19) != 0);
  postRestoreSync = (hingeOn != 0);
  eleLoadDirty = false;

  int dmpTag = (int)idData(11);
  if (dmpTag) {
    if (theDamping == 0 || theDamping->getClassTag() != dmpTag) {
      if (theDamping) delete theDamping;
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) return -1;
    }
    theDamping->setDbTag((int)idData(12));
    if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0)
      return -1;
  } else if (theDamping) {
    delete theDamping;
    theDamping = 0;
  }

  return 0;
}

int
LadrunoForceBeamColumn2d::getInitialFlexibility(Matrix &fe, bool includeHinges)
{
  fe.Zero();
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();

    Matrix fb(workArea, order, NEBD);
    double xL = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    const Matrix &fSec = sections[i]->getInitialFlexibility();
    fb.Zero();
    double tmp;
    for (int ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        if (!axialRelease) {
          for (int jj = 0; jj < order; jj++)
            fb(jj, 0) += fSec(jj, ii) * wtL;
        }
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < order; jj++) {
          tmp = fSec(jj, ii) * wtL;
          fb(jj, 1) += xL1 * tmp;
          fb(jj, 2) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < order; jj++) {
          tmp = oneOverL * fSec(jj, ii) * wtL;
          fb(jj, 1) += tmp;
          fb(jj, 2) += tmp;
        }
        break;
      default:
        break;
      }
    }
    for (int ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        if (axialRelease) {
          fe(0, 0) += fSec(ii, ii) * wtL;
        } else {
          for (int jj = 0; jj < NEBD; jj++)
            fe(0, jj) += fb(ii, jj);
        }
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < NEBD; jj++) {
          tmp = fb(ii, jj);
          fe(1, jj) += xL1 * tmp;
          fe(2, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < NEBD; jj++) {
          tmp = oneOverL * fb(ii, jj);
          fe(1, jj) += tmp;
          fe(2, jj) += tmp;
        }
        break;
      default:
        break;
      }
    }
  }
  if (GAvy > 0.0) {
    double fsY = oneOverL / GAvy;
    fe(1, 1) += fsY;
    fe(1, 2) += fsY;
    fe(2, 1) += fsY;
    fe(2, 2) += fsY;
  }
  if (includeHinges && hingeOn) {
    double betaI = 1.0 - hingeXiI, gammaI = -hingeXiI;
    double betaJ = hingeXiJ - 1.0, gammaJ = hingeXiJ;
    double dummyDeriv = 0.0;
    double fPM = (pmP0 > 0.0) ? lfbc2dEvalPM(Se(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dummyDeriv) : 1.0;
    if (theHinge != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHinge) == 0) {
      double k0I = fPM * theHinge->getInitialTangent();
      if (k0I > 0.0) {
        double f0I = 1.0 / k0I;
        fe(1, 1) += betaI * betaI * f0I;
        fe(1, 2) += betaI * gammaI * f0I;
        fe(2, 1) += gammaI * betaI * f0I;
        fe(2, 2) += gammaI * gammaI * f0I;
      }
    }
    if (theHingeJ != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) == 0) {
      double k0J = fPM * theHingeJ->getInitialTangent();
      if (k0J > 0.0) {
        double f0J = 1.0 / k0J;
        fe(1, 1) += betaJ * betaJ * f0J;
        fe(1, 2) += betaJ * gammaJ * f0J;
        fe(2, 1) += gammaJ * betaJ * f0J;
        fe(2, 2) += gammaJ * gammaJ * f0J;
      }
    }
  }
  return 0;
}

int
LadrunoForceBeamColumn2d::getInitialDeformations(Vector &v0)
{
  v0.Zero();
  if (numEleLoads < 1)
    return 0;

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();

    double xL = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    static thread_local Vector sp;
    sp.setData(workArea, order);
    sp.Zero();
    this->computeSectionForces(sp, i);

    const Matrix &fse = sections[i]->getInitialFlexibility();
    static thread_local Vector e;
    e.setData(&workArea[order], order);
    e.addMatrixVector(0.0, fse, sp, 1.0);

    double dei, tmp;
    for (int ii = 0; ii < order; ii++) {
      dei = e(ii) * wtL;
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        v0(0) += dei;
        break;
      case SECTION_RESPONSE_MZ:
        v0(1) += xL1 * dei;
        v0(2) += xL * dei;
        break;
      case SECTION_RESPONSE_VY:
        tmp = oneOverL * dei;
        v0(1) += tmp;
        v0(2) += tmp;
        break;
      default:
        break;
      }
    }
  }
  return 0;
}

void
LadrunoForceBeamColumn2d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
  static thread_local Vector ub(NEBD);
  ub = crdTransf->getBasicTrialDisp();
  double L = crdTransf->getInitialLength();

  static thread_local double xi_pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi_pts);

  Matrix ls(numSections, numSections);
  getCBDIinfluenceMatrix(numSections, xi_pts, L, ls);

  Vector kappa(numSections);
  static thread_local Vector vs_loc;
  for (int i = 0; i < numSections; i++) {
    int sectionKey = 0;
    const ID &code = sections[i]->getType();
    int ii;
    for (ii = 0; ii < code.Size(); ii++) {
      if (code(ii) == SECTION_RESPONSE_MZ) {
        sectionKey = ii;
        break;
      }
    }
    if (ii == code.Size()) {
      opserr << "FATAL LadrunoForceBeamColumn2d::compSectionDispls - section does not provide Mz response\n";
      exit(-1);
    }
    vs_loc = sections[i]->getSectionDeformation();
    kappa(i) = vs_loc(sectionKey);
  }

  Vector w(numSections);
  static thread_local Vector xl(NDM), uxb(NDM);
  w.addMatrixVector(0.0, ls, kappa, 1.0);

  for (int i = 0; i < numSections; i++) {
    double xi = xi_pts[i];
    xl(0) = xi * L;
    xl(1) = 0;
    sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);
    uxb(0) = xi * ub(0);
    uxb(1) = w(i);
    sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
  }
}

void
LadrunoForceBeamColumn2d::Print(OPS_Stream &s, int flag)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nElement: " << this->getTag() << " Type: LadrunoForceBeamColumn2d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
    beamIntegr->Print(s, flag);
    double P = Secommit(0);
    double M1 = Secommit(1);
    double M2 = Secommit(2);
    double L = crdTransf->getInitialLength();
    double V = (M1 + M2) / L;
    double p0[3] = {0.0, 0.0, 0.0};
    if (numEleLoads > 0)
      this->computeReactions(p0);
    s << "\tEnd 1 Forces (P V M): " << -P + p0[0] << " " << V + p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P << " " << -V + p0[2] << " " << M2 << endln;
  }
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoForceBeamColumn2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << "\"" << sections[i]->getTag() << "\", ";
    s << "\"" << sections[numSections - 1]->getTag() << "\"], ";
    s << "\"integration\": ";
    beamIntegr->Print(s, flag);
    s << ", \"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
  }
}

OPS_Stream &operator<<(OPS_Stream &s, LadrunoForceBeamColumn2d &E)
{
  E.Print(s);
  return s;
}

void
LadrunoForceBeamColumn2d::setSectionPointers(int numSec, SectionForceDeformation **secPtrs)
{
  if (numSec > maxNumSections)
    opserr << "Error: LadrunoForceBeamColumn2d::setSectionPointers -- max number of sections exceeded\n";

  numSections = numSec;
  if (secPtrs == 0)
    opserr << "Error: LadrunoForceBeamColumn2d::setSectionPointers -- invalid section pointer\n";

  sections = new SectionForceDeformation *[numSections];
  for (int i = 0; i < numSections; i++) {
    sections[i] = secPtrs[i]->getCopy();
  }

  fs = new Matrix[numSections];
  vs = new Vector[numSections];
  Ssr = new Vector[numSections];
  vscommit = new Vector[numSections];
  fscommit = new Matrix[numSections];
  Ssrcommit = new Vector[numSections];
}

int
LadrunoForceBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
  static thread_local Vector v1(3);
  static thread_local Vector v2(3);
  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);
  return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response *
LadrunoForceBeamColumn2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoForceBeamColumn2d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  if (LadrunoResp::isGlobalForce(argv[0])) {
    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Mz_2");
    theResponse = new ElementResponse(this, 1, theVector);
  }
  else if (LadrunoResp::isLocalForce(argv[0])) {
    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "V_1");
    output.tag("ResponseType", "M_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "V_2");
    output.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 2, theVector);
  }
  else if (LadrunoResp::isBasicForce(argv[0])) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "M_1");
    output.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 7, Vector(3));
  }
  else if (theDamping && (strcmp(argv[0], "globalDampingForce") == 0 || strcmp(argv[0], "globalDampingForces") == 0 ||
                          strcmp(argv[0], "dampingForce") == 0)) {
    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Mz_2");
    theResponse = new ElementResponse(this, 21, theVector);
  }
  else if (theDamping && (strcmp(argv[0], "localDampingForce") == 0 || strcmp(argv[0], "localDampingForces") == 0)) {
    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "V_1");
    output.tag("ResponseType", "M_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "V_2");
    output.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 22, theVector);
  }
  else if (theDamping && (strcmp(argv[0], "basicDampingForce") == 0 || strcmp(argv[0], "basicDampingForces") == 0)) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "M_1");
    output.tag("ResponseType", "M_2");
    theResponse = new ElementResponse(this, 23, Vector(3));
  }
  else if (LadrunoResp::is(argv[0], "basicStiffness")) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "M1");
    output.tag("ResponseType", "M2");
    theResponse = new ElementResponse(this, 19, Matrix(3, 3));
  }
  else if (LadrunoResp::isBasicDef(argv[0])) {
    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "theta_1");
    output.tag("ResponseType", "theta_2");
    theResponse = new ElementResponse(this, 3, Vector(3));
  }
  else if (LadrunoResp::isPlasticDef(argv[0])) {
    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaP_1");
    output.tag("ResponseType", "thetaP_2");
    theResponse = new ElementResponse(this, 4, Vector(3));
  }
  else if (strcmp(argv[0], "inflectionPoint") == 0) {
    output.tag("ResponseType", "inflectionPoint");
    theResponse = new ElementResponse(this, 5, 0.0);
  }
  else if (strcmp(argv[0], "tangentDrift") == 0) {
    theResponse = new ElementResponse(this, 6, Vector(2));
  }
  else if (strcmp(argv[0], "integrationPoints") == 0) {
    theResponse = new ElementResponse(this, 10, Vector(numSections));
  }
  else if (strcmp(argv[0], "integrationWeights") == 0) {
    theResponse = new ElementResponse(this, 11, Vector(numSections));
  }
  else if (strcmp(argv[0], "sectionTags") == 0) {
    theResponse = new ElementResponse(this, 110, ID(numSections));
  }
  else if (strcmp(argv[0], "sectionDisplacements") == 0) {
    theResponse = new ElementResponse(this, 111, Matrix(numSections, 3));
  }
  else if (strcmp(argv[0], "cbdiDisplacements") == 0) {
    theResponse = new ElementResponse(this, 112, Matrix(1, 3));
  }
  else if (LadrunoResp::isRayleighForce(argv[0]) || strcmp(argv[0], "dampingForces") == 0) {
    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Mz_2");
    theResponse = new ElementResponse(this, 13, theVector);
  }
  else if (strcmp(argv[0], "sectionX") == 0) {
    if (argc > 2) {
      float sectionLoc = (float)atof(argv[1]);
      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      sectionLoc /= (float)L;
      float minDistance = (float)fabs(xi[0] - sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = (float)fabs(xi[i] - sectionLoc);
          sectionNum = i;
        }
      }
      output.tag("GaussPointOutput");
      output.attr("number", sectionNum + 1);
      output.attr("eta", xi[sectionNum] * L);
      if (strcmp(argv[2], "dsdh") != 0) {
        theResponse = sections[sectionNum]->setResponse(&argv[2], argc - 2, output);
      } else {
        int order = sections[sectionNum]->getOrder();
        theResponse = new ElementResponse(this, 76, Vector(order));
        Information &info = theResponse->getInformation();
        info.theInt = sectionNum;
      }
      output.endTag();
    }
  }
  else if (strcmp(argv[0], "section") == 0) {
    if (argc > 1) {
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {
        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);
        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", xi[sectionNum - 1] * L);
        if (strcmp(argv[2], "dsdh") != 0) {
          theResponse = sections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);
        } else {
          int order = sections[sectionNum - 1]->getOrder();
          theResponse = new ElementResponse(this, 76, Vector(order));
          Information &info = theResponse->getInformation();
          info.theInt = sectionNum;
        }
        output.endTag();
      } else if (sectionNum == 0) {
        CompositeResponse *theCResponse = new CompositeResponse();
        int numResponse = 0;
        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);
        for (int i = 0; i < numSections; i++) {
          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", xi[i] * L);
          Response *theSectionResponse = sections[i]->setResponse(&argv[1], argc - 1, output);
          if (theSectionResponse != 0)
            numResponse = theCResponse->addResponse(theSectionResponse);
          output.endTag();
        }
        if (numResponse == 0)
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }
  else if (LadrunoResp::isEnergy(argv[0])) {
    output.tag("ResponseType", "Energy");
    theResponse = new ElementResponse(this, 14, 0.0);
  }
  else if (LadrunoResp::isLocalAxes(argv[0])) {
    theResponse = new ElementResponse(this, 30, Vector(9));
  }
  else if (hingeOn && argc > 1 &&
           (strcmp(argv[0], "hinge") == 0 || strcmp(argv[0], "hingeZ") == 0)) {
    UniaxialMaterial *hm = (theHinge != 0) ? theHinge : theHingeJ;
    if (hm != 0)
      theResponse = hm->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHinge != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hingeZI") == 0 || strcmp(argv[0], "hinge1") == 0)) {
    theResponse = theHinge->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeJ != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hingeZJ") == 0 || strcmp(argv[0], "hinge2") == 0)) {
    theResponse = theHingeJ->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && (theHinge != 0 || theHingeJ != 0) && LadrunoResp::is(argv[0], "hingeSummary")) {
    output.tag("ResponseType", "theta_i");
    output.tag("ResponseType", "M_i");
    output.tag("ResponseType", "theta_j");
    output.tag("ResponseType", "M_j");
    theResponse = new ElementResponse(this, 31, Vector(4));
  }
  else if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
    output.tag("ResponseType", "implexError");
    theResponse = new ElementResponse(this, 32, 0.0);
  }
  else if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
    output.tag("ResponseType", "avgImplexError");
    theResponse = new ElementResponse(this, 33, 0.0);
  }
  else if (strcmp(argv[0], "hingeJump") == 0 || strcmp(argv[0], "hingeJumpI") == 0) {
    output.tag("ResponseType", "hingeJumpI");
    theResponse = new ElementResponse(this, 34, 0.0);
  }
  else if (strcmp(argv[0], "hingeJumpJ") == 0) {
    output.tag("ResponseType", "hingeJumpJ");
    theResponse = new ElementResponse(this, 35, 0.0);
  }

  if (theResponse == 0)
    theResponse = crdTransf->setResponse(argv, argc, output);

  output.endTag(); // ElementOutput

  if (theResponse == 0)
    return Element::setResponse(argv, argc, output);
  return theResponse;
}

int
LadrunoForceBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  static thread_local Vector vp(3);
  static thread_local Matrix fe(3, 3);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 30) {
    static thread_local Vector la(9);
    static thread_local Vector vx(3), vy(3), vz(3);
    crdTransf->getLocalAxes(vx, vy, vz);
    for (int i = 0; i < 3; i++) {
      la(i) = vx(i);
      la(i + 3) = vy(i);
      la(i + 6) = vz(i);
    }
    return eleInfo.setVector(la);
  }

  else if (responseID == 31) {
    static thread_local Vector hs(4);
    hs.Zero();
    double dummyDeriv = 0.0;
    if (hingeOn && theHinge != 0) {
      double fPMI = (dynamic_cast<LadrunoCohesiveHinge*>(theHinge) == 0 && pmP0 > 0.0)
                        ? lfbc2dEvalPM(Se(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dummyDeriv) : 1.0;
      hs(0) = hingeJump;
      hs(1) = fPMI * theHinge->getStress();
    }
    if (hingeOn && theHingeJ != 0) {
      double fPMJ = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) == 0 && pmP0 > 0.0)
                        ? lfbc2dEvalPM(Se(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dummyDeriv) : 1.0;
      hs(2) = hingeJumpJ;
      hs(3) = fPMJ * theHingeJ->getStress();
    }
    return eleInfo.setVector(hs);
  }

  else if (responseID == 32 || responseID == 33) {
    double e1 = 0.0, e2 = 0.0;
    int cnt = 0;
    if (lfbc2dGetHingeImplexError(theHinge, e1))
      cnt++;
    if (lfbc2dGetHingeImplexError(theHingeJ, e2))
      cnt++;
    if (cnt > 0) {
      double val = (responseID == 32) ? ((e1 > e2) ? e1 : e2) : ((e1 + e2) / cnt);
      return eleInfo.setDouble(val);
    }
    return eleInfo.setDouble(useImplex ? eleImplexErr : 0.0);
  }

  else if (responseID == 34) {
    return eleInfo.setDouble(hingeOn ? hingeJump : 0.0);
  }

  else if (responseID == 35) {
    return eleInfo.setDouble(hingeOn ? hingeJumpJ : 0.0);
  }

  else if (responseID == 13) {
    static thread_local Vector damp(6);
    damp.Zero();
    damp.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    return eleInfo.setVector(damp);
  }

  else if (responseID == 2) {
    double p0[3] = {0.0, 0.0, 0.0};
    if (numEleLoads > 0)
      this->computeReactions(p0);
    theVector(3) = Se(0);
    theVector(0) = -Se(0) + p0[0];
    theVector(2) = Se(1);
    theVector(5) = Se(2);
    double V = (Se(1) + Se(2)) / crdTransf->getInitialLength();
    theVector(1) = V + p0[1];
    theVector(4) = -V + p0[2];
    return eleInfo.setVector(theVector);
  }

  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    static thread_local Vector v0(3);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return eleInfo.setVector(vp);
  }

  else if (responseID == 5) {
    double LI = 0.0;
    if (fabs(Se(1) + Se(2)) > DBL_EPSILON) {
      double L = crdTransf->getInitialLength();
      LI = Se(1) / (Se(1) + Se(2)) * L;
    }
    return eleInfo.setDouble(LI);
  }

  else if (responseID == 6) {
    double d2 = 0.0, d3 = 0.0;
    double L = crdTransf->getInitialLength();
    double LI = 0.0;
    if (fabs(Se(1) + Se(2)) > DBL_EPSILON)
      LI = Se(1) / (Se(1) + Se(2)) * L;

    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    for (int i = 0; i < numSections; i++) {
      double x = pts[i] * L;
      if (x > LI) continue;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      for (int j = 0; j < order; j++)
        if (type(j) == SECTION_RESPONSE_MZ)
          kappa += vs[i](j);
      double b = -LI + x;
      d2 += (wts[i] * L) * kappa * b;
    }
    for (int i = numSections - 1; i >= 0; i--) {
      double x = pts[i] * L;
      if (x < LI) continue;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      for (int j = 0; j < order; j++)
        if (type(j) == SECTION_RESPONSE_MZ)
          kappa += vs[i](j);
      double b = x - LI;
      d3 += (wts[i] * L) * kappa * b;
    }
    static thread_local Vector d(2);
    d(0) = d2;
    d(1) = d3;
    return eleInfo.setVector(d);
  }

  else if (responseID == 7)
    return eleInfo.setVector(Se);

  else if (responseID == 21)
    return eleInfo.setVector(this->getDampingForce());

  else if (responseID == 22) {
    Vector Sd(NEBD);
    Sd = theDamping->getDampingForce();
    theVector(3) = Sd(0);
    theVector(0) = -Sd(0);
    theVector(2) = Sd(1);
    theVector(5) = Sd(2);
    double V = (Sd(1) + Sd(2)) / crdTransf->getInitialLength();
    theVector(1) = V;
    theVector(4) = -V;
    return eleInfo.setVector(theVector);
  }

  else if (responseID == 23)
    return eleInfo.setVector(theDamping->getDampingForce());

  else if (responseID == 19)
    return eleInfo.setMatrix(kv);

  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i] * L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i] * L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = sections[i]->getTag();
    return eleInfo.setID(tags);
  }

  else if (responseID == 111) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    Vector kappa(numSections);
    for (int i = 0; i < numSections; i++) {
      const ID &code = sections[i]->getType();
      const Vector &e = sections[i]->getSectionDeformation();
      int order = sections[i]->getOrder();
      for (int j = 0; j < order; j++)
        if (code(j) == SECTION_RESPONSE_MZ)
          kappa(i) += e(j);
    }
    Vector dispsy(numSections);
    dispsy.addMatrixVector(0.0, ls, kappa, 1.0);
    static thread_local Vector uxb(2);
    static thread_local Vector uxg(2);
    Matrix disps(numSections, 3);
    vp = crdTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i] * vp(0);
      uxb(1) = dispsy(i);
      uxg = crdTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = 0.0;
    }
    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 112) {
    double L = crdTransf->getInitialLength();
    double ipts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, ipts);
    double pts[1];
    pts[0] = eleInfo.theDouble;
    Matrix ls(1, numSections);
    getCBDIinfluenceMatrix(1, pts, numSections, ipts, L, ls);
    Vector kappa(numSections);
    for (int i = 0; i < numSections; i++) {
      const ID &code = sections[i]->getType();
      const Vector &e = sections[i]->getSectionDeformation();
      int order = sections[i]->getOrder();
      for (int j = 0; j < order; j++)
        if (code(j) == SECTION_RESPONSE_MZ)
          kappa(i) += e(j);
    }
    Vector dispsy(1);
    dispsy.addMatrixVector(0.0, ls, kappa, 1.0);
    static thread_local Vector uxb(2);
    static thread_local Vector uxg(2);
    Matrix disps(1, 3);
    vp = crdTransf->getBasicTrialDisp();
    uxb(0) = pts[0] * vp(0);
    uxb(1) = dispsy(0);
    uxg = crdTransf->getPointGlobalDisplFromBasic(pts[0], uxb);
    disps(0, 0) = uxg(0);
    disps(0, 1) = uxg(1);
    disps(0, 2) = 0.0;
    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 14) {
    double wt[maxNumSections];
    double L = crdTransf->getInitialLength();
    beamIntegr->getSectionWeights(numSections, L, wt);
    double energy = 0.0;
    for (int i = 0; i < numSections; i++)
      energy += sections[i]->getEnergy() * wt[i] * L;
    double dummyDeriv = 0.0;
    if (hingeOn && theHinge != 0) {
      double fPMI = (dynamic_cast<LadrunoCohesiveHinge*>(theHinge) == 0 && pmP0 > 0.0)
                        ? lfbc2dEvalPM(Se(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dummyDeriv) : 1.0;
      energy += fPMI * theHinge->getEnergy();
    }
    if (hingeOn && theHingeJ != 0) {
      double fPMJ = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeJ) == 0 && pmP0 > 0.0)
                        ? lfbc2dEvalPM(Se(0), pmP0, pmAlpha, pmP0t, pmPbal, pmMuBal, dummyDeriv) : 1.0;
      energy += fPMJ * theHingeJ->getEnergy();
    }
    return eleInfo.setDouble(energy);
  }

  else
    return Element::getResponse(responseID, eleInfo);
}

double
LadrunoForceBeamColumn2d::getCharacteristicLength(void)
{
  if (lchMode == 1)
    return crdTransf ? crdTransf->getInitialLength() : Element::getCharacteristicLength();
  if (lchMode == 2 && userLch > 0.0)
    return userLch;
  if (current_section_lch > 0.0)
    return current_section_lch;
  return Element::getCharacteristicLength();
}

int
LadrunoForceBeamColumn2d::getResponseSensitivity(int responseID, int gradNumber, Information &eleInfo)
{
  if (responseID == 3) {
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }
  else if (responseID == 7) {
    static thread_local Vector dqdh(3);
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);
    return eleInfo.setVector(dqdh);
  }
  else if (responseID == 76) {
    int sectionNum = eleInfo.theInt;
    int order = sections[sectionNum - 1]->getOrder();
    Vector dsdh(order);
    dsdh.Zero();
    if (numEleLoads > 0)
      this->computeSectionForceSensitivity(dsdh, sectionNum - 1, gradNumber);

    static thread_local Vector dqdh(3);
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    const ID &code = sections[sectionNum - 1]->getType();
    double xL = pts[sectionNum - 1];
    double xL1 = xL - 1.0;

    for (int ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        dsdh(ii) += dqdh(0);
        break;
      case SECTION_RESPONSE_MZ:
        dsdh(ii) += xL1 * dqdh(1) + xL * dqdh(2);
        break;
      case SECTION_RESPONSE_VY:
        dsdh(ii) += oneOverL * (dqdh(1) + dqdh(2));
        break;
      default:
        break;
      }
    }

    double dLdh = crdTransf->getdLdh();
    double d1oLdh = crdTransf->getd1overLdh();
    double dptsdh[maxNumSections];
    beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);
    double dxLdh = dptsdh[sectionNum - 1];

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        dsdh(j) += dxLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_VY:
        dsdh(j) += d1oLdh * (Se(1) + Se(2));
        break;
      default:
        break;
      }
    }
    return eleInfo.setVector(dsdh);
  }
  else if (responseID == 4) {
    static thread_local Vector dvpdh(3);
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    dvpdh = dvdh;
    static thread_local Matrix fe(3, 3);
    this->getInitialFlexibility(fe);
    const Vector &dqdh = this->computedqdh(gradNumber);
    dvpdh.addMatrixVector(1.0, fe, dqdh, -1.0);
    static thread_local Matrix fek(3, 3);
    fek.addMatrixProduct(0.0, fe, kv, 1.0);
    dvpdh.addMatrixVector(1.0, fek, dvdh, -1.0);
    const Matrix &dfedh = this->computedfedh(gradNumber);
    dvpdh.addMatrixVector(1.0, dfedh, Se, -1.0);
    return eleInfo.setVector(dvpdh);
  }
  return -1;
}

int
LadrunoForceBeamColumn2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "lch") == 0 || strcmp(argv[0], "Lch") == 0 ||
      strcmp(argv[0], "characteristicLength") == 0) {
    param.setValue(this->getCharacteristicLength());
    return param.addObject(3, this);
  }
  if (strcmp(argv[0], "GAy") == 0 || strcmp(argv[0], "GAvy") == 0) {
    param.setValue(GAvy);
    return param.addObject(4, this);
  }

  if (strstr(argv[0], "damp") != 0) {
    if (argc < 2 || !theDamping)
      return -1;
    return theDamping->setParameter(&argv[1], argc - 1, param);
  }

  if (hingeOn && argc > 1) {
    if (strcmp(argv[0], "hinge") == 0 || strcmp(argv[0], "hingeZ") == 0 ||
        strcmp(argv[0], "hingeIJ") == 0 || strcmp(argv[0], "hingeZIJ") == 0) {
      int rI = (theHinge != 0) ? theHinge->setParameter(&argv[1], argc - 1, param) : -1;
      int rJ = (theHingeJ != 0) ? theHingeJ->setParameter(&argv[1], argc - 1, param) : -1;
      int res = (rI != -1) ? rI : rJ;
      if (res != -1) param.addObject(2, this);
      return res;
    }
    if (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hingeZI") == 0 || strcmp(argv[0], "hinge1") == 0) {
      int res = (theHinge != 0) ? theHinge->setParameter(&argv[1], argc - 1, param) : -1;
      if (res != -1) param.addObject(2, this);
      return res;
    }
    if (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hingeZJ") == 0 || strcmp(argv[0], "hinge2") == 0) {
      int res = (theHingeJ != 0) ? theHingeJ->setParameter(&argv[1], argc - 1, param) : -1;
      if (res != -1) param.addObject(2, this);
      return res;
    }
  }

  if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = (float)atof(argv[1]);
      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      sectionLoc /= (float)L;
      float minDistance = (float)fabs(xi[0] - sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = (float)fabs(xi[i] - sectionLoc);
          sectionNum = i;
        }
      }
      int res = sections[sectionNum]->setParameter(&argv[2], argc - 2, param);
      if (res != -1)
        param.addObject(2, this);
      return res;
    }
  }

  if (strcmp(argv[0], "allSections") == 0) {
    if (argc < 2)
      return -1;
    for (int i = 0; i < numSections; i++) {
      int ok = sections[i]->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }
    if (result != -1)
      param.addObject(2, this);
    return result;
  }

  if (strstr(argv[0], "section") != 0) {
    if (argc < 3)
      return -1;
    int sectionNum = atoi(argv[1]);
    if (sectionNum > 0 && sectionNum <= numSections) {
      int res = sections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);
      if (res != -1)
        param.addObject(2, this);
      return res;
    }
    return -1;
  }

  if (strstr(argv[0], "integration") != 0) {
    if (argc < 2)
      return -1;
    int res = beamIntegr->setParameter(&argv[1], argc - 1, param);
    if (res != -1)
      param.addObject(2, this);
    return res;
  }

  for (int i = 0; i < numSections; i++) {
    int ok = sections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  int ok = beamIntegr->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;
  if (hingeOn) {
    if (theHinge != 0) {
      ok = theHinge->setParameter(argv, argc, param);
      if (ok != -1) result = ok;
    }
    if (theHingeJ != 0) {
      ok = theHingeJ->setParameter(argv, argc, param);
      if (ok != -1) result = ok;
    }
  }
  if (result != -1)
    param.addObject(2, this);
  return result;
}

int
LadrunoForceBeamColumn2d::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else if (parameterID == 2 || parameterID == 3 || parameterID == 4) {
    if (parameterID == 3 && info.theDouble > 0.0) {
      userLch = info.theDouble;
      lchMode = 2;
      current_section_lch = userLch;
    }
    if (parameterID == 4 && info.theDouble >= 0.0) {
      GAvy = info.theDouble;
    }
    if (Ki != 0) {
      delete Ki;
      Ki = 0;
    }
    if (theNodes[0] != 0 && theNodes[1] != 0) {
      ops_TheActiveElement = this;
      if (initialFlag != 0) {
        double L = crdTransf->getInitialLength();
        double wt[maxNumSections];
        beamIntegr->getSectionWeights(numSections, L, wt);
        for (int i = 0; i < numSections; i++) {
          if (lchMode == 1)
            current_section_lch = L;
          else if (lchMode == 2)
            current_section_lch = userLch;
          else
            current_section_lch = wt[i] * L;
          sections[i]->setTrialSectionDeformation(vs[i]);
          Ssr[i] = sections[i]->getStressResultant();
          fs[i] = sections[i]->getSectionFlexibility();
        }
        current_section_lch = (lchMode == 2) ? userLch : L;
        initialFlag = 3;
      }
      int rc = this->update();
      ops_TheActiveElement = 0;
      return rc;
    }
    return 0;
  }
  return -1;
}

int
LadrunoForceBeamColumn2d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}

const Matrix &
LadrunoForceBeamColumn2d::getInitialStiffSensitivity(int gradNumber)
{
  theMatrix.Zero();
  return theMatrix;
}

const Matrix &
LadrunoForceBeamColumn2d::getMassSensitivity(int gradNumber)
{
  theMatrix.Zero();
  if (rho == 0.0 || parameterID != 1)
    return theMatrix;

  double L = crdTransf->getInitialLength();
  if (cMass == 0) {
    theMatrix(0, 0) = theMatrix(1, 1) = theMatrix(3, 3) = theMatrix(4, 4) = 0.5 * L;
  } else {
    static thread_local Matrix ml(6, 6);
    ml.Zero();
    double m = L / 420.0;
    ml(0, 0) = ml(3, 3) = m * 140.0;
    ml(0, 3) = ml(3, 0) = m * 70.0;
    ml(1, 1) = ml(4, 4) = m * 156.0;
    ml(1, 4) = ml(4, 1) = m * 54.0;
    ml(2, 2) = ml(5, 5) = m * 4.0 * L * L;
    ml(2, 5) = ml(5, 2) = -m * 3.0 * L * L;
    ml(1, 2) = ml(2, 1) = m * 22.0 * L;
    ml(4, 5) = ml(5, 4) = -ml(1, 2);
    ml(1, 5) = ml(5, 1) = -m * 13.0 * L;
    ml(2, 4) = ml(4, 2) = -ml(1, 5);
    theMatrix = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  return theMatrix;
}

const Vector &
LadrunoForceBeamColumn2d::getResistingForceSensitivity(int gradNumber)
{
  static thread_local Vector dqdh(3);
  dqdh = this->computedqdh(gradNumber);

  double dp0dh[3] = {0.0, 0.0, 0.0};
  this->computeReactionSensitivity(dp0dh, gradNumber);
  Vector dp0dhVec(dp0dh, 3);

  static thread_local Vector P(6);
  P.Zero();

  if (crdTransf->isShapeSensitivity()) {
    P = crdTransf->getGlobalResistingForceShapeSensitivity(Se, dp0dhVec, gradNumber);
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);
  }
  P += crdTransf->getGlobalResistingForce(dqdh, dp0dhVec);
  return P;
}

int
LadrunoForceBeamColumn2d::commitSensitivity(int gradNumber, int numGrads)
{
  int err = 0;
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wts);

  double dLdh = crdTransf->getdLdh();
  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);
  double d1oLdh = crdTransf->getd1overLdh();

  static thread_local Vector dqdh(3);
  dqdh = this->computedqdh(gradNumber);

  const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  dqdh.addMatrixVector(1.0, kv, dvdh, 1.0);

  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    double xL = pts[i];
    double xL1 = xL - 1.0;
    double dxLdh = dptsdh[i];

    Vector ds(workArea, order);
    ds.Zero();
    if (numEleLoads > 0)
      this->computeSectionForceSensitivity(ds, i, gradNumber);

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
        ds(j) += dqdh(0);
        break;
      case SECTION_RESPONSE_MZ:
        ds(j) += xL1 * dqdh(1) + xL * dqdh(2);
        break;
      case SECTION_RESPONSE_VY:
        ds(j) += oneOverL * (dqdh(1) + dqdh(2));
        break;
      default:
        break;
      }
    }

    const Vector &dsdh = sections[i]->getStressResultantSensitivity(gradNumber, true);
    ds -= dsdh;

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        ds(j) += dxLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_VY:
        ds(j) += d1oLdh * (Se(1) + Se(2));
        break;
      default:
        break;
      }
    }

    Vector de(&workArea[order], order);
    const Matrix &fSec = sections[i]->getSectionFlexibility();
    de.addMatrixVector(0.0, fSec, ds, 1.0);
    err += sections[i]->commitSensitivity(de, gradNumber, numGrads);
  }
  return err;
}

const Vector &
LadrunoForceBeamColumn2d::computedqdh(int gradNumber)
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, pts);
  double wts[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wts);

  double dLdh = crdTransf->getdLdh();
  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);
  double dwtsdh[maxNumSections];
  beamIntegr->getWeightsDeriv(numSections, L, dLdh, dwtsdh);
  double d1oLdh = crdTransf->getd1overLdh();

  static thread_local Vector dvdh(3);
  dvdh.Zero();

  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    double xL = pts[i];
    double xL1 = xL - 1.0;
    double wtL = wts[i] * L;
    double dxLdh = dptsdh[i];
    double dwtLdh = wts[i] * dLdh + dwtsdh[i] * L;

    Vector dsdh(&workArea[order], order);
    dsdh = sections[i]->getStressResultantSensitivity(gradNumber, true);

    Vector dspdh(&workArea[2 * order], order);
    dspdh.Zero();
    if (numEleLoads > 0)
      this->computeSectionForceSensitivity(dspdh, i, gradNumber);
    dsdh.addVector(1.0, dspdh, -1.0);

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        dsdh(j) -= dxLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_VY:
        dsdh(j) -= d1oLdh * (Se(1) + Se(2));
        break;
      default:
        break;
      }
    }

    Vector dedh(workArea, order);
    const Matrix &fSec = sections[i]->getSectionFlexibility();
    dedh.addMatrixVector(0.0, fSec, dsdh, 1.0);

    for (int j = 0; j < order; j++) {
      double dei = dedh(j) * wtL;
      switch (code(j)) {
      case SECTION_RESPONSE_P:
        dvdh(0) += dei;
        break;
      case SECTION_RESPONSE_MZ:
        dvdh(1) += xL1 * dei;
        dvdh(2) += xL * dei;
        break;
      case SECTION_RESPONSE_VY:
        dei = oneOverL * dei;
        dvdh(1) += dei;
        dvdh(2) += dei;
        break;
      default:
        break;
      }
    }

    const Vector &e = vs[i];
    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
        dvdh(0) -= e(j) * dwtLdh;
        break;
      case SECTION_RESPONSE_MZ:
        dvdh(1) -= xL1 * e(j) * dwtLdh;
        dvdh(2) -= xL * e(j) * dwtLdh;
        dvdh(1) -= dxLdh * e(j) * wtL;
        dvdh(2) -= dxLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_VY:
        dvdh(1) -= oneOverL * e(j) * dwtLdh;
        dvdh(2) -= oneOverL * e(j) * dwtLdh;
        dvdh(1) -= d1oLdh * e(j) * wtL;
        dvdh(2) -= d1oLdh * e(j) * wtL;
        break;
      default:
        break;
      }
    }
  }

  static thread_local Vector dqdh(3);
  dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
  return dqdh;
}

const Matrix &
LadrunoForceBeamColumn2d::computedfedh(int gradNumber)
{
  static thread_local Matrix dfedh(3, 3);
  dfedh.Zero();

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    Matrix fb(workArea, order, NEBD);
    double xL = xi[i];
    double xL1 = xL - 1.0;
    double wtL = wt[i] * L;

    const Matrix &dfsdh = sections[i]->getInitialFlexibilitySensitivity(gradNumber);
    fb.Zero();
    double tmp;
    for (int ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < order; jj++)
          fb(jj, 0) += dfsdh(jj, ii) * wtL;
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < order; jj++) {
          tmp = dfsdh(jj, ii) * wtL;
          fb(jj, 1) += xL1 * tmp;
          fb(jj, 2) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < order; jj++) {
          tmp = oneOverL * dfsdh(jj, ii) * wtL;
          fb(jj, 1) += tmp;
          fb(jj, 2) += tmp;
        }
        break;
      default:
        break;
      }
    }
    for (int ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        for (int jj = 0; jj < NEBD; jj++)
          dfedh(0, jj) += fb(ii, jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (int jj = 0; jj < NEBD; jj++) {
          tmp = fb(ii, jj);
          dfedh(1, jj) += xL1 * tmp;
          dfedh(2, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int jj = 0; jj < NEBD; jj++) {
          tmp = oneOverL * fb(ii, jj);
          dfedh(1, jj) += tmp;
          dfedh(2, jj) += tmp;
        }
        break;
      default:
        break;
      }
    }
  }
  return dfedh;
}

bool
LadrunoForceBeamColumn2d::ladrunoThreadSafeUpdate(void) const
{
  if (theDamping != 0 || crdTransf == 0 || beamIntegr == 0 || sections == 0)
    return false;
  int ctTag = crdTransf->getClassTag();
  if (ctTag != CRDTR_TAG_LinearCrdTransf2d &&
      ctTag != CRDTR_TAG_PDeltaCrdTransf2d &&
      ctTag != CRDTR_TAG_CorotCrdTransf2d)
    return false;
  int biTag = beamIntegr->getClassTag();
  if (biTag != BEAM_INTEGRATION_TAG_Lobatto &&
      biTag != BEAM_INTEGRATION_TAG_Legendre &&
      biTag != BEAM_INTEGRATION_TAG_Radau &&
      biTag != BEAM_INTEGRATION_TAG_HingeRadau &&
      biTag != BEAM_INTEGRATION_TAG_HingeRadauTwo &&
      biTag != BEAM_INTEGRATION_TAG_UserDefined &&
      biTag != BEAM_INTEGRATION_TAG_FixedLocation)
    return false;
  for (int i = 0; i < numSections; i++) {
    if (sections[i] == 0 || !sections[i]->ladrunoThreadSafeUpdate())
      return false;
  }
  if (theHinge != 0 && !theHinge->ladrunoThreadSafeUpdate())
    return false;
  if (theHingeJ != 0 && !theHingeJ->ladrunoThreadSafeUpdate())
    return false;
  return true;
}

bool
LadrunoForceBeamColumn2d::ladrunoThreadSafeForm(void) const
{
  if (numEleLoads > 0 || !this->ladrunoThreadSafeUpdate())
    return false;
  for (int i = 0; i < numSections; i++) {
    if (!sections[i]->ladrunoThreadSafeForm())
      return false;
  }
  if (theHinge != 0 && !theHinge->ladrunoThreadSafeForm())
    return false;
  if (theHingeJ != 0 && !theHingeJ->ladrunoThreadSafeForm())
    return false;
  return true;
}
