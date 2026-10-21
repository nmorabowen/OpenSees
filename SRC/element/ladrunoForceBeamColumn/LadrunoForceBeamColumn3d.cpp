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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <cmath>
#include <vector>

#include <Information.h>
#include <Parameter.h>
#include <LadrunoForceBeamColumn3d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <elementAPI.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <ElementIter.h>
#include <LadrunoThreads.h>
#include <LadrunoResponseTokens.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <LadrunoCohesiveHinge.h>
#include <LadrunoCohesiveHingeBiaxial.h>
#include <Bilin.h>
#include <ModIMKPeakOriented.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <CompositeSimpsonBeamIntegration.h>

#define DefaultLoverGJ 1.0e-10

thread_local Matrix LadrunoForceBeamColumn3d::theMatrix(12, 12);
thread_local Vector LadrunoForceBeamColumn3d::theVector(12);
thread_local double LadrunoForceBeamColumn3d::workArea[400];

thread_local Vector LadrunoForceBeamColumn3d::vsSubdivide[maxNumSections];
thread_local Matrix LadrunoForceBeamColumn3d::fsSubdivide[maxNumSections];
thread_local Vector LadrunoForceBeamColumn3d::SsrSubdivide[maxNumSections];

static inline double
lfbc3dScalePM(double P, double P0, double alpha)
{
  if (P0 <= 0.0)
    return 1.0;
  double n = P / P0;
  double f = 1.0 - alpha * n - (1.0 - alpha) * n * n;
  if (f < 0.05) f = 0.05;
  else if (f > 2.5) f = 2.5;
  return f;
}

static inline bool
lfbc3dHingeIsImplex(UniaxialMaterial *m)
{
  if (m == 0)
    return false;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(m))
    return ch->isImplex();
  if (Bilin *b = dynamic_cast<Bilin*>(m))
    return b->isImplex();
  if (ModIMKPeakOriented *mk = dynamic_cast<ModIMKPeakOriented*>(m))
    return mk->isImplex();
  return false;
}

static inline bool
lfbc3dHingeBiaxialIsImplex(NDMaterial *m)
{
  if (m == 0)
    return false;
  if (LadrunoCohesiveHingeBiaxial *ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(m))
    return ch->isImplex();
  return false;
}

static inline bool
lfbc3dGetHingeImplexError(UniaxialMaterial *m, double &errOut)
{
  errOut = 0.0;
  if (m == 0)
    return false;
  if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(m)) {
    errOut = ch->getImplexError();
    return true;
  }
  if (Bilin *b = dynamic_cast<Bilin*>(m)) {
    errOut = b->getImplexError();
    return true;
  }
  if (ModIMKPeakOriented *mk = dynamic_cast<ModIMKPeakOriented*>(m)) {
    errOut = mk->getImplexError();
    return true;
  }
  return false;
}

static void
lfbc3dFlooredInv4x4(const double Kaa[4][4], const bool active[4], double floor,
                    double invKaa[4][4])
{
  double A[4][4];
  double V[4][4];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      V[i][j] = (i == j) ? 1.0 : 0.0;
      if (active[i] && active[j])
        A[i][j] = 0.5 * (Kaa[i][j] + Kaa[j][i]);
      else
        A[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int sweep = 0; sweep < 15; sweep++) {
    double off = 0.0;
    for (int p = 0; p < 3; p++)
      for (int q = p + 1; q < 4; q++)
        off += fabs(A[p][q]);
    if (off < 1.0e-30)
      break;
    for (int p = 0; p < 3; p++) {
      for (int q = p + 1; q < 4; q++) {
        double apq = A[p][q];
        if (fabs(apq) < 1.0e-30)
          continue;
        double app = A[p][p];
        double aqq = A[q][q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0) ? (1.0 / (tau + sqrt(1.0 + tau * tau)))
                                : (-1.0 / (-tau + sqrt(1.0 + tau * tau)));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = A[q][p] = 0.0;
        for (int r = 0; r < 4; r++) {
          if (r != p && r != q) {
            double arp = A[r][p];
            double arq = A[r][q];
            A[r][p] = A[p][r] = c * arp - s * arq;
            A[r][q] = A[q][r] = s * arp + c * arq;
          }
          double vrp = V[r][p];
          double vrq = V[r][q];
          V[r][p] = c * vrp - s * vrq;
          V[r][q] = s * vrp + c * vrq;
        }
      }
    }
  }

  double invLam[4];
  for (int k = 0; k < 4; k++) {
    if (!active[k]) {
      invLam[k] = 0.0;
    } else {
      double lam = A[k][k];
      if (fabs(lam) < floor)
        lam = (lam < 0.0) ? -floor : floor;
      invLam[k] = 1.0 / lam;
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (!active[i] || !active[j]) {
        invKaa[i][j] = 0.0;
      } else {
        double sum = 0.0;
        for (int k = 0; k < 4; k++)
          sum += V[i][k] * invLam[k] * V[j][k];
        invKaa[i][j] = sum;
      }
    }
  }
}

static bool
lfbc3dTryOptionalInt(int &valOut)
{
  if (OPS_GetNumRemainingInputArgs() < 1)
    return false;
  char buf[64];
  const char *peek = OPS_GetStringFromAll(buf, sizeof(buf));
  if (peek == 0 || peek[0] == '\0')
    return false;
  char *endptr = 0;
  long v = strtol(peek, &endptr, 10);
  if (endptr == peek || *endptr != '\0') {
    OPS_ResetCurrentInputArg(-1);
    return false;
  }
  valOut = (int)v;
  return true;
}

static bool
lfbc3dTryOptionalDouble(double &valOut)
{
  if (OPS_GetNumRemainingInputArgs() < 1)
    return false;
  char buf[64];
  const char *peek = OPS_GetStringFromAll(buf, sizeof(buf));
  if (peek == 0 || peek[0] == '\0' || (peek[0] == '-' && (peek[1] < '0' || peek[1] > '9') && peek[1] != '.')) {
    OPS_ResetCurrentInputArg(-1);
    return false;
  }
  char *endptr = 0;
  double v = strtod(peek, &endptr);
  if (endptr == peek || *endptr != '\0' || !std::isfinite(v)) {
    OPS_ResetCurrentInputArg(-1);
    return false;
  }
  valOut = v;
  return true;
}

void *
OPS_LadrunoForceBeamColumn3d()
{
  int dampingTag = 0;
  Damping *theDamping = 0;
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING insufficient arguments: eleTag iNode jNode transfTag integrationTag\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm != 3 || ndf != 6) {
    opserr << "WARNING LadrunoForceBeamColumn3d: ndm must be 3 and ndf must be 6\n";
    return 0;
  }

  int iData[3];
  int numData = 3;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING LadrunoForceBeamColumn3d: invalid int inputs\n";
    return 0;
  }

  int arg4 = 0, arg5 = 0;
  numData = 1;
  if (OPS_GetIntInput(&numData, &arg4) < 0) {
    opserr << "WARNING LadrunoForceBeamColumn3d: invalid 4th int input\n";
    return 0;
  }
  if (OPS_GetIntInput(&numData, &arg5) < 0) {
    opserr << "WARNING LadrunoForceBeamColumn3d: invalid 5th int input\n";
    return 0;
  }

  bool legacyForm = false;
  int transfTag = arg4;
  int beamIntTag = arg5;
  int legacyNumSecs = 0;
  std::vector<int> legacySecTags;
  BeamIntegration *legacyBi = 0;

  if (OPS_getCrdTransf(arg4) == 0 || OPS_getBeamIntegrationRule(arg5) == 0) {
    if (arg4 > 0 && arg4 <= 20 && OPS_GetNumRemainingInputArgs() >= arg4) {
      legacyForm = true;
      legacyNumSecs = arg4;
      legacySecTags.resize(legacyNumSecs);
      legacySecTags[0] = arg5;
      for (int i = 1; i < legacyNumSecs; i++) {
        int sTag = 0;
        numData = 1;
        if (OPS_GetIntInput(&numData, &sTag) < 0) {
          opserr << "WARNING LadrunoForceBeamColumn3d: invalid legacy section tag\n";
          return 0;
        }
        legacySecTags[i] = sTag;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &transfTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid legacy transfTag\n";
        return 0;
      }
    }
  }

  double mass = 0.0, tol = 1e-12, subFac = 10.0, userLch = 0.0;
  int maxIter = 10, numSub = 4;
  int cMass = 0, lchMode = 0;
  UniaxialMaterial *hingeMatZ = 0;
  UniaxialMaterial *hingeMatY = 0;
  NDMaterial *hingeMatC = 0;
  UniaxialMaterial *hingeMatZJ = 0;
  UniaxialMaterial *hingeMatYJ = 0;
  NDMaterial *hingeMatCJ = 0;
  int dualHingeMode = 0;
  double hingeXiI = 0.0, hingeXiJ = 1.0;
  double elePmP0 = 0.0, elePmAlpha = 0.5;
  bool eleImplex = false;
  double gavy = 0.0, gavz = 0.0;

  numData = 1;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-Integration") == 0 || strcmp(type, "-integration") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: insufficient args for -Integration\n";
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
        opserr << "WARNING LadrunoForceBeamColumn3d: unknown -Integration type '" << integType << "'\n";
        return 0;
      }
    }
    else if (strcmp(type, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING insufficient args for -iter in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &maxIter) < 0 || maxIter <= 0) {
        opserr << "WARNING invalid maxIter in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      if (OPS_GetDoubleInput(&numData, &tol) < 0 || tol <= 0.0) {
        opserr << "WARNING invalid tol in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-subdivide") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING insufficient args for -subdivide in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &numSub) < 0 || numSub <= 0) {
        opserr << "WARNING invalid numSubdivide in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      if (OPS_GetDoubleInput(&numData, &subFac) < 0 || subFac <= 1.0) {
        opserr << "WARNING invalid subdivideFactor in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING insufficient args for -mass in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &mass) < 0 || mass < 0.0 || !std::isfinite(mass)) {
        opserr << "WARNING invalid mass in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-cMass") == 0 || strcmp(type, "-cmass") == 0) {
      cMass = 1;
    }
    else if (strcmp(type, "-lMass") == 0 || strcmp(type, "-lmass") == 0) {
      cMass = 0;
    }
    else if (strcmp(type, "-lch") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING insufficient args for -lch in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      char lchBuf[64];
      const char *lchOpt = OPS_GetStringFromAll(lchBuf, sizeof(lchBuf));
      if (lchOpt == 0) {
        opserr << "WARNING failed to read -lch argument in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      if (strcmp(lchOpt, "ip") == 0 || strcmp(lchOpt, "IP") == 0 || strcmp(lchOpt, "wtL") == 0) {
        lchMode = 0;
        userLch = 0.0;
      } else if (strcmp(lchOpt, "element") == 0 || strcmp(lchOpt, "L") == 0) {
        lchMode = 1;
        userLch = 0.0;
      } else {
        char *endptr = 0;
        double val = strtod(lchOpt, &endptr);
        if (endptr == lchOpt || *endptr != '\0' || !std::isfinite(val) || !(val > 0.0)) {
          opserr << "WARNING invalid -lch value '" << lchOpt << "' in LadrunoForceBeamColumn3d\n";
          delete legacyBi;
          return 0;
        }
        lchMode = 2;
        userLch = val;
      }
    }
    else if (strcmp(type, "-totalElementLch") == 0) {
      lchMode = 1;
      userLch = 0.0;
    }
    else if (strcmp(type, "-damp") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING insufficient args for -damp in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &dampingTag) < 0) {
        delete legacyBi;
        return 0;
      }
      theDamping = OPS_getDamping(dampingTag);
      if (theDamping == 0) {
        opserr << "WARNING damping not found in LadrunoForceBeamColumn3d\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-hinge") == 0 || strcmp(type, "-hingeZ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hinge needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hinge material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatZ = OPS_getUniaxialMaterial(hTagI);
      if (hingeMatZ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hinge uniaxialMaterial " << hTagI << " not found\n";
        delete legacyBi;
        return 0;
      }
      int hTagJ = 0;
      if (lfbc3dTryOptionalInt(hTagJ)) {
        hingeMatZJ = OPS_getUniaxialMaterial(hTagJ);
        if (hingeMatZJ == 0) {
          opserr << "WARNING LadrunoForceBeamColumn3d: -hinge second uniaxialMaterial " << hTagJ << " not found\n";
          delete legacyBi;
          return 0;
        }
        dualHingeMode = 1;
      }
    }
    else if (strcmp(type, "-hingeI") == 0 || strcmp(type, "-hinge1") == 0 || strcmp(type, "-hingeZI") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeI needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeI material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatZ = OPS_getUniaxialMaterial(hTag);
      if (hingeMatZ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeI uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeJ") == 0 || strcmp(type, "-hinge2") == 0 || strcmp(type, "-hingeZJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeJ needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatZJ = OPS_getUniaxialMaterial(hTag);
      if (hingeMatZJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeJ uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeIJ") == 0 || strcmp(type, "-hingeZIJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeIJ needs at least one uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0, hTagJ = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeIJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hTagJ = hTagI;
      int optTagJ = 0;
      if (lfbc3dTryOptionalInt(optTagJ))
        hTagJ = optTagJ;
      hingeMatZ = OPS_getUniaxialMaterial(hTagI);
      hingeMatZJ = OPS_getUniaxialMaterial(hTagJ);
      if (hingeMatZ == 0 || hingeMatZJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeIJ uniaxialMaterial not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeY") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeY needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeY material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatY = OPS_getUniaxialMaterial(hTagI);
      if (hingeMatY == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeY uniaxialMaterial " << hTagI << " not found\n";
        delete legacyBi;
        return 0;
      }
      int hTagJ = 0;
      if (lfbc3dTryOptionalInt(hTagJ)) {
        hingeMatYJ = OPS_getUniaxialMaterial(hTagJ);
        if (hingeMatYJ == 0) {
          opserr << "WARNING LadrunoForceBeamColumn3d: -hingeY second uniaxialMaterial " << hTagJ << " not found\n";
          delete legacyBi;
          return 0;
        }
        dualHingeMode = 1;
      }
    }
    else if (strcmp(type, "-hingeYI") == 0 || strcmp(type, "-hingeY1") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeYI needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeYI material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatY = OPS_getUniaxialMaterial(hTag);
      if (hingeMatY == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeYI uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeYJ") == 0 || strcmp(type, "-hingeY2") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeYJ needs a uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeYJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatYJ = OPS_getUniaxialMaterial(hTag);
      if (hingeMatYJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeYJ uniaxialMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeYIJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeYIJ needs at least one uniaxialMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0, hTagJ = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeYIJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hTagJ = hTagI;
      int optTagJ = 0;
      if (lfbc3dTryOptionalInt(optTagJ))
        hTagJ = optTagJ;
      hingeMatY = OPS_getUniaxialMaterial(hTagI);
      hingeMatYJ = OPS_getUniaxialMaterial(hTagJ);
      if (hingeMatY == 0 || hingeMatYJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeYIJ uniaxialMaterial not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeBiaxial") == 0 || strcmp(type, "-hingeC") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxial needs an NDMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeBiaxial material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatC = OPS_getNDMaterial(hTagI);
      if (hingeMatC == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxial NDMaterial " << hTagI << " not found\n";
        delete legacyBi;
        return 0;
      }
      int hTagJ = 0;
      if (lfbc3dTryOptionalInt(hTagJ)) {
        hingeMatCJ = OPS_getNDMaterial(hTagJ);
        if (hingeMatCJ == 0) {
          opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxial second NDMaterial " << hTagJ << " not found\n";
          delete legacyBi;
          return 0;
        }
        dualHingeMode = 1;
      }
    }
    else if (strcmp(type, "-hingeBiaxialI") == 0 || strcmp(type, "-hingeBiaxial1") == 0 || strcmp(type, "-hingeCI") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxialI needs an NDMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeBiaxialI material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatC = OPS_getNDMaterial(hTag);
      if (hingeMatC == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxialI NDMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeBiaxialJ") == 0 || strcmp(type, "-hingeBiaxial2") == 0 || strcmp(type, "-hingeCJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxialJ needs an NDMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTag = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTag) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeBiaxialJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hingeMatCJ = OPS_getNDMaterial(hTag);
      if (hingeMatCJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxialJ NDMaterial " << hTag << " not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeBiaxialIJ") == 0 || strcmp(type, "-hingeCIJ") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxialIJ needs at least one NDMaterial tag\n";
        delete legacyBi;
        return 0;
      }
      int hTagI = 0, hTagJ = 0;
      numData = 1;
      if (OPS_GetIntInput(&numData, &hTagI) < 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: invalid -hingeBiaxialIJ material tag\n";
        delete legacyBi;
        return 0;
      }
      hTagJ = hTagI;
      int optTagJ = 0;
      if (lfbc3dTryOptionalInt(optTagJ))
        hTagJ = optTagJ;
      hingeMatC = OPS_getNDMaterial(hTagI);
      hingeMatCJ = OPS_getNDMaterial(hTagJ);
      if (hingeMatC == 0 || hingeMatCJ == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxialIJ NDMaterial not found\n";
        delete legacyBi;
        return 0;
      }
      dualHingeMode = 1;
    }
    else if (strcmp(type, "-hingeLoc") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeLoc needs xiI [xiJ]\n";
        delete legacyBi;
        return 0;
      }
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &hingeXiI) < 0 || !std::isfinite(hingeXiI) ||
          hingeXiI < 0.0 || hingeXiI > 1.0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -hingeLoc xiI must be in [0, 1]\n";
        delete legacyBi;
        return 0;
      }
      double optXiJ = 0.0;
      if (lfbc3dTryOptionalDouble(optXiJ)) {
        if (!std::isfinite(optXiJ) || optXiJ < 0.0 || optXiJ > 1.0) {
          opserr << "WARNING LadrunoForceBeamColumn3d: -hingeLoc xiJ must be in [0, 1]\n";
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
        opserr << "WARNING LadrunoForceBeamColumn3d: -pm needs P0 > 0\n";
        delete legacyBi;
        return 0;
      }
      double optAlpha = elePmAlpha;
      if (lfbc3dTryOptionalDouble(optAlpha)) {
        if (!std::isfinite(optAlpha)) {
          opserr << "WARNING LadrunoForceBeamColumn3d: invalid -pm alpha\n";
          delete legacyBi;
          return 0;
        }
        elePmAlpha = optAlpha;
      }
    }
    else if (strcmp(type, "-shear") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          OPS_GetDoubleInput(&numData, &gavy) < 0 || !std::isfinite(gavy) || gavy <= 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -shear needs GAvy > 0 [GAvz > 0]\n";
        delete legacyBi;
        return 0;
      }
      double optGavz = 0.0;
      if (lfbc3dTryOptionalDouble(optGavz)) {
        if (!std::isfinite(optGavz) || optGavz <= 0.0) {
          opserr << "WARNING LadrunoForceBeamColumn3d: -shear GAvz must be > 0\n";
          delete legacyBi;
          return 0;
        }
        gavz = optGavz;
      } else {
        gavz = gavy;
      }
    }
    else if (strcmp(type, "-shearY") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          OPS_GetDoubleInput(&numData, &gavy) < 0 || !std::isfinite(gavy) || gavy <= 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -shearY needs GAvy > 0\n";
        delete legacyBi;
        return 0;
      }
    }
    else if (strcmp(type, "-shearZ") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          OPS_GetDoubleInput(&numData, &gavz) < 0 || !std::isfinite(gavz) || gavz <= 0.0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: -shearZ needs GAvz > 0\n";
        delete legacyBi;
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoForceBeamColumn3d: unknown option '" << type << "'\n";
      delete legacyBi;
      return 0;
    }
  }

  if ((hingeMatC != 0 && (hingeMatZ != 0 || hingeMatY != 0)) ||
      (hingeMatCJ != 0 && (hingeMatZJ != 0 || hingeMatYJ != 0))) {
    opserr << "WARNING LadrunoForceBeamColumn3d: -hingeBiaxial is mutually exclusive with -hinge/-hingeY at the same end\n";
    delete legacyBi;
    return 0;
  }

  if (hingeMatZJ != 0 || hingeMatYJ != 0 || hingeMatCJ != 0)
    dualHingeMode = 1;

  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING LadrunoForceBeamColumn3d: coord transformation " << transfTag << " not found\n";
    delete legacyBi;
    return 0;
  }

  BeamIntegration *bi = 0;
  SectionForceDeformation **sections = 0;
  int nSecs = 0;

  if (!legacyForm) {
    BeamIntegrationRule *theRule = OPS_getBeamIntegrationRule(beamIntTag);
    if (theRule == 0) {
      opserr << "WARNING LadrunoForceBeamColumn3d: beam integration " << beamIntTag << " not found\n";
      return 0;
    }
    bi = theRule->getBeamIntegration();
    if (bi == 0) {
      opserr << "WARNING LadrunoForceBeamColumn3d: beam integration is null\n";
      return 0;
    }
    const ID &secTags = theRule->getSectionTags();
    nSecs = secTags.Size();
    if (nSecs <= 0 || nSecs > 20) {
      opserr << "WARNING LadrunoForceBeamColumn3d: number of sections (" << nSecs << ") must be in 1..20\n";
      return 0;
    }
    sections = new SectionForceDeformation *[nSecs];
    for (int i = 0; i < nSecs; i++) {
      sections[i] = OPS_getSectionForceDeformation(secTags(i));
      if (sections[i] == 0) {
        opserr << "WARNING LadrunoForceBeamColumn3d: section " << secTags(i) << " not found\n";
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
        opserr << "WARNING LadrunoForceBeamColumn3d: section " << legacySecTags[i] << " not found\n";
        delete [] sections;
        delete legacyBi;
        return 0;
      }
    }
  }

  UniaxialMaterial *hZ_tmp = 0, *hY_tmp = 0, *hZJ_tmp = 0, *hYJ_tmp = 0;
  NDMaterial *hC_tmp = 0, *hCJ_tmp = 0;
  if (eleImplex || elePmP0 > 0.0) {
    auto cfgUni = [&](UniaxialMaterial *&hPtr, UniaxialMaterial *&tmpPtr) {
      if (hPtr != 0) {
        tmpPtr = hPtr->getCopy();
        if (LadrunoCohesiveHinge *ch = dynamic_cast<LadrunoCohesiveHinge*>(tmpPtr)) {
          if (eleImplex) ch->setImplex(true);
          if (elePmP0 > 0.0) ch->setPM(elePmP0, elePmAlpha);
        } else if (Bilin *b = dynamic_cast<Bilin*>(tmpPtr)) {
          if (eleImplex) b->setImplex(true);
        } else if (ModIMKPeakOriented *m = dynamic_cast<ModIMKPeakOriented*>(tmpPtr)) {
          if (eleImplex) m->setImplex(true);
        }
        hPtr = tmpPtr;
      }
    };
    auto cfgND = [&](NDMaterial *&hPtr, NDMaterial *&tmpPtr) {
      if (hPtr != 0) {
        tmpPtr = hPtr->getCopy("BeamFiber");
        if (tmpPtr == 0) tmpPtr = hPtr->getCopy();
        if (LadrunoCohesiveHingeBiaxial *ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(tmpPtr)) {
          if (eleImplex) ch->setImplex(true);
          if (elePmP0 > 0.0) ch->setPM(elePmP0, elePmAlpha);
        }
        hPtr = tmpPtr;
      }
    };
    cfgUni(hingeMatZ, hZ_tmp);
    cfgUni(hingeMatY, hY_tmp);
    cfgUni(hingeMatZJ, hZJ_tmp);
    cfgUni(hingeMatYJ, hYJ_tmp);
    cfgND(hingeMatC, hC_tmp);
    cfgND(hingeMatCJ, hCJ_tmp);
  }

  Element *theEle = new LadrunoForceBeamColumn3d(
      iData[0], iData[1], iData[2], nSecs, sections,
      *bi, *theTransf, mass, maxIter, tol, numSub, subFac,
      theDamping, cMass, lchMode, userLch,
      hingeMatZ, hingeMatY, hingeMatC,
      hingeMatZJ, hingeMatYJ, hingeMatCJ,
      dualHingeMode, hingeXiI, hingeXiJ,
      elePmP0, elePmAlpha, eleImplex,
      gavy, gavz);

  delete hZ_tmp;
  delete hY_tmp;
  delete hZJ_tmp;
  delete hYJ_tmp;
  delete hC_tmp;
  delete hCJ_tmp;
  delete [] sections;
  delete legacyBi;
  return theEle;
}

LadrunoForceBeamColumn3d::LadrunoForceBeamColumn3d()
  : Element(0, ELE_TAG_LadrunoForceBeamColumn3d),
    connectedExternalNodes(2),
    beamIntegr(0), numSections(0), sections(0), crdTransf(0),
    rho(0.0), cMass(0), maxIters(10), tol(1.0e-12),
    initialFlag(0),
    kv(NEBD, NEBD), Se(NEBD),
    kvcommit(NEBD, NEBD), Secommit(NEBD),
    kvbulk(NEBD, NEBD), kvbulkCommit(NEBD, NEBD),
    fs(0), vs(0), Ssr(0), vscommit(0), fscommit(0), Ssrcommit(0),
    numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0), load(12), eleLoadDirty(false),
    Ki(0), isTorsion(false), theDamping(0),
    current_section_lch(0.0), lchMode(0), userLch(0.0), GAvy(0.0), GAvz(0.0),
    hingeOn(0), dualHingeMode(0), hingeXiI(0.0), hingeXiJ(1.0),
    theHingeZ(0), theHingeY(0), theHingeC(0),
    theHingeZJ(0), theHingeYJ(0), theHingeCJ(0),
    hingeJumpZ(0.0), hingeJumpCommitZ(0.0),
    hingeJumpY(0.0), hingeJumpCommitY(0.0),
    hingeJumpZJ(0.0), hingeJumpCommitZJ(0.0),
    hingeJumpYJ(0.0), hingeJumpCommitYJ(0.0),
    hingeMscaleZ(0.0), hingeMscaleCommitZ(0.0),
    hingeMscaleY(0.0), hingeMscaleCommitY(0.0),
    hingeMscaleZJ(0.0), hingeMscaleCommitZJ(0.0),
    hingeMscaleYJ(0.0), hingeMscaleCommitYJ(0.0),
    pmP0(0.0), pmAlpha(0.5),
    useImplex(false), eleImplexErr(0.0), eleImplexErrCommit(0.0),
    postRestoreSync(false), vRestoreValid(false),
    vCommitBulk(NEBD),
    maxSubdivisions(4), subdivideFactor(10.0),
    parameterID(0)
{
  load.Zero();
  vCommitBulk.Zero();
  theNodes[0] = 0;
  theNodes[1] = 0;
  for (int a = 0; a < 6; a++) {
    vRestore[a] = 0.0;
    for (int m = 0; m < 4; m++) {
      hingeDualKv[a][m] = 0.0;
      hingeDualKvCommit[a][m] = 0.0;
    }
  }
  for (int m = 0; m < 4; m++) {
    for (int n = 0; n < 4; n++) {
      hingeDualKaaInv[m][n] = 0.0;
      hingeDualKaaInvCommit[m][n] = 0.0;
    }
  }
}

LadrunoForceBeamColumn3d::LadrunoForceBeamColumn3d(
    int tag, int nodeI, int nodeJ,
    int numSec, SectionForceDeformation **sec,
    BeamIntegration &bi,
    CrdTransf &coordTransf, double massDensPerUnitLength,
    int maxNumIters, double tolerance,
    int maxNumSub, double subFac,
    Damping *damping,
    int cm, int lchM, double uLch,
    UniaxialMaterial *hingeMatZ,
    UniaxialMaterial *hingeMatY,
    NDMaterial *hingeMatC,
    UniaxialMaterial *hingeMatZJ,
    UniaxialMaterial *hingeMatYJ,
    NDMaterial *hingeMatCJ,
    int dualMode, double xiI, double xiJ,
    double inPmP0, double inPmAlpha,
    bool eleImplex,
    double gavy_, double gavz_)
  : Element(tag, ELE_TAG_LadrunoForceBeamColumn3d),
    connectedExternalNodes(2),
    beamIntegr(0), numSections(0), sections(0), crdTransf(0),
    rho(massDensPerUnitLength), cMass(cm), maxIters(maxNumIters), tol(tolerance),
    initialFlag(0),
    kv(NEBD, NEBD), Se(NEBD),
    kvcommit(NEBD, NEBD), Secommit(NEBD),
    kvbulk(NEBD, NEBD), kvbulkCommit(NEBD, NEBD),
    fs(0), vs(0), Ssr(0), vscommit(0), fscommit(0), Ssrcommit(0),
    numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0), load(12), eleLoadDirty(false),
    Ki(0), isTorsion(false), theDamping(0),
    current_section_lch(uLch > 0.0 ? uLch : 0.0), lchMode(lchM), userLch(uLch), GAvy(gavy_), GAvz(gavz_),
    hingeOn((hingeMatZ != 0 || hingeMatY != 0 || hingeMatC != 0 ||
             hingeMatZJ != 0 || hingeMatYJ != 0 || hingeMatCJ != 0) ? 1 : 0),
    dualHingeMode((dualMode != 0 || hingeMatZJ != 0 || hingeMatYJ != 0 || hingeMatCJ != 0) ? 1 : 0),
    hingeXiI(xiI), hingeXiJ(xiJ),
    theHingeZ(0), theHingeY(0), theHingeC(0),
    theHingeZJ(0), theHingeYJ(0), theHingeCJ(0),
    hingeJumpZ(0.0), hingeJumpCommitZ(0.0),
    hingeJumpY(0.0), hingeJumpCommitY(0.0),
    hingeJumpZJ(0.0), hingeJumpCommitZJ(0.0),
    hingeJumpYJ(0.0), hingeJumpCommitYJ(0.0),
    hingeMscaleZ(0.0), hingeMscaleCommitZ(0.0),
    hingeMscaleY(0.0), hingeMscaleCommitY(0.0),
    hingeMscaleZJ(0.0), hingeMscaleCommitZJ(0.0),
    hingeMscaleYJ(0.0), hingeMscaleCommitYJ(0.0),
    pmP0(inPmP0), pmAlpha(inPmAlpha),
    useImplex(eleImplex), eleImplexErr(0.0), eleImplexErrCommit(0.0),
    postRestoreSync(false), vRestoreValid(false),
    vCommitBulk(NEBD),
    maxSubdivisions(maxNumSub), subdivideFactor(subFac),
    parameterID(0)
{
  if (maxSubdivisions < 1)
    maxSubdivisions = 1;
  if (subdivideFactor < 1.0)
    subdivideFactor = 1.0;

  load.Zero();
  vCommitBulk.Zero();
  theNodes[0] = 0;
  theNodes[1] = 0;
  for (int a = 0; a < 6; a++) {
    vRestore[a] = 0.0;
    for (int m = 0; m < 4; m++) {
      hingeDualKv[a][m] = 0.0;
      hingeDualKvCommit[a][m] = 0.0;
    }
  }
  for (int m = 0; m < 4; m++) {
    for (int n = 0; n < 4; n++) {
      hingeDualKaaInv[m][n] = 0.0;
      hingeDualKaaInvCommit[m][n] = 0.0;
    }
  }

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;

  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: LadrunoForceBeamColumn3d::LadrunoForceBeamColumn3d: could not create copy of beam integration object" << endln;
    exit(-1);
  }

  crdTransf = coordTransf.getCopy3d();
  if (crdTransf == 0) {
    opserr << "Error: LadrunoForceBeamColumn3d::LadrunoForceBeamColumn3d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }

  if (damping) {
    theDamping = (*damping).getCopy();
    if (!theDamping) {
      opserr << "Error: LadrunoForceBeamColumn3d::LadrunoForceBeamColumn3d: could not create copy of damping object\n";
      exit(-1);
    }
  }

  if (hingeMatZ != 0)  theHingeZ  = hingeMatZ->getCopy();
  if (hingeMatY != 0)  theHingeY  = hingeMatY->getCopy();
  if (hingeMatC != 0) {
    theHingeC = hingeMatC->getCopy("BeamFiber");
    if (theHingeC == 0) theHingeC = hingeMatC->getCopy();
  }
  if (hingeMatZJ != 0) theHingeZJ = hingeMatZJ->getCopy();
  if (hingeMatYJ != 0) theHingeYJ = hingeMatYJ->getCopy();
  if (hingeMatCJ != 0) {
    theHingeCJ = hingeMatCJ->getCopy("BeamFiber");
    if (theHingeCJ == 0) theHingeCJ = hingeMatCJ->getCopy();
  }

  this->setSectionPointers(numSec, sec);
}

LadrunoForceBeamColumn3d::~LadrunoForceBeamColumn3d()
{
  if (sections != 0) {
    for (int i = 0; i < numSections; i++)
      if (sections[i] != 0)
        delete sections[i];
    delete [] sections;
  }

  if (sizeEleLoads != 0) {
    if (eleLoads != 0)
      delete [] eleLoads;
    if (eleLoadFactors != 0)
      delete [] eleLoadFactors;
  }

  if (fs != 0)
    delete [] fs;
  if (vs != 0)
    delete [] vs;
  if (Ssr != 0)
    delete [] Ssr;
  if (vscommit != 0)
    delete [] vscommit;
  if (fscommit != 0)
    delete [] fscommit;
  if (Ssrcommit != 0)
    delete [] Ssrcommit;

  if (crdTransf != 0)
    delete crdTransf;
  if (beamIntegr != 0)
    delete beamIntegr;
  if (Ki != 0)
    delete Ki;
  if (theDamping)
    delete theDamping;
  if (theHingeZ != 0)
    delete theHingeZ;
  if (theHingeY != 0)
    delete theHingeY;
  if (theHingeC != 0)
    delete theHingeC;
  if (theHingeZJ != 0)
    delete theHingeZJ;
  if (theHingeYJ != 0)
    delete theHingeYJ;
  if (theHingeCJ != 0)
    delete theHingeCJ;
}

int
LadrunoForceBeamColumn3d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
LadrunoForceBeamColumn3d::getExternalNodes(void)
{
  return connectedExternalNodes;
}

Node **
LadrunoForceBeamColumn3d::getNodePtrs(void)
{
  return theNodes;
}

int
LadrunoForceBeamColumn3d::getNumDOF(void)
{
  return NEGD;
}

void
LadrunoForceBeamColumn3d::setDomain(Domain *theDomain)
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
    opserr << "LadrunoForceBeamColumn3d::setDomain: node does not exist in model\n";
    return;
  }

  this->DomainComponent::setDomain(theDomain);

  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();

  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "LadrunoForceBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof\n";
    return;
  }

  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "LadrunoForceBeamColumn3d::setDomain(): Error initializing coordinate transformation\n";
    return;
  }

  if (theDamping && theDamping->setDomain(theDomain, NEBD)) {
    opserr << "LadrunoForceBeamColumn3d::setDomain(): Error initializing damping\n";
    return;
  }

  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "LadrunoForceBeamColumn3d::setDomain() -- zero length for element with tag: " << this->getTag() << endln;
    return;
  }

  if (initialFlag == 0)
    this->initializeSectionHistoryVariables();

  ops_TheActiveElement = this;
  this->update();
  ops_TheActiveElement = 0;
}

int
LadrunoForceBeamColumn3d::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping) {
    if (theDamping)
      delete theDamping;
    theDamping = (*damping).getCopy();
    if (!theDamping) {
      opserr << "LadrunoForceBeamColumn3d::setDamping -- failed to get copy of damping\n";
      return -1;
    }
    if (theDamping->setDomain(theDomain, NEBD)) {
      opserr << "LadrunoForceBeamColumn3d::setDamping -- Error initializing damping\n";
      return -2;
    }
  }
  return 0;
}

int
LadrunoForceBeamColumn3d::commitState()
{
  int err = 0;
  int i = 0;

  if ((err = this->Element::commitState()) != 0) {
    opserr << "LadrunoForceBeamColumn3d::commitState () - failed in base class\n";
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

  kvcommit = kv;
  kvbulkCommit = kvbulk;
  Secommit = Se;

  const Vector &v = crdTransf->getBasicTrialDisp();
  for (int m = 0; m < 6; m++)
    vRestore[m] = v(m);
  vRestoreValid = true;
  postRestoreSync = false;

  if (hingeOn) {
    const double betaI  = 1.0 - hingeXiI;
    const double gammaI = hingeXiI;
    const double betaJ  = 1.0 - hingeXiJ;
    const double gammaJ = hingeXiJ;
    vCommitBulk(0) = v(0);
    vCommitBulk(1) = v(1) - betaI * hingeJumpZ  - betaJ * hingeJumpZJ;
    vCommitBulk(2) = v(2) - gammaI * hingeJumpZ - gammaJ * hingeJumpZJ;
    vCommitBulk(3) = v(3) - betaI * hingeJumpY  - betaJ * hingeJumpYJ;
    vCommitBulk(4) = v(4) - gammaI * hingeJumpY - gammaJ * hingeJumpYJ;
    vCommitBulk(5) = v(5);

    hingeJumpCommitZ  = hingeJumpZ;
    hingeJumpCommitY  = hingeJumpY;
    hingeJumpCommitZJ = hingeJumpZJ;
    hingeJumpCommitYJ = hingeJumpYJ;
    hingeMscaleCommitZ  = hingeMscaleZ;
    hingeMscaleCommitY  = hingeMscaleY;
    hingeMscaleCommitZJ = hingeMscaleZJ;
    hingeMscaleCommitYJ = hingeMscaleYJ;
    for (int a = 0; a < 6; a++)
      for (int m = 0; m < 4; m++)
        hingeDualKvCommit[a][m] = hingeDualKv[a][m];
    for (int m = 0; m < 4; m++)
      for (int n = 0; n < 4; n++)
        hingeDualKaaInvCommit[m][n] = hingeDualKaaInv[m][n];

    if (theHingeZ  != 0) err += theHingeZ->commitState();
    if (theHingeY  != 0) err += theHingeY->commitState();
    if (theHingeC  != 0) err += theHingeC->commitState();
    if (theHingeZJ != 0) err += theHingeZJ->commitState();
    if (theHingeYJ != 0) err += theHingeYJ->commitState();
    if (theHingeCJ != 0) err += theHingeCJ->commitState();
    if (err != 0)
      return err;
  } else {
    vCommitBulk = v;
  }

  eleImplexErrCommit = eleImplexErr;

  if (theDamping && (err = theDamping->commitState()))
    return err;

  return err;
}

int
LadrunoForceBeamColumn3d::revertToLastCommit()
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

  Se = Secommit;
  kv = kvcommit;
  kvbulk = kvbulkCommit;

  if (hingeOn) {
    hingeJumpZ  = hingeJumpCommitZ;
    hingeJumpY  = hingeJumpCommitY;
    hingeJumpZJ = hingeJumpCommitZJ;
    hingeJumpYJ = hingeJumpCommitYJ;
    hingeMscaleZ  = hingeMscaleCommitZ;
    hingeMscaleY  = hingeMscaleCommitY;
    hingeMscaleZJ = hingeMscaleCommitZJ;
    hingeMscaleYJ = hingeMscaleCommitYJ;
    for (int a = 0; a < 6; a++)
      for (int m = 0; m < 4; m++)
        hingeDualKv[a][m] = hingeDualKvCommit[a][m];
    for (int m = 0; m < 4; m++)
      for (int n = 0; n < 4; n++)
        hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n];
    if (theHingeZ  != 0) err += theHingeZ->revertToLastCommit();
    if (theHingeY  != 0) err += theHingeY->revertToLastCommit();
    if (theHingeC  != 0) err += theHingeC->revertToLastCommit();
    if (theHingeZJ != 0) err += theHingeZJ->revertToLastCommit();
    if (theHingeYJ != 0) err += theHingeYJ->revertToLastCommit();
    if (theHingeCJ != 0) err += theHingeCJ->revertToLastCommit();
    if (err != 0)
      return err;
    postRestoreSync = true;
  }

  const double betaI  = 1.0 - hingeXiI;
  const double gammaI = hingeXiI;
  const double betaJ  = 1.0 - hingeXiJ;
  const double gammaJ = hingeXiJ;
  vCommitBulk(0) = vRestore[0];
  vCommitBulk(1) = vRestore[1] - (hingeOn ? (betaI * hingeJumpCommitZ  + betaJ * hingeJumpCommitZJ) : 0.0);
  vCommitBulk(2) = vRestore[2] - (hingeOn ? (gammaI * hingeJumpCommitZ + gammaJ * hingeJumpCommitZJ) : 0.0);
  vCommitBulk(3) = vRestore[3] - (hingeOn ? (betaI * hingeJumpCommitY  + betaJ * hingeJumpCommitYJ) : 0.0);
  vCommitBulk(4) = vRestore[4] - (hingeOn ? (gammaI * hingeJumpCommitY + gammaJ * hingeJumpCommitYJ) : 0.0);
  vCommitBulk(5) = vRestore[5];

  eleImplexErr = eleImplexErrCommit;
  eleLoadDirty = (numEleLoads > 0);
  initialFlag = 1;

  if (theDamping && (err = theDamping->revertToLastCommit()))
    return err;

  return err;
}

int
LadrunoForceBeamColumn3d::revertToStart()
{
  int err = 0;
  int i = 0;

  if (Ki != 0) {
    delete Ki;
    Ki = 0;
  }

  do {
    fs[i].Zero();
    fscommit[i].Zero();
    vs[i].Zero();
    Ssr[i].Zero();
    Ssrcommit[i].Zero();
    vscommit[i].Zero();
    err = sections[i++]->revertToStart();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  this->initializeSectionHistoryVariables();

  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  Secommit.Zero();
  kvcommit.Zero();
  kvbulkCommit.Zero();
  Se.Zero();
  kv.Zero();
  kvbulk.Zero();
  vCommitBulk.Zero();

  if (hingeOn) {
    hingeJumpZ = hingeJumpCommitZ = 0.0;
    hingeJumpY = hingeJumpCommitY = 0.0;
    hingeJumpZJ = hingeJumpCommitZJ = 0.0;
    hingeJumpYJ = hingeJumpCommitYJ = 0.0;
    hingeMscaleZ = hingeMscaleCommitZ = 0.0;
    hingeMscaleY = hingeMscaleCommitY = 0.0;
    hingeMscaleZJ = hingeMscaleCommitZJ = 0.0;
    hingeMscaleYJ = hingeMscaleCommitYJ = 0.0;
    for (int a = 0; a < 6; a++) {
      for (int m = 0; m < 4; m++) {
        hingeDualKv[a][m] = 0.0;
        hingeDualKvCommit[a][m] = 0.0;
      }
    }
    for (int m = 0; m < 4; m++) {
      for (int n = 0; n < 4; n++) {
        hingeDualKaaInv[m][n] = 0.0;
        hingeDualKaaInvCommit[m][n] = 0.0;
      }
    }
    if (theHingeZ  != 0) err += theHingeZ->revertToStart();
    if (theHingeY  != 0) err += theHingeY->revertToStart();
    if (theHingeC  != 0) err += theHingeC->revertToStart();
    if (theHingeZJ != 0) err += theHingeZJ->revertToStart();
    if (theHingeYJ != 0) err += theHingeYJ->revertToStart();
    if (theHingeCJ != 0) err += theHingeCJ->revertToStart();
    if (err != 0)
      return err;
  }

  for (int m = 0; m < 6; m++)
    vRestore[m] = 0.0;
  vRestoreValid = false;
  postRestoreSync = false;
  eleImplexErr = eleImplexErrCommit = 0.0;
  eleLoadDirty = false;
  initialFlag = 0;

  if (theDamping && (err = theDamping->revertToStart()))
    return err;

  if (this->getDomain() != 0 && theNodes[0] != 0 && theNodes[1] != 0) {
    ops_TheActiveElement = this;
    err += this->update();
    ops_TheActiveElement = 0;
  }

  return err;
}

const Matrix &
LadrunoForceBeamColumn3d::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  static thread_local Matrix f(NEBD, NEBD);
  this->getInitialFlexibility(f, false);

  static thread_local Matrix I(NEBD, NEBD);
  I.Zero();
  for (int i = 0; i < NEBD; i++)
    I(i, i) = 1.0;

  static thread_local Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "LadrunoForceBeamColumn3d::getInitialStiff() -- could not invert flexibility for element with tag: " << this->getTag() << endln;

  if (hingeOn) {
    const double betaI  = 1.0 - hingeXiI;
    const double gammaI = hingeXiI;
    const double betaJ  = 1.0 - hingeXiJ;
    const double gammaJ = hingeXiJ;
    double W[6][4] = {{0.0}};
    W[1][0] = betaI;  W[2][0] = gammaI;
    W[3][1] = betaI;  W[4][1] = gammaI;
    W[1][2] = betaJ;  W[2][2] = gammaJ;
    W[3][3] = betaJ;  W[4][3] = gammaJ;

    bool active[4];
    active[0] = (theHingeZ  != 0 || theHingeC  != 0);
    active[1] = (theHingeY  != 0 || theHingeC  != 0);
    active[2] = (theHingeZJ != 0 || theHingeCJ != 0);
    active[3] = (theHingeYJ != 0 || theHingeCJ != 0);

    double Kva0[6][4] = {{0.0}};
    for (int m = 0; m < 4; m++) {
      if (active[m]) {
        for (int a = 0; a < 6; a++) {
          double kv_am = 0.0;
          for (int b = 0; b < 6; b++)
            kv_am -= kvInit(a, b) * W[b][m];
          Kva0[a][m] = kv_am;
        }
      }
    }

    double Kaa0[4][4] = {{0.0}};
    for (int m = 0; m < 4; m++) {
      for (int n = 0; n < 4; n++) {
        if (active[m] && active[n]) {
          double sum = 0.0;
          for (int a = 0; a < 6; a++)
            sum -= W[a][m] * Kva0[a][n];
          Kaa0[m][n] = sum;
        }
      }
    }

    double CdiagSum = 0.0;
    if (theHingeC != 0) {
      const Matrix &Kc0 = theHingeC->getInitialTangent();
      double c00 = Kc0(0, 0), c11 = Kc0(1, 1), c01 = 0.5 * (Kc0(0, 1) + Kc0(1, 0));
      Kaa0[0][0] += c00; Kaa0[1][1] += c11;
      Kaa0[0][1] += c01; Kaa0[1][0] += c01;
      CdiagSum += fabs(c00) + fabs(c11);
    } else {
      if (theHingeZ != 0) {
        double c00 = theHingeZ->getInitialTangent();
        Kaa0[0][0] += c00;
        CdiagSum += fabs(c00);
      }
      if (theHingeY != 0) {
        double c11 = theHingeY->getInitialTangent();
        Kaa0[1][1] += c11;
        CdiagSum += fabs(c11);
      }
    }
    if (theHingeCJ != 0) {
      const Matrix &Kc0J = theHingeCJ->getInitialTangent();
      double c22 = Kc0J(0, 0), c33 = Kc0J(1, 1), c23 = 0.5 * (Kc0J(0, 1) + Kc0J(1, 0));
      Kaa0[2][2] += c22; Kaa0[3][3] += c33;
      Kaa0[2][3] += c23; Kaa0[3][2] += c23;
      CdiagSum += fabs(c22) + fabs(c33);
    } else {
      if (theHingeZJ != 0) {
        double c22 = theHingeZJ->getInitialTangent();
        Kaa0[2][2] += c22;
        CdiagSum += fabs(c22);
      }
      if (theHingeYJ != 0) {
        double c33 = theHingeYJ->getInitialTangent();
        Kaa0[3][3] += c33;
        CdiagSum += fabs(c33);
      }
    }

    double bulkDiagSum = 0.0;
    for (int a = 1; a <= 4; a++) bulkDiagSum += fabs(kvInit(a, a));
    double floor = 1.0e-8 * (bulkDiagSum + CdiagSum) + 1.0e-300;
    double invKaa0[4][4] = {{0.0}};
    lfbc3dFlooredInv4x4(Kaa0, active, floor, invKaa0);

    for (int a = 0; a < 6; a++) {
      for (int b = 0; b < 6; b++) {
        double sub = 0.0;
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            sub += Kva0[a][m] * invKaa0[m][n] * Kva0[b][n];
        kvInit(a, b) -= sub;
      }
    }
  }

  if (theDamping)
    kvInit *= theDamping->getStiffnessMultiplier();

  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));
  return *Ki;
}

const Matrix &
LadrunoForceBeamColumn3d::getTangentStiff(void)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  crdTransf->update();
  return crdTransf->getGlobalStiffMatrix(kv, Se);
}

void
LadrunoForceBeamColumn3d::computeReactions(double *p0)
{
  int type;
  double L = crdTransf->getInitialLength();

  for (int i = 0; i < numEleLoads; i++) {
    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0) * loadFactor;
      double wz = data(1) * loadFactor;
      double wa = data(2) * loadFactor;

      p0[0] -= wa * L;
      double V = 0.5 * wy * L;
      p0[1] -= V;
      p0[2] -= V;
      V = 0.5 * wz * L;
      p0[3] -= V;
      p0[4] -= V;
    }
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wy = data(0) * loadFactor;
      double wz = data(1) * loadFactor;
      double wa = data(2) * loadFactor;
      double a = data(3) * L;
      double b = data(4) * L;
      if (b <= a || a < 0.0 || b > L)
        continue;
      double wyb = data(5) * loadFactor;
      double wzb = data(6) * loadFactor;
      double wab = data(7) * loadFactor;
      p0[0] -= wa * (b - a) + 0.5 * (wab - wa) * (b - a);
      double c = a + 0.5 * (b - a);
      double Fy = wy * (b - a);
      p0[1] -= Fy * (1 - c / L);
      p0[2] -= Fy * c / L;
      double Fz = wz * (b - a);
      p0[3] -= Fz * (1 - c / L);
      p0[4] -= Fz * c / L;
      c = a + 2.0 / 3.0 * (b - a);
      Fy = 0.5 * (wyb - wy) * (b - a);
      p0[1] -= Fy * (1 - c / L);
      p0[2] -= Fy * c / L;
      Fz = 0.5 * (wzb - wz) * (b - a);
      p0[3] -= Fz * (1 - c / L);
      p0[4] -= Fz * c / L;
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0) * loadFactor;
      double Pz = data(1) * loadFactor;
      double N = data(2) * loadFactor;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double V1 = Py * (1.0 - aOverL);
      double V2 = Py * aOverL;
      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
      V1 = Pz * (1.0 - aOverL);
      V2 = Pz * aOverL;
      p0[3] -= V1;
      p0[4] -= V2;
    }
  }
}

void
LadrunoForceBeamColumn3d::computeReactionSensitivity(double *dp0dh, int gradNumber)
{
  int type;
  double L = crdTransf->getInitialLength();
  double dLdh = crdTransf->getdLdh();

  for (int i = 0; i < numEleLoads; i++) {
    const Vector &data = eleLoads[i]->getData(type, 1.0);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0);
      double wz = data(1);
      double wa = data(2);

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwzdh = sens(1);
      double dwadh = sens(2);

      dp0dh[0] -= wa * dLdh + dwadh * L;
      double dVdh = 0.5 * (wy * dLdh + dwydh * L);
      dp0dh[1] -= dVdh;
      dp0dh[2] -= dVdh;
      dVdh = 0.5 * (wz * dLdh + dwzdh * L);
      dp0dh[3] -= dVdh;
      dp0dh[4] -= dVdh;
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0);
      double Pz = data(1);
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPydh = sens(0);
      double dPzdh = sens(1);
      double dNdh = sens(2);
      double daLdh = sens(3);

      double dV1dh = Py * (0.0 - daLdh) + dPydh * (1.0 - aOverL);
      double dV2dh = Py * daLdh + dPydh * aOverL;

      dp0dh[0] -= dNdh;
      dp0dh[1] -= dV1dh;
      dp0dh[2] -= dV2dh;

      dV1dh = Pz * (0.0 - daLdh) + dPzdh * (1.0 - aOverL);
      dV2dh = Pz * daLdh + dPzdh * aOverL;
      dp0dh[3] -= dV1dh;
      dp0dh[4] -= dV2dh;
    }
  }
}

const Vector &
LadrunoForceBeamColumn3d::getResistingForce(void)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  crdTransf->update();

  double p0[5];
  Vector p0Vec(p0, 5);
  p0Vec.Zero();

  if (numEleLoads > 0)
    this->computeReactions(p0);

  theVector = crdTransf->getGlobalResistingForce(Se, p0Vec);

  if (rho != 0)
    theVector.addVector(1.0, load, -1.0);

  return theVector;
}

const Vector &
LadrunoForceBeamColumn3d::getDampingForce(void)
{
  crdTransf->update();
  return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(5));
}

void
LadrunoForceBeamColumn3d::initializeSectionHistoryVariables(void)
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
LadrunoForceBeamColumn3d::evalBulkState(const Vector &v_tilde, const Vector &v_tilde_commit,
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

  SeCur = (initialFlag == 0) ? Se : Secommit;
  kvCur = (initialFlag == 0) ? kvbulk : kvbulkCommit;
  for (int i = 0; i < numSections; i++) {
    vsCur[i] = (initialFlag == 0) ? vs[i] : vscommit[i];
    fsCur[i] = (initialFlag == 0) ? fs[i] : fscommit[i];
    SsrCur[i] = (initialFlag == 0) ? Ssr[i] : Ssrcommit[i];
  }

  vin = v_tilde_commit;
  dvToDo = v_tilde;
  dvToDo -= v_tilde_commit;
  dvTrial = dvToDo;

  double factor = subdivideFactor;
  double dW = 0.0, dW0 = 0.0;
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
                Ss(ii) = SeTrial(0);
                break;
              case SECTION_RESPONSE_MZ:
                Ss(ii) = xL1 * SeTrial(1) + xL * SeTrial(2);
                break;
              case SECTION_RESPONSE_VY:
                Ss(ii) = oneOverL * (SeTrial(1) + SeTrial(2));
                break;
              case SECTION_RESPONSE_MY:
                Ss(ii) = xL1 * SeTrial(3) + xL * SeTrial(4);
                break;
              case SECTION_RESPONSE_VZ:
                Ss(ii) = oneOverL * (SeTrial(3) + SeTrial(4));
                break;
              case SECTION_RESPONSE_T:
                Ss(ii) = SeTrial(5);
                break;
              default:
                Ss(ii) = 0.0;
                break;
              }
            }

            if (numEleLoads > 0)
              this->computeSectionForces(Ss, i);

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
              { LadrunoOutputGuard ladrunoOut; opserr << "LadrunoForceBeamColumn3d::evalBulkState() - section failed in setTrial\n"; }
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
                for (int jj = 0; jj < order; jj++)
                  fb(jj, 0) += fSec(jj, ii) * wtL;
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
              case SECTION_RESPONSE_MY:
                for (int jj = 0; jj < order; jj++) {
                  tmp = fSec(jj, ii) * wtL;
                  fb(jj, 3) += xL1 * tmp;
                  fb(jj, 4) += xL * tmp;
                }
                break;
              case SECTION_RESPONSE_VZ:
                for (int jj = 0; jj < order; jj++) {
                  tmp = oneOverL * fSec(jj, ii) * wtL;
                  fb(jj, 3) += tmp;
                  fb(jj, 4) += tmp;
                }
                break;
              case SECTION_RESPONSE_T:
                for (int jj = 0; jj < order; jj++)
                  fb(jj, 5) += fSec(jj, ii) * wtL;
                break;
              default:
                break;
              }
            }

            for (int ii = 0; ii < order; ii++) {
              switch (code(ii)) {
              case SECTION_RESPONSE_P:
                for (int jj = 0; jj < NEBD; jj++)
                  f(0, jj) += fb(ii, jj);
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
              case SECTION_RESPONSE_MY:
                for (int jj = 0; jj < NEBD; jj++) {
                  tmp = fb(ii, jj);
                  f(3, jj) += xL1 * tmp;
                  f(4, jj) += xL * tmp;
                }
                break;
              case SECTION_RESPONSE_VZ:
                for (int jj = 0; jj < NEBD; jj++) {
                  tmp = oneOverL * fb(ii, jj);
                  f(3, jj) += tmp;
                  f(4, jj) += tmp;
                }
                break;
              case SECTION_RESPONSE_T:
                for (int jj = 0; jj < NEBD; jj++)
                  f(5, jj) += fb(ii, jj);
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
                vr(0) += dei;
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
              case SECTION_RESPONSE_MY:
                vr(3) += xL1 * dei;
                vr(4) += xL * dei;
                break;
              case SECTION_RESPONSE_VZ:
                tmp = oneOverL * dei;
                vr(3) += tmp;
                vr(4) += tmp;
                break;
              case SECTION_RESPONSE_T:
                vr(5) += dei;
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
            double dvShearY = fsY * (SeTrial(1) + SeTrial(2));
            vr(1) += dvShearY;
            vr(2) += dvShearY;
          }
          if (GAvz > 0.0) {
            double fsZ = oneOverL / GAvz;
            f(3, 3) += fsZ;
            f(3, 4) += fsZ;
            f(4, 3) += fsZ;
            f(4, 4) += fsZ;
            double dvShearZ = fsZ * (SeTrial(3) + SeTrial(4));
            vr(3) += dvShearZ;
            vr(4) += dvShearZ;
          }

          if (!isTorsion) {
            f(5, 5) = DefaultLoverGJ;
            vr(5) = SeTrial(5) * DefaultLoverGJ;
          }

          if (f.Solve(I, kvTrial) < 0) {
            { LadrunoOutputGuard ladrunoOut; opserr << "LadrunoForceBeamColumn3d::evalBulkState() -- could not invert flexibility for element with tag: " << this->getTag() << endln; }
            return -1;
          }

          dv = vin;
          dv += dvTrial;
          dv -= vr;

          dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);

          dW = dv ^ dSe;
          if (dW0 == 0.0)
            dW0 = dW;

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
  }

  return 0;
}

int
LadrunoForceBeamColumn3d::solveHingeState(const Vector &v)
{
  const double betaI  = 1.0 - hingeXiI;
  const double gammaI = hingeXiI;
  const double betaJ  = 1.0 - hingeXiJ;
  const double gammaJ = hingeXiJ;

  double W[6][4] = {{0.0}};
  W[1][0] = betaI;  W[2][0] = gammaI;
  W[3][1] = betaI;  W[4][1] = gammaI;
  W[1][2] = betaJ;  W[2][2] = gammaJ;
  W[3][3] = betaJ;  W[4][3] = gammaJ;

  bool active[4];
  active[0] = (theHingeZ  != 0 || theHingeC  != 0);
  active[1] = (theHingeY  != 0 || theHingeC  != 0);
  active[2] = (theHingeZJ != 0 || theHingeCJ != 0);
  active[3] = (theHingeYJ != 0 || theHingeCJ != 0);

  double alpha[4];
  alpha[0] = active[0] ? hingeJumpCommitZ  : 0.0;
  alpha[1] = active[1] ? hingeJumpCommitY  : 0.0;
  alpha[2] = active[2] ? hingeJumpCommitZJ : 0.0;
  alpha[3] = active[3] ? hingeJumpCommitYJ : 0.0;

  LadrunoCohesiveHinge *chZI = dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ);
  LadrunoCohesiveHinge *chYI = dynamic_cast<LadrunoCohesiveHinge*>(theHingeY);
  LadrunoCohesiveHingeBiaxial *chCI = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC);
  LadrunoCohesiveHinge *chZJ = dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ);
  LadrunoCohesiveHinge *chYJ = dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ);
  LadrunoCohesiveHingeBiaxial *chCJ = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeCJ);

  const bool freezeCI = postRestoreSync && lfbc3dHingeBiaxialIsImplex(theHingeC);
  const bool freezeCJ = postRestoreSync && lfbc3dHingeBiaxialIsImplex(theHingeCJ);
  const bool freezeZI = postRestoreSync && lfbc3dHingeIsImplex(theHingeZ);
  const bool freezeYI = postRestoreSync && lfbc3dHingeIsImplex(theHingeY);
  const bool freezeZJ = postRestoreSync && lfbc3dHingeIsImplex(theHingeZJ);
  const bool freezeYJ = postRestoreSync && lfbc3dHingeIsImplex(theHingeYJ);

  static thread_local Vector v_tilde(NEBD);
  static thread_local Vector q_bulk(NEBD);
  static thread_local Matrix Kbulk(NEBD, NEBD);

  double Kva[6][4] = {{0.0}};
  double invKaa[4][4] = {{0.0}};

  const int maxOuter = 20;
  for (int outer = 0; outer < maxOuter; outer++) {
    v_tilde(0) = v(0);
    v_tilde(1) = v(1) - betaI * alpha[0] - betaJ * alpha[2];
    v_tilde(2) = v(2) - gammaI * alpha[0] - gammaJ * alpha[2];
    v_tilde(3) = v(3) - betaI * alpha[1] - betaJ * alpha[3];
    v_tilde(4) = v(4) - gammaI * alpha[1] - gammaJ * alpha[3];
    v_tilde(5) = v(5);

    int err = this->evalBulkState(v_tilde, vCommitBulk, q_bulk, Kbulk, false);
    if (err != 0)
      return err;

    for (int m = 0; m < 4; m++) {
      for (int a = 0; a < 6; a++) {
        double kv_am = 0.0;
        if (active[m]) {
          for (int b = 0; b < 6; b++)
            kv_am -= Kbulk(a, b) * W[b][m];
        }
        Kva[a][m] = kv_am;
      }
    }

    auto evalCohesive = [&](const double aVec[4], double q0_ax,
                            double Mcoh[4], double Kaa[4][4], double &CdiagSum) {
      for (int m = 0; m < 4; m++) {
        Mcoh[m] = 0.0;
        for (int n = 0; n < 4; n++) {
          if (active[m] && active[n]) {
            double sum = 0.0;
            for (int a = 0; a < 6; a++)
              sum -= W[a][m] * Kva[a][n];
            Kaa[m][n] = sum;
          } else {
            Kaa[m][n] = (m == n) ? 1.0 : 0.0;
          }
        }
      }
      CdiagSum = 0.0;

      if (theHingeC != 0) {
        if (chCI != 0) chCI->setAxialForce(q0_ax);
        double fPMCI = (chCI == 0 && pmP0 > 0.0) ? lfbc3dScalePM(q0_ax, pmP0, pmAlpha) : 1.0;
        if (!freezeCI) {
          static thread_local Vector jumpI(2);
          jumpI(0) = aVec[0];
          jumpI(1) = aVec[1];
          theHingeC->setTrialStrain(jumpI);
        }
        const Vector &MI = theHingeC->getStress();
        const Matrix &KcI = theHingeC->getTangent();
        Mcoh[0] = fPMCI * MI(0);
        Mcoh[1] = fPMCI * MI(1);
        double c00 = fPMCI * KcI(0, 0), c11 = fPMCI * KcI(1, 1), c01 = 0.5 * fPMCI * (KcI(0, 1) + KcI(1, 0));
        Kaa[0][0] += c00; Kaa[1][1] += c11;
        Kaa[0][1] += c01; Kaa[1][0] += c01;
        CdiagSum += fabs(c00) + fabs(c11);
      } else {
        if (theHingeZ != 0) {
          if (chZI != 0) chZI->setAxialForce(q0_ax);
          double fPMZI = (chZI == 0 && pmP0 > 0.0) ? lfbc3dScalePM(q0_ax, pmP0, pmAlpha) : 1.0;
          if (!freezeZI)
            theHingeZ->setTrialStrain(aVec[0]);
          Mcoh[0] = fPMZI * theHingeZ->getStress();
          double c00 = fPMZI * theHingeZ->getTangent();
          Kaa[0][0] += c00;
          CdiagSum += fabs(c00);
        }
        if (theHingeY != 0) {
          if (chYI != 0) chYI->setAxialForce(q0_ax);
          double fPMYI = (chYI == 0 && pmP0 > 0.0) ? lfbc3dScalePM(q0_ax, pmP0, pmAlpha) : 1.0;
          if (!freezeYI)
            theHingeY->setTrialStrain(aVec[1]);
          Mcoh[1] = fPMYI * theHingeY->getStress();
          double c11 = fPMYI * theHingeY->getTangent();
          Kaa[1][1] += c11;
          CdiagSum += fabs(c11);
        }
      }

      if (theHingeCJ != 0) {
        if (chCJ != 0) chCJ->setAxialForce(q0_ax);
        double fPMCJ = (chCJ == 0 && pmP0 > 0.0) ? lfbc3dScalePM(q0_ax, pmP0, pmAlpha) : 1.0;
        if (!freezeCJ) {
          static thread_local Vector jumpJ(2);
          jumpJ(0) = aVec[2];
          jumpJ(1) = aVec[3];
          theHingeCJ->setTrialStrain(jumpJ);
        }
        const Vector &MJ = theHingeCJ->getStress();
        const Matrix &KcJ = theHingeCJ->getTangent();
        Mcoh[2] = fPMCJ * MJ(0);
        Mcoh[3] = fPMCJ * MJ(1);
        double c22 = fPMCJ * KcJ(0, 0), c33 = fPMCJ * KcJ(1, 1), c23 = 0.5 * fPMCJ * (KcJ(0, 1) + KcJ(1, 0));
        Kaa[2][2] += c22; Kaa[3][3] += c33;
        Kaa[2][3] += c23; Kaa[3][2] += c23;
        CdiagSum += fabs(c22) + fabs(c33);
      } else {
        if (theHingeZJ != 0) {
          if (chZJ != 0) chZJ->setAxialForce(q0_ax);
          double fPMZJ = (chZJ == 0 && pmP0 > 0.0) ? lfbc3dScalePM(q0_ax, pmP0, pmAlpha) : 1.0;
          if (!freezeZJ)
            theHingeZJ->setTrialStrain(aVec[2]);
          Mcoh[2] = fPMZJ * theHingeZJ->getStress();
          double c22 = fPMZJ * theHingeZJ->getTangent();
          Kaa[2][2] += c22;
          CdiagSum += fabs(c22);
        }
        if (theHingeYJ != 0) {
          if (chYJ != 0) chYJ->setAxialForce(q0_ax);
          double fPMYJ = (chYJ == 0 && pmP0 > 0.0) ? lfbc3dScalePM(q0_ax, pmP0, pmAlpha) : 1.0;
          if (!freezeYJ)
            theHingeYJ->setTrialStrain(aVec[3]);
          Mcoh[3] = fPMYJ * theHingeYJ->getStress();
          double c33 = fPMYJ * theHingeYJ->getTangent();
          Kaa[3][3] += c33;
          CdiagSum += fabs(c33);
        }
      }
    };

    double Mcoh[4], Kaa[4][4], CdiagSum = 0.0;
    evalCohesive(alpha, q_bulk(0), Mcoh, Kaa, CdiagSum);

    double Mbulk[4] = {0.0, 0.0, 0.0, 0.0};
    for (int m = 0; m < 4; m++)
      for (int a = 0; a < 6; a++)
        Mbulk[m] += W[a][m] * q_bulk(a);

    bool freezeDOF[4] = {
      freezeCI || freezeZI,
      freezeCI || freezeYI,
      freezeCJ || freezeZJ,
      freezeCJ || freezeYJ
    };

    double h[4] = {0.0, 0.0, 0.0, 0.0};
    for (int m = 0; m < 4; m++)
      if (active[m] && !freezeDOF[m] && initialFlag != 0)
        h[m] = -Mbulk[m] + Mcoh[m];

    double bulkDiagSum = 0.0;
    for (int a = 1; a <= 4; a++) bulkDiagSum += fabs(Kbulk(a, a));
    double floor = 1.0e-8 * (bulkDiagSum + CdiagSum) + 1.0e-300;
    lfbc3dFlooredInv4x4(Kaa, active, floor, invKaa);

    double *mscales[4] = {&hingeMscaleZ, &hingeMscaleY, &hingeMscaleZJ, &hingeMscaleYJ};
    bool allOk = true;
    bool allZero = true;
    for (int m = 0; m < 4; m++) {
      if (active[m]) {
        double Mhere = fabs(Mcoh[m]) + fabs(Mbulk[m]);
        if (Mhere > *mscales[m]) *mscales[m] = Mhere;
        if (h[m] != 0.0)
          allZero = false;
        if (fabs(h[m]) > 1.0e-12 * (*mscales[m]) + 1.0e-14)
          allOk = false;
      }
    }

    if (allZero || (outer > 0 && allOk))
      break;

    // Inner 4-DOF Newton loop on the linearized bulk response
    double a0[4] = {alpha[0], alpha[1], alpha[2], alpha[3]};
    double damp = 1.0;
    double hnormPrev = 1.0e300;
    for (int inner = 0; inner < 35; inner++) {
      double q_lin[6];
      for (int a = 0; a < 6; a++) {
        q_lin[a] = q_bulk(a);
        for (int m = 0; m < 4; m++)
          q_lin[a] += Kva[a][m] * (alpha[m] - a0[m]);
      }

      evalCohesive(alpha, q_lin[0], Mcoh, Kaa, CdiagSum);

      for (int m = 0; m < 4; m++) {
        Mbulk[m] = 0.0;
        for (int a = 0; a < 6; a++)
          Mbulk[m] += W[a][m] * q_lin[a];
        h[m] = (active[m] && !freezeDOF[m] && initialFlag != 0) ? (-Mbulk[m] + Mcoh[m]) : 0.0;
      }

      floor = 1.0e-8 * (bulkDiagSum + CdiagSum) + 1.0e-300;
      lfbc3dFlooredInv4x4(Kaa, active, floor, invKaa);

      bool innerOk = true;
      bool innerAllZero = true;
      double hnormSq = 0.0;
      for (int m = 0; m < 4; m++) {
        if (active[m]) {
          double Mhere = fabs(Mcoh[m]) + fabs(Mbulk[m]);
          if (Mhere > *mscales[m]) *mscales[m] = Mhere;
          if (h[m] != 0.0)
            innerAllZero = false;
          if (fabs(h[m]) > 1.0e-13 * (*mscales[m]) + 1.0e-15)
            innerOk = false;
          double hn = h[m] / (*mscales[m] + 1.0e-12);
          hnormSq += hn * hn;
        }
      }
      if (innerAllZero || (inner > 0 && innerOk))
        break;

      double hnorm = sqrt(hnormSq);
      if (hnorm > hnormPrev)
        damp *= 0.5;
      else
        damp = (damp * 1.5 < 1.0) ? damp * 1.5 : 1.0;
      hnormPrev = hnorm;

      for (int m = 0; m < 4; m++) {
        if (active[m] && !freezeDOF[m]) {
          double da = 0.0;
          for (int n = 0; n < 4; n++)
            da += invKaa[m][n] * h[n];
          alpha[m] -= damp * da;
        }
      }
    }
  }

  v_tilde(0) = v(0);
  v_tilde(1) = v(1) - betaI * alpha[0] - betaJ * alpha[2];
  v_tilde(2) = v(2) - gammaI * alpha[0] - gammaJ * alpha[2];
  v_tilde(3) = v(3) - betaI * alpha[1] - betaJ * alpha[3];
  v_tilde(4) = v(4) - gammaI * alpha[1] - gammaJ * alpha[3];
  v_tilde(5) = v(5);

  int err = this->evalBulkState(v_tilde, vCommitBulk, Se, kvbulk, true);
  if (err != 0)
    return err;

  if (theHingeC != 0) {
    if (chCI != 0) chCI->setAxialForce(Se(0));
    if (!freezeCI) {
      static thread_local Vector jumpI(2);
      jumpI(0) = alpha[0];
      jumpI(1) = alpha[1];
      theHingeC->setTrialStrain(jumpI);
    }
  } else {
    if (theHingeZ != 0) {
      if (chZI != 0) chZI->setAxialForce(Se(0));
      if (!freezeZI) theHingeZ->setTrialStrain(alpha[0]);
    }
    if (theHingeY != 0) {
      if (chYI != 0) chYI->setAxialForce(Se(0));
      if (!freezeYI) theHingeY->setTrialStrain(alpha[1]);
    }
  }
  if (theHingeCJ != 0) {
    if (chCJ != 0) chCJ->setAxialForce(Se(0));
    if (!freezeCJ) {
      static thread_local Vector jumpJ(2);
      jumpJ(0) = alpha[2];
      jumpJ(1) = alpha[3];
      theHingeCJ->setTrialStrain(jumpJ);
    }
  } else {
    if (theHingeZJ != 0) {
      if (chZJ != 0) chZJ->setAxialForce(Se(0));
      if (!freezeZJ) theHingeZJ->setTrialStrain(alpha[2]);
    }
    if (theHingeYJ != 0) {
      if (chYJ != 0) chYJ->setAxialForce(Se(0));
      if (!freezeYJ) theHingeYJ->setTrialStrain(alpha[3]);
    }
  }

  hingeJumpZ  = alpha[0];
  hingeJumpY  = alpha[1];
  hingeJumpZJ = alpha[2];
  hingeJumpYJ = alpha[3];
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      hingeDualKv[a][m] = Kva[a][m];
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      hingeDualKaaInv[m][n] = invKaa[m][n];

  kv = kvbulk;
  for (int a = 0; a < 6; a++) {
    for (int b = 0; b < 6; b++) {
      double sub = 0.0;
      for (int m = 0; m < 4; m++)
        for (int n = 0; n < 4; n++)
          sub += Kva[a][m] * invKaa[m][n] * Kva[b][n];
      kv(a, b) -= sub;
    }
  }

  if (postRestoreSync) {
    for (int a = 0; a < 6; a++)
      for (int m = 0; m < 4; m++)
        hingeDualKvCommit[a][m] = hingeDualKv[a][m];
    for (int m = 0; m < 4; m++)
      for (int n = 0; n < 4; n++)
        hingeDualKaaInvCommit[m][n] = hingeDualKaaInv[m][n];
  }

  return 0;
}

int
LadrunoForceBeamColumn3d::update()
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

  if (postRestoreSync && vRestoreValid) {
    bool moved = false;
    for (int m = 0; m < 6; m++)
      if (fabs(v(m) - vRestore[m]) > 1.0e-14)
        moved = true;
    if (moved)
      postRestoreSync = false;
  }

  static thread_local Vector dvPred(NEBD);
  static thread_local Vector SePred(NEBD);
  for (int m = 0; m < NEBD; m++)
    dvPred(m) = v(m) - vRestore[m];
  SePred = Secommit;
  SePred.addMatrixVector(1.0, kvcommit, dvPred, 1.0);

  int err = 0;
  if (hingeOn) {
    err = this->solveHingeState(v);
  } else {
    err = this->evalBulkState(v, vCommitBulk, Se, kvbulk, true);
    kv = kvbulk;
  }

  if (err != 0) {
    LadrunoOutputGuard ladrunoOut;
    opserr << "WARNING - LadrunoForceBeamColumn3d::update - failed to get compatible "
           << "element forces & deformations for element: " << this->getTag() << endln;
    ops_TheActiveElement = 0;
    return -1;
  }

  if (initialFlag == 0 || (initialFlag == 3 && !vRestoreValid)) {
    kvcommit = kv;
    kvbulkCommit = kvbulk;
    Secommit = Se;
    for (int k = 0; k < numSections; k++) {
      vscommit[k] = vs[k];
      fscommit[k] = fs[k];
      Ssrcommit[k] = Ssr[k];
    }
    const double betaI  = 1.0 - hingeXiI;
    const double gammaI = hingeXiI;
    const double betaJ  = 1.0 - hingeXiJ;
    const double gammaJ = hingeXiJ;
    vCommitBulk(0) = v(0);
    vCommitBulk(1) = v(1) - (hingeOn ? (betaI * hingeJumpZ  + betaJ * hingeJumpZJ) : 0.0);
    vCommitBulk(2) = v(2) - (hingeOn ? (gammaI * hingeJumpZ + gammaJ * hingeJumpZJ) : 0.0);
    vCommitBulk(3) = v(3) - (hingeOn ? (betaI * hingeJumpY  + betaJ * hingeJumpYJ) : 0.0);
    vCommitBulk(4) = v(4) - (hingeOn ? (gammaI * hingeJumpY + gammaJ * hingeJumpYJ) : 0.0);
    vCommitBulk(5) = v(5);
  }

  if (theDamping) {
    kv *= theDamping->getStiffnessMultiplier();
    theDamping->update(Se);
  }

  eleImplexErr = 0.0;
  if (useImplex && initialFlag != 0 && initialFlag != 3) {
    static thread_local Vector dSeErr(NEBD);
    dSeErr = Se;
    dSeErr -= SePred;
    eleImplexErr = dSeErr.Norm() / (Se.Norm() + 1.0e-12);
  }
  if (hingeOn) {
    UniaxialMaterial *uniList[4] = {theHingeZ, theHingeY, theHingeZJ, theHingeYJ};
    for (int k = 0; k < 4; k++) {
      double e = 0.0;
      if (lfbc3dGetHingeImplexError(uniList[k], e))
        if (e > eleImplexErr) eleImplexErr = e;
    }
    NDMaterial *ndList[2] = {theHingeC, theHingeCJ};
    for (int k = 0; k < 2; k++) {
      if (LadrunoCohesiveHingeBiaxial *ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(ndList[k])) {
        double e = ch->getImplexError();
        if (e > eleImplexErr) eleImplexErr = e;
      }
    }
  }

  initialFlag = 1;
  eleLoadDirty = false;
  ops_TheActiveElement = 0;
  return 0;
}

void
LadrunoForceBeamColumn3d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
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
    case SECTION_RESPONSE_MY:
      b(i, 3) = xi - 1.0;
      b(i, 4) = xi;
      break;
    case SECTION_RESPONSE_VZ:
      b(i, 3) = b(i, 4) = 1.0 / L;
      break;
    case SECTION_RESPONSE_T:
      b(i, 5) = 1.0;
      break;
    default:
      break;
    }
  }
}

void
LadrunoForceBeamColumn3d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
{
  bp.Zero();
  double L = crdTransf->getInitialLength();
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:
      bp(i, 1) = xi * (xi - 1) * L * L / 2;
      break;
    case SECTION_RESPONSE_P:
      bp(i, 0) = (1 - xi) * L;
      break;
    case SECTION_RESPONSE_VY:
      bp(i, 1) = (xi - 0.5) * L;
      break;
    case SECTION_RESPONSE_MY:
      bp(i, 2) = xi * (1 - xi) * L * L / 2;
      break;
    case SECTION_RESPONSE_VZ:
      bp(i, 2) = (0.5 - xi) * L;
      break;
    default:
      break;
    }
  }
}

const Matrix &
LadrunoForceBeamColumn3d::getMass(void)
{
  theMatrix.Zero();
  if (rho == 0.0)
    return theMatrix;

  double L = crdTransf->getInitialLength();
  if (cMass == 0) {
    theMatrix(0, 0) = theMatrix(1, 1) = theMatrix(2, 2) =
      theMatrix(6, 6) = theMatrix(7, 7) = theMatrix(8, 8) = 0.5 * L * rho;
  } else {
    static thread_local Matrix ml(12, 12);
    ml.Zero();
    double m = rho * L / 420.0;
    ml(0, 0) = ml(6, 6) = m * 140.0;
    ml(0, 6) = ml(6, 0) = m * 70.0;

    ml(2, 2) = ml(8, 8) = m * 156.0;
    ml(2, 8) = ml(8, 2) = m * 54.0;
    ml(4, 4) = ml(10, 10) = m * 4.0 * L * L;
    ml(4, 10) = ml(10, 4) = -m * 3.0 * L * L;
    ml(2, 4) = ml(4, 2) = -m * 22.0 * L;
    ml(8, 10) = ml(10, 8) = -ml(2, 4);
    ml(2, 10) = ml(10, 2) = m * 13.0 * L;
    ml(4, 8) = ml(8, 4) = -ml(2, 10);

    ml(1, 1) = ml(7, 7) = m * 156.0;
    ml(1, 7) = ml(7, 1) = m * 54.0;
    ml(5, 5) = ml(11, 11) = m * 4.0 * L * L;
    ml(5, 11) = ml(11, 5) = -m * 3.0 * L * L;
    ml(1, 5) = ml(5, 1) = m * 22.0 * L;
    ml(7, 11) = ml(11, 7) = -ml(1, 5);
    ml(1, 11) = ml(11, 1) = -m * 13.0 * L;
    ml(5, 7) = ml(7, 5) = -ml(1, 11);

    theMatrix = crdTransf->getGlobalMatrixFromLocal(ml);
  }

  return theMatrix;
}

void
LadrunoForceBeamColumn3d::zeroLoad(void)
{
  load.Zero();
  if (numEleLoads > 0)
    eleLoadDirty = true;
  numEleLoads = 0;
}

int
LadrunoForceBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  if (numEleLoads == sizeEleLoads) {
    ElementalLoad **theNextEleLoads = new ElementalLoad *[sizeEleLoads + 1];
    double *theNextEleLoadFactors = new double[sizeEleLoads + 1];
    for (int i = 0; i < numEleLoads; i++) {
      theNextEleLoads[i] = eleLoads[i];
      theNextEleLoadFactors[i] = eleLoadFactors[i];
    }
    delete [] eleLoads;
    delete [] eleLoadFactors;
    eleLoads = theNextEleLoads;
    eleLoadFactors = theNextEleLoadFactors;
    sizeEleLoads += 1;
  }

  eleLoadFactors[numEleLoads] = loadFactor;
  eleLoads[numEleLoads] = theLoad;
  numEleLoads++;
  postRestoreSync = false;
  eleLoadDirty = true;

  return 0;
}

void
LadrunoForceBeamColumn3d::computeSectionForces(Vector &sp, int isec)
{
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

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0) * loadFactor;
      double wz = data(1) * loadFactor;
      double wa = data(2) * loadFactor;

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
        case SECTION_RESPONSE_MY:
          sp(ii) += wz * 0.5 * x * (L - x);
          break;
        case SECTION_RESPONSE_VZ:
          sp(ii) += wz * (0.5 * L - x);
          break;
        default:
          break;
        }
      }
    }
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wy = data(0) * loadFactor;
      double wz = data(1) * loadFactor;
      double wa = data(2) * loadFactor;
      double a = data(3) * L;
      double b = data(4) * L;
      if (b <= a || a < 0.0 || b > L)
        continue;
      double wyb = data(5) * loadFactor;
      double wzb = data(6) * loadFactor;
      double wab = data(7) * loadFactor;
      double Fa = wa * (b - a) + 0.5 * (wab - wa) * (b - a);
      double Fy = wy * (b - a);
      double Fz = wz * (b - a);
      double c = a + 0.5 * (b - a);
      double VyI = Fy * (1 - c / L);
      double VyJ = Fy * c / L;
      double VzI = Fz * (1 - c / L);
      double VzJ = Fz * c / L;
      Fy = 0.5 * (wyb - wy) * (b - a);
      Fz = 0.5 * (wzb - wz) * (b - a);
      c = a + 2.0 / 3.0 * (b - a);
      VyI += Fy * (1 - c / L);
      VyJ += Fy * c / L;
      VzI += Fz * (1 - c / L);
      VzJ += Fz * c / L;

      for (int ii = 0; ii < order; ii++) {
        if (x <= a) {
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            sp(ii) += Fa;
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) -= VyI * x;
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) += VzI * x;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) -= VyI;
            break;
          case SECTION_RESPONSE_VZ:
            sp(ii) += VzI;
            break;
          default:
            break;
          }
        }
        else if (x >= b) {
          switch (code(ii)) {
          case SECTION_RESPONSE_MZ:
            sp(ii) += VyJ * (x - L);
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) -= VzJ * (x - L);
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) += VyJ;
            break;
          case SECTION_RESPONSE_VZ:
            sp(ii) -= VzJ;
            break;
          default:
            break;
          }
        }
        else {
          double wyy = wy + (wyb - wy) / (b - a) * (x - a);
          double wzz = wz + (wzb - wz) / (b - a) * (x - a);
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            sp(ii) += Fa - wa * (x - a) - 0.5 * (wab - wa) / (b - a) * (x - a) * (x - a);
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) += -VyI * x + 0.5 * wy * (x - a) * (x - a) + 0.5 * (wyy - wy) * (x - a) * (x - a) / 3.0;
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) += VzI * x - 0.5 * wz * (x - a) * (x - a) - 0.5 * (wzz - wz) * (x - a) * (x - a) / 3.0;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) += -VyI + wy * (x - a) + 0.5 * (wyy - wy) * (x - a);
            break;
          case SECTION_RESPONSE_VZ:
            sp(ii) -= -VzI + wz * (x - a) + 0.5 * (wzz - wz) * (x - a);
            break;
          default:
            break;
          }
        }
      }
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0) * loadFactor;
      double Pz = data(1) * loadFactor;
      double N = data(2) * loadFactor;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      double a = aOverL * L;
      double Vy1 = Py * (1.0 - aOverL);
      double Vy2 = Py * aOverL;
      double Vz1 = Pz * (1.0 - aOverL);
      double Vz2 = Pz * aOverL;

      for (int ii = 0; ii < order; ii++) {
        if (x <= a) {
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            sp(ii) += N;
            break;
          case SECTION_RESPONSE_MZ:
            sp(ii) -= x * Vy1;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) -= Vy1;
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) += x * Vz1;
            break;
          case SECTION_RESPONSE_VZ:
            sp(ii) += Vz1;
            break;
          default:
            break;
          }
        }
        else {
          switch (code(ii)) {
          case SECTION_RESPONSE_MZ:
            sp(ii) -= (L - x) * Vy2;
            break;
          case SECTION_RESPONSE_VY:
            sp(ii) += Vy2;
            break;
          case SECTION_RESPONSE_MY:
            sp(ii) += (L - x) * Vz2;
            break;
          case SECTION_RESPONSE_VZ:
            sp(ii) -= Vz2;
            break;
          default:
            break;
          }
        }
      }
    }
  }
}

void
LadrunoForceBeamColumn3d::computeSectionForceSensitivity(Vector &dspdh, int isec,
                                                         int gradNumber)
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

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0);
      double wz = data(1);
      double wa = data(2);

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwzdh = sens(1);
      double dwadh = sens(2);

      for (int ii = 0; ii < order; ii++) {
        switch (code(ii)) {
        case SECTION_RESPONSE_P:
          dspdh(ii) += dwadh * (L - x) + wa * (dLdh - dxdh);
          break;
        case SECTION_RESPONSE_MZ:
          dspdh(ii) += 0.5 * (dwydh * x * (x - L) + wy * (dxdh * (2 * x - L) - x * dLdh));
          break;
        case SECTION_RESPONSE_VY:
          dspdh(ii) += dwydh * (x - 0.5 * L) + wy * (dxdh - 0.5 * dLdh);
          break;
        case SECTION_RESPONSE_MY:
          dspdh(ii) += 0.5 * (dwzdh * x * (L - x) + wz * (dxdh * (L - 2 * x) + x * dLdh));
          break;
        case SECTION_RESPONSE_VZ:
          dspdh(ii) += dwzdh * (0.5 * L - x) + wz * (0.5 * dLdh - dxdh);
          break;
        default:
          break;
        }
      }
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0);
      double Pz = data(1);
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
        continue;

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPydh = sens(0);
      double dPzdh = sens(1);
      double dNdh = sens(2);
      double daLdh = sens(3);

      double a = aOverL * L;

      double Vy1 = Py * (1.0 - aOverL);
      double Vy2 = Py * aOverL;
      double dVy1dh = Py * (0.0 - daLdh) + dPydh * (1.0 - aOverL);
      double dVy2dh = Py * daLdh + dPydh * aOverL;

      double Vz1 = Pz * (1.0 - aOverL);
      double Vz2 = Pz * aOverL;
      double dVz1dh = Pz * (0.0 - daLdh) + dPzdh * (1.0 - aOverL);
      double dVz2dh = Pz * daLdh + dPzdh * aOverL;

      for (int ii = 0; ii < order; ii++) {
        if (x <= a) {
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
            dspdh(ii) += dNdh;
            break;
          case SECTION_RESPONSE_MZ:
            dspdh(ii) -= (dxdh * Vy1 + x * dVy1dh);
            break;
          case SECTION_RESPONSE_VY:
            dspdh(ii) -= dVy1dh;
            break;
          case SECTION_RESPONSE_MY:
            dspdh(ii) += (dxdh * Vz1 + x * dVz1dh);
            break;
          case SECTION_RESPONSE_VZ:
            dspdh(ii) -= dVz1dh;
            break;
          default:
            break;
          }
        }
        else {
          switch (code(ii)) {
          case SECTION_RESPONSE_MZ:
            dspdh(ii) -= (dLdh - dxdh) * Vy2 + (L - x) * dVy2dh;
            break;
          case SECTION_RESPONSE_VY:
            dspdh(ii) += dVy2dh;
            break;
          case SECTION_RESPONSE_MY:
            dspdh(ii) += (dLdh - dxdh) * Vz2 + (L - x) * dVz2dh;
            break;
          case SECTION_RESPONSE_VZ:
            dspdh(ii) += dVz2dh;
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
LadrunoForceBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
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
    load(2) -= m * Raccel1(2);
    load(6) -= m * Raccel2(0);
    load(7) -= m * Raccel2(1);
    load(8) -= m * Raccel2(2);
  } else {
    static thread_local Vector Raccel(12);
    for (int i = 0; i < 6; i++) {
      Raccel(i) = Raccel1(i);
      Raccel(i + 6) = Raccel2(i);
    }
    load.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

const Vector &
LadrunoForceBeamColumn3d::getResistingForceIncInertia()
{
  static thread_local Vector res(12);

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
      res(2) += m * accel1(2);
      res(6) += m * accel2(0);
      res(7) += m * accel2(1);
      res(8) += m * accel2(2);
    } else {
      static thread_local Vector accel(12);
      for (int i = 0; i < 6; i++) {
        accel(i) = accel1(i);
        accel(i + 6) = accel2(i);
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
LadrunoForceBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  int i, j, k;
  int loc = 0;

  static thread_local ID idData(29);
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;
  idData(6) = (isTorsion) ? 1 : 0;
  idData(13) = maxSubdivisions;
  idData(14) = cMass + 2 * lchMode;

  idData(7) = crdTransf->getClassTag();
  int crdTransfDbTag = crdTransf->getDbTag();
  if (crdTransfDbTag == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag != 0)
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(8) = crdTransfDbTag;

  idData(9) = beamIntegr->getClassTag();
  int beamIntegrDbTag = beamIntegr->getDbTag();
  if (beamIntegrDbTag == 0) {
    beamIntegrDbTag = theChannel.getDbTag();
    if (beamIntegrDbTag != 0)
      beamIntegr->setDbTag(beamIntegrDbTag);
  }
  idData(10) = beamIntegrDbTag;

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

  idData(15) = hingeOn + 2 * dualHingeMode + (useImplex ? 4 : 0);
  UniaxialMaterial *uniList[4] = {theHingeZ, theHingeY, theHingeZJ, theHingeYJ};
  for (int m = 0; m < 4; m++) {
    if (uniList[m] != 0) {
      idData(16 + 2 * m) = uniList[m]->getClassTag();
      int hDb = uniList[m]->getDbTag();
      if (hDb == 0) {
        hDb = theChannel.getDbTag();
        if (hDb != 0)
          uniList[m]->setDbTag(hDb);
      }
      idData(17 + 2 * m) = hDb;
    } else {
      idData(16 + 2 * m) = 0;
      idData(17 + 2 * m) = 0;
    }
  }
  NDMaterial *ndList[2] = {theHingeC, theHingeCJ};
  for (int m = 0; m < 2; m++) {
    if (ndList[m] != 0) {
      idData(24 + 2 * m) = ndList[m]->getClassTag();
      int hDb = ndList[m]->getDbTag();
      if (hDb == 0) {
        hDb = theChannel.getDbTag();
        if (hDb != 0)
          ndList[m]->setDbTag(hDb);
      }
      idData(25 + 2 * m) = hDb;
    } else {
      idData(24 + 2 * m) = 0;
      idData(25 + 2 * m) = 0;
    }
  }
  idData(28) = vRestoreValid ? 1 : 0;

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "LadrunoForceBeamColumn3d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn3d::sendSelf() - failed to send crdTransf\n";
    return -1;
  }

  if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn3d::sendSelf() - failed to send beamIntegr\n";
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
    opserr << "LadrunoForceBeamColumn3d::sendSelf() - failed to send section IDs\n";
    return -1;
  }

  for (j = 0; j < numSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoForceBeamColumn3d::sendSelf() - section " << j << " failed to send itself\n";
      return -1;
    }
  }

  int secDefSize = 0;
  for (i = 0; i < numSections; i++)
    secDefSize += sections[i]->getOrder();

  // Extra doubles: kvbulkCommit (36), vCommitBulk (6), vRestore (6),
  // hingeXiI, hingeXiJ, 4x jumpCommit, 4x mscaleCommit, 24x hingeDualKvCommit,
  // 16x hingeDualKaaInvCommit, pmP0, pmAlpha, eleImplexErrCommit, GAvy, GAvz -> 36+6+6+2+4+4+24+16+5 = 103
  Vector dData(3 + NEBD + NEBD * NEBD + secDefSize + 5 + 103);
  loc = 0;

  dData(loc++) = rho;
  dData(loc++) = tol;
  dData(loc++) = subdivideFactor;

  for (i = 0; i < NEBD; i++)
    dData(loc++) = Secommit(i);

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      dData(loc++) = kvcommit(i, j);

  for (k = 0; k < numSections; k++)
    for (i = 0; i < sections[k]->getOrder(); i++)
      dData(loc++) = (vscommit[k])(i);

  dData(loc++) = alphaM;
  dData(loc++) = betaK;
  dData(loc++) = betaK0;
  dData(loc++) = betaKc;
  dData(loc++) = userLch;

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      dData(loc++) = kvbulkCommit(i, j);
  for (i = 0; i < NEBD; i++)
    dData(loc++) = vCommitBulk(i);
  for (i = 0; i < NEBD; i++)
    dData(loc++) = vRestore[i];
  dData(loc++) = hingeXiI;
  dData(loc++) = hingeXiJ;
  dData(loc++) = hingeJumpCommitZ;
  dData(loc++) = hingeJumpCommitY;
  dData(loc++) = hingeJumpCommitZJ;
  dData(loc++) = hingeJumpCommitYJ;
  dData(loc++) = hingeMscaleCommitZ;
  dData(loc++) = hingeMscaleCommitY;
  dData(loc++) = hingeMscaleCommitZJ;
  dData(loc++) = hingeMscaleCommitYJ;
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      dData(loc++) = hingeDualKvCommit[a][m];
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      dData(loc++) = hingeDualKaaInvCommit[m][n];
  dData(loc++) = pmP0;
  dData(loc++) = pmAlpha;
  dData(loc++) = eleImplexErrCommit;
  dData(loc++) = GAvy;
  dData(loc++) = GAvz;

  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
    opserr << "LadrunoForceBeamColumn3d::sendSelf() - failed to send Vector data\n";
    return -1;
  }

  for (int m = 0; m < 4; m++) {
    if (uniList[m] != 0 && uniList[m]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoForceBeamColumn3d::sendSelf -- could not send uniaxial hinge " << m << "\n";
      return -1;
    }
  }
  for (int m = 0; m < 2; m++) {
    if (ndList[m] != 0 && ndList[m]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoForceBeamColumn3d::sendSelf -- could not send biaxial hinge " << m << "\n";
      return -1;
    }
  }

  if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoForceBeamColumn3d::sendSelf -- could not send Damping\n";
    return -1;
  }

  return 0;
}

int
LadrunoForceBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  int i, j, k;

  static thread_local ID idData(29);
  if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
    opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters = idData(4);
  initialFlag = idData(5);
  isTorsion = (idData(6) == 1) ? true : false;
  maxSubdivisions = idData(13);
  cMass = idData(14) & 1;
  lchMode = idData(14) >> 1;
  hingeOn = idData(15) & 1;
  dualHingeMode = (idData(15) >> 1) & 1;
  useImplex = ((idData(15) >> 2) & 1) != 0;
  vRestoreValid = (idData(28) != 0);

  int crdTransfClassTag = idData(7);
  int crdTransfDbTag = idData(8);
  int beamIntegrClassTag = idData(9);
  int beamIntegrDbTag = idData(10);

  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
    if (crdTransf != 0)
      delete crdTransf;
    crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);
    if (crdTransf == 0) {
      opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to obtain CrdTransf\n";
      return -2;
    }
  }
  crdTransf->setDbTag(crdTransfDbTag);
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to recv crdTransf\n";
    return -3;
  }

  if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntegrClassTag) {
    if (beamIntegr != 0)
      delete beamIntegr;
    beamIntegr = theBroker.getNewBeamIntegration(beamIntegrClassTag);
    if (beamIntegr == 0) {
      opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to obtain BeamIntegration\n";
      return -2;
    }
  }
  beamIntegr->setDbTag(beamIntegrDbTag);
  if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to recv beamIntegr\n";
    return -3;
  }

  ID idSections(2 * idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
    opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to recv section IDs\n";
    return -1;
  }

  if (numSections != idData(3)) {
    if (numSections != 0) {
      for (int ii = 0; ii < numSections; ii++)
        delete sections[ii];
      delete [] sections;
    }
    numSections = idData(3);
    if (vscommit != 0) delete [] vscommit;
    if (fs != 0) delete [] fs;
    if (vs != 0) delete [] vs;
    if (Ssr != 0) delete [] Ssr;
    if (fscommit != 0) delete [] fscommit;
    if (Ssrcommit != 0) delete [] Ssrcommit;

    vscommit = new Vector[numSections];
    fs = new Matrix[numSections];
    vs = new Vector[numSections];
    Ssr = new Vector[numSections];
    fscommit = new Matrix[numSections];
    Ssrcommit = new Vector[numSections];
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

  int secDefSize = 0;
  for (int ii = 0; ii < numSections; ii++)
    secDefSize += sections[ii]->getOrder();

  Vector dData(3 + NEBD + NEBD * NEBD + secDefSize + 5 + 103);
  if (theChannel.recvVector(dbTag, commitTag, dData) < 0) {
    opserr << "LadrunoForceBeamColumn3d::recvSelf() - failed to recv Vector data\n";
    return -1;
  }

  loc = 0;
  rho = dData(loc++);
  tol = dData(loc++);
  subdivideFactor = dData(loc++);

  for (i = 0; i < NEBD; i++)
    Secommit(i) = dData(loc++);

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      kvcommit(i, j) = dData(loc++);

  kv = kvcommit;
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

  for (i = 0; i < NEBD; i++)
    for (j = 0; j < NEBD; j++)
      kvbulkCommit(i, j) = dData(loc++);
  kvbulk = kvbulkCommit;

  for (i = 0; i < NEBD; i++)
    vCommitBulk(i) = dData(loc++);
  for (i = 0; i < NEBD; i++)
    vRestore[i] = dData(loc++);
  hingeXiI = dData(loc++);
  hingeXiJ = dData(loc++);
  hingeJumpZ  = hingeJumpCommitZ  = dData(loc++);
  hingeJumpY  = hingeJumpCommitY  = dData(loc++);
  hingeJumpZJ = hingeJumpCommitZJ = dData(loc++);
  hingeJumpYJ = hingeJumpCommitYJ = dData(loc++);
  hingeMscaleZ  = hingeMscaleCommitZ  = dData(loc++);
  hingeMscaleY  = hingeMscaleCommitY  = dData(loc++);
  hingeMscaleZJ = hingeMscaleCommitZJ = dData(loc++);
  hingeMscaleYJ = hingeMscaleCommitYJ = dData(loc++);
  for (int a = 0; a < 6; a++)
    for (int m = 0; m < 4; m++)
      hingeDualKv[a][m] = hingeDualKvCommit[a][m] = dData(loc++);
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++)
      hingeDualKaaInv[m][n] = hingeDualKaaInvCommit[m][n] = dData(loc++);
  pmP0 = dData(loc++);
  pmAlpha = dData(loc++);
  eleImplexErr = eleImplexErrCommit = dData(loc++);
  GAvy = dData(loc++);
  GAvz = dData(loc++);

  if (Ki != 0) {
    delete Ki;
    Ki = 0;
  }

  UniaxialMaterial **uniPtrs[4] = {&theHingeZ, &theHingeY, &theHingeZJ, &theHingeYJ};
  for (int m = 0; m < 4; m++) {
    int hClassTag = idData(16 + 2 * m);
    int hDbTag = idData(17 + 2 * m);
    if (hClassTag != 0) {
      if (*uniPtrs[m] == 0 || (*uniPtrs[m])->getClassTag() != hClassTag) {
        if (*uniPtrs[m] != 0) delete *uniPtrs[m];
        *uniPtrs[m] = theBroker.getNewUniaxialMaterial(hClassTag);
        if (*uniPtrs[m] == 0) return -1;
      }
      (*uniPtrs[m])->setDbTag(hDbTag);
      if ((*uniPtrs[m])->recvSelf(commitTag, theChannel, theBroker) < 0)
        return -1;
    } else if (*uniPtrs[m] != 0) {
      delete *uniPtrs[m];
      *uniPtrs[m] = 0;
    }
  }

  NDMaterial **ndPtrs[2] = {&theHingeC, &theHingeCJ};
  for (int m = 0; m < 2; m++) {
    int hClassTag = idData(24 + 2 * m);
    int hDbTag = idData(25 + 2 * m);
    if (hClassTag != 0) {
      if (*ndPtrs[m] == 0 || (*ndPtrs[m])->getClassTag() != hClassTag) {
        if (*ndPtrs[m] != 0) delete *ndPtrs[m];
        *ndPtrs[m] = theBroker.getNewNDMaterial(hClassTag);
        if (*ndPtrs[m] == 0) return -1;
      }
      (*ndPtrs[m])->setDbTag(hDbTag);
      if ((*ndPtrs[m])->recvSelf(commitTag, theChannel, theBroker) < 0)
        return -1;
    } else if (*ndPtrs[m] != 0) {
      delete *ndPtrs[m];
      *ndPtrs[m] = 0;
    }
  }

  initialFlag = 2;
  postRestoreSync = (hingeOn != 0);
  eleLoadDirty = false;

  int dmpTag = (int)idData(11);
  if (dmpTag) {
    if (theDamping == 0 || theDamping->getClassTag() != dmpTag) {
      if (theDamping != 0) delete theDamping;
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) return -1;
    }
    theDamping->setDbTag((int)idData(12));
    if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0)
      return -1;
  } else {
    if (theDamping) {
      delete theDamping;
      theDamping = 0;
    }
  }

  return 0;
}

int
LadrunoForceBeamColumn3d::getInitialFlexibility(Matrix &fe, bool includeHinges)
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
    int ii, jj;
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        for (jj = 0; jj < order; jj++)
          fb(jj, 0) += fSec(jj, ii) * wtL;
        break;
      case SECTION_RESPONSE_MZ:
        for (jj = 0; jj < order; jj++) {
          tmp = fSec(jj, ii) * wtL;
          fb(jj, 1) += xL1 * tmp;
          fb(jj, 2) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (jj = 0; jj < order; jj++) {
          tmp = oneOverL * fSec(jj, ii) * wtL;
          fb(jj, 1) += tmp;
          fb(jj, 2) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (jj = 0; jj < order; jj++) {
          tmp = fSec(jj, ii) * wtL;
          fb(jj, 3) += xL1 * tmp;
          fb(jj, 4) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (jj = 0; jj < order; jj++) {
          tmp = oneOverL * fSec(jj, ii) * wtL;
          fb(jj, 3) += tmp;
          fb(jj, 4) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (jj = 0; jj < order; jj++)
          fb(jj, 5) += fSec(jj, ii) * wtL;
        break;
      default:
        break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        for (jj = 0; jj < NEBD; jj++)
          fe(0, jj) += fb(ii, jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = fb(ii, jj);
          fe(1, jj) += xL1 * tmp;
          fe(2, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = oneOverL * fb(ii, jj);
          fe(1, jj) += tmp;
          fe(2, jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = fb(ii, jj);
          fe(3, jj) += xL1 * tmp;
          fe(4, jj) += xL * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = oneOverL * fb(ii, jj);
          fe(3, jj) += tmp;
          fe(4, jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (jj = 0; jj < NEBD; jj++)
          fe(5, jj) += fb(ii, jj);
        break;
      default:
        break;
      }
    }
  }

  if (GAvy > 0.0) {
    double fShearY = oneOverL / GAvy;
    fe(1, 1) += fShearY;
    fe(1, 2) += fShearY;
    fe(2, 1) += fShearY;
    fe(2, 2) += fShearY;
  }
  if (GAvz > 0.0) {
    double fShearZ = oneOverL / GAvz;
    fe(3, 3) += fShearZ;
    fe(3, 4) += fShearZ;
    fe(4, 3) += fShearZ;
    fe(4, 4) += fShearZ;
  }

  if (!isTorsion)
    fe(5, 5) = DefaultLoverGJ;

  if (includeHinges && hingeOn) {
    const double betaI  = 1.0 - hingeXiI;
    const double gammaI = hingeXiI;
    const double betaJ  = 1.0 - hingeXiJ;
    const double gammaJ = hingeXiJ;
    const double fPM = (pmP0 > 0.0) ? lfbc3dScalePM(Se(0), pmP0, pmAlpha) : 1.0;
    if (theHingeC != 0) {
      if (dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC) == 0) {
        const Matrix &C0 = theHingeC->getInitialTangent();
        if (C0.noRows() >= 2 && C0.noCols() >= 2) {
          double det = fPM * (C0(0, 0) * C0(1, 1) - C0(0, 1) * C0(1, 0));
          if (fabs(det) > 1e-30) {
            double invC00 =  C0(1, 1) / det;
            double invC01 = -C0(0, 1) / det;
            double invC10 = -C0(1, 0) / det;
            double invC11 =  C0(0, 0) / det;
            fe(1, 1) += betaI * betaI * invC00;
            fe(1, 2) += betaI * gammaI * invC00;
            fe(2, 1) += gammaI * betaI * invC00;
            fe(2, 2) += gammaI * gammaI * invC00;
            fe(1, 3) += betaI * betaI * invC01;
            fe(1, 4) += betaI * gammaI * invC01;
            fe(2, 3) += gammaI * betaI * invC01;
            fe(2, 4) += gammaI * gammaI * invC01;
            fe(3, 1) += betaI * betaI * invC10;
            fe(3, 2) += betaI * gammaI * invC10;
            fe(4, 1) += gammaI * betaI * invC10;
            fe(4, 2) += gammaI * gammaI * invC10;
            fe(3, 3) += betaI * betaI * invC11;
            fe(3, 4) += betaI * gammaI * invC11;
            fe(4, 3) += gammaI * betaI * invC11;
            fe(4, 4) += gammaI * gammaI * invC11;
          }
        }
      }
    } else {
      if (theHingeZ != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ) == 0) {
        double k0 = fPM * theHingeZ->getInitialTangent();
        if (k0 > 0.0) {
          double f0 = 1.0 / k0;
          fe(1, 1) += betaI * betaI * f0;
          fe(1, 2) += betaI * gammaI * f0;
          fe(2, 1) += gammaI * betaI * f0;
          fe(2, 2) += gammaI * gammaI * f0;
        }
      }
      if (theHingeY != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHingeY) == 0) {
        double k0 = fPM * theHingeY->getInitialTangent();
        if (k0 > 0.0) {
          double f0 = 1.0 / k0;
          fe(3, 3) += betaI * betaI * f0;
          fe(3, 4) += betaI * gammaI * f0;
          fe(4, 3) += gammaI * betaI * f0;
          fe(4, 4) += gammaI * gammaI * f0;
        }
      }
    }
    if (theHingeCJ != 0) {
      if (dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeCJ) == 0) {
        const Matrix &C0 = theHingeCJ->getInitialTangent();
        if (C0.noRows() >= 2 && C0.noCols() >= 2) {
          double det = fPM * (C0(0, 0) * C0(1, 1) - C0(0, 1) * C0(1, 0));
          if (fabs(det) > 1e-30) {
            double invC00 =  C0(1, 1) / det;
            double invC01 = -C0(0, 1) / det;
            double invC10 = -C0(1, 0) / det;
            double invC11 =  C0(0, 0) / det;
            fe(1, 1) += betaJ * betaJ * invC00;
            fe(1, 2) += betaJ * gammaJ * invC00;
            fe(2, 1) += gammaJ * betaJ * invC00;
            fe(2, 2) += gammaJ * gammaJ * invC00;
            fe(1, 3) += betaJ * betaJ * invC01;
            fe(1, 4) += betaJ * gammaJ * invC01;
            fe(2, 3) += gammaJ * betaJ * invC01;
            fe(2, 4) += gammaJ * gammaJ * invC01;
            fe(3, 1) += betaJ * betaJ * invC10;
            fe(3, 2) += betaJ * gammaJ * invC10;
            fe(4, 1) += gammaJ * betaJ * invC10;
            fe(4, 2) += gammaJ * gammaJ * invC10;
            fe(3, 3) += betaJ * betaJ * invC11;
            fe(3, 4) += betaJ * gammaJ * invC11;
            fe(4, 3) += gammaJ * betaJ * invC11;
            fe(4, 4) += gammaJ * gammaJ * invC11;
          }
        }
      }
    } else {
      if (theHingeZJ != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ) == 0) {
        double k0 = fPM * theHingeZJ->getInitialTangent();
        if (k0 > 0.0) {
          double f0 = 1.0 / k0;
          fe(1, 1) += betaJ * betaJ * f0;
          fe(1, 2) += betaJ * gammaJ * f0;
          fe(2, 1) += gammaJ * betaJ * f0;
          fe(2, 2) += gammaJ * gammaJ * f0;
        }
      }
      if (theHingeYJ != 0 && dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ) == 0) {
        double k0 = fPM * theHingeYJ->getInitialTangent();
        if (k0 > 0.0) {
          double f0 = 1.0 / k0;
          fe(3, 3) += betaJ * betaJ * f0;
          fe(3, 4) += betaJ * gammaJ * f0;
          fe(4, 3) += gammaJ * betaJ * f0;
          fe(4, 4) += gammaJ * gammaJ * f0;
        }
      }
    }
  }

  return 0;
}

int
LadrunoForceBeamColumn3d::getInitialDeformations(Vector &v0)
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
      case SECTION_RESPONSE_MY:
        v0(3) += xL1 * dei;
        v0(4) += xL * dei;
        break;
      case SECTION_RESPONSE_VZ:
        tmp = oneOverL * dei;
        v0(3) += tmp;
        v0(4) += tmp;
        break;
      default:
        break;
      }
    }
  }

  return 0;
}

void
LadrunoForceBeamColumn3d::compSectionDisplacements(Vector sectionCoords[],
                                                   Vector sectionDispls[]) const
{
  static thread_local Vector ub(NEBD);
  ub = crdTransf->getBasicTrialDisp();

  double L = crdTransf->getInitialLength();

  static thread_local double pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, pts);

  int i;
  double xi;

  Matrix ls(numSections, numSections);
  getCBDIinfluenceMatrix(numSections, pts, L, ls);

  Vector kappa_y(numSections);
  Vector kappa_z(numSections);
  static thread_local Vector vs_loc;

  for (i = 0; i < numSections; i++) {
    int sectionKey1 = 0;
    int sectionKey2 = 0;
    const ID &code = sections[i]->getType();
    for (int j = 0; j < code.Size(); j++) {
      if (code(j) == SECTION_RESPONSE_MZ)
        sectionKey1 = j;
      if (code(j) == SECTION_RESPONSE_MY)
        sectionKey2 = j;
    }
    vs_loc = sections[i]->getSectionDeformation();
    kappa_z(i) = vs_loc(sectionKey1);
    kappa_y(i) = vs_loc(sectionKey2);
  }

  Vector v(numSections), w(numSections);
  static thread_local Vector xl(NDM), uxb(NDM);

  v.addMatrixVector(0.0, ls, kappa_z, 1.0);
  w.addMatrixVector(0.0, ls, kappa_y, -1.0);

  for (i = 0; i < numSections; i++) {
    xi = pts[i];
    xl(0) = xi * L;
    xl(1) = 0;
    xl(2) = 0;

    sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

    uxb(0) = xi * ub(0);
    uxb(1) = v(i);
    uxb(2) = w(i);

    sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
  }
}

void
LadrunoForceBeamColumn3d::Print(OPS_Stream &s, int flag)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nElement: " << this->getTag() << " Type: LadrunoForceBeamColumn3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
    beamIntegr->Print(s, flag);
    double P = Secommit(0);
    double MZ1 = Secommit(1);
    double MZ2 = Secommit(2);
    double MY1 = Secommit(3);
    double MY2 = Secommit(4);
    double L = crdTransf->getInitialLength();
    double VY = (MZ1 + MZ2) / L;
    double VZ = -(MY1 + MY2) / L;
    double T = Secommit(5);

    double p0[5]; p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);

    s << "\tEnd 1 Forces (P MZ VY MY VZ T): "
      << -P + p0[0] << " " << MZ1 << " " << VY + p0[1] << " "
      << MY1 << " " << VZ + p0[3] << " " << -T << endln;
    s << "\tEnd 2 Forces (P MZ VY MY VZ T): "
      << P << " " << MZ2 << " " << -VY + p0[2] << " "
      << MY2 << " " << -VZ + p0[4] << " " << T << endln;

    for (int i = 0; i < numSections; i++)
      sections[i]->Print(s, flag);
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoForceBeamColumn3d\", ";
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

OPS_Stream &
operator<<(OPS_Stream &s, LadrunoForceBeamColumn3d &E)
{
  E.Print(s);
  return s;
}

int
LadrunoForceBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
  static thread_local Vector v1(3);
  static thread_local Vector v2(3);

  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);

  return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response *
LadrunoForceBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoForceBeamColumn3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  if (LadrunoResp::isGlobalForce(argv[0])) {
    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 1, theVector);
  }
  else if (LadrunoResp::isLocalForce(argv[0])) {
    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 2, theVector);
  }
  else if (LadrunoResp::isBasicForce(argv[0])) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 7, Vector(6));
  }
  else if (LadrunoResp::is(argv[0], "basicStiffness")) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Mz_2");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "T");

    theResponse = new ElementResponse(this, 19, Matrix(6, 6));
  }
  else if (theDamping && (strcmp(argv[0], "globalDampingForce") == 0 || strcmp(argv[0], "globalDampingForces") == 0)) {
    theResponse = new ElementResponse(this, 21, theVector);
  }
  else if (theDamping && (strcmp(argv[0], "localDampingForce") == 0 || strcmp(argv[0], "localDampingForces") == 0)) {
    theResponse = new ElementResponse(this, 22, theVector);
  }
  else if (theDamping && (strcmp(argv[0], "basicDampingForce") == 0 || strcmp(argv[0], "basicDampingForces") == 0)) {
    theResponse = new ElementResponse(this, 23, Vector(6));
  }
  else if (LadrunoResp::isBasicDef(argv[0])) {
    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));
  }
  else if (LadrunoResp::isPlasticDef(argv[0])) {
    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));
  }
  else if (strcmp(argv[0], "inflectionPoint") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(2));
  }
  else if (strcmp(argv[0], "tangentDrift") == 0) {
    theResponse = new ElementResponse(this, 6, Vector(4));
  }
  else if (LadrunoResp::isRayleighForce(argv[0])) {
    theResponse = new ElementResponse(this, 12, theVector);
  }
  else if (strcmp(argv[0], "sections") == 0) {
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
    if (argc > 1 && strcmp(argv[1], "local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(numSections, 3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(numSections, 3));
  }
  else if (strcmp(argv[0], "cbdiDisplacements") == 0) {
    theResponse = new ElementResponse(this, 112, Matrix(1, 3));
  }
  else if (strcmp(argv[0], "sectionX") == 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);
      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = fabs(xi[i] - sectionLoc);
          sectionNum = i;
        }
      }
      output.tag("GaussPointOutput");
      output.attr("number", sectionNum + 1);
      output.attr("eta", xi[sectionNum] * L);
      theResponse = sections[sectionNum]->setResponse(&argv[2], argc - 2, output);
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
    theResponse = new ElementResponse(this, 2000, 0.0);
  }
  else if (LadrunoResp::isLocalAxes(argv[0])) {
    theResponse = new ElementResponse(this, 30, Vector(9));
  }
  else if (hingeOn && argc > 1 &&
           (strcmp(argv[0], "hinge") == 0 || strcmp(argv[0], "hingeZ") == 0)) {
    UniaxialMaterial *hm = (theHingeZ != 0) ? theHingeZ : theHingeZJ;
    if (hm != 0)
      theResponse = hm->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeZ != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hingeZI") == 0 || strcmp(argv[0], "hinge1") == 0)) {
    theResponse = theHingeZ->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeZJ != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hingeZJ") == 0 || strcmp(argv[0], "hinge2") == 0)) {
    theResponse = theHingeZJ->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && argc > 1 && strcmp(argv[0], "hingeY") == 0) {
    UniaxialMaterial *hm = (theHingeY != 0) ? theHingeY : theHingeYJ;
    if (hm != 0)
      theResponse = hm->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeY != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeYI") == 0 || strcmp(argv[0], "hingeY1") == 0)) {
    theResponse = theHingeY->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeYJ != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeYJ") == 0 || strcmp(argv[0], "hingeY2") == 0)) {
    theResponse = theHingeYJ->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && argc > 1 &&
           (strcmp(argv[0], "hingeBiaxial") == 0 || strcmp(argv[0], "hingeC") == 0)) {
    NDMaterial *hm = (theHingeC != 0) ? theHingeC : theHingeCJ;
    if (hm != 0)
      theResponse = hm->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeC != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeBiaxialI") == 0 || strcmp(argv[0], "hingeBiaxial1") == 0 || strcmp(argv[0], "hingeCI") == 0)) {
    theResponse = theHingeC->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && theHingeCJ != 0 && argc > 1 &&
           (strcmp(argv[0], "hingeBiaxialJ") == 0 || strcmp(argv[0], "hingeBiaxial2") == 0 || strcmp(argv[0], "hingeCJ") == 0)) {
    theResponse = theHingeCJ->setResponse(&argv[1], argc - 1, output);
  }
  else if (hingeOn && LadrunoResp::is(argv[0], "hingeSummary")) {
    output.tag("ResponseType", "thetaZ_i");
    output.tag("ResponseType", "Mz_i");
    output.tag("ResponseType", "thetaZ_j");
    output.tag("ResponseType", "Mz_j");
    output.tag("ResponseType", "thetaY_i");
    output.tag("ResponseType", "My_i");
    output.tag("ResponseType", "thetaY_j");
    output.tag("ResponseType", "My_j");
    theResponse = new ElementResponse(this, 31, Vector(8));
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
    output.tag("ResponseType", "thetaZ_i");
    output.tag("ResponseType", "thetaY_i");
    theResponse = new ElementResponse(this, 34, Vector(2));
  }
  else if (strcmp(argv[0], "hingeJumpJ") == 0) {
    output.tag("ResponseType", "thetaZ_j");
    output.tag("ResponseType", "thetaY_j");
    theResponse = new ElementResponse(this, 35, Vector(2));
  }

  if (theResponse == 0)
    theResponse = crdTransf->setResponse(argv, argc, output);

  output.endTag();

  if (theResponse == 0)
    return this->Element::setResponse(argv, argc, output);
  return theResponse;
}

int
LadrunoForceBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  if (eleLoadDirty && theNodes[0] != 0 && theNodes[1] != 0)
    this->update();
  static thread_local Vector vp(6);
  static thread_local Matrix fe(6, 6);

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
    static thread_local Vector hs(8);
    hs.Zero();
    if (hingeOn) {
      hs(0) = hingeJumpZ;
      hs(2) = hingeJumpZJ;
      hs(4) = hingeJumpY;
      hs(6) = hingeJumpYJ;
      double fPM = (pmP0 > 0.0) ? lfbc3dScalePM(Se(0), pmP0, pmAlpha) : 1.0;
      if (theHingeC != 0) {
        double f = (dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC) == 0) ? fPM : 1.0;
        const Vector &s = theHingeC->getStress();
        if (s.Size() >= 2) {
          hs(1) = f * s(0);
          hs(5) = f * s(1);
        }
      } else {
        if (theHingeZ != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ) == 0) ? fPM : 1.0;
          hs(1) = f * theHingeZ->getStress();
        }
        if (theHingeY != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeY) == 0) ? fPM : 1.0;
          hs(5) = f * theHingeY->getStress();
        }
      }
      if (theHingeCJ != 0) {
        double f = (dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeCJ) == 0) ? fPM : 1.0;
        const Vector &sJ = theHingeCJ->getStress();
        if (sJ.Size() >= 2) {
          hs(3) = f * sJ(0);
          hs(7) = f * sJ(1);
        }
      } else {
        if (theHingeZJ != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ) == 0) ? fPM : 1.0;
          hs(3) = f * theHingeZJ->getStress();
        }
        if (theHingeYJ != 0) {
          double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ) == 0) ? fPM : 1.0;
          hs(7) = f * theHingeYJ->getStress();
        }
      }
    }
    return eleInfo.setVector(hs);
  }

  else if (responseID == 32 || responseID == 33) {
    double maxErr = 0.0, sumErr = 0.0;
    int cnt = 0;
    UniaxialMaterial *uniList[4] = {theHingeZ, theHingeY, theHingeZJ, theHingeYJ};
    for (int k = 0; k < 4; k++) {
      double e = 0.0;
      if (lfbc3dGetHingeImplexError(uniList[k], e)) {
        if (e > maxErr) maxErr = e;
        sumErr += e;
        cnt++;
      }
    }
    NDMaterial *ndList[2] = {theHingeC, theHingeCJ};
    for (int k = 0; k < 2; k++) {
      if (LadrunoCohesiveHingeBiaxial *ch = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(ndList[k])) {
        double e = ch->getImplexError();
        if (e > maxErr) maxErr = e;
        sumErr += e;
        cnt++;
      }
    }
    if (cnt > 0) {
      double val = (responseID == 32) ? maxErr : (sumErr / cnt);
      return eleInfo.setDouble(val);
    }
    return eleInfo.setDouble(useImplex ? eleImplexErr : 0.0);
  }

  else if (responseID == 34) {
    static thread_local Vector hjI(2);
    hjI(0) = hingeOn ? hingeJumpZ : 0.0;
    hjI(1) = hingeOn ? hingeJumpY : 0.0;
    return eleInfo.setVector(hjI);
  }

  else if (responseID == 35) {
    static thread_local Vector hjJ(2);
    hjJ(0) = hingeOn ? hingeJumpZJ : 0.0;
    hjJ(1) = hingeOn ? hingeJumpYJ : 0.0;
    return eleInfo.setVector(hjJ);
  }

  else if (responseID == 2) {
    double p0[5]; p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);
    double N = Se(0);
    theVector(6) = N;
    theVector(0) = -N + p0[0];

    double T = Se(5);
    theVector(9) = T;
    theVector(3) = -T;

    double M1 = Se(1);
    double M2 = Se(2);
    theVector(5) = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1 + M2) / L;
    theVector(1) = V + p0[1];
    theVector(7) = -V + p0[2];

    M1 = Se(3);
    M2 = Se(4);
    theVector(4) = M1;
    theVector(10) = M2;
    V = (M1 + M2) / L;
    theVector(2) = -V + p0[3];
    theVector(8) = V + p0[4];

    return eleInfo.setVector(theVector);
  }

  else if (responseID == 21)
    return eleInfo.setVector(this->getDampingForce());

  else if (responseID == 22) {
    Vector Sd(NEBD);
    Sd = theDamping->getDampingForce();
    double N = Sd(0);
    theVector(6) = N;
    theVector(0) = -N;

    double T = Sd(5);
    theVector(9) = T;
    theVector(3) = -T;

    double M1 = Sd(1);
    double M2 = Sd(2);
    theVector(5) = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1 + M2) / L;
    theVector(1) = V;
    theVector(7) = -V;

    M1 = Sd(3);
    M2 = Sd(4);
    theVector(4) = M1;
    theVector(10) = M2;
    V = (M1 + M2) / L;
    theVector(2) = -V;
    theVector(8) = V;

    return eleInfo.setVector(theVector);
  }

  else if (responseID == 23)
    return eleInfo.setVector(theDamping->getDampingForce());

  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  else if (responseID == 7)
    return eleInfo.setVector(Se);

  else if (responseID == 19)
    return eleInfo.setMatrix(kv);

  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    Vector v0(6);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return eleInfo.setVector(vp);
  }

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

  else if (responseID == 111 || responseID == 1111) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    Vector kappaz(numSections);
    Vector kappay(numSections);
    for (int i = 0; i < numSections; i++) {
      const ID &code = sections[i]->getType();
      const Vector &e = sections[i]->getSectionDeformation();
      int order = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz(i) += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay(i) += e(j);
      }
    }
    Vector dispsy(numSections);
    Vector dispsz(numSections);
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    beamIntegr->getSectionLocations(numSections, L, pts);
    static thread_local Vector uxb(3);
    static thread_local Vector uxg(3);
    Matrix disps(numSections, 3);
    vp = crdTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i] * vp(0);
      uxb(1) = dispsy(i);
      uxb(2) = dispsz(i);
      if (responseID == 111)
        uxg = crdTransf->getPointGlobalDisplFromBasic(pts[i], uxb);
      else
        uxg = crdTransf->getPointLocalDisplFromBasic(pts[i], uxb);
      disps(i, 0) = uxg(0);
      disps(i, 1) = uxg(1);
      disps(i, 2) = uxg(2);
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
    Vector kappaz(numSections);
    Vector kappay(numSections);
    for (int i = 0; i < numSections; i++) {
      const ID &code = sections[i]->getType();
      const Vector &e = sections[i]->getSectionDeformation();
      int order = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
        if (code(j) == SECTION_RESPONSE_MZ)
          kappaz(i) += e(j);
        if (code(j) == SECTION_RESPONSE_MY)
          kappay(i) += e(j);
      }
    }
    Vector dispsy(1);
    Vector dispsz(1);
    dispsy.addMatrixVector(0.0, ls, kappaz, 1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    static thread_local Vector uxb(3);
    static thread_local Vector uxg(3);
    Matrix disps(1, 3);
    vp = crdTransf->getBasicTrialDisp();
    uxb(0) = pts[0] * vp(0);
    uxb(1) = dispsy(0);
    uxb(2) = dispsz(0);
    uxg = crdTransf->getPointGlobalDisplFromBasic(pts[0], uxb);
    disps(0, 0) = uxg(0);
    disps(0, 1) = uxg(1);
    disps(0, 2) = uxg(2);
    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

  else if (responseID == 5) {
    static thread_local Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;
    double L = crdTransf->getInitialLength();
    if (fabs(Se(1) + Se(2)) > DBL_EPSILON)
      LI(0) = Se(1) / (Se(1) + Se(2)) * L;
    if (fabs(Se(3) + Se(4)) > DBL_EPSILON)
      LI(1) = Se(3) / (Se(3) + Se(4)) * L;
    return eleInfo.setVector(LI);
  }

  else if (responseID == 6) {
    double d2z = 0.0, d2y = 0.0, d3z = 0.0, d3y = 0.0;
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    double LIz = 0.0;
    if (fabs(Se(1) + Se(2)) > DBL_EPSILON)
      LIz = Se(1) / (Se(1) + Se(2)) * L;

    double LIy = 0.0;
    if (fabs(Se(3) + Se(4)) > DBL_EPSILON)
      LIy = Se(3) / (Se(3) + Se(4)) * L;

    for (int i = 0; i < numSections; i++) {
      double x = pts[i] * L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x < LIz) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += vs[i](j);
        d2z += (wts[i] * L) * kappa * (-LIz + x);
      }
      kappa = 0.0;
      if (x < LIy) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += vs[i](j);
        d2y += (wts[i] * L) * kappa * (-LIy + x);
      }
    }

    for (int i = numSections - 1; i >= 0; i--) {
      double x = pts[i] * L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x > LIz) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += vs[i](j);
        d3z += (wts[i] * L) * kappa * (x - LIz);
      }
      kappa = 0.0;
      if (x > LIy) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += vs[i](j);
        d3y += (wts[i] * L) * kappa * (x - LIy);
      }
    }

    static thread_local Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;
    return eleInfo.setVector(d);
  }

  else if (responseID == 2000) {
    double xi[maxNumSections];
    double L = crdTransf->getInitialLength();
    beamIntegr->getSectionWeights(numSections, L, xi);
    double energy = 0.0;
    for (int i = 0; i < numSections; i++)
      energy += sections[i]->getEnergy() * xi[i] * L;
    if (hingeOn) {
      double fPM = (pmP0 > 0.0) ? lfbc3dScalePM(Se(0), pmP0, pmAlpha) : 1.0;
      if (theHingeZ != 0) {
        double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZ) == 0) ? fPM : 1.0;
        energy += f * theHingeZ->getEnergy();
      }
      if (theHingeY != 0) {
        double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeY) == 0) ? fPM : 1.0;
        energy += f * theHingeY->getEnergy();
      }
      if (theHingeC != 0) {
        if (LadrunoCohesiveHingeBiaxial *hc = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeC))
          energy += hc->getEnergy();
      }
      if (theHingeZJ != 0) {
        double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeZJ) == 0) ? fPM : 1.0;
        energy += f * theHingeZJ->getEnergy();
      }
      if (theHingeYJ != 0) {
        double f = (dynamic_cast<LadrunoCohesiveHinge*>(theHingeYJ) == 0) ? fPM : 1.0;
        energy += f * theHingeYJ->getEnergy();
      }
      if (theHingeCJ != 0) {
        if (LadrunoCohesiveHingeBiaxial *hcJ = dynamic_cast<LadrunoCohesiveHingeBiaxial*>(theHingeCJ))
          energy += hcJ->getEnergy();
      }
    }
    return eleInfo.setDouble(energy);
  }

  return this->Element::getResponse(responseID, eleInfo);
}

double
LadrunoForceBeamColumn3d::getCharacteristicLength(void)
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
LadrunoForceBeamColumn3d::getResponseSensitivity(int responseID, int gradNumber,
                                                 Information &eleInfo)
{
  if (responseID == 3) {
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }
  else if (responseID == 7) {
    static thread_local Vector dqdh(6);
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);
    return eleInfo.setVector(dqdh);
  }
  return -1;
}

int
LadrunoForceBeamColumn3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "lch") == 0 || strcmp(argv[0], "Lch") == 0 ||
      strcmp(argv[0], "characteristicLength") == 0) {
    param.setValue(this->getCharacteristicLength());
    return param.addObject(2, this);
  }
  if (strcmp(argv[0], "GAy") == 0 || strcmp(argv[0], "GAvy") == 0 || strcmp(argv[0], "shearY") == 0) {
    param.setValue(GAvy);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0], "GAz") == 0 || strcmp(argv[0], "GAvz") == 0 || strcmp(argv[0], "shearZ") == 0) {
    param.setValue(GAvz);
    return param.addObject(4, this);
  }

  if (hingeOn && argc >= 2) {
    if (strcmp(argv[0], "hingeI") == 0 || strcmp(argv[0], "hingeZI") == 0 || strcmp(argv[0], "hinge1") == 0) {
      int res = -1;
      if (theHingeZ != 0) res = theHingeZ->setParameter(&argv[1], argc - 1, param);
      else if (theHingeC != 0) res = theHingeC->setParameter(&argv[1], argc - 1, param);
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeJ") == 0 || strcmp(argv[0], "hingeZJ") == 0 || strcmp(argv[0], "hinge2") == 0) {
      int res = -1;
      if (theHingeZJ != 0) res = theHingeZJ->setParameter(&argv[1], argc - 1, param);
      else if (theHingeCJ != 0) res = theHingeCJ->setParameter(&argv[1], argc - 1, param);
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hinge") == 0 || strcmp(argv[0], "hingeZ") == 0 ||
        strcmp(argv[0], "hingeIJ") == 0 || strcmp(argv[0], "hingeZIJ") == 0) {
      int res = -1;
      if (theHingeZ != 0) {
        int ok = theHingeZ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeZJ != 0) {
        int ok = theHingeZJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeC != 0) {
        int ok = theHingeC->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeCJ != 0) {
        int ok = theHingeCJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeY") == 0 || strcmp(argv[0], "hingeYIJ") == 0) {
      int res = -1;
      if (theHingeY != 0) {
        int ok = theHingeY->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeYJ != 0) {
        int ok = theHingeYJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeYI") == 0 || strcmp(argv[0], "hingeY1") == 0) {
      int res = (theHingeY != 0) ? theHingeY->setParameter(&argv[1], argc - 1, param) : -1;
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeYJ") == 0 || strcmp(argv[0], "hingeY2") == 0) {
      int res = (theHingeYJ != 0) ? theHingeYJ->setParameter(&argv[1], argc - 1, param) : -1;
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeBiaxial") == 0 || strcmp(argv[0], "hingeC") == 0 ||
        strcmp(argv[0], "hingeBiaxialIJ") == 0 || strcmp(argv[0], "hingeCIJ") == 0) {
      int res = -1;
      if (theHingeC != 0) {
        int ok = theHingeC->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (theHingeCJ != 0) {
        int ok = theHingeCJ->setParameter(&argv[1], argc - 1, param);
        if (ok != -1) res = ok;
      }
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeBiaxialI") == 0 || strcmp(argv[0], "hingeBiaxial1") == 0 ||
        strcmp(argv[0], "hingeCI") == 0) {
      int res = (theHingeC != 0) ? theHingeC->setParameter(&argv[1], argc - 1, param) : -1;
      if (res != -1) param.addObject(5, this);
      return res;
    }
    if (strcmp(argv[0], "hingeBiaxialJ") == 0 || strcmp(argv[0], "hingeBiaxial2") == 0 ||
        strcmp(argv[0], "hingeCJ") == 0) {
      int res = (theHingeCJ != 0) ? theHingeCJ->setParameter(&argv[1], argc - 1, param) : -1;
      if (res != -1) param.addObject(5, this);
      return res;
    }
  }

  if (strstr(argv[0], "damp") != 0) {
    if (argc < 2 || !theDamping)
      return -1;
    return theDamping->setParameter(&argv[1], argc - 1, param);
  }

  if (strstr(argv[0], "sectionX") != 0) {
    if (argc < 3)
      return -1;
    float sectionLoc = atof(argv[1]);
    double xi[maxNumSections];
    double L = crdTransf->getInitialLength();
    beamIntegr->getSectionLocations(numSections, L, xi);
    sectionLoc /= L;

    float minDistance = fabs(xi[0] - sectionLoc);
    int sectionNum = 0;
    for (int i = 1; i < numSections; i++) {
      if (fabs(xi[i] - sectionLoc) < minDistance) {
        minDistance = fabs(xi[i] - sectionLoc);
        sectionNum = i;
      }
    }
    int res = sections[sectionNum]->setParameter(&argv[2], argc - 2, param);
    if (res != -1)
      param.addObject(5, this);
    return res;
  }

  if (strcmp(argv[0], "allSections") == 0) {
    if (argc < 2)
      return -1;
    int result = -1;
    for (int i = 0; i < numSections; i++) {
      int ok = sections[i]->setParameter(&argv[1], argc - 1, param);
      if (ok != -1)
        result = ok;
    }
    if (result != -1)
      param.addObject(5, this);
    return result;
  }

  if (strstr(argv[0], "section") != 0) {
    if (argc < 3)
      return -1;
    int sectionNum = atoi(argv[1]);
    if (sectionNum > 0 && sectionNum <= numSections) {
      int res = sections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);
      if (res != -1)
        param.addObject(5, this);
      return res;
    }
    return -1;
  }

  if (strstr(argv[0], "integration") != 0) {
    if (argc < 2)
      return -1;
    int res = beamIntegr->setParameter(&argv[1], argc - 1, param);
    if (res != -1)
      param.addObject(5, this);
    return res;
  }

  int result = -1;
  for (int i = 0; i < numSections; i++) {
    int ok = sections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  int ok = beamIntegr->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  if (hingeOn) {
    UniaxialMaterial *uniList[4] = {theHingeZ, theHingeY, theHingeZJ, theHingeYJ};
    for (int k = 0; k < 4; k++) {
      if (uniList[k] != 0) {
        ok = uniList[k]->setParameter(argv, argc, param);
        if (ok != -1) result = ok;
      }
    }
    NDMaterial *ndList[2] = {theHingeC, theHingeCJ};
    for (int k = 0; k < 2; k++) {
      if (ndList[k] != 0) {
        ok = ndList[k]->setParameter(argv, argc, param);
        if (ok != -1) result = ok;
      }
    }
  }

  if (result != -1)
    param.addObject(5, this);
  return result;
}

int
LadrunoForceBeamColumn3d::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  if (parameterID == 2 || parameterID == 3 || parameterID == 4 || parameterID == 5) {
    if (parameterID == 2) {
      if (info.theDouble <= 0.0)
        return -1;
      lchMode = 2;
      userLch = info.theDouble;
      current_section_lch = userLch;
    }
    if (parameterID == 3)
      GAvy = (info.theDouble > 0.0) ? info.theDouble : 0.0;
    if (parameterID == 4)
      GAvz = (info.theDouble > 0.0) ? info.theDouble : 0.0;

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
LadrunoForceBeamColumn3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}

const Matrix &
LadrunoForceBeamColumn3d::getKiSensitivity(int gradNumber)
{
  theMatrix.Zero();
  return theMatrix;
}

const Matrix &
LadrunoForceBeamColumn3d::getMassSensitivity(int gradNumber)
{
  theMatrix.Zero();
  if (rho == 0.0 || parameterID != 1)
    return theMatrix;

  double L = crdTransf->getInitialLength();
  if (cMass == 0) {
    theMatrix(0, 0) = theMatrix(1, 1) = theMatrix(2, 2) =
      theMatrix(6, 6) = theMatrix(7, 7) = theMatrix(8, 8) = 0.5 * L;
  } else {
    static thread_local Matrix ml(12, 12);
    ml.Zero();
    double m = L / 420.0;
    ml(0, 0) = ml(6, 6) = m * 140.0;
    ml(0, 6) = ml(6, 0) = m * 70.0;

    ml(2, 2) = ml(8, 8) = m * 156.0;
    ml(2, 8) = ml(8, 2) = m * 54.0;
    ml(4, 4) = ml(10, 10) = m * 4.0 * L * L;
    ml(4, 10) = ml(10, 4) = -m * 3.0 * L * L;
    ml(2, 4) = ml(4, 2) = -m * 22.0 * L;
    ml(8, 10) = ml(10, 8) = -ml(2, 4);
    ml(2, 10) = ml(10, 2) = m * 13.0 * L;
    ml(4, 8) = ml(8, 4) = -ml(2, 10);

    ml(1, 1) = ml(7, 7) = m * 156.0;
    ml(1, 7) = ml(7, 1) = m * 54.0;
    ml(5, 5) = ml(11, 11) = m * 4.0 * L * L;
    ml(5, 11) = ml(11, 5) = -m * 3.0 * L * L;
    ml(1, 5) = ml(5, 1) = m * 22.0 * L;
    ml(7, 11) = ml(11, 7) = -ml(1, 5);
    ml(1, 11) = ml(11, 1) = -m * 13.0 * L;
    ml(5, 7) = ml(7, 5) = -ml(1, 11);

    theMatrix = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  return theMatrix;
}

const Vector &
LadrunoForceBeamColumn3d::getResistingForceSensitivity(int gradNumber)
{
  static thread_local Vector dqdh(6);
  dqdh = this->computedqdh(gradNumber);

  double dp0dh[6];
  for (int k = 0; k < 6; k++)
    dp0dh[k] = 0.0;
  if (numEleLoads > 0)
    this->computeReactionSensitivity(dp0dh, gradNumber);
  Vector dp0dhVec(dp0dh, 6);

  static thread_local Vector P(12);
  P.Zero();

  if (crdTransf->isShapeSensitivity()) {
    double p0[6];
    for (int k = 0; k < 6; k++)
      p0[k] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);
    Vector p0Vec(p0, 6);
    P += crdTransf->getGlobalResistingForceShapeSensitivity(Se, p0Vec, gradNumber);
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);
  }

  P += crdTransf->getGlobalResistingForce(dqdh, dp0dhVec);
  return P;
}

int
LadrunoForceBeamColumn3d::commitSensitivity(int gradNumber, int numGrads)
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

  static thread_local Vector dvdh(6);
  dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  static thread_local Vector dqdh(6);
  dqdh = this->computedqdh(gradNumber);
  dqdh.addMatrixVector(1.0, kv, dvdh, 1.0);

  if (crdTransf->isShapeSensitivity()) {
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);
  }

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

    int j;
    for (j = 0; j < order; j++) {
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
      case SECTION_RESPONSE_MY:
        ds(j) += xL1 * dqdh(3) + xL * dqdh(4);
        break;
      case SECTION_RESPONSE_VZ:
        ds(j) += oneOverL * (dqdh(3) + dqdh(4));
        break;
      case SECTION_RESPONSE_T:
        ds(j) += dqdh(5);
        break;
      default:
        ds(j) += 0.0;
        break;
      }
    }

    const Vector &dsdh = sections[i]->getStressResultantSensitivity(gradNumber, true);
    ds -= dsdh;

    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        ds(j) += dxLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_VY:
        ds(j) += d1oLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_MY:
        ds(j) += dxLdh * (Se(3) + Se(4));
        break;
      case SECTION_RESPONSE_VZ:
        ds(j) += d1oLdh * (Se(3) + Se(4));
        break;
      default:
        break;
      }
    }

    Vector de(&workArea[order], order);
    const Matrix &fsSec = sections[i]->getSectionFlexibility();
    de.addMatrixVector(0.0, fsSec, ds, 1.0);

    err += sections[i]->commitSensitivity(de, gradNumber, numGrads);
  }

  return err;
}

const Vector &
LadrunoForceBeamColumn3d::computedqdh(int gradNumber)
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

  static thread_local Vector dvdh(6);
  dvdh.Zero();

  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();

    double xL = pts[i];
    double xL1 = xL - 1.0;
    double wtL = wts[i] * L;

    double dxLdh = dptsdh[i];
    double dwtLdh = wts[i] * dLdh + dwtsdh[i] * L;

    Vector ds(workArea, order);
    ds.Zero();

    if (numEleLoads > 0)
      this->computeSectionForceSensitivity(ds, i, gradNumber);

    const Vector &dsdh = sections[i]->getStressResultantSensitivity(gradNumber, true);
    ds -= dsdh;

    int j;
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
        ds(j) += dxLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_VY:
        ds(j) += d1oLdh * (Se(1) + Se(2));
        break;
      case SECTION_RESPONSE_MY:
        ds(j) += dxLdh * (Se(3) + Se(4));
        break;
      case SECTION_RESPONSE_VZ:
        ds(j) += d1oLdh * (Se(3) + Se(4));
        break;
      default:
        break;
      }
    }

    Vector de(&workArea[order], order);
    const Matrix &fsSec = sections[i]->getSectionFlexibility();
    de.addMatrixVector(0.0, fsSec, ds, 1.0);

    for (j = 0; j < order; j++) {
      double dei = -de(j) * wtL;
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
      case SECTION_RESPONSE_MY:
        dvdh(3) += xL1 * dei;
        dvdh(4) += xL * dei;
        break;
      case SECTION_RESPONSE_VZ:
        dei = oneOverL * dei;
        dvdh(3) += dei;
        dvdh(4) += dei;
        break;
      case SECTION_RESPONSE_T:
        dvdh(5) += dei;
        break;
      default:
        break;
      }
    }

    const Vector &e = vs[i];
    for (j = 0; j < order; j++) {
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
      case SECTION_RESPONSE_MY:
        dvdh(3) -= xL1 * e(j) * dwtLdh;
        dvdh(4) -= xL * e(j) * dwtLdh;
        dvdh(3) -= dxLdh * e(j) * wtL;
        dvdh(4) -= dxLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_VZ:
        dvdh(3) -= oneOverL * e(j) * dwtLdh;
        dvdh(4) -= oneOverL * e(j) * dwtLdh;
        dvdh(3) -= d1oLdh * e(j) * wtL;
        dvdh(4) -= d1oLdh * e(j) * wtL;
        break;
      case SECTION_RESPONSE_T:
        dvdh(5) -= e(j) * dwtLdh;
        break;
      default:
        break;
      }
    }
  }

  static thread_local Vector dqdh(6);
  dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
  return dqdh;
}

const Matrix &
LadrunoForceBeamColumn3d::computedfedh(int gradNumber)
{
  static thread_local Matrix dfedh(6, 6);
  dfedh.Zero();
  return dfedh;
}

void
LadrunoForceBeamColumn3d::setSectionPointers(int numSec, SectionForceDeformation **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: LadrunoForceBeamColumn3d::setSectionPointers -- max number of sections exceeded\n";
  }

  numSections = numSec;

  if (secPtrs == 0) {
    opserr << "Error: LadrunoForceBeamColumn3d::setSectionPointers -- invalid section pointer\n";
  }

  sections = new SectionForceDeformation *[numSections];
  for (int i = 0; i < numSections; i++) {
    if (secPtrs[i] == 0) {
      opserr << "Error: LadrunoForceBeamColumn3d::setSectionPointers -- null section pointer " << i << endln;
    }
    sections[i] = secPtrs[i]->getCopy();
    if (sections[i] == 0) {
      opserr << "Error: LadrunoForceBeamColumn3d::setSectionPointers -- could not create copy of section " << i << endln;
    }
    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    for (int j = 0; j < order; j++) {
      if (code(j) == SECTION_RESPONSE_T)
        isTorsion = true;
    }
  }

  fs = new Matrix[numSections];
  vs = new Vector[numSections];
  Ssr = new Vector[numSections];
  vscommit = new Vector[numSections];
  fscommit = new Matrix[numSections];
  Ssrcommit = new Vector[numSections];
}

bool
LadrunoForceBeamColumn3d::ladrunoThreadSafeUpdate(void) const
{
  if (theDamping != 0 || crdTransf == 0 || beamIntegr == 0 || sections == 0)
    return false;
  int ctTag = crdTransf->getClassTag();
  if (ctTag != CRDTR_TAG_LinearCrdTransf3d &&
      ctTag != CRDTR_TAG_PDeltaCrdTransf3d &&
      ctTag != CRDTR_TAG_CorotCrdTransf3d)
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
  for (int i = 0; i < numSections; i++)
    if (sections[i] == 0 || !sections[i]->ladrunoThreadSafeUpdate())
      return false;
  if (hingeOn) {
    if (theHingeZ  != 0 && !theHingeZ->ladrunoThreadSafeUpdate())  return false;
    if (theHingeY  != 0 && !theHingeY->ladrunoThreadSafeUpdate())  return false;
    if (theHingeC  != 0 && !theHingeC->ladrunoThreadSafeUpdate())  return false;
    if (theHingeZJ != 0 && !theHingeZJ->ladrunoThreadSafeUpdate()) return false;
    if (theHingeYJ != 0 && !theHingeYJ->ladrunoThreadSafeUpdate()) return false;
    if (theHingeCJ != 0 && !theHingeCJ->ladrunoThreadSafeUpdate()) return false;
  }
  return true;
}

bool
LadrunoForceBeamColumn3d::ladrunoThreadSafeForm(void) const
{
  if (!this->ladrunoThreadSafeUpdate())
    return false;
  if (numEleLoads != 0 || eleLoadDirty)
    return false;
  for (int i = 0; i < numSections; i++)
    if (!sections[i]->ladrunoThreadSafeForm())
      return false;
  if (hingeOn) {
    if (theHingeZ  != 0 && !theHingeZ->ladrunoThreadSafeForm())  return false;
    if (theHingeY  != 0 && !theHingeY->ladrunoThreadSafeForm())  return false;
    if (theHingeC  != 0 && !theHingeC->ladrunoThreadSafeForm())  return false;
    if (theHingeZJ != 0 && !theHingeZJ->ladrunoThreadSafeForm()) return false;
    if (theHingeYJ != 0 && !theHingeYJ->ladrunoThreadSafeForm()) return false;
    if (theHingeCJ != 0 && !theHingeCJ->ladrunoThreadSafeForm()) return false;
  }
  return true;
}
