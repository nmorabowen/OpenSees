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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 10/2026
//
// Parser/factory for LadrunoCPML — unified 2D/3D Convolutional Perfectly
// Matched Layer (CFS-PML / C-PML, Matzen 2011 / Xie et al. 2014) wave
// absorbing element (ELE 33023).
//
//   element LadrunoCPML $tag $n1 .. $nk
//       (-E $E -nu $nu -rho $rho | -mat $matTag <-rho $rho>)
//       (-box $xmin $xmax $ymin $ymax <$zmin $zmax> <$Lx <$Ly <$Lz>>>
//        | -plane $x0 $y0 <$z0> $nx $ny <$nz> <$Lpml>
//        | -plane -x0 $x0 $y0 <$z0> -n $nx $ny <$nz> -L $Lpml)
//       <-L $L1 <$L2 <$L3>>> <-Lx $Lx> <-Ly $Ly> <-Lz $Lz>
//       <-m $mDeg> <-R0 $R0> <-kappaMax $k> <-alphaMax $a | -f0 $f0>
//       <-thick $t> <-lumped> <-stage 0|1>
//       <-shape T3|Q4|T6|BT6|bezier|lagrange|TET4|H8|TET10|BTET10|H20>
//       <-body $bx $by <$bz>>
//       <-up -Kf $Kf -poro $n -rhoF $rhoF
//        (-perm $kx <$ky <$kz>> | -permH $khx <$khy <$khz>> -gammaW $gw)
//        <-alphaB $alphaB> <-Ks $Ks> <-stab off|auto <$a0>|manual $a|$a>
//        <-pOrder equal|linear> <-formulation std|bbar>
//        <-fluidBody $fbx $fby <$fbz>>>

#include <LadrunoCPML.h>
#include <NDMaterial.h>
#include <Domain.h>
#include <Node.h>
#include <elementAPI.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>

static const char *cpmlGetTok(char *buf, int len)
{
  const char *s = OPS_GetStringFromAll(buf, len);
  if (s == 0)
    return 0;
  if (s != buf) {
    strncpy(buf, s, (size_t)len - 1);
    buf[len - 1] = '\0';
  }
  return buf;
}

static bool cpmlTryReadDouble(double &val)
{
  if (OPS_GetNumRemainingInputArgs() <= 0)
    return false;
  char tok[64];
  if (cpmlGetTok(tok, sizeof(tok)) == 0)
    return false;
  if (tok[0] == '-' && (tok[1] < '0' || tok[1] > '9') && tok[1] != '.') {
    OPS_ResetCurrentInputArg(-1);
    return false;
  }
  char *endp = 0;
  double v = strtod(tok, &endp);
  if (endp == tok || *endp != '\0' || !std::isfinite(v)) {
    OPS_ResetCurrentInputArg(-1);
    return false;
  }
  val = v;
  return true;
}

void *OPS_LadrunoCPML()
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "ERROR LadrunoCPML -- insufficient arguments\n"
              "Want: element LadrunoCPML tag n1 .. nk (-E E -nu nu -rho rho | -mat matTag) "
              "(-box ... | -plane ...) <options>\n";
    return 0;
  }

  // ---- Greedy read of leading integers: tag + n1..nk (up to 20 nodes) ------
  const int maxLead = 22;   // 1 tag + up to 20 node tags + 1 slack
  int lead[maxLead];
  int nLead = 0;
  while (OPS_GetNumRemainingInputArgs() > 0 && nLead < maxLead) {
    char tok[64];
    if (cpmlGetTok(tok, sizeof(tok)) == 0)
      break;
    if (tok[0] == '-' && (tok[1] < '0' || tok[1] > '9')) {
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    char *endp = 0;
    long v = strtol(tok, &endp, 10);
    if (endp == tok || *endp != '\0') {
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    lead[nLead++] = (int)v;
  }

  if (nLead < 4) {
    opserr << "ERROR LadrunoCPML -- need tag + >=3 node tags before flags "
              "(got " << nLead << " leading integers)\n";
    return 0;
  }
  if (nLead >= maxLead) {
    opserr << "ERROR LadrunoCPML -- too many leading integer arguments "
              "(max 20 node tags supported)\n";
    return 0;
  }

  int tag = lead[0];
  int nNodes = nLead - 1;
  ID nodeTags(nNodes);
  for (int i = 0; i < nNodes; i++)
    nodeTags(i) = lead[1 + i];

  // ---- Determine ndm from domain nodes first, fallback to OPS_GetNDM() -----
  int ndm = 0;
  Domain *dom = OPS_GetDomain();
  if (dom != 0) {
    Node *n0 = dom->getNode(nodeTags(0));
    if (n0 != 0) {
      ndm = n0->getCrds().Size();
    }
  }
  if (ndm != 2 && ndm != 3) {
    if (nNodes == 3 || nNodes == 6)
      ndm = 2;
    else if (nNodes == 8 || nNodes == 10 || nNodes == 20)
      ndm = 3;
    else
      ndm = OPS_GetNDM();
  }

  if (ndm != 2 && ndm != 3) {
    opserr << "ERROR LadrunoCPML " << tag << " -- model ndm must be 2 or 3, got "
           << ndm << "\n";
    return 0;
  }

  // ---- Validate (ndm, nNodes) legality -------------------------------------
  bool validTopo = (ndm == 2 && (nNodes == 3 || nNodes == 4 || nNodes == 6)) ||
                   (ndm == 3 && (nNodes == 4 || nNodes == 8 || nNodes == 10 || nNodes == 20));
  if (!validTopo) {
    opserr << "ERROR LadrunoCPML " << tag << " -- invalid topology for ndm="
           << ndm << " with " << nNodes
           << " nodes (want 2D: 3/4/6 nodes; 3D: 4/8/10/20 nodes)\n";
    return 0;
  }

  // ---- Option defaults -----------------------------------------------------
  int shapeOverride = 0;    // 0=auto, 1=T6, 2=BT6, 3=TET10, 4=BTET10
  double E = 0.0;           bool EGiven = false;
  double nu = 0.0;          bool nuGiven = false;
  double rho = 0.0;         bool rhoGiven = false;
  bool matGiven = false;
  double matE = 0.0;
  double matNu = 0.0;
  double matRho = 0.0;

  double thickness = 1.0;   bool thickGiven = false;
  int geomMode = -1;        // 0 = box, 1 = plane
  bool boxGiven = false;
  bool planeGiven = false;

  double boxBounds[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Lx = 0.0, Ly = 0.0, Lz = 0.0;
  bool LxGiven = false, LyGiven = false, LzGiven = false;
  bool hasLx = false, hasLy = false, hasLz = false, hasMultiL = false;

  double planeX0[3] = {0.0, 0.0, 0.0};
  double planeN[3]  = {0.0, 0.0, 0.0};
  double Lpml = 0.0;
  bool x0Given = false, nGiven = false, LpmlGiven = false;
  bool hasX0 = false, hasN = false;

  double mDeg = 2.0;
  double R0 = 1.0e-4;
  double kappaMax = 5.0;
  double alphaMax = 31.41592653589793; // 10 * pi
  bool alphaGiven = false;
  bool f0Given = false;
  bool lumped = false;
  int initialStage = 1;
  Vector bodyForces(ndm);

  // Coupled Biot u-p mode options
  bool upMode = false;
  double Kf = 0.0;          bool KfGiven = false;
  double poro = 0.0;        bool poroGiven = false;
  double rhoF = 0.0;        bool rhoFGiven = false;
  double perm[3]  = {0.0, 0.0, 0.0}; bool permGiven = false;
  double permH[3] = {0.0, 0.0, 0.0}; bool permHGiven = false;
  double gammaW = 0.0;      bool gammaWGiven = false;
  double biotAlpha = 1.0;   bool alphaBGiven = false;
  double Ks = 0.0;          bool KsGiven = false;
  int stabMode = 1;         // 0 = off, 1 = auto, 2 = manual
  double stabValue = 0.25;  bool stabGiven = false;
  int pOrder = 0;           bool pOrderGiven = false;
  int formulation = 0;      bool formGiven = false;
  Vector fluidBodyForces(ndm);
  bool fluidBodyGiven = false;

  int num = 1;

  // ---- Flag loop (Fail-closed: unknown flags are fatal) --------------------
  while (OPS_GetNumRemainingInputArgs() > 0) {
    char opt[64];
    if (cpmlGetTok(opt, sizeof(opt)) == 0)
      break;

    if (strcmp(opt, "-E") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &E) < 0 || E <= 0.0 || !std::isfinite(E)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -E (want E > 0)\n";
        return 0;
      }
      EGiven = true;

    } else if (strcmp(opt, "-nu") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &nu) < 0 || nu <= -1.0 || nu >= 0.499999 || !std::isfinite(nu)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad -nu (want -1.0 < nu < 0.499999)\n";
        return 0;
      }
      nuGiven = true;

    } else if (strcmp(opt, "-rho") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rho) < 0 || rho <= 0.0 || !std::isfinite(rho)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -rho (want rho > 0)\n";
        return 0;
      }
      rhoGiven = true;

    } else if (strcmp(opt, "-mat") == 0) {
      int matTag = 0;
      num = 1;
      if (OPS_GetIntInput(&num, &matTag) < 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -mat tag\n";
        return 0;
      }
      NDMaterial *matBase = OPS_getNDMaterial(matTag);
      if (matBase == nullptr) {
        opserr << "WARNING element LadrunoCPML (eleTag " << tag
               << "): NDMaterial with tag " << matTag << " not found.\n";
        return nullptr;
      }
      const char *copyType = (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
      NDMaterial *mat = matBase->getCopy(copyType);
      if (mat == nullptr) {
        mat = matBase->getCopy();
      }
      if (mat == nullptr) {
        opserr << "WARNING element LadrunoCPML (eleTag " << tag
               << "): failed to obtain " << copyType << " copy of NDMaterial " << matTag << ".\n";
        return nullptr;
      }
      const char *mtype = mat->getType();
      if (mtype == nullptr ||
          (strstr(mtype, "ElasticIsotropic") == nullptr &&
           strstr(mtype, "PlaneStrain") == nullptr &&
           strstr(mtype, "ThreeDimensional") == nullptr)) {
        opserr << "WARNING element LadrunoCPML (eleTag " << tag
               << "): NDMaterial " << matTag << " has type '"
               << (mtype ? mtype : "null")
               << "'; only ElasticIsotropic (PlaneStrain / ThreeDimensional) is supported.\n";
        delete mat;
        return nullptr;
      }
      double rhoFromMat = mat->getRho();
      if (rhoFromMat <= 0.0) {
        rhoFromMat = matBase->getRho();
      }
      const Matrix &D = mat->getInitialTangent();
      if (D.noRows() < 3 || D.noCols() < 3) {
        opserr << "WARNING element LadrunoCPML (eleTag " << tag
               << "): NDMaterial " << matTag
               << " initial tangent is too small (" << D.noRows() << "x" << D.noCols() << ").\n";
        delete mat;
        return nullptr;
      }
      double c11 = D(0, 0);
      double c12 = D(0, 1);
      delete mat;
      mat = nullptr;
      double denom = c11 + c12;
      if (!std::isfinite(c11) || !std::isfinite(c12) || std::fabs(denom) < 1.0e-15) {
        opserr << "WARNING element LadrunoCPML (eleTag " << tag
               << "): failed to recover valid (E, nu) from NDMaterial " << matTag << ".\n";
        return nullptr;
      }
      matNu = c12 / denom;
      matE  = (c11 - c12) * (1.0 + matNu);
      if (!std::isfinite(matE) || !std::isfinite(matNu) ||
          matE <= 0.0 || matNu <= -1.0 || matNu >= 0.499999) {
        opserr << "WARNING element LadrunoCPML (eleTag " << tag
               << "): NDMaterial " << matTag
               << " extracted (E=" << matE << ", nu=" << matNu
               << ") out of isotropic elastic bounds.\n";
        return nullptr;
      }
      matRho = rhoFromMat;
      matGiven = true;

    } else if (strcmp(opt, "-box") == 0) {
      if (planeGiven) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- cannot specify both -box and -plane\n";
        return 0;
      }
      num = 2 * ndm;
      if (OPS_GetDoubleInput(&num, boxBounds) < 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- -box requires " << (2 * ndm)
               << " coordinate bounds\n";
        return 0;
      }
      for (int d = 0; d < 2 * ndm; d++) {
        if (!std::isfinite(boxBounds[d])) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -box bounds must be finite numbers\n";
          return 0;
        }
      }
      for (int d = 0; d < ndm; d++) {
        if (boxBounds[2 * d + 1] < boxBounds[2 * d]) {
          opserr << "ERROR LadrunoCPML " << tag << " -- -box bound max < min on axis "
                 << d << " (" << boxBounds[2 * d + 1] << " < " << boxBounds[2 * d] << ")\n";
          return 0;
        }
      }
      // Optional inline Lx <Ly <Lz>> right after box bounds
      double l0 = 0.0;
      if (cpmlTryReadDouble(l0)) {
        if (l0 <= 0.0 || !std::isfinite(l0)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- inline -box PML thickness must be > 0\n";
          return 0;
        }
        double l1 = 0.0;
        if (cpmlTryReadDouble(l1)) {
          if (l1 <= 0.0 || !std::isfinite(l1)) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- inline -box PML thickness Ly must be > 0\n";
            return 0;
          }
          Lx = l0; LxGiven = true;
          Ly = l1; LyGiven = true;
          if (ndm == 3) {
            double l2 = 0.0;
            if (!cpmlTryReadDouble(l2) || l2 <= 0.0 || !std::isfinite(l2)) {
              opserr << "ERROR LadrunoCPML " << tag
                     << " -- inline 3D -box PML thickness requires 1 or 3 positive values\n";
              return 0;
            }
            Lz = l2; LzGiven = true;
          }
        } else {
          Lx = l0; LxGiven = true;
          Ly = l0; LyGiven = true;
          Lz = l0; LzGiven = true;
        }
      }
      boxGiven = true;
      geomMode = 0;

    } else if (strcmp(opt, "-plane") == 0) {
      if (boxGiven) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- cannot specify both -box and -plane\n";
        return 0;
      }
      // Check if inline positional numbers follow -plane:
      // -plane $x0 $y0 <$z0> $nx $ny <$nz> <$Lpml>
      double firstVal = 0.0;
      if (cpmlTryReadDouble(firstVal)) {
        planeX0[0] = firstVal;
        for (int d = 1; d < ndm; d++) {
          if (!cpmlTryReadDouble(planeX0[d]) || !std::isfinite(planeX0[d])) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- -plane inline form needs " << ndm << " x0 coords and "
                   << ndm << " normal components\n";
            return 0;
          }
        }
        for (int d = 0; d < ndm; d++) {
          if (!cpmlTryReadDouble(planeN[d]) || !std::isfinite(planeN[d])) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- -plane inline form needs " << ndm << " normal components\n";
            return 0;
          }
        }
        x0Given = true;
        nGiven = true;
        double lVal = 0.0;
        if (cpmlTryReadDouble(lVal)) {
          if (lVal <= 0.0 || !std::isfinite(lVal)) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- -plane inline Lpml must be > 0\n";
            return 0;
          }
          Lpml = lVal;
          LpmlGiven = true;
        }
      }
      planeGiven = true;
      geomMode = 1;

    } else if (strcmp(opt, "-x0") == 0) {
      num = ndm;
      if (OPS_GetDoubleInput(&num, planeX0) < 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- -x0 requires " << ndm
               << " coordinates\n";
        return 0;
      }
      for (int d = 0; d < ndm; d++) {
        if (!std::isfinite(planeX0[d])) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -x0 coordinates must be finite\n";
          return 0;
        }
      }
      x0Given = true;
      hasX0 = true;

    } else if (strcmp(opt, "-n") == 0 || strcmp(opt, "-normal") == 0) {
      num = ndm;
      if (OPS_GetDoubleInput(&num, planeN) < 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- " << opt << " requires "
               << ndm << " normal components\n";
        return 0;
      }
      for (int d = 0; d < ndm; d++) {
        if (!std::isfinite(planeN[d])) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- normal components must be finite\n";
          return 0;
        }
      }
      nGiven = true;
      hasN = true;

    } else if (strcmp(opt, "-L") == 0) {
      double l0 = 0.0;
      num = 1;
      if (OPS_GetDoubleInput(&num, &l0) < 0 || l0 <= 0.0 || !std::isfinite(l0)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -L (want L > 0)\n";
        return 0;
      }
      double l1 = 0.0;
      if (cpmlTryReadDouble(l1)) {
        if (l1 <= 0.0 || !std::isfinite(l1)) {
          opserr << "ERROR LadrunoCPML " << tag << " -- bad -L Ly (want Ly > 0)\n";
          return 0;
        }
        Lx = l0; LxGiven = true;
        Ly = l1; LyGiven = true;
        hasMultiL = true;
        if (ndm == 3) {
          double l2 = 0.0;
          if (!cpmlTryReadDouble(l2) || l2 <= 0.0 || !std::isfinite(l2)) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- 3D -L requires 1 or 3 positive values\n";
            return 0;
          }
          Lz = l2; LzGiven = true;
        }
      } else {
        Lx = l0; LxGiven = true;
        Ly = l0; LyGiven = true;
        Lz = l0; LzGiven = true;
        Lpml = l0; LpmlGiven = true;
      }

    } else if (strcmp(opt, "-Lx") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Lx) < 0 || Lx <= 0.0 || !std::isfinite(Lx)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -Lx (want Lx > 0)\n";
        return 0;
      }
      LxGiven = true;
      hasLx = true;

    } else if (strcmp(opt, "-Ly") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Ly) < 0 || Ly <= 0.0 || !std::isfinite(Ly)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -Ly (want Ly > 0)\n";
        return 0;
      }
      LyGiven = true;
      hasLy = true;

    } else if (strcmp(opt, "-Lz") == 0) {
      if (ndm == 2) {
        opserr << "ERROR LadrunoCPML " << tag << " -- -Lz is 3D only\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &Lz) < 0 || Lz <= 0.0 || !std::isfinite(Lz)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -Lz (want Lz > 0)\n";
        return 0;
      }
      LzGiven = true;
      hasLz = true;

    } else if (strcmp(opt, "-m") == 0 || strcmp(opt, "-mDeg") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &mDeg) < 0 || mDeg < 0.0 || !std::isfinite(mDeg)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad " << opt
               << " (want m >= 0)\n";
        return 0;
      }

    } else if (strcmp(opt, "-R0") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &R0) < 0 || R0 <= 0.0 || R0 >= 1.0 || !std::isfinite(R0)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad -R0 (want 0 < R0 < 1)\n";
        return 0;
      }

    } else if (strcmp(opt, "-kappaMax") == 0 || strcmp(opt, "-kappa") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &kappaMax) < 0 || kappaMax < 1.0 || !std::isfinite(kappaMax)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad " << opt
               << " (want kappaMax >= 1.0)\n";
        return 0;
      }

    } else if (strcmp(opt, "-alphaMax") == 0 || strcmp(opt, "-alpha") == 0) {
      if (f0Given) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- specify either -alphaMax or -f0, not both\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &alphaMax) < 0 || alphaMax < 0.0 || !std::isfinite(alphaMax)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad " << opt
               << " (want alphaMax >= 0.0)\n";
        return 0;
      }
      alphaGiven = true;

    } else if (strcmp(opt, "-f0") == 0) {
      if (alphaGiven) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- specify either -alphaMax or -f0, not both\n";
        return 0;
      }
      double f0 = 0.0;
      num = 1;
      if (OPS_GetDoubleInput(&num, &f0) < 0 || f0 < 0.0 || !std::isfinite(f0)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -f0 (want f0 >= 0.0)\n";
        return 0;
      }
      alphaMax = 3.14159265358979323846 * f0;
      f0Given = true;

    } else if (strcmp(opt, "-thick") == 0 || strcmp(opt, "-thickness") == 0) {
      if (ndm == 3) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- -thick is 2D only (plane strain); remove it in 3D\n";
        return 0;
      }
      num = 1;
      if (OPS_GetDoubleInput(&num, &thickness) < 0 || thickness <= 0.0 || !std::isfinite(thickness)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -thick (want t > 0)\n";
        return 0;
      }
      thickGiven = true;

    } else if (strcmp(opt, "-lumped") == 0) {
      lumped = true;

    } else if (strcmp(opt, "-stage") == 0) {
      num = 1;
      if (OPS_GetIntInput(&num, &initialStage) < 0 ||
          (initialStage != 0 && initialStage != 1)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad -stage (want 0 or 1)\n";
        return 0;
      }

    } else if (strcmp(opt, "-shape") == 0) {
      char sh[32];
      if (cpmlGetTok(sh, sizeof(sh)) == 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- -shape needs a shape name\n";
        return 0;
      }
      if (strcmp(sh, "T3") == 0 || strcmp(sh, "t3") == 0) {
        if (!(ndm == 2 && nNodes == 3)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape T3 requires ndm=2 and 3 nodes\n";
          return 0;
        }
        shapeOverride = 0;
      } else if (strcmp(sh, "Q4") == 0 || strcmp(sh, "q4") == 0) {
        if (!(ndm == 2 && nNodes == 4)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape Q4 requires ndm=2 and 4 nodes\n";
          return 0;
        }
        shapeOverride = 0;
      } else if (strcmp(sh, "T6") == 0 || strcmp(sh, "t6") == 0) {
        if (!(ndm == 2 && nNodes == 6)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape T6 requires ndm=2 and 6 nodes\n";
          return 0;
        }
        shapeOverride = 1;
      } else if (strcmp(sh, "BT6") == 0 || strcmp(sh, "bt6") == 0 ||
                 strcmp(sh, "BezierTri6") == 0) {
        if (!(ndm == 2 && nNodes == 6)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape BT6 requires ndm=2 and 6 nodes\n";
          return 0;
        }
        shapeOverride = 2;
      } else if (strcmp(sh, "TET4") == 0 || strcmp(sh, "tet4") == 0 ||
                 strcmp(sh, "Tet4") == 0) {
        if (!(ndm == 3 && nNodes == 4)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape TET4 requires ndm=3 and 4 nodes\n";
          return 0;
        }
        shapeOverride = 0;
      } else if (strcmp(sh, "H8") == 0 || strcmp(sh, "h8") == 0) {
        if (!(ndm == 3 && nNodes == 8)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape H8 requires ndm=3 and 8 nodes\n";
          return 0;
        }
        shapeOverride = 0;
      } else if (strcmp(sh, "TET10") == 0 || strcmp(sh, "tet10") == 0 ||
                 strcmp(sh, "Tet10") == 0) {
        if (!(ndm == 3 && nNodes == 10)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape TET10 requires ndm=3 and 10 nodes\n";
          return 0;
        }
        shapeOverride = 3;
      } else if (strcmp(sh, "BTET10") == 0 || strcmp(sh, "btet10") == 0 ||
                 strcmp(sh, "BezierTet10") == 0) {
        if (!(ndm == 3 && nNodes == 10)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape BTET10 requires ndm=3 and 10 nodes\n";
          return 0;
        }
        shapeOverride = 4;
      } else if (strcmp(sh, "H20") == 0 || strcmp(sh, "h20") == 0) {
        if (!(ndm == 3 && nNodes == 20)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape H20 requires ndm=3 and 20 nodes\n";
          return 0;
        }
        shapeOverride = 0;
      } else if (strcmp(sh, "bezier") == 0 || strcmp(sh, "Bezier") == 0) {
        if (ndm == 2 && nNodes == 6)
          shapeOverride = 2;
        else if (ndm == 3 && nNodes == 10)
          shapeOverride = 4;
        else {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -shape bezier requires (ndm=2, 6 nodes) or (ndm=3, 10 nodes)\n";
          return 0;
        }
      } else if (strcmp(sh, "lagrange") == 0 || strcmp(sh, "Lagrange") == 0) {
        if (ndm == 2 && nNodes == 6)
          shapeOverride = 1;
        else if (ndm == 3 && nNodes == 10)
          shapeOverride = 3;
        else
          shapeOverride = 0;
      } else {
        opserr << "ERROR LadrunoCPML " << tag << " -- unknown -shape '" << sh
               << "' (want T3|Q4|T6|BT6|bezier|lagrange|TET4|H8|TET10|BTET10|H20)\n";
        return 0;
      }

    } else if (strcmp(opt, "-body") == 0 || strcmp(opt, "-b") == 0) {
      double bb[3] = {0.0, 0.0, 0.0};
      num = ndm;
      if (OPS_GetDoubleInput(&num, bb) < 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- " << opt << " requires "
               << ndm << " values\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++) {
        if (!std::isfinite(bb[i])) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- body force values must be finite\n";
          return 0;
        }
        bodyForces(i) = bb[i];
      }

    } else if (strcmp(opt, "-up") == 0) {
      upMode = true;

    } else if (strcmp(opt, "-Kf") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Kf) < 0 || Kf <= 0.0 || !std::isfinite(Kf)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad -Kf (want Kf > 0)\n";
        return 0;
      }
      KfGiven = true;

    } else if (strcmp(opt, "-poro") == 0 || strcmp(opt, "-n_poro") == 0 ||
               strcmp(opt, "-porosity") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &poro) < 0 || poro < 0.0 || poro >= 1.0 || !std::isfinite(poro)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad " << opt << " (want 0 <= n < 1)\n";
        return 0;
      }
      poroGiven = true;

    } else if (strcmp(opt, "-rhoF") == 0 || strcmp(opt, "-rhof") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &rhoF) < 0 || rhoF < 0.0 || !std::isfinite(rhoF)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad " << opt << " (want rhoF >= 0)\n";
        return 0;
      }
      rhoFGiven = true;

    } else if (strcmp(opt, "-perm") == 0 || strcmp(opt, "-kbar") == 0) {
      double k0 = 0.0;
      num = 1;
      if (OPS_GetDoubleInput(&num, &k0) < 0 || k0 < 0.0 || !std::isfinite(k0)) {
        opserr << "ERROR LadrunoCPML " << tag << " -- bad " << opt
               << " (want k >= 0)\n";
        return 0;
      }
      double k1 = 0.0;
      if (cpmlTryReadDouble(k1)) {
        if (k1 < 0.0 || !std::isfinite(k1)) {
          opserr << "ERROR LadrunoCPML " << tag << " -- bad " << opt
                 << " ky (want ky >= 0)\n";
          return 0;
        }
        perm[0] = k0;
        perm[1] = k1;
        if (ndm == 3) {
          double k2 = 0.0;
          if (!cpmlTryReadDouble(k2) || k2 < 0.0 || !std::isfinite(k2)) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- 3D " << opt << " requires 1 or 3 non-negative values\n";
            return 0;
          }
          perm[2] = k2;
        }
      } else {
        perm[0] = k0;
        perm[1] = k0;
        perm[2] = k0;
      }
      permGiven = true;

    } else if (strcmp(opt, "-permH") == 0) {
      double kh0 = 0.0;
      num = 1;
      if (OPS_GetDoubleInput(&num, &kh0) < 0 || kh0 < 0.0 || !std::isfinite(kh0)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad -permH (want kh >= 0)\n";
        return 0;
      }
      double kh1 = 0.0;
      if (cpmlTryReadDouble(kh1)) {
        if (kh1 < 0.0 || !std::isfinite(kh1)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- bad -permH khy (want khy >= 0)\n";
          return 0;
        }
        permH[0] = kh0;
        permH[1] = kh1;
        if (ndm == 3) {
          double kh2 = 0.0;
          if (!cpmlTryReadDouble(kh2) || kh2 < 0.0 || !std::isfinite(kh2)) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- 3D -permH requires 1 or 3 non-negative values\n";
            return 0;
          }
          permH[2] = kh2;
        }
      } else {
        permH[0] = kh0;
        permH[1] = kh0;
        permH[2] = kh0;
      }
      permHGiven = true;

    } else if (strcmp(opt, "-gammaW") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &gammaW) < 0 || gammaW <= 0.0 || !std::isfinite(gammaW)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad -gammaW (want gammaW > 0)\n";
        return 0;
      }
      gammaWGiven = true;

    } else if (strcmp(opt, "-alphaB") == 0 || strcmp(opt, "-biotAlpha") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &biotAlpha) < 0 ||
          biotAlpha < 0.0 || biotAlpha > 1.0 || !std::isfinite(biotAlpha)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad " << opt << " (want n <= alphaB <= 1.0)\n";
        return 0;
      }
      alphaBGiven = true;

    } else if (strcmp(opt, "-Ks") == 0) {
      num = 1;
      if (OPS_GetDoubleInput(&num, &Ks) < 0 || Ks <= 0.0 || !std::isfinite(Ks)) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- bad -Ks (want Ks > 0)\n";
        return 0;
      }
      KsGiven = true;

    } else if (strcmp(opt, "-stab") == 0) {
      char stabTok[64];
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          cpmlGetTok(stabTok, sizeof(stabTok)) == 0) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- -stab requires off|none|auto <alpha0>|manual <alpha>|<alpha>\n";
        return 0;
      }
      if (strcmp(stabTok, "off") == 0 || strcmp(stabTok, "none") == 0) {
        stabMode = 0;
        stabValue = 0.0;
        stabGiven = true;
      } else if (strcmp(stabTok, "auto") == 0) {
        double alpha0 = 0.25;
        if (cpmlTryReadDouble(alpha0)) {
          if (alpha0 < 0.0 || !std::isfinite(alpha0)) {
            opserr << "ERROR LadrunoCPML " << tag
                   << " -- -stab auto alpha0 must be >= 0\n";
            return 0;
          }
          stabMode = 1;
          stabValue = alpha0;
        } else {
          stabMode = 1;
          stabValue = 0.25;
        }
        stabGiven = true;
      } else if (strcmp(stabTok, "manual") == 0) {
        double alpha = 0.0;
        num = 1;
        if (OPS_GetDoubleInput(&num, &alpha) < 0 || alpha < 0.0 || !std::isfinite(alpha)) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -stab manual requires numeric alpha >= 0\n";
          return 0;
        }
        stabMode = 2;
        stabValue = alpha;
        stabGiven = true;
      } else {
        char *endp = 0;
        double v = strtod(stabTok, &endp);
        if (endp == stabTok || *endp != '\0' || !std::isfinite(v) || v < 0.0) {
          opserr << "ERROR LadrunoCPML " << tag << " -- bad -stab value '"
                 << stabTok << "' (want off|none|auto|manual|numeric >= 0)\n";
          return 0;
        }
        if (v == 0.0) {
          stabMode = 0;
          stabValue = 0.0;
        } else {
          stabMode = 1;
          stabValue = v;
        }
        stabGiven = true;
      }

    } else if (strcmp(opt, "-pOrder") == 0) {
      char ord[32];
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          cpmlGetTok(ord, sizeof(ord)) == 0) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- -pOrder requires equal|linear\n";
        return 0;
      }
      if (strcmp(ord, "equal") == 0 || strcmp(ord, "0") == 0) {
        pOrder = 0;
        pOrderGiven = true;
      } else if (strcmp(ord, "linear") == 0 || strcmp(ord, "1") == 0 ||
                 strcmp(ord, "TH") == 0 || strcmp(ord, "taylor-hood") == 0) {
        pOrder = 1;
        pOrderGiven = true;
      } else {
        opserr << "ERROR LadrunoCPML " << tag << " -- unknown -pOrder '" << ord
               << "' (want equal|linear)\n";
        return 0;
      }

    } else if (strcmp(opt, "-formulation") == 0 || strcmp(opt, "-form") == 0) {
      char form[32];
      if (OPS_GetNumRemainingInputArgs() < 1 ||
          cpmlGetTok(form, sizeof(form)) == 0) {
        opserr << "ERROR LadrunoCPML " << tag
               << " -- " << opt << " requires std|bbar\n";
        return 0;
      }
      if (strcmp(form, "std") == 0 || strcmp(form, "standard") == 0 ||
          strcmp(form, "0") == 0) {
        formulation = 0;
        formGiven = true;
      } else if (strcmp(form, "bbar") == 0 || strcmp(form, "Bbar") == 0 ||
                 strcmp(form, "1") == 0) {
        formulation = 1;
        formGiven = true;
      } else {
        opserr << "ERROR LadrunoCPML " << tag << " -- unknown " << opt << " '"
               << form << "' (want std|bbar)\n";
        return 0;
      }

    } else if (strcmp(opt, "-fluidBody") == 0) {
      double fb[3] = {0.0, 0.0, 0.0};
      num = ndm;
      if (OPS_GetDoubleInput(&num, fb) < 0) {
        opserr << "ERROR LadrunoCPML " << tag << " -- -fluidBody requires "
               << ndm << " values\n";
        return 0;
      }
      for (int i = 0; i < ndm; i++) {
        if (!std::isfinite(fb[i])) {
          opserr << "ERROR LadrunoCPML " << tag
                 << " -- -fluidBody values must be finite\n";
          return 0;
        }
        fluidBodyForces(i) = fb[i];
      }
      fluidBodyGiven = true;

    } else {
      opserr << "ERROR LadrunoCPML " << tag << " -- unknown flag '" << opt
             << "' (unknown flags are fatal for this element)\n";
      return 0;
    }
  }

  (void)thickGiven;

  // ---- Material resolution and validation ----------------------------------
  if (matGiven) {
    if (EGiven || nuGiven) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- specify either -mat or (-E, -nu), not both\n";
      return 0;
    }
    E = matE;
    nu = matNu;
    EGiven = true;
    nuGiven = true;
    if (!rhoGiven && matRho > 0.0) {
      rho = matRho;
      rhoGiven = true;
    }
  }

  if (!EGiven || !nuGiven || !rhoGiven) {
    opserr << "ERROR LadrunoCPML " << tag << " -- missing required material property:";
    if (!EGiven)   opserr << " -E";
    if (!nuGiven)  opserr << " -nu";
    if (!rhoGiven) opserr << " -rho";
    opserr << " (or -mat with positive density)\n";
    return 0;
  }

  // ---- Geometry mode validation & geomParams packing -----------------------
  if (!boxGiven && !planeGiven) {
    opserr << "ERROR LadrunoCPML " << tag
           << " -- must specify either -box or -plane geometry mode\n";
    return 0;
  }

  Vector geomParams;
  if (geomMode == 0) {
    // -box mode
    if (hasX0 || hasN) {
      opserr << "WARNING element LadrunoCPML (eleTag " << tag
             << "): -x0 and -n/-normal flags are only valid in -plane mode, not -box mode.\n";
      return nullptr;
    }
    if (!LxGiven || !LyGiven || (ndm == 3 && !LzGiven)) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -box requires positive PML layer thickness (-L or -Lx/-Ly"
             << (ndm == 3 ? "/-Lz" : "") << ")\n";
      return 0;
    }
    if (ndm == 2) {
      geomParams.resize(6);
      geomParams(0) = boxBounds[0];
      geomParams(1) = boxBounds[1];
      geomParams(2) = boxBounds[2];
      geomParams(3) = boxBounds[3];
      geomParams(4) = Lx;
      geomParams(5) = Ly;
    } else {
      geomParams.resize(9);
      geomParams(0) = boxBounds[0];
      geomParams(1) = boxBounds[1];
      geomParams(2) = boxBounds[2];
      geomParams(3) = boxBounds[3];
      geomParams(4) = boxBounds[4];
      geomParams(5) = boxBounds[5];
      geomParams(6) = Lx;
      geomParams(7) = Ly;
      geomParams(8) = Lz;
    }
  } else {
    // -plane mode
    if (hasLx || hasLy || hasLz || hasMultiL) {
      opserr << "WARNING element LadrunoCPML (eleTag " << tag
             << "): -Lx/-Ly/-Lz and multi-value -L are only valid in -box mode, not -plane mode.\n";
      return nullptr;
    }
    if (!x0Given || !nGiven || !LpmlGiven || Lpml <= 0.0) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -plane requires x0, normal vector n, and positive thickness Lpml (-L)\n";
      return 0;
    }
    double nNorm2 = 0.0;
    for (int d = 0; d < ndm; d++)
      nNorm2 += planeN[d] * planeN[d];
    if (std::sqrt(nNorm2) <= 1.0e-12) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -plane normal vector has zero length\n";
      return 0;
    }
    if (ndm == 2) {
      geomParams.resize(5);
      geomParams(0) = planeX0[0];
      geomParams(1) = planeX0[1];
      geomParams(2) = planeN[0];
      geomParams(3) = planeN[1];
      geomParams(4) = Lpml;
    } else {
      geomParams.resize(7);
      geomParams(0) = planeX0[0];
      geomParams(1) = planeX0[1];
      geomParams(2) = planeX0[2];
      geomParams(3) = planeN[0];
      geomParams(4) = planeN[1];
      geomParams(5) = planeN[2];
      geomParams(6) = Lpml;
    }
  }

  // ---- Coupled Biot u-p mode validation ------------------------------------
  const bool anyBiotFlag = KfGiven || poroGiven || rhoFGiven || permGiven ||
                           permHGiven || gammaWGiven || alphaBGiven || KsGiven ||
                           stabGiven || pOrderGiven || formGiven || fluidBodyGiven;
  if (!upMode) {
    if (anyBiotFlag) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- Biot u-p flag specified without -up\n";
      return 0;
    }
  }

  ladruno_cpml::CPMLBiotParams biot{};
  if (upMode) {
    if ((permGiven && permHGiven) || (permGiven && gammaWGiven)) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- specify either -perm or (-permH -gammaW), not both\n";
      return 0;
    }
    if (permHGiven != gammaWGiven) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -permH and -gammaW must be specified together\n";
      return 0;
    }
    if (!KfGiven || !poroGiven || !rhoFGiven || (!permGiven && !permHGiven)) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -up mode missing required Biot argument(s):";
      if (!KfGiven)   opserr << " -Kf";
      if (!poroGiven) opserr << " -poro";
      if (!rhoFGiven) opserr << " -rhoF";
      if (!permGiven && !permHGiven) opserr << " -perm (or -permH -gammaW)";
      opserr << "\n";
      return 0;
    }
    if (biotAlpha < poro || biotAlpha > 1.0) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- Biot alphaB=" << biotAlpha
             << " must satisfy poro <= alphaB <= 1.0 (poro=" << poro << ")\n";
      return 0;
    }
    if (KsGiven && !(Ks > 0.0)) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -Ks must be > 0 when specified\n";
      return 0;
    }

    const bool isQuadraticTopo =
        (ndm == 2 && nNodes == 6) ||
        (ndm == 3 && (nNodes == 10 || nNodes == 20));
    if (!pOrderGiven) {
      pOrder = 0;
      if (isQuadraticTopo && dom != 0) {
        Node *midNode = dom->getNode(nodeTags(nNodes - 1));
        if (midNode != 0 && midNode->getNumberDOF() == ndm) {
          pOrder = 1;
        }
      }
    }
    if (pOrder == 1 && !isQuadraticTopo) {
      opserr << "ERROR LadrunoCPML " << tag
             << " -- -pOrder linear (Taylor-Hood) is only valid on quadratic "
                "topologies (6, 10, or 20 nodes; got " << nNodes << " nodes)\n";
      return 0;
    }

    if (!fluidBodyGiven) {
      fluidBodyForces = bodyForces;
    }

    biot.enabled = true;
    biot.Kf = Kf;
    biot.poro = poro;
    biot.rhoF = rhoF;
    biot.biotAlpha = biotAlpha;
    biot.Ks = Ks;
    biot.pOrder = pOrder;
    biot.formulation = formulation;
    biot.stabMode = stabMode;
    biot.stabValue = stabValue;
    for (int d = 0; d < 3; d++) {
      if (d < ndm) {
        biot.kbar[d] = permHGiven ? (permH[d] / gammaW) : perm[d];
      } else {
        biot.kbar[d] = 0.0;
      }
    }

    static bool upSolverNoticePrinted = false;
    if (!upSolverNoticePrinted) {
      opserr << "NOTE LadrunoCPML (-up): honest-p Biot C-PML tangent is unsymmetric "
                "(-Q_eff in K, +Q_eff^T in C); use a general solver "
                "(system UmfPack / SparseGEN / FullGeneral, or MUMPS -ICNTL14 with SYM=0).\n";
      upSolverNoticePrinted = true;
    }
  }

  LadrunoCPML *ele = nullptr;
  if (!upMode) {
    ele = new LadrunoCPML(tag,
                          ndm,
                          nodeTags,
                          shapeOverride,
                          E,
                          nu,
                          rho,
                          thickness,
                          geomMode,
                          geomParams,
                          mDeg,
                          R0,
                          kappaMax,
                          alphaMax,
                          lumped,
                          initialStage,
                          bodyForces);
  } else {
    ele = new LadrunoCPML(tag,
                          ndm,
                          nodeTags,
                          shapeOverride,
                          E,
                          nu,
                          rho,
                          thickness,
                          geomMode,
                          geomParams,
                          mDeg,
                          R0,
                          kappaMax,
                          alphaMax,
                          lumped,
                          initialStage,
                          bodyForces,
                          biot,
                          fluidBodyForces);
  }
  if (ele == 0) {
    opserr << "ERROR LadrunoCPML " << tag << " -- out of memory\n";
    return 0;
  }
  if (!ele->isValid()) {
    delete ele;
    opserr << "ERROR LadrunoCPML " << tag
           << " -- element failed constructor validation\n";
    return 0;
  }
  return ele;
}
