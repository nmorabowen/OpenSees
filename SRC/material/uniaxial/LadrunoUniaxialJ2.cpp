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

// Ladruno: uniaxial combined isotropic (Voce+linear) + Chaboche Armstrong-
// Frederick kinematic von Mises (J2) UniaxialMaterial. See LadrunoUniaxialJ2.h
// and Ladruno_implementation/13_ladruno_uniaxial_j2_adr.md for the derivation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoUniaxialJ2.h>
#include <LadrunoHardening.h>     // shared sig_y(pbar) -- oracle contract w/ LadrunoJ2
#include <LadrunoDamage.h>        // shared Lemaitre damage law -- oracle contract w/ LadrunoJ2
#include "LadrunoUniaxialJ2Kernel.h"  // WP-262: templated return map (dual-number DDM)
#include <Element.h>
#include <OPS_Globals.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoUniaxialJ2 tag E -iso voce s0 Qinf b Hiso
//                                              -kin N C1 g1 C2 g2 ...
// ===========================================================================
void* OPS_LadrunoUniaxialJ2(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoUniaxialJ2 tag? E? "
           << "-iso voce s0? Qinf? b? Hiso? -kin N? C1? g1? ...\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoUniaxialJ2 tag\n";
    return 0;
  }

  double E;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &E) < 0) {
    opserr << "WARNING invalid LadrunoUniaxialJ2 E\n";
    return 0;
  }

  // defaults: perfectly plastic, no kinematic hardening, no damage
  double s0 = 0.0, Qinf = 0.0, bIso = 0.0, Hiso = 0.0;
  int nBack = 0;
  double C[LadrunoUniaxialJ2::MAXBACK], gam[LadrunoUniaxialJ2::MAXBACK];
  for (int i = 0; i < LadrunoUniaxialJ2::MAXBACK; i++) { C[i] = 0.0; gam[i] = 0.0; }
  bool dmgOn = false;
  double dmgR = 0.0, dmgS = 1.0, dmgPD = 0.0, dmgDc = 1.0;
  bool dmgImplex = false;
  double dmgLchRef = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();

    if (strcmp(flag, "-iso") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "voce") == 0 || strcmp(law, "Voce") == 0) {
        double d[4];
        numData = 4;
        if (OPS_GetDoubleInput(&numData, d) < 0) {
          opserr << "WARNING LadrunoUniaxialJ2: -iso voce wants s0 Qinf b Hiso\n";
          return 0;
        }
        s0 = d[0]; Qinf = d[1]; bIso = d[2]; Hiso = d[3];
      } else {
        opserr << "WARNING LadrunoUniaxialJ2: unknown -iso law '" << law
               << "' (only 'voce' in v1)\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-kin") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &nBack) < 0) {
        opserr << "WARNING LadrunoUniaxialJ2: -kin wants N\n";
        return 0;
      }
      if (nBack < 0 || nBack > LadrunoUniaxialJ2::MAXBACK) {
        opserr << "WARNING LadrunoUniaxialJ2: -kin N out of range [0,"
               << LadrunoUniaxialJ2::MAXBACK << "]\n";
        return 0;
      }
      for (int k = 0; k < nBack; k++) {
        double cg[2];
        numData = 2;
        if (OPS_GetDoubleInput(&numData, cg) < 0) {
          opserr << "WARNING LadrunoUniaxialJ2: -kin term " << k+1 << " wants C gamma\n";
          return 0;
        }
        C[k] = cg[0]; gam[k] = cg[1];
      }
    }
    else if (strcmp(flag, "-damage") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "lemaitre") == 0 || strcmp(law, "Lemaitre") == 0) {
        double d[4];
        numData = 4;
        if (OPS_GetDoubleInput(&numData, d) < 0) {
          opserr << "WARNING LadrunoUniaxialJ2: -damage lemaitre wants r s pD Dc\n";
          return 0;
        }
        dmgOn = true; dmgR = d[0]; dmgS = d[1]; dmgPD = d[2]; dmgDc = d[3];
        if (dmgR <= 0.0 || dmgDc <= 0.0 || dmgDc > 1.0) {
          opserr << "WARNING LadrunoUniaxialJ2: -damage lemaitre needs r>0, 0<Dc<=1\n";
          return 0;
        }
      } else if (strcmp(law, "none") == 0) {
        dmgOn = false;
      } else {
        opserr << "WARNING LadrunoUniaxialJ2: unknown -damage law '" << law
               << "' (only 'lemaitre' in v1)\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-implex") == 0) {
      dmgImplex = true;
    }
    else if (strcmp(flag, "-lchRef") == 0 || strcmp(flag, "-lch_ref") == 0 ||
             strcmp(flag, "-autoRegularization") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &dmgLchRef) < 0 || dmgLchRef < 0.0) {
        opserr << "WARNING LadrunoUniaxialJ2: -lchRef wants non-negative lchRef\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoUniaxialJ2: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  UniaxialMaterial* mat =
    new LadrunoUniaxialJ2(tag, E, s0, Qinf, bIso, Hiso, nBack, C, gam,
                          dmgOn, dmgR, dmgS, dmgPD, dmgDc, dmgImplex, dmgLchRef);
  if (mat == 0) {
    opserr << "WARNING LadrunoUniaxialJ2: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoUniaxialJ2::LadrunoUniaxialJ2(int tag, double e,
                                     double s0, double Qi, double b, double Hi,
                                     int nb, const double* C, const double* gam,
                                     bool dOn, double dR, double dS, double dPD, double dDc,
                                     bool dImplex, double dLchRef)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoUniaxialJ2),
    E(e), sig0(s0), Qinf(Qi), bIso(b), Hiso(Hi), nBack(nb),
    dmgOn(dOn), dmgR(dR), dmgS(dS), dmgPD(dPD), dmgDc(dDc), dmgImplex(dImplex),
    dmgLchRef(dLchRef),
    parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) {
    Ckin[k] = (k < nb) ? C[k]   : 0.0;
    gKin[k] = (k < nb) ? gam[k] : 0.0;
  }
  this->revertToStart();
}

LadrunoUniaxialJ2::LadrunoUniaxialJ2()
  : UniaxialMaterial(0, MAT_TAG_LadrunoUniaxialJ2),
    E(0.0), sig0(0.0), Qinf(0.0), bIso(0.0), Hiso(0.0), nBack(0),
    dmgOn(false), dmgR(0.0), dmgS(1.0), dmgPD(0.0), dmgDc(1.0), dmgImplex(false),
    dmgLchRef(0.0),
    parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) { Ckin[k] = 0.0; gKin[k] = 0.0; }
  this->revertToStart();
}

LadrunoUniaxialJ2::~LadrunoUniaxialJ2() {}

// ===========================================================================
//  hardening (shared Ladruno law -- byte-identical with 3D LadrunoJ2)
// ===========================================================================
double LadrunoUniaxialJ2::yieldStress(double p) const
{
  return Ladruno::yieldStressVoceLinear(p, sig0, Qinf, bIso, Hiso);
}

double LadrunoUniaxialJ2::yieldSlope(double p) const
{
  return Ladruno::yieldSlopeVoceLinear(p, sig0, Qinf, bIso, Hiso);
}

// ===========================================================================
//  the return map: fixed direction s = sign(sigma_tr - X_n) + scalar Newton on
//  dp (= accumulated plastic strain increment). No M(dGamma) tensor solve.
// ===========================================================================
// WP-262: pack the instance's constants / committed history for the kernel.
template <class P>
static void fillKernelParams(ladruno_uniaxial_j2::ParamsT<P>& p, double E, double sig0,
                             double Qinf, double bIso, double Hiso, int nBack,
                             const double* Ckin, const double* gKin,
                             bool dmgOn, double dmgR, double dmgS, double dmgPD, double dmgDc,
                             bool dmgImplex, double lchScale)
{
  p.E = E; p.sig0 = sig0; p.Qinf = Qinf; p.bIso = bIso; p.Hiso = Hiso;
  p.nBack = nBack;
  for (int k = 0; k < ladruno_uniaxial_j2::MAXBACK; k++) { p.C[k] = Ckin[k]; p.gam[k] = gKin[k]; }
  p.dmgOn = dmgOn; p.dmgR = dmgR; p.dmgS = dmgS; p.dmgPD = dmgPD; p.dmgDc = dmgDc;
  p.dmgImplex = dmgImplex; p.lchScale = lchScale;
}

template <class T>
static void committedState(ladruno_uniaxial_j2::StateT<T>& C, double plasticStrain, double ebarp,
                           double dGamma, double wp, double D, const double* X, int nBack)
{
  C.plasticStrain = plasticStrain; C.ebarp = ebarp; C.dGamma = dGamma; C.wp = wp; C.D = D;
  for (int k = 0; k < ladruno_uniaxial_j2::MAXBACK; k++) C.X[k] = (k < nBack) ? X[k] : 0.0;
}

int LadrunoUniaxialJ2::setTrialStrain(double strain, double)
{
  if (fabs(Tstrain - strain) < DBL_EPSILON)
    return 0;

  Tstrain = strain;

  // WP-262: the return map lives in the templated kernel
  // (LadrunoUniaxialJ2Kernel.h); this is its double instantiation, which performs
  // exactly the operations of the former in-class code (bit-identity gate:
  // tests/ladruno_uniaxialj2_dual_check.cpp). The crack-band scale lch/lch_ref is
  // read from the active element only when the regularization is on (it is 1
  // otherwise, exactly as before), and LATCHED: the DDM pass replays this return
  // map outside Domain::update, where ops_TheActiveElement no longer names this
  // material's element (WP-262 review finding 1).
  ladruno_uniaxial_j2::ParamsT<double> p;
  fillKernelParams(p, E, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin,
                   dmgOn, dmgR, dmgS, dmgPD, dmgDc, dmgImplex, 1.0);
  ladruno_uniaxial_j2::StateT<double> C, Tr;
  committedState(C, CplasticStrain, Cebarp, CdGamma, Cwp, CD, CX, nBack);
  if (dmgOn && dmgLchRef > 0.0) p.lchScale = this->lchScaleFromElement();
  lchScaleTrial = p.lchScale;

  double resid = 0.0;
  int status = ladruno_uniaxial_j2::returnMap(p, Tstrain, C, CdD, Tr, Tstress, Ttangent, &resid);

  TplasticStrain = Tr.plasticStrain;
  Tebarp = Tr.ebarp;
  TdGamma = Tr.dGamma;
  Twp = Tr.wp;
  TD = Tr.D;
  for (int k = 0; k < nBack; k++) TX[k] = Tr.X[k];

  if (status == ladruno_uniaxial_j2::STATUS_SINGULAR)
    opserr << "WARNING LadrunoUniaxialJ2: singular local Jacobian (tag "
           << this->getTag() << ")\n";
  else if (status == ladruno_uniaxial_j2::STATUS_NO_CONVERGE)
    opserr << "WARNING LadrunoUniaxialJ2: local Newton did not converge (tag "
           << this->getTag() << ", |R|=" << resid << ")\n";

  return 0;
}

double LadrunoUniaxialJ2::lchScaleFromElement(void) const
{
  const double lch = (ops_TheActiveElement != 0) ? ops_TheActiveElement->getCharacteristicLength() : 0.0;
  return (dmgLchRef > 0.0 && lch > 0.0) ? (lch / dmgLchRef) : 1.0;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoUniaxialJ2::commitState(void)
{
  CplasticStrain = TplasticStrain;
  Cebarp = Tebarp;
  CdGamma = TdGamma;
  Cwp = Twp;
  CdD = TD - CD;                 // implicit increment (for the next IMPL-EX extrapolation)
  CD = TD;
  for (int k = 0; k < nBack; k++) CX[k] = TX[k];
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  return 0;
}

int LadrunoUniaxialJ2::revertToLastCommit(void)
{
  TplasticStrain = CplasticStrain;
  Tebarp = Cebarp;
  TdGamma = CdGamma;
  Twp = Cwp;
  TD = CD;
  for (int k = 0; k < nBack; k++) TX[k] = CX[k];
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  return 0;
}

int LadrunoUniaxialJ2::revertToStart(void)
{
  SHV.clear();                                   // WP-262: history sensitivities
  CplasticStrain = 0.0; Cebarp = 0.0; CdGamma = 0.0; Cwp = 0.0; CD = 0.0; CdD = 0.0;
  TplasticStrain = 0.0; Tebarp = 0.0; TdGamma = 0.0; Twp = 0.0; TD = 0.0;
  for (int k = 0; k < MAXBACK; k++) { CX[k] = 0.0; TX[k] = 0.0; }
  Cstrain = 0.0;
  Cstress = 0.0;
  Ctangent = E;
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = E;
  return 0;
}

double LadrunoUniaxialJ2::getEnergy(void)
{
  return Twp + 0.5 * Tstress * (Tstrain - TplasticStrain);
}

// ===========================================================================
//  copy
// ===========================================================================
UniaxialMaterial* LadrunoUniaxialJ2::getCopy(void)
{
  LadrunoUniaxialJ2* theCopy =
    new LadrunoUniaxialJ2(this->getTag(), E, sig0, Qinf, bIso, Hiso,
                          nBack, Ckin, gKin, dmgOn, dmgR, dmgS, dmgPD, dmgDc,
                          dmgImplex, dmgLchRef);

  theCopy->CplasticStrain = CplasticStrain;
  theCopy->Cebarp = Cebarp;
  theCopy->CdGamma = CdGamma;
  theCopy->Cwp = Cwp;
  theCopy->CD = CD;
  theCopy->CdD = CdD;
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ctangent = Ctangent;
  theCopy->TplasticStrain = TplasticStrain;
  theCopy->Tebarp = Tebarp;
  theCopy->TdGamma = TdGamma;
  theCopy->Twp = Twp;
  theCopy->TD = TD;
  for (int k = 0; k < MAXBACK; k++) { theCopy->CX[k] = CX[k]; theCopy->TX[k] = TX[k]; }
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  theCopy->parameterID = parameterID;
  theCopy->SHV = SHV;                            // WP-262: history sensitivities
  theCopy->lchScaleTrial = lchScaleTrial;

  return theCopy;
}

// ===========================================================================
//  parallel
// ===========================================================================
int LadrunoUniaxialJ2::sendSelf(int cTag, Channel& theChannel)
{
  // tag + 5 params + nBack + 7 damage params + 2*MAXBACK (C,gamma) + committed
  //   history (CplasticStrain, Cebarp, CdGamma, Cwp, CD, CdD + MAXBACK CX) + Cstrain, Cstress, Ctangent
  static thread_local Vector data(1 + 5 + 1 + 7 + 2*MAXBACK + 4 + 1 + 1 + MAXBACK + 3);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = E;
  data(c++) = sig0;
  data(c++) = Qinf;
  data(c++) = bIso;
  data(c++) = Hiso;
  data(c++) = nBack;
  data(c++) = dmgOn ? 1.0 : 0.0;
  data(c++) = dmgR;
  data(c++) = dmgS;
  data(c++) = dmgPD;
  data(c++) = dmgDc;
  data(c++) = dmgImplex ? 1.0 : 0.0;
  data(c++) = dmgLchRef;
  for (int k = 0; k < MAXBACK; k++) data(c++) = Ckin[k];
  for (int k = 0; k < MAXBACK; k++) data(c++) = gKin[k];
  data(c++) = CplasticStrain;
  data(c++) = Cebarp;
  data(c++) = CdGamma;
  data(c++) = Cwp;
  data(c++) = CD;
  data(c++) = CdD;
  for (int k = 0; k < MAXBACK; k++) data(c++) = CX[k];
  data(c++) = Cstrain;
  data(c++) = Cstress;
  data(c++) = Ctangent;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoUniaxialJ2::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoUniaxialJ2::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static thread_local Vector data(1 + 5 + 1 + 7 + 2*MAXBACK + 4 + 1 + 1 + MAXBACK + 3);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoUniaxialJ2::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  E    = data(c++);
  sig0 = data(c++);
  Qinf = data(c++);
  bIso = data(c++);
  Hiso = data(c++);
  nBack = (int)data(c++);
  dmgOn = (data(c++) != 0.0);
  dmgR  = data(c++);
  dmgS  = data(c++);
  dmgPD = data(c++);
  dmgDc = data(c++);
  dmgImplex = (data(c++) != 0.0);
  dmgLchRef = data(c++);
  for (int k = 0; k < MAXBACK; k++) Ckin[k] = data(c++);
  for (int k = 0; k < MAXBACK; k++) gKin[k] = data(c++);
  CplasticStrain = data(c++);
  Cebarp = data(c++);
  CdGamma = data(c++);
  Cwp = data(c++);
  CD = data(c++);
  CdD = data(c++);
  for (int k = 0; k < MAXBACK; k++) CX[k] = data(c++);
  Cstrain = data(c++);
  Cstress = data(c++);
  Ctangent = data(c++);

  // sync trial to committed
  this->revertToLastCommit();
  return 0;
}

// ===========================================================================
//  parameters / print
// ===========================================================================
int LadrunoUniaxialJ2::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  // Scalar elastic / isotropic parameters (FE sensitivity + updateMaterialStage).
  // NOTE: the dispatch is intentionally a flat strcmp ladder so future layers
  // (e.g. a local-buckling slenderness key, ADR 13 sec 6.2) can be added here
  // without reshaping it -- setTrialStrain stays scalar.
  if (strcmp(argv[0], "E") == 0)                                  { param.setValue(E);    return param.addObject(1, this); }
  if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "sig0") == 0 ||
      strcmp(argv[0], "fy") == 0 || strcmp(argv[0], "Fy") == 0)   { param.setValue(sig0); return param.addObject(2, this); }
  if (strcmp(argv[0], "Qinf") == 0)                               { param.setValue(Qinf); return param.addObject(3, this); }
  if (strcmp(argv[0], "b") == 0)                                  { param.setValue(bIso); return param.addObject(4, this); }
  if (strcmp(argv[0], "Hiso") == 0 || strcmp(argv[0], "H_iso") == 0) { param.setValue(Hiso); return param.addObject(5, this); }
  if (strcmp(argv[0], "lchRef") == 0 || strcmp(argv[0], "lch_ref") == 0) { param.setValue(dmgLchRef); return param.addObject(6, this); }
  return -1;
}

int LadrunoUniaxialJ2::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case 1: E         = info.theDouble; return 0;
    case 2: sig0      = info.theDouble; return 0;
    case 3: Qinf      = info.theDouble; return 0;
    case 4: bIso      = info.theDouble; return 0;
    case 5: Hiso      = info.theDouble; return 0;
    case 6: dmgLchRef = info.theDouble; return 0;
    default: return -1;
  }
}

int LadrunoUniaxialJ2::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

// ===========================================================================
//  WP-262: DDM sensitivities by dual numbers
//
//  The OpenSees DDM contract (StaticAnalysis: computeSensitivities BEFORE
//  commit): getStressSensitivity(g) is d(sigma)/d(theta) at the CURRENT trial
//  strain held fixed, given the committed history and ITS sensitivities;
//  commitSensitivity(d eps/d theta, g, n) then advances the history
//  sensitivities with the converged total strain derivative.
//
//  Both are one pass of the templated kernel in Dual<1>: the active parameter
//  (parameterID, set by Parameter::activate; 0 when this material is not the
//  one being differentiated) is seeded with 1, the committed history with its
//  stored sensitivities, the strain with d eps/d theta (0 for the conditional
//  stress sensitivity). The local Newton is closed with the IFT step, so the
//  result is the exact derivative of the discrete return map, history included.
//  Not differentiated: the damage constants and lchRef (parameter 6), which
//  return 0 here like every material without a derivative for them.
// ===========================================================================
int LadrunoUniaxialJ2::sensitivityPass(int gradIndex, double strainGradient,
                                       double& dStress, double& dTangent, double* dTrialState)
{
  typedef ladruno_dual::Dual<1> D1;
  ladruno_uniaxial_j2::ParamsT<D1> p;
  fillKernelParams(p, E, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin,
                   dmgOn, dmgR, dmgS, dmgPD, dmgDc, dmgImplex, 1.0);
  // the scale setTrialStrain used for THIS trial state: never re-read the active
  // element here (outside Domain::update it is the last element updated)
  p.lchScale = lchScaleTrial;
  switch (parameterID) {
    case 1: p.E.d[0]    = 1.0; break;
    case 2: p.sig0.d[0] = 1.0; break;
    case 3: p.Qinf.d[0] = 1.0; break;
    case 4: p.bIso.d[0] = 1.0; break;
    case 5: p.Hiso.d[0] = 1.0; break;
    default: break;                        // not this material's parameter (or lchRef)
  }

  ladruno_uniaxial_j2::StateT<D1> C, Tr;
  committedState(C, CplasticStrain, Cebarp, CdGamma, Cwp, CD, CX, nBack);
  D1 CdDd(CdD);
  const size_t base = (size_t)gradIndex * SHV_N;
  if (gradIndex >= 0 && SHV.size() >= base + SHV_N) {
    const double* sh = &SHV[base];
    C.plasticStrain.d[0] = sh[0];
    C.ebarp.d[0]         = sh[1];
    C.D.d[0]             = sh[2];
    CdDd.d[0]            = sh[3];
    for (int k = 0; k < nBack; k++) C.X[k].d[0] = sh[4 + k];
  }

  D1 strain(Tstrain);
  strain.d[0] = strainGradient;
  D1 stress, tangent;
  int status = ladruno_uniaxial_j2::returnMap(p, strain, C, CdDd, Tr, stress, tangent);

  dStress = stress.d[0];
  dTangent = tangent.d[0];
  if (dTrialState) {
    dTrialState[0] = Tr.plasticStrain.d[0];
    dTrialState[1] = Tr.ebarp.d[0];
    dTrialState[2] = Tr.D.d[0];
    dTrialState[3] = Tr.D.d[0] - C.D.d[0];          // commitState sets CdD = TD - CD
    for (int k = 0; k < MAXBACK; k++) dTrialState[4 + k] = (k < nBack) ? Tr.X[k].d[0] : 0.0;
  }
  return status;
}

double LadrunoUniaxialJ2::getStressSensitivity(int gradIndex, bool)
{
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, 0.0, dS, dT, 0);
  return dS;
}

double LadrunoUniaxialJ2::getTangentSensitivity(int gradIndex)
{
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, 0.0, dS, dT, 0);
  return dT;
}

double LadrunoUniaxialJ2::getInitialTangentSensitivity(int)
{
  return (parameterID == 1) ? 1.0 : 0.0;    // getInitialTangent() == E
}

int LadrunoUniaxialJ2::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
  if (gradIndex < 0 || numGrads <= gradIndex) return -1;
  if (SHV.size() < (size_t)numGrads * SHV_N) SHV.resize((size_t)numGrads * SHV_N, 0.0);
  double dTr[SHV_N];
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, strainGradient, dS, dT, dTr);
  for (int i = 0; i < SHV_N; i++) SHV[(size_t)gradIndex * SHV_N + i] = dTr[i];
  return 0;
}

void LadrunoUniaxialJ2::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoUniaxialJ2\", ";
    s << "\"E\": " << E << ", ";
    s << "\"sig0\": " << sig0 << ", \"Qinf\": " << Qinf << ", \"b\": " << bIso
      << ", \"Hiso\": " << Hiso << ", ";
    s << "\"nBack\": " << nBack << "}";
    return;
  }
  s << "LadrunoUniaxialJ2 (uniaxial combined iso + Chaboche AF kinematic J2)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  E      : " << E << endln;
  s << "  iso    : sig0=" << sig0 << " Qinf=" << Qinf << " b=" << bIso << " Hiso=" << Hiso << endln;
  s << "  nBack  : " << nBack << endln;
  for (int k = 0; k < nBack; k++)
    s << "    term " << k+1 << ": C=" << Ckin[k] << " gamma=" << gKin[k] << endln;
  if (dmgOn)
    s << "  damage : Lemaitre r=" << dmgR << " s=" << dmgS
      << " pD=" << dmgPD << " Dc=" << dmgDc << "  (D=" << CD << ")" << endln;
}

// ===========================================================================
//  recordable responses (stress/strain/tangent + internal state hooks)
// ===========================================================================
Response* LadrunoUniaxialJ2::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return UniaxialMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "backStress") == 0 || strcmp(a, "backstress") == 0 || strcmp(a, "alpha") == 0)
    return new MaterialResponse(this, 100, 0.0);
  if (strcmp(a, "plasticStrain") == 0 || strcmp(a, "plasticStrains") == 0)
    return new MaterialResponse(this, 101, 0.0);
  if (strcmp(a, "equivalentPlasticStrain") == 0 || strcmp(a, "plasticStrainEq") == 0 ||
      strcmp(a, "ebarP") == 0)
    return new MaterialResponse(this, 102, 0.0);
  if (strcmp(a, "plasticWork") == 0 || strcmp(a, "plasticEnergy") == 0)
    return new MaterialResponse(this, 103, 0.0);
  if (strcmp(a, "damage") == 0 || strcmp(a, "D") == 0)
    return new MaterialResponse(this, 104, 0.0);
  if (strcmp(a, "failed") == 0 || strcmp(a, "rupture") == 0)
    return new MaterialResponse(this, 105, 0.0);
  if (strcmp(a, "Y") == 0 || strcmp(a, "energyReleaseRate") == 0)
    return new MaterialResponse(this, 106, 0.0);

  return UniaxialMaterial::setResponse(argv, argc, s);
}

int LadrunoUniaxialJ2::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 100: {                              // total backstress X = sum_k X_k
      double tot = 0.0;
      for (int k = 0; k < nBack; k++) tot += TX[k];
      matInfo.theDouble = tot;
      return 0;
    }
    case 101: matInfo.theDouble = TplasticStrain; return 0;
    case 102: matInfo.theDouble = Tebarp;         return 0;
    case 103: matInfo.theDouble = Twp;            return 0;
    case 104: matInfo.theDouble = TD;             return 0;   // Lemaitre damage
    case 105: matInfo.theDouble = (dmgOn && TD >= dmgDc) ? 1.0 : 0.0; return 0; // failed
    case 106: {                                                // damage energy release rate Y
      const double sigEff = E * (Tstrain - TplasticStrain);
      matInfo.theDouble = (E > 0.0) ? sigEff*sigEff/(2.0*E) : 0.0;  // R_v = 1 at uniaxial
      return 0;
    }
    default:
      return UniaxialMaterial::getResponse(responseID, matInfo);
  }
}
