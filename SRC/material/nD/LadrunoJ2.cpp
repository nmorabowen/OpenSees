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

// Ladruno: combined isotropic (Voce+linear) + Chaboche Armstrong-Frederick
// kinematic von Mises (J2) nDMaterial. See LadrunoJ2.h and
// Ladruno_implementation/10_ladruno_j2_plasticity.md for the full derivation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoJ2.h>
#include <LadrunoJ2Kernel.h>    // the return map (which consumes LadrunoHardening.h)
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

// The combined-hardening von Mises return map and the elastic tangent live in the
// header-only, OpenSees-free LadrunoJ2Kernel.h so the SAME verified map serves
// both this small-strain material AND the finite-strain path (LadrunoJ2 wrapped
// by LogStrainNDMaterial). See LadrunoJ2Kernel.h.
using ladruno_j2_kernel::Params;

// ===========================================================================
//  OPS parser
//   nDMaterial LadrunoJ2 tag K G -iso voce s0 Qinf b Hiso -kin N C1 g1 ... [-rho r]
// ===========================================================================
void* OPS_LadrunoJ2(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoJ2 tag? K? G? "
           << "-iso voce s0? Qinf? b? Hiso? -kin N? C1? g1? ... <-rho rho?>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoJ2 tag\n";
    return 0;
  }

  double KG[2];
  numData = 2;
  if (OPS_GetDoubleInput(&numData, KG) < 0) {
    opserr << "WARNING invalid LadrunoJ2 K G\n";
    return 0;
  }
  double K = KG[0], G = KG[1];

  // defaults: perfectly plastic, no kinematic hardening
  double s0 = 0.0, Qinf = 0.0, bIso = 0.0, Hiso = 0.0, rho = 0.0;
  int nBack = 0;
  double C[LadrunoJ2::MAXBACK], gam[LadrunoJ2::MAXBACK];
  for (int i = 0; i < LadrunoJ2::MAXBACK; i++) { C[i] = 0.0; gam[i] = 0.0; }

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();

    if (strcmp(flag, "-iso") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "voce") == 0 || strcmp(law, "Voce") == 0) {
        double d[4];
        numData = 4;
        if (OPS_GetDoubleInput(&numData, d) < 0) {
          opserr << "WARNING LadrunoJ2: -iso voce wants s0 Qinf b Hiso\n";
          return 0;
        }
        s0 = d[0]; Qinf = d[1]; bIso = d[2]; Hiso = d[3];
      } else {
        opserr << "WARNING LadrunoJ2: unknown -iso law '" << law
               << "' (only 'voce' in v1)\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-kin") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &nBack) < 0) {
        opserr << "WARNING LadrunoJ2: -kin wants N\n";
        return 0;
      }
      if (nBack < 0 || nBack > LadrunoJ2::MAXBACK) {
        opserr << "WARNING LadrunoJ2: -kin N out of range [0," << LadrunoJ2::MAXBACK << "]\n";
        return 0;
      }
      for (int k = 0; k < nBack; k++) {
        double cg[2];
        numData = 2;
        if (OPS_GetDoubleInput(&numData, cg) < 0) {
          opserr << "WARNING LadrunoJ2: -kin term " << k+1 << " wants C gamma\n";
          return 0;
        }
        C[k] = cg[0]; gam[k] = cg[1];
      }
    }
    else if (strcmp(flag, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0) {
        opserr << "WARNING LadrunoJ2: -rho wants rho\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoJ2: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  NDMaterial* mat = new LadrunoJ2(tag, K, G, s0, Qinf, bIso, Hiso, nBack, C, gam, rho);
  if (mat == 0) {
    opserr << "WARNING LadrunoJ2: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoJ2::LadrunoJ2()
  : NDMaterial(0, ND_TAG_LadrunoJ2),
    bulk(0.0), shear(0.0), sig0(0.0), Qinf(0.0), bIso(0.0), Hiso(0.0), rho(0.0),
    nBack(0), dim(DIM_3D), ncomp(6), condense(false), cEps22(0.0),
    dGamma_n(0.0), ebarP_n(0.0), ebarP(0.0), dGammaTrial(0.0), parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) { Ckin[k] = 0.0; gKin[k] = 0.0; }
  this->setupDim();
  this->revertToStart();
}

LadrunoJ2::LadrunoJ2(int tag, double K, double G,
                     double s0, double Qi, double b, double Hi,
                     int nb, const double* C, const double* gam, double r, int dimMode)
  : NDMaterial(tag, ND_TAG_LadrunoJ2),
    bulk(K), shear(G), sig0(s0), Qinf(Qi), bIso(b), Hiso(Hi), rho(r),
    nBack(nb), dim(dimMode), ncomp(6), condense(false), cEps22(0.0),
    dGamma_n(0.0), ebarP_n(0.0), ebarP(0.0), dGammaTrial(0.0), parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) {
    Ckin[k] = (k < nb) ? C[k]   : 0.0;
    gKin[k] = (k < nb) ? gam[k] : 0.0;
  }
  this->setupDim();
  this->revertToStart();
}

// ---------------------------------------------------------------------------
//  dimensional view setup: reduced-vector order -> full 6-comp index map
// ---------------------------------------------------------------------------
void LadrunoJ2::setupDim(void)
{
  // full ordering: 0:00 1:11 2:22 3:01 4:12 5:02
  switch (dim) {
    case DIM_PSTRAIN:    ncomp = 3; { int m[3]={0,1,3};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=false; break;
    case DIM_AXISYM:     ncomp = 4; { int m[4]={0,1,2,3};     for(int a=0;a<4;a++) vmap[a]=m[a]; } condense=false; break;
    case DIM_PLATEFIBER: ncomp = 5; { int m[5]={0,1,3,4,5};   for(int a=0;a<5;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_PSTRESS:    ncomp = 3; { int m[3]={0,1,3};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_3D:
    default:             ncomp = 6; { int m[6]={0,1,2,3,4,5}; for(int a=0;a<6;a++) vmap[a]=m[a]; } condense=false; break;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
}

// Static condensation of the out-of-plane (33, index 2) dof so sigma_22 = 0:
//   Dtan[I][J] -= Dtan[I][2] Dtan[2][J] / Dtan[2][2].
void LadrunoJ2::condenseTangent(void)
{
  double d22 = Dtan[2][2];
  if (fabs(d22) < 1.0e-300) return;
  double col[6], row[6];
  for (int i = 0; i < 6; i++) { col[i] = Dtan[i][2]; row[i] = Dtan[2][i]; }
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan[I][J] -= col[I]*row[J]/d22;
}

LadrunoJ2::~LadrunoJ2() {}

// ===========================================================================
//  kernel parameter pack (copy the instance's material constants)
// ===========================================================================
static void fillParams(Params& p, double bulk, double shear,
                       double sig0, double Qinf, double bIso, double Hiso,
                       int nBack, const double* Ckin, const double* gKin)
{
  p.K = bulk; p.G = shear;
  p.sig0 = sig0; p.Qinf = Qinf; p.bIso = bIso; p.Hiso = Hiso;
  p.nBack = nBack;
  for (int k = 0; k < ladruno_j2_kernel::MAXBACK; k++) {
    p.C[k]   = Ckin[k];
    p.gam[k] = gKin[k];
  }
}

// ===========================================================================
//  strain interface
// ===========================================================================
int LadrunoJ2::setTrialStrain(const Vector& e)
{
  // map the element's reduced (engineering-shear) vector into the full
  // 6-comp tensor strain. For condensed modes (PlaneStress/PlateFiber) the
  // out-of-plane eps_22 (index 2) is an internal unknown — preserve its
  // running value as the starting guess.
  double eps22 = strain6[2];
  for (int i = 0; i < 6; i++) strain6[i] = 0.0;
  if (condense) strain6[2] = eps22;

  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double val = e(a);
    if (full >= 3) val *= 0.5;     // engineering -> tensor shear
    strain6[full] = val;
  }

  if (!condense) {
    this->integrate();
    return 0;
  }

  // enforce sigma_22 = 0: Newton on eps_22 (= strain6[2]); dSNPO sec 9.2.3.
  const int maxIt = 25;
  for (int it = 0; it < maxIt; it++) {
    this->integrate();
    double d22 = Dtan[2][2];
    double smag = 0.0;
    for (int i = 0; i < 6; i++) smag += stress6[i]*stress6[i];
    smag = sqrt(smag);
    double tol22 = 1.0e-10 * (smag > 1.0 ? smag : 1.0);
    if (fabs(stress6[2]) <= tol22) break;
    if (fabs(d22) < 1.0e-300) break;
    strain6[2] -= stress6[2] / d22;
    if (it == maxIt - 1)
      opserr << "WARNING LadrunoJ2: sigma_22 condensation did not converge (tag "
             << this->getTag() << ", |s22|=" << fabs(stress6[2]) << ")\n";
  }
  this->condenseTangent();
  return 0;
}

int LadrunoJ2::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoJ2::setTrialStrainIncr(const Vector& v)
{
  // current reduced strain + increment
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double cur = (full >= 3) ? 2.0*strain6[full] : strain6[full];
    ne(a) = cur + v(a);
  }
  return this->setTrialStrain(ne);
}

int LadrunoJ2::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  the return map (scalar Newton, Kobayashi-Ohno) + consistent tangent
//
//  Thin wrapper: pack the instance's parameters + committed history, run the
//  header-only kernel (LadrunoJ2Kernel.h — the SINGLE source of the return map,
//  shared with the finite-strain path via LogStrainNDMaterial), then surface its
//  diagnostics as the warnings the inline version used to print.
// ===========================================================================
void LadrunoJ2::integrate(void)
{
  Params p;
  fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);

  double resid = 0.0;
  int status = ladruno_j2_kernel::returnMap(
      p, strain6, epsP_n, ebarP_n, alpha_n,
      stress6, Dtan, epsP, ebarP, alpha, dGammaTrial, &resid);

  if (status == ladruno_j2_kernel::STATUS_SINGULAR)
    opserr << "WARNING LadrunoJ2: singular local Jacobian (tag "
           << this->getTag() << ")\n";
  else if (status == ladruno_j2_kernel::STATUS_NO_CONVERGE)
    opserr << "WARNING LadrunoJ2: local Newton did not converge (tag "
           << this->getTag() << ", |R|=" << resid << ")\n";
}

void LadrunoJ2::buildElasticTangent(double Kt[6][6]) const
{
  Params p;
  fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);
  ladruno_j2_kernel::elasticTangent(p, Kt);
}

void LadrunoJ2::setStateToCommitted(void)
{
  for (int i = 0; i < 6; i++) epsP[i] = epsP_n[i];
  ebarP = ebarP_n;
  dGammaTrial = 0.0;
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++) alpha[k][i] = alpha_n[k][i];
}

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoJ2::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = stress6[vmap[a]];  // shear = true stress
  return stressOut;
}

const Vector& LadrunoJ2::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    strainOut(a) = (full >= 3) ? 2.0*strain6[full] : strain6[full];  // tensor -> engineering
  }
  return strainOut;
}

const Matrix& LadrunoJ2::getTangent(void)
{
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Dtan[vmap[a]][vmap[b]];
  return tangentOut;
}

const Matrix& LadrunoJ2::getInitialTangent(void)
{
  // elastic tangent mapped to the reduced view (uncondensed, matching the
  // J2Plasticity specializations' getInitialTangent convention).
  double Kt[6][6];
  this->buildElasticTangent(Kt);
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Kt[vmap[a]][vmap[b]];
  return tangentOut;
}

const char* LadrunoJ2::getType(void) const
{
  switch (dim) {
    case DIM_PSTRAIN:    return "PlaneStrain";
    case DIM_AXISYM:     return "AxiSymmetric";
    case DIM_PLATEFIBER: return "PlateFiber";
    case DIM_PSTRESS:    return "PlaneStress";
    default:             return "ThreeDimensional";
  }
}

int LadrunoJ2::getOrder(void) const { return ncomp; }

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoJ2::commitState(void)
{
  for (int i = 0; i < 6; i++) epsP_n[i] = epsP[i];
  ebarP_n  = ebarP;
  dGamma_n = dGammaTrial;
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++) alpha_n[k][i] = alpha[k][i];
  cEps22 = strain6[2];   // converged out-of-plane strain (condensed modes)
  return 0;
}

int LadrunoJ2::revertToLastCommit(void)
{
  this->setStateToCommitted();
  if (condense) strain6[2] = cEps22;
  return 0;
}

int LadrunoJ2::revertToStart(void)
{
  for (int i = 0; i < 6; i++) {
    strain6[i] = 0.0;
    stress6[i] = 0.0;
    epsP_n[i]  = 0.0;
    epsP[i]    = 0.0;
  }
  ebarP_n = 0.0; ebarP = 0.0;
  dGamma_n = 0.0; dGammaTrial = 0.0;
  cEps22 = 0.0;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) { alpha_n[k][i] = 0.0; alpha[k][i] = 0.0; }
  this->buildElasticTangent(Dtan);
  return 0;
}

// ===========================================================================
//  copies
// ===========================================================================
NDMaterial* LadrunoJ2::getCopy(void)
{
  return new LadrunoJ2(this->getTag(), bulk, shear, sig0, Qinf, bIso, Hiso,
                       nBack, Ckin, gKin, rho, dim);
}

NDMaterial* LadrunoJ2::getCopy(const char* type)
{
  int d = -1;
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) d = DIM_3D;
  else if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0) d = DIM_PSTRAIN;
  else if (strcmp(type, "AxiSymmetric") == 0 || strcmp(type, "AxiSymmetric2D") == 0) d = DIM_AXISYM;
  else if (strcmp(type, "PlateFiber") == 0) d = DIM_PLATEFIBER;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) d = DIM_PSTRESS;

  if (d < 0)
    return NDMaterial::getCopy(type);   // let the base report the unsupported type

  return new LadrunoJ2(this->getTag(), bulk, shear, sig0, Qinf, bIso, Hiso,
                       nBack, Ckin, gKin, rho, d);
}

// ===========================================================================
//  parallel
// ===========================================================================
int LadrunoJ2::sendSelf(int commitTag, Channel& theChannel)
{
  // tag + 7 params + nBack + dim + 2*MAXBACK (C,gamma) + epsP(6) + ebarP(1)
  //     + alpha(MAXBACK*6) + dGamma(1) + cEps22(1)
  static Vector data(1 + 7 + 1 + 1 + 2*MAXBACK + 6 + 1 + MAXBACK*6 + 1 + 1);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = bulk;
  data(c++) = shear;
  data(c++) = sig0;
  data(c++) = Qinf;
  data(c++) = bIso;
  data(c++) = Hiso;
  data(c++) = rho;
  data(c++) = nBack;
  data(c++) = dim;
  for (int k = 0; k < MAXBACK; k++) data(c++) = Ckin[k];
  for (int k = 0; k < MAXBACK; k++) data(c++) = gKin[k];
  for (int i = 0; i < 6; i++) data(c++) = epsP_n[i];
  data(c++) = ebarP_n;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) data(c++) = alpha_n[k][i];
  data(c++) = dGamma_n;
  data(c++) = cEps22;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoJ2::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static Vector data(1 + 7 + 1 + 1 + 2*MAXBACK + 6 + 1 + MAXBACK*6 + 1 + 1);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  bulk  = data(c++);
  shear = data(c++);
  sig0  = data(c++);
  Qinf  = data(c++);
  bIso  = data(c++);
  Hiso  = data(c++);
  rho   = data(c++);
  nBack = (int)data(c++);
  dim   = (int)data(c++);
  for (int k = 0; k < MAXBACK; k++) Ckin[k] = data(c++);
  for (int k = 0; k < MAXBACK; k++) gKin[k] = data(c++);
  for (int i = 0; i < 6; i++) epsP_n[i] = data(c++);
  ebarP_n = data(c++);
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) alpha_n[k][i] = data(c++);
  dGamma_n = data(c++);
  cEps22   = data(c++);

  this->setupDim();             // rebuild vmap/ncomp/condense + resize buffers
  // sync trial to committed
  this->setStateToCommitted();
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  if (condense) strain6[2] = cEps22;
  this->buildElasticTangent(Dtan);
  return 0;
}

// ===========================================================================
//  parameters / print
// ===========================================================================
int LadrunoJ2::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "K") == 0)   return param.addObject(1, this);
  if (strcmp(argv[0], "G") == 0 || strcmp(argv[0], "mu") == 0) return param.addObject(2, this);
  if (strcmp(argv[0], "rho") == 0) return param.addObject(3, this);
  if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "sig0") == 0) return param.addObject(4, this);
  if (strcmp(argv[0], "Hiso") == 0) return param.addObject(5, this);
  if (strcmp(argv[0], "Qinf") == 0) return param.addObject(6, this);
  if (strcmp(argv[0], "b") == 0)    return param.addObject(7, this);
  return -1;
}

int LadrunoJ2::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case 1: bulk  = info.theDouble; return 0;
    case 2: shear = info.theDouble; return 0;
    case 3: rho   = info.theDouble; return 0;
    case 4: sig0  = info.theDouble; return 0;
    case 5: Hiso  = info.theDouble; return 0;
    case 6: Qinf  = info.theDouble; return 0;
    case 7: bIso  = info.theDouble; return 0;
    default: return -1;
  }
}

int LadrunoJ2::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

void LadrunoJ2::Print(OPS_Stream& s, int flag)
{
  s << endln;
  s << "LadrunoJ2 (combined isotropic + Chaboche AF kinematic J2)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  K, G   : " << bulk << ", " << shear << endln;
  s << "  iso    : sig0=" << sig0 << " Qinf=" << Qinf << " b=" << bIso << " Hiso=" << Hiso << endln;
  s << "  nBack  : " << nBack << endln;
  for (int k = 0; k < nBack; k++)
    s << "    term " << k+1 << ": C=" << Ckin[k] << " gamma=" << gKin[k] << endln;
  s << "  rho    : " << rho << endln;
}

// ===========================================================================
//  recordable responses (stress/strain/tangent + internal state)
// ===========================================================================
Response* LadrunoJ2::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "backStress") == 0 || strcmp(a, "backstress") == 0 || strcmp(a, "alpha") == 0)
    return new MaterialResponse(this, 4, Vector(ncomp));
  if (strcmp(a, "plasticStrain") == 0 || strcmp(a, "plasticStrains") == 0)
    return new MaterialResponse(this, 5, Vector(ncomp));
  if (strcmp(a, "equivalentPlasticStrain") == 0 || strcmp(a, "plasticStrainEq") == 0 ||
      strcmp(a, "ebarP") == 0)
    return new MaterialResponse(this, 6, Vector(1));

  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoJ2::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 1:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStress();
      return 0;
    case 2:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();
      return 0;
    case 3:
      if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent();
      return 0;
    case 4:                                   // total backstress, reduced view (stress convention)
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int a = 0; a < ncomp; a++) {
          double tot = 0.0;
          for (int k = 0; k < nBack; k++) tot += alpha[k][vmap[a]];
          v(a) = tot;
        }
      }
      return 0;
    case 5:                                   // plastic strain, reduced view (engineering shear)
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int a = 0; a < ncomp; a++) {
          int full = vmap[a];
          v(a) = (full >= 3) ? 2.0*epsP[full] : epsP[full];
        }
      }
      return 0;
    case 6:                                   // equivalent (accumulated) plastic strain
      if (matInfo.theVector) (*(matInfo.theVector))(0) = ebarP;
      return 0;
    default:
      return -1;
  }
}
