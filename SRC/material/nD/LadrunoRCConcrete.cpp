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

// Ladruno RC plastic-damage nDMaterial — see LadrunoRCConcrete.h. classTag 33015.

#include "LadrunoRCConcrete.h"
#include <elementAPI.h>
#include <Channel.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace ladruno_rc_kernel;

// ===========================================================================
//  parser:  nDMaterial LadrunoRCConcrete $tag $E $nu
//             -Ce {..} -Cs {..} <-Cd {..}>  (compression backbone: strain/nom-stress/damage)
//             -Te {..} -Ts {..} <-Td {..}>  (tension backbone)
//             <-Kc $Kc> <-beta> <-betaFloor $f> <-lublinerReduced>
//             <-rho $rho> <-secant | -numericalTangent>
// ===========================================================================
void* OPS_LadrunoRCConcrete(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "nDMaterial LadrunoRCConcrete error: too few args.\n"
           << "  want: tag E nu -Ce {..} -Cs {..} <-Cd {..}> -Te {..} -Ts {..} <-Td {..}>"
           << " <-Kc Kc> <-beta> <-betaFloor f> <-lublinerReduced> <-rho rho> <-secant|-numericalTangent>\n";
    return 0;
  }

  int numData = 1;
  int tag = 0;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "nDMaterial LadrunoRCConcrete error: invalid tag.\n";
    return 0;
  }
  double E = 0.0, nu = 0.0;
  if (OPS_GetDoubleInput(&numData, &E) < 0)  { opserr << "LadrunoRCConcrete: invalid E.\n";  return 0; }
  if (OPS_GetDoubleInput(&numData, &nu) < 0) { opserr << "LadrunoRCConcrete: invalid nu.\n"; return 0; }

  std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
  double Kc = 2.0 / 3.0, betaFloor = 0.1, rho = 0.0;
  bool betaOn = false, lubRed = false;
  int  tanMode = 0;

  auto readList = [](std::vector<double>& v) {
    v.clear();
    int nd = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      double item;
      int oldn = OPS_GetNumRemainingInputArgs();
      if (OPS_GetDoubleInput(&nd, &item) < 0) {
        int newn = OPS_GetNumRemainingInputArgs();
        if (newn < oldn) OPS_ResetCurrentInputArg(-1);   // un-consume the flag/string
        break;
      }
      v.push_back(item);
    }
  };

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if      (strcmp(opt, "-Ce") == 0) readList(Ce);
    else if (strcmp(opt, "-Cs") == 0) readList(Cs);
    else if (strcmp(opt, "-Cd") == 0) readList(Cd);
    else if (strcmp(opt, "-Te") == 0) readList(Te);
    else if (strcmp(opt, "-Ts") == 0) readList(Ts);
    else if (strcmp(opt, "-Td") == 0) readList(Td);
    else if (strcmp(opt, "-Kc") == 0)        { int nd = 1; OPS_GetDoubleInput(&nd, &Kc); }
    else if (strcmp(opt, "-betaFloor") == 0) { int nd = 1; OPS_GetDoubleInput(&nd, &betaFloor); }
    else if (strcmp(opt, "-rho") == 0)       { int nd = 1; OPS_GetDoubleInput(&nd, &rho); }
    else if (strcmp(opt, "-beta") == 0)              betaOn = true;
    else if (strcmp(opt, "-lublinerReduced") == 0)   lubRed = true;
    else if (strcmp(opt, "-secant") == 0)            tanMode = 1;
    else if (strcmp(opt, "-numericalTangent") == 0)  tanMode = 2;
    // unknown tokens are ignored (forward-compat)
  }

  if (Ce.size() < 2 || Cs.size() != Ce.size()) {
    opserr << "nDMaterial LadrunoRCConcrete error: need -Ce and -Cs of equal length (>=2).\n";
    return 0;
  }
  if (Te.size() < 2 || Ts.size() != Te.size()) {
    opserr << "nDMaterial LadrunoRCConcrete error: need -Te and -Ts of equal length (>=2).\n";
    return 0;
  }

  Params P;
  P.E = E; P.nu = nu; P.Kc = Kc; P.betaFloor = betaFloor;
  P.cdf = 0.0; P.eta = 0.0;
  P.betaOn = betaOn; P.lublinerTCReduced = lubRed; P.tangentMode = tanMode;

  // build backbones via the faithful ASDConcrete3D HardeningLaw c-tor + adjust()
  // (elastic-consistent q). -Cd/-Td are optional -> pad with zeros to match length.
  std::vector<double> Cdf(Ce.size(), 0.0), Tdf(Te.size(), 0.0);
  for (size_t i = 0; i < Cd.size() && i < Cdf.size(); ++i) Cdf[i] = Cd[i];
  for (size_t i = 0; i < Td.size() && i < Tdf.size(); ++i) Tdf[i] = Td[i];
  buildBackbone(P.hc, E, Ce.data(), Cs.data(), Cdf.data(), (int)Ce.size());
  buildBackbone(P.ht, E, Te.data(), Ts.data(), Tdf.data(), (int)Te.size());
  setupParams(P);

  return new LadrunoRCConcrete(tag, P, rho, LadrunoRCConcrete::DIM_3D);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoRCConcrete::LadrunoRCConcrete()
  : NDMaterial(0, ND_TAG_LadrunoRCConcrete),
    rho(0.0), dim(DIM_3D), ncomp(6), condense(false), cEps33(0.0), status(STATUS_OK),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  // safe defaults until recvSelf populates P
  P.E = 1.0; P.nu = 0.0; P.Kc = 2.0/3.0; P.fcft_ratio = 5.0; P.betaFloor = 0.1;
  P.cdf = 0.0; P.eta = 0.0; P.betaOn = false; P.lublinerTCReduced = false; P.tangentMode = 0;
  P.ht.n = 0; P.hc.n = 0;
  this->setupDim();
  this->revertToStart();
}

LadrunoRCConcrete::LadrunoRCConcrete(int tag, const Params& P_, double rho_, int dimMode)
  : NDMaterial(tag, ND_TAG_LadrunoRCConcrete),
    P(P_), rho(rho_), dim(dimMode), ncomp(6), condense(false), cEps33(0.0), status(STATUS_OK),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->setupDim();
  this->revertToStart();
}

LadrunoRCConcrete::~LadrunoRCConcrete() {}

// ===========================================================================
//  dimensional view setup
// ===========================================================================
void LadrunoRCConcrete::setupDim(void)
{
  switch (dim) {
    case DIM_PSTRESS:    ncomp = 3; { int m[3]={0,1,3};     for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_PLATEFIBER: ncomp = 5; { int m[5]={0,1,3,4,5}; for(int a=0;a<5;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_3D:
    default:             ncomp = 6; { int m[6]={0,1,2,3,4,5}; for(int a=0;a<6;a++) vmap[a]=m[a]; } condense=false; break;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
}

void LadrunoRCConcrete::condenseTangent(void)
{
  double d22 = Dtan[2][2];
  if (fabs(d22) < 1.0e-300) return;
  double col[6], row[6];
  for (int i = 0; i < 6; i++) { col[i] = Dtan[i][2]; row[i] = Dtan[2][i]; }
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan[I][J] -= col[I]*row[J]/d22;
}

// ===========================================================================
//  strain interface  (ENGINEERING shear — no tensor conversion)
// ===========================================================================
void LadrunoRCConcrete::integrate(void)
{
  status = returnMap3D(P, strain6, histN, stress6, Dtan, histTr, true, 0);
}

int LadrunoRCConcrete::setTrialStrain(const Vector& e)
{
  double eps33 = strain6[2];
  for (int i = 0; i < 6; i++) strain6[i] = 0.0;
  if (condense) strain6[2] = eps33;
  for (int a = 0; a < ncomp; a++) strain6[vmap[a]] = e(a);   // engineering shear, no factor

  if (!condense) { this->integrate(); return (status == STATUS_NO_CONVERGE) ? -1 : 0; }

  // enforce sigma_33 = 0: guarded Newton on eps_33 (= strain6[2]); dSNPO sec 9.4
  const int maxIt = 25;
  double prevAbs = 1.0e300;
  bool converged = false;
  for (int it = 0; it < maxIt; it++) {
    this->integrate();
    double smag = 0.0; for (int i = 0; i < 6; i++) smag += stress6[i]*stress6[i];
    smag = sqrt(smag);
    double tol22 = 1.0e-10 * (smag > 1.0 ? smag : 1.0);
    double abs33 = fabs(stress6[2]);
    if (abs33 <= tol22) { converged = true; break; }
    double d22 = Dtan[2][2];
    double step;
    if (fabs(d22) < 1.0e-14 * P.E) step = (stress6[2] > 0.0 ? 1.0 : -1.0) * 1.0e-6;  // bisection nudge
    else                          step = stress6[2] / d22;
    double lambda = (abs33 > prevAbs) ? 0.5 : 1.0;   // damp if residual grew
    strain6[2] -= lambda * step;
    prevAbs = abs33;
  }
  if (!converged) {
    status = STATUS_NO_CONVERGE;
    opserr << "WARNING LadrunoRCConcrete: sigma_33 condensation did not converge (tag "
           << this->getTag() << ", |s33|=" << fabs(stress6[2]) << ")\n";
  }
  this->condenseTangent();
  return (status == STATUS_NO_CONVERGE) ? -1 : 0;
}

int LadrunoRCConcrete::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoRCConcrete::setTrialStrainIncr(const Vector& v)
{
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) ne(a) = strain6[vmap[a]] + v(a);
  return this->setTrialStrain(ne);
}
int LadrunoRCConcrete::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoRCConcrete::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = stress6[vmap[a]];   // true stress
  return stressOut;
}
const Vector& LadrunoRCConcrete::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) strainOut(a) = strain6[vmap[a]];   // engineering shear
  return strainOut;
}
const Matrix& LadrunoRCConcrete::getTangent(void)
{
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Dtan[vmap[a]][vmap[b]];
  return tangentOut;
}
const Matrix& LadrunoRCConcrete::getInitialTangent(void)
{
  double C0[6][6];
  elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = C0[vmap[a]][vmap[b]];
  return tangentOut;
}

const char* LadrunoRCConcrete::getType(void) const
{
  switch (dim) {
    case DIM_PSTRESS:    return "PlaneStress";
    case DIM_PLATEFIBER: return "PlateFiber";
    default:             return "ThreeDimensional";
  }
}
int LadrunoRCConcrete::getOrder(void) const { return ncomp; }

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRCConcrete::commitState(void)
{
  histN = histTr;
  cEps33 = strain6[2];
  return 0;
}
int LadrunoRCConcrete::revertToLastCommit(void)
{
  histTr = histN;
  if (condense) strain6[2] = cEps33;
  return 0;
}
int LadrunoRCConcrete::revertToStart(void)
{
  histZero(histN);
  histZero(histTr);
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  cEps33 = 0.0;
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan[a][b] = C0[a][b];
  return 0;
}

// ===========================================================================
//  copies
// ===========================================================================
NDMaterial* LadrunoRCConcrete::getCopy(void)
{
  return new LadrunoRCConcrete(this->getTag(), P, rho, dim);
}
NDMaterial* LadrunoRCConcrete::getCopy(const char* type)
{
  int d = -1;
  if      (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) d = DIM_3D;
  else if (strcmp(type, "PlateFiber") == 0)                                  d = DIM_PLATEFIBER;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) d = DIM_PSTRESS;
  if (d < 0) return NDMaterial::getCopy(type);
  return new LadrunoRCConcrete(this->getTag(), P, rho, d);
}

// ===========================================================================
//  parallel  (serialize params + backbones + committed history)
// ===========================================================================
static const int RC_NSCALAR = 3 /*tag,dim,rho*/ + 10 /*E,nu,Kc,fcft,betaFloor,cdf,eta,betaOn,lubRed,tanMode*/;
static const int RC_BACK = 1 + 3*MAXPTS;   // n + x[]+y[]+q[]
static const int RC_HIST = 6 + 6 + 6;      // stress_eff, strain, (xt,xc,dt_bar,dc_bar,beta,eps1)
static const int RC_DATA = RC_NSCALAR + 2*RC_BACK + RC_HIST + 1 /*cEps33*/;

int LadrunoRCConcrete::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(RC_DATA);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = dim;
  data(c++) = rho;
  data(c++) = P.E; data(c++) = P.nu; data(c++) = P.Kc; data(c++) = P.fcft_ratio;
  data(c++) = P.betaFloor; data(c++) = P.cdf; data(c++) = P.eta;
  data(c++) = P.betaOn ? 1.0 : 0.0;
  data(c++) = P.lublinerTCReduced ? 1.0 : 0.0;
  data(c++) = P.tangentMode;
  data(c++) = P.ht.n;
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.x[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.y[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.q[i];
  data(c++) = P.hc.n;
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.x[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.y[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.q[i];
  for (int i = 0; i < 6; i++) data(c++) = histN.stress_eff[i];
  for (int i = 0; i < 6; i++) data(c++) = histN.strain[i];
  data(c++) = histN.xt; data(c++) = histN.xc;
  data(c++) = histN.dt_bar; data(c++) = histN.dc_bar;
  data(c++) = histN.beta; data(c++) = histN.eps1;
  data(c++) = cEps33;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCConcrete::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoRCConcrete::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(RC_DATA);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCConcrete::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  dim = (int)data(c++);
  rho = data(c++);
  P.E = data(c++); P.nu = data(c++); P.Kc = data(c++); P.fcft_ratio = data(c++);
  P.betaFloor = data(c++); P.cdf = data(c++); P.eta = data(c++);
  P.betaOn = (data(c++) != 0.0);
  P.lublinerTCReduced = (data(c++) != 0.0);
  P.tangentMode = (int)data(c++);
  P.ht.n = (int)data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.x[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.y[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.q[i] = data(c++);
  P.hc.n = (int)data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.x[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.y[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.q[i] = data(c++);
  for (int i = 0; i < 6; i++) histN.stress_eff[i] = data(c++);
  for (int i = 0; i < 6; i++) histN.strain[i] = data(c++);
  histN.xt = data(c++); histN.xc = data(c++);
  histN.dt_bar = data(c++); histN.dc_bar = data(c++);
  histN.beta = data(c++); histN.eps1 = data(c++);
  cEps33 = data(c++);

  this->setupDim();
  histTr = histN;
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  if (condense) strain6[2] = cEps33;
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan[a][b] = C0[a][b];
  return 0;
}

void LadrunoRCConcrete::Print(OPS_Stream& s, int)
{
  s << endln;
  s << "LadrunoRCConcrete (RC plastic-damage + MCFT compression softening)" << endln;
  s << "  tag   : " << this->getTag() << endln;
  s << "  E, nu : " << P.E << ", " << P.nu << endln;
  s << "  Kc    : " << P.Kc << "   fc/ft ratio: " << P.fcft_ratio << endln;
  s << "  beta  : " << (P.betaOn ? "ON" : "off") << "  floor=" << P.betaFloor
    << "  lublinerReduced=" << (P.lublinerTCReduced ? 1 : 0) << endln;
  s << "  view  : " << this->getType() << " (order " << ncomp << ")" << endln;
}

// ===========================================================================
//  recordable responses
// ===========================================================================
Response* LadrunoRCConcrete::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];
  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "damage") == 0 || strcmp(a, "Damage") == 0)
    return new MaterialResponse(this, 4, Vector(2));   // (dt_bar, dc_bar)
  if (strcmp(a, "beta") == 0 || strcmp(a, "compressionSoftening") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (strcmp(a, "eps1") == 0 || strcmp(a, "transverseTensileStrain") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoRCConcrete::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 1: if (matInfo.theVector) *(matInfo.theVector) = this->getStress();  return 0;
    case 2: if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();  return 0;
    case 3: if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent(); return 0;
    case 4: if (matInfo.theVector) { Vector& v = *(matInfo.theVector); v(0) = histTr.dt_bar; v(1) = histTr.dc_bar; } return 0;
    case 5: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.beta;   return 0;
    case 6: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.eps1;   return 0;
    default: return -1;
  }
}
